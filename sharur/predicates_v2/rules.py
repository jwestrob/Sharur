"""
Rule loader: reads V1 mapping dicts + YAML overrides to produce typed rules.

This module bridges V1's flat predicate mappings with V2's typed atoms by:
1. Importing V1 mapping dicts (PFAM_TO_PREDICATES, KEGG_TO_PREDICATES, etc.)
2. Loading facet assignments from YAML config
3. Loading relation overrides from YAML config
4. Providing lookup functions for facet + relation per (predicate_id, source, accession)
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Optional

import yaml

from sharur.predicates_v2.model import ClaimRelation, SemanticFacet

# ---------------------------------------------------------------------------
# Prefer wheel-packaged config, with the source-tree path for editable work.
# ---------------------------------------------------------------------------
_PACKAGED_CONFIG_DIR = Path(__file__).parent / "config"
_SOURCE_CONFIG_DIR = Path(__file__).parent.parent.parent / "config" / "predicates_v2"
_CONFIG_DIR = (
    _PACKAGED_CONFIG_DIR
    if _PACKAGED_CONFIG_DIR.is_dir()
    else _SOURCE_CONFIG_DIR
)


def _load_yaml(filename: str) -> dict[str, Any]:
    """Load a YAML config file from the predicates_v2 config directory."""
    path = _CONFIG_DIR / filename
    if not path.exists():
        return {}
    with open(path) as f:
        return yaml.safe_load(f) or {}


# ---------------------------------------------------------------------------
# Facet assignment lookup
# ---------------------------------------------------------------------------
_facet_cache: Optional[dict[str, SemanticFacet]] = None

# Category -> facet fallback mapping (matches spec §3.3)
CATEGORY_TO_FACET: dict[str, SemanticFacet] = {
    "enzyme": SemanticFacet.activity,
    "transport": SemanticFacet.role,
    "regulation": SemanticFacet.role,
    "metabolism": SemanticFacet.activity,
    "cazy": SemanticFacet.activity,
    "structure": SemanticFacet.architecture,
    "binding": SemanticFacet.role,
    "envelope": SemanticFacet.localization,
    "mobile": SemanticFacet.role,
    "stress": SemanticFacet.role,
    "size": SemanticFacet.size_class,
    "annotation": SemanticFacet.quality_flag,
    "composition": SemanticFacet.quality_flag,
    "info_processing": SemanticFacet.role,
    "division": SemanticFacet.role,
    "topology": SemanticFacet.topology,
    "viral": SemanticFacet.role,
}


def _load_facet_assignments() -> dict[str, SemanticFacet]:
    """Load facet assignments from YAML, with V1 vocabulary fallback."""
    global _facet_cache
    if _facet_cache is not None:
        return _facet_cache

    raw = _load_yaml("facet_assignments.yaml")
    result: dict[str, SemanticFacet] = {}

    for pred_id, facet_str in raw.items():
        # Strip comments from facet value
        facet_str = str(facet_str).strip()
        try:
            result[pred_id] = SemanticFacet(facet_str)
        except ValueError:
            # Skip malformed entries
            pass

    _facet_cache = result
    return result


def get_facet(predicate_id: str) -> SemanticFacet:
    """Look up the facet for a predicate ID.

    Priority:
    1. YAML facet_assignments.yaml override
    2. V1 vocabulary category -> facet default
    3. Fallback to 'role'
    """
    assignments = _load_facet_assignments()
    if predicate_id in assignments:
        return assignments[predicate_id]

    # Try V1 vocabulary lookup
    try:
        from sharur.predicates.vocabulary import PREDICATE_BY_ID
        pred = PREDICATE_BY_ID.get(predicate_id)
        if pred:
            return CATEGORY_TO_FACET.get(pred.category, SemanticFacet.role)
    except ImportError:
        pass

    return SemanticFacet.role


# ---------------------------------------------------------------------------
# Relation override lookup
# ---------------------------------------------------------------------------
_relation_cache: Optional[dict[str, Any]] = None


def _load_relation_overrides() -> dict[str, Any]:
    """Load relation overrides from YAML."""
    global _relation_cache
    if _relation_cache is not None:
        return _relation_cache

    _relation_cache = _load_yaml("relation_overrides.yaml")
    return _relation_cache


# Source name normalization
_SOURCE_ALIASES: dict[str, str] = {
    "vog": "vogdb",
    "vogdb": "vogdb",
    "pfam": "pfam",
    "kegg": "kegg",
    "kofam": "kegg",
    "cazy": "cazy",
    "hyddb": "hyddb",
    "hyddb_subgroup": "hyddb_subgroup",
    "defensefinder": "defensefinder",
    "defensefinder_system": "defensefinder_system",
    "txsscan": "txsscan",
    "txsscan_system": "txsscan_system",
    "cant_hyd": "cant_hyd",
    "_property": "_property",
    "_computed": "_computed",
}


def get_relation(
    source_db: str,
    accession: Optional[str] = None,
) -> ClaimRelation:
    """Look up the claim relation for a source + accession.

    Priority:
    1. Accession-level override in YAML
    2. Source-level default in YAML
    3. Hardcoded default (supports)
    """
    overrides = _load_relation_overrides()

    # Check accession-level override first
    if accession:
        acc_overrides = overrides.get("accession_overrides", {})
        if accession in acc_overrides:
            return ClaimRelation(acc_overrides[accession])

    # Check source-level default
    src = _SOURCE_ALIASES.get(source_db.lower(), source_db.lower())
    src_defaults = overrides.get("source_defaults", {})

    # Try normalized source, then raw
    for key in (src, source_db.lower()):
        if key in src_defaults:
            return ClaimRelation(src_defaults[key])

    # Hardcoded fallback
    return ClaimRelation.supports


# ---------------------------------------------------------------------------
# Cache invalidation (for testing)
# ---------------------------------------------------------------------------

def clear_caches() -> None:
    """Clear all cached config data. Used in tests."""
    global _facet_cache, _relation_cache
    _facet_cache = None
    _relation_cache = None


__all__ = [
    "CATEGORY_TO_FACET",
    "clear_caches",
    "get_facet",
    "get_relation",
]
