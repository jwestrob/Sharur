"""
Composite predicate evaluator.

Loads composite definitions from YAML and evaluates them against
a protein's atom set to derive higher-order predicates.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

import yaml

from sharur.predicates_v2.model import ClaimRelation, SemanticAtom, SemanticFacet

# ---------------------------------------------------------------------------
# Config path
# ---------------------------------------------------------------------------
_CONFIG_DIR = Path(__file__).parent.parent.parent / "config" / "predicates_v2"
_COMPOSITES_FILE = _CONFIG_DIR / "composites.yaml"


@dataclass
class CompositeDefinition:
    """A composite predicate definition loaded from YAML."""

    name: str
    description: str
    facet: str
    requires: dict[str, Any]


# ---------------------------------------------------------------------------
# YAML loader
# ---------------------------------------------------------------------------
_composites_cache: Optional[list[CompositeDefinition]] = None


def load_composites(path: Optional[Path] = None) -> list[CompositeDefinition]:
    """Load composite definitions from YAML.

    Args:
        path: Override path for testing. Uses default if None.

    Returns:
        List of CompositeDefinition objects.
    """
    global _composites_cache
    if _composites_cache is not None and path is None:
        return _composites_cache

    target = path or _COMPOSITES_FILE
    if not target.exists():
        return []

    with open(target) as f:
        raw = yaml.safe_load(f) or {}

    defs = []
    for name, config in raw.items():
        if not isinstance(config, dict):
            continue
        defs.append(CompositeDefinition(
            name=name,
            description=config.get("description", ""),
            facet=config.get("facet", "quality_flag"),
            requires=config.get("requires", {}),
        ))

    if path is None:
        _composites_cache = defs

    return defs


def clear_composites_cache() -> None:
    """Clear cached composite definitions. Used in tests."""
    global _composites_cache
    _composites_cache = None


# ---------------------------------------------------------------------------
# Evaluator
# ---------------------------------------------------------------------------

def evaluate_composites(
    atoms: list[SemanticAtom],
    composites: Optional[list[CompositeDefinition]] = None,
    topology: Optional[dict[str, Any]] = None,
) -> list[str]:
    """Evaluate composite predicates against a protein's atom set.

    Args:
        atoms: All atoms for the protein.
        composites: Composite definitions. Loaded from YAML if None.
        topology: Topology dict from SemanticState (for property_* conditions).

    Returns:
        List of composite predicate names that matched.
    """
    if composites is None:
        composites = load_composites()

    # Build lookup structures for efficient evaluation
    atom_index = _build_atom_index(atoms)

    matched = []
    for comp in composites:
        if _evaluate_condition(comp.requires, atom_index, topology or {}):
            matched.append(comp.name)

    return sorted(matched)


def _build_atom_index(atoms: list[SemanticAtom]) -> dict[str, list[SemanticAtom]]:
    """Build index: atom_id -> list of atoms with that ID."""
    index: dict[str, list[SemanticAtom]] = {}
    for atom in atoms:
        if atom.atom_id not in index:
            index[atom.atom_id] = []
        index[atom.atom_id].append(atom)
    return index


def _evaluate_condition(
    condition: dict[str, Any],
    atom_index: dict[str, list[SemanticAtom]],
    topology: dict[str, Any],
) -> bool:
    """Evaluate a single condition or compound condition.

    Handles: all_of, any_of, none_of, has_atom (with optional filters).

    When a condition dict contains multiple compound keys (e.g., all_of + none_of),
    ALL of them must be satisfied.
    """
    has_compound = False
    result = True

    if "all_of" in condition:
        has_compound = True
        if not all(
            _evaluate_condition(c, atom_index, topology)
            for c in condition["all_of"]
        ):
            result = False

    if "any_of" in condition:
        has_compound = True
        if not any(
            _evaluate_condition(c, atom_index, topology)
            for c in condition["any_of"]
        ):
            result = False

    if "none_of" in condition:
        has_compound = True
        if any(
            _evaluate_condition(c, atom_index, topology)
            for c in condition["none_of"]
        ):
            result = False

    if has_compound:
        return result

    # Leaf condition: has_atom with optional filters
    if "has_atom" in condition:
        return _evaluate_has_atom(condition, atom_index)

    # Property conditions (operate on topology dict)
    if "property_equals" in condition:
        for key, value in condition["property_equals"].items():
            if topology.get(key) != value:
                return False
        return True

    if "property_gte" in condition:
        for key, value in condition["property_gte"].items():
            actual = topology.get(key)
            if actual is None or actual < value:
                return False
        return True

    if "property_lte" in condition:
        for key, value in condition["property_lte"].items():
            actual = topology.get(key)
            if actual is None or actual > value:
                return False
        return True

    # Unknown condition type
    return False


def _evaluate_has_atom(
    condition: dict[str, Any],
    atom_index: dict[str, list[SemanticAtom]],
) -> bool:
    """Evaluate a has_atom condition with optional relation/source filters.

    Condition keys:
        has_atom: atom_id to check
        relation: exact relation match
        relation_in: list of acceptable relations
        source_db: required source database
    """
    atom_id = condition["has_atom"]
    candidates = atom_index.get(atom_id, [])

    if not candidates:
        return False

    # Apply filters
    relation_filter = condition.get("relation")
    relation_in_filter = condition.get("relation_in")
    source_db_filter = condition.get("source_db")

    for atom in candidates:
        # Check relation filter
        if relation_filter and atom.relation.value != relation_filter:
            continue
        if relation_in_filter and atom.relation.value not in relation_in_filter:
            continue
        # Check source_db filter
        if source_db_filter and atom.source_db != source_db_filter:
            continue
        # Atom matches all filters
        return True

    return False


__all__ = [
    "CompositeDefinition",
    "clear_composites_cache",
    "evaluate_composites",
    "load_composites",
]
