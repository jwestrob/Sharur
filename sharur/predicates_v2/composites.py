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

from sharur.predicates_v2.model import SemanticAtom

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


def explain_composites(
    atoms: list[SemanticAtom],
    composites: Optional[list[CompositeDefinition]] = None,
    topology: Optional[dict[str, Any]] = None,
    only: Optional[list[str]] = None,
) -> dict[str, list[dict[str, Any]]]:
    """Return atom witnesses explaining matched composite predicates.

    Args:
        atoms: All atoms for the protein.
        composites: Composite definitions. Loaded from YAML if None.
        topology: Topology dict from SemanticState.
        only: Optional composite names to explain. When omitted, all matched
            composites are returned.

    Returns:
        Mapping of composite name to the atom dictionaries that satisfied the
        positive parts of its rule. ``none_of`` conditions contribute no
        witnesses because their successful state is absence.
    """
    if composites is None:
        composites = load_composites()

    only_set = set(only) if only is not None else None
    atom_index = _build_atom_index(atoms)
    explanations: dict[str, list[dict[str, Any]]] = {}

    for comp in composites:
        if only_set is not None and comp.name not in only_set:
            continue
        matched, witnesses = _explain_condition(
            comp.requires,
            atom_index,
            topology or {},
        )
        if matched:
            explanations[comp.name] = _dedupe_witnesses(witnesses)

    return dict(sorted(explanations.items()))


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


def _explain_condition(
    condition: dict[str, Any],
    atom_index: dict[str, list[SemanticAtom]],
    topology: dict[str, Any],
) -> tuple[bool, list[dict[str, Any]]]:
    """Evaluate a condition and return positive atom witnesses."""
    witnesses: list[dict[str, Any]] = []

    if "all_of" in condition:
        for child in condition["all_of"]:
            matched, child_witnesses = _explain_condition(child, atom_index, topology)
            if not matched:
                return False, []
            witnesses.extend(child_witnesses)

    if "any_of" in condition:
        any_witnesses: list[dict[str, Any]] | None = None
        for child in condition["any_of"]:
            matched, child_witnesses = _explain_condition(child, atom_index, topology)
            if matched:
                any_witnesses = child_witnesses
                break
        if any_witnesses is None:
            return False, []
        witnesses.extend(any_witnesses)

    if "none_of" in condition:
        for child in condition["none_of"]:
            matched, _child_witnesses = _explain_condition(child, atom_index, topology)
            if matched:
                return False, []

    if any(key in condition for key in ("all_of", "any_of", "none_of")):
        return True, witnesses

    if "has_atom" in condition:
        atom = _matching_atom(condition, atom_index)
        if atom is None:
            return False, []
        return True, [atom.to_dict()]

    if "property_equals" in condition:
        for key, value in condition["property_equals"].items():
            if topology.get(key) != value:
                return False, []
        return True, []

    if "property_gte" in condition:
        for key, value in condition["property_gte"].items():
            actual = topology.get(key)
            if actual is None or actual < value:
                return False, []
        return True, []

    if "property_lte" in condition:
        for key, value in condition["property_lte"].items():
            actual = topology.get(key)
            if actual is None or actual > value:
                return False, []
        return True, []

    return False, []


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
    return _matching_atom(condition, atom_index) is not None


def _matching_atom(
    condition: dict[str, Any],
    atom_index: dict[str, list[SemanticAtom]],
) -> SemanticAtom | None:
    """Return the first atom matching a has_atom condition."""
    candidates = atom_index.get(condition["has_atom"], [])
    relation_filter = condition.get("relation")
    relation_in_filter = condition.get("relation_in")
    source_db_filter = condition.get("source_db")

    for atom in candidates:
        if relation_filter and atom.relation.value != relation_filter:
            continue
        if relation_in_filter and atom.relation.value not in relation_in_filter:
            continue
        if source_db_filter and atom.source_db != source_db_filter:
            continue
        return atom

    return None


def _dedupe_witnesses(witnesses: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Deduplicate witness atoms while preserving rule order."""
    seen: set[tuple[object, ...]] = set()
    deduped = []
    for witness in witnesses:
        key = (
            witness.get("protein_id"),
            witness.get("atom_id"),
            witness.get("relation"),
            witness.get("source_db"),
            witness.get("source_accession"),
        )
        if key in seen:
            continue
        seen.add(key)
        deduped.append(witness)
    return deduped


__all__ = [
    "CompositeDefinition",
    "clear_composites_cache",
    "evaluate_composites",
    "explain_composites",
    "load_composites",
]
