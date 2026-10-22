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
# Prefer wheel-packaged config, with the source-tree path for editable work.
# ---------------------------------------------------------------------------
_PACKAGED_CONFIG_DIR = Path(__file__).parent / "config"
_SOURCE_CONFIG_DIR = Path(__file__).parent.parent.parent / "config" / "predicates_v2"
_CONFIG_DIR = (
    _PACKAGED_CONFIG_DIR
    if _PACKAGED_CONFIG_DIR.is_dir()
    else _SOURCE_CONFIG_DIR
)
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

    problems = composite_problems(defs)
    if problems:
        raise ValueError(f"{target}: invalid composite rules: " + "; ".join(problems))

    if path is None:
        _composites_cache = defs

    return defs


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------
_COMPOUND_KEYS = ("all_of", "any_of", "none_of")
_ATOM_KEYS = ("has_atom", "relation", "relation_in", "source_db")
_PROPERTY_KEYS = ("property_equals", "property_gte", "property_lte")


def composite_problems(
    composites: list[CompositeDefinition],
    vocabulary: Optional[set[str]] = None,
) -> list[str]:
    """Rule errors that would make a composite clause inert or misread.

    Every ``has_atom`` must name a vocabulary predicate (an unknown name never
    matches, so an ``all_of`` clause never passes and a ``none_of`` clause never
    excludes); relations, facets and operator keys must be ones the evaluator
    knows. Returns one message per problem; empty when every rule is valid.
    """
    if vocabulary is None:
        from sharur.predicates.vocabulary import PREDICATE_BY_ID

        vocabulary = set(PREDICATE_BY_ID)
    relations = {r.value for r in ClaimRelation}
    facets = {f.value for f in SemanticFacet}
    problems: list[str] = []

    def check(condition: Any, where: str) -> None:
        if not isinstance(condition, dict) or not condition:
            problems.append(f"{where}: condition must be a non-empty mapping")
            return
        unknown = set(condition) - {*_COMPOUND_KEYS, *_ATOM_KEYS, *_PROPERTY_KEYS}
        if unknown:
            problems.append(f"{where}: unknown operator(s) {sorted(unknown)}")
        # The evaluator reads one kind per mapping: compound keys, else has_atom, else properties.
        kinds = [keys for keys in (_COMPOUND_KEYS, _ATOM_KEYS, _PROPERTY_KEYS)
                 if any(key in condition for key in keys)]
        if len(kinds) > 1:
            problems.append(f"{where}: mixes compound, has_atom and property conditions; "
                            "nest them under all_of")
        for key in _COMPOUND_KEYS:
            if key in condition:
                children = condition[key]
                if not isinstance(children, list) or not children:
                    problems.append(f"{where}.{key}: must be a non-empty list")
                    continue
                for i, child in enumerate(children):
                    check(child, f"{where}.{key}[{i}]")
        if "has_atom" in condition:
            atom = condition["has_atom"]
            if atom not in vocabulary:
                problems.append(f"{where}: has_atom {atom!r} is not a vocabulary predicate")
        elif any(key in condition for key in _ATOM_KEYS):
            problems.append(f"{where}: relation/source_db filters need has_atom")
        if "relation" in condition and condition["relation"] not in relations:
            problems.append(f"{where}: unknown relation {condition['relation']!r}")
        if "relation_in" in condition:
            listed = condition["relation_in"]
            if not isinstance(listed, list) or not listed or set(listed) - relations:
                problems.append(f"{where}: relation_in must list relations from {sorted(relations)}")
        properties = [key for key in _PROPERTY_KEYS if key in condition]
        if len(properties) > 1:
            # The evaluator reads the first property operator only.
            problems.append(f"{where}: one property operator per condition; nest {properties} under all_of")
        for key in properties:
            if not isinstance(condition[key], dict) or not condition[key]:
                problems.append(f"{where}.{key}: must map property names to values")

    for comp in composites:
        if comp.facet not in facets:
            problems.append(f"{comp.name}: unknown facet {comp.facet!r}")
        check(comp.requires, comp.name)
    return problems


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
    "composite_problems",
    "evaluate_composites",
    "explain_composites",
    "load_composites",
]
