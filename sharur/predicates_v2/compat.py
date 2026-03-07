"""
V1 compatibility layer: SemanticState -> flat predicate list.

Converts V2 semantic state back to flat predicate lists for downstream
code that expects V1 format (e.g., search_by_predicates).
"""

from __future__ import annotations

from sharur.predicates_v2.model import SemanticState


def semantic_state_to_predicates(state: SemanticState) -> list[str]:
    """Convert V2 state to V1-compatible flat predicate list.

    This ensures V2 can replace V1 without breaking any downstream
    operator or agent that calls search_by_predicates().

    Args:
        state: The V2 semantic state for a protein.

    Returns:
        Sorted list of unique predicate IDs matching V1 format.
    """
    predicates: set[str] = set()

    predicates.update(state.activities)
    predicates.update(state.roles)
    predicates.update(state.architecture)
    predicates.update(state.localization)
    predicates.update(state.quality_flags)
    predicates.update(state.composite_predicates)

    if state.size_class:
        predicates.add(state.size_class)

    # Topology -> predicates
    if state.topology.get("tm_count", 0) > 0:
        predicates.add("transmembrane_predicted")
    if state.topology.get("signal_peptide"):
        predicates.add("signal_peptide")

    return sorted(predicates)


__all__ = ["semantic_state_to_predicates"]
