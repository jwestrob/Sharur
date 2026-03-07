"""
V1 compatibility layer: SemanticState -> flat predicate list.

Converts V2 semantic state back to flat predicate lists for downstream
code that expects V1 format (e.g., search_by_predicates).
"""

from __future__ import annotations

from typing import Optional

from sharur.predicates_v2.model import SemanticAtom, SemanticState


def semantic_state_to_predicates(
    state: SemanticState,
    atoms: Optional[list[SemanticAtom]] = None,
) -> list[str]:
    """Convert V2 state to V1-compatible flat predicate list.

    This ensures V2 can replace V1 without breaking any downstream
    operator or agent that calls search_by_predicates().

    Args:
        state: The V2 semantic state for a protein.
        atoms: Optional raw atoms for reconstructing V1 bookkeeping
            predicates (source flags, confidence, hypothetical).

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
        # V1 emits both "giant" and "massive" for >2000aa proteins
        if state.size_class == "massive":
            predicates.add("giant")

    # Topology -> predicates
    if state.topology.get("tm_count", 0) > 0:
        predicates.add("transmembrane_predicted")
    if state.topology.get("signal_peptide"):
        predicates.add("signal_peptide")

    # Reconstruct V1 bookkeeping predicates from atom metadata
    if atoms:
        predicates.update(_bookkeeping_predicates_from_atoms(atoms))

    return sorted(predicates)


# Source DB names that map to V1 "{source}_annotated" predicates
_SOURCE_TO_ANNOTATED = {
    "pfam": "pfam_annotated",
    "kegg": "kegg_annotated",
    "cazy": "cazy_annotated",
    "vog": "vog_annotated",
    "vogdb": "vog_annotated",
    "hyddb": "hyddb_annotated",
    "defensefinder": "defensefinder_annotated",
    "defensefinder_system": "defensefinder_annotated",
    "txsscan_system": "txsscan_annotated",
    "cant_hyd": "cant_hyd_annotated",
}


def _bookkeeping_predicates_from_atoms(atoms: list[SemanticAtom]) -> set[str]:
    """Reconstruct V1 bookkeeping predicates from V2 atom metadata.

    V1 emits these per-annotation predicates that V2 captures structurally:
    - pfam_annotated, kegg_annotated, etc. (from source_db)
    - confident_hit (evalue < 1e-10)
    - weak_hit (evalue > 1e-5)
    - hypothetical (from atom_id)
    """
    preds: set[str] = set()

    has_confident = False
    has_weak = False

    for atom in atoms:
        # Skip computed/property atoms
        if atom.source_db.startswith("_"):
            continue

        # Source-annotated flags (including from witness atoms)
        annotated_pred = _SOURCE_TO_ANNOTATED.get(atom.source_db)
        if annotated_pred:
            preds.add(annotated_pred)

        # Skip witness atoms for confidence/hypothetical checks
        if atom.atom_id.startswith("_source_witness:"):
            continue

        # Confidence from e-value
        if atom.evidence_evalue is not None:
            if atom.evidence_evalue < 1e-10:
                has_confident = True
            elif atom.evidence_evalue > 1e-5:
                has_weak = True

        # Hypothetical status
        if atom.atom_id == "hypothetical":
            preds.add("hypothetical")

    if has_confident:
        preds.add("confident_hit")
    if has_weak:
        preds.add("weak_hit")

    return preds


__all__ = ["semantic_state_to_predicates"]
