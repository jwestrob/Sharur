"""
V1 compatibility layer: SemanticState -> flat predicate list.

Converts V2 semantic state back to flat predicate lists for downstream
code that expects V1 format (e.g., search_by_predicates).
"""

from __future__ import annotations

from typing import TYPE_CHECKING


if TYPE_CHECKING:
    from sharur.predicates_v2.model import SemanticAtom, SemanticState


def semantic_state_to_predicates(
    state: SemanticState,
    atoms: list[SemanticAtom] | None = None,
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
        predicates.update(_direct_access_predicates_from_atoms(atoms))

    return sorted(predicates)


# Source DB names that map to V1 "{source}_annotated" predicates
_SOURCE_TO_ANNOTATED = {
    "pfam": "pfam_annotated",
    "kegg": "kegg_annotated",
    "kofam": "kegg_annotated",
    "cazy": "cazy_annotated",
    "vog": "vog_annotated",
    "vogdb": "vog_annotated",
    "hyddb": "hyddb_annotated",
    "defensefinder": "defensefinder_annotated",
    "defensefinder_system": "defensefinder_annotated",
    "txsscan_system": "txsscan_annotated",
    "cant_hyd": "cant_hyd_annotated",
}


# Source DB names that map to V1 "{source}:{accession}" direct-access predicates.
SOURCE_TO_DIRECT_PREFIX = {
    "pfam": "pfam",
    "kegg": "kegg",
    "kofam": "kegg",
    "cazy": "cazy",
    "vog": "vog",
    "vogdb": "vog",
    "hyddb": "hyddb",
    "hyddb_subgroup": "hyddb_subgroup",
    "defensefinder": "defensefinder",
    "defensefinder_system": "defensefinder",
    "txsscan": "txsscan",
    "txsscan_system": "txsscan",
    "cant_hyd": "cant_hyd",
}


def direct_access_predicate_from_atom(atom: SemanticAtom) -> str | None:
    """Return a V1-style direct accession predicate for an atom, if valid."""
    prefix = SOURCE_TO_DIRECT_PREFIX.get(atom.source_db)
    if not prefix:
        return None
    if atom.atom_id.startswith("_source_witness:"):
        return None
    if not atom.source_accession or atom.source_accession.startswith("_"):
        return None
    return f"{prefix}:{atom.source_accession}"


def direct_access_prefix_case_sql(source_column: str = "source_db") -> str:
    """Return a DuckDB CASE expression for source -> direct-access prefix."""
    cases = " ".join(
        f"WHEN '{source}' THEN '{prefix}'"
        for source, prefix in SOURCE_TO_DIRECT_PREFIX.items()
    )
    return f"CASE {source_column} {cases} ELSE NULL END"


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


def _direct_access_predicates_from_atoms(atoms: list[SemanticAtom]) -> set[str]:
    """Reconstruct V1 direct-access predicates from V2 evidence fields."""
    preds: set[str] = set()

    for atom in atoms:
        if direct_access := direct_access_predicate_from_atom(atom):
            preds.add(direct_access)

    return preds


__all__ = [
    "SOURCE_TO_DIRECT_PREFIX",
    "direct_access_prefix_case_sql",
    "direct_access_predicate_from_atom",
    "semantic_state_to_predicates",
]
