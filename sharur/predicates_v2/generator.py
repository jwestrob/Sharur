"""
Atom generator: annotation -> SemanticAtom list.

Replaces V1's PredicateGenerator by producing typed SemanticAtom objects
instead of flat predicate strings. Uses V1 mapping dicts under the hood
but wraps each emitted predicate with facet + relation metadata.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Optional

from sharur.predicates.generator import (
    AnnotationRecord,
    PredicateGenerator,
    ProteinRecord,
)
from sharur.predicates_v2.model import ClaimRelation, SemanticAtom, SemanticFacet
from sharur.predicates_v2.rules import get_facet, get_relation

if TYPE_CHECKING:
    pass


class AtomGenerator:
    """Generate SemanticAtom objects for proteins from their annotations.

    This wraps V1's PredicateGenerator and converts its flat predicate output
    into typed SemanticAtom objects with facet + relation metadata.
    """

    def __init__(
        self,
        expand_hierarchy: bool = True,
        predict_topology: bool = True,
        vog_reference_path: Optional[str] = None,
    ):
        """Initialize the atom generator.

        Args:
            expand_hierarchy: Whether to expand parent predicates.
            predict_topology: Whether to predict TM topology (requires pyTMHMM).
            vog_reference_path: Path to VOGdb annotations TSV.
        """
        # We use V1's PredicateGenerator internally with direct_access=True
        # so we can track which accession produced each predicate.
        self._v1_gen = PredicateGenerator(
            expand_hierarchy=expand_hierarchy,
            include_direct_access=True,
            predict_topology=predict_topology,
            vog_reference_path=vog_reference_path,
        )

    def generate_atoms(
        self,
        protein: ProteinRecord,
        annotations: list[AnnotationRecord],
    ) -> list[SemanticAtom]:
        """Generate all semantic atoms for a protein.

        This runs the V1 predicate pipeline and wraps each result
        with facet + relation metadata. Additionally, annotations
        that produce NO predicates are emitted as unresolved atoms.

        Args:
            protein: Protein record with properties.
            annotations: List of annotation records.

        Returns:
            List of SemanticAtom objects.
        """
        atoms: list[SemanticAtom] = []

        # Track which accessions produced predicates
        accessions_that_mapped: set[str] = set()

        # --- Per-annotation atoms ---
        for ann in annotations:
            ann_atoms = self._atoms_from_annotation(protein.protein_id, ann)
            if ann_atoms:
                accessions_that_mapped.add(ann.accession)
            atoms.extend(ann_atoms)

        # --- Property-based atoms (size, GC) ---
        atoms.extend(self._atoms_from_properties(protein))

        # --- Annotation status atoms ---
        atoms.extend(self._annotation_status_atoms(protein.protein_id, annotations))

        # --- Topology atoms ---
        if self._v1_gen._topology_available and protein.sequence:
            atoms.extend(self._topology_atoms(protein))

        # --- Unresolved accessions ---
        for ann in annotations:
            if ann.accession not in accessions_that_mapped:
                # This annotation had no mapping rule
                atoms.append(SemanticAtom(
                    protein_id=protein.protein_id,
                    atom_id=f"unresolved:{ann.accession}",
                    facet=SemanticFacet.quality_flag,
                    relation=ClaimRelation.unresolved,
                    source_accession=ann.accession,
                    source_db=ann.source.lower(),
                    evidence_evalue=ann.evalue,
                    evidence_score=ann.score,
                ))

        # --- Source witness atoms ---
        # Record which annotation sources were seen, even if all semantic
        # atoms from that source are later stripped by validation.
        # This lets the compat layer reconstruct *_annotated flags.
        seen_sources = {ann.source.lower() for ann in annotations}
        for source_db in seen_sources:
            atoms.append(SemanticAtom(
                protein_id=protein.protein_id,
                atom_id=f"_source_witness:{source_db}",
                facet=SemanticFacet.quality_flag,
                relation=ClaimRelation.implies,
                source_accession="_witness",
                source_db=source_db,
            ))

        # --- Post-processing: hydrogenase validation ---
        # Replicate V1's _validate_hydrogenase_calls: NiFe hydrogenase
        # hits from HydDB are false positives when Complex I domains
        # are present without PF00374 (NiFeSe_Hases).
        atoms = self._validate_hydrogenase_atoms(atoms)

        return atoms

    def _atoms_from_annotation(
        self,
        protein_id: str,
        ann: AnnotationRecord,
    ) -> list[SemanticAtom]:
        """Generate atoms from a single annotation.

        Uses V1's _predicates_from_annotation internally, then wraps
        each predicate with the appropriate facet and relation.
        """
        # Get V1 predicates for this annotation
        v1_preds = self._v1_gen._predicates_from_annotation(ann)

        # Also run hierarchy expansion on them
        if self._v1_gen.expand_hierarchy:
            v1_preds = self._v1_gen._expand_hierarchy(v1_preds)

        atoms = []
        source_db = ann.source.lower()

        # Get the base relation for this source + accession
        base_relation = get_relation(source_db, ann.accession)

        # Meta-predicates that V1 emits per-annotation but are not
        # semantic mappings — they indicate annotation quality, not function.
        # (hypothetical is kept — it flows through as a quality_flag atom)
        _META_PREDICATES = {
            "confident_hit", "weak_hit",
        }

        for pred_id in v1_preds:
            # Skip direct-access predicates (pfam:PF00005, kegg:K00001, etc.)
            # These are not semantic atoms — they're index keys
            if ":" in pred_id and pred_id.split(":")[0] in (
                "pfam", "kegg", "kofam", "cazy", "vog", "hyddb",
                "hyddb_subgroup", "defensefinder", "txsscan",
            ):
                continue

            # Skip source-annotated flags that are quality indicators
            # (pfam_annotated, kegg_annotated, etc.) — handled separately
            if pred_id.endswith("_annotated"):
                continue

            # Skip meta-predicates (confidence, hypothetical status)
            # These are generated per-annotation by V1 but belong in
            # annotation status atoms, not semantic mapping atoms.
            if pred_id in _META_PREDICATES:
                continue

            facet = get_facet(pred_id)
            relation = base_relation

            atoms.append(SemanticAtom(
                protein_id=protein_id,
                atom_id=pred_id,
                facet=facet,
                relation=relation,
                source_accession=ann.accession,
                source_db=source_db,
                evidence_evalue=ann.evalue,
                evidence_score=ann.score,
            ))

        return atoms

    def _atoms_from_properties(
        self, protein: ProteinRecord,
    ) -> list[SemanticAtom]:
        """Generate atoms from protein properties (size, GC)."""
        v1_preds = self._v1_gen._predicates_from_properties(protein)
        atoms = []

        for pred_id in v1_preds:
            facet = get_facet(pred_id)
            atoms.append(SemanticAtom(
                protein_id=protein.protein_id,
                atom_id=pred_id,
                facet=facet,
                relation=ClaimRelation.implies,  # Properties are deterministic
                source_accession="_computed",
                source_db="_property",
            ))

        return atoms

    def _annotation_status_atoms(
        self,
        protein_id: str,
        annotations: list[AnnotationRecord],
    ) -> list[SemanticAtom]:
        """Generate annotation status atoms (unannotated, well_annotated, etc.)."""
        v1_preds = self._v1_gen._annotation_status_predicates(annotations)
        atoms = []

        for pred_id in v1_preds:
            facet = get_facet(pred_id)
            atoms.append(SemanticAtom(
                protein_id=protein_id,
                atom_id=pred_id,
                facet=facet,
                relation=ClaimRelation.implies,  # Status is deterministic
                source_accession="_computed",
                source_db="_computed",
            ))

        return atoms

    def _validate_hydrogenase_atoms(
        self, atoms: list[SemanticAtom],
    ) -> list[SemanticAtom]:
        """Validate hydrogenase atoms against Complex I false positives.

        When HydDB NiFe hit exists but PF00374 (NiFeSe_Hases) is absent
        and Complex I domains are present, emit `excludes` atoms. The
        aggregator will surface implies+excludes as a conflict for review.
        Use the composite `nife_hydrogenase_validated` to find true positives.
        """
        atom_ids = {a.atom_id for a in atoms}
        source_accessions = {a.source_accession for a in atoms}

        has_hyddb_nife = "nife_hydrogenase" in atom_ids
        if not has_hyddb_nife:
            return atoms

        has_nifese_hases = "PF00374" in source_accessions
        has_complex1 = "PF00346" in source_accessions or "PF00329" in source_accessions

        if not has_nifese_hases and has_complex1:
            # Emit excludes atoms — aggregator will surface as conflict
            pid = atoms[0].protein_id
            complex1_acc = "PF00346" if "PF00346" in source_accessions else "PF00329"
            for excluded_id in ("nife_hydrogenase", "hydrogenase", "hydrogen_metabolism"):
                atoms.append(SemanticAtom(
                    protein_id=pid,
                    atom_id=excluded_id,
                    facet=get_facet(excluded_id),
                    relation=ClaimRelation.excludes,
                    source_accession=complex1_acc,
                    source_db="_validation",
                ))

        return atoms

    def _topology_atoms(
        self, protein: ProteinRecord,
    ) -> list[SemanticAtom]:
        """Generate topology atoms from TM prediction."""
        if not protein.sequence:
            return []

        v1_preds = self._v1_gen._predicates_from_topology(protein.sequence)
        atoms = []

        for pred_id in v1_preds:
            facet = get_facet(pred_id)
            atoms.append(SemanticAtom(
                protein_id=protein.protein_id,
                atom_id=pred_id,
                facet=facet,
                relation=ClaimRelation.implies,  # Topology is deterministic
                source_accession="_computed",
                source_db="_topology",
            ))

        return atoms


__all__ = [
    "AtomGenerator",
    "AnnotationRecord",
    "ProteinRecord",
]
