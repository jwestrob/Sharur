"""
Predicate generator for computing predicates from annotations.

This module provides the main interface for generating predicates from
protein annotations. It combines:
1. Annotation-derived predicates (from PFAM, KEGG, CAZy mappings)
2. Computed predicates (from protein properties like size)
3. Hierarchical predicate expansion

Usage:
    from bennu.predicates.generator import PredicateGenerator

    gen = PredicateGenerator()
    predicates = gen.generate_for_protein(protein_data, annotations)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Optional
import re

from bennu.predicates.mappings.pfam_map import get_predicates_for_pfam
from bennu.predicates.mappings.kegg_map import get_predicates_for_kegg
from bennu.predicates.mappings.cazy_map import get_predicates_for_cazy
from bennu.predicates.mappings.vog_map import get_vog_predicates
from bennu.predicates.vocabulary import get_hierarchy, PREDICATE_BY_ID

if TYPE_CHECKING:
    from bennu.storage.duckdb_store import DuckDBStore


@dataclass
class AnnotationRecord:
    """Annotation data for predicate generation."""
    source: str  # 'pfam', 'kegg', 'cazy'
    accession: str
    name: Optional[str] = None
    description: Optional[str] = None
    evalue: Optional[float] = None
    score: Optional[float] = None


@dataclass
class ProteinRecord:
    """Protein data for predicate generation."""
    protein_id: str
    sequence_length: Optional[int] = None
    gc_content: Optional[float] = None
    contig_gc_mean: Optional[float] = None
    contig_gc_std: Optional[float] = None
    sequence: Optional[str] = None  # For topology prediction (requires pyTMHMM)


class PredicateGenerator:
    """
    Generates predicates for proteins based on annotations and properties.

    This is the main class for computing predicates. It:
    1. Maps annotations to predicates using PFAM/KEGG/CAZy mappings
    2. Computes property-based predicates (size, GC, etc.)
    3. Expands predicates through the hierarchy (child -> parent)
    4. Adds annotation status predicates
    5. Optionally computes topology predicates (requires pyTMHMM)
    """

    # Class-level cache for VOGdb reference data
    _vog_reference: Optional[dict[str, tuple[str, str]]] = None

    def __init__(
        self,
        expand_hierarchy: bool = True,
        include_direct_access: bool = True,
        predict_topology: bool = True,
        vog_reference_path: Optional[str] = None,
    ):
        """
        Initialize the predicate generator.

        Args:
            expand_hierarchy: Whether to add parent predicates automatically
            include_direct_access: Whether to add pfam:PF*, kegg:K*, etc. predicates
            predict_topology: Whether to predict TM topology when sequences available
            vog_reference_path: Path to VOGdb annotations TSV (auto-detected if not provided)
        """
        self.expand_hierarchy = expand_hierarchy
        self.include_direct_access = include_direct_access
        self.predict_topology = predict_topology

        # Check if topology prediction is available
        self._topology_available = False
        if self.predict_topology:
            try:
                from bennu.predicates.topology import is_available
                self._topology_available = is_available()
            except ImportError:
                pass

        # Load VOGdb reference for category/description lookup
        self._load_vog_reference(vog_reference_path)

    def _load_vog_reference(self, path: Optional[str] = None) -> None:
        """Load VOGdb annotations reference for functional category lookup."""
        # Use class-level cache
        if PredicateGenerator._vog_reference is not None:
            return

        import os
        from pathlib import Path

        # Auto-detect reference file location
        if path is None:
            # Check common locations
            candidates = [
                Path("data/reference/vogdb/vog.annotations.tsv"),
                Path(__file__).parent.parent.parent.parent / "data/reference/vogdb/vog.annotations.tsv",
                Path(os.path.expanduser("~/.bennu/vogdb/vog.annotations.tsv")),
            ]
            for candidate in candidates:
                if candidate.exists():
                    path = str(candidate)
                    break

        if path is None or not Path(path).exists():
            # No reference file found - VOGdb will use patterns only
            PredicateGenerator._vog_reference = {}
            return

        # Load the reference file
        try:
            vog_ref: dict[str, tuple[str, str]] = {}
            with open(path) as f:
                next(f)  # Skip header
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) >= 5:
                        vog_id = parts[0]
                        category = parts[3]
                        description = parts[4]
                        vog_ref[vog_id] = (category, description)
            PredicateGenerator._vog_reference = vog_ref
        except Exception:
            PredicateGenerator._vog_reference = {}

    def generate_for_protein(
        self,
        protein: ProteinRecord,
        annotations: list[AnnotationRecord],
    ) -> list[str]:
        """
        Generate all predicates for a protein.

        Args:
            protein: Protein record with properties
            annotations: List of annotation records

        Returns:
            Sorted list of unique predicate IDs
        """
        predicates = set()

        # Annotation-derived predicates
        for ann in annotations:
            predicates.update(self._predicates_from_annotation(ann))

        # Property-based predicates
        predicates.update(self._predicates_from_properties(protein))

        # Annotation status predicates
        predicates.update(self._annotation_status_predicates(annotations))

        # Topology predicates (if sequence available and pyTMHMM installed)
        if self._topology_available and protein.sequence:
            predicates.update(self._predicates_from_topology(protein.sequence))

        # Validate HydDB hydrogenase calls against PFAM domains
        # NiFe hydrogenase requires NiFeSe_Hases (PF00374) to distinguish from Complex I
        predicates = self._validate_hydrogenase_calls(predicates)

        # Expand hierarchy (add parent predicates)
        if self.expand_hierarchy:
            predicates = self._expand_hierarchy(predicates)

        return sorted(predicates)

    def _predicates_from_topology(self, sequence: str) -> set[str]:
        """Get predicates from transmembrane topology prediction."""
        predicates = set()

        try:
            from bennu.predicates.topology import predict_topology, get_topology_predicates

            prediction = predict_topology(sequence)
            if prediction:
                predicates.update(get_topology_predicates(prediction))
        except Exception:
            # Silently ignore topology prediction failures
            pass

        return predicates

    def _validate_hydrogenase_calls(self, predicates: set[str]) -> set[str]:
        """
        Validate HydDB hydrogenase classifications against PFAM domains.

        HydDB NiFe hits have ~44% false positive rate from Complex I (NADH dehydrogenase)
        because the [NiFe] binding site is similar. We require the NiFeSe_Hases domain
        (PF00374) to confirm a true NiFe hydrogenase.

        Similarly, FeFe hydrogenases should have Fe_hyd_lg_C (PF02906) or Fe_hyd_SSU (PF02256).
        """
        # Check for HydDB NiFe hit
        has_hyddb_nife = "hyddb:NiFe" in predicates

        if has_hyddb_nife:
            # Check for validating PFAM domain
            has_nifese_hases = "pfam:PF00374" in predicates  # NiFeSe_Hases

            # Check for Complex I domains (indicates this is NADH dehydrogenase, not hydrogenase)
            has_complex1 = (
                "pfam:PF00346" in predicates or  # Complex1_49kDa
                "pfam:PF00329" in predicates     # Complex1_30kDa
            )

            if not has_nifese_hases and has_complex1:
                # This is Complex I (NADH dehydrogenase), not a hydrogenase
                # Don't assign hydrogenase predicates - let PFAM provide correct annotation
                predicates.discard("nife_hydrogenase")
                predicates.discard("hydrogenase")
                predicates.discard("hydrogen_metabolism")

        # Check for HydDB FeFe hit
        has_hyddb_fefe = "hyddb:FeFe" in predicates

        if has_hyddb_fefe:
            # FeFe should have Fe_hyd domains
            has_fe_hyd = (
                "pfam:PF02906" in predicates or  # Fe_hyd_lg_C
                "pfam:PF02256" in predicates     # Fe_hyd_SSU
            )

            if not has_fe_hyd:
                # Uncertain FeFe call - keep the annotation but note it's unvalidated
                # We don't remove the predicate since FeFe false positives are less common
                pass

        return predicates

    def _predicates_from_annotation(self, ann: AnnotationRecord) -> set[str]:
        """Get predicates from a single annotation."""
        predicates = set()
        source = ann.source.lower()

        if source == "pfam":
            # PFAM domain mapping
            preds = get_predicates_for_pfam(
                ann.accession,
                ann.name or "",
                ann.description or "",
            )
            predicates.update(preds)

            # Direct access predicate
            if self.include_direct_access:
                predicates.add(f"pfam:{ann.accession}")

            # Add pfam_annotated
            predicates.add("pfam_annotated")

        elif source == "kegg":
            # KEGG ortholog mapping
            preds = get_predicates_for_kegg(
                ann.accession,
                ann.description or "",
            )
            predicates.update(preds)

            # Direct access predicate
            if self.include_direct_access:
                predicates.add(f"kegg:{ann.accession}")

            # Add kegg_annotated
            predicates.add("kegg_annotated")

        elif source == "cazy":
            # CAZy family mapping
            preds = get_predicates_for_cazy(ann.accession)
            predicates.update(preds)

            # Direct access predicate
            if self.include_direct_access:
                predicates.add(f"cazy:{ann.accession}")

            # Add cazy_annotated
            predicates.add("cazy_annotated")

        elif source in ("vog", "vogdb"):
            # VOGdb orthologous group mapping
            # Look up category and description from reference if not in annotation
            category = ann.name if ann.name and ann.name.startswith("X") else None
            description = ann.description or ""

            # Try reference lookup if we don't have category/description
            if (not category or not description) and PredicateGenerator._vog_reference:
                ref_data = PredicateGenerator._vog_reference.get(ann.accession)
                if ref_data:
                    if not category:
                        category = ref_data[0]
                    if not description:
                        description = ref_data[1]

            preds = get_vog_predicates(
                ann.accession,
                category=category,
                description=description,
            )
            predicates.update(preds)

            # Direct access predicate
            if self.include_direct_access:
                predicates.add(f"vog:{ann.accession}")

            # Add vog_annotated
            predicates.add("vog_annotated")

        elif source == "hyddb":
            # HydDB hydrogenase classification
            # Maps to specific hydrogenase predicates
            hmm_name = ann.accession.lower() if ann.accession else ""
            if hmm_name == "nife":
                predicates.add("nife_hydrogenase")
                predicates.add("hydrogenase")
            elif hmm_name == "fefe":
                predicates.add("fefe_hydrogenase")
                predicates.add("hydrogenase")
            elif hmm_name in ("fe_only", "fe-only", "feonly"):
                predicates.add("fe_only_hydrogenase")
                predicates.add("hydrogenase")

            # Direct access predicate
            if self.include_direct_access:
                predicates.add(f"hyddb:{ann.accession}")

            # Add hyddb_annotated
            predicates.add("hyddb_annotated")

        elif source == "defensefinder":
            # DefenseFinder defense system predictions
            # The accession contains the system/component name
            hmm_name = ann.accession or ""

            # Add direct access predicate
            if self.include_direct_access:
                predicates.add(f"defensefinder:{hmm_name}")

            # Add defense-related semantic predicates based on HMM name patterns
            hmm_lower = hmm_name.lower()

            # CRISPR-Cas systems
            if "cas" in hmm_lower or "crispr" in hmm_lower:
                predicates.add("crispr_associated")
                predicates.add("defense_system")

            # Restriction-modification systems
            if "rm_" in hmm_lower or "restriction" in hmm_lower or "methylase" in hmm_lower:
                predicates.add("restriction_modification")
                predicates.add("defense_system")

            # Toxin-antitoxin systems
            if "ta_" in hmm_lower or "toxin" in hmm_lower or "antitoxin" in hmm_lower:
                predicates.add("toxin_antitoxin")
                predicates.add("defense_system")

            # Abortive infection
            if "abi" in hmm_lower:
                predicates.add("abortive_infection")
                predicates.add("defense_system")

            # CBASS
            if "cbass" in hmm_lower or "cyclic" in hmm_lower:
                predicates.add("cbass")
                predicates.add("defense_system")

            # General defense marker
            predicates.add("defensefinder_annotated")

        # Check for hypothetical annotations
        if self._is_hypothetical(ann):
            predicates.add("hypothetical")

        # Check annotation confidence
        if ann.evalue is not None:
            if ann.evalue < 1e-10:
                predicates.add("confident_hit")
            elif ann.evalue > 1e-5:
                predicates.add("weak_hit")

        return predicates

    def _predicates_from_properties(self, protein: ProteinRecord) -> set[str]:
        """Get predicates from protein properties."""
        predicates = set()

        # Size predicates
        if protein.sequence_length is not None:
            length = protein.sequence_length

            if length < 50:
                predicates.add("tiny")
            elif length < 150:
                predicates.add("small")
            elif length < 400:
                predicates.add("medium")
            elif length < 1000:
                predicates.add("large")
            elif length < 2000:
                predicates.add("giant")
            else:
                predicates.add("giant")
                predicates.add("massive")

        # GC outlier
        if (
            protein.gc_content is not None
            and protein.contig_gc_mean is not None
            and protein.contig_gc_std is not None
            and protein.contig_gc_std > 0
        ):
            deviation = abs(protein.gc_content - protein.contig_gc_mean)
            if deviation > 2 * protein.contig_gc_std:
                predicates.add("gc_outlier")

            if protein.gc_content > 0.6:
                predicates.add("gc_high")
            elif protein.gc_content < 0.4:
                predicates.add("gc_low")

        return predicates

    def _annotation_status_predicates(
        self, annotations: list[AnnotationRecord]
    ) -> set[str]:
        """Get annotation status predicates."""
        predicates = set()

        if not annotations:
            predicates.add("unannotated")
            return predicates

        # Count confident hits
        confident_count = sum(
            1 for a in annotations if a.evalue is not None and a.evalue < 1e-10
        )

        if confident_count >= 3:
            predicates.add("well_annotated")

        # Check source diversity
        sources = set(a.source.lower() for a in annotations)
        if len(sources) == 1:
            predicates.add("single_source")
        elif len(sources) >= 2:
            predicates.add("multi_source")

        # Count domain annotations for multi_domain
        pfam_domains = set(
            a.accession for a in annotations if a.source.lower() == "pfam"
        )
        if len(pfam_domains) >= 3:
            predicates.add("multi_domain")
        elif len(pfam_domains) == 1:
            predicates.add("single_domain")

        return predicates

    def _is_hypothetical(self, ann: AnnotationRecord) -> bool:
        """Check if annotation indicates hypothetical/unknown function."""
        text = f"{ann.name or ''} {ann.description or ''}".lower()
        patterns = [
            r"hypothetical",
            r"uncharacterized",
            r"unknown function",
            r"duf\d+",  # Domain of Unknown Function
            r"predicted protein",
            r"conserved protein",
        ]
        return any(re.search(p, text) for p in patterns)

    def _expand_hierarchy(self, predicates: set[str]) -> set[str]:
        """Expand predicates by adding parent predicates."""
        expanded = set(predicates)

        for pred_id in list(predicates):
            # Skip direct access predicates (pfam:*, kegg:*, etc.)
            if ":" in pred_id:
                continue

            # Get hierarchy and add parents
            hierarchy = get_hierarchy(pred_id)
            expanded.update(hierarchy)

        return expanded


def generate_predicates_for_proteins(
    store: "DuckDBStore",
    protein_ids: Optional[list[str]] = None,
    batch_size: int = 1000,
) -> dict[str, list[str]]:
    """
    Generate predicates for multiple proteins from database.

    Args:
        store: DuckDB store
        protein_ids: Optional subset of proteins (None = all)
        batch_size: Batch size for processing

    Returns:
        Dict mapping protein_id -> list of predicates
    """
    gen = PredicateGenerator()

    # Get GC statistics per contig for GC outlier calculation
    gc_stats = _get_contig_gc_stats(store)

    # Get proteins
    if protein_ids:
        placeholders = ",".join(["?"] * len(protein_ids))
        query = f"""
            SELECT protein_id, sequence_length, gc_content, contig_id
            FROM proteins
            WHERE protein_id IN ({placeholders})
        """
        proteins = store.execute(query, protein_ids)
    else:
        proteins = store.execute("""
            SELECT protein_id, sequence_length, gc_content, contig_id
            FROM proteins
        """)

    # Get all annotations
    annotations_query = """
        SELECT protein_id, source, accession, name, description, evalue, score
        FROM annotations
    """
    if protein_ids:
        placeholders = ",".join(["?"] * len(protein_ids))
        annotations_query = f"""
            SELECT protein_id, source, accession, name, description, evalue, score
            FROM annotations
            WHERE protein_id IN ({placeholders})
        """
        all_annotations = store.execute(annotations_query, protein_ids)
    else:
        all_annotations = store.execute(annotations_query)

    # Group annotations by protein
    annotations_by_protein: dict[str, list[AnnotationRecord]] = {}
    for row in all_annotations:
        pid = row[0]
        if pid not in annotations_by_protein:
            annotations_by_protein[pid] = []
        annotations_by_protein[pid].append(
            AnnotationRecord(
                source=row[1],
                accession=row[2],
                name=row[3],
                description=row[4],
                evalue=row[5],
                score=row[6],
            )
        )

    # Generate predicates
    result = {}
    for row in proteins:
        pid = row[0]
        length = row[1]
        gc = row[2]
        contig = row[3]

        # Get contig GC stats
        gc_mean, gc_std = gc_stats.get(contig, (None, None))

        protein = ProteinRecord(
            protein_id=pid,
            sequence_length=length,
            gc_content=gc,
            contig_gc_mean=gc_mean,
            contig_gc_std=gc_std,
        )

        annotations = annotations_by_protein.get(pid, [])
        predicates = gen.generate_for_protein(protein, annotations)
        result[pid] = predicates

    return result


def _get_contig_gc_stats(store: "DuckDBStore") -> dict[str, tuple[float, float]]:
    """Get GC content mean and std per contig."""
    rows = store.execute("""
        SELECT contig_id, AVG(gc_content), STDDEV(gc_content)
        FROM proteins
        WHERE gc_content IS NOT NULL
        GROUP BY contig_id
    """)
    return {
        row[0]: (row[1], row[2] if row[2] else 0.0)
        for row in rows
    }


def persist_generated_predicates(
    store: "DuckDBStore",
    protein_predicates: Optional[dict[str, list[str]]] = None,
) -> int:
    """
    Persist generated predicates to the protein_predicates table.

    Args:
        store: DuckDB store
        protein_predicates: Pre-computed predicates, or None to generate fresh

    Returns:
        Number of proteins updated
    """
    if protein_predicates is None:
        protein_predicates = generate_predicates_for_proteins(store)

    # Clear existing and insert new
    store.execute("DELETE FROM protein_predicates")

    count = 0
    for protein_id, predicates in protein_predicates.items():
        store.execute(
            """
            INSERT INTO protein_predicates (protein_id, predicates, updated_at)
            VALUES (?, ?, CURRENT_TIMESTAMP)
            """,
            [protein_id, predicates],
        )
        count += 1

    return count


__all__ = [
    "AnnotationRecord",
    "ProteinRecord",
    "PredicateGenerator",
    "generate_predicates_for_proteins",
    "persist_generated_predicates",
]
