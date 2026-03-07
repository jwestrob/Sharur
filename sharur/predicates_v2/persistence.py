"""
Persistence: generate atoms, aggregate, evaluate composites, write to DuckDB.

Orchestrates the full V2 pipeline and persists results to the
semantic_atoms and semantic_state tables.
"""

from __future__ import annotations

import json
from typing import TYPE_CHECKING, Optional

from sharur.predicates.generator import AnnotationRecord, ProteinRecord
from sharur.predicates_v2.aggregator import aggregate_atoms
from sharur.predicates_v2.composites import evaluate_composites, load_composites
from sharur.predicates_v2.generator import AtomGenerator
from sharur.predicates_v2.model import SemanticAtom, SemanticState, create_v2_tables
from sharur.predicates_v2.review_queue import (
    build_review_queue,
    format_review_queue_tsv,
)

if TYPE_CHECKING:
    from sharur.storage.duckdb_store import DuckDBStore


def generate_and_persist_v2(
    store: "DuckDBStore",
    protein_ids: Optional[list[str]] = None,
    output_review_queue: Optional[str] = None,
) -> dict[str, SemanticState]:
    """Generate V2 atoms + state and persist to DuckDB.

    Full pipeline:
    1. Read proteins + annotations from DuckDB
    2. Generate atoms via AtomGenerator
    3. Aggregate atoms -> SemanticState per protein
    4. Evaluate composites
    5. Write to semantic_atoms and semantic_state tables
    6. Optionally generate review queue TSV

    Args:
        store: DuckDB store with proteins and annotations tables.
        protein_ids: Optional subset of proteins (None = all).
        output_review_queue: Optional path to write review queue TSV.

    Returns:
        Dict mapping protein_id -> SemanticState.
    """
    # Ensure V2 tables exist
    create_v2_tables(store)

    # Load composites once
    composites = load_composites()

    # Initialize generator
    gen = AtomGenerator()

    # Get GC stats for outlier detection
    gc_stats = _get_contig_gc_stats(store)

    # Get proteins
    if protein_ids:
        placeholders = ",".join(["?"] * len(protein_ids))
        proteins = store.execute(
            f"SELECT protein_id, sequence_length, gc_content, contig_id "
            f"FROM proteins WHERE protein_id IN ({placeholders})",
            protein_ids,
        )
    else:
        proteins = store.execute(
            "SELECT protein_id, sequence_length, gc_content, contig_id FROM proteins"
        )

    # Get all annotations
    if protein_ids:
        placeholders = ",".join(["?"] * len(protein_ids))
        all_annotations = store.execute(
            f"SELECT protein_id, source, accession, name, description, evalue, score "
            f"FROM annotations WHERE protein_id IN ({placeholders})",
            protein_ids,
        )
    else:
        all_annotations = store.execute(
            "SELECT protein_id, source, accession, name, description, evalue, score "
            "FROM annotations"
        )

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

    # Get protein -> genome mapping for review queue
    protein_genomes: dict[str, str] = {}
    for row in proteins:
        if row[3]:  # contig_id
            protein_genomes[row[0]] = row[3]  # Use contig_id as proxy

    # Try to get bin_id for actual genome mapping
    try:
        if protein_ids:
            placeholders = ",".join(["?"] * len(protein_ids))
            bin_rows = store.execute(
                f"SELECT protein_id, bin_id FROM proteins "
                f"WHERE protein_id IN ({placeholders}) AND bin_id IS NOT NULL",
                protein_ids,
            )
        else:
            bin_rows = store.execute(
                "SELECT protein_id, bin_id FROM proteins WHERE bin_id IS NOT NULL"
            )
        for row in bin_rows:
            protein_genomes[row[0]] = row[1]
    except Exception:
        pass

    # Process each protein
    all_atoms: list[SemanticAtom] = []
    results: dict[str, SemanticState] = {}

    for row in proteins:
        pid = row[0]
        length = row[1]
        gc = row[2]
        contig = row[3]

        gc_mean, gc_std = gc_stats.get(contig, (None, None)) if contig else (None, None)

        protein = ProteinRecord(
            protein_id=pid,
            sequence_length=length,
            gc_content=gc,
            contig_gc_mean=gc_mean,
            contig_gc_std=gc_std,
        )

        annotations = annotations_by_protein.get(pid, [])

        # Generate atoms
        atoms = gen.generate_atoms(protein, annotations)
        all_atoms.extend(atoms)

        # Aggregate
        state = aggregate_atoms(pid, atoms)

        # Evaluate composites
        composite_hits = evaluate_composites(atoms, composites, state.topology)
        state.composite_predicates = composite_hits

        results[pid] = state

    # Persist to DuckDB
    _persist_atoms(store, all_atoms)
    _persist_states(store, results)

    # Generate review queue
    if output_review_queue:
        queue = build_review_queue(
            all_atoms,
            protein_genomes=protein_genomes,
            total_proteins=len(proteins),
        )
        tsv = format_review_queue_tsv(queue)
        with open(output_review_queue, "w") as f:
            f.write(tsv)

    return results


def _persist_atoms(store: "DuckDBStore", atoms: list[SemanticAtom]) -> None:
    """Write atoms to semantic_atoms table."""
    # Clear existing
    store.execute("DELETE FROM semantic_atoms;")

    for atom in atoms:
        store.execute(
            """
            INSERT OR REPLACE INTO semantic_atoms
            (protein_id, atom_id, facet, relation, source_accession,
             source_db, evidence_evalue, evidence_score)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """,
            [
                atom.protein_id,
                atom.atom_id,
                atom.facet.value,
                atom.relation.value,
                atom.source_accession,
                atom.source_db,
                atom.evidence_evalue,
                atom.evidence_score,
            ],
        )


def _persist_states(
    store: "DuckDBStore",
    states: dict[str, SemanticState],
) -> None:
    """Write aggregated states to semantic_state table."""
    # Clear existing
    store.execute("DELETE FROM semantic_state;")

    for pid, state in states.items():
        store.execute(
            """
            INSERT OR REPLACE INTO semantic_state
            (protein_id, activities, roles, architecture, localization,
             topology, size_class, quality_flags, composite_predicates,
             unresolved_count, updated_at)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, CURRENT_TIMESTAMP)
            """,
            [
                pid,
                state.activities,
                state.roles,
                state.architecture,
                state.localization,
                json.dumps(state.topology),
                state.size_class,
                state.quality_flags,
                state.composite_predicates,
                len(state.unresolved_accessions),
            ],
        )


def _get_contig_gc_stats(
    store: "DuckDBStore",
) -> dict[str, tuple[float, float]]:
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


__all__ = ["generate_and_persist_v2"]
