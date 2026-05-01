"""
Persistence: generate atoms, aggregate, evaluate composites, write to DuckDB.

Orchestrates the full V2 pipeline and persists results to the semantic_atoms,
semantic_state, and semantic_terms tables.
"""

from __future__ import annotations

import json
from typing import TYPE_CHECKING

import pandas as pd

from sharur.predicates.generator import AnnotationRecord, ProteinRecord
from sharur.predicates_v2.aggregator import aggregate_atoms
from sharur.predicates_v2.compat import (
    direct_access_prefix_case_sql,
    direct_access_predicate_from_atom,
    semantic_state_to_predicates,
)
from sharur.predicates_v2.composites import evaluate_composites, load_composites
from sharur.predicates_v2.generator import AtomGenerator
from sharur.predicates_v2.model import (
    ClaimRelation,
    SemanticAtom,
    SemanticFacet,
    SemanticState,
    create_v2_tables,
)
from sharur.predicates_v2.review_queue import (
    build_review_queue,
    format_review_queue_tsv,
)
from sharur.predicates_v2.validated_systems import (
    fetch_validated_system_annotations,
    materialize_system_proteins,
)


if TYPE_CHECKING:
    from collections.abc import Iterable

    from sharur.storage.duckdb_store import DuckDBStore


def generate_and_persist_v2(
    store: DuckDBStore,
    protein_ids: list[str] | None = None,
    output_review_queue: str | None = None,
    chunk_size: int = 25_000,
    update_legacy_predicates: bool = False,
    return_states: bool = True,
    predict_topology: bool = False,
) -> dict[str, SemanticState]:
    """Generate V2 atoms + state and persist to DuckDB.

    Full pipeline:
    1. Read proteins + annotations from DuckDB in chunks
    2. Generate atoms via AtomGenerator
    3. Aggregate atoms -> SemanticState per protein
    4. Evaluate composites
    5. Write to semantic_atoms and semantic_state tables
    6. Optionally generate review queue TSV

    Args:
        store: DuckDB store with proteins and annotations tables.
        protein_ids: Optional subset of proteins (None = all).
        output_review_queue: Optional path to write review queue TSV.
        chunk_size: Number of proteins to process per database batch.
        update_legacy_predicates: Also materialize V2-compatible flat
            predicates into protein_predicates for legacy operators.
        return_states: Return generated SemanticState objects in memory.
            Disable for full-dataset runs to avoid retaining millions of
            states after persistence.
        predict_topology: Run sequence-based topology prediction. Defaults
            off for full-dataset V2 generation because topology prediction is
            much slower than annotation-derived atom generation.

    Returns:
        Dict mapping protein_id -> SemanticState when return_states is True,
        otherwise an empty dict.
    """
    # Ensure V2 tables exist
    create_v2_tables(store)

    if chunk_size <= 0:
        raise ValueError("chunk_size must be positive")

    if protein_ids is not None:
        protein_ids = list(dict.fromkeys(protein_ids))
        if not protein_ids:
            return {}
        _backfill_semantic_terms_if_empty(store)

    # Load composites once
    composites = load_composites()

    # Initialize generator
    gen = AtomGenerator(predict_topology=predict_topology)

    materialize_system_proteins(store)
    validated_system_annotations = fetch_validated_system_annotations(
        store,
        protein_ids=set(protein_ids) if protein_ids is not None else None,
    )

    # Get GC stats for outlier detection
    gc_stats = _get_contig_gc_stats(store)

    # Full refresh replaces the V2 tables. Subset refresh must only replace
    # rows for the requested proteins; otherwise a targeted rerun destroys
    # previously generated state for the rest of the dataset.
    if protein_ids is None:
        _clear_v2_tables(store, update_legacy_predicates=update_legacy_predicates)
    else:
        _delete_for_proteins(
            store,
            protein_ids,
            update_legacy_predicates=update_legacy_predicates,
        )

    results: dict[str, SemanticState] = {}
    unresolved_atoms: list[SemanticAtom] = []
    protein_genomes: dict[str, str] = {}
    total_processed = 0

    for protein_rows in _iter_protein_chunks(store, protein_ids, chunk_size):
        if not protein_rows:
            continue

        chunk_protein_ids = [row[0] for row in protein_rows]
        annotations_by_protein = _fetch_annotations_by_protein(store, chunk_protein_ids)
        context_atoms_by_protein = _fetch_context_atoms_by_protein(
            store,
            chunk_protein_ids,
        )

        chunk_atoms: list[SemanticAtom] = []
        chunk_states: dict[str, SemanticState] = {}
        legacy_rows: list[tuple[str, list[str]]] = []

        for row in protein_rows:
            pid = row[0]
            length = row[1]
            gc = row[2]
            contig = row[3]
            genome = row[4] or contig
            if genome:
                protein_genomes[pid] = genome

            gc_mean, gc_std = (
                gc_stats.get(contig, (None, None)) if contig else (None, None)
            )

            protein = ProteinRecord(
                protein_id=pid,
                sequence_length=length,
                gc_content=gc,
                contig_gc_mean=gc_mean,
                contig_gc_std=gc_std,
            )

            annotations = annotations_by_protein.get(pid, [])
            if pid in validated_system_annotations:
                annotations = annotations + validated_system_annotations[pid]
            atoms = gen.generate_atoms(protein, annotations)
            atoms.extend(context_atoms_by_protein.get(pid, []))
            chunk_atoms.extend(atoms)

            state = aggregate_atoms(pid, atoms)
            state.composite_predicates = evaluate_composites(
                atoms, composites, state.topology
            )
            chunk_states[pid] = state

            if update_legacy_predicates:
                legacy_rows.append((pid, semantic_state_to_predicates(state, atoms=atoms)))

            if return_states:
                results[pid] = state

        if output_review_queue:
            unresolved_atoms.extend(
                atom for atom in chunk_atoms if atom.relation.value == "unresolved"
            )

        _persist_chunk(
            store,
            chunk_atoms,
            chunk_states,
            legacy_rows=legacy_rows,
            update_legacy_predicates=update_legacy_predicates,
        )

        total_processed += len(protein_rows)

    # Generate review queue
    if output_review_queue:
        queue = build_review_queue(
            unresolved_atoms,
            protein_genomes=protein_genomes,
            total_proteins=total_processed,
        )
        tsv = format_review_queue_tsv(queue)
        with open(output_review_queue, "w") as f:
            f.write(tsv)

    return results


def materialize_legacy_predicates_from_v2(
    store: DuckDBStore,
    protein_ids: list[str] | None = None,
    chunk_size: int = 100_000,
) -> int:
    """Materialize protein_predicates from persisted V2 state.

    This is the repair path for stale compatibility caches. It does not
    recompute atoms; it reads semantic_state and semantic_atoms and converts
    them through the V2 compatibility layer.

    Args:
        store: DuckDB store with populated semantic_state and semantic_atoms.
        protein_ids: Optional subset to refresh.
        chunk_size: Number of semantic_state rows per batch.

    Returns:
        Number of protein_predicates rows written.
    """
    if chunk_size <= 0:
        raise ValueError("chunk_size must be positive")

    if protein_ids is not None:
        protein_ids = list(dict.fromkeys(protein_ids))
        if not protein_ids:
            return 0
        _delete_legacy_for_proteins(store, protein_ids)
    else:
        store.execute("DELETE FROM protein_predicates;")

    count = 0
    for state_rows in _iter_semantic_state_chunks(store, protein_ids, chunk_size):
        if not state_rows:
            continue

        chunk_protein_ids = [row[0] for row in state_rows]
        atoms_by_protein = _fetch_atoms_by_protein(store, chunk_protein_ids)
        legacy_rows = []

        for row in state_rows:
            state = _semantic_state_from_row(row)
            atoms = atoms_by_protein.get(state.protein_id, [])
            legacy_rows.append((
                state.protein_id,
                semantic_state_to_predicates(state, atoms=atoms),
            ))

        _persist_legacy_predicates(store, legacy_rows)
        count += len(legacy_rows)

    return count


def materialize_semantic_terms_from_v2(
    store: DuckDBStore,
    protein_ids: list[str] | None = None,
    chunk_size: int = 100_000,
) -> int:
    """Materialize the unified V2 search table from persisted atoms/state.

    This is the repair/backfill path for databases generated before
    ``semantic_terms`` existed. It does not recompute atoms; it reads
    semantic_atoms and semantic_state and writes the unified terms used by
    search_by_atoms().

    Args:
        store: DuckDB store with populated semantic_state and semantic_atoms.
        protein_ids: Optional subset to refresh.
        chunk_size: Number of semantic_state rows per batch.

    Returns:
        Number of semantic_terms rows written.
    """
    create_v2_tables(store)

    if chunk_size <= 0:
        raise ValueError("chunk_size must be positive")

    if protein_ids is not None:
        protein_ids = list(dict.fromkeys(protein_ids))
        if not protein_ids:
            return 0
        _delete_semantic_terms_for_proteins(store, protein_ids)
    else:
        store.execute("DELETE FROM semantic_terms;")

    if protein_ids is None:
        _insert_semantic_terms_from_v2_sql(store)
        return store.execute("SELECT COUNT(*) FROM semantic_terms")[0][0]

    count = 0
    for id_chunk in _chunks(protein_ids, chunk_size):
        _insert_semantic_terms_from_v2_sql(store, id_chunk)
        count += _count_semantic_terms_for_proteins(store, id_chunk)

    return count


def refresh_composite_predicates_from_v2(
    store: DuckDBStore,
    protein_ids: list[str] | None = None,
    chunk_size: int = 100_000,
    update_semantic_terms: bool = True,
    update_legacy_predicates: bool = False,
) -> int:
    """Recompute composite predicates from persisted atoms.

    Use this after editing ``config/predicates_v2/composites.yaml`` when the
    raw atoms are still valid. This avoids a full atom-generation pass while
    keeping semantic_state, semantic_terms, and optionally protein_predicates
    consistent.

    Args:
        store: DuckDB store with populated semantic_state and semantic_atoms.
        protein_ids: Optional subset of proteins to refresh.
        chunk_size: Number of semantic_state rows per batch.
        update_semantic_terms: Rebuild materialized search terms after updating
            semantic_state.composite_predicates.
        update_legacy_predicates: Rebuild the V2-derived compatibility cache.

    Returns:
        Number of semantic_state rows whose composite list was refreshed.
    """
    if chunk_size <= 0:
        raise ValueError("chunk_size must be positive")

    if protein_ids is not None:
        protein_ids = list(dict.fromkeys(protein_ids))
        if not protein_ids:
            return 0

    composites = load_composites()
    count = 0

    for state_rows in _iter_semantic_state_chunks(store, protein_ids, chunk_size):
        if not state_rows:
            continue

        chunk_protein_ids = [row[0] for row in state_rows]
        atoms_by_protein = _fetch_atoms_by_protein(store, chunk_protein_ids)
        rows = []

        for row in state_rows:
            state = _semantic_state_from_row(row)
            atoms = atoms_by_protein.get(state.protein_id, [])
            rows.append((
                state.protein_id,
                evaluate_composites(atoms, composites, state.topology),
            ))

        _persist_composite_predicate_updates(store, rows)
        count += len(rows)

    if update_semantic_terms:
        materialize_semantic_terms_from_v2(
            store,
            protein_ids=protein_ids,
            chunk_size=chunk_size,
        )

    if update_legacy_predicates:
        materialize_legacy_predicates_from_v2(
            store,
            protein_ids=protein_ids,
            chunk_size=chunk_size,
        )

    return count


def _iter_protein_chunks(
    store: DuckDBStore,
    protein_ids: list[str] | None,
    chunk_size: int,
) -> Iterable[list[tuple]]:
    """Yield protein rows in deterministic chunks."""
    select_cols = "protein_id, sequence_length, gc_content, contig_id, bin_id"

    if protein_ids is not None:
        for id_chunk in _chunks(protein_ids, chunk_size):
            placeholders = ",".join(["?"] * len(id_chunk))
            yield store.execute(
                f"SELECT {select_cols} FROM proteins "
                f"WHERE protein_id IN ({placeholders}) "
                f"ORDER BY protein_id",
                id_chunk,
            )
        return

    last_protein_id: str | None = None
    while True:
        if last_protein_id is None:
            rows = store.execute(
                f"SELECT {select_cols} FROM proteins "
                "ORDER BY protein_id LIMIT ?",
                [chunk_size],
            )
        else:
            rows = store.execute(
                f"SELECT {select_cols} FROM proteins "
                "WHERE protein_id > ? "
                "ORDER BY protein_id LIMIT ?",
                [last_protein_id, chunk_size],
            )
        if not rows:
            break
        yield rows
        last_protein_id = rows[-1][0]


def _fetch_annotations_by_protein(
    store: DuckDBStore,
    protein_ids: list[str],
) -> dict[str, list[AnnotationRecord]]:
    """Load annotations for one protein chunk and group them by protein."""
    if not protein_ids:
        return {}

    placeholders = ",".join(["?"] * len(protein_ids))
    rows = store.execute(
        f"SELECT protein_id, source, accession, name, description, evalue, score "
        f"FROM annotations WHERE protein_id IN ({placeholders})",
        protein_ids,
    )

    annotations_by_protein: dict[str, list[AnnotationRecord]] = {}
    for row in rows:
        pid = row[0]
        annotations_by_protein.setdefault(pid, []).append(
            AnnotationRecord(
                source=row[1],
                accession=row[2],
                name=row[3],
                description=row[4],
                evalue=row[5],
                score=row[6],
            )
        )
    return annotations_by_protein


def _iter_semantic_state_chunks(
    store: DuckDBStore,
    protein_ids: list[str] | None,
    chunk_size: int,
) -> Iterable[list[tuple]]:
    """Yield semantic_state rows in deterministic chunks."""
    select_cols = (
        "protein_id, activities, roles, architecture, localization, topology, "
        "size_class, quality_flags, composite_predicates"
    )

    if protein_ids is not None:
        for id_chunk in _chunks(protein_ids, chunk_size):
            placeholders = ",".join(["?"] * len(id_chunk))
            yield store.execute(
                f"SELECT {select_cols} FROM semantic_state "
                f"WHERE protein_id IN ({placeholders}) "
                "ORDER BY protein_id",
                id_chunk,
            )
        return

    last_protein_id: str | None = None
    while True:
        if last_protein_id is None:
            rows = store.execute(
                f"SELECT {select_cols} FROM semantic_state "
                "ORDER BY protein_id LIMIT ?",
                [chunk_size],
            )
        else:
            rows = store.execute(
                f"SELECT {select_cols} FROM semantic_state "
                "WHERE protein_id > ? "
                "ORDER BY protein_id LIMIT ?",
                [last_protein_id, chunk_size],
            )
        if not rows:
            break
        yield rows
        last_protein_id = rows[-1][0]


def _semantic_state_from_row(row: tuple) -> SemanticState:
    """Deserialize one semantic_state row."""
    topology = row[5]
    if isinstance(topology, str):
        topology = json.loads(topology) if topology else {}
    elif topology is None:
        topology = {}

    return SemanticState(
        protein_id=row[0],
        activities=_as_list(row[1]),
        roles=_as_list(row[2]),
        architecture=_as_list(row[3]),
        localization=_as_list(row[4]),
        topology=topology,
        size_class=row[6] or "",
        quality_flags=_as_list(row[7]),
        composite_predicates=_as_list(row[8]),
    )


def _fetch_atoms_by_protein(
    store: DuckDBStore,
    protein_ids: list[str],
) -> dict[str, list[SemanticAtom]]:
    """Load persisted atoms for one protein chunk."""
    if not protein_ids:
        return {}

    placeholders = ",".join(["?"] * len(protein_ids))
    rows = store.execute(
        f"""
        SELECT protein_id, atom_id, facet, relation, source_accession,
               source_db, evidence_evalue, evidence_score
        FROM semantic_atoms
        WHERE protein_id IN ({placeholders})
        ORDER BY protein_id, atom_id, source_db, source_accession
        """,
        protein_ids,
    )

    atoms_by_protein: dict[str, list[SemanticAtom]] = {}
    for row in rows:
        atoms_by_protein.setdefault(row[0], []).append(
            SemanticAtom(
                protein_id=row[0],
                atom_id=row[1],
                facet=SemanticFacet(row[2]),
                relation=ClaimRelation(row[3]),
                source_accession=row[4],
                source_db=row[5],
                evidence_evalue=row[6],
                evidence_score=row[7],
            )
        )

    return atoms_by_protein


def _fetch_context_atoms_by_protein(
    store: DuckDBStore,
    protein_ids: list[str],
) -> dict[str, list[SemanticAtom]]:
    """Load context-derived atoms that are not annotation rows."""
    if not protein_ids:
        return {}

    placeholders = ",".join(["?"] * len(protein_ids))
    rows = store.execute(
        f"""
        SELECT DISTINCT p.protein_id, l.locus_id
        FROM proteins p
        JOIN loci l ON p.contig_id = l.contig_id
        WHERE p.protein_id IN ({placeholders})
          AND l.locus_type = 'crispr'
          AND p.start < l.end_coord
          AND p.end_coord > l.start
        ORDER BY p.protein_id, l.locus_id
        """,
        protein_ids,
    )

    atoms_by_protein: dict[str, list[SemanticAtom]] = {}
    for protein_id, locus_id in rows:
        atoms_by_protein.setdefault(protein_id, []).append(
            SemanticAtom(
                protein_id=protein_id,
                atom_id="in_crispr_array",
                facet=SemanticFacet.quality_flag,
                relation=ClaimRelation.implies,
                source_accession=locus_id,
                source_db="_locus",
            )
        )

    return atoms_by_protein


def _clear_v2_tables(
    store: DuckDBStore,
    *,
    update_legacy_predicates: bool,
) -> None:
    """Clear V2 tables before a full-dataset refresh."""
    store.execute("DELETE FROM semantic_atoms;")
    store.execute("DELETE FROM semantic_state;")
    store.execute("DELETE FROM semantic_terms;")
    if update_legacy_predicates:
        store.execute("DELETE FROM protein_predicates;")


def _delete_for_proteins(
    store: DuckDBStore,
    protein_ids: list[str],
    *,
    update_legacy_predicates: bool,
) -> None:
    """Delete persisted V2 state for a subset without touching other proteins."""
    for id_chunk in _chunks(protein_ids, 10_000):
        placeholders = ",".join(["?"] * len(id_chunk))
        store.execute(
            f"DELETE FROM semantic_atoms WHERE protein_id IN ({placeholders})",
            id_chunk,
        )
        store.execute(
            f"DELETE FROM semantic_state WHERE protein_id IN ({placeholders})",
            id_chunk,
        )
        store.execute(
            f"DELETE FROM semantic_terms WHERE protein_id IN ({placeholders})",
            id_chunk,
        )
        if update_legacy_predicates:
            store.execute(
                f"DELETE FROM protein_predicates WHERE protein_id IN ({placeholders})",
                id_chunk,
            )


def _delete_semantic_terms_for_proteins(
    store: DuckDBStore,
    protein_ids: list[str],
) -> None:
    """Delete materialized semantic terms for a subset."""
    for id_chunk in _chunks(protein_ids, 10_000):
        placeholders = ",".join(["?"] * len(id_chunk))
        store.execute(
            f"DELETE FROM semantic_terms WHERE protein_id IN ({placeholders})",
            id_chunk,
        )


def _backfill_semantic_terms_if_empty(store: DuckDBStore) -> None:
    """Backfill old V2 databases before a subset refresh makes terms partial."""
    term_count = store.execute("SELECT COUNT(*) FROM semantic_terms")[0][0]
    if term_count:
        return

    state_count = store.execute("SELECT COUNT(*) FROM semantic_state")[0][0]
    if state_count:
        materialize_semantic_terms_from_v2(store)


def _protein_filter_sql(protein_ids: list[str] | None, table_alias: str = "") -> str:
    """Build an IN-filter fragment for a non-empty protein ID chunk."""
    if not protein_ids:
        return ""
    column = f"{table_alias}.protein_id" if table_alias else "protein_id"
    placeholders = ",".join(["?"] * len(protein_ids))
    return f" AND {column} IN ({placeholders})"


def _insert_semantic_terms_from_v2_sql(
    store: DuckDBStore,
    protein_ids: list[str] | None = None,
) -> None:
    """Backfill semantic_terms using set-oriented DuckDB INSERTs."""
    protein_filter = _protein_filter_sql(protein_ids)
    params = protein_ids or []

    store.execute(
        f"""
        INSERT INTO semantic_terms (
            protein_id, term_id, term_kind, facet, relation,
            source_db, source_accession
        )
        SELECT DISTINCT
            protein_id,
            atom_id AS term_id,
            'atom' AS term_kind,
            facet,
            relation,
            COALESCE(source_db, '') AS source_db,
            COALESCE(source_accession, '') AS source_accession
        FROM semantic_atoms
        WHERE atom_id NOT LIKE '_source_witness:%'
        {protein_filter}
        """,
        params,
    )

    source_prefix_case = direct_access_prefix_case_sql("source_db")
    store.execute(
        f"""
        INSERT INTO semantic_terms (
            protein_id, term_id, term_kind, facet, relation,
            source_db, source_accession
        )
        WITH direct_terms AS (
            SELECT
                protein_id,
                {source_prefix_case} AS source_prefix,
                facet,
                relation,
                COALESCE(source_db, '') AS source_db,
                source_accession
            FROM semantic_atoms
            WHERE atom_id NOT LIKE '_source_witness:%'
              AND source_accession IS NOT NULL
              AND LEFT(source_accession, 1) != '_'
              {protein_filter}
        )
        SELECT DISTINCT
            protein_id,
            source_prefix || ':' || source_accession AS term_id,
            'direct_access' AS term_kind,
            facet,
            relation,
            source_db,
            source_accession
        FROM direct_terms
        WHERE source_prefix IS NOT NULL
        """,
        params,
    )

    state_filter = _protein_filter_sql(protein_ids)
    store.execute(
        f"""
        INSERT INTO semantic_terms (
            protein_id, term_id, term_kind, facet, relation,
            source_db, source_accession
        )
        SELECT DISTINCT
            protein_id,
            composite AS term_id,
            'composite' AS term_kind,
            NULL AS facet,
            'implies' AS relation,
            '_composite' AS source_db,
            composite AS source_accession
        FROM (
            SELECT protein_id, UNNEST(composite_predicates) AS composite
            FROM semantic_state
            WHERE LEN(composite_predicates) > 0
            {state_filter}
        ) composites
        """,
        params,
    )


def _count_semantic_terms_for_proteins(
    store: DuckDBStore,
    protein_ids: list[str],
) -> int:
    """Count semantic_terms rows for one protein ID chunk."""
    if not protein_ids:
        return 0
    placeholders = ",".join(["?"] * len(protein_ids))
    return store.execute(
        f"SELECT COUNT(*) FROM semantic_terms WHERE protein_id IN ({placeholders})",
        protein_ids,
    )[0][0]


def _delete_legacy_for_proteins(
    store: DuckDBStore,
    protein_ids: list[str],
) -> None:
    """Delete compatibility cache rows for a subset."""
    for id_chunk in _chunks(protein_ids, 10_000):
        placeholders = ",".join(["?"] * len(id_chunk))
        store.execute(
            f"DELETE FROM protein_predicates WHERE protein_id IN ({placeholders})",
            id_chunk,
        )


def _as_list(value: object) -> list:
    """Normalize DuckDB list values to Python lists."""
    if value is None:
        return []
    if isinstance(value, list):
        return value
    if isinstance(value, tuple):
        return list(value)
    if isinstance(value, str):
        try:
            parsed = json.loads(value)
            return parsed if isinstance(parsed, list) else [value]
        except json.JSONDecodeError:
            return [value]
    return list(value)


def _chunks(items: list[str], chunk_size: int) -> Iterable[list[str]]:
    """Yield fixed-size chunks from a list."""
    for start in range(0, len(items), chunk_size):
        yield items[start:start + chunk_size]


def _persist_chunk(
    store: DuckDBStore,
    atoms: list[SemanticAtom],
    states: dict[str, SemanticState],
    *,
    legacy_rows: list[tuple[str, list[str]]],
    update_legacy_predicates: bool,
) -> None:
    """Persist one generation chunk inside one DuckDB transaction."""
    store.conn.execute("BEGIN TRANSACTION;")
    try:
        _persist_atoms(store, atoms)
        _persist_states(store, states)
        _persist_semantic_terms(store, _build_semantic_terms(atoms, states))
        if update_legacy_predicates:
            _persist_legacy_predicates(store, legacy_rows)
    except BaseException:
        store.conn.execute("ROLLBACK;")
        raise
    else:
        store.conn.execute("COMMIT;")


def _persist_atoms(store: DuckDBStore, atoms: list[SemanticAtom]) -> None:
    """Write atoms to semantic_atoms table using batch insert."""
    if not atoms:
        return

    # Batch insert via executemany for ~10-50x speedup over per-row inserts
    rows_by_key: dict[tuple[str, str, str], tuple] = {}
    for atom in atoms:
        rows_by_key[(atom.protein_id, atom.atom_id, atom.source_accession)] = (
            atom.protein_id,
            atom.atom_id,
            atom.facet.value,
            atom.relation.value,
            atom.source_accession,
            atom.source_db,
            atom.evidence_evalue,
            atom.evidence_score,
        )
    rows = list(rows_by_key.values())
    _insert_dataframe(
        store,
        view_name="_semantic_atoms_batch",
        rows=rows,
        columns=[
            "protein_id",
            "atom_id",
            "facet",
            "relation",
            "source_accession",
            "source_db",
            "evidence_evalue",
            "evidence_score",
        ],
        sql="""
            INSERT INTO semantic_atoms
            (protein_id, atom_id, facet, relation, source_accession,
             source_db, evidence_evalue, evidence_score)
            SELECT protein_id, atom_id, facet, relation, source_accession,
                   source_db, evidence_evalue, evidence_score
            FROM _semantic_atoms_batch
        """,
    )


def _persist_states(
    store: DuckDBStore,
    states: dict[str, SemanticState],
) -> None:
    """Write aggregated states to semantic_state table using batch insert."""
    if not states:
        return

    rows = [
        (
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
        )
        for pid, state in states.items()
    ]
    _insert_dataframe(
        store,
        view_name="_semantic_state_batch",
        rows=rows,
        columns=[
            "protein_id",
            "activities",
            "roles",
            "architecture",
            "localization",
            "topology",
            "size_class",
            "quality_flags",
            "composite_predicates",
            "unresolved_count",
        ],
        sql="""
            INSERT INTO semantic_state
            (protein_id, activities, roles, architecture, localization,
             topology, size_class, quality_flags, composite_predicates,
             unresolved_count, updated_at)
            SELECT protein_id, activities, roles, architecture, localization,
                   topology, size_class, quality_flags, composite_predicates,
                   unresolved_count, CURRENT_TIMESTAMP
            FROM _semantic_state_batch
        """,
    )


def _build_semantic_terms(
    atoms: list[SemanticAtom],
    states: dict[str, SemanticState],
) -> list[tuple[str, str, str, str | None, str | None, str, str]]:
    """Build unified search terms from raw atoms plus resolved composites."""
    rows_by_key: dict[
        tuple[str, str, str, str | None, str | None, str, str],
        tuple[str, str, str, str | None, str | None, str, str],
    ] = {}

    def add(
        protein_id: str,
        term_id: str,
        term_kind: str,
        facet: str | None,
        relation: str | None,
        source_db: str | None,
        source_accession: str | None,
    ) -> None:
        row = (
            protein_id,
            term_id,
            term_kind,
            facet,
            relation,
            source_db or "",
            source_accession or "",
        )
        rows_by_key[row] = row

    for atom in atoms:
        if not atom.atom_id.startswith("_source_witness:"):
            add(
                atom.protein_id,
                atom.atom_id,
                "atom",
                atom.facet.value,
                atom.relation.value,
                atom.source_db,
                atom.source_accession,
            )

        if direct_access := direct_access_predicate_from_atom(atom):
            add(
                atom.protein_id,
                direct_access,
                "direct_access",
                atom.facet.value,
                atom.relation.value,
                atom.source_db,
                atom.source_accession,
            )

    for state in states.values():
        for composite in state.composite_predicates:
            add(
                state.protein_id,
                composite,
                "composite",
                None,
                "implies",
                "_composite",
                composite,
            )

    return list(rows_by_key.values())


def _persist_semantic_terms(
    store: DuckDBStore,
    rows: list[tuple[str, str, str, str | None, str | None, str, str]],
) -> None:
    """Write materialized V2 search terms using batch insert."""
    if not rows:
        return

    _insert_dataframe(
        store,
        view_name="_semantic_terms_batch",
        rows=rows,
        columns=[
            "protein_id",
            "term_id",
            "term_kind",
            "facet",
            "relation",
            "source_db",
            "source_accession",
        ],
        sql="""
            INSERT INTO semantic_terms
            (protein_id, term_id, term_kind, facet, relation,
             source_db, source_accession)
            SELECT protein_id, term_id, term_kind, facet, relation,
                   source_db, source_accession
            FROM _semantic_terms_batch
        """,
    )


def _persist_legacy_predicates(
    store: DuckDBStore,
    rows: list[tuple[str, list[str]]],
) -> None:
    """Write V2-derived flat predicates into the V1 compatibility cache."""
    if not rows:
        return

    _insert_dataframe(
        store,
        view_name="_protein_predicates_batch",
        rows=rows,
        columns=["protein_id", "predicates"],
        sql="""
            INSERT INTO protein_predicates
            (protein_id, predicates, updated_at)
            SELECT protein_id, predicates, CURRENT_TIMESTAMP
            FROM _protein_predicates_batch
        """,
    )


def _persist_composite_predicate_updates(
    store: DuckDBStore,
    rows: list[tuple[str, list[str]]],
) -> None:
    """Update semantic_state composite_predicates using a DataFrame batch."""
    if not rows:
        return

    _insert_dataframe(
        store,
        view_name="_composite_predicates_batch",
        rows=rows,
        columns=["protein_id", "composite_predicates"],
        sql="""
            UPDATE semantic_state
            SET
                composite_predicates = batch.composite_predicates,
                updated_at = CURRENT_TIMESTAMP
            FROM _composite_predicates_batch AS batch
            WHERE semantic_state.protein_id = batch.protein_id
        """,
    )


def _insert_dataframe(
    store: DuckDBStore,
    *,
    view_name: str,
    rows: list[tuple],
    columns: list[str],
    sql: str,
) -> None:
    """Insert rows through a registered DataFrame instead of row-wise binds."""
    if not rows:
        return

    frame = pd.DataFrame.from_records(rows, columns=columns)
    store.conn.register(view_name, frame)
    try:
        store.conn.execute(sql)
    finally:
        store.conn.unregister(view_name)


def _get_contig_gc_stats(
    store: DuckDBStore,
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


__all__ = [
    "generate_and_persist_v2",
    "materialize_legacy_predicates_from_v2",
    "materialize_semantic_terms_from_v2",
    "refresh_composite_predicates_from_v2",
]
