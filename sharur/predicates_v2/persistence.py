"""
Persistence: generate atoms, aggregate, evaluate composites, write to DuckDB.

Orchestrates the full V2 pipeline and persists results to the semantic_atoms,
semantic_state, and semantic_terms tables.
"""

from __future__ import annotations

import hashlib
import json
import logging
import math
import multiprocessing
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pandas as pd

from sharur.predicates.generator import AnnotationRecord, ProteinRecord
from sharur.predicates_v2.aggregator import aggregate_atoms
from sharur.predicates_v2.compat import (
    direct_access_predicate_from_atom,
    direct_access_prefix_case_sql,
    semantic_state_to_predicates,
)
from sharur.predicates_v2.composites import evaluate_composites, load_composites
from sharur.predicates_v2.generator import AtomGenerator
from sharur.predicates_v2.model import (
    ClaimRelation,
    SemanticAtom,
    SemanticFacet,
    SemanticState,
    create_v2_indexes,
    create_v2_tables,
    drop_v2_indexes,
)
from sharur.predicates_v2.review_queue import (
    format_review_queue_tsv,
    suggest_facet,
)
from sharur.predicates_v2.validated_systems import (
    fetch_validated_system_annotations,
    materialize_system_proteins,
)


if TYPE_CHECKING:
    from collections.abc import Iterable

    from sharur.storage.duckdb_store import DuckDBStore


logger = logging.getLogger(__name__)

_GENERATION_KEY = "full_v2"
_GENERATION_CONTRACT_VERSION = 2
_GENERATION_TABLES = (
    "v2_generation_atoms",
    "v2_generation_state",
    "v2_generation_terms",
    "v2_generation_legacy",
)


@dataclass
class _ProteinTransformInput:
    """Pickle-safe input for one protein's pure semantic transform."""

    protein_id: str
    sequence_length: int | None
    gc_content: float | None
    contig_gc_mean: float | None
    contig_gc_std: float | None
    annotations: list[AnnotationRecord]
    context_atoms: list[SemanticAtom]


@dataclass
class _TransformBatchResult:
    """Ordered output from one worker microbatch."""

    atoms: list[SemanticAtom]
    states: dict[str, SemanticState]
    legacy_rows: list[tuple[str, list[str]]]


@dataclass
class _TransformWorkerState:
    """Process-local immutable rule state initialized once per worker."""

    generator: AtomGenerator | None = None
    composites: list[Any] | None = None
    update_legacy_predicates: bool = False


_WORKER_STATE = _TransformWorkerState()


def generate_and_persist_v2(
    store: DuckDBStore,
    protein_ids: list[str] | None = None,
    output_review_queue: str | None = None,
    chunk_size: int = 25_000,
    workers: int = 1,
    worker_batch_size: int | None = None,
    resume: bool = False,
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
        workers: Number of persistent transform processes. DuckDB reads and
            writes remain in this parent process.
        worker_batch_size: Proteins per worker task. By default, each database
            chunk is split into about two tasks per worker, capped at 5,000.
        resume: Resume the last full-dataset generation checkpoint. Resume is
            supported for full refreshes with ``return_states=False``.
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
    if chunk_size <= 0:
        raise ValueError("chunk_size must be positive")
    if workers <= 0:
        raise ValueError("workers must be positive")
    if worker_batch_size is not None and worker_batch_size <= 0:
        raise ValueError("worker_batch_size must be positive")
    if resume and protein_ids is not None:
        raise ValueError("resume is supported only for full-dataset generation")
    if resume and return_states:
        raise ValueError("resume requires return_states=False")

    # Subset refreshes keep the canonical indexes online. Full refreshes write
    # into constraint-free generation tables and restore canonical indexes
    # once, after bulk promotion.
    create_v2_tables(store, create_indexes=protein_ids is not None)

    if protein_ids is not None:
        protein_ids = list(dict.fromkeys(protein_ids))
        if not protein_ids:
            return {}
        _backfill_semantic_terms_if_empty(store)

    materialize_system_proteins(store)

    start_after: str | None = None
    total_processed = 0
    checkpoint_complete = False
    checkpoint_enabled = protein_ids is None
    semantic_fingerprint = ""
    input_signature = ""

    # Full refreshes use a transaction-backed checkpoint. Subset refreshes
    # retain their existing replacement semantics and do not participate in
    # the full-dataset checkpoint.
    if protein_ids is None:
        semantic_fingerprint = _semantic_generation_fingerprint(
            predict_topology=predict_topology,
            update_legacy_predicates=update_legacy_predicates,
        )
        input_signature = _generation_input_signature(store)
        if resume:
            start_after, total_processed, checkpoint_complete = _resume_full_generation(
                store,
                semantic_fingerprint=semantic_fingerprint,
                input_signature=input_signature,
            )
        else:
            _initialize_full_generation(
                store,
                semantic_fingerprint=semantic_fingerprint,
                input_signature=input_signature,
                total_count=_protein_count_from_signature(input_signature),
                update_legacy_predicates=update_legacy_predicates,
            )
    else:
        _delete_for_proteins(
            store,
            protein_ids,
            update_legacy_predicates=update_legacy_predicates,
        )

    total_count = (
        _protein_count_from_signature(input_signature)
        if protein_ids is None
        else len(protein_ids)
    )
    initial_processed = total_processed
    results: dict[str, SemanticState] = {}
    started_at = time.monotonic()
    executor: ProcessPoolExecutor | None = None
    serial_generator: AtomGenerator | None = None
    serial_composites: list[Any] | None = None

    try:
        validated_system_annotations = fetch_validated_system_annotations(
            store,
            protein_ids=set(protein_ids) if protein_ids is not None else None,
        )

        # Load GC statistics after checkpoint validation so incompatible
        # resume attempts remain cheap.
        gc_stats = _get_contig_gc_stats(store)

        if workers == 1:
            serial_generator = AtomGenerator(predict_topology=predict_topology)
            serial_composites = load_composites()
        else:
            executor = ProcessPoolExecutor(
                max_workers=workers,
                mp_context=multiprocessing.get_context("spawn"),
                initializer=_initialize_transform_worker,
                initargs=(predict_topology, update_legacy_predicates),
            )

        for protein_rows in _iter_protein_chunks(
            store,
            protein_ids,
            chunk_size,
            start_after=start_after,
        ):
            if not protein_rows:
                continue

            chunk_started_at = time.monotonic()
            chunk_protein_ids = [row[0] for row in protein_rows]
            annotations_by_protein = _fetch_annotations_by_protein(
                store,
                chunk_protein_ids,
            )
            context_atoms_by_protein = _fetch_context_atoms_by_protein(
                store,
                chunk_protein_ids,
            )
            inputs = _prepare_transform_inputs(
                protein_rows,
                annotations_by_protein=annotations_by_protein,
                validated_system_annotations=validated_system_annotations,
                context_atoms_by_protein=context_atoms_by_protein,
                gc_stats=gc_stats,
            )
            prepared_at = time.monotonic()
            batch_size = _resolve_worker_batch_size(
                len(inputs),
                workers=workers,
                requested=worker_batch_size,
            )
            batches = list(_chunks(inputs, batch_size))

            if executor is None:
                assert serial_generator is not None
                assert serial_composites is not None
                batch_results = [
                    _transform_batch(
                        batch,
                        generator=serial_generator,
                        composites=serial_composites,
                        update_legacy_predicates=update_legacy_predicates,
                    )
                    for batch in batches
                ]
            else:
                # map() preserves input order; the parent merges and commits
                # every database chunk deterministically through one writer.
                batch_results = list(executor.map(_transform_batch_worker, batches))

            chunk_atoms: list[SemanticAtom] = []
            chunk_states: dict[str, SemanticState] = {}
            legacy_rows: list[tuple[str, list[str]]] = []
            for batch_result in batch_results:
                chunk_atoms.extend(batch_result.atoms)
                chunk_states.update(batch_result.states)
                legacy_rows.extend(batch_result.legacy_rows)
            transformed_at = time.monotonic()

            chunk_count = len(protein_rows)
            chunk_atom_count = len(chunk_atoms)
            checkpoint_last_protein_id = protein_rows[-1][0]
            del (
                annotations_by_protein,
                batch_results,
                batches,
                context_atoms_by_protein,
                inputs,
                protein_rows,
            )

            next_processed = total_processed + chunk_count
            _persist_chunk(
                store,
                chunk_atoms,
                chunk_states,
                legacy_rows=legacy_rows,
                update_legacy_predicates=update_legacy_predicates,
                use_generation_tables=checkpoint_enabled,
                checkpoint_processed=(next_processed if checkpoint_enabled else None),
                checkpoint_last_protein_id=(
                    checkpoint_last_protein_id if checkpoint_enabled else None
                ),
            )
            persisted_at = time.monotonic()
            total_processed = next_processed
            if return_states:
                results.update(chunk_states)

            elapsed = max(time.monotonic() - started_at, 1e-9)
            run_processed = total_processed - initial_processed
            rate = run_processed / elapsed if run_processed > 0 else 0.0
            remaining = max(total_count - total_processed, 0)
            eta_seconds = remaining / rate if rate > 0 else 0.0
            logger.info(
                "V2 progress: %s/%s proteins (%.2f%%), %.1f proteins/s, ETA %.1fh; "
                "chunk=%s proteins/%s atoms prepare=%.1fs transform=%.1fs persist=%.1fs total=%.1fs",
                f"{total_processed:,}",
                f"{total_count:,}",
                (100.0 * total_processed / total_count) if total_count else 100.0,
                rate,
                eta_seconds / 3600.0,
                f"{chunk_count:,}",
                f"{chunk_atom_count:,}",
                prepared_at - chunk_started_at,
                transformed_at - prepared_at,
                persisted_at - transformed_at,
                persisted_at - chunk_started_at,
            )

        if checkpoint_enabled:
            state_table = (
                "semantic_state" if checkpoint_complete else "v2_generation_state"
            )
            persisted_states = int(
                store.execute(f"SELECT COUNT(*) FROM {state_table}")[0][0]
            )
            if total_processed != total_count or persisted_states != total_count:
                raise RuntimeError(
                    "V2 full refresh ended before every protein state was committed"
                )

            if not checkpoint_complete:
                _promote_full_generation(
                    store,
                    update_legacy_predicates=update_legacy_predicates,
                    expected_count=total_count,
                )

            index_started = time.monotonic()
            create_v2_indexes(store)
            logger.info(
                "V2 canonical indexes ready in %.1fs",
                time.monotonic() - index_started,
            )

        if output_review_queue:
            queue = _build_review_queue_from_store(
                store,
                total_proteins=(total_count if protein_ids is None else total_processed),
                protein_ids=protein_ids,
            )
            review_path = Path(output_review_queue)
            review_path.parent.mkdir(parents=True, exist_ok=True)
            review_path.write_text(format_review_queue_tsv(queue))

        if checkpoint_enabled:
            _mark_full_generation_complete(store)
            _clear_generation_tables_best_effort(store)
    except BaseException:
        if checkpoint_enabled:
            _mark_full_generation_failed(store)
        raise
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)

    return results


def _initialize_transform_worker(
    predict_topology: bool,
    update_legacy_predicates: bool,
) -> None:
    """Initialize immutable rule state once in each spawned process."""
    _WORKER_STATE.generator = AtomGenerator(predict_topology=predict_topology)
    _WORKER_STATE.composites = load_composites()
    _WORKER_STATE.update_legacy_predicates = update_legacy_predicates


def _transform_batch_worker(
    inputs: list[_ProteinTransformInput],
) -> _TransformBatchResult:
    """Transform one ordered microbatch in a spawned process."""
    if _WORKER_STATE.generator is None or _WORKER_STATE.composites is None:
        raise RuntimeError("V2 transform worker was not initialized")
    return _transform_batch(
        inputs,
        generator=_WORKER_STATE.generator,
        composites=_WORKER_STATE.composites,
        update_legacy_predicates=_WORKER_STATE.update_legacy_predicates,
    )


def _transform_batch(
    inputs: list[_ProteinTransformInput],
    *,
    generator: AtomGenerator,
    composites: list[Any],
    update_legacy_predicates: bool,
) -> _TransformBatchResult:
    """Run the CPU-only semantic transformation for an ordered microbatch."""
    atoms_out: list[SemanticAtom] = []
    states_out: dict[str, SemanticState] = {}
    legacy_rows: list[tuple[str, list[str]]] = []

    for item in inputs:
        protein = ProteinRecord(
            protein_id=item.protein_id,
            sequence_length=item.sequence_length,
            gc_content=item.gc_content,
            contig_gc_mean=item.contig_gc_mean,
            contig_gc_std=item.contig_gc_std,
        )
        atoms = generator.generate_atoms(protein, item.annotations)
        atoms.extend(item.context_atoms)
        state = aggregate_atoms(item.protein_id, atoms)
        state.composite_predicates = evaluate_composites(
            atoms,
            composites,
            state.topology,
        )

        atoms_out.extend(atoms)
        states_out[item.protein_id] = state
        if update_legacy_predicates:
            legacy_rows.append((
                item.protein_id,
                semantic_state_to_predicates(state, atoms=atoms),
            ))

    return _TransformBatchResult(
        atoms=atoms_out,
        states=states_out,
        legacy_rows=legacy_rows,
    )


def _prepare_transform_inputs(
    protein_rows: list[tuple],
    *,
    annotations_by_protein: dict[str, list[AnnotationRecord]],
    validated_system_annotations: dict[str, list[AnnotationRecord]],
    context_atoms_by_protein: dict[str, list[SemanticAtom]],
    gc_stats: dict[str, tuple[float, float]],
) -> list[_ProteinTransformInput]:
    """Prepare one bounded database chunk for worker serialization."""
    inputs: list[_ProteinTransformInput] = []
    for row in protein_rows:
        protein_id, sequence_length, gc_content, contig_id = row[:4]
        gc_mean, gc_std = (
            gc_stats.get(contig_id, (None, None))
            if contig_id
            else (None, None)
        )
        annotations = list(annotations_by_protein.get(protein_id, ()))
        annotations.extend(validated_system_annotations.get(protein_id, ()))
        inputs.append(_ProteinTransformInput(
            protein_id=protein_id,
            sequence_length=sequence_length,
            gc_content=gc_content,
            contig_gc_mean=gc_mean,
            contig_gc_std=gc_std,
            annotations=annotations,
            context_atoms=list(context_atoms_by_protein.get(protein_id, ())),
        ))
    return inputs


def _resolve_worker_batch_size(
    input_count: int,
    *,
    workers: int,
    requested: int | None,
) -> int:
    """Choose bounded worker shards while retaining enough scheduling slack."""
    if requested is not None:
        return requested
    return max(1, min(5_000, math.ceil(input_count / max(workers * 2, 1))))


def _semantic_generation_fingerprint(
    *,
    predict_topology: bool,
    update_legacy_predicates: bool,
) -> str:
    """Hash semantic code/config plus output-affecting generation options."""
    digest = hashlib.sha256()
    options = {
        "contract_version": _GENERATION_CONTRACT_VERSION,
        "predict_topology": predict_topology,
        "update_legacy_predicates": update_legacy_predicates,
    }
    digest.update(json.dumps(options, sort_keys=True).encode())

    package_root = Path(__file__).resolve().parent
    trees = [
        ("predicates_v2", package_root),
        ("predicates_v1", package_root.parent / "predicates"),
        ("source_config", package_root.parent.parent / "config" / "predicates_v2"),
    ]
    allowed_suffixes = {".py", ".yaml", ".yml", ".json", ".tsv"}
    for label, root in trees:
        if not root.exists():
            continue
        paths = [root] if root.is_file() else sorted(root.rglob("*"))
        for path in paths:
            if not path.is_file() or path.suffix not in allowed_suffixes:
                continue
            digest.update(f"{label}/{path.relative_to(root)}\0".encode())
            digest.update(path.read_bytes())
            digest.update(b"\0")
    return digest.hexdigest()


def _generation_input_signature(store: DuckDBStore) -> str:
    """Capture a lightweight source-table identity for safe resume."""
    protein = store.execute("""
        SELECT COUNT(*), MIN(protein_id), MAX(protein_id),
               COALESCE(SUM(sequence_length), 0)
        FROM proteins
    """)[0]
    annotations = store.execute("""
        SELECT COUNT(*), COALESCE(MAX(annotation_id), 0)
        FROM annotations
    """)[0]
    signature = {
        "protein_count": int(protein[0]),
        "min_protein_id": protein[1],
        "max_protein_id": protein[2],
        "sequence_length_sum": int(protein[3]),
        "annotation_count": int(annotations[0]),
        "max_annotation_id": int(annotations[1]),
        "locus_count": int(store.execute("SELECT COUNT(*) FROM loci")[0][0]),
        "system_protein_count": int(
            store.execute("SELECT COUNT(*) FROM system_proteins")[0][0]
        ),
    }
    return json.dumps(signature, sort_keys=True, separators=(",", ":"))


def _protein_count_from_signature(input_signature: str) -> int:
    """Read the protein count from a checkpoint input signature."""
    return int(json.loads(input_signature)["protein_count"])


def _initialize_full_generation(
    store: DuckDBStore,
    *,
    semantic_fingerprint: str,
    input_signature: str,
    total_count: int,
    update_legacy_predicates: bool,
) -> None:
    """Atomically clear full-refresh products and create its checkpoint."""
    store.conn.execute("BEGIN TRANSACTION;")
    try:
        drop_v2_indexes(store)
        _clear_v2_tables(
            store,
            update_legacy_predicates=update_legacy_predicates,
        )
        _clear_generation_tables(store)
        store.execute(
            "DELETE FROM v2_generation_checkpoint WHERE generation_key = ?",
            [_GENERATION_KEY],
        )
        store.execute(
            """
            INSERT INTO v2_generation_checkpoint (
                generation_key, semantic_fingerprint, input_signature,
                status, last_protein_id, processed_count, total_count,
                started_at, updated_at, completed_at
            ) VALUES (?, ?, ?, 'running', NULL, 0, ?,
                      CURRENT_TIMESTAMP, CURRENT_TIMESTAMP, NULL)
            """,
            [_GENERATION_KEY, semantic_fingerprint, input_signature, total_count],
        )
    except BaseException:
        store.conn.execute("ROLLBACK;")
        raise
    else:
        store.conn.execute("COMMIT;")


def _resume_full_generation(
    store: DuckDBStore,
    *,
    semantic_fingerprint: str,
    input_signature: str,
) -> tuple[str | None, int, bool]:
    """Validate and reopen a full-refresh checkpoint."""
    rows = store.execute(
        """
        SELECT semantic_fingerprint, input_signature, status,
               last_protein_id, processed_count, total_count
        FROM v2_generation_checkpoint
        WHERE generation_key = ?
        """,
        [_GENERATION_KEY],
    )
    if not rows:
        raise ValueError("No full V2 generation checkpoint is available to resume")

    stored_fingerprint, stored_input, status, last_protein_id, processed, total = rows[0]
    if stored_fingerprint != semantic_fingerprint:
        raise ValueError(
            "V2 semantic code/config changed since the checkpoint; start a fresh run"
        )
    if stored_input != input_signature:
        raise ValueError(
            "V2 source tables changed since the checkpoint; start a fresh run"
        )
    if status not in {"running", "failed", "complete"}:
        raise ValueError(f"Unsupported V2 checkpoint status: {status!r}")
    if int(total) != _protein_count_from_signature(input_signature):
        raise ValueError("V2 checkpoint protein total disagrees with source tables")

    processed = int(processed)
    complete = status == "complete"
    state_table = "semantic_state" if complete else "v2_generation_state"
    persisted_count, persisted_max = store.execute(
        f"SELECT COUNT(*), MAX(protein_id) FROM {state_table}"
    )[0]
    if int(persisted_count) != processed or persisted_max != last_protein_id:
        raise ValueError(
            f"V2 checkpoint disagrees with persisted {state_table}; start a fresh run"
        )
    if last_protein_id is not None:
        prefix_count = store.execute(
            "SELECT COUNT(*) FROM proteins WHERE protein_id <= ?",
            [last_protein_id],
        )[0][0]
        if int(prefix_count) != processed:
            raise ValueError(
                "V2 checkpoint boundary disagrees with protein ordering; start a fresh run"
            )

    if not complete:
        drop_v2_indexes(store)
        store.execute(
            """
            UPDATE v2_generation_checkpoint
            SET status = 'running', updated_at = CURRENT_TIMESTAMP,
                completed_at = NULL
            WHERE generation_key = ?
            """,
            [_GENERATION_KEY],
        )
    logger.info(
        "Resuming V2 generation after %s (%s/%s proteins committed)",
        last_protein_id,
        f"{processed:,}",
        f"{int(total):,}",
    )
    return last_protein_id, processed, complete


def _mark_full_generation_complete(store: DuckDBStore) -> None:
    """Mark the full refresh complete after all derived artifacts exist."""
    store.execute(
        """
        UPDATE v2_generation_checkpoint
        SET status = 'complete', updated_at = CURRENT_TIMESTAMP,
            completed_at = CURRENT_TIMESTAMP
        WHERE generation_key = ?
        """,
        [_GENERATION_KEY],
    )


def _mark_full_generation_failed(store: DuckDBStore) -> None:
    """Record a recoverable failure while preserving the committed boundary."""
    try:
        store.execute(
            """
            UPDATE v2_generation_checkpoint
            SET status = 'failed', updated_at = CURRENT_TIMESTAMP
            WHERE generation_key = ? AND status != 'complete'
            """,
            [_GENERATION_KEY],
        )
    except Exception:
        logger.exception("Failed to record V2 checkpoint failure state")


def _build_review_queue_from_store(
    store: DuckDBStore,
    *,
    total_proteins: int,
    protein_ids: list[str] | None,
) -> list[dict[str, Any]]:
    """Aggregate unresolved accessions set-wise from persisted atoms."""
    join_filter = ""
    registered_filter = False
    if protein_ids is not None:
        frame = pd.DataFrame({"protein_id": protein_ids})
        store.conn.register("_v2_review_proteins", frame)
        registered_filter = True
        join_filter = (
            "JOIN _v2_review_proteins selected "
            "ON selected.protein_id = atoms.protein_id"
        )

    try:
        rows = store.execute(f"""
            WITH unresolved AS (
                SELECT DISTINCT
                    atoms.source_db,
                    atoms.source_accession,
                    atoms.protein_id,
                    COALESCE(proteins.bin_id, proteins.contig_id, '') AS genome_id
                FROM semantic_atoms AS atoms
                JOIN proteins ON proteins.protein_id = atoms.protein_id
                {join_filter}
                WHERE atoms.relation = 'unresolved'
            ), ranked AS (
                SELECT *, ROW_NUMBER() OVER (
                    PARTITION BY source_db, source_accession
                    ORDER BY protein_id
                ) AS example_rank
                FROM unresolved
            )
            SELECT
                source_db,
                source_accession,
                COUNT(DISTINCT protein_id) AS n_proteins,
                COUNT(DISTINCT NULLIF(genome_id, '')) AS n_genomes,
                STRING_AGG(protein_id, ';' ORDER BY protein_id)
                    FILTER (WHERE example_rank <= 5) AS examples
            FROM ranked
            GROUP BY source_db, source_accession
            ORDER BY n_proteins * GREATEST(n_genomes, 1) DESC,
                     source_db, source_accession
        """)
    finally:
        if registered_filter:
            store.conn.unregister("_v2_review_proteins")

    queue = []
    for source_db, accession, n_proteins, n_genomes, examples in rows:
        queue.append({
            "accession": accession,
            "source_db": source_db,
            "n_proteins": int(n_proteins),
            "n_genomes": int(n_genomes),
            "pct_proteome": round(
                (100.0 * int(n_proteins) / total_proteins)
                if total_proteins > 0
                else 0.0,
                2,
            ),
            "example_protein_ids": examples or "",
            "suggested_facet": suggest_facet(accession or ""),
        })
    return queue


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
    *,
    start_after: str | None = None,
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

    last_protein_id = start_after
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


def _clear_generation_tables(store: DuckDBStore) -> None:
    """Clear constraint-free full-refresh append tables."""
    for table_name in _GENERATION_TABLES:
        store.execute(f"DELETE FROM {table_name};")


def _clear_generation_tables_best_effort(store: DuckDBStore) -> None:
    """Release promoted build rows while preserving a completed generation."""
    try:
        for table_name in _GENERATION_TABLES:
            store.execute(f"DROP TABLE IF EXISTS {table_name};")
    except Exception:
        logger.exception("Failed to clear promoted V2 generation tables")


def _promote_full_generation(
    store: DuckDBStore,
    *,
    update_legacy_predicates: bool,
    expected_count: int,
) -> None:
    """Bulk-copy validated generation tables into canonical constrained tables."""
    generation_states = int(
        store.execute("SELECT COUNT(*) FROM v2_generation_state")[0][0]
    )
    if generation_states != expected_count:
        raise RuntimeError(
            "V2 generation-state count changed before canonical promotion"
        )
    if update_legacy_predicates:
        generation_legacy = int(
            store.execute("SELECT COUNT(*) FROM v2_generation_legacy")[0][0]
        )
        if generation_legacy != expected_count:
            raise RuntimeError(
                "V2 generation legacy-cache count changed before promotion"
            )

    started_at = time.monotonic()
    store.conn.execute("BEGIN TRANSACTION;")
    try:
        _clear_v2_tables(
            store,
            update_legacy_predicates=update_legacy_predicates,
        )
        store.execute("""
            INSERT INTO semantic_atoms
            SELECT * FROM v2_generation_atoms
        """)
        store.execute("""
            INSERT INTO semantic_state
            SELECT * FROM v2_generation_state
        """)
        store.execute("""
            INSERT INTO semantic_terms
            SELECT * FROM v2_generation_terms
        """)
        if update_legacy_predicates:
            store.execute("""
                INSERT INTO protein_predicates
                SELECT * FROM v2_generation_legacy
            """)
    except BaseException:
        store.conn.execute("ROLLBACK;")
        raise
    else:
        store.conn.execute("COMMIT;")

    canonical_states = int(
        store.execute("SELECT COUNT(*) FROM semantic_state")[0][0]
    )
    if canonical_states != expected_count:
        raise RuntimeError("Canonical V2 state count changed during promotion")
    logger.info(
        "V2 generation tables promoted in %.1fs",
        time.monotonic() - started_at,
    )


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


def _chunks(items: list[Any], chunk_size: int) -> Iterable[list[Any]]:
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
    use_generation_tables: bool = False,
    checkpoint_processed: int | None = None,
    checkpoint_last_protein_id: str | None = None,
) -> None:
    """Persist one generation chunk inside one DuckDB transaction."""
    if use_generation_tables:
        _persist_generation_chunk(
            store,
            atoms,
            states,
            legacy_rows=legacy_rows,
            update_legacy_predicates=update_legacy_predicates,
            checkpoint_processed=checkpoint_processed,
            checkpoint_last_protein_id=checkpoint_last_protein_id,
        )
        return

    store.conn.execute("BEGIN TRANSACTION;")
    try:
        _persist_atoms(store, atoms)
        _persist_states(store, states)
        _persist_semantic_terms(
            store,
            _build_semantic_terms(atoms, states),
        )
        if update_legacy_predicates:
            _persist_legacy_predicates(store, legacy_rows)
        if checkpoint_processed is not None:
            updated = store.execute(
                """
                UPDATE v2_generation_checkpoint
                SET last_protein_id = ?, processed_count = ?,
                    updated_at = CURRENT_TIMESTAMP
                WHERE generation_key = ?
                RETURNING generation_key
                """,
                [
                    checkpoint_last_protein_id,
                    checkpoint_processed,
                    _GENERATION_KEY,
                ],
            )
            if updated != [(_GENERATION_KEY,)]:
                raise RuntimeError("V2 checkpoint row disappeared during generation")
    except BaseException:
        store.conn.execute("ROLLBACK;")
        raise
    else:
        store.conn.execute("COMMIT;")


def _persist_generation_chunk(
    store: DuckDBStore,
    atoms: list[SemanticAtom],
    states: dict[str, SemanticState],
    *,
    legacy_rows: list[tuple[str, list[str]]],
    update_legacy_predicates: bool,
    checkpoint_processed: int | None,
    checkpoint_last_protein_id: str | None,
) -> None:
    """Append one full-refresh chunk through reusable vectorized batch views."""
    raw_atom_rows = [
        (
            position,
            atom.protein_id,
            atom.atom_id,
            atom.facet.value,
            atom.relation.value,
            atom.source_accession,
            atom.source_db,
            atom.evidence_evalue,
            atom.evidence_score,
        )
        for position, atom in enumerate(atoms)
    ]
    raw_atom_frame = pd.DataFrame.from_records(
        raw_atom_rows,
        columns=[
            "input_position",
            "protein_id",
            "atom_id",
            "facet",
            "relation",
            "source_accession",
            "source_db",
            "evidence_evalue",
            "evidence_score",
        ],
    )
    state_frame = pd.DataFrame.from_records(
        [
            (
                protein_id,
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
            for protein_id, state in states.items()
        ],
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
    )
    legacy_frame = pd.DataFrame.from_records(
        legacy_rows,
        columns=["protein_id", "predicates"],
    )

    registered_views = [
        ("_v2_raw_atoms_batch", raw_atom_frame),
        ("_v2_states_batch", state_frame),
    ]
    if update_legacy_predicates:
        registered_views.append(("_v2_legacy_batch", legacy_frame))
    for view_name, frame in registered_views:
        store.conn.register(view_name, frame)

    source_prefix_case = direct_access_prefix_case_sql("source_db")
    try:
        store.conn.execute("BEGIN TRANSACTION;")
        try:
            # The canonical atom key keeps the last occurrence. ROW_NUMBER
            # reproduces the prior ordered Python-dict behavior in DuckDB.
            store.execute("""
                INSERT INTO v2_generation_atoms
                SELECT
                    protein_id, atom_id, facet, relation, source_accession,
                    source_db, evidence_evalue, evidence_score
                FROM _v2_raw_atoms_batch
                QUALIFY ROW_NUMBER() OVER (
                    PARTITION BY protein_id, atom_id, source_accession
                    ORDER BY input_position DESC
                ) = 1
            """)
            store.execute("""
                INSERT INTO v2_generation_state
                SELECT
                    protein_id, activities, roles, architecture, localization,
                    topology, size_class, quality_flags, composite_predicates,
                    unresolved_count, CURRENT_TIMESTAMP
                FROM _v2_states_batch
            """)
            store.execute(f"""
                INSERT INTO v2_generation_terms
                SELECT DISTINCT
                    protein_id,
                    atom_id AS term_id,
                    'atom' AS term_kind,
                    facet,
                    relation,
                    COALESCE(source_db, '') AS source_db,
                    COALESCE(source_accession, '') AS source_accession
                FROM _v2_raw_atoms_batch
                WHERE NOT starts_with(atom_id, '_source_witness:')

                UNION ALL

                SELECT DISTINCT
                    protein_id,
                    source_prefix || ':' || source_accession AS term_id,
                    'direct_access' AS term_kind,
                    facet,
                    relation,
                    COALESCE(source_db, '') AS source_db,
                    source_accession
                FROM (
                    SELECT
                        protein_id,
                        {source_prefix_case} AS source_prefix,
                        facet,
                        relation,
                        source_db,
                        source_accession
                    FROM _v2_raw_atoms_batch
                    WHERE NOT starts_with(atom_id, '_source_witness:')
                      AND source_accession IS NOT NULL
                      AND NOT starts_with(source_accession, '_')
                ) direct_terms
                WHERE source_prefix IS NOT NULL

                UNION ALL

                SELECT DISTINCT
                    protein_id,
                    composite AS term_id,
                    'composite' AS term_kind,
                    NULL AS facet,
                    'implies' AS relation,
                    '_composite' AS source_db,
                    composite AS source_accession
                FROM (
                    SELECT
                        protein_id,
                        UNNEST(composite_predicates) AS composite
                    FROM _v2_states_batch
                    WHERE LEN(composite_predicates) > 0
                ) composites
            """)
            if update_legacy_predicates:
                store.execute("""
                    INSERT INTO v2_generation_legacy
                    SELECT protein_id, predicates, CURRENT_TIMESTAMP
                    FROM _v2_legacy_batch
                """)
            if checkpoint_processed is not None:
                updated = store.execute(
                    """
                    UPDATE v2_generation_checkpoint
                    SET last_protein_id = ?, processed_count = ?,
                        updated_at = CURRENT_TIMESTAMP
                    WHERE generation_key = ?
                    RETURNING generation_key
                    """,
                    [
                        checkpoint_last_protein_id,
                        checkpoint_processed,
                        _GENERATION_KEY,
                    ],
                )
                if updated != [(_GENERATION_KEY,)]:
                    raise RuntimeError(
                        "V2 checkpoint row disappeared during generation"
                    )
        except BaseException:
            store.conn.execute("ROLLBACK;")
            raise
        else:
            store.conn.execute("COMMIT;")
    finally:
        for view_name, _ in reversed(registered_views):
            store.conn.unregister(view_name)


def _persist_atoms(
    store: DuckDBStore,
    atoms: list[SemanticAtom],
) -> None:
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
