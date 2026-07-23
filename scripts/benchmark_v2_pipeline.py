#!/usr/bin/env python3
"""Benchmark the bounded Stage 07 V2 pipeline on a synthetic DuckDB corpus."""

from __future__ import annotations

import argparse
import logging
import shutil
import tempfile
import time
from pathlib import Path

import sharur.predicates_v2.persistence as persistence
from sharur.predicates_v2.persistence import generate_and_persist_v2
from sharur.storage.duckdb_store import DuckDBStore


def _seed_database(path: Path, protein_count: int, *, threads: int) -> None:
    with DuckDBStore(path, threads=threads) as store:
        store.conn.execute(
            """
            INSERT INTO bins (
                bin_id, completeness, contamination, taxonomy,
                n_contigs, total_length
            )
            VALUES ('benchmark_bin', 100.0, 0.0, 'synthetic', 1, ?)
            """,
            [protein_count * 300],
        )
        store.conn.execute(
            """
            INSERT INTO contigs (contig_id, bin_id, length, gc_content)
            VALUES ('benchmark_contig', 'benchmark_bin', ?, 0.5)
            """,
            [protein_count * 300],
        )
        store.conn.execute(
            """
            INSERT INTO proteins (
                protein_id, contig_id, bin_id, start, end_coord, strand,
                gene_index, sequence_length, gc_content
            )
            SELECT
                printf('protein_%09d', protein_index),
                'benchmark_contig',
                'benchmark_bin',
                protein_index * 300 + 1,
                protein_index * 300 + 300,
                '+',
                protein_index,
                100,
                0.5
            FROM range(?) AS source(protein_index)
            """,
            [protein_count],
        )
        store.conn.execute(
            """
            INSERT INTO annotations (
                annotation_id, protein_id, source, accession, name,
                description, evalue, score
            )
            SELECT
                protein_index + 1,
                printf('protein_%09d', protein_index),
                'pfam',
                'PF00005',
                'ABC_tran',
                'synthetic benchmark annotation',
                1e-30,
                100.0
            FROM range(?) AS source(protein_index)
            """,
            [protein_count],
        )


def _run_once(
    seed_path: Path,
    run_path: Path,
    *,
    workers: int,
    chunk_size: int,
    worker_batch_size: int | None,
    pipeline_depth: int,
    object_worker_output: bool,
) -> dict[str, float | int]:
    shutil.copy2(seed_path, run_path)
    started_at = time.monotonic()
    original_frame_worker = persistence._transform_batch_frames_worker
    if object_worker_output:
        persistence._transform_batch_frames_worker = (
            persistence._transform_batch_worker
        )
    try:
        with DuckDBStore(run_path, threads=workers) as store:
            generate_and_persist_v2(
                store,
                chunk_size=chunk_size,
                workers=workers,
                worker_batch_size=worker_batch_size,
                pipeline_depth=pipeline_depth,
                update_legacy_predicates=True,
                return_states=False,
            )
            elapsed = time.monotonic() - started_at
            states = int(
                store.execute("SELECT COUNT(*) FROM semantic_state")[0][0]
            )
            atoms = int(store.execute("SELECT COUNT(*) FROM semantic_atoms")[0][0])
    finally:
        persistence._transform_batch_frames_worker = original_frame_worker
    return {
        "pipeline_depth": pipeline_depth,
        "seconds": elapsed,
        "proteins_per_second": states / elapsed,
        "states": states,
        "atoms": atoms,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--proteins", type=int, default=100_000)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--chunk-size", type=int, default=10_000)
    parser.add_argument("--worker-batch-size", type=int)
    parser.add_argument(
        "--pipeline-depth",
        type=int,
        nargs="+",
        default=[1, 2, 3],
        dest="pipeline_depths",
    )
    parser.add_argument("--work-dir", type=Path)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument(
        "--object-worker-output",
        action="store_true",
        help="Diagnostic baseline using Python object graphs across worker pipes",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
    )

    if args.proteins < 1:
        parser.error("--proteins must be positive")
    if args.workers < 1:
        parser.error("--workers must be positive")
    if args.chunk_size < 1:
        parser.error("--chunk-size must be positive")
    if args.worker_batch_size is not None and args.worker_batch_size < 1:
        parser.error("--worker-batch-size must be positive")
    if any(depth < 1 for depth in args.pipeline_depths):
        parser.error("--pipeline-depth values must be positive")

    temporary = None
    if args.work_dir is None:
        temporary = tempfile.TemporaryDirectory(prefix="sharur-v2-benchmark-")
        work_dir = Path(temporary.name)
    else:
        work_dir = args.work_dir.expanduser().resolve()
        work_dir.mkdir(parents=True, exist_ok=True)

    try:
        seed_path = work_dir / "seed.duckdb"
        _seed_database(seed_path, args.proteins, threads=args.workers)

        print(
            "depth\tseconds\tproteins_per_second\tstates\tatoms",
            flush=True,
        )
        for run_index, depth in enumerate(args.pipeline_depths, start=1):
            result = _run_once(
                seed_path,
                work_dir / f"run_{run_index}_depth_{depth}.duckdb",
                workers=args.workers,
                chunk_size=args.chunk_size,
                worker_batch_size=args.worker_batch_size,
                pipeline_depth=depth,
                object_worker_output=args.object_worker_output,
            )
            print(
                f"{result['pipeline_depth']}\t"
                f"{result['seconds']:.3f}\t"
                f"{result['proteins_per_second']:.1f}\t"
                f"{result['states']}\t"
                f"{result['atoms']}",
                flush=True,
            )
    finally:
        if temporary is not None:
            temporary.cleanup()


if __name__ == "__main__":
    main()
