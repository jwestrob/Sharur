#!/usr/bin/env python3
"""
Stage 10 (optional): CDD per protein with RPS-BLAST ("resembles" rows).

Reads every protein of the Stage 07 database (read-only), writes each unique sequence once
(named by its MD5), runs RPS-BLAST against the CDD database of ``--config`` with every
allocated CPU (workers of four threads), and expands hits at E <= ``--max-evalue`` (default
0.01) to every protein carrying that sequence. The staging file keeps every hit and marks the
best model of each region (``--max-overlap``); activation writes those best rows to
``annotations`` by default (``--cdd-rows all`` writes every hit).

Outputs in ``stage10_cdd/``:

- ``remote_homology_cdd.duckdb``: staging DuckDB with ``cdd_hits`` (per unique sequence MD5:
  model id, E-value, bits, ranges), ``cdd_models`` and ``cdd_runs`` (database name and version,
  E-value, counts); activation expands the hits to every protein with that sequence;
- ``cdd/export.json``, ``cdd/compact/`` and ``cdd/out/chunk_*.tsv`` (raw RPS-BLAST tables);
- ``processing_manifest.json``.

The dataset database stays unchanged. Load the rows (``annotations`` with ``source =
'cdd'``) with ``sharur remote-homology activate --db DB --staging
stage10_cdd/remote_homology_cdd.duckdb --publish``. CDD rows are observed similarities: the
SQL surface labels them ``resembles`` and predicate generation leaves them out.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

from sharur.remote_homology import runner
from sharur.remote_homology.store import RemoteHomologyError, export_cdd, load_cdd


def main() -> int:
    parser = argparse.ArgumentParser(description="CDD per protein with RPS-BLAST (remote homology).")
    parser.add_argument("--db", type=Path, required=True, help="Stage 07 sharur.duckdb (read-only)")
    parser.add_argument("--out-dir", type=Path, required=True, help="Stage directory (stage10_cdd)")
    parser.add_argument("--config", type=Path, required=True, help="Remote-homology configuration (JSON)")
    parser.add_argument("--threads", type=int, default=None, help="Default: $SLURM_CPUS_ON_NODE or all CPUs")
    parser.add_argument("--max-evalue", type=float, default=0.01)
    parser.add_argument("--max-overlap", type=float, default=0.5,
                        help="Best model per region: overlap fraction of the shorter range that makes one region")
    parser.add_argument("--chunk-size", type=int, default=2000)
    parser.add_argument("--keep-chunks", action="store_true", help="Keep the exported sequence chunks")
    parser.add_argument("--force", action="store_true", help="Replace an existing stage directory")
    args = parser.parse_args()

    if not args.db.is_file():
        raise SystemExit(f"Database not found: {args.db}")
    manifest_path = args.out_dir / "processing_manifest.json"
    if manifest_path.exists() and not args.force:
        raise SystemExit(f"{args.out_dir} already holds a finished run; pass --force to replace it")
    if args.force and args.out_dir.exists():
        shutil.rmtree(args.out_dir)
    threads = args.threads or runner.allocated_cpus()
    try:
        runner.check_database(runner.load_config(args.config), "cdd")
        export = export_cdd(args.db, args.out_dir, chunk_size=args.chunk_size, threads=min(threads, 8))
        code = runner.main(["cdd", "--work", str(args.out_dir), "--config", str(args.config),
                            "--threads", str(threads)])
        if code:
            print(f"Stage 10 failed: RPS-BLAST reported errors (see {args.out_dir}/logs)", file=sys.stderr)
            return 1
        staging = args.out_dir / "remote_homology_cdd.duckdb"
        result = load_cdd(args.out_dir, staging, args.config, max_evalue=args.max_evalue,
                          max_overlap=args.max_overlap, threads=min(threads, 8))
    except (RemoteHomologyError, runner.ConfigError, FileNotFoundError) as exc:
        print(f"Stage 10 failed: {exc}", file=sys.stderr)
        return 1
    if not args.keep_chunks:
        shutil.rmtree(args.out_dir / "cdd" / "chunks", ignore_errors=True)
    manifest = {"stage": "stage10_cdd", "db": str(args.db.resolve()), "export": export, **result}
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"CDD: {result['hit_rows']:,} hits on {result['sequences_with_hits']:,} of {result['unique_sequences']:,} "
          f"unique sequences ({result['proteins']:,} proteins)")
    print(f"Activate: sharur remote-homology activate --db {args.db} --staging {staging} --publish")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
