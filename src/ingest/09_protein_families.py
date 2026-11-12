#!/usr/bin/env python3
"""
Stage 9 (optional): annotation-independent protein families with MMseqs2.

Reads every protein of the Stage 07 database (read-only) and clusters it twice,
hierarchically: subfamilies at ``--subfamily-min-seq-id`` over all proteins, then
families at ``--family-min-seq-id`` over the subfamily representatives, so each
subfamily lies inside one family. Both levels require coverage of both sequences
(``-c 0.8 --cov-mode 0``), sensitivity 6 and ``--cluster-reassign``.

Outputs in ``stage09_families/``:

- ``protein_families.duckdb``: staging DuckDB with ``protein_families``,
  ``family_stats``, ``protein_family_runs`` (a family run and a subfamily run)
  and the per-genome presence views;
- ``clustering/`` (``level0.tsv.gz`` subfamilies, ``level1.tsv.gz`` families,
  ``clustering.json``, ``export.json``) for provenance;
- ``processing_manifest.json``.

With ``--rescue`` (the default follows ``store.RESCUE_BY_DEFAULT``) the staging file also
gets the truncated-gene rescue layer: gene calls that run off a contig end and ended up
family singletons are placed into existing families when an alignment covers at least
``--rescue-min-query-coverage`` of them at ``--rescue-min-seq-id`` identity. The clustering
runs stay unchanged; the layer is a pair of derived runs (``docs/guides/protein-families.md``).

The dataset database is left unchanged, which keeps the ingest ledger's Stage 07
snapshot valid. Activate the families with
``sharur families activate --db DB --staging stage09_families/protein_families.duckdb --publish``.

MMseqs2 and DuckDB use every allocated CPU (``$SLURM_CPUS_ON_NODE`` when set).
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

from sharur.protein_families import mmseqs
from sharur.protein_families.store import RESCUE_BY_DEFAULT, FamilyError, RescueParams, run_pipeline


def main() -> int:
    parser = argparse.ArgumentParser(description="Annotation-independent protein families (MMseqs2).")
    parser.add_argument("--db", type=Path, required=True, help="Stage 07 sharur.duckdb (read-only)")
    parser.add_argument("--out-dir", type=Path, required=True, help="Stage directory (stage09_families)")
    parser.add_argument("--threads", type=int, default=None,
                        help="Default: $SLURM_CPUS_ON_NODE or the host CPU count")
    parser.add_argument("--subfamily-min-seq-id", type=float, default=mmseqs.SUBFAMILY_MIN_SEQ_ID)
    parser.add_argument("--family-min-seq-id", type=float, default=mmseqs.FAMILY_MIN_SEQ_ID)
    parser.add_argument("--coverage", type=float, default=mmseqs.COVERAGE)
    parser.add_argument("--sensitivity", type=float, default=mmseqs.SENSITIVITY)
    parser.add_argument("--memory-limit", default=None, help="DuckDB memory limit (default: 25%% of the allocation)")
    parser.add_argument("--split-memory-limit", default=None, help="MMseqs2 --split-memory-limit")
    parser.add_argument("--mmseqs", default=None, help="MMseqs2 executable (default: $SHARUR_MMSEQS or mmseqs)")
    parser.add_argument("--force", action="store_true", help="Replace an existing stage directory")
    parser.add_argument("--rescue", action=argparse.BooleanOptionalAction, default=RESCUE_BY_DEFAULT,
                        help="Add the truncated-gene rescue layer")
    parser.add_argument("--rescue-min-query-coverage", type=float, default=RescueParams.min_query_coverage)
    parser.add_argument("--rescue-min-seq-id", type=float, default=None,
                        help="Default: the family identity floor")
    parser.add_argument("--rescue-min-family-size", type=int, default=RescueParams.min_family_size)
    parser.add_argument("--rescue-include", action=argparse.BooleanOptionalAction, default=RESCUE_BY_DEFAULT,
                        help="Make the rescued runs the defaults (views include rescued members)")
    args = parser.parse_args()

    if not args.db.is_file():
        raise SystemExit(f"Database not found: {args.db}")
    manifest_path = args.out_dir / "processing_manifest.json"
    if manifest_path.exists() and not args.force:
        raise SystemExit(f"{args.out_dir} already holds a finished run; pass --force to replace it")
    if args.force and args.out_dir.exists():
        shutil.rmtree(args.out_dir)
    try:
        levels = mmseqs.default_levels(args.subfamily_min_seq_id, args.family_min_seq_id,
                                       coverage=args.coverage, sensitivity=args.sensitivity)
        rescue = (RescueParams(min_query_coverage=args.rescue_min_query_coverage,
                               min_seq_id=args.rescue_min_seq_id, min_family_size=args.rescue_min_family_size)
                  if args.rescue else None)
        result = run_pipeline(args.db, args.out_dir, levels=levels, threads=args.threads,
                              memory_limit=args.memory_limit, split_memory_limit=args.split_memory_limit,
                              binary=args.mmseqs, rescue_params=rescue, rescue_include=args.rescue_include)
    except (FamilyError, ValueError, FileNotFoundError) as exc:
        print(f"Stage 09 failed: {exc}", file=sys.stderr)
        return 1
    manifest = {"stage": "stage09_families", "db": str(args.db.resolve()), **result}
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Protein families: {result['families']:,} families and {result['subfamilies']:,} subfamilies over "
          f"{result['proteins']:,} proteins (build {result['build_id']})")
    print(f"Activate: sharur families activate --db {args.db} --staging {result['staging']} --publish")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
