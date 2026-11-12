#!/usr/bin/env python3
"""
Stage 11 (optional): HMM-HMM remote homology per protein family ("resembles" rows).

Needs a protein-family clustering (Stage 09's staging file, or a run imported into it with
``sharur families import``). Selects the families of ``--run`` (default: the default family
run) found in ``--min-genomes`` or more genomes whose members carry no Pfam hit, writes up to
``--max-sequences`` unique member sequences per family (representative first), and for each
family runs ``mafft --auto``, optionally HHblits enrichment against UniRef30
(``msa_kinds: ["enriched", ...]`` in the configuration), and HHsearch against the Pfam and
PDB profile databases of ``--config``, with every allocated CPU (workers of four threads).

Outputs in ``stage11_family_homology/``:

- ``remote_homology_families.duckdb``: staging DuckDB with ``family_annotations`` (hits at
  probability >= 50%, with database name and version), ``family_annotation_status`` and
  ``family_annotation_runs``, and the ``family_resemblance`` view (the "resembles" rule:
  probability >= 90% over >= 30 aligned columns);
- ``out/`` (gzipped ``.hhr`` files, ``hits.tsv`` and ``done.json`` per family);
- ``processing_manifest.json``.

A rerun with the same configuration skips finished families. The dataset database stays
unchanged; load with ``sharur remote-homology activate --db DB --staging
stage11_family_homology/remote_homology_families.duckdb --families-staging
stage09_families/protein_families.duckdb --publish``.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

from sharur.remote_homology import runner
from sharur.remote_homology.store import (
    MAX_SEQUENCES,
    MIN_GENOMES,
    RemoteHomologyError,
    export_families,
    load_families,
)


def main() -> int:
    parser = argparse.ArgumentParser(description="HMM-HMM remote homology per protein family.")
    parser.add_argument("--db", type=Path, required=True, help="Stage 07 sharur.duckdb (read-only)")
    parser.add_argument("--families", type=Path, required=True,
                        help="Protein-family staging DuckDB (stage09_families/protein_families.duckdb)")
    parser.add_argument("--run", default=None, help="Clustering run (default: the default family run)")
    parser.add_argument("--out-dir", type=Path, required=True, help="Stage directory (stage11_family_homology)")
    parser.add_argument("--config", type=Path, required=True, help="Remote-homology configuration (JSON)")
    parser.add_argument("--threads", type=int, default=None, help="Default: $SLURM_CPUS_ON_NODE or all CPUs")
    parser.add_argument("--min-genomes", type=int, default=MIN_GENOMES)
    parser.add_argument("--include-pfam-annotated", action="store_true",
                        help="Also search families whose members carry Pfam hits")
    parser.add_argument("--max-sequences", type=int, default=MAX_SEQUENCES)
    parser.add_argument("--force", action="store_true", help="Replace an existing stage directory")
    args = parser.parse_args()

    if not args.db.is_file():
        raise SystemExit(f"Database not found: {args.db}")
    if not args.families.is_file():
        raise SystemExit(f"Protein-family staging file not found: {args.families}")
    manifest_path = args.out_dir / "processing_manifest.json"
    if manifest_path.exists() and not args.force:
        raise SystemExit(f"{args.out_dir} already holds a finished run; pass --force to replace it")
    if args.force and args.out_dir.exists():
        shutil.rmtree(args.out_dir)
    threads = args.threads or runner.allocated_cpus()
    staging = args.out_dir / "remote_homology_families.duckdb"
    try:
        cfg = runner.load_config(args.config)
        for key in cfg["family_targets"] + (["uniref"] if "enriched" in cfg["msa_kinds"] else []):
            runner.check_database(cfg, key)
        export = export_families(args.db, args.out_dir, families_db=args.families, run=args.run,
                                 min_genomes=args.min_genomes, require_pfam_dark=not args.include_pfam_annotated,
                                 max_sequences=args.max_sequences, threads=min(threads, 8))
        code = runner.main(["families", "--work", str(args.out_dir), "--config", str(args.config),
                            "--threads", str(threads)])
        if code:
            print(f"Stage 11 failed: some families failed (see {args.out_dir}/out/*/*/error.json); "
                  "rerun with the same arguments minus --force to retry them", file=sys.stderr)
            return 1
        result = load_families(args.out_dir, staging, args.config, threads=min(threads, 8))
    except (RemoteHomologyError, runner.ConfigError, FileNotFoundError) as exc:
        print(f"Stage 11 failed: {exc}", file=sys.stderr)
        return 1
    shutil.rmtree(args.out_dir / "fasta", ignore_errors=True)
    manifest = {"stage": "stage11_family_homology", "db": str(args.db.resolve()),
                "families_db": str(args.families.resolve()), "export": export, **result}
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Family tier: {result['ok']:,} families searched, {result['hit_rows']:,} hit rows (run "
          f"{result['clustering_run']})")
    print(f"Activate: sharur remote-homology activate --db {args.db} --staging {staging} "
          f"--families-staging {args.families} --publish")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
