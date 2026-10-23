#!/usr/bin/env python3
"""
Reload DefenseFinder annotations and systems from a per-domain hits table.

Rebuilds the `defensefinder` annotation rows of an existing database from an
Aksha/Astra hits table (sequence_id, hmm_name, bitscore, evalue, c_evalue,
i_evalue, env_from, env_to, dom_bitscore, ali_from, ali_to, hmm_from, hmm_to)
with MacSyFinder/DefenseFinder hit selection: GA bit thresholds where the
profile has them, then i-Evalue <= 1e-3 and profile coverage >= 0.4
(DefenseFinder 2.0.1's --coverage default; MacSyFinder's own default is 0.5).
Profile lengths and GA values come from the defense-finder-models directory
given with --models, whose version is recorded in `annotation_source_versions`.
It then reruns the co-location caller (`defense_systems`, `defensefinder_system`
annotations, `system_proteins`) and regenerates V2 semantics and the V1
compatibility cache for every protein whose DefenseFinder evidence changes.

The command stages a copy of the database, validates it, and reports the
difference: annotation rows, systems per model and per model family (the
definitions folder, comparable across releases that rename models), genomes
with any system, and system and profile families gained or lost. Without --publish the copy
is discarded. With --publish it atomically replaces the database, keeps the
prior file as a backup beside it, and re-seals when a dataset seal exists.

Usage:
    # dry run; the staged copy goes to --work-dir
    python scripts/reload_defensefinder.py --db data/DATASET/sharur.duckdb \\
        --hits DefenseFinder_hits_df.tsv --models defense-finder-models \\
        --work-dir /tmp/df_reload --report-dir reports/
    # publish (staging beside the database, for an atomic rename)
    python scripts/reload_defensefinder.py --db data/DATASET/sharur.duckdb \\
        --hits DefenseFinder_hits_df.tsv --models defense-finder-models --publish
"""

import argparse
import logging
import os
import sys
from pathlib import Path


sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sharur.defensefinder.refresh import reload_defensefinder
from sharur.defensefinder.selection import (
    DEFAULT_COVERAGE_PROFILE,
    DEFAULT_I_EVALUE_SEL,
    MACSYFINDER_E_VALUE_SEARCH,
    SelectionParams,
)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Reload DefenseFinder annotations and systems from a hits table",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--db", required=True, type=Path, help="Path to sharur.duckdb")
    parser.add_argument("--hits", required=True, type=Path,
                        help="Per-domain DefenseFinder hits table (Aksha/Astra *_hits_df.tsv)")
    parser.add_argument("--models", required=True, type=Path,
                        help="defense-finder-models directory (metadata.yml, definitions/, profiles/) "
                             "used for the search")
    parser.add_argument("--previous-models", type=Path,
                        help="defense-finder-models directory that produced the current systems; "
                             "maps their model names to families in the report")
    parser.add_argument("--publish", action="store_true",
                        help="Replace the database after validation (default: dry run)")
    parser.add_argument("--work-dir", type=Path,
                        help="Directory for the staged copy (default: beside the database; "
                             "publishing requires the same filesystem)")
    parser.add_argument("--report-dir", type=Path, help="Write summary, per-type and per-profile tables here")
    parser.add_argument("--label", help="Report file prefix (default: the dataset directory name)")
    parser.add_argument("--i-evalue", type=float, default=DEFAULT_I_EVALUE_SEL,
                        help=f"Maximum domain i-Evalue (default {DEFAULT_I_EVALUE_SEL:g})")
    parser.add_argument("--coverage", type=float, default=DEFAULT_COVERAGE_PROFILE,
                        help=f"Minimum profile coverage (default {DEFAULT_COVERAGE_PROFILE:g})")
    parser.add_argument("--threads", type=int, default=2, help="DuckDB threads (default 2)")
    parser.add_argument("--memory-limit", default="4GB", help="DuckDB memory limit (default 4GB)")
    parser.add_argument("--no-reseal", action="store_true",
                        help="Leave an existing dataset seal untouched after publishing")
    args = parser.parse_args()
    logging.basicConfig(level=logging.WARNING, format="%(levelname)s: %(message)s")

    params = SelectionParams(e_value_search=MACSYFINDER_E_VALUE_SEARCH,
                             i_evalue_sel=args.i_evalue, coverage_profile=args.coverage)
    report = reload_defensefinder(
        args.db,
        args.hits,
        args.models,
        dry_run=not args.publish,
        work_dir=args.work_dir,
        params=params,
        threads=args.threads,
        memory_limit=args.memory_limit,
        reseal=not args.no_reseal,
        report_dir=args.report_dir,
        label=args.label,
        previous_models=args.previous_models,
    )
    print(report.to_text())
    return 0


if __name__ == "__main__":
    sys.exit(main())
