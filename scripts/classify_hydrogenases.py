#!/usr/bin/env python3
"""
Hydrogenase subgroup assignment against HydDB references.

Proteins with a HydDB HMM hit (Stage 04, Astra) receive the subgroup of their
best DIAMOND match among the installed HydDB reference sequences. This is a
Sharur nearest-reference assignment: the published HydDB classifier also votes
over k=4 neighbors, screens homologous non-hydrogenase families, and uses the
downstream gene to separate FeFe Group A subtypes. Subgroup interpretations
follow Søndergaard et al. (2016) Table 1; see sharur/hydrogenase/subgroups.py.

Each protein gets one reconciled record in `hydrogenase_classifications`
(discovery HMM class, nearest reference, identity/e-value/bit score, Pfam domain
observations, curation reason, reference checksum, classifier version) and its
derived labels as `hyddb_subgroup` annotations. Raw `hyddb` rows are unchanged.

The command stages a copy of the database, reclassifies, regenerates V2 and the
V1 compatibility cache for affected proteins, and validates the copy. Without
--publish it reports what would change and discards the copy. With --publish it
atomically replaces the database, keeps the prior file as a backup beside it,
and re-seals when a dataset seal exists.

Usage:
    python scripts/classify_hydrogenases.py --db data/DATASET/sharur.duckdb
    python scripts/classify_hydrogenases.py --db data/DATASET/sharur.duckdb --publish

References:
    Søndergaard D, Pedersen CNS, Greening C (2016) HydDB: A web tool for hydrogenase
    classification and analysis. Sci Rep 6:34212. doi:10.1038/srep34212
"""

import argparse
import os
import sys
from pathlib import Path


sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sharur.hydrogenase.refresh import refresh_hydrogenases


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Refresh HydDB hydrogenase subgroup assignments",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--db", required=True, help="Path to sharur.duckdb")
    parser.add_argument("--threads", type=int, default=os.cpu_count() or 4,
                        help="DIAMOND threads (default: all CPUs)")
    parser.add_argument("--publish", action="store_true",
                        help="Replace the database after validation (default: dry run)")
    parser.add_argument("--transitions", type=Path,
                        help="Write per-protein label/curation transitions to this TSV")
    parser.add_argument("--reference-dir", type=Path, help="HydDB reference directory")
    parser.add_argument("--no-reseal", action="store_true",
                        help="Leave an existing dataset seal untouched after publishing")
    args = parser.parse_args()

    report = refresh_hydrogenases(
        args.db,
        threads=args.threads,
        dry_run=not args.publish,
        reference_dir=args.reference_dir,
        transitions_path=args.transitions,
        reseal=not args.no_reseal,
    )
    print(report.to_text())
    return 0


if __name__ == "__main__":
    sys.exit(main())
