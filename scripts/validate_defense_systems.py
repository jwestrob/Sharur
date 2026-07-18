#!/usr/bin/env python3
"""Validate DefenseFinder hits with Sharur's in-process co-location caller.

This module remains as a compatibility entry point for older ingest scripts.
The authoritative implementation is :mod:`sharur.colocation`, which reads
protein positions and HMM hits directly from DuckDB and scopes every call to a
single ``(genome_id, contig_id)`` replicon.

Usage:
    python scripts/validate_defense_systems.py data/susan_genomes
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import TYPE_CHECKING


if TYPE_CHECKING:
    import pandas as pd

from sharur.colocation import integrate_defense_results, validate_systems


logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def validate_defense_systems(
    db_path: str | Path,
    data_dir: str | Path | None = None,
    workers: int = 4,
    verbose: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return curated DefenseFinder systems and their member assignments.

    ``workers`` is retained for API compatibility. The in-process caller is
    position-indexed and does not launch MacSyFinder subprocesses.
    """
    del workers
    systems_df, genes_df = validate_systems(
        db_path=db_path,
        source="defensefinder",
        model_family="defense-finder-models",
        verbose=verbose,
    )

    if data_dir is not None:
        results_dir = Path(data_dir) / "stage04_astra" / "defensefinder_results"
        results_dir.mkdir(parents=True, exist_ok=True)
        systems_df.to_csv(
            results_dir / "validated_defense_systems.tsv",
            sep="\t",
            index=False,
        )
        genes_df.to_csv(
            results_dir / "validated_defense_genes.tsv",
            sep="\t",
            index=False,
        )

    return systems_df, genes_df


def integrate_results(
    db_path: str | Path,
    all_systems: pd.DataFrame,
    all_genes: pd.DataFrame,
) -> set[str]:
    """Compatibility wrapper for replacement-based database integration."""
    return integrate_defense_results(db_path, all_systems, all_genes)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Validate DefenseFinder HMM hits by replicon-local co-location"
    )
    parser.add_argument(
        "data_dir",
        type=Path,
        help="Dataset directory containing sharur.duckdb",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Retained for compatibility; the in-process caller does not use it",
    )
    parser.add_argument(
        "--no-integrate",
        action="store_true",
        help="Write validation TSVs without replacing database caller tables",
    )
    args = parser.parse_args()

    db_path = args.data_dir / "sharur.duckdb"
    systems_df, genes_df = validate_defense_systems(
        db_path=db_path,
        data_dir=args.data_dir,
        workers=args.workers,
    )

    if not args.no_integrate and db_path.exists():
        integrate_results(db_path, systems_df, genes_df)

    logger.info(
        "Done: %d systems and %d member assignments",
        len(systems_df),
        len(genes_df),
    )


if __name__ == "__main__":
    main()
