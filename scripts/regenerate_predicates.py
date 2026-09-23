#!/usr/bin/env python3
"""Regenerate V2 predicates and the legacy compatibility cache.

Usage:
    python scripts/regenerate_predicates.py --db data/DATASET/sharur.duckdb
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sharur.predicates_v2.persistence import generate_and_persist_v2
from sharur.storage.duckdb_store import DuckDBStore


def _count_rows(store: DuckDBStore, table: str) -> int:
    try:
        return store.execute(f"SELECT COUNT(*) FROM {table}")[0][0]
    except Exception:
        return 0


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Regenerate V2 predicates for a Sharur database",
    )
    parser.add_argument("--db", required=True, help="Path to sharur.duckdb")
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=100_000,
        help="Number of proteins per V2 generation batch",
    )
    parser.add_argument(
        "--protein",
        help="Optional single protein_id to refresh",
    )
    parser.add_argument(
        "--review-queue",
        help="Optional TSV path for unresolved-accession review queue",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help="Deprecated compatibility option; V2 generation is database-chunked.",
    )
    args = parser.parse_args()

    store = DuckDBStore(Path(args.db))
    protein_ids = [args.protein] if args.protein else None

    if args.workers is not None:
        print("--workers is ignored; V2 generation uses chunked DuckDB batches.")

    before = {
        "protein_predicates": _count_rows(store, "protein_predicates"),
        "semantic_atoms": _count_rows(store, "semantic_atoms"),
        "semantic_state": _count_rows(store, "semantic_state"),
    }
    print("Current rows:")
    for table, count in before.items():
        print(f"  {table}: {count:,}")

    generate_and_persist_v2(
        store,
        protein_ids=protein_ids,
        output_review_queue=args.review_queue,
        chunk_size=args.chunk_size,
        update_legacy_predicates=True,
        return_states=False,
        predict_topology=False,
    )

    after = {
        "protein_predicates": store.execute("SELECT COUNT(*) FROM protein_predicates")[0][0],
        "semantic_atoms": store.execute("SELECT COUNT(*) FROM semantic_atoms")[0][0],
        "semantic_state": store.execute("SELECT COUNT(*) FROM semantic_state")[0][0],
        "unique_atoms": store.execute("SELECT COUNT(DISTINCT atom_id) FROM semantic_atoms")[0][0],
    }

    print("\nRegenerated V2 predicates:")
    for table, count in after.items():
        print(f"  {table}: {count:,}")

    if args.review_queue:
        print(f"  review_queue: {args.review_queue}")


if __name__ == "__main__":
    main()
