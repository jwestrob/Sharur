#!/usr/bin/env python3
"""Load additional Astra annotation TSVs into an existing Bennu DuckDB.

Usage:
    python scripts/load_additional_annotations.py \
        --db data/omni_production/bennu.duckdb \
        --tsv data/omni_production/stage04_astra/hyddb_results/HydDB_hits_df.tsv \
        --source hyddb

    python scripts/load_additional_annotations.py \
        --db data/omni_production/bennu.duckdb \
        --tsv data/omni_production/stage04_astra/defensefinder_results/DefenseFinder_hits_df.tsv \
        --source defensefinder
"""

import argparse
import duckdb
import pandas as pd


def load_annotations(db_path: str, tsv_path: str, source: str) -> None:
    conn = duckdb.connect(db_path)

    source = source.lower()  # Enforce lowercase source names
    print(f"Loading {tsv_path} as source={source}...")
    df = pd.read_csv(tsv_path, sep="\t")
    print(f"  Raw rows: {len(df):,}")

    # Normalize column names to match Bennu schema
    df = df.rename(columns={
        "sequence_id": "protein_id",
        "hmm_name": "accession",
        "bitscore": "score",
        "dom_bitscore": "dom_score",
        "env_from": "start_aa",
        "env_to": "end_aa",
    })

    # For DefenseFinder, split system__component into name
    if source == "defensefinder" and "accession" in df.columns:
        df["name"] = df["accession"].apply(
            lambda x: x.split("__")[-1] if isinstance(x, str) and "__" in x else x
        )
        df["description"] = df["accession"].apply(
            lambda x: x.split("__")[0] if isinstance(x, str) and "__" in x else ""
        )
    else:
        if "name" not in df.columns:
            df["name"] = df["accession"]
        if "description" not in df.columns:
            df["description"] = ""

    df["source"] = source

    # Only keep proteins that exist in the DB
    existing = set(r[0] for r in conn.execute(
        "SELECT DISTINCT protein_id FROM proteins"
    ).fetchall())
    before = len(df)
    df = df[df["protein_id"].isin(existing)]
    print(f"  Matched to existing proteins: {len(df):,} ({before - len(df):,} dropped)")

    # Get next annotation_id
    next_id = conn.execute(
        "SELECT COALESCE(MAX(annotation_id), 0) + 1 FROM annotations"
    ).fetchone()[0]

    # Build final DataFrame
    keep = ["protein_id", "source", "accession", "name", "description",
            "evalue", "score", "start_aa", "end_aa"]
    df = df.reindex(columns=keep, fill_value=None)
    df["annotation_id"] = range(next_id, next_id + len(df))

    # Bulk insert via DataFrame register
    conn.register("tmp_new_annots", df)
    conn.execute("""
        INSERT INTO annotations (
            annotation_id, protein_id, source, accession, name,
            description, evalue, score, start_aa, end_aa
        )
        SELECT annotation_id, protein_id, source, accession, name,
               description, evalue, score, start_aa, end_aa
        FROM tmp_new_annots
    """)
    conn.unregister("tmp_new_annots")

    final_count = conn.execute(
        f"SELECT COUNT(*) FROM annotations WHERE source = '{source}'"
    ).fetchone()[0]
    print(f"  Inserted {len(df):,} annotations (total {source} in DB: {final_count:,})")

    conn.close()
    print("Done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", required=True)
    parser.add_argument("--tsv", required=True)
    parser.add_argument("--source", required=True)
    args = parser.parse_args()
    load_annotations(args.db, args.tsv, args.source)
