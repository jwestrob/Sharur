#!/usr/bin/env python3
"""Load additional Astra annotation TSVs into an existing Sharur DuckDB.

Usage:
    python scripts/load_additional_annotations.py \
        --db data/omni_production/sharur.duckdb \
        --tsv data/omni_production/stage04_astra/hyddb_results/HydDB_hits_df.tsv \
        --source hyddb

    python scripts/load_additional_annotations.py \
        --db data/omni_production/sharur.duckdb \
        --tsv data/omni_production/vogdb_results/VOGdb_hits_df.tsv \
        --source vogdb --evalue 1e-10
"""

import argparse
import re

import duckdb
import pandas as pd

# Default e-value thresholds per source (same as 07_build_knowledge_base.py)
DEFAULT_EVALUE_THRESHOLDS = {
    "pfam": 1e-5,
    "kegg": 1e-5,
    "hyddb": 1e-5,
    "defensefinder": 1e-15,
    "vogdb": 1e-15,
    "cant_hyd": 1e-10,
    "cazy": 1e-5,
}

# Regex for DefenseFinder accessions with numbered FAM variants.
# Matches trailing _FAM_N or _FAMN at end of accession string.
# See 07_build_knowledge_base.py for full explanation.
_DEFENSEFINDER_FAM_RE = re.compile(r"_FAM_?\d+$")


def _dedup_defensefinder_fam_variants(df: pd.DataFrame) -> pd.DataFrame:
    """Deduplicate DefenseFinder hits that differ only by FAM variant number.

    DefenseFinder has ~40+ Type_II_MTases_FAM_* HMM profiles that all detect
    variations of the same Rossmann methyltransferase fold.  When a single
    protein is scanned, it matches many FAMs simultaneously, inflating hit
    counts.  This keeps only the best-scoring hit per (protein_id, component)
    when multiple _FAM_N profiles match.
    """
    if df.empty:
        return df
    df = df.copy()
    df["_dedup_group"] = df["accession"].apply(
        lambda x: _DEFENSEFINDER_FAM_RE.sub("", x) if isinstance(x, str) else x
    )
    has_fam = df["_dedup_group"] != df["accession"]
    if not has_fam.any():
        df.drop(columns=["_dedup_group"], inplace=True)
        return df
    no_fam_df = df[~has_fam].drop(columns=["_dedup_group"])
    fam_df = df[has_fam].copy()
    fam_df = fam_df.sort_values(
        ["score", "evalue"], ascending=[False, True], na_position="last"
    )
    fam_deduped = fam_df.drop_duplicates(
        subset=["protein_id", "_dedup_group"], keep="first"
    ).drop(columns=["_dedup_group"])
    n_removed = len(fam_df) - len(fam_deduped)
    if n_removed > 0:
        print(
            f"  DefenseFinder FAM dedup: removed {n_removed:,} redundant variant hits "
            f"({len(fam_df):,} → {len(fam_deduped):,} FAM rows)"
        )
    return pd.concat([no_fam_df, fam_deduped], ignore_index=True)


def load_annotations(db_path: str, tsv_path: str, source: str,
                     evalue_threshold: float | None = None) -> None:
    conn = duckdb.connect(db_path)

    source = source.lower()  # Enforce lowercase source names
    print(f"Loading {tsv_path} as source={source}...")
    df = pd.read_csv(tsv_path, sep="\t")
    print(f"  Raw rows: {len(df):,}")

    # Normalize column names to match Sharur schema
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

    # Apply e-value filtering
    if evalue_threshold is None:
        evalue_threshold = DEFAULT_EVALUE_THRESHOLDS.get(source)
    if evalue_threshold is not None and "evalue" in df.columns:
        before_filter = len(df)
        df = df[df["evalue"].isna() | (df["evalue"] <= evalue_threshold)]
        dropped = before_filter - len(df)
        print(f"  E-value filter (e <= {evalue_threshold:.0e}): kept {len(df):,}, dropped {dropped:,}")

    # DefenseFinder FAM variant dedup: keep best hit per
    # (protein_id, system_component) when multiple _FAM_N profiles match.
    if source == "defensefinder":
        df = _dedup_defensefinder_fam_variants(df)

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
    parser.add_argument("--evalue", type=float, default=None,
                        help="E-value threshold (default: source-specific, see DEFAULT_EVALUE_THRESHOLDS)")
    args = parser.parse_args()
    load_annotations(args.db, args.tsv, args.source, evalue_threshold=args.evalue)
