#!/usr/bin/env python3
"""
Targeted reload of DefenseFinder annotations with FAM variant deduplication.

This script:
1. Reads the existing DefenseFinder TSV from stage04_astra
2. Applies the same FAM dedup logic as 07_build_knowledge_base.py
3. Deletes existing defensefinder annotations from the DB
4. Re-inserts the deduped annotations
5. Regenerates predicates for proteins that had defensefinder annotations

Usage:
    python scripts/reload_defensefinder_dedup.py data/susan_genomes/
"""

from __future__ import annotations

import logging
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import duckdb
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Match _FAM_N or _FAMN at end of accession (same regex as 07_build_knowledge_base.py)
_DEFENSEFINDER_FAM_RE = re.compile(r"_FAM_?\d+$")


def _dedup_defensefinder_fam_variants(df: pd.DataFrame) -> pd.DataFrame:
    """Deduplicate DefenseFinder hits that differ only by FAM variant number."""
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
        logger.info(
            f"  defensefinder: deduplicated {n_removed:,} redundant FAM variant hits "
            f"({len(fam_df):,} -> {len(fam_deduped):,} FAM rows)"
        )

    result = pd.concat([no_fam_df, fam_deduped], ignore_index=True)
    return result


def reload_defensefinder(data_dir: Path) -> None:
    db_path = data_dir / "sharur.duckdb"
    if not db_path.exists():
        logger.error(f"Database not found: {db_path}")
        sys.exit(1)

    # Find ALL DefenseFinder TSVs — stage07 uses rglob("*_hits_df.tsv") and
    # loads every TSV under stage04_astra that contains "defense" in its path.
    # We must do the same to avoid losing proteins from secondary TSVs.
    stage04_dir = data_dir / "stage04_astra"
    tsv_candidates = [
        p for p in stage04_dir.rglob("*_hits_df.tsv")
        if "defense" in p.as_posix().lower() and "dbcan" not in p.as_posix().lower()
    ]
    if not tsv_candidates:
        logger.error(f"No DefenseFinder TSV found under {stage04_dir}")
        sys.exit(1)

    logger.info(f"Found {len(tsv_candidates)} DefenseFinder TSV file(s):")
    dfs = []
    for tsv_path in tsv_candidates:
        logger.info(f"  Loading: {tsv_path.relative_to(data_dir)}")
        df_part = pd.read_csv(tsv_path, sep="\t")
        logger.info(f"    {len(df_part):,} rows, {df_part['sequence_id'].nunique():,} unique proteins")
        dfs.append(df_part)

    df = pd.concat(dfs, ignore_index=True)
    logger.info(f"Combined: {len(df):,} total rows, {df['sequence_id'].nunique():,} unique proteins")

    # Normalize columns to match what 07_build_knowledge_base.py does.
    # The TSV may use 'bit_score' or 'bitscore'; handle both.
    df = df.rename(columns={
        "sequence_id": "protein_id",
        "hmm_name": "accession",
        "bit_score": "score",       # Astra newer output
        "bitscore": "score_primary", # Astra older output (will be used as fallback below)
        "dom_bitscore": "dom_bitscore",
        "i_evalue": "i_evalue",
        "c_evalue": "c_evalue",
        "ali_from": "start_aa",
        "ali_to": "end_aa",
        "env_from": "env_from",
        "env_to": "env_to",
    })

    # Fill e-value / score from alternate columns
    if "evalue" not in df.columns:
        df["evalue"] = None
    if "score" not in df.columns:
        df["score"] = None
    if "i_evalue" in df.columns:
        df["evalue"] = df["evalue"].fillna(df["i_evalue"])
    if "c_evalue" in df.columns:
        df["evalue"] = df["evalue"].fillna(df["c_evalue"])
    if "dom_bitscore" in df.columns:
        df["score"] = df["score"].fillna(df["dom_bitscore"])
    if "score_primary" in df.columns:
        df["score"] = df["score"].fillna(df["score_primary"])
    if "start_aa" not in df.columns:
        df["start_aa"] = None
    if "end_aa" not in df.columns:
        df["end_aa"] = None
    if "env_from" in df.columns:
        df["start_aa"] = df["start_aa"].fillna(df["env_from"])
    if "env_to" in df.columns:
        df["end_aa"] = df["end_aa"].fillna(df["env_to"])

    # Apply e-value threshold (1e-15 for defensefinder, same as 07_build_knowledge_base.py)
    evalue_threshold = 1e-15
    if "evalue" in df.columns:
        before_filter = len(df)
        df = df[df["evalue"].isna() | (df["evalue"] <= evalue_threshold)]
        dropped = before_filter - len(df)
        if dropped > 0:
            logger.info(f"  Filtered {dropped:,} hits with e-value > {evalue_threshold:.0e} ({len(df):,} retained)")

    # Apply FAM dedup
    df_deduped = _dedup_defensefinder_fam_variants(df)
    logger.info(f"After dedup: {len(df_deduped):,} rows ({len(df):,} before dedup)")

    # Fill name/description if not present
    if "name" not in df_deduped.columns:
        df_deduped["name"] = df_deduped.get("accession")
    if "description" not in df_deduped.columns:
        df_deduped["description"] = ""
    df_deduped["description"] = df_deduped["description"].fillna("")

    # Set source
    df_deduped["source"] = "defensefinder"

    # Connect to database
    conn = duckdb.connect(str(db_path))

    # Record before state
    before_rows = conn.execute("SELECT COUNT(*) FROM annotations WHERE source = 'defensefinder'").fetchone()[0]
    before_unique = conn.execute("SELECT COUNT(DISTINCT protein_id) FROM annotations WHERE source = 'defensefinder'").fetchone()[0]
    logger.info(f"BEFORE: {before_rows:,} rows, {before_unique:,} unique proteins")

    # Delete existing defensefinder annotations
    conn.execute("DELETE FROM annotations WHERE source = 'defensefinder'")
    logger.info("Deleted existing defensefinder annotations")

    # Get next annotation_id
    next_id = conn.execute("SELECT COALESCE(MAX(annotation_id), 0) FROM annotations").fetchone()[0]

    # Select only the columns we need
    keep = ["annotation_id", "protein_id", "source", "accession", "name",
            "description", "evalue", "score", "start_aa", "end_aa"]
    df_deduped = df_deduped.reindex(columns=keep, fill_value=None)
    df_deduped["annotation_id"] = range(next_id + 1, next_id + 1 + len(df_deduped))

    # Filter to only proteins that exist in DB (safety check)
    existing_proteins = set(r[0] for r in conn.execute("SELECT DISTINCT protein_id FROM proteins").fetchall())
    before_filter_prot = len(df_deduped)
    df_deduped = df_deduped[df_deduped["protein_id"].isin(existing_proteins)]
    if len(df_deduped) < before_filter_prot:
        logger.info(f"  Filtered {before_filter_prot - len(df_deduped):,} rows for proteins not in DB")

    # Insert
    conn.register("tmp_df_ann", df_deduped)
    conn.execute("""
        INSERT INTO annotations (annotation_id, protein_id, source, accession, name,
                                 description, evalue, score, start_aa, end_aa)
        SELECT * FROM tmp_df_ann
    """)
    conn.unregister("tmp_df_ann")
    logger.info(f"Inserted {len(df_deduped):,} deduped rows")

    # Regenerate predicates for affected proteins
    affected_proteins = set(df_deduped["protein_id"].dropna().tolist())
    logger.info(f"Regenerating predicates for {len(affected_proteins):,} affected proteins...")

    # Import predicate generator
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
    from sharur.predicates.generator import PredicateGenerator, AnnotationRecord, ProteinRecord

    gen = PredicateGenerator(predict_topology=False)

    # Get protein data for affected proteins
    affected_list = list(affected_proteins)
    placeholders = ",".join(["?"] * len(affected_list))
    proteins = conn.execute(
        f"SELECT protein_id, sequence_length, gc_content, contig_id FROM proteins WHERE protein_id IN ({placeholders})",
        affected_list
    ).fetchall()

    # Get GC stats per contig
    gc_stats = {}
    gc_rows = conn.execute("""
        SELECT contig_id, AVG(gc_content), STDDEV(gc_content)
        FROM proteins
        WHERE gc_content IS NOT NULL
        GROUP BY contig_id
    """).fetchall()
    for row in gc_rows:
        gc_stats[row[0]] = (row[1], row[2] if row[2] else 0.0)

    # Get all annotations for affected proteins (full set, not just defensefinder)
    ann_rows = conn.execute(
        f"SELECT protein_id, source, accession, name, description, evalue, score FROM annotations WHERE protein_id IN ({placeholders})",
        affected_list
    ).fetchall()
    annotations_by_protein: dict = {}
    for row in ann_rows:
        pid = row[0]
        if pid not in annotations_by_protein:
            annotations_by_protein[pid] = []
        annotations_by_protein[pid].append((row[1], row[2], row[3], row[4], row[5], row[6]))

    # Generate new predicates
    updated = 0
    for pid, length, gc, contig in proteins:
        gc_mean, gc_std = gc_stats.get(contig, (None, None))
        protein = ProteinRecord(
            protein_id=pid,
            sequence_length=length,
            gc_content=gc,
            contig_gc_mean=gc_mean,
            contig_gc_std=gc_std,
        )
        anns = [
            AnnotationRecord(source=t[0], accession=t[1], name=t[2],
                             description=t[3], evalue=t[4], score=t[5])
            for t in annotations_by_protein.get(pid, [])
        ]
        predicates = gen.generate_for_protein(protein, anns)

        conn.execute("""
            UPDATE protein_predicates
            SET predicates = ?, updated_at = ?
            WHERE protein_id = ?
        """, [predicates, datetime.now(timezone.utc), pid])
        updated += 1

    logger.info(f"Updated predicates for {updated:,} proteins")

    # Verify after state
    after_rows = conn.execute("SELECT COUNT(*) FROM annotations WHERE source = 'defensefinder'").fetchone()[0]
    after_unique = conn.execute("SELECT COUNT(DISTINCT protein_id) FROM annotations WHERE source = 'defensefinder'").fetchone()[0]
    logger.info(f"AFTER: {after_rows:,} rows, {after_unique:,} unique proteins")
    logger.info(f"Row reduction: {before_rows:,} -> {after_rows:,} ({before_rows - after_rows:,} removed)")
    logger.info(f"Unique proteins preserved: {before_unique:,} -> {after_unique:,}")

    conn.commit()
    conn.close()
    logger.info("Done. Database updated successfully.")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <data_dir>")
        sys.exit(1)
    reload_defensefinder(Path(sys.argv[1]))
