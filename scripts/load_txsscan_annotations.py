#!/usr/bin/env python3
"""Load MacSyFinder HMM hits into `annotations` as source='txsscan'.

The co-location engine validates hits that already live in DuckDB, so it can
only run for a system family whose profile hits have been ingested. Astra's
standard database set does not include the secretion-system profiles, which is
why `validate_systems(source='txsscan')` has never had input to work on.

This ingests the per-profile hit tables MacSyFinder writes alongside its system
calls (`hmmer_results/*.res_hmm_extract`), giving the engine the same hits the
reference implementation used. That makes a like-for-like concordance comparison
possible: identical input, so any difference in output is attributable to
clustering, quorum or best-solution logic rather than to search settings.

The extract header records the thresholds the reference run applied
(`i_evalue threshold`, `coverage threshold`). They are parsed and returned so a
later astra-based ingest can be configured to match rather than guessed at --
the profiles carry no GA/TC/NC cutoffs, so an astra run defaults to different
filtering and would silently diverge at the hit level.

Usage:
    python scripts/load_txsscan_annotations.py RESULTS_DIR --out-dir DIR --shard i/N
    python scripts/load_txsscan_annotations.py RESULTS_DIR --out-dir DIR --load DATASET_DIR
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import duckdb
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

SOURCE = "txsscan"
COLUMNS = [
    "hit_id", "replicon_name", "position_hit", "hit_sequence_length",
    "gene_name", "gene_system", "i_eval", "score",
    "profile_coverage", "sequence_coverage", "begin", "end",
]


def parse_extract(path: Path) -> tuple[list[dict], dict]:
    """Parse one .res_hmm_extract into (rows, thresholds)."""
    rows: list[dict] = []
    thresholds: dict = {}
    try:
        text = path.read_text()
    except OSError:
        return rows, thresholds
    for line in text.splitlines():
        line = line.rstrip("\n")
        if not line.strip():
            continue
        if line.startswith("#"):
            low = line.lower()
            if "i_evalue threshold" in low:
                thresholds["i_evalue"] = line.split("=")[-1].strip()
            elif "coverage threshold" in low:
                thresholds["coverage"] = line.split("=")[-1].strip()
            continue
        parts = line.split()
        if len(parts) < len(COLUMNS):
            continue
        rec = dict(zip(COLUMNS, parts[: len(COLUMNS)]))
        rows.append(rec)
    return rows, thresholds


def collect(results_dir: Path, shard: tuple[int, int] | None) -> tuple[pd.DataFrame, dict]:
    genomes = sorted(p for p in results_dir.iterdir() if p.is_dir())
    if shard is not None:
        i, n = shard
        genomes = [g for k, g in enumerate(genomes) if k % n == i]
    logger.info("shard covers %d genome dirs", len(genomes))

    all_rows: list[dict] = []
    thresholds: dict = {}
    for gi, gdir in enumerate(genomes, start=1):
        hdir = gdir / "hmmer_results"
        if not hdir.is_dir():
            continue
        for ext in hdir.glob("*.res_hmm_extract"):
            rows, th = parse_extract(ext)
            thresholds.update(th)
            for r in rows:
                r["genome_id"] = gdir.name
            all_rows.extend(rows)
        if gi % 50 == 0:
            logger.info("  %d/%d genomes, %d hits", gi, len(genomes), len(all_rows))
    return pd.DataFrame(all_rows), thresholds


def load_into_db(data_dir: Path, frames: list[Path]) -> dict:
    """Replace source='txsscan' rows in `annotations` with the parsed hits."""
    db_path = data_dir / "sharur.duckdb"
    dfs = [pd.read_parquet(f) for f in frames if f.exists()]
    df = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()
    if df.empty:
        logger.warning("no txsscan hits parsed; nothing to load")
        return {"rows": 0}

    df = df.rename(columns={"hit_id": "protein_id", "gene_name": "accession"})
    df["source"] = SOURCE
    df["evalue"] = pd.to_numeric(df["i_eval"], errors="coerce")
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    df["start_aa"] = pd.to_numeric(df["begin"], errors="coerce").astype("Int64")
    df["end_aa"] = pd.to_numeric(df["end"], errors="coerce").astype("Int64")
    df["name"] = df["accession"]
    df["description"] = df["gene_system"]
    keep = ["protein_id", "source", "accession", "name", "description",
            "evalue", "score", "start_aa", "end_aa"]
    df = df[keep].drop_duplicates(subset=["protein_id", "accession", "start_aa", "end_aa"])

    con = duckdb.connect(str(db_path))
    cols = [r[0] for r in con.execute("DESCRIBE annotations").fetchall()]
    con.execute("DELETE FROM annotations WHERE source = ?", [SOURCE])
    con.register("incoming", df)
    # annotation_id may be generated; insert only the columns that exist.
    shared = [c for c in keep if c in cols]
    if "annotation_id" in cols:
        con.execute(
            f"""INSERT INTO annotations (annotation_id, {','.join(shared)})
                SELECT 'txsscan_' || CAST(row_number() OVER () AS VARCHAR),
                       {','.join(shared)} FROM incoming"""
        )
    else:
        con.execute(f"INSERT INTO annotations ({','.join(shared)}) SELECT {','.join(shared)} FROM incoming")
    n = con.execute("SELECT COUNT(*) FROM annotations WHERE source = ?", [SOURCE]).fetchone()[0]
    prot = con.execute(
        "SELECT COUNT(DISTINCT protein_id) FROM annotations WHERE source = ?", [SOURCE]
    ).fetchone()[0]
    con.close()
    logger.info("loaded %d txsscan annotations over %d proteins", n, prot)
    return {"rows": int(n), "proteins": int(prot)}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("results_dir", type=Path, help="MacSyFinder per-genome results root")
    ap.add_argument("--out-dir", type=Path, required=True, help="Where shard parquets are written")
    ap.add_argument("--shard", type=str, default=None, help="i/N")
    ap.add_argument("--load", type=Path, default=None, help="Dataset dir; load all shard parquets into its DuckDB")
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    if args.load is not None:
        stats = load_into_db(args.load, sorted(args.out_dir.glob("shard_*.parquet")))
        (args.out_dir / "load_manifest.json").write_text(json.dumps(stats, indent=2))
        return 0

    shard = None
    tag = "all"
    if args.shard:
        i, n = args.shard.split("/")
        shard = (int(i), int(n))
        tag = f"{int(i):03d}"
    df, thresholds = collect(args.results_dir, shard)
    out = args.out_dir / f"shard_{tag}.parquet"
    df.to_parquet(out, index=False)
    if thresholds:
        (args.out_dir / f"thresholds_{tag}.json").write_text(json.dumps(thresholds, indent=2))
    logger.info("wrote %d hits to %s (thresholds: %s)", len(df), out, thresholds)
    return 0


if __name__ == "__main__":
    sys.exit(main())
