#!/usr/bin/env python3
"""
Run DefenseFinder proper (MacSyFinder system-level validation) on all genomes
in a dataset and integrate system-level annotations into the Sharur database.

DefenseFinder proper goes beyond raw HMM hits (which is all Astra provides):
- Checks gene co-localization within inter_gene_max_space
- Validates mandatory/accessory component requirements
- Produces system-level calls (e.g., "RM_Type_II locus at contig X, genes 85-86")

This eliminates superfamily false positives (e.g., standalone Mokosh kinases,
scattered RM MTase Rossmann folds) by requiring co-localization of system components.

Usage:
    python scripts/run_defensefinder_systems.py data/susan_genomes/ [--workers 8]
"""

from __future__ import annotations

import argparse
import logging
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import duckdb
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def find_genome_faas(data_dir: Path) -> list[Path]:
    """Find per-genome FAA files from stage03_prodigal."""
    prodigal_dir = data_dir / "stage03_prodigal" / "genomes"
    if not prodigal_dir.exists():
        logger.error(f"Prodigal output not found: {prodigal_dir}")
        sys.exit(1)

    faa_files = []
    for genome_dir in sorted(prodigal_dir.iterdir()):
        if not genome_dir.is_dir() or genome_dir.name == "all_protein_symlinks":
            continue
        faa = genome_dir / f"{genome_dir.name}.faa"
        if faa.exists():
            faa_files.append(faa)
        else:
            logger.warning(f"  Missing FAA for {genome_dir.name}")

    return faa_files


def run_defensefinder_on_genome(
    faa_path: Path, output_dir: Path, workers: int = 8, db_type: str = "gembase"
) -> tuple[pd.DataFrame, pd.DataFrame] | None:
    """Run defense-finder on a single genome FAA, return (systems_df, genes_df)."""
    genome_name = faa_path.stem
    genome_out = output_dir / genome_name

    if genome_out.exists():
        shutil.rmtree(genome_out)

    cmd = [
        "defense-finder", "run",
        str(faa_path),
        "-o", str(genome_out),
        "-w", str(workers),
        "--preserve-raw",
        "--db-type", db_type,
    ]

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        # Check if we got partial results (DefenseFinder + RM may succeed even if CasFinder fails)
        stderr = result.stderr or ""
        if "Error" in stderr or "error" in stderr.lower():
            logger.warning(f"  defense-finder returned non-zero for {genome_name}: {stderr[-200:]}")

    # Read output TSVs
    systems_tsv = genome_out / f"{genome_name}_defense_finder_systems.tsv"
    genes_tsv = genome_out / f"{genome_name}_defense_finder_genes.tsv"

    if not genes_tsv.exists():
        # Try to read from raw MacSyFinder output if post-treatment failed
        logger.warning(f"  No post-treatment output for {genome_name}, checking raw results...")
        return _read_raw_macsyfinder_results(genome_out / "defense-finder-tmp", genome_name)

    systems_df = pd.read_csv(systems_tsv, sep="\t") if systems_tsv.exists() else pd.DataFrame()
    genes_df = pd.read_csv(genes_tsv, sep="\t") if genes_tsv.exists() else pd.DataFrame()

    return systems_df, genes_df


def _read_raw_macsyfinder_results(
    tmp_dir: Path, genome_name: str
) -> tuple[pd.DataFrame, pd.DataFrame] | None:
    """Fallback: read raw MacSyFinder best_solution.tsv files when post-treatment fails."""
    all_genes = []

    for family_dir in sorted(tmp_dir.iterdir()):
        if not family_dir.is_dir():
            continue
        bs_file = family_dir / "best_solution.tsv"
        if bs_file.exists() and bs_file.stat().st_size > 0:
            try:
                df = pd.read_csv(bs_file, sep="\t", comment="#")
                if not df.empty:
                    all_genes.append(df)
            except Exception as e:
                logger.warning(f"  Failed to read {bs_file}: {e}")

    if not all_genes:
        logger.info(f"  No systems found in raw results for {genome_name}")
        return pd.DataFrame(), pd.DataFrame()

    genes_df = pd.concat(all_genes, ignore_index=True)
    # Build minimal systems summary from genes
    if "sys_id" in genes_df.columns:
        systems_df = (
            genes_df.groupby("sys_id")
            .agg(
                genes_count=("hit_id", "count"),
                protein_in_syst=("hit_id", lambda x: ",".join(x)),
            )
            .reset_index()
        )
    else:
        systems_df = pd.DataFrame()

    return systems_df, genes_df


def integrate_results(
    data_dir: Path,
    all_systems: pd.DataFrame,
    all_genes: pd.DataFrame,
) -> None:
    """Load defense-finder system-level results into the Sharur database."""
    db_path = data_dir / "sharur.duckdb"
    if not db_path.exists():
        logger.error(f"Database not found: {db_path}")
        sys.exit(1)

    conn = duckdb.connect(str(db_path))

    # Check if defense_systems table exists, create if not
    tables = [r[0] for r in conn.execute("SHOW TABLES").fetchall()]
    if "defense_systems" not in tables:
        conn.execute("""
            CREATE TABLE defense_systems (
                system_id VARCHAR PRIMARY KEY,
                genome_id VARCHAR,
                system_type VARCHAR,
                system_subtype VARCHAR,
                activity VARCHAR,
                genes_count INTEGER,
                protein_ids VARCHAR,
                profile_names VARCHAR,
                sys_beg VARCHAR,
                sys_end VARCHAR,
                created_at TIMESTAMP
            )
        """)
        logger.info("Created defense_systems table")
    else:
        before_count = conn.execute("SELECT COUNT(*) FROM defense_systems").fetchone()[0]
        conn.execute("DELETE FROM defense_systems")
        logger.info(f"Cleared {before_count} existing defense_systems rows")

    # Insert systems
    if not all_systems.empty:
        now = datetime.now(timezone.utc)
        rows = []
        seen_ids: set[str] = set()
        for _, row in all_systems.iterrows():
            sys_id = row.get("sys_id", "")
            # MacSyFinder can produce duplicate sys_ids across model runs (DefenseFinder vs RM)
            # Make unique by appending a suffix if needed
            orig_id = sys_id
            suffix = 1
            while sys_id in seen_ids:
                sys_id = f"{orig_id}_{suffix}"
                suffix += 1
            seen_ids.add(sys_id)
            rows.append({
                "system_id": sys_id,
                "genome_id": row.get("genome_id", ""),
                "system_type": row.get("type", ""),
                "system_subtype": row.get("subtype", ""),
                "activity": row.get("activity", "Defense"),
                "genes_count": int(row.get("genes_count", 0)),
                "protein_ids": row.get("protein_in_syst", ""),
                "profile_names": row.get("name_of_profiles_in_sys", ""),
                "sys_beg": row.get("sys_beg", ""),
                "sys_end": row.get("sys_end", ""),
                "created_at": now,
            })

        systems_insert = pd.DataFrame(rows)
        conn.register("tmp_systems", systems_insert)
        conn.execute("INSERT INTO defense_systems SELECT * FROM tmp_systems")
        conn.unregister("tmp_systems")
        logger.info(f"Inserted {len(systems_insert)} defense systems")

    # Also add system membership to annotations table
    # This creates source='defensefinder_system' annotations that link proteins to validated systems
    if not all_genes.empty and "hit_id" in all_genes.columns:
        # Check existing defensefinder_system annotations
        existing = conn.execute(
            "SELECT COUNT(*) FROM annotations WHERE source = 'defensefinder_system'"
        ).fetchone()[0]
        if existing > 0:
            conn.execute("DELETE FROM annotations WHERE source = 'defensefinder_system'")
            logger.info(f"Cleared {existing} existing defensefinder_system annotations")

        # Get next annotation_id
        next_id = conn.execute("SELECT COALESCE(MAX(annotation_id), 0) FROM annotations").fetchone()[0]

        # Filter to proteins that exist in DB
        existing_proteins = set(
            r[0] for r in conn.execute("SELECT DISTINCT protein_id FROM proteins").fetchall()
        )

        ann_rows = []
        for _, row in all_genes.iterrows():
            pid = row.get("hit_id", "")
            if pid not in existing_proteins:
                continue
            next_id += 1
            sys_type = row.get("type", "")
            sys_subtype = row.get("subtype", "")
            sys_id = row.get("sys_id", "")
            gene_name = row.get("gene_name", "")

            ann_rows.append({
                "annotation_id": next_id,
                "protein_id": pid,
                "source": "defensefinder_system",
                "accession": gene_name,
                "name": f"{sys_type}/{sys_subtype}",
                "description": f"System: {sys_id}",
                "evalue": row.get("hit_i_eval"),
                "score": row.get("hit_score"),
                "start_aa": row.get("hit_begin_match"),
                "end_aa": row.get("hit_end_match"),
            })

        if ann_rows:
            ann_df = pd.DataFrame(ann_rows)
            keep_cols = [
                "annotation_id", "protein_id", "source", "accession", "name",
                "description", "evalue", "score", "start_aa", "end_aa",
            ]
            ann_df = ann_df.reindex(columns=keep_cols, fill_value=None)
            conn.register("tmp_ann", ann_df)
            conn.execute("""
                INSERT INTO annotations (annotation_id, protein_id, source, accession,
                                         name, description, evalue, score, start_aa, end_aa)
                SELECT * FROM tmp_ann
            """)
            conn.unregister("tmp_ann")
            logger.info(f"Inserted {len(ann_df)} defensefinder_system annotations")

    # Summary comparison
    astra_count = conn.execute(
        "SELECT COUNT(DISTINCT protein_id) FROM annotations WHERE source = 'defensefinder'"
    ).fetchone()[0]
    system_count = conn.execute(
        "SELECT COUNT(DISTINCT protein_id) FROM annotations WHERE source = 'defensefinder_system'"
    ).fetchone()[0]
    logger.info(f"\nComparison:")
    logger.info(f"  Astra HMM-only (source='defensefinder'): {astra_count} unique proteins")
    logger.info(f"  System-validated (source='defensefinder_system'): {system_count} unique proteins")
    if astra_count > 0:
        logger.info(f"  Reduction: {astra_count - system_count} proteins were HMM hits without system context")

    # System type summary
    if not all_systems.empty:
        type_counts = all_systems.groupby("type").size().sort_values(ascending=False)
        logger.info("\nSystem type summary:")
        for stype, count in type_counts.items():
            logger.info(f"  {stype}: {count} systems")

    conn.commit()
    conn.close()


def main():
    parser = argparse.ArgumentParser(description="Run DefenseFinder system-level validation")
    parser.add_argument("data_dir", type=Path, help="Dataset directory (e.g., data/susan_genomes/)")
    parser.add_argument("--workers", type=int, default=8, help="Number of HMM search workers")
    parser.add_argument("--keep-raw", action="store_true", help="Keep raw defense-finder output")
    parser.add_argument(
        "--db-type",
        default="gembase",
        choices=["gembase", "ordered_replicon", "unordered"],
        help=(
            "gembase scopes each contig as its own replicon and is correct for MAGs. "
            "ordered_replicon declares a whole bin one replicon; that is only safe for "
            "model sets that are strictly single-locus, and silently permits "
            "cross-contig assembly for any model allowing multiple loci."
        ),
    )
    args = parser.parse_args()

    data_dir = args.data_dir
    faa_files = find_genome_faas(data_dir)
    logger.info(f"Found {len(faa_files)} genome FAA files")

    # Create output directory
    df_output_dir = data_dir / "stage04_astra" / "defensefinder_systems"
    df_output_dir.mkdir(parents=True, exist_ok=True)

    all_systems_list = []
    all_genes_list = []
    total_genomes = len(faa_files)

    for i, faa_path in enumerate(faa_files, 1):
        genome_name = faa_path.stem
        logger.info(f"\n[{i}/{total_genomes}] Processing {genome_name}...")

        result = run_defensefinder_on_genome(
            faa_path, df_output_dir, workers=args.workers, db_type=args.db_type
        )
        if result is None:
            continue

        systems_df, genes_df = result

        # Tag with genome_id
        bin_id = genome_name
        if not systems_df.empty:
            systems_df["genome_id"] = bin_id
            all_systems_list.append(systems_df)
            logger.info(f"  {len(systems_df)} systems found")

        if not genes_df.empty:
            genes_df["genome_id"] = bin_id
            all_genes_list.append(genes_df)
            logger.info(f"  {len(genes_df) - 1 if 'replicon' in genes_df.columns else len(genes_df)} genes in systems")

    # Merge results
    all_systems = pd.concat(all_systems_list, ignore_index=True) if all_systems_list else pd.DataFrame()
    all_genes = pd.concat(all_genes_list, ignore_index=True) if all_genes_list else pd.DataFrame()

    logger.info(f"\n{'='*60}")
    logger.info(f"TOTAL: {len(all_systems)} systems, {len(all_genes)} gene assignments across {total_genomes} genomes")

    # Save merged results
    if not all_systems.empty:
        all_systems.to_csv(df_output_dir / "all_defense_systems.tsv", sep="\t", index=False)
    if not all_genes.empty:
        all_genes.to_csv(df_output_dir / "all_defense_genes.tsv", sep="\t", index=False)

    # Integrate into database
    integrate_results(data_dir, all_systems, all_genes)

    # Cleanup raw output if not keeping
    if not args.keep_raw:
        for genome_dir in df_output_dir.iterdir():
            if genome_dir.is_dir() and genome_dir.name != "all_protein_symlinks":
                raw_dir = genome_dir / "defense-finder-tmp"
                if raw_dir.exists():
                    shutil.rmtree(raw_dir)

    logger.info("\nDone.")


if __name__ == "__main__":
    main()
