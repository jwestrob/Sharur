#!/usr/bin/env python3
"""
Validate DefenseFinder HMM hits via MacSyFinder co-localization.

Takes the global Astra DefenseFinder output (with --write_macsyfinder), splits
the hmmsearch results per-genome, and runs MacSyFinder --previous-run on each
genome to check gene co-localization.  This eliminates ~72% of HMM-only false
positives (standalone kinases, scattered Rossmann MTases, etc.).

Works with both binned MAGs (per-genome dirs in stage03) and unbinned
metagenomes (falls back to unordered mode on the full set).

Usage:
    # Standalone
    python scripts/validate_defense_systems.py data/susan_genomes/ [--workers 4]

    # As library (called from 07_build_knowledge_base.py)
    from validate_defense_systems import validate_defense_systems
    systems_df, genes_df = validate_defense_systems(db_path, data_dir)
"""

from __future__ import annotations

import argparse
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import duckdb
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


# ------------------------------------------------------------------ #
# Protein → genome mapping from stage03 directory structure
# ------------------------------------------------------------------ #

def _build_protein_to_genome(genomes_dir: Path) -> dict[str, str]:
    """Read per-genome FAA headers to build protein_id → genome_name map."""
    mapping: dict[str, str] = {}
    for genome_dir in sorted(genomes_dir.iterdir()):
        if not genome_dir.is_dir() or genome_dir.name == "all_protein_symlinks":
            continue
        genome_name = genome_dir.name
        faa = genome_dir / f"{genome_name}.faa"
        if not faa.exists():
            # Try any .faa in the directory
            faas = list(genome_dir.glob("*.faa"))
            if not faas:
                continue
            faa = faas[0]
        with open(faa) as f:
            for line in f:
                if line.startswith(">"):
                    pid = line[1:].split()[0].split(" # ")[0]
                    mapping[pid] = genome_name
    return mapping


def _build_protein_to_genome_from_db(db_path: Path) -> dict[str, str]:
    """Fallback: build protein_id → bin_id mapping from the database."""
    conn = duckdb.connect(str(db_path), read_only=True)
    rows = conn.execute("SELECT protein_id, bin_id FROM proteins").fetchall()
    conn.close()
    return {r[0]: r[1] for r in rows}


# ------------------------------------------------------------------ #
# Filter hmmsearch output files per-genome
# ------------------------------------------------------------------ #

def _filter_hmmsearch_file(
    src_path: Path, dst_path: Path, protein_ids: set[str]
) -> int:
    """Filter a .search_hmm.out file to keep only hits for given protein IDs.

    Returns the number of kept hit blocks.
    """
    kept = 0
    with open(src_path) as fin, open(dst_path, "w") as fout:
        in_hit_block = False
        keep_block = False
        block_lines: list[str] = []

        for line in fin:
            if line.startswith(">> "):
                # Flush previous block
                if in_hit_block and keep_block:
                    fout.writelines(block_lines)
                    kept += 1
                # Start new block
                pid = line[3:].split()[0]
                keep_block = pid in protein_ids
                block_lines = [line]
                in_hit_block = True
            elif line.startswith("//"):
                # End-of-query marker: flush last block, write marker
                if in_hit_block and keep_block:
                    fout.writelines(block_lines)
                    kept += 1
                fout.write(line)
                in_hit_block = False
                block_lines = []
            elif in_hit_block:
                block_lines.append(line)
            else:
                # Header lines before first hit block
                fout.write(line)

    return kept


def _prepare_genome_macsyfinder_dir(
    genome_name: str,
    genome_faa: Path,
    global_hmmer_dir: Path,
    protein_ids: set[str],
    work_dir: Path,
) -> Optional[Path]:
    """Create a per-genome macsyfinder_compat directory with filtered hmmsearch output.

    Returns the path to the created directory, or None if no hits were found.
    """
    genome_dir = work_dir / genome_name
    hmmer_out = genome_dir / "hmmer_results"
    hmmer_out.mkdir(parents=True, exist_ok=True)

    total_kept = 0
    for src_file in sorted(global_hmmer_dir.glob("*.search_hmm.out")):
        dst_file = hmmer_out / src_file.name
        n = _filter_hmmsearch_file(src_file, dst_file, protein_ids)
        total_kept += n
        # Remove empty files (MacSyFinder is fine without them)
        if n == 0:
            dst_file.unlink()

    if total_kept == 0:
        shutil.rmtree(genome_dir, ignore_errors=True)
        return None

    # Write macsyfinder.conf pointing to this genome's FAA
    conf_path = genome_dir / "macsyfinder.conf"
    with open(conf_path, "w") as fh:
        fh.write("[base]\n")
        fh.write(f"sequence_db = {genome_faa.resolve()}\n")
        fh.write("db_type = ordered_replicon\n")
        fh.write("hmmer = hmmsearch\n\n")
        fh.write("[hmmer]\n")
        fh.write("e_value_search = 0.1\n")

    return genome_dir


# ------------------------------------------------------------------ #
# Run MacSyFinder --previous-run
# ------------------------------------------------------------------ #

def _find_models_dir() -> Optional[Path]:
    """Locate the defense-finder-models directory."""
    candidates = [
        Path.home() / ".macsyfinder" / "models",
        Path.home() / ".mdmlab" / "macsyfinder" / "models",
    ]
    for c in candidates:
        if (c / "defense-finder-models").is_dir():
            return c
    return None


def _run_macsyfinder_on_dir(
    genome_dir: Path, models_dir: Path, workers: int = 4
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run MacSyFinder --previous-run on a prepared per-genome directory.

    Returns (systems_df, genes_df).
    """
    # MacSyFinder writes output into a results directory inside the working dir
    results_dir = genome_dir / "macsyfinder_results"
    results_dir.mkdir(exist_ok=True)

    cmd = [
        "macsyfinder",
        "--previous-run", str(genome_dir),
        "--models-dir", str(models_dir),
        "--models", "defense-finder-models", "all",
        "-w", str(workers),
        "-o", str(results_dir),
    ]

    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True, timeout=300
        )
    except subprocess.TimeoutExpired:
        logger.warning(f"  MacSyFinder timed out for {genome_dir.name}")
        return pd.DataFrame(), pd.DataFrame()
    except FileNotFoundError:
        logger.error("macsyfinder not found on PATH")
        return pd.DataFrame(), pd.DataFrame()

    if proc.returncode != 0:
        stderr = proc.stderr or ""
        # Some warnings are non-fatal (e.g., no hits for some models)
        if "error" in stderr.lower() and "no systems found" not in stderr.lower():
            logger.warning(
                f"  MacSyFinder non-zero exit for {genome_dir.name}: "
                f"{stderr[-300:]}"
            )

    # Parse best_solution.tsv
    best_solution = results_dir / "best_solution.tsv"
    if not best_solution.exists():
        # Check inside subdirectories (MacSyFinder output structure varies)
        for bs in results_dir.rglob("best_solution.tsv"):
            best_solution = bs
            break

    if not best_solution.exists() or best_solution.stat().st_size == 0:
        return pd.DataFrame(), pd.DataFrame()

    try:
        df = pd.read_csv(best_solution, sep="\t", comment="#")
    except Exception as e:
        logger.warning(f"  Failed to parse {best_solution}: {e}")
        return pd.DataFrame(), pd.DataFrame()

    if df.empty:
        return pd.DataFrame(), pd.DataFrame()

    # Build systems summary from gene-level output
    genes_df = df.copy()

    # Extract system type from model_fqn (e.g. "defense-finder-models/DefenseFinder/Pycsar/Pycsar")
    # The second-to-last path component is the system type.
    if "model_fqn" in genes_df.columns:
        genes_df["type"] = genes_df["model_fqn"].apply(_extract_system_type)
        genes_df["subtype"] = genes_df["model_fqn"].apply(_extract_system_subtype)

    if "sys_id" in genes_df.columns:
        agg_dict = {
            "genes_count": ("hit_id", "count"),
            "protein_in_syst": ("hit_id", lambda x: ",".join(str(v) for v in x)),
            "name_of_profiles_in_sys": (
                "gene_name",
                lambda x: ",".join(str(v) for v in x),
            ),
        }
        if "type" in genes_df.columns:
            agg_dict["type"] = ("type", "first")
            agg_dict["subtype"] = ("subtype", "first")

        systems_df = (
            genes_df.groupby("sys_id")
            .agg(**agg_dict)
            .reset_index()
        )
    else:
        systems_df = pd.DataFrame()

    return systems_df, genes_df


def _extract_system_type(model_fqn: str) -> str:
    """Extract system type from model_fqn path.

    Example: 'defense-finder-models/DefenseFinder/Pycsar/Pycsar' → 'Pycsar'
             'defense-finder-models/DefenseFinder/RM/RM_Type_II' → 'RM_Type_II'
    """
    if not isinstance(model_fqn, str):
        return str(model_fqn)
    parts = model_fqn.rstrip("/").split("/")
    # The model name is the last component; the system family is second-to-last
    return parts[-1] if len(parts) >= 2 else model_fqn


def _extract_system_subtype(model_fqn: str) -> str:
    """Extract system subtype (same as type for DefenseFinder models)."""
    return _extract_system_type(model_fqn)


# ------------------------------------------------------------------ #
# Main validation entry point
# ------------------------------------------------------------------ #

def validate_defense_systems(
    db_path: str | Path,
    data_dir: str | Path | None = None,
    workers: int = 4,
    verbose: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run MacSyFinder co-localization validation on Astra DefenseFinder output.

    Args:
        db_path: Path to sharur.duckdb (used for protein→genome mapping fallback
                 and for integrating results).
        data_dir: Dataset root directory (contains stage03_prodigal/, stage04_astra/).
                  If None, inferred from db_path parent.
        workers: MacSyFinder worker threads per genome.
        verbose: Print progress to console.

    Returns:
        (systems_df, genes_df) — merged DataFrames across all genomes.
    """
    db_path = Path(db_path)
    if data_dir is None:
        data_dir = db_path.parent
    else:
        data_dir = Path(data_dir)

    # Locate Astra macsyfinder_compat output
    macsyfinder_compat = None
    for candidate in [
        data_dir / "stage04_astra" / "defensefinder_results" / "macsyfinder_compat",
    ]:
        if candidate.is_dir() and list(candidate.glob("hmmer_results/*.search_hmm.out")):
            macsyfinder_compat = candidate
            break

    if macsyfinder_compat is None:
        logger.error(
            "No macsyfinder_compat directory found. "
            "Re-run stage 04 Astra DefenseFinder with --write_macsyfinder."
        )
        return pd.DataFrame(), pd.DataFrame()

    global_hmmer_dir = macsyfinder_compat / "hmmer_results"
    n_hmm_files = len(list(global_hmmer_dir.glob("*.search_hmm.out")))
    if verbose:
        logger.info(f"Found {n_hmm_files} HMM result files in {global_hmmer_dir}")

    # Find MacSyFinder models
    models_dir = _find_models_dir()
    if models_dir is None:
        logger.error(
            "defense-finder-models not found. Install with: "
            "macsyfinder --install-models defense-finder-models"
        )
        return pd.DataFrame(), pd.DataFrame()

    # Build protein → genome mapping
    genomes_dir = data_dir / "stage03_prodigal" / "genomes"
    genome_dirs_exist = (
        genomes_dir.exists()
        and any(
            d.is_dir() and d.name != "all_protein_symlinks"
            for d in genomes_dir.iterdir()
        )
    )

    if genome_dirs_exist:
        if verbose:
            logger.info("Building protein→genome mapping from stage03 FAA files...")
        prot_to_genome = _build_protein_to_genome(genomes_dir)
    else:
        if verbose:
            logger.info("No per-genome dirs found; using DB for protein→genome mapping...")
        prot_to_genome = _build_protein_to_genome_from_db(db_path)

    if not prot_to_genome:
        logger.error("Could not build protein→genome mapping")
        return pd.DataFrame(), pd.DataFrame()

    # Group proteins by genome
    genome_proteins: dict[str, set[str]] = {}
    for pid, genome in prot_to_genome.items():
        genome_proteins.setdefault(genome, set()).add(pid)

    if verbose:
        logger.info(
            f"  {len(prot_to_genome):,} proteins across {len(genome_proteins)} genomes"
        )

    # Map genome names to their FAA files
    genome_faas: dict[str, Path] = {}
    if genome_dirs_exist:
        for genome_dir in genomes_dir.iterdir():
            if not genome_dir.is_dir() or genome_dir.name == "all_protein_symlinks":
                continue
            faa = genome_dir / f"{genome_dir.name}.faa"
            if faa.exists():
                genome_faas[genome_dir.name] = faa
            else:
                faas = list(genome_dir.glob("*.faa"))
                if faas:
                    genome_faas[genome_dir.name] = faas[0]

    # Create temp working directory for per-genome macsyfinder dirs
    work_dir = data_dir / "stage04_astra" / "defensefinder_results" / "_macsyfinder_work"
    if work_dir.exists():
        shutil.rmtree(work_dir)
    work_dir.mkdir(parents=True)

    # Process each genome
    all_systems: list[pd.DataFrame] = []
    all_genes: list[pd.DataFrame] = []
    total_genomes = len(genome_proteins)

    for i, (genome_name, protein_ids) in enumerate(
        sorted(genome_proteins.items()), 1
    ):
        if verbose:
            logger.info(f"  [{i}/{total_genomes}] {genome_name} ({len(protein_ids)} proteins)")

        genome_faa = genome_faas.get(genome_name)
        if genome_faa is None:
            logger.warning(f"    No FAA found for {genome_name}, skipping")
            continue

        # Filter hmmsearch output for this genome
        genome_msf_dir = _prepare_genome_macsyfinder_dir(
            genome_name, genome_faa, global_hmmer_dir, protein_ids, work_dir
        )
        if genome_msf_dir is None:
            if verbose:
                logger.info(f"    No DefenseFinder hits for {genome_name}")
            continue

        # Run MacSyFinder
        systems_df, genes_df = _run_macsyfinder_on_dir(
            genome_msf_dir, models_dir, workers=workers
        )

        if not systems_df.empty:
            systems_df["genome_id"] = genome_name
            all_systems.append(systems_df)
        if not genes_df.empty:
            genes_df["genome_id"] = genome_name
            all_genes.append(genes_df)

        n_sys = len(systems_df) if not systems_df.empty else 0
        if verbose and n_sys > 0:
            logger.info(f"    {n_sys} systems found")

    # Merge results
    merged_systems = (
        pd.concat(all_systems, ignore_index=True) if all_systems else pd.DataFrame()
    )
    merged_genes = (
        pd.concat(all_genes, ignore_index=True) if all_genes else pd.DataFrame()
    )

    if verbose:
        logger.info(
            f"\nTotal: {len(merged_systems)} systems, "
            f"{len(merged_genes)} gene assignments across {total_genomes} genomes"
        )

    # Save merged results
    results_out = data_dir / "stage04_astra" / "defensefinder_results"
    if not merged_systems.empty:
        merged_systems.to_csv(
            results_out / "validated_defense_systems.tsv", sep="\t", index=False
        )
    if not merged_genes.empty:
        merged_genes.to_csv(
            results_out / "validated_defense_genes.tsv", sep="\t", index=False
        )

    # Cleanup work dir
    shutil.rmtree(work_dir, ignore_errors=True)

    return merged_systems, merged_genes


# ------------------------------------------------------------------ #
# DB integration (reused from run_defensefinder_systems.py)
# ------------------------------------------------------------------ #

def integrate_results(
    db_path: Path,
    all_systems: pd.DataFrame,
    all_genes: pd.DataFrame,
) -> None:
    """Load defense system validation results into the Sharur database."""
    conn = duckdb.connect(str(db_path))

    # Create defense_systems table if needed
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
            sys_id = str(row.get("sys_id", ""))
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
                "activity": "Defense",
                "genes_count": int(row.get("genes_count", 0)),
                "protein_ids": row.get("protein_in_syst", ""),
                "profile_names": row.get("name_of_profiles_in_sys", ""),
                "sys_beg": str(row.get("sys_beg", "")),
                "sys_end": str(row.get("sys_end", "")),
                "created_at": now,
            })

        systems_insert = pd.DataFrame(rows)
        conn.register("tmp_systems", systems_insert)
        conn.execute("INSERT INTO defense_systems SELECT * FROM tmp_systems")
        conn.unregister("tmp_systems")
        logger.info(f"Inserted {len(systems_insert)} defense systems")

    # Add source='defensefinder_system' annotations
    if not all_genes.empty and "hit_id" in all_genes.columns:
        existing = conn.execute(
            "SELECT COUNT(*) FROM annotations WHERE source = 'defensefinder_system'"
        ).fetchone()[0]
        if existing > 0:
            conn.execute("DELETE FROM annotations WHERE source = 'defensefinder_system'")
            logger.info(f"Cleared {existing} existing defensefinder_system annotations")

        next_id = conn.execute(
            "SELECT COALESCE(MAX(annotation_id), 0) FROM annotations"
        ).fetchone()[0]

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
                "evalue": row.get("hit_i_eval") if "hit_i_eval" in row.index else row.get("i_eval"),
                "score": row.get("hit_score") if "hit_score" in row.index else row.get("score"),
                "start_aa": None,
                "end_aa": None,
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
        fp_rate = (astra_count - system_count) / astra_count * 100
        logger.info(f"  FP reduction: {astra_count - system_count} proteins ({fp_rate:.1f}%)")

    if not all_systems.empty and "type" in all_systems.columns:
        type_counts = all_systems.groupby("type").size().sort_values(ascending=False)
        logger.info("\nSystem type summary:")
        for stype, count in type_counts.items():
            logger.info(f"  {stype}: {count} systems")

    conn.commit()
    conn.close()


# ------------------------------------------------------------------ #
# CLI
# ------------------------------------------------------------------ #

def main():
    parser = argparse.ArgumentParser(
        description="Validate DefenseFinder HMM hits via MacSyFinder co-localization"
    )
    parser.add_argument(
        "data_dir", type=Path,
        help="Dataset directory (e.g., data/susan_genomes/)"
    )
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument(
        "--no-integrate", action="store_true",
        help="Skip database integration (just produce TSV output)"
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

    logger.info("Done.")


if __name__ == "__main__":
    main()
