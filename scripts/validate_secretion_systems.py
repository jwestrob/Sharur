#!/usr/bin/env python3
"""
Validate TXSScan HMM hits via MacSyFinder co-localization.

Takes the global Astra TXSScan output (with --write_macsyfinder), splits
the hmmsearch results per-genome, and runs MacSyFinder --previous-run on each
genome to check gene co-localization for secretion system identification.

Usage:
    # Standalone
    python scripts/validate_secretion_systems.py data/DATASET/ [--workers 4]

    # As library (called from 07_build_knowledge_base.py)
    from validate_secretion_systems import validate_secretion_systems
    systems_df, genes_df = validate_secretion_systems(db_path, data_dir)
"""

from __future__ import annotations

import argparse
import logging
import shutil
import subprocess
import time as _time
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import duckdb
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


# ------------------------------------------------------------------ #
# Reuse protein→genome mapping from validate_defense_systems
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
# Global MacSyFinder run: one invocation, all genomes
# ------------------------------------------------------------------ #

def _prepare_global_macsyfinder_dir(
    macsyfinder_compat: Path,
    genome_faas: dict[str, Path],
    work_dir: Path,
) -> Path:
    """Create a MacSyFinder --previous-run directory for the whole dataset.

    Creates a combined FASTA and conf file. The existing hmmsearch results
    from Astra are symlinked in — no copying or splitting needed.

    Returns the prepared directory path.
    """
    t0 = _time.time()

    run_dir = work_dir / "global_run"
    run_dir.mkdir(parents=True, exist_ok=True)

    # Symlink hmmer_results from the Astra output (already global)
    hmmer_link = run_dir / "hmmer_results"
    if not hmmer_link.exists():
        hmmer_link.symlink_to((macsyfinder_compat / "hmmer_results").resolve())

    # Create combined FASTA preserving per-genome protein order
    combined_faa = run_dir / "combined_proteins.faa"
    with open(combined_faa, "w") as fout:
        for genome_name in sorted(genome_faas):
            faa = genome_faas[genome_name]
            with open(faa) as fin:
                for line in fin:
                    fout.write(line)

    # Write macsyfinder.conf — ordered_replicon treats each contig independently
    conf_path = run_dir / "macsyfinder.conf"
    with open(conf_path, "w") as fh:
        fh.write("[base]\n")
        fh.write(f"sequence_db = {combined_faa.resolve()}\n")
        fh.write("db_type = ordered_replicon\n")
        fh.write("hmmer = hmmsearch\n\n")
        fh.write("[hmmer]\n")
        fh.write("e_value_search = 0.1\n")

    logger.info(f"  Prepared global MacSyFinder dir in {_time.time() - t0:.1f}s")
    return run_dir


# ------------------------------------------------------------------ #
# Run MacSyFinder --previous-run with TXSScan models
# ------------------------------------------------------------------ #

def _find_models_dir() -> Optional[Path]:
    """Locate the MacSyFinder models directory containing TXSScan."""
    candidates = [
        Path.home() / ".macsyfinder" / "models",
        Path.home() / ".mdmlab" / "macsyfinder" / "models",
    ]
    for c in candidates:
        if (c / "TXSScan").is_dir():
            return c
    return None


def _run_macsyfinder_on_dir(
    run_dir: Path, models_dir: Path, workers: int = 4
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run MacSyFinder --previous-run on a prepared directory.

    CRITICAL: Do NOT pass -o flag. With --previous-run, MacSyFinder only reuses
    existing hmmsearch results if -o is omitted. With -o, it starts a fresh run
    and re-runs all hmmsearch from scratch (hours on large datasets).

    Returns (systems_df, genes_df).
    """
    run_dir_abs = run_dir.resolve()
    cmd = [
        "macsyfinder",
        "--previous-run", str(run_dir_abs),
        "--models-dir", str(models_dir),
        "--models", "TXSScan", "all",
        "-w", str(workers),
        # NO -o flag — critical for reusing hmmsearch results
    ]

    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True, timeout=1800,
            cwd=str(run_dir_abs),
        )
    except subprocess.TimeoutExpired:
        logger.warning(f"  MacSyFinder timed out for {run_dir.name}")
        return pd.DataFrame(), pd.DataFrame()
    except FileNotFoundError:
        logger.error("macsyfinder not found on PATH")
        return pd.DataFrame(), pd.DataFrame()

    if proc.returncode != 0:
        stderr = proc.stderr or ""
        if "error" in stderr.lower() and "no systems found" not in stderr.lower():
            logger.warning(
                f"  MacSyFinder non-zero exit for {run_dir.name}: "
                f"{stderr[-300:]}"
            )

    # Find the output directory — MacSyFinder creates macsyfinder-YYYYMMDD_HHMMSS/
    results_dir = None
    for bs in sorted(run_dir.rglob("best_solution.tsv"), key=lambda p: p.stat().st_mtime, reverse=True):
        results_dir = bs.parent
        break

    if results_dir is None:
        for candidate in sorted(run_dir.glob("macsyfinder-*"), reverse=True):
            if candidate.is_dir():
                results_dir = candidate
                break

    if results_dir is None:
        return pd.DataFrame(), pd.DataFrame()

    best_solution = results_dir / "best_solution.tsv"
    if not best_solution.exists():
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

    genes_df = df.copy()

    # Extract system type from model_fqn (e.g. "TXSScan/bacteria/diderm/T2SS")
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

    Example: 'TXSScan/bacteria/diderm/T2SS' → 'T2SS'
             'TXSScan/bacteria/monoderm/ComM' → 'ComM'
    """
    if not isinstance(model_fqn, str):
        return str(model_fqn)
    parts = model_fqn.rstrip("/").split("/")
    return parts[-1] if len(parts) >= 2 else model_fqn


def _extract_system_subtype(model_fqn: str) -> str:
    """Extract system subtype — includes membrane type context.

    Example: 'TXSScan/bacteria/diderm/T2SS' → 'diderm/T2SS'
    """
    if not isinstance(model_fqn, str):
        return str(model_fqn)
    parts = model_fqn.rstrip("/").split("/")
    if len(parts) >= 3:
        return "/".join(parts[-2:])
    return parts[-1] if parts else model_fqn


# ------------------------------------------------------------------ #
# Main validation entry point
# ------------------------------------------------------------------ #

def validate_secretion_systems(
    db_path: str | Path,
    data_dir: str | Path | None = None,
    workers: int = 4,
    verbose: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run MacSyFinder co-localization validation on Astra TXSScan output.

    Args:
        db_path: Path to sharur.duckdb.
        data_dir: Dataset root directory (contains stage03_prodigal/, stage04_astra/).
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

    # Locate Astra TXSScan macsyfinder_compat output
    macsyfinder_compat = data_dir / "stage04_astra" / "txsscan_results" / "macsyfinder_compat"
    if not macsyfinder_compat.is_dir() or not list(macsyfinder_compat.glob("hmmer_results/*.search_hmm.out")):
        logger.error(
            "No TXSScan macsyfinder_compat directory found. "
            "Re-run stage 04 with TXSScan --write_macsyfinder."
        )
        return pd.DataFrame(), pd.DataFrame()

    global_hmmer_dir = macsyfinder_compat / "hmmer_results"
    n_hmm_files = len(list(global_hmmer_dir.glob("*.search_hmm.out")))
    if verbose:
        logger.info(f"Found {n_hmm_files} TXSScan HMM result files in {global_hmmer_dir}")

    # Find MacSyFinder models
    models_dir = _find_models_dir()
    if models_dir is None:
        logger.error(
            "TXSScan models not found. Install with: "
            "macsyfinder --install-models TXSScan"
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

    # Create temp working directory
    work_dir = data_dir / "stage04_astra" / "txsscan_results" / "_macsyfinder_work"
    if work_dir.exists():
        shutil.rmtree(work_dir)
    work_dir.mkdir(parents=True)

    # Prepare a single MacSyFinder run directory with combined FASTA + symlinked hmmer_results
    run_dir = _prepare_global_macsyfinder_dir(
        macsyfinder_compat, genome_faas, work_dir
    )

    # Run MacSyFinder ONCE on the entire dataset
    if verbose:
        logger.info(
            f"  Running MacSyFinder once on {len(genome_faas)} genomes "
            f"(ordered_replicon mode, {workers} workers)..."
        )

    t0 = _time.time()
    merged_systems, merged_genes = _run_macsyfinder_on_dir(
        run_dir, models_dir, workers=workers
    )
    elapsed = _time.time() - t0
    if verbose:
        logger.info(f"  MacSyFinder completed in {elapsed:.1f}s")

    # Post-process: add genome_id column by mapping protein IDs back to genomes
    if not merged_genes.empty and "hit_id" in merged_genes.columns:
        merged_genes["genome_id"] = merged_genes["hit_id"].map(prot_to_genome)
    if not merged_systems.empty:
        if "protein_in_syst" in merged_systems.columns:
            merged_systems["genome_id"] = merged_systems["protein_in_syst"].apply(
                lambda pids: prot_to_genome.get(
                    str(pids).split(",")[0].strip(), ""
                ) if pd.notna(pids) else ""
            )

    n_genomes_with_systems = (
        len(merged_systems["genome_id"].unique()) if not merged_systems.empty else 0
    )
    if verbose:
        logger.info(
            f"\nTotal: {len(merged_systems)} secretion systems, "
            f"{len(merged_genes)} gene assignments across {n_genomes_with_systems} genomes"
        )

    # Save merged results
    results_out = data_dir / "stage04_astra" / "txsscan_results"
    if not merged_systems.empty:
        merged_systems.to_csv(
            results_out / "validated_secretion_systems.tsv", sep="\t", index=False
        )
    if not merged_genes.empty:
        merged_genes.to_csv(
            results_out / "validated_secretion_genes.tsv", sep="\t", index=False
        )

    # Cleanup work dir
    shutil.rmtree(work_dir, ignore_errors=True)

    return merged_systems, merged_genes


# ------------------------------------------------------------------ #
# DB integration
# ------------------------------------------------------------------ #

def integrate_results(
    db_path: Path,
    all_systems: pd.DataFrame,
    all_genes: pd.DataFrame,
) -> None:
    """Load secretion system validation results into the Sharur database."""
    conn = duckdb.connect(str(db_path))

    # Create secretion_systems table if needed
    tables = [r[0] for r in conn.execute("SHOW TABLES").fetchall()]
    if "secretion_systems" not in tables:
        conn.execute("""
            CREATE TABLE secretion_systems (
                system_id VARCHAR PRIMARY KEY,
                genome_id VARCHAR,
                system_type VARCHAR,
                system_subtype VARCHAR,
                genes_count INTEGER,
                protein_ids VARCHAR,
                profile_names VARCHAR,
                sys_beg VARCHAR,
                sys_end VARCHAR,
                created_at TIMESTAMP
            )
        """)
        logger.info("Created secretion_systems table")
    else:
        before_count = conn.execute("SELECT COUNT(*) FROM secretion_systems").fetchone()[0]
        conn.execute("DELETE FROM secretion_systems")
        logger.info(f"Cleared {before_count} existing secretion_systems rows")

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
                "genes_count": int(row.get("genes_count", 0)),
                "protein_ids": row.get("protein_in_syst", ""),
                "profile_names": row.get("name_of_profiles_in_sys", ""),
                "sys_beg": str(row.get("sys_beg", "")),
                "sys_end": str(row.get("sys_end", "")),
                "created_at": now,
            })

        systems_insert = pd.DataFrame(rows)
        conn.register("tmp_systems", systems_insert)
        conn.execute("INSERT INTO secretion_systems SELECT * FROM tmp_systems")
        conn.unregister("tmp_systems")
        logger.info(f"Inserted {len(systems_insert)} secretion systems")

    # Add source='txsscan_system' annotations
    if not all_genes.empty and "hit_id" in all_genes.columns:
        existing = conn.execute(
            "SELECT COUNT(*) FROM annotations WHERE source = 'txsscan_system'"
        ).fetchone()[0]
        if existing > 0:
            conn.execute("DELETE FROM annotations WHERE source = 'txsscan_system'")
            logger.info(f"Cleared {existing} existing txsscan_system annotations")

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
                "source": "txsscan_system",
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
            logger.info(f"Inserted {len(ann_df)} txsscan_system annotations")

    # Summary comparison
    astra_count = conn.execute(
        "SELECT COUNT(DISTINCT protein_id) FROM annotations WHERE source = 'txsscan'"
    ).fetchone()[0]
    system_count = conn.execute(
        "SELECT COUNT(DISTINCT protein_id) FROM annotations WHERE source = 'txsscan_system'"
    ).fetchone()[0]
    logger.info(f"\nComparison:")
    logger.info(f"  Astra HMM-only (source='txsscan'): {astra_count} unique proteins")
    logger.info(f"  System-validated (source='txsscan_system'): {system_count} unique proteins")
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
        description="Validate TXSScan HMM hits via MacSyFinder co-localization"
    )
    parser.add_argument(
        "data_dir", type=Path,
        help="Dataset directory (e.g., data/DATASET/)"
    )
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument(
        "--no-integrate", action="store_true",
        help="Skip database integration (just produce TSV output)"
    )
    args = parser.parse_args()

    db_path = args.data_dir / "sharur.duckdb"

    systems_df, genes_df = validate_secretion_systems(
        db_path=db_path,
        data_dir=args.data_dir,
        workers=args.workers,
    )

    if not args.no_integrate and db_path.exists():
        integrate_results(db_path, systems_df, genes_df)

    logger.info("Done.")


if __name__ == "__main__":
    main()
