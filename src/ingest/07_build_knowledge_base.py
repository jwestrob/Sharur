#!/usr/bin/env python3
"""
Stage 07: Build Sharur knowledge base (DuckDB) from ingest outputs.

Lightweight implementation:
- Uses Sharur DuckDB schema from sharur.storage.schema.SCHEMA
- Ingests bins (stage02), proteins/contigs (stage03), annotations (stage04/dbCAN),
  loci (stage05a/05c when present), and basic feature_store metrics.
- Skips ELSA/vector store work per current deferral.
"""

from __future__ import annotations

import json
import logging
import os
import re
import shutil
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import duckdb
import pandas as pd
import typer
from rich.console import Console

from sharur.storage.migrations import run_migrations
from sharur.storage.schema import SCHEMA

console = Console()
logger = logging.getLogger(__name__)

# --------------------------------------------------------------------------- #
# Default e-value thresholds applied at load time.
# Databases with --cut_ga (PFAM, KOFAM, HydDB) are already clean; these
# thresholds act as a safety net.  Databases WITHOUT --cut_ga
# (DefenseFinder, VOGdb, CANT-HYD) rely on this filter.
# --------------------------------------------------------------------------- #
DEFAULT_EVALUE_THRESHOLDS: Dict[str, float] = {
    "pfam": 1e-5,            # Already clean via --cut_ga (max ~3e-4)
    "kegg": 1e-5,            # --cut_ga --cascade; this is a safety net only
    "hyddb": 1e-5,           # Already clean via --cut_ga
    "defensefinder": 1e-15,  # Has GA but values are permissive; superfamily HMMs need strict e-value
    "vogdb": 1e-15,          # No GA thresholds (0/48,439 profiles); e-value filter required
    "cant_hyd": 1e-10,       # No --cut_ga
    "cazy": 1e-5,            # dbCAN has own thresholds
    "txsscan": 1e-10,        # No GA thresholds; secretion system HMMs
}


def _resolve_thread_count(requested: Optional[int] = None) -> int:
    """Resolve the worker count, honoring a Slurm allocation when present."""
    if requested is not None:
        threads = requested
        source = "--threads"
    elif os.environ.get("SLURM_CPUS_ON_NODE"):
        source = "SLURM_CPUS_ON_NODE"
        try:
            threads = int(os.environ[source])
        except ValueError as exc:
            raise ValueError(f"{source} must be a positive integer") from exc
    else:
        threads = os.cpu_count() or 1
        source = "os.cpu_count()"

    if threads < 1:
        raise ValueError(f"Thread count from {source} must be positive (got {threads})")
    return threads

# Regex for DefenseFinder accessions with numbered FAM variants.
# Matches trailing _FAM_N or _FAMN (where N is one or more digits) at the
# end of an accession string.  Example matches:
#   RM_Type_II__Type_II_MTases_FAM_1  → group prefix = RM_Type_II__Type_II_MTases
#   RM_Type_III__Type_III_MTases_FAM10 → group prefix = RM_Type_III__Type_III_MTases
# Does NOT match component-level suffixes like MkoA_B (no leading _FAM).
_DEFENSEFINDER_FAM_RE = re.compile(r"_FAM_?\d+$")


def _dedup_defensefinder_fam_variants(df: "pd.DataFrame") -> "pd.DataFrame":
    """Deduplicate DefenseFinder hits that differ only by FAM variant number.

    DefenseFinder has ~40+ Type_II_MTases_FAM_* HMM profiles that all detect
    variations of the same Rossmann methyltransferase fold.  When a single
    protein is scanned, it matches many FAMs simultaneously, inflating hit
    counts (e.g. 27 unique MTase proteins → 259 annotation rows).

    This function groups hits by (protein_id, accession-with-FAM-suffix-stripped)
    and keeps only the best-scoring hit per group.  Accessions that do NOT end
    in a _FAM_N pattern are left untouched.

    Args:
        df: DataFrame with at least columns 'protein_id', 'accession', 'score', 'evalue'.

    Returns:
        DataFrame with redundant FAM variant rows removed.
    """
    if df.empty:
        return df

    # Compute the dedup group key by stripping trailing _FAM_N
    df = df.copy()
    df["_dedup_group"] = df["accession"].apply(
        lambda x: _DEFENSEFINDER_FAM_RE.sub("", x) if isinstance(x, str) else x
    )

    # Identify rows that actually have FAM variants (group key differs from accession)
    has_fam = df["_dedup_group"] != df["accession"]

    if not has_fam.any():
        df.drop(columns=["_dedup_group"], inplace=True)
        return df

    # Split: rows without FAM suffix pass through; rows with FAM suffix get deduped
    no_fam_df = df[~has_fam].drop(columns=["_dedup_group"])

    fam_df = df[has_fam].copy()

    # For each (protein_id, dedup_group), keep the row with the highest score.
    # Break ties by lowest e-value.  NaN scores sort last (ascending=False puts
    # NaN at the bottom), so real scores are preferred.
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
            f"({len(fam_df):,} → {len(fam_deduped):,} FAM rows)"
        )

    result = pd.concat([no_fam_df, fam_deduped], ignore_index=True)
    return result


def _dedup_kegg_best_hit(df: "pd.DataFrame") -> "pd.DataFrame":
    """Keep only the best-scoring KEGG KO per protein.

    KOFAM assigns individual score thresholds per KO profile.  Proteins
    matching broad enzyme superfamilies (PKS, P450, etc.) can hit dozens
    of variant-specific KO models simultaneously.  Only the top-scoring
    KO is informative; the rest describe the same fold.

    Args:
        df: DataFrame with at least columns 'protein_id', 'score', 'evalue'.

    Returns:
        DataFrame with one row per protein (the best-scoring KO hit).
    """
    if df.empty or len(df) <= 1:
        return df

    before = len(df)
    df = df.sort_values(
        ["score", "evalue"], ascending=[False, True], na_position="last"
    )
    df = df.drop_duplicates(subset=["protein_id"], keep="first")
    n_removed = before - len(df)
    if n_removed > 0:
        logger.info(
            f"  kegg: deduplicated {n_removed:,} secondary KO hits "
            f"({before:,} → {len(df):,} rows, best-hit-per-protein)"
        )
    return df


# --------------------------------------------------------------------------- #
# Store adapter for predicate generation (matches DuckDBStore.execute interface)
# --------------------------------------------------------------------------- #
class _StoreAdapter:
    """Minimal adapter to make raw DuckDB connection work with predicate generator."""

    def __init__(self, conn: duckdb.DuckDBPyConnection):
        self.conn = conn

    def execute(self, query: str, params=None) -> list:
        if params:
            return self.conn.execute(query, params).fetchall()
        return self.conn.execute(query).fetchall()


# --------------------------------------------------------------------------- #
# Data containers
# --------------------------------------------------------------------------- #
@dataclass
class PipelineOutputs:
    stage00_dir: Path
    stage01_dir: Path
    stage02_dir: Path
    stage03_dir: Path
    stage04_dir: Path
    stage05a_dir: Path
    stage05b_dir: Path
    stage05c_dir: Path
    stage06_dir: Path

    def validate(self) -> List[str]:
        present = []
        for name, path in self.__dict__.items():
            if path.exists():
                present.append(name)
        return present


# --------------------------------------------------------------------------- #
# Builder
# --------------------------------------------------------------------------- #
class KnowledgeBaseBuilder:
    def __init__(self, outputs: PipelineOutputs, db_path: Path, force: bool = False,
                 enable_cazymes: bool = False, threads: Optional[int] = None):
        self.outputs = outputs
        self.db_path = db_path
        self.force = force
        self.enable_cazymes = enable_cazymes
        self.threads = _resolve_thread_count(threads)
        self.conn: Optional[duckdb.DuckDBPyConnection] = None
        self.embeddings_path: Optional[str] = None
        self.stats: Dict[str, int] = {
            "bins": 0,
            "contigs": 0,
            "proteins": 0,
            "annotations": 0,
            "loci": 0,
            "embeddings": 0,
            "predicates": 0,
            "semantic_atoms": 0,
            "semantic_state": 0,
            "review_queue": 0,
        }
        # reference maps for annotation names
        self.ref_dir = Path(__file__).resolve().parents[2] / "data" / "reference"
        self.pfam_id_to_acc: Dict[str, str] = {}
        self.pfam_acc_to_name: Dict[str, str] = {}
        self.pfam_acc_to_desc: Dict[str, str] = {}
        self.kegg_def: Dict[str, str] = {}
        self.kegg_def_short: Dict[str, str] = {}
        self._load_reference_maps()

    def _load_reference_maps(self) -> None:
        # PFAM mapping: id -> accession, accession -> name/description
        pfam_path = self.ref_dir / "pfam_id_desc.tsv"
        if pfam_path.exists():
            try:
                pdf = pd.read_csv(pfam_path, sep="\t", header=None, names=["accession", "id", "description"])
                self.pfam_id_to_acc = {row.id: row.accession for row in pdf.itertuples()}
                self.pfam_acc_to_name = {row.accession: row.id for row in pdf.itertuples()}
                self.pfam_acc_to_desc = {row.accession: row.description for row in pdf.itertuples()}
            except Exception as exc:
                logger.warning(f"Failed to load PFAM reference ({pfam_path}): {exc}")
        # KEGG KOFAM mapping
        ko_path = self.ref_dir / "ko_list"
        if ko_path.exists():
            try:
                kdf = pd.read_csv(ko_path, sep="\t")
                self.kegg_def = dict(zip(kdf["knum"], kdf["definition"]))
                self.kegg_def_short = dict(zip(kdf["knum"], kdf.get("simplified_definition", kdf["definition"])))
            except Exception as exc:
                logger.warning(f"Failed to load KO reference ({ko_path}): {exc}")

    # --- orchestrate ---------------------------------------------------- #
    def _run_step(self, name: str, func, step_num: int, total: int) -> None:
        """Run a build step with timing and logging."""
        console.print(f"\n[bold cyan]━━━ Step {step_num}/{total}: {name} ━━━[/bold cyan]")
        t0 = time.time()
        try:
            func()
            elapsed = time.time() - t0
            console.print(f"  [green]✓ {name} completed in {elapsed:.1f}s[/green]")
        except Exception as e:
            elapsed = time.time() - t0
            console.print(f"  [red]✗ {name} FAILED after {elapsed:.1f}s: {e}[/red]")
            raise

    def build(self) -> Dict[str, int]:
        build_start = time.time()
        self._init_db()

        steps = [
            ("Load bins", self._load_bins),
            ("Load proteins/contigs", self._load_proteins_and_contigs),
            ("Load annotations", self._load_annotations),
            ("Load loci", self._load_loci),
            ("Load embeddings info", self._load_embeddings_info),
            ("Compute metrics", self._compute_basic_metrics),
            ("Classify hydrogenases", self._classify_hydrogenases),
            *([("Classify CAZymes", self._classify_cazymes)] if self.enable_cazymes else []),
            ("Validate defense systems", self._validate_defense_systems),
            ("Validate secretion systems", self._validate_secretion_systems),
            ("Generate predicates", self._generate_predicates),
            ("Finalize indexes", self._finalize),
        ]

        for i, (name, func) in enumerate(steps, 1):
            self._run_step(name, func, i, len(steps))

        total_elapsed = time.time() - build_start
        console.print(f"\n[bold green]Build complete in {total_elapsed:.0f}s ({total_elapsed/60:.1f}m)[/bold green]")
        console.print(f"Stats: {self.stats}")
        return self.stats

    def resume_v2(self) -> Dict[str, int]:
        """Continue an interrupted V2 generation against an existing DB."""
        if not self.db_path.exists():
            raise FileNotFoundError(f"{self.db_path} is required for --resume-v2")
        self._reacquire_db()
        resume_start = time.time()
        self._run_step(
            "Resume predicates",
            lambda: self._generate_predicates(resume=True),
            1,
            2,
        )
        self._run_step("Finalize indexes", self._finalize, 2, 2)
        elapsed = time.time() - resume_start
        console.print(
            f"\n[bold green]V2 resume complete in {elapsed:.0f}s "
            f"({elapsed / 60:.1f}m)[/bold green]"
        )
        console.print(f"Stats: {self.stats}")
        return self.stats

    def _finalize(self) -> None:
        self._create_indexes()
        self._update_stats()

    # --- connection management ------------------------------------------ #
    def _release_db(self) -> None:
        """Close the DB connection so external scripts can open their own."""
        if self.conn is not None:
            self.conn.commit()
            self.conn.close()
            self.conn = None

    def _reacquire_db(self) -> None:
        """Re-open the DB connection after external scripts finish."""
        self.conn = duckdb.connect(str(self.db_path))
        self.conn.execute(f"SET threads = {self.threads}")

    # --- init ----------------------------------------------------------- #
    def _init_db(self) -> None:
        if self.db_path.exists():
            if self.force:
                self.db_path.unlink()
            else:
                raise FileExistsError(f"{self.db_path} exists (use --force to overwrite)")

        self.conn = duckdb.connect(str(self.db_path))
        self.conn.execute(f"SET threads = {self.threads}")
        self.conn.execute(SCHEMA)
        run_migrations(self.conn)
        console.print(f"[blue]Created DuckDB at {self.db_path}[/blue]")

    # --- bins ----------------------------------------------------------- #
    def _load_bins(self) -> None:
        manifest = self._load_manifest(self.outputs.stage02_dir / "processing_manifest.json")
        records = []
        for genome in manifest.get("genomes", []) if manifest else []:
            tax = genome.get("taxonomy", {})
            records.append(
                {
                    "bin_id": genome.get("genome_id"),
                    "completeness": tax.get("completeness"),
                    "contamination": tax.get("contamination"),
                    "taxonomy": tax.get("name") or tax.get("taxonomy") or "unknown",
                    "n_contigs": genome.get("n_contigs") or 0,
                    "total_length": genome.get("total_length") or 0,
                }
            )
        if records:
            df = pd.DataFrame(records)
            for col in ["n_contigs", "total_length"]:
                if col not in df.columns:
                    df[col] = 0
            self.conn.execute("INSERT INTO bins SELECT * FROM df")
            self.stats["bins"] = len(records)

    # --- proteins/contigs ---------------------------------------------- #
    def _load_proteins_and_contigs(self) -> None:
        genomes_dir = self.outputs.stage03_dir / "genomes"
        protein_rows = []
        contig_lengths: Dict[str, int] = {}
        if not genomes_dir.exists():
            logger.warning("stage03 genomes directory missing; skipping proteins")
            return

        faa_files = [
            faa for faa in genomes_dir.glob("**/*.faa")
            if faa.parent.name != "all_protein_symlinks"
        ]
        console.print(
            f"  Parsing {len(faa_files):,} FAA files with {self.threads} threads..."
        )
        t0 = time.time()

        def parse_faa(faa: Path):
            local_contig_lengths: Dict[tuple, int] = {}
            bin_id = faa.parent.name
            rows = self._parse_prodigal_faa(faa, bin_id, local_contig_lengths)
            return rows, local_contig_lengths

        # Bound outstanding results to one full worker wave. This keeps file
        # order deterministic and prevents completed parse results from
        # accumulating behind a slower early file.
        with ThreadPoolExecutor(max_workers=self.threads) as executor:
            for batch_start in range(0, len(faa_files), self.threads):
                batch = faa_files[batch_start : batch_start + self.threads]
                futures = [executor.submit(parse_faa, faa) for faa in batch]
                for future in futures:
                    rows, local_contig_lengths = future.result()
                    protein_rows.extend(rows)
                    for key, length in local_contig_lengths.items():
                        contig_lengths[key] = max(contig_lengths.get(key, 0), length)

                processed = batch_start + len(batch)
                if processed % 500 == 0 or processed == len(faa_files):
                    console.print(
                        f"    {processed:,}/{len(faa_files):,} files, "
                        f"{len(protein_rows):,} proteins so far"
                    )
        console.print(f"  Parsed {len(protein_rows):,} proteins from {len(faa_files):,} files ({time.time()-t0:.1f}s)")

        # Ensure bins exist even if stage02 was skipped
        if contig_lengths:
            bin_stats: Dict[str, Dict[str, int]] = {}
            for (contig_id, bid), length in contig_lengths.items():
                stats = bin_stats.setdefault(bid, {"n_contigs": 0, "total_length": 0})
                stats["n_contigs"] += 1
                stats["total_length"] += length
            missing_bins = []
            for bid, stats in bin_stats.items():
                exists = self.conn.execute("SELECT 1 FROM bins WHERE bin_id = ?", [bid]).fetchone()
                if not exists:
                    missing_bins.append(
                        {
                            "bin_id": bid,
                            "completeness": None,
                            "contamination": None,
                            "taxonomy": "unknown",
                            "n_contigs": stats["n_contigs"],
                            "total_length": stats["total_length"],
                        }
                    )
            if missing_bins:
                bdf = pd.DataFrame(missing_bins)
                self.conn.register("tmp_bins", bdf)
                self.conn.execute(
                    "INSERT INTO bins (bin_id, completeness, contamination, taxonomy, n_contigs, total_length) SELECT * FROM tmp_bins"
                )
                self.stats["bins"] += len(missing_bins)

        if contig_lengths:
            seen: set[str] = set()
            contigs = []
            for (cid, bid), length in contig_lengths.items():
                if cid in seen:
                    continue
                seen.add(cid)
                contigs.append(
                    {
                        "contig_id": cid,
                        "bin_id": bid,
                        "length": length,
                        "gc_content": None,
                        "is_circular": False,
                        "taxonomy": None,
                    }
                )
            cdf = pd.DataFrame(contigs)
            cdf = cdf.reindex(
                columns=["contig_id", "bin_id", "length", "gc_content", "is_circular", "taxonomy"],
                fill_value=None,
            )
            self.conn.register("tmp_cdf", cdf)
            self.conn.execute(
                """
                INSERT INTO contigs (contig_id, bin_id, length, gc_content, is_circular, taxonomy)
                SELECT * FROM tmp_cdf
                """,
            )
            self.stats["contigs"] = len(contigs)

        if protein_rows:
            pdf = pd.DataFrame(protein_rows)
            # Deduplicate proteins that appear in multiple bins (same contig
            # assigned to different bins by different binners). Keep first
            # occurrence; the bin_id of the first-seen bin is retained.
            before_dedup = len(pdf)
            pdf = pdf.drop_duplicates(subset=["protein_id"], keep="first")
            if len(pdf) < before_dedup:
                console.print(
                    f"  [yellow]Deduplicated {before_dedup - len(pdf):,} "
                    f"shared proteins across bins ({len(pdf):,} unique)[/yellow]"
                )
            # Ensure column order matches schema
            pdf = pdf.reindex(
                columns=[
                    "protein_id",
                    "contig_id",
                    "bin_id",
                    "start",
                    "end_coord",
                    "strand",
                    "gene_index",
                    "sequence",
                    "sequence_length",
                    "gc_content",
                ],
                fill_value=None,
            )
            self.conn.register("tmp_pdf", pdf)
            self.conn.execute(
                """
                INSERT INTO proteins (
                    protein_id, contig_id, bin_id, start, end_coord, strand,
                    gene_index, sequence, sequence_length, gc_content
                )
                SELECT * FROM tmp_pdf
                """,
            )
            self.stats["proteins"] = len(protein_rows)

    def _parse_prodigal_faa(
        self, path: Path, bin_id: str, contig_lengths: Dict[tuple, int]
    ) -> List[Dict]:
        rows: List[Dict] = []
        contig_gene_idx: Dict[str, int] = {}
        with open(path) as f:
            current_id = None
            seq_parts: List[str] = []
            header_meta = None
            for line in f:
                if line.startswith(">"):
                    # flush previous
                    if current_id is not None:
                        seq = "".join(seq_parts)
                        if not seq.rstrip("*"):
                            raise ValueError(
                                f"Empty protein sequence for {current_id} in {path}"
                            )
                        row = self._row_from_header(current_id, header_meta, seq, bin_id)
                        contig_gene_idx[row["contig_id"]] = contig_gene_idx.get(row["contig_id"], 0)
                        row["gene_index"] = contig_gene_idx[row["contig_id"]]
                        contig_gene_idx[row["contig_id"]] += 1
                        rows.append(row)
                        seq_parts = []
                    header = line[1:].strip()
                    parts = header.split(" # ")
                    current_id = parts[0].split()[0]
                    header_meta = parts[1:] if len(parts) > 3 else None
                else:
                    seq_parts.append(line.strip())
            if current_id is not None:
                seq = "".join(seq_parts)
                if not seq.rstrip("*"):
                    raise ValueError(f"Empty protein sequence for {current_id} in {path}")
                row = self._row_from_header(current_id, header_meta, seq, bin_id)
                contig_gene_idx[row["contig_id"]] = contig_gene_idx.get(row["contig_id"], 0)
                row["gene_index"] = contig_gene_idx[row["contig_id"]]
                contig_gene_idx[row["contig_id"]] += 1
                rows.append(row)

        # track contig length
        for r in rows:
            contig_id = r["contig_id"]
            key = (contig_id, r["bin_id"])
            contig_lengths[key] = max(contig_lengths.get(key, 0), r["end_coord"])
        return rows

    def _row_from_header(self, protein_id: str, meta: Optional[List[str]], sequence: str, bin_id: str) -> Dict:
        # Meta fields: start, end, strand_flag (1/-1)
        start, end, strand = 0, len(sequence) * 3, "+"
        if meta and len(meta) >= 3:
            try:
                start = int(meta[0])
                end = int(meta[1])
                strand = "+" if meta[2].strip() == "1" else "-"
            except Exception:
                pass
        contig_id = protein_id.rsplit("_", 1)[0] if "_" in protein_id else protein_id
        gene_index = None
        return {
            "protein_id": protein_id,
            "contig_id": contig_id,
            "bin_id": bin_id,
            "start": start,
            "end_coord": end,
            "strand": strand,
            "sequence": sequence,
            "sequence_length": len(sequence) if sequence else None,
            "gene_index": gene_index,
            "gc_content": None,
        }

    # --- annotations ---------------------------------------------------- #
    def _load_annotations(self) -> None:
        # Astra PFAM/KOFAM TSV
        tsv_files = list(self.outputs.stage04_dir.rglob("*_hits_df.tsv")) if self.outputs.stage04_dir.exists() else []
        console.print(f"  Found {len(tsv_files)} annotation TSV files")
        for tsv in tsv_files:
            # Skip dbCAN raw HMM hits — _classify_cazymes() handles CAZy
            # annotation loading via 3-tool consensus pipeline
            if "dbcan" in tsv.as_posix().lower():
                logger.info(f"Skipping {tsv.name} (handled by CAZyme consensus pipeline)")
                continue
            try:
                t_ann = time.time()
                df = pd.read_csv(tsv, sep="\t")
                # Derive source name from filename/path (always lowercase)
                tsv_lower = tsv.as_posix().lower()
                if "pfam" in tsv_lower:
                    source = "pfam"
                elif "kofam" in tsv_lower or "kegg" in tsv_lower:
                    source = "kegg"
                elif "hyddb" in tsv_lower:
                    source = "hyddb"
                elif "defense" in tsv_lower:
                    source = "defensefinder"
                elif "vog" in tsv_lower:
                    source = "vogdb"
                elif "cant" in tsv_lower:
                    source = "cant_hyd"
                elif "txss" in tsv_lower:
                    source = "txsscan"
                else:
                    source = tsv.stem.split("_")[0].lower()
                original_name = df.get("hmm_name", None)
                # Normalize columns and avoid duplicate labels
                df = df.rename(
                    columns={
                        "sequence_id": "protein_id",
                        "hmm_name": "accession",
                        "evalue": "evalue",
                        "i_evalue": "i_evalue",
                        "c_evalue": "c_evalue",
                        "bit_score": "score",
                        "bitscore": "score_primary",
                        "dom_bitscore": "dom_bitscore",
                        "ali_from": "start_aa",
                        "ali_to": "end_aa",
                        "env_from": "env_from",
                        "env_to": "env_to",
                    }
                )
                # Map names/descriptions from reference
                if source == "pfam" and original_name is not None:
                    df["accession"] = df["accession"].map(self.pfam_id_to_acc).fillna(df["accession"])
                    df["name"] = df["accession"].map(self.pfam_acc_to_name)
                    df["description"] = df["accession"].map(self.pfam_acc_to_desc)
                elif source == "kegg":
                    df["name"] = df["accession"].map(self.kegg_def_short)
                    df["description"] = df["accession"].map(self.kegg_def)
                # Fallbacks
                df["name"] = df.get("name")
                df["description"] = df.get("description")
                if original_name is not None and "name" in df:
                    df["name"] = df["name"].fillna(original_name)
                df["description"] = df["description"].fillna("")

                # Prefer inner hit e-value/score; fill from c/i/dom if missing
                if "evalue" not in df:
                    df["evalue"] = None
                if "score" not in df:
                    df["score"] = None
                if "i_evalue" in df:
                    df["evalue"] = df["evalue"].fillna(df["i_evalue"])
                if "c_evalue" in df:
                    df["evalue"] = df["evalue"].fillna(df["c_evalue"])
                if "score" in df and "dom_bitscore" in df:
                    df["score"] = df["score"].fillna(df["dom_bitscore"])
                if "score" in df and "score_primary" in df:
                    df["score"] = df["score"].fillna(df["score_primary"])
                # Coordinates: prefer ali_, fallback to env_
                if "start_aa" not in df:
                    df["start_aa"] = None
                if "end_aa" not in df:
                    df["end_aa"] = None
                if "env_from" in df:
                    df["start_aa"] = df["start_aa"].fillna(df["env_from"])
                if "env_to" in df:
                    df["end_aa"] = df["end_aa"].fillna(df["env_to"])

                df["source"] = source

                # Apply e-value threshold filtering
                evalue_threshold = DEFAULT_EVALUE_THRESHOLDS.get(source)
                if evalue_threshold is not None and "evalue" in df.columns:
                    before_filter = len(df)
                    df = df[df["evalue"].isna() | (df["evalue"] <= evalue_threshold)]
                    dropped = before_filter - len(df)
                    if dropped > 0:
                        logger.info(
                            f"  {source}: filtered {dropped:,} hits with e-value > {evalue_threshold:.0e} "
                            f"({len(df):,} retained)"
                        )

                # DefenseFinder FAM variant dedup: keep best hit per
                # (protein_id, system_component) when multiple _FAM_N
                # profiles match the same protein.
                if source == "defensefinder":
                    df = _dedup_defensefinder_fam_variants(df)

                # KEGG best-hit dedup: keep only the top-scoring KO per
                # protein.  KOFAM profiles for enzyme superfamilies (PKS,
                # P450, etc.) produce dozens of hits per protein from
                # variant-specific models that all match the same fold.
                if source == "kegg":
                    df = _dedup_kegg_best_hit(df)

                keep = [
                    "annotation_id",
                    "protein_id",
                    "source",
                    "accession",
                    "name",
                    "description",
                    "evalue",
                    "score",
                    "start_aa",
                    "end_aa",
                ]
                next_id = self.conn.execute("SELECT COALESCE(MAX(annotation_id), 0) FROM annotations").fetchone()[0]
                df = df.reindex(columns=keep, fill_value=None)
                df["annotation_id"] = range(next_id + 1, next_id + 1 + len(df))
                self.conn.register("tmp_adf", df)
                self.conn.execute(
                    """
                    INSERT INTO annotations (
                        annotation_id, protein_id, source, accession, name, description, evalue,
                        score, start_aa, end_aa
                    ) SELECT * FROM tmp_adf
                    """,
                )
                self.stats["annotations"] += len(df)
                console.print(f"    {source}: {len(df):,} annotations from {tsv.name} ({time.time()-t_ann:.1f}s)")
            except Exception as exc:
                logger.exception(f"Failed to load annotations from {tsv}")
                console.print(f"    [red]FAILED: {tsv.name}: {exc}[/red]")
                raise

        # dbCAN JSON — backward compatibility for datasets that ran dbcan_cazyme.py.
        # For new ingestions, _classify_cazymes() runs DIAMOND directly and takes precedence.
        for jf in (self.outputs.stage05b_dir.glob("*_cazyme_results.json") if self.outputs.stage05b_dir.exists() else []):
            try:
                data = json.loads(jf.read_text())
                ann = data.get("annotations", [])
                if not ann:
                    continue
                df = pd.DataFrame(ann)
                df = df.rename(columns={"cazyme_family": "accession"})
                df["source"] = "cazy"
                try:
                    existing = set(self.conn.execute("SELECT protein_id FROM proteins").fetchdf()["protein_id"])
                    df = df[df["protein_id"].isin(existing)]
                except Exception:
                    pass
                keep = [
                    "annotation_id",
                    "protein_id",
                    "source",
                    "accession",
                    "name",
                    "description",
                    "evalue",
                    "score",
                    "start_aa",
                    "end_aa",
                ]
                next_id = self.conn.execute("SELECT COALESCE(MAX(annotation_id), 0) FROM annotations").fetchone()[0]
                df = df.reindex(columns=keep, fill_value=None)
                df["annotation_id"] = range(next_id + 1, next_id + 1 + len(df))
                self.conn.register("tmp_adf2", df)
                self.conn.execute(
                    """
                    INSERT INTO annotations (
                        annotation_id, protein_id, source, accession, name, description, evalue,
                        score, start_aa, end_aa
                    ) SELECT * FROM tmp_adf2
                    """,
                )
                self.stats["annotations"] += len(df)
            except Exception as exc:
                logger.warning(f"Failed to load CAZy annotations from {jf}: {exc}")

    # --- loci ----------------------------------------------------------- #
    def _load_loci(self) -> None:
        gecco_json = self.outputs.stage05a_dir / "combined_bgc_data.json"
        if gecco_json.exists():
            try:
                clusters = json.loads(gecco_json.read_text()).get("clusters", [])
                loci_rows = []
                lp_rows = []
                for cl in clusters:
                    loci_rows.append(
                        {
                            "locus_id": cl.get("cluster_id"),
                            "locus_type": "bgc",
                            "contig_id": cl.get("contig"),
                            "start": cl.get("start", 0),
                            "end_coord": cl.get("end", 0),
                            "confidence": 1.0,
                            "metadata": json.dumps({"bgc_type": cl.get("bgc_type")}),
                        }
                    )
                    for idx, pid in enumerate(cl.get("protein_list") or []):
                        lp_rows.append({"locus_id": cl.get("cluster_id"), "protein_id": pid, "position": idx})
                if loci_rows:
                    ldf = pd.DataFrame(loci_rows).reindex(
                        columns=[
                            "locus_id",
                            "locus_type",
                            "contig_id",
                            "start",
                            "end_coord",
                            "confidence",
                            "metadata",
                        ],
                        fill_value=None,
                    )
                    self.conn.register("tmp_ldf", ldf)
                    self.conn.execute(
                        """
                        INSERT INTO loci (
                            locus_id, locus_type, contig_id, start, end_coord, confidence,
                            metadata
                        ) SELECT * FROM tmp_ldf
                        """,
                    )
                    self.stats["loci"] += len(loci_rows)
                if lp_rows:
                    lpdf = pd.DataFrame(lp_rows).reindex(
                        columns=["locus_id", "protein_id", "position"], fill_value=None
                    )
                    self.conn.register("tmp_lpdf", lpdf)
                    self.conn.execute("INSERT INTO locus_proteins SELECT * FROM tmp_lpdf")
            except Exception as exc:
                logger.warning(f"Failed to load GECCO loci: {exc}")

        for crispr_json in (self.outputs.stage05c_dir.glob("*_crispr_arrays.json") if self.outputs.stage05c_dir.exists() else []):
            try:
                arrays = json.loads(crispr_json.read_text()).get("arrays", [])
                loci_rows = []
                for arr in arrays:
                    loci_rows.append(
                        {
                            "locus_id": arr.get("id"),
                            "locus_type": "crispr",
                            "contig_id": arr.get("contig"),
                            "start": arr.get("startCoordinate", 0),
                            "end_coord": arr.get("endCoordinate", 0),
                            "confidence": 1.0,
                            "metadata": json.dumps(arr),
                        }
                    )
                if loci_rows:
                    ldf = pd.DataFrame(loci_rows).reindex(
                        columns=[
                            "locus_id",
                            "locus_type",
                            "contig_id",
                            "start",
                            "end_coord",
                            "confidence",
                            "metadata",
                        ],
                        fill_value=None,
                    )
                    self.conn.register("tmp_ldf_cr", ldf)
                    # Filter to contigs that exist in the DB (CRISPR arrays
                    # may land on contigs with no predicted genes)
                    n_before = len(ldf)
                    self.conn.execute(
                        """
                        INSERT INTO loci (
                            locus_id, locus_type, contig_id, start, end_coord, confidence,
                            metadata
                        ) SELECT * FROM tmp_ldf_cr
                        WHERE contig_id IN (SELECT contig_id FROM contigs)
                        """,
                    )
                    actual = self.conn.execute(
                        "SELECT COUNT(*) FROM loci WHERE locus_type = 'crispr'"
                    ).fetchone()[0]
                    skipped = n_before - actual if actual < n_before else 0
                    if skipped > 0:
                        logger.info(f"  crispr: skipped {skipped}/{n_before} arrays on contigs without genes")
                    self.stats["loci"] += actual
            except Exception as exc:
                logger.warning(f"Failed to load CRISPR arrays from {crispr_json}: {exc}")

    # --- metrics -------------------------------------------------------- #
    def _compute_basic_metrics(self) -> None:
        console.print("  Computing length z-scores...")
        self.conn.execute(
            """
            INSERT INTO feature_store (protein_id, metric_name, metric_value)
            SELECT protein_id, 'length_zscore',
                   CASE WHEN std > 0 THEN (len - mean)/std ELSE 0 END
            FROM (
                SELECT p.protein_id,
                       (p.end_coord - p.start) AS len,
                       stats.mean,
                       stats.std
                FROM proteins p
                LEFT JOIN (
                    SELECT bin_id, AVG(end_coord - start) AS mean, STDDEV(end_coord - start) AS std
                    FROM proteins
                    GROUP BY bin_id
                ) stats ON p.bin_id = stats.bin_id
            );
            """
        )

    # --- predicates ----------------------------------------------------- #
    def _generate_predicates(self, *, resume: bool = False) -> None:
        """Generate semantic atoms and V2-derived legacy predicates.

        Stage 07 treats V2 as the normal predicate product. The compatibility
        table is still materialized so existing operators and analysis scripts
        that read protein_predicates continue to work.
        """
        from sharur.predicates_v2.persistence import generate_and_persist_v2

        protein_count = self.conn.execute("SELECT COUNT(*) FROM proteins").fetchone()[0]
        if protein_count == 0:
            console.print("  No proteins found, skipping predicate generation")
            return

        reports_dir = self.db_path.parent / "reports"
        reports_dir.mkdir(parents=True, exist_ok=True)
        review_queue_path = reports_dir / "predicates_v2_review_queue.tsv"
        chunk_size = int(os.environ.get("SHARUR_V2_CHUNK_SIZE", "100000"))

        console.print(
            f"  Generating V2 atoms/states for {protein_count:,} proteins "
            f"(chunk_size={chunk_size:,}, workers={self.threads}, "
            f"resume={resume})"
        )
        generate_and_persist_v2(
            _StoreAdapter(self.conn),
            output_review_queue=str(review_queue_path),
            chunk_size=chunk_size,
            workers=self.threads,
            resume=resume,
            update_legacy_predicates=True,
            return_states=False,
            predict_topology=False,
        )

        self.stats["predicates"] = self.conn.execute(
            "SELECT COUNT(*) FROM protein_predicates"
        ).fetchone()[0]
        self.stats["semantic_atoms"] = self.conn.execute(
            "SELECT COUNT(*) FROM semantic_atoms"
        ).fetchone()[0]
        self.stats["semantic_state"] = self.conn.execute(
            "SELECT COUNT(*) FROM semantic_state"
        ).fetchone()[0]

        if review_queue_path.exists():
            # Header-only TSV is one line; queue entries are line_count - 1.
            line_count = sum(1 for _ in review_queue_path.open())
            self.stats["review_queue"] = max(0, line_count - 1)

        console.print(
            "  V2 predicates: "
            f"{self.stats['semantic_state']:,} states, "
            f"{self.stats['semantic_atoms']:,} atoms, "
            f"{self.stats['review_queue']:,} review entries"
        )

    # --- hydrogenase classification --------------------------------------- #
    def _classify_hydrogenases(self) -> None:
        """Run HydDB subgroup classification if HydDB annotations exist."""
        # Check if HydDB annotations exist
        hyddb_count = self.conn.execute(
            "SELECT COUNT(*) FROM annotations WHERE LOWER(source) = 'hyddb'"
        ).fetchone()[0]

        if hyddb_count == 0:
            console.print("  No HydDB annotations found, skipping hydrogenase classification")
            return

        console.print(f"  Found {hyddb_count} HydDB annotations, running subgroup classification...")

        # Release DB so the external script can open its own connection
        self._release_db()
        try:
            import sys
            scripts_dir = Path(__file__).resolve().parents[2] / "scripts"
            sys.path.insert(0, str(scripts_dir))

            from classify_hydrogenases import classify_hydrogenases as run_classification

            results = run_classification(
                db_path=str(self.db_path),
                threads=self.threads,
                update_predicates=True,
                verbose=False,
            )

            if not results.empty:
                n_curation = results['needs_curation'].sum() if 'needs_curation' in results.columns else 0
                console.print(f"  Classified {len(results)} hydrogenases, {n_curation} need curation")
            else:
                console.print("  No hydrogenases classified")

        except ImportError as e:
            logger.warning(f"Could not import hydrogenase classification: {e}")
            console.print("  [yellow]Hydrogenase classification skipped (missing dependencies)[/yellow]")
        except FileNotFoundError as e:
            logger.warning(f"HydDB reference not found: {e}")
            console.print("  [yellow]Hydrogenase classification skipped (HydDB reference not found)[/yellow]")
        except Exception as e:
            logger.warning(f"Hydrogenase classification failed: {e}")
            console.print(f"  [yellow]Hydrogenase classification failed: {e}[/yellow]")
        finally:
            self._reacquire_db()

    # --- cazyme classification --------------------------------------------- #
    def _classify_cazymes(self) -> None:
        """Run dbCAN 3-tool consensus CAZyme classification.

        Uses DIAMOND + HMMER (dbCAN.hmm) + HMMER (dbCAN-sub.hmm) with ≥2-tool
        consensus filter. Falls back to 2-tool mode if dbCAN-sub.hmm is absent.
        """
        protein_count = self.conn.execute("SELECT COUNT(*) FROM proteins").fetchone()[0]
        if protein_count == 0:
            console.print("  No proteins found, skipping CAZyme classification")
            return

        # Copy stage04 dbCAN HMM results to annotations dir for classify_cazymes cache
        data_dir = self.db_path.parent
        ann_dir = data_dir / "annotations"
        ann_dir.mkdir(parents=True, exist_ok=True)
        cached_tsv = ann_dir / "dbCAN_hits_df.tsv"
        if not cached_tsv.exists():
            for tsv in (self.outputs.stage04_dir.rglob("*dbCAN*_hits_df.tsv")
                        if self.outputs.stage04_dir.exists() else []):
                shutil.copy2(tsv, cached_tsv)
                console.print(f"  Cached dbCAN HMM results → {cached_tsv.name}")
                break

        console.print(f"  Running dbCAN 3-tool consensus CAZyme classification for {protein_count:,} proteins...")

        # Release DB so the external script can open its own connection
        self._release_db()
        try:
            import sys
            scripts_dir = Path(__file__).resolve().parents[2] / "scripts"
            sys.path.insert(0, str(scripts_dir))

            from classify_cazymes import classify_cazymes as run_classification

            results = run_classification(
                db_path=str(self.db_path),
                threads=self.threads,
                update_predicates=True,
                verbose=False,
            )

            if not results.empty:
                n_proteins = results['protein_id'].nunique()
                n_families = results['family_class'].nunique()
                console.print(f"  Classified {n_proteins:,} CAZymes across {n_families} families")
            else:
                console.print("  No CAZymes classified")

        except ImportError as e:
            logger.warning(f"Could not import CAZyme classification: {e}")
            console.print("  [yellow]CAZyme classification skipped (missing dependencies)[/yellow]")
        except FileNotFoundError as e:
            logger.warning(f"CAZyDB reference not found: {e}")
            console.print("  [yellow]CAZyme classification skipped (CAZyDB not found)[/yellow]")
        except Exception as e:
            logger.warning(f"CAZyme classification failed: {e}")
            console.print(f"  [yellow]CAZyme classification failed: {e}[/yellow]")
        finally:
            self._reacquire_db()

    # --- defense system validation ---------------------------------------- #
    def _validate_defense_systems(self) -> None:
        """Run co-location validation on DefenseFinder HMM hits.

        Uses the in-process co-location engine (sharur.colocation) which reads
        HMM hits directly from DuckDB annotations. No MacSyFinder subprocess
        or macsyfinder_compat directory needed.

        Skipped if no defensefinder annotations exist in the database.
        """
        self._run_colocation("defensefinder", "defense")

    # --- secretion system validation -------------------------------------- #
    def _validate_secretion_systems(self) -> None:
        """Run co-location validation on TXSScan HMM hits.

        Uses the in-process co-location engine (sharur.colocation) which reads
        HMM hits directly from DuckDB annotations. No MacSyFinder subprocess
        or macsyfinder_compat directory needed.

        Skipped if no txsscan annotations exist in the database.
        """
        self._run_colocation("txsscan", "secretion")

    def _run_colocation(self, source: str, system_type: str) -> None:
        """Generic co-location validation for any MacSyFinder-model-based source.

        Args:
            source: Annotation source in DuckDB ('defensefinder' or 'txsscan').
            system_type: Human label for logging ('defense' or 'secretion').
        """
        # Check if annotations exist for this source
        n_hits = self.conn.execute(
            "SELECT COUNT(*) FROM annotations WHERE source = ?", [source]
        ).fetchone()[0]
        if n_hits == 0:
            console.print(f"  No {source} annotations in DB, skipping {system_type} system validation")
            return

        console.print(f"  Running co-location validation on {n_hits:,} {source} hits...")

        # Release DB so the co-location engine can open its own connection
        self._release_db()
        try:
            from sharur.colocation import (
                integrate_defense_results,
                integrate_secretion_results,
                validate_systems,
            )

            systems_df, genes_df = validate_systems(
                db_path=str(self.db_path),
                source=source,
                verbose=True,
            )

            # Integration is replacement-based and must also run for an empty
            # result so stale structured calls are removed.
            if source == "defensefinder":
                integrate_defense_results(self.db_path, systems_df, genes_df)
            elif source == "txsscan":
                integrate_secretion_results(self.db_path, systems_df, genes_df)

            if not systems_df.empty or not genes_df.empty:
                n_systems = len(systems_df)
                n_genes = len(genes_df)
                console.print(f"  Validated {n_systems} {system_type} systems ({n_genes} genes)")
            else:
                console.print(f"  No {system_type} systems validated")

        except Exception as e:
            logger.exception(f"{system_type.title()} system validation failed")
            console.print(
                f"  [red]{system_type.title()} system validation failed; "
                f"aborting build: {e}[/red]"
            )
            raise
        finally:
            self._reacquire_db()

    # --- indexes/stats -------------------------------------------------- #
    def _create_indexes(self) -> None:
        idx_sql = [
            "CREATE INDEX IF NOT EXISTS idx_proteins_spatial ON proteins (contig_id, start, end_coord)",
            "CREATE INDEX IF NOT EXISTS idx_annotations_src_acc ON annotations (source, accession)",
            "CREATE INDEX IF NOT EXISTS idx_loci_type ON loci (locus_type)",
            "CREATE OR REPLACE VIEW bgc_loci AS SELECT * FROM loci WHERE locus_type = 'bgc'",
        ]
        for sql in idx_sql:
            try:
                self.conn.execute(sql)
            except Exception as exc:
                logger.warning(f"Index creation warning: {exc}")

    def _update_stats(self) -> None:
        for table in [
            "bins",
            "contigs",
            "proteins",
            "annotations",
            "loci",
            "protein_predicates",
            "semantic_atoms",
            "semantic_state",
        ]:
            try:
                count = self.conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
                # Map protein_predicates -> predicates for stats key.
                key = "predicates" if table == "protein_predicates" else table
                self.stats[key] = count
            except Exception:
                key = "predicates" if table == "protein_predicates" else table
                self.stats[key] = 0

    def _load_embeddings_info(self) -> None:
        """Record embeddings count/path if stage06 outputs are present."""
        manifest = self.outputs.stage06_dir / "embedding_manifest.json"
        if manifest.exists():
            try:
                data = json.loads(manifest.read_text())
                self.stats["embeddings"] = int(data.get("total_proteins", 0))
                h5_path = data.get("output_files", {}).get("embeddings")
                if h5_path:
                    self.embeddings_path = h5_path
            except Exception as exc:
                logger.warning(f"Failed to read embeddings manifest: {exc}")
        else:
            self.stats["embeddings"] = 0

    # --- helpers -------------------------------------------------------- #
    @staticmethod
    def _load_manifest(path: Path) -> Dict:
        if not path.exists():
            return {}
        try:
            return json.loads(path.read_text())
        except Exception:
            return {}


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def main(
    data_dir: Path = typer.Option(Path("data"), "--data-dir", "-d"),
    output: Path = typer.Option(Path("data/sharur.duckdb"), "--output", "-o"),
    force: bool = typer.Option(False, "--force"),
    resume_v2: bool = typer.Option(
        False,
        "--resume-v2",
        help="Continue V2 generation from the latest committed database chunk",
    ),
    enable_cazymes: bool = typer.Option(False, "--enable-cazymes", help="Run dbCAN CAZyme classification (slow, off by default)"),
    threads: Optional[int] = typer.Option(
        None,
        "--threads",
        "-t",
        help="Worker threads (defaults to SLURM_CPUS_ON_NODE, then host CPU count)",
    ),
) -> None:
    logging.basicConfig(level=logging.INFO)
    if resume_v2 and force:
        raise typer.BadParameter("--resume-v2 and --force select different build modes")
    canonical_embeddings_dir = data_dir / "embeddings"
    legacy_embeddings_dir = data_dir / "stage06_embeddings"
    embeddings_dir = (
        canonical_embeddings_dir
        if canonical_embeddings_dir.exists() or not legacy_embeddings_dir.exists()
        else legacy_embeddings_dir
    )
    outputs = PipelineOutputs(
        stage00_dir=data_dir / "stage00_prepared",
        stage01_dir=data_dir / "stage01_quast",
        stage02_dir=data_dir / "stage02_dfast_qc",
        stage03_dir=data_dir / "stage03_prodigal",
        stage04_dir=data_dir / "stage04_astra",
        stage05a_dir=data_dir / "stage05a_gecco",
        stage05b_dir=data_dir / "stage05b_dbcan",
        stage05c_dir=data_dir / "stage05c_crispr",
        stage06_dir=embeddings_dir,
    )
    console.print(f"Detected stage outputs: {outputs.validate()}")
    builder = KnowledgeBaseBuilder(
        outputs,
        output,
        force=force,
        enable_cazymes=enable_cazymes,
        threads=threads,
    )
    stats = builder.resume_v2() if resume_v2 else builder.build()
    console.print(f"[green]Build complete[/green]: {stats}")


if __name__ == "__main__":
    typer.run(main)
