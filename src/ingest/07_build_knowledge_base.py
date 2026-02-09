#!/usr/bin/env python3
"""
Stage 07: Build Sharur knowledge base (DuckDB) from ingest outputs.

Lightweight implementation aligned to BENNU_PROJECT_SEED and INGEST_PIPELINE_SEED:
- Uses Sharur DuckDB schema from sharur.storage.schema.SCHEMA
- Ingests bins (stage02), proteins/contigs (stage03), annotations (stage04/dbCAN),
  loci (stage05a/05c when present), and basic feature_store metrics.
- Skips ELSA/vector store work per current deferral.
"""

from __future__ import annotations

import json
import logging
import multiprocessing
import os
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import duckdb
import pandas as pd
import typer
from rich.console import Console
from rich.progress import Progress

from sharur.storage.schema import SCHEMA
from sharur.predicates.generator import (
    PredicateGenerator,
    AnnotationRecord,
    ProteinRecord,
)

console = Console()
logger = logging.getLogger(__name__)


# --------------------------------------------------------------------------- #
# Top-level worker function for parallel predicate generation (must be
# picklable, so it lives at module scope rather than as a method).
# --------------------------------------------------------------------------- #

def _generate_predicates_chunk(item: tuple) -> tuple:
    """Generate predicates for a single protein.  Called by Pool.imap_unordered.

    Args:
        item: (pid, length, gc, gc_mean, gc_std, ann_tuples)
              where ann_tuples is a list of (source, accession, name, desc, evalue, score)

    Returns:
        (protein_id, predicate_list)
    """
    pid, length, gc, gc_mean, gc_std, ann_tuples = item

    # Each worker lazily initialises its own generator (cached on the function object)
    gen = getattr(_generate_predicates_chunk, "_gen", None)
    if gen is None:
        gen = PredicateGenerator(predict_topology=False)
        _generate_predicates_chunk._gen = gen

    protein = ProteinRecord(
        protein_id=pid,
        sequence_length=length,
        gc_content=gc,
        contig_gc_mean=gc_mean,
        contig_gc_std=gc_std,
    )
    annotations = [
        AnnotationRecord(source=t[0], accession=t[1], name=t[2],
                         description=t[3], evalue=t[4], score=t[5])
        for t in ann_tuples
    ]
    predicates = gen.generate_for_protein(protein, annotations)
    return (pid, predicates)


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
    def __init__(self, outputs: PipelineOutputs, db_path: Path, force: bool = False):
        self.outputs = outputs
        self.db_path = db_path
        self.force = force
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
    def build(self) -> Dict[str, int]:
        self._init_db()

        with Progress() as progress:
            task = progress.add_task("Building knowledge base...", total=9)

            progress.update(task, description="Loading bins")
            self._load_bins()
            progress.advance(task)

            progress.update(task, description="Loading proteins/contigs")
            self._load_proteins_and_contigs()
            progress.advance(task)

            progress.update(task, description="Loading annotations")
            self._load_annotations()
            progress.advance(task)

            progress.update(task, description="Loading loci")
            self._load_loci()
            progress.advance(task)

            progress.update(task, description="Loading embeddings info")
            self._load_embeddings_info()
            progress.advance(task)

            progress.update(task, description="Computing metrics")
            self._compute_basic_metrics()
            progress.advance(task)

            progress.update(task, description="Generating predicates")
            self._generate_predicates()
            progress.advance(task)

            progress.update(task, description="Classifying hydrogenases")
            self._classify_hydrogenases()
            progress.advance(task)

            progress.update(task, description="Finalizing")
            self._create_indexes()
            self._update_stats()
            progress.advance(task)

        return self.stats

    # --- init ----------------------------------------------------------- #
    def _init_db(self) -> None:
        if self.db_path.exists():
            if self.force:
                self.db_path.unlink()
            else:
                raise FileExistsError(f"{self.db_path} exists (use --force to overwrite)")

        self.conn = duckdb.connect(str(self.db_path))
        self.conn.execute(SCHEMA)
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

        for faa in genomes_dir.glob("**/*.faa"):
            if faa.parent.name == "all_protein_symlinks":
                continue  # skip symlink aggregation dir to avoid duplicates
            bin_id = faa.parent.name
            protein_rows.extend(self._parse_prodigal_faa(faa, bin_id, contig_lengths))

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
        for tsv in (self.outputs.stage04_dir.rglob("*_hits_df.tsv") if self.outputs.stage04_dir.exists() else []):
            try:
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
            except Exception as exc:
                logger.warning(f"Failed to load annotations from {tsv}: {exc}")

        # dbCAN JSON
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
                    self.conn.execute(
                        """
                        INSERT INTO loci (
                            locus_id, locus_type, contig_id, start, end_coord, confidence,
                            metadata
                        ) SELECT * FROM tmp_ldf_cr
                        """,
                    )
                    self.stats["loci"] += len(loci_rows)
            except Exception as exc:
                logger.warning(f"Failed to load CRISPR arrays from {crispr_json}: {exc}")

    # --- metrics -------------------------------------------------------- #
    def _compute_basic_metrics(self) -> None:
        # length_zscore within bins
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
    def _generate_predicates(self) -> None:
        """Generate semantic predicates for all proteins based on annotations.

        Uses multiprocessing to parallelize the pure-function predicate
        generation across CPU cores, then bulk-inserts results.
        """
        console.print("[cyan]Loading protein and annotation data for predicate generation...[/cyan]")

        # Get GC stats per contig
        gc_stats = {}
        gc_rows = self.conn.execute("""
            SELECT contig_id, AVG(gc_content), STDDEV(gc_content)
            FROM proteins
            WHERE gc_content IS NOT NULL
            GROUP BY contig_id
        """).fetchall()
        for row in gc_rows:
            gc_stats[row[0]] = (row[1], row[2] if row[2] else 0.0)

        # Get all proteins
        proteins = self.conn.execute("""
            SELECT protein_id, sequence_length, gc_content, contig_id
            FROM proteins
        """).fetchall()

        console.print(f"  {len(proteins):,} proteins loaded")

        # Get all annotations grouped by protein
        annotations_by_protein: Dict[str, List[tuple]] = {}
        ann_rows = self.conn.execute("""
            SELECT protein_id, source, accession, name, description, evalue, score
            FROM annotations
        """).fetchall()
        for row in ann_rows:
            pid = row[0]
            if pid not in annotations_by_protein:
                annotations_by_protein[pid] = []
            # Store as plain tuples for pickle-ability across processes
            annotations_by_protein[pid].append(
                (row[1], row[2], row[3], row[4], row[5], row[6])
            )

        console.print(f"  {len(ann_rows):,} annotations loaded")

        # Build work items: list of (pid, length, gc, contig_gc_mean, contig_gc_std, ann_tuples)
        work_items = []
        for row in proteins:
            pid, length, gc, contig = row[0], row[1], row[2], row[3]
            gc_mean, gc_std = gc_stats.get(contig, (None, None))
            anns = annotations_by_protein.get(pid, [])
            work_items.append((pid, length, gc, gc_mean, gc_std, anns))

        # Determine worker count
        n_workers = min(os.cpu_count() or 4, 12)
        chunk_size = max(1, len(work_items) // (n_workers * 4))

        console.print(f"[cyan]Generating predicates with {n_workers} workers (chunk_size={chunk_size:,})...[/cyan]")

        # Parallel predicate generation
        results = []
        with multiprocessing.Pool(n_workers) as pool:
            for batch_result in pool.imap_unordered(
                _generate_predicates_chunk, work_items, chunksize=chunk_size
            ):
                results.append(batch_result)
                if len(results) % 100_000 == 0:
                    console.print(f"  {len(results):,}/{len(work_items):,} proteins processed")

        console.print(f"  {len(results):,} proteins processed, inserting into DB...")

        # Bulk insert via DataFrame register
        insert_batch_size = 50_000
        for i in range(0, len(results), insert_batch_size):
            batch = results[i:i + insert_batch_size]
            df = pd.DataFrame(batch, columns=["protein_id", "predicates"])
            df["updated_at"] = datetime.now(timezone.utc)
            self.conn.register("tmp_pred_batch", df)
            self.conn.execute("""
                INSERT INTO protein_predicates (protein_id, predicates, updated_at)
                SELECT protein_id, predicates, updated_at FROM tmp_pred_batch
            """)
            self.conn.unregister("tmp_pred_batch")

        self.stats["predicates"] = len(proteins)
        console.print(f"[green]  Predicates generated for {len(proteins):,} proteins[/green]")

        # Flag proteins overlapping CRISPR arrays
        self._flag_crispr_array_overlaps()

    def _flag_crispr_array_overlaps(self) -> None:
        """Add 'in_crispr_array' predicate to proteins overlapping CRISPR arrays."""
        # Find proteins that overlap CRISPR arrays
        overlaps = self.conn.execute("""
            SELECT DISTINCT p.protein_id
            FROM proteins p
            JOIN loci l ON p.contig_id = l.contig_id
            WHERE l.locus_type = 'crispr_array'
              AND p.start < l.end_coord
              AND p.end_coord > l.start
        """).fetchall()

        if not overlaps:
            return

        # Add predicate to each overlapping protein
        for (protein_id,) in overlaps:
            self.conn.execute("""
                UPDATE protein_predicates
                SET predicates = list_append(predicates, 'in_crispr_array'),
                    updated_at = CURRENT_TIMESTAMP
                WHERE protein_id = ?
                  AND NOT list_contains(predicates, 'in_crispr_array')
            """, [protein_id])

        console.print(f"  Flagged {len(overlaps)} proteins overlapping CRISPR arrays")

    def _insert_predicate_batch(self, batch: List[tuple]) -> None:
        """Insert batch of predicates into protein_predicates table.
        Legacy row-by-row method kept for compatibility; _generate_predicates
        now uses bulk DataFrame insert instead.
        """
        for protein_id, predicates in batch:
            self.conn.execute(
                """
                INSERT INTO protein_predicates (protein_id, predicates, updated_at)
                VALUES (?, ?, CURRENT_TIMESTAMP)
                """,
                [protein_id, predicates],
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

        # Commit current changes so classification script can read them
        self.conn.commit()

        console.print(f"  Found {hyddb_count} HydDB annotations, running subgroup classification...")

        try:
            # Import the classification function
            import sys
            scripts_dir = Path(__file__).resolve().parents[2] / "scripts"
            sys.path.insert(0, str(scripts_dir))

            from classify_hydrogenases import classify_hydrogenases as run_classification

            # Run classification (updates predicates in place)
            results = run_classification(
                db_path=str(self.db_path),
                threads=4,
                update_predicates=True,
                verbose=False,
            )

            if not results.empty:
                validated = results['validated'].sum()
                console.print(f"  Classified {len(results)} hydrogenases, {validated} validated")
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
        for table in ["bins", "contigs", "proteins", "annotations", "loci", "protein_predicates"]:
            try:
                count = self.conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
                # Map protein_predicates -> predicates for stats key
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
                lancedb_path = data.get("output_files", {}).get("lancedb")
                if lancedb_path:
                    self.embeddings_path = lancedb_path
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
) -> None:
    logging.basicConfig(level=logging.INFO)
    outputs = PipelineOutputs(
        stage00_dir=data_dir / "stage00_prepared",
        stage01_dir=data_dir / "stage01_quast",
        stage02_dir=data_dir / "stage02_dfast_qc",
        stage03_dir=data_dir / "stage03_prodigal",
        stage04_dir=data_dir / "stage04_astra",
        stage05a_dir=data_dir / "stage05a_gecco",
        stage05b_dir=data_dir / "stage05b_dbcan",
        stage05c_dir=data_dir / "stage05c_crispr",
        stage06_dir=data_dir / "stage06_embeddings",
    )
    console.print(f"Detected stage outputs: {outputs.validate()}")
    builder = KnowledgeBaseBuilder(outputs, output, force=force)
    stats = builder.build()
    console.print(f"[green]Build complete[/green]: {stats}")


if __name__ == "__main__":
    typer.run(main)
