#!/usr/bin/env python3
"""
Stage 6: ESM2 Protein Embeddings

Generate protein embeddings using ESM2 and stream directly to HDF5.
No in-memory accumulation — embeddings are written to disk in chunks
as they are generated, so memory usage stays constant regardless of
dataset size.

Output format (ELSA-compatible):
    protein_embeddings.h5
    ├── protein_ids   # byte string array, shape (N,)
    └── embeddings    # float32 matrix, shape (N, D)
"""

import json
import logging
import os
import subprocess
import sys
from datetime import datetime, timezone
from itertools import chain
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import torch
from rich.console import Console
from rich.progress import Progress
from transformers import AutoModel, AutoTokenizer

from sharur.ingest.embedding_stream import (
    SequenceStatistics,
    iter_stage03_proteins,
    stream_embeddings_to_h5,
)


console = Console()
logger = logging.getLogger(__name__)


def _write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        temporary.write_text(json.dumps(payload, indent=2))
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _build_index_isolated(
    h5_path: Path,
    *,
    force: bool,
    threads: int | None,
) -> dict[str, Any]:
    """Build FAISS sidecars without importing FAISS into the Torch process."""
    command = [
        sys.executable,
        "-m",
        "sharur.ingest.vector_index_runner",
        "--embeddings",
        str(h5_path),
    ]
    if threads is not None:
        command.extend(["--threads", str(threads)])
    if force:
        command.append("--force")
    process = subprocess.run(
        command,
        check=True,
        capture_output=True,
        text=True,
    )
    try:
        return json.loads(process.stdout)
    except json.JSONDecodeError as exc:
        raise RuntimeError(
            "Vector-index subprocess returned invalid machine-readable output"
        ) from exc


class ESM2EmbeddingGenerator:
    """Generate protein embeddings using ESM2, streaming to HDF5."""

    def __init__(self, model_name: str = "facebook/esm2_t6_8M_UR50D", device: str | None = None):
        self.model_name = model_name

        if device is None:
            if torch.backends.mps.is_available() and torch.backends.mps.is_built():
                self.device = "mps"
            elif torch.cuda.is_available():
                self.device = "cuda"
            else:
                self.device = "cpu"
        else:
            self.device = device

        console.print(f"Loading ESM2 model: {model_name}")
        console.print(f"Using device: {self.device}")

        self.tokenizer = AutoTokenizer.from_pretrained(model_name)
        self.model = AutoModel.from_pretrained(model_name)
        self.model.to(self.device)
        self.model.eval()

        self.embedding_dim = self.model.config.hidden_size
        reported_max_length = getattr(self.tokenizer, "model_max_length", 1024)
        if not isinstance(reported_max_length, int) or not (1 <= reported_max_length < 1_000_000):
            reported_max_length = 1024
        self.max_token_length = min(1024, reported_max_length)
        self.max_residue_length = self.max_token_length - int(
            self.tokenizer.num_special_tokens_to_add(pair=False)
        )
        console.print(f"Embedding dimension: {self.embedding_dim}")
        console.print(
            f"Maximum residues per embedding: {self.max_residue_length} "
            "(longer proteins are truncated)"
        )

    def embed_records_to_h5(
        self,
        records,
        h5_path: Path,
        batch_size: int = 16,
    ) -> SequenceStatistics:
        """Stream FASTA records through the model without retaining the proteome."""
        with Progress(console=console) as progress:
            task = progress.add_task("Embedding proteins...", total=None)
            with torch.no_grad():
                return stream_embeddings_to_h5(
                    records,
                    h5_path,
                    embedding_dim=self.embedding_dim,
                    model_name=self.model_name,
                    batch_size=batch_size,
                    embed_batch=self._embed_batch,
                    truncation_residue_limit=self.max_residue_length,
                    on_progress=lambda count: progress.advance(task, count),
                )

    def _embed_batch(self, sequences: list[str]) -> np.ndarray:
        """Embed a batch of sequences, returning (batch_size, D) float32 array."""
        try:
            inputs = self.tokenizer(
                sequences,
                return_tensors="pt",
                padding=True,
                truncation=True,
                max_length=self.max_token_length,
                return_special_tokens_mask=True,
            ).to(self.device)

            special_tokens_mask = inputs.pop("special_tokens_mask")
            residue_mask = inputs["attention_mask"].bool() & ~special_tokens_mask.bool()
            outputs = self.model(**inputs)
            weights = residue_mask.unsqueeze(-1).to(outputs.last_hidden_state.dtype)
            embeddings = (outputs.last_hidden_state * weights).sum(dim=1)
            embeddings = embeddings / weights.sum(dim=1).clamp_min(1)
            result = embeddings.cpu().numpy().astype(np.float32)

            del inputs, outputs, embeddings, residue_mask, special_tokens_mask, weights
            if self.device == "cuda":
                torch.cuda.empty_cache()
            elif self.device == "mps":
                torch.mps.empty_cache()

            return result

        except Exception as e:
            logger.exception("Embedding batch failed")
            if self.device == "cuda":
                torch.cuda.empty_cache()
            elif self.device == "mps":
                torch.mps.empty_cache()
            raise RuntimeError("ESM2 inference failed; refusing to write zero embeddings") from e


def write_manifest(
    h5_path: Path,
    stats: SequenceStatistics,
    model_name: str,
) -> Path:
    """Write embedding_manifest.json alongside the H5 file."""
    summary = stats.to_dict()

    manifest = {
        "version": "0.2.0",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "model_name": model_name,
        "embedding_dim": None,  # filled below
        "total_proteins": summary["total_proteins"],
        "genomes_processed": summary["genomes_processed"],
        "output_files": {
            "embeddings": str(h5_path),
        },
        "statistics": {
            "min_protein_length": summary["min_protein_length"],
            "max_protein_length": summary["max_protein_length"],
            "mean_protein_length": summary["mean_protein_length"],
            "truncated_sequences": summary["truncated_sequences"],
            "genome_protein_counts": summary["genome_protein_counts"],
        },
    }

    # Read back embedding dim from the H5 we just wrote
    with h5py.File(h5_path, "r") as f:
        manifest["embedding_dim"] = int(f.attrs["embedding_dim"])

    manifest_path = h5_path.parent / "embedding_manifest.json"
    _write_json_atomic(manifest_path, manifest)

    console.print(f"Manifest saved: {manifest_path}")
    return manifest_path


def run_esm2_embeddings(
    stage03_dir: Path,
    output_dir: Path,
    model_name: str = "facebook/esm2_t6_8M_UR50D",
    batch_size: int | None = None,
    force: bool = False,
    build_index: bool = True,
    index_threads: int | None = None,
    device: str | None = None,
) -> dict[str, Any]:
    """Generate ESM2 embeddings and their persistent similarity sidecars."""
    console.print("[bold blue]Stage 6: ESM2 Protein Embeddings[/bold blue]")

    h5_path = output_dir / "protein_embeddings.h5"
    manifest_path = output_dir / "embedding_manifest.json"

    if manifest_path.exists() and h5_path.is_file() and not force:
        console.print(f"[yellow]Embeddings already exist: {manifest_path}[/yellow]")
        with open(manifest_path) as f:
            existing = json.load(f)
        if build_index and h5_path.is_file():
            inspection = _build_index_isolated(
                h5_path,
                force=False,
                threads=index_threads,
            )
            existing.setdefault("output_files", {}).update(
                {
                    "faiss_index": inspection["index_path"],
                    "protein_id_map": inspection["id_map_path"],
                    "index_manifest": inspection["manifest_path"],
                }
            )
            existing["similarity_index"] = inspection
            _write_json_atomic(manifest_path, existing)
        else:
            console.print("[yellow]Use --force to overwrite[/yellow]")
        return existing
    if manifest_path.exists() and not h5_path.is_file():
        console.print(
            "[yellow]Embedding manifest exists but its H5 payload is missing; "
            "regenerating instead of reusing it.[/yellow]"
        )

    records = iter_stage03_proteins(stage03_dir)
    first_record = next(records, None)
    if first_record is None:
        raise ValueError("No protein sequences found in stage03 output")
    records = chain((first_record,), records)

    # Initialize model
    generator = ESM2EmbeddingGenerator(model_name=model_name, device=device)

    if batch_size is None:
        batch_size = {"mps": 16, "cuda": 32}.get(generator.device, 4)

    console.print(f"Using batch size: {batch_size}")
    console.print(f"Streaming embeddings to: {h5_path}")

    # Stream embeddings directly to H5 — constant memory
    temporary_h5 = h5_path.with_name(f".{h5_path.name}.{os.getpid()}.tmp")
    temporary_h5.unlink(missing_ok=True)
    try:
        stats = generator.embed_records_to_h5(
            records,
            temporary_h5,
            batch_size,
        )
        os.replace(temporary_h5, h5_path)
    finally:
        temporary_h5.unlink(missing_ok=True)
    embedding_dim = generator.embedding_dim

    # Write manifest
    write_manifest(h5_path, stats, model_name)

    n_embedded = stats.total_proteins
    console.print(
        f"[green]Done. {n_embedded} proteins from "
        f"{len(stats.genome_protein_counts)} genomes embedded to {h5_path}[/green]"
    )

    index_result = None
    if build_index:
        # Model inference is complete before the CPU FAISS build begins. FAISS
        # runs in a Torch-free subprocess to isolate native OpenMP runtimes.
        del generator
        if torch.backends.mps.is_available():
            torch.mps.empty_cache()

        console.print("[blue]Building persistent FAISS similarity sidecars...[/blue]")
        index_result = _build_index_isolated(
            h5_path,
            force=force,
            threads=index_threads,
        )
        manifest = json.loads(manifest_path.read_text())
        manifest.setdefault("output_files", {}).update(
            {
                "faiss_index": index_result["index_path"],
                "protein_id_map": index_result["id_map_path"],
                "index_manifest": index_result["manifest_path"],
            }
        )
        manifest["similarity_index"] = index_result
        _write_json_atomic(manifest_path, manifest)
        console.print(f"[green]Similarity index ready: {index_result['index_path']}[/green]")

    result = {
        "embeddings_file": str(h5_path),
        "total_proteins": n_embedded,
        "embedding_dim": embedding_dim,
    }
    if index_result is not None:
        result["similarity_index"] = index_result
    return result


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate ESM2 H5 embeddings and persistent FAISS sidecars."
    )
    parser.add_argument("stage03_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--device",
        choices=("cpu", "mps", "cuda"),
        default=None,
        help="Execution device (default: auto-detect MPS, CUDA, then CPU).",
    )
    parser.add_argument(
        "--skip-index",
        action="store_true",
        help="Do not build persistent FAISS and protein-ID sidecars.",
    )
    parser.add_argument(
        "--index-threads",
        type=int,
        default=None,
        help="FAISS CPU threads for sidecar construction.",
    )
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)

    result = run_esm2_embeddings(
        args.stage03_dir,
        args.output_dir,
        force=args.force,
        build_index=not args.skip_index,
        index_threads=args.index_threads,
        device=args.device,
    )
    print(f"ESM2 embeddings completed: {result['total_proteins']} proteins")
