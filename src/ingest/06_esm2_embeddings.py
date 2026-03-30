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

import logging
import json
import h5py
import numpy as np
from pathlib import Path
from typing import List, Optional, Dict, Any
from datetime import datetime

import torch
from transformers import AutoTokenizer, AutoModel
from rich.console import Console
from rich.progress import Progress

console = Console()
logger = logging.getLogger(__name__)


class ESM2EmbeddingGenerator:
    """Generate protein embeddings using ESM2, streaming to HDF5."""

    def __init__(self, model_name: str = "facebook/esm2_t6_8M_UR50D",
                 device: Optional[str] = None):
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
        console.print(f"Embedding dimension: {self.embedding_dim}")

    def embed_to_h5(self, protein_ids: List[str], sequences: List[str],
                    h5_path: Path, batch_size: int = 16) -> int:
        """
        Generate embeddings and stream directly to an HDF5 file.

        Preallocates the output datasets and writes each batch in-place,
        so peak memory is O(batch_size * D) not O(N * D).

        Returns the number of proteins embedded.
        """
        n = len(protein_ids)
        h5_path.parent.mkdir(parents=True, exist_ok=True)

        with h5py.File(h5_path, "w") as f:
            # Write protein IDs upfront
            f.create_dataset("protein_ids",
                             data=np.array(protein_ids, dtype="S"))

            # Preallocate embedding matrix — never fully in memory
            emb_ds = f.create_dataset("embeddings",
                                      shape=(n, self.embedding_dim),
                                      dtype="float32")

            f.attrs["model_name"] = self.model_name
            f.attrs["embedding_dim"] = self.embedding_dim
            f.attrs["num_proteins"] = n
            f.attrs["created_at"] = datetime.now().isoformat()

            # Stream batches directly to disk
            offset = 0
            with Progress(console=console) as progress:
                task = progress.add_task("Embedding proteins...", total=n)

                with torch.no_grad():
                    for i in range(0, n, batch_size):
                        batch_seqs = sequences[i:i + batch_size]
                        batch_embs = self._embed_batch(batch_seqs)

                        bs = batch_embs.shape[0]
                        emb_ds[offset:offset + bs] = batch_embs
                        offset += bs

                        progress.advance(task, bs)

        return offset

    def _embed_batch(self, sequences: List[str]) -> np.ndarray:
        """Embed a batch of sequences, returning (batch_size, D) float32 array."""
        try:
            inputs = self.tokenizer(
                sequences,
                return_tensors="pt",
                padding=True,
                truncation=True,
                max_length=1024,
            ).to(self.device)

            outputs = self.model(**inputs)
            embeddings = outputs.last_hidden_state.mean(dim=1)
            result = embeddings.cpu().numpy().astype(np.float32)

            del inputs, outputs, embeddings
            if self.device == "cuda":
                torch.cuda.empty_cache()
            elif self.device == "mps":
                torch.mps.empty_cache()

            return result

        except Exception as e:
            logger.error(f"Error embedding batch: {e}")
            if self.device == "cuda":
                torch.cuda.empty_cache()
            elif self.device == "mps":
                torch.mps.empty_cache()
            return np.zeros((len(sequences), self.embedding_dim), dtype=np.float32)


def load_protein_sequences(stage03_dir: Path) -> tuple[list[str], list[str], dict[str, str]]:
    """
    Load protein IDs and sequences from stage03 prodigal output.

    Returns:
        (protein_ids, sequences, genome_map)
        where genome_map maps protein_id -> genome_id
    """
    genomes_dir = stage03_dir / "genomes"
    if not genomes_dir.exists():
        raise FileNotFoundError(f"Genomes directory not found: {genomes_dir}")

    console.print(f"Loading protein sequences from: {genomes_dir}")

    protein_ids: list[str] = []
    sequences: list[str] = []
    genome_map: dict[str, str] = {}
    genome_count = 0

    for genome_dir in sorted(genomes_dir.iterdir()):
        if not genome_dir.is_dir() or genome_dir.name == "all_protein_symlinks":
            continue

        genome_id = genome_dir.name
        faa_files = list(genome_dir.glob("*.faa"))
        if not faa_files:
            logger.warning(f"No .faa file for genome: {genome_id}")
            continue

        genome_count += 1
        _read_fasta(faa_files[0], genome_id, protein_ids, sequences, genome_map)

    console.print(f"Loaded {len(protein_ids)} protein sequences from {genome_count} genomes")
    return protein_ids, sequences, genome_map


def _read_fasta(path: Path, genome_id: str,
                protein_ids: list[str], sequences: list[str],
                genome_map: dict[str, str]):
    """Parse a FASTA file, appending to the provided lists in-place."""
    with open(path) as f:
        current_id = None
        parts: list[str] = []

        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id and parts:
                    seq = "".join(parts).rstrip("*")
                    protein_ids.append(current_id)
                    sequences.append(seq)
                    genome_map[current_id] = genome_id
                current_id = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)

        if current_id and parts:
            seq = "".join(parts).rstrip("*")
            protein_ids.append(current_id)
            sequences.append(seq)
            genome_map[current_id] = genome_id


def write_manifest(h5_path: Path, protein_ids: List[str], sequences: List[str],
                   genome_map: dict[str, str], model_name: str) -> Path:
    """Write embedding_manifest.json alongside the H5 file."""
    from collections import Counter

    lengths = [len(s) for s in sequences]
    genome_counts = Counter(genome_map.values())

    manifest = {
        "version": "0.2.0",
        "created_at": datetime.now().isoformat(),
        "model_name": model_name,
        "embedding_dim": None,  # filled below
        "total_proteins": len(protein_ids),
        "genomes_processed": len(genome_counts),
        "output_files": {
            "embeddings": str(h5_path),
        },
        "statistics": {
            "min_protein_length": min(lengths),
            "max_protein_length": max(lengths),
            "mean_protein_length": float(np.mean(lengths)),
            "genome_protein_counts": dict(genome_counts),
        },
    }

    # Read back embedding dim from the H5 we just wrote
    with h5py.File(h5_path, "r") as f:
        manifest["embedding_dim"] = int(f.attrs["embedding_dim"])

    manifest_path = h5_path.parent / "embedding_manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)

    console.print(f"Manifest saved: {manifest_path}")
    return manifest_path


def run_esm2_embeddings(stage03_dir: Path, output_dir: Path,
                        model_name: str = "facebook/esm2_t6_8M_UR50D",
                        batch_size: Optional[int] = None,
                        force: bool = False) -> Dict[str, Any]:
    """Generate ESM2 embeddings for all proteins in stage03 output."""
    console.print("[bold blue]Stage 6: ESM2 Protein Embeddings[/bold blue]")

    h5_path = output_dir / "protein_embeddings.h5"
    manifest_path = output_dir / "embedding_manifest.json"

    if manifest_path.exists() and not force:
        console.print(f"[yellow]Embeddings already exist: {manifest_path}[/yellow]")
        console.print("[yellow]Use --force to overwrite[/yellow]")
        with open(manifest_path) as f:
            return json.load(f)

    # Load sequences (just IDs + strings, no heavy objects)
    protein_ids, sequences, genome_map = load_protein_sequences(stage03_dir)
    if not protein_ids:
        raise ValueError("No protein sequences found in stage03 output")

    # Initialize model
    generator = ESM2EmbeddingGenerator(model_name=model_name)

    if batch_size is None:
        batch_size = {"mps": 16, "cuda": 32}.get(generator.device, 4)

    console.print(f"Using batch size: {batch_size}")
    console.print(f"Streaming {len(protein_ids)} embeddings to: {h5_path}")

    # Stream embeddings directly to H5 — constant memory
    n_embedded = generator.embed_to_h5(protein_ids, sequences, h5_path, batch_size)

    # Write manifest
    write_manifest(h5_path, protein_ids, sequences, genome_map, model_name)

    console.print(f"[green]Done. {n_embedded} proteins embedded to {h5_path}[/green]")

    return {
        "embeddings_file": str(h5_path),
        "total_proteins": n_embedded,
        "embedding_dim": generator.embedding_dim,
    }


if __name__ == "__main__":
    import sys

    if len(sys.argv) not in (3, 4):
        print("Usage: python 06_esm2_embeddings.py <stage03_dir> <output_dir> [--force]")
        sys.exit(1)

    stage03_dir = Path(sys.argv[1])
    output_dir = Path(sys.argv[2])
    force = "--force" in sys.argv

    logging.basicConfig(level=logging.INFO)

    result = run_esm2_embeddings(stage03_dir, output_dir, force=force)
    print(f"ESM2 embeddings completed: {result['total_proteins']} proteins")
