"""Memory-bounded protein FASTA to HDF5 embedding stream helpers."""

from __future__ import annotations

import logging
from collections import Counter
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import TYPE_CHECKING

import h5py
import numpy as np


if TYPE_CHECKING:
    from collections.abc import Callable, Iterable, Iterator
    from pathlib import Path

logger = logging.getLogger(__name__)


@dataclass
class SequenceStatistics:
    total_proteins: int = 0
    min_protein_length: int | None = None
    max_protein_length: int = 0
    total_residues: int = 0
    truncated_sequences: int = 0
    genome_protein_counts: Counter[str] = field(default_factory=Counter)

    def observe(
        self,
        sequence: str,
        genome_id: str,
        *,
        truncation_residue_limit: int | None = None,
    ) -> None:
        length = len(sequence)
        self.total_proteins += 1
        self.total_residues += length
        self.min_protein_length = (
            length if self.min_protein_length is None else min(self.min_protein_length, length)
        )
        self.max_protein_length = max(self.max_protein_length, length)
        if truncation_residue_limit is not None and length > truncation_residue_limit:
            self.truncated_sequences += 1
        self.genome_protein_counts[genome_id] += 1

    def to_dict(self) -> dict:
        mean = self.total_residues / self.total_proteins if self.total_proteins else 0.0
        return {
            "total_proteins": self.total_proteins,
            "genomes_processed": len(self.genome_protein_counts),
            "min_protein_length": self.min_protein_length,
            "max_protein_length": self.max_protein_length,
            "mean_protein_length": float(mean),
            "truncated_sequences": self.truncated_sequences,
            "genome_protein_counts": dict(self.genome_protein_counts),
        }


def _iter_fasta(path: Path, genome_id: str) -> Iterator[tuple[str, str, str]]:
    with open(path) as handle:
        protein_id: str | None = None
        parts: list[str] = []
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if protein_id is not None:
                    yield protein_id, "".join(parts).rstrip("*"), genome_id
                protein_id = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
        if protein_id is not None:
            yield protein_id, "".join(parts).rstrip("*"), genome_id


def iter_stage03_proteins(stage03_dir: Path) -> Iterator[tuple[str, str, str]]:
    """Yield `(protein_id, sequence, genome_id)` without materializing the proteome."""
    genomes_dir = stage03_dir / "genomes"
    if not genomes_dir.is_dir():
        raise FileNotFoundError(f"Genomes directory not found: {genomes_dir}")

    for genome_dir in sorted(genomes_dir.iterdir()):
        if not genome_dir.is_dir() or genome_dir.name == "all_protein_symlinks":
            continue
        faa_files = sorted(genome_dir.glob("*.faa"))
        if not faa_files:
            logger.warning("No .faa file for genome: %s", genome_dir.name)
            continue
        yield from _iter_fasta(faa_files[0], genome_dir.name)


def _validate_batch(
    values: np.ndarray,
    *,
    expected_rows: int,
    embedding_dim: int,
) -> np.ndarray:
    embeddings = np.asarray(values, dtype=np.float32)
    expected_shape = (expected_rows, embedding_dim)
    if embeddings.shape != expected_shape:
        raise ValueError(f"Embedding batch has shape {embeddings.shape}; expected {expected_shape}")
    if not np.isfinite(embeddings).all():
        raise ValueError("Embedding batch contains NaN or infinite values")
    if np.any(np.linalg.norm(embeddings, axis=1) == 0):
        raise ValueError("Embedding batch contains a zero vector")
    return embeddings


def stream_embeddings_to_h5(
    records: Iterable[tuple[str, str, str]],
    h5_path: Path,
    *,
    embedding_dim: int,
    model_name: str,
    batch_size: int,
    embed_batch: Callable[[list[str]], np.ndarray],
    truncation_residue_limit: int | None = None,
    on_progress: Callable[[int], None] | None = None,
) -> SequenceStatistics:
    """Stream records through an embedding callback into extendable HDF5 datasets."""
    if batch_size < 1:
        raise ValueError("batch_size must be positive")
    if embedding_dim < 1:
        raise ValueError("embedding_dim must be positive")

    h5_path.parent.mkdir(parents=True, exist_ok=True)
    stats = SequenceStatistics()
    batch_ids: list[str] = []
    batch_sequences: list[str] = []
    offset = 0

    with h5py.File(h5_path, "w") as handle:
        string_dtype = h5py.string_dtype(encoding="utf-8")
        ids_ds = handle.create_dataset(
            "protein_ids",
            shape=(0,),
            maxshape=(None,),
            chunks=(max(1, min(batch_size * 64, 65_536)),),
            dtype=string_dtype,
        )
        embeddings_ds = handle.create_dataset(
            "embeddings",
            shape=(0, embedding_dim),
            maxshape=(None, embedding_dim),
            chunks=(max(1, min(batch_size, 1_024)), embedding_dim),
            dtype="float32",
        )
        handle.attrs["model_name"] = model_name
        handle.attrs["embedding_dim"] = embedding_dim
        handle.attrs["created_at"] = datetime.now(timezone.utc).isoformat()

        def flush() -> None:
            nonlocal offset
            if not batch_ids:
                return
            embedded = _validate_batch(
                embed_batch(batch_sequences),
                expected_rows=len(batch_ids),
                embedding_dim=embedding_dim,
            )
            next_offset = offset + len(batch_ids)
            ids_ds.resize((next_offset,))
            embeddings_ds.resize((next_offset, embedding_dim))
            ids_ds[offset:next_offset] = batch_ids
            embeddings_ds[offset:next_offset] = embedded
            if on_progress is not None:
                on_progress(len(batch_ids))
            offset = next_offset
            batch_ids.clear()
            batch_sequences.clear()

        for protein_id, sequence, genome_id in records:
            if not protein_id:
                raise ValueError("Protein ID must be non-empty")
            if not sequence:
                raise ValueError(f"Protein {protein_id} has an empty sequence")
            batch_ids.append(protein_id)
            batch_sequences.append(sequence)
            stats.observe(
                sequence,
                genome_id,
                truncation_residue_limit=truncation_residue_limit,
            )
            if len(batch_ids) == batch_size:
                flush()
        flush()

        if stats.total_proteins == 0:
            raise ValueError("No protein sequences found in stage03 output")
        handle.attrs["num_proteins"] = stats.total_proteins

    return stats


__all__ = [
    "SequenceStatistics",
    "iter_stage03_proteins",
    "stream_embeddings_to_h5",
]
