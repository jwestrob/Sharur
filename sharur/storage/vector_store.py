"""Vector store interfaces for protein embedding similarity search."""

from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Optional

import numpy as np


class VectorStore(ABC):
    """Abstract interface for vector similarity search."""

    @abstractmethod
    def add(self, id: str, embedding: np.ndarray, metadata: dict | None = None) -> None:
        """Add embedding to store."""

    @abstractmethod
    def query(
        self,
        query_id: str,
        k: int = 10,
        threshold: Optional[float] = None,
        include_distances: bool = False,
    ) -> list[str] | list[tuple[str, float]]:
        """Find k nearest neighbors."""

    @abstractmethod
    def query_vector(
        self,
        vector: np.ndarray,
        k: int = 10,
        threshold: Optional[float] = None,
    ) -> list[tuple[str, float]]:
        """Find nearest neighbors to raw vector."""


def _cosine_similarity(a: np.ndarray, b: np.ndarray) -> float:
    """Compute cosine similarity with small epsilon for stability."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    denom = (np.linalg.norm(a) * np.linalg.norm(b)) + 1e-12
    return float(np.dot(a, b) / denom)


class FAISSStore(VectorStore):
    """FAISS-backed vector store. Loads from H5 embeddings on init.

    Builds an inner-product index on L2-normalized vectors so that
    FAISS scores equal cosine similarity directly.
    """

    def __init__(
        self,
        h5_path: str | Path,
        id_column: str = "protein_id",
        nprobe: int = 32,
    ):
        self._id_to_idx: dict[str, int] = {}
        self._idx_to_id: list[str] = []
        self._index = None
        self._vectors: np.ndarray | None = None

        h5_path = Path(h5_path)
        if not h5_path.exists():
            return

        try:
            import faiss
            import h5py
        except ImportError:
            return

        with h5py.File(h5_path, "r") as f:
            ids_raw = f["protein_ids"][:]
            ids = [x.decode() if isinstance(x, bytes) else x for x in ids_raw]
            vectors = np.asarray(f["embeddings"][:], dtype=np.float32)

        # L2-normalize so inner product == cosine similarity
        norms = np.linalg.norm(vectors, axis=1, keepdims=True)
        norms[norms == 0] = 1.0
        vectors = vectors / norms

        self._idx_to_id = ids
        self._id_to_idx = {pid: i for i, pid in enumerate(ids)}
        self._vectors = vectors

        dim = vectors.shape[1]
        n = vectors.shape[0]

        # Use IVF index for large datasets, flat for small
        if n > 100_000:
            nlist = min(int(np.sqrt(n)), 4096)
            quantizer = faiss.IndexFlatIP(dim)
            self._index = faiss.IndexIVFFlat(quantizer, dim, nlist, faiss.METRIC_INNER_PRODUCT)
            self._index.train(vectors)
            self._index.add(vectors)
            self._index.nprobe = nprobe
        else:
            self._index = faiss.IndexFlatIP(dim)
            self._index.add(vectors)

    def add(self, id: str, embedding: np.ndarray, metadata: dict | None = None) -> None:
        # FAISS index is built from H5 at init; runtime adds go to memory only
        vec = np.asarray(embedding, dtype=np.float32).reshape(1, -1)
        norm = np.linalg.norm(vec)
        if norm > 0:
            vec = vec / norm
        idx = len(self._idx_to_id)
        self._idx_to_id.append(id)
        self._id_to_idx[id] = idx
        if self._index is not None:
            self._index.add(vec)

    def query(
        self,
        query_id: str,
        k: int = 10,
        threshold: Optional[float] = None,
        include_distances: bool = False,
    ) -> list[str] | list[tuple[str, float]]:
        if self._index is None or query_id not in self._id_to_idx:
            return []

        idx = self._id_to_idx[query_id]
        vec = self._vectors[idx].reshape(1, -1) if self._vectors is not None else None
        if vec is None:
            return []

        # k+1 to account for self-match
        scores, indices = self._index.search(vec, k + 1)

        pairs = []
        for score, i in zip(scores[0], indices[0]):
            if i < 0 or i >= len(self._idx_to_id):
                continue
            rid = self._idx_to_id[i]
            if rid == query_id:
                continue
            sim = float(score)  # already cosine similarity (inner product of normalized vecs)
            if threshold is not None and sim < threshold:
                continue
            pairs.append((rid, sim))

        pairs = pairs[:k]

        if include_distances:
            return pairs
        return [rid for rid, _ in pairs]

    def query_vector(
        self,
        vector: np.ndarray,
        k: int = 10,
        threshold: Optional[float] = None,
    ) -> list[tuple[str, float]]:
        if self._index is None:
            return []

        vec = np.asarray(vector, dtype=np.float32).reshape(1, -1)
        norm = np.linalg.norm(vec)
        if norm > 0:
            vec = vec / norm

        scores, indices = self._index.search(vec, k)

        pairs = []
        for score, i in zip(scores[0], indices[0]):
            if i < 0 or i >= len(self._idx_to_id):
                continue
            rid = self._idx_to_id[i]
            sim = float(score)
            if threshold is not None and sim < threshold:
                continue
            pairs.append((rid, sim))

        return pairs


__all__ = ["VectorStore", "FAISSStore"]
