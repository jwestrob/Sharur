"""Vector store interfaces (Part 8.3)."""

from __future__ import annotations

from abc import ABC, abstractmethod
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


class LanceDBStore(VectorStore):
    """LanceDB implementation (with in-memory fallback if lancedb unavailable)."""

    def __init__(self, uri: str, table_name: str = "embeddings", id_column: str = "id"):
        self.table_name = table_name
        self.id_column = id_column
        self._mem: dict[str, np.ndarray] = {}

        try:
            import lancedb
            self._db = lancedb.connect(uri)
            if table_name in self._db.table_names():
                self._table = self._db.open_table(table_name)
            else:
                # Avoid creating schema when database is empty; fall back to in-memory only.
                self._table = None
        except ModuleNotFoundError:
            # Graceful degradation: keep everything in memory.
            self._db = None
            self._table = None

    def add(self, id: str, embedding: np.ndarray, metadata: dict | None = None) -> None:
        metadata = metadata or {}
        vec = np.asarray(embedding, dtype=float)
        if self._table is not None:
            self._table.add([{"id": id, "vector": vec.tolist(), "metadata": metadata}])
        # Always keep in-memory cache for quick fallback queries
        self._mem[id] = vec

    def query(
        self,
        query_id: str,
        k: int = 10,
        threshold: Optional[float] = None,
        include_distances: bool = False,
    ) -> list[str] | list[tuple[str, float]]:
        # Try to use in-memory cache first (fast, avoids extra round trips)
        if query_id not in self._mem:
            if self._table is None:
                return []
            # Pull vector from LanceDB if not cached
            try:
                row = (
                    self._table.search(None)
                    .where(f"{self.id_column} = '{query_id}'")
                    .limit(1)
                    .to_list()
                )
            except Exception:
                row = []
            if not row:
                return []
            self._mem[query_id] = np.asarray(row[0]["vector"], dtype=float)

        anchor = self._mem[query_id]

        # If table is available, delegate search; otherwise use in-memory cosine.
        if self._table is not None:
            try:
                results = self._table.search(anchor).limit(k * 2).to_list()
                pairs = []
                for r in results:
                    rid = r.get(self.id_column) or r.get("id")
                    if rid == query_id:
                        continue  # skip self
                    dist = r.get("_distance")
                    if dist is not None:
                        # LanceDB returns L2 distance by default; convert to similarity
                        # Using 1/(1+d) gives range (0, 1] where 1 is identical
                        score = 1.0 / (1.0 + float(dist))
                    else:
                        # Estimate similarity manually if only vector present
                        score = _cosine_similarity(anchor, np.asarray(r["vector"], dtype=float))
                    pairs.append((rid, float(score)))
            except Exception:
                pairs = []
        else:
            pairs = []

        # Merge with in-memory cache for robustness
        for rid, vec in self._mem.items():
            if rid == query_id:
                continue
            pairs.append((rid, _cosine_similarity(anchor, vec)))

        # Deduplicate keeping max score
        merged: dict[str, float] = {}
        for rid, score in pairs:
            merged[rid] = max(score, merged.get(rid, -1.0))

        # Apply threshold and sort
        filtered = [
            (rid, score)
            for rid, score in merged.items()
            if threshold is None or score >= threshold
        ]
        filtered.sort(key=lambda x: x[1], reverse=True)
        filtered = filtered[:k]

        if include_distances:
            return filtered
        return [rid for rid, _ in filtered]

    def query_vector(
        self, vector: np.ndarray, k: int = 10, threshold: Optional[float] = None
    ) -> list[tuple[str, float]]:
        vec = np.asarray(vector, dtype=float)
        pairs: list[tuple[str, float]] = []

        if self._table is not None:
            try:
                results = self._table.search(vec).limit(k * 2).to_list()
                for r in results:
                    rid = r.get(self.id_column) or r.get("id")
                    dist = r.get("_distance")
                    if dist is not None:
                        # L2 distance to similarity: 1/(1+d)
                        score = 1.0 / (1.0 + float(dist))
                    else:
                        score = _cosine_similarity(vec, np.asarray(r["vector"], dtype=float))
                    pairs.append((rid, score))
            except Exception:
                pairs = []

        # Add in-memory cache
        for rid, cached in self._mem.items():
            pairs.append((rid, _cosine_similarity(vec, cached)))

        merged: dict[str, float] = {}
        for rid, score in pairs:
            merged[rid] = max(score, merged.get(rid, -1.0))

        filtered = [
            (rid, score)
            for rid, score in merged.items()
            if threshold is None or score >= threshold
        ]
        filtered.sort(key=lambda x: x[1], reverse=True)
        return filtered[:k]


class PgVectorStore(VectorStore):
    """pgvector implementation (best-effort, simple cosine similarity)."""

    def __init__(self, connection_string: str, table_name: str = "embeddings"):
        self.table_name = table_name
        try:
            import psycopg2

            self._pg = psycopg2
            self.conn = psycopg2.connect(connection_string)
            self.conn.autocommit = True
        except ModuleNotFoundError:
            self._pg = None
            self.conn = None

        # In-memory fallback to keep interface usable without Postgres
        self._mem: dict[str, np.ndarray] = {}

    def add(self, id: str, embedding: np.ndarray, metadata: dict | None = None) -> None:
        vec = np.asarray(embedding, dtype=float)
        self._mem[id] = vec

        if self.conn is None:
            return

        import json

        with self.conn.cursor() as cur:
            cur.execute(
                f"""
                CREATE TABLE IF NOT EXISTS {self.table_name} (
                    id TEXT PRIMARY KEY,
                    embedding VECTOR,
                    metadata JSONB
                );
                """,
            )
            cur.execute(
                f"""
                INSERT INTO {self.table_name} (id, embedding, metadata)
                VALUES (%s, %s, %s)
                ON CONFLICT (id) DO UPDATE
                  SET embedding = EXCLUDED.embedding,
                      metadata = EXCLUDED.metadata;
                """,
                (id, vec.tolist(), json.dumps(metadata or {})),
            )

    def query(
        self,
        query_id: str,
        k: int = 10,
        threshold: Optional[float] = None,
        include_distances: bool = False,
    ) -> list[str] | list[tuple[str, float]]:
        if query_id not in self._mem:
            return []
        anchor = self._mem[query_id]

        pairs = [
            (rid, _cosine_similarity(anchor, vec))
            for rid, vec in self._mem.items()
            if rid != query_id
        ]

        if threshold is not None:
            pairs = [(rid, s) for rid, s in pairs if s >= threshold]
        pairs.sort(key=lambda x: x[1], reverse=True)
        pairs = pairs[:k]

        if include_distances:
            return pairs
        return [rid for rid, _ in pairs]

    def query_vector(
        self, vector: np.ndarray, k: int = 10, threshold: Optional[float] = None
    ) -> list[tuple[str, float]]:
        vec = np.asarray(vector, dtype=float)
        pairs = [(rid, _cosine_similarity(vec, emb)) for rid, emb in self._mem.items()]
        if threshold is not None:
            pairs = [(rid, s) for rid, s in pairs if s >= threshold]
        pairs.sort(key=lambda x: x[1], reverse=True)
        return pairs[:k]


__all__ = ["VectorStore", "LanceDBStore", "PgVectorStore"]
