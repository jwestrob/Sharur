"""Persistent FAISS-backed protein embedding similarity search.

The HDF5 file remains the canonical embedding artifact.  Sharur builds a
validated, generation-scoped FAISS sidecar plus a SQLite row-to-protein map.
The sidecars make ordinary startup cheap and let large IVF indexes be opened
with FAISS' read-only mmap mode instead of rebuilding the index and two Python
ID maps in every process.
"""

from __future__ import annotations

import contextlib
import hashlib
import json
import math
import os
import sqlite3
import threading
import uuid
from abc import ABC, abstractmethod
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np


if TYPE_CHECKING:
    from collections.abc import Iterator

INDEX_MANIFEST_VERSION = 1
DEFAULT_CHUNK_SIZE = 50_000
FLAT_INDEX_LIMIT = 100_000
MAX_TRAINING_VECTORS = 200_000


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
        threshold: float | None = None,
        include_distances: bool = False,
    ) -> list[str] | list[tuple[str, float]]:
        """Find k nearest neighbors."""

    @abstractmethod
    def query_vector(
        self,
        vector: np.ndarray,
        k: int = 10,
        threshold: float | None = None,
    ) -> list[tuple[str, float]]:
        """Find nearest neighbors to a raw vector."""


@dataclass(frozen=True)
class VectorIndexInspection:
    """Non-loading status for a persistent vector index."""

    state: str  # available | unavailable | stale | failed
    summary: str
    h5_path: str
    manifest_path: str
    index_path: str | None = None
    id_map_path: str | None = None
    count: int | None = None
    dimension: int | None = None
    index_kind: str | None = None
    source_signature: str | None = None
    detail: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class _IndexPaths:
    h5: Path
    manifest: Path
    lock: Path

    @classmethod
    def from_h5(cls, h5_path: str | Path) -> _IndexPaths:
        h5 = Path(h5_path).expanduser().resolve()
        stem = h5.stem
        return cls(
            h5=h5,
            manifest=h5.with_name(f"{stem}.index.json"),
            lock=h5.with_name(f"{stem}.index.lock"),
        )


def _normalize(vectors: np.ndarray, *, copy: bool = True) -> np.ndarray:
    values = np.asarray(vectors, dtype=np.float32)
    if values.ndim != 2:
        raise ValueError(f"Embedding matrix must be 2-D, got shape {values.shape}")
    if copy:
        values = values.copy()
    if not np.isfinite(values).all():
        raise ValueError("Embedding matrix contains NaN or infinite values")
    norms = np.linalg.norm(values, axis=1, keepdims=True)
    zero_rows = np.flatnonzero(norms[:, 0] == 0)
    if zero_rows.size:
        raise ValueError(
            f"Embedding matrix contains zero vectors; first invalid row is {int(zero_rows[0])}"
        )
    values /= norms
    return values


def _decode_id(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def _source_description(h5_path: Path) -> tuple[dict[str, Any], str]:
    """Describe the canonical H5 artifact without reading its vector matrix."""
    import h5py  # noqa: PLC0415 - optional heavy dependency

    stat = h5_path.stat()
    with h5py.File(h5_path, "r") as handle:
        if "protein_ids" not in handle or "embeddings" not in handle:
            raise ValueError("H5 must contain protein_ids and embeddings datasets")
        ids = handle["protein_ids"]
        embeddings = handle["embeddings"]
        if len(ids.shape) != 1 or len(embeddings.shape) != 2:
            raise ValueError("protein_ids must be 1-D and embeddings must be 2-D")
        if ids.shape[0] != embeddings.shape[0]:
            raise ValueError(
                "protein_ids and embeddings row counts differ: "
                f"{ids.shape[0]} != {embeddings.shape[0]}"
            )
        if embeddings.shape[0] == 0 or embeddings.shape[1] == 0:
            raise ValueError("Embedding H5 contains an empty embedding matrix")
        if not np.issubdtype(embeddings.dtype, np.number):
            raise ValueError(f"Embedding dataset must be numeric, got {embeddings.dtype}")
        description = {
            "path": str(h5_path),
            "size": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
            "count": int(embeddings.shape[0]),
            "dimension": int(embeddings.shape[1]),
            "model_name": _decode_id(handle.attrs.get("model_name", "")),
            "created_at": _decode_id(handle.attrs.get("created_at", "")),
        }
    encoded = json.dumps(
        description,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return description, hashlib.sha256(encoded).hexdigest()


def _safe_generation_path(parent: Path, name: Any) -> Path | None:
    if not isinstance(name, str) or Path(name).name != name:
        return None
    path = parent / name
    return path if path.is_file() and path.stat().st_size > 0 else None


def inspect_vector_index(h5_path: str | Path) -> VectorIndexInspection:
    """Inspect sidecar readiness without opening the FAISS index."""
    paths = _IndexPaths.from_h5(h5_path)
    base = {
        "h5_path": str(paths.h5),
        "manifest_path": str(paths.manifest),
    }
    if not paths.h5.is_file():
        return VectorIndexInspection(
            state="unavailable",
            summary="Canonical embedding H5 is absent.",
            **base,
        )

    try:
        source, source_signature = _source_description(paths.h5)
    except Exception as exc:
        return VectorIndexInspection(
            state="failed",
            summary="Canonical embedding H5 could not be inspected.",
            detail=f"{type(exc).__name__}: {exc}",
            **base,
        )

    common = {
        **base,
        "count": source["count"],
        "dimension": source["dimension"],
        "source_signature": source_signature,
    }
    if not paths.manifest.is_file():
        return VectorIndexInspection(
            state="unavailable",
            summary="Persistent FAISS sidecars have not been built.",
            **common,
        )

    try:
        manifest = json.loads(paths.manifest.read_text())
    except Exception as exc:
        return VectorIndexInspection(
            state="failed",
            summary="Vector index manifest is unreadable.",
            detail=f"{type(exc).__name__}: {exc}",
            **common,
        )

    if manifest.get("manifest_version") != INDEX_MANIFEST_VERSION:
        return VectorIndexInspection(
            state="stale",
            summary="Vector index manifest version is obsolete.",
            detail=(
                f"found {manifest.get('manifest_version')!r}; expected {INDEX_MANIFEST_VERSION}"
            ),
            **common,
        )
    if manifest.get("source_signature") != source_signature:
        return VectorIndexInspection(
            state="stale",
            summary="Vector sidecars do not match the current embedding H5.",
            index_kind=manifest.get("index_kind"),
            detail="Rebuild the persistent vector index.",
            **common,
        )

    index_path = _safe_generation_path(paths.h5.parent, manifest.get("index_file"))
    id_map_path = _safe_generation_path(paths.h5.parent, manifest.get("id_map_file"))
    if index_path is None or id_map_path is None:
        return VectorIndexInspection(
            state="failed",
            summary="Vector index manifest references missing sidecar files.",
            index_kind=manifest.get("index_kind"),
            detail="Rebuild the persistent vector index.",
            **common,
        )

    try:
        uri = f"{id_map_path.resolve().as_uri()}?mode=ro&immutable=1"
        conn = sqlite3.connect(uri, uri=True)
        try:
            row_count = conn.execute("SELECT COUNT(*) FROM vector_ids").fetchone()[0]
            db_meta = dict(conn.execute("SELECT key, value FROM metadata").fetchall())
        finally:
            conn.close()
        if row_count != source["count"]:
            raise ValueError(f"ID map has {row_count} rows; embeddings have {source['count']}")
        if db_meta.get("source_signature") != source_signature:
            raise ValueError("ID map source signature does not match")
        if db_meta.get("build_id") != manifest.get("build_id"):
            raise ValueError("ID map build ID does not match manifest")
        if db_meta.get("id_sha256") != manifest.get("id_sha256"):
            raise ValueError("ID map digest does not match manifest")
    except Exception as exc:
        return VectorIndexInspection(
            state="failed",
            summary="Persistent protein-ID map failed validation.",
            index_path=str(index_path),
            id_map_path=str(id_map_path),
            index_kind=manifest.get("index_kind"),
            detail=f"{type(exc).__name__}: {exc}",
            **common,
        )

    return VectorIndexInspection(
        state="available",
        summary="Persistent FAISS index and stable protein-ID map are ready.",
        index_path=str(index_path),
        id_map_path=str(id_map_path),
        index_kind=manifest.get("index_kind"),
        **common,
    )


@contextlib.contextmanager
def _exclusive_build_lock(lock_path: Path) -> Iterator[None]:
    """Serialize sidecar construction across local processes."""
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with open(lock_path, "a+b") as handle:
        try:
            import fcntl  # noqa: PLC0415 - unavailable on some platforms
        except ImportError:  # pragma: no cover - Sharur targets macOS/Linux
            yield
            return
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def _build_id_map(
    h5_path: Path,
    sqlite_path: Path,
    *,
    source_signature: str,
    build_id: str,
    chunk_size: int,
) -> str:
    import h5py  # noqa: PLC0415 - optional heavy dependency

    digest = hashlib.sha256()
    conn = sqlite3.connect(str(sqlite_path))
    try:
        conn.executescript(
            """
            PRAGMA journal_mode=OFF;
            PRAGMA synchronous=OFF;
            PRAGMA temp_store=MEMORY;
            CREATE TABLE vector_ids (
                row_index INTEGER PRIMARY KEY,
                protein_id TEXT NOT NULL UNIQUE
            );
            CREATE TABLE metadata (
                key TEXT PRIMARY KEY,
                value TEXT NOT NULL
            );
            """
        )
        with h5py.File(h5_path, "r") as handle:
            ids_ds = handle["protein_ids"]
            count = int(ids_ds.shape[0])
            for start in range(0, count, chunk_size):
                decoded = [_decode_id(value) for value in ids_ds[start : start + chunk_size]]
                rows = []
                for offset, protein_id in enumerate(decoded):
                    encoded = protein_id.encode("utf-8")
                    digest.update(len(encoded).to_bytes(8, "big"))
                    digest.update(encoded)
                    rows.append((start + offset, protein_id))
                conn.executemany(
                    "INSERT INTO vector_ids(row_index, protein_id) VALUES (?, ?)",
                    rows,
                )
        conn.executemany(
            "INSERT INTO metadata(key, value) VALUES (?, ?)",
            [
                ("source_signature", source_signature),
                ("build_id", build_id),
                ("id_sha256", digest.hexdigest()),
            ],
        )
        conn.commit()
        conn.execute("PRAGMA optimize")
    finally:
        conn.close()
    return digest.hexdigest()


def _training_rows(count: int, nlist: int) -> np.ndarray:
    training_count = min(count, MAX_TRAINING_VECTORS)
    training_count = min(
        count,
        max(training_count if count <= MAX_TRAINING_VECTORS else 0, nlist * 40),
    )
    return np.linspace(0, count - 1, num=training_count, dtype=np.int64)


def _build_faiss_index(
    h5_path: Path,
    index_path: Path,
    *,
    chunk_size: int,
    nprobe: int,
    threads: int | None,
) -> tuple[str, int]:
    import faiss  # noqa: PLC0415 - optional heavy dependency
    import h5py  # noqa: PLC0415 - optional heavy dependency

    if threads is not None:
        faiss.omp_set_num_threads(max(1, threads))

    with h5py.File(h5_path, "r") as handle:
        vectors_ds = handle["embeddings"]
        count, dimension = map(int, vectors_ds.shape)
        if count <= FLAT_INDEX_LIMIT:
            index = faiss.IndexFlatIP(dimension)
            index_kind = "flat_ip"
            nlist = 0
        else:
            nlist = min(max(1, int(math.sqrt(count))), 4096)
            quantizer = faiss.IndexFlatIP(dimension)
            index = faiss.IndexIVFFlat(
                quantizer,
                dimension,
                nlist,
                faiss.METRIC_INNER_PRODUCT,
            )
            rows = _training_rows(count, nlist)
            training = _normalize(
                np.asarray(vectors_ds[rows], dtype=np.float32),
                copy=False,
            )
            index.train(training)
            del training
            index_kind = "ivf_flat_ip"

        for start in range(0, count, chunk_size):
            chunk = _normalize(
                np.asarray(vectors_ds[start : start + chunk_size], dtype=np.float32),
                copy=False,
            )
            index.add(chunk)

        if hasattr(index, "nprobe"):
            index.nprobe = min(max(1, nprobe), nlist)
        faiss.write_index(index, str(index_path))
    return index_kind, nlist


def build_vector_index(
    h5_path: str | Path,
    *,
    force: bool = False,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
    nprobe: int = 32,
    threads: int | None = None,
) -> VectorIndexInspection:
    """Build persistent sidecars and atomically publish their manifest."""
    if chunk_size < 1:
        raise ValueError("chunk_size must be positive")
    paths = _IndexPaths.from_h5(h5_path)
    if not paths.h5.is_file():
        raise FileNotFoundError(paths.h5)

    with _exclusive_build_lock(paths.lock):
        existing = inspect_vector_index(paths.h5)
        if existing.state == "available" and not force:
            return existing

        source, source_signature = _source_description(paths.h5)
        build_id = (
            f"{source_signature[:12]}-"
            f"{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S')}-"
            f"{uuid.uuid4().hex[:8]}"
        )
        stem = paths.h5.stem
        index_name = f"{stem}.{build_id}.faiss"
        id_map_name = f"{stem}.{build_id}.ids.sqlite"
        index_final = paths.h5.parent / index_name
        id_map_final = paths.h5.parent / id_map_name
        index_tmp = paths.h5.parent / f".{index_name}.tmp"
        id_map_tmp = paths.h5.parent / f".{id_map_name}.tmp"
        manifest_tmp = paths.h5.parent / f".{paths.manifest.name}.{build_id}.tmp"

        for temporary in (index_tmp, id_map_tmp, manifest_tmp):
            temporary.unlink(missing_ok=True)

        try:
            id_sha256 = _build_id_map(
                paths.h5,
                id_map_tmp,
                source_signature=source_signature,
                build_id=build_id,
                chunk_size=chunk_size,
            )
            index_kind, nlist = _build_faiss_index(
                paths.h5,
                index_tmp,
                chunk_size=chunk_size,
                nprobe=nprobe,
                threads=threads,
            )
            manifest = {
                "manifest_version": INDEX_MANIFEST_VERSION,
                "build_id": build_id,
                "created_at": datetime.now(timezone.utc).isoformat(),
                "source": source,
                "source_signature": source_signature,
                "id_sha256": id_sha256,
                "count": source["count"],
                "dimension": source["dimension"],
                "index_kind": index_kind,
                "nlist": nlist,
                "nprobe": nprobe,
                "index_file": index_name,
                "id_map_file": id_map_name,
            }
            manifest_tmp.write_text(json.dumps(manifest, indent=2, sort_keys=True))

            os.replace(index_tmp, index_final)
            os.replace(id_map_tmp, id_map_final)
            os.replace(manifest_tmp, paths.manifest)
        finally:
            for temporary in (index_tmp, id_map_tmp, manifest_tmp):
                temporary.unlink(missing_ok=True)

        inspection = inspect_vector_index(paths.h5)
        if inspection.state != "available":
            raise RuntimeError(
                f"Built vector sidecars did not validate: "
                f"{inspection.state}: {inspection.summary} {inspection.detail or ''}"
            )
        # The manifest is the atomic commit point. Older generation files can
        # now be unlinked; already-open mmap readers retain their file handle.
        keep = {
            Path(inspection.index_path).resolve(),
            Path(inspection.id_map_path).resolve(),
        }
        patterns = (
            f"{paths.h5.stem}.*.faiss",
            f"{paths.h5.stem}.*.ids.sqlite",
            f"{paths.h5.stem}.*.ids.sqlite-*",
        )
        for pattern in patterns:
            for candidate in paths.h5.parent.glob(pattern):
                if candidate.resolve() not in keep:
                    with contextlib.suppress(OSError):
                        candidate.unlink()
        return inspection


class FAISSStore(VectorStore):
    """Read-only persistent FAISS store backed by a canonical embedding H5."""

    def __init__(
        self,
        h5_path: str | Path,
        id_column: str = "protein_id",
        nprobe: int = 32,
        *,
        build_if_missing: bool = True,
    ):
        del id_column  # retained for API compatibility
        self.h5_path = Path(h5_path).expanduser().resolve()
        self.nprobe = nprobe
        self._index = None
        self._id_conn: sqlite3.Connection | None = None
        self._h5 = None
        self._lock = threading.RLock()
        self.inspection = inspect_vector_index(self.h5_path)

        if self.inspection.state != "available" and build_if_missing:
            self.inspection = build_vector_index(
                self.h5_path,
                force=self.inspection.state in {"stale", "failed"},
                nprobe=nprobe,
            )
        if self.inspection.state == "available":
            self._open()

    @classmethod
    def inspect(cls, h5_path: str | Path) -> VectorIndexInspection:
        return inspect_vector_index(h5_path)

    @classmethod
    def build_persistent_index(
        cls,
        h5_path: str | Path,
        *,
        force: bool = False,
        chunk_size: int = DEFAULT_CHUNK_SIZE,
        nprobe: int = 32,
        threads: int | None = None,
    ) -> VectorIndexInspection:
        return build_vector_index(
            h5_path,
            force=force,
            chunk_size=chunk_size,
            nprobe=nprobe,
            threads=threads,
        )

    @property
    def ready(self) -> bool:
        return self._index is not None and self._id_conn is not None

    def _open(self) -> None:
        import faiss  # noqa: PLC0415 - optional heavy dependency
        import h5py  # noqa: PLC0415 - optional heavy dependency

        assert self.inspection.index_path is not None
        assert self.inspection.id_map_path is not None
        flags = getattr(faiss, "IO_FLAG_MMAP", 0) | getattr(
            faiss,
            "IO_FLAG_READ_ONLY",
            0,
        )
        if not flags:
            raise RuntimeError("Installed FAISS does not expose read-only mmap flags")
        self._index = faiss.read_index(self.inspection.index_path, flags)
        expected = self.inspection.count
        if expected is not None and int(self._index.ntotal) != expected:
            raise RuntimeError(f"FAISS index has {self._index.ntotal} vectors; expected {expected}")
        if hasattr(self._index, "nprobe"):
            nlist = int(getattr(self._index, "nlist", self.nprobe))
            self._index.nprobe = min(max(1, self.nprobe), nlist)

        id_path = Path(self.inspection.id_map_path).resolve()
        uri = f"{id_path.as_uri()}?mode=ro&immutable=1"
        self._id_conn = sqlite3.connect(
            uri,
            uri=True,
            check_same_thread=False,
        )
        self._h5 = h5py.File(self.h5_path, "r")

    def close(self) -> None:
        with self._lock:
            if self._id_conn is not None:
                self._id_conn.close()
                self._id_conn = None
            if self._h5 is not None:
                self._h5.close()
                self._h5 = None
            self._index = None

    def __del__(self):  # pragma: no cover - deterministic callers can use close()
        with contextlib.suppress(Exception):
            self.close()

    def add(self, id: str, embedding: np.ndarray, metadata: dict | None = None) -> None:
        raise RuntimeError(
            "Persistent FAISSStore is read-only; update the canonical H5 and rebuild "
            "the sidecars instead."
        )

    def _row_for_id(self, protein_id: str) -> int | None:
        assert self._id_conn is not None
        row = self._id_conn.execute(
            "SELECT row_index FROM vector_ids WHERE protein_id = ?",
            (protein_id,),
        ).fetchone()
        return int(row[0]) if row is not None else None

    def _ids_for_rows(self, row_indices: list[int]) -> dict[int, str]:
        assert self._id_conn is not None
        result: dict[int, str] = {}
        unique = list(dict.fromkeys(row_indices))
        for start in range(0, len(unique), 900):
            chunk = unique[start : start + 900]
            placeholders = ",".join("?" for _ in chunk)
            rows = self._id_conn.execute(
                f"SELECT row_index, protein_id FROM vector_ids WHERE row_index IN ({placeholders})",
                chunk,
            ).fetchall()
            result.update((int(row), protein_id) for row, protein_id in rows)
        return result

    def _vector_for_row(self, row_index: int) -> np.ndarray:
        assert self._h5 is not None
        vector = np.asarray(
            self._h5["embeddings"][row_index],
            dtype=np.float32,
        ).reshape(1, -1)
        return _normalize(vector)

    def query(
        self,
        query_id: str,
        k: int = 10,
        threshold: float | None = None,
        include_distances: bool = False,
    ) -> list[str] | list[tuple[str, float]]:
        if not self.ready or k < 1:
            return []

        with self._lock:
            row_index = self._row_for_id(query_id)
            if row_index is None:
                return []
            vector = self._vector_for_row(row_index)
            scores, indices = self._index.search(vector, k + 1)
            valid_rows = [int(value) for value in indices[0] if value >= 0]
            ids = self._ids_for_rows(valid_rows)

        pairs: list[tuple[str, float]] = []
        for score, raw_index in zip(scores[0], indices[0], strict=True):
            result_index = int(raw_index)
            result_id = ids.get(result_index)
            if result_id is None or result_id == query_id:
                continue
            similarity = float(score)
            if threshold is not None and similarity < threshold:
                continue
            pairs.append((result_id, similarity))
            if len(pairs) == k:
                break

        if include_distances:
            return pairs
        return [protein_id for protein_id, _ in pairs]

    def query_vector(
        self,
        vector: np.ndarray,
        k: int = 10,
        threshold: float | None = None,
    ) -> list[tuple[str, float]]:
        if not self.ready or k < 1:
            return []
        query = np.asarray(vector, dtype=np.float32).reshape(1, -1)
        if self.inspection.dimension is not None and query.shape[1] != self.inspection.dimension:
            raise ValueError(
                f"Query dimension {query.shape[1]} does not match index dimension "
                f"{self.inspection.dimension}"
            )
        query = _normalize(query)

        with self._lock:
            scores, indices = self._index.search(query, k)
            valid_rows = [int(value) for value in indices[0] if value >= 0]
            ids = self._ids_for_rows(valid_rows)

        pairs: list[tuple[str, float]] = []
        for score, result_index in zip(scores[0], indices[0], strict=True):
            result_id = ids.get(int(result_index))
            if result_id is None:
                continue
            similarity = float(score)
            if threshold is not None and similarity < threshold:
                continue
            pairs.append((result_id, similarity))
        return pairs


__all__ = [
    "FAISSStore",
    "VectorIndexInspection",
    "VectorStore",
    "build_vector_index",
    "inspect_vector_index",
]
