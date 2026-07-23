"""Atomic campaign-local staging for sealed Sharur DuckDB databases."""

from __future__ import annotations

import contextlib
import fcntl
import hashlib
import json
import os
import shutil
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from sharur.dataset_seal import (
    DEFAULT_SEAL_NAME,
    DatasetSealError,
    verify_artifact_digest,
    verify_dataset_seal,
)
from sharur.ingest.input_integrity import write_json_atomic


DEFAULT_RESERVE_BYTES = 10 * 1024**3
COPY_CHUNK_BYTES = 16 * 1024**2


@dataclass(frozen=True)
class StagedDatabase:
    """Identity and location of one immutable campaign-local replica."""

    path: Path
    source_path: Path
    seal_path: Path
    dataset_id: str
    seal_strength: str
    artifact_digest: dict[str, Any]
    reused: bool
    staged_at: str

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        for key in ("path", "source_path", "seal_path"):
            payload[key] = str(payload[key])
        return payload


def _read_seal(seal_path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(seal_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise DatasetSealError(f"Could not read dataset seal {seal_path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise DatasetSealError(f"Dataset seal must contain an object: {seal_path}")
    return payload


def _database_artifact(
    seal: dict[str, Any],
    source_path: Path,
) -> dict[str, Any]:
    dataset = seal.get("dataset")
    identity = seal.get("identity")
    if not isinstance(dataset, dict) or not isinstance(identity, dict):
        raise DatasetSealError("Dataset seal lacks dataset identity fields")
    database_name = dataset.get("database")
    artifacts = identity.get("canonical_artifacts")
    if not isinstance(database_name, str) or not isinstance(artifacts, list):
        raise DatasetSealError("Dataset seal lacks a canonical database artifact")
    if source_path.name != Path(database_name).name:
        raise DatasetSealError(
            f"Seal records database {database_name!r}, received {source_path.name!r}"
        )
    for artifact in artifacts:
        if isinstance(artifact, dict) and artifact.get("path") == database_name:
            if not isinstance(artifact.get("size"), int) or not isinstance(
                artifact.get("digest"), dict
            ):
                break
            return artifact
    raise DatasetSealError(
        f"Canonical artifact record for {database_name!r} is absent or malformed"
    )


def _metadata_path(database_path: Path) -> Path:
    return database_path.with_suffix(f"{database_path.suffix}.stage.json")


def _lock_path(database_path: Path) -> Path:
    return database_path.with_suffix(f"{database_path.suffix}.stage.lock")


def _existing_replica_matches(
    database_path: Path,
    artifact: dict[str, Any],
) -> bool:
    return verify_artifact_digest(
        database_path,
        expected_size=int(artifact["size"]),
        expected_digest=artifact["digest"],
    )


def _fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def stage_database(
    source_path: str | Path,
    stage_dir: str | Path,
    *,
    seal_path: str | Path | None = None,
    reserve_bytes: int = DEFAULT_RESERVE_BYTES,
    verify_source: bool = True,
) -> StagedDatabase:
    """Verify and atomically stage one sealed DuckDB for a query campaign.

    A dataset-ID-qualified target is shared by all query-service launches that
    use the same staging directory. The copy is published with ``os.replace``
    only after its sealed artifact digest matches.
    """

    source = Path(source_path).expanduser().resolve()
    destination_dir = Path(stage_dir).expanduser().resolve()
    resolved_seal = (
        Path(seal_path).expanduser().resolve()
        if seal_path is not None
        else source.parent / DEFAULT_SEAL_NAME
    )
    if reserve_bytes < 0:
        raise ValueError("reserve_bytes must be non-negative")
    if not source.is_file():
        raise FileNotFoundError(f"DuckDB file does not exist: {source}")
    if not resolved_seal.is_file():
        raise DatasetSealError(
            f"Completed datasets require {DEFAULT_SEAL_NAME} before staging: "
            f"{resolved_seal}"
        )

    seal = _read_seal(resolved_seal)
    dataset_id = seal.get("dataset_id")
    if not isinstance(dataset_id, str) or len(dataset_id) < 12:
        raise DatasetSealError("Dataset seal lacks a usable dataset_id")
    artifact = _database_artifact(seal, source)
    if verify_source:
        verification = verify_dataset_seal(resolved_seal, db_path=source)
        if not verification.valid:
            changed = ", ".join(verification.changed_sections) or "dataset identity"
            raise DatasetSealError(
                f"Source dataset seal verification failed; changed: {changed}"
            )
    elif not _existing_replica_matches(source, artifact):
        raise DatasetSealError("Source DuckDB does not match its sealed artifact digest")

    destination_dir.mkdir(parents=True, exist_ok=True)
    target = destination_dir / f"{source.stem}-{dataset_id[:12]}.duckdb"
    lock_path = _lock_path(target)
    staged_at = datetime.now(timezone.utc).isoformat()

    with lock_path.open("a+b") as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        try:
            if _existing_replica_matches(target, artifact):
                metadata_path = _metadata_path(target)
                if metadata_path.is_file():
                    with contextlib.suppress(OSError, json.JSONDecodeError):
                        metadata = json.loads(metadata_path.read_text())
                        staged_at = str(metadata.get("staged_at") or staged_at)
                return StagedDatabase(
                    path=target,
                    source_path=source,
                    seal_path=resolved_seal,
                    dataset_id=dataset_id,
                    seal_strength=str(seal.get("seal_strength", "unknown")),
                    artifact_digest=artifact["digest"],
                    reused=True,
                    staged_at=staged_at,
                )

            required = int(artifact["size"]) + reserve_bytes
            available = shutil.disk_usage(destination_dir).free
            if available < required:
                raise OSError(
                    f"Staging requires {required:,} free bytes including reserve; "
                    f"{available:,} bytes are available in {destination_dir}"
                )

            partial = target.with_name(
                f".{target.name}.partial-{os.getpid()}-{hashlib.sha256(os.urandom(16)).hexdigest()[:8]}"
            )
            source_before = source.stat()
            try:
                with source.open("rb") as reader, partial.open("xb") as writer:
                    shutil.copyfileobj(reader, writer, length=COPY_CHUNK_BYTES)
                    writer.flush()
                    os.fsync(writer.fileno())
                source_after = source.stat()
                if (
                    source_before.st_size,
                    source_before.st_mtime_ns,
                ) != (
                    source_after.st_size,
                    source_after.st_mtime_ns,
                ):
                    raise DatasetSealError(
                        f"Source DuckDB changed during replica staging: {source}"
                    )
                if not _existing_replica_matches(partial, artifact):
                    raise DatasetSealError(
                        "Staged DuckDB failed its sealed artifact digest check"
                    )
                partial.chmod(0o444)
                os.replace(partial, target)
                _fsync_directory(destination_dir)
            finally:
                with contextlib.suppress(FileNotFoundError):
                    partial.unlink()

            result = StagedDatabase(
                path=target,
                source_path=source,
                seal_path=resolved_seal,
                dataset_id=dataset_id,
                seal_strength=str(seal.get("seal_strength", "unknown")),
                artifact_digest=artifact["digest"],
                reused=False,
                staged_at=staged_at,
            )
            write_json_atomic(result.to_dict(), _metadata_path(target))
            return result
        finally:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)


__all__ = [
    "COPY_CHUNK_BYTES",
    "DEFAULT_RESERVE_BYTES",
    "StagedDatabase",
    "stage_database",
]
