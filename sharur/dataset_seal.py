"""Portable, disk-light integrity seals for completed Sharur datasets."""

from __future__ import annotations

import contextlib
import hashlib
import json
import os
import platform
import sqlite3
import subprocess
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

from sharur import __version__
from sharur.capabilities import CapabilityBrief, build_capability_brief
from sharur.ingest.input_integrity import write_json_atomic
from sharur.storage.duckdb_store import DuckDBStore


SEAL_SCHEMA_VERSION = 1
DEFAULT_SEAL_NAME = "dataset.seal.json"
DEFAULT_MAX_HASH_BYTES = 10 * 1024 * 1024
LARGE_FILE_SAMPLE_BYTES = 64 * 1024
LARGE_FILE_SAMPLE_COUNT = 9

_CANONICAL_PATTERNS = (
    "*/processing_manifest.json",
    "*/input_integrity.json",
    "*/summary_stats.json",
    "*/embedding_manifest.json",
    "*/protein_embeddings.h5",
    "*/protein_embeddings.index.json",
    "*/findings.jsonl",
)
_JSON_FIELDS = {
    "command",
    "config",
    "inputs",
    "outputs",
    "resource_profile",
    "result",
}


class DatasetSealError(RuntimeError):
    """Raised when a stable dataset seal cannot be constructed or inspected."""


@dataclass(frozen=True)
class SealVerification:
    """Comparison between a stored seal and the dataset currently on disk."""

    valid: bool
    expected_dataset_id: str
    observed_dataset_id: str
    seal_strength: str
    database_path: str
    changed_sections: tuple[str, ...]
    changed_artifacts: tuple[str, ...]
    provenance_changed: bool

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def to_markdown(self) -> str:
        status = "verified" if self.valid else "MISMATCH"
        lines = [
            "# Sharur Dataset Seal Verification",
            "",
            f"- Status: **{status}**",
            f"- Database: `{self.database_path}`",
            f"- Strength: `{self.seal_strength}`",
            f"- Expected ID: `{self.expected_dataset_id}`",
            f"- Observed ID: `{self.observed_dataset_id}`",
        ]
        if self.changed_sections:
            lines.append(f"- Changed sections: `{', '.join(self.changed_sections)}`")
        if self.changed_artifacts:
            lines.extend(
                [
                    "",
                    "## Changed canonical artifacts",
                    *[f"- `{path}`" for path in self.changed_artifacts],
                ]
            )
        if self.provenance_changed:
            lines.extend(
                [
                    "",
                    "> Runtime/software provenance differs from the sealing environment; "
                    "this does not by itself invalidate the dataset identity.",
                ]
            )
        return "\n".join(lines)


def _canonical_bytes(payload: Any) -> bytes:
    return json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        default=str,
    ).encode("utf-8")


def _payload_sha256(payload: Any) -> str:
    return hashlib.sha256(_canonical_bytes(payload)).hexdigest()


def _normalized_absolute(path: Path) -> Path:
    return Path(os.path.abspath(path.expanduser()))


def _relative_path(path: Path, root: Path) -> str:
    normalized = _normalized_absolute(path)
    try:
        return normalized.relative_to(root).as_posix()
    except ValueError as exc:
        raise DatasetSealError(
            f"Canonical artifact escapes the dataset directory: {normalized}"
        ) from exc


def _stable_file_digest(path: Path, *, full: bool) -> dict[str, Any]:
    before = path.stat()
    digest = hashlib.sha256()
    if full:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        result = {
            "algorithm": "sha256",
            "scope": "full",
            "value": digest.hexdigest(),
            "bytes_read": before.st_size,
        }
    else:
        sample_size = min(LARGE_FILE_SAMPLE_BYTES, before.st_size)
        span = max(0, before.st_size - sample_size)
        if LARGE_FILE_SAMPLE_COUNT == 1:
            offsets = [0]
        else:
            offsets = sorted(
                {
                    round(index * span / (LARGE_FILE_SAMPLE_COUNT - 1))
                    for index in range(LARGE_FILE_SAMPLE_COUNT)
                }
            )
        digest.update(f"size:{before.st_size}\n".encode())
        bytes_read = 0
        with path.open("rb") as handle:
            for offset in offsets:
                handle.seek(offset)
                chunk = handle.read(sample_size)
                digest.update(f"offset:{offset}:length:{len(chunk)}\n".encode())
                digest.update(chunk)
                bytes_read += len(chunk)
        result = {
            "algorithm": "sha256",
            "scope": "sampled",
            "value": digest.hexdigest(),
            "bytes_read": bytes_read,
            "sample_offsets": offsets,
            "sample_bytes": sample_size,
        }
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise DatasetSealError(f"Artifact changed while it was being sealed: {path}")
    return result


def verify_artifact_digest(
    path: str | Path,
    *,
    expected_size: int,
    expected_digest: dict[str, Any],
) -> bool:
    """Verify one canonical artifact against its stored seal record.

    This is intentionally narrower than :func:`verify_dataset_seal`: replica
    staging can validate the copied DuckDB artifact without reconstructing a
    second dataset directory around it.
    """

    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file() or resolved.stat().st_size != expected_size:
        return False
    if expected_digest.get("algorithm") != "sha256":
        raise DatasetSealError(
            f"Unsupported artifact digest algorithm: {expected_digest.get('algorithm')!r}"
        )
    scope = expected_digest.get("scope")
    if scope not in {"full", "sampled"}:
        raise DatasetSealError(f"Unsupported artifact digest scope: {scope!r}")
    observed = _stable_file_digest(resolved, full=scope == "full")
    return observed == expected_digest


def _artifact_snapshot(
    path: Path,
    *,
    root: Path,
    kind: str,
    hash_large_files: bool,
    max_hash_bytes: int,
) -> dict[str, Any]:
    size = path.stat().st_size
    digest = _stable_file_digest(
        path,
        full=hash_large_files or size <= max_hash_bytes,
    )
    return {
        "path": _relative_path(path, root),
        "kind": kind,
        "size": size,
        "digest": digest,
    }


def _read_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise DatasetSealError(
            f"Could not read canonical JSON artifact {path}: {type(exc).__name__}: {exc}"
        ) from exc


def _normalize_dataset_paths(value: Any, root: Path) -> Any:
    if isinstance(value, dict):
        return {
            str(key): _normalize_dataset_paths(item, root)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, list):
        return [_normalize_dataset_paths(item, root) for item in value]
    if isinstance(value, str):
        return value.replace(str(root), "<DATASET>")
    return value


def _stage00_manifest(root: Path) -> Path | None:
    preferred = root / "stage00_prepared" / "processing_manifest.json"
    if preferred.is_file():
        return preferred
    candidates = sorted(root.glob("stage00*/processing_manifest.json"))
    return candidates[0] if candidates else None


def _input_contract(root: Path) -> tuple[dict[str, Any], list[Path]]:
    manifest_path = _stage00_manifest(root)
    if manifest_path is None:
        return {"status": "absent", "assemblies": []}, []

    payload = _read_json(manifest_path)
    if not isinstance(payload, dict):
        raise DatasetSealError(f"Stage-00 manifest must contain a JSON object: {manifest_path}")

    report_path = manifest_path.parent / "input_integrity.json"
    report: dict[str, Any] = {}
    if report_path.is_file():
        loaded_report = _read_json(report_path)
        if not isinstance(loaded_report, dict):
            raise DatasetSealError(
                f"Stage-00 integrity report must contain a JSON object: {report_path}"
            )
        report = loaded_report

    report_files = {
        entry.get("filename"): entry
        for entry in report.get("files", [])
        if isinstance(entry, dict) and isinstance(entry.get("filename"), str)
    }
    assemblies: list[dict[str, Any]] = []
    prepared_paths: list[Path] = []
    genomes = payload.get("genomes", [])
    if not isinstance(genomes, list):
        raise DatasetSealError(f"Stage-00 manifest has a non-list genomes field: {manifest_path}")
    for entry in genomes:
        if not isinstance(entry, dict):
            continue
        filename = entry.get("filename")
        report_entry = report_files.get(filename, {})
        checksums = report_entry.get("checksums", {})
        sha256 = entry.get("sha256")
        if sha256 is None and isinstance(checksums, dict):
            sha256 = checksums.get("sha256")
        assemblies.append(
            {
                "filename": filename,
                "genome_id": entry.get("genome_id"),
                "file_size": entry.get("file_size"),
                "sha256": sha256,
                "legacy_checksum": entry.get("checksum"),
                "legacy_checksum_algorithm": entry.get("checksum_algorithm", "md5"),
                "format_valid": entry.get("format_valid"),
                "sequence_count": entry.get("sequence_count"),
                "total_length": entry.get("total_length"),
            }
        )
        if not isinstance(filename, str):
            continue
        output_path = entry.get("output_path")
        candidates: list[Path] = [manifest_path.parent / filename]
        if isinstance(output_path, str):
            candidates.insert(0, Path(output_path))
        prepared = next(
            (
                candidate
                for candidate in candidates
                if candidate.is_file()
                and _normalized_absolute(candidate).is_relative_to(root)
            ),
            None,
        )
        if prepared is not None:
            prepared_paths.append(_normalized_absolute(prepared))

    assemblies.sort(key=lambda item: (str(item["genome_id"]), str(item["filename"])))
    contract = {
        "status": payload.get("integrity_status", report.get("status", "legacy")),
        "manifest": _relative_path(manifest_path, root),
        "manifest_version": payload.get("version"),
        "source_set_sha256": payload.get("source_set_sha256"),
        "summary": payload.get("summary", report.get("summary", {})),
        "assemblies": assemblies,
    }
    if report_path.is_file():
        contract["integrity_report"] = _relative_path(report_path, root)
        contract["integrity_schema_version"] = report.get("schema_version")
    return contract, prepared_paths


def _database_snapshot(db_path: Path) -> dict[str, Any]:
    try:
        with DuckDBStore(db_path, read_only=True) as store:
            table_rows = store.execute(
                """
                SELECT table_name
                FROM information_schema.tables
                WHERE table_schema = 'main' AND table_type = 'BASE TABLE'
                ORDER BY table_name
                """
            )
            tables = [row[0] for row in table_rows]
            columns = [
                {
                    "table": table,
                    "ordinal": int(ordinal),
                    "column": column,
                    "data_type": data_type,
                    "nullable": nullable == "YES",
                    "default": default,
                }
                for table, ordinal, column, data_type, nullable, default in store.execute(
                    """
                    SELECT table_name, ordinal_position, column_name, data_type,
                           is_nullable, column_default
                    FROM information_schema.columns
                    WHERE table_schema = 'main'
                    ORDER BY table_name, ordinal_position
                    """
                )
            ]
            counts = []
            for table in tables:
                escaped = table.replace('"', '""')
                count = int(store.execute(f'SELECT COUNT(*) FROM "{escaped}"')[0][0])
                counts.append({"table": table, "rows": count})

            schema_version = None
            if "schema_version" in tables:
                rows = store.execute("SELECT MAX(version) FROM schema_version")
                if rows and rows[0][0] is not None:
                    schema_version = int(rows[0][0])
    except Exception as exc:
        raise DatasetSealError(
            f"DuckDB could not be sealed: {type(exc).__name__}: {exc}"
        ) from exc

    return {
        "database_file": db_path.name,
        "schema_version": schema_version,
        "tables": counts,
        "columns": columns,
    }


def _capability_identity(brief: CapabilityBrief) -> dict[str, Any]:
    capabilities = {
        item.capability_id: item
        for item in brief.capabilities
    }
    annotation = capabilities.get("annotation_sources")
    callers = capabilities.get("structured_callers")
    return {
        "annotation_sources": (
            annotation.evidence.get("sources", []) if annotation is not None else []
        ),
        "structured_callers": (
            callers.evidence.get("resources", []) if callers is not None else []
        ),
    }


def _stage_manifests(root: Path) -> list[dict[str, Any]]:
    summaries: list[dict[str, Any]] = []
    for path in sorted(root.glob("*/processing_manifest.json")):
        payload = _read_json(path)
        if not isinstance(payload, dict):
            raise DatasetSealError(f"Processing manifest must contain an object: {path}")
        summaries.append(
            {
                "path": _relative_path(path, root),
                "stage": payload.get("stage", path.parent.name),
                "version": payload.get("version"),
                "source_set_sha256": payload.get("source_set_sha256"),
                "integrity_status": payload.get("integrity_status"),
            }
        )
    return summaries


def _coerce_local_artifact(value: Any, *, manifest: Path, root: Path) -> Path | None:
    if not isinstance(value, str) or not value:
        return None
    raw = Path(value).expanduser()
    candidates = (
        [raw]
        if raw.is_absolute()
        else [manifest.parent / raw, root / raw]
    )
    for candidate in candidates:
        normalized = _normalized_absolute(candidate)
        if normalized.is_relative_to(root) and normalized.is_file():
            return normalized
    return None


def _discover_artifacts(
    root: Path,
    db_path: Path,
    prepared_inputs: list[Path],
) -> dict[Path, str]:
    artifacts: dict[Path, str] = {}

    def add(path: Path, kind: str) -> None:
        normalized = _normalized_absolute(path)
        if normalized.is_relative_to(root) and normalized.is_file():
            artifacts.setdefault(normalized, kind)

    add(db_path, "duckdb")
    add(Path(f"{db_path}.wal"), "duckdb_wal")
    add(root / "manifest.json", "analysis_manifest")
    for pattern in _CANONICAL_PATTERNS:
        for path in root.glob(pattern):
            kind = {
                "findings.jsonl": "canonical_findings",
                "processing_manifest.json": "stage_manifest",
                "input_integrity.json": "input_integrity",
                "summary_stats.json": "stage_summary",
                "embedding_manifest.json": "embedding_manifest",
                "protein_embeddings.h5": "protein_embeddings",
                "protein_embeddings.index.json": "vector_index_manifest",
            }.get(path.name, "canonical_artifact")
            add(path, kind)
    for path in prepared_inputs:
        add(path, "prepared_input")

    json_manifests = [
        path
        for path in artifacts
        if path.name in {"embedding_manifest.json", "protein_embeddings.index.json"}
    ]
    for manifest in json_manifests:
        payload = _read_json(manifest)
        if not isinstance(payload, dict):
            continue
        output_files = payload.get("output_files", {})
        if isinstance(output_files, dict):
            for key, value in output_files.items():
                path = _coerce_local_artifact(value, manifest=manifest, root=root)
                if path is not None:
                    add(path, f"embedding_{key}")
        for key, kind in (
            ("index_file", "vector_index"),
            ("id_map_file", "vector_id_map"),
        ):
            path = _coerce_local_artifact(payload.get(key), manifest=manifest, root=root)
            if path is not None:
                add(path, kind)
    return artifacts


def _decode_sqlite_row(row: sqlite3.Row) -> dict[str, Any]:
    decoded = dict(row)
    for field in _JSON_FIELDS:
        value = decoded.get(field)
        if not isinstance(value, str):
            continue
        with contextlib.suppress(json.JSONDecodeError):
            decoded[field] = json.loads(value)
    return decoded


def _ledger_snapshot(root: Path) -> tuple[dict[str, Any], dict[str, Any] | None]:
    ledger_path = root / "sharur_ops.db"
    if not ledger_path.is_file():
        return {"status": "absent"}, None

    try:
        uri = f"{ledger_path.resolve().as_uri()}?mode=ro"
        connection = sqlite3.connect(uri, uri=True)
        connection.row_factory = sqlite3.Row
        try:
            tables = {
                row[0]
                for row in connection.execute(
                    "SELECT name FROM sqlite_master WHERE type = 'table'"
                )
            }
            run_count = (
                int(connection.execute("SELECT COUNT(*) FROM runs").fetchone()[0])
                if "runs" in tables
                else 0
            )
            latest = None
            completed = None
            active = None
            if "runs" in tables:
                row = connection.execute(
                    """
                    SELECT * FROM runs
                    WHERE run_type = 'ingest'
                    ORDER BY created_ts DESC
                    LIMIT 1
                    """
                ).fetchone()
                latest = _decode_sqlite_row(row) if row is not None else None
                row = connection.execute(
                    """
                    SELECT * FROM runs
                    WHERE run_type = 'ingest'
                      AND status IN ('pending', 'submitted', 'running')
                    ORDER BY created_ts DESC
                    LIMIT 1
                    """
                ).fetchone()
                active = _decode_sqlite_row(row) if row is not None else None
                row = connection.execute(
                    """
                    SELECT * FROM runs
                    WHERE run_type = 'ingest' AND status = 'complete'
                    ORDER BY completed_ts DESC, created_ts DESC
                    LIMIT 1
                    """
                ).fetchone()
                completed = _decode_sqlite_row(row) if row is not None else None

            identity_stages: list[dict[str, Any]] = []
            provenance_stages: list[dict[str, Any]] = []
            if completed is not None and "run_stages" in tables:
                rows = connection.execute(
                    """
                    SELECT stage_id, attempt, signature, status, resource_profile,
                           reused_from_run_id, reused_from_attempt
                    FROM run_stages
                    WHERE run_id = ?
                    ORDER BY stage_id, attempt
                    """,
                    (completed["id"],),
                )
                decoded_stages = [_decode_sqlite_row(row) for row in rows]
                latest_by_stage: dict[str, dict[str, Any]] = {}
                for stage in decoded_stages:
                    latest_by_stage[stage["stage_id"]] = stage
                provenance_stages = [
                    _normalize_dataset_paths(stage, root)
                    for stage in latest_by_stage.values()
                ]
                identity_stages = [
                    {
                        "stage_id": stage["stage_id"],
                        "signature": stage["signature"],
                        "status": stage["status"],
                    }
                    for stage in latest_by_stage.values()
                ]
        finally:
            connection.close()
    except Exception as exc:
        raise DatasetSealError(
            f"Operational ledger could not be inspected: {type(exc).__name__}: {exc}"
        ) from exc

    if active is not None:
        raise DatasetSealError(
            "Refusing to seal while an ingest run is active "
            f"({active['id']}: {active['status']})"
        )

    operational = {
        "status": "available",
        "path": ledger_path.name,
        "tables": sorted(tables),
        "runs": run_count,
        "latest_ingest": (
            {
                "id": latest.get("id"),
                "status": latest.get("status"),
                "created_ts": latest.get("created_ts"),
                "completed_ts": latest.get("completed_ts"),
            }
            if latest is not None
            else None
        ),
        "latest_completed_ingest": (
            {
                "id": completed.get("id"),
                "completed_ts": completed.get("completed_ts"),
                "config": _normalize_dataset_paths(completed.get("config", {}), root),
                "stages": provenance_stages,
            }
            if completed is not None
            else None
        ),
    }
    if completed is None:
        return operational, None
    stable_run = {
        "stages": identity_stages,
    }
    return operational, stable_run


def _config_fingerprints() -> list[dict[str, Any]]:
    package_root = Path(__file__).resolve().parent
    project_root = package_root.parent
    roots = (
        (
            "predicates_v2",
            next(
                (
                    path
                    for path in (
                        package_root / "predicates_v2" / "config",
                        project_root / "config" / "predicates_v2",
                    )
                    if path.is_dir()
                ),
                None,
            ),
        ),
        ("predicate_mappings", package_root / "predicates" / "mappings" / "data"),
    )
    records: list[dict[str, Any]] = []
    for scope, directory in roots:
        if directory is None or not directory.is_dir():
            continue
        for path in sorted(directory.rglob("*")):
            if not path.is_file() or path.suffix.lower() not in {
                ".json",
                ".tsv",
                ".yaml",
                ".yml",
            }:
                continue
            records.append(
                {
                    "path": f"{scope}/{path.relative_to(directory).as_posix()}",
                    "sha256": _stable_file_digest(path, full=True)["value"],
                }
            )
    return records


def _git_provenance() -> dict[str, Any]:
    project_root = Path(__file__).resolve().parent.parent
    try:
        revision = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=project_root,
            check=True,
            capture_output=True,
            text=True,
            timeout=5,
        ).stdout.strip()
        status = subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=no"],
            cwd=project_root,
            check=True,
            capture_output=True,
            text=True,
            timeout=5,
        ).stdout
    except (OSError, subprocess.SubprocessError):
        return {"available": False}
    return {
        "available": True,
        "revision": revision,
        "tracked_worktree_dirty": bool(status.strip()),
    }


def _package_versions() -> dict[str, str]:
    packages = ("duckdb", "numpy", "pandas", "pydantic", "sharur")
    result: dict[str, str] = {}
    for package in packages:
        try:
            result[package] = version(package)
        except PackageNotFoundError:
            result[package] = "not-installed"
    return result


def _capability_provenance(brief: CapabilityBrief, root: Path) -> dict[str, Any]:
    return {
        "overall_state": brief.overall_state.value,
        "capabilities": [
            _normalize_dataset_paths(capability.to_dict(), root)
            for capability in brief.capabilities
        ],
    }


def _software_provenance() -> dict[str, Any]:
    return {
        "sharur_version": __version__,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "packages": _package_versions(),
        "git": _git_provenance(),
        "configuration_files": _config_fingerprints(),
    }


def _database_file_state(db_path: Path) -> dict[str, tuple[int, int]]:
    states: dict[str, tuple[int, int]] = {}
    for path in (db_path, Path(f"{db_path}.wal")):
        if path.is_file():
            stat = path.stat()
            states[str(path)] = (stat.st_size, stat.st_mtime_ns)
    return states


def build_dataset_seal(
    db_path: str | Path,
    *,
    hash_large_files: bool = False,
    max_hash_bytes: int = DEFAULT_MAX_HASH_BYTES,
    include_tools: bool = False,
) -> dict[str, Any]:
    """Build a read-only seal over one dataset's canonical scientific state."""
    if max_hash_bytes < 0:
        raise ValueError("max_hash_bytes must be non-negative")
    resolved_db = Path(db_path).expanduser().resolve()
    if not resolved_db.is_file():
        raise DatasetSealError(f"DuckDB file does not exist: {resolved_db}")
    root = resolved_db.parent
    initial_db_state = _database_file_state(resolved_db)

    inputs, prepared_inputs = _input_contract(root)
    database = _database_snapshot(resolved_db)
    brief = build_capability_brief(
        resolved_db,
        include_tools=include_tools,
        include_execution=False,
    )
    operations, completed_ingest = _ledger_snapshot(root)
    artifact_paths = _discover_artifacts(root, resolved_db, prepared_inputs)
    artifacts = [
        _artifact_snapshot(
            path,
            root=root,
            kind=kind,
            hash_large_files=hash_large_files,
            max_hash_bytes=max_hash_bytes,
        )
        for path, kind in sorted(
            artifact_paths.items(),
            key=lambda pair: _relative_path(pair[0], root),
        )
    ]
    if initial_db_state != _database_file_state(resolved_db):
        raise DatasetSealError("DuckDB or its WAL changed while the seal was being built")

    identity = {
        "inputs": inputs,
        "database": database,
        "annotation_contract": _capability_identity(brief),
        "stage_manifests": _stage_manifests(root),
        "completed_ingest": completed_ingest,
        "canonical_artifacts": artifacts,
    }
    dataset_id = _payload_sha256(identity)
    policy = {
        "artifact_scope": "canonical-v1",
        "hash_large_files": hash_large_files,
        "max_hash_bytes": max_hash_bytes,
        "large_file_sample_bytes": LARGE_FILE_SAMPLE_BYTES,
        "large_file_sample_count": LARGE_FILE_SAMPLE_COUNT,
        "include_tools": include_tools,
    }
    return {
        "schema_version": SEAL_SCHEMA_VERSION,
        "seal_type": "sharur_dataset",
        "seal_strength": "content" if hash_large_files else "structural",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "dataset_id": dataset_id,
        "dataset": {
            "root": str(root),
            "database": resolved_db.name,
        },
        "policy": policy,
        "identity": identity,
        "provenance": {
            "software": _software_provenance(),
            "capability_brief": _capability_provenance(brief, root),
            "operations": operations,
        },
    }


def write_dataset_seal(
    seal: dict[str, Any],
    path: str | Path,
    *,
    force: bool = False,
) -> Path:
    """Atomically write a seal, refusing replacement unless explicitly forced."""
    output = Path(path).expanduser().resolve()
    if output.exists() and not force:
        raise FileExistsError(f"Seal already exists: {output}")
    write_json_atomic(seal, output)
    return output


def seal_dataset(
    db_path: str | Path,
    *,
    output_path: str | Path | None = None,
    hash_large_files: bool = False,
    max_hash_bytes: int = DEFAULT_MAX_HASH_BYTES,
    include_tools: bool = False,
    force: bool = False,
) -> tuple[Path, dict[str, Any]]:
    """Build and atomically save one dataset seal."""
    resolved_db = Path(db_path).expanduser().resolve()
    output = (
        Path(output_path).expanduser().resolve()
        if output_path is not None
        else resolved_db.parent / DEFAULT_SEAL_NAME
    )
    if output.exists() and not force:
        raise FileExistsError(f"Seal already exists: {output}")
    seal = build_dataset_seal(
        resolved_db,
        hash_large_files=hash_large_files,
        max_hash_bytes=max_hash_bytes,
        include_tools=include_tools,
    )
    return write_dataset_seal(seal, output, force=force), seal


def _verification_database(
    stored: dict[str, Any],
    seal_path: Path,
    explicit_db: str | Path | None,
) -> Path:
    if explicit_db is not None:
        return Path(explicit_db).expanduser().resolve()
    dataset = stored.get("dataset")
    if not isinstance(dataset, dict) or not isinstance(dataset.get("database"), str):
        raise DatasetSealError("Seal does not contain a usable dataset.database field")
    relative_database = Path(dataset["database"])
    beside_seal = seal_path.parent / relative_database
    if beside_seal.is_file():
        return beside_seal.resolve()
    root = dataset.get("root")
    if isinstance(root, str):
        original = Path(root).expanduser() / relative_database
        if original.is_file():
            return original.resolve()
    raise DatasetSealError(
        "Could not locate the sealed DuckDB; pass --db with its current path"
    )


def _changed_artifacts(expected: Any, observed: Any) -> tuple[str, ...]:
    def by_path(value: Any) -> dict[str, Any]:
        if not isinstance(value, list):
            return {}
        return {
            item["path"]: item
            for item in value
            if isinstance(item, dict) and isinstance(item.get("path"), str)
        }

    expected_paths = by_path(expected)
    observed_paths = by_path(observed)
    return tuple(
        path
        for path in sorted(expected_paths.keys() | observed_paths.keys())
        if expected_paths.get(path) != observed_paths.get(path)
    )


def verify_dataset_seal(
    seal_path: str | Path,
    *,
    db_path: str | Path | None = None,
) -> SealVerification:
    """Recompute a stored seal policy and compare the canonical dataset identity."""
    resolved_seal = Path(seal_path).expanduser().resolve()
    if not resolved_seal.is_file():
        raise DatasetSealError(f"Seal file does not exist: {resolved_seal}")
    stored = _read_json(resolved_seal)
    if not isinstance(stored, dict):
        raise DatasetSealError("Seal must contain a JSON object")
    if stored.get("schema_version") != SEAL_SCHEMA_VERSION:
        raise DatasetSealError(
            f"Unsupported seal schema version: {stored.get('schema_version')!r}"
        )
    if stored.get("seal_type") != "sharur_dataset":
        raise DatasetSealError(f"Unsupported seal type: {stored.get('seal_type')!r}")
    policy = stored.get("policy")
    if not isinstance(policy, dict):
        raise DatasetSealError("Seal does not contain a usable policy")
    database = _verification_database(stored, resolved_seal, db_path)
    observed = build_dataset_seal(
        database,
        hash_large_files=bool(policy.get("hash_large_files", False)),
        max_hash_bytes=int(policy.get("max_hash_bytes", DEFAULT_MAX_HASH_BYTES)),
        include_tools=bool(policy.get("include_tools", False)),
    )

    expected_identity = stored.get("identity")
    observed_identity = observed["identity"]
    if not isinstance(expected_identity, dict):
        raise DatasetSealError("Seal does not contain a usable identity object")
    changed_sections = tuple(
        section
        for section in sorted(expected_identity.keys() | observed_identity.keys())
        if expected_identity.get(section) != observed_identity.get(section)
    )
    artifacts = _changed_artifacts(
        expected_identity.get("canonical_artifacts"),
        observed_identity.get("canonical_artifacts"),
    )
    expected_id = str(stored.get("dataset_id", ""))
    observed_id = str(observed["dataset_id"])
    return SealVerification(
        valid=expected_id == observed_id and not changed_sections,
        expected_dataset_id=expected_id,
        observed_dataset_id=observed_id,
        seal_strength=str(stored.get("seal_strength", "unknown")),
        database_path=str(database),
        changed_sections=changed_sections,
        changed_artifacts=artifacts,
        provenance_changed=stored.get("provenance") != observed.get("provenance"),
    )


__all__ = [
    "DEFAULT_MAX_HASH_BYTES",
    "DEFAULT_SEAL_NAME",
    "SEAL_SCHEMA_VERSION",
    "DatasetSealError",
    "SealVerification",
    "build_dataset_seal",
    "seal_dataset",
    "verify_artifact_digest",
    "verify_dataset_seal",
    "write_dataset_seal",
]
