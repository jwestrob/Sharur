"""
Helpers for writing canonical findings records to JSONL files.

Writers should use this module instead of manually serializing findings so
schema_version, phase, and provenance normalization are applied consistently.
"""

from __future__ import annotations

import fcntl
import json
import os
import uuid
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping

from sharur.core.analysis_record_compat import (
    FindingNormalizationResult,
    normalize_finding,
)
from sharur.core.analysis_records import canonicalize_finding_phase


class FindingValidationError(ValueError):
    """Raised when a new canonical finding violates the write contract."""


def _count_existing_records(path: Path) -> int:
    if not path.exists():
        return 0
    with path.open() as handle:
        return sum(1 for line in handle if line.strip())


@contextmanager
def _exclusive_findings_lock(path: Path) -> Iterator[None]:
    """Serialize JSONL writers across agent processes."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lock_path = path.with_name(f".{path.name}.lock")
    with lock_path.open("a+") as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)


def _prepare_with_ordinal(
    finding: dict[str, Any],
    path: Path,
    *,
    phase: str | None,
    ordinal: int,
) -> FindingNormalizationResult:
    raw = dict(finding)
    normalized_phase = canonicalize_finding_phase(phase)
    if normalized_phase and "phase" not in raw:
        raw["phase"] = normalized_phase
    result = normalize_finding(raw, source_path=path, ordinal=ordinal)
    return _apply_write_contract(raw, result)


def _apply_write_contract(
    raw: Mapping[str, Any],
    result: FindingNormalizationResult,
) -> FindingNormalizationResult:
    """Add strict new-write requirements without breaking legacy reads."""
    issues = list(result.issues)
    for field in ("title", "category", "description"):
        value = raw.get(field)
        if not isinstance(value, str) or not value.strip():
            issues.append(f"missing {field}")

    if "evidence" not in raw:
        issues.append("missing evidence")
    elif not isinstance(raw.get("evidence"), Mapping):
        issues.append("evidence must be an object")

    novelty = result.finding.get("novelty", 0)
    if novelty >= 2:
        falsification = raw.get("falsification")
        if not isinstance(falsification, str) or not falsification.strip():
            issues.append("missing falsification for novelty >= 2")

    return FindingNormalizationResult(
        finding=result.finding,
        issues=tuple(dict.fromkeys(issues)),
        is_key_finding=result.is_key_finding,
        source_path=result.source_path,
    )


def _raise_for_issues(result: FindingNormalizationResult) -> None:
    if result.issues:
        raise FindingValidationError(
            f"Finding {result.finding.get('id')} failed validation: "
            + "; ".join(result.issues)
        )


def _read_existing_ids(handle) -> tuple[int, set[str]]:
    handle.seek(0)
    count = 0
    finding_ids: set[str] = set()
    for line_number, line in enumerate(handle, start=1):
        if not line.strip():
            continue
        count += 1
        try:
            record = json.loads(line)
        except json.JSONDecodeError as exc:
            raise FindingValidationError(
                f"Cannot append to malformed JSONL at line {line_number}: {exc}"
            ) from exc
        finding_id = record.get("id") if isinstance(record, dict) else None
        if finding_id:
            finding_ids.add(str(finding_id))
    return count, finding_ids


def prepare_finding_for_write(
    finding: dict[str, Any],
    path: str | Path,
    *,
    phase: str | None = None,
    ordinal: int | None = None,
) -> FindingNormalizationResult:
    """Normalize a finding before writing it to disk."""
    path = Path(path)
    next_ordinal = ordinal or (_count_existing_records(path) + 1)
    return _prepare_with_ordinal(
        finding,
        path,
        phase=phase,
        ordinal=next_ordinal,
    )


def append_finding_record(
    path: str | Path,
    finding: dict[str, Any],
    *,
    phase: str | None = None,
    strict: bool = True,
) -> FindingNormalizationResult:
    """Append one normalized finding under an inter-process file lock.

    Strict mode is the default for canonical archives: validation issues abort
    the write. Set ``strict=False`` only for explicitly non-canonical draft
    spools that will pass through a strict merge step later.
    """
    path = Path(path)
    with _exclusive_findings_lock(path):
        with path.open("a+", encoding="utf-8") as handle:
            record_count, existing_ids = _read_existing_ids(handle)
            ordinal = record_count + 1
            result = _prepare_with_ordinal(
                finding,
                path,
                phase=phase,
                ordinal=ordinal,
            )

            if not finding.get("id"):
                while str(result.finding["id"]) in existing_ids:
                    ordinal += 1
                    result = _prepare_with_ordinal(
                        finding,
                        path,
                        phase=phase,
                        ordinal=ordinal,
                    )

            finding_id = str(result.finding["id"])
            if finding_id in existing_ids:
                raise FindingValidationError(
                    f"Finding ID {finding_id!r} already exists in {path}"
                )
            if strict:
                _raise_for_issues(result)

            handle.seek(0, os.SEEK_END)
            handle.write(json.dumps(result.finding, default=str) + "\n")
            handle.flush()
            os.fsync(handle.fileno())
    return result


def replace_finding_record(
    path: str | Path,
    finding_id: str,
    finding: dict[str, Any],
    *,
    phase: str | None = None,
    strict: bool = True,
) -> FindingNormalizationResult:
    """Atomically replace one canonical finding while preserving line order.

    This is the safe correction path for an existing JSONL archive. Unrelated
    records, including legacy records, are retained byte-for-byte.
    """
    path = Path(path)
    target_id = str(finding_id)
    with _exclusive_findings_lock(path):
        if not path.exists():
            raise FindingValidationError(f"Findings archive does not exist: {path}")

        lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
        record_ordinal = 0
        target_line_indexes: list[tuple[int, int]] = []
        for line_index, line in enumerate(lines):
            if not line.strip():
                continue
            record_ordinal += 1
            try:
                record = json.loads(line)
            except json.JSONDecodeError as exc:
                raise FindingValidationError(
                    f"Cannot replace record in malformed JSONL at "
                    f"line {line_index + 1}: {exc}"
                ) from exc
            record_id = record.get("id") if isinstance(record, dict) else None
            if record_id is not None and str(record_id) == target_id:
                target_line_indexes.append((line_index, record_ordinal))

        if not target_line_indexes:
            raise FindingValidationError(
                f"Finding ID {target_id!r} does not exist in {path}"
            )
        if len(target_line_indexes) > 1:
            raise FindingValidationError(
                f"Finding ID {target_id!r} occurs multiple times in {path}"
            )

        raw = dict(finding)
        supplied_id = raw.get("id")
        if supplied_id is not None and str(supplied_id) != target_id:
            raise FindingValidationError(
                f"Replacement ID {supplied_id!r} does not match target "
                f"{target_id!r}"
            )
        raw["id"] = target_id

        line_index, ordinal = target_line_indexes[0]
        result = _prepare_with_ordinal(
            raw,
            path,
            phase=phase,
            ordinal=ordinal,
        )
        if strict:
            _raise_for_issues(result)

        newline = "\n" if lines[line_index].endswith("\n") else ""
        lines[line_index] = (
            json.dumps(result.finding, default=str) + newline
        )

        temporary_path = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
        with temporary_path.open("w", encoding="utf-8") as handle:
            handle.writelines(lines)
            handle.flush()
            os.fsync(handle.fileno())
        temporary_path.replace(path)

    return result


def write_findings_records(
    path: str | Path,
    findings: Iterable[dict[str, Any]],
    *,
    phase: str | None = None,
    strict: bool = True,
) -> list[FindingNormalizationResult]:
    """Atomically replace a JSONL file with normalized finding records."""
    path = Path(path)
    materialized_findings = list(findings)
    results: list[FindingNormalizationResult] = []
    seen_ids: set[str] = set()
    for ordinal, finding in enumerate(materialized_findings, start=1):
        result = _prepare_with_ordinal(
            finding,
            path,
            phase=phase,
            ordinal=ordinal,
        )
        finding_id = str(result.finding["id"])
        if finding_id in seen_ids:
            raise FindingValidationError(
                f"Duplicate finding ID {finding_id!r} in replacement payload"
            )
        if strict:
            _raise_for_issues(result)
        seen_ids.add(finding_id)
        results.append(result)

    with _exclusive_findings_lock(path):
        temporary_path = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
        with temporary_path.open("w", encoding="utf-8") as handle:
            for result in results:
                handle.write(json.dumps(result.finding, default=str) + "\n")
            handle.flush()
            os.fsync(handle.fileno())
        temporary_path.replace(path)

    return results


__all__ = [
    "append_finding_record",
    "FindingValidationError",
    "prepare_finding_for_write",
    "replace_finding_record",
    "write_findings_records",
]
