"""
Helpers for writing canonical findings records to JSONL files.

Writers should use this module instead of manually serializing findings so
schema_version, phase, and provenance normalization are applied consistently.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Iterable

from sharur.core.analysis_record_compat import (
    FindingNormalizationResult,
    normalize_finding,
)
from sharur.core.analysis_records import canonicalize_finding_phase


def _count_existing_records(path: Path) -> int:
    if not path.exists():
        return 0
    with path.open() as handle:
        return sum(1 for line in handle if line.strip())


def prepare_finding_for_write(
    finding: dict[str, Any],
    path: str | Path,
    *,
    phase: str | None = None,
    ordinal: int | None = None,
) -> FindingNormalizationResult:
    """Normalize a finding before writing it to disk."""
    path = Path(path)
    raw = dict(finding)

    normalized_phase = canonicalize_finding_phase(phase)
    if normalized_phase and "phase" not in raw:
        raw["phase"] = normalized_phase

    next_ordinal = ordinal or (_count_existing_records(path) + 1)
    return normalize_finding(raw, source_path=path, ordinal=next_ordinal)


def append_finding_record(
    path: str | Path,
    finding: dict[str, Any],
    *,
    phase: str | None = None,
) -> FindingNormalizationResult:
    """Append a normalized finding record to a JSONL file."""
    path = Path(path)
    result = prepare_finding_for_write(finding, path, phase=phase)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as handle:
        handle.write(json.dumps(result.finding, default=str) + "\n")
    return result


def write_findings_records(
    path: str | Path,
    findings: Iterable[dict[str, Any]],
    *,
    phase: str | None = None,
) -> list[FindingNormalizationResult]:
    """Write normalized findings records to a JSONL file, replacing contents."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    results: list[FindingNormalizationResult] = []
    with path.open("w") as handle:
        for ordinal, finding in enumerate(findings, start=1):
            result = prepare_finding_for_write(
                finding,
                path,
                phase=phase,
                ordinal=ordinal,
            )
            handle.write(json.dumps(result.finding, default=str) + "\n")
            results.append(result)
    return results


__all__ = [
    "append_finding_record",
    "prepare_finding_for_write",
    "write_findings_records",
]
