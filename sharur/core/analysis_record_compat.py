"""
Compatibility helpers for canonical findings records.

This is the only layer that should interpret legacy top-level findings fields
such as `genes`, `priority`, or `provenance.query`.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

from sharur.core.analysis_records import (
    CANONICAL_FINDING_FILE_PHASES,
    CANONICAL_FINDING_SCHEMA_VERSION,
    canonical_finding_id,
    canonicalize_finding_phase,
    infer_finding_phase_from_path,
)


@dataclass(frozen=True)
class FindingNormalizationResult:
    """Normalized finding plus compatibility metadata."""

    finding: dict[str, Any]
    issues: tuple[str, ...] = ()
    is_key_finding: bool = False
    source_path: str | None = None


def _stringify_identifier(value: Any) -> str | None:
    if value is None:
        return None
    if isinstance(value, (str, int, float)):
        text = str(value).strip()
        return text or None
    return None


def _unique_in_order(values: Sequence[str]) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    for value in values:
        if value not in seen:
            seen.add(value)
            ordered.append(value)
    return ordered


def _normalize_string_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, list):
        values = []
        for item in value:
            text = _stringify_identifier(item)
            if text:
                values.append(text)
        return _unique_in_order(values)

    text = _stringify_identifier(value)
    return [text] if text else []


def _iter_protein_refs(value: Any) -> list[str]:
    if not isinstance(value, list):
        return []

    protein_ids: list[str] = []
    for item in value:
        if isinstance(item, Mapping):
            for key in ("protein_id", "gene_id", "id"):
                text = _stringify_identifier(item.get(key))
                if text:
                    protein_ids.append(text)
                    break
        else:
            text = _stringify_identifier(item)
            if text:
                protein_ids.append(text)
    return protein_ids


def _collect_protein_ids(raw: Mapping[str, Any]) -> list[str]:
    protein_ids: list[str] = []

    protein_ids.extend(_normalize_string_list(raw.get("protein_ids")))
    protein_ids.extend(_iter_protein_refs(raw.get("genes")))

    evidence = raw.get("evidence")
    if isinstance(evidence, Mapping):
        protein_id = _stringify_identifier(evidence.get("protein_id"))
        if protein_id:
            protein_ids.append(protein_id)

        protein_ids.extend(_normalize_string_list(evidence.get("protein_ids")))
        protein_ids.extend(_iter_protein_refs(evidence.get("genes")))

    return _unique_in_order(protein_ids)


def _collect_contigs(raw: Mapping[str, Any]) -> list[str]:
    contigs: list[str] = []

    contig_id = _stringify_identifier(raw.get("contig_id"))
    if contig_id:
        contigs.append(contig_id)
    contigs.extend(_normalize_string_list(raw.get("contigs")))

    evidence = raw.get("evidence")
    if isinstance(evidence, Mapping):
        evidence_contig = _stringify_identifier(evidence.get("contig_id"))
        if evidence_contig:
            contigs.append(evidence_contig)
        contigs.extend(_normalize_string_list(evidence.get("contig_ids")))

    return _unique_in_order(contigs)


def _normalize_provenance(value: Any) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        return {}

    provenance = dict(value)
    legacy_query = provenance.pop("query", None)
    if legacy_query is not None and "query_used" not in provenance:
        provenance["query_used"] = legacy_query
    return provenance


def _normalize_verification(value: Any, issues: list[str]) -> list[dict[str, Any]]:
    if value is None:
        issues.append("missing verification")
        return []

    if not isinstance(value, list):
        issues.append("verification must be a list")
        return []

    normalized: list[dict[str, Any]] = []
    for index, item in enumerate(value, start=1):
        if not isinstance(item, Mapping):
            issues.append(f"verification[{index}] must be an object")
            continue

        entry = dict(item)
        missing_fields = [
            field for field in ("claim", "query", "expected") if field not in entry
        ]
        if missing_fields:
            issues.append(
                f"verification[{index}] missing {', '.join(missing_fields)}"
            )
        normalized.append(entry)

    if not normalized:
        issues.append("missing verification")

    return normalized


def _coerce_int(value: Any, default: int) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def _coerce_float(value: Any, default: float) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _derive_key_finding(raw: Mapping[str, Any], novelty: int) -> bool:
    priority = raw.get("priority")
    if isinstance(priority, (int, float)):
        if int(priority) >= 2:
            return True
    else:
        text = _stringify_identifier(priority)
        if text and text.lower() in {"high", "urgent"}:
            return True

    return novelty >= 2


def normalize_finding(
    raw: Mapping[str, Any],
    *,
    source_path: str | Path | None = None,
    ordinal: int | None = None,
) -> FindingNormalizationResult:
    """
    Normalize a finding into the canonical shape used by reports and manifests.

    This function preserves existing IDs on read, including legacy IDs, but new
    synthesized IDs always use the canonical `{phase}-{NNN}` pattern.
    """
    issues: list[str] = []
    raw_dict = dict(raw)

    phase = (
        canonicalize_finding_phase(_stringify_identifier(raw_dict.get("phase")))
        or infer_finding_phase_from_path(source_path)
        or "unknown"
    )

    finding_id = _stringify_identifier(raw_dict.get("id"))
    if not finding_id:
        finding_id = canonical_finding_id(phase, ordinal or 1)

    evidence = raw_dict.get("evidence")
    if evidence is None:
        evidence = {}
    elif isinstance(evidence, Mapping):
        evidence = dict(evidence)
    else:
        issues.append("evidence should be an object")

    novelty = _coerce_int(raw_dict.get("novelty"), default=0)

    normalized = {
        key: value
        for key, value in raw_dict.items()
        if key not in {"genes", "priority"}
    }
    normalized.update(
        {
            "id": finding_id,
            "schema_version": str(
                raw_dict.get("schema_version")
                or CANONICAL_FINDING_SCHEMA_VERSION
            ),
            "phase": phase,
            "title": str(raw_dict.get("title") or "Untitled"),
            "category": str(raw_dict.get("category") or "general"),
            "description": str(raw_dict.get("description") or ""),
            "evidence": evidence,
            "verification": _normalize_verification(
                raw_dict.get("verification"), issues
            ),
            "protein_ids": _collect_protein_ids(raw_dict),
            "contigs": _collect_contigs(raw_dict),
            "provenance": _normalize_provenance(raw_dict.get("provenance")),
            "related_findings": _normalize_string_list(
                raw_dict.get("related_findings")
            ),
            "novelty": novelty,
            "confidence": _coerce_float(raw_dict.get("confidence"), default=0.5),
        }
    )

    if phase == "unknown":
        issues.append("missing phase")

    return FindingNormalizationResult(
        finding=normalized,
        issues=tuple(issues),
        is_key_finding=_derive_key_finding(raw_dict, novelty),
        source_path=str(source_path) if source_path else None,
    )


def load_findings_file(path: str | Path) -> list[FindingNormalizationResult]:
    """Load and normalize findings from a single JSONL file."""
    path = Path(path)
    if not path.exists():
        return []

    findings: list[FindingNormalizationResult] = []
    record_index = 0
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            record_index += 1
            findings.append(
                normalize_finding(
                    json.loads(line),
                    source_path=path,
                    ordinal=record_index,
                )
            )
    return findings


def load_dataset_findings(
    dataset_dir: str | Path,
    *,
    phases: Sequence[str] = CANONICAL_FINDING_FILE_PHASES,
) -> list[FindingNormalizationResult]:
    """Load canonical findings from the standard dataset findings files."""
    dataset_dir = Path(dataset_dir)
    findings: list[FindingNormalizationResult] = []
    for phase in phases:
        findings_path = dataset_dir / phase / "findings.jsonl"
        findings.extend(load_findings_file(findings_path))
    return findings


def build_findings_summary(
    records: Sequence[FindingNormalizationResult],
) -> tuple[dict[str, Any], int]:
    """Build a manifest-friendly summary from normalized findings."""
    findings_info: dict[str, Any] = {
        "count": 0,
        "by_category": {},
        "by_phase": {},
        "key_findings": [],
        "validation_issues": [],
    }
    proteins_with_findings: set[str] = set()

    for record in records:
        finding = record.finding
        findings_info["count"] += 1

        category = finding.get("category") or "general"
        phase = finding.get("phase") or "unknown"
        findings_info["by_category"][category] = (
            findings_info["by_category"].get(category, 0) + 1
        )
        findings_info["by_phase"][phase] = (
            findings_info["by_phase"].get(phase, 0) + 1
        )

        for protein_id in finding.get("protein_ids", []):
            proteins_with_findings.add(str(protein_id))

        if record.is_key_finding:
            findings_info["key_findings"].append(
                {
                    "id": finding.get("id"),
                    "title": finding.get("title"),
                    "phase": phase,
                    "category": category,
                    "novelty": finding.get("novelty", 0),
                    "confidence": finding.get("confidence"),
                }
            )

        if record.issues:
            findings_info["validation_issues"].append(
                {
                    "id": finding.get("id"),
                    "issues": list(record.issues),
                }
            )

    findings_info["key_findings"].sort(
        key=lambda item: (
            -_coerce_int(item.get("novelty"), 0),
            -_coerce_float(item.get("confidence"), 0.0),
            str(item.get("id") or ""),
        )
    )

    return findings_info, len(proteins_with_findings)


__all__ = [
    "FindingNormalizationResult",
    "build_findings_summary",
    "load_dataset_findings",
    "load_findings_file",
    "normalize_finding",
]
