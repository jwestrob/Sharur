"""
Audit and safe-migration helpers for findings JSONL archives.

This module rewrites only fields that can be normalized without inventing
scientific content. It preserves legacy evidence strings and missing
verification blocks, but reports those as remaining issues.
"""

from __future__ import annotations

import json
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

from sharur.core.analysis_record_compat import normalize_finding

_SENTENCE_SPLIT_RE = re.compile(r"(?<=[.!?])\s+(?=[A-Z(])")


@dataclass(frozen=True)
class FindingAuditRecord:
    """Audit result for a single finding record."""

    source_path: str
    record_index: int
    finding_id: str | None
    migrated: dict[str, Any]
    changes: tuple[str, ...]
    issues: tuple[str, ...]


def _safe_string_list(value: Any) -> list[str]:
    if isinstance(value, list):
        out: list[str] = []
        for item in value:
            if isinstance(item, (str, int, float)) and str(item).strip():
                out.append(str(item).strip())
        return out
    return []


def _migrate_provenance(
    raw_provenance: Any,
    normalized_provenance: Mapping[str, Any],
) -> tuple[Any, list[str]]:
    if not isinstance(raw_provenance, Mapping):
        return raw_provenance, []

    migrated = dict(raw_provenance)
    changes: list[str] = []
    query_used = normalized_provenance.get("query_used")
    if query_used is not None and migrated.get("query_used") != query_used:
        migrated["query_used"] = query_used
        changes.append("normalized provenance.query_used")
    if "query" in migrated:
        migrated.pop("query", None)
        changes.append("removed legacy provenance.query")
    return migrated, changes


def _structure_legacy_evidence(raw_evidence: str) -> dict[str, Any]:
    """Convert legacy free-text evidence into a minimally structured object."""
    text = raw_evidence.strip()
    sentences = [
        sentence.strip()
        for sentence in _SENTENCE_SPLIT_RE.split(text)
        if sentence.strip()
    ]

    statements: list[dict[str, str]] = []
    for sentence in sentences:
        clauses = [clause.strip() for clause in sentence.split(";") if clause.strip()]
        if not clauses:
            clauses = [sentence]

        for clause in clauses:
            entry: dict[str, str] = {}
            if ":" in clause:
                label, body = clause.split(":", 1)
                label = label.strip()
                body = body.strip()
                if label and body and len(label) <= 120:
                    entry["label"] = label
                    entry["text"] = body
                else:
                    entry["text"] = clause
            else:
                entry["text"] = clause
            statements.append(entry)

    evidence = {
        "source_format": "legacy_free_text",
        "legacy_text": text,
        "statements": statements or [{"text": text}],
    }
    return evidence


def safe_migrate_finding(
    raw: Mapping[str, Any],
    *,
    source_path: str | Path,
    ordinal: int,
    convert_legacy_evidence: bool = False,
) -> FindingAuditRecord:
    """Apply only safe structural rewrites to a finding record."""
    raw_for_normalization = dict(raw)
    converted_evidence: dict[str, Any] | None = None
    if convert_legacy_evidence and isinstance(raw_for_normalization.get("evidence"), str):
        converted_evidence = _structure_legacy_evidence(raw_for_normalization["evidence"])
        raw_for_normalization["evidence"] = converted_evidence

    normalized = normalize_finding(
        raw_for_normalization,
        source_path=source_path,
        ordinal=ordinal,
    )
    migrated = dict(raw)
    changes: list[str] = []

    if migrated.get("id") != normalized.finding.get("id") and "id" not in migrated:
        migrated["id"] = normalized.finding["id"]
        changes.append("added canonical id")

    if migrated.get("schema_version") != normalized.finding.get("schema_version"):
        migrated["schema_version"] = normalized.finding["schema_version"]
        changes.append("set schema_version")

    normalized_phase = normalized.finding.get("phase")
    if normalized_phase and migrated.get("phase") != normalized_phase:
        migrated["phase"] = normalized_phase
        changes.append("normalized phase")

    if converted_evidence is not None:
        migrated["evidence"] = converted_evidence
        changes.append("converted legacy evidence to structured object")

    if normalized.finding.get("protein_ids"):
        migrated_protein_ids = _safe_string_list(migrated.get("protein_ids"))
        if migrated_protein_ids != normalized.finding["protein_ids"]:
            migrated["protein_ids"] = list(normalized.finding["protein_ids"])
            changes.append("normalized protein_ids")

    if normalized.finding.get("contigs"):
        migrated_contigs = _safe_string_list(migrated.get("contigs"))
        if migrated_contigs != normalized.finding["contigs"]:
            migrated["contigs"] = list(normalized.finding["contigs"])
            changes.append("normalized contigs")

    raw_provenance = migrated.get("provenance")
    migrated_provenance, provenance_changes = _migrate_provenance(
        raw_provenance,
        normalized.finding.get("provenance", {}),
    )
    if provenance_changes:
        migrated["provenance"] = migrated_provenance
        changes.extend(provenance_changes)

    if "related_findings" in migrated:
        related = normalized.finding.get("related_findings", [])
        if related and migrated.get("related_findings") != related:
            migrated["related_findings"] = related
            changes.append("normalized related_findings")

    return FindingAuditRecord(
        source_path=str(source_path),
        record_index=ordinal,
        finding_id=str(migrated.get("id")) if migrated.get("id") else None,
        migrated=migrated,
        changes=tuple(changes),
        issues=normalized.issues,
    )


def audit_findings_file(
    path: str | Path,
    *,
    apply_changes: bool = False,
    convert_legacy_evidence: bool = False,
) -> dict[str, Any]:
    """Audit a single findings JSONL file."""
    path = Path(path)
    records: list[FindingAuditRecord] = []
    issue_counts: Counter[str] = Counter()
    change_counts: Counter[str] = Counter()

    if not path.exists():
        return {
            "path": str(path),
            "exists": False,
            "total": 0,
            "changed": 0,
            "issues": {},
            "changes": {},
            "records": [],
        }

    with path.open() as handle:
        for ordinal, line in enumerate(handle, start=1):
            line = line.strip()
            if not line:
                continue
            record = safe_migrate_finding(
                json.loads(line),
                source_path=path,
                ordinal=ordinal,
                convert_legacy_evidence=convert_legacy_evidence,
            )
            records.append(record)
            issue_counts.update(record.issues)
            change_counts.update(record.changes)

    if apply_changes and records:
        with path.open("w") as handle:
            for record in records:
                handle.write(json.dumps(record.migrated, default=str) + "\n")

    return {
        "path": str(path),
        "exists": True,
        "total": len(records),
        "changed": sum(1 for record in records if record.changes),
        "issues": dict(issue_counts),
        "changes": dict(change_counts),
        "records": [
            {
                "record_index": record.record_index,
                "id": record.finding_id,
                "changes": list(record.changes),
                "issues": list(record.issues),
            }
            for record in records
        ],
    }


def audit_dataset_findings(
    dataset_dir: str | Path,
    *,
    phases: Sequence[str] = ("survey", "exploration"),
    apply_changes: bool = False,
    convert_legacy_evidence: bool = False,
) -> dict[str, Any]:
    """Audit the canonical findings files for a dataset."""
    dataset_dir = Path(dataset_dir)
    file_summaries = []
    aggregate_issues: Counter[str] = Counter()
    aggregate_changes: Counter[str] = Counter()
    total = 0
    changed = 0

    for phase in phases:
        summary = audit_findings_file(
            dataset_dir / phase / "findings.jsonl",
            apply_changes=apply_changes,
            convert_legacy_evidence=convert_legacy_evidence,
        )
        file_summaries.append(summary)
        total += summary["total"]
        changed += summary["changed"]
        aggregate_issues.update(summary["issues"])
        aggregate_changes.update(summary["changes"])

    return {
        "dataset": dataset_dir.name,
        "dataset_dir": str(dataset_dir),
        "apply_changes": apply_changes,
        "convert_legacy_evidence": convert_legacy_evidence,
        "total": total,
        "changed": changed,
        "issues": dict(aggregate_issues),
        "changes": dict(aggregate_changes),
        "files": file_summaries,
    }


def render_findings_audit_markdown(summary: Mapping[str, Any]) -> str:
    """Render a dataset findings audit summary as Markdown."""
    lines = [
        f"# Findings Schema Audit: {summary['dataset']}",
        "",
        f"- Dataset: `{summary['dataset_dir']}`",
        f"- Findings scanned: {summary['total']}",
        f"- Records changed: {summary['changed']}",
        f"- Apply mode: `{summary['apply_changes']}`",
        f"- Convert legacy evidence: `{summary.get('convert_legacy_evidence', False)}`",
        "",
        "## Aggregate Issues",
        "",
    ]

    if summary["issues"]:
        for issue, count in sorted(
            summary["issues"].items(),
            key=lambda item: (-item[1], item[0]),
        ):
            lines.append(f"- `{issue}`: {count}")
    else:
        lines.append("- None")

    lines.extend(["", "## Safe Rewrites", ""])
    if summary["changes"]:
        for change, count in sorted(
            summary["changes"].items(),
            key=lambda item: (-item[1], item[0]),
        ):
            lines.append(f"- `{change}`: {count}")
    else:
        lines.append("- None")

    for file_summary in summary["files"]:
        lines.extend(
            [
                "",
                f"## {Path(file_summary['path']).relative_to(summary['dataset_dir'])}",
                "",
                f"- Exists: `{file_summary['exists']}`",
                f"- Findings: {file_summary['total']}",
                f"- Changed: {file_summary['changed']}",
            ]
        )
        if file_summary["issues"]:
            lines.append("- Top issues:")
            for issue, count in sorted(
                file_summary["issues"].items(),
                key=lambda item: (-item[1], item[0]),
            )[:5]:
                lines.append(f"  - `{issue}`: {count}")

    return "\n".join(lines)


__all__ = [
    "FindingAuditRecord",
    "audit_dataset_findings",
    "audit_findings_file",
    "render_findings_audit_markdown",
    "safe_migrate_finding",
]
