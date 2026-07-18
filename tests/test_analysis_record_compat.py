from __future__ import annotations

import importlib.util
import json
import multiprocessing
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pytest

from sharur.core.analysis_record_compat import normalize_finding
from sharur.core.analysis_record_io import (
    FindingValidationError,
    append_finding_record,
    replace_finding_record,
)

REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_module(relative_path: str, module_name: str):
    module_path = REPO_ROOT / relative_path
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    assert spec is not None
    assert spec.loader is not None

    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


MANIFEST_MODULE = _load_module(
    "sharur/operators/manifest.py",
    "test_manifest_module",
)
REPORT_MODULE = _load_module(
    "scripts/compile_comprehensive_report.py",
    "test_compile_comprehensive_report_module",
)


def _append_concurrent_finding(args: tuple[str, int]) -> str:
    path, worker_id = args
    result = append_finding_record(
        path,
        {
            "title": f"Concurrent finding {worker_id}",
            "category": "general",
            "description": "Lock allocation test.",
            "evidence": {},
            "verification": [
                {
                    "claim": f"worker {worker_id}",
                    "query": f"SELECT {worker_id}",
                    "expected": worker_id,
                }
            ],
        },
        phase="exploration",
    )
    return str(result.finding["id"])


def _write_jsonl(path: Path, records: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for record in records:
            handle.write(json.dumps(record) + "\n")


def test_normalize_finding_maps_legacy_fields() -> None:
    result = normalize_finding(
        {
            "id": "E001",
            "title": "Legacy finding",
            "category": "general",
            "description": "Legacy finding description.",
            "evidence": {
                "genes": [
                    {"protein_id": "evidence_gene_1"},
                    {"gene_id": "evidence_gene_2"},
                ],
                "contig_id": "contig_1",
            },
            "genes": ["legacy_gene"],
            "priority": "HIGH",
            "provenance": {
                "query": "SELECT 1",
                "source_report": "exploration/report.md",
            },
        },
        source_path=Path("/tmp/dataset/exploration/findings.jsonl"),
        ordinal=1,
    )

    assert result.finding["id"] == "E001"
    assert result.finding["phase"] == "exploration"
    assert result.finding["schema_version"] == "2.0"
    assert result.finding["protein_ids"] == [
        "legacy_gene",
        "evidence_gene_1",
        "evidence_gene_2",
    ]
    assert result.finding["contigs"] == ["contig_1"]
    assert result.finding["provenance"] == {
        "source_report": "exploration/report.md",
        "query_used": "SELECT 1",
    }
    assert result.is_key_finding is True
    assert "missing verification" in result.issues


def test_manifest_update_findings_aggregates_survey_and_exploration(tmp_path: Path) -> None:
    dataset_dir = tmp_path / "dataset"
    dataset_dir.mkdir()
    (dataset_dir / "sharur.duckdb").write_text("")

    _write_jsonl(
        dataset_dir / "survey" / "findings.jsonl",
        [
            {
                "id": "survey-001",
                "title": "Survey finding",
                "category": "energy_metabolism",
                "description": "Survey description.",
                "evidence": {"protein_id": "survey_gene"},
                "verification": [
                    {
                        "claim": "Survey finding count",
                        "query": "SELECT 1",
                        "expected": 1,
                    }
                ],
                "protein_ids": ["survey_gene"],
            }
        ],
    )
    _write_jsonl(
        dataset_dir / "exploration" / "findings.jsonl",
        [
            {
                "title": "Exploration finding",
                "category": "defense_systems",
                "description": "Exploration description.",
                "evidence": {
                    "genes": [
                        {"protein_id": "explore_gene_1"},
                        {"protein_id": "explore_gene_2"},
                    ]
                },
                "verification": [
                    {
                        "claim": "Exploration finding count",
                        "query": "SELECT 2",
                        "expected": 2,
                    }
                ],
                "priority": "HIGH",
            }
        ],
    )

    manifest = MANIFEST_MODULE.AnalysisManifest(dataset_dir / "sharur.duckdb")
    manifest.update_findings()

    assert manifest.data["findings"]["count"] == 2
    assert manifest.data["findings"]["by_phase"] == {
        "survey": 1,
        "exploration": 1,
    }
    assert manifest.data["exploration"]["coverage"]["proteins_with_findings"] == 3
    assert manifest.data["findings"]["key_findings"] == [
        {
            "id": "exploration-001",
            "title": "Exploration finding",
            "phase": "exploration",
            "category": "defense_systems",
            "novelty": 0,
            "confidence": 0.5,
        }
    ]
    assert (
        manifest.data["findings"]["high_priority"]
        == manifest.data["findings"]["key_findings"]
    )


def test_compile_report_includes_verification_and_normalized_provenance(
    tmp_path: Path,
) -> None:
    dataset_dir = tmp_path / "dataset"
    dataset_dir.mkdir()

    _write_jsonl(
        dataset_dir / "exploration" / "findings.jsonl",
        [
            {
                "title": "Report finding",
                "category": "general",
                "description": "Report description.",
                "evidence": {"protein_id": "gene_1", "contig_id": "contig_1"},
                "verification": [
                    {
                        "claim": "42 proteins support the claim",
                        "query": "SELECT COUNT(DISTINCT protein_id) FROM annotations",
                        "expected": 42,
                    }
                ],
                "provenance": {
                    "query": "SELECT protein_id FROM annotations LIMIT 5",
                },
            }
        ],
    )

    report = REPORT_MODULE.compile_report(dataset_dir)

    assert "### [exploration-001] Report finding" in report
    assert "**Verification:**" in report
    assert "42 proteins support the claim" in report
    assert "SELECT COUNT(DISTINCT protein_id) FROM annotations" in report
    assert "**Provenance:**" in report
    assert "SELECT protein_id FROM annotations LIMIT 5" in report


def test_append_finding_record_emits_canonical_defaults(tmp_path: Path) -> None:
    findings_path = tmp_path / "exploration" / "findings.jsonl"

    result = append_finding_record(
        findings_path,
        {
            "title": "Writer helper finding",
            "category": "general",
            "description": "Helper output should be canonicalized.",
            "evidence": {"protein_id": "protein_1"},
            "verification": [
                {"claim": "one protein", "query": "SELECT 1", "expected": 1}
            ],
            "provenance": {"query": "SELECT protein_id FROM proteins LIMIT 1"},
        },
    )

    assert result.finding["id"] == "exploration-001"
    assert result.finding["schema_version"] == "2.0"
    assert result.finding["phase"] == "exploration"
    assert result.finding["provenance"]["query_used"] == (
        "SELECT protein_id FROM proteins LIMIT 1"
    )
    assert findings_path.read_text().strip()


def test_append_rejects_validation_issues_by_default(tmp_path: Path) -> None:
    findings_path = tmp_path / "survey" / "findings.jsonl"

    with pytest.raises(FindingValidationError, match="missing verification"):
        append_finding_record(
            findings_path,
            {
                "title": "Invalid new finding",
                "description": "No verification supplied.",
                "evidence": {},
            },
        )

    assert not findings_path.exists() or not findings_path.read_text()


def test_append_enforces_lean_canonical_fields(tmp_path: Path) -> None:
    findings_path = tmp_path / "survey" / "findings.jsonl"

    with pytest.raises(FindingValidationError, match="missing category"):
        append_finding_record(
            findings_path,
            {
                "title": "Incomplete finding",
                "description": "Category is required.",
                "evidence": {},
                "verification": [
                    {"claim": "one", "query": "SELECT 1", "expected": 1}
                ],
            },
        )


def test_append_requires_falsification_for_novel_finding(tmp_path: Path) -> None:
    findings_path = tmp_path / "exploration" / "findings.jsonl"

    with pytest.raises(FindingValidationError, match="missing falsification"):
        append_finding_record(
            findings_path,
            {
                "title": "Novel finding",
                "category": "general",
                "description": "Novel claims need a falsification condition.",
                "evidence": {},
                "verification": [
                    {"claim": "one", "query": "SELECT 1", "expected": 1}
                ],
                "novelty": 2,
            },
        )


def test_append_can_write_explicit_draft_spool(tmp_path: Path) -> None:
    findings_path = tmp_path / "drafts" / "findings.jsonl"

    result = append_finding_record(
        findings_path,
        {
            "title": "Draft finding",
            "description": "Drafts may retain validation issues before merge.",
            "evidence": {},
        },
        phase="exploration",
        strict=False,
    )

    assert "missing verification" in result.issues
    assert findings_path.read_text().strip()


def test_append_rejects_duplicate_explicit_id(tmp_path: Path) -> None:
    findings_path = tmp_path / "survey" / "findings.jsonl"
    finding = {
        "id": "survey-001",
        "title": "Finding",
        "category": "general",
        "description": "Duplicate ID test.",
        "evidence": {},
        "verification": [{"claim": "one", "query": "SELECT 1", "expected": 1}],
    }
    append_finding_record(findings_path, finding)

    with pytest.raises(FindingValidationError, match="already exists"):
        append_finding_record(findings_path, finding)


def test_replace_finding_record_preserves_other_lines(tmp_path: Path) -> None:
    findings_path = tmp_path / "exploration" / "findings.jsonl"
    first = {
        "id": "legacy-A001",
        "title": "Legacy record",
        "category": "general",
        "description": "Must remain byte-for-byte stable.",
        "evidence": {},
        "verification": [
            {"claim": "one", "query": "SELECT 1", "expected": 1}
        ],
    }
    second = {
        "id": "exploration-002",
        "title": "Replaceable record",
        "category": "general",
        "description": "Before.",
        "evidence": {},
        "verification": [
            {"claim": "one", "query": "SELECT 1", "expected": 1}
        ],
    }
    append_finding_record(findings_path, first)
    append_finding_record(findings_path, second)
    original_first_line = findings_path.read_text().splitlines(keepends=True)[0]

    second["description"] = "After."
    result = replace_finding_record(
        findings_path,
        "exploration-002",
        second,
    )
    lines = findings_path.read_text().splitlines(keepends=True)

    assert result.finding["description"] == "After."
    assert result.finding["id"] == "exploration-002"
    assert len(lines) == 2
    assert lines[0] == original_first_line


def test_replace_finding_record_rejects_missing_id(tmp_path: Path) -> None:
    findings_path = tmp_path / "exploration" / "findings.jsonl"
    finding = {
        "id": "exploration-001",
        "title": "Existing record",
        "category": "general",
        "description": "Existing.",
        "evidence": {},
        "verification": [
            {"claim": "one", "query": "SELECT 1", "expected": 1}
        ],
    }
    append_finding_record(findings_path, finding)

    with pytest.raises(FindingValidationError, match="does not exist"):
        replace_finding_record(
            findings_path,
            "exploration-999",
            finding,
        )


def test_parallel_appends_allocate_unique_ids(tmp_path: Path) -> None:
    findings_path = tmp_path / "exploration" / "findings.jsonl"
    work = [(str(findings_path), worker_id) for worker_id in range(12)]

    with ProcessPoolExecutor(
        max_workers=4,
        mp_context=multiprocessing.get_context("spawn"),
    ) as executor:
        finding_ids = list(executor.map(_append_concurrent_finding, work))

    records = [
        json.loads(line)
        for line in findings_path.read_text().splitlines()
        if line.strip()
    ]
    assert len(records) == 12
    assert len(set(finding_ids)) == 12
    assert {record["id"] for record in records} == set(finding_ids)
