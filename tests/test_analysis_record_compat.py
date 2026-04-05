from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

from sharur.core.analysis_record_compat import normalize_finding
from sharur.core.analysis_record_io import append_finding_record

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
