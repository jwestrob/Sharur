from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

from sharur import Hypothesis as PublicHypothesis
from sharur import LegacyHypothesis
from sharur.core.analysis_record_audit import (
    audit_dataset_findings,
    render_findings_audit_markdown,
)
from sharur.core.models import Hypothesis as ModelsHypothesis
from sharur.core.types import Hypothesis as CanonicalHypothesis
from sharur.reports.findings_pdf import generate_phase_report

REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_script(relative_path: str, module_name: str):
    module_path = REPO_ROOT / relative_path
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


AUDIT_SCRIPT = _load_script(
    "scripts/audit_findings_schema.py",
    "test_audit_findings_schema_script",
)


def _write_jsonl(path: Path, records: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for record in records:
            handle.write(json.dumps(record) + "\n")


def test_audit_dataset_applies_safe_rewrites(tmp_path: Path) -> None:
    dataset_dir = tmp_path / "dataset"
    _write_jsonl(
        dataset_dir / "survey" / "findings.jsonl",
        [
            {
                "title": "Legacy survey finding",
                "category": "general",
                "description": "Old schema finding.",
                "evidence": "legacy evidence string",
                "provenance": {"query": "SELECT 1"},
            }
        ],
    )

    summary = audit_dataset_findings(dataset_dir, apply_changes=True)
    report = render_findings_audit_markdown(summary)

    assert summary["total"] == 1
    assert summary["changed"] == 1
    assert summary["issues"]["missing verification"] == 1
    assert summary["issues"]["evidence should be an object"] == 1
    assert summary["changes"]["set schema_version"] == 1
    assert summary["changes"]["normalized phase"] == 1
    assert "Findings Schema Audit" in report

    migrated = json.loads((dataset_dir / "survey" / "findings.jsonl").read_text().strip())
    assert migrated["schema_version"] == "2.0"
    assert migrated["phase"] == "survey"
    assert migrated["evidence"] == "legacy evidence string"
    assert migrated["provenance"]["query_used"] == "SELECT 1"
    assert "query" not in migrated["provenance"]


def test_audit_dataset_can_convert_legacy_evidence_strings(tmp_path: Path) -> None:
    dataset_dir = tmp_path / "dataset"
    _write_jsonl(
        dataset_dir / "survey" / "findings.jsonl",
        [
            {
                "title": "Legacy survey finding",
                "category": "general",
                "description": "Old schema finding.",
                "evidence": "Rnf subunits: RnfC (997 proteins, 911 genomes); RnfD (871, 799). Complete (>=5): 493 genomes.",
                "verification": [
                    {"claim": "one", "query": "SELECT 1", "expected": 1}
                ],
            }
        ],
    )

    summary = audit_dataset_findings(
        dataset_dir,
        apply_changes=True,
        convert_legacy_evidence=True,
    )

    assert summary["changed"] == 1
    assert (
        summary["changes"]["converted legacy evidence to structured object"] == 1
    )
    assert "evidence should be an object" not in summary["issues"]

    migrated = json.loads((dataset_dir / "survey" / "findings.jsonl").read_text().strip())
    assert migrated["evidence"]["source_format"] == "legacy_free_text"
    assert migrated["evidence"]["legacy_text"].startswith("Rnf subunits:")
    assert migrated["evidence"]["statements"][0]["label"] == "Rnf subunits"
    assert migrated["evidence"]["statements"][0]["text"].startswith("RnfC")


def test_audit_script_main_writes_reports(tmp_path: Path) -> None:
    dataset_dir = tmp_path / "dataset"
    _write_jsonl(
        dataset_dir / "exploration" / "findings.jsonl",
        [
            {
                "title": "Exploration finding",
                "category": "general",
                "description": "Needs verification.",
                "evidence": {"protein_id": "p1"},
                "verification": [
                    {"claim": "one", "query": "SELECT 1", "expected": 1}
                ],
            }
        ],
    )
    report_path = dataset_dir / "reports" / "audit.md"
    json_path = dataset_dir / "reports" / "audit.json"

    old_argv = sys.argv
    try:
        sys.argv = [
            "audit_findings_schema.py",
            "--dataset",
            str(dataset_dir),
            "--convert-legacy-evidence",
            "--report",
            str(report_path),
            "--json",
            str(json_path),
        ]
        assert AUDIT_SCRIPT.main() == 0
    finally:
        sys.argv = old_argv

    assert report_path.exists()
    assert json_path.exists()


def test_generate_phase_report_creates_pdf(tmp_path: Path) -> None:
    dataset_dir = tmp_path / "dataset"
    _write_jsonl(
        dataset_dir / "survey" / "findings.jsonl",
        [
            {
                "title": "Survey report finding",
                "category": "energy_metabolism",
                "description": "Report body.",
                "evidence": {"protein_id": "p1", "count": 1},
                "verification": [
                    {"claim": "one", "query": "SELECT 1", "expected": 1}
                ],
            }
        ],
    )

    output_path = Path(generate_phase_report(dataset_dir, "survey"))
    assert output_path.exists()
    assert output_path.stat().st_size > 0


def test_public_hypothesis_export_is_canonical() -> None:
    assert PublicHypothesis is CanonicalHypothesis
    assert LegacyHypothesis is ModelsHypothesis
