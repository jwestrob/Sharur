"""Tests for manifest normalization and live resume behavior."""

import json

import pytest

from sharur.core.analysis_record_io import append_finding_record
from sharur.operators import Sharur
from sharur.operators.manifest import AnalysisManifest


def test_legacy_manifest_is_normalized_without_crashing(tmp_path):
    db_path = tmp_path / "dpann" / "sharur.duckdb"
    db_path.parent.mkdir()
    (db_path.parent / "manifest.json").write_text(
        json.dumps(
            {
                "dataset": "dpann",
                "description": "legacy metadata",
                "stats": {"total_genomes": 10},
                "sessions": [
                    {
                        "phase": "setup",
                        "note": "download started",
                        "timestamp": "2026-01-01T00:00:00+00:00",
                    }
                ],
            }
        )
    )

    manifest = AnalysisManifest(db_path)
    summary = manifest.get_status_summary()

    assert "Analysis State: dpann" in summary
    assert "legacy manifest shape" in summary
    assert manifest.data["dataset"]["legacy_metadata"]["stats"] == {"total_genomes": 10}
    assert manifest.data["session_log"][0]["action"] == "setup"


def test_resume_refreshes_findings_from_canonical_archive(tmp_path):
    dataset_dir = tmp_path / "dataset"
    dataset_dir.mkdir()
    db_path = dataset_dir / "sharur.duckdb"
    b = Sharur(db_path)
    b.store.execute("SELECT 1")

    append_finding_record(
        dataset_dir / "survey" / "findings.jsonl",
        {
            "title": "Live finding",
            "category": "general",
            "description": "Present in the canonical archive.",
            "evidence": {},
            "verification": [{"claim": "one", "query": "SELECT 1", "expected": 1}],
        },
    )

    summary = b.resume()

    assert "**Status:** in_progress" in summary
    assert "**Findings:** 1" in summary
    assert b.manifest.data["findings"]["count"] == 1


def test_corrupt_manifest_surfaces_warning_without_overwriting_file(tmp_path):
    db_path = tmp_path / "dataset" / "sharur.duckdb"
    db_path.parent.mkdir()
    manifest_path = db_path.parent / "manifest.json"
    manifest_path.write_text("{not json")

    with pytest.warns(RuntimeWarning, match="Could not load manifest"):
        manifest = AnalysisManifest(db_path)

    assert manifest_path.read_text() == "{not json"
    assert "Manifest warnings" in manifest.get_status_summary()


def test_manifest_save_is_atomic_and_cleans_temporary_file(tmp_path):
    db_path = tmp_path / "dataset" / "sharur.duckdb"
    manifest = AnalysisManifest(db_path)

    manifest.save()

    assert json.loads((db_path.parent / "manifest.json").read_text())
    assert list(db_path.parent.glob(".manifest.json.*.tmp")) == []
