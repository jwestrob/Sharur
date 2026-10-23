"""Staged DefenseFinder reload: dry run, publish, idempotence and failure handling."""

from __future__ import annotations

import hashlib
from typing import TYPE_CHECKING

import duckdb
import pytest

from sharur.defensefinder.refresh import ReloadValidationError, reload_defensefinder
from tests.defensefinder_fixtures import VERSION, build_toy_db, hit, hits_frame, write_model_set


if TYPE_CHECKING:
    from pathlib import Path


def _query(path, sql, params=()):
    conn = duckdb.connect(str(path), read_only=True)
    try:
        return conn.execute(sql, list(params)).fetchall()
    finally:
        conn.close()


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture
def family(tmp_path):
    return write_model_set(tmp_path / "models")


@pytest.fixture
def db_path(tmp_path, family):
    return build_toy_db(tmp_path / "data" / "toy" / "sharur.duckdb", family)


@pytest.fixture
def hits_path(tmp_path):
    path = tmp_path / "hits" / "defensefinder_hits_df.tsv"
    path.parent.mkdir()
    hits_frame(
        hit("p1", "Toy__ToyA", score=40.0, evalue=1e-12),
        hit("p2", "ToyB_internal", score=35.0, evalue=1e-10, hmm_from=1, hmm_to=180),
        hit("p3", "Solo__SoloA", score=60.0),
        hit("p4", "Anti__AntiX", score=80.0),
        hit("p5", "Solo__SoloA", score=60.0, hmm_from=1, hmm_to=20),       # coverage 0.2
        hit("p7", "Toy__ToyA", score=50.0),                                # lone mandatory gene
    ).to_csv(path, sep="\t", index=False)
    return path


class TestReload:
    def test_dry_run_reports_and_leaves_database_unchanged(self, db_path, hits_path, family, tmp_path):
        digest = _sha(db_path)
        work = tmp_path / "work"
        reports = tmp_path / "reports"
        report = reload_defensefinder(db_path, hits_path, family, dry_run=True, work_dir=work,
                                      report_dir=reports, label="toy")
        assert _sha(db_path) == digest
        assert not list(work.glob("*.staging*"))
        assert not list(db_path.parent.glob("*.staging*"))
        assert report.checks
        assert all(v.startswith("PASS") for v in report.checks.values())
        assert report.model_version == VERSION
        assert (report.before.annotation_rows, report.after.annotation_rows) == (1, 5)
        assert report.before.systems_by_type == {"Solo": 1}
        assert report.after.systems_by_type == {"Solo": 1, "Toy": 1}
        assert report.stats.failed_coverage == 1
        assert report.sensitivity["applied"] == 5
        assert report.new_tables == ["annotation_source_versions"]
        text = report.to_text()
        assert "System types newly called (1): Toy (1)" in text
        assert {p.name for p in reports.iterdir()} == {
            "toy_reload_summary.txt", "toy_systems_by_type.tsv", "toy_profiles.tsv", "toy_reload.json",
        }

    def test_publish_replaces_rows_systems_and_semantics(self, db_path, hits_path, family):
        unrelated = _query(db_path, "SELECT * FROM semantic_atoms WHERE protein_id = 'p6' ORDER BY ALL")
        other_rows = _query(db_path, "SELECT * FROM annotations WHERE source NOT LIKE 'defensefinder%' "
                                     "ORDER BY ALL")
        assert _query(db_path, "SELECT COUNT(*) FROM semantic_atoms WHERE protein_id = 'p5' "
                               "AND source_db LIKE 'defensefinder%'")[0][0] > 0

        report = reload_defensefinder(db_path, hits_path, family, dry_run=False)

        assert report.backup is not None
        assert report.backup.exists()
        assert _query(db_path, "SELECT protein_id, accession FROM annotations "
                               "WHERE source = 'defensefinder' ORDER BY 1") == [
            ("p1", "Toy__ToyA"), ("p2", "ToyB_internal"), ("p3", "Solo__SoloA"),
            ("p4", "Anti__AntiX"), ("p7", "Toy__ToyA"),
        ]
        systems = _query(db_path, "SELECT system_type, protein_ids FROM defense_systems ORDER BY 1")
        assert systems == [("Solo", "p3"), ("Toy", "p1,p2")]
        assert _query(db_path, "SELECT model_version, loaded_rows FROM annotation_source_versions") == [
            (VERSION, 5)]
        # The former Solo member lost every DefenseFinder-derived atom; others are untouched.
        assert _query(db_path, "SELECT COUNT(*) FROM semantic_atoms WHERE protein_id = 'p5' "
                               "AND source_db LIKE 'defensefinder%'")[0][0] == 0
        assert _query(db_path, "SELECT COUNT(*) FROM semantic_atoms WHERE protein_id = 'p2' "
                               "AND source_db = 'defensefinder_system'")[0][0] > 0
        assert _query(db_path, "SELECT * FROM semantic_atoms WHERE protein_id = 'p6' "
                               "ORDER BY ALL") == unrelated
        assert _query(db_path, "SELECT * FROM annotations WHERE source NOT LIKE 'defensefinder%' "
                               "ORDER BY ALL") == other_rows

    def test_reload_is_idempotent(self, db_path, hits_path, family):
        reload_defensefinder(db_path, hits_path, family, dry_run=False)
        systems = _query(db_path, "SELECT system_type, protein_ids FROM defense_systems ORDER BY ALL")
        atoms = _query(db_path, "SELECT * FROM semantic_atoms ORDER BY ALL")
        second = reload_defensefinder(db_path, hits_path, family, dry_run=False)
        assert second.before.systems_by_type == second.after.systems_by_type
        assert _query(db_path, "SELECT system_type, protein_ids FROM defense_systems ORDER BY ALL") == systems
        assert _query(db_path, "SELECT * FROM semantic_atoms ORDER BY ALL") == atoms

    def test_failure_before_publication_leaves_database_intact(self, db_path, hits_path, family):
        digest = _sha(db_path)

        def explode(staging):
            assert staging.exists()
            raise RuntimeError("injected publication failure")

        with pytest.raises(RuntimeError, match="injected"):
            reload_defensefinder(db_path, hits_path, family, dry_run=False, _before_publish=explode)
        assert _sha(db_path) == digest
        assert not list(db_path.parent.glob("*.staging*"))
        assert not list(db_path.parent.glob("*.pre-defensefinder-reload-*"))

    def test_hits_for_proteins_outside_the_database_fail_closed(self, db_path, family, tmp_path):
        bad = tmp_path / "bad_hits.tsv"
        hits_frame(hit("p_missing", "Solo__SoloA")).to_csv(bad, sep="\t", index=False)
        digest = _sha(db_path)
        with pytest.raises(ReloadValidationError, match="absent from proteins"):
            reload_defensefinder(db_path, bad, family, dry_run=False)
        assert _sha(db_path) == digest
