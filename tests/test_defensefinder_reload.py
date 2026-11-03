"""Staged DefenseFinder reload: dry run, publish, idempotence and failure handling."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from typing import TYPE_CHECKING

import duckdb
import pytest

import sharur.colocation as colocation
import sharur.dataset_seal as dataset_seal
from sharur.defensefinder.refresh import ReloadValidationError, reload_defensefinder
from sharur.defensefinder.selection import UNRESOLVED_COMPONENTS_KEY
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
        assert (report.before.annotation_rows, report.after.annotation_rows) == (1, 4)
        assert report.stats.out_of_scope_profile == 1      # the anti-defense hit on p4
        assert report.before.systems_by_type == {"Solo": 1}
        assert report.after.systems_by_type == {"Solo": 1, "Toy": 1}
        assert report.stats.failed_coverage == 1
        assert report.sensitivity["applied"] == 4
        assert report.new_tables == ["annotation_source_versions"]
        text = report.to_text()
        assert "System families newly called (1): Toy (1 systems, 1 genomes)" in text
        assert report.after.systems_by_family == {"Solo": 1, "Toy": 1}
        assert report.after.proteins_by_profile_family == {"Solo": 1, "Toy": 3}
        assert {p.name for p in reports.iterdir()} == {
            "toy_reload_summary.txt", "toy_systems_by_type.tsv", "toy_systems_by_family.tsv",
            "toy_profile_families.tsv", "toy_reload.json",
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
            ("p1", "Toy__ToyA"), ("p2", "ToyB_internal"), ("p3", "Solo__SoloA"), ("p7", "Toy__ToyA"),
        ]
        systems = _query(db_path, "SELECT system_type, protein_ids FROM defense_systems ORDER BY 1")
        assert systems == [("Solo", "p3"), ("Toy", "p1,p2")]
        assert _query(db_path, "SELECT model_version, loaded_rows FROM annotation_source_versions") == [
            (VERSION, 4)]
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


def test_model_families_stay_comparable_across_release_renames(tmp_path):
    from sharur.defensefinder.refresh import align_families, model_families
    from sharur.defensefinder.selection import load_model_set

    old = write_model_set(tmp_path / "v2", version="2.0.2")
    new = write_model_set(tmp_path / "v3", version="3.1.0")
    rm_xml = '<model inter_gene_max_space="5"><gene name="Toy__ToyA" presence="mandatory"/></model>'
    for path in (old / "definitions" / "RM" / "RM" / "RM_Type_I.xml",
                 old / "definitions" / "RM" / "RM" / "RM_Type_IIG_2.xml",
                 old / "definitions" / "RM" / "RM" / "RM_Type_IV_1.xml",
                 new / "definitions" / "RM" / "RM_Type_I" / "N6_RM_Type_I.xml",
                 new / "definitions" / "RM" / "RM_Type_IIG" / "N6_RM_Type_IIG.xml",
                 new / "definitions" / "RM" / "McrBC" / "McrBC.xml",
                 new / "definitions" / "DefenseFinder" / "Other_Type_IV" / "Other_Type_IV.xml",
                 old / "definitions" / "DefenseFinder" / "All_UG" / "UG13a.xml",
                 new / "definitions" / "DefenseFinder" / "DRT-like_UG" / "UG13a.xml",
                 old / "definitions" / "DefenseFinder" / "GasderMIN" / "GasderMIN.xml",
                 new / "definitions" / "DefenseFinder" / "Gasdermin" / "bGSDM.xml"):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(rm_xml)

    before, after = model_families(load_model_set(old)), model_families(load_model_set(new))
    assert before["RM_Type_I"] == after["N6_RM_Type_I"] == "RM_Type_I"
    assert before["RM_Type_IIG_2"] == after["N6_RM_Type_IIG"] == "RM_Type_IIG"
    assert before["RM_Type_IV_1"] == after["McrBC"] == after["Other_Type_IV"] == "RM_Type_IV"
    assert after["Toy"] == "Toy"
    aligned = align_families(before, after)
    assert aligned["UG13a"] == "DRT-like_UG"
    assert aligned["GasderMIN"] == "Gasdermin"
    assert aligned["RM_Type_I"] == "RM_Type_I"


def test_stamp_fingerprint_ignores_columns_added_after_the_snapshot():
    """Regeneration may add predicate_provenance columns (reconcile support); stamps that existed before
    keep their values, so the unchanged-stamps check must not flag the new nullable columns."""
    import duckdb
    from sharur.defensefinder.refresh import _columns, _fingerprint

    conn = duckdb.connect()
    conn.execute("CREATE TABLE predicate_provenance (generation_id BIGINT, scope VARCHAR, pfam_map_sha256 VARCHAR)")
    conn.execute("INSERT INTO predicate_provenance VALUES (1, 'full', 'a'), (2, 'subset', 'a')")
    cols = _columns(conn, "predicate_provenance")
    before = _fingerprint(conn, "predicate_provenance", "WHERE generation_id <= 2", columns=cols)
    conn.execute("ALTER TABLE predicate_provenance ADD COLUMN reconciled_from BIGINT")
    conn.execute("INSERT INTO predicate_provenance VALUES (3, 'subset', 'a', NULL)")
    assert _fingerprint(conn, "predicate_provenance", "WHERE generation_id <= 2", columns=cols) == before
    assert _fingerprint(conn, "predicate_provenance", "WHERE generation_id <= 2") != before  # whole-row hash shifts
    conn.execute("UPDATE predicate_provenance SET pfam_map_sha256 = 'b' WHERE generation_id = 1")
    assert _fingerprint(conn, "predicate_provenance", "WHERE generation_id <= 2", columns=cols) != before


def _hold_open(db_path, *, write: bool):
    """Another process that opens the source database and keeps it open until released.

    ``write`` commits one pfam row first; DuckDB keeps the commit in the WAL, so the
    database file's size and mtime do not change.
    """
    statement = (
        "c.execute(\"INSERT INTO annotations (annotation_id, protein_id, source, accession, name) "
        "VALUES (999, 'p6', 'pfam', 'PF99999', 'concurrent')\")" if write else "pass"
    )
    code = (
        "import duckdb, sys\n"
        f"c = duckdb.connect({str(db_path)!r}, read_only={not write})\n"
        f"{statement}\n"
        "print('ready', flush=True)\n"
        "sys.stdin.readline()\n"
        "c.close()\n"
    )
    process = subprocess.Popen([sys.executable, "-c", code], stdin=subprocess.PIPE,
                               stdout=subprocess.PIPE, text=True)
    assert process.stdout.readline().strip() == "ready"
    return process


def _release(process) -> None:
    process.stdin.write("\n")
    process.stdin.flush()
    process.wait(timeout=60)


class TestConcurrentUse:
    """The source may be opened between the start-of-reload checkpoint and the swap."""

    @pytest.fixture
    def during_integration(self, monkeypatch):
        started = {}
        original = colocation.integrate_defense_results

        def install(db_path, *, write):
            def integrate(*args, **kwargs):
                started["process"] = _hold_open(db_path, write=write)
                return original(*args, **kwargs)
            monkeypatch.setattr(colocation, "integrate_defense_results", integrate)
            return started

        yield install
        if "process" in started and started["process"].poll() is None:
            _release(started["process"])

    def test_write_committed_during_the_reload_blocks_publication(
        self, db_path, hits_path, family, during_integration,
    ):
        started = during_integration(db_path, write=True)
        with pytest.raises(ReloadValidationError, match="changed during the reload"):
            reload_defensefinder(db_path, hits_path, family, dry_run=False)
        _release(started["process"])
        # The concurrent commit stays in the source; nothing was swapped or backed up.
        assert _query(db_path, "SELECT COUNT(*) FROM annotations WHERE accession = 'PF99999'")[0][0] == 1
        assert not list(db_path.parent.glob("*.pre-defensefinder-reload-*"))
        assert not list(db_path.parent.glob("*.staging*"))

    def test_another_open_connection_blocks_publication(
        self, db_path, hits_path, family, during_integration,
    ):
        digest = _sha(db_path)
        started = during_integration(db_path, write=False)
        with pytest.raises(ReloadValidationError, match="open in another process"):
            reload_defensefinder(db_path, hits_path, family, dry_run=False)
        _release(started["process"])
        assert _sha(db_path) == digest
        assert not list(db_path.parent.glob("*.pre-defensefinder-reload-*"))


class TestReportStatus:
    def test_summary_reads_staged_until_publication_succeeds(self, db_path, hits_path, family, tmp_path):
        reports = tmp_path / "reports"

        def explode(staging):
            raise RuntimeError("injected publication failure")

        with pytest.raises(RuntimeError, match="injected"):
            reload_defensefinder(db_path, hits_path, family, dry_run=False, report_dir=reports,
                                 label="toy", _before_publish=explode)
        first = (reports / "toy_reload_summary.txt").read_text().splitlines()[0]
        assert first.startswith("DefenseFinder reload STAGED (not published)")
        report = reload_defensefinder(db_path, hits_path, family, dry_run=False, report_dir=reports,
                                      label="toy")
        assert report.published
        assert (reports / "toy_reload_summary.txt").read_text().startswith("DefenseFinder reload PUBLISHED:")
        assert '"status": "PUBLISHED"' in (reports / "toy_reload.json").read_text()

    def test_failed_reseal_after_the_swap_is_reported(self, db_path, hits_path, family, tmp_path,
                                                      monkeypatch):
        (db_path.parent / "dataset.seal.json").write_text("{}\n")

        def fail(*args, **kwargs):
            raise RuntimeError("seal writer failed")

        monkeypatch.setattr(dataset_seal, "seal_dataset", fail)
        reports = tmp_path / "reports"
        with pytest.raises(RuntimeError, match="seal writer failed"):
            reload_defensefinder(db_path, hits_path, family, dry_run=False, report_dir=reports,
                                 label="toy")
        text = (reports / "toy_reload_summary.txt").read_text()
        assert text.startswith("DefenseFinder reload PUBLISHED, SEAL STALE:")
        assert "Seal update failed after the database was replaced" in text
        assert _query(db_path, "SELECT COUNT(*) FROM defense_systems")[0][0] == 2


def test_reload_reports_and_records_components_beyond_the_solver_budget(
    db_path, hits_path, family, tmp_path, monkeypatch,
):
    monkeypatch.setattr(colocation, "MAX_COMPONENT_CANDIDATES", 0)
    report = reload_defensefinder(db_path, hits_path, family, dry_run=False)
    assert [(c["first_gene_index"], c["models"], c["pass"]) for c in report.unresolved_components] == [
        (1, ["Toy"], "DefenseFinder"), (10, ["Solo"], "DefenseFinder"),
    ]
    assert report.checks["caller_skips_recorded"].startswith("PASS")
    text = report.to_text()
    assert "Conflict components skipped by the conflict-solver budget (2;" in text
    assert "bin_1 / c1 genes 1-2: 1 candidates over 2 proteins" in text
    parameters = json.loads(_query(
        db_path, "SELECT parameters FROM annotation_source_versions WHERE source = 'defensefinder'"
    )[0][0])
    assert [c["first_gene_index"] for c in parameters[UNRESOLVED_COMPONENTS_KEY]] == [1, 10]
    assert parameters["coverage_profile"] == 0.4
