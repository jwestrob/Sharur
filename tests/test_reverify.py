"""sharur reverify: replay findings' verification triples."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import duckdb
import pytest

from sharur.reverify import classify, match, reverify


def _db(path: Path, n: int) -> Path:
    con = duckdb.connect(str(path))
    con.execute("CREATE TABLE defense_systems (system_id INTEGER, genome_id VARCHAR)")
    con.execute("INSERT INTO defense_systems SELECT i, 'g' || (i % 3) FROM range(?) t(i)", [n])
    con.execute("CREATE TABLE bins (bin_id VARCHAR)")
    con.execute("INSERT INTO bins VALUES ('g0'), ('g1'), ('g2')")
    con.close()
    return path


@pytest.fixture()
def repo(tmp_path: Path) -> Path:
    root = tmp_path / "repo"
    ds = root / "data" / "ds"
    (ds / "exploration").mkdir(parents=True)
    _db(ds / "sharur.duckdb", 10)
    _db(tmp_path / "staged.duckdb", 12)
    lane = root / "data" / "ds" / "exploration" / "lane"
    lane.mkdir()
    # helper module that opens the dataset database by its own path
    (lane / "helper.py").write_text(
        "import duckdb\n"
        "def n_systems():\n"
        "    con = duckdb.connect('data/ds/sharur.duckdb', read_only=True)\n"
        "    return con.execute('SELECT COUNT(*) FROM defense_systems').fetchone()[0]\n")
    (lane / "count.py").write_text(
        "import duckdb, json\n"
        "con = duckdb.connect('data/ds/sharur.duckdb')\n"
        "print(json.dumps(con.execute('SELECT COUNT(*) FROM defense_systems').fetchone()[0]))\n")
    findings = [
        {"id": "F1", "title": "systems", "verification": [
            {"claim": "10 systems", "query": "SELECT COUNT(*) FROM defense_systems", "expected": 10},
            {"claim": "3 genomes", "query": "SELECT COUNT(*) FROM bins", "expected": 3},
        ]},
        {"id": "F2", "title": "attached", "verification": [
            {"claim": "10 systems (attach)", "query": "SELECT COUNT(*) FROM d.defense_systems",
             "expected": 10, "attach": {"d": "data/ds/sharur.duckdb"}},
            {"claim": "ratio", "query": "python: round(10/3, 3)", "expected": 3.333},
        ]},
        {"id": "F3", "title": "helper", "verification": [
            {"claim": "helper count", "expected": 10,
             "query": "python: (__import__('sys').path.insert(0, 'data/ds/exploration/lane') "
                      "or __import__('helper')).n_systems()"},
            {"claim": "script count", "expected": 10,
             "query": f"shell: {sys.executable} data/ds/exploration/lane/count.py"},
        ]},
        {"id": "F4", "title": "unrelated", "verification": [
            {"claim": "echo", "query": "shell: echo 7", "expected": 7},
            {"claim": "wrong", "query": "SELECT 2", "expected": 3},
        ]},
    ]
    (ds / "exploration" / "findings.jsonl").write_text("".join(json.dumps(f) + "\n" for f in findings))
    reports = root / "data" / "ds" / "reports"
    reports.mkdir()
    (reports / "R.md").write_text("Defense counts (F1, F3).\nUnrelated F4.\n")
    return root


def test_classify_query_forms():
    assert classify("SELECT 1") == "sql"
    assert classify("  WITH x AS (SELECT 1) SELECT * FROM x") == "sql"
    assert classify("shell: wc -l f") == "shell"
    assert classify("python: 1+1") == "python"
    assert classify("python -c \"print(1)\"") == "shell"
    assert classify("Sharur('x', read_only=True).card('p')") == "python"
    assert classify("anything", {"a": "x"}) == "sql"


def test_match_rules():
    assert match(10, 10) == "OK"
    assert match(10.0, 10) == "OK"
    assert match(10.5, 10) == "FAIL"
    assert match(3.33333, 3.333) == "OK"  # rounded to the precision written
    assert match(3.3346, 3.333) == "FAIL"
    assert match([1, 2.004], [1, 2.0]) == "OK"
    assert match({"a": 1}, {"a": 1}) == "OK"
    assert match({"a": 1, "b": 2}, {"a": 1}) == "FAIL"
    assert match(True, 1) == "FAIL"
    assert match("7", 7) == "OK"


def test_replays_every_form(repo: Path):
    report = reverify([repo / "data" / "ds"], workers=2)
    status = {(r.triple.finding_id, r.triple.index): r.replay.status for r in report.results}
    assert status == {("F1", 1): "OK", ("F1", 2): "OK", ("F2", 1): "OK", ("F2", 2): "OK",
                      ("F3", 1): "OK", ("F3", 2): "OK", ("F4", 1): "OK", ("F4", 2): "FAIL"}


def test_touching_selects_whole_findings_including_helper_scripts(repo: Path):
    report = reverify([repo / "data" / "ds"], touching="defense_systems", workers=2)
    # F3 mentions the table only inside its helper scripts; F1's bins triple rides along
    assert sorted({r.triple.finding_id for r in report.results}) == ["F1", "F2", "F3"]
    assert len(report.results) == 6


def test_compare_redirects_sql_attach_helpers_and_shell_children(repo: Path, tmp_path: Path):
    staged = tmp_path / "staged.duckdb"
    report = reverify([repo / "data" / "ds"], touching="defense_systems", workers=2, compare=True,
                      substitutions=[("data/ds/sharur.duckdb", str(staged))],
                      report_globs=["data/ds/reports/*.md"])
    changed = {(r.triple.finding_id, r.triple.index) for r in report.changed()}
    # every triple that reads defense_systems sees the staged 12, through every route
    assert changed == {("F1", 1), ("F2", 1), ("F3", 1), ("F3", 2)}
    for r in report.results:
        if (r.triple.finding_id, r.triple.index) in changed:
            assert r.baseline.observed == 10 and r.replay.observed == 12
            assert r.replay.status == "FAIL"
    assert report.citations["F1"] == ["data/ds/reports/R.md:1"]
    assert report.citations["F3"] == ["data/ds/reports/R.md:1"]
    assert "F4" not in report.citations
    text = report.to_markdown()
    assert "4 triples change" in text and "cited in `data/ds/reports/R.md:1`" in text


def test_baseline_substitution_compares_backup_to_live(repo: Path, tmp_path: Path):
    backup = tmp_path / "staged.duckdb"  # stands in for a pre-publish backup
    report = reverify([repo / "data" / "ds"], ids=["F1"], workers=1, compare=True,
                      baseline_substitutions=[("data/ds/sharur.duckdb", str(backup))])
    by_index = {r.triple.index: r for r in report.results}
    assert by_index[1].baseline.observed == 12 and by_index[1].replay.observed == 10
    assert by_index[1].changed and not by_index[2].changed


def test_replays_open_databases_read_only(repo: Path):
    findings = repo / "data" / "ds" / "exploration" / "findings.jsonl"
    findings.write_text(json.dumps({"id": "W", "verification": [
        {"claim": "write attempt", "expected": 1,
         "query": "python: __import__('duckdb').connect('data/ds/sharur.duckdb')"
                  ".execute('CREATE TABLE t AS SELECT 1 AS x').fetchall()"},
    ]}) + "\n")
    report = reverify([findings], workers=1)
    assert report.results[0].replay.status == "ERROR"
    assert "read-only" in (report.results[0].replay.error or "").lower()
    con = duckdb.connect(str(repo / "data" / "ds" / "sharur.duckdb"), read_only=True)
    assert con.execute("SELECT COUNT(*) FROM information_schema.tables WHERE table_name = 't'").fetchone()[0] == 0


def test_compare_needs_a_substitution(repo: Path):
    with pytest.raises(ValueError):
        reverify([repo / "data" / "ds"], compare=True)
