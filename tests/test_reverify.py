"""sharur reverify: replay findings' verification triples."""

from __future__ import annotations

import contextlib
import json
import os
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

import duckdb
import pytest

from sharur.reverify import (
    WORKER_OVERHEAD_BYTES,
    available_memory,
    classify,
    match,
    memory_bytes,
    judge,
    plan_workers,
    reverify,
)


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


def test_one_row_result_matches_a_recorded_row():
    """SQL returns [[row]]; findings record the row. 2026-11-02: 195 of 198 'FAIL on both' were this."""
    assert match([[908, 714, 2490, 134]], [908, 714, 2490, 134]) == "FAIL"
    assert judge([[908, 714, 2490, 134]], [908, 714, 2490, 134]) == "OK"
    assert judge([[908, 718, 2608, 135]], [908, 714, 2490, 134]) == "FAIL"
    assert judge([[1, 2], [3, 4]], [1, 2]) == "FAIL"  # two rows are a table, never a row
    assert judge([[[1, 2], 5]], [[1, 2], 5]) == "FAIL"  # nested expected compares exactly
    assert judge([[25, 44]], [[25, 44]]) == "OK"


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


def test_substitution_to_a_name_that_extends_the_original(repo: Path):
    """A backup named <db>.pre-... must not be rewritten a second time (attach + helper ATTACH)."""
    live = repo / "data" / "ds" / "sharur.duckdb"
    backup = live.with_name("sharur.duckdb.pre-reload-1")
    _db(backup, 12)
    report = reverify([repo / "data" / "ds"], ids=["F1", "F2", "F3"], workers=2, compare=True,
                      baseline_substitutions=[("data/ds/sharur.duckdb", str(backup))])
    errors = [r for r in report.results if r.baseline.status == "ERROR" or r.replay.status == "ERROR"]
    assert errors == []
    changed = {(r.triple.finding_id, r.triple.index) for r in report.changed()}
    assert changed == {("F1", 1), ("F2", 1), ("F3", 1), ("F3", 2)}


def test_substitution_text_is_idempotent_and_whole_path(tmp_path: Path):
    from sharur.reverify import Substitution

    root = tmp_path
    sub = Substitution([("data/ds/sharur.duckdb", str(root / "data/ds/sharur.duckdb.pre-x"))], root)
    once = sub.text("ATTACH 'data/ds/sharur.duckdb' AS a; -- data/ds/sharur.duckdb.wal stays")
    assert once.count("sharur.duckdb.pre-x") == 1 and "sharur.duckdb.wal" in once
    assert sub.text(once) == once
    assert sub.text(f"'{root}/data/ds/sharur.duckdb'") == f"'{root}/data/ds/sharur.duckdb.pre-x'"


# --------------------------------------------------------------------------- #
# Notes, labelled commands, command output and read-only files
# --------------------------------------------------------------------------- #
def test_interpret_tells_notes_from_statements():
    from sharur.reverify import interpret

    def kind(q, **kw):
        return interpret(q, **kw).kind

    # notes: methods references, prose, tool provenance, commands annotated with their output file
    assert kind("Python: see the methods section for the neighborhood scan") == "narrative"
    assert kind("see reports/methods.md") == "narrative"
    assert kind("blastn -query a.fna -db b") == "narrative"
    assert kind("python a.py -> tables/out.json") == "narrative"  # sh would truncate out.json
    assert kind("python a.py  # then read out.json key n") == "narrative"
    assert kind("python a.py (table out.tsv, column n)") == "narrative"
    assert kind("SELECT 1", check="narrative") == "narrative"
    # statements, including a labelled command that says what it prints
    assert kind("SELECT 1") == "sql" and kind("CALL pragma_version()") == "sql"
    assert kind("grep -c '>' x.faa | head -1") == "shell"  # a quoted '>' is no redirect
    assert kind("python a.py > /dev/null 2>&1 && python a.py --summary") == "shell"
    reading = interpret("Python: cd lane && python facts.py n_units  (prints JSON; key 'n_units')")
    assert (reading.kind, reading.command, reading.select) == ("shell", "cd lane && python facts.py n_units", "n_units")
    assert interpret("Python/shell (not SQL): python f.py (prints JSON)").select is None
    assert kind("shell: echo $(date)") == "shell" and kind("echo $(date)") == "narrative"


def test_read_output_reads_printed_python_values_and_number_lists():
    from sharur.reverify import read_output

    assert read_output('{"a": 1}', {"a": 1}) == {"a": 1}
    assert read_output("{'a': np.float64(0.25), 'b': [1, 2]}", {"a": 0.25}) == {"a": 0.25, "b": [1, 2]}
    assert read_output("12 7", [12, 7]) == [12, 7]
    assert read_output("True", "True") == "True"  # a string expectation keeps the text
    assert read_output("12 genomes", [12]) == "12 genomes"
    assert match(read_output("{'H': 5.5, 'p': 0.01}", {"H": 5.5, "p": 0.01}), {"H": 5.5, "p": 0.01}) == "OK"


def _write_findings(repo: Path, items: list[dict]) -> Path:
    path = repo / "data" / "ds" / "exploration" / "findings.jsonl"
    path.write_text(json.dumps({"id": "N", "title": "notes", "verification": items}) + "\n")
    return path


def test_narrative_triples_are_reported_and_never_run(repo: Path):
    lane = repo / "data" / "ds" / "exploration" / "lane"
    (lane / "facts.py").write_text("import json, sys\nprint(json.dumps({'n_units': 9, 'other': 1}))\n")
    findings = _write_findings(repo, [
        {"claim": "methods", "query": "Python: see the methods section", "expected": 36},
        {"claim": "arrow", "query": f"{sys.executable} data/ds/exploration/lane/facts.py -> data/ds/out.json",
         "expected": 9},
        {"claim": "labelled", "expected": 9,
         "query": f"Python: cd data/ds/exploration/lane && {sys.executable} facts.py  (prints JSON; key 'n_units')"},
        {"claim": "explicit", "query": "SELECT 1", "expected": 1, "check": "narrative"},
    ])
    report = reverify([findings], workers=1)
    status = [r.replay.status for r in report.results]
    assert status == ["NARRATIVE", "NARRATIVE", "OK", "NARRATIVE"]
    assert not (repo / "data" / "ds" / "out.json").exists()
    assert report.counts()["NARRATIVE"] == 3 and len(report.attention()) == 0
    assert "1 OK, 0 FAIL, 0 ERROR, 3 NARRATIVE" in report.to_markdown()
    assert "output file" in report.results[1].replay.error


def test_replays_refuse_to_write_files(repo: Path):
    findings = _write_findings(repo, [
        {"claim": "python write", "expected": 1, "query": "python: open('data/ds/a.txt', 'w').write('x')"},
        {"claim": "child write", "expected": 1,
         "query": f"shell: {sys.executable} -c \"open('data/ds/b.txt', 'w'); print(1)\""},
        {"claim": "redirect", "expected": 1, "query": "shell: echo 1 > data/ds/c.txt"},
        {"claim": "outside", "expected": 1, "query": "python: len(open(__import__('os').devnull).read()) + 1"},
    ])
    report = reverify([findings], workers=1)
    status = [r.replay.status for r in report.results]
    assert status == ["ERROR", "ERROR", "ERROR", "OK"]
    assert "read-only" in report.results[0].replay.error and "read-only" in report.results[2].replay.error
    assert not any((repo / "data" / "ds" / name).exists() for name in ("a.txt", "b.txt", "c.txt"))


@pytest.mark.skipif(sys.platform != "darwin" or not shutil.which("sandbox-exec"), reason="macOS sandbox")
def test_shell_triples_run_in_a_write_denying_sandbox(repo: Path):
    findings = _write_findings(repo, [
        {"claim": "native write", "expected": 1, "query": "shell: touch data/ds/d.txt; echo 1"}])
    report = reverify([findings], workers=1)
    if report.results[0].replay.status == "OK":  # the sandbox applied: touch failed, echo ran
        assert not (repo / "data" / "ds" / "d.txt").exists()


def test_compare_baseline_keeps_dataframe_scans_in_helper_scripts(repo: Path, tmp_path: Path):
    """A helper that queries a pandas frame by name (FROM df) works under a substitution too."""
    lane = repo / "data" / "ds" / "exploration" / "lane"
    (lane / "frame.py").write_text(
        "import duckdb, json\nimport pandas as pd\n"
        "con = duckdb.connect('data/ds/sharur.duckdb', read_only=True)\n"
        "df = pd.DataFrame({'x': [1, 2, 3]})\n"
        "n = con.execute('SELECT COUNT(*) FROM df').fetchone()[0]\n"
        "print(json.dumps(n + con.execute('SELECT COUNT(*) FROM defense_systems').fetchone()[0]))\n")
    findings = _write_findings(repo, [
        {"claim": "frame + systems", "expected": 13,
         "query": f"shell: {sys.executable} data/ds/exploration/lane/frame.py"}])
    report = reverify([findings], workers=1, compare=True,
                      baseline_substitutions=[("data/ds/sharur.duckdb", str(tmp_path / "staged.duckdb"))])
    result = report.results[0]
    assert (result.baseline.status, result.baseline.observed) == ("FAIL", 15), result.baseline.error
    assert (result.replay.status, result.replay.observed) == ("OK", 13)


# --------------------------------------------------------------------------- #
# Pool sizing and process lifetime
# --------------------------------------------------------------------------- #
GiB = 1 << 30


def test_memory_bytes_parses_duckdb_sizes():
    assert memory_bytes("4GB") == 4 * GiB
    assert memory_bytes("512MiB") == 512 << 20
    assert memory_bytes("1.5 G") == int(1.5 * GiB)
    with pytest.raises(ValueError, match="not understood"):
        memory_bytes("lots")


def test_plan_workers_fits_the_pool_to_available_memory(monkeypatch):
    monkeypatch.setattr("os.cpu_count", lambda: 16)
    per_worker = 4 * GiB + WORKER_OVERHEAD_BYTES
    # 2026-11-02: eight 4 GB workers on a loaded 48 GB laptop
    n, note = plan_workers(8, 2, "4GB", 1000, available=12 * GiB)
    assert n == 1
    assert "8 workers lowered to 1" in note
    n, note = plan_workers(None, 2, "4GB", 1000, available=64 * GiB)
    assert n == (32 * GiB) // per_worker
    assert n < 8
    assert note
    assert plan_workers(2, 2, "4GB", 1000, available=64 * GiB) == (2, None)
    assert plan_workers(8, 2, "4GB", 3, available=512 * GiB) == (3, None)
    # one worker always runs, however little memory is free
    assert plan_workers(None, 2, "1GB", 1000, available=1 * GiB)[0] == 1
    # memory unknown: the CPU-based count stands
    monkeypatch.setattr("sharur.reverify.available_memory", lambda: None)
    assert plan_workers(8, 2, "4GB", 1000) == (8, None)


def test_available_memory_respects_a_slurm_allocation(monkeypatch):
    monkeypatch.setenv("SLURM_MEM_PER_NODE", "1024")
    assert available_memory() <= 1 * GiB


def _alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def _wait_gone(pids, seconds: float) -> list[int]:
    deadline = time.time() + seconds
    while time.time() < deadline:
        if not any(_alive(p) for p in pids):
            return []
        time.sleep(0.1)
    return [p for p in pids if _alive(p)]


def _kill_all(pids) -> None:
    for pid in pids:
        with contextlib.suppress(ProcessLookupError, PermissionError):
            os.kill(pid, signal.SIGKILL)


def _sleeper_dataset(tmp_path: Path, pidfile: Path, seconds: int) -> Path:
    """A dataset whose one triple records its shell, worker and child PIDs, then waits."""
    ds = tmp_path / "repo" / "data" / "ds" / "exploration"
    ds.mkdir(parents=True)
    query = (f"shell: sleep {seconds} & echo $$ $PPID $! > {pidfile}.tmp; "
             f"mv {pidfile}.tmp {pidfile}; wait; echo 1")
    (ds / "findings.jsonl").write_text(json.dumps(
        {"id": "S1", "title": "sleeper", "verification": [{"claim": "c", "query": query, "expected": 1}]}) + "\n")
    return ds.parent


def test_shell_timeout_kills_every_process_the_command_started(tmp_path: Path):
    pidfile = tmp_path / "pids"
    ds = _sleeper_dataset(tmp_path, pidfile, 60)
    report = reverify([ds], workers=1, timeout=3)
    shell_pid, _worker_pid, sleep_pid = (int(x) for x in pidfile.read_text().split())
    try:
        assert report.results[0].replay.status == "ERROR"
        assert "timeout" in report.results[0].replay.error
        assert _wait_gone([shell_pid, sleep_pid], 5) == []
    finally:
        _kill_all([shell_pid, sleep_pid])


@pytest.mark.parametrize("sig", ["SIGKILL", "SIGTERM"])
def test_workers_and_their_shell_children_exit_with_the_reverify_process(tmp_path: Path, sig: str):
    """Killing `sharur reverify` must not leave replay workers behind (2026-11-02 hang)."""
    pidfile = tmp_path / "pids"
    ds = _sleeper_dataset(tmp_path, pidfile, 120)
    driver = tmp_path / "driver.py"
    driver.write_text("import sys\nfrom pathlib import Path\nfrom sharur.reverify import reverify\n"
                      "if __name__ == '__main__':\n"
                      "    reverify([Path(sys.argv[1])], workers=1, timeout=300)\n")
    env = dict(os.environ, PYTHONPATH=os.pathsep.join(
        filter(None, [str(Path(__file__).resolve().parents[1]), os.environ.get("PYTHONPATH")])))
    parent = subprocess.Popen([sys.executable, str(driver), str(ds)], env=env,
                              stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    pids: list[int] = []
    try:
        deadline = time.time() + 90
        while not pidfile.exists():
            assert parent.poll() is None, parent.stderr.read().decode()[-2000:]
            assert time.time() < deadline, "replay never started"
            time.sleep(0.1)
        pids = [int(x) for x in pidfile.read_text().split()]
        assert all(_alive(p) for p in pids)
        parent.send_signal(getattr(signal, sig))
        parent.wait(timeout=30)
        assert _wait_gone(pids, 15) == [], "orphaned replay processes"
    finally:
        if parent.poll() is None:
            parent.kill()
            parent.wait()
        _kill_all(pids)
