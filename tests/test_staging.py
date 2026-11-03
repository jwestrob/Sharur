"""sharur.staging: clone, diff, publish."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import duckdb
import pytest

from sharur.staging import (
    Stage,
    StagingError,
    diff_databases,
    discard,
    publish,
    stage,
    touching_regex,
)


def _make(path: Path) -> Path:
    con = duckdb.connect(str(path))
    con.execute("CREATE TABLE proteins (protein_id VARCHAR PRIMARY KEY, bin_id VARCHAR, sequence VARCHAR)")
    con.execute("INSERT INTO proteins VALUES ('p1','b1','MKV'), ('p2','b1','MAA'), ('p3','b2','MGG')")
    con.execute("CREATE TABLE annotations (protein_id VARCHAR, source VARCHAR, name VARCHAR)")
    con.execute("INSERT INTO annotations VALUES ('p1','pfam','A'), ('p2','pfam','B'), ('p2','cazy','GT2')")
    con.execute("CREATE TABLE notes (k INTEGER)")
    con.execute("INSERT INTO notes VALUES (1), (1)")
    con.close()
    return path


def test_diff_reports_rows_keys_sources_and_schema(tmp_path: Path):
    a = _make(tmp_path / "a.duckdb")
    b = tmp_path / "b.duckdb"
    b.write_bytes(a.read_bytes())
    con = duckdb.connect(str(b))
    con.execute("UPDATE proteins SET bin_id = 'b9' WHERE protein_id = 'p3'")
    con.execute("INSERT INTO proteins VALUES ('p4','b2','MCC')")
    con.execute("DELETE FROM proteins WHERE protein_id = 'p1'")
    con.execute("DELETE FROM annotations WHERE source = 'cazy'")
    con.execute("INSERT INTO annotations VALUES ('p3','cazy','GH5'), ('p3','cazy','CBM1')")
    con.execute("ALTER TABLE notes ADD COLUMN extra INTEGER")
    con.execute("CREATE TABLE fresh (x INTEGER)")
    con.close()

    diff = diff_databases(a, b, examples=5)
    by = {t.table: t for t in diff.tables}
    p = by["proteins"]
    assert p.status == "changed" and p.key == ["protein_id"]
    assert (p.keys_added, p.keys_removed, p.keys_changed) == (1, 1, 1)
    assert p.example_keys == {"added": ["p4"], "removed": ["p1"]}
    ann = by["annotations"]
    assert (ann.rows_only_before, ann.rows_only_after) == (1, 2)
    assert ann.changed_sources == ["cazy"]
    assert by["notes"].status == "changed" and by["notes"].columns_added == ["extra"]
    assert by["notes"].rows_only_before is None  # same rows on the common columns
    assert by["fresh"].status == "added"
    regex = touching_regex(diff)
    assert "cazy" in regex and "annotations" not in regex and "proteins" in regex
    text = diff.to_markdown()
    assert "MKV" not in text and "MCC" not in text  # sequence content stays out


def test_diff_identical(tmp_path: Path):
    a = _make(tmp_path / "a.duckdb")
    b = tmp_path / "b.duckdb"
    b.write_bytes(a.read_bytes())
    diff = diff_databases(a, b)
    assert diff.changed() == [] and touching_regex(diff) is None
    assert "identical" in diff.to_markdown()


def test_stage_publish_swaps_keeps_backup_and_reseals(tmp_path: Path, monkeypatch):
    (tmp_path / "ds").mkdir()
    db = _make(tmp_path / "ds" / "sharur.duckdb")
    (tmp_path / "ds" / "dataset.seal.json").write_text(json.dumps({"policy": {"hash_large_files": False}}))
    sealed = []
    import sharur.dataset_seal as ds_mod

    monkeypatch.setattr(ds_mod, "seal_dataset", lambda db, **kw: sealed.append((Path(db), kw)))
    st = stage(db, label="cazy-test")
    assert Path(st.staged).is_file() and st.manifest.is_file()
    con = duckdb.connect(st.staged)
    con.execute("DELETE FROM annotations WHERE source = 'cazy'")
    con.close()
    # the canonical file is untouched until publish
    assert duckdb.connect(str(db), read_only=True).execute(
        "SELECT COUNT(*) FROM annotations").fetchone()[0] == 3
    result = publish(Stage.load(st.staged))
    assert duckdb.connect(str(db), read_only=True).execute(
        "SELECT COUNT(*) FROM annotations").fetchone()[0] == 2
    assert Path(result.backup).is_file() and ".pre-cazy-test-" in result.backup
    assert duckdb.connect(result.backup, read_only=True).execute(
        "SELECT COUNT(*) FROM annotations").fetchone()[0] == 3
    assert result.resealed and sealed[0][1]["hash_large_files"] is False
    assert not Path(st.staged).exists() and not st.manifest.exists()


def test_publish_refuses_when_source_changed_after_staging(tmp_path: Path):
    db = _make(tmp_path / "sharur.duckdb")
    st = stage(db, label="t")
    con = duckdb.connect(str(db))
    con.execute("INSERT INTO notes VALUES (2)")
    con.execute("CHECKPOINT")
    con.close()
    with pytest.raises(StagingError, match="changed after it was staged"):
        publish(st)
    assert duckdb.connect(str(db), read_only=True).execute("SELECT COUNT(*) FROM notes").fetchone()[0] == 3
    assert Path(st.staged).exists()  # left for inspection
    assert not list(tmp_path.glob("sharur.duckdb.pre-*"))
    discard(st)
    assert not Path(st.staged).exists()


def test_publish_refuses_while_another_process_holds_the_source(tmp_path: Path):
    db = _make(tmp_path / "sharur.duckdb")
    st = stage(db, label="t")
    holder = subprocess.Popen(
        [sys.executable, "-c",
         "import duckdb, sys, time; c = duckdb.connect(sys.argv[1], read_only=True); "
         "print('open', flush=True); time.sleep(30)", str(db)],
        stdout=subprocess.PIPE, text=True)
    try:
        assert holder.stdout.readline().strip() == "open"
        with pytest.raises(StagingError, match="open in another process"):
            publish(st)
    finally:
        holder.kill()
        holder.wait()
    assert not list(tmp_path.glob("sharur.duckdb.pre-*"))


def test_stage_run_cli_end_to_end(tmp_path: Path):
    root = tmp_path / "repo"
    ds = root / "data" / "ds"
    (ds / "exploration").mkdir(parents=True)
    db = _make(ds / "sharur.duckdb")
    (ds / "exploration" / "findings.jsonl").write_text(json.dumps({"id": "F1", "verification": [
        {"claim": "1 cazy row", "query": "SELECT COUNT(*) FROM annotations WHERE source = 'cazy'",
         "expected": 1},
        {"claim": "untouched", "query": "SELECT COUNT(*) FROM notes", "expected": 2}]}) + "\n" +
        json.dumps({"id": "F2", "verification": [
            {"claim": "pfam", "query": "SELECT COUNT(*) FROM annotations WHERE source = 'pfam'",
             "expected": 2}]}) + "\n")
    writer = ("import duckdb, os; c = duckdb.connect(os.environ['SHARUR_STAGED_DB']); "
              "c.execute(\"INSERT INTO annotations VALUES ('p3','cazy','GH5')\")")
    out = subprocess.run(
        [sys.executable, "-m", "sharur.cli", "stage", "run", "--db", str(db), "--label", "t1",
         "--findings", str(ds), "--report-dir", str(tmp_path / "rep"), "-j", "1", "--",
         sys.executable, "-c", writer],
        capture_output=True, text=True, cwd=root)
    assert out.returncode == 0, out.stderr
    assert "annotations" in out.stdout and "1 triples change" in out.stdout
    assert "F2" not in (tmp_path / "rep" / "t1.reverify.md").read_text()  # pfam finding not replayed
    staged = ds / "sharur.duckdb.t1.staging"
    assert staged.is_file()
    pub = subprocess.run([sys.executable, "-m", "sharur.cli", "stage", "publish", str(staged)],
                         capture_output=True, text=True, cwd=root)
    assert pub.returncode == 0, pub.stderr
    assert duckdb.connect(str(db), read_only=True).execute(
        "SELECT COUNT(*) FROM annotations WHERE source = 'cazy'").fetchone()[0] == 2
