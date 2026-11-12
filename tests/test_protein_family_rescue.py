"""Truncated-gene rescue: a removable layer of derived runs over a protein-family build."""

from __future__ import annotations

import json
import os
import random
import shutil
import subprocess
import sys
from pathlib import Path

import duckdb
import pytest
from typer.testing import CliRunner

from sharur.protein_families import mmseqs, store
from sharur.storage.duckdb_store import DuckDBStore
from tests.test_protein_families import AA, _proteins, _rows, build_db

REPO_ROOT = Path(__file__).resolve().parents[1]
STAGE09 = REPO_ROOT / "src" / "ingest" / "09_protein_families.py"
MMSEQS_SCRIPT = REPO_ROOT / "sharur" / "protein_families" / "mmseqs.py"
requires_mmseqs = pytest.mark.skipif(shutil.which("mmseqs") is None, reason="MMseqs2 is not installed")
QUIET = {"log": lambda _m: None}


def fragment_db(path: Path, *, include_g4: bool = False) -> Path:
    """The family test database plus gene calls cut by contig ends.

    - frag_g3_xa: N-terminal 60% of a family-X protein, partial '10' (rescued into X and subfamily Xa)
    - frag_g2_y: C-terminal 110 aa of a family-Y protein, partial '01' (rescued into Y)
    - novel_g2: an unrelated 100-aa call, partial '01' (stays a singleton)
    - part_g1_xa: residues 100-250 of the same X protein, complete call '00' (no candidate)
    - frag_g4_xb (with g4): C-terminal 70% of an Xb protein, partial '11' (a new family by assignment,
      then rescued). Fragments of related proteins that share an end can instead form a family of
      their own; those are no singletons and stay where clustering or assignment put them.
    """
    db = build_db(path, include_g4=include_g4)
    seq = {pid: s for pid, _, s in _proteins()}
    rng = random.Random(11)
    novel = "M" + "".join(rng.choice(AA) for _ in range(99))
    rows = [("frag_g3_xa", "g3", seq["xa_g1_1"][:156], "10"), ("frag_g2_y", "g2", seq["y_g1_1"][-110:], "01"),
            ("novel_g2", "g2", novel, "01"), ("part_g1_xa", "g1", seq["xa_g1_1"][100:250], "00")]
    if include_g4:
        rows.append(("frag_g4_xb", "g4", seq["xb_g2_1"][-182:], "11"))
    con = duckdb.connect(str(db))
    con.execute("UPDATE proteins SET partial = '00'")
    con.executemany(
        "INSERT INTO proteins(protein_id, contig_id, bin_id, start, end_coord, strand, gene_index, sequence, "
        "sequence_length, partial) VALUES (?, ?, ?, ?, ?, '+', ?, ?, ?, ?)",
        [(pid, f"{b}_c1", b, 200001 + i * 1000, 200001 + i * 1000 + 3 * len(s), 200 + i, s, len(s), flag)
         for i, (pid, b, s, flag) in enumerate(rows)])
    con.close()
    return db


@pytest.fixture(scope="module")
def rescued(tmp_path_factory):
    root = tmp_path_factory.mktemp("rescue")
    db = fragment_db(root / "sharur.duckdb")
    plain = store.run_pipeline(db, root / "plain", threads=2, duckdb_threads=2, memory_limit="1GB",
                               split_memory_limit="4G", **QUIET)
    result = store.run_pipeline(db, root / "rescued", threads=2, duckdb_threads=2, memory_limit="1GB",
                                split_memory_limit="4G", rescue_params=store.RescueParams(),
                                rescue_include=True, **QUIET)
    yield {"root": root, "db": db, "plain": Path(plain["staging"]), "staging": Path(result["staging"]),
           "result": result, "plain_result": plain}
    shutil.rmtree(root, ignore_errors=True)


def _levels(path: Path) -> dict[str, tuple]:
    return {r[0]: r[1:] for r in _rows(path, "SELECT protein_id, family_id, subfamily_id, rescued "
                                             "FROM protein_family_levels")}


@requires_mmseqs
def test_rescue_places_truncated_singletons_in_their_family(rescued):
    staging, result = rescued["staging"], rescued["result"]
    assert result["rescue"]["candidates"] == 3 and result["rescue"]["rescued"] == 2
    levels = _levels(staging)
    assert levels["frag_g3_xa"][:2] == levels["xa_g1_1"][:2] and levels["frag_g3_xa"][2] is True
    assert levels["frag_g2_y"][0] == levels["y_g1_1"][0] and levels["frag_g2_y"][2] is True
    assert levels["novel_g2"][0] == "fam:" + store.build_tag(result["build_id"]) + ":novel_g2"
    assert levels["novel_g2"][2] is False and levels["part_g1_xa"][2] is False
    evidence = {r[0]: r[1:] for r in _rows(staging, """
        SELECT protein_id, family_id, joined_subfamily, identity, query_coverage, partial, original_family_id,
               hit_protein_id FROM protein_family_rescues""")}
    assert set(evidence) == {"frag_g3_xa", "frag_g2_y"}
    fam, joined, identity, qcov, partial, original, hit = evidence["frag_g3_xa"]
    assert fam == levels["xa_g1_1"][0] and joined is True and identity == pytest.approx(1.0)
    assert qcov == pytest.approx(1.0) and partial == "10" and original.endswith(":frag_g3_xa")
    assert levels[hit][0] == fam
    # the clustering partition and build id are untouched: base runs equal a build without the layer
    sql = ("SELECT clustering_run, protein_id, family_id, is_representative, assignment FROM protein_families "
           "WHERE clustering_run NOT LIKE '%-rescued' ORDER BY 1, 2")
    assert _rows(staging, sql) == _rows(rescued["plain"], sql)
    assert result["build_id"] == rescued["plain_result"]["build_id"]
    runs = dict(_rows(staging, "SELECT clustering_run, origin || ':' || is_default FROM protein_family_runs"))
    assert runs == {f"{result['build_id']}-family": "sharur:false", f"{result['build_id']}-subfamily": "sharur:false",
                    f"{result['build_id']}-family-rescued": "rescue:true",
                    f"{result['build_id']}-subfamily-rescued": "rescue:true"}


@requires_mmseqs
def test_views_and_statistics_include_or_exclude_rescued_members(rescued):
    staging = rescued["staging"]
    fam_x = _levels(staging)["xa_g1_1"][0]
    rescued_run = f"{rescued['result']['build_id']}-family-rescued"
    base_run = f"{rescued['result']['build_id']}-family"
    assert _rows(staging, f"SELECT genome_id, n_proteins, n_rescued FROM genome_protein_families "
                          f"WHERE feature_id = '{fam_x}' ORDER BY 1") == [("g1", 3, 0), ("g2", 2, 0), ("g3", 3, 1)]
    assert _rows(staging, f"SELECT genome_id, n_proteins FROM genome_protein_families_clustered "
                          f"WHERE feature_id = '{fam_x}' ORDER BY 1") == [("g1", 3), ("g2", 2), ("g3", 2)]
    singleton = "fam:" + store.build_tag(rescued["result"]["build_id"]) + ":frag_g3_xa"
    assert _rows(staging, f"SELECT COUNT(*) FROM genome_protein_families WHERE feature_id = '{singleton}'") == [(0,)]
    assert _rows(staging, f"SELECT COUNT(*) FROM genome_protein_families_clustered "
                          f"WHERE feature_id = '{singleton}'") == [(1,)]
    stats = dict(_rows(staging, f"SELECT clustering_run, n_members || '/' || n_rescued FROM family_stats "
                                f"WHERE family_id = '{fam_x}'"))
    assert stats == {base_run: "7/0", rescued_run: "8/1"}
    assert _rows(staging, f"SELECT COUNT(DISTINCT clustering_run) FROM genome_family_presence") == [(4,)]


@requires_mmseqs
def test_parameters_gate_the_rescue(rescued, tmp_path):
    work = tmp_path / "s.duckdb"
    shutil.copyfile(rescued["plain"], work)
    out = store.rescue_staging(rescued["db"], work, params=store.RescueParams(min_family_size=4),
                               set_default=True, mmseqs_threads=2, **QUIET)
    assert out["rescued"] == 1  # family Y has 3 members, family X 7
    assert {r[0] for r in _rows(work, "SELECT protein_id FROM protein_family_rescues")} == {"frag_g3_xa"}
    out = store.rescue_staging(rescued["db"], work, params=store.RescueParams(min_family_size=100),
                               mmseqs_threads=2, **QUIET)
    assert out["rescued"] == 0 and out["candidates"] == 3
    with pytest.raises(ValueError):
        store.RescueParams(min_query_coverage=1.01)
    out = store.rescue_staging(rescued["db"], work, params=store.RescueParams(min_seq_id=0.3), set_default=False,
                               mmseqs_threads=2, **QUIET)
    assert out["rescued"] == 2 and out["included_by_default"] is False
    assert _rows(work, "SELECT clustering_run FROM protein_family_runs WHERE is_default AND level = 'family'") == \
        [(f"{rescued['result']['build_id']}-family",)]
    dropped = store.rescue_staging(rescued["db"], work, drop=True, **QUIET)
    assert sorted(dropped["dropped"]) == sorted(r for r in (f"{rescued['result']['build_id']}-family-rescued",
                                                           f"{rescued['result']['build_id']}-subfamily-rescued"))
    assert _rows(work, "SELECT COUNT(*) FROM protein_family_rescues") == [(0,)]
    assert _rows(work, "SELECT COUNT(*) FROM protein_family_runs") == [(2,)]


@requires_mmseqs
def test_rescue_from_hits_computed_elsewhere_matches(rescued, tmp_path):
    con = store.connect(rescued["plain"], read_only=True)
    con.execute(f"ATTACH '{rescued['db']}' AS __pf_src (READ_ONLY)")
    inputs = store.rescue_inputs(con, tmp_path / "inputs", prefix=store.SRC)
    con.close()
    assert inputs["candidates"] == 3
    proc = subprocess.run([sys.executable, "-I", str(MMSEQS_SCRIPT), "search", "--query",
                           str(tmp_path / "inputs" / "candidates.faa"), "--target",
                           str(tmp_path / "inputs" / "targets.faa"), "--out", str(tmp_path / "hits.tsv"),
                           "--min-seq-id", "0.3", "--coverage", "0.8", "--cov-mode", "2", "--top", "5",
                           "--threads", "2", "--split-memory-limit", "4G"], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr[-1500:]
    work = tmp_path / "s.duckdb"
    shutil.copyfile(rescued["plain"], work)
    from sharur.cli import app

    result = CliRunner().invoke(app, ["families", "rescue", "--db", str(rescued["db"]), "--staging", str(work),
                                      "--hits", str(tmp_path / "hits.tsv"), "--tool-version", "18-8cc5c"])
    assert result.exit_code == 0, result.output
    assert _rows(work, "SELECT COUNT(*) FROM protein_family_runs WHERE is_default AND origin = 'rescue'") == \
        [(2 if store.RESCUE_BY_DEFAULT else 0,)]
    result = CliRunner().invoke(app, ["families", "rescue", "--db", str(rescued["db"]), "--staging", str(work),
                                      "--hits", str(tmp_path / "hits.tsv"), "--tool-version", "18-8cc5c",
                                      "--include"])
    assert result.exit_code == 0, result.output
    sql = "SELECT protein_id, family_id, subfamily_id, hit_protein_id FROM protein_family_rescues ORDER BY 1"
    assert _rows(work, sql) == _rows(rescued["staging"], sql)
    assert not list(tmp_path.glob("s.duckdb.rescue.tmp*"))


@requires_mmseqs
def test_rescue_layer_activates_and_follows_assignment(rescued, tmp_path):
    db = fragment_db(tmp_path / "full.duckdb", include_g4=True)
    con = store.connect(db, threads=2, memory_limit="1GB")
    status = store.load_staging(con, rescued["staging"], **QUIET)
    assert status["rescue"]["rescued"] == 2 and status["rescue"]["included_by_default"]
    assert status["proteins_without_family"] == 6
    result = store.assign(con, tmp_path / "work", threads=2, **QUIET)
    assert result.new_proteins == 6 and result.rescue["rescued"] == 3
    levels = {r[0]: r[1:] for r in con.execute(
        "SELECT protein_id, family_id, rescued FROM protein_family_levels").fetchall()}
    assert levels["frag_g4_xb"] == (levels["xb_g2_1"][0], True)
    base = con.execute("SELECT assignment FROM protein_families WHERE protein_id = 'frag_g4_xb' "
                       "AND clustering_run NOT LIKE '%-rescued' AND clustering_run LIKE '%-family'").fetchone()
    assert base == ("new_family",)
    assert store.status(con)["proteins_without_family"] == 0
    con.close()
    from sharur.capabilities import build_capability_brief

    assert build_capability_brief(db, include_execution=False).get("protein_families").state.value == "available"


@requires_mmseqs
def test_sql_macros_include_or_exclude_rescued(rescued, tmp_path):
    db = tmp_path / "sql.duckdb"
    shutil.copyfile(rescued["db"], db)
    con = store.connect(db, threads=2, memory_limit="1GB")
    store.load_staging(con, rescued["staging"], **QUIET)
    con.close()
    with DuckDBStore(db, read_only=True) as ro:
        fam = ro.conn.execute("SELECT family_id FROM protein_family_levels WHERE protein_id = 'xa_g1_1'").fetchone()[0]
        rows = dict(ro.conn.execute(f"SELECT protein_id, evidence FROM family_proteins('{fam}')").fetchall())
        assert rows["frag_g3_xa"] == "rescued" and len(rows) == 8
        assert "frag_g3_xa" not in dict(ro.conn.execute(
            f"SELECT protein_id, evidence FROM family_proteins('{fam}', include_rescued := false)").fetchall())
        g3 = ro.conn.execute(f"SELECT n_proteins, evidence FROM family_presence('{fam}') WHERE bin_id = 'g3'").fetchone()
        assert g3 == (3, ["observed", "rescued"])
        runs = [r[0] for r in ro.conn.execute("SELECT clustering_run FROM families_of('frag_g3_xa')").fetchall()]
        assert all(r.endswith("-rescued") for r in runs) and len(runs) == 2


@requires_mmseqs
def test_older_staging_files_load(rescued, tmp_path):
    """A staging file written before the rescue layer (no n_rescued column, no rescues table) still loads."""
    old = tmp_path / "old.duckdb"
    shutil.copyfile(rescued["plain"], old)
    con = duckdb.connect(str(old))
    for view in store.VIEWS:
        con.execute(f"DROP VIEW IF EXISTS {view}")
    con.execute("DROP TABLE protein_family_rescues")
    con.execute("ALTER TABLE family_stats DROP COLUMN n_rescued")
    con.close()
    db = tmp_path / "db.duckdb"
    shutil.copyfile(rescued["db"], db)
    con = store.connect(db, threads=2, memory_limit="1GB")
    status = store.load_staging(con, old, **QUIET)
    assert status["proteins_without_family"] == 0 and "rescue" not in status
    assert con.execute("SELECT COUNT(*) FROM family_stats WHERE n_rescued IS NULL").fetchone()[0] > 0
    con.close()
    # the layer goes onto an older file too; its base runs stay as they were (fingerprint check inside)
    out = store.rescue_staging(rescued["db"], old, set_default=False, mmseqs_threads=2, **QUIET)
    assert out["rescued"] == 2 and out["included_by_default"] is False
    assert _rows(old, "SELECT COUNT(*) FROM protein_family_rescues") == [(2,)]


@requires_mmseqs
def test_stage09_rescue_flag(tmp_path):
    db = fragment_db(tmp_path / "sharur.duckdb")
    env = dict(os.environ, PYTHONPATH=str(REPO_ROOT))
    out = tmp_path / "stage09_families"
    proc = subprocess.run([sys.executable, str(STAGE09), "--db", str(db), "--out-dir", str(out), "--threads", "2",
                           "--memory-limit", "1GB", "--split-memory-limit", "4G", "--rescue"],
                          capture_output=True, text=True, env=env)
    assert proc.returncode == 0, proc.stderr[-2000:]
    manifest = json.loads((out / "processing_manifest.json").read_text())
    assert manifest["rescue"]["rescued"] == 2
    assert manifest["rescue"]["included_by_default"] is store.RESCUE_BY_DEFAULT
    from scripts import ingest

    def stage09(**options):
        plan = ingest.run(input_dir=Path("dummy_dataset"), data_dir=tmp_path / "d", output=tmp_path / "d" / "x.duckdb",
                          mode="tools", dry_run=True, skip_embeddings=True, with_families=True, **options)
        return next(c for c in plan if c[0].endswith("09_protein_families.py"))

    assert "--rescue" in stage09(family_rescue=True) and "--no-rescue" in stage09(family_rescue=False)
    assert ("--rescue" if store.RESCUE_BY_DEFAULT else "--no-rescue") in stage09()


@requires_mmseqs
def test_compaction_keeps_rows_and_drops_legacy_indexes(rescued, tmp_path):
    path = tmp_path / "s.duckdb"
    shutil.copyfile(rescued["staging"], path)
    con = duckdb.connect(str(path))
    for view in store.VIEWS:
        con.execute(f"DROP VIEW IF EXISTS {view}")
    con.execute("CREATE INDEX idx_protein_families_protein ON protein_families(protein_id)")
    con.close()
    tables = ("protein_families", "family_stats", "protein_family_runs", "protein_family_rescues")
    before = {t: _rows(path, f"SELECT * FROM {t} ORDER BY ALL") for t in tables}
    out = store.compact_family_file(path, threads=2, memory_limit="1GB")
    assert out["rows"]["protein_families"] == len(before["protein_families"])
    assert {t: _rows(path, f"SELECT * FROM {t} ORDER BY ALL") for t in tables} == before
    assert _rows(path, "SELECT COUNT(*) FROM duckdb_indexes() WHERE index_name LIKE 'idx_protein_families_%'") == \
        [(0,)]
    assert _rows(path, "SELECT COUNT(*) FROM genome_protein_families WHERE n_rescued > 0") == [(2,)]
    assert _rows(path, "SELECT COUNT(*) FROM genome_family_presence WHERE n_rescued > 0") == [(4,)]
