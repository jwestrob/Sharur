"""Protein families: MMseqs2 stage, run-scoped schema, id stability, import, activation, incremental assignment."""

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

from sharur.protein_families import mmseqs
from sharur.protein_families import store
from sharur.storage.duckdb_store import DuckDBStore

REPO_ROOT = Path(__file__).resolve().parents[1]
STAGE09 = REPO_ROOT / "src" / "ingest" / "09_protein_families.py"
DUMMY = Path(os.environ.get("SHARUR_DUMMY_DATASET", REPO_ROOT / "dummy_dataset"))

requires_mmseqs = pytest.mark.skipif(shutil.which("mmseqs") is None, reason="MMseqs2 is not installed")

AA = "ACDEFGHIKLMNPQRSTVWY"
GENOMES = {  # bin_id -> taxonomy
    "g1": "d__Archaea;p__P1;c__C1;o__OrderA;f__F1;g__G1;s__S1",
    "g2": "d__Archaea;p__P1;c__C1;o__OrderA;f__F1;g__G2;s__S2",
    "g3": "d__Archaea;p__P1;c__C1;o__OrderB;f__F2;g__G3;s__S3",
    "g4": "d__Archaea;p__P1;c__C1;o__OrderC;f__F3;g__G4;s__S4",
}
PERMISSIVE = "mmseqs cluster -s 7.5 -e 1e-3 --min-seq-id 0 -c 0.5 --cov-mode 0 --cluster-mode 0"


def _random(rng: random.Random, n: int = 260) -> str:
    return "M" + "".join(rng.choice(AA) for _ in range(n - 1))


def _mutate(rng: random.Random, seq: str, identity: float) -> str:
    """Substitute exactly round((1 - identity) * L) positions (never the first residue)."""
    positions = rng.sample(range(1, len(seq)), round((1 - identity) * len(seq)))
    out = list(seq)
    for i in positions:
        out[i] = rng.choice([a for a in AA if a != seq[i]])
    return "".join(out)


def _proteins() -> list[tuple[str, str, str]]:
    """(protein_id, bin_id, sequence) with known structure.

    Family X: subfamily Xa (base A and 95% variants) and subfamily Xb (B, 42% identical
    to A, with 95% variants). Family Y: one subfamily (C and 90% variants), Pfam-dark.
    Three unrelated singletons. Genome g4 holds the proteins added later: a 98% copy of
    A (joins Xa), a 42% relative of C (new subfamily in Y), a pair of novel homologs
    (one new family) and a novel singleton.
    """
    rng = random.Random(20261112)
    a, c = _random(rng), _random(rng)
    b = _mutate(rng, a, 0.42)
    novel = _random(rng)
    rows = [("xa_g1_1", "g1", a), ("xa_g2_1", "g2", _mutate(rng, a, 0.95)),
            ("xa_g3_1", "g3", _mutate(rng, a, 0.95)), ("xa_g1_2", "g1", _mutate(rng, a, 0.95)),
            ("xb_g2_1", "g2", b), ("xb_g3_1", "g3", _mutate(rng, b, 0.95)),
            ("xb_g1_3", "g1", _mutate(rng, b, 0.95)),
            ("y_g1_1", "g1", c), ("y_g2_1", "g2", _mutate(rng, c, 0.9)), ("y_g3_1", "g3", _mutate(rng, c, 0.9)),
            ("s_g1_1", "g1", _random(rng)), ("s_g2_1", "g2", _random(rng)), ("s_g3_1", "g3", _random(rng) + "*"),
            ("new_g4_xa", "g4", _mutate(rng, a, 0.98)), ("new_g4_y", "g4", _mutate(rng, c, 0.42)),
            ("new_g4_n1", "g4", novel), ("new_g4_n2", "g4", _mutate(rng, novel, 0.9)),
            ("new_g4_m", "g4", _random(rng))]
    return rows


def build_db(path: Path, *, include_g4: bool = True) -> Path:
    rows = [r for r in _proteins() if include_g4 or r[1] != "g4"]
    bins = sorted({r[1] for r in rows})
    with DuckDBStore(path) as db:
        con = db.conn
        con.executemany("INSERT INTO bins(bin_id, taxonomy, n_contigs, total_length) VALUES (?, ?, 1, 100000)",
                        [(b, GENOMES[b]) for b in bins])
        con.executemany("INSERT INTO contigs(contig_id, bin_id, length) VALUES (?, ?, 100000)",
                        [(f"{b}_c1", b) for b in bins])
        con.executemany(
            "INSERT INTO proteins(protein_id, contig_id, bin_id, start, end_coord, strand, gene_index, sequence, "
            "sequence_length) VALUES (?, ?, ?, ?, ?, '+', ?, ?, ?)",
            [(pid, f"{b}_c1", b, 1 + i * 1000, 1 + i * 1000 + 3 * len(seq), i, seq, len(seq.rstrip("*")))
             for i, (pid, b, seq) in enumerate(rows)])
        annotations = []
        for pid in ("xa_g1_1", "xa_g2_1", "xa_g3_1", "xa_g1_2"):
            annotations += [(pid, "pfam", "PF00001", "DomA", 5, 120, 1e-30),
                            (pid, "pfam", "PF00002", "DomB", 130, 250, 1e-20)]
        annotations += [("xb_g2_1", "pfam", "PF00001", "DomA", 5, 120, 1e-25),
                        ("s_g1_1", "pfam", "PF00009", "DomZ", 10, 200, 1e-10)]
        annotations += [(pid, "kofam", "K00001", "K00001", None, None, 1e-50)
                        for pid in ("xa_g1_1", "xa_g2_1", "xa_g3_1")]
        annotations += [("xb_g2_1", "kofam", "K00002", "K00002", None, None, 1e-40)]
        con.executemany(
            "INSERT INTO annotations(annotation_id, protein_id, source, accession, name, start_aa, end_aa, evalue) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
            [(i, *row) for i, row in enumerate(annotations, start=1)])
    return path


@pytest.fixture(scope="module")
def families(tmp_path_factory):
    """Families built on genomes g1-g3; the full database (with g4) for activation and assignment."""
    root = tmp_path_factory.mktemp("families")
    small = build_db(root / "small.duckdb", include_g4=False)
    full = build_db(root / "full.duckdb", include_g4=True)
    result = store.run_pipeline(small, root / "stage09", threads=2, duckdb_threads=2, memory_limit="1GB",
                                split_memory_limit="4G", log=lambda _m: None)
    return {"root": root, "small": small, "full": full, "result": result,
            "staging": Path(result["staging"])}


def _rows(path: Path, sql: str):
    con = duckdb.connect(str(path), read_only=True)
    try:
        return con.execute(sql).fetchall()
    finally:
        con.close()


def _levels(path: Path) -> dict[str, tuple[str, str]]:
    return {pid: (fam, sub) for pid, fam, sub in _rows(
        path, "SELECT protein_id, family_id, subfamily_id FROM protein_family_levels")}


def _write_tsv(path: Path, clusters: list[list[str]]) -> Path:
    """createtsv layout: representative<TAB>member, the representative first in its own cluster."""
    path.write_text("".join(f"{c[0]}\t{m}\n" for c in clusters for m in c))
    return path


def test_levels_require_family_below_subfamily():
    assert [lv.min_seq_id for lv in mmseqs.default_levels(0.5, 0.3)] == [0.5, 0.3]
    with pytest.raises(ValueError):
        mmseqs.default_levels(0.3, 0.5)
    level = mmseqs.default_levels()[1]
    args = level.cluster_args()
    assert args[args.index("--min-seq-id") + 1] == "0.3"
    assert args[args.index("-c") + 1] == "0.8" and args[args.index("--cov-mode") + 1] == "0"
    assert args[args.index("--cluster-reassign") + 1] == "1"
    assert mmseqs._covered(mmseqs.LevelParams(0.3, cov_mode=2), 0.9, 0.4) == 0.9
    assert mmseqs._covered(mmseqs.LevelParams(0.3, cov_mode=0), 0.9, 0.4) == 0.4


def test_mmseqs_module_runs_without_sharur(tmp_path):
    """The cluster-node entry point imports only the standard library."""
    script = REPO_ROOT / "sharur" / "protein_families" / "mmseqs.py"
    for command in ("cluster", "search"):
        out = subprocess.run([sys.executable, "-I", str(script), command, "--help"], cwd=tmp_path,
                             capture_output=True, text=True)
        assert out.returncode == 0, out.stderr


def test_ids_are_documented_and_parse_back():
    assert store.family_id("abc", "P.1_2") == "fam:abc:P.1_2"
    assert store.subfamily_id("abc", "P:odd") == "sub:abc:P:odd"
    assert store.representative_of("sub:abc:P:odd") == "P:odd"
    assert store.build_tag("pfr-0123456789") == "0123456789"


@requires_mmseqs
def test_stage_builds_nested_families_and_stats(families):
    staging, result = families["staging"], families["result"]
    assert result["proteins"] == 13
    tag = store.build_tag(result["build_id"])
    runs = _rows(staging, "SELECT clustering_run, level, parent_run, build_id, origin, is_default, families, "
                          "input_proteins FROM protein_family_runs ORDER BY level")
    assert runs == [(result["family_run"], "family", None, result["build_id"], "sharur", True, result["families"], 13),
                    (result["subfamily_run"], "subfamily", result["family_run"], result["build_id"], "sharur", True,
                     result["subfamilies"], 13)]
    members = _levels(staging)
    assert set(members) == {pid for pid, b, _ in _proteins() if b != "g4"}
    # subfamilies: Xa and Xb separate, Y one; family X holds both X subfamilies
    assert len({members[p][1] for p in ("xa_g1_1", "xa_g2_1", "xa_g3_1", "xa_g1_2")}) == 1
    assert len({members[p][1] for p in ("xb_g2_1", "xb_g3_1", "xb_g1_3")}) == 1
    assert members["xa_g1_1"][1] != members["xb_g2_1"][1]
    assert members["xa_g1_1"][0] == members["xb_g2_1"][0]
    assert members["y_g1_1"][0] != members["xa_g1_1"][0]
    assert len({members[p][0] for p in ("s_g1_1", "s_g2_1", "s_g3_1")}) == 3
    for fam, sub in members.values():
        assert fam.startswith(f"fam:{tag}:") and sub.startswith(f"sub:{tag}:")
        assert members[store.representative_of(fam)][0] == fam
        assert members[store.representative_of(sub)][1] == sub
    reps = _rows(staging, "SELECT clustering_run, COUNT(*) FILTER (WHERE is_representative), "
                          "COUNT(DISTINCT family_id) FROM protein_families GROUP BY 1")
    assert all(n_rep == n_fam for _, n_rep, n_fam in reps)
    assert _rows(staging, "SELECT COUNT(*) FROM (SELECT subfamily_id FROM protein_family_levels "
                          "GROUP BY 1 HAVING COUNT(DISTINCT family_id) > 1)") == [(0,)]
    fx = members["xa_g1_1"][0]
    cols = ["n_members", "n_subfamilies", "n_genomes", "n_orders", "pfam_fraction", "dominant_pfam_architecture",
            "dominant_pfam_architecture_fraction", "ko_fraction", "dominant_ko", "dominant_ko_fraction", "pfam_dark",
            "rep_protein_id", "rep_length", "parent_family_id"]
    stats = dict(zip(cols, _rows(staging, f"SELECT {', '.join(cols)} FROM family_stats "
                                          f"WHERE family_id = '{fx}'")[0]))
    assert stats["n_members"] == 7 and stats["n_subfamilies"] == 2
    assert stats["n_genomes"] == 3 and stats["n_orders"] == 2
    assert stats["pfam_fraction"] == pytest.approx(5 / 7)
    assert stats["dominant_pfam_architecture"] == "DomA - DomB"
    assert stats["dominant_pfam_architecture_fraction"] == pytest.approx(4 / 7)
    assert stats["dominant_ko"] == "K00001" and stats["dominant_ko_fraction"] == pytest.approx(3 / 7)
    assert stats["ko_fraction"] == pytest.approx(4 / 7)
    assert stats["pfam_dark"] is False and stats["parent_family_id"] is None
    assert stats["rep_protein_id"] == store.representative_of(fx) and stats["rep_length"] == 260
    sub_xb = members["xb_g2_1"][1]
    assert _rows(staging, f"SELECT parent_family_id, n_subfamilies FROM family_stats "
                          f"WHERE family_id = '{sub_xb}'") == [(fx, None)]
    fy = members["y_g1_1"][0]
    assert _rows(staging, f"SELECT pfam_dark, dominant_ko, n_orders FROM family_stats "
                          f"WHERE family_id = '{fy}'") == [(True, None, 2)]
    # per-genome presence: (genome_id, feature_id, n_proteins) for the default runs, and per run
    assert _rows(staging, f"SELECT genome_id, n_proteins FROM genome_protein_families WHERE feature_id = '{fx}' "
                          "ORDER BY 1") == [("g1", 3), ("g2", 2), ("g3", 2)]
    for view in ("genome_protein_families", "genome_protein_subfamilies"):
        columns = [r[0] for r in _rows(staging, "SELECT column_name FROM information_schema.columns "
                                                f"WHERE table_name = '{view}' ORDER BY ordinal_position")]
        assert columns == ["genome_id", "feature_id", "n_proteins"]
    assert _rows(staging, "SELECT COUNT(DISTINCT clustering_run) FROM genome_family_presence") == [(2,)]
    params = json.loads(_rows(staging, f"SELECT parameters FROM protein_family_runs "
                                       f"WHERE clustering_run = '{result['family_run']}'")[0][0])
    assert params["min_seq_id"] == 0.3 and params["coverage"] == 0.8 and params["cov_mode"] == 0


@requires_mmseqs
def test_ids_are_stable_across_reruns_and_change_with_parameters(families, tmp_path):
    again = store.run_pipeline(families["small"], tmp_path / "again", threads=2, duckdb_threads=2,
                               memory_limit="1GB", split_memory_limit="4G", log=lambda _m: None)
    assert again["build_id"] == families["result"]["build_id"]
    sql = "SELECT clustering_run, protein_id, family_id, is_representative FROM protein_families ORDER BY 1, 2"
    assert _rows(Path(again["staging"]), sql) == _rows(families["staging"], sql)
    other = store.run_pipeline(families["small"], tmp_path / "other", threads=2, duckdb_threads=2,
                               memory_limit="1GB", split_memory_limit="4G", log=lambda _m: None,
                               levels=mmseqs.default_levels(0.7, 0.3))
    assert other["build_id"] != families["result"]["build_id"]


@requires_mmseqs
def test_build_refuses_changed_sequences(families, tmp_path):
    db = tmp_path / "changed.duckdb"
    shutil.copyfile(families["small"], db)
    con = duckdb.connect(str(db))
    con.execute("UPDATE proteins SET sequence = 'MKV' || sequence WHERE protein_id = 'y_g2_1'")
    con.close()
    with pytest.raises(store.FamilyError, match="differ from the clustered export"):
        store.build_staging(db, families["root"] / "stage09" / "clustering", tmp_path / "x.duckdb",
                            threads=2, memory_limit="1GB", log=lambda _m: None)


@requires_mmseqs
def test_activation_loads_runs_and_checks_identity(families, tmp_path):
    db = tmp_path / "scratch.duckdb"
    shutil.copyfile(families["full"], db)
    con = store.connect(db, threads=2, memory_limit="1GB")
    status = store.load_staging(con, families["staging"], log=lambda _m: None)
    con.close()
    assert status["present"] and status["proteins_with_family"] == 13
    assert status["proteins_without_family"] == 5  # genome g4, added after clustering
    assert status["default_family_run"] == families["result"]["family_run"]
    assert _rows(db, "SELECT COUNT(*) FROM genome_protein_families") == \
        _rows(families["staging"], "SELECT COUNT(*) FROM genome_protein_families")
    from sharur.capabilities import build_capability_brief

    capability = build_capability_brief(db, include_execution=False).get("protein_families")
    assert capability.state.value == "stale" and capability.evidence["proteins_without_family"] == 5
    assert build_capability_brief(families["full"], include_execution=False).get(
        "protein_families").state.value == "unavailable"
    # a protein moved to another genome, or a sequence edited after clustering, is refused
    bad = tmp_path / "edited.duckdb"
    shutil.copyfile(families["full"], bad)
    con = duckdb.connect(str(bad))
    con.execute("UPDATE proteins SET sequence = sequence || 'W' WHERE protein_id = 'xa_g2_1'")
    con.execute("UPDATE proteins SET bin_id = 'g2' WHERE protein_id = 'y_g1_1'")
    con.close()
    con = store.connect(bad, threads=2, memory_limit="1GB")
    with pytest.raises(store.FamilyError, match="another genome"):
        store.load_staging(con, families["staging"], log=lambda _m: None)
    con.execute("UPDATE proteins SET bin_id = 'g1' WHERE protein_id = 'y_g1_1'")
    with pytest.raises(store.FamilyError, match="changed since clustering"):
        store.load_staging(con, families["staging"], log=lambda _m: None)
    assert not store.has_families(con)
    con.close()


@requires_mmseqs
def test_activate_command_stages_and_publishes(families, tmp_path):
    from sharur.cli import app

    db = tmp_path / "dataset" / "sharur.duckdb"
    db.parent.mkdir()
    shutil.copyfile(families["full"], db)
    result = CliRunner().invoke(app, ["families", "activate", "--db", str(db), "--staging",
                                      str(families["staging"]), "--publish"])
    assert result.exit_code == 0, result.output
    assert "Published" in result.output
    assert list(db.parent.glob("sharur.duckdb.pre-protein-families-*"))
    assert _rows(db, "SELECT COUNT(*) FROM protein_families") == [(26,)]
    result = CliRunner().invoke(app, ["families", "status", "--db", str(db)])
    assert json.loads(result.stdout)["proteins_without_family"] == 5


@requires_mmseqs
def test_import_registers_an_external_run_beside_the_build(families, tmp_path):
    staging = tmp_path / "staging.duckdb"
    shutil.copyfile(families["staging"], staging)
    proteins = [pid for pid, b, _ in _proteins() if b != "g4"]
    x = ["xa_g1_1", "xa_g2_1", "xa_g3_1", "xa_g1_2", "xb_g2_1", "xb_g3_1", "xb_g1_3"]
    y = ["y_g1_1", "y_g2_1", "y_g3_1"]
    singles = [[p] for p in proteins if p not in x + y]
    tsv = _write_tsv(tmp_path / "clu.tsv", [x, y, *singles])
    result = store.import_clustering(families["small"], tsv, staging, name="ext-permissive", level="permissive",
                                     parameters=PERMISSIVE, tool_version="18-8cc5c", log=lambda _m: None)
    assert result["families"] == 2 + len(singles) and result["proteins_not_covered"] == 0
    runs = dict(_rows(staging, "SELECT clustering_run, level FROM protein_family_runs"))
    assert runs["ext-permissive"] == "permissive" and len(runs) == 3
    row = _rows(staging, "SELECT origin, is_default, parameters, tool_version FROM protein_family_runs "
                         "WHERE clustering_run = 'ext-permissive'")[0]
    assert row[:2] == ("import", False) and json.loads(row[2]) == PERMISSIVE and row[3] == "18-8cc5c"
    fam = _rows(staging, "SELECT family_id, is_representative FROM protein_families "
                         "WHERE clustering_run = 'ext-permissive' AND protein_id = 'xa_g1_1'")[0]
    assert fam[0].startswith("fam:") and fam[0].endswith(":xa_g1_1") and fam[1] is True
    assert _rows(staging, "SELECT n_members, n_genomes, rep_protein_id FROM family_stats "
                          f"WHERE family_id = '{fam[0]}'") == [(7, 3, "xa_g1_1")]
    # the default views still follow the Sharur build; the presence view lists every run
    assert _rows(staging, "SELECT COUNT(DISTINCT clustering_run) FROM genome_family_presence") == [(3,)]
    assert _rows(staging, "SELECT COUNT(*) FROM genome_protein_families") == \
        _rows(families["staging"], "SELECT COUNT(*) FROM genome_protein_families")
    # a re-import replaces the run; --default takes over its level
    store.import_clustering(families["small"], tsv, staging, name="ext-permissive", level="family",
                            parameters=PERMISSIVE, set_default=True, log=lambda _m: None)
    assert _rows(staging, "SELECT clustering_run FROM protein_family_runs WHERE is_default AND level = 'family'") \
        == [("ext-permissive",)]
    assert _rows(staging, "SELECT COUNT(*) FROM protein_families WHERE clustering_run = 'ext-permissive'") == [(13,)]
    # activation into a database that already holds the build keeps it and adds the imported run
    db = tmp_path / "db.duckdb"
    shutil.copyfile(families["full"], db)
    con = store.connect(db, threads=2, memory_limit="1GB")
    store.load_staging(con, families["staging"], log=lambda _m: None)
    status = store.load_staging(con, staging, runs=["ext-permissive"], log=lambda _m: None)
    con.close()
    assert {r["clustering_run"] for r in status["runs"]} == set(runs)
    assert status["default_family_run"] == "ext-permissive"


@requires_mmseqs
def test_import_rejects_malformed_tables(families, tmp_path):
    staging = tmp_path / "s.duckdb"
    cases = {
        "absent from": "xa_g1_1\txa_g1_1\nxa_g1_1\tnot_a_protein\n",
        "more than one cluster": "xa_g1_1\txa_g1_1\nxa_g1_1\txa_g2_1\ny_g1_1\ty_g1_1\ny_g1_1\txa_g2_1\n",
        "missing from their own cluster": "xa_g1_1\txa_g2_1\n",
    }
    for message, text in cases.items():
        tsv = tmp_path / "bad.tsv"
        tsv.write_text(text)
        with pytest.raises(store.FamilyError, match=message):
            store.import_clustering(families["small"], tsv, staging, name="bad", parameters="x", log=lambda _m: None)
    with pytest.raises(store.FamilyError, match="run name"):
        store.import_clustering(families["small"], tsv, staging, name="pfr-taken", parameters="x")


@requires_mmseqs
def test_import_command(families, tmp_path):
    from sharur.cli import app

    tsv = _write_tsv(tmp_path / "clu.tsv", [[pid] for pid, b, _ in _proteins() if b != "g4"])
    result = CliRunner().invoke(app, ["families", "import", "--db", str(families["small"]), "--tsv", str(tsv),
                                      "--staging", str(tmp_path / "s.duckdb"), "--name", "singletons",
                                      "--parameters", "none"])
    assert result.exit_code == 0, result.output
    assert json.loads(result.stdout[result.stdout.index("{"):])["families"] == 13


@requires_mmseqs
def test_assign_places_new_genomes_without_reclustering(families, tmp_path):
    db = tmp_path / "assign.duckdb"
    shutil.copyfile(families["full"], db)
    con = store.connect(db, threads=2, memory_limit="1GB")
    store.load_staging(con, families["staging"], log=lambda _m: None)
    before = con.execute("SELECT clustering_run, protein_id, family_id FROM protein_families ORDER BY 1, 2").fetchall()
    result = store.assign(con, tmp_path / "work", threads=2, log=lambda _m: None)
    assert result.new_proteins == 5
    levels = {pid: (fam, sub) for pid, fam, sub in con.execute(
        "SELECT protein_id, family_id, subfamily_id FROM protein_family_levels").fetchall()}
    how = dict(con.execute("SELECT protein_id, assignment FROM protein_families WHERE clustering_run = ?",
                           [families["result"]["family_run"]]).fetchall())
    assert levels["new_g4_xa"] == levels["xa_g1_1"] and how["new_g4_xa"] == "search"
    assert levels["new_g4_y"][0] == levels["y_g1_1"][0] and levels["new_g4_y"][1] != levels["y_g1_1"][1]
    assert how["new_g4_y"] == "new_subfamily"
    assert levels["new_g4_n1"][0] == levels["new_g4_n2"][0] and how["new_g4_n1"] == "new_family"
    assert how["new_g4_m"] == "new_family" and levels["new_g4_m"][0] != levels["new_g4_n1"][0]
    tag = store.build_tag(families["result"]["build_id"])
    assert all(fam.startswith(f"fam:{tag}:") and sub.startswith(f"sub:{tag}:") for fam, sub in levels.values())
    # a new subfamily's representative carries the flag; existing rows are untouched
    assert con.execute("SELECT is_representative FROM protein_families WHERE protein_id = 'new_g4_y' "
                       "AND clustering_run = ?", [families["result"]["subfamily_run"]]).fetchone() == (True,)
    assert con.execute("SELECT clustering_run, protein_id, family_id FROM protein_families "
                       "WHERE assignment = 'cluster' ORDER BY 1, 2").fetchall() == before
    assert con.execute(f"SELECT n_members, n_genomes, n_orders FROM family_stats "
                       f"WHERE family_id = '{levels['xa_g1_1'][0]}'").fetchone() == (8, 4, 3)
    assert con.execute("SELECT COUNT(*) FROM (SELECT clustering_run, family_id FROM protein_families GROUP BY ALL "
                       "HAVING COUNT(*) FILTER (WHERE is_representative) <> 1)").fetchone() == (0,)
    batch = con.execute("SELECT batch_id, build_id, input_proteins FROM protein_family_assignments").fetchall()
    assert batch == [(result.batch_id, families["result"]["build_id"], 5)]
    status = store.status(con)
    assert status["proteins_without_family"] == 0
    assert status["assigned_fraction"] == pytest.approx(5 / 18) and status["recluster_recommended"]
    assert store.assign(con, tmp_path / "work2", threads=2, log=lambda _m: None).new_proteins == 0
    con.close()
    # the extended layer, assignment batch included, copies into another database with its digests intact
    copy = tmp_path / "copy.duckdb"
    shutil.copyfile(families["full"], copy)
    target = store.connect(copy, threads=2, memory_limit="1GB")
    store.load_staging(target, db, log=lambda _m: None)
    assert store.status(target)["proteins_without_family"] == 0
    target.close()


@requires_mmseqs
def test_stage09_script_writes_manifest_and_staging(tmp_path):
    db = build_db(tmp_path / "sharur.duckdb", include_g4=True)
    out = tmp_path / "stage09_families"
    env = dict(os.environ, PYTHONPATH=str(REPO_ROOT))
    proc = subprocess.run([sys.executable, str(STAGE09), "--db", str(db), "--out-dir", str(out), "--threads", "2",
                           "--memory-limit", "1GB", "--split-memory-limit", "4G"],
                          capture_output=True, text=True, env=env)
    assert proc.returncode == 0, proc.stderr[-2000:]
    manifest = json.loads((out / "processing_manifest.json").read_text())
    assert manifest["stage"] == "stage09_families" and manifest["proteins"] == 18
    assert (out / "protein_families.duckdb").is_file() and (out / "clustering" / "level1.tsv.gz").is_file()
    assert not (out / "proteins.faa").exists()
    proc = subprocess.run([sys.executable, str(STAGE09), "--db", str(db), "--out-dir", str(out)],
                          capture_output=True, text=True, env=env)
    assert proc.returncode != 0 and "--force" in proc.stderr


def test_ingest_plan_adds_stage09_only_on_request(tmp_path):
    from scripts import ingest

    def plan(**options):
        return ingest.run(input_dir=Path("dummy_dataset"), data_dir=tmp_path / "d", output=tmp_path / "d" / "x.duckdb",
                          mode="tools", dry_run=True, skip_embeddings=True, **options)

    assert not any(c[0].endswith("09_protein_families.py") for c in plan())
    commands = plan(with_families=True, subfamily_min_seq_id=0.6)
    names = [c[0].rsplit("/", 1)[-1] for c in commands]
    assert names.index("09_protein_families.py") > names.index("07_build_knowledge_base.py")
    stage09 = commands[names.index("09_protein_families.py")]
    assert stage09[stage09.index("--subfamily-min-seq-id") + 1] == "0.6"
    assert stage09[stage09.index("--family-min-seq-id") + 1] == "0.3"
    assert stage09[stage09.index("--out-dir") + 1] == str(tmp_path / "d" / "stage09_families")
    from sharur.ingest_cli import app

    result = CliRunner().invoke(app, ["--input-dir", "dummy_dataset", "--data-dir", str(tmp_path / "e"),
                                      "--dry-run", "--with-families", "--subfamily-min-seq-id", "0.2"])
    assert result.exit_code != 0


@requires_mmseqs
@pytest.mark.slow
@pytest.mark.skipif(not DUMMY.is_dir() or shutil.which("prodigal") is None,
                    reason="dummy_dataset or Prodigal unavailable (set SHARUR_DUMMY_DATASET)")
def test_stage09_and_import_on_the_dummy_dataset(tmp_path):
    """Stages 00, 03 and 07, then Stage 09, on the repo's four-genome test dataset; then an external
    clustering with permissive parameters imported as a further run, and everything activated."""
    env = dict(os.environ, PYTHONPATH=str(REPO_ROOT))
    data = tmp_path / "dummy"
    steps = [
        [REPO_ROOT / "src/ingest/00_prepare_inputs.py", "-i", DUMMY, "-o", data / "stage00_prepared", "--copy"],
        [REPO_ROOT / "src/ingest/03_prodigal.py", "-i", data / "stage00_prepared", "-o", data / "stage03_prodigal",
         "--max-workers", "4"],
        [REPO_ROOT / "src/ingest/07_build_knowledge_base.py", "-d", data, "-o", data / "sharur.duckdb", "-t", "4"],
        [STAGE09, "--db", data / "sharur.duckdb", "--out-dir", data / "stage09_families", "--threads", "4",
         "--memory-limit", "2GB", "--split-memory-limit", "4G"],
    ]
    for step in steps:
        proc = subprocess.run([sys.executable, *map(str, step)], capture_output=True, text=True, env=env)
        assert proc.returncode == 0, (step[0], proc.stderr[-2000:])
    db = data / "sharur.duckdb"
    n_proteins = _rows(db, "SELECT COUNT(*) FROM proteins")[0][0]
    staging = data / "stage09_families" / "protein_families.duckdb"
    counts = _rows(staging, "SELECT COUNT(*), COUNT(DISTINCT family_id), COUNT(DISTINCT subfamily_id) "
                            "FROM protein_family_levels")[0]
    assert counts[0] == n_proteins and 0 < counts[1] <= counts[2] <= n_proteins
    # an externally produced clustering of the same proteins, registered without re-clustering
    ext = tmp_path / "ext"
    ext.mkdir()
    con = store.connect(db, read_only=True)
    store.export_proteins(con, ext / "p.faa")
    con.close()
    mm = shutil.which("mmseqs")
    for args in (["createdb", ext / "p.faa", ext / "db"],
                 ["cluster", ext / "db", ext / "clu", ext / "tmp", "-s", "7.5", "-e", "1e-3", "--min-seq-id", "0",
                  "-c", "0.5", "--cov-mode", "0", "--cluster-mode", "0", "--threads", "4",
                  "--split-memory-limit", "4G"],
                 ["createtsv", ext / "db", ext / "db", ext / "clu", ext / "clu.tsv"]):
        subprocess.run([mm, *map(str, args)], check=True, capture_output=True)
    n_ext = len({line.split("\t")[0] for line in (ext / "clu.tsv").read_text().splitlines()})
    imported = store.import_clustering(db, ext / "clu.tsv", staging, name="permissive-e3c50s75",
                                       level="permissive", parameters=PERMISSIVE, log=lambda _m: None)
    assert imported["families"] == n_ext and imported["proteins"] == n_proteins
    con = store.connect(db, threads=4, memory_limit="2GB")
    status = store.load_staging(con, staging, log=lambda _m: None)
    con.close()
    assert status["proteins_without_family"] == 0 and status["families"] == counts[1]
    assert {r["level"] for r in status["runs"]} == {"family", "subfamily", "permissive"}
