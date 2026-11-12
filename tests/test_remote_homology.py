"""Remote homology: .hhr parsing, family selection and export, the runner, staging, activation, SQL surface.

Everything runs on synthetic proteins and synthetic profiles. Tests that call HH-suite, MAFFT or
BLAST+ skip when those programs are missing; the others replace them with small stand-in
scripts that write HH-suite-format output.
"""

from __future__ import annotations

import gzip
import json
import os
import random
import re
import shutil
import stat
import subprocess
import sys
from pathlib import Path

import duckdb
import pytest
from typer.testing import CliRunner

from sharur.protein_families import store as pf
from sharur.remote_homology import RESEMBLES_SOURCES, hhr, runner, store
from sharur.storage.duckdb_store import DuckDBStore
from tests.remote_homology_fixtures import PDB_HITS, PFAM_HITS, TITLES, UNIREF_HITS, render_hhr


REPO_ROOT = Path(__file__).resolve().parents[1]
AA = "ACDEFGHIKLMNPQRSTVWY"
GENOMES = {f"g{i}": f"d__Archaea;p__P1;c__C1;o__Order{'AB'[i % 2]};f__F{i};g__G{i};s__S{i}" for i in range(1, 6)}
HH_TOOLS = ("hhsearch", "hhmake", "hhconsensus", "ffindex_build", "mafft")
requires_hhsuite = pytest.mark.skipif(any(shutil.which(t) is None for t in HH_TOOLS),
                                      reason="HH-suite (hhsearch, hhmake, hhconsensus, ffindex_build) or MAFFT is not installed")
requires_blast = pytest.mark.skipif(any(shutil.which(t) is None for t in ("rpsblast", "psiblast", "makeprofiledb")),
                                    reason="BLAST+ (rpsblast, psiblast, makeprofiledb) is not installed")


# --------------------------------------------------------------------------- #
# Synthetic proteins, families and .hhr text
# --------------------------------------------------------------------------- #
def _random(rng: random.Random, n: int) -> str:
    return "M" + "".join(rng.choice(AA) for _ in range(n - 1))


def _mutate(rng: random.Random, seq: str, identity: float) -> str:
    out = list(seq)
    for i in rng.sample(range(1, len(seq)), round((1 - identity) * len(seq))):
        out[i] = rng.choice([a for a in AA if a != seq[i]])
    return "".join(out)


def _dataset():
    """Proteins (id, bin, sequence) and families (representative -> members).

    fa  Pfam-dark, 4 genomes, one duplicated sequence      -> selected
    fb  Pfam-annotated, 3 genomes                          -> not selected (Pfam hit)
    fc  Pfam-dark, 2 genomes                               -> not selected (< 3 genomes)
    fd  Pfam-dark, 3 genomes, 8 unique sequences            -> selected (cap test)
    s1, s2  singletons; s2 also carries the synthetic CDD profile's sequence family
    """
    rng = random.Random(20261112)
    a, b, c, d, cdd = (_random(rng, n) for n in (120, 140, 90, 110, 120))
    rows = [("fa_g1", "g1", a), ("fa_g2", "g2", _mutate(rng, a, 0.85)), ("fa_g3", "g3", _mutate(rng, a, 0.85)),
            ("fa_g4", "g4", _mutate(rng, a, 0.85))]
    rows.append(("fa_g4b", "g4", rows[-1][2]))  # identical to fa_g4
    rows += [("fb_g1", "g1", b), ("fb_g2", "g2", _mutate(rng, b, 0.9)), ("fb_g3", "g3", _mutate(rng, b, 0.9))]
    rows += [("fc_g1", "g1", c), ("fc_g2", "g2", _mutate(rng, c, 0.9))]
    rows += [(f"fd_{i}", f"g{1 + i % 3}", d if i == 0 else _mutate(rng, d, 0.8)) for i in range(8)]
    rows += [("s1", "g5", _random(rng, 100)), ("s2", "g5", _mutate(rng, cdd, 0.8))]
    families = {"fa_g1": ["fa_g1", "fa_g2", "fa_g3", "fa_g4", "fa_g4b"], "fb_g1": ["fb_g1", "fb_g2", "fb_g3"],
                "fc_g1": ["fc_g1", "fc_g2"], "fd_0": [f"fd_{i}" for i in range(8)], "s1": ["s1"], "s2": ["s2"]}
    return rows, families, {"a": a, "cdd": cdd}


def build_db(path: Path) -> Path:
    rows, _, _ = _dataset()
    with DuckDBStore(path) as db:
        con = db.conn
        con.executemany("INSERT INTO bins(bin_id, taxonomy, n_contigs, total_length, completeness, contamination) "
                        "VALUES (?, ?, 1, 100000, 95, 1)", sorted(GENOMES.items()))
        con.executemany("INSERT INTO contigs(contig_id, bin_id, length) VALUES (?, ?, 100000)",
                        [(f"{b}_c1", b) for b in GENOMES])
        con.executemany(
            "INSERT INTO proteins(protein_id, contig_id, bin_id, start, end_coord, strand, gene_index, sequence, "
            "sequence_length) VALUES (?, ?, ?, ?, ?, '+', ?, ?, ?)",
            [(pid, f"{b}_c1", b, 1 + i * 1000, 1 + i * 1000 + 3 * len(seq), i, seq + "*", len(seq) + 1)
             for i, (pid, b, seq) in enumerate(rows)])
        con.executemany(
            "INSERT INTO annotations(annotation_id, protein_id, source, accession, name, evalue) "
            "VALUES (?, ?, 'pfam', 'PF00001', 'DomA', 1e-30)", [(1, "fb_g1"), (2, "fb_g2")])
    return path


def build_families(db: Path, staging: Path, tmp: Path) -> str:
    _, families, _ = _dataset()
    tsv = tmp / "clu.tsv"
    tsv.write_text("".join(f"{rep}\t{m}\n" for rep, members in families.items() for m in members))
    return pf.import_clustering(db, tsv, staging, name="synthetic-permissive", level="permissive",
                                parameters="synthetic", log=lambda _m: None)["tag"]


# --------------------------------------------------------------------------- #
# Parser
# --------------------------------------------------------------------------- #
def test_parse_hhr_reads_header_summary_and_alignment_blocks():
    text = render_hhr("fa_g1", 120, PFAM_HITS, alignments=3)
    result = hhr.parse_hhr_text(text)
    assert (result.query, result.match_columns, result.searched_hmms) == ("fa_g1", 120, 42)
    assert [h.rank for h in result.hits] == [1, 2, 3, 4]
    top = result.hits[0]
    assert top.label == "PF90001.3 SynthA; synthetic family A"
    assert (top.probability, top.aligned_cols, top.q_start, top.q_end, top.t_start, top.t_end, top.t_length) == \
        (99.6, 100, 5, 110, 1, 105, 108)
    assert top.identities == 25.0 and top.evalue == pytest.approx(1e-25)
    # a four-digit template length glued to its parenthesis, and a hit without an alignment block
    assert result.hits[2].t_length == 1162
    assert result.hits[3].identities is None and result.hits[3].label.startswith("PF90005.1 NoiseDom")


def test_summary_names_with_digits_and_truncation():
    hits = [dict(PDB_HITS[2], name="6zzz_y 1-233 mol:protein length:233  LONG NAME THAT IS TRUNCATED")]
    result = hhr.parse_hhr_text(render_hhr("q", 80, hits, alignments=0))
    assert result.hits[0].label.startswith("6zzz_y 1-233") and result.hits[0].t_length == 233
    with pytest.raises(ValueError):
        hhr.parse_hhr_text(render_hhr("q", 80, hits).replace(" 88.0 ", " xx.x "))


def test_split_label_by_database_kind():
    assert hhr.split_label("PF90001.3 SynthA; synthetic family A", "pfam") == ("PF90001", "SynthA",
                                                                              "synthetic family A")
    titles = {"9yyy": "Crystal structure of a protein of unknown function"}
    assert hhr.split_label("9zzz_A mol:protein length:205  SYNTHETIC ENZYME", "pdb") == ("9zzz_A", "9zzz_A",
                                                                                         "SYNTHETIC ENZYME")
    assert hhr.split_label("9yyy_B", "pdb", titles)[2] == titles["9yyy"]
    assert hhr.split_label("9xxx_C 1-233", "pdb", {"9xxx": "T"}) == ("9xxx_C", "9xxx_C", "T")
    acc, _, desc = hhr.split_label(UNIREF_HITS[0]["name"], "uniref")
    assert (acc, desc) == ("UniRef100_SYN000", "Uncharacterized protein")


def test_hit_rows_threshold_cap_and_unknown_function():
    result = hhr.parse_hhr_text(render_hhr("q", 120, PFAM_HITS))
    rows = [dict(zip(hhr.HIT_COLUMNS, r)) for r in hhr.hit_rows(
        result, source=hhr.SOURCE_PFAM, msa_kind=hhr.MSA_FAMILY_ONLY, kind="pfam", min_probability=50)]
    assert [r["hit_accession"] for r in rows] == ["PF90001", "PF90003", "PF90004"]
    assert [r["unknown_function"] for r in rows] == [False, True, False]
    uniref = hhr.parse_hhr_text(render_hhr("q", 120, UNIREF_HITS))
    assert len(list(hhr.hit_rows(uniref, source=hhr.SOURCE_UNIREF, msa_kind=hhr.MSA_ENRICHED, kind="uniref",
                                 min_probability=50, max_hits=20))) == 20
    assert hhr.is_unknown_function("Uncharacterized protein") and not hhr.is_unknown_function("ATPase family")


def test_prelude_lists_the_same_resembles_sources():
    text = (REPO_ROOT / "sharur" / "sql" / "prelude.sql").read_text()
    literal = "(" + ", ".join(f"'{s}'" for s in sorted(RESEMBLES_SOURCES)) + ")"
    assert text.count(f"NOT IN {literal}") >= 2 and text.count(f"source IN {literal}") >= 4
    for match in re.finditer(r"IN \(([^()]*'cdd'[^()]*)\)", text):
        assert {s.strip(" '") for s in match.group(1).split(",")} - {"pfam", "kofam", "kegg"} == RESEMBLES_SOURCES


# --------------------------------------------------------------------------- #
# Selection and export
# --------------------------------------------------------------------------- #
@pytest.fixture(scope="module")
def dataset(tmp_path_factory):
    root = tmp_path_factory.mktemp("remote_homology")
    db = build_db(root / "sharur.duckdb")
    fam = root / "protein_families.duckdb"
    tag = build_families(db, fam, root)
    return {"root": root, "db": db, "families": fam, "tag": tag}


def test_export_selects_pfam_dark_families_in_three_genomes(dataset, tmp_path):
    work = tmp_path / "work"
    manifest = store.export_families(dataset["db"], work, families_db=dataset["families"], max_sequences=5,
                                     log=lambda _m: None)
    fams = runner.read_families(work)
    tag = dataset["tag"]
    assert [f["family_id"] for f in fams] == [f"fam:{tag}:fa_g1", f"fam:{tag}:fd_0"]
    assert manifest["clustering_run"] == "synthetic-permissive" and manifest["families"] == 2
    fa, fd = fams
    assert (fa["n_members"], fa["n_genomes"], fa["n_unique"], fa["n_sequences"]) == (5, 4, 4, 4)
    assert (fd["n_unique"], fd["n_sequences"]) == (8, 5)
    records = runner.read_fasta(work / fa["fasta"])
    assert records[0][0] == "fa_g1" and fa["rep_length"] == len(records[0][1]) == 120
    assert "fa_g4b" not in {r[0] for r in records}  # identical sequence written once, under the lower id
    assert not any(s.endswith("*") for _, s in records)
    again = store.export_families(dataset["db"], tmp_path / "again", families_db=dataset["families"],
                                  max_sequences=5, log=lambda _m: None)
    assert again["export_digest"] == manifest["export_digest"]


def test_export_named_families_and_errors(dataset, tmp_path):
    tag = dataset["tag"]
    m = store.export_families(dataset["db"], tmp_path / "w", families_db=dataset["families"],
                              family_ids=[f"fam:{tag}:fb_g1"], log=lambda _m: None)
    assert m["families"] == 1
    with pytest.raises(store.RemoteHomologyError, match="not in run"):
        store.export_families(dataset["db"], tmp_path / "x", families_db=dataset["families"],
                              family_ids=["fam:nope:x"], log=lambda _m: None)
    with pytest.raises(store.RemoteHomologyError, match="no clustering run"):
        store.export_families(dataset["db"], tmp_path / "y", families_db=dataset["families"], run="other",
                              log=lambda _m: None)
    with pytest.raises(store.RemoteHomologyError, match="no protein-family tables"):
        store.export_families(dataset["db"], tmp_path / "z", log=lambda _m: None)


def test_export_cdd_writes_unique_sequences_by_md5(dataset, tmp_path):
    m = store.export_cdd(dataset["db"], tmp_path, chunk_size=7, log=lambda _m: None)
    rows, _, _ = _dataset()
    assert m["proteins"] == len(rows) and m["unique_sequences"] == len(rows) - 1
    chunks = sorted((tmp_path / "cdd" / "chunks").glob("chunk_*.faa"))
    assert len(chunks) == m["chunks"] == -(-m["unique_sequences"] // 7)
    names = [h for c in chunks for h, _ in runner.read_fasta(c)]
    assert names == sorted(names) and all(re.fullmatch(r"[0-9a-f]{32}", h) for h in names)


# --------------------------------------------------------------------------- #
# Runner with stand-in tools
# --------------------------------------------------------------------------- #
FAKE_TOOL = r'''#!{python}
import json, os, shutil, sys
args = sys.argv[1:]
name = os.path.basename(sys.argv[0])
log = os.environ.get("FAKE_TOOL_LOG")
if {"-h", "--version", "-version"} & set(args):
    print(name + " (stand-in)")
    sys.exit(0)
if log:
    with open(log, "a") as f:
        f.write(name + "\n")
if name == "mafft":
    sys.stdout.write(open(args[-1]).read())
    sys.exit(0)
def opt(flag):
    return args[args.index(flag) + 1] if flag in args else None
inp, out, db = opt("-i"), opt("-o"), opt("-d") or ""
seqs, cur = [], None
for line in open(inp):
    if line.startswith(">"):
        seqs.append("")
    elif seqs:
        seqs[-1] += line.strip()
first = seqs[0]
ncols = sum(1 for c in first if c.isalpha() and c.isupper())
query = open(inp).readline()[1:].split()[0]
sys.path.insert(0, os.environ["FAKE_TOOL_TESTS"])
from remote_homology_fixtures import PDB_HITS, PFAM_HITS, UNIREF_HITS, render_hhr
hits = UNIREF_HITS if name == "hhblits" else PDB_HITS if "pdb" in db else PFAM_HITS
open(out, "w").write(render_hhr(query, ncols, hits))
if name == "hhblits":
    shutil.copyfile(inp, opt("-oa3m"))
'''


@pytest.fixture
def fake_tools(tmp_path, monkeypatch):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    for name in ("hhsearch", "hhblits", "mafft"):
        path = bin_dir / name
        path.write_text(FAKE_TOOL.replace("{python}", sys.executable))
        path.chmod(path.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setenv("FAKE_TOOL_TESTS", str(Path(__file__).parent))
    monkeypatch.setenv("FAKE_TOOL_LOG", str(tmp_path / "calls.log"))
    for part in ("pfam_hhm", "pdb_hhm", "uniref_hhm"):
        (tmp_path / "db").mkdir(exist_ok=True)
        (tmp_path / "db" / f"{part}.ffindex").write_text("")
    (tmp_path / "titles.tsv").write_text(TITLES)
    config = {
        "hhsuite_bin": str(bin_dir), "mafft": str(bin_dir / "mafft"),
        "databases": {
            "pfam": {"prefix": str(tmp_path / "db" / "pfam"), "name": "Pfam", "version": "99.0", "kind": "pfam"},
            "pdb": {"prefix": str(tmp_path / "db" / "pdb"), "name": "PDB100", "version": "2099-01-01",
                    "kind": "pdb", "titles": str(tmp_path / "titles.tsv")},
            "uniref": {"prefix": str(tmp_path / "db" / "uniref"), "name": "UniRef30", "version": "2099_01",
                       "kind": "uniref"}},
        "msa_kinds": ["enriched", "family_only"], "family_targets": ["pfam", "pdb"]}
    path = tmp_path / "config.json"
    path.write_text(json.dumps(config))
    return {"config": path, "calls": tmp_path / "calls.log", "dict": config}


@pytest.fixture
def searched(dataset, fake_tools, tmp_path):
    work = tmp_path / "work"
    store.export_families(dataset["db"], work, families_db=dataset["families"], log=lambda _m: None)
    code = runner.main(["families", "--work", str(work), "--config", str(fake_tools["config"]), "--threads", "8"])
    assert code == 0
    return {"work": work, **fake_tools}


def test_runner_searches_parses_compresses_and_resumes(searched):
    work = searched["work"]
    fams = runner.read_families(work)
    d = runner.family_dir(work, fams[0])
    done = json.loads((d / "done.json").read_text())
    assert done["status"] == "ok" and done["verified"] and set(done["steps"]) == {
        "mafft", "hhblits", "hhsearch_pfam_enriched", "hhsearch_pdb_enriched", "hhsearch_pfam_family_only",
        "hhsearch_pdb_family_only"}
    assert sorted(p.name for p in d.iterdir()) == sorted(
        ["done.json", "hits.tsv", "pdb_enriched.hhr.gz", "pdb_family_only.hhr.gz", "pfam_enriched.hhr.gz",
         "pfam_family_only.hhr.gz", "uniref.hhr.gz"])
    with gzip.open(d / "pfam_enriched.hhr.gz", "rt") as handle:
        assert hhr.parse_hhr_text(handle.read()).match_columns == fams[0]["rep_length"]
    lines = (d / "hits.tsv").read_text().splitlines()
    # (3 Pfam + 3 PDB) per kind, and the first 20 UniRef30 hits
    assert len(lines) - 1 == 2 * (3 + 3) + 20
    calls = searched["calls"].read_text().count("\n")
    assert runner.main(["families", "--work", str(work), "--config", str(searched["config"]), "--threads", "4"]) == 0
    assert searched["calls"].read_text().count("\n") == calls  # every family skipped


def test_runner_redoes_families_when_the_configuration_changes(searched, tmp_path):
    cfg = dict(searched["dict"], msa_kinds=["family_only"])
    path = tmp_path / "family_only.json"
    path.write_text(json.dumps(cfg))
    assert runner.family_digest(runner.load_config(path)) != runner.family_digest(runner.load_config(searched["config"]))
    before = searched["calls"].read_text().count("hhsearch")
    assert runner.main(["families", "--work", str(searched["work"]), "--config", str(path), "--threads", "4"]) == 0
    # outputs present (gzipped) are reused: no new search runs, only re-parsing
    assert searched["calls"].read_text().count("hhsearch") == before


def test_runner_records_errors_and_load_refuses_partial(dataset, fake_tools, tmp_path):
    work = tmp_path / "work"
    store.export_families(dataset["db"], work, families_db=dataset["families"], log=lambda _m: None)
    fams = runner.read_families(work)
    (work / fams[1]["fasta"]).unlink()
    assert runner.main(["families", "--work", str(work), "--config", str(fake_tools["config"])]) == 1
    assert json.loads((runner.family_dir(work, fams[1]) / "error.json").read_text())["status"] == "error"
    with pytest.raises(store.RemoteHomologyError, match="1 errors"):
        store.load_families(work, tmp_path / "s.duckdb", fake_tools["config"], log=lambda _m: None)
    result = store.load_families(work, tmp_path / "s.duckdb", fake_tools["config"], allow_partial=True,
                                 log=lambda _m: None)
    assert (result["ok"], result["error"]) == (1, 1)


def test_select_shard_and_config_errors(tmp_path):
    assert runner.select_shard(list(range(7)), 1, 3) == [1, 4]
    with pytest.raises(ValueError):
        runner.select_shard([1], 3, 3)
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps({"msa_kinds": ["nope"]}))
    with pytest.raises(runner.ConfigError):
        runner.load_config(bad)
    ok = tmp_path / "ok.json"
    ok.write_text(json.dumps({"databases": {"pfam": {"prefix": "/nonexistent/pfam"}}}))
    with pytest.raises(runner.ConfigError, match="name is required"):
        runner.database(runner.load_config(ok), "pfam")


# --------------------------------------------------------------------------- #
# Staging, activation, SQL surface
# --------------------------------------------------------------------------- #
CDD_ROWS = [  # qseqid replaced by the MD5 of a protein's sequence
    ("s2", "gnl|CDD|990001", 120, 3, 118, 120, 1, 117, 1e-40, 150.0, "cdT0001, SynthA, synthetic profile A"),
    ("s2", "gnl|CDD|990002", 120, 40, 90, 60, 5, 55, 5e-3, 30.0, "cdT0002, DUF9999, domain of unknown function"),
    ("s2", "gnl|CDD|990003", 120, 40, 90, 60, 5, 55, 0.5, 20.0, "cdT0003, Loose, above the E-value cutoff"),
    ("fa_g4", "gnl|CDD|990001", 120, 1, 50, 120, 1, 50, 1e-5, 45.0, "cdT0001, SynthA, synthetic profile A"),
]


def _md5(seq: str) -> str:
    import hashlib

    return hashlib.md5(seq.encode()).hexdigest()


def _fake_cdd(db: Path, work: Path, config: Path) -> None:
    store.export_cdd(db, work, chunk_size=1000, log=lambda _m: None)
    seqs = {pid: s for pid, _, s in _dataset()[0]}
    out = work / "cdd" / "out"
    out.mkdir(parents=True)
    (out / "chunk_00000.tsv").write_text("".join(
        "\t".join(map(str, (_md5(seqs[r[0]]), *r[1:]))) + "\n" for r in CDD_ROWS))
    (out / "chunk_00000.done").write_text(json.dumps({"status": "ok", "evalue": 0.01}))


@pytest.fixture
def staged(dataset, searched, tmp_path):
    cfg = dict(searched["dict"])
    cfg["databases"] = dict(cfg["databases"], cdd={"prefix": str(tmp_path / "Cdd"), "name": "CDD", "version": "9.99"})
    config = tmp_path / "with_cdd.json"
    config.write_text(json.dumps(cfg))
    staging = tmp_path / "remote_homology.duckdb"
    fam_result = store.load_families(searched["work"], staging, searched["config"], log=lambda _m: None)
    _fake_cdd(dataset["db"], tmp_path / "cddwork", config)
    cdd_result = store.load_cdd(tmp_path / "cddwork", staging, config, log=lambda _m: None)
    db = tmp_path / "dataset" / "sharur.duckdb"
    db.parent.mkdir()
    shutil.copyfile(dataset["db"], db)
    return {"staging": staging, "db": db, "fam": fam_result, "cdd": cdd_result, "families": dataset["families"],
            "tag": dataset["tag"], "cddwork": tmp_path / "cddwork", "config": config, "source_db": dataset["db"]}


def _rows(path: Path, sql: str):
    con = duckdb.connect(str(path), read_only=True)
    try:
        return con.execute(sql).fetchall()
    finally:
        con.close()


def test_staging_tables_carry_database_names_versions_and_the_rule(staged):
    s = staged["staging"]
    assert staged["fam"]["ok"] == 2 and staged["fam"]["hit_rows"] == 2 * (2 * 6 + 20)
    assert set(_rows(s, "SELECT DISTINCT source, db_name, db_version FROM family_annotations")) == {
        ("hhsearch_pfam", "Pfam", "99.0"), ("hhsearch_pdb", "PDB100", "2099-01-01"),
        ("hhblits_uniref", "UniRef30", "2099_01")}
    assert _rows(s, "SELECT MIN(probability) FROM family_annotations")[0][0] >= 50
    # rule: Pfam/PDB, probability >= 90 over >= 30 columns; DUF and unknown-function titles flagged
    rule = _rows(s, "SELECT DISTINCT hit_accession, hit_class FROM family_resemblance ORDER BY 1")
    assert rule == [("9yyy_B", "unknown_function"), ("9zzz_A", "characterized"), ("PF90001", "characterized"),
                    ("PF90003", "unknown_function")]
    assert _rows(s, "SELECT COUNT(*) FROM family_annotations WHERE hit_accession = '9xxx_C'")[0][0] == 4
    assert _rows(s, "SELECT DISTINCT hit_description FROM family_annotations WHERE hit_accession = '9xxx_C'") == [
        ("Structure of a synthetic transporter",)]
    # E <= 0.01 kept; fa_g4 and fa_g4b share one sequence, so its hit lands on both
    # staged per unique sequence: s2 (two hits) and the sequence fa_g4 and fa_g4b share (one hit)
    assert staged["cdd"]["hit_rows"] == 3 and staged["cdd"]["sequences_with_hits"] == 2
    seqs = {pid: seq for pid, _, seq in _dataset()[0]}
    assert _rows(s, f"""SELECT m.accession, m.name FROM cdd_hits h JOIN cdd_models m USING (pssm_id)
                        WHERE h.sequence_md5 = '{_md5(seqs["s2"])}' ORDER BY 1""") == [
        ("cdT0001", "SynthA"), ("cdT0002", "DUF9999")]
    assert _rows(s, "SELECT pssm_id, model_length FROM cdd_models ORDER BY 1") == [(990001, 120), (990002, 60)]
    compact = json.loads((staged["cddwork"] / "cdd" / "compact" / "compact.json").read_text())
    assert (compact["raw_rows"], compact["rows"]) == (4, 3)
    # s2's DUF9999 hit lies inside its higher-scoring SynthA hit: one region, SynthA kept
    assert _rows(s, "SELECT m.accession, h.best_in_region FROM cdd_hits h JOIN cdd_models m USING (pssm_id) "
                    "ORDER BY h.score DESC") == [("cdT0001", True), ("cdT0001", True), ("cdT0002", False)]
    assert _rows(s, "SELECT hit_rows, best_rows, max_overlap FROM cdd_runs") == [(3, 2, 0.5)]
    assert staged["cdd"]["best_in_region_rows"] == 2


def test_best_per_region_keeps_the_top_bit_score_of_each_region():
    from sharur.architecture import MAX_OVERLAP

    assert store.MAX_OVERLAP == MAX_OVERLAP  # one overlap default across Pfam architectures and CDD regions
    hits = [(150.0, 1e-40, 3, 3, 118),   # SynthA over the whole protein
            (30.0, 5e-3, 2, 40, 90),     # inside it
            (80.0, 1e-20, 1, 130, 200),  # a second region
            (80.0, 1e-25, 4, 140, 205),  # same region, same bits, lower E-value -> wins the tie
            (20.0, 1e-3, 5, 100, 135)]   # shares 19 of its 36 residues with the first region
    assert store.best_per_region(hits) == [True, False, False, True, False]
    assert store.best_per_region(hits, max_overlap=0.6) == [True, False, False, True, True]
    assert store.best_per_region([(10.0, 1.0, 9, 5, 50)]) == [True]
    assert store.best_per_region([]) == []


def test_activation_needs_the_clustering_run(staged):
    con = store.connect(staged["db"])
    try:
        with pytest.raises(store.RemoteHomologyError, match="no protein-family tables"):
            store.load_staging(con, [staged["staging"]], log=lambda _m: None)
    finally:
        con.close()


def test_activation_loads_rows_and_the_sql_surface_reads_resembles(staged):
    con = store.connect(staged["db"])
    try:
        out = store.load_staging(con, [staged["staging"]], families_staging=staged["families"], log=lambda _m: None)
        again = store.load_staging(con, [staged["staging"]], families_staging=staged["families"],
                                   log=lambda _m: None)
    finally:
        con.close()
    assert out["family_tier"]["hit_rows"] == again["family_tier"]["hit_rows"] == 64
    assert again["cdd"]["rows_removed"] == again["cdd"]["rows_loaded"] == 3 and again["cdd"]["rows_all"] == 4
    db = staged["db"]
    assert _rows(db, "SELECT COUNT(*) FROM family_annotations") == [(64,)]
    assert _rows(db, "SELECT model_set, model_version FROM annotation_source_versions WHERE source = 'cdd'") == [
        ("CDD", "9.99")]
    # best model per region by default: s2 keeps SynthA, the shared sequence of fa_g4 and fa_g4b its one hit
    assert _rows(db, "SELECT COUNT(*), COUNT(DISTINCT protein_id) FROM annotations WHERE source = 'cdd'") == [(3, 3)]
    assert _rows(db, "SELECT COUNT(DISTINCT annotation_id) FROM annotations") == [(5,)]
    assert status_cdd(db) == {"rows": 3, "rows_stored": "best_per_region", "rows_all_hits": 4, "max_overlap": 0.5}
    status = again["status"]
    assert status["family_tier"]["families_resembling_characterized"] == 2 and status["cdd"]["proteins"] == 3
    assert status["family_tier"]["families_resembling_unknown_function_only"] == 0
    assert set(status["family_tier"]["by_clustering_run"]) == {"synthetic-permissive"}

    from sharur.sql import connect as sql_connect

    con = sql_connect(db)
    try:
        q = lambda sql: con.execute(sql).fetchall()  # noqa: E731
        assert q("SELECT DISTINCT evidence FROM feature_proteins('cdd:cdT0001')") == [("resembles",)]
        assert q("SELECT COUNT(*) FROM feature_proteins('SynthA')") == [(0,)]  # unqualified: curated sources only
        assert q("SELECT DISTINCT evidence FROM feature_proteins('DomA')") == [("observed",)]
        summary = q("SELECT pfam, other, n_annotations, resembles FROM protein_summary('s2')")
        assert summary == [([], [], 0, ["cdd:cdT0001"])]
        tag = staged["tag"]
        rows = q("SELECT tier, family_id, source, hit_accession, hit_class, evidence FROM resembles_of('fa_g2')")
        assert {r[3] for r in rows} == {"PF90001", "PF90003", "9zzz_A", "9yyy_B"}
        assert {r[1] for r in rows} == {f"fam:{tag}:fa_g1"} and {r[5] for r in rows} == {"resembles"}
        cdd_rows = q("SELECT tier, db_name, db_version, hit_accession, unknown_function FROM resembles_of('s2') "
                     "ORDER BY hit_accession")
        assert cdd_rows == [("protein", "CDD", "9.99", "cdT0001", False)]
        strict = q("SELECT COUNT(*) FROM resembles_of('fa_g2', min_probability := 99)")
        assert strict == [(2,)]  # PF90001 at 99.6, once per MSA kind
        members = q("SELECT DISTINCT protein_id FROM resembles_proteins('PF90001') ORDER BY 1")
        assert [m[0] for m in members] == ["fa_g1", "fa_g2", "fa_g3", "fa_g4", "fa_g4b", "fd_0", "fd_1", "fd_2",
                                           "fd_3", "fd_4", "fd_5", "fd_6", "fd_7"]
        grid = q("SELECT bin_id, present, n_proteins FROM resembles_presence('cdT0001') ORDER BY bin_id")
        assert grid == [("g1", False, 0), ("g2", False, 0), ("g3", False, 0), ("g4", True, 2), ("g5", True, 1)]
        nb = q("SELECT resembles FROM neighbors('s2', k := 0)")
        assert nb == [(["cdd:cdT0001"],)]
    finally:
        con.close()


def status_cdd(db: Path) -> dict:
    con = store.connect(db, read_only=True)
    try:
        info = store.status(con)["cdd"]
    finally:
        con.close()
    return {k: info[k] for k in ("rows", "rows_stored", "rows_all_hits", "max_overlap")}


def test_activation_can_store_every_cdd_hit(staged):
    con = store.connect(staged["db"])
    try:
        out = store.load_staging(con, [staged["staging"]], families_staging=staged["families"], cdd_rows="all",
                                 log=lambda _m: None)
        with pytest.raises(store.RemoteHomologyError, match="cdd_rows must be one of"):
            store.load_staging(con, [staged["staging"]], cdd_rows="some", log=lambda _m: None)
    finally:
        con.close()
    assert out["cdd"]["rows_loaded"] == out["cdd"]["rows_all"] == 4
    assert status_cdd(staged["db"]) == {"rows": 4, "rows_stored": "all", "rows_all_hits": 4, "max_overlap": 0.5}
    from sharur.sql import connect as sql_connect

    con = sql_connect(staged["db"])
    try:
        rows = con.execute("SELECT hit_accession, unknown_function FROM resembles_of('s2') ORDER BY 1").fetchall()
    finally:
        con.close()
    assert rows == [("cdT0001", False), ("cdT0002", True)]


def test_activation_refuses_changed_sequences(staged, tmp_path):
    cdd_only = tmp_path / "cdd_only.duckdb"
    store.load_cdd(staged["cddwork"], cdd_only, staged["config"], log=lambda _m: None)
    con = store.connect(staged["db"])
    try:
        con.execute("UPDATE proteins SET sequence = 'MKV' WHERE protein_id = 's2'")
        with pytest.raises(store.RemoteHomologyError, match="differ from the 19"):
            store.load_staging(con, [cdd_only], log=lambda _m: None)
        with pytest.raises(pf.FamilyError, match="changed since clustering"):
            store.load_staging(con, [staged["staging"]], families_staging=staged["families"], log=lambda _m: None)
        assert con.execute("SELECT COUNT(*) FROM annotations WHERE source = 'cdd'").fetchone()[0] == 0
    finally:
        con.close()


def test_activate_command_stages_publishes_and_reports(staged):
    from sharur.cli import app

    db = staged["db"]
    result = CliRunner().invoke(app, ["remote-homology", "activate", "--db", str(db), "--staging",
                                      str(staged["staging"]), "--families-staging", str(staged["families"])])
    assert result.exit_code == 0, result.output
    assert "family_annotations" in result.output and "Publish with" in result.output
    assert _rows(db, "SELECT COUNT(*) FROM annotations WHERE source = 'cdd'") == [(0,)]  # untouched until publish
    result = CliRunner().invoke(app, ["remote-homology", "activate", "--db", str(db), "--staging",
                                      str(staged["staging"]), "--families-staging", str(staged["families"]),
                                      "--publish"])
    assert result.exit_code == 0, result.output and "Published" in result.output
    assert list(db.parent.glob("sharur.duckdb.pre-remote-homology-*"))
    result = CliRunner().invoke(app, ["remote-homology", "status", "--db", str(db)])
    info = json.loads(result.stdout)
    assert info["cdd"]["rows"] == 3 and info["cdd"]["rows_all_hits"] == 4
    assert info["family_tier"]["families_resembling_any"] == 2

    from sharur.capabilities import _remote_homology_check
    from sharur.storage.duckdb_store import DuckDBStore as Store

    with Store(db, read_only=True) as s:
        cap = _remote_homology_check(s)
    assert cap.state.value == "available" and "resembles" in cap.summary


def test_predicates_leave_cdd_rows_out():
    from sharur.predicates.generator import AnnotationRecord, PredicateGenerator, ProteinRecord
    from sharur.predicates_v2.generator import AtomGenerator

    protein = ProteinRecord("p", sequence_length=200)
    cdd = [AnnotationRecord("cdd", "cdT0001", "SynthA", "synthetic profile A", 1e-30, 100.0)]
    assert "unannotated" in PredicateGenerator().generate_for_protein(protein, cdd)
    atoms = AtomGenerator(predict_topology=False).generate_atoms(protein, cdd)
    ids = {a.atom_id for a in atoms}
    assert "unannotated" in ids and not any("cdT0001" in i or "cdd" in i for i in ids)


def test_ingest_dag_adds_stages_10_and_11(tmp_path):
    from sharur.ingest_cli import app

    cfg = tmp_path / "rh.json"
    cfg.write_text("{}")
    args = ["--input-dir", str(tmp_path), "--data-dir", str(tmp_path / "d"), "--dry-run"]
    result = CliRunner().invoke(app, [*args, "--with-families", "--with-cdd", "--with-family-homology",
                                      "--remote-homology-config", str(cfg), "--family-homology-run", "runX"])
    assert result.exit_code == 0, result.output
    assert "10_remote_homology_cdd.py" in result.output and "11_remote_homology_families.py" in result.output
    assert "runX" in result.output
    assert CliRunner().invoke(app, [*args, "--with-cdd"]).exit_code != 0
    assert CliRunner().invoke(app, [*args, "--with-family-homology", "--remote-homology-config",
                                    str(cfg)]).exit_code != 0


# --------------------------------------------------------------------------- #
# Real tools (skipped when not installed)
# --------------------------------------------------------------------------- #
def _alignment(rng: random.Random, base: str, n: int, identity: float) -> str:
    return "".join(f">m{i}\n{_mutate(rng, base, identity)}\n" for i in range(n))


def _hhsuite_data() -> Path | None:
    """Directory with HH-suite's context libraries (cs219.lib), from $HHLIB or the install prefix."""
    candidates = [Path(os.environ["HHLIB"]) / "data"] if os.environ.get("HHLIB") else []
    exe = shutil.which("hhsearch")
    if exe:
        real = Path(exe).resolve()
        try:  # a pyenv shim resolves to itself; ask pyenv for the real binary
            out = subprocess.run(["pyenv", "which", "hhsearch"], check=False, capture_output=True, text=True)
            if out.returncode == 0 and out.stdout.strip():
                real = Path(out.stdout.strip()).resolve()
        except OSError:
            pass
        candidates += [real.parent.parent / "data", real.parent.parent / "share" / "hhsuite" / "data"]
    return next((c for c in candidates if (c / "cs219.lib").exists()), None)


def _build_hh_db(db_dir: Path, prefix: str, entries: list[tuple[str, str]]) -> None:
    """HH-suite database ``prefix`` (_hhm, _a3m and, when the context library is found, _cs219)."""
    hhms, a3ms = [], []
    for i, (label, aln) in enumerate(entries):
        fa = db_dir / f"e{i}.fa"
        records = aln.split(">")[1:]
        fa.write_text(f">{label}\n" + records[0].split("\n", 1)[1] + "".join(">" + r for r in records[1:]))
        subprocess.run(["hhmake", "-i", str(fa), "-o", str(db_dir / f"e{i}.hhm"), "-name", label, "-M", "50",
                        "-v", "0"], check=True, capture_output=True)
        subprocess.run(["hhconsensus", "-i", str(fa), "-oa3m", str(db_dir / f"e{i}.a3m"), "-M", "50", "-v", "0"],
                       check=True, capture_output=True)
        hhms.append(str(db_dir / f"e{i}.hhm"))
        a3ms.append(str(db_dir / f"e{i}.a3m"))
    for part, files in (("hhm", hhms), ("a3m", a3ms)):
        (db_dir / f"{part}.txt").write_text("\n".join(files) + "\n")
        subprocess.run(["ffindex_build", "-s", "-f", str(db_dir / f"{part}.txt"),
                        str(db_dir / f"{prefix}_{part}.ffdata"), str(db_dir / f"{prefix}_{part}.ffindex")],
                       check=True, capture_output=True)
    data = _hhsuite_data()
    if data is not None and shutil.which("cstranslate"):
        subprocess.run(["cstranslate", "-f", "-x", "0.3", "-c", "4", "-I", "a3m", "-i", str(db_dir / f"{prefix}_a3m"),
                        "-o", str(db_dir / f"{prefix}_cs219"), "-A", str(data / "cs219.lib"),
                        "-D", str(data / "context_data.lib")], check=True, capture_output=True)


@requires_hhsuite
def test_family_tier_with_hhsuite_finds_the_matching_profile(dataset, tmp_path):
    rng = random.Random(7)
    a = _dataset()[2]["a"]
    db_dir = tmp_path / "hhdb"
    db_dir.mkdir()
    _build_hh_db(db_dir, "synth", [("PF90001.1 SynthA; synthetic", _alignment(rng, a, 12, 0.7)),
                                   ("PF90002.1 Unrelated; synthetic", _alignment(rng, _random(rng, 130), 12, 0.7))])
    config = tmp_path / "config.json"
    config.write_text(json.dumps({
        "databases": {"pfam": {"prefix": str(db_dir / "synth"), "name": "Pfam", "version": "test", "kind": "pfam"}},
        "family_targets": ["pfam"], "msa_kinds": ["family_only"]}))
    work = tmp_path / "work"
    store.export_families(dataset["db"], work, families_db=dataset["families"], log=lambda _m: None)
    code = runner.main(["families", "--work", str(work), "--config", str(config), "--threads", "2"])
    errors = list(work.glob("out/*/*/error.json"))
    if code and errors and "cs219" in errors[0].read_text():
        pytest.skip("this HH-suite needs a _cs219 index and its context library was not found (set HHLIB)")
    assert code == 0
    staging = tmp_path / "s.duckdb"
    store.load_families(work, staging, config, log=lambda _m: None)
    best = _rows(staging, f"""SELECT hit_accession, probability FROM family_resemblance
                               WHERE family_id = 'fam:{dataset["tag"]}:fa_g1' ORDER BY probability DESC LIMIT 1""")
    assert best and best[0][0] == "PF90001" and best[0][1] >= 90
    assert not _rows(staging, f"""SELECT 1 FROM family_resemblance
                                  WHERE family_id = 'fam:{dataset["tag"]}:fd_0' AND hit_accession = 'PF90001'""")
    done = json.loads((runner.family_dir(work, runner.read_families(work)[0]) / "done.json").read_text())
    assert done["verified"]


@requires_blast
def test_cdd_with_rpsblast_finds_the_synthetic_profile(dataset, tmp_path):
    rng = random.Random(11)
    cdd = _dataset()[2]["cdd"]
    aln = tmp_path / "aln.fa"
    aln.write_text(_alignment(rng, cdd, 10, 0.8))
    (tmp_path / "q.fa").write_text(f">q\n{cdd}\n")
    smp = tmp_path / "p.smp"
    subprocess.run(["psiblast", "-in_msa", str(aln), "-subject", str(tmp_path / "q.fa"), "-out_pssm", str(smp)],
                   check=True, capture_output=True)
    text = smp.read_text()
    text = text.replace("      },\n      inst {",
                        '      },\n      descr {\n        title "cdT0001, SynthA, synthetic profile A"\n      },'
                        "\n      inst {", 1)
    smp.write_text(re.sub(r'local str "[^"]*"', 'general { db "CDD", tag id 990001 }', text, count=1))
    (tmp_path / "list.pn").write_text(str(smp) + "\n")
    subprocess.run(["makeprofiledb", "-in", str(tmp_path / "list.pn"), "-out", str(tmp_path / "Cdd"),
                    "-dbtype", "rps", "-title", "synthCDD"], check=True, capture_output=True)
    config = tmp_path / "config.json"
    config.write_text(json.dumps({"databases": {"cdd": {"prefix": str(tmp_path / "Cdd"), "name": "CDD",
                                                        "version": "test"}}}))
    work = tmp_path / "work"
    store.export_cdd(dataset["db"], work, chunk_size=10, log=lambda _m: None)
    assert runner.main(["cdd", "--work", str(work), "--config", str(config), "--threads", "2"]) == 0
    result = store.load_cdd(work, tmp_path / "s.duckdb", config, log=lambda _m: None)
    assert result["sequences_with_hits"] >= 1
    seqs = {pid: seq for pid, _, seq in _dataset()[0]}
    hits = _rows(tmp_path / "s.duckdb", "SELECT h.sequence_md5, m.accession, m.name FROM cdd_hits h "
                                       "JOIN cdd_models m USING (pssm_id) ORDER BY 1")
    assert (_md5(seqs["s2"]), "cdT0001", "SynthA") in hits
    assert all(r[1] == "cdT0001" for r in hits)
