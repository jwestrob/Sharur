"""The Sharur SQL surface (sharur/sql/prelude.sql): TEMP macros on read-only connections."""

from __future__ import annotations

import json
from pathlib import Path

import duckdb
import pytest
from typer.testing import CliRunner

from sharur.sql import (
    available_tables,
    install_prelude,
    looks_like_sequence,
    parse_prelude,
    prelude_text,
    sequence_leak,
    write_fasta,
)
from sharur.storage.duckdb_store import DuckDBStore


SEQ = "M" + "KLVAEGRT" * 6  # 49 residues, upper case: what a leak looks like
NEIGHBOR = "pfam:PF90002,pfam:PF90003"
TAX = "d__Archaea;p__Testarchaeota;c__Testarchaeia;o__{o};f__F{o};g__G{o};s__G{o} sp{n}"

# Gene layout (gene_index order per contig). All feature names are invented:
# A = anchor family (KO K90001 or Pfam PF90001 AnchorFam), N = neighbor family (PF90002 NeighborFam /
# PF90003 NeighborFam2), S = split feature (K90003), P = partner (K90002).
#   g1 (OrdA, 95%):  c1  p0+  p1+ N  p2+ A  p3-  p4+ P  p5+
#   g2 (OrdA, 92%):  c1  q0- P  q1- A(+AnchorFam)  q2- N  q3+  q4-
#   g3 (OrdB, 70%):  c1  r0+ A(AnchorFam only)  r1+ S
#   g4 (OrdB, 96%):  c1  s0- (no annotation)  s1+ A  s2+ S ;  c2  s3+ (CBASS system)  s4+ (hydrogenase predicate)
GENES = {
    "g1": {"g1_c1": [("p0", "+"), ("p1", "+"), ("p2", "+"), ("p3", "-"), ("p4", "+"), ("p5", "+")]},
    "g2": {"g2_c1": [("q0", "-"), ("q1", "-"), ("q2", "-"), ("q3", "+"), ("q4", "-")]},
    "g3": {"g3_c1": [("r0", "+"), ("r1", "+")]},
    "g4": {"g4_c1": [("s0", "-"), ("s1", "+"), ("s2", "+")], "g4_c2": [("s3", "+"), ("s4", "+")]},
}
GENOMES = [("g1", 95.0, 1.0, "OrdA", 1), ("g2", 92.0, 2.0, "OrdA", 2), ("g3", 70.0, 1.0, "OrdB", 3),
           ("g4", 96.0, 0.0, "OrdB", 4)]
ANNOTATIONS = [
    ("p1", "pfam", "PF90002", "NeighborFam", 5), ("p1", "pfam", "PF90004", "NeighborFam_C", 120),
    ("p1", "kofam", "K90004", "K90004", None),
    ("p2", "kofam", "K90001", "K90001", None),
    ("p2", "kegg", "K90001", "anchor family enzyme", None),
    ("p4", "kofam", "K90002", "K90002", None),
    ("q0", "kofam", "K90002", "K90002", None),
    ("q1", "kofam", "K90001", "K90001", None), ("q1", "pfam", "PF90001", "AnchorFam", 3),
    ("q2", "pfam", "PF90003", "NeighborFam2", 2),
    ("r0", "pfam", "PF90001", "AnchorFam", 1),
    ("r1", "kofam", "K90003", "K90003", None),
    ("s1", "kofam", "K90001", "K90001", None),
    ("s2", "kofam", "K90003", "K90003", None),
    ("s3", "defensefinder", "CBASS__Cyclase", "CBASS__Cyclase", None),
    ("s3", "defensefinder_system", "CBASS__Cyclase", "CBASS/CBASS_I", None),
]


def build_db(path: Path) -> Path:
    store = DuckDBStore(path)
    con = store.conn
    for bin_id, comp, cont, order, n in GENOMES:
        con.execute("INSERT INTO bins VALUES (?, ?, ?, ?, ?, 10000)",
                    [bin_id, comp, cont, TAX.format(o=order, n=n), len(GENES[bin_id])])
    for bin_id, contigs in GENES.items():
        for contig, genes in contigs.items():
            con.execute("INSERT INTO contigs (contig_id, bin_id, length) VALUES (?, ?, 10000)", [contig, bin_id])
            for i, (pid, strand) in enumerate(genes):
                con.execute(
                    "INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, gene_index, "
                    "sequence, sequence_length) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                    [pid, contig, bin_id, 100 + 1000 * i, 900 + 1000 * i, strand, i, SEQ + "*", len(SEQ)])
                preds = ["hydrogenase"] if pid == "s4" else ["enzyme"]
                con.execute("INSERT INTO protein_predicates (protein_id, predicates) VALUES (?, ?)", [pid, preds])
    for n, (pid, source, acc, name, start) in enumerate(ANNOTATIONS, 1):
        con.execute("INSERT INTO annotations (annotation_id, protein_id, source, accession, name, start_aa) "
                    "VALUES (?, ?, ?, ?, ?, ?)", [n, pid, source, acc, name, start])
    con.execute("INSERT INTO defense_systems (system_id, genome_id, contig_id, system_type, system_subtype) "
                "VALUES ('sysA', 'g4', 'g4_c2', 'CBASS', 'CBASS_I')")
    con.execute("INSERT INTO system_proteins (system_id, protein_id, system_source) "
                "VALUES ('sysA', 's3', 'defensefinder_system')")
    store.close()
    return path


@pytest.fixture
def db(tmp_path: Path) -> Path:
    return build_db(tmp_path / "sharur.duckdb")


@pytest.fixture
def con(db: Path):
    store = DuckDBStore(db, read_only=True)
    yield store.conn
    store.close()


def rows(con, sql: str) -> list[dict]:
    cur = con.execute(sql)
    cols = [d[0] for d in cur.description]
    return [dict(zip(cols, r)) for r in cur.fetchall()]


# --------------------------------------------------------------------------- installation
def test_prelude_file_parses_into_documented_objects():
    parsed = parse_prelude()
    names = [o.name for o in parsed.objects]
    assert names[:5] == ["strand_sign", "id_list", "rank_of", "genes", "genomes"]
    assert {"feature_proteins", "feature_presence", "neighbors", "slot", "slot_tally",
            "protein_summary", "anchor_proteins", "neighbor_features"} <= set(names)
    assert parsed.version == "1" and parsed.requires == ("proteins", "annotations", "bins")
    for obj in parsed.objects:
        assert obj.doc and "CREATE OR REPLACE TEMP MACRO" in obj.statement
        assert "CREATE VIEW" not in obj.statement.upper() and "CREATE TABLE" not in obj.statement.upper()
    assert next(o for o in parsed.objects if o.name == "neighbors").signature == "neighbors(anchor, k := 5)"


def test_store_installs_prelude_on_read_only_connection(con):
    objects = {r["name"]: r for r in rows(con, "FROM sharur_sql()")}
    assert {"genes", "neighbors", "slot", "slot_tally", "feature_presence", "write_fasta", "sharur_sql"} <= set(objects)
    assert objects["slot"]["version"] == "1"


def test_prelude_leaves_catalog_and_seal_snapshot_unchanged(db, monkeypatch):
    from sharur.dataset_seal import _database_snapshot

    plain = duckdb.connect(str(db), read_only=True)
    before = (plain.execute("FROM information_schema.tables ORDER BY ALL").fetchall(),
              plain.execute("FROM information_schema.columns ORDER BY ALL").fetchall(),
              plain.execute("SHOW TABLES").fetchall())
    report = install_prelude(plain)
    after = (plain.execute("FROM information_schema.tables ORDER BY ALL").fetchall(),
             plain.execute("FROM information_schema.columns ORDER BY ALL").fetchall(),
             plain.execute("SHOW TABLES").fetchall())
    plain.close()
    assert report.ok and before == after
    with_prelude = _database_snapshot(db)
    monkeypatch.setenv("SHARUR_SQL_PRELUDE", "0")
    assert _database_snapshot(db) == with_prelude


def test_opt_out_and_missing_tables(db, tmp_path):
    store = DuckDBStore(db, read_only=True, sql_prelude=False)
    with pytest.raises(duckdb.Error):
        store.conn.execute("FROM genes()")
    store.close()
    empty = duckdb.connect(str(tmp_path / "empty.duckdb"))
    report = install_prelude(empty)
    assert report.skipped == "missing tables: proteins, annotations, bins" and not report.ok
    assert report.installed == []
    with pytest.raises(RuntimeError, match="missing tables"):
        install_prelude(empty, strict=True)


def test_generic_client_runs_the_file_after_attach_and_use(db):
    mem = duckdb.connect()
    mem.execute(f"ATTACH '{db}' AS ds (READ_ONLY)")
    mem.execute("USE ds")
    # `sharur sql prelude --db`: the script for this database, without the macros of absent
    # optional layers (DuckDB 1.3+ binds a table macro's tables when it is created).
    mem.execute(prelude_text(available_tables(mem)))
    assert mem.execute("SELECT count(*) FROM genes()").fetchone()[0] == 18
    names = {r[0] for r in mem.execute("SELECT name FROM sharur_sql()").fetchall()}
    assert "slot_tally" in names and "family_proteins" not in names
    assert mem.execute("SELECT count(*) FROM slot('K90001')").fetchone()[0] == 3


def test_prelude_script_leaves_out_only_absent_layers():
    core = {"proteins", "annotations", "bins"}
    families = core | {"protein_families", "family_stats", "protein_family_runs"}
    assert "MACRO family_proteins" not in prelude_text(core) and "MACRO neighbors" in prelude_text(core)
    assert "MACRO family_proteins" in prelude_text(families)
    assert "MACRO resembles_of" not in prelude_text(families)
    assert "MACRO resembles_of" in prelude_text(families | {"family_annotations"})
    assert "MACRO resembles_presence" in prelude_text()


# --------------------------------------------------------------------------- tables and features
def test_genes_has_no_sequence_and_genomes_split_ranks(con):
    cols = [d[0] for d in con.execute("FROM genes() LIMIT 0").description]
    assert "sequence" not in cols and {"protein_id", "strand_sign", "sequence_length"} <= set(cols)
    assert con.execute("SELECT strand_sign('-1'), strand_sign('+'), strand_sign('?')").fetchone() == (-1, 1, None)
    g = rows(con, "SELECT bin_id, phylum, class, \"order\", species FROM genomes() ORDER BY bin_id")
    assert g[0] == {"bin_id": "g1", "phylum": "Testarchaeota", "class": "Testarchaeia", "order": "OrdA",
                    "species": "GOrdA sp1"}
    assert con.execute("SELECT id_list(' a, b ,,c'), id_list(['x', 'y'])").fetchone() == (["a", "b", "c"], ["x", "y"])


def test_feature_specs_and_evidence(con):
    def hits(spec):
        return {(r["protein_id"], r["evidence"]) for r in rows(con, f"FROM feature_proteins('{spec}')")}

    assert hits("K90001") == {("p2", "observed"), ("q1", "observed"), ("s1", "observed")}
    assert hits("ko:K90001") == hits("K90001")
    assert hits("kofam:K90001") == hits("K90001")
    assert hits("pfam:AnchorFam") == {("q1", "observed"), ("r0", "observed")}
    assert hits("PF90001") == hits("pfam:AnchorFam")
    assert {p for p, _ in hits("pfam:*")} == {"p1", "q1", "q2", "r0"}
    assert hits("CBASS__Cyclase") == {("s3", "observed"), ("s3", "caller")}  # raw profile hit vs projection
    assert hits("defense:CBASS") == {("s3", "caller")}
    assert hits("defense:CBASS_I") == {("s3", "caller")}
    assert hits("predicate:hydrogenase") == {("s4", "predicate")}
    assert hits("K99999") == set()
    both = rows(con, "SELECT feature, count(DISTINCT protein_id) AS n FROM feature_proteins(['K90001', 'K90003']) "
                     "GROUP BY 1 ORDER BY 1")
    assert both == [{"feature": "K90001", "n": 3}, {"feature": "K90003", "n": 2}]


def test_feature_presence_is_a_full_genome_by_feature_grid(con):
    grid = rows(con, "SELECT bin_id, feature, present, n_proteins, \"order\" FROM feature_presence('K90003,defense:CBASS') "
                     "ORDER BY bin_id, feature")
    assert len(grid) == 8
    present = {(r["bin_id"], r["feature"]) for r in grid if r["present"]}
    assert present == {("g3", "K90003"), ("g4", "K90003"), ("g4", "defense:CBASS")}
    hq = rows(con, "SELECT DISTINCT bin_id FROM feature_presence('K90003', min_completeness := 90) ORDER BY 1")
    assert [r["bin_id"] for r in hq] == ["g1", "g2", "g4"]


# --------------------------------------------------------------------------- neighbors
def test_neighbors_orient_offsets_to_the_anchor_strand(con):
    plus = {r["protein_id"]: r for r in rows(con, "FROM neighbors('p2', k := 2)")}
    assert {p: r["offset"] for p, r in plus.items()} == {"p0": -2, "p1": -1, "p2": 0, "p3": 1, "p4": 2}
    assert plus["p3"]["same_strand"] is False and plus["p1"]["same_strand"] is True
    assert plus["p1"]["pfam"] == ["NeighborFam", "NeighborFam_C"] and plus["p1"]["ko"] == ["K90004"]
    assert plus["p2"]["ko_name"] == "anchor family enzyme"
    assert (plus["p2"]["genes_upstream_to_edge"], plus["p2"]["genes_downstream_to_edge"]) == (2, 3)
    minus = {r["protein_id"]: r for r in rows(con, "FROM neighbors('q1', k := 1)")}
    assert {p: r["offset"] for p, r in minus.items()} == {"q2": -1, "q1": 0, "q0": 1}
    assert (minus["q1"]["genes_upstream_to_edge"], minus["q1"]["genes_downstream_to_edge"]) == (3, 1)
    sysrow = rows(con, "SELECT systems, other FROM neighbors('s4', k := 1) WHERE protein_id = 's3'")[0]
    assert sysrow == {"systems": ["defense:CBASS_I"], "other": ["defensefinder:CBASS__Cyclase"]}


def test_neighbors_of_a_feature_cover_every_anchor(con):
    anchors = rows(con, "SELECT anchor_id, anchor_features FROM neighbors('K90001,pfam:AnchorFam', k := 1) "
                        "WHERE \"offset\" = 0 ORDER BY anchor_id")
    assert [(r["anchor_id"], sorted(r["anchor_features"])) for r in anchors] == [
        ("p2", ["K90001"]), ("q1", ["K90001", "pfam:AnchorFam"]), ("r0", ["pfam:AnchorFam"]), ("s1", ["K90001"])]
    counts = rows(con, "SELECT \"offset\", feature, n_genomes FROM neighbor_features('K90001', k := 1) ORDER BY 1, 2")
    assert {"offset": -1, "feature": "pfam:PF90002", "n_genomes": 1} in counts
    assert {"offset": -1, "feature": "pfam:PF90003", "n_genomes": 1} in counts


# --------------------------------------------------------------------------- slot
def test_slot_single_anchor_upstream_and_downstream(con):
    up = {r["anchor_id"]: r for r in rows(con, "FROM slot('K90001,pfam:AnchorFam')")}
    assert {a: (r["occupant_id"], r["occupant_offset"], r["occupant_same_strand"]) for a, r in up.items()} == {
        "p2": ("p1", -1, True), "q1": ("q2", -1, True), "r0": (None, None, None), "s1": ("s0", -1, False)}
    assert up["r0"]["at_contig_edge"] is True and up["s1"]["pfam"] == []
    assert up["p2"]["order"] == "OrdA" and up["p2"]["side"] == "upstream"
    down = {r["anchor_id"]: r["occupant_id"] for r in rows(con, "FROM slot('K90001', side := 'downstream')")}
    assert down == {"p2": "p3", "q1": "q0", "s1": "s2"}
    hq = {r["anchor_id"] for r in rows(con, "FROM slot('K90001,pfam:AnchorFam', min_completeness := 90)")}
    assert hq == {"p2", "q1", "s1"}


def test_slot_between_anchor_and_partner(con):
    res = rows(con, "FROM slot('K90001', partner := 'K90002', side := 'either') ORDER BY anchor_id, slot_position")
    by = {}
    for r in res:
        by.setdefault(r["anchor_id"], []).append(r)
    assert [(r["partner_id"], r["partner_offset"], r["n_between"], r["occupant_id"]) for r in by["p2"]] == [
        ("p4", 2, 1, "p3")]
    assert by["p2"][0]["occupant_same_strand"] is False and by["p2"][0]["slot_position"] == 1
    assert [(r["partner_id"], r["partner_offset"], r["n_between"], r["occupant_id"]) for r in by["q1"]] == [
        ("q0", 1, 0, None)]
    assert [(r["partner_id"], r["occupant_id"]) for r in by["s1"]] == [(None, None)]
    upstream_only = {r["anchor_id"]: r["partner_id"] for r in rows(
        con, "FROM slot('K90001', partner := 'K90002', side := 'upstream')")}
    assert upstream_only == {"p2": None, "q1": None, "s1": None}


def test_slot_tally_counts_genomes_per_taxon_and_split(con):
    t = rows(con, f"FROM slot_tally('K90001,pfam:AnchorFam', '{NEIGHBOR}', rank := 'order', split := 'K90003')")
    key = {(r["taxon"], r["genome_has_split"]): r for r in t}
    a = key[("OrdA", False)]
    assert (a["genomes_with_anchor"], a["genomes_occupant_same_strand"], a["anchor_proteins"]) == (2, 2, 2)
    b = key[("OrdB", True)]
    assert (b["genomes_with_anchor"], b["genomes_occupant_same_strand"], b["genomes_occupant_no_pfam"],
            b["genomes_at_contig_edge"]) == (2, 0, 1, 1)
    assert b["phylum"] == "Testarchaeota" and b["split_feature"] == "K90003"
    hq = rows(con, f"FROM slot_tally('K90001,pfam:AnchorFam', '{NEIGHBOR}', min_completeness := 90)")
    assert {(r["taxon"], r["genomes_with_anchor"]) for r in hq} == {("OrdA", 2), ("OrdB", 1)}
    assert all(r["genome_has_split"] is None for r in hq)
    pair = rows(con, "FROM slot_tally('K90001', 'pfam:*', partner := 'K90002', side := 'either', rank := 'phylum')")
    assert pair[0]["genomes_with_anchor"] == 3 and pair[0]["genomes_no_partner_in_span"] == 1


# --------------------------------------------------------------------------- sequences stay out
def test_write_fasta_reports_path_and_counts_only(con, db, tmp_path):
    out = tmp_path / "fa" / "slot.faa"
    res = con.execute(f"SELECT write_fasta('{out}', (SELECT list(occupant_id ORDER BY occupant_id) "
                      f"FROM slot('K90001') WHERE occupant_id IS NOT NULL))").fetchone()[0]
    assert res == {"path": str(out.resolve()), "n_sequences": 3, "n_requested": 3, "n_missing": 0}
    assert not looks_like_sequence(res)
    text = out.read_text()
    assert text.count(">") == 3 and text.splitlines()[0] == ">p1" and "*" not in text
    store = DuckDBStore(db, read_only=True, sql_prelude=False)
    direct = write_fasta(store, "s3, nope, s3", tmp_path / "d.faa")
    assert direct["n_sequences"] == 1 and direct["n_missing"] == 1
    store.close()


def test_sequence_guard_and_cli(db):
    assert sequence_leak(["protein_id", "sequence"], [("p1", "x")]) == "result has sequence column(s): sequence"
    assert sequence_leak(["s"], [(SEQ,)]) == "column 's' holds sequence-like text"
    assert sequence_leak(["d"], [("alcohol dehydrogenase [EC:1.1.1.1]",)]) is None
    from sharur.cli import app

    runner = CliRunner()
    ok = runner.invoke(app, ["sql", "query", "SELECT taxon, genomes_with_anchor FROM slot_tally('K90001', "
                             "'pfam:*') ORDER BY 1", "--db", str(db), "--format", "json"])
    assert ok.exit_code == 0, ok.output
    assert json.loads(ok.stdout)["rows"] == [{"taxon": "OrdA", "genomes_with_anchor": 2},
                                             {"taxon": "OrdB", "genomes_with_anchor": 1}]
    leak = runner.invoke(app, ["sql", "query", "SELECT protein_id, upper(sequence) AS s FROM proteins", "--db", str(db)])
    assert leak.exit_code == 2 and SEQ not in leak.stdout
    script = runner.invoke(app, ["sql", "prelude"])
    assert script.exit_code == 0 and "CREATE OR REPLACE TEMP MACRO slot_tally" in script.stdout


def test_reverify_sql_triples_see_the_surface(tmp_path):
    from sharur.reverify import reverify

    root = tmp_path / "repo"
    ds = root / "data" / "ds"
    (ds / "exploration").mkdir(parents=True)
    build_db(ds / "sharur.duckdb")
    finding = {"id": "F1", "title": "slot", "verification": [
        {"claim": "2 OrdA genomes carry a neighbor-family gene upstream of the anchor",
         "query": f"SELECT genomes_occupant_same_strand FROM slot_tally('K90001,pfam:AnchorFam', '{NEIGHBOR}') "
                  "WHERE taxon = 'OrdA'", "expected": 2}]}
    (ds / "exploration" / "findings.jsonl").write_text(json.dumps(finding) + "\n")
    report = reverify([ds], workers=1)
    assert [r.replay.status for r in report.results] == ["OK"], report.results[0].replay
