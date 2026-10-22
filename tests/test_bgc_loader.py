"""GECCO BGC loading into an existing database (sharur.bgc, `sharur load-bgcs`)."""

import json

import pytest
from typer.testing import CliRunner

from sharur.bgc import load_gecco_bgcs
from sharur.storage.duckdb_store import DuckDBStore

HEADER = ("sequence_id\tbgc_id\tstart\tend\taverage_p\tmax_p\ttype\talkaloid_probability\t"
          "nrp_probability\tpolyketide_probability\tripp_probability\tsaccharide_probability\t"
          "terpene_probability\tproteins\tdomains\n")


def _shard(root, name, records, rows):
    shard = root / name
    shard.mkdir(parents=True)
    (shard / "record_map.tsv").write_text(
        "record\tbin_id\tcontig_id\n" + "".join(f"{r}\t{b}\t{c}\n" for r, (b, c) in records.items()))
    (shard / "gecco_version.txt").write_text("gecco 0.9.10\n")
    (shard / "input.clusters.tsv").write_text(HEADER + "".join(rows))
    return shard


@pytest.fixture
def store():
    store = DuckDBStore()
    store.conn.execute("""
        INSERT INTO bins (bin_id) VALUES ('b1'), ('b2');
        INSERT INTO contigs (contig_id, bin_id, length) VALUES ('c1', 'b1', 50000), ('c2', 'b2', 40000);
        INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand) VALUES
            ('c1_1', 'c1', 'b1', 100, 1000, '+'), ('c1_2', 'c1', 'b1', 1200, 2000, '+'),
            ('c1_3', 'c1', 'b1', 30000, 31000, '-'), ('c2_1', 'c2', 'b2', 10, 900, '+');
        INSERT INTO loci (locus_id, locus_type, contig_id, start, end_coord) VALUES
            ('b1::CRISPR1', 'crispr', 'c1', 40000, 40500);
    """)
    return store


def test_load_maps_records_and_numbers_loci_per_genome(store, tmp_path):
    _shard(tmp_path, "shard_000", {"r0": ("b1", "c1"), "r1": ("b2", "c2")}, [
        "r0\tr0_cluster_2\t29000\t31500\t0.9\t0.99\tTerpene\t0\t0\t0\t0\t0\t0.9\tc1_3\tPF01397\n",
        "r0\tr0_cluster_1\t50\t2100\t0.8\t0.95\tPolyketide\t0\t0.1\t0.8\t0\t0\t0\tc1_1;c1_2\tPF00109;PF02801\n",
        "r1\tr1_cluster_1\t5\t950\t0.7\t0.9\tUnknown\t0\t0\t0\t0\t0\t0\tc2_1;missing_9\t\n",
    ])
    stats = load_gecco_bgcs(store.conn, tmp_path)

    assert stats["clusters"] == 3 and stats["genomes"] == 2
    assert stats["dropped_members"] == 1  # missing_9 is not a stored protein
    loci = store.conn.execute(
        "SELECT locus_id, contig_id, start, end_coord FROM loci WHERE locus_type = 'bgc' ORDER BY 1").fetchall()
    assert loci == [("b1::bgc1", "c1", 50, 2100), ("b1::bgc2", "c1", 29000, 31500),
                    ("b2::bgc1", "c2", 5, 950)]
    members = store.conn.execute(
        "SELECT protein_id, position FROM locus_proteins WHERE locus_id = 'b1::bgc1' ORDER BY 2").fetchall()
    assert members == [("c1_1", 0), ("c1_2", 1)]
    meta = json.loads(store.conn.execute(
        "SELECT metadata FROM loci WHERE locus_id = 'b1::bgc1'").fetchone()[0])
    assert meta["caller"] == "gecco" and meta["bgc_type"] == "Polyketide"
    assert meta["gecco_version"] == "gecco 0.9.10"
    assert store.conn.execute("SELECT COUNT(*) FROM bgc_loci").fetchone()[0] == 3
    # other locus types are untouched
    assert store.conn.execute("SELECT COUNT(*) FROM loci WHERE locus_type = 'crispr'").fetchone()[0] == 1


def test_reload_replaces_and_drops_members_from_another_genome(store, tmp_path):
    first = tmp_path / "first"
    _shard(first, "s0", {"r0": ("b1", "c1")}, [
        "r0\tx\t50\t2100\t0.8\t0.95\tNRP\t0\t0.8\t0\t0\t0\t0\tc1_1;c1_2\t\n"])
    load_gecco_bgcs(store.conn, first)
    second = tmp_path / "second"
    # a member stored under another genome is dropped; the cluster keeps its own proteins
    _shard(second, "s0", {"r0": ("b1", "c1")}, [
        "r0\tx\t50\t2100\t0.8\t0.95\tRiPP\t0\t0\t0\t0.8\t0\t0\tc1_1;c2_1\t\n"])
    stats = load_gecco_bgcs(store.conn, second)
    assert stats["replaced"] == 1 and stats["clusters"] == 1 and stats["dropped_members"] == 1
    assert store.conn.execute(
        "SELECT protein_id FROM locus_proteins WHERE locus_id = 'b1::bgc1'").fetchall() == [("c1_1",)]


def test_unmapped_record_or_unfinished_shard_is_an_error(store, tmp_path):
    _shard(tmp_path, "s0", {"r0": ("b1", "c1")}, [
        "r9\tx\t50\t2100\t0.8\t0.95\tNRP\t0\t0.8\t0\t0\t0\t0\tc1_1\t\n"])
    with pytest.raises(ValueError, match="absent from"):
        load_gecco_bgcs(store.conn, tmp_path)
    unfinished = tmp_path / "other"
    (unfinished / "s1").mkdir(parents=True)
    (unfinished / "s1" / "record_map.tsv").write_text("record\tbin_id\tcontig_id\n")
    with pytest.raises(ValueError, match="without a clusters table"):
        load_gecco_bgcs(store.conn, unfinished)


def test_load_bgcs_cli(tmp_path):
    from sharur.cli import app

    db = tmp_path / "sharur.duckdb"
    store = DuckDBStore(db)
    store.conn.execute("""
        INSERT INTO bins (bin_id) VALUES ('b1');
        INSERT INTO contigs (contig_id, bin_id, length) VALUES ('c1', 'b1', 5000);
        INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand)
        VALUES ('c1_1', 'c1', 'b1', 100, 1000, '+');
    """)
    store.close()
    gecco = tmp_path / "gecco"
    _shard(gecco, "s0", {"r0": ("b1", "c1")}, [
        "r0\tx\t50\t1100\t0.8\t0.95\tSaccharide\t0\t0\t0\t0\t0.8\t0\tc1_1\t\n"])
    result = CliRunner().invoke(app, ["load-bgcs", "--db", str(db), "--gecco-dir", str(gecco)])
    assert result.exit_code == 0, result.output
    assert "BGC loci: 1 in 1 genomes" in result.output
    assert "Saccharide: 1" in result.output
