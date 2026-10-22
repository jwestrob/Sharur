"""rebuild_database: restructure a foreign-keyed DuckDB file by copying it table by table."""

import duckdb
import pytest

from sharur.storage.rebuild import rebuild_database, table_order


def _source(path):
    con = duckdb.connect(str(path))
    con.execute("""
        CREATE SEQUENCE seq_ann START 1;
        CREATE TABLE bins (bin_id VARCHAR PRIMARY KEY, n_contigs INTEGER);
        CREATE TABLE contigs (contig_id VARCHAR PRIMARY KEY, bin_id VARCHAR,
                              FOREIGN KEY (bin_id) REFERENCES main.bins(bin_id));
        CREATE TABLE proteins (protein_id VARCHAR PRIMARY KEY, contig_id VARCHAR NOT NULL, bin_id VARCHAR,
                               sequence VARCHAR, FOREIGN KEY (contig_id) REFERENCES main.contigs(contig_id));
        CREATE TABLE annotations (annotation_id INTEGER PRIMARY KEY DEFAULT nextval('seq_ann'),
                                  protein_id VARCHAR NOT NULL, accession VARCHAR,
                                  FOREIGN KEY (protein_id) REFERENCES main.proteins(protein_id));
        CREATE INDEX idx_proteins_bin ON proteins(bin_id);
        CREATE VIEW pfam_only AS SELECT * FROM annotations WHERE accession LIKE 'PF%';
        INSERT INTO bins VALUES ('keep', 1), ('drop', 1);
        INSERT INTO contigs VALUES ('c1', 'keep'), ('c2', 'drop');
        INSERT INTO proteins VALUES ('p1', 'c1', 'keep', 'MK'), ('p2', 'c2', 'drop', 'MA'), ('p3', 'c2', 'drop', 'MV');
        INSERT INTO annotations (protein_id, accession) VALUES ('p1', 'PF00001'), ('p2', 'PF00002'), ('p3', 'K00001');
    """)
    # In-place restructuring is what DuckDB refuses: the reason this helper exists.
    with pytest.raises(duckdb.ConstraintException):
        con.execute("UPDATE proteins SET bin_id = 'keep' WHERE protein_id = 'p3'")
    con.close()


def test_table_order_follows_foreign_keys():
    ddl = {"annotations": "CREATE TABLE annotations(x INT, FOREIGN KEY (x) REFERENCES main.proteins(x))",
           "proteins": "CREATE TABLE proteins(x INT, FOREIGN KEY (x) REFERENCES main.contigs(x))",
           "contigs": "CREATE TABLE contigs(x INT)"}
    assert table_order(ddl) == ["contigs", "proteins", "annotations"]
    with pytest.raises(ValueError, match="Cyclic"):
        table_order({"a": "REFERENCES b(x)", "b": "REFERENCES a(x)"})


def test_rebuild_drops_a_genome_and_moves_a_protein(tmp_path):
    src, dst = tmp_path / "src.duckdb", tmp_path / "dst.duckdb"
    _source(src)
    counts = rebuild_database(src, dst, {
        "bins": "SELECT * FROM src.bins WHERE bin_id = 'keep'",
        "contigs": "SELECT contig_id, 'keep' AS bin_id FROM src.contigs",
        "proteins": "SELECT * REPLACE ('keep' AS bin_id) FROM src.proteins WHERE protein_id <> 'p2'",
        "annotations": "SELECT * FROM src.annotations WHERE protein_id <> 'p2'",
    }, setup_sql=["CREATE TEMP TABLE _unused AS SELECT 1 AS x"])
    assert counts == {"bins": 1, "contigs": 2, "proteins": 2, "annotations": 2}
    con = duckdb.connect(str(dst))
    assert con.execute("SELECT protein_id, bin_id FROM proteins ORDER BY 1").fetchall() == [
        ("p1", "keep"), ("p3", "keep")]
    # constraints, index, view and sequence survive
    with pytest.raises(duckdb.ConstraintException):
        con.execute("INSERT INTO proteins VALUES ('p9', 'missing', 'keep', 'M')")
    assert {r[0] for r in con.execute("SELECT index_name FROM duckdb_indexes()").fetchall()} == {"idx_proteins_bin"}
    assert con.execute("SELECT COUNT(*) FROM pfam_only").fetchone()[0] == 1
    con.execute("INSERT INTO annotations (protein_id, accession) VALUES ('p1', 'PF00009')")
    # the sequence continues past the copied IDs (a gap is possible: DuckDB may draw once during the copy)
    new_id = con.execute("SELECT annotation_id FROM annotations WHERE accession = 'PF00009'").fetchone()[0]
    assert new_id > 3
    assert {r[0] for r in con.execute("SELECT table_name FROM duckdb_tables()").fetchall()} == {
        "bins", "contigs", "proteins", "annotations"}


def test_rebuild_refuses_unknown_tables_and_cleans_up(tmp_path):
    src, dst = tmp_path / "src.duckdb", tmp_path / "dst.duckdb"
    _source(src)
    with pytest.raises(KeyError, match="absent"):
        rebuild_database(src, dst, {"nope": "SELECT 1"})
    assert not dst.exists()
    # a select that breaks a foreign key fails and leaves no partial file
    with pytest.raises(duckdb.ConstraintException):
        rebuild_database(src, dst, {"bins": "SELECT * FROM src.bins WHERE bin_id = 'keep'"})
    assert not dst.exists()
    dst.touch()
    with pytest.raises(FileExistsError):
        rebuild_database(src, dst)
