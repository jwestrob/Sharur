"""Contig-edge context: assembly lengths, Prodigal truncation flags, edge status."""

import json
from importlib.machinery import SourceFileLoader
from pathlib import Path

import duckdb
import pytest

from sharur.contig_context import (
    backfill_contig_context,
    describe_edge,
    edge_context,
    fasta_lengths,
    parse_partial,
    prodigal_partials,
)
from sharur.modules import genome_modules, modules_markdown
from sharur.operators.cards import card, card_markdown
from sharur.operators.navigation import get_neighborhood
from sharur.storage.duckdb_store import DuckDBStore


REPO_ROOT = Path(__file__).resolve().parents[1]


def _gene(pid, contig, start, end, index, partial=None, bin_id="b1"):
    return f"('{pid}', '{contig}', '{bin_id}', {start}, {end}, '+', {index}, 300, " + (
        f"'{partial}')" if partial else "NULL)")


@pytest.fixture
def store():
    store = DuckDBStore()
    genes = [
        # c1: 20 kb assembly contig, ten genes; g0 runs off the contig start
        *(_gene(f"g{i}", "c1", 1 + 2000 * i, 1800 + 2000 * i, i, "10" if i == 0 else "00") for i in range(10)),
        # c2: gene-span length only (legacy), three genes
        *(_gene(f"h{i}", "c2", 1 + 1000 * i, 900 + 1000 * i, i) for i in range(3)),
        # c3: circular
        _gene("k0", "c3", 1, 900, 0),
        # protein-only record without coordinates
        _gene("x0", "x0", 0, 900, 0),
    ]
    store.conn.execute(f"""
        INSERT INTO bins (bin_id, completeness, contamination) VALUES ('b1', 80.0, 1.0);
        INSERT INTO contigs (contig_id, bin_id, length, is_circular, length_source) VALUES
            ('c1', 'b1', 20000, FALSE, 'assembly'), ('c2', 'b1', 2900, FALSE, 'gene_span'),
            ('c3', 'b1', 5000, TRUE, 'assembly'), ('x0', 'b1', 900, FALSE, 'gene_span');
        INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, gene_index,
                              sequence_length, partial)
        VALUES {", ".join(genes)};
    """)
    return store


def test_edge_status_by_position(store):
    ctx = edge_context(store, ["g0", "g1", "g2", "g3", "g5", "g9", "h1", "k0", "x0", "missing"])
    assert "missing" not in ctx
    assert ctx["g0"].edge_status == "truncated" and ctx["g0"].truncated_start
    assert ctx["g1"].edge_status == "contig_edge"
    assert (ctx["g1"].genes_to_start, ctx["g1"].bp_to_start) == (1, 2000)
    assert ctx["g3"].edge_status == "interior"
    assert ctx["g5"].edge_status == "interior"
    assert ctx["g9"].edge_status == "contig_edge"
    assert (ctx["g9"].genes_to_end, ctx["g9"].bp_to_end) == (0, 200)
    assert ctx["h1"].bp_to_end is None and ctx["h1"].contig_length is None
    assert ctx["k0"].edge_status == "circular"
    assert ctx["x0"].edge_status == "no_coordinates"
    assert "contig length from gene span" in describe_edge(ctx["h1"])
    assert "runs off contig start" in describe_edge(ctx["g0"])


def test_edge_thresholds_are_parameters(store):
    assert edge_context(store, ["g3"], edge_genes=4, edge_bp=10000)["g3"].edge_status == "contig_edge"


def test_card_reports_contig_position(store):
    result = card(store, "g9")
    assert result["contig_edge"]["edge_status"] == "contig_edge"
    text = card_markdown(result)
    assert "of a 20,000 bp contig" in text and "Contig position: contig_edge" in text
    assert "spanning at least 2,900 bp" in card_markdown(card(store, "h1"))


def test_neighborhood_marks_contig_ends(store):
    hood = get_neighborhood(store, "g1", window=2)
    assert hood.raw["contig_start_in_window"] and not hood.raw["contig_end_in_window"]
    assert "runs off contig end: g0" in hood.data
    assert {p["protein_id"]: p["edge_status"] for p in hood.raw["proteins"]}["g0"] == "truncated"


def test_module_steps_flag_mates_at_contig_edges(store, tmp_path):
    (tmp_path / "kegg_modules.tsv").write_text(
        "# module\tname\tclass\tdefinition\n"
        "M90001\tedge module\tPathway modules; Test\tK90001 K90002 K90009\n"
        "M90002\tinterior module\tPathway modules; Test\tK90003 K90004 K90009\n"
        "M90003\tlone module\tPathway modules; Test\tK90005 K90009\n")
    store.conn.execute("""
        INSERT INTO annotations (annotation_id, protein_id, source, accession, name, evalue)
        VALUES (1, 'g8', 'kofam', 'K90001', 'a', 1e-50), (2, 'g9', 'kofam', 'K90002', 'a', 1e-50),
               (3, 'g4', 'kofam', 'K90003', 'b', 1e-50), (4, 'g5', 'kofam', 'K90004', 'b', 1e-50),
               (5, 'h2', 'kofam', 'K90005', 'c', 1e-50)""")
    rows = {r["module"]: r for r in genome_modules(store, directory=tmp_path)}
    assert rows["M90001"]["found_at_contig_edge"] == ["g8", "g9"]   # operon-like pair reaching the contig end
    assert rows["M90002"]["found_at_contig_edge"] == []             # interior pair
    assert rows["M90003"]["found_at_contig_edge"] == []             # a lone gene at an edge is not a cluster
    assert "could lie beyond the contig end" in modules_markdown(list(rows.values()), "t")


def test_fasta_and_header_parsing(tmp_path):
    fasta = tmp_path / "bin1.fna"
    fasta.write_text(">ctg1 desc\nACGT\nAC\n>ctg2\nA\n")
    assert fasta_lengths(fasta) == {"ctg1": 6, "ctg2": 1}
    assert parse_partial(["ID=1_1;partial=01;start_type=ATG"]) == "01"
    assert parse_partial(["ID=1_1;start_type=ATG"]) is None
    faa = tmp_path / "bin1.faa"
    faa.write_text(">ctg1_1 # 1 # 300 # 1 # ID=1_1;partial=10;start_type=Edge\nMAAA\n>ncbi_protein desc\nMA\n")
    assert list(prodigal_partials(faa)) == [("ctg1_1", "10")]


def _legacy_db(path: Path):
    conn = duckdb.connect(str(path))
    conn.execute("""
        CREATE TABLE schema_version (version INTEGER PRIMARY KEY, description VARCHAR,
                                     applied_at TIMESTAMP DEFAULT current_timestamp);
        CREATE TABLE contigs (contig_id VARCHAR PRIMARY KEY, bin_id VARCHAR, length INTEGER NOT NULL,
                              gc_content FLOAT, is_circular BOOLEAN DEFAULT FALSE, taxonomy VARCHAR);
        CREATE TABLE proteins (protein_id VARCHAR PRIMARY KEY, contig_id VARCHAR NOT NULL, bin_id VARCHAR,
                               start INTEGER NOT NULL, end_coord INTEGER NOT NULL, strand VARCHAR(1) NOT NULL,
                               gene_index INTEGER, sequence TEXT, sequence_length INTEGER, gc_content FLOAT);
        INSERT INTO schema_version (version, description) VALUES (7, 'legacy');
        INSERT INTO contigs (contig_id, bin_id, length) VALUES ('ctg1', 'bin1', 300), ('ctg9', 'bin1', 90);
        INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand)
        VALUES ('ctg1_1', 'ctg1', 'bin1', 1, 300, '+'), ('ctg9_1', 'ctg9', 'bin1', 1, 90, '+');
    """)
    return conn


def test_backfill_records_assembly_lengths_and_partials(tmp_path):
    conn = _legacy_db(tmp_path / "legacy.duckdb")
    fasta = tmp_path / "bin1.fna"
    fasta.write_text(">ctg1\n" + "A" * 1000 + "\n")
    faa = tmp_path / "bin1.faa"
    faa.write_text(">ctg1_1 # 1 # 300 # 1 # ID=1_1;partial=10\nMAAA\n")
    stats = backfill_contig_context(conn, assemblies={"bin1": fasta}, protein_faas=[faa], threads=2)
    assert stats == {"contigs_updated": 1, "contigs_unmatched": 1, "proteins_flagged": 1}
    assert conn.execute("SELECT contig_id, length, length_source FROM contigs ORDER BY 1").fetchall() == [
        ("ctg1", 1000, "assembly"), ("ctg9", 90, "gene_span")]
    assert conn.execute("SELECT partial FROM proteins WHERE protein_id = 'ctg1_1'").fetchone()[0] == "10"
    assert conn.execute("SELECT MAX(version) FROM schema_version").fetchone()[0] >= 8


def test_backfill_rejects_a_mismatched_assembly(tmp_path):
    conn = _legacy_db(tmp_path / "legacy.duckdb")
    fasta = tmp_path / "bin1.fna"
    fasta.write_text(">ctg1\nACGT\n")
    with pytest.raises(ValueError, match="wrong assembly"):
        backfill_contig_context(conn, assemblies={"bin1": fasta})
    assert conn.execute("SELECT length FROM contigs WHERE contig_id = 'ctg1'").fetchone()[0] == 300


def test_stage07_reads_assembly_lengths_and_partial_flags(tmp_path):
    module = SourceFileLoader("kb_build_edges", str(REPO_ROOT / "src/ingest/07_build_knowledge_base.py")).load_module()
    data = tmp_path / "data"
    assembly = data / "stage00_prepared" / "genomes" / "bin1.fna"
    assembly.parent.mkdir(parents=True)
    assembly.write_text(">bin1_contig1\n" + "A" * 5000 + "\n")
    (data / "stage00_prepared" / "processing_manifest.json").write_text(json.dumps(
        {"genomes": [{"genome_id": "bin1", "filename": "bin1.fna", "output_path": str(assembly)}]}))
    faa = data / "stage03_prodigal" / "genomes" / "bin1" / "bin1.faa"
    faa.parent.mkdir(parents=True)
    faa.write_text(">bin1_contig1_1 # 1 # 300 # 1 # ID=1_1;partial=10;start_type=Edge\nMAAA\n"
                   ">bin1_contig1_2 # 400 # 900 # 1 # ID=1_2;partial=00;start_type=ATG\nMAAA\n")
    outputs = module.PipelineOutputs(
        **{f"stage{s}_dir": data / name for s, name in [
            ("00", "stage00_prepared"), ("01", "stage01_quast"), ("02", "stage02_dfast_qc"),
            ("03", "stage03_prodigal"), ("04", "stage04_astra"), ("05a", "stage05a_gecco"),
            ("05b", "stage05b_dbcan"), ("05c", "stage05c_crispr"), ("06", "stage06_embeddings")]})
    module.KnowledgeBaseBuilder(outputs, data / "sharur.duckdb", force=True, threads=1).build()
    conn = duckdb.connect(str(data / "sharur.duckdb"), read_only=True)
    assert conn.execute("SELECT length, length_source FROM contigs").fetchall() == [(5000, "assembly")]
    assert conn.execute("SELECT protein_id, partial FROM proteins ORDER BY 1").fetchall() == [
        ("bin1_contig1_1", "10"), ("bin1_contig1_2", "00")]
