"""Protein-page context: same architecture elsewhere, copies in the genome, usual neighbours."""

import pytest
from fastapi.testclient import TestClient

from sharur.browser import create_app
from sharur.browser import protein_context as pc
from sharur.browser.routes_protein_context import exact_pattern
from sharur.storage.duckdb_store import DuckDBStore

A, B, C = ("PF00001", "A"), ("PF00002", "B"), ("PF00003", "C")


@pytest.fixture
def db(tmp_path):
    """Six genomes, six genes each. Gene 2 is A-B with KO K00001 (B-A in g5); g0 has an adjacent copy
    at gene 3; domain C sits two genes downstream of each carrier."""
    path = tmp_path / "d.duckdb"
    s = DuckDBStore(str(path))
    proteins, ann = [], []

    def hit(pid, dom, lo, ev=1e-20):
        ann.append((len(ann) + 1, pid, "pfam", dom[0], dom[1], ev, lo, lo + 90))

    for i in range(6):
        b = f"g{i}"
        s.conn.execute("INSERT INTO bins (bin_id, taxonomy) VALUES (?, ?)",
                       [b, f"d__Archaea;p__P;c__C;o__O;f__F{i % 2};g__G{i};s__S{i}"])
        s.conn.execute("INSERT INTO contigs (contig_id, bin_id, length) VALUES (?, ?, 7000)", [f"{b}_c1", b])
        for k in range(1, 7):
            proteins.append((f"{b}_c1_{k}", f"{b}_c1", b, k * 1000, k * 1000 + 600, "+", k - 1, 200, "00"))
        carriers = [2, 3] if i == 0 else [2]
        for k in carriers:
            pid = f"{b}_c1_{k}"
            first, second = (B, A) if i == 5 else (A, B)
            hit(pid, first, 5)
            hit(pid, second, 105)
            ann.append((len(ann) + 1, pid, "kofam", "K00001", "K00001", 1e-30, None, None))
        hit(f"{b}_c1_{max(carriers) + 2}", C, 5)
    s.conn.executemany("INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, gene_index, "
                       "sequence_length, partial) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)", proteins)
    s.conn.executemany("INSERT INTO annotations (annotation_id, protein_id, source, accession, name, evalue, "
                       "start_aa, end_aa) VALUES (?, ?, ?, ?, ?, ?, ?, ?)", ann)
    s.close()
    return path


@pytest.fixture
def store(db):
    s = DuckDBStore(str(db), read_only=True)
    yield s
    s.close()


def _domains(store, pid):
    from sharur.architecture import architecture

    return [d.to_dict() for d in architecture(store, pid)]


def test_same_architecture_counts_exact_order(store):
    from sharur.browser.catalog import load_catalog

    catalog = load_catalog(store)
    out = pc.same_architecture(store, catalog, "g1_c1_2", _domains(store, "g1_c1_2"), 200, ko="K00001")
    # A-B in g0 (two copies) and g1-g4; g5 carries B-A. Genomes count the other proteins' genomes.
    assert (out["kind"], out["label"], out["others"], out["genomes"], out["same_domains"]) == \
        ("architecture", "A - B", 5, 4, 7)
    assert out["estimated"] is None and out["lengths"]["median"] == 200


def test_paralogs_flag_adjacent_copies(store):
    rows = pc.paralogs(store, "g0_c1_2", "g0", ko="K00001", domains=_domains(store, "g0_c1_2"))
    assert [(r["protein_id"], r["tandem"]) for r in rows] == [("g0_c1_3", "adjacent")]
    assert rows[0]["why"] == ["same KO (K00001)", "same Pfam domains"]
    assert pc.paralogs(store, "g1_c1_2", "g1", ko="K00001", domains=_domains(store, "g1_c1_2")) == []


def test_neighbourhood_ranks_downstream_family(store):
    family = pc.family_of(store, "g1_c1_2", _domains(store, "g1_c1_2"))
    assert family == {"kind": "ko", "id": "K00001", "name": "K00001", "accession": "K00001"}
    hood = pc.neighbourhood(store, family)
    assert (hood["carriers"], hood["sampled"], hood["vocabulary"]) == (7, 6, "Pfam domain")
    top = hood["rows"][0]
    assert (top["name"], top["share"], top["side"]) == ("C", 1.0, "downstream")


def test_contig_position_and_pattern(store):
    pos = pc.contig_position(store, "g1_c1_2", "g1_c1")
    assert (pos["index"], pos["count"]) == (2, 6)
    assert pc.contig_position(store, "x", "x") is None
    assert exact_pattern([{"name": "A", "accession": "PF1"}, {"name": "A", "accession": "PF1"},
                          {"name": "odd name", "accession": "PF00009.3"}]) == "^ A{2} PF00009 $"


def test_context_fragment_and_ko_stack(db):
    client = TestClient(create_app(db, background=False))
    page = client.get("/protein/g0_c1_2")
    assert page.status_code == 200 and "/protein/g0_c1_2/context" in page.text
    frag = client.get("/protein/g0_c1_2/context")
    assert frag.status_code == 200
    for text in ("Same architecture elsewhere", "adjacent", "Usual neighbours", "Position on contig"):
        assert text in frag.text, text
    assert client.get("/protein/nope/context").status_code == 404
    assert client.get("/stack/ko/K00001").status_code == 200
