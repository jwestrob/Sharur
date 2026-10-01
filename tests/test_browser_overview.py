"""Overview (home) page: quality tiers, the quality scatter, labels and graceful gaps."""

import threading
from types import SimpleNamespace

import pytest
from fastapi.testclient import TestClient

from sharur.browser import charts, create_app, overview
from sharur.browser.catalog import Catalog, Genome
from sharur.storage.duckdb_store import DuckDBStore


def _catalog(quality):
    catalog = Catalog()
    phyla = ["A", "A", "A", "B", "B", "C", "D"]
    for i, ((comp, cont), phylum) in enumerate(zip(quality, phyla)):
        g = Genome(i, f"g{i}", {"phylum": phylum}, comp, cont, 1, 1_000_000, 1000, 900, 600, 300)
        catalog.genomes.append(g)
        catalog.by_bin[g.bin_id] = g
    return catalog


def test_quality_tiers_and_scatter_groups():
    catalog = _catalog([(95, 1), (92, 6), (60, 2), (45, 1), (99, 0.5), (80, 9), (None, None)])
    q = overview.quality(catalog)
    assert (q["scored"], q["total"]) == (6, 7)
    # high: >90% and <5%; medium: >=50% and <10% (g1 at 6% contamination, g2, g5); low: g3 at 45% complete
    assert dict(q["tiers"]) == {"high": 2, "medium": 3, "low": 1}
    assert q["legend"] == [("A", 3), ("B", 2), ("C", 1)] and q["other"] == 0
    assert {p["group"] for p in q["points"]} == {0, 1, 2}
    svg = charts.quality_scatter(q["points"])
    assert svg.count("<circle") == 6 and 'data-g="g0"' in svg and "q-hq" in svg


def test_quality_without_completeness():
    q = overview.quality(_catalog([(None, None)] * 7))
    assert q["scored"] == 0 and q["points"] == [] and q["median"] is None
    assert charts.quality_scatter([]) == ""


def test_histogram_decimals():
    assert "median 1.18" in charts.histogram([1.18, 1.18, 1.2], "Mb", decimals=2)
    assert "median 921" in charts.histogram([921.0, 921.0, 930.0], "proteins")


def test_highlights_follow_discover_feeds():
    catalog = _catalog([(95, 1)] * 7)
    assert overview.highlights(catalog) == []
    catalog.notable = {"giants": [{"protein_id": "p", "bin_id": "g0", "length": 9000, "domains": [],
                                   "architecture": "A x3", "partial": True}],
                       "fusions": [], "islands": [], "giant_clades": {"rows": []}}
    (h,) = overview.highlights(catalog)
    assert h["title"] == "Longest protein" and h["value"] == "9,000 aa" and h["href"] == "/discover/giants"


def test_length_basis_label():
    ctx = SimpleNamespace(lock=threading.Lock(), store=None)
    assert overview.length_basis(ctx) == "assembly"     # unreadable store falls back to "assembled"


@pytest.fixture
def db(tmp_path):
    path = tmp_path / "sharur.duckdb"
    store = DuckDBStore(str(path))
    store.conn.execute("""
        INSERT INTO bins (bin_id, completeness, contamination, taxonomy)
            VALUES ('b1', 95.0, 1.0, 'd__Archaea;p__P1'), ('b2', 70.0, 3.0, 'd__Archaea;p__P2');
        INSERT INTO contigs (contig_id, bin_id, length, length_source) VALUES ('b1_c1', 'b1', 5000, 'gene_span'),
                                                                              ('b2_c1', 'b2', 5000, 'gene_span');
        INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, gene_index, sequence_length)
            VALUES ('b1_c1_1', 'b1_c1', 'b1', 1, 900, '+', 0, 300), ('b2_c1_1', 'b2_c1', 'b2', 1, 900, '+', 0, 300);
        INSERT INTO annotations (annotation_id, protein_id, source, accession, name, evalue, start_aa, end_aa)
            VALUES (1, 'b1_c1_1', 'pfam', 'PF00005', 'ABC_tran', 1e-30, 5, 200);
    """)
    store.close()
    return path


def test_home_page(db):
    client = TestClient(create_app(db, background=False))
    page = client.get("/", headers={"Accept-Encoding": "gzip"})
    assert page.status_code == 200 and page.headers.get("content-encoding") == "gzip"
    text = page.text
    assert "Genome quality" in text and 'data-g="b1"' in text and "high-quality genomes" in text
    assert "sequence spanned by gene calls" in text          # contig lengths stored as gene spans
    for href in ("/tree", "/landscape", "/health", "/matrix", "/discover"):
        assert f'href="{href}"' in text
    assert "50%" in text                                      # pfam: 1 of 2 genomes


def test_home_without_completeness(tmp_path):
    path = tmp_path / "sharur.duckdb"
    store = DuckDBStore(str(path))
    store.conn.execute("INSERT INTO bins (bin_id, taxonomy) VALUES ('b1', 'd__Bacteria;p__P1')")
    store.close()
    text = TestClient(create_app(path, background=False)).get("/").text
    assert "No completeness estimates yet" in text and 'class="qscatter"' not in text
