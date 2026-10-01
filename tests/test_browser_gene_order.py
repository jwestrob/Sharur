"""Gene order between two genomes: matching, collinear blocks, drawing, and the Compare page."""

import pytest
from fastapi.testclient import TestClient

from sharur.browser import create_app
from sharur.browser import gene_order as go
from sharur.predicates_v2.persistence import generate_and_persist_v2
from sharur.storage.duckdb_store import DuckDBStore


def _genome(bin_id, contig, families, strands=None, reverse=False):
    genes = []
    order = list(range(len(families)))
    for pos, i in enumerate(reversed(order) if reverse else order):
        genes.append({"protein_id": f"{bin_id}_{i}", "contig": contig, "start": 1 + pos * 1000, "end": 900 + pos * 1000,
                      "strand": (strands or [1] * len(families))[i], "family": families[i], "pos": pos})
    genes.sort(key=lambda g: g["pos"])
    return {"bin_id": bin_id, "genes": genes, "by_contig": {contig: genes}, "contigs": [(contig, len(genes) * 1000)],
            "offset": {contig: 0}, "total": len(genes) * 1000}


def test_forward_block_with_a_gap_and_multicopy_cap():
    a = _genome("a", "ca", ["K1", "K2", "x", "K3", "K4", "M", "M", "M", "M", "M"])
    b = _genome("b", "cb", ["K1", "K2", "K3", "y", "K4", "M", "M", "M", "M", "M"])
    m = go.match(a, b)
    assert m["shared"] == 5 and m["skipped_families"] == 1          # 5 x 5 copies of M exceed the cap
    blocks = go.blocks(m["pairs"])
    assert len(blocks) == 1 and blocks[0]["size"] == 4 and not blocks[0]["reverse"]


def test_inverted_block():
    a = _genome("a", "ca", ["K1", "K2", "K3", "K4"])
    b = _genome("b", "cb", ["K1", "K2", "K3", "K4"], reverse=True)
    blocks = go.blocks(go.match(a, b)["pairs"])
    assert blocks[0]["size"] == 4 and blocks[0]["reverse"]


def test_compare_order_summary_and_drawings():
    a = _genome("a", "ca", ["K1", "K2", "K3", None, "K9"], strands=[1, 1, -1, 1, 1])
    b = _genome("b", "cb", ["K1", "K2", "K3", "K7"])
    r = go.compare_order(a, b)
    assert r["longest"] == 3 and r["in_blocks"] == 3 and r["share_a"] == pytest.approx(0.6)
    svg = go.dotplot(a, b, r, "A", "B")
    assert svg.count('class="dp-dot same') == 2 and svg.count('class="dp-dot flip') == 1   # K3 sits on opposite strands
    assert 'data-a="a_0"' in svg
    rib = go.ribbon(a, b, r["blocks"][0])
    assert rib.count('class="rb-link"') == 3 and 'href="/protein/a_1"' in rib


@pytest.fixture
def client(tmp_path):
    store = DuckDBStore(str(tmp_path / "d.duckdb"))
    for g in ("g1", "g2"):
        store.conn.execute("INSERT INTO bins (bin_id, taxonomy) VALUES (?, 'd__Archaea;p__P')", [g])
        store.conn.execute("INSERT INTO contigs (contig_id, bin_id, length) VALUES (?, ?, 6000)", [f"{g}_c", g])
        for i in range(5):
            pid = f"{g}_c_{i}"
            store.conn.execute("INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, gene_index, "
                               "sequence_length) VALUES (?, ?, ?, ?, ?, '+', ?, 300)", [pid, f"{g}_c", g, 1 + i * 1000, 900 + i * 1000, i])
            store.conn.execute("INSERT INTO annotations (annotation_id, protein_id, source, accession, name, evalue) "
                               "VALUES (?, ?, 'kofam', ?, ?, 1e-30)", [hash(pid) % 10**9, pid, f"K0000{i}", f"K0000{i}"])
    generate_and_persist_v2(store, chunk_size=10, return_states=False, update_legacy_predicates=True)
    store.close()
    return TestClient(create_app(str(tmp_path / "d.duckdb"), background=False))


def test_compare_page_shows_gene_order_for_two_genomes(client):
    page = client.get("/compare", params={"a": "g1", "b": "g2"}).text
    assert 'id="gene-order"' in page and "svg class=\"dotplot\"" in page and "5 genes" in page
    clades = client.get("/compare", params={"a": "phylum:P", "b": "g2"}).text
    assert 'id="gene-order"' not in clades
