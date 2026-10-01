"""Taxonomy tree: clade counts, feature shares, feature search and the /tree page."""

import numpy as np
import pytest
from fastapi.testclient import TestClient

from sharur.browser import create_app
from sharur.browser import routes_tree as tree
from sharur.browser.catalog import Genome
from sharur.storage.duckdb_store import DuckDBStore


def _genome(i, family, genus):
    taxonomy = {"domain": "Archaea", "phylum": "P", "class": "C", "order": "O", "family": family, "genus": genus}
    return Genome(i, f"g{i}", taxonomy, 90.0, 1.0, 1, 1000, 1000, 10, 5, 100)


def test_build_tree_counts_carriers_per_clade():
    genomes = [_genome(0, "F1", "G1"), _genome(1, "F1", "G1"), _genome(2, "F1", "G2"), _genome(3, "F2", None)]
    mask = np.array([True, False, True, True])
    root = tree.build_tree(genomes, None, None, [mask])
    assert (root["n"], root["k"]) == (4, [3])
    shown, passed = tree.effective_root(root)
    # domain, phylum, class and order each hold one clade: the tree opens at the order
    assert (shown["rank"], shown["name"]) == ("order", "O") and len(passed) == 4
    families = {c["name"]: c for c in shown["children"]}
    assert (families["F1"]["n"], families["F1"]["k"]) == (3, [2])
    assert [c["name"] for c in families["F1"]["children"]] == ["G1", "G2"]
    # unclassified at genus: the genome stops under an Unclassified genus
    assert families["F2"]["children"][0]["name"] == "Unclassified" and families["F2"]["children"][0]["children"] == []


def test_build_tree_from_a_clade_root():
    genomes = [_genome(0, "F1", "G1"), _genome(1, "F1", "G2")]
    root = tree.build_tree(genomes, "family", "F1", [])
    assert (root["rank"], root["name"], root["n"]) == ("family", "F1", 2)
    assert [c["rank"] for c in root["children"]] == ["genus", "genus"]


def test_parse_features_keeps_order_and_known_kinds():
    parsed = tree.parse_features("ko:K00001, pfam:PF00005,bogus:x,system:defense:CBASS,ko:K00001")
    assert parsed == [("ko", "K00001"), ("pfam", "PF00005"), ("system", "defense:CBASS")]
    assert len(tree.parse_features(",".join(f"ko:K{i:05d}" for i in range(12)))) == tree.MAX_FEATURES


def test_search_features_ranks_exact_then_prefix_then_substring():
    options = [{"token": "pfam:PF00005", "kind": "pfam", "id": "PF00005", "label": "ABC_tran", "genomes": 9},
               {"token": "function:abc", "kind": "function", "id": "abc_transporter", "label": "ABC transporter", "genomes": 5},
               {"token": "ko:K1", "kind": "ko", "id": "K00001", "label": "adh · alcohol dehydrogenase", "genomes": None}]
    assert [o["id"] for o in tree.search_features(options, "abc")] == ["PF00005", "abc_transporter"]
    assert tree.search_features(options, "alcohol")[0]["id"] == "K00001"
    assert tree.search_features(options, "a") == []


@pytest.fixture
def client(tmp_path):
    path = tmp_path / "d.duckdb"
    store = DuckDBStore(str(path))
    for i, family in enumerate(["F1", "F1", "F2"]):
        b = f"g{i}"
        store.conn.execute("INSERT INTO bins (bin_id, taxonomy) VALUES (?, ?)",
                           [b, f"d__Archaea;p__P;c__C;o__O;f__{family};g__G{i};s__S{i}"])
        store.conn.execute("INSERT INTO contigs (contig_id, bin_id, length) VALUES (?, ?, 5000)", [f"{b}_c", b])
        store.conn.execute("INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, gene_index, "
                           "sequence_length) VALUES (?, ?, ?, 1, 900, '+', 0, 300)", [f"{b}_1", f"{b}_c", b])
        if i < 2:
            store.conn.execute("INSERT INTO annotations (annotation_id, protein_id, source, accession, name, evalue, "
                               "start_aa, end_aa) VALUES (?, ?, 'pfam', 'PF00005', 'ABC_tran', 1e-20, 5, 200)",
                               [i + 1, f"{b}_1"])
    store.close()
    return TestClient(create_app(str(path), background=False))


def test_tree_pages(client):
    page = client.get("/tree", params={"features": "pfam:PF00005"})
    assert page.status_code == 200 and 'id="tree-data"' in page.text and "ABC_tran" in page.text
    assert '"k": [2]' in page.text or '"k":[2]' in page.text
    assert client.get("/tree", params={"root": "family:F1", "layout": "rect"}).status_code == 200
    assert client.get("/tree", params={"root": "family:nope"}).status_code == 404
    assert client.get("/tree", params={"root": "bogus"}).status_code == 404
    missing = client.get("/tree", params={"features": "pfam:PF99999"})
    assert missing.status_code == 200 and "PF99999 is not in this dataset" in missing.text
    found = client.get("/api/tree/features", params={"q": "abc_tr"}).json()
    assert found and found[0]["token"] == "pfam:PF00005"
