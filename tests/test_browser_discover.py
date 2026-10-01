"""Discovery feeds: fusions, unannotated stretches, giants, repeats, rare systems."""

import pytest
from fastapi.testclient import TestClient

from sharur.browser import create_app, discover
from sharur.browser.catalog import Catalog, Genome
from sharur.storage.duckdb_store import DuckDBStore

FAMILIES = ["F1"] * 4 + ["F2"] * 24


@pytest.fixture
def store(tmp_path):
    """28 genomes in two families. Domains A and B are common everywhere, fused only in three F1 genomes;
    g0 carries an unannotated stretch and a fully unannotated contig, and a giant with an A repeat."""
    s = DuckDBStore(str(tmp_path / "d.duckdb"))
    rows, ann = [], []
    aid = 0
    for i, fam in enumerate(FAMILIES):
        b = f"g{i}"
        s.conn.execute("INSERT INTO bins (bin_id, taxonomy) VALUES (?, ?)",
                       [b, f"d__Archaea;p__P;c__C;o__O;f__{fam};g__G{i};s__S{i}"])
        s.conn.execute("INSERT INTO contigs (contig_id, bin_id, length) VALUES (?, ?, 100000)", [f"{b}_c1", b])
        for k in range(10):
            pid = f"{b}_c1_{k + 1}"
            length = 9000 if (i == 0 and k == 0) else 300
            rows.append((pid, f"{b}_c1", b, 1 + k * 1000, 900 + k * 1000, "+", k, length, "10" if length > 5000 else "00"))
        # every genome: a separate A protein and B protein (genes 9 and 10)
        for pid, acc, name in ((f"{b}_c1_9", "PF00001", "A"), (f"{b}_c1_10", "PF00002", "B")):
            aid += 1
            ann.append((aid, pid, "pfam", acc, name, 1e-20, 10, 200))
        if fam == "F1" and i < 3:   # the fusion: A and B on gene 1
            for acc, name, lo in (("PF00001", "A", 10), ("PF00002", "B", 220)):
                aid += 1
                ann.append((aid, f"{b}_c1_1", "pfam", acc, name, 1e-20, lo, lo + 180))
        if i != 0:                   # genes 2-8 annotated outside g0
            for k in range(2, 9):
                aid += 1
                ann.append((aid, f"{b}_c1_{k}", "pfam", "PF00003", "C", 1e-10, 5, 100))
    # g0: a contig whose six genes are all unannotated
    s.conn.execute("INSERT INTO contigs (contig_id, bin_id, length) VALUES ('g0_c2', 'g0', 9000)")
    for k in range(6):
        rows.append((f"g0_c2_{k + 1}", "g0_c2", "g0", 1 + k * 1000, 900 + k * 1000, "-", k, 300, "00"))
    # g0's giant: eight A repeats in a row
    for r in range(8):
        aid += 1
        ann.append((aid, "g0_c1_1", "pfam", "PF00001", "A", 1e-5, 2000 + r * 100, 2090 + r * 100))
    s.conn.executemany("INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, gene_index, "
                       "sequence_length, partial) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)", rows)
    s.conn.executemany("INSERT INTO annotations (annotation_id, protein_id, source, accession, name, evalue, "
                       "start_aa, end_aa) VALUES (?, ?, ?, ?, ?, ?, ?, ?)", ann)
    yield s
    s.close()


def _catalog():
    catalog = Catalog()
    for i, fam in enumerate(FAMILIES):
        g = Genome(i, f"g{i}", {"domain": "Archaea", "phylum": "P", "class": "C", "order": "O", "family": fam,
                                "genus": f"G{i}"}, 90.0, 1.0, 1, 10000, 1000, 10, 9, 300)
        catalog.genomes.append(g)
        catalog.by_bin[g.bin_id] = g
    return catalog


def test_fusion_restricted_to_one_family(store):
    fused, total = discover.fusions(store, _catalog(), outside=10)
    assert total == 1
    f = fused[0]
    assert {name for _, name in f["domains"]} == {"A", "B"}
    assert (f["rank"], f["clade"], f["genomes"], f["clade_size"]) == ("family", "F1", 3, 4)
    # common domains need carriers outside the clade
    assert discover.fusions(store, _catalog(), outside=100)[1] == 0


def test_unannotated_stretches(store):
    found = {i["contig_id"]: i for i in discover.islands(store, min_genes=6)}
    assert found["g0_c1"]["genes"] == 7 and not found["g0_c1"]["at_edge"]   # genes 2-8
    whole = found["g0_c2"]
    assert whole["whole_contig"] and whole["genes"] == 6 and whole["strands"] == "------"
    assert len(found) == 2


def test_giants_and_repeats(store):
    giant = discover.giants(store)[0]
    assert giant["protein_id"] == "g0_c1_1" and giant["length"] == 9000 and giant["partial"]
    rep = discover.repeats(store)[0]
    assert (rep["protein_id"], rep["domain"], rep["run"], rep["lo"], rep["hi"]) == ("g0_c1_1", "A", 8, 2000, 2790)


def test_rare_systems():
    catalog = _catalog()
    catalog.systems = [{"kind": "defense", "type": "Rare", "bin_id": "g1"},
                       *({"kind": "defense", "type": "Common", "bin_id": f"g{i}"} for i in range(5))]
    rows = discover.rare_systems(catalog, lambda: [{"confident": True, "prediction": "I-B", "bin_id": "g2"},
                                                   {"confident": False, "prediction": "False", "bin_id": "g3"}])
    assert [(r["kind"], r["type"]) for r in rows] == [("crispr", "I-B"), ("defense", "Rare")]


def test_discover_pages(store, tmp_path):
    store.close()
    client = TestClient(create_app(str(tmp_path / "d.duckdb"), background=False))
    page = client.get("/discover")
    assert page.status_code == 200 and "Giant proteins" in page.text and "Unannotated stretches" in page.text
    for feed in ("giants", "clades", "fusions", "islands", "dark", "repeats", "systems"):
        assert client.get(f"/discover/{feed}").status_code == 200, feed
    assert client.get("/discover/nope").status_code == 404
    assert client.get("/discover/random", follow_redirects=False).status_code == 303
