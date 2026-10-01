"""Functional landscape: PCA placement, closest genomes, genome selections as groups."""

import threading
from types import SimpleNamespace

import numpy as np
import pytest

from sharur.browser import routes_landscape as ls
from sharur.browser import routes_matrix as mx
from sharur.browser.catalog import Catalog, Genome
from sharur.browser.routes_compare import clade_filter, genome_set, resolve_side


def _catalog():
    """Twelve genomes in two classes; class A carries features f0-f9, class B f10-f19, all carry f20-f24."""
    catalog = Catalog()
    for i in range(12):
        cls = "A" if i < 6 else "B"
        g = Genome(i, f"g{i}", {"domain": "Archaea", "class": cls, "order": f"O{cls}"},
                   90.0 - i, 1.0, 1, 10000, 1000, 30, 20, 300)
        catalog.genomes.append(g)
        catalog.by_bin[g.bin_id] = g
    return catalog


def _features(catalog):
    rng = np.random.default_rng(1)
    pairs = []
    for g in catalog.genomes:
        own = range(10) if g.taxonomy["class"] == "A" else range(10, 20)
        for f in own:
            if rng.random() < 0.85:
                pairs.append((g.bin_id, f"f{f:02d}", 1.0))
        for f in range(20, 25):
            pairs.append((g.bin_id, f"f{f:02d}", 1.0))
    # one class-B genome shares most of class A's repertoire
    pairs += [("g11", f"f{f:02d}", 1.0) for f in range(10)]
    ids, bins, feats, values = mx._sparse(catalog, pairs)
    return mx.FeatureSet("ko", ids, ids, [f"name {i}" for i in ids], bins, feats, values)


def test_first_axis_separates_the_two_repertoires():
    catalog = _catalog()
    emb = ls.embed(_features(catalog), catalog)
    x = emb.coords[:, 0]
    a, b = x[:6], x[6:11]
    assert (a.max() < b.min()) or (b.max() < a.min())
    assert 0 < emb.explained[1] < emb.explained[0] < 1
    assert emb.used_features == 20   # f20-f24 are in every genome and carry no contrast
    assert {f["id"] for f in emb.loadings[0]["high"] + emb.loadings[0]["low"]} <= {f"f{i:02d}" for i in range(20)}


def test_closest_genomes_by_jaccard_and_cross_lineage_flag():
    catalog = _catalog()
    emb = ls.embed(_features(catalog), catalog)
    near = ls.neighbours(emb, catalog, 0, top=11)
    assert all(n["genome"].taxonomy["class"] == "A" for n in near[:5])
    g11 = next(n for n in near if n["genome"].bin_id == "g11")
    assert g11["differs_at"] == "class" and 0 < g11["similarity"] < 1
    assert near == sorted(near, key=lambda n: -n["similarity"])


def test_too_few_genomes_returns_a_note():
    catalog = _catalog()
    catalog.genomes, catalog.by_bin = catalog.genomes[:2], {g.bin_id: g for g in catalog.genomes[:2]}
    fs = mx.FeatureSet("ko", ["f00"], ["f00"], ["n"], np.array([0, 1], dtype=np.int32),
                       np.array([0, 0], dtype=np.int32), np.array([1.0, 1.0], dtype=np.float32))
    emb = ls.embed(fs, catalog)
    assert emb.note and np.isnan(emb.coords).all()


def test_genome_selection_token():
    catalog = _catalog()
    picked = genome_set(catalog, "genomes:g1,g2,nope,g1")
    assert [g.bin_id for g in picked] == ["g1", "g2"]
    assert genome_set(catalog, "class:A") is None
    side = resolve_side(catalog, "genomes:g1,g2", lambda *p: "/" + "/".join(p))
    assert side.kind == "selection" and len(side.genomes) == 2 and side.label == "2 selected genomes"
    assert resolve_side(catalog, "genomes:g3", lambda *p: "/" + "/".join(p)).kind == "genome"
    assert resolve_side(catalog, "genomes:nope", lambda *p: "") is None
    assert clade_filter(catalog, "genomes:g1,g7") == ("2 selected genomes", {"g1", "g7"})
    catalog.selections["selection:abc"] = ["g4", "g5", "gone"]
    assert [g.bin_id for g in genome_set(catalog, "selection:abc")] == ["g4", "g5"]
    assert genome_set(catalog, "selection:unknown") == []
    ctx = SimpleNamespace(catalog=catalog, lock=threading.Lock(), store=None)
    group = mx.resolve_group(ctx, "genomes:g0,g6")
    assert group.label == "2 selected genomes" and [g.bin_id for g in group.genomes] == ["g0", "g6"]


def test_landscapes_cache_and_errors():
    catalog = _catalog()
    fs = _features(catalog)
    calls = []

    class Sets:
        def get(self, kind):
            calls.append(kind)
            if kind == "function":
                raise LookupError("no labels")
            return fs

    maps = ls.Landscapes(SimpleNamespace(catalog=catalog, feature_sets=Sets()), background=False)
    first = maps.get("ko")
    assert first is maps.get("ko") and calls == ["ko"]
    with pytest.raises(LookupError):
        maps.get("function")
    data = ls.payload(first, catalog)
    assert len(data["ids"]) == 12 and data["default_rank"] in ("class", "order") and data["feature_url"]


def test_landscape_routes(tmp_path):
    from fastapi.testclient import TestClient

    from sharur.browser import create_app
    from sharur.storage.duckdb_store import DuckDBStore

    store = DuckDBStore(str(tmp_path / "d.duckdb"))
    for i in range(8):
        b = f"g{i}"
        store.conn.execute("INSERT INTO bins (bin_id, taxonomy, completeness) VALUES (?, ?, 90)",
                           [b, f"d__Archaea;p__P;c__C{i % 2};o__O{i % 2};f__F;g__G;s__S{i}"])
        store.conn.execute("INSERT INTO contigs (contig_id, bin_id, length) VALUES (?, ?, 5000)", [f"{b}_c", b])
        for k in range(4):
            pid = f"{b}_c_{k}"
            store.conn.execute("INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, "
                               "gene_index, sequence_length) VALUES (?, ?, ?, ?, ?, '+', ?, 300)",
                               [pid, f"{b}_c", b, 1 + k * 1000, 900 + k * 1000, k])
            ko = f"K{(i % 2) * 10 + k + (i % 3):05d}"
            store.conn.execute("INSERT INTO annotations (annotation_id, protein_id, source, accession, name, evalue) "
                               "VALUES (?, ?, 'kofam', ?, ?, 1e-20)", [i * 10 + k, pid, ko, ko])
    store.close()
    client = TestClient(create_app(str(tmp_path / "d.duckdb"), background=False))
    page = client.get("/landscape")
    assert page.status_code == 200 and "landscape.js" in page.text
    data = client.get("/api/landscape").json()
    assert data["status"] == "ready" and len(data["ids"]) == 8 and len(data["x"]) == 8
    assert client.get("/api/landscape", params={"kind": "nope"}).status_code == 404
    near = client.get("/api/neighbors/g0")
    assert near.status_code == 200 and "Jaccard" in near.text and "g2" in near.text
    assert client.get("/api/neighbors/nope").status_code == 404
    assert "Functionally closest genomes" in client.get("/genome/g0").text
    matrix = client.get("/matrix", params={"a": "genomes:g0,g1,g2", "kind": "ko", "features": "K00000,K00001"})
    assert matrix.status_code == 200 and "3 selected genomes" in matrix.text
    stored = client.post("/api/selection", json={"genomes": ["g3", "g1", "nope", "g1"]}).json()
    assert stored["genomes"] == 2 and stored["token"].startswith("selection:")
    assert stored == client.post("/api/selection", json={"genomes": ["g1", "g3"]}).json()   # same set, same key
    compare = client.get("/compare", params={"a": stored["token"], "b": "class:C0"})
    assert compare.status_code == 200 and "2 selected genomes" in compare.text
    assert client.post("/api/selection", json={"genomes": ["nope"]}).status_code == 400
