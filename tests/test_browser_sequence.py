"""Protein sequence properties: hydropathy segments, low complexity, N-terminus, self-similarity."""

import base64
import random

from fastapi.testclient import TestClient

from sharur.browser import create_app
from sharur.browser import sequence_properties as sp
from sharur.storage.duckdb_store import DuckDBStore

AMINO = "ACDEFGHIKLMNPQRSTVWY"


def _random(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice(AMINO) for _ in range(n))


def test_hydrophobic_segment_and_n_terminus():
    seq = "MKR" + "LLIVALLIVALLIVALLIVAL" + "DEKRDEKRSG" * 6 + "LIVLAVLLIFVALLIVLAG" + "DEKRSTNQ" * 5
    a = sp.analyse(seq)
    assert len(a["segments"]) == 2
    first = a["segments"][0]
    assert first["start"] <= 4 and first["end"] >= 24
    assert a["n_terminal"] and a["n_terminal"]["start"] <= 10


def test_low_complexity_found_only_where_entropy_drops():
    seq = _random(80, 1) + "Q" * 40 + _random(80, 2)
    regions = sp.low_complexity(seq)
    assert len(regions) == 1 and regions[0]["top"] == "Q"
    assert regions[0]["start"] <= 82 and regions[0]["end"] >= 119
    assert sp.low_complexity(_random(300, 3)) == []


def test_tandem_repeat_period_and_composition():
    unit = _random(50, 4)
    seq = _random(100, 5) + unit * 8 + _random(100, 6)
    s = sp.self_similarity(seq)
    assert s["period"]["period"] == 50
    assert s["period"]["start"] <= 110 and s["period"]["end"] >= 480
    assert s["cells"].shape == (len(seq), len(seq)) or s["bins"] <= sp.GRID
    assert sp.self_similarity(_random(600, 7))["period"] is None
    comp = sp.composition("KKKRDDE" + "A" * 3)
    assert comp["net_charge"] == 1 and comp["length"] == 10


def test_large_proteins_bin_onto_a_bounded_grid():
    seq = (_random(400, 8) * 30)[:12000]
    s = sp.self_similarity(seq)
    assert s["bins"] == sp.GRID and 4 <= s["k"] <= 8 and s["period"]["period"] == 400
    assert sp.word_length("ACDEFGHIKLMNPQRSTVWY" * 5, 100) == 4               # one residue per cell: the floor
    assert sp.word_length("ACDEFGHIKLMNPQRSTVWY" * 4000, 600) > 4             # 133 residues per cell
    png = base64.b64decode(sp.dotplot_png(s["cells"]))
    assert png.startswith(b"\x89PNG\r\n\x1a\n")


def test_properties_fragment(tmp_path):
    path = tmp_path / "sharur.duckdb"
    store = DuckDBStore(str(path))
    seq = "MKR" + "LLIVALLIVALLIVALLIVAL" + _random(200, 9)
    store.conn.execute("INSERT INTO bins (bin_id, taxonomy) VALUES ('b', 'd__Archaea')")
    store.conn.execute("INSERT INTO contigs (contig_id, bin_id, length) VALUES ('b_c1', 'b', 5000)")
    store.conn.execute("INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, gene_index, "
                       "sequence, sequence_length) VALUES ('p1', 'b_c1', 'b', 1, 700, '+', 0, ?, ?)", [seq, len(seq)])
    store.close()
    client = TestClient(create_app(path, background=False))
    page = client.get("/protein/p1/properties")
    assert page.status_code == 200
    assert "Sequence properties" in page.text and "Self-similarity" in page.text and "seqp-fig" in page.text
    assert seq not in page.text                      # the fragment carries drawings and numbers only
    assert 'id="sequence-properties-slot"' in client.get("/protein/p1").text
    assert client.get("/protein/nope/properties").status_code == 404
