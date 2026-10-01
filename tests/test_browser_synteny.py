"""Synteny in the browser: present with a sidecar, absent without one."""

import json
import time

from fastapi.testclient import TestClient

from sharur.browser import create_app
from sharur.storage.duckdb_store import DuckDBStore
from tests.test_synteny import _build_synteny_sidecar

GENOMES = {
    "genome_a": ("contig_a", "d__Archaea;p__P1;c__C1;o__O1;f__F1;g__G1;s__S1"),
    "genome_b": ("contig_b", "d__Archaea;p__P2;c__C2;o__O2;f__F2;g__G2;s__S2"),
    "genome_c": ("contig_c", "d__Archaea;p__P2;c__C2;o__O2;f__F3;g__G3;s__S3"),
}
PROTEINS = [("p1", "genome_a", 1, 270, "+"), ("p10", "genome_a", 301, 570, "+"),
            ("q1", "genome_b", 301, 570, "-"), ("q2", "genome_b", 601, 870, "-")]


def _core(directory, proteins=PROTEINS):
    store = DuckDBStore(str(directory / "sharur.duckdb"))
    for bin_id, (contig, taxonomy) in GENOMES.items():
        store.conn.execute("INSERT INTO bins (bin_id, taxonomy) VALUES (?, ?)", [bin_id, taxonomy])
        store.conn.execute("INSERT INTO contigs (contig_id, bin_id, length) VALUES (?, ?, 5000)", [contig, bin_id])
    for i, (pid, bin_id, start, end, strand) in enumerate(proteins):
        store.conn.execute(
            "INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, gene_index, "
            "sequence_length) VALUES (?, ?, ?, ?, ?, ?, ?, 90)", [pid, GENOMES[bin_id][0], bin_id, start, end, strand, i])
    store.conn.execute("INSERT INTO annotations (annotation_id, protein_id, source, accession, name, evalue, start_aa, "
                       "end_aa) VALUES (1, 'p1', 'pfam', 'PF00001', 'DomA', 1e-20, 1, 80)")
    store.close()
    return directory / "sharur.duckdb"


def _seal(directory, dataset_id):
    (directory / "dataset.seal.json").write_text(json.dumps({"dataset_id": dataset_id}))


def _client(db):
    app = create_app(db, background=False)
    view = getattr(app.state.ctx, "synteny", None)
    for _ in range(100):           # summaries load in a background thread
        if view is None or view.summary is not None or view.summary_error:
            break
        time.sleep(0.02)
    return TestClient(app)


def test_without_a_sidecar_nothing_about_synteny_renders(tmp_path):
    client = _client(_core(tmp_path))
    for path in ("/synteny", "/synteny/cluster/7", "/synteny/protein/p1", "/discover/synteny"):
        assert client.get(path).status_code == 404, path
    for path in ("/", "/protein/p1", "/discover", "/genome/genome_a"):
        page = client.get(path)
        assert page.status_code == 200 and "ynteny" not in page.text, path


def test_synteny_pages_with_a_sidecar(tmp_path):
    db = _core(tmp_path)
    _build_synteny_sidecar(tmp_path / "synteny.duckdb")
    client = _client(db)
    overview = client.get("/synteny")
    assert overview.status_code == 200 and "/synteny/cluster/7" in overview.text
    assert "earlier version" not in overview.text            # no seal: nothing to compare
    cluster = client.get("/synteny/cluster/7")
    assert cluster.status_code == 200 and cluster.text.count('class="stack-item"') == 2
    assert client.get("/synteny/cluster/cluster%3A7").status_code == 200
    assert client.get("/synteny/cluster/99").status_code == 404
    panel = client.get("/synteny/protein/p1").text
    assert "Cluster 7" in panel and "Pair 8" in panel and "anchor" in panel and "context" in panel
    assert 'data-synteny-url="/synteny/protein/p1"' in client.get("/protein/p1").text
    found = client.get("/synteny", params={"q": "p1"}, follow_redirects=False)
    assert found.status_code == 303 and found.headers["location"].endswith("/protein/p1#synteny")
    assert "/synteny/cluster/7" in client.get("/synteny", params={"q": "genome_b"}).text
    for page in (overview.text, cluster.text, panel, client.get("/discover").text):
        assert "elsa" not in page.lower()


def test_sidecar_from_an_earlier_dataset_version_shows_a_note(tmp_path):
    db = _core(tmp_path)
    _build_synteny_sidecar(tmp_path / "synteny.duckdb", dataset_id="dataset-old")
    _seal(tmp_path, "dataset-new")
    client = _client(db)
    assert "earlier version of this dataset" in client.get("/synteny").text


def test_sidecar_whose_genes_no_longer_resolve_stays_hidden(tmp_path):
    db = _core(tmp_path, proteins=PROTEINS[:2])                # q1, q2 are gone
    _build_synteny_sidecar(tmp_path / "synteny.duckdb", dataset_id="dataset-old")
    _seal(tmp_path, "dataset-new")
    client = _client(db)
    assert client.get("/synteny").status_code == 404
    assert "ynteny" not in client.get("/protein/p1").text
