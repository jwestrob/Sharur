"""Read-only dataset browser."""

import pytest
from fastapi.testclient import TestClient

from sharur.browser import create_app
from sharur.predicates_v2.persistence import generate_and_persist_v2
from sharur.storage.duckdb_store import DuckDBStore


SEQ = "MKTAYIAKQRQISFVKSHFSRQ" * 10


@pytest.fixture
def db(tmp_path):
    path = tmp_path / "sharur.duckdb"
    store = DuckDBStore(str(path))
    store.conn.execute(f"""
        INSERT INTO bins (bin_id, completeness, contamination, taxonomy) VALUES ('bin|1', 91.5, 1.2, 'd__Archaea');
        INSERT INTO contigs (contig_id, bin_id, length, length_source) VALUES ('bin|1_c1', 'bin|1', 9000, 'assembly');
        INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, gene_index,
                              sequence, sequence_length, partial)
        VALUES ('bin|1_c1_1', 'bin|1_c1', 'bin|1', 1, 600, '+', 0, '{SEQ}', 220, '10'),
               ('bin|1_c1_2', 'bin|1_c1', 'bin|1', 700, 1400, '-', 1, '{SEQ}', 220, '00'),
               ('bin|1_c1_3', 'bin|1_c1', 'bin|1', 1500, 2100, '+', 2, '{SEQ}', 220, '00');
        INSERT INTO annotations (annotation_id, protein_id, source, accession, name, description, evalue, score,
                                 start_aa, end_aa)
        VALUES (1, 'bin|1_c1_2', 'pfam', 'PF00005', 'ABC_tran', 'ABC <transporter>', 1e-40, 150.0, 20, 160);
    """)
    generate_and_persist_v2(store, chunk_size=1, return_states=False, update_legacy_predicates=True)
    store.close()
    return path


@pytest.fixture
def client(db):
    return TestClient(create_app(db))


def test_home_and_navigation(client):
    home = client.get("/")
    assert home.status_code == 200 and "3 proteins in 1 genomes" in home.text
    protein = client.get("/go", params={"q": "bin|1_c1_2"}, follow_redirects=False)
    assert protein.status_code == 303 and protein.headers["location"] == "/protein/bin%7C1_c1_2"
    assert client.get("/go", params={"q": "bin|1"}, follow_redirects=False).headers["location"] == "/genome/bin%7C1"
    assert client.get("/go", params={"q": "ABC_tran"}, follow_redirects=False).headers["location"].startswith(
        "/architecture?pattern=")


def test_protein_page_escapes_and_omits_sequences(client):
    page = client.get("/protein/bin%7C1_c1_2")
    assert page.status_code == 200
    assert "ABC &lt;transporter&gt;" in page.text and "<transporter>" not in page.text
    assert "MKTAYIAKQ" not in page.text
    assert 'aria-label="Domain architecture"' in page.text and 'aria-label="Gene neighborhood"' in page.text
    assert "/protein/bin%7C1_c1_2/why/abc_transporter" in page.text
    assert client.get("/protein/missing").status_code == 404


def test_why_genome_predicate_and_pattern_pages(client):
    why = client.get("/protein/bin%7C1_c1_2/why/abc_transporter")
    assert why.status_code == 200 and "PF00005" in why.text
    genome = client.get("/genome/bin%7C1")
    assert genome.status_code == 200 and "3 proteins" in genome.text and "91.5" in genome.text
    predicate = client.get("/predicate/abc_transporter")
    assert "bin|1_c1_2" in predicate.text
    pattern = client.get("/architecture", params={"pattern": "ABC_tran"})
    assert "1 proteins match" in pattern.text
    assert "Unbalanced" in client.get("/architecture", params={"pattern": "( ABC_tran"}).text


def test_token_is_exchanged_for_a_cookie(db):
    client = TestClient(create_app(db, token="s3cret"))
    assert client.get("/").status_code == 401
    exchanged = client.get("/?token=s3cret", follow_redirects=False)
    assert exchanged.status_code == 303 and "token" not in exchanged.headers["location"]
    assert client.get("/").status_code == 200
