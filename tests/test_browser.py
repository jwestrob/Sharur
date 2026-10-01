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
    return TestClient(create_app(db, background=False))


def test_overview_and_search(client):
    home = client.get("/")
    assert home.status_code == 200 and "3 proteins" in home.text
    protein = client.get("/search", params={"q": "bin|1_c1_2"}, follow_redirects=False)
    assert protein.status_code == 303 and protein.headers["location"] == "/protein/bin%7C1_c1_2"
    assert client.get("/search", params={"q": "bin|1"}, follow_redirects=False).headers["location"] == "/genome/bin%7C1"
    suggestions = client.get("/api/suggest", params={"q": "Archaea"}).json()
    assert {"kind": "taxon", "label": "Archaea", "sub": "domain", "url": "/taxa/domain/Archaea"} in suggestions
    assert "Nothing matched" in client.get("/search", params={"q": "zzzz"}).text


def test_protein_page_escapes_html_and_offers_the_sequence(client):
    page = client.get("/protein/bin%7C1_c1_2")
    assert page.status_code == 200
    assert "ABC &lt;transporter&gt;" in page.text and "<transporter>" not in page.text
    # the sequence is shown in 10-residue blocks and held once, unbroken, for the copy button
    assert SEQ[:10] in page.text and f'id="seq-raw" class="visually-hidden" readonly aria-hidden="true" tabindex="-1">{SEQ}<' in page.text
    assert "Copy FASTA" in page.text and "/fasta/bin%7C1_c1_2" in page.text
    fasta = client.get("/fasta/bin%7C1_c1_2")
    lines = fasta.text.splitlines()
    assert fasta.status_code == 200 and lines[0] == ">bin|1_c1_2 genome=bin|1"
    assert "".join(lines[1:]) == SEQ and all(len(line) <= 60 for line in lines[1:])
    assert client.get("/fasta/missing").status_code == 404
    assert 'aria-label="Domain architecture"' in page.text and 'aria-label="Gene neighborhood"' in page.text
    assert "ABC transporter" in page.text  # functional label shown by name
    assert "/protein/bin%7C1_c1_2/why/abc_transporter" in page.text
    assert client.get("/protein/missing").status_code == 404


def test_why_genome_predicate_and_pattern_pages(client):
    why = client.get("/protein/bin%7C1_c1_2/why/abc_transporter")
    assert why.status_code == 200 and "PF00005" in why.text
    genome = client.get("/genome/bin%7C1")
    assert genome.status_code == 200 and "Contig landscape" in genome.text and "91.5" in genome.text
    assert client.get("/genome/bin%7C1/proteins").status_code == 200
    function = client.get("/function/abc_transporter")
    assert function.status_code == 200 and "bin|1_c1_2" in function.text
    for page in ("/taxa", "/taxa/domain/Archaea", "/genomes", "/functions", "/systems", "/discover", "/pathways",
                 "/domains", "/domains?q=abc&sort=rare", "/genome/bin%7C1/contigs"):
        assert client.get(page).status_code == 200, page
    pattern = client.get("/architecture", params={"pattern": "ABC_tran"})
    assert "1 proteins match" in pattern.text
    assert "Unbalanced" in client.get("/architecture", params={"pattern": "( ABC_tran"}).text


def test_token_is_exchanged_for_a_cookie(db):
    client = TestClient(create_app(db, token="s3cret", background=False))
    assert client.get("/").status_code == 401
    exchanged = client.get("/?token=s3cret", follow_redirects=False)
    assert exchanged.status_code == 303 and "token" not in exchanged.headers["location"]
    assert client.get("/").status_code == 200


def test_domain_and_contig_pages(client):
    domain = client.get("/domain/PF00005")
    assert domain.status_code == 200 and "ABC_tran" in domain.text and "Functional labels from this family" in domain.text
    assert "abc_transporter" in domain.text or "ABC transporter" in domain.text
    assert client.get("/domain/PF99999").status_code == 404
    assert client.get("/search", params={"q": "PF00005"}, follow_redirects=False).headers["location"] == "/domain/PF00005"
    contig = client.get("/contig/bin%7C1_c1")
    assert contig.status_code == 200 and 'aria-label="Contig genes"' in contig.text and "gene 2" in contig.text
    assert client.get("/contig/bin%7C1_c1", params={"start": 500, "span": 2000}).status_code == 200
    assert client.get("/contig/missing").status_code == 404
    genome = client.get("/genome/bin%7C1")
    assert "/contig/bin%7C1_c1" in genome.text  # contig landscape links to the viewer
