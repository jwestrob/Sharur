"""Per-call system pages: members in genomic context, overlays, links."""

import pytest
from fastapi.testclient import TestClient

from sharur.browser import create_app
from sharur.predicates_v2.persistence import generate_and_persist_v2
from sharur.storage.duckdb_store import DuckDBStore


@pytest.fixture
def client(tmp_path):
    path = tmp_path / "ds" / "sharur.duckdb"
    path.parent.mkdir()
    store = DuckDBStore(str(path))
    genes = ", ".join(f"('c_{i}', 'c', 'g1', {1 + 1000 * i}, {900 + 1000 * i}, '+', {i}, 300)" for i in range(10))
    store.conn.execute(f"""
        INSERT INTO bins (bin_id, taxonomy) VALUES ('g1', 'd__Bacteria;p__P;c__C');
        INSERT INTO contigs (contig_id, bin_id, length) VALUES ('c', 'g1', 10000);
        INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, gene_index, sequence_length)
        VALUES {genes};
        INSERT INTO annotations (annotation_id, protein_id, source, accession, name, description, evalue, start_aa, end_aa)
        VALUES (1, 'c_7', 'pfam', 'PF00005', 'ABC_tran', 'ABC transporter', 1e-40, 10, 150);
        INSERT INTO defense_systems (system_id, genome_id, contig_id, system_type, system_subtype, genes_count,
                                     protein_ids, profile_names)
        VALUES ('g1_Toy_1', 'g1', 'c', 'Toy', 'Toy_I', 2, 'c_4,c_5', 'Toy__ToyA,Toy__ToyB'),
               ('g1_Other_1', 'g1', 'c', 'Other', 'Other', 1, 'c_6', 'Other__OthA');
        INSERT INTO loci (locus_id, locus_type, contig_id, start, end_coord, confidence, metadata)
        VALUES ('g1_prophage_1', 'prophage', 'c', 6000, 9000, 1.0, '{{}}');
    """)
    generate_and_persist_v2(store, chunk_size=10, return_states=False, update_legacy_predicates=True)
    store.close()
    return TestClient(create_app(path, background=False, notes_path=tmp_path / "notes.sqlite"))


def test_call_page_shows_members_in_context(client):
    page = client.get("/call/g1_Toy_1", params={"flank": 2000})
    assert page.status_code == 200
    text = page.text
    assert 'aria-label="Contig genes"' in text and 'class="gene member"' in text
    assert "ToyA" in text and "ToyB" in text                      # profile labels
    assert "/call/g1_Other_1" in text and "prophage" in text      # overlapping call and locus
    assert "Flanking genes (4)" in text                           # c_2, c_3, c_6, c_7 within ±2 kb
    assert 'data-notes-kind="system"' in text
    assert client.get("/call/missing").status_code == 404


def test_call_window_reaches_contig_ends(client):
    text = client.get("/call/g1_Toy_1", params={"flank": 40000}).text
    assert text.count('class="contig-end"') == 2


def test_calls_are_linked_from_system_pages_and_proteins(client):
    assert "/call/g1_Toy_1" in client.get("/system/defense/Toy").text
    assert "/call/g1_Toy_1" in client.get("/system/defense/Toy/loci").text
    assert "/call/g1_Toy_1" in client.get("/protein/c_4").text
    assert "/call/g1_Toy_1" in client.get("/genome/g1").text
    client.post("/api/notes", json={"kind": "system", "id": "g1_Toy_1", "flag": "verified"})
    assert "/call/g1_Toy_1" in client.get("/flags").text
