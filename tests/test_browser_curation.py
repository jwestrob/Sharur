"""Browser curation: stacked loci, notes and flags, triage, collection, previews."""

import pytest
from fastapi.testclient import TestClient

from sharur.browser import create_app
from sharur.browser.notes import NotesStore
from sharur.predicates_v2.persistence import generate_and_persist_v2
from sharur.storage.duckdb_store import DuckDBStore


SEQ = "MKTAYIAKQRQISFVKSHFSRQ" * 10


@pytest.fixture
def db(tmp_path):
    path = tmp_path / "data" / "sharur.duckdb"
    path.parent.mkdir()
    store = DuckDBStore(str(path))
    genes = []
    for b, strand in (("g1", "+"), ("g2", "-")):
        for i in range(6):
            genes.append(f"('{b}_c_{i}', '{b}_c', '{b}', {1 + 1000 * i}, {900 + 1000 * i}, '{strand}', {i}, '{SEQ}', 220)")
    store.conn.execute(f"""
        INSERT INTO bins (bin_id, taxonomy) VALUES ('g1', 'd__Bacteria;p__P;c__C1'), ('g2', 'd__Bacteria;p__P;c__C2');
        INSERT INTO contigs (contig_id, bin_id, length) VALUES ('g1_c', 'g1', 6000), ('g2_c', 'g2', 6000);
        INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, gene_index, sequence,
                              sequence_length) VALUES {", ".join(genes)};
        INSERT INTO annotations (annotation_id, protein_id, source, accession, name, description, evalue, start_aa, end_aa)
        VALUES (1, 'g1_c_2', 'pfam', 'PF00005', 'ABC_tran', 'ABC transporter', 1e-40, 10, 150),
               (2, 'g2_c_3', 'pfam', 'PF00005', 'ABC_tran', 'ABC transporter', 1e-40, 10, 150),
               (3, 'g1_c_3', 'pfam', 'PF00004', 'AAA', 'ATPase', 1e-30, 10, 150),
               (4, 'g2_c_2', 'pfam', 'PF00004', 'AAA', 'ATPase', 1e-30, 10, 150);
        INSERT INTO defense_systems (system_id, genome_id, contig_id, system_type, system_subtype, genes_count,
                                     protein_ids, profile_names)
        VALUES ('g1_Toy_1', 'g1', 'g1_c', 'Toy', 'Toy_I', 2, 'g1_c_2,g1_c_3', 'Toy__ToyA,Toy__ToyB'),
               ('g2_Toy_1', 'g2', 'g2_c', 'Toy', 'Toy_I', 2, 'g2_c_3,g2_c_2', 'Toy__ToyA,Toy__ToyB');
    """)
    generate_and_persist_v2(store, chunk_size=2, return_states=False, update_legacy_predicates=True)
    store.close()
    return path


@pytest.fixture
def client(db, tmp_path):
    c = TestClient(create_app(db, background=False, notes_path=tmp_path / "notes.sqlite"))
    c.cookies.set("sharur_user", "ada")
    return c


def test_system_loci_align_on_the_core_component(client):
    page = client.get("/system/defense/Toy/loci")
    assert page.status_code == 200
    assert page.text.count('class="stack-item"') == 2
    assert "<b>ToyA</b>" in page.text            # core component named in the lede
    assert "reversed" in page.text               # the minus-strand call is turned around
    assert "ToyB" in page.text and "Triage these calls" in page.text
    assert client.get("/system/defense/Toy/loci", params={"subtype": "Toy_II"}).text.count('class="stack-item"') == 0
    assert client.get("/system/defense/Nope/loci").status_code == 404


def test_family_stacks(client):
    page = client.get("/stack/domain/PF00005", params={"flank": 2})
    assert page.status_code == 200 and page.text.count('class="stack-item"') == 2
    assert "AAA" in page.text  # the neighbor family appears in the legend
    assert client.get("/stack/domain/PF00005", params={"rank": "class", "clade": "C1"}).text.count(
        'class="stack-item"') == 1
    assert client.get("/stack/bogus/x").status_code == 404


def test_flags_and_notes_round_trip(client, tmp_path):
    state = client.post("/api/notes", json={"kind": "protein", "id": "g1_c_2", "flag": "verified"}).json()
    assert state["flags"] == {"verified": ["ada"]}
    assert client.post("/api/notes", json={"kind": "protein", "id": "g1_c_2", "flag": "verified"}).json()["flags"] == {}
    state = client.post("/api/notes", json={"kind": "protein", "id": "g1_c_2", "text": "operon with AAA"}).json()
    note_id = state["notes"][0]["id"]
    assert state["notes"][0]["author"] == "ada"
    client.post("/api/notes", json={"kind": "system", "id": "g2_Toy_1", "flag": "suspicious"})
    assert client.post("/api/notes", json={"kind": "nope", "id": "x", "flag": "verified"}).status_code == 400
    other = TestClient(client.app)
    other.cookies.set("sharur_user", "bob")
    assert other.post(f"/api/notes/{note_id}/delete").status_code == 403
    flags = client.get("/flags")
    assert "operon with AAA" in flags.text and "Toy in g2" in flags.text
    tsv = client.get("/flags.tsv").text.splitlines()
    assert tsv[0].startswith("id\tkind\tentity") and len(tsv) == 3  # header, note, system flag (toggled-off flag hidden)
    assert client.post(f"/api/notes/{note_id}/delete").status_code == 200
    assert NotesStore(tmp_path / "notes.sqlite").state("protein", "g1_c_2")["notes"] == []


def test_triage_walks_items(client):
    first = client.get("/triage", params={"source": "system", "kind": "defense", "id": "Toy"})
    assert first.status_code == 200 and "1 / 2" in first.text and 'class="stack-row"' in first.text
    second = client.get("/triage", params={"source": "system", "kind": "defense", "id": "Toy", "i": 5})
    assert "2 / 2" in second.text  # clamped to the last item
    carriers = client.get("/triage", params={"source": "family", "kind": "domain", "id": "PF00005"})
    assert "1 / 2" in carriers.text and "ABC_tran" in carriers.text
    client.post("/api/notes", json={"kind": "protein", "id": "g2_c_3", "flag": "follow_up"})
    flagged = client.get("/triage", params={"source": "flags", "flag": "follow_up"})
    assert "1 / 1" in flagged.text and "g2_c_3" in flagged.text
    assert "Nothing to review" in client.get("/triage", params={"source": "flags", "flag": "verified"}).text


def test_collection_previews_and_gene_navigation(client):
    rows = client.post("/api/collection", json={"proteins": ["g1_c_2"], "genomes": ["g2"]}).json()
    assert [r["kind"] for r in rows] == ["protein", "genome"]
    fasta = client.post("/api/fasta", data={"ids": "g1_c_2\ng2_c_3"}).text.splitlines()
    assert fasta[0] == ">g1_c_2 genome=g1" and sum(line.startswith(">") for line in fasta) == 2
    preview = client.get("/api/preview", params={"href": "/protein/g1_c_2"})
    assert preview.status_code == 200 and "ABC_tran" in preview.text and SEQ[:10] not in preview.text
    assert "C2" in client.get("/api/preview", params={"href": "/genome/g2"}).text
    page = client.get("/protein/g1_c_2").text
    assert 'data-key="prev-gene"' in page and 'data-key="next-gene"' in page and 'data-collect="protein"' in page
    assert 'data-notes-kind="protein"' in page
