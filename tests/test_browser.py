"""Read-only dataset browser."""

import pytest
from fastapi.testclient import TestClient

from sharur.browser import charts, create_app
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


def test_protein_neighborhood_draws_crispr_arrays_in_the_plotted_stretch(db):
    store = DuckDBStore(str(db))
    store.conn.execute("""
        INSERT INTO contigs (contig_id, bin_id, length, length_source) VALUES ('bin|1_c2', 'bin|1', 5000, 'assembly');
        INSERT INTO loci (locus_id, locus_type, contig_id, start, end_coord, confidence) VALUES
            ('bin|1_CRISPR1', 'crispr', 'bin|1_c1', 2300, 2700, 1.0),
            ('bin|1_CRISPR2', 'crispr', 'bin|1_c2', 100, 400, 1.0);
    """)
    store.close()
    page = TestClient(create_app(db, background=False)).get("/protein/bin%7C1_c1_2")
    assert page.status_code == 200
    # the array past the last gene sits between it and the contig end, so the plot widens to show it
    assert '<a href="/crispr/bin%7C1_CRISPR1">' in page.text and "bin|1_CRISPR1 · 2,300–2,700" in page.text
    assert 'viewBox="0 0 1000 88"' in page.text and "CRISPR array</span>" in page.text
    assert "bin|1_CRISPR2" not in page.text  # another contig's array


def test_neighborhood_features_clip_overlaps_and_widen_for_flanking_ones():
    genes = [{"protein_id": "g1", "start": 1000, "end": 2000, "strand": "+", "number": 1},
             {"protein_id": "g2", "start": 3000, "end": 4000, "strand": "-", "number": 2}]
    plain = charts.neighborhood(genes, start_edge=False, end_edge=False)
    assert 'viewBox="0 0 1000 72"' in plain and "<rect" not in plain and "976.0,28" in plain
    inner = charts.neighborhood(genes, start_edge=False, end_edge=False,
                                features=[{"start": 500, "end": 2500, "kind": "prophage", "label": "Prophage region"}])
    # overlapping feature is clipped at the first gene, so the genes keep their positions
    assert 'points="24.0,28' in inner and '<rect x="24.0" y="60"' in inner
    flank = charts.neighborhood(genes, start_edge=False, end_edge=True,
                                features=[{"start": 4100, "end": 4300, "kind": "crispr", "label": "CRISPR array"}])
    # flanking array widens the range: the last gene now ends before the right padding
    assert '<rect x="' in flank and 'points="24.0,28' in flank and "976.0,28" not in flank
    # a narrow band at the right edge carries its label on the left, inside the plot
    assert 'text-anchor="end">CRISPR array</text>' in flank


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


def test_vog_families_use_the_vogdb_annotation_table(db, tmp_path, monkeypatch):
    import duckdb

    from sharur.predicates.mappings import vog_map

    table = tmp_path / "vog.annotations.tsv"
    table.write_text("#GroupName\tProteinCount\tSpeciesCount\tFunctionalCategory\tConsensusFunctionalDescription\n"
                     "VOG00042\t40\t31\tXr\tsp|P0A000|INT_LAMBD Integrase\n")
    monkeypatch.setenv("SHARUR_VOG_ANNOTATIONS", str(table))
    vog_map.load_vog_annotations.cache_clear()
    conn = duckdb.connect(str(db))
    conn.execute("""INSERT INTO annotations (annotation_id, protein_id, source, accession, name, description, evalue,
                                             start_aa, end_aa)
                    VALUES (99, 'bin|1_c1_3', 'vogdb', 'VOG00042', 'VOG00042', '', 1e-20, 5, 180)""")
    conn.close()
    try:
        client = TestClient(create_app(db, background=False))
        catalog = client.get("/vogs")
        assert catalog.status_code == 200 and "Integrase" in catalog.text and "sp|P0A000" not in catalog.text
        page = client.get("/vog/VOG00042")
        assert page.status_code == 200 and "Replication" in page.text and "40 viral proteins from 31" in page.text
        assert "Integrase (VOG00042)" in client.get("/contig/bin%7C1_c1").text  # hit labels use the description
        assert client.get("/search", params={"q": "VOG00042"}, follow_redirects=False).headers["location"] == "/vog/VOG00042"
        assert client.get("/vog/VOG99999").status_code == 404
    finally:
        vog_map.load_vog_annotations.cache_clear()


def test_scoped_search(client):
    page = client.get("/search", params={"q": "abc in bin|1"})
    assert page.status_code == 200 and "bin|1_c1_2" in page.text and "1 proteins in 1 of 1 genome" in page.text
    clade = client.get("/search", params={"q": "transporter in Archaea"})
    assert "bin|1_c1_2" in clade.text and "Search within domain" in clade.text
    assert "No protein here matches" in client.get("/search", params={"q": "rubisco in bin|1"}).text
    assert "By genome" in page.text and "/matrix?a=bin%7C1&kind=pfam&features=PF00005" in page.text
    # terms match where a word starts ("tran" in ABC_tran), never mid-word ("BC_tran")
    assert "bin|1_c1_2" in client.get("/search", params={"q": "tran in bin|1"}).text
    assert "No protein here matches" in client.get("/search", params={"q": "BC_tran in bin|1"}).text


def test_search_page_groups_matches_and_offers_near_misses(client):
    page = client.get("/search", params={"q": "ABC_tr"}).text
    assert 'class="sr-chips"' in page and "<mark>ABC_tr</mark>" in page and "best match" in page
    assert "Did you mean" in client.get("/search", params={"q": "ABC_trn"}).text


def test_matrix_pages(client):
    assert client.get("/matrix").status_code == 200
    page = client.get("/matrix", params={"a": "bin|1", "kind": "pfam", "features": "PF00005"})
    assert page.status_code == 200 and 'id="matrix-data"' in page.text and "ABC_tran" in page.text
    tsv = client.get("/matrix.tsv", params={"a": "domain:Archaea", "kind": "pfam", "features": "PF00005"})
    assert tsv.status_code == 200 and tsv.text.splitlines()[1].startswith("bin|1\t") and tsv.text.rstrip().endswith("1")
    assert client.get("/matrix", params={"a": "no such clade"}).status_code == 404
    # one genome: no variable features, a clear message
    assert "No features to show" in client.get("/matrix", params={"a": "bin|1", "kind": "pfam"}).text
    # category (curated function label) and keyword column searches
    by_category = client.get("/matrix", params={"a": "domain:Archaea", "kind": "pfam", "category": "ABC transporter"})
    assert 'id="matrix-data"' in by_category.text and "ABC_tran" in by_category.text and "category ABC transporter" in by_category.text
    by_keyword = client.get("/api/v1/matrix", params={"a": "domain:Archaea", "kind": "pfam", "q": "abc_tran"}).json()
    assert [f["id"] for f in by_keyword["features"]] == ["PF00005"]
    assert "Unknown category" in client.get("/matrix", params={"a": "bin|1", "category": "nonsense label"}).text
    assert 'id="matrix-data"' not in client.get("/matrix", params={"a": "bin|1", "kind": "system",
                                                                   "category": "ribosomal_protein"}).text


@pytest.fixture
def ribosomal_client(tmp_path):
    """One genome with uL2 (its N-terminal profile stored by name), an L11 methyltransferase and ABC_tran."""
    path = tmp_path / "ribo.duckdb"
    store = DuckDBStore(str(path))
    store.conn.execute(f"""
        INSERT INTO bins (bin_id, completeness, contamination, taxonomy) VALUES ('bin|1', 91.5, 1.2, 'd__Archaea');
        INSERT INTO contigs (contig_id, bin_id, length, length_source) VALUES ('bin|1_c1', 'bin|1', 9000, 'assembly');
        INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, gene_index,
                              sequence, sequence_length, partial)
        VALUES ('bin|1_c1_1', 'bin|1_c1', 'bin|1', 1, 600, '+', 0, '{SEQ}', 220, '00'),
               ('bin|1_c1_2', 'bin|1_c1', 'bin|1', 700, 1400, '-', 1, '{SEQ}', 220, '00');
        INSERT INTO annotations (annotation_id, protein_id, source, accession, name, description, evalue, score,
                                 start_aa, end_aa)
        VALUES (1, 'bin|1_c1_1', 'pfam', 'Ribosomal_L2', 'Ribosomal_L2', '', 1e-30, 110.0, 5, 80),
               (2, 'bin|1_c1_1', 'pfam', 'PF03947.19', 'Ribosomal_L2_C', 'Ribosomal Proteins L2, C-terminal domain',
                1e-50, 180.0, 90, 210),
               (3, 'bin|1_c1_2', 'pfam', 'PF06325.18', 'PrmA', 'Ribosomal protein L11 methyltransferase (PrmA)',
                1e-60, 200.0, 10, 200),
               (4, 'bin|1_c1_2', 'pfam', 'PF00005.31', 'ABC_tran', 'ABC transporter', 1e-20, 80.0, 1, 9),
               (5, 'bin|1_c1_2', 'kofam', 'K02986', 'K02986', 'small subunit ribosomal protein S4', 1e-30, 90.0,
                NULL, NULL);
    """)
    store.close()
    return TestClient(create_app(path, background=False))


def test_matrix_ribosomal_pfam_panel(ribosomal_client):
    client = ribosomal_client
    page = client.get("/matrix", params={"a": "all", "kind": "ko", "category": "ribosomal_proteins_pfam"})
    assert page.status_code == 200 and "panel Ribosomal proteins (Pfam)" in page.text
    assert '<option value="pfam" selected' in page.text                 # the panel sets the feature kind
    api = client.get("/api/v1/matrix", params={"a": "all", "category": "Ribosomal proteins (Pfam)"}).json()
    assert api["kind"] == "mixed"                   # Pfam families, and the KO for archaeal uS4
    cols = {f["id"]: f for f in api["features"]}
    assert set(cols) == {"PF03947", "Ribosomal_L2", "K02986"}
    assert cols["PF03947"]["label"] == "uL2 (L2p) \u00b7 L2_C" and cols["Ribosomal_L2"]["label"] == "uL2 (L2p) \u00b7 L2"
    assert cols["K02986"]["kind"] == "ko" and cols["K02986"]["label"] == "uS4 (S4p)"
    assert cols["PF03947"]["name"].startswith("uL2 (L2p) \u00b7 ") and "KOfam call" in cols["K02986"]["name"]
    assert "Ribosomal proteins (Pfam) <span class=\"count\">3</span>" in page.text      # the preset chip
    # the broader category also reaches the L11 methyltransferase, and points to the panel
    broad = client.get("/api/v1/matrix", params={"a": "all", "kind": "pfam", "category": "ribosomal_protein"}).json()
    assert {"PF03947", "Ribosomal_L2", "PF06325"} <= {f["id"] for f in broad["features"]}
    assert any("curated panel Ribosomal proteins (Pfam)" in note for note in broad["notes"])
    tsv = client.get("/matrix.tsv", params={"a": "all", "category": "ribosomal_proteins_pfam"})
    assert tsv.status_code == 200 and "K02986" in tsv.text.splitlines()[0]


def test_json_api(client):
    index = client.get("/api/v1").json()
    assert any(e["path"] == "/api/v1/genomes" for e in index["endpoints"])
    assert client.get("/api").status_code == 200
    genomes = client.get("/api/v1/genomes").json()
    assert genomes["total"] == 1 and genomes["rows"][0]["bin_id"] == "bin|1"
    assert client.get("/api/v1/genome/bin|1").json()["completeness"] == pytest.approx(91.5)
    protein = client.get("/api/v1/protein/bin|1_c1_2").json()
    assert protein["found"] and "sequence" not in protein
    assert client.get("/api/v1/protein/bin|1_c1_2", params={"sequence": 1}).json()["sequence"]
    assert client.get("/api/v1/protein/nope").status_code == 404
    matrix = client.get("/api/v1/matrix", params={"a": "bin|1", "kind": "pfam", "features": "PF00005"}).json()
    assert matrix["values"] == [[1.0]]
    assert client.get("/api/v1/search", params={"q": "abc in bin|1"}).json()["total"] == 1
    assert client.get("/api/v1/discover/nope").status_code == 404
