"""CRISPR array pages: repeats and spacers recovered from the assembly, Cas-domain context."""

import json

import pytest
from fastapi.testclient import TestClient

from sharur.browser import create_app
from sharur.browser.routes_crispr import find_repeats, spacers_between
from sharur.predicates_v2.persistence import generate_and_persist_v2
from sharur.storage.duckdb_store import DuckDBStore

REPEAT = "GTTTCAATCCCATAAAGGTAGTTTAGGAAC"  # 30 bp
SPACERS = ["ACGTTGCAAGTCTGATCGGATCCATGCAATGC", "TTGACCGGTACGATTAGCCATGGCTTAACGGA",
           "CATGCCTTAGGACCTAGTTGGCAATCGTAACG", "GGATCCAATGCTAGCTTGACCATCGGTACCAT"]


def _array() -> str:
    degenerate = REPEAT[:20] + "AAAACCCCAA"  # terminal repeat, 10 mismatches (33%)? keep 7 for the 30% budget
    degenerate = REPEAT[:23] + "CCCCCCC"
    parts = []
    for s in SPACERS:
        parts += [REPEAT, s]
    return "".join(parts) + degenerate


@pytest.fixture
def db(tmp_path):
    data = tmp_path / "ds"
    (data / "genomes_fna").mkdir(parents=True)
    left, right = "A" * 1500, "T" * 1500
    contig = left + _array() + right
    (data / "genomes_fna" / "g1.fna").write_text(">c1 test\n" + "\n".join(contig[i:i + 70] for i in range(0, len(contig), 70)) + "\n")
    start, end = len(left) + 1, len(left) + len(_array())
    path = data / "sharur.duckdb"
    store = DuckDBStore(str(path))
    meta = json.dumps({"id": "CRISPR1", "metadata": {"rpt_unit_seq": REPEAT}})
    store.conn.execute(f"""
        INSERT INTO bins (bin_id, taxonomy) VALUES ('g1', 'd__Archaea;p__P;c__C');
        INSERT INTO contigs (contig_id, bin_id, length) VALUES ('c1', 'g1', {len(contig)});
        INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, gene_index, sequence_length)
        VALUES ('c1_1', 'c1', 'g1', 200, 1200, '+', 0, 333), ('c1_2', 'c1', 'g1', {end + 300}, {end + 1200}, '-', 1, 300);
        INSERT INTO annotations (annotation_id, protein_id, source, accession, name, description, evalue, start_aa, end_aa)
        VALUES (1, 'c1_1', 'pfam', 'PF01867', 'Cas_Cas1', 'CRISPR-associated protein Cas1', 1e-50, 5, 300),
               (2, 'c1_2', 'pfam', 'PF00005', 'ABC_tran', 'ABC transporter', 1e-40, 10, 150);
        INSERT INTO loci (locus_id, locus_type, contig_id, start, end_coord, confidence, metadata)
        VALUES ('g1_CRISPR1', 'crispr', 'c1', {start}, {end}, 1.0, '{meta}');
    """)
    generate_and_persist_v2(store, chunk_size=10, return_states=False, update_legacy_predicates=True)
    store.close()
    return path


def test_repeat_recovery_includes_a_diverged_terminal_repeat():
    seq = "A" * 100 + _array() + "T" * 100
    repeats = find_repeats(seq, REPEAT, 101, 100 + len(_array()))
    assert len(repeats) == 5
    assert [r["diverged"] for r in repeats] == [False, False, False, False, True]
    assert [s["seq"] for s in spacers_between(seq, repeats)] == SPACERS


def test_crispr_pages(db, tmp_path):
    client = TestClient(create_app(db, background=False, notes_path=tmp_path / "notes.sqlite"))
    catalog = client.get("/crispr")
    assert catalog.status_code == 200 and "g1_CRISPR1" in catalog.text
    page = client.get("/crispr/g1_CRISPR1")
    assert page.status_code == 200
    assert "<b>5</b><span>repeats" in page.text and "<b>4</b><span>spacers" in page.text
    assert 'aria-label="CRISPR array"' in page.text and 'aria-label="Genomic context"' in page.text
    assert "Cas_Cas1" in page.text and "<b>1</b><span>Cas-domain genes" in page.text
    assert 'class="repeat variant"' in page.text
    fasta = client.get("/crispr/g1_CRISPR1/spacers.fasta").text.splitlines()
    assert fasta[0].startswith(">g1_CRISPR1_spacer1") and fasta[1::2] == SPACERS
    assert "CRISPR arrays" in client.get("/genome/g1").text
    assert client.post("/api/notes", json={"kind": "crispr", "id": "g1_CRISPR1", "flag": "interesting"}).status_code == 200
    assert client.get("/crispr/missing").status_code == 404


def test_minced_report_is_the_primary_source(db, tmp_path):
    from sharur.crispr import parse_minced_text

    data = db.parent
    contig = (data / "genomes_fna" / "g1.fna").read_text().split("\n", 1)[1].replace("\n", "")
    start = 1501
    rows, pos = [], start
    for s in SPACERS:
        rows.append(f"{pos}\t\t{REPEAT}\t{s}\t[ {len(REPEAT)}, {len(s)} ]")
        pos += len(REPEAT) + len(s)
    rows.append(f"{pos}\t\t{REPEAT}")
    report = data / "stage05c_crispr" / "g1_crispr.txt"
    report.parent.mkdir()
    report.write_text(f"Sequence 'c1' ({len(contig)} bp)\n\nCRISPR 1   Range: {start} - {pos + len(REPEAT) - 1}\n"
                      "POSITION\tREPEAT\t\t\t\tSPACER\n--------\t----\t----\n" + "\n".join(rows) +
                      "\n--------\t----\t----\nRepeats: 5\tAverage Length: 30\t\tAverage Length: 32\n")
    (array,) = parse_minced_text(report)
    assert (array["contig"], len(array["repeats"]), [s["seq"] for s in array["spacers"]]) == ("c1", 5, SPACERS)
    page = TestClient(create_app(db, background=False, notes_path=tmp_path / "n.sqlite")).get("/crispr/g1_CRISPR1")
    assert "from MinCED" in page.text and "<b>5</b><span>repeats" in page.text
