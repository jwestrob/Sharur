"""Browser views of CRISPRCasTyper subtype calls (crispr_cas_systems, crispr_array_types)."""

import duckdb
from fastapi.testclient import TestClient

from sharur.browser import create_app
from tests.test_browser_crispr import db  # noqa: F401  (fixture: one genome, one array, a Cas1 gene)


def _add_calls(path):
    conn = duckdb.connect(str(path))
    end = conn.execute("SELECT end_coord FROM loci WHERE locus_id = 'g1_CRISPR1'").fetchone()[0]
    conn.execute("""INSERT INTO crispr_cas_systems (system_id, genome_id, contig_id, start, end_coord, status,
                        prediction, prediction_cas, best_type, best_score, complete_interference,
                        complete_adaptation, strand_interference, strand_adaptation, genes_count, protein_ids,
                        profile_names, crispr_locus_ids, crispr_distances, caller)
                    VALUES ('cctyper:c1:200-1200', 'g1', 'c1', 200, 1200, 'crispr_cas', 'I-B', 'I-B', 'I-B', 9.0,
                            '60%', '100%', '1', '1', 1, 'c1_1', 'Cas1_0_IB', 'g1_CRISPR1', '300', 'cctyper-port 1.9.0'),
                           ('cctyper:c1:9000-9900', 'g1', 'c1', ?, ?, 'cas_putative', 'False', 'False', 'I-A', 1.0,
                            '0%', '0%', 'NA', 'NA', 1, 'c1_2', 'Cas6_0_IA', '', '', 'cctyper-port 1.9.0')""",
                 [end + 300, end + 1200])
    conn.execute("""INSERT INTO crispr_array_types VALUES ('g1_CRISPR1', 'g1', 'c1', 30, 5, 'I-B', 0.93, 'I-B',
                    100.0, 30.0, 0.1, TRUE, TRUE, 'cctyper-port 1.9.0')""")
    conn.execute("""INSERT INTO system_proteins (system_id, protein_id, system_source, position, profile_name, score)
                    VALUES ('cctyper:c1:200-1200', 'c1_1', 'cctyper', 0, 'Cas1_0_IB', 210.5)""")
    conn.close()


def test_cctyper_calls_across_pages(db, tmp_path):  # noqa: F811
    _add_calls(db)
    client = TestClient(create_app(db, background=False, notes_path=tmp_path / "n.sqlite"))
    systems = client.get("/systems").text
    assert "CRISPR-Cas subtypes" in systems and "/crispr-cas/calls?subtype=I-B" in systems
    assert "/crispr-cas/calls?status=cas_putative" in systems and "Cas-domain loci" in systems

    listing = client.get("/crispr-cas/calls", params={"subtype": "I-B"})
    assert listing.status_code == 200 and listing.text.count('class="stack-item"') == 1
    assert ">Cas1<" in listing.text  # member labeled by its profile's gene
    putative = client.get("/crispr-cas/calls", params={"status": "cas_putative"}).text
    assert putative.count('class="stack-item"') == 1 and "CRISPR-Cas candidates" in putative

    call = client.get("/cas-system/cctyper:c1:200-1200")
    assert call.status_code == 200
    assert "Cas1_0_IB" in call.text and "210.5" in call.text and 'aria-label="Genomic context"' in call.text
    assert "/crispr/g1_CRISPR1" in call.text and "0.93" in call.text
    candidate = client.get("/cas-system/cctyper:c1:9000-9900").text
    assert "Cas operon candidate" in candidate and "Unnamed (False)" in candidate
    assert client.get("/call/cctyper:c1:200-1200", follow_redirects=False).status_code == 307
    assert client.get("/cas-system/cctyper:none").status_code == 404

    array = client.get("/crispr/g1_CRISPR1").text
    assert "Repeat subtype <b>I-B</b>" in array and "trusted array" in array
    genome = client.get("/genome/g1").text
    assert "CRISPR-Cas I-B" in genome and "Cas candidate (putative)" in genome


def test_without_calls_pages_hide_the_panels(db, tmp_path):  # noqa: F811
    conn = duckdb.connect(str(db))
    conn.execute("DROP TABLE crispr_cas_systems")
    conn.execute("DROP TABLE crispr_array_types")
    conn.close()
    client = TestClient(create_app(db, background=False, notes_path=tmp_path / "n.sqlite"))
    assert "CRISPR-Cas subtypes" not in client.get("/systems").text
    assert client.get("/crispr/g1_CRISPR1").status_code == 200
    assert client.get("/genome/g1").status_code == 200
    assert client.get("/crispr-cas/calls").status_code == 200
