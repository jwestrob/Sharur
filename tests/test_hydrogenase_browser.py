"""Hydrogenase neighborhood evidence and the /hydrogenases pages."""

import pytest
from fastapi.testclient import TestClient

from sharur.browser import create_app
from sharur.hydrogenase.classifier import CLASSIFICATION_COLUMNS
from sharur.hydrogenase.neighborhood import neighborhood_contexts
from sharur.storage.duckdb_store import DuckDBStore

COLUMNS = [c.split()[0] for c in CLASSIFICATION_COLUMNS]

# focal proteins at gene 1 of their contig; neighbor at gene 2 carries the marker
CASES = {
    "maturation": ("K04651", "K04651"),        # hypA: KEGG names it a hydrogenase -> supported
    "nuo": ("K00333", "K00333"),               # nuoD -> Complex I context
    "small_subunit": ("PF01058", "Oxidored_q6"),  # shared with [NiFe] small subunits -> no marker
}


@pytest.fixture
def db(tmp_path):
    path = tmp_path / "sharur.duckdb"
    s = DuckDBStore(str(path))
    c = s.conn
    c.execute("INSERT INTO bins (bin_id, taxonomy) VALUES ('g1', 'd__Archaea;p__P;c__C')")
    ann, rows, n = [], [], 0
    for case, (acc, name) in CASES.items():
        c.execute("INSERT INTO contigs (contig_id, bin_id, length) VALUES (?, 'g1', 9000)", [case])
        for i in range(3):
            c.execute("INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, gene_index, "
                      "sequence_length) VALUES (?, ?, 'g1', ?, ?, '+', ?, 300)", [f"{case}_{i}", case, 1000 * i + 1, 1000 * i + 900, i])
        n += 1
        ann.append((n, f"{case}_2", "kofam" if acc.startswith("K") else "pfam", acc, name))
        rows.append(dict.fromkeys(COLUMNS) | {
            "protein_id": f"{case}_1", "outcome": "assigned", "reference_label": "[NiFe]_Group_4e",
            "reference_class": "NiFe", "reference_subgroup": "Group_4e", "interpretation_status": "characterized",
            "reference_role": "Ferredoxin-coupled, Ech-type", "pident": 60.0, "has_nifese_hases": case == "maturation",
            "has_fe_hyd": False, "has_complex1": True, "has_hmd": False,
            "curation_status": "domain_check_cleared" if case == "maturation" else "needs_curation",
            "ko_support": "none", "reference_release": "MM2022", "classifier_version": "test"})
    c.execute("INSERT INTO contigs (contig_id, bin_id, length) VALUES ('loose', 'g1', 900)")
    c.execute("INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand) "
              "VALUES ('loose', 'loose', 'g1', 0, 900, '+')")
    c.executemany("INSERT INTO annotations (annotation_id, protein_id, source, accession, name) VALUES (?, ?, ?, ?, ?)", ann)
    c.execute(f"CREATE TABLE hydrogenase_classifications ({', '.join(CLASSIFICATION_COLUMNS)})")
    c.executemany(f"INSERT INTO hydrogenase_classifications VALUES ({','.join('?' * len(COLUMNS))})",
                  [[r[k] for k in COLUMNS] for r in rows])
    s.close()
    return path


def test_neighborhood_verdicts_read_kegg_hydrogenase_and_nuo_kos(db):
    store = DuckDBStore(str(db), read_only=True)
    try:
        contexts = neighborhood_contexts(store, ["maturation_1", "nuo_1", "small_subunit_1", "loose"])
    finally:
        store.close()
    assert {k: v.verdict for k, v in contexts.items()} == {
        "maturation_1": "supported", "nuo_1": "complex_i_context", "small_subunit_1": "no_markers",
        "loose": "no_position"}
    assert contexts["maturation_1"].hydrogenase_genes == {"maturation_2"}


def test_hydrogenase_pages(db):
    client = TestClient(create_app(db, background=False))

    page = client.get("/hydrogenases")
    assert page.status_code == 200 and "need curation" in page.text and "[NiFe] 4e" in page.text
    calls = client.get("/hydrogenases/calls", params={"class": "NiFe", "evidence": "flagged", "verdict": "complex_i_context"})
    assert calls.status_code == 200 and "nuo_1" in calls.text and "maturation_1" not in calls.text
    assert "Hydrogenases" in client.get("/systems").text
    assert "by nearest HydDB reference" in client.get("/protein/maturation_1").text
