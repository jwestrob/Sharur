"""Browser enrichments: similar proteins, structures, agent findings."""

import json

import numpy as np
import pytest
from fastapi.testclient import TestClient

from sharur.browser import create_app
from sharur.browser.routes_insight import _compare, is_safe_select, safe_name
from sharur.predicates_v2.persistence import generate_and_persist_v2
from sharur.storage.duckdb_store import DuckDBStore


SEQ = "MKTAYIAKQRQISFVKSHFSRQ" * 10
PROTEINS = ["bin|1_c1_1", "bin|1_c1_2", "bin|1_c1_3", "bin|2_c1_1"]
PDB = (
    "ATOM      1  CA  MET A   1      11.104   6.134  -6.504  1.00  0.91           C\n"
    "ATOM      2  CA  LYS A   2      12.560   9.580  -6.123  1.00  0.88           C\n"
    "ATOM      3  CA  THR A   3      15.872  10.210  -4.300  1.00  0.62           C\nEND\n"
)


def _dataset(root, *, embeddings=True, structures=True, findings=True):
    root.mkdir(parents=True, exist_ok=True)
    db = root / "sharur.duckdb"
    store = DuckDBStore(str(db))
    store.conn.execute(f"""
        INSERT INTO bins (bin_id, taxonomy) VALUES ('bin|1', 'd__Archaea;p__Nanoarchaeota;c__Nanoarchaeia'),
                                                   ('bin|2', 'd__Archaea;p__Nanoarchaeota;c__Nanoarchaeia');
        INSERT INTO contigs (contig_id, bin_id, length) VALUES ('bin|1_c1', 'bin|1', 9000), ('bin|2_c1', 'bin|2', 9000);
        INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, gene_index, sequence,
                              sequence_length)
        VALUES ('bin|1_c1_1', 'bin|1_c1', 'bin|1', 1, 600, '+', 0, '{SEQ}', 220),
               ('bin|1_c1_2', 'bin|1_c1', 'bin|1', 700, 1400, '-', 1, '{SEQ}', 220),
               ('bin|1_c1_3', 'bin|1_c1', 'bin|1', 1500, 2100, '+', 2, '{SEQ}', 220),
               ('bin|2_c1_1', 'bin|2_c1', 'bin|2', 1, 600, '+', 0, '{SEQ}', 220);
        INSERT INTO annotations (annotation_id, protein_id, source, accession, name, description, evalue, start_aa, end_aa)
        VALUES (1, 'bin|1_c1_2', 'pfam', 'PF00005', 'ABC_tran', 'ABC transporter', 1e-40, 20, 160),
               (2, 'bin|2_c1_1', 'pfam', 'PF00005', 'ABC_tran', 'ABC transporter', 1e-35, 15, 150);
        INSERT INTO defense_systems (system_id, genome_id, system_type, system_subtype, genes_count, protein_ids)
        VALUES ('sys1', 'bin|1', 'Viperin', 'pVip', 1, 'bin|1_c1_3');
    """)
    generate_and_persist_v2(store, chunk_size=1, return_states=False, update_legacy_predicates=True)
    store.close()
    if embeddings:
        import h5py

        from sharur.storage.vector_store import build_vector_index

        (root / "embeddings").mkdir()
        h5 = root / "embeddings" / "protein_embeddings.h5"
        vectors = np.array([[1, 0, 0, 0], [0.9, 0.1, 0, 0], [0, 1, 0, 0], [0.8, 0.2, 0.1, 0]], dtype=np.float32)
        with h5py.File(h5, "w") as handle:
            handle.create_dataset("protein_ids", data=np.array(PROTEINS, dtype="S"))
            handle.create_dataset("embeddings", data=vectors)
            handle.attrs["model_name"] = "test"
        assert build_vector_index(h5).state == "available"
    if structures:
        sdir = root / "structures"
        sdir.mkdir()
        (sdir / "tip_01.pdb").write_text(PDB)
        (sdir / "tip_results.json").write_text(json.dumps({"tips": [
            {"label": "tip_01", "protein_id": "bin|1_c1_2", "pdb_path": "tip_01.pdb", "plddt": 0.81, "length": 120,
             "foldseek_hits": [{"target": "AF-Q9XYZ1-F1", "description": "ABC transporter ATP-binding", "evalue": 1e-12,
                                "prob": 1.0}]}]}))
        (sdir / f"{safe_name('bin|1_c1_3')}.pdb").write_text(PDB)  # canonical name, no record
        (sdir / "orphan_7.pdb").write_text(PDB)                     # no record, no matching protein: ignored
    if findings:
        rows = [
            {"id": "F1", "title": "Viperin flanks ABC transporters in bin|1", "category": "defense",
             "protein_ids": ["bin|1_c1_2", "bin|1_c1_3"], "n_genomes": 1,
             "evidence": {"genomes": ["bin|1"], "note": "two loci"},
             "verification": [
                 {"claim": "4 proteins", "query": "SELECT COUNT(*) FROM proteins", "expected": 4},
                 {"claim": "2 ABC_tran proteins", "query": "SELECT COUNT(DISTINCT protein_id) FROM annotations "
                  "WHERE name = 'ABC_tran'", "expected": 3},
                 {"claim": "python check", "query": "Sharur('x').card('y')", "expected": {"ok": True}},
                 {"claim": "malicious", "query": "COPY proteins TO 'leak.csv'", "expected": 0},
                 {"claim": "file read", "query": "SELECT * FROM read_csv('/etc/hosts')", "expected": 0}]},
            {"id": "F2", "title": "Contig-level observation", "category": "mobile", "contigs": ["bin|2_c1"]},
        ]
        (root / "exploration").mkdir()
        (root / "exploration" / "findings.jsonl").write_text("\n".join(json.dumps(r) for r in rows) + "\n")
    return db


@pytest.fixture
def client(tmp_path):
    return TestClient(create_app(_dataset(tmp_path / "ds"), background=False))


def test_similar_proteins_from_the_persistent_index(client):
    data = client.get("/api/similar/bin%7C1_c1_1").json()
    assert data["state"] == "available"
    ids = [n["protein_id"] for n in data["neighbors"]]
    assert ids[0] == "bin|1_c1_2" and "bin|1_c1_1" not in ids
    top = data["neighbors"][0]
    assert top["architecture"] == "ABC_tran" and top["lineage"] == "Nanoarchaeia" and top["length"] == 220
    assert data["neighbors"][0]["similarity"] >= data["neighbors"][-1]["similarity"]
    page = client.get("/protein/bin%7C1_c1_1")
    assert 'data-similar-url="/api/similar/bin%7C1_c1_1"' in page.text and "Loading" in page.text


def test_similar_panel_explains_a_missing_index(tmp_path):
    client = TestClient(create_app(_dataset(tmp_path / "ds", embeddings=False), background=False))
    assert client.get("/api/similar/bin%7C1_c1_1").json()["state"] == "no_embeddings"
    assert "no protein embeddings" in client.get("/protein/bin%7C1_c1_1").text


def test_structures_come_from_explicit_records_only(client):
    insight = client.app.state.insight
    assert set(insight.structures) == {"bin|1_c1_2", "bin|1_c1_3"}       # orphan_7.pdb is not attached
    record = insight.structures["bin|1_c1_2"][0]
    assert record.plddt == pytest.approx(81.0)                            # fractional pLDDT normalized
    assert record.hits[0]["target"] == "AF-Q9XYZ1-F1"
    page = client.get("/protein/bin%7C1_c1_2")
    assert 'id="mol"' in page.text and "AF-Q9XYZ1-F1" in page.text and "partial model: 120 of 220 aa" in page.text
    assert client.get(f"/structure-file/{record.key}").text.startswith("ATOM")
    assert client.get("/structure-file/999").status_code == 404
    assert 'id="mol"' not in client.get("/protein/bin%7C1_c1_1").text


def test_findings_link_to_proteins_genomes_clades_and_systems(client):
    insight = client.app.state.insight
    assert [f["id"] for f in insight.findings_for_protein("bin|1_c1_2")] == ["F1"]
    assert {f["id"] for f in insight.findings_for_genome("bin|1")} == {"F1"}
    assert {f["id"] for f in insight.findings_for_genome("bin|2")} == {"F2"}      # via its contig
    assert [f["id"] for f in insight.findings_for_system("Viperin")] == ["F1"]
    assert {f["id"] for f in insight.findings_for_clade("class", "Nanoarchaeia")} == {"F1", "F2"}
    assert "Agent findings" in client.get("/genome/bin%7C1").text
    assert "Agent findings" in client.get("/system/defense/Viperin").text
    assert "Agent findings" in client.get("/taxa/class/Nanoarchaeia").text
    index = client.get("/findings")
    assert index.status_code == 200 and "Viperin flanks" in index.text
    assert "Contig-level" not in client.get("/findings", params={"category": "defense"}).text
    detail = client.get("/finding/F1")
    assert detail.status_code == 200 and "Re-run checks" in detail.text and "bin|1_c1_3" in detail.text
    assert client.get("/finding/NOPE").status_code == 404


def test_verification_runs_only_safe_selects(client):
    results = client.post("/api/finding/F1/verify").json()["results"]
    statuses = [r["status"] for r in results]
    assert statuses == ["pass", "fail", "skipped", "skipped", "skipped"]
    assert results[1]["observed"] == 2
    assert "SELECT" in results[3]["detail"] and "not allowed" in results[4]["detail"]
    assert client.post("/api/finding/NOPE/verify").status_code == 404


def test_sql_guard_and_comparison():
    assert is_safe_select("WITH a AS (SELECT 1) SELECT * FROM a")[0]
    for bad in ["SELECT 1; SELECT 2", "PRAGMA version", "ATTACH 'x.db'", "SELECT getenv('HOME')",
                "SELECT * FROM glob('*')", "INSERT INTO bins VALUES ('x')", ""]:
        assert not is_safe_select(bad)[0], bad
    assert _compare(3, [(3,)]) == ("pass", 3)
    assert _compare(0.5, [(0.5000000001,)])[0] == "pass"
    assert _compare(["a", "b"], [("a",), ("b",)])[0] == "pass"
    assert _compare(None, [(1,)])[0] == "unchecked"


def test_pages_degrade_without_structures_or_findings(tmp_path):
    client = TestClient(create_app(_dataset(tmp_path / "ds", embeddings=False, structures=False, findings=False),
                                   background=False))
    assert 'id="mol"' not in client.get("/protein/bin%7C1_c1_2").text
    assert "No <code>findings.jsonl</code>" in client.get("/findings").text
    assert "Agent findings" not in client.get("/genome/bin%7C1").text
