"""CAZyme consensus table: compute/load split and cache separation (scripts/classify_cazymes.py)."""

import importlib.util
from pathlib import Path

import pandas as pd
import pytest

from sharur.storage.duckdb_store import DuckDBStore

REPO_ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("classify_cazymes", REPO_ROOT / "scripts/classify_cazymes.py")
cc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cc)


def _db(tmp_path):
    path = tmp_path / "sharur.duckdb"
    store = DuckDBStore(path)
    store.conn.execute("""
        INSERT INTO bins (bin_id) VALUES ('b1');
        INSERT INTO contigs (contig_id, bin_id, length) VALUES ('c1', 'b1', 9000);
        INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand) VALUES
            ('c1_1', 'c1', 'b1', 1, 900, '+'), ('c1_2', 'c1', 'b1', 1000, 2000, '-');
        INSERT INTO annotations (annotation_id, protein_id, source, accession, name)
        VALUES (7, 'c1_1', 'cazy', 'GH999', 'GH999'), (8, 'c1_1', 'pfam', 'PF00150', 'Cellulase');
    """)
    store.close()
    return path


def test_load_replaces_cazy_rows_and_drops_unknown_proteins(tmp_path):
    db = _db(tmp_path)
    results = pd.DataFrame([
        {"protein_id": "c1_1", "family_class": "GH5", "cazy_class": "GH", "evalue": 1e-40,
         "score": 150.0, "n_tools": 3, "tools": "dbcan_hmm,dbcan_sub,diamond"},
        {"protein_id": "c1_2", "family_class": "GT2", "cazy_class": "GT", "evalue": 1e-30,
         "score": 90.0, "n_tools": 2, "tools": "dbcan_hmm,diamond"},
        {"protein_id": "elsewhere_9", "family_class": "CE4", "cazy_class": "CE", "evalue": 1e-20,
         "score": 60.0, "n_tools": 2, "tools": "dbcan_hmm,diamond"},
    ])
    stats = cc.load_cazyme_results(str(db), results, update_predicates=False, verbose=False)
    assert stats == {"annotations": 2, "proteins": 2, "dropped_unknown": 1}

    import duckdb
    con = duckdb.connect(str(db), read_only=True)
    rows = con.execute("SELECT protein_id, accession, annotation_id FROM annotations "
                       "WHERE source = 'cazy' ORDER BY 1").fetchall()
    assert [r[:2] for r in rows] == [("c1_1", "GH5"), ("c1_2", "GT2")]
    assert min(r[2] for r in rows) > 8  # new ids follow the existing maximum
    assert con.execute("SELECT COUNT(*) FROM annotations WHERE source = 'pfam'").fetchone()[0] == 1


def test_diamond_cache_is_separate_from_the_consensus_table(tmp_path):
    db = tmp_path / "sharur.duckdb"
    # a consensus table beside the database is never read back as DIAMOND evidence
    pd.DataFrame([{"protein_id": "p1", "family_class": "GH13", "evalue": 1e-50, "score": 200.0}]).to_csv(
        tmp_path / "cazyme_classification.tsv", sep="\t", index=False)
    assert cc.load_diamond_cache(str(db), verbose=False) is None
    cc._save_diamond_cache({"p1": {"GH13": (1e-50, 210.0)}}, tmp_path, verbose=False)
    cached = cc.load_diamond_cache(str(db), verbose=False)
    assert list(cached) == ["p1"] and list(cached["p1"]) == ["GH13"]
    assert cached["p1"]["GH13"][0] == pytest.approx(1e-50) and cached["p1"]["GH13"][1] == 210.0


def test_explicit_dbcan_directory(tmp_path, monkeypatch):
    (tmp_path / "CAZy.dmnd").write_text("")
    (tmp_path / "dbCAN.hmm").write_text("")
    assert cc.get_dbcan_path(tmp_path) == tmp_path
    monkeypatch.setenv("SHARUR_DBCAN_DIR", str(tmp_path))
    assert cc.get_dbcan_path() == tmp_path
