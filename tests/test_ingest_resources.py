from pathlib import Path


def test_dbcan_database_present():
    db_root = Path("data/dbcan_db")
    required = ["dbCAN.hmm", "CAZy.dmnd", "dbCAN-sub.hmm", "fam-substrate-mapping.tsv"]
    for fname in required:
        assert (db_root / fname).exists(), f"Missing dbCAN asset: {fname}"
