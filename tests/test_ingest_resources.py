from pathlib import Path

import pytest


DB_ROOT = Path("data/dbcan_db")


@pytest.mark.skipif(
    not DB_ROOT.is_dir(),
    reason="optional local dbCAN reference bundle is not installed",
)
def test_dbcan_database_present():
    required = ["dbCAN.hmm", "CAZy.dmnd", "dbCAN-sub.hmm", "fam-substrate-mapping.tsv"]
    for fname in required:
        assert (DB_ROOT / fname).exists(), f"Missing dbCAN asset: {fname}"
