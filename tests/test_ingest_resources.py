from pathlib import Path

import pytest


DB_ROOT = Path("data/dbcan_db")


@pytest.mark.skipif(
    not DB_ROOT.is_dir(),
    reason="optional local dbCAN reference bundle is not installed",
)
def test_dbcan_database_present():
    required = ["dbCAN.hmm", "CAZy.dmnd", "fam-substrate-mapping.tsv"]
    for fname in required:
        assert (DB_ROOT / fname).exists(), f"Missing dbCAN asset: {fname}"
    # dbCAN releases name the sub-family profiles dbCAN_sub.hmm; older bundles used dbCAN-sub.hmm.
    assert any((DB_ROOT / n).exists() for n in ("dbCAN-sub.hmm", "dbCAN_sub.hmm")), "Missing dbCAN sub-family HMMs"
