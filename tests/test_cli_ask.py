import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
DB_PATH = REPO_ROOT / "data" / "sharur.duckdb"


@pytest.mark.integration
def test_cli_ask_runs_with_db():
    if not DB_PATH.exists():
        pytest.skip("DuckDB not built; run ingest first.")
    cmd = [
        "python",
        "-m",
        "sharur.cli",
        "ask",
        "--db",
        str(DB_PATH),
        "find proteins with PFAM",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "Query:" in result.stdout
    assert "find_proteins" in result.stdout.lower()
