"""Test that CLI search command works (replaces stale ask command test)."""

import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
DB_PATH = REPO_ROOT / "data" / "sharur.duckdb"


@pytest.mark.integration
def test_cli_search_runs_with_db():
    if not DB_PATH.exists():
        pytest.skip("DuckDB not built; run ingest first.")
    cmd = [
        "python",
        "-m",
        "sharur.cli",
        "search",
        "--has",
        "giant",
        "--db",
        str(DB_PATH),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
