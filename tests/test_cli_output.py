"""Machine-readable output tests for read-only CLI commands."""

import json

from typer.testing import CliRunner

from sharur.cli import app
from sharur.storage.duckdb_store import DuckDBStore


def _make_cli_database(tmp_path):
    db_path = tmp_path / "sharur.duckdb"
    store = DuckDBStore(db_path)
    store.execute(
        """
        INSERT INTO bins (
            bin_id, completeness, contamination, taxonomy, n_contigs, total_length
        )
        VALUES ('bin_1', 90, 1, 'd__Bacteria;p__Test', 1, 1000)
        """
    )
    store.conn.close()
    return db_path


def test_overview_json_output_is_structured(tmp_path):
    db_path = _make_cli_database(tmp_path)

    result = CliRunner().invoke(
        app,
        ["overview", "--db", str(db_path), "--format", "json"],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["status"] == "ok"
    assert payload["raw"]["genome_count"] == 1
    assert payload["trace"]["schema_version"] == "4"


def test_genomes_jsonl_output_emits_one_record_per_row(tmp_path):
    db_path = _make_cli_database(tmp_path)

    result = CliRunner().invoke(
        app,
        ["genomes", "--db", str(db_path), "--format", "jsonl"],
    )

    assert result.exit_code == 0, result.output
    records = [json.loads(line) for line in result.output.splitlines()]
    assert records == [
        {
            "bin_id": "bin_1",
            "completeness": 90.0,
            "contamination": 1.0,
            "taxonomy": "d__Bacteria;p__Test",
            "n_contigs": 1,
            "total_length": 1000,
        }
    ]
