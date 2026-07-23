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
    assert payload["trace"]["schema_version"] == "6"


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


def test_inspect_json_preserves_asymmetric_context(case_database):
    result = CliRunner().invoke(
        app,
        [
            "inspect",
            "target_minus",
            "--type",
            "system",
            "--window",
            "8",
            "--upstream",
            "2",
            "--downstream",
            "1",
            "--db",
            str(case_database),
            "--format",
            "json",
        ],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["entity"]["source_table"] == "defense_systems"
    assert payload["context_window"]["default_orfs"] == 8
    assert payload["context_window"]["upstream_orfs"] == 2
    assert payload["context_window"]["downstream_orfs"] == 1
    assert payload["context_window"]["orientation"] == "-"


def test_compare_context_json_exposes_denominators_and_recipe(case_database):
    result = CliRunner().invoke(
        app,
        [
            "compare-context",
            "target_plus",
            "--type",
            "system",
            "--feature",
            "pfam:PFTEST",
            "--window",
            "2",
            "--min-components",
            "2",
            "--db",
            str(case_database),
            "--format",
            "json",
        ],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["composite"]["foreground_positive"] == 2
    assert payload["composite"]["foreground_total"] == 2
    assert payload["composite"]["background_positive"] == 1
    assert payload["composite"]["background_total"] == 2
    assert payload["recipe"]["upstream_orfs"] == 2
    assert payload["recipe"]["downstream_orfs"] == 2
