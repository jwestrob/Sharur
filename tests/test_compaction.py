"""Exact, reversible physical-copy controls using tiny local fixtures."""

import hashlib
import json
import os
from pathlib import Path

import duckdb
import pytest
from typer.testing import CliRunner

from sharur.cli import app
from sharur.storage import compaction


def file_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture
def source(tmp_path):
    path = tmp_path / "source.duckdb"
    connection = duckdb.connect(str(path), config={"threads": 2, "memory_limit": "128MB"})
    connection.execute('CREATE TABLE "z parent"(id INTEGER PRIMARY KEY)')
    connection.execute(
        'CREATE TABLE "a child"(pid INTEGER REFERENCES "z parent"(id), value VARCHAR DEFAULT \'fixture\')'
    )
    connection.execute('INSERT INTO "z parent" VALUES (1)')
    connection.execute('INSERT INTO "a child"(pid) VALUES (1)')
    connection.execute("CREATE SCHEMA extra")
    connection.execute("CREATE TABLE extra.empty(i INTEGER)")
    connection.execute(
        "CREATE TABLE observations(label VARCHAR, payload JSON, items VARCHAR[], sequence VARCHAR, n DECIMAL(9,3), f DOUBLE)"
    )
    row = (
        "fixture",
        '{ "number": 1.00, "nested": [null, 2] }',
        ["fixture", None],
        "X" * 64,
        "2.500",
        -0.0,
    )
    connection.executemany("INSERT INTO observations VALUES (?,?,?,?,?,?)", [row, row])
    connection.execute("INSERT INTO observations VALUES (NULL,NULL,NULL,NULL,NULL,NULL)")
    connection.execute("CREATE INDEX observed_label ON observations(label)")
    connection.execute(
        "CREATE VIEW z_first AS SELECT label,payload,items,length(sequence) AS residues FROM observations"
    )
    connection.execute("CREATE VIEW a_second AS SELECT * FROM z_first")
    connection.execute("COMMENT ON TABLE observations IS 'table fixture'")
    connection.execute("COMMENT ON COLUMN observations.payload IS 'column fixture'")
    connection.execute("COMMENT ON INDEX observed_label IS 'index fixture'")
    connection.execute("COMMENT ON VIEW z_first IS 'view fixture'")
    connection.close()
    return path


def test_exact_copy_preserves_fk_rows_json_duplicates_comments_and_sidecars(source, tmp_path):
    original_hash = file_hash(source)
    sidecar = tmp_path / "dataset.seal.json"
    sidecar.write_text('{"existing":true}\n')
    sidecar_hash = file_hash(sidecar)
    output, receipt = tmp_path / "compacted.duckdb", tmp_path / "receipt.json"
    report = compaction.compact_database(
        source, output, threads=2, memory_limit="128MB", receipt=receipt
    )
    assert report["status"] == "PASS"
    assert report["source_sha256"] == original_hash == file_hash(source)
    assert file_hash(sidecar) == sidecar_hash
    assert report["output_sha256"] == file_hash(output)
    assert report["table_count"] == 4
    assert report["explicit_index_count"] == 1
    assert report["view_count"] == 2
    assert report["reseal_required"] is True
    assert sum(row["rows"] for row in report["tables"]) == 5
    evidence = receipt.read_text()
    assert str(tmp_path) not in evidence
    assert "number" not in evidence
    assert "X" * 64 not in evidence
    connection = duckdb.connect(str(output), read_only=True)
    assert (
        connection.execute('SELECT COUNT(*) FROM "a child" JOIN "z parent" ON pid=id').fetchone()[0]
        == 1
    )
    assert connection.execute("SELECT COUNT(*) FROM a_second").fetchone()[0] == 3
    assert (
        connection.execute(
            "SELECT CAST(payload AS VARCHAR) FROM observations WHERE payload IS NOT NULL"
        ).fetchall()
        == [('{ "number": 1.00, "nested": [null, 2] }',)] * 2
    )
    assert (
        connection.execute("SELECT items FROM observations WHERE items IS NOT NULL").fetchall()
        == [(["fixture", None],)] * 2
    )
    connection.close()
    assert not list(tmp_path.glob("*.building"))


def test_dependency_ordered_copy_succeeds_where_native_copy_fails(tmp_path):
    source = tmp_path / "native_fk.duckdb"
    setup = duckdb.connect(str(source))
    setup.execute("CREATE TABLE p(id INTEGER PRIMARY KEY)")
    setup.execute("CREATE TABLE c(pid INTEGER REFERENCES p(id))")
    setup.execute("INSERT INTO p VALUES (1)")
    setup.execute("INSERT INTO c VALUES (1)")
    setup.close()
    connection = duckdb.connect(":memory:", config={"threads": 2, "memory_limit": "128MB"})
    connection.execute("ATTACH " + compaction._literal(str(source)) + " AS original (READ_ONLY)")
    connection.execute("ATTACH ':memory:' AS target")
    with pytest.raises(duckdb.ConstraintException):
        connection.execute("COPY FROM DATABASE original TO target")
    connection.close()
    assert (
        compaction.compact_database(source, tmp_path / "ordered.duckdb", memory_limit="128MB")[
            "status"
        ]
        == "PASS"
    )


@pytest.mark.parametrize(
    "mutation",
    [
        'UPDATE observations SET payload=json(\'{"number":1.0,"nested":[null,2]}\') WHERE payload IS NOT NULL',
        "DELETE FROM observations WHERE rowid=0",
        "UPDATE observations SET sequence=repeat('Y',64) WHERE sequence IS NOT NULL",
        "DROP INDEX observed_label",
        "COMMENT ON COLUMN observations.payload IS NULL",
        "DROP VIEW a_second",
    ],
)
def test_corrupted_candidate_never_published(source, tmp_path, monkeypatch, mutation):
    original_copy = compaction._copy
    source_hash = file_hash(source)

    def corrupted(connection, catalog, ordered):
        original_copy(connection, catalog, ordered)
        connection.execute(mutation)

    monkeypatch.setattr(compaction, "_copy", corrupted)
    output, receipt = tmp_path / "rejected.duckdb", tmp_path / "failure.json"
    with pytest.raises(compaction.CompactionError):
        compaction.compact_database(source, output, memory_limit="128MB", receipt=receipt)
    assert not output.exists()
    assert file_hash(source) == source_hash
    evidence = json.loads(receipt.read_text())
    assert evidence["status"] == "FAILED"
    assert evidence["failure_phase"] == "validation"
    assert "X" * 64 not in receipt.read_text()
    assert "Y" * 64 not in receipt.read_text()


@pytest.mark.parametrize(
    "definition",
    [
        "CREATE SEQUENCE local_sequence",
        "CREATE TYPE local_type AS ENUM ('fixture')",
        "CREATE MACRO local_function(x) AS x+1",
    ],
)
def test_unsupported_catalog_fails_before_candidate_creation(source, tmp_path, definition):
    connection = duckdb.connect(str(source))
    connection.execute(definition)
    connection.close()
    source_hash = file_hash(source)
    output = tmp_path / "rejected.duckdb"
    with pytest.raises(compaction.CompactionError, match="preservation adapter"):
        compaction.compact_database(source, output, memory_limit="128MB")
    assert not output.exists()
    assert file_hash(source) == source_hash
    assert not list(tmp_path.glob("*.building"))


@pytest.mark.parametrize(
    "case",
    [
        "same_source",
        "source_wal",
        "existing_output",
        "existing_receipt",
        "output_symlink",
        "receipt_output_wal",
    ],
)
def test_refuses_overwrite_and_aliases(source, tmp_path, case):
    source_hash = file_hash(source)
    output, receipt = tmp_path / "output.duckdb", None
    if case == "same_source":
        output = source
    if case == "source_wal":
        output = Path(str(source) + ".wal")
    if case == "existing_output":
        output.write_text("keep")
    if case == "existing_receipt":
        receipt = tmp_path / "receipt.json"
        receipt.write_text("keep")
    if case == "output_symlink":
        output.symlink_to(source)
    if case == "receipt_output_wal":
        receipt = tmp_path / "output.duckdb.wal"
    with pytest.raises(compaction.CompactionError):
        compaction.compact_database(source, output, receipt=receipt)
    assert file_hash(source) == source_hash
    if case == "existing_output":
        assert output.read_text() == "keep"
    if case == "existing_receipt":
        assert receipt.read_text() == "keep"


@pytest.mark.parametrize(
    "options",
    [
        {"threads": 0},
        {"threads": (os.cpu_count() or 1) + 1},
        {"memory_limit": "0GB"},
        {"memory_limit": "unlimited"},
        {"memory_limit": "1GB';SELECT 1"},
        {"storage_version": "latest"},
    ],
)
def test_invalid_resource_or_storage_setting_creates_nothing(source, tmp_path, options):
    output = tmp_path / "rejected.duckdb"
    with pytest.raises(compaction.CompactionError):
        compaction.compact_database(source, output, **options)
    assert not output.exists()
    assert not list(tmp_path.glob("*.building"))


def test_exclusive_publication_preserves_destination_created_during_copy(
    source, tmp_path, monkeypatch
):
    output = tmp_path / "occupied.duckdb"
    original_copy = compaction._copy

    def concurrent_creation(connection, catalog, ordered):
        original_copy(connection, catalog, ordered)
        output.write_text("other writer")

    monkeypatch.setattr(compaction, "_copy", concurrent_creation)
    with pytest.raises(compaction.CompactionError, match="publication"):
        compaction.compact_database(source, output, memory_limit="128MB")
    assert output.read_text() == "other writer"


def test_cli_requires_fresh_output_and_emits_reseal_contract(source, tmp_path):
    runner = CliRunner()
    assert runner.invoke(app, ["compact", "--db", str(source)]).exit_code == 2
    output = tmp_path / "cli.duckdb"
    result = runner.invoke(
        app,
        [
            "compact",
            "--db",
            str(source),
            "--output",
            str(output),
            "--memory-limit",
            "128MB",
            "--storage-version",
            "v1.0.0",
        ],
    )
    assert result.exit_code == 0, result.output
    assert output.is_file()
    assert "Rebuild the content seal" in result.output


def test_cleanup_failure_preserves_successful_publication(source, tmp_path, monkeypatch):
    original_rmdir = Path.rmdir

    def fail_owned_cleanup(path):
        if path.name.endswith(".building"):
            raise OSError("fixture cleanup failure")
        return original_rmdir(path)

    monkeypatch.setattr(Path, "rmdir", fail_owned_cleanup)
    output = tmp_path / "validated.duckdb"
    report = compaction.compact_database(source, output, memory_limit="128MB")
    assert report["status"] == "PASS"
    assert output.is_file()
    assert report["output_sha256"] == file_hash(output)
    assert report["cleanup_warning"]["task_owned_temporary_artifacts_retained"] is True


def test_interrupted_copy_keeps_source_and_records_failure(source, tmp_path, monkeypatch):
    source_hash = file_hash(source)

    def interrupted(*args):
        raise KeyboardInterrupt

    monkeypatch.setattr(compaction, "_copy", interrupted)
    output, receipt = tmp_path / "interrupted.duckdb", tmp_path / "interrupted.json"
    with pytest.raises(KeyboardInterrupt):
        compaction.compact_database(source, output, memory_limit="128MB", receipt=receipt)
    assert not output.exists()
    assert file_hash(source) == source_hash
    assert json.loads(receipt.read_text())["error_type"] == "KeyboardInterrupt"


def test_source_change_prevents_publication(source, tmp_path, monkeypatch):
    original_copy = compaction._copy

    def source_touched(connection, catalog, ordered):
        original_copy(connection, catalog, ordered)
        stat = source.stat()
        os.utime(source, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1))

    monkeypatch.setattr(compaction, "_copy", source_touched)
    output = tmp_path / "changed_source.duckdb"
    with pytest.raises(compaction.CompactionError, match="Source changed"):
        compaction.compact_database(source, output, memory_limit="128MB")
    assert not output.exists()


def test_column_collation_mutation_is_rejected(source, tmp_path, monkeypatch):
    original_ddl = compaction._table_ddl

    def changed_collation(schema, entry):
        sql = original_ddl(schema, entry)
        if '"label" VARCHAR' in sql:
            return sql.replace('"label" VARCHAR', '"label" VARCHAR COLLATE nocase')
        return sql

    monkeypatch.setattr(compaction, "_table_ddl", changed_collation)
    output = tmp_path / "wrong_collation.duckdb"
    with pytest.raises(compaction.CompactionError, match="catalog definitions"):
        compaction.compact_database(source, output, memory_limit="128MB")
    assert not output.exists()
