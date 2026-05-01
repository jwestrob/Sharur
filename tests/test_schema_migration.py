"""Tests for schema versioning and migrations."""

import duckdb
import pytest

from sharur.storage.migrations import get_current_version, run_migrations
from sharur.storage.schema import SCHEMA, SCHEMA_VERSION


@pytest.fixture
def fresh_conn():
    """In-memory DuckDB connection with schema applied."""
    conn = duckdb.connect(":memory:")
    conn.execute(SCHEMA)
    return conn


@pytest.fixture
def bare_conn():
    """In-memory DuckDB connection with schema but no migrations."""
    conn = duckdb.connect(":memory:")
    conn.execute(SCHEMA)
    return conn


def test_fresh_db_gets_latest_version(fresh_conn):
    """A fresh database should reach the latest version after migrations."""
    assert get_current_version(fresh_conn) == 0
    applied = run_migrations(fresh_conn)
    assert applied == SCHEMA_VERSION
    assert get_current_version(fresh_conn) == SCHEMA_VERSION


def test_db_without_version_table_returns_0():
    """A database without schema_version table returns version 0."""
    conn = duckdb.connect(":memory:")
    assert get_current_version(conn) == 0


def test_migrations_idempotent(fresh_conn):
    """Running migrations twice applies nothing the second time."""
    run_migrations(fresh_conn)
    assert get_current_version(fresh_conn) == SCHEMA_VERSION

    applied = run_migrations(fresh_conn)
    assert applied == 0
    assert get_current_version(fresh_conn) == SCHEMA_VERSION


def test_schema_version_table_has_correct_data(fresh_conn):
    """After migration, schema_version table has the expected row."""
    run_migrations(fresh_conn)
    row = fresh_conn.execute(
        "SELECT version, description FROM schema_version WHERE version = 1"
    ).fetchone()
    assert row is not None
    assert row[0] == 1
    assert "Initial" in row[1]


def test_schema_version_constant():
    """SCHEMA_VERSION constant matches latest migration."""
    assert SCHEMA_VERSION == 4


def test_duckdb_store_runs_migrations():
    """DuckDBStore._initialize_schema runs migrations on fresh DB."""
    from sharur.storage.duckdb_store import DuckDBStore

    store = DuckDBStore()  # in-memory
    # Accessing .conn triggers _initialize_schema
    _ = store.conn
    version = get_current_version(store._conn)
    assert version == SCHEMA_VERSION


def test_duckdb_store_read_only_skips_schema_initialization(tmp_path):
    """Read-only stores should not require a writer lock for schema setup."""
    from sharur.storage.duckdb_store import DuckDBStore

    db_path = tmp_path / "readonly.duckdb"
    writer = DuckDBStore(db_path)
    writer.execute("INSERT INTO bins (bin_id) VALUES ('bin1')")
    writer.conn.close()

    store = DuckDBStore(db_path, read_only=True)

    assert store.read_only is True
    assert store.execute("SELECT COUNT(*) FROM bins")[0][0] == 1
    with pytest.raises(duckdb.InvalidInputException, match="read-only"):
        store.execute("INSERT INTO bins (bin_id) VALUES ('bin2')")


def test_existing_v1_db_gets_current_semantic_tables():
    """Pending migrations should run on DBs already marked version 1."""
    conn = duckdb.connect(":memory:")
    conn.execute("""
        CREATE TABLE schema_version (
            version INTEGER PRIMARY KEY,
            applied_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            description TEXT
        );
        INSERT INTO schema_version (version, description)
            VALUES (1, 'Initial schema version tracking');
        CREATE TABLE proteins (protein_id VARCHAR PRIMARY KEY);
    """)

    applied = run_migrations(conn)
    assert applied == SCHEMA_VERSION - 1
    assert get_current_version(conn) == SCHEMA_VERSION
    assert conn.execute(
        "SELECT COUNT(*) FROM information_schema.tables "
        "WHERE table_name IN ("
        "  'semantic_atoms', 'semantic_state', 'semantic_terms', 'system_proteins'"
        ")"
    ).fetchone()[0] == 4
