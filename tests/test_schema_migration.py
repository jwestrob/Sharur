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


def test_fresh_db_gets_version_1(fresh_conn):
    """A fresh database should reach version 1 after migrations."""
    assert get_current_version(fresh_conn) == 0
    applied = run_migrations(fresh_conn)
    assert applied == 1
    assert get_current_version(fresh_conn) == 1


def test_db_without_version_table_returns_0():
    """A database without schema_version table returns version 0."""
    conn = duckdb.connect(":memory:")
    assert get_current_version(conn) == 0


def test_migrations_idempotent(fresh_conn):
    """Running migrations twice applies nothing the second time."""
    run_migrations(fresh_conn)
    assert get_current_version(fresh_conn) == 1

    applied = run_migrations(fresh_conn)
    assert applied == 0
    assert get_current_version(fresh_conn) == 1


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
    assert SCHEMA_VERSION == 1


def test_duckdb_store_runs_migrations():
    """DuckDBStore._initialize_schema runs migrations on fresh DB."""
    from sharur.storage.duckdb_store import DuckDBStore

    store = DuckDBStore()  # in-memory
    # Accessing .conn triggers _initialize_schema
    _ = store.conn
    version = get_current_version(store._conn)
    assert version == 1
