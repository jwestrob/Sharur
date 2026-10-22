"""Tests for schema versioning and migrations."""

import duckdb
import pytest

from sharur.storage.duckdb_store import DuckDBStore
from sharur.storage.migrations import MIGRATIONS, get_current_version, run_migrations
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
    assert fresh_conn.execute(
        """
        SELECT COUNT(*)
        FROM duckdb_indexes()
        WHERE index_name = 'idx_contigs_bin'
        """
    ).fetchone()[0] == 1


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
    assert SCHEMA_VERSION == max(version for version, _, _ in MIGRATIONS) == 9


def test_duckdb_store_runs_migrations():
    """DuckDBStore._initialize_schema runs migrations on fresh DB."""
    store = DuckDBStore()  # in-memory
    # Accessing .conn triggers _initialize_schema
    _ = store.conn
    version = get_current_version(store._conn)
    assert version == SCHEMA_VERSION


def test_duckdb_store_read_only_skips_schema_initialization(tmp_path):
    """Read-only stores should not require a writer lock for schema setup."""
    db_path = tmp_path / "readonly.duckdb"
    writer = DuckDBStore(db_path)
    writer.execute("INSERT INTO bins (bin_id) VALUES ('bin1')")
    writer.conn.close()

    store = DuckDBStore(db_path, read_only=True)

    assert store.read_only is True
    assert store.execute("SELECT COUNT(*) FROM bins")[0][0] == 1
    with pytest.raises(duckdb.InvalidInputException, match="read-only"):
        store.execute("INSERT INTO bins (bin_id) VALUES ('bin2')")


def _write_pre_v5_db(db_path) -> None:
    """A file shaped like a v4 build: system tables without contig_id."""
    conn = duckdb.connect(str(db_path))
    conn.execute(SCHEMA)
    run_migrations(conn, target_version=4)
    for table in ("defense_systems", "secretion_systems"):
        for suffix in ("genome", "replicon", "type"):
            conn.execute(f"DROP INDEX IF EXISTS idx_{table}_{suffix}")
        conn.execute(f"ALTER TABLE {table} DROP COLUMN contig_id")
        conn.execute(f"CREATE INDEX idx_{table}_genome ON {table}(genome_id)")
        conn.execute(f"CREATE INDEX idx_{table}_type ON {table}(system_type)")
    conn.execute("""
        INSERT INTO bins (bin_id) VALUES ('bin1');
        INSERT INTO contigs (contig_id, bin_id, length)
            VALUES ('contig1', 'bin1', 300);
        INSERT INTO proteins (
            protein_id, contig_id, bin_id, start, end_coord, strand
        ) VALUES ('p1', 'contig1', 'bin1', 1, 300, '+');
        INSERT INTO annotations (
            annotation_id, protein_id, source, accession, name
        ) VALUES
            (1, 'p1', 'defensefinder', 'raw-profile', 'raw-profile'),
            (2, 'p1', 'defensefinder_system', 'legacy-profile', 'legacy-system');
        INSERT INTO defense_systems (
            system_id, genome_id, system_type, system_subtype,
            genes_count, protein_ids, profile_names
        ) VALUES (
            'legacy-system', 'bin1', 'RM', 'RM_Type_IV',
            1, 'p1', 'legacy-profile'
        );
    """)
    conn.close()


def test_duckdb_store_opens_pre_v5_file(tmp_path):
    """A writable store migrates a v4 file before applying SCHEMA's indexes."""
    db_path = tmp_path / "pre_v5.duckdb"
    _write_pre_v5_db(db_path)

    store = DuckDBStore(db_path)
    conn = store.conn

    assert get_current_version(conn) == SCHEMA_VERSION
    for table in ("defense_systems", "secretion_systems"):
        columns = {row[0] for row in conn.execute(f"DESCRIBE {table}").fetchall()}
        assert "contig_id" in columns
        assert conn.execute(
            "SELECT COUNT(*) FROM duckdb_indexes() WHERE index_name = ?",
            [f"idx_{table}_replicon"],
        ).fetchone()[0] == 1
    # v5 quarantine: legacy named calls removed, raw HMM observations kept.
    assert conn.execute("SELECT COUNT(*) FROM defense_systems").fetchone()[0] == 0
    assert conn.execute(
        "SELECT source FROM annotations ORDER BY annotation_id"
    ).fetchall() == [("defensefinder",)]
    store.close()


def test_duckdb_store_reopens_current_file_without_reapplying(tmp_path):
    """Reopening a current-version file applies no further migrations."""
    db_path = tmp_path / "current.duckdb"
    DuckDBStore(db_path).conn.execute("INSERT INTO bins (bin_id) VALUES ('b')")
    store = DuckDBStore(db_path)
    conn = store.conn
    assert get_current_version(conn) == SCHEMA_VERSION
    assert conn.execute(
        "SELECT COUNT(*) FROM schema_version"
    ).fetchone()[0] == SCHEMA_VERSION
    assert run_migrations(conn) == 0
    store.close()


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


def test_v5_migration_quarantines_legacy_named_system_calls():
    """Pre-v5 caller output must not survive as validation-level evidence."""
    conn = duckdb.connect(":memory:")
    conn.execute(SCHEMA)
    run_migrations(conn, target_version=4)
    conn.execute("""
        INSERT INTO bins (bin_id) VALUES ('bin1');
        INSERT INTO contigs (contig_id, bin_id, length)
            VALUES ('contig1', 'bin1', 300);
        INSERT INTO proteins (
            protein_id, contig_id, bin_id, start, end_coord, strand
        ) VALUES ('p1', 'contig1', 'bin1', 1, 300, '+');
        INSERT INTO annotations (
            annotation_id, protein_id, source, accession, name
        ) VALUES
            (1, 'p1', 'txsscan', 'raw-profile', 'raw-profile'),
            (2, 'p1', 'txsscan_system', 'legacy-profile', 'legacy-system');
        INSERT INTO secretion_systems (
            system_id, genome_id, system_type, system_subtype,
            genes_count, protein_ids, profile_names
        ) VALUES (
            'legacy-system', 'bin1', 'T4P', 'T4P',
            1, 'p1', 'legacy-profile'
        );
        INSERT INTO system_proteins (
            system_id, protein_id, system_source, position, profile_name
        ) VALUES (
            'legacy-system', 'p1', 'txsscan_system', 1, 'legacy-profile'
        );
        INSERT INTO semantic_atoms (
            protein_id, atom_id, facet, relation, source_accession, source_db
        ) VALUES
            ('p1', 'secretion_system', 'role', 'implies',
             'legacy-system', 'txsscan_system'),
            ('p1', 'raw_observation', 'role', 'supports',
             'raw-profile', 'txsscan');
        INSERT INTO semantic_state (
            protein_id, activities, roles, architecture, localization,
            topology, size_class, quality_flags, composite_predicates,
            unresolved_count
        ) VALUES (
            'p1', [], ['secretion_system'], [], [], '{}', 'standard',
            [], ['secretion_system_validated'], 0
        );
        INSERT INTO semantic_terms (
            protein_id, term_id, term_kind, source_db, source_accession
        ) VALUES
            ('p1', 'secretion_system_validated', 'composite', '', ''),
            ('p1', 'secretion_system', 'atom', 'txsscan_system', 'legacy-system');
        INSERT INTO protein_predicates (protein_id, predicates)
            VALUES ('p1', ['secretion_system_validated']);
    """)

    applied = run_migrations(conn)

    assert applied == SCHEMA_VERSION - 4
    assert get_current_version(conn) == SCHEMA_VERSION
    assert conn.execute("SELECT COUNT(*) FROM secretion_systems").fetchone()[0] == 0
    assert conn.execute(
        "SELECT source FROM annotations ORDER BY annotation_id"
    ).fetchall() == [("txsscan",)]
    assert conn.execute("SELECT COUNT(*) FROM system_proteins").fetchone()[0] == 0
    assert conn.execute(
        "SELECT atom_id, source_db FROM semantic_atoms"
    ).fetchall() == [("raw_observation", "txsscan")]
    assert conn.execute("SELECT COUNT(*) FROM semantic_state").fetchone()[0] == 0
    assert conn.execute("SELECT COUNT(*) FROM semantic_terms").fetchone()[0] == 0
    assert conn.execute("SELECT COUNT(*) FROM protein_predicates").fetchone()[0] == 0
    assert {
        row[0] for row in conn.execute("DESCRIBE secretion_systems").fetchall()
    } >= {"contig_id"}
    assert conn.execute(
        """
        SELECT COUNT(*)
        FROM duckdb_indexes()
        WHERE index_name = 'idx_contigs_bin'
        """
    ).fetchone()[0] == 1
