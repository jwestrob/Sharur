"""Schema versioning and migrations for Sharur databases."""

from __future__ import annotations

from contextlib import suppress
from dataclasses import dataclass
from typing import TYPE_CHECKING, NamedTuple

from sharur.predicates.provenance import CREATE_SQL as PROVENANCE_SQL


if TYPE_CHECKING:
    from collections.abc import Callable

    import duckdb


@dataclass(frozen=True)
class Destructive:
    """Marks a migration that removes or rewrites rows already in a database.

    ``run_migrations`` applies such a migration only with
    ``allow_destructive=True``, which ``sharur migrate`` passes; every other
    caller (a writable ``DuckDBStore`` open, Stage 07, the system callers, the
    contig backfill) stops with :class:`DestructiveMigrationError` instead.
    ``removes_rows`` reports whether applying the migration to this database
    would remove anything, so a fresh or empty database still migrates on
    open; without one, the migration counts as destructive for every database.
    """

    effect: str
    removes_rows: Callable[[duckdb.DuckDBPyConnection], bool] | None = None


class Migration(NamedTuple):
    version: int
    description: str
    sql: str
    # None for an additive migration (new tables, columns, indexes). Mark every
    # migration that deletes or rewrites existing rows.
    destructive: Destructive | None = None


MIGRATIONS: list[Migration] = [
    Migration(
        1,
        "Initial schema version tracking",
        """
        CREATE TABLE IF NOT EXISTS schema_version (
            version INTEGER PRIMARY KEY,
            applied_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            description TEXT
        );
        INSERT INTO schema_version (version, description)
            VALUES (1, 'Initial schema version tracking');
        """,
    ),
    Migration(
        2,
        "Predicate V2 semantic tables and indexes",
        """
        CREATE TABLE IF NOT EXISTS semantic_atoms (
            protein_id VARCHAR NOT NULL,
            atom_id VARCHAR NOT NULL,
            facet VARCHAR NOT NULL,
            relation VARCHAR NOT NULL,
            source_accession VARCHAR,
            source_db VARCHAR,
            evidence_evalue DOUBLE,
            evidence_score DOUBLE,

            PRIMARY KEY (protein_id, atom_id, source_accession)
        );
        CREATE INDEX IF NOT EXISTS idx_semantic_atoms_protein
            ON semantic_atoms(protein_id);
        CREATE INDEX IF NOT EXISTS idx_semantic_atoms_atom
            ON semantic_atoms(atom_id);
        CREATE INDEX IF NOT EXISTS idx_semantic_atoms_facet_atom
            ON semantic_atoms(facet, atom_id);
        CREATE INDEX IF NOT EXISTS idx_semantic_atoms_relation
            ON semantic_atoms(relation);
        CREATE INDEX IF NOT EXISTS idx_semantic_atoms_source
            ON semantic_atoms(source_db, source_accession);

        CREATE TABLE IF NOT EXISTS semantic_state (
            protein_id VARCHAR PRIMARY KEY,
            activities VARCHAR[],
            roles VARCHAR[],
            architecture VARCHAR[],
            localization VARCHAR[],
            topology JSON,
            size_class VARCHAR,
            quality_flags VARCHAR[],
            composite_predicates VARCHAR[],
            unresolved_count INTEGER,
            updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        );
        CREATE INDEX IF NOT EXISTS idx_semantic_state_size
            ON semantic_state(size_class);
        CREATE INDEX IF NOT EXISTS idx_semantic_state_updated
            ON semantic_state(updated_at);

        INSERT INTO schema_version (version, description)
            VALUES (2, 'Predicate V2 semantic tables and indexes');
        """,
    ),
    Migration(
        3,
        "Predicate V2 materialized search terms",
        """
        CREATE TABLE IF NOT EXISTS semantic_terms (
            protein_id VARCHAR NOT NULL,
            term_id VARCHAR NOT NULL,
            term_kind VARCHAR NOT NULL,
            facet VARCHAR,
            relation VARCHAR,
            source_db VARCHAR NOT NULL DEFAULT '',
            source_accession VARCHAR NOT NULL DEFAULT ''
        );
        CREATE INDEX IF NOT EXISTS idx_semantic_terms_term
            ON semantic_terms(term_id, protein_id);
        CREATE INDEX IF NOT EXISTS idx_semantic_terms_protein
            ON semantic_terms(protein_id);
        CREATE INDEX IF NOT EXISTS idx_semantic_terms_facet_term
            ON semantic_terms(facet, term_id);
        CREATE INDEX IF NOT EXISTS idx_semantic_terms_source
            ON semantic_terms(source_db, source_accession);
        CREATE INDEX IF NOT EXISTS idx_semantic_terms_kind
            ON semantic_terms(term_kind);

        INSERT INTO schema_version (version, description)
            VALUES (3, 'Predicate V2 materialized search terms');
        """,
    ),
    Migration(
        4,
        "Validated system membership join table",
        """
        CREATE TABLE IF NOT EXISTS system_proteins (
            system_id VARCHAR NOT NULL,
            protein_id VARCHAR NOT NULL,
            system_source VARCHAR NOT NULL,
            position INTEGER,
            profile_name VARCHAR,
            score DOUBLE,

            PRIMARY KEY (system_id, protein_id, system_source),
            FOREIGN KEY (protein_id) REFERENCES proteins(protein_id)
        );
        CREATE INDEX IF NOT EXISTS idx_system_proteins_protein
            ON system_proteins(protein_id);
        CREATE INDEX IF NOT EXISTS idx_system_proteins_system
            ON system_proteins(system_id);
        CREATE INDEX IF NOT EXISTS idx_system_proteins_source
            ON system_proteins(system_source);

        INSERT INTO schema_version (version, description)
            VALUES (4, 'Validated system membership join table');
        """,
    ),
    Migration(
        5,
        "Replicon provenance for validated systems",
        """
        INSERT INTO schema_version (version, description)
            VALUES (5, 'Replicon provenance for validated systems');
        """,
        Destructive(
            effect=(
                "removes every DefenseFinder/TXSScan system call made before replicon-local "
                "validation (defense_systems, secretion_systems, *_system annotations and "
                "system_proteins rows) and the V2 rows of their member proteins; raw HMM hits "
                "stay, and Stage 07 (or the callers plus a V2 refresh) rebuilds the calls"
            ),
            removes_rows=lambda conn: _pre_v5_system_calls_present(conn),
        ),
    ),
    Migration(
        6,
        "Contig navigation index",
        """
        INSERT INTO schema_version (version, description)
            VALUES (6, 'Contig navigation index');
        """,
    ),
    Migration(
        7,
        "Predicate map provenance stamps",
        f"""
        {PROVENANCE_SQL};
        INSERT INTO schema_version (version, description)
            VALUES (7, 'Predicate map provenance stamps');
        """,
    ),
    Migration(
        8,
        "Contig length source and gene truncation flags",
        """
        INSERT INTO schema_version (version, description)
            VALUES (8, 'Contig length source and gene truncation flags');
        """,
    ),
    Migration(
        9,
        "CRISPR-Cas subtype calls",
        """
        CREATE TABLE IF NOT EXISTS crispr_cas_systems (
            system_id VARCHAR PRIMARY KEY,
            genome_id VARCHAR,
            contig_id VARCHAR,
            start INTEGER,
            end_coord INTEGER,
            status VARCHAR,                 -- crispr_cas | crispr_cas_putative | cas | cas_putative
            prediction VARCHAR,             -- joint operon + array call
            prediction_cas VARCHAR,         -- operon call (subtype, Ambiguous, Hybrid(...), False)
            best_type VARCHAR,              -- top-scoring subtype(s), comma-joined
            best_score DOUBLE,
            complete_interference VARCHAR,
            complete_adaptation VARCHAR,
            strand_interference VARCHAR,
            strand_adaptation VARCHAR,
            genes_count INTEGER,
            protein_ids VARCHAR,
            profile_names VARCHAR,
            crispr_locus_ids VARCHAR,
            crispr_distances VARCHAR,
            caller VARCHAR,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        );
        CREATE INDEX IF NOT EXISTS idx_crispr_cas_systems_genome ON crispr_cas_systems(genome_id);

        -- Repeat-subtype predictions and trust statistics for CRISPR array loci
        CREATE TABLE IF NOT EXISTS crispr_array_types (
            locus_id VARCHAR PRIMARY KEY,
            genome_id VARCHAR,
            contig_id VARCHAR,
            repeat_length INTEGER,
            n_repeats INTEGER,
            subtype VARCHAR,                -- most probable repeat subtype
            probability DOUBLE,
            prediction VARCHAR,             -- subtype, or Unknown below the probability cutoff
            repeat_identity DOUBLE,
            spacer_identity DOUBLE,
            spacer_sem DOUBLE,
            trusted BOOLEAN,
            near_cas BOOLEAN,
            caller VARCHAR
        );
        INSERT INTO schema_version (version, description)
            VALUES (9, 'CRISPR-Cas subtype calls');
        """,
    ),
]


def _column_names(conn: duckdb.DuckDBPyConnection, table_name: str) -> set[str]:
    """Return columns for an existing table."""
    return {
        str(row[1])
        for row in conn.execute(f"PRAGMA table_info('{table_name}')").fetchall()
    }


def _add_v6_contig_index(conn: duckdb.DuckDBPyConnection) -> None:
    """Add the navigation index when the legacy database has contigs."""
    tables = {str(row[0]) for row in conn.execute("SHOW TABLES").fetchall()}
    if "contigs" in tables:
        conn.execute(
            """
            CREATE INDEX IF NOT EXISTS idx_contigs_bin
                ON contigs(bin_id, contig_id)
            """
        )


def _add_v8_contig_context_columns(conn: duckdb.DuckDBPyConnection) -> None:
    """Add contig length provenance and Prodigal truncation flags."""
    tables = {str(row[0]) for row in conn.execute("SHOW TABLES").fetchall()}
    if "contigs" in tables and "length_source" not in _column_names(conn, "contigs"):
        conn.execute("ALTER TABLE contigs ADD COLUMN length_source VARCHAR")
    if "proteins" in tables and "partial" not in _column_names(conn, "proteins"):
        conn.execute("ALTER TABLE proteins ADD COLUMN partial VARCHAR(2)")


_PRE_V5_SYSTEM_SOURCES = ("defensefinder_system", "txsscan_system")
# Tables v5 empties entirely, and (table, column) pairs it filters on the sources above.
_PRE_V5_SYSTEM_TABLES = ("defense_systems", "secretion_systems")
_PRE_V5_SOURCE_COLUMNS = (
    ("annotations", "source"),
    ("system_proteins", "system_source"),
    ("semantic_atoms", "source_db"),
    ("semantic_terms", "source_db"),
)


def _pre_v5_system_calls_present(conn: duckdb.DuckDBPyConnection) -> bool:
    """True when :func:`_quarantine_pre_v5_system_calls` would delete any row.

    It deletes all rows of the two system tables and rows carrying a
    ``*_system`` source; the semantic_state, semantic_terms and
    protein_predicates rows it drops belong to proteins those rows name, so
    they exist only when one of these probes finds a row.
    """
    tables = {str(row[0]) for row in conn.execute("SHOW TABLES").fetchall()}
    for table_name in _PRE_V5_SYSTEM_TABLES:
        if table_name in tables and conn.execute(
            f"SELECT 1 FROM {table_name} LIMIT 1"
        ).fetchone():
            return True
    for table_name, column in _PRE_V5_SOURCE_COLUMNS:
        if (
            table_name in tables
            and column in _column_names(conn, table_name)
            and conn.execute(
                f"SELECT 1 FROM {table_name} WHERE {column} IN (?, ?) LIMIT 1",
                list(_PRE_V5_SYSTEM_SOURCES),
            ).fetchone()
        ):
            return True
    return False


def _quarantine_pre_v5_system_calls(conn: duckdb.DuckDBPyConnection) -> None:
    """Remove caller output that predates replicon-local validation.

    The old caller cannot be repaired by deriving a contig from its serialized
    members: global profile competition and quorum evaluation must be rerun.
    Preserve raw HMM observations, but remove named calls and invalidate V2
    rows for proteins that inherited those calls. This makes a migrated
    database fail closed until Stage 07 reruns the current caller and refreshes
    the affected semantic subset.
    """
    tables = {str(row[0]) for row in conn.execute("SHOW TABLES").fetchall()}
    sources = _PRE_V5_SYSTEM_SOURCES

    conn.execute("""
        CREATE TEMP TABLE sharur_stale_system_members (
            protein_id VARCHAR
        )
    """)

    if "system_proteins" in tables:
        columns = _column_names(conn, "system_proteins")
        if {"protein_id", "system_source"} <= columns:
            conn.execute(
                """
                INSERT INTO sharur_stale_system_members
                SELECT DISTINCT protein_id
                FROM system_proteins
                WHERE system_source IN (?, ?)
                """,
                list(sources),
            )

    if "annotations" in tables:
        columns = _column_names(conn, "annotations")
        if {"protein_id", "source"} <= columns:
            conn.execute(
                """
                INSERT INTO sharur_stale_system_members
                SELECT DISTINCT protein_id
                FROM annotations
                WHERE source IN (?, ?)
                """,
                list(sources),
            )

    if "semantic_atoms" in tables:
        columns = _column_names(conn, "semantic_atoms")
        if {"protein_id", "source_db"} <= columns:
            conn.execute(
                """
                INSERT INTO sharur_stale_system_members
                SELECT DISTINCT protein_id
                FROM semantic_atoms
                WHERE source_db IN (?, ?)
                """,
                list(sources),
            )

    if "semantic_terms" in tables:
        columns = _column_names(conn, "semantic_terms")
        if {"protein_id", "source_db"} <= columns:
            conn.execute(
                """
                INSERT INTO sharur_stale_system_members
                SELECT DISTINCT protein_id
                FROM semantic_terms
                WHERE source_db IN (?, ?)
                """,
                list(sources),
            )

    if "semantic_atoms" in tables:
        conn.execute(
            "DELETE FROM semantic_atoms WHERE source_db IN (?, ?)",
            list(sources),
        )
    if "semantic_state" in tables:
        conn.execute(
            """
            DELETE FROM semantic_state
            WHERE protein_id IN (
                SELECT protein_id FROM sharur_stale_system_members
            )
            """
        )
    if "semantic_terms" in tables:
        conn.execute(
            """
            DELETE FROM semantic_terms
            WHERE protein_id IN (
                SELECT protein_id FROM sharur_stale_system_members
            )
            """
        )
    if "protein_predicates" in tables:
        conn.execute(
            """
            DELETE FROM protein_predicates
            WHERE protein_id IN (
                SELECT protein_id FROM sharur_stale_system_members
            )
            """
        )
    if "annotations" in tables:
        conn.execute(
            "DELETE FROM annotations WHERE source IN (?, ?)",
            list(sources),
        )
    if "system_proteins" in tables:
        conn.execute(
            "DELETE FROM system_proteins WHERE system_source IN (?, ?)",
            list(sources),
        )
    if "defense_systems" in tables:
        conn.execute("DELETE FROM defense_systems")
    if "secretion_systems" in tables:
        conn.execute("DELETE FROM secretion_systems")

    conn.execute("DROP TABLE sharur_stale_system_members")


def _add_v5_replicon_columns(conn: duckdb.DuckDBPyConnection) -> None:
    """Idempotently add caller provenance columns before data migration.

    Current DuckDB releases reject altering and deleting from the same table
    in one transaction. Adding a nullable column is safe as a pre-step: if the
    later quarantine transaction fails, the database remains at schema version
    4 and the whole migration can be retried.
    """
    tables = {str(row[0]) for row in conn.execute("SHOW TABLES").fetchall()}
    for table_name in ("defense_systems", "secretion_systems"):
        if (
            table_name in tables
            and "contig_id" not in _column_names(conn, table_name)
        ):
            conn.execute(f"ALTER TABLE {table_name} ADD COLUMN contig_id VARCHAR")


def get_current_version(conn: duckdb.DuckDBPyConnection) -> int:
    """Return current schema version, or 0 if no version table exists."""
    try:
        result = conn.execute(
            "SELECT MAX(version) FROM schema_version"
        ).fetchone()
        return result[0] if result and result[0] is not None else 0
    except Exception:
        # Table doesn't exist
        return 0


class DestructiveMigrationError(RuntimeError):
    """A pending migration would remove rows; only ``sharur migrate`` applies it."""

    def __init__(self, path: str | None, current: int, migrations: list[Migration]):
        self.path = path
        self.current = current
        self.migrations = migrations
        target = path or "PATH"
        steps = "; ".join(
            f"v{m.version} ({m.description}) {m.destructive.effect}" for m in migrations
        )
        super().__init__(
            f"{target} is at schema {current} and needs a destructive migration: {steps}. "
            "Opening or writing a database applies only additive migrations; nothing was "
            "changed. Copy the file if its current contents matter, then run "
            f"`sharur migrate --db {target}` in a maintenance window."
        )


def _pending(current: int, target_version: int | None) -> list[Migration]:
    if target_version is None:
        target_version = max(m.version for m in MIGRATIONS) if MIGRATIONS else 0
    return [
        m for m in sorted(MIGRATIONS, key=lambda m: m.version)
        if current < m.version <= target_version
    ]


def blocking_migrations(
    conn: duckdb.DuckDBPyConnection, target_version: int | None = None
) -> list[Migration]:
    """Pending destructive migrations that would remove rows from this database.

    Read-only: safe on a read-only connection.
    """
    return [
        m for m in _pending(get_current_version(conn), target_version)
        if m.destructive is not None
        and (m.destructive.removes_rows is None or m.destructive.removes_rows(conn))
    ]


def _database_path(conn: duckdb.DuckDBPyConnection) -> str | None:
    with suppress(Exception):
        row = conn.execute(
            "SELECT path FROM duckdb_databases() WHERE database_name = current_database()"
        ).fetchone()
        return row[0] if row else None
    return None


def run_migrations(
    conn: duckdb.DuckDBPyConnection,
    target_version: int | None = None,
    *,
    allow_destructive: bool = False,
) -> int:
    """Apply pending migrations up to target_version (default: latest).

    Additive migrations always apply. When a pending destructive migration
    would remove rows, nothing is applied and :class:`DestructiveMigrationError`
    names ``sharur migrate``, unless ``allow_destructive`` is set (that command
    sets it). Returns count of migrations applied.
    """
    current = get_current_version(conn)
    if not allow_destructive:
        blocking = blocking_migrations(conn, target_version)
        if blocking:
            raise DestructiveMigrationError(_database_path(conn), current, blocking)

    applied = 0
    for version, _description, sql, _destructive in _pending(current, target_version):
        if version == 5:
            _add_v5_replicon_columns(conn)
            conn.execute("BEGIN TRANSACTION")
            try:
                _quarantine_pre_v5_system_calls(conn)
                conn.execute(sql)
                conn.commit()
            except Exception:
                # DuckDB can close an aborted transaction automatically.
                # Preserve the migration error rather than masking it.
                with suppress(Exception):
                    conn.rollback()
                raise
        elif version == 6:
            _add_v6_contig_index(conn)
            conn.execute(sql)
        elif version == 8:
            _add_v8_contig_context_columns(conn)
            conn.execute(sql)
        else:
            conn.execute(sql)
        applied += 1

    return applied


__all__ = [
    "MIGRATIONS",
    "Destructive",
    "DestructiveMigrationError",
    "Migration",
    "blocking_migrations",
    "get_current_version",
    "run_migrations",
]
