"""Schema versioning and migrations for Sharur databases."""

from __future__ import annotations

from typing import TYPE_CHECKING


if TYPE_CHECKING:
    import duckdb

# Each migration: (version, description, sql)
MIGRATIONS: list[tuple[int, str, str]] = [
    (
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
    (
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
    (
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
    (
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
]


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


def run_migrations(
    conn: duckdb.DuckDBPyConnection, target_version: int | None = None
) -> int:
    """Apply pending migrations up to target_version (default: latest).

    Returns count of migrations applied.
    """
    current = get_current_version(conn)
    if target_version is None:
        target_version = max(v for v, _, _ in MIGRATIONS) if MIGRATIONS else 0

    applied = 0
    for version, _description, sql in sorted(MIGRATIONS, key=lambda m: m[0]):
        if version <= current:
            continue
        if version > target_version:
            break

        conn.execute(sql)
        applied += 1

    return applied


__all__ = ["MIGRATIONS", "get_current_version", "run_migrations"]
