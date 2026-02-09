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
]


def get_current_version(conn: "duckdb.DuckDBPyConnection") -> int:
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
    conn: "duckdb.DuckDBPyConnection", target_version: int | None = None
) -> int:
    """Apply pending migrations up to target_version (default: latest).

    Returns count of migrations applied.
    """
    current = get_current_version(conn)
    if target_version is None:
        target_version = max(v for v, _, _ in MIGRATIONS) if MIGRATIONS else 0

    applied = 0
    for version, description, sql in sorted(MIGRATIONS, key=lambda m: m[0]):
        if version <= current:
            continue
        if version > target_version:
            break

        conn.execute(sql)
        applied += 1

    return applied


__all__ = ["MIGRATIONS", "get_current_version", "run_migrations"]
