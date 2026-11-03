"""Dataset-level settings recorded in the database itself.

``dataset_meta`` is a key-value table with the shape of the sidecar stores'
``abundance_meta`` and ``expression_meta``. It holds facts about how a
dataset's inputs are named or produced that analysis code must not guess
globally, for example which contig-name topology convention its assembler
used (:mod:`sharur.contig_context`). A database without the table has no
settings recorded; the table is created on the first write.
"""

from __future__ import annotations

from typing import TYPE_CHECKING


if TYPE_CHECKING:
    import duckdb


CREATE_SQL = "CREATE TABLE IF NOT EXISTS dataset_meta (key VARCHAR PRIMARY KEY, value VARCHAR NOT NULL)"


def _has_table(conn: duckdb.DuckDBPyConnection) -> bool:
    return bool(conn.execute(
        "SELECT COUNT(*) FROM duckdb_tables() "
        "WHERE database_name = current_database() AND schema_name = 'main' AND table_name = 'dataset_meta'"
    ).fetchone()[0])


def get_dataset_meta(conn: duckdb.DuckDBPyConnection, key: str) -> str | None:
    """The recorded value for ``key``, or None. Works on read-only connections."""
    if not _has_table(conn):
        return None
    row = conn.execute("SELECT value FROM dataset_meta WHERE key = ?", [key]).fetchone()
    return row[0] if row else None


def set_dataset_meta(conn: duckdb.DuckDBPyConnection, key: str, value: str | None) -> None:
    """Record ``value`` for ``key``; None removes the key."""
    if value is None:
        if _has_table(conn):
            conn.execute("DELETE FROM dataset_meta WHERE key = ?", [key])
        return
    conn.execute(CREATE_SQL)
    conn.execute("INSERT OR REPLACE INTO dataset_meta (key, value) VALUES (?, ?)", [key, value])


__all__ = ["CREATE_SQL", "get_dataset_meta", "set_dataset_meta"]
