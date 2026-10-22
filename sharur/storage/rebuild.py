"""Rebuild a DuckDB database into a new file, rewriting rows table by table.

DuckDB checks foreign keys eagerly. A row that another table references cannot
be deleted, and an indexed column of such a row cannot be updated, even inside a
transaction that also removes the referencing rows (DuckDB 1.2). Restructuring a
canonical database, such as dropping genomes or moving proteins between
genomes, therefore copies every table into a fresh file with the original DDL,
filtering or rewriting rows on the way. The caller swaps the new file in once it
validates, and the source stays untouched as the backup.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping
from pathlib import Path

_REFERENCES = re.compile(r'REFERENCES\s+(?:"?\w+"?\.)?"?(\w+)"?\s*\(', re.IGNORECASE)


def table_order(ddl: Mapping[str, str]) -> list[str]:
    """Tables ordered so every table follows the tables its foreign keys reference."""
    deps = {t: {r for r in _REFERENCES.findall(sql) if r != t and r in ddl} for t, sql in ddl.items()}
    order: list[str] = []
    done: set[str] = set()
    while len(order) < len(ddl):
        ready = sorted(t for t in ddl if t not in done and deps[t] <= done)
        if not ready:
            raise ValueError(f"Cyclic foreign keys among {sorted(set(ddl) - done)}")
        order.extend(ready)
        done.update(ready)
    return order


def rebuild_database(
    src: Path,
    dst: Path,
    selects: Mapping[str, str] | None = None,
    *,
    setup_sql: Iterable[str] = (),
    threads: int | None = None,
) -> dict[str, int]:
    """Copy ``src`` into a new ``dst`` with the same tables, constraints, indexes, views and sequences.

    ``selects`` maps a table name to a SELECT that produces its rows, matched to
    the table's columns by name; a table without an entry is copied whole. The
    source is attached read-only as ``src``, so a select reads ``src.<table>``.
    ``setup_sql`` runs first in the destination connection, e.g. to create
    TEMP staging tables that the selects join against; temporary objects are
    not persisted. Indexes are created after the rows load. Returns the row
    count of every table. ``dst`` must not exist, and it is removed again if the
    rebuild fails.
    """
    import duckdb  # noqa: PLC0415

    src, dst = Path(src), Path(dst)
    selects = dict(selects or {})
    if dst.exists():
        raise FileExistsError(dst)
    con = duckdb.connect(str(dst), config={"threads": threads} if threads else {})
    try:
        con.execute(f"ATTACH '{src}' AS src (READ_ONLY)")
        ddl = dict(con.execute(
            "SELECT table_name, sql FROM duckdb_tables() "
            "WHERE database_name = 'src' AND schema_name = 'main' AND NOT temporary").fetchall())
        unknown = sorted(set(selects) - set(ddl))
        if unknown:
            raise KeyError(f"selects name tables absent from the source: {unknown}")
        for statement in setup_sql:
            con.execute(statement)
        for name, start, increment, minimum, maximum, cycle, last in con.execute(
                "SELECT sequence_name, start_value, increment_by, min_value, max_value, cycle, last_value "
                "FROM duckdb_sequences() WHERE database_name = 'src' AND NOT temporary").fetchall():
            nxt = start if last is None else last + increment
            con.execute(f'CREATE SEQUENCE "{name}" START WITH {nxt} INCREMENT BY {increment} '
                        f"MINVALUE {minimum} MAXVALUE {maximum} {'CYCLE' if cycle else 'NO CYCLE'}")
        order = table_order(ddl)
        for table in order:
            con.execute(ddl[table])
        counts: dict[str, int] = {}
        for table in order:
            query = selects.get(table, f'SELECT * FROM src.main."{table}"')
            con.execute(f'INSERT INTO main."{table}" BY NAME {query}')
            counts[table] = con.execute(f'SELECT COUNT(*) FROM main."{table}"').fetchone()[0]
        for (sql,) in con.execute(
                "SELECT sql FROM duckdb_indexes() WHERE database_name = 'src' AND sql IS NOT NULL").fetchall():
            con.execute(sql)
        for (sql,) in con.execute(
                "SELECT sql FROM duckdb_views() WHERE database_name = 'src' AND NOT internal "
                "AND schema_name = 'main'").fetchall():
            con.execute(sql)
        con.execute("DETACH src")
        con.execute("CHECKPOINT")
    except BaseException:
        con.close()
        dst.unlink(missing_ok=True)
        Path(f"{dst}.wal").unlink(missing_ok=True)
        raise
    con.close()
    return counts
