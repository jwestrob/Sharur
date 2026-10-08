"""Validated physical copies of DuckDB files, with immutable source access.

Compaction preserves logical rows and catalog objects. The resulting physical
file needs a new content seal after the owner chooses its final dataset location.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import tempfile
import time
from pathlib import Path
from typing import Any

import duckdb


class CompactionError(RuntimeError):
    """A copy failed its resource, preservation or publication contract."""


def _identifier(value: str) -> str:
    return '"' + value.replace('"', '""') + '"'


def _literal(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def _table(catalog: str, schema: str, table: str) -> str:
    return ".".join(_identifier(value) for value in (catalog, schema, table))


def _state(path: Path) -> tuple[int, int, int, int]:
    stat = path.stat()
    return stat.st_size, stat.st_mtime_ns, stat.st_ino, stat.st_dev


def _file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _catalog(connection: Any, database: str) -> dict[str, Any]:
    tables = {}
    for schema, table, sql, comment in connection.execute(
        "SELECT schema_name,table_name,sql,comment FROM duckdb_tables() "
        "WHERE database_name=? AND NOT internal AND NOT temporary ORDER BY 1,2",
        [database],
    ).fetchall():
        columns = connection.execute(
            "SELECT column_name,data_type,is_nullable,column_default,comment "
            "FROM duckdb_columns() WHERE database_name=? AND schema_name=? "
            "AND table_name=? ORDER BY column_index",
            [database, schema, table],
        ).fetchall()
        constraints = connection.execute(
            "SELECT constraint_type,constraint_column_names,referenced_table,"
            "referenced_column_names,constraint_text FROM duckdb_constraints() "
            "WHERE database_name=? AND schema_name=? AND table_name=?",
            [database, schema, table],
        ).fetchall()
        tables[(schema, table)] = {
            "sql": sql,
            "columns": columns,
            "comment": comment,
            "constraints": sorted(constraints, key=lambda row: json.dumps(row)),
        }
    queries = {
        "schemas": "SELECT schema_name,comment FROM duckdb_schemas() WHERE database_name=? AND NOT internal ORDER BY 1",
        "indexes": "SELECT schema_name,index_name,table_name,is_unique,is_primary,expressions,sql,comment FROM duckdb_indexes() WHERE database_name=? ORDER BY 1,2",
        "views": "SELECT schema_name,view_name,sql,comment FROM duckdb_views() WHERE database_name=? AND NOT internal ORDER BY 1,2",
    }
    return {
        "tables": tables,
        **{name: connection.execute(sql, [database]).fetchall() for name, sql in queries.items()},
    }


def _preflight(connection: Any, database: str, catalog: dict[str, Any]) -> list[tuple[str, str]]:
    for function in ("duckdb_sequences", "duckdb_types", "duckdb_functions"):
        restriction = "" if function == "duckdb_sequences" else " AND NOT internal"
        count = connection.execute(
            f"SELECT COUNT(*) FROM {function}() WHERE database_name=?" + restriction,
            [database],
        ).fetchone()[0]
        if count:
            raise CompactionError(
                "Custom sequences, types and functions require a preservation adapter."
            )
    # Tags lack a supported reconstruction adapter. Refuse to silently discard them.
    for function in (
        "duckdb_schemas",
        "duckdb_tables",
        "duckdb_columns",
        "duckdb_indexes",
        "duckdb_views",
    ):
        fields = {
            column[0]
            for column in connection.execute(f"SELECT * FROM {function}() LIMIT 0").description
        }
        if (
            "tags" in fields
            and connection.execute(
                f"SELECT COUNT(*) FROM {function}() WHERE database_name=? AND cardinality(tags)>0",
                [database],
            ).fetchone()[0]
        ):
            raise CompactionError("Catalog tags require a preservation adapter.")
    tables = catalog["tables"]
    dependencies = {key: set() for key in tables}
    for key, entry in tables.items():
        for kind, _columns, parent, _parent_columns, _text in entry["constraints"]:
            if kind == "FOREIGN KEY":
                parent_key = (key[0], parent)
                if parent_key not in tables:
                    raise CompactionError(
                        "Foreign-key dependency lies outside the copied table set."
                    )
                if parent_key == key:
                    raise CompactionError(
                        "Self-referencing foreign keys require a preservation adapter."
                    )
                dependencies[key].add(parent_key)
    ordered = []
    remaining = set(tables)
    while remaining:
        ready = sorted(key for key in remaining if dependencies[key].isdisjoint(remaining))
        if not ready:
            raise CompactionError("Cyclic foreign keys require a preservation adapter.")
        ordered.extend(ready)
        remaining.difference_update(ready)
    return ordered


def _row_digest(
    connection: Any, database: str, schema: str, table: str, columns: list
) -> dict[str, Any]:
    # Preserve exact JSON storage text, including whitespace/number formatting.
    # Stream sorted SHA-256/multiplicity pairs; sorting can spill within the
    # configured budget, unlike a large string_agg state outside the buffer pool.
    expressions = ",".join(
        _identifier(name)
        + ":="
        + (f"CAST({_identifier(name)} AS VARCHAR)" if "JSON" in dtype else _identifier(name))
        for name, dtype, *_rest in columns
    )
    rows = connection.execute(
        f"SELECT sha256(to_json(struct_pack({expressions}))) AS row_hash, COUNT(*) AS n "
        f"FROM {_table(database, schema, table)} GROUP BY 1 ORDER BY 1"
    )
    digest = hashlib.sha256()
    count = groups = 0
    while batch := rows.fetchmany(4096):
        for row_hash, multiplicity in batch:
            digest.update(f"{row_hash}:{multiplicity}\n".encode("ascii"))
            count += multiplicity
            groups += 1
    return {"rows": count, "distinct_row_hashes": groups, "sha256": digest.hexdigest()}


def _same_catalog(left: dict[str, Any], right: dict[str, Any]) -> bool:
    # DuckDB can reorder equivalent constraint declarations in serialized DDL.
    # Preserve complete constraint multisets plus every column property instead.
    def comparable(catalog):
        return {
            **catalog,
            "tables": {
                key: {**entry, "sql": _logical_ddl(entry)}
                for key, entry in catalog["tables"].items()
            },
        }

    return comparable(left) == comparable(right)


def _ddl_spans(sql: str) -> list[tuple[int, int]]:
    """Locate serialized top-level table clauses, preserving quoted literals."""
    depth = 0
    quote = None
    start = None
    clauses = []
    index = 0
    while index < len(sql):
        char = sql[index]
        if quote:
            if char == quote:
                if index + 1 < len(sql) and sql[index + 1] == quote:
                    index += 2
                    continue
                quote = None
        elif char in ("'", '"'):
            quote = char
        elif char == "(":
            depth += 1
            if depth == 1:
                start = index + 1
        elif char == ")":
            if depth == 1 and start is not None:
                clauses.append((start, index))
            depth -= 1
        elif char == "," and depth == 1:
            clauses.append((start, index))
            start = index + 1
        index += 1
    return clauses


def _logical_ddl(entry: dict[str, Any]) -> dict[str, Any]:
    # Complete DDL retains properties absent from duckdb_columns (collations,
    # generated expressions, etc.); only constraint declaration order is relaxed.
    declared = {row[4] for row in entry["constraints"]}
    clauses = [entry["sql"][first:last].strip() for first, last in _ddl_spans(entry["sql"])]
    return {
        "columns": [clause for clause in clauses if clause not in declared],
        "constraints": sorted(clause for clause in clauses if clause in declared),
    }


def _table_ddl(schema: str, entry: dict[str, Any]) -> str:
    """Requote only serialized top-level FK clauses (DuckDB 1.2 omits quotes)."""
    sql = entry["sql"]
    replacements = []
    foreign_keys = [row for row in entry["constraints"] if row[0] == "FOREIGN KEY"]
    for first, last in _ddl_spans(sql):
        clause = sql[first:last].strip()
        if not clause.startswith("FOREIGN KEY"):
            continue
        matches = [row for row in foreign_keys if row[4] == clause]
        if len(matches) != 1:
            raise CompactionError("Serialized foreign-key clause requires a preservation adapter.")
        _kind, columns, parent, parent_columns, _text = matches[0]
        key = ", ".join(_identifier(name) for name in columns)
        referenced = ", ".join(_identifier(name) for name in parent_columns)
        reference = clause.split("REFERENCES ", 1)[1].rsplit("(", 1)[0].strip()
        if reference == parent:
            parent_name = _identifier(parent)
        elif reference == schema + "." + parent:
            parent_name = ".".join(_identifier(name) for name in (schema, parent))
        else:
            raise CompactionError(
                "Qualified foreign-key reference requires a preservation adapter."
            )
        replacements.append(
            (first, last, f"FOREIGN KEY ({key}) REFERENCES {parent_name} ({referenced})")
        )
    if len(replacements) != len(foreign_keys):
        raise CompactionError("Serialized foreign keys require a preservation adapter.")
    for first, last, replacement in reversed(replacements):
        sql = sql[:first] + replacement + sql[last:]
    return sql


def _copy(connection: Any, catalog: dict[str, Any], ordered: list[tuple[str, str]]) -> None:
    connection.execute("USE compacted")
    for schema, _comment in catalog["schemas"]:
        connection.execute("CREATE SCHEMA " + _identifier(schema))
    for schema, table in ordered:
        connection.execute(_table_ddl(schema, catalog["tables"][(schema, table)]))
        connection.execute(
            f"INSERT INTO {_table('compacted', schema, table)} "
            f"SELECT * FROM {_table('original', schema, table)}"
        )
    for _schema, _name, _table_name, _unique, _primary, _expressions, sql, _comment in catalog[
        "indexes"
    ]:
        connection.execute(sql)
    pending = list(catalog["views"])
    while pending:
        retry = []
        for view in pending:
            try:
                connection.execute(view[2])
            except (duckdb.BinderException, duckdb.CatalogException):
                retry.append(view)
        if len(retry) == len(pending):
            raise CompactionError("View dependencies require a preservation adapter.")
        pending = retry
    for schema, comment in catalog["schemas"]:
        if comment is not None:
            connection.execute(f"COMMENT ON SCHEMA {_identifier(schema)} IS {_literal(comment)}")
    for (schema, table), entry in catalog["tables"].items():
        if entry["comment"] is not None:
            connection.execute(
                f"COMMENT ON TABLE {_table('compacted', schema, table)} IS {_literal(entry['comment'])}"
            )
        for column, _dtype, _nullable, _default, comment in entry["columns"]:
            if comment is not None:
                name = ".".join(_identifier(value) for value in (schema, table, column))
                connection.execute(f"COMMENT ON COLUMN {name} IS {_literal(comment)}")
    for schema, name, *_fields in catalog["indexes"]:
        if _fields[-1] is not None:
            connection.execute(
                f"COMMENT ON INDEX {_identifier(schema)}.{_identifier(name)} IS {_literal(_fields[-1])}"
            )
    for schema, name, _sql, comment in catalog["views"]:
        if comment is not None:
            connection.execute(
                f"COMMENT ON VIEW {_identifier(schema)}.{_identifier(name)} IS {_literal(comment)}"
            )


def compact_database(
    source: str | Path,
    output: str | Path,
    *,
    threads: int = 2,
    memory_limit: str = "8GB",
    storage_version: str | None = None,
    receipt: str | Path | None = None,
) -> dict[str, Any]:
    """Publish one validated fresh copy, preserving the source and sidecars.

    Existing outputs and receipts are refused. Only candidate-local temporary
    files are written. A failure retains the building file for local diagnosis;
    exception messages and receipts contain no copied row payloads.
    """
    source_path = Path(source).expanduser().resolve()
    requested_output = Path(output).expanduser()
    if requested_output.is_symlink():
        raise CompactionError("Output must be a new regular file.")
    output_path = requested_output.resolve()
    requested_receipt = Path(receipt).expanduser() if receipt is not None else None
    if requested_receipt is not None and requested_receipt.is_symlink():
        raise CompactionError("Receipt must be a new regular file.")
    receipt_path = requested_receipt.resolve() if requested_receipt is not None else None
    if not source_path.is_file():
        raise CompactionError("Source database is missing.")
    if (
        output_path in (source_path, Path(str(source_path) + ".wal"))
        or output_path.exists()
        or Path(str(output_path) + ".wal").exists()
    ):
        raise CompactionError("Output must be a fresh path distinct from the source.")
    if receipt_path is not None and (
        receipt_path in (source_path, output_path)
        or receipt_path.exists()
        or receipt_path in (Path(str(source_path) + ".wal"), Path(str(output_path) + ".wal"))
    ):
        raise CompactionError("Receipt must be a fresh path distinct from database files.")
    if Path(str(source_path) + ".wal").exists():
        raise CompactionError("Source must be checkpointed with its writer closed.")
    if (
        isinstance(threads, bool)
        or not isinstance(threads, int)
        or not 1 <= threads <= (os.cpu_count() or 1)
    ):
        raise CompactionError("Threads must be a positive count within available CPUs.")
    if (
        not re.fullmatch(
            r"[0-9]+(?:\.[0-9]+)?\s*(?:KB|MB|GB|TB|KiB|MiB|GiB|TiB)", memory_limit, re.IGNORECASE
        )
        or float(re.match(r"[0-9.]+", memory_limit)[0]) <= 0
    ):
        raise CompactionError("Memory limit must specify a positive bounded size.")
    if storage_version is not None and not re.fullmatch(r"v\d+\.\d+\.\d+", storage_version):
        raise CompactionError("Storage version must be an explicit version such as v1.0.0.")
    initial_state = _state(source_path)
    started = time.monotonic()
    report: dict[str, Any] = {
        "status": "BUILDING",
        "source_filename": source_path.name,
        "output_filename": output_path.name,
        "duckdb_version": duckdb.__version__,
        "threads": threads,
        "memory_limit": memory_limit,
        "storage_version": storage_version or "runtime_default",
        "source_read_only": True,
        "tables": [],
        "reseal_required": True,
    }
    building = None
    work_directory = None
    receipt_handle = None
    connection = None
    phase = "preflight"
    try:
        connection = duckdb.connect(
            ":memory:", config={"threads": threads, "memory_limit": memory_limit}
        )
        connection.execute(f"ATTACH {_literal(str(source_path))} AS original (READ_ONLY)")
        catalog = _catalog(connection, "original")
        ordered = _preflight(connection, "original", catalog)
        source_hash = _file_hash(source_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        work_directory = Path(
            tempfile.mkdtemp(
                prefix=f".{output_path.name}.", suffix=".building", dir=output_path.parent
            )
        )
        building = work_directory / "database.duckdb"
        spill_directory = work_directory / "spill"
        report["building_filename"] = str(building.relative_to(output_path.parent))
        if receipt_path is not None:
            receipt_path.parent.mkdir(parents=True, exist_ok=True)
            receipt_handle = receipt_path.open("x")
        connection.execute("SET temp_directory=" + _literal(str(spill_directory)))
        option = f" (STORAGE_VERSION {_literal(storage_version)})" if storage_version else ""
        connection.execute(f"ATTACH {_literal(str(building))} AS compacted" + option)
        phase = "copy"
        _copy(connection, catalog, ordered)
        connection.execute("CHECKPOINT compacted")
        connection.close()
        connection = duckdb.connect(
            ":memory:",
            config={
                "threads": threads,
                "memory_limit": memory_limit,
                "temp_directory": str(spill_directory),
            },
        )
        connection.execute(f"ATTACH {_literal(str(source_path))} AS original (READ_ONLY)")
        connection.execute(f"ATTACH {_literal(str(building))} AS compacted (READ_ONLY)")
        phase = "validation"
        if not _same_catalog(catalog, _catalog(connection, "compacted")):
            raise CompactionError("Copied catalog definitions differ from the source.")
        for schema, table in ordered:
            columns = catalog["tables"][(schema, table)]["columns"]
            before = _row_digest(connection, "original", schema, table, columns)
            after = _row_digest(connection, "compacted", schema, table, columns)
            if before != after:
                raise CompactionError("Copied full-row digest differs from the source.")
            report["tables"].append(
                {"schema": schema, "table": table, **after, "full_row_parity": True}
            )
        # A copied view must remain usable with the source catalog detached.
        connection.execute("DETACH original")
        for schema, name, _sql, _comment in catalog["views"]:
            connection.execute(f"SELECT * FROM {_table('compacted', schema, name)} LIMIT 0")
        # Retain a source read lock through publication while releasing every
        # candidate handle. Reattachment detects writers in the detach interval.
        connection.execute(f"ATTACH {_literal(str(source_path))} AS original (READ_ONLY)")
        connection.execute("DETACH compacted")
        if (
            Path(str(source_path) + ".wal").exists()
            or _state(source_path) != initial_state
            or _file_hash(source_path) != source_hash
        ):
            raise CompactionError("Source changed while compacting; candidate remains unpublished.")
        if Path(str(building) + ".wal").exists():
            raise CompactionError("Candidate retains a WAL after checkpoint and close.")
        report.update(
            source_sha256=source_hash,
            output_sha256=_file_hash(building),
            source_bytes=initial_state[0],
            output_bytes=building.stat().st_size,
            source_file_preserved=True,
            catalog_parity=True,
            table_count=len(ordered),
            explicit_index_count=len(catalog["indexes"]),
            view_count=len(catalog["views"]),
            elapsed_seconds=round(time.monotonic() - started, 3),
        )
        phase = "publication"
        if Path(str(output_path) + ".wal").exists():
            raise CompactionError(
                "Destination WAL appeared during validation; candidate remains unpublished."
            )
        # Atomic exclusive publication on the same filesystem. rename() alone
        # overwrites an output created by another process during validation.
        os.link(building, output_path)
        report["status"] = "PASS"
        try:
            connection.close()
            connection = None
            building.unlink()
            if spill_directory.exists():
                spill_directory.rmdir()
            work_directory.rmdir()
        except OSError as exc:
            report["cleanup_warning"] = {
                "error_type": type(exc).__name__,
                "task_owned_temporary_artifacts_retained": True,
            }
        else:
            report.pop("building_filename", None)
        if receipt_handle is not None:
            json.dump(report, receipt_handle, indent=2)
            receipt_handle.write("\n")
        return report
    except BaseException as exc:
        report.update(status="FAILED", failure_phase=phase, error_type=type(exc).__name__)
        if receipt_handle is not None:
            receipt_handle.seek(0)
            receipt_handle.truncate()
            json.dump(report, receipt_handle, indent=2)
            receipt_handle.write("\n")
        if isinstance(exc, (CompactionError, KeyboardInterrupt, SystemExit)):
            raise
        raise CompactionError(
            f"Compaction failed during {phase} ({type(exc).__name__}); source remains unchanged."
        ) from None
    finally:
        if connection is not None:
            connection.close()
        if receipt_handle is not None:
            receipt_handle.close()
