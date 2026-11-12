"""Sharur's SQL surface: per-connection macros that work on read-only connections.

The prelude (``prelude.sql`` beside this module) defines TEMP macros only: no views, no
tables. A read-only connection can create them, they live as long as the connection, and
they never appear in ``information_schema`` or ``SHOW TABLES``, so table enumerations,
capability briefs and dataset seals read the same database with or without them.

:class:`~sharur.storage.duckdb_store.DuckDBStore` installs the prelude on every connection
it opens (``sql_prelude=False`` or ``SHARUR_SQL_PRELUDE=0`` turns that off), so
``Sharur(path, read_only=True).store.execute("FROM neighbors('PF00069', k := 3)")`` works.
``sharur reverify`` and the review verifier install it on their replay connections. A
generic DuckDB client runs the file itself (or :func:`prelude_text`, ``sharur sql prelude``)
after ``ATTACH ... (READ_ONLY); USE alias;``.

File format. Header lines::

    -- @sharur-sql-version: 1
    -- @requires: proteins, annotations, bins

then one statement per ``-- @object NAME KIND: DOC`` marker. The installer runs each
statement on its own and records what installed; when a ``@requires`` table is missing from
the default database, nothing is installed and the report says why.

Sequences stay out of query results by construction: ``genes()`` is the protein table
without its sequence column, no macro returns a sequence, and :func:`write_fasta` (also
registered as the SQL function ``write_fasta(path, ids)``) writes FASTA to disk and returns
only the path and the counts.
"""

from __future__ import annotations

import functools
import logging
import os
import re
from collections.abc import Iterable, Sequence
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any


logger = logging.getLogger(__name__)

PRELUDE_PATH = Path(__file__).with_name("prelude.sql")
ENV_DISABLE = "SHARUR_SQL_PRELUDE"
REGISTRY_MACRO = "sharur_sql"

_OBJECT = re.compile(r"^--\s*@object\s+([A-Za-z_][\w]*)\s+([\w]+)\s*:\s*(.*)$")
_HEADER = re.compile(r"^--\s*@(sharur-sql-version|requires)\s*:\s*(.*)$")
_SIGNATURE = re.compile(r"\bMACRO\s+(?:[\w.]+\.)?([A-Za-z_]\w*)\s*\(([^)]*)\)", re.I)


# --------------------------------------------------------------------------- #
# Parsing
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class SqlObject:
    """One installable statement of the prelude."""

    name: str
    kind: str
    doc: str
    signature: str
    statement: str


@dataclass(frozen=True)
class SqlFile:
    path: str
    version: str | None
    requires: tuple[str, ...]
    objects: tuple[SqlObject, ...]


def _signature(name: str, statement: str) -> str:
    m = _SIGNATURE.search(statement)
    if not m or m.group(1) != name:
        return name
    return f"{name}({' '.join(m.group(2).split())})"


@functools.lru_cache(maxsize=8)
def _parse_cached(path: str, mtime_ns: int) -> SqlFile:
    version: str | None = None
    requires: tuple[str, ...] = ()
    objects: list[SqlObject] = []
    current: tuple[str, str, str] | None = None
    body: list[str] = []

    def flush() -> None:
        statement = "\n".join(body).strip()
        if current is not None and statement:
            name, kind, doc = current
            objects.append(SqlObject(name=name, kind=kind, doc=doc.strip(),
                                     signature=_signature(name, statement), statement=statement))

    for line in Path(path).read_text().splitlines():
        m = _OBJECT.match(line.strip())
        if m:
            flush()
            current, body = (m.group(1), m.group(2), m.group(3)), []
            continue
        if current is None:
            h = _HEADER.match(line.strip())
            if h and h.group(1) == "sharur-sql-version":
                version = h.group(2).strip()
            elif h:
                requires = tuple(t.strip() for t in h.group(2).split(",") if t.strip())
            continue
        body.append(line)
    flush()
    return SqlFile(path=path, version=version, requires=requires, objects=tuple(objects))


def parse_prelude(path: str | Path = PRELUDE_PATH) -> SqlFile:
    """Objects, version and table requirements of the prelude file."""
    p = Path(path).expanduser().resolve()
    return _parse_cached(str(p), p.stat().st_mtime_ns)


# --------------------------------------------------------------------------- #
# Installation
# --------------------------------------------------------------------------- #
@dataclass
class PreludeReport:
    """What one installation did."""

    version: str | None = None
    installed: list[str] = field(default_factory=list)
    failed: list[dict[str, str]] = field(default_factory=list)
    skipped: str | None = None
    udfs: list[str] = field(default_factory=list)
    # Macros left out because their optional layer is absent: {"name", "layer"}.
    unavailable: list[dict[str, str]] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.failed and self.skipped is None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


# Tables of the optional layers. DuckDB 1.5 binds the tables of a table macro when the macro is
# created, so a macro over an absent layer cannot be installed; it is reported as unavailable,
# as are macros built on such a macro, and the rest of the surface installs.
LAYER_TABLES = {
    "protein_families": "protein families",
    "family_stats": "protein families",
    "protein_family_runs": "protein families",
    "family_annotations": "remote homology",
}


def _absent_layer(error: Exception, unavailable: dict[str, str]) -> str | None:
    """The optional layer whose absence explains ``error``, or None for a genuine failure."""
    message = str(error)
    match = re.search(r"(Table|Table Function) with name (\w+) does not exist", message)
    if not match:
        return None
    name = match.group(2)
    return LAYER_TABLES.get(name) or unavailable.get(name)


def object_layers(objects: Sequence[SqlObject]) -> dict[str, set[str]]:
    """Object name -> the optional layers it reads, directly or through other objects."""
    layers: dict[str, set[str]] = {}
    for obj in objects:  # objects are defined before they are used
        needed = {layer for table, layer in LAYER_TABLES.items() if re.search(rf"\b{table}\b", obj.statement)}
        for name, inherited in layers.items():
            if re.search(rf"\b{name}\s*\(", obj.statement):
                needed |= inherited
        if needed:
            layers[obj.name] = needed
    return layers


def available_tables(con: Any) -> set[str]:
    """Tables and views of the connection's default database."""
    rows = con.execute(
        "SELECT table_name FROM information_schema.tables WHERE table_catalog = current_database()"
    ).fetchall()
    return {str(r[0]) for r in rows}


def registry_statement(objects: Sequence[SqlObject], udfs: Sequence[str] = (),
                       version: str | None = None) -> str:
    """``sharur_sql()``: one row per object of the surface (name, kind, signature, doc, version)."""
    rows = [(o.name, o.kind, o.signature, o.doc) for o in objects]
    if "write_fasta" in udfs:
        rows.append(("write_fasta", "python_function", "write_fasta(path, ids)",
                     "write the sequences of ids (a list) to a FASTA file; returns {path, n_sequences, "
                     "n_requested, n_missing}, never a sequence"))
    rows.append((REGISTRY_MACRO, "table_macro", f"{REGISTRY_MACRO}()",
                 "the objects of the Sharur SQL surface on this connection"))

    def lit(value: Any) -> str:
        return "NULL" if value is None else "'" + str(value).replace("'", "''") + "'"

    values = ",\n        ".join("(" + ", ".join(lit(v) for v in (*row, version)) + ")" for row in rows)
    return (f"CREATE OR REPLACE TEMP MACRO {REGISTRY_MACRO}() AS TABLE\n"
            f"    SELECT * FROM (VALUES\n        {values}\n    ) t(name, kind, signature, doc, version)")


def install_prelude(con: Any, *, udfs: bool = True, strict: bool = False) -> PreludeReport:
    """Install the prelude on one DuckDB connection (read-only connections included).

    ``strict`` raises on the first failing statement; otherwise failures are recorded in
    the report and the remaining objects still install.
    """
    sql_file = parse_prelude()
    report = PreludeReport(version=sql_file.version)
    missing = [t for t in sql_file.requires if t not in available_tables(con)]
    if missing:
        report.skipped = f"missing tables: {', '.join(missing)}"
        if strict:
            raise RuntimeError(f"Sharur SQL prelude needs {report.skipped}")
        return report
    installed: list[SqlObject] = []
    unavailable: dict[str, str] = {}
    for obj in sql_file.objects:
        try:
            con.execute(obj.statement)
        except Exception as exc:
            layer = _absent_layer(exc, unavailable)
            if layer is not None:
                unavailable[obj.name] = layer
                report.unavailable.append({"name": obj.name, "layer": layer})
                continue
            if strict:
                raise
            report.failed.append({"name": obj.name, "error": f"{type(exc).__name__}: {exc}"[:500]})
            continue
        installed.append(obj)
        report.installed.append(obj.name)
    if udfs:
        try:
            register_udfs(con)
            report.udfs.append("write_fasta")
        except Exception as exc:
            if strict:
                raise
            report.failed.append({"name": "write_fasta", "error": f"{type(exc).__name__}: {exc}"[:500]})
    con.execute(registry_statement(installed, report.udfs, sql_file.version))
    report.installed.append(REGISTRY_MACRO)
    return report


def prelude_enabled(flag: bool | None = None) -> bool:
    """``flag`` when given; otherwise on unless ``SHARUR_SQL_PRELUDE`` is 0/false/off/no."""
    if flag is not None:
        return bool(flag)
    return os.environ.get(ENV_DISABLE, "1").strip().lower() not in {"0", "false", "off", "no"}


def install_prelude_quietly(con: Any) -> PreludeReport | None:
    """Install for a general-purpose connection; never raises (problems are logged)."""
    try:
        report = install_prelude(con)
    except Exception as exc:
        logger.warning("Sharur SQL prelude not installed: %s: %s", type(exc).__name__, exc)
        return None
    if report.skipped:
        logger.debug("Sharur SQL prelude skipped: %s", report.skipped)
    for item in report.failed:
        logger.debug("Sharur SQL prelude object %s failed: %s", item["name"], item["error"])
    return report


def prelude_text(tables: set[str] | None = None) -> str:
    """The prelude as one SQL script for a generic DuckDB client, with the registry macro.

    ``tables``: the tables of the database the script is for. Macros of an optional layer whose
    tables are absent are then left out, because DuckDB 1.3 and later bind a table macro's tables
    when the macro is created and stop a one-shot script at the first missing one.
    """
    sql_file = parse_prelude()
    layers = object_layers(sql_file.objects)
    absent = {layer for table, layer in LAYER_TABLES.items() if tables is not None and table not in tables}
    objects = [o for o in sql_file.objects if not layers.get(o.name, set()) & absent]
    parts = [f"-- Sharur SQL surface v{sql_file.version}, generated by `sharur sql prelude`. Run after:",
             "--   ATTACH 'PATH/sharur.duckdb' AS ds (READ_ONLY); USE ds;",
             f"-- Needs tables: {', '.join(sql_file.requires)}"]
    if absent:
        parts.append(f"-- Left out (layer not in the database): {', '.join(sorted(absent))}")
    parts.append("")
    for obj in objects:
        parts.append(f"-- @object {obj.name} {obj.kind}: {obj.doc}")
        parts.append(obj.statement.rstrip().rstrip(";") + ";")
        parts.append("")
    parts.append(registry_statement(objects, version=sql_file.version) + ";")
    return "\n".join(parts) + "\n"


# --------------------------------------------------------------------------- #
# FASTA out of the database, never into a result
# --------------------------------------------------------------------------- #
def _id_list(ids: str | Iterable[str]) -> list[str]:
    items = ids.split(",") if isinstance(ids, str) else [str(i) for i in ids if i is not None]
    seen: dict[str, None] = {}
    for item in items:
        item = item.strip()
        if item:
            seen.setdefault(item, None)
    return list(seen)


def write_fasta(
    con: Any,
    ids: str | Iterable[str],
    path: str | Path,
    *,
    strip_stop: bool = True,
) -> dict[str, Any]:
    """Write the sequences of ``ids`` to ``path`` (atomic replace); return path and counts only.

    ``con`` is a DuckDB connection or a DuckDBStore. Records follow the order of ``ids``;
    IDs absent from ``proteins`` are counted in ``n_missing``. A trailing stop ``*`` is
    removed unless ``strip_stop=False``.
    """
    conn = getattr(con, "conn", con)
    wanted = _id_list(ids)
    out = Path(path).expanduser().resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    rows: list[tuple] = []
    if wanted:
        cursor = conn.cursor()
        try:
            rows = cursor.execute(
                "SELECT protein_id, sequence FROM proteins "
                "WHERE protein_id IN (SELECT unnest(?::VARCHAR[]))", [wanted]).fetchall()
        finally:
            cursor.close()
    by_id = {str(pid): seq for pid, seq in rows if seq is not None}
    tmp = out.with_name(f".{out.name}.{os.getpid()}.tmp")
    written = 0
    with tmp.open("w") as handle:
        for pid in wanted:
            seq = by_id.get(pid)
            if seq is None:
                continue
            seq = "".join(str(seq).split())
            if strip_stop:
                seq = seq.rstrip("*")
            handle.write(f">{pid}\n{seq}\n")
            written += 1
    os.replace(tmp, out)
    return {"path": str(out), "n_sequences": written, "n_requested": len(wanted),
            "n_missing": len(wanted) - written}


def register_udfs(con: Any) -> None:
    """Register ``write_fasta(path VARCHAR, ids VARCHAR[])`` on a DuckDB connection."""
    import duckdb  # noqa: PLC0415
    try:  # DuckDB 1.4 added duckdb.sqltypes; 1.5 removed duckdb.typing
        import duckdb.sqltypes as T  # noqa: PLC0415
    except ImportError:  # DuckDB before 1.4
        import duckdb.typing as T  # noqa: PLC0415

    def _udf(path: str, ids: list[str]) -> dict[str, Any]:
        return write_fasta(con, ids or [], path)

    try:
        con.remove_function("write_fasta")
    except Exception:
        pass
    con.create_function(
        "write_fasta", _udf, [T.VARCHAR, duckdb.list_type(T.VARCHAR)],
        duckdb.struct_type({"path": T.VARCHAR, "n_sequences": T.BIGINT,
                            "n_requested": T.BIGINT, "n_missing": T.BIGINT}),
        side_effects=True,
    )


# --------------------------------------------------------------------------- #
# Connections and sequence-free results
# --------------------------------------------------------------------------- #
SEQUENCE_COLUMNS = {"sequence", "protein_sequence", "amino_acid_sequence", "nucleotide_sequence",
                    "aa_sequence", "nt_sequence", "dna_sequence", "seq"}
_SEQUENCE_TEXT = re.compile(r"[ACDEFGHIKLMNPQRSTVWYBXZUO*\-]{30,}")


def connect(db_path: str | Path, *, threads: int | None = 4, memory_limit: str | None = None) -> Any:
    """A read-only DuckDB connection to a dataset database with the prelude installed."""
    import duckdb  # noqa: PLC0415

    con = duckdb.connect(str(Path(db_path).expanduser()), read_only=True)
    if threads:
        con.execute(f"SET threads = {int(threads)}")
    if memory_limit:
        con.execute("SET memory_limit = ?", [memory_limit])
    install_prelude(con)
    return con


def looks_like_sequence(value: Any) -> bool:
    """True for an unbroken run of 30+ upper-case residue letters (in nested values too).

    Descriptions and names carry spaces, digits or lower case, so they pass.
    """
    if isinstance(value, str):
        return bool(_SEQUENCE_TEXT.search(value))
    if isinstance(value, (list, tuple)):
        return any(looks_like_sequence(v) for v in value)
    if isinstance(value, dict):
        return any(looks_like_sequence(v) for v in value.values())
    return False


def sequence_leak(columns: Sequence[str], rows: Iterable[Sequence[Any]]) -> str | None:
    """Why a result would expose raw sequence, or None when it is sequence-free."""
    named = [c for c in columns if str(c).strip().lower() in SEQUENCE_COLUMNS]
    if named:
        return f"result has sequence column(s): {', '.join(named)}"
    for row in rows:
        for column, value in zip(columns, row):
            if looks_like_sequence(value):
                return f"column {column!r} holds sequence-like text"
    return None


__all__ = [
    "ENV_DISABLE",
    "PRELUDE_PATH",
    "PreludeReport",
    "SqlFile",
    "SqlObject",
    "available_tables",
    "connect",
    "install_prelude",
    "install_prelude_quietly",
    "looks_like_sequence",
    "parse_prelude",
    "prelude_enabled",
    "prelude_text",
    "object_layers",
    "register_udfs",
    "registry_statement",
    "sequence_leak",
    "write_fasta",
]
