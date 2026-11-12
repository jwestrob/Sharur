"""``sharur sql``: the SQL surface from the command line (prelude script, read-only queries)."""

from __future__ import annotations

import json
import sys
from enum import Enum
from pathlib import Path
from typing import Any

import typer


sql_app = typer.Typer(help="Sharur's SQL surface: per-connection macros (neighbors, slot, feature_presence, ...).",
                      no_args_is_help=True)


class QueryFormat(str, Enum):
    markdown = "markdown"
    tsv = "tsv"
    json = "json"


def _cell(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, (list, tuple, dict)):
        return json.dumps(value, default=str)
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def render(columns: list[str], rows: list[tuple], fmt: QueryFormat, *, total: int | None = None) -> str:
    """Rows as a markdown table, TSV or JSON records."""
    if fmt == QueryFormat.json:
        return json.dumps({"columns": columns, "rows": [dict(zip(columns, r)) for r in rows],
                           "total_rows": total if total is not None else len(rows)}, default=str, indent=1)
    if fmt == QueryFormat.tsv:
        lines = ["\t".join(columns)]
        lines += ["\t".join(_cell(v).replace("\t", " ").replace("\n", " ") for v in r) for r in rows]
        return "\n".join(lines)
    lines = ["| " + " | ".join(columns) + " |", "|" + "---|" * len(columns)]
    lines += ["| " + " | ".join(_cell(v).replace("|", "\\|").replace("\n", " ") for v in r) + " |" for r in rows]
    if total is not None and total > len(rows):
        lines.append(f"\n{len(rows)} of {total} rows shown (--limit).")
    return "\n".join(lines)


@sql_app.command("prelude")
def prelude_command(
    out: Path | None = typer.Option(None, "--out", "-o", help="Write the script here instead of stdout."),
    db: Path | None = typer.Option(
        None, "--db", "-d",
        help="Dataset DuckDB the script is for: leaves out macros of optional layers it lacks "
             "(needed with DuckDB 1.3 and later, which stop the script at the first such macro).",
    ),
):
    """Print the prelude as one SQL script for a generic DuckDB client.

    Run it after `ATTACH 'data/DATASET/sharur.duckdb' AS ds (READ_ONLY); USE ds;`.
    """
    from sharur.sql import available_tables, prelude_text  # noqa: PLC0415

    tables = None
    if db is not None:
        import duckdb  # noqa: PLC0415

        con = duckdb.connect(str(db), read_only=True)
        try:
            tables = available_tables(con)
        finally:
            con.close()
    text = prelude_text(tables)
    if out is not None:
        out.write_text(text)
        typer.echo(f"Wrote {out}", err=True)
    else:
        typer.echo(text)


@sql_app.command("objects")
def objects_command(
    db: Path = typer.Option(..., "--db", "-d", help="Dataset DuckDB."),
    output_format: QueryFormat = typer.Option(QueryFormat.markdown, "--format", "-f"),
):
    """List the SQL objects installed on a connection to this database, and any that failed."""
    import duckdb  # noqa: PLC0415

    from sharur.sql import install_prelude  # noqa: PLC0415

    con = duckdb.connect(str(db), read_only=True)
    try:
        report = install_prelude(con)
        cur = con.execute("FROM sharur_sql()")
        columns = [d[0] for d in cur.description]
        rows = cur.fetchall()
    finally:
        con.close()
    typer.echo(render(columns, rows, output_format))
    if report.skipped:
        typer.echo(f"prelude skipped: {report.skipped}", err=True)
    for item in report.unavailable:
        typer.echo(f"not installed {item['name']}: the {item['layer']} layer is not in this database", err=True)
    for item in report.failed:
        typer.echo(f"failed {item['name']}: {item['error']}", err=True)


@sql_app.command("query")
def query_command(
    sql: str = typer.Argument(..., help="SQL to run; '-' reads it from stdin."),
    db: Path = typer.Option(..., "--db", "-d", help="Dataset DuckDB (opened read-only)."),
    output_format: QueryFormat = typer.Option(QueryFormat.markdown, "--format", "-f"),
    limit: int = typer.Option(200, "--limit", "-n", help="Rows to print (0 = all)."),
    out: Path | None = typer.Option(None, "--out", "-o", help="Also write every row as TSV here."),
    threads: int = typer.Option(4, "--threads", help="DuckDB threads."),
    memory_limit: str = typer.Option("8GB", "--memory-limit", help="DuckDB memory limit."),
):
    """Run one read-only query with the prelude installed; results never carry sequences.

    Example: sharur sql query "FROM slot_tally('pfam:PF00069', 'pfam:PF00072', split := 'K00001')" -d DB
    """
    from sharur.sql import sequence_leak  # noqa: PLC0415
    from sharur.storage.duckdb_store import DuckDBStore  # noqa: PLC0415

    text = sys.stdin.read() if sql == "-" else sql
    store = DuckDBStore(db, read_only=True, threads=threads, memory_limit=memory_limit, sql_prelude=True)
    try:
        try:
            cur = store.conn.execute(text)
        except Exception as exc:
            typer.echo(f"query failed: {type(exc).__name__}: {exc}", err=True)
            raise typer.Exit(1) from exc
        if cur.description is None:
            typer.echo("(no result)")
            return
        columns = [d[0] for d in cur.description]
        rows = cur.fetchall()
    finally:
        store.close()
    leak = sequence_leak(columns, rows)
    if leak:
        typer.echo(f"Refusing to print: {leak}. Select genes() columns, or write FASTA with "
                   "write_fasta(path, list(protein_id)).", err=True)
        raise typer.Exit(2)
    if out is not None:
        out.write_text(render(columns, rows, QueryFormat.tsv) + "\n")
    shown = rows if limit <= 0 else rows[:limit]
    typer.echo(render(columns, shown, output_format, total=len(rows)))
