"""Staged writes to a canonical DuckDB database: clone, run, diff, publish.

Any writer (a loader, a caller rerun, a migration) can run against a staged
clone of a dataset database instead of the canonical file:

1. :func:`stage` checkpoints the source and clones it beside itself (an APFS
   clone or reflink where the filesystem supports one, so staging a 15 GB
   database takes milliseconds and no space until pages diverge), and records
   the source's file state in ``<staged>.stage.json``;
2. the writer runs against the staged path;
3. :func:`diff_databases` compares every table of the source and the staged
   copy: schema, row counts, an order-free content hash, and for each changed
   table the rows only in one side (multiset ``EXCEPT ALL``), split into
   added/removed/changed keys when the table has a key;
4. :func:`publish` replaces the source with the staged copy when the source is
   unchanged since staging, holding DuckDB's exclusive lock on the source
   through ``os.replace``, keeps the prior file as a hard-linked backup, and
   re-seals the dataset when a seal exists.

The publish protocol is the one the DefenseFinder reload uses
(:mod:`sharur.defensefinder.refresh`): DuckDB commits land in ``<db>.wal``
before a checkpoint touches the database file, so the change check covers the
size and mtime of both files; DuckDB locks with POSIX ``fcntl`` record locks,
which stay on the original inode through the rename; the backup link is made
before the lock because closing any descriptor of a file drops the process's
``fcntl`` locks on it.
"""

from __future__ import annotations

import json
import os
import platform
import re
import shutil
import subprocess
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Sequence

import duckdb

STAGE_SUFFIX = ".stage.json"


class StagingError(RuntimeError):
    """A staging or publish precondition failed; the canonical file is untouched."""


# --------------------------------------------------------------------------- #
# File state, cloning, locking
# --------------------------------------------------------------------------- #
def wal_path(path: Path) -> Path:
    return path.with_name(path.name + ".wal")


def db_state(path: Path) -> tuple[int, int, int, int]:
    """Size and mtime of the database file and of its WAL ((0, 0) when absent)."""
    stat = path.stat()
    try:
        wal = wal_path(path).stat()
        wal_state = (wal.st_size, wal.st_mtime_ns)
    except FileNotFoundError:
        wal_state = (0, 0)
    return stat.st_size, stat.st_mtime_ns, *wal_state


def clone_file(src: Path, dst: Path) -> str:
    """Copy ``src`` to ``dst`` as a copy-on-write clone when possible; return the method used."""
    dst.unlink(missing_ok=True)
    system = platform.system()
    if system == "Darwin":
        if subprocess.run(["cp", "-c", "-p", str(src), str(dst)], capture_output=True).returncode == 0:
            return "clone"
    elif system == "Linux":
        if subprocess.run(["cp", "--reflink=always", "-p", str(src), str(dst)],
                          capture_output=True).returncode == 0:
            return "reflink"
    shutil.copy2(src, dst)
    return "copy"


def exclusive_lock(db: Path, config: dict | None = None) -> duckdb.DuckDBPyConnection:
    """Open ``db`` read-write, which takes DuckDB's exclusive lock; the caller closes it."""
    try:
        return duckdb.connect(str(db), config=config or {})
    except duckdb.IOException as exc:
        raise StagingError(f"{db} is open in another process; close every connection to it: {exc}") from exc


def checkpoint(db: Path, config: dict | None = None) -> None:
    """Fold a pending WAL into ``db`` (needs the exclusive lock)."""
    con = exclusive_lock(db, config)
    try:
        con.execute("CHECKPOINT")
    finally:
        con.close()
    if wal_path(db).exists() and wal_path(db).stat().st_size > 0:
        raise StagingError(f"{wal_path(db)} persists after CHECKPOINT; resolve it first")


# --------------------------------------------------------------------------- #
# Diff
# --------------------------------------------------------------------------- #
@dataclass
class TableDiff:
    table: str
    status: str  # unchanged | changed | added | removed
    rows_before: int | None = None
    rows_after: int | None = None
    columns_added: list[str] = field(default_factory=list)
    columns_removed: list[str] = field(default_factory=list)
    columns_retyped: list[str] = field(default_factory=list)
    rows_only_before: int | None = None  # multiset rows (common columns) only in the source
    rows_only_after: int | None = None
    key: list[str] | None = None
    keys_added: int | None = None
    keys_removed: int | None = None
    keys_changed: int | None = None
    example_keys: dict[str, list[Any]] = field(default_factory=dict)
    changed_sources: list[str] | None = None  # tables with a `source` column: sources of differing rows


@dataclass
class DatabaseDiff:
    before: str
    after: str
    tables: list[TableDiff]
    seconds: float = 0.0

    def changed(self) -> list[TableDiff]:
        return [t for t in self.tables if t.status != "unchanged"]

    def to_dict(self) -> dict[str, Any]:
        return {"before": self.before, "after": self.after, "seconds": round(self.seconds, 1),
                "changed_tables": [t.table for t in self.changed()],
                "tables": [asdict(t) for t in self.tables]}

    def to_markdown(self) -> str:
        changed = self.changed()
        lines = [f"# Database diff: {len(changed)} of {len(self.tables)} tables differ ({self.seconds:.0f} s)",
                 "", f"- before: `{self.before}`", f"- after: `{self.after}`", ""]
        if not changed:
            lines.append("Every table is identical (schema, rows and content).")
            return "\n".join(lines)
        lines += ["| table | status | rows before | rows after | only before | only after | keys +/-/~ | schema |",
                  "|---|---|---:|---:|---:|---:|---|---|"]
        for t in changed:
            keys = (f"{t.keys_added:,}/{t.keys_removed:,}/{t.keys_changed:,} on {','.join(t.key)}"
                    if t.keys_added is not None else "")
            schema = "; ".join(filter(None, [
                f"+{','.join(t.columns_added)}" if t.columns_added else "",
                f"-{','.join(t.columns_removed)}" if t.columns_removed else "",
                f"~{','.join(t.columns_retyped)}" if t.columns_retyped else ""]))
            lines.append(f"| {t.table} | {t.status} | {_n(t.rows_before)} | {_n(t.rows_after)} | "
                         f"{_n(t.rows_only_before)} | {_n(t.rows_only_after)} | {keys} | {schema} |")
        unchanged = [t.table for t in self.tables if t.status == "unchanged"]
        sourced = [(t.table, t.changed_sources) for t in changed if t.changed_sources]
        if sourced:
            lines += [""] + [f"- {table}: differing rows come from source(s) {', '.join(s)}"
                             for table, s in sourced]
        lines += ["", f"Unchanged ({len(unchanged)}): " + ", ".join(unchanged)]
        examples = [(t.table, k, v) for t in changed for k, v in t.example_keys.items() if v]
        if examples:
            lines += ["", "Example keys:"]
            lines += [f"- {table} {kind}: {', '.join(map(str, vals))}" for table, kind, vals in examples]
        return "\n".join(lines)


def _n(value: int | None) -> str:
    return "" if value is None else f"{value:,}"


def _q(name: str) -> str:
    return '"' + name.replace('"', '""') + '"'


def _schema(con, alias: str) -> dict[str, dict[str, str]]:
    rows = con.execute(
        "SELECT table_name, column_name, data_type FROM information_schema.columns "
        "WHERE table_catalog = ? AND table_schema = 'main' AND table_name IN "
        "(SELECT table_name FROM information_schema.tables WHERE table_catalog = ? "
        " AND table_schema = 'main' AND table_type = 'BASE TABLE') ORDER BY table_name, ordinal_position",
        [alias, alias]).fetchall()
    out: dict[str, dict[str, str]] = {}
    for table, column, dtype in rows:
        out.setdefault(table, {})[column] = dtype
    return out


def _primary_key(con, alias: str, table: str) -> list[str] | None:
    rows = con.execute(
        "SELECT constraint_column_names FROM duckdb_constraints() WHERE database_name = ? "
        "AND schema_name = 'main' AND table_name = ? AND constraint_type IN ('PRIMARY KEY', 'UNIQUE') "
        "ORDER BY constraint_type = 'PRIMARY KEY' DESC, len(constraint_column_names)",
        [alias, table]).fetchall()
    return list(rows[0][0]) if rows else None


# Columns that never appear in example output (sequence content stays out of reports).
_PRIVATE_COLUMNS = re.compile(r"sequence|seq$|_seq\b|translation|residues", re.I)


def diff_databases(
    before: Path | str,
    after: Path | str,
    *,
    tables: Sequence[str] | None = None,
    keys: dict[str, list[str]] | None = None,
    examples: int = 0,
    threads: int = 4,
    memory_limit: str = "8GB",
) -> DatabaseDiff:
    """Compare every base table of two DuckDB files (opened read-only)."""
    started = time.time()
    con = duckdb.connect(":memory:")
    con.execute(f"SET threads = {int(threads)}")
    con.execute(f"SET memory_limit = '{memory_limit}'")
    con.execute(f"ATTACH '{str(before).replace(chr(39), chr(39) * 2)}' AS b (READ_ONLY)")
    con.execute(f"ATTACH '{str(after).replace(chr(39), chr(39) * 2)}' AS a (READ_ONLY)")
    try:
        sb, sa = _schema(con, "b"), _schema(con, "a")
        names = sorted(set(sb) | set(sa))
        if tables:
            names = [t for t in names if t in set(tables)]
        out: list[TableDiff] = []
        for table in names:
            if table not in sa:
                rows = con.execute(f"SELECT COUNT(*) FROM b.{_q(table)}").fetchone()[0]
                out.append(TableDiff(table, "removed", rows_before=rows))
                continue
            if table not in sb:
                rows = con.execute(f"SELECT COUNT(*) FROM a.{_q(table)}").fetchone()[0]
                out.append(TableDiff(table, "added", rows_after=rows))
                continue
            cb, ca = sb[table], sa[table]
            common = [c for c in cb if c in ca]
            diff = TableDiff(
                table, "unchanged",
                columns_added=[c for c in ca if c not in cb],
                columns_removed=[c for c in cb if c not in ca],
                columns_retyped=[c for c in common if cb[c] != ca[c]],
            )
            cols = ", ".join(_q(c) for c in common)
            packed = "struct_pack(" + ", ".join(f"{_q(c)} := t.{_q(c)}" for c in common) + ")"
            fp = {}
            for alias in ("b", "a"):
                fp[alias] = con.execute(
                    f"SELECT COUNT(*), COALESCE(SUM(hash({packed})::HUGEINT), 0) FROM {alias}.{_q(table)} t"
                ).fetchone()
            diff.rows_before, diff.rows_after = int(fp["b"][0]), int(fp["a"][0])
            schema_changed = bool(diff.columns_added or diff.columns_removed or diff.columns_retyped)
            if fp["b"] == fp["a"] and not schema_changed:
                out.append(diff)
                continue
            diff.status = "changed"
            if fp["b"] != fp["a"]:
                diff.rows_only_before = int(con.execute(
                    f"SELECT COUNT(*) FROM (SELECT {cols} FROM b.{_q(table)} EXCEPT ALL "
                    f"SELECT {cols} FROM a.{_q(table)})").fetchone()[0])
                diff.rows_only_after = int(con.execute(
                    f"SELECT COUNT(*) FROM (SELECT {cols} FROM a.{_q(table)} EXCEPT ALL "
                    f"SELECT {cols} FROM b.{_q(table)})").fetchone()[0])
                if "source" in common:
                    diff.changed_sources = sorted(str(r[0]) for r in con.execute(
                        f"SELECT DISTINCT source FROM ((SELECT {cols} FROM b.{_q(table)} EXCEPT ALL "
                        f"SELECT {cols} FROM a.{_q(table)}) UNION ALL (SELECT {cols} FROM a.{_q(table)} "
                        f"EXCEPT ALL SELECT {cols} FROM b.{_q(table)}))").fetchall())
                key = (keys or {}).get(table) or _primary_key(con, "a", table)
                if key and all(k in common for k in key):
                    diff.key = key
                    kc = ", ".join(_q(k) for k in key)
                    join = " AND ".join(f"x.{_q(k)} IS NOT DISTINCT FROM y.{_q(k)}" for k in key)
                    diff.keys_added = int(con.execute(
                        f"SELECT COUNT(*) FROM (SELECT DISTINCT {kc} FROM a.{_q(table)} EXCEPT "
                        f"SELECT DISTINCT {kc} FROM b.{_q(table)})").fetchone()[0])
                    diff.keys_removed = int(con.execute(
                        f"SELECT COUNT(*) FROM (SELECT DISTINCT {kc} FROM b.{_q(table)} EXCEPT "
                        f"SELECT DISTINCT {kc} FROM a.{_q(table)})").fetchone()[0])
                    rowhash = "hash(" + packed.replace("t.", "x.") + ")"
                    rowhash_y = "hash(" + packed.replace("t.", "y.") + ")"
                    diff.keys_changed = int(con.execute(
                        f"SELECT COUNT(*) FROM (SELECT {kc}, {rowhash} h FROM b.{_q(table)} x) x "
                        f"JOIN (SELECT {kc}, {rowhash_y} h FROM a.{_q(table)} y) y ON {join} "
                        f"WHERE x.h <> y.h").fetchone()[0])
                    if examples:
                        safe = [k for k in key if not _PRIVATE_COLUMNS.search(k)]
                        if safe:
                            sk = ", ".join(_q(k) for k in safe)
                            for kind, (left, right) in {"added": ("a", "b"), "removed": ("b", "a")}.items():
                                vals = con.execute(
                                    f"SELECT {sk} FROM (SELECT DISTINCT {kc} FROM {left}.{_q(table)} EXCEPT "
                                    f"SELECT DISTINCT {kc} FROM {right}.{_q(table)}) ORDER BY ALL LIMIT {int(examples)}"
                                ).fetchall()
                                diff.example_keys[kind] = [v[0] if len(v) == 1 else list(v) for v in vals]
            out.append(diff)
        return DatabaseDiff(str(before), str(after), out, time.time() - started)
    finally:
        con.close()


# --------------------------------------------------------------------------- #
# Stage and publish
# --------------------------------------------------------------------------- #
@dataclass
class Stage:
    source: str
    staged: str
    source_state: list[int]
    method: str
    created: str
    label: str

    @property
    def manifest(self) -> Path:
        return Path(self.staged + STAGE_SUFFIX)

    def save(self) -> None:
        self.manifest.write_text(json.dumps(asdict(self), indent=2) + "\n")

    @classmethod
    def load(cls, staged: Path | str) -> "Stage":
        path = Path(str(staged) + STAGE_SUFFIX)
        if not path.is_file():
            raise StagingError(f"{path} not found; stage the database with `sharur stage run` first")
        return cls(**json.loads(path.read_text()))


def stage(db: Path | str, *, work_dir: Path | str | None = None, label: str = "staged") -> Stage:
    """Checkpoint ``db`` and clone it to ``<work_dir>/<db name>.<label>.staging``."""
    db = Path(db).expanduser().resolve()
    if not db.is_file():
        raise StagingError(f"{db} does not exist")
    if not re.fullmatch(r"[\w.-]+", label):
        raise StagingError(f"label {label!r}: use letters, digits, '.', '_' or '-'")
    checkpoint(db)
    stage_dir = Path(work_dir).expanduser().resolve() if work_dir else db.parent
    stage_dir.mkdir(parents=True, exist_ok=True)
    staged = stage_dir / f"{db.name}.{label}.staging"
    for stale in (staged, wal_path(staged), Path(str(staged) + STAGE_SUFFIX)):
        stale.unlink(missing_ok=True)
    state = db_state(db)
    method = clone_file(db, staged)
    if method == "copy" and shutil.disk_usage(stage_dir).free < db.stat().st_size * 0.2:
        staged.unlink(missing_ok=True)
        raise StagingError(f"Less than 20% headroom in {stage_dir} after copying {db}")
    if db_state(db) != state:
        staged.unlink(missing_ok=True)
        raise StagingError(f"{db} changed while it was being staged")
    result = Stage(source=str(db), staged=str(staged), source_state=list(state), method=method,
                   created=time.strftime("%Y-%m-%dT%H:%M:%S%z"), label=label)
    result.save()
    return result


@dataclass
class PublishResult:
    source: str
    backup: str
    resealed: bool
    seal_error: str | None = None


def discard(stage_: Stage) -> None:
    staged = Path(stage_.staged)
    for path in (staged, wal_path(staged), stage_.manifest):
        path.unlink(missing_ok=True)


def publish(stage_: Stage, *, reseal: bool = True, backup: bool = True) -> PublishResult:
    """Swap the staged copy into place; see the module docstring for the protocol."""
    db, staged = Path(stage_.source), Path(stage_.staged)
    if not staged.is_file():
        raise StagingError(f"staged copy {staged} is gone")
    if os.stat(staged.parent).st_dev != os.stat(db.parent).st_dev:
        raise StagingError(f"{staged} is on another filesystem than {db}; publishing needs an atomic rename")
    if tuple(db_state(db)) != tuple(stage_.source_state):
        raise StagingError(f"{db} changed after it was staged; re-stage and rerun the writer")
    checkpoint(staged)
    backup_path = db.with_name(f"{db.name}.pre-{stage_.label}-{time.strftime('%Y%m%dT%H%M%S')}")
    if backup:
        try:
            os.link(db, backup_path)
        except OSError:
            clone_file(db, backup_path)
    swapped = False
    try:
        lock = exclusive_lock(db)
        try:
            if wal_path(db).exists() or tuple(db_state(db)) != tuple(stage_.source_state):
                raise StagingError(f"{db} changed after it was staged; staged copy left in place")
            os.replace(staged, db)
            swapped = True
        finally:
            lock.close()
    finally:
        if not swapped and backup:
            backup_path.unlink(missing_ok=True)
    stage_.manifest.unlink(missing_ok=True)
    result = PublishResult(source=str(db), backup=str(backup_path) if backup else "", resealed=False)
    seal = db.parent / "dataset.seal.json"
    if reseal and seal.exists():
        from sharur.dataset_seal import seal_dataset  # noqa: PLC0415

        try:
            previous = json.loads(seal.read_text()).get("policy", {})
            seal_dataset(db, output_path=seal, force=True,
                         hash_large_files=bool(previous.get("hash_large_files", False)))
            result.resealed = True
        except Exception as exc:  # noqa: BLE001 - the swap happened; report the seal failure
            result.seal_error = f"{type(exc).__name__}: {exc}"
    return result


def touching_regex(diff: DatabaseDiff) -> str | None:
    """A ``sharur reverify --touching`` pattern for what changed.

    A table with a ``source`` column (``annotations``) contributes the source values of its
    differing rows, since nearly every finding reads that table; other tables contribute
    their names.
    """
    terms = []
    for t in diff.changed():
        if t.changed_sources:
            terms += [re.escape(s) for s in t.changed_sources]
        else:
            terms.append(re.escape(t.table))
    if not terms:
        return None
    return r"\b(" + "|".join(sorted(set(terms))) + r")\b"
