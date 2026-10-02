"""Browser startup summaries, computed once per database file in a short-lived worker process.

The catalog's startup aggregations (per-genome function categories, annotated-protein counts, Pfam
domain, discovery, Compare-page and synteny summaries) return small tables but build large DuckDB hash tables, and the process
keeps that memory: about 5 GB for a 3.6M-protein dataset. A worker process runs them and saves their
result rows; the server replays the rows and its working set stays near 1 GB.

Rows are keyed by query text and parameters, and the file by the database's size and modification time
plus the Sharur code that issues the queries. A query absent from the file runs live, so a stale or
partial file costs memory, never correctness.
"""

from __future__ import annotations

import hashlib
import logging
import os
import pickle
import subprocess
import sys
import threading
import time
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

VERSION = 1
PACKAGE = Path(__file__).resolve().parents[1]
CODE_SUFFIXES = {".py", ".tsv", ".json", ".yaml", ".yml"}


def default_dir() -> Path:
    return Path(os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache") / "sharur" / "browser"


def cache_file(db_path: str | Path, cache_dir: str | Path) -> Path:
    db = str(Path(db_path).resolve())
    return Path(cache_dir) / f"{Path(db).stem}-{hashlib.sha1(db.encode()).hexdigest()[:12]}.pkl"


def fingerprint(db_path: str | Path) -> dict[str, Any]:
    st = Path(db_path).resolve().stat()
    code = hashlib.sha1()
    for f in sorted(PACKAGE.rglob("*")):
        if f.suffix in CODE_SUFFIXES and "__pycache__" not in f.parts:
            s = f.stat()
            code.update(f"{f.relative_to(PACKAGE)}:{s.st_size}:{s.st_mtime_ns}\n".encode())
    return {"version": VERSION, "db_size": st.st_size, "db_mtime_ns": st.st_mtime_ns, "code": code.hexdigest()}


def _interned(value: Any) -> Any:
    if type(value) is str:
        return sys.intern(value)
    if type(value) is list:
        return [sys.intern(v) if type(v) is str else v for v in value]
    return value


def _key(query: str, params: Any) -> tuple[str, str]:
    return query, repr(params)


class Recorder:
    """Store proxy that keeps every query's rows, plus named summaries computed elsewhere (``artifacts``)."""

    def __init__(self, store) -> None:
        self._store, self.rows, self.artifacts = store, {}, {}

    def execute(self, query: str, params: Any = None) -> list[tuple]:
        rows = self._store.execute(query, params)
        # A few thousand genome ids, KOs and labels repeat across millions of rows; interned, the file
        # stores each once and the server's rows share them.
        self.rows[_key(query, params)] = [tuple(_interned(v) for v in row) for row in rows]
        return rows

    def __getattr__(self, name: str) -> Any:
        return getattr(self._store, name)


class Replayer:
    """Store proxy that answers recorded queries from the cache and runs the rest live."""

    def __init__(self, store, rows: dict, artifacts: dict | None = None) -> None:
        self._store, self._rows, self.misses = store, rows, 0
        self._artifacts = artifacts or {}

    def artifact(self, key: str) -> Any:
        """A named summary the worker computed (handed over once), or None."""
        return self._artifacts.pop(key, None)

    def execute(self, query: str, params: Any = None) -> list[tuple]:
        rows = self._rows.pop(_key(query, params), None)   # each startup query runs once; free rows as served
        if rows is None:
            self.misses += 1
            return self._store.execute(query, params)
        return rows

    def __getattr__(self, name: str) -> Any:
        return getattr(self._store, name)


def build(db_path: str | Path, out: str | Path) -> int:
    """Run the startup summaries against ``db_path`` and write their rows to ``out``; returns the query count."""
    from sharur.browser.catalog import load_background, load_catalog  # noqa: PLC0415
    from sharur.storage.duckdb_store import DuckDBStore  # noqa: PLC0415

    stamp = fingerprint(db_path)
    store = DuckDBStore(str(db_path), read_only=True)
    try:
        from sharur.browser.routes_compare import PfamPairs  # noqa: PLC0415

        recorder, lock = Recorder(store), threading.Lock()
        catalog = load_catalog(recorder)
        load_background(recorder, catalog, lock)
        if catalog.status != "ready":
            raise RuntimeError(f"startup summaries did not finish: {catalog.status}")
        PfamPairs().load(recorder, catalog, lock)
        _synteny_summary(db_path, store, catalog, lock, recorder.artifacts)
    finally:
        store.close()
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(out.name + f".{os.getpid()}.tmp")
    with open(tmp, "wb") as fh:
        pickle.dump({"fingerprint": stamp, "rows": recorder.rows, "artifacts": recorder.artifacts}, fh,
                    protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tmp, out)
    return len(recorder.rows)


def _synteny_summary(db_path, store, catalog, lock, artifacts: dict) -> None:
    """The synteny sidecar's overview summary, when the dataset has one (its key carries the sidecar's state)."""
    from types import SimpleNamespace  # noqa: PLC0415

    from sharur.browser.routes_synteny import open_view  # noqa: PLC0415

    view = open_view(SimpleNamespace(db_path=Path(db_path), store=store, lock=lock, catalog=catalog))
    if view is not None:
        artifacts[view.summary_key()] = view._summarize()


def load(db_path: str | Path, cache_dir: str | Path) -> tuple[dict, dict] | None:
    """Recorded rows and named summaries for this database and code, or None."""
    path = cache_file(db_path, cache_dir)
    try:
        with open(path, "rb") as fh:
            data = pickle.load(fh)
    except (OSError, EOFError, pickle.UnpicklingError, AttributeError, ImportError):
        return None
    if data.get("fingerprint") != fingerprint(db_path):
        return None
    return data["rows"], data.get("artifacts", {})


def ensure(db_path: str | Path, cache_dir: str | Path) -> tuple[dict, dict] | None:
    """Recorded rows and summaries, built in a worker process first when absent or stale; None if that fails."""
    cached = load(db_path, cache_dir)
    if cached is not None:
        return cached
    out = cache_file(db_path, cache_dir)
    logger.warning("Summarizing %s for the browser (once per database version)...", Path(db_path).name)
    t0 = time.time()
    done = subprocess.run([sys.executable, "-m", "sharur.browser.startup_cache", str(db_path), str(out)],
                          capture_output=True, text=True)
    if done.returncode != 0:
        logger.warning("Startup summary worker failed; summarizing in the server instead.\n%s", done.stderr[-2000:])
        return None
    logger.warning("Summaries ready in %.0fs.", time.time() - t0)
    return load(db_path, cache_dir)


if __name__ == "__main__":
    build(sys.argv[1], sys.argv[2])
