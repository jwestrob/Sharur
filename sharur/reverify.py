"""Replay the verification triples recorded in findings files.

Every finding carries ``verification``: a list of ``{claim, query, expected}``
triples (``docs/findings_spec.md``). This module re-executes them and compares
each result with ``expected``, so a refreshed database, a new caller or a
corrected map shows exactly which recorded numbers still hold.

Query forms, as findings record them:

- SQL. With ``attach`` (``{alias: path}``) the query runs in an in-memory
  connection with every alias attached read-only; otherwise it runs against the
  dataset database (``data/<dataset>/sharur.duckdb`` for a findings file under
  ``data/<dataset>/``, or ``--db``).
- ``python: EXPR``, and bare Python expressions that open ``Sharur(...)``,
  evaluated from the repository root with ``json``, ``pd``, ``np``, ``duckdb``,
  ``Sharur``, ``os``, ``re``, ``math`` and ``Path`` in scope.
- ``shell: CMD`` and ``python -c ...`` commands, run by ``/bin/sh`` from the
  repository root; stdout is parsed as JSON when it parses.

Every DuckDB connection a replay opens is read-only and carries the configured
thread, memory and spill limits. That holds for Python evaluated in-process and
for Python processes a shell triple starts, through a ``sitecustomize`` shim.

Path substitution (``OLD=NEW``) redirects every database or file path a replay
touches: query text, ``attach`` maps, ``duckdb.connect`` arguments and
``ATTACH`` statements issued by helper scripts. Compare mode runs each triple
twice, once under the baseline substitutions (none by default: the paths as
recorded) and once under the replay substitutions, and reports every triple
whose result changes. A staged refresh, or a backup taken before a publish, can then
be checked against the live database before anyone reads the new numbers.

Status per triple: ``OK``, ``FAIL`` (result differs from ``expected``),
``ERROR`` (the replay raised) or ``NARRATIVE`` (the recorded query is a note,
such as a methods reference or a command annotated with the file its value was
read from, and is reported without being run). Integers, strings, booleans and
lists compare exactly; a float matches when the observed value, rounded to the
precision written in ``expected``, equals it.

Replays are read-only for files too: Python code a replay runs may not write
inside the repository (an audit hook refuses it), and on macOS shell triples run
under ``sandbox-exec`` with writes inside the repository denied.
"""

from __future__ import annotations

import ast
import concurrent.futures as cf
import contextlib
import decimal
import json
import logging
import multiprocessing
import os
import re
import shlex
import shutil
import signal
import subprocess
import sys
import tempfile
import threading
import time
from dataclasses import asdict, dataclass, field
from decimal import ROUND_HALF_EVEN, ROUND_HALF_UP, Decimal
from pathlib import Path
from typing import Any, Iterable, Sequence

STATUSES = ("OK", "FAIL", "ERROR", "NARRATIVE")
logger = logging.getLogger(__name__)
_SQL_START = re.compile(r"(?is)^\s*(?:--[^\n]*\n\s*)*\(?\s*(select|with|from|pragma|describe|summarize|show|values)\b")
_PY_CMD = re.compile(r"^\s*(?:\S*/)?python3?\s+(-c|-m|\S+\.py)\b")
_IMPORT_NAME = re.compile(r"__import__\(\s*['\"]([A-Za-z_][\w.]*)['\"]")
_PATH_INSERT = re.compile(r"sys\.path\.insert\(\s*\d+\s*,\s*['\"]([^'\"]+)['\"]")
_PY_PATH = re.compile(r"[\w./\-]+\.py\b")
_LOCAL_IMPORT = re.compile(r"^\s*(?:from\s+([A-Za-z_]\w*)\s+import|import\s+([A-Za-z_]\w*))", re.M)


# --------------------------------------------------------------------------- #
# Triples
# --------------------------------------------------------------------------- #
@dataclass
class Triple:
    """One verification triple with its finding context."""

    source: str  # findings file, relative to root when possible
    line: int  # 1-based line in the findings file
    finding_id: str
    title: str
    index: int  # 1-based position in the finding's verification list
    claim: str
    query: str
    expected: Any
    attach: dict[str, str] | None
    default_db: str | None
    kind: str = ""
    # statistical checks (sharur/stat_checks.py): per-value tolerances and the recorded claim
    tolerance: dict | None = None
    stat: dict | None = None
    # what a replay runs (the query without its form prefix or annotation), the JSON key it
    # reads from a command's output, and why a narrative triple is not replayed
    command: str = ""
    select: str | None = None
    note: str | None = None

    @property
    def key(self) -> str:
        return f"{self.source}:{self.line}#{self.index}"


@dataclass
class Outcome:
    """Result of one replay."""

    status: str
    observed: Any = None
    error: str | None = None
    seconds: float = 0.0


@dataclass
class TripleResult:
    triple: Triple
    flags: list[str]
    baseline: Outcome | None = None  # compare mode: run under the baseline substitutions
    replay: Outcome = field(default_factory=lambda: Outcome("ERROR"))

    @property
    def changed(self) -> bool:
        if self.baseline is None:
            return False
        return (self.baseline.status, _canonical(self.baseline.observed)) != (
            self.replay.status, _canonical(self.replay.observed))


@dataclass
class Reading:
    """How a replay treats one recorded query."""

    kind: str  # sql, python, shell or narrative
    command: str  # what runs: the query without its form prefix or annotation
    select: str | None = None  # JSON key read from a command's output
    note: str | None = None  # why a narrative query is not run


# "Label: text" where the label is not one of the replayable prefixes (``shell:``, ``python:``).
_LABEL = re.compile(r"^(?P<label>[A-Za-z][\w/ ().-]{0,40}?)\s*:\s+(?P<body>.+)$", re.S)
# A command annotated with what it prints: "(prints JSON)", "(prints JSON; key 'K')".
_PRINTS_JSON = re.compile(
    r"\s*\(\s*prints JSON\s*(?:;\s*(?:key\s+'(?P<key>[^']+)'|then index\s+(?P<index>[^\s)]+)))?\s*\)?\s*$")
# Commands a replay runs without a ``shell:`` prefix: interpreters and read-only text tools.
_TEXT_TOOLS = frozenset({
    "awk", "bc", "cat", "column", "comm", "cut", "echo", "egrep", "expr", "grep", "gunzip", "gzip",
    "head", "jq", "join", "paste", "printf", "sed", "sort", "tail", "test", "tr", "true", "uniq", "wc",
    "zcat", "zgrep"})
_PY_EXE = re.compile(r"^(?:\S*/)?python(?:3(?:\.\d+)?)?$")
_ASSIGN = re.compile(r"^[A-Za-z_]\w*=")
_SEPARATORS = frozenset({";", "&&", "||", "|"})
_SHELL_PUNCT = set("();<>|&")


def _shell_tokens(text: str) -> list[str]:
    """Shell words and operators; quoted words keep their quotes (``'>'`` is not a redirect)."""
    lex = shlex.shlex(text, posix=False, punctuation_chars=True)
    lex.whitespace_split = True
    lex.commenters = ""
    return list(lex)


def _redirect_targets(tokens: Sequence[str]) -> list[str]:
    """Files an output redirection writes (``>``, ``>>``, ``&>``); fd copies and /dev/null excluded."""
    out = []
    for i, tok in enumerate(tokens):
        if set(tok) <= _SHELL_PUNCT and ">" in tok:
            target = tokens[i + 1].strip("'\"") if i + 1 < len(tokens) else ""
            if (tok.endswith("&") and target.isdigit()) or target == "/dev/null":
                continue
            out.append(target)
    return out


def _command_problem(text: str) -> str | None:
    """``None`` when ``text`` is a plain command line a replay may run without a ``shell:`` prefix.

    Such a line chains (``;``, ``&&``, ``||``, ``|``) ``cd``, Python interpreters and read-only
    text tools, and writes no file. Anything else (parenthesized remarks, ``#`` comments,
    ``-> file`` arrows, which the shell reads as redirections, other programs) marks a note.
    """
    try:
        tokens = _shell_tokens(text)
    except ValueError as exc:
        return f"unbalanced quoting ({exc})"
    if not tokens:
        return "empty"
    if any(tok.startswith("#") for tok in tokens):
        return "a remark follows the command (#)"
    targets = _redirect_targets(tokens)
    if targets:
        return f"it names an output file ({targets[0]})"
    start = True
    for tok in tokens:
        if set(tok) <= _SHELL_PUNCT:
            if "(" in tok or ")" in tok:
                return "a parenthesized remark follows the command"
            if tok in _SEPARATORS:
                start = True
            elif not ("<" in tok or ">" in tok):
                return f"unsupported shell operator {tok!r}"
            continue
        if start:
            if _ASSIGN.match(tok):
                continue
            word = tok.strip("'\"")
            if not (word == "cd" or os.path.basename(word) in _TEXT_TOOLS or _PY_EXE.match(word)):
                return f"{word[:40]!r} is not an interpreter or a read-only text tool"
            start = False
    return None


_PARSER: Any = None


def _parses_as_sql(query: str) -> bool:
    global _PARSER
    import duckdb  # noqa: PLC0415

    if _PARSER is None:
        _PARSER = duckdb.connect(":memory:")
    try:
        return bool(_PARSER.extract_statements(query))
    except Exception:  # noqa: BLE001 - any parser error means "not SQL"
        return False


def interpret(query: str, attach: dict | None = None, check: str | None = None) -> Reading:
    """The replay form of a recorded query, or ``narrative`` for a note.

    Order: an explicit ``check: "narrative"``; the ``shell:`` and ``python:`` prefixes; a Python
    command line; SQL (``attach``, or a leading SQL keyword); Python opening ``Sharur(...)``;
    a labelled note (``Python: ...``), run when its text is a plain command line, optionally
    annotated ``(prints JSON; key 'K')``; text DuckDB parses as SQL; a plain command line of
    text tools. Everything else is a note.
    """
    q = query.strip()
    if check == "narrative":
        return Reading("narrative", q, note="recorded as narrative (check: narrative)")
    if q.startswith("shell:"):
        return Reading("shell", q[len("shell:"):].strip())
    if q.startswith("python:"):
        return Reading("python", q[len("python:"):].strip())
    if _PY_CMD.match(q):
        problem = _command_problem(q)
        return Reading("shell", q) if problem is None else Reading("narrative", q, note=f"command note: {problem}")
    if attach or _SQL_START.match(q):
        return Reading("sql", q)
    if "Sharur(" in q or "__import__" in q:
        return Reading("python", q)
    label = _LABEL.match(q)
    if label:
        body, select = label.group("body").strip(), None
        annotated = _PRINTS_JSON.search(body)
        if annotated:
            body = body[:annotated.start()].strip()
            select = annotated.group("key") or annotated.group("index")
        problem = _command_problem(body)
        if problem is None:
            return Reading("shell", body, select=select)
        return Reading("narrative", q, note=f"labelled note ({label.group('label')}:): {problem}")
    if _parses_as_sql(q):
        return Reading("sql", q)
    problem = _command_problem(q)
    if problem is None:
        return Reading("shell", q)
    return Reading("narrative", q, note=f"note: neither SQL nor a command line ({problem})")


def classify(query: str, attach: dict | None = None) -> str:
    """``sql``, ``python``, ``shell`` or ``narrative`` for a recorded query string."""
    return interpret(query, attach).kind


def find_root(path: Path) -> Path:
    """Repository root for a findings file: the parent of its ``data`` ancestor."""
    path = path.resolve()
    for parent in path.parents:
        if parent.name == "data":
            return parent.parent
    return Path.cwd().resolve()


def dataset_of(path: Path, root: Path) -> str | None:
    try:
        rel = path.resolve().relative_to(root / "data")
    except ValueError:
        return None
    return rel.parts[0] if rel.parts else None


def resolve_findings(paths: Sequence[Path]) -> list[Path]:
    """Findings files from file or dataset-directory arguments."""
    out: list[Path] = []
    for p in paths:
        p = Path(p)
        if p.is_dir():
            candidate = p / "exploration" / "findings.jsonl"
            if not candidate.is_file():
                raise FileNotFoundError(f"{p} has no exploration/findings.jsonl")
            out.append(candidate)
        elif p.is_file():
            out.append(p)
        else:
            raise FileNotFoundError(p)
    return out


def load_triples(findings: Sequence[Path], root: Path, db: Path | None = None) -> list[Triple]:
    triples: list[Triple] = []
    for path in findings:
        ds = dataset_of(path, root)
        default_db = None
        if db is not None:
            default_db = str(Path(db).resolve())
        elif ds and (root / "data" / ds / "sharur.duckdb").is_file():
            default_db = str(root / "data" / ds / "sharur.duckdb")
        try:
            source = str(path.resolve().relative_to(root))
        except ValueError:
            source = str(path)
        with path.open() as handle:
            for line_no, line in enumerate(handle, 1):
                if not line.strip():
                    continue
                record = json.loads(line)
                for index, item in enumerate(record.get("verification") or [], 1):
                    if not isinstance(item, dict) or not isinstance(item.get("query"), str):
                        continue
                    attach = item.get("attach") or None
                    reading = interpret(item["query"], attach, item.get("check"))
                    triples.append(Triple(
                        source=source, line=line_no,
                        finding_id=str(record.get("id") or record.get("finding_id") or ""),
                        title=str(record.get("title") or ""), index=index,
                        claim=str(item.get("claim") or ""), query=item["query"],
                        expected=item.get("expected"), attach=attach, default_db=default_db,
                        kind=reading.kind,
                        tolerance=item.get("tolerance"), stat=item.get("statistical_claim"),
                        command=reading.command, select=reading.select, note=reading.note,
                    ))
    return triples


# --------------------------------------------------------------------------- #
# Dependency text: what a triple reads, for --touching selection
# --------------------------------------------------------------------------- #
def dependency_text(triple: Triple, root: Path, _cache: dict[Path, str] | None = None) -> str:
    """The query plus the source of local scripts it runs or imports (two levels)."""
    cache = _cache if _cache is not None else {}
    parts = [triple.query, json.dumps(triple.attach or {})]
    dirs = [root / d if not Path(d).is_absolute() else Path(d) for d in _PATH_INSERT.findall(triple.query)]
    files: list[Path] = []
    for match in _PY_PATH.findall(triple.query):
        candidate = Path(match) if Path(match).is_absolute() else root / match
        if candidate.is_file():
            files.append(candidate)
            dirs.append(candidate.parent)
    for name in _IMPORT_NAME.findall(triple.query):
        for d in dirs:
            candidate = d / (name.replace(".", "/") + ".py")
            if candidate.is_file():
                files.append(candidate)
                break
    seen: set[Path] = set()
    frontier = files
    for _depth in range(2):
        nxt: list[Path] = []
        for f in frontier:
            f = f.resolve()
            if f in seen:
                continue
            seen.add(f)
            if f not in cache:
                try:
                    cache[f] = f.read_text(errors="replace")
                except OSError:
                    cache[f] = ""
            text = cache[f]
            parts.append(text)
            for a, b in _LOCAL_IMPORT.findall(text):
                mod = a or b
                for d in {f.parent, *dirs, f.parent / "scripts"}:
                    candidate = d / f"{mod}.py"
                    if candidate.is_file():
                        nxt.append(candidate)
        frontier = nxt
    return "\n".join(parts)


# --------------------------------------------------------------------------- #
# Path substitution
# --------------------------------------------------------------------------- #
# A path matches only as a whole path: never inside a longer name (``sharur.duckdb`` inside
# ``sharur.duckdb.pre-reload``) and a relative path never inside an absolute one. That keeps
# a rewrite idempotent when the replacement path starts with the original.
_BOUNDARY_BEFORE_ABS = r"(?<![\w.-])"
_BOUNDARY_BEFORE_REL = r"(?<![\w./-])"
_BOUNDARY_AFTER = r"(?![\w.-])"


def _boundary_pattern(variants: Sequence[str]) -> re.Pattern | None:
    if not variants:
        return None
    return re.compile("|".join(
        (_BOUNDARY_BEFORE_ABS if v.startswith("/") else _BOUNDARY_BEFORE_REL) + re.escape(v) + _BOUNDARY_AFTER
        for v in variants))


class Substitution:
    """Rewrites paths in text and resolves ``duckdb.connect`` targets."""

    def __init__(self, pairs: Sequence[tuple[str, str]], root: Path):
        self.root = root
        self.map: dict[str, str] = {}
        variants: dict[str, str] = {}
        for old, new in pairs:
            old_abs = (root / old).resolve() if not Path(old).is_absolute() else Path(old).resolve()
            new_abs = str((root / new).resolve() if not Path(new).is_absolute() else Path(new).resolve())
            self.map[str(old_abs)] = new_abs
            variants[str(old_abs)] = new_abs
            try:
                variants[str(old_abs.relative_to(root))] = new_abs
            except ValueError:
                pass
            variants[old] = new_abs
        ordered = sorted(variants, key=len, reverse=True)
        self._variants = variants
        self._pattern = _boundary_pattern(ordered)

    def __bool__(self) -> bool:
        return bool(self.map)

    def text(self, value: str) -> str:
        if not self._pattern or not isinstance(value, str):
            return value
        return self._pattern.sub(lambda m: self._variants[m.group(0)], value)

    def path(self, value: Any) -> Any:
        if not self.map or not isinstance(value, (str, os.PathLike)):
            return value
        raw = os.fspath(value)
        if raw in ("", ":memory:") or raw.startswith(("md:", "motherduck:")):
            return value
        p = Path(raw)
        resolved = str((self.root / p).resolve() if not p.is_absolute() else p.resolve())
        return self.map.get(resolved, value)

    def env(self) -> str:
        return json.dumps({"root": str(self.root), "map": self.map, "variants": self._variants})


# --------------------------------------------------------------------------- #
# Read-only, limited duckdb.connect (in-process and for shell children)
# --------------------------------------------------------------------------- #
_SHIM = r'''
"""sharur reverify: read-only, limited DuckDB connections for replayed commands."""
import json as _json, os as _os, re as _re, sys as _sys
_cfg = _json.loads(_os.environ.get("SHARUR_REVERIFY_SHIM") or "{}")
if _cfg.get("guard") and not getattr(_sys, "_sharur_reverify_guard", None):
    # Replays are read-only for files too: refuse Python-level writes inside the guarded paths.
    _sys.dont_write_bytecode = True
    _guard = tuple({f(p) for p in _cfg["guard"] for f in (_os.path.abspath, _os.path.realpath)})
    _WFLAGS = _os.O_WRONLY | _os.O_RDWR | _os.O_CREAT | _os.O_TRUNC | _os.O_APPEND
    _TARGETS = {"os.remove": (0,), "os.rmdir": (0,), "os.mkdir": (0,), "os.truncate": (0,),
                "os.rename": (0, 1), "os.symlink": (1,), "os.link": (1,), "shutil.rmtree": (0,),
                "shutil.copyfile": (1,), "shutil.copytree": (1,), "shutil.move": (0, 1)}

    def _guarded(p):
        if not isinstance(p, (str, bytes, _os.PathLike)):
            return None
        try:
            p = _os.path.abspath(_os.fsdecode(p))
        except (TypeError, ValueError):
            return None
        return p if any(p == g or p.startswith(g.rstrip(_os.sep) + _os.sep) for g in _guard) else None

    def _reverify_guard(event, args):
        if event == "open":
            path, mode, flags = (tuple(args) + (None, None, None))[:3]
            if (isinstance(mode, str) and any(c in mode for c in "wax+")) or (
                    isinstance(flags, int) and flags & _WFLAGS):
                hit = _guarded(path)
                if hit:
                    raise PermissionError(13, "sharur reverify: replays are read-only; refused to write", hit)
        elif event in _TARGETS:
            for i in _TARGETS[event]:
                hit = _guarded(args[i]) if i < len(args) else None
                if hit:
                    raise PermissionError(13, f"sharur reverify: replays are read-only; refused {event}", hit)

    _sys.addaudithook(_reverify_guard)
    _sys._sharur_reverify_guard = _guard
if _cfg:
    try:
        import duckdb as _duckdb
    except Exception:  # noqa: BLE001
        _duckdb = None
    if _duckdb is not None and not getattr(_duckdb.connect, "_sharur_reverify", False):
        _real = _duckdb.connect
        _map = _cfg.get("map", {})
        _variants = _cfg.get("variants", {})
        _root = _cfg.get("root", ".")
        _pat = _re.compile("|".join(
            (r"(?<![\w.-])" if v.startswith("/") else r"(?<![\w./-])") + _re.escape(v) + r"(?![\w.-])"
            for v in sorted(_variants, key=len, reverse=True))) if _variants else None

        def _rw(sql):
            return _pat.sub(lambda m: _variants[m.group(0)], sql) if (_pat and isinstance(sql, str)) else sql

        def _mapped(db):
            raw = _os.fspath(db)
            if raw in ("", ":memory:") or raw.startswith(("md:", "motherduck:")):
                return raw, False
            p = raw if _os.path.isabs(raw) else _os.path.join(_root, raw)
            p = _os.path.realpath(p)
            return _map.get(p, raw), True

        # The wrapper's frame sits between the caller and DuckDB, so the connection scans every
        # frame for Python objects named in a query (a DataFrame ``df`` in ``FROM df``), and the
        # wrapper's own locals carry names no query uses.
        class _Conn:
            def __init__(self, con):
                object.__setattr__(self, "_con", con)

            def _wrap(self, r):
                return self if r is self._con else r

            def execute(self, _sharur_q, *_sharur_a, **_sharur_k):
                return self._wrap(self._con.execute(_rw(_sharur_q), *_sharur_a, **_sharur_k))

            def sql(self, _sharur_q, *_sharur_a, **_sharur_k):
                return self._con.sql(_rw(_sharur_q), *_sharur_a, **_sharur_k)

            query = sql

            def executemany(self, _sharur_q, *_sharur_a, **_sharur_k):
                return self._wrap(self._con.executemany(_rw(_sharur_q), *_sharur_a, **_sharur_k))

            def __getattr__(self, name):
                return getattr(self._con, name)

            def __setattr__(self, name, value):
                setattr(self._con, name, value)

            def __enter__(self):
                self._con.__enter__()
                return self

            def __exit__(self, *exc):
                return self._con.__exit__(*exc)

        def connect(database=":memory:", read_only=False, config=None, **kw):
            db, is_file = _mapped(database)
            cfg = dict(config or {})
            cfg.pop("access_mode", None)
            con = _real(db, read_only=True if is_file else read_only, config=cfg, **kw)
            con.execute("SET enable_progress_bar = false")
            con.execute(f"SET threads = {int(_cfg.get('threads', 2))}")
            con.execute(f"SET memory_limit = '{_cfg.get('memory_limit', '4GB')}'")
            if _cfg.get("temp_directory"):
                con.execute(f"SET temp_directory = '{_cfg['temp_directory']}'")
            if not _variants:
                return con
            try:
                con.execute("SET python_scan_all_frames = true")
            except Exception:  # noqa: BLE001 - a DuckDB without the setting scans the calling frame
                pass
            return _Conn(con)

        connect._sharur_reverify = True
        _duckdb.connect = connect
'''


def _shim_config(sub: Substitution, threads: int, memory_limit: str, temp_dir: str,
                 guard: Sequence[str] = ()) -> str:
    cfg = json.loads(sub.env())
    cfg.update({"threads": threads, "memory_limit": memory_limit, "temp_directory": temp_dir,
                "guard": list(guard)})
    return json.dumps(cfg)


def _sbpl(path: str) -> str:
    return path.replace("\\", "\\\\").replace('"', '\\"')


def sandbox_prefix(paths: Sequence[str]) -> list[str] | None:
    """``sandbox-exec`` argv that denies file writes under ``paths`` (macOS), or ``None``.

    Probed once: where the profile cannot be applied (no ``sandbox-exec``, or a parent sandbox
    forbids it), shell triples run without it and the Python audit hook alone guards writes.
    """
    exe = shutil.which("sandbox-exec") if sys.platform == "darwin" else None
    if not exe:
        return None
    rules = "".join(f'(deny file-write* (subpath "{_sbpl(os.path.realpath(p))}"))' for p in paths)
    argv = [exe, "-p", "(version 1)(allow default)" + rules]
    try:
        ok = subprocess.run([*argv, "/usr/bin/true"], capture_output=True, timeout=30).returncode == 0
    except (OSError, subprocess.SubprocessError):
        ok = False
    return argv if ok else None


# --------------------------------------------------------------------------- #
# Pool sizing
# --------------------------------------------------------------------------- #
# Resident memory a worker holds beyond its DuckDB limit: the interpreter, pandas and numpy,
# fetched results, and helper modules a replay imports. Workers with a 4 GB limit reached
# 1.7-5.5 GB resident replaying a set of defense findings (2026-11-02).
WORKER_OVERHEAD_BYTES = 3 << 29  # 1.5 GiB
# Share of currently available memory the pool may commit.
MEMORY_FRACTION = 0.5

_UNITS = {"": 1, "B": 1, "K": 1 << 10, "KB": 1 << 10, "KIB": 1 << 10, "M": 1 << 20, "MB": 1 << 20,
          "MIB": 1 << 20, "G": 1 << 30, "GB": 1 << 30, "GIB": 1 << 30, "T": 1 << 40, "TB": 1 << 40,
          "TIB": 1 << 40}


def memory_bytes(text: str) -> int:
    """Bytes in a DuckDB-style size (``4GB``, ``512MiB``, ``1.5 G``); units count in 1024s."""
    m = re.fullmatch(r"\s*([\d.]+)\s*([A-Za-z]*)\s*", str(text))
    if not m or m.group(2).upper() not in _UNITS:
        raise ValueError(f"memory size not understood: {text!r}")
    return int(float(m.group(1)) * _UNITS[m.group(2).upper()])


def available_memory() -> int | None:
    """Bytes this process can still commit: free memory, capped by a Slurm allocation."""
    avail: int | None = None
    try:
        import psutil  # noqa: PLC0415

        avail = int(psutil.virtual_memory().available)
    except ImportError:
        try:
            meminfo = Path("/proc/meminfo").read_text()
            m = re.search(r"^MemAvailable:\s+(\d+) kB", meminfo, re.M)
            avail = int(m.group(1)) * 1024 if m else None
        except OSError:
            try:  # macOS: free, inactive, speculative and purgeable pages can be reclaimed
                out = subprocess.run(["vm_stat"], capture_output=True, text=True, timeout=5, check=False).stdout
                page = int(re.search(r"page size of (\d+) bytes", out).group(1))
                pages = sum(int(n) for n in re.findall(
                    r"Pages (?:free|inactive|speculative|purgeable):\s+(\d+)", out))
                avail = pages * page
            except (OSError, AttributeError, ValueError, subprocess.SubprocessError):
                avail = None
    slurm_mb = os.environ.get("SLURM_MEM_PER_NODE")
    if not slurm_mb and os.environ.get("SLURM_MEM_PER_CPU"):
        cpus = os.environ.get("SLURM_CPUS_ON_NODE") or os.environ.get("SLURM_CPUS_PER_TASK") or "1"
        slurm_mb = str(int(os.environ["SLURM_MEM_PER_CPU"]) * int(cpus))
    if slurm_mb and slurm_mb.isdigit():
        cap = int(slurm_mb) << 20
        avail = cap if avail is None else min(avail, cap)
    return avail


def plan_workers(requested: int | None, threads: int, memory_limit: str, n_jobs: int,
                 available: int | None = None) -> tuple[int, str | None]:
    """Worker count for a replay pool, and a note when free memory lowered it.

    The default is one worker per ``threads`` CPUs. Any count, requested or default, is
    capped so that workers * (``memory_limit`` + ``WORKER_OVERHEAD_BYTES``) stays within
    ``MEMORY_FRACTION`` of available memory, with one worker always allowed.
    """
    cpu_default = max(1, (os.cpu_count() or 2) // max(1, threads))
    want = max(1, min(requested or cpu_default, n_jobs or 1))
    avail = available_memory() if available is None else available
    if avail is None:
        return want, None
    per_worker = memory_bytes(memory_limit) + WORKER_OVERHEAD_BYTES
    budget = int(avail * MEMORY_FRACTION)
    fits = max(1, budget // per_worker)
    if want <= fits:
        return want, None
    gib = 1 << 30
    return fits, (f"{want} workers lowered to {fits}: {avail / gib:.1f} GiB available, "
                  f"{MEMORY_FRACTION:.0%} budget at {per_worker / gib:.1f} GiB per worker "
                  f"(--memory-limit {memory_limit} + {WORKER_OVERHEAD_BYTES / gib:.1f} GiB)")


# --------------------------------------------------------------------------- #
# Worker
# --------------------------------------------------------------------------- #
_WORKER: dict[str, Any] = {}
# Seconds between a worker's checks that the process that started it is still alive.
PARENT_POLL_SECONDS = 1.0


def _kill_shell_group() -> None:
    """SIGKILL the process group of the shell triple this worker is running, if any.

    The group outlives /bin/sh while anything the command started is still in it, so the
    kill goes to the group whether or not the shell itself has exited.
    """
    proc = _WORKER.get("shell")
    if proc is not None:
        with contextlib.suppress(ProcessLookupError, PermissionError):
            os.killpg(proc.pid, signal.SIGKILL)


def _exit_worker(*_: Any) -> None:
    _kill_shell_group()
    os._exit(1)


def _watch_parent(parent_pid: int) -> None:
    """Exit when the reverify process goes away; a worker reparented to init is an orphan."""
    while True:
        time.sleep(PARENT_POLL_SECONDS)
        if os.getppid() != parent_pid:
            _exit_worker()


def _worker_init(root: str, shim_dir: str, threads: int, memory_limit: str, temp_dir: str,
                 parent_pid: int | None = None, sandbox: list[str] | None = None) -> None:
    _WORKER.update(root=root, shim_dir=shim_dir, threads=threads,
                   memory_limit=memory_limit, temp_dir=temp_dir, sandbox=sandbox)
    # A worker lives only as long as the reverify process: SIGTERM from the pool, or the
    # parent's death by any signal, ends the worker and the shell triple it is running.
    signal.signal(signal.SIGTERM, _exit_worker)
    if parent_pid is not None:
        threading.Thread(target=_watch_parent, args=(parent_pid,), name="reverify-parent-watch",
                         daemon=True).start()
    # Replayed code and DuckDB's native progress bar write to fds 1 and 2; the parent's
    # stdout carries the report, so workers write to /dev/null.
    devnull = os.open(os.devnull, os.O_WRONLY)
    os.dup2(devnull, 1)
    os.dup2(devnull, 2)
    # The shim patches duckdb.connect for this process too; the configuration
    # is swapped per task, so it reads the environment at call time.
    os.environ.setdefault("SHARUR_REVERIFY_SHIM", json.dumps({"root": root, "map": {}, "variants": {}}))
    sys.path.insert(0, shim_dir)
    import importlib  # noqa: PLC0415

    import duckdb  # noqa: PLC0415

    _WORKER["real_connect"] = duckdb.connect
    importlib.import_module("sitecustomize")


def _install(shim_cfg: str) -> None:
    """Re-patch duckdb.connect with this task's substitution."""
    import importlib  # noqa: PLC0415

    import duckdb  # noqa: PLC0415

    os.environ["SHARUR_REVERIFY_SHIM"] = shim_cfg
    duckdb.connect = _WORKER["real_connect"]
    sys.modules.pop("sitecustomize", None)
    importlib.import_module("sitecustomize")


def normalize_rows(rows: Any) -> Any:
    if isinstance(rows, list):
        if len(rows) == 1 and len(rows[0]) == 1:
            return rows[0][0]
        if rows and all(len(r) == 1 for r in rows):
            return [r[0] for r in rows]
        return [list(r) for r in rows]
    return rows


def to_plain(value: Any) -> Any:
    """JSON-comparable Python values (numpy scalars, Decimals, tuples, frames)."""
    if hasattr(value, "item") and not isinstance(value, (list, dict, str, bytes)):
        try:
            return to_plain(value.item())
        except (ValueError, TypeError, AttributeError):
            pass
    if hasattr(value, "tolist") and not isinstance(value, (list, dict, str, bytes)):
        try:
            return to_plain(value.tolist())
        except (ValueError, TypeError, AttributeError):
            pass
    if isinstance(value, dict):
        return {str(k): to_plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        seq = sorted(value, key=repr) if isinstance(value, (set, frozenset)) else value
        return [to_plain(v) for v in seq]
    if isinstance(value, decimal.Decimal):
        return float(value)
    if hasattr(value, "to_json") and callable(value.to_json):
        try:
            return json.loads(value.to_json())
        except (TypeError, ValueError):
            pass
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def _run_one(triple: Triple, shim_cfg: str, timeout: float) -> Outcome:
    root = Path(_WORKER["root"])
    sub_cfg = json.loads(shim_cfg)
    sub = Substitution([], root)
    if sub_cfg.get("map"):
        sub = Substitution([(k, v) for k, v in sub_cfg["map"].items()], root)
    _install(shim_cfg)
    start = time.time()
    query = sub.text(triple.command or triple.query.strip())
    saved_path = list(sys.path)
    saved_modules = set(sys.modules)
    cwd = os.getcwd()
    try:
        os.chdir(root)
        if triple.kind == "shell":
            command = query[len("shell:"):] if query.startswith("shell:") else query
            written = _writes_inside(command, root, sub_cfg.get("guard") or [str(root)])
            if written:
                raise PermissionError(f"refused: the command writes {written}; replays are read-only")
            env = dict(os.environ)
            env["PYTHONPATH"] = os.pathsep.join(filter(None, [_WORKER["shim_dir"], env.get("PYTHONPATH")]))
            env["PYTHONDONTWRITEBYTECODE"] = "1"
            # Its own session, so a timeout or the worker's exit kills every process the
            # command starts, not only /bin/sh.
            proc = subprocess.Popen([*(_WORKER.get("sandbox") or []), "/bin/sh", "-c", command], cwd=root,
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env,
                                    start_new_session=True)
            _WORKER["shell"] = proc
            try:
                stdout, stderr = proc.communicate(timeout=timeout)
            except subprocess.TimeoutExpired:
                _kill_shell_group()
                try:  # a descendant that left the group could hold the pipes open
                    proc.communicate(timeout=10)
                except subprocess.TimeoutExpired:
                    proc.kill()
                raise
            finally:
                _kill_shell_group()  # background children the command left behind
                _WORKER.pop("shell", None)
            out = stdout.strip()
            if proc.returncode != 0 and not out:
                raise RuntimeError(f"exit {proc.returncode}: {stderr.strip()[-300:]}")
            observed = read_output(out, triple.expected)
            if triple.select is not None:
                observed = select_key(observed, triple.select)
        elif triple.kind == "python":
            import math  # noqa: PLC0415

            import duckdb  # noqa: PLC0415
            import numpy as np  # noqa: PLC0415
            import pandas as pd  # noqa: PLC0415

            from sharur.operators import Sharur  # noqa: PLC0415

            expr = query[len("python:"):] if query.startswith("python:") else query
            scope = {"json": json, "pd": pd, "np": np, "duckdb": duckdb, "Sharur": Sharur,
                     "os": os, "re": re, "math": math, "Path": Path, "sys": sys}
            try:
                code = compile(expr.strip(), "<verification>", "eval")
                observed = eval(code, scope)  # noqa: S307 - replaying recorded verification
            except SyntaxError:
                exec(compile(expr.strip(), "<verification>", "exec"), scope)  # noqa: S102
                observed = scope.get("result", scope.get("_"))
        else:
            import duckdb  # noqa: PLC0415

            if triple.attach:
                con = duckdb.connect(":memory:")
                for alias, target in triple.attach.items():
                    target = sub.path(target)
                    target_abs = target if Path(target).is_absolute() else str(root / target)
                    con.execute(f"ATTACH '{_sql_quote(target_abs)}' AS {alias} (READ_ONLY)")
            else:
                if not triple.default_db:
                    raise RuntimeError("SQL without attach and no default database (pass --db)")
                con = duckdb.connect(sub.path(triple.default_db), read_only=True)
                _install_sql_prelude(con)  # the SQL surface (neighbors, slot, ...), TEMP macros
            try:
                observed = normalize_rows(con.execute(query).fetchall())
            finally:
                con.close()
        return Outcome(match(to_plain(observed), triple.expected), to_plain(observed),
                       seconds=time.time() - start)
    except subprocess.TimeoutExpired:
        return Outcome("ERROR", error=f"timeout after {timeout:.0f}s", seconds=time.time() - start)
    except (Exception, SystemExit) as exc:  # noqa: BLE001 - every failure is a reportable ERROR
        return Outcome("ERROR", error=f"{type(exc).__name__}: {exc}"[:500], seconds=time.time() - start)
    finally:
        os.chdir(cwd)
        sys.path[:] = saved_path
        # helper modules imported by this triple may cache connections or results
        for name in set(sys.modules) - saved_modules:
            module = sys.modules.get(name)
            path = getattr(module, "__file__", None) or ""
            if path and str(root) in path and "site-packages" not in path:
                sys.modules.pop(name, None)


def _sql_quote(value: str) -> str:
    return value.replace("'", "''")


def _writes_inside(command: str, root: Path, guard: Sequence[str]) -> str | None:
    """The first file inside ``guard`` that ``command`` redirects output to, if any."""
    try:
        targets = _redirect_targets(_shell_tokens(command))
    except ValueError:
        return None
    roots = [os.path.realpath(g) for g in guard]
    for target in targets:
        path = os.path.realpath(target if os.path.isabs(target) else os.path.join(root, target))
        if any(path == g or path.startswith(g + os.sep) for g in roots):
            return target
    return None


_NP_SCALAR = re.compile(r"\bnp\.(?:float|int|uint|bool_?|str_)\d*\(([^()]*)\)")


def read_output(text: str, expected: Any) -> Any:
    """A command's stdout as a value.

    JSON when it parses. Otherwise, when ``expected`` is not a string (so the raw text could
    never match it): a Python literal such as a printed dict or list, with numpy scalar wrappers
    (``np.float64(0.5)``) read as their numbers, or two or more whitespace-separated numbers as
    a list. Anything else stays text.
    """
    try:
        return json.loads(text)
    except (json.JSONDecodeError, ValueError):
        pass
    if isinstance(expected, str):
        return text
    with contextlib.suppress(ValueError, SyntaxError, TypeError, MemoryError, RecursionError):
        return ast.literal_eval(_NP_SCALAR.sub(r"\1", text))
    numbers = []
    for token in text.split():
        try:
            value = json.loads(token)
        except (json.JSONDecodeError, ValueError):
            return text
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            return text
        numbers.append(value)
    return numbers if len(numbers) > 1 else text


def select_key(observed: Any, key: str) -> Any:
    """``observed[key]`` for a command annotated ``(prints JSON; key 'K')``."""
    if not isinstance(observed, dict):
        raise TypeError(f"output is not a JSON object, so key {key!r} cannot be read")
    if key not in observed:
        raise KeyError(f"output has no key {key!r}")
    return observed[key]


def _install_sql_prelude(con: Any) -> None:
    """The Sharur SQL surface (TEMP macros, no FASTA function) on a replay connection."""
    from sharur.sql import install_prelude  # noqa: PLC0415

    with contextlib.suppress(Exception):
        install_prelude(con, udfs=False)


def _statistical_flags(triple: Triple, root: Path) -> list[str]:
    """Input, seal and pre-registration flags of a statistical check (none for count checks)."""
    if triple.stat is None:
        return []
    from sharur.stat_checks import static_flags  # noqa: PLC0415

    return static_flags(triple, root)


def _statistical_rejudge(result: TripleResult, outcome: Outcome, root: Path) -> Outcome:
    """Tolerance and verdict rules of a statistical check; count checks pass through unchanged."""
    if result.triple.tolerance is None and result.triple.stat is None:
        return outcome
    from sharur.stat_checks import rejudge  # noqa: PLC0415

    outcome, flags = rejudge(result.triple, outcome, root)
    result.flags.extend(f for f in flags if f not in result.flags)
    return outcome


# --------------------------------------------------------------------------- #
# Comparison
# --------------------------------------------------------------------------- #
def _quantum(expected: float) -> Decimal:
    d = Decimal(repr(float(expected)))
    return Decimal(1).scaleb(d.as_tuple().exponent)


def _float_eq(observed: float, expected: float) -> bool:
    q = _quantum(expected)
    o = Decimal(repr(float(observed)))
    e = Decimal(repr(float(expected))).quantize(q)
    return o.quantize(q, rounding=ROUND_HALF_UP) == e or o.quantize(q, rounding=ROUND_HALF_EVEN) == e


def match(observed: Any, expected: Any) -> str:
    """``OK`` when ``observed`` reproduces ``expected`` under the findings rules."""
    if isinstance(expected, bool) or isinstance(observed, bool):
        same = isinstance(expected, bool) and isinstance(observed, bool) and observed == expected
        return "OK" if same else "FAIL"
    if isinstance(expected, int) and isinstance(observed, (int, float)):
        return "OK" if observed == expected and float(observed).is_integer() else "FAIL"
    if isinstance(expected, float) and isinstance(observed, (int, float)):
        try:
            return "OK" if _float_eq(observed, expected) else "FAIL"
        except (decimal.InvalidOperation, ValueError, OverflowError):
            return "OK" if observed == expected else "FAIL"
    if isinstance(expected, dict) and isinstance(observed, dict):
        if set(map(str, expected)) != set(map(str, observed)):
            return "FAIL"
        return "FAIL" if any(match(observed[str(k)], v) == "FAIL" for k, v in expected.items()) else "OK"
    if isinstance(expected, list) and isinstance(observed, list):
        if len(expected) != len(observed):
            return "FAIL"
        return "FAIL" if any(match(o, e) == "FAIL" for o, e in zip(observed, expected)) else "OK"
    if isinstance(expected, list) and not isinstance(observed, list) and len(expected) == 1:
        return match(observed, expected[0])
    if isinstance(expected, str) and not isinstance(observed, str):
        return "OK" if str(observed) == expected else "FAIL"
    if isinstance(expected, (int, float)) and isinstance(observed, str):
        try:
            return match(json.loads(observed), expected)
        except (json.JSONDecodeError, ValueError):
            return "FAIL"
    return "OK" if observed == expected else "FAIL"


def _canonical(value: Any) -> str:
    return json.dumps(value, sort_keys=True, default=str)


def flags_for(triple: Triple) -> list[str]:
    """Durability notes for a triple (scratch paths, empty expectations)."""
    out = []
    e, q = triple.expected, triple.query
    if e in ("", None, []) or (isinstance(e, list) and any(x in ("", None) for x in e)):
        out.append("EMPTY_EXPECTED")
    if re.search(r"/private/tmp|scratchpad|(?<![\w/])/tmp/", q):
        out.append("SCRATCH_PATH")
    if re.search(r"grep\s+(-\w*P|--perl)", q):
        out.append("GREP_P")
    if triple.kind == "sql" and not triple.attach and not triple.default_db:
        out.append("NO_DEFAULT_DB")
    return out


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
@dataclass
class ReverifyReport:
    root: str
    findings: list[str]
    substitutions: dict[str, str]
    baseline_substitutions: dict[str, str]
    compare: bool
    results: list[TripleResult]
    selected_findings: int
    seconds: float
    citations: dict[str, list[str]] = field(default_factory=dict)
    workers: int = 0
    worker_note: str | None = None

    def counts(self, which: str = "replay") -> dict[str, int]:
        out = dict.fromkeys(STATUSES, 0)
        for r in self.results:
            outcome = r.replay if which == "replay" else r.baseline
            if outcome is not None:
                out[outcome.status] += 1
        return out

    def changed(self) -> list[TripleResult]:
        return [r for r in self.results if r.changed]

    def attention(self) -> list[TripleResult]:
        if self.compare:
            return self.changed()
        return [r for r in self.results if r.replay.status not in ("OK", "NARRATIVE")]

    def to_dict(self) -> dict[str, Any]:
        def row(r: TripleResult) -> dict[str, Any]:
            d = {"key": r.triple.key, "finding_id": r.triple.finding_id, "title": r.triple.title,
                 "claim": r.triple.claim, "kind": r.triple.kind, "expected": r.triple.expected,
                 "flags": r.flags, "replay": asdict(r.replay)}
            if r.baseline is not None:
                d["baseline"] = asdict(r.baseline)
                d["changed"] = r.changed
            return d

        return {
            "root": self.root, "findings": self.findings, "substitutions": self.substitutions,
            "baseline_substitutions": self.baseline_substitutions,
            "compare": self.compare, "selected_findings": self.selected_findings,
            "triples": len(self.results), "seconds": round(self.seconds, 1),
            "workers": self.workers, "worker_note": self.worker_note,
            "replay_counts": self.counts("replay"),
            "baseline_counts": self.counts("baseline") if self.compare else None,
            "changed": len(self.changed()) if self.compare else None,
            "citations": self.citations,
            "results": [row(r) for r in self.results],
        }

    def to_markdown(self, limit: int = 200) -> str:
        lines = [f"# Reverify: {len(self.results)} triples from {self.selected_findings} findings "
                 f"({self.seconds:.0f} s, {self.workers} worker{'' if self.workers == 1 else 's'})", ""]
        if self.worker_note:
            lines.append(f"Workers: {self.worker_note}.")
        if self.baseline_substitutions:
            lines.append("Baseline substitutions: " + "; ".join(
                f"`{k}` → `{v}`" for k, v in self.baseline_substitutions.items()))
        if self.substitutions:
            lines.append("Substitutions: " + "; ".join(f"`{k}` → `{v}`" for k, v in self.substitutions.items()))
        def tally(c: dict[str, int]) -> str:
            return f"{c['OK']} OK, {c['FAIL']} FAIL, {c['ERROR']} ERROR, {c['NARRATIVE']} NARRATIVE"

        if self.compare:
            lines.append(f"Baseline: {tally(self.counts('baseline'))}. Replay: {tally(self.counts('replay'))}. "
                         f"**{len(self.changed())} triples change.**")
        else:
            lines.append(f"{tally(self.counts())}.")
        if self.counts()["NARRATIVE"]:
            lines.append("NARRATIVE triples record a note (a methods reference, or a command annotated "
                         "with the file its value was read from) and are reported without being run.")
        lines.append("")
        by_finding: dict[tuple[str, str], list[TripleResult]] = {}
        for r in self.attention():
            by_finding.setdefault((r.triple.source, r.triple.finding_id), []).append(r)
        shown = 0
        for (source, fid), items in by_finding.items():
            if shown >= limit:
                lines.append(f"... {len(by_finding) - shown} more findings (see --out)")
                break
            shown += 1
            lines.append(f"## {fid} ({source}:{items[0].triple.line}) {items[0].triple.title[:120]}")
            for cite in self.citations.get(fid, [])[:8]:
                lines.append(f"- cited in `{cite}`")
            for r in items:
                exp = _short(r.triple.expected)
                if self.compare:
                    before = _outcome_text(r.baseline)
                    after = _outcome_text(r.replay)
                    lines.append(f"- #{r.triple.index} {r.triple.claim[:160]} | expected {exp} | "
                                 f"{before} → {after}")
                else:
                    lines.append(f"- #{r.triple.index} {r.replay.status} {r.triple.claim[:160]} | "
                                 f"expected {exp} | {_outcome_text(r.replay)}")
                if r.triple.stat is not None and r.flags:
                    lines.append(f"  - statistical check: {'; '.join(r.flags)}")
            lines.append("")
        return "\n".join(lines)


def _short(value: Any, n: int = 120) -> str:
    text = json.dumps(value, default=str)
    return text if len(text) <= n else text[: n - 1] + "…"


def _outcome_text(outcome: Outcome | None) -> str:
    if outcome is None:
        return "-"
    if outcome.status in ("ERROR", "NARRATIVE"):
        return f"{outcome.status} {outcome.error or ''}"[:200]
    return f"{outcome.status} {_short(outcome.observed)}"


def write_tsv(report: ReverifyReport, path: Path) -> None:
    cols = ["key", "finding_id", "triple", "kind", "status", "baseline_status", "changed", "flags",
            "seconds", "claim", "expected", "observed", "baseline_observed", "error", "baseline_error"]
    with path.open("w") as handle:
        handle.write("\t".join(cols) + "\n")
        for r in report.results:
            rec = r.baseline
            row = [r.triple.key, r.triple.finding_id, r.triple.index, r.triple.kind, r.replay.status,
                   rec.status if rec else "", r.changed if rec else "", ",".join(r.flags),
                   f"{r.replay.seconds:.1f}", r.triple.claim, _short(r.triple.expected, 400),
                   _short(r.replay.observed, 400), _short(rec.observed, 400) if rec else "",
                   (r.replay.error or ""), (rec.error or "") if rec else ""]
            handle.write("\t".join(str(x).replace("\t", " ").replace("\n", " ") for x in row) + "\n")


def find_citations(finding_ids: Iterable[str], globs: Sequence[str], root: Path) -> dict[str, list[str]]:
    """Report lines that mention each finding ID (``file:line``)."""
    ids = sorted({i for i in finding_ids if i}, key=len, reverse=True)
    if not ids or not globs:
        return {}
    pattern = re.compile(r"(?<![\w-])(" + "|".join(re.escape(i) for i in ids) + r")(?![\w-])")
    out: dict[str, list[str]] = {}
    for g in globs:
        for path in sorted(root.glob(g)):
            if not path.is_file():
                continue
            try:
                text = path.read_text(errors="replace")
            except OSError:
                continue
            rel = path.relative_to(root) if path.is_relative_to(root) else path
            for n, line in enumerate(text.splitlines(), 1):
                for m in pattern.finditer(line):
                    out.setdefault(m.group(1), []).append(f"{rel}:{n}")
    return out


def _terminate_pool(pool: cf.ProcessPoolExecutor) -> None:
    """Stop a replay pool now: drop queued replays, SIGTERM every worker, reap them."""
    procs = list((getattr(pool, "_processes", None) or {}).values())
    pool.shutdown(wait=False, cancel_futures=True)
    for proc in procs:
        if proc.is_alive():
            proc.terminate()  # the worker's SIGTERM handler kills its shell triple's group
    deadline = time.time() + 10
    for proc in procs:
        proc.join(max(0.0, deadline - time.time()))
        if proc.is_alive():
            proc.kill()
            proc.join(5)


def reverify(
    findings: Sequence[Path],
    *,
    root: Path | None = None,
    db: Path | None = None,
    touching: str | None = None,
    ids: Sequence[str] | None = None,
    substitutions: Sequence[tuple[str, str]] = (),
    baseline_substitutions: Sequence[tuple[str, str]] = (),
    compare: bool = False,
    workers: int | None = None,
    threads: int = 2,
    memory_limit: str = "4GB",
    timeout: float = 1800.0,
    report_globs: Sequence[str] = (),
    progress: Any = None,
) -> ReverifyReport:
    """Replay verification triples; see the module docstring."""
    files = resolve_findings(findings)
    root = (root or find_root(files[0])).resolve()
    triples = load_triples(files, root, db)
    if ids:
        wanted = set(ids)
        triples = [t for t in triples if t.finding_id in wanted]
    if touching:
        rx = re.compile(touching, re.I)
        cache: dict[Path, str] = {}
        hit_findings = {(t.source, t.finding_id) for t in triples if rx.search(dependency_text(t, root, cache))}
        # a finding that touches the table is replayed whole: its other numbers often derive from it
        triples = [t for t in triples if (t.source, t.finding_id) in hit_findings]
    if compare and not (substitutions or baseline_substitutions):
        raise ValueError("compare mode needs a substitution (OLD=NEW) for the replay or the baseline")
    sub = Substitution(substitutions, root)
    base = Substitution(baseline_substitutions, root)
    for old, new in [*sub.map.items(), *base.map.items()]:
        if not Path(new).exists():
            raise FileNotFoundError(f"substitution target does not exist: {new} (for {old})")

    tmp = Path(tempfile.mkdtemp(prefix="sharur-reverify-"))
    shim_dir = tmp / "shim"
    shim_dir.mkdir()
    (shim_dir / "sitecustomize.py").write_text(_SHIM)
    spill = tmp / "spill"
    spill.mkdir()
    # nothing a replay runs may write inside the repository or to a substituted database
    guard = [str(root), *sub.map.values(), *base.map.values()]
    plain_cfg = _shim_config(base, threads, memory_limit, str(spill), guard)
    sub_cfg = _shim_config(sub, threads, memory_limit, str(spill), guard)
    sandbox = sandbox_prefix(guard) if any(t.kind == "shell" for t in triples) else None

    results = [TripleResult(t, flags_for(t) + _statistical_flags(t, root)) for t in triples]
    jobs: list[tuple[int, str, str]] = []
    for i, result in enumerate(results):
        if result.triple.kind == "narrative":  # a note: reported, never run
            result.replay = Outcome("NARRATIVE", error=result.triple.note)
            if compare:
                result.baseline = Outcome("NARRATIVE", error=result.triple.note)
            continue
        if compare:
            jobs.append((i, "baseline", plain_cfg))
        jobs.append((i, "replay", sub_cfg))
    n_workers, worker_note = plan_workers(workers, threads, memory_limit, len(jobs))
    if worker_note:
        logger.warning("reverify: %s", worker_note)
    started = time.time()
    # SIGTERM (pkill, a scheduler) becomes KeyboardInterrupt so the pool is torn down below.
    previous_term: Any = None
    term_installed = threading.current_thread() is threading.main_thread()
    if term_installed:
        def _on_term(signum: int, frame: Any) -> None:
            raise KeyboardInterrupt(f"signal {signum}")

        previous_term = signal.signal(signal.SIGTERM, _on_term)
    pool: cf.ProcessPoolExecutor | None = None
    try:
        ctx = multiprocessing.get_context("spawn")
        # No max_tasks_per_child: CPython 3.12's pool stops replacing retired workers and
        # hangs. Each task drops the helper modules it imported, which gives the same isolation.
        pool = cf.ProcessPoolExecutor(
            max_workers=n_workers, mp_context=ctx, initializer=_worker_init,
            initargs=(str(root), str(shim_dir), threads, memory_limit, str(spill), os.getpid(), sandbox))
        futures = {pool.submit(_run_one, results[i].triple, cfg, timeout): (i, slot)
                   for i, slot, cfg in jobs}
        done = 0
        for fut in cf.as_completed(futures):
            i, slot = futures[fut]
            try:
                outcome = fut.result()
            except Exception as exc:  # noqa: BLE001 - a crashed worker is an ERROR
                outcome = Outcome("ERROR", error=f"worker: {type(exc).__name__}: {exc}"[:500])
            outcome = _statistical_rejudge(results[i], outcome, root)
            setattr(results[i], slot, outcome)
            done += 1
            if progress is not None:
                progress(done, len(jobs), results[i], slot)
        pool.shutdown(wait=True)
    except BaseException:
        if pool is not None:
            _terminate_pool(pool)
        raise
    finally:
        if term_installed:
            signal.signal(signal.SIGTERM, previous_term if previous_term is not None else signal.SIG_DFL)
        shutil.rmtree(tmp, ignore_errors=True)
    report = ReverifyReport(
        root=str(root), findings=[str(f) for f in files], substitutions=dict(sub.map),
        baseline_substitutions=dict(base.map),
        compare=compare, results=results,
        selected_findings=len({(t.source, t.finding_id) for t in triples}),
        seconds=time.time() - started, workers=n_workers, worker_note=worker_note,
    )
    if report_globs:
        report.citations = find_citations([r.triple.finding_id for r in report.attention()], report_globs, root)
    return report
