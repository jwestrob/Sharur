"""DuckDB 1.2.1 returns wrong multithreaded window results; Sharur detects and rejects it."""

import warnings

import duckdb
import pytest

from sharur import diagnostics
from sharur.storage import duckdb_store
from sharur.storage.duckdb_store import KNOWN_BAD_DUCKDB, DuckDBStore, duckdb_version_problem

BAD = duckdb.__version__ in KNOWN_BAD_DUCKDB

# An ordered window plus a whole-partition SUM/MAX over the same partition key: the shape that
# DuckDB 1.2.1 answers wrongly (NULL or wrong values) when run on more than one thread.
QUERY = """
SELECT hash(list(row(id, running, total, biggest) ORDER BY id)) FROM (
    SELECT id,
           SUM(v) OVER (PARTITION BY g ORDER BY v DESC, id) AS running,
           SUM(v) OVER (PARTITION BY g) AS total,
           MAX(v) OVER (PARTITION BY g) AS biggest
    FROM t)
"""


@pytest.mark.xfail(BAD, reason=f"DuckDB {duckdb.__version__} is a known-bad release (see KNOWN_BAD_DUCKDB)",
                   strict=False)
def test_multithreaded_window_results_match_single_thread():
    conn = duckdb.connect()
    conn.execute("""CREATE TABLE t AS SELECT (i % 2000)::VARCHAR AS g, (hash(i) % 100000)::INTEGER AS v, i AS id
                    FROM range(400000) r(i)""")
    conn.execute("SET threads = 1")
    truth = conn.execute(QUERY).fetchone()[0]
    conn.execute("SET threads = 8")
    assert [conn.execute(QUERY).fetchone()[0] for _ in range(10)] == [truth] * 10


def test_known_bad_release_is_reported():
    assert duckdb_version_problem("1.2.1").startswith("DuckDB 1.2.1:")
    assert "pip install -U duckdb" in duckdb_version_problem("1.2.1")
    assert duckdb_version_problem("1.2.2") is None
    assert duckdb_version_problem("1.1.3") is None


def test_doctor_fails_the_core_check_on_a_bad_release(monkeypatch):
    monkeypatch.setattr(duckdb, "__version__", "1.2.1")
    check = diagnostics.check_duckdb()
    assert (check.status, check.core) == (diagnostics.MISSING, True)
    assert diagnostics.has_core_failure([check])
    monkeypatch.setattr(duckdb, "__version__", "1.2.2")
    assert diagnostics.check_duckdb().status == diagnostics.OK


def test_store_warns_once_on_a_bad_release(monkeypatch):
    monkeypatch.setattr(duckdb, "__version__", "1.2.1")
    monkeypatch.setattr(duckdb_store, "_warned_bad_duckdb", False)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        DuckDBStore().conn
        DuckDBStore().conn
    messages = [str(w.message) for w in caught if issubclass(w.category, RuntimeWarning)]
    assert len(messages) == 1 and "DuckDB 1.2.1" in messages[0]


def test_dependency_pin_excludes_bad_releases():
    from pathlib import Path

    text = (Path(__file__).resolve().parents[1] / "pyproject.toml").read_text()
    pin = next(line for line in text.splitlines() if line.strip().startswith('"duckdb'))
    assert all(f"!={version}" in pin for version in KNOWN_BAD_DUCKDB)
