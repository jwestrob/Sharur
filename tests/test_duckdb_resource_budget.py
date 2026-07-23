"""Per-agent DuckDB resource budget contracts."""

from __future__ import annotations

import pytest

from sharur.core.session import ExplorationSession
from sharur.operators import Sharur
from sharur.storage.duckdb_store import DuckDBStore


def test_duckdb_store_applies_explicit_thread_memory_and_spill_limits(tmp_path):
    spill = tmp_path / "spill"
    store = DuckDBStore(
        threads=2,
        memory_limit="256MB",
        temp_directory=spill,
    )

    threads, memory, temp_directory = store.conn.execute(
        """
        SELECT current_setting('threads'),
               current_setting('memory_limit'),
               current_setting('temp_directory')
        """
    ).fetchone()
    assert threads == 2
    assert "MiB" in memory
    assert temp_directory == str(spill.resolve())
    assert store.resource_budget == {
        "threads": 2,
        "memory_limit": "256MB",
        "temp_directory": str(spill.resolve()),
    }


def test_duckdb_budget_propagates_through_session_and_facade(tmp_path):
    session = ExplorationSession(
        duckdb_threads=1,
        duckdb_memory_limit="128MB",
        duckdb_temp_directory=tmp_path / "session-spill",
    )
    assert session.db.resource_budget["threads"] == 1

    facade = Sharur(
        duckdb_threads=2,
        duckdb_memory_limit="256MB",
        duckdb_temp_directory=tmp_path / "facade-spill",
    )
    assert facade.store.resource_budget["threads"] == 2
    assert facade.store.conn.execute(
        "SELECT current_setting('threads')"
    ).fetchone()[0] == 2


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"threads": 0}, "positive integer"),
        ({"memory_limit": "lots"}, "supported unit"),
    ],
)
def test_duckdb_budget_rejects_ambiguous_limits(kwargs, message):
    with pytest.raises(ValueError, match=message):
        DuckDBStore(**kwargs)
