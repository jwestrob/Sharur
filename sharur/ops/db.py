"""SQLite connection ownership and pooling for the Sharur Ops service."""

from __future__ import annotations

import contextlib
import fcntl
import os
import queue
import socket
import sqlite3
import threading
import time
from pathlib import Path
from typing import TYPE_CHECKING

from sharur.ops.schema import ensure_ops_schema


if TYPE_CHECKING:
    from collections.abc import Iterator


DEFAULT_BUSY_TIMEOUT_MS = 15_000


def open_ops_connection(
    db_path: str | Path,
    *,
    initialize: bool = False,
    busy_timeout_ms: int = DEFAULT_BUSY_TIMEOUT_MS,
) -> sqlite3.Connection:
    """Open one consistently configured coordination connection."""

    path = Path(db_path).expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(
        str(path),
        check_same_thread=False,
        isolation_level="DEFERRED",
    )
    conn.row_factory = sqlite3.Row
    conn.execute(f"PRAGMA busy_timeout={int(busy_timeout_ms)}")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA wal_autocheckpoint=1000")
    if initialize:
        ensure_ops_schema(conn)
    return conn


class SQLiteConnectionPool:
    """Small bounded pool; each checked-out connection has one request owner."""

    def __init__(
        self,
        db_path: str | Path,
        *,
        size: int = 4,
        busy_timeout_ms: int = DEFAULT_BUSY_TIMEOUT_MS,
    ):
        if size < 1:
            raise ValueError("SQLite pool size must be positive")
        self.db_path = Path(db_path).expanduser().resolve()
        self.size = int(size)
        self._queue: queue.LifoQueue[sqlite3.Connection] = queue.LifoQueue(size)
        self._lock = threading.Lock()
        self._closed = False
        self._checked_out = 0
        self._high_watermark = 0
        self._wait_count = 0
        self._wait_seconds_total = 0.0
        self._wait_seconds_max = 0.0
        self._sqlite_write_transactions = 0
        self._sqlite_write_wait_seconds = 0.0
        self._sqlite_write_wait_max_seconds = 0.0

        first = open_ops_connection(
            self.db_path,
            initialize=True,
            busy_timeout_ms=busy_timeout_ms,
        )
        self._queue.put(first)
        for _ in range(1, size):
            self._queue.put(
                open_ops_connection(
                    self.db_path,
                    initialize=False,
                    busy_timeout_ms=busy_timeout_ms,
                )
            )

    def acquire(self, *, timeout: float = 30.0) -> sqlite3.Connection:
        started = time.perf_counter()
        with self._lock:
            if self._closed:
                raise RuntimeError("SQLite connection pool is closed")
            if self._queue.empty():
                self._wait_count += 1
        try:
            conn = self._queue.get(timeout=timeout)
        except queue.Empty as exc:
            raise TimeoutError("Timed out waiting for a Sharur Ops database connection") from exc
        waited = time.perf_counter() - started
        with self._lock:
            self._checked_out += 1
            self._high_watermark = max(self._high_watermark, self._checked_out)
            self._wait_seconds_total += waited
            self._wait_seconds_max = max(self._wait_seconds_max, waited)
        return conn

    def release(self, conn: sqlite3.Connection) -> None:
        with contextlib.suppress(sqlite3.Error):
            if conn.in_transaction:
                conn.rollback()
        with self._lock:
            self._checked_out -= 1
            closed = self._closed
        if closed:
            conn.close()
        else:
            self._queue.put(conn)

    @contextlib.contextmanager
    def connection(self, *, timeout: float = 30.0) -> Iterator[sqlite3.Connection]:
        conn = self.acquire(timeout=timeout)
        try:
            yield conn
        finally:
            self.release(conn)

    def observe_sqlite_write_wait(self, duration_seconds: float) -> None:
        """Record time spent acquiring SQLite's serialized writer slot."""

        with self._lock:
            self._sqlite_write_transactions += 1
            self._sqlite_write_wait_seconds += duration_seconds
            self._sqlite_write_wait_max_seconds = max(
                self._sqlite_write_wait_max_seconds,
                duration_seconds,
            )

    def stats(self) -> dict[str, int | float]:
        with self._lock:
            return {
                "size": self.size,
                "available": self._queue.qsize(),
                "checked_out": self._checked_out,
                "high_watermark": self._high_watermark,
                "wait_count": self._wait_count,
                "wait_seconds_total": self._wait_seconds_total,
                "wait_seconds_max": self._wait_seconds_max,
                "sqlite_write_transactions": self._sqlite_write_transactions,
                "sqlite_write_wait_seconds": self._sqlite_write_wait_seconds,
                "sqlite_write_wait_max_seconds": self._sqlite_write_wait_max_seconds,
            }

    def close(self) -> None:
        with self._lock:
            if self._closed:
                return
            self._closed = True
        while True:
            try:
                self._queue.get_nowait().close()
            except queue.Empty:
                break


class SQLiteServerLock:
    """Advisory process lock enforcing one HTTP owner for a SQLite file."""

    def __init__(self, db_path: str | Path):
        self.db_path = Path(db_path).expanduser().resolve()
        self.path = self.db_path.with_name(f".{self.db_path.name}.server.lock")
        self._handle = None

    def acquire(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        handle = self.path.open("a+", encoding="utf-8")
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            handle.seek(0)
            owner = handle.read().strip() or "unknown owner"
            handle.close()
            raise RuntimeError(
                f"Sharur Ops database already has an HTTP owner: {owner}"
            ) from exc
        handle.seek(0)
        handle.truncate()
        handle.write(
            f"pid={os.getpid()} host={socket.gethostname()} acquired_ts={time.time():.6f}\n"
        )
        handle.flush()
        os.fsync(handle.fileno())
        self._handle = handle

    def release(self) -> None:
        if self._handle is None:
            return
        with contextlib.suppress(OSError):
            fcntl.flock(self._handle.fileno(), fcntl.LOCK_UN)
        self._handle.close()
        self._handle = None

    def __enter__(self) -> SQLiteServerLock:
        self.acquire()
        return self

    def __exit__(self, *_exc_info: object) -> None:
        self.release()


class SQLiteDirectAccessLock:
    """Shared lock that excludes the HTTP server's exclusive owner lock."""

    def __init__(self, db_path: str | Path):
        path = Path(db_path).expanduser().resolve()
        self.path = path.with_name(f".{path.name}.server.lock")
        self._handle = None

    def acquire(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        handle = self.path.open("a+", encoding="utf-8")
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_SH | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            handle.seek(0)
            owner = handle.read().strip() or "active HTTP server"
            handle.close()
            raise RuntimeError(
                "Direct SQLite access is unavailable while Sharur Ops owns "
                f"the database: {owner}"
            ) from exc
        self._handle = handle

    def release(self) -> None:
        if self._handle is None:
            return
        with contextlib.suppress(OSError):
            fcntl.flock(self._handle.fileno(), fcntl.LOCK_UN)
        self._handle.close()
        self._handle = None


__all__ = [
    "DEFAULT_BUSY_TIMEOUT_MS",
    "SQLiteConnectionPool",
    "SQLiteDirectAccessLock",
    "SQLiteServerLock",
    "open_ops_connection",
]
