"""Single-owner DuckDB runtime, admission control, cancellation, and telemetry."""

from __future__ import annotations

import asyncio
import contextlib
import fcntl
import hashlib
import re
import tempfile
import threading
import time
from collections import defaultdict, deque
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

import duckdb


if TYPE_CHECKING:
    from collections.abc import Callable

    from sharur.query.auth import QueryPrincipal


_SIZE_PATTERN = re.compile(
    r"(?:\d+(?:\.\d+)?)\s*(?:B|KB|MB|GB|TB|KiB|MiB|GiB|TiB)",
    flags=re.IGNORECASE,
)


class AdmissionRejectedError(RuntimeError):
    """Raised when the bounded query queue has reached capacity."""


class AdmissionTimeoutError(RuntimeError):
    """Raised when a query cannot acquire its weighted capacity in time."""


class DuplicateQueryIDError(RuntimeError):
    """Raised when two active requests use the same query identifier."""


class QueryCancelledError(RuntimeError):
    """Raised when an authenticated caller interrupts an active query."""


class QueryTimeoutError(RuntimeError):
    """Raised when a query exceeds its server-side execution limit."""


def _validate_size(value: str, name: str) -> str:
    normalized = value.strip()
    if not _SIZE_PATTERN.fullmatch(normalized):
        raise ValueError(f"{name} must include a supported unit, for example '16GB'")
    return normalized


@dataclass(eq=False)
class _AdmissionTicket:
    weight: int


class WeightedAdmissionController:
    """FIFO weighted admission with a bounded wait queue."""

    def __init__(self, *, capacity: int, max_queue: int):
        if capacity < 1:
            raise ValueError("Admission capacity must be positive")
        if max_queue < 0:
            raise ValueError("Admission queue bound must be non-negative")
        self.capacity = capacity
        self.max_queue = max_queue
        self._condition = asyncio.Condition()
        self._queue: deque[_AdmissionTicket] = deque()
        self._active_units = 0
        self._active_queries = 0
        self._admitted = 0
        self._rejected = 0
        self._timed_out = 0
        self._max_active_units = 0
        self._max_queued = 0
        self._wait_seconds = 0.0

    async def acquire(self, weight: int, *, timeout_seconds: float) -> float:
        if weight < 1 or weight > self.capacity:
            raise ValueError("Admission weight must fit within total capacity")
        if timeout_seconds <= 0:
            raise ValueError("Admission timeout must be positive")
        started = time.perf_counter()
        ticket = _AdmissionTicket(weight)
        async with self._condition:
            if not self._queue and self._active_units + weight <= self.capacity:
                self._admit(weight, 0.0)
                return 0.0
            if len(self._queue) >= self.max_queue:
                self._rejected += 1
                raise AdmissionRejectedError("Sharur query admission queue is full")
            self._queue.append(ticket)
            self._max_queued = max(self._max_queued, len(self._queue))
            try:
                await asyncio.wait_for(
                    self._wait_until_ready(ticket),
                    timeout=timeout_seconds,
                )
            except TimeoutError as exc:
                with contextlib.suppress(ValueError):
                    self._queue.remove(ticket)
                self._timed_out += 1
                self._condition.notify_all()
                raise AdmissionTimeoutError(
                    "Timed out waiting for Sharur query capacity"
                ) from exc
            except asyncio.CancelledError:
                with contextlib.suppress(ValueError):
                    self._queue.remove(ticket)
                self._condition.notify_all()
                raise
            wait_seconds = time.perf_counter() - started
            self._queue.popleft()
            self._admit(weight, wait_seconds)
            self._condition.notify_all()
            return wait_seconds

    async def _wait_until_ready(self, ticket: _AdmissionTicket) -> None:
        while (
            not self._queue
            or self._queue[0] is not ticket
            or self._active_units + ticket.weight > self.capacity
        ):
            await self._condition.wait()

    def _admit(self, weight: int, wait_seconds: float) -> None:
        self._active_units += weight
        self._active_queries += 1
        self._admitted += 1
        self._wait_seconds += wait_seconds
        self._max_active_units = max(self._max_active_units, self._active_units)

    async def release(self, weight: int) -> None:
        async with self._condition:
            self._active_units -= weight
            self._active_queries -= 1
            if self._active_units < 0 or self._active_queries < 0:
                raise RuntimeError("Sharur query admission accounting underflow")
            self._condition.notify_all()

    def snapshot(self) -> dict[str, int | float]:
        return {
            "capacity_units": self.capacity,
            "active_units": self._active_units,
            "active_queries": self._active_queries,
            "queued_queries": len(self._queue),
            "admitted_total": self._admitted,
            "rejected_total": self._rejected,
            "queue_timeouts_total": self._timed_out,
            "max_active_units": self._max_active_units,
            "max_queued_queries": self._max_queued,
            "wait_seconds_total": self._wait_seconds,
        }


class CursorStore:
    """DuckDBStore-compatible read facade over one thread-local cursor."""

    def __init__(
        self,
        cursor: duckdb.DuckDBPyConnection,
        *,
        db_path: Path,
        resource_budget: dict[str, Any],
        dataset_id: str | None = None,
    ):
        self._cursor = cursor
        self.db_path = db_path
        self.dataset_id = dataset_id
        self.read_only = True
        self._resource_budget = resource_budget

    @property
    def conn(self) -> duckdb.DuckDBPyConnection:
        return self._cursor

    @property
    def resource_budget(self) -> dict[str, Any]:
        return dict(self._resource_budget)

    def execute(
        self,
        query: str,
        params: dict | list | tuple | None = None,
    ) -> list[tuple]:
        if params:
            return self._cursor.execute(query, params).fetchall()
        return self._cursor.execute(query).fetchall()

    def execute_df(
        self,
        query: str,
        params: dict | list | tuple | None = None,
    ):
        if params:
            return self._cursor.execute(query, params).df()
        return self._cursor.execute(query).df()

    def interrupt(self) -> None:
        self._cursor.interrupt()


@dataclass
class _ActiveQuery:
    query_id: str
    principal: QueryPrincipal
    operator: str
    store: CursorStore
    started_at: float
    cancel_reason: str | None = None


@dataclass(frozen=True)
class QueryExecution:
    value: Any
    query_id: str
    operator: str
    execution_seconds: float


@dataclass
class _OperatorMetrics:
    queries: int = 0
    errors: int = 0
    cancellations: int = 0
    timeouts: int = 0
    duration_seconds: float = 0.0
    max_duration_seconds: float = 0.0


class ReadOnlyDuckDBRuntime:
    """Own one read-only DuckDB instance and share its cache across agents."""

    def __init__(
        self,
        db_path: str | Path,
        *,
        threads: int = 8,
        memory_limit: str = "16GB",
        temp_directory: str | Path,
        max_temp_directory_size: str = "128GB",
        dataset_id: str | None = None,
    ):
        if isinstance(threads, bool) or not isinstance(threads, int) or threads < 1:
            raise ValueError("threads must be a positive integer")
        self.db_path = Path(db_path).expanduser().resolve()
        self.threads = threads
        self.memory_limit = _validate_size(memory_limit, "memory_limit")
        self.temp_directory = Path(temp_directory).expanduser().resolve()
        self.max_temp_directory_size = _validate_size(
            max_temp_directory_size,
            "max_temp_directory_size",
        )
        self.dataset_id = dataset_id
        self._root: duckdb.DuckDBPyConnection | None = None
        self._thread_local = threading.local()
        self._state_lock = threading.RLock()
        self._cursors: list[duckdb.DuckDBPyConnection] = []
        self._active: dict[str, _ActiveQuery] = {}
        self._operator_metrics: dict[str, _OperatorMetrics] = defaultdict(
            _OperatorMetrics
        )
        self._opened_at: float | None = None
        self._settings: dict[str, str] = {}
        self._generation = 0
        self._owner_lock_handle: Any = None
        lock_root = Path(tempfile.gettempdir()) / "sharur-query-locks"
        lock_key = hashlib.sha256(str(self.db_path).encode("utf-8")).hexdigest()[:24]
        self._owner_lock_path = lock_root / f"{lock_key}.lock"

    @property
    def is_open(self) -> bool:
        return self._root is not None

    @property
    def resource_budget(self) -> dict[str, Any]:
        return {
            "threads": self.threads,
            "memory_limit": self.memory_limit,
            "temp_directory": str(self.temp_directory),
            "max_temp_directory_size": self.max_temp_directory_size,
        }

    def open(self) -> None:
        with self._state_lock:
            if self._root is not None:
                return
            if not self.db_path.is_file():
                raise FileNotFoundError(f"DuckDB file does not exist: {self.db_path}")
            self.temp_directory.mkdir(parents=True, exist_ok=True)
            self._owner_lock_path.parent.mkdir(parents=True, exist_ok=True)
            owner_lock = self._owner_lock_path.open("a+b")
            try:
                fcntl.flock(
                    owner_lock.fileno(),
                    fcntl.LOCK_EX | fcntl.LOCK_NB,
                )
            except BlockingIOError as exc:
                owner_lock.close()
                raise RuntimeError(
                    f"Another Sharur Query process already owns {self.db_path}"
                ) from exc
            root: duckdb.DuckDBPyConnection | None = None
            try:
                # DuckDB eagerly spawns cpu_count-1 scheduler threads at import
                # for the module default connection; cap it so a 224-core host
                # does not carry ~223 idle threads in the query-server process.
                duckdb.execute("SET threads TO 1")
                root = duckdb.connect(
                    str(self.db_path),
                    read_only=True,
                    config={"threads": self.threads},
                )
                root.execute("SET threads = ?", [self.threads])
                root.execute("SET memory_limit = ?", [self.memory_limit])
                root.execute("SET temp_directory = ?", [str(self.temp_directory)])
                root.execute(
                    "SET max_temp_directory_size = ?",
                    [self.max_temp_directory_size],
                )
                root.execute("SET enable_external_access = false")
                root.execute("SET autoload_known_extensions = false")
                root.execute("SET autoinstall_known_extensions = false")
                root.execute("SET lock_configuration = true")
                settings = root.execute(
                    """
                    SELECT name, value
                    FROM duckdb_settings()
                    WHERE name IN (
                        'threads', 'memory_limit', 'temp_directory',
                        'max_temp_directory_size', 'enable_external_access',
                        'autoload_known_extensions',
                        'autoinstall_known_extensions', 'lock_configuration'
                    )
                    ORDER BY name
                    """
                ).fetchall()
            except Exception:
                if root is not None:
                    root.close()
                fcntl.flock(owner_lock.fileno(), fcntl.LOCK_UN)
                owner_lock.close()
                raise
            assert root is not None
            self._owner_lock_handle = owner_lock
            self._root = root
            self._settings = {str(name): str(value) for name, value in settings}
            self._opened_at = time.time()
            self._generation += 1

    def close(self) -> None:
        with self._state_lock:
            active = list(self._active.values())
            cursors = list(self._cursors)
            root = self._root
            owner_lock = self._owner_lock_handle
            self._active.clear()
            self._cursors.clear()
            self._root = None
            self._owner_lock_handle = None
        for query in active:
            with contextlib.suppress(Exception):
                query.store.interrupt()
        for cursor in cursors:
            with contextlib.suppress(Exception):
                cursor.close()
        if root is not None:
            with contextlib.suppress(Exception):
                root.close()
        if owner_lock is not None:
            with contextlib.suppress(Exception):
                fcntl.flock(owner_lock.fileno(), fcntl.LOCK_UN)
                owner_lock.close()

    def _store_for_thread(self) -> CursorStore:
        existing = getattr(self._thread_local, "store", None)
        generation = getattr(self._thread_local, "generation", None)
        if existing is not None and generation == self._generation:
            return existing
        with self._state_lock:
            if self._root is None:
                raise RuntimeError("Sharur query runtime is closed")
            cursor = self._root.cursor()
            self._cursors.append(cursor)
        store = CursorStore(
            cursor,
            db_path=self.db_path,
            resource_budget=self.resource_budget,
            dataset_id=self.dataset_id,
        )
        self._thread_local.store = store
        self._thread_local.generation = self._generation
        return store

    def cursor_store(self) -> CursorStore:
        """Return this thread's cursor over the shared read-only DuckDB instance.

        Offline bounded workflows use this surface to parallelize typed
        operators while retaining one database owner, buffer cache, memory
        ceiling, thread pool, and spill budget. HTTP execution continues
        through :meth:`execute` for cancellation and telemetry.
        """
        return self._store_for_thread()

    def _discard_store(self, store: CursorStore) -> None:
        with self._state_lock, contextlib.suppress(ValueError):
            self._cursors.remove(store.conn)
        if getattr(self._thread_local, "store", None) is store:
            self._thread_local.store = None
            self._thread_local.generation = None
        with contextlib.suppress(Exception):
            store.conn.close()

    def execute(
        self,
        *,
        query_id: str,
        principal: QueryPrincipal,
        operator: str,
        timeout_seconds: float,
        operation: Callable[[CursorStore], Any],
    ) -> QueryExecution:
        if timeout_seconds <= 0:
            raise ValueError("Query timeout must be positive")
        store = self._store_for_thread()
        active = _ActiveQuery(
            query_id=query_id,
            principal=principal,
            operator=operator,
            store=store,
            started_at=time.perf_counter(),
        )
        with self._state_lock:
            if query_id in self._active:
                raise DuplicateQueryIDError(
                    f"Query ID is already active: {query_id}"
                )
            self._active[query_id] = active

        completed = threading.Event()

        def expire() -> None:
            if not completed.is_set():
                self._interrupt(query_id, reason="timeout")

        timer = threading.Timer(timeout_seconds, expire)
        timer.daemon = True
        timer.start()
        error: BaseException | None = None
        try:
            value = operation(store)
            if active.cancel_reason == "timeout":
                raise QueryTimeoutError(
                    f"Query exceeded its {timeout_seconds:g}s execution limit"
                )
            if active.cancel_reason is not None:
                raise QueryCancelledError(f"Query {query_id} was cancelled")
            return QueryExecution(
                value=value,
                query_id=query_id,
                operator=operator,
                execution_seconds=time.perf_counter() - active.started_at,
            )
        except BaseException as exc:
            error = exc
            if isinstance(exc, (QueryTimeoutError, QueryCancelledError)):
                raise
            if active.cancel_reason == "timeout":
                raise QueryTimeoutError(
                    f"Query exceeded its {timeout_seconds:g}s execution limit"
                ) from exc
            if active.cancel_reason is not None:
                raise QueryCancelledError(
                    f"Query {query_id} was cancelled"
                ) from exc
            raise
        finally:
            completed.set()
            timer.cancel()
            duration = time.perf_counter() - active.started_at
            interrupted = active.cancel_reason is not None
            with self._state_lock:
                self._active.pop(query_id, None)
                metrics = self._operator_metrics[operator]
                metrics.queries += 1
                metrics.duration_seconds += duration
                metrics.max_duration_seconds = max(
                    metrics.max_duration_seconds,
                    duration,
                )
                if error is not None:
                    metrics.errors += 1
                if active.cancel_reason == "timeout":
                    metrics.timeouts += 1
                elif active.cancel_reason is not None:
                    metrics.cancellations += 1
            if interrupted:
                self._discard_store(store)

    def _interrupt(self, query_id: str, *, reason: str) -> bool:
        with self._state_lock:
            active = self._active.get(query_id)
            if active is None:
                return False
            if active.cancel_reason is None:
                active.cancel_reason = reason
            store = active.store
        store.interrupt()
        return True

    def cancel(
        self,
        query_id: str,
        *,
        principal: QueryPrincipal,
        reason: str = "client",
    ) -> bool:
        with self._state_lock:
            active = self._active.get(query_id)
            if active is None:
                return False
            if (
                active.principal.agent_id != principal.agent_id
                and principal.role != "operator"
            ):
                raise PermissionError(
                    f"Credential {principal.agent_id} cannot cancel query {query_id}"
                )
        return self._interrupt(query_id, reason=reason)

    def snapshot(self) -> dict[str, Any]:
        with self._state_lock:
            operators = {
                name: {
                    "queries": metrics.queries,
                    "errors": metrics.errors,
                    "cancellations": metrics.cancellations,
                    "timeouts": metrics.timeouts,
                    "duration_seconds": metrics.duration_seconds,
                    "max_duration_seconds": metrics.max_duration_seconds,
                }
                for name, metrics in sorted(self._operator_metrics.items())
            }
            active_by_operator: dict[str, int] = defaultdict(int)
            for query in self._active.values():
                active_by_operator[query.operator] += 1
            return {
                "open": self._root is not None,
                "database": str(self.db_path),
                "dataset_id": self.dataset_id,
                "owner_lock": str(self._owner_lock_path),
                "opened_at": self._opened_at,
                "settings": dict(self._settings),
                "cursor_count": len(self._cursors),
                "active_queries": len(self._active),
                "active_by_operator": dict(sorted(active_by_operator.items())),
                "operators": operators,
            }


__all__ = [
    "AdmissionRejectedError",
    "AdmissionTimeoutError",
    "CursorStore",
    "DuplicateQueryIDError",
    "QueryCancelledError",
    "QueryExecution",
    "QueryTimeoutError",
    "ReadOnlyDuckDBRuntime",
    "WeightedAdmissionController",
]
