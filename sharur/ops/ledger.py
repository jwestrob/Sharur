"""Durable run and stage-attempt ledger stored in the Sharur ops SQLite file."""

from __future__ import annotations

import contextlib
import json
import sqlite3
import threading
import time
import uuid
from pathlib import Path
from typing import TYPE_CHECKING, Any

from sharur.ops.db import SQLiteDirectAccessLock, open_ops_connection
from sharur.ops.schema import ensure_ops_schema


if TYPE_CHECKING:
    from collections.abc import Callable


_JSON_FIELDS = {
    "config",
    "result",
    "payload",
    "command",
    "inputs",
    "outputs",
    "resource_profile",
}


def _canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)


def _decode_row(row: sqlite3.Row | None) -> dict[str, Any] | None:
    if row is None:
        return None
    result = dict(row)
    for field in _JSON_FIELDS:
        value = result.get(field)
        if isinstance(value, str):
            with contextlib.suppress(json.JSONDecodeError):
                result[field] = json.loads(value)
    return result


class StageAlreadyCompleteError(RuntimeError):
    """Raised when a duplicate launcher reaches an already successful stage."""

    def __init__(self, stage: dict[str, Any]):
        self.stage = stage
        super().__init__(f"Stage {stage['stage_id']} already completed in run {stage['run_id']}")


class RunLedger:
    """Append-oriented record of runs, events, and stage attempts."""

    def __init__(
        self,
        db_path: str | Path,
        *,
        connection: sqlite3.Connection | None = None,
        lock: threading.RLock | None = None,
        initialize: bool = True,
        close_callback: Callable[[sqlite3.Connection], None] | None = None,
        agent_id: str = "operator",
        transaction_wait_observer: Callable[[float], None] | None = None,
    ):
        self.db_path = Path(db_path)
        self.agent_id = agent_id
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._lock = lock or threading.RLock()
        self._close_callback = close_callback
        self._transaction_wait_observer = transaction_wait_observer
        self._owns_connection = connection is None
        self._direct_lock: SQLiteDirectAccessLock | None = None
        if connection is None:
            self._direct_lock = SQLiteDirectAccessLock(self.db_path)
            self._direct_lock.acquire()
            try:
                self._conn = open_ops_connection(
                    self.db_path,
                    initialize=initialize,
                )
            except Exception:
                self._direct_lock.release()
                self._direct_lock = None
                raise
        else:
            self._conn = connection
        if connection is not None:
            self._conn.row_factory = sqlite3.Row
            self._conn.execute("PRAGMA foreign_keys=ON")
            if initialize:
                with self._lock:
                    ensure_ops_schema(self._conn)

    def _begin_immediate(self) -> None:
        started = time.perf_counter()
        self._conn.execute("BEGIN IMMEDIATE")
        if self._transaction_wait_observer is not None:
            self._transaction_wait_observer(time.perf_counter() - started)

    def create_run(
        self,
        run_type: str,
        dataset_path: str | Path,
        *,
        created_by: str = "operator",
        config: dict[str, Any] | None = None,
        idempotency_key: str | None = None,
        parent_run_id: str | None = None,
        campaign_id: str | None = None,
    ) -> str:
        if not run_type.strip():
            raise ValueError("run_type must be non-empty")
        if idempotency_key is not None and not idempotency_key.strip():
            raise ValueError("idempotency_key must be non-empty when provided")
        normalized_path = str(Path(dataset_path).expanduser().resolve())
        config_json = _canonical_json(config or {})
        with self._lock:
            self._begin_immediate()
            try:
                if idempotency_key:
                    existing = self._conn.execute(
                        "SELECT * FROM runs WHERE idempotency_key = ?",
                        (idempotency_key,),
                    ).fetchone()
                    if existing is not None:
                        if (
                            existing["run_type"] != run_type
                            or existing["dataset_path"] != normalized_path
                            or existing["config"] != config_json
                            or existing["created_by"] != created_by
                            or existing["parent_run_id"] != parent_run_id
                            or existing["campaign_id"] != campaign_id
                        ):
                            raise ValueError(
                                "Run idempotency key already exists with a different payload"
                            )
                        self._conn.commit()
                        return str(existing["id"])

                run_id = str(uuid.uuid4())
                now = time.time()
                self._conn.execute(
                    """
                    INSERT INTO runs(
                        id, idempotency_key, run_type, dataset_path, created_by,
                        status, created_ts, parent_run_id, config, campaign_id
                    ) VALUES (?, ?, ?, ?, ?, 'pending', ?, ?, ?, ?)
                    """,
                    (
                        run_id,
                        idempotency_key,
                        run_type,
                        normalized_path,
                        created_by,
                        now,
                        parent_run_id,
                        config_json,
                        campaign_id,
                    ),
                )
                self._append_event_locked(
                    run_id,
                    "run_created",
                    payload=config or {},
                    actor_agent_id=created_by,
                )
                self._conn.commit()
                return run_id
            except Exception:
                self._conn.rollback()
                raise

    def start_run(self, run_id: str) -> dict[str, Any]:
        now = time.time()
        with self._lock:
            cursor = self._conn.execute(
                """
                UPDATE runs
                SET status = 'running',
                    started_ts = COALESCE(started_ts, ?),
                    heartbeat_ts = ?
                WHERE id = ? AND status = 'pending'
                """,
                (now, now, run_id),
            )
            if cursor.rowcount == 0:
                raise ValueError(f"Run {run_id} is missing or not startable")
            self._append_event_locked(run_id, "run_started")
            self._conn.commit()
            return self.get_run(run_id)

    def submit_run(self, run_id: str) -> dict[str, Any]:
        """Atomically mark a scheduler bundle submitted without starting its work."""
        now = time.time()
        with self._lock:
            cursor = self._conn.execute(
                """
                UPDATE runs
                SET status = 'submitted', heartbeat_ts = ?
                WHERE id = ? AND status = 'pending'
                """,
                (now, run_id),
            )
            if cursor.rowcount == 0:
                raise ValueError(f"Run {run_id} is missing or not submittable")
            self._append_event_locked(run_id, "run_submitted")
            self._conn.commit()
            return self.get_run(run_id)

    def heartbeat_run(self, run_id: str) -> dict[str, Any]:
        with self._lock:
            now = time.time()
            cursor = self._conn.execute(
                "UPDATE runs SET heartbeat_ts = ? WHERE id = ? AND status = 'running'",
                (now, run_id),
            )
            if cursor.rowcount == 0:
                raise ValueError(f"Run {run_id} is not running")
            self._append_event_locked(run_id, "run_heartbeat")
            self._conn.commit()
            return self.get_run(run_id)

    def wait_for_scheduler(self, run_id: str) -> dict[str, Any]:
        """Return an active run to submitted state when no stage is executing."""
        now = time.time()
        with self._lock:
            self._begin_immediate()
            try:
                run = self._conn.execute(
                    "SELECT status FROM runs WHERE id = ?",
                    (run_id,),
                ).fetchone()
                if run is None:
                    raise ValueError(f"Run {run_id} does not exist")
                if run["status"] in {"complete", "failed"}:
                    self._conn.commit()
                    return self.get_run(run_id)
                if run["status"] not in {"submitted", "running"}:
                    raise ValueError(f"Run {run_id} is not scheduler-managed and active")
                running_stage = self._conn.execute(
                    """
                    SELECT stage_id FROM run_stages
                    WHERE run_id = ? AND status = 'running'
                    ORDER BY started_ts DESC
                    LIMIT 1
                    """,
                    (run_id,),
                ).fetchone()
                if running_stage is None:
                    self._conn.execute(
                        """
                        UPDATE runs
                        SET status = 'submitted', heartbeat_ts = ?, current_stage = NULL
                        WHERE id = ?
                        """,
                        (now, run_id),
                    )
                    self._append_event_locked(run_id, "run_waiting_for_scheduler")
                else:
                    self._conn.execute(
                        """
                        UPDATE runs
                        SET status = 'running', heartbeat_ts = ?, current_stage = ?
                        WHERE id = ?
                        """,
                        (now, running_stage["stage_id"], run_id),
                    )
                    self._append_event_locked(
                        run_id,
                        "run_scheduler_state_refreshed",
                        stage_id=running_stage["stage_id"],
                    )
                self._conn.commit()
                return self.get_run(run_id)
            except Exception:
                self._conn.rollback()
                raise

    def complete_run(
        self,
        run_id: str,
        *,
        result: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        now = time.time()
        with self._lock:
            cursor = self._conn.execute(
                """
                UPDATE runs
                SET status = 'complete', heartbeat_ts = ?, completed_ts = ?,
                    current_stage = NULL, result = ?, error = NULL
                WHERE id = ? AND status IN ('pending', 'submitted', 'running')
                """,
                (now, now, _canonical_json(result or {}), run_id),
            )
            if cursor.rowcount == 0:
                raise ValueError(f"Run {run_id} is missing or not completable")
            self._append_event_locked(run_id, "run_completed", payload=result or {})
            self._conn.commit()
            return self.get_run(run_id)

    def fail_run(self, run_id: str, error: str) -> dict[str, Any]:
        now = time.time()
        with self._lock:
            cursor = self._conn.execute(
                """
                UPDATE runs
                SET status = 'failed', heartbeat_ts = ?, completed_ts = ?,
                    current_stage = NULL, error = ?
                WHERE id = ? AND status IN ('pending', 'submitted', 'running')
                """,
                (now, now, error, run_id),
            )
            if cursor.rowcount == 0:
                raise ValueError(f"Run {run_id} is missing or not fail-able")
            self._append_event_locked(
                run_id,
                "run_failed",
                payload={"error": error},
            )
            self._conn.commit()
            return self.get_run(run_id)

    def get_run(self, run_id: str) -> dict[str, Any]:
        row = self._conn.execute(
            "SELECT * FROM runs WHERE id = ?",
            (run_id,),
        ).fetchone()
        decoded = _decode_row(row)
        if decoded is None:
            raise KeyError(f"Run {run_id} not found")
        return decoded

    def list_runs(
        self,
        *,
        dataset_path: str | Path | None = None,
        run_type: str | None = None,
        status: str | None = None,
        campaign_id: str | None = None,
        before_ts: float | None = None,
        limit: int = 50,
    ) -> list[dict[str, Any]]:
        clauses: list[str] = []
        params: list[Any] = []
        if dataset_path is not None:
            clauses.append("dataset_path = ?")
            params.append(str(Path(dataset_path).expanduser().resolve()))
        if run_type is not None:
            clauses.append("run_type = ?")
            params.append(run_type)
        if status is not None:
            clauses.append("status = ?")
            params.append(status)
        if campaign_id is not None:
            clauses.append("campaign_id = ?")
            params.append(campaign_id)
        if before_ts is not None:
            clauses.append("created_ts < ?")
            params.append(before_ts)
        where = f"WHERE {' AND '.join(clauses)}" if clauses else ""
        params.append(limit)
        rows = self._conn.execute(
            f"SELECT * FROM runs {where} ORDER BY created_ts DESC LIMIT ?",
            params,
        ).fetchall()
        return [_decode_row(row) for row in rows]  # type: ignore[misc]

    def list_stages(
        self,
        run_id: str,
        *,
        stage_id: str | None = None,
    ) -> list[dict[str, Any]]:
        clauses = ["run_id = ?"]
        params: list[Any] = [run_id]
        if stage_id is not None:
            clauses.append("stage_id = ?")
            params.append(stage_id)
        rows = self._conn.execute(
            f"""
            SELECT * FROM run_stages
            WHERE {" AND ".join(clauses)}
            ORDER BY started_ts, stage_id, attempt
            """,
            params,
        ).fetchall()
        return [_decode_row(row) for row in rows]  # type: ignore[misc]

    def latest_stage(
        self,
        run_id: str,
        stage_id: str,
    ) -> dict[str, Any] | None:
        row = self._conn.execute(
            """
            SELECT * FROM run_stages
            WHERE run_id = ? AND stage_id = ?
            ORDER BY started_ts DESC, attempt DESC
            LIMIT 1
            """,
            (run_id, stage_id),
        ).fetchone()
        return _decode_row(row)

    def append_event(
        self,
        run_id: str,
        event_type: str,
        *,
        stage_id: str | None = None,
        attempt: int | None = None,
        payload: dict[str, Any] | None = None,
    ) -> int:
        with self._lock:
            event_id = self._append_event_locked(
                run_id,
                event_type,
                stage_id=stage_id,
                attempt=attempt,
                payload=payload,
            )
            self._conn.commit()
            return event_id

    def _append_event_locked(
        self,
        run_id: str,
        event_type: str,
        *,
        stage_id: str | None = None,
        attempt: int | None = None,
        payload: dict[str, Any] | None = None,
        actor_agent_id: str | None = None,
    ) -> int:
        now = time.time()
        cursor = self._conn.execute(
            """
            INSERT INTO run_events(run_id, ts, event_type, stage_id, attempt, payload)
            VALUES (?, ?, ?, ?, ?, ?)
            """,
            (
                run_id,
                now,
                event_type,
                stage_id,
                attempt,
                _canonical_json(payload or {}),
            ),
        )
        run = self._conn.execute(
            "SELECT campaign_id FROM runs WHERE id = ?",
            (run_id,),
        ).fetchone()
        self._conn.execute(
            """
            INSERT INTO coordination_events(
                ts, event_type, actor_agent_id, campaign_id, entity_type,
                entity_id, run_id, payload
            ) VALUES (?, ?, ?, ?, 'run', ?, ?, ?)
            """,
            (
                now,
                event_type,
                actor_agent_id or self.agent_id,
                run["campaign_id"] if run is not None else None,
                run_id,
                run_id,
                _canonical_json(
                    {
                        "stage_id": stage_id,
                        "attempt": attempt,
                        **(payload or {}),
                    }
                ),
            ),
        )
        return int(cursor.lastrowid)

    def start_stage(
        self,
        run_id: str,
        stage_id: str,
        signature: str,
        *,
        command: list[str] | None = None,
        inputs: dict[str, Any] | None = None,
        outputs: dict[str, Any] | None = None,
        resource_profile: dict[str, Any] | None = None,
    ) -> int:
        now = time.time()
        with self._lock:
            self._begin_immediate()
            try:
                run = self._conn.execute(
                    "SELECT status FROM runs WHERE id = ?",
                    (run_id,),
                ).fetchone()
                if run is None:
                    raise ValueError(f"Run {run_id} does not exist")
                latest = self._conn.execute(
                    """
                    SELECT * FROM run_stages
                    WHERE run_id = ? AND stage_id = ?
                    ORDER BY started_ts DESC, attempt DESC
                    LIMIT 1
                    """,
                    (run_id, stage_id),
                ).fetchone()
                if latest is not None and latest["status"] == "running":
                    raise ValueError(
                        f"Stage {stage_id} already has a running attempt in run {run_id}"
                    )
                if latest is not None and latest["status"] in {"complete", "reused"}:
                    decoded = _decode_row(latest)
                    assert decoded is not None
                    if latest["signature"] != signature:
                        raise ValueError(
                            f"Stage {stage_id} already completed with a different signature"
                        )
                    raise StageAlreadyCompleteError(decoded)
                if run["status"] not in {"pending", "submitted", "running"}:
                    raise ValueError(f"Run {run_id} is not active")

                attempt = int(
                    self._conn.execute(
                        """
                        SELECT COALESCE(MAX(attempt), 0) + 1
                        FROM run_stages WHERE run_id = ? AND stage_id = ?
                        """,
                        (run_id, stage_id),
                    ).fetchone()[0]
                )
                self._conn.execute(
                    """
                    INSERT INTO run_stages(
                        run_id, stage_id, attempt, signature, status, started_ts,
                        heartbeat_ts, command, inputs, outputs, resource_profile,
                        owner_agent_id
                    ) VALUES (?, ?, ?, ?, 'running', ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        run_id,
                        stage_id,
                        attempt,
                        signature,
                        now,
                        now,
                        _canonical_json(command or []),
                        _canonical_json(inputs or {}),
                        _canonical_json(outputs or {}),
                        _canonical_json(resource_profile or {}),
                        self.agent_id,
                    ),
                )
                self._conn.execute(
                    """
                    UPDATE runs
                    SET status = 'running', started_ts = COALESCE(started_ts, ?),
                        heartbeat_ts = ?, current_stage = ?
                    WHERE id = ?
                    """,
                    (now, now, stage_id, run_id),
                )
                self._append_event_locked(
                    run_id,
                    "stage_started",
                    stage_id=stage_id,
                    attempt=attempt,
                    payload={"signature": signature},
                )
                self._conn.commit()
                return attempt
            except Exception:
                self._conn.rollback()
                raise

    def heartbeat_stage(self, run_id: str, stage_id: str, attempt: int) -> None:
        now = time.time()
        with self._lock:
            cursor = self._conn.execute(
                """
                UPDATE run_stages SET heartbeat_ts = ?
                WHERE run_id = ? AND stage_id = ? AND attempt = ?
                  AND status = 'running'
                """,
                (now, run_id, stage_id, attempt),
            )
            if cursor.rowcount == 0:
                raise ValueError("Stage attempt is not running")
            self._conn.execute(
                "UPDATE runs SET heartbeat_ts = ? WHERE id = ?",
                (now, run_id),
            )
            self._append_event_locked(
                run_id,
                "stage_heartbeat",
                stage_id=stage_id,
                attempt=attempt,
            )
            self._conn.commit()

    def complete_stage(
        self,
        run_id: str,
        stage_id: str,
        attempt: int,
        *,
        outputs: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        now = time.time()
        with self._lock:
            cursor = self._conn.execute(
                """
                UPDATE run_stages
                SET status = 'complete', heartbeat_ts = ?, completed_ts = ?,
                    outputs = ?, error = NULL
                WHERE run_id = ? AND stage_id = ? AND attempt = ?
                  AND status = 'running'
                """,
                (
                    now,
                    now,
                    _canonical_json(outputs or {}),
                    run_id,
                    stage_id,
                    attempt,
                ),
            )
            if cursor.rowcount == 0:
                raise ValueError("Stage attempt is not running")
            self._conn.execute(
                "UPDATE runs SET heartbeat_ts = ? WHERE id = ?",
                (now, run_id),
            )
            self._refresh_run_stage_state_locked(run_id, now)
            self._append_event_locked(
                run_id,
                "stage_completed",
                stage_id=stage_id,
                attempt=attempt,
                payload=outputs or {},
            )
            self._conn.commit()
            row = self._conn.execute(
                """
                SELECT * FROM run_stages
                WHERE run_id = ? AND stage_id = ? AND attempt = ?
                """,
                (run_id, stage_id, attempt),
            ).fetchone()
            return _decode_row(row)  # type: ignore[return-value]

    def fail_stage(
        self,
        run_id: str,
        stage_id: str,
        attempt: int,
        error: str,
    ) -> dict[str, Any]:
        now = time.time()
        with self._lock:
            cursor = self._conn.execute(
                """
                UPDATE run_stages
                SET status = 'failed', heartbeat_ts = ?, completed_ts = ?, error = ?
                WHERE run_id = ? AND stage_id = ? AND attempt = ?
                  AND status = 'running'
                """,
                (now, now, error, run_id, stage_id, attempt),
            )
            if cursor.rowcount == 0:
                raise ValueError("Stage attempt is not running")
            self._conn.execute(
                "UPDATE runs SET heartbeat_ts = ? WHERE id = ?",
                (now, run_id),
            )
            self._refresh_run_stage_state_locked(run_id, now)
            self._append_event_locked(
                run_id,
                "stage_failed",
                stage_id=stage_id,
                attempt=attempt,
                payload={"error": error},
            )
            self._conn.commit()
            row = self._conn.execute(
                """
                SELECT * FROM run_stages
                WHERE run_id = ? AND stage_id = ? AND attempt = ?
                """,
                (run_id, stage_id, attempt),
            ).fetchone()
            return _decode_row(row)  # type: ignore[return-value]

    def reuse_stage(
        self,
        run_id: str,
        stage_id: str,
        signature: str,
        source: dict[str, Any],
    ) -> None:
        now = time.time()
        with self._lock:
            self._begin_immediate()
            try:
                run = self._conn.execute(
                    "SELECT status FROM runs WHERE id = ?",
                    (run_id,),
                ).fetchone()
                if run is None or run["status"] not in {"pending", "submitted", "running"}:
                    raise ValueError(f"Run {run_id} is not active")
                existing = self._conn.execute(
                    """
                    SELECT * FROM run_stages
                    WHERE run_id = ? AND stage_id = ? AND attempt = 0
                    """,
                    (run_id, stage_id),
                ).fetchone()
                if existing is not None:
                    if (
                        existing["status"] != "reused"
                        or existing["signature"] != signature
                        or existing["reused_from_run_id"] != source["run_id"]
                        or existing["reused_from_attempt"] != source["attempt"]
                    ):
                        raise ValueError(
                            f"Conflicting reuse record for stage {stage_id} in run {run_id}"
                        )
                    self._conn.commit()
                    return

                self._conn.execute(
                    """
                    INSERT INTO run_stages(
                        run_id, stage_id, attempt, signature, status, started_ts,
                        heartbeat_ts, completed_ts, command, inputs, outputs,
                        resource_profile, reused_from_run_id, reused_from_attempt
                    ) VALUES (?, ?, 0, ?, 'reused', ?, ?, ?, '[]', '{}', ?, '{}', ?, ?)
                    """,
                    (
                        run_id,
                        stage_id,
                        signature,
                        now,
                        now,
                        now,
                        _canonical_json(source.get("outputs") or {}),
                        source["run_id"],
                        source["attempt"],
                    ),
                )
                self._append_event_locked(
                    run_id,
                    "stage_reused",
                    stage_id=stage_id,
                    attempt=0,
                    payload={
                        "source_run_id": source["run_id"],
                        "source_attempt": source["attempt"],
                        "signature": signature,
                    },
                )
                self._conn.commit()
            except Exception:
                self._conn.rollback()
                raise

    def find_reusable_stage(
        self,
        dataset_path: str | Path,
        run_type: str,
        stage_id: str,
        signature: str,
    ) -> dict[str, Any] | None:
        normalized = str(Path(dataset_path).expanduser().resolve())
        row = self._conn.execute(
            """
            SELECT rs.*
            FROM run_stages rs
            JOIN runs r ON r.id = rs.run_id
            WHERE r.dataset_path = ? AND r.run_type = ?
              AND rs.stage_id = ? AND rs.signature = ?
              AND rs.status IN ('complete', 'reused')
            ORDER BY rs.completed_ts DESC
            LIMIT 1
            """,
            (normalized, run_type, stage_id, signature),
        ).fetchone()
        return _decode_row(row)

    def recover_stale_runs(
        self,
        *,
        stale_after_seconds: int = 300,
        now: float | None = None,
    ) -> dict[str, list[Any]]:
        """Fail abandoned running stages/runs whose durable heartbeat expired."""
        if stale_after_seconds < 1:
            raise ValueError("stale_after_seconds must be positive")
        current = time.time() if now is None else now
        cutoff = current - stale_after_seconds
        recovered_stages: list[dict[str, Any]] = []
        recovered_runs: set[str] = set()
        with self._lock:
            self._begin_immediate()
            try:
                rows = self._conn.execute(
                    """
                    SELECT run_id, stage_id, attempt
                    FROM run_stages
                    WHERE status = 'running'
                      AND COALESCE(heartbeat_ts, started_ts, 0) <= ?
                    """,
                    (cutoff,),
                ).fetchall()
                for row in rows:
                    error = f"stage heartbeat expired; recovered after {stale_after_seconds}s"
                    self._conn.execute(
                        """
                        UPDATE run_stages
                        SET status = 'failed', completed_ts = ?, error = ?
                        WHERE run_id = ? AND stage_id = ? AND attempt = ?
                          AND status = 'running'
                        """,
                        (
                            current,
                            error,
                            row["run_id"],
                            row["stage_id"],
                            row["attempt"],
                        ),
                    )
                    self._append_event_locked(
                        row["run_id"],
                        "stage_recovered_after_crash",
                        stage_id=row["stage_id"],
                        attempt=row["attempt"],
                        payload={"error": error},
                    )
                    recovered_stages.append(
                        {
                            "run_id": row["run_id"],
                            "stage_id": row["stage_id"],
                            "attempt": row["attempt"],
                        }
                    )
                    recovered_runs.add(row["run_id"])
                    self._refresh_run_stage_state_locked(row["run_id"], current)

                stale_runs = self._conn.execute(
                    """
                    SELECT id FROM runs
                    WHERE status = 'running'
                      AND COALESCE(heartbeat_ts, started_ts, created_ts) <= ?
                    """,
                    (cutoff,),
                ).fetchall()
                recovered_runs.update(row["id"] for row in stale_runs)
                for run_id in recovered_runs:
                    error = f"run heartbeat expired; recovered after {stale_after_seconds}s"
                    self._conn.execute(
                        """
                        UPDATE runs
                        SET status = 'failed', completed_ts = ?, current_stage = NULL,
                            error = ?
                        WHERE id = ? AND status = 'running'
                        """,
                        (current, error, run_id),
                    )
                    self._append_event_locked(
                        run_id,
                        "run_recovered_after_crash",
                        payload={"error": error},
                    )
                self._conn.commit()
            except Exception:
                self._conn.rollback()
                raise
        return {
            "runs": sorted(recovered_runs),
            "stages": recovered_stages,
        }

    def _refresh_run_stage_state_locked(self, run_id: str, now: float) -> None:
        """Keep the run summary coherent when stages execute concurrently."""

        active = self._conn.execute(
            """
            SELECT stage_id
            FROM run_stages
            WHERE run_id = ? AND status = 'running'
            ORDER BY started_ts DESC, stage_id
            LIMIT 1
            """,
            (run_id,),
        ).fetchone()
        if active is None:
            self._conn.execute(
                "UPDATE runs SET current_stage = NULL, heartbeat_ts = ? WHERE id = ?",
                (now, run_id),
            )
        else:
            self._conn.execute(
                """
                UPDATE runs
                SET status = 'running', current_stage = ?, heartbeat_ts = ?
                WHERE id = ? AND status IN ('pending', 'submitted', 'running')
                """,
                (active["stage_id"], now, run_id),
            )

    def events(
        self,
        run_id: str,
        *,
        after_id: int = 0,
        limit: int = 1_000,
    ) -> list[dict[str, Any]]:
        rows = self._conn.execute(
            "SELECT * FROM run_events WHERE run_id = ? AND id > ? ORDER BY id LIMIT ?",
            (run_id, after_id, limit),
        ).fetchall()
        return [_decode_row(row) for row in rows]  # type: ignore[misc]

    def close(self) -> None:
        if self._close_callback is not None:
            self._close_callback(self._conn)
            self._close_callback = None
        elif self._owns_connection:
            self._conn.close()
        if self._direct_lock is not None:
            self._direct_lock.release()
            self._direct_lock = None


__all__ = ["RunLedger", "StageAlreadyCompleteError"]
