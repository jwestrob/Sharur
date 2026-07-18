# ruff: noqa: RUF013
"""
Sharur Ops Store — direct SQLite backend for agent coordination.

Usage:
    from sharur.ops.store import OpsStore
    ops = OpsStore("sharur_ops.db", agent_id="my_agent")

    fid = ops.finding(finding_type="observation", domain="cross_domain",
                      summary="Found something", confidence=0.8, novelty=2)
    results = ops.recent_findings(min_novelty=1)

Same API as SharurOps (HTTP client) but hits SQLite directly.
No server needed.
"""

import contextlib
import json
import sqlite3
import threading
import time
import uuid
from pathlib import Path

from sharur.ops.ledger import RunLedger
from sharur.ops.schema import (
    DEFAULT_LEASE_SECONDS,
    OPS_SCHEMA,
    ensure_ops_schema,
)


_SCHEMA = OPS_SCHEMA


def _row_to_dict(row: sqlite3.Row) -> dict:
    d = dict(row)
    for k in (
        "evidence",
        "source_finding_ids",
        "domains_relevant",
        "evidence_for",
        "evidence_against",
        "params",
        "result_finding_ids",
        "dependency_ids",
        "referenced_finding_ids",
        "referenced_hypothesis_ids",
        "decisions_made",
    ):
        if k in d and isinstance(d[k], str):
            with contextlib.suppress(json.JSONDecodeError, TypeError):
                d[k] = json.loads(d[k])
    return d


class OpsStore:
    """Direct SQLite ops store. Same interface as the HTTP client."""

    def __init__(self, db_path: str = "sharur_ops.db", agent_id: str = "unnamed"):
        self.db_path = Path(db_path)
        self.agent_id = agent_id
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._lock = threading.RLock()
        self._conn = sqlite3.connect(str(self.db_path), check_same_thread=False)
        self._conn.execute("PRAGMA journal_mode=WAL")
        self._conn.execute("PRAGMA busy_timeout=5000")
        self._conn.execute("PRAGMA synchronous=NORMAL")
        self._conn.row_factory = sqlite3.Row
        with self._lock:
            ensure_ops_schema(self._conn)
        self.ledger = RunLedger(self.db_path)

    # ------------------------------------------------------------------
    # Findings
    # ------------------------------------------------------------------

    def finding(
        self,
        finding_type: str,
        domain: str,
        summary: str,
        evidence: dict = None,
        confidence: float = 0.5,
        novelty: int = 0,
        parent_finding_id: str = None,
        reasoning: str = "",
    ) -> str:
        fid = str(uuid.uuid4())
        ts = time.time()
        with self._lock:
            self._conn.execute(
                """INSERT INTO findings
                   (id, agent_id, ts, finding_type, domain, summary,
                    evidence, confidence, novelty, parent_finding_id, reasoning)
                   VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                (
                    fid,
                    self.agent_id,
                    ts,
                    finding_type,
                    domain,
                    summary,
                    json.dumps(evidence or {}),
                    confidence,
                    novelty,
                    parent_finding_id,
                    reasoning,
                ),
            )
            self._conn.commit()
        return fid

    def recent_findings(
        self,
        since: float = 0,
        min_novelty: int = 0,
        finding_type: str = None,
        domain: str = None,
        agent_id: str = None,
        limit: int = 50,
    ) -> list[dict]:
        sql = "SELECT * FROM findings WHERE ts > ? AND novelty >= ?"
        params: list = [since, min_novelty]
        if finding_type:
            sql += " AND finding_type = ?"
            params.append(finding_type)
        if domain:
            sql += " AND domain = ?"
            params.append(domain)
        if agent_id:
            sql += " AND agent_id = ?"
            params.append(agent_id)
        sql += " ORDER BY ts DESC LIMIT ?"
        params.append(limit)
        with self._lock:
            rows = self._conn.execute(sql, params).fetchall()
        return [_row_to_dict(r) for r in rows]

    def get_finding(self, finding_id: str) -> dict:
        with self._lock:
            row = self._conn.execute(
                "SELECT * FROM findings WHERE id = ?",
                (finding_id,),
            ).fetchone()
        if not row:
            raise KeyError(f"Finding {finding_id} not found")
        return _row_to_dict(row)

    def search_findings(self, text: str, limit: int = 20) -> list[dict]:
        pattern = f"%{text}%"
        with self._lock:
            rows = self._conn.execute(
                """SELECT * FROM findings
                   WHERE summary LIKE ? OR reasoning LIKE ? OR evidence LIKE ?
                   ORDER BY novelty DESC, ts DESC LIMIT ?""",
                (pattern, pattern, pattern, limit),
            ).fetchall()
        return [_row_to_dict(r) for r in rows]

    # ------------------------------------------------------------------
    # Hypotheses
    # ------------------------------------------------------------------

    def hypothesis(
        self,
        hypothesis: str,
        source_finding_ids: list[str] = None,
        domains_relevant: list[str] = None,
    ) -> str:
        hid = str(uuid.uuid4())
        ts = time.time()
        with self._lock:
            self._conn.execute(
                """INSERT INTO hypotheses
                   (id, source_agent_id, source_finding_ids, ts, hypothesis, domains_relevant)
                   VALUES (?, ?, ?, ?, ?, ?)""",
                (
                    hid,
                    self.agent_id,
                    json.dumps(source_finding_ids or []),
                    ts,
                    hypothesis,
                    json.dumps(domains_relevant or []),
                ),
            )
            self._conn.commit()
        return hid

    def list_hypotheses(
        self,
        *,
        status: str = None,
        unassigned: bool = False,
        limit: int = 50,
    ) -> list[dict]:
        clauses = []
        params: list = []
        if status:
            clauses.append("status = ?")
            params.append(status)
        if unassigned:
            clauses.append("assigned_agent_id IS NULL")
        where = (" WHERE " + " AND ".join(clauses)) if clauses else ""
        params.append(limit)
        with self._lock:
            rows = self._conn.execute(
                f"SELECT * FROM hypotheses{where} ORDER BY ts DESC LIMIT ?",
                params,
            ).fetchall()
        return [_row_to_dict(r) for r in rows]

    def open_hypotheses(self, unassigned: bool = True) -> list[dict]:
        return self.list_hypotheses(
            status="proposed",
            unassigned=unassigned,
        )

    def get_hypothesis(self, hypothesis_id: str) -> dict:
        with self._lock:
            row = self._conn.execute(
                "SELECT * FROM hypotheses WHERE id = ?",
                (hypothesis_id,),
            ).fetchone()
        if row is None:
            raise KeyError(f"Hypothesis {hypothesis_id} not found")
        return _row_to_dict(row)

    def update_hypothesis(self, hyp_id: str, **kwargs) -> dict:
        allowed = {
            "status",
            "assigned_agent_id",
            "domains_relevant",
            "evidence_for",
            "evidence_against",
            "resolution_summary",
        }
        unknown = sorted(set(kwargs) - allowed)
        if unknown:
            raise ValueError(f"Unsupported hypothesis fields: {', '.join(unknown)}")
        if not kwargs:
            raise ValueError("No hypothesis updates supplied")

        with self._lock:
            existing_row = self._conn.execute(
                "SELECT * FROM hypotheses WHERE id = ?",
                (hyp_id,),
            ).fetchone()
            if existing_row is None:
                raise KeyError(f"Hypothesis {hyp_id} not found")

            existing = _row_to_dict(existing_row)
            sets = []
            params = []
            for key, value in kwargs.items():
                updated_value = value
                if key in ("evidence_for", "evidence_against"):
                    combined = list(existing.get(key) or [])
                    combined.extend(updated_value or [])
                    updated_value = list(dict.fromkeys(combined))
                if isinstance(updated_value, (list, dict)):
                    updated_value = json.dumps(updated_value)
                sets.append(f"{key} = ?")
                params.append(updated_value)
            params.append(hyp_id)

            cur = self._conn.execute(
                f"UPDATE hypotheses SET {', '.join(sets)} WHERE id = ?", params
            )
            self._conn.commit()
            if cur.rowcount == 0:
                raise KeyError(f"Hypothesis {hyp_id} not found")
            updated = self._conn.execute(
                "SELECT * FROM hypotheses WHERE id = ?",
                (hyp_id,),
            ).fetchone()
        return _row_to_dict(updated)

    # ------------------------------------------------------------------
    # Tasks
    # ------------------------------------------------------------------

    def create_task(
        self,
        task_type: str,
        description: str,
        params: dict = None,
        priority: int = 1,
        domain_hint: str = None,
        assigned_to: str = None,
        *,
        run_id: str = None,
        idempotency_key: str = None,
        depends_on: list[str] = None,
        max_attempts: int = 3,
        lease_seconds: int = DEFAULT_LEASE_SECONDS,
    ) -> str:
        if max_attempts < 1:
            raise ValueError("max_attempts must be positive")
        if lease_seconds < 1:
            raise ValueError("lease_seconds must be positive")
        if idempotency_key is not None and not idempotency_key.strip():
            raise ValueError("idempotency_key must be non-empty when provided")
        dependencies = list(dict.fromkeys(depends_on or []))
        params_json = json.dumps(params or {}, sort_keys=True)
        dependencies_json = json.dumps(dependencies)
        ts = time.time()
        status = "claimed" if assigned_to else "pending"
        claimed_ts = ts if assigned_to else None
        attempt_count = 1 if assigned_to else 0
        lease_expires_ts = ts + lease_seconds if assigned_to else None
        with self._lock:
            self._conn.execute("BEGIN IMMEDIATE")
            try:
                if idempotency_key:
                    existing = self._conn.execute(
                        "SELECT * FROM tasks WHERE idempotency_key = ?",
                        (idempotency_key,),
                    ).fetchone()
                    if existing is not None:
                        immutable = {
                            "task_type": task_type,
                            "description": description,
                            "params": params_json,
                            "priority": priority,
                            "domain_hint": domain_hint,
                            "run_id": run_id,
                            "dependency_ids": dependencies_json,
                            "max_attempts": max_attempts,
                            "lease_seconds": lease_seconds,
                        }
                        conflicts = [
                            key for key, value in immutable.items() if existing[key] != value
                        ]
                        if conflicts:
                            raise ValueError(
                                "Task idempotency key already exists with a "
                                f"different payload: {', '.join(conflicts)}"
                            )
                        self._conn.commit()
                        return str(existing["id"])

                self._validate_dependencies_locked(dependencies)
                if assigned_to and not self._dependencies_complete_locked(dependencies):
                    raise ValueError("Cannot preassign a task until all dependencies are complete")
                tid = str(uuid.uuid4())
                self._conn.execute(
                    """INSERT INTO tasks
                       (id, created_by, assigned_to, ts, claimed_ts, status, priority,
                        task_type, description, params, domain_hint, run_id,
                        idempotency_key, dependency_ids, attempt_count, max_attempts,
                        lease_seconds, lease_expires_ts, heartbeat_ts)
                       VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                    (
                        tid,
                        self.agent_id,
                        assigned_to,
                        ts,
                        claimed_ts,
                        status,
                        priority,
                        task_type,
                        description,
                        params_json,
                        domain_hint,
                        run_id,
                        idempotency_key,
                        dependencies_json,
                        attempt_count,
                        max_attempts,
                        lease_seconds,
                        lease_expires_ts,
                        claimed_ts,
                    ),
                )
                self._conn.commit()
            except Exception:
                self._conn.rollback()
                raise
        return tid

    def my_tasks(self, status: str = None) -> list[dict]:
        return self.list_tasks(status=status, assigned_to=self.agent_id)

    def list_tasks(
        self,
        *,
        status: str = None,
        assigned_to: str = None,
        unassigned: bool = False,
        limit: int = 50,
    ) -> list[dict]:
        if unassigned:
            return self.available_tasks()[:limit]
        clauses = []
        params: list = []
        if status:
            clauses.append("status = ?")
            params.append(status)
        if assigned_to:
            clauses.append("assigned_to = ?")
            params.append(assigned_to)
        where = (" WHERE " + " AND ".join(clauses)) if clauses else ""
        params.append(limit)
        with self._lock:
            rows = self._conn.execute(
                f"SELECT * FROM tasks{where} ORDER BY priority DESC, ts ASC LIMIT ?",
                params,
            ).fetchall()
        return [_row_to_dict(row) for row in rows]

    def available_tasks(self) -> list[dict]:
        with self._lock:
            self._conn.execute("BEGIN IMMEDIATE")
            try:
                now = time.time()
                self._recover_expired_locked(now)
                self._promote_retries_locked(now)
                rows = self._conn.execute(
                    """
                    SELECT * FROM tasks
                    WHERE status = 'pending' AND assigned_to IS NULL
                      AND attempt_count < max_attempts
                    ORDER BY priority DESC, ts ASC
                    """
                ).fetchall()
                available = [
                    row
                    for row in rows
                    if self._dependencies_complete_locked(json.loads(row["dependency_ids"] or "[]"))
                ]
                self._conn.commit()
            except Exception:
                self._conn.rollback()
                raise
        return [_row_to_dict(row) for row in available]

    def claim_task(
        self,
        task_id: str,
        *,
        lease_seconds: int = None,
    ) -> dict:
        with self._lock:
            self._conn.execute("BEGIN IMMEDIATE")
            try:
                now = time.time()
                self._recover_expired_locked(now)
                self._promote_retries_locked(now)
                row = self._conn.execute(
                    "SELECT * FROM tasks WHERE id = ?",
                    (task_id,),
                ).fetchone()
                if row is None:
                    raise ValueError(f"Task {task_id} not found")
                dependencies = json.loads(row["dependency_ids"] or "[]")
                if not self._dependencies_complete_locked(dependencies):
                    raise ValueError(f"Task {task_id} has incomplete dependencies")
                duration = int(lease_seconds or row["lease_seconds"])
                if duration < 1:
                    raise ValueError("lease_seconds must be positive")
                cursor = self._conn.execute(
                    """
                    UPDATE tasks
                    SET status = 'claimed', assigned_to = ?, claimed_ts = ?,
                        heartbeat_ts = ?, lease_expires_ts = ?,
                        attempt_count = attempt_count + 1,
                        retry_after_ts = NULL
                    WHERE id = ? AND status = 'pending' AND assigned_to IS NULL
                      AND attempt_count < max_attempts
                    """,
                    (
                        self.agent_id,
                        now,
                        now,
                        now + duration,
                        task_id,
                    ),
                )
                if cursor.rowcount == 0:
                    raise ValueError(f"Task {task_id} already claimed, exhausted, or not pending")
                self._conn.commit()
            except Exception:
                self._conn.rollback()
                raise
        return self._get_task(task_id)

    def heartbeat_task(
        self,
        task_id: str,
        *,
        lease_seconds: int = None,
        in_progress: bool = True,
    ) -> dict:
        with self._lock:
            now = time.time()
            row = self._conn.execute(
                "SELECT lease_seconds FROM tasks WHERE id = ?",
                (task_id,),
            ).fetchone()
            if row is None:
                raise ValueError(f"Task {task_id} not found")
            duration = int(lease_seconds or row["lease_seconds"])
            status = "in_progress" if in_progress else "claimed"
            cursor = self._conn.execute(
                """
                UPDATE tasks
                SET status = ?, heartbeat_ts = ?, lease_expires_ts = ?
                WHERE id = ? AND assigned_to = ?
                  AND status IN ('claimed', 'in_progress')
                  AND (lease_expires_ts IS NULL OR lease_expires_ts > ?)
                """,
                (
                    status,
                    now,
                    now + duration,
                    task_id,
                    self.agent_id,
                    now,
                ),
            )
            self._conn.commit()
            if cursor.rowcount == 0:
                raise ValueError(f"Task {task_id} has no live lease owned by {self.agent_id}")
        return self._get_task(task_id)

    def complete_task(self, task_id: str, result_finding_ids: list[str] = None) -> dict:
        with self._lock:
            now = time.time()
            cur = self._conn.execute(
                "UPDATE tasks SET status = 'complete', completed_ts = ?, "
                "result_finding_ids = ?, heartbeat_ts = ?, lease_expires_ts = NULL "
                "WHERE id = ? AND assigned_to = ? "
                "AND status IN ('claimed', 'in_progress') "
                "AND (lease_expires_ts IS NULL OR lease_expires_ts > ?)",
                (
                    now,
                    json.dumps(result_finding_ids or []),
                    now,
                    task_id,
                    self.agent_id,
                    now,
                ),
            )
            self._conn.commit()
            if cur.rowcount == 0:
                raise ValueError(f"Task {task_id} is not active and assigned to {self.agent_id}")
            row = self._conn.execute(
                "SELECT * FROM tasks WHERE id = ?",
                (task_id,),
            ).fetchone()
        return _row_to_dict(row)

    def fail_task(
        self,
        task_id: str,
        *,
        error: str = None,
        retryable: bool = False,
        retry_delay_seconds: int = 0,
    ) -> dict:
        with self._lock:
            self._conn.execute("BEGIN IMMEDIATE")
            try:
                now = time.time()
                row = self._conn.execute(
                    "SELECT * FROM tasks WHERE id = ?",
                    (task_id,),
                ).fetchone()
                if (
                    row is None
                    or row["assigned_to"] != self.agent_id
                    or row["status"] not in {"claimed", "in_progress"}
                    or (row["lease_expires_ts"] is not None and row["lease_expires_ts"] <= now)
                ):
                    raise ValueError(
                        f"Task {task_id} is not active and assigned to {self.agent_id}"
                    )
                can_retry = retryable and row["attempt_count"] < row["max_attempts"]
                status = "retry_wait" if can_retry else "failed"
                completed_ts = None if can_retry else now
                retry_after = now + max(0, retry_delay_seconds) if can_retry else None
                self._conn.execute(
                    """
                    UPDATE tasks
                    SET status = ?, assigned_to = NULL, completed_ts = ?,
                        heartbeat_ts = ?, lease_expires_ts = NULL,
                        retry_after_ts = ?, last_error = ?
                    WHERE id = ?
                    """,
                    (
                        status,
                        completed_ts,
                        now,
                        retry_after,
                        error,
                        task_id,
                    ),
                )
                self._conn.commit()
            except Exception:
                self._conn.rollback()
                raise
        return self._get_task(task_id)

    def recover_expired_tasks(self, *, now: float = None) -> dict:
        with self._lock:
            self._conn.execute("BEGIN IMMEDIATE")
            try:
                recovered, exhausted = self._recover_expired_locked(
                    time.time() if now is None else now
                )
                self._conn.commit()
            except Exception:
                self._conn.rollback()
                raise
        return {"recovered": recovered, "exhausted": exhausted}

    def _recover_expired_locked(self, now: float) -> tuple[list[str], list[str]]:
        rows = self._conn.execute(
            """
            SELECT id, attempt_count, max_attempts
            FROM tasks
            WHERE status IN ('claimed', 'in_progress')
              AND lease_expires_ts IS NOT NULL
              AND lease_expires_ts <= ?
            """,
            (now,),
        ).fetchall()
        recovered: list[str] = []
        exhausted: list[str] = []
        for row in rows:
            if row["attempt_count"] < row["max_attempts"]:
                status = "pending"
                completed_ts = None
                recovered.append(row["id"])
            else:
                status = "failed"
                completed_ts = now
                exhausted.append(row["id"])
            self._conn.execute(
                """
                UPDATE tasks
                SET status = ?, assigned_to = NULL, completed_ts = ?,
                    lease_expires_ts = NULL, retry_after_ts = NULL,
                    last_error = COALESCE(last_error, 'worker lease expired')
                WHERE id = ?
                """,
                (status, completed_ts, row["id"]),
            )
        return recovered, exhausted

    def _promote_retries_locked(self, now: float) -> None:
        self._conn.execute(
            """
            UPDATE tasks
            SET status = 'pending', retry_after_ts = NULL
            WHERE status = 'retry_wait'
              AND (retry_after_ts IS NULL OR retry_after_ts <= ?)
              AND attempt_count < max_attempts
            """,
            (now,),
        )

    def _validate_dependencies_locked(self, dependency_ids: list[str]) -> None:
        if not dependency_ids:
            return
        placeholders = ",".join("?" for _ in dependency_ids)
        found = {
            row[0]
            for row in self._conn.execute(
                f"SELECT id FROM tasks WHERE id IN ({placeholders})",
                dependency_ids,
            )
        }
        missing = [task_id for task_id in dependency_ids if task_id not in found]
        if missing:
            raise ValueError(f"Unknown task dependencies: {', '.join(missing)}")

    def _dependencies_complete_locked(self, dependency_ids: list[str]) -> bool:
        if not dependency_ids:
            return True
        placeholders = ",".join("?" for _ in dependency_ids)
        rows = self._conn.execute(
            f"SELECT id, status FROM tasks WHERE id IN ({placeholders})",
            dependency_ids,
        ).fetchall()
        return len(rows) == len(dependency_ids) and all(row["status"] == "complete" for row in rows)

    def _get_task(self, task_id: str) -> dict:
        row = self._conn.execute(
            "SELECT * FROM tasks WHERE id = ?",
            (task_id,),
        ).fetchone()
        if row is None:
            raise KeyError(f"Task {task_id} not found")
        return _row_to_dict(row)

    def get_task(self, task_id: str) -> dict:
        with self._lock:
            return self._get_task(task_id)

    # ------------------------------------------------------------------
    # Runs
    # ------------------------------------------------------------------

    def create_run(
        self,
        run_type: str,
        dataset_path: str | Path,
        *,
        config: dict = None,
        idempotency_key: str = None,
        parent_run_id: str = None,
    ) -> str:
        return self.ledger.create_run(
            run_type,
            dataset_path,
            created_by=self.agent_id,
            config=config,
            idempotency_key=idempotency_key,
            parent_run_id=parent_run_id,
        )

    def start_run(self, run_id: str) -> dict:
        return self.ledger.start_run(run_id)

    def submit_run(self, run_id: str) -> dict:
        return self.ledger.submit_run(run_id)

    def wait_for_scheduler(self, run_id: str) -> dict:
        return self.ledger.wait_for_scheduler(run_id)

    def heartbeat_run(self, run_id: str) -> dict:
        return self.ledger.heartbeat_run(run_id)

    def complete_run(self, run_id: str, result: dict = None) -> dict:
        return self.ledger.complete_run(run_id, result=result)

    def fail_run(self, run_id: str, error: str) -> dict:
        return self.ledger.fail_run(run_id, error)

    def get_run(self, run_id: str) -> dict:
        return self.ledger.get_run(run_id)

    def list_runs(self, **filters) -> list[dict]:
        return self.ledger.list_runs(**filters)

    def list_run_stages(self, run_id: str, *, stage_id: str = None) -> list[dict]:
        return self.ledger.list_stages(run_id, stage_id=stage_id)

    def run_events(self, run_id: str) -> list[dict]:
        self.ledger.get_run(run_id)
        return self.ledger.events(run_id)

    def recover_stale_runs(
        self,
        *,
        stale_after_seconds: int = 300,
        now: float = None,
    ) -> dict:
        return self.ledger.recover_stale_runs(
            stale_after_seconds=stale_after_seconds,
            now=now,
        )

    # ------------------------------------------------------------------
    # Coordinator log
    # ------------------------------------------------------------------

    def log(
        self,
        action_type: str,
        reasoning: str,
        referenced_finding_ids: list[str] = None,
        referenced_hypothesis_ids: list[str] = None,
        decisions_made: dict = None,
    ) -> str:
        lid = str(uuid.uuid4())
        ts = time.time()
        with self._lock:
            self._conn.execute(
                """INSERT INTO coordinator_log
                   (id, ts, action_type, reasoning,
                    referenced_finding_ids, referenced_hypothesis_ids, decisions_made)
                   VALUES (?, ?, ?, ?, ?, ?, ?)""",
                (
                    lid,
                    ts,
                    action_type,
                    reasoning,
                    json.dumps(referenced_finding_ids or []),
                    json.dumps(referenced_hypothesis_ids or []),
                    json.dumps(decisions_made or {}),
                ),
            )
            self._conn.commit()
        return lid

    def recent_log(self, limit: int = 20, since: float = 0) -> list[dict]:
        rows = self._conn.execute(
            "SELECT * FROM coordinator_log WHERE ts > ? ORDER BY ts DESC LIMIT ?",
            (since, limit),
        ).fetchall()
        return [_row_to_dict(r) for r in rows]

    def get_log_entry(self, entry_id: str) -> dict:
        with self._lock:
            row = self._conn.execute(
                "SELECT * FROM coordinator_log WHERE id = ?",
                (entry_id,),
            ).fetchone()
        if row is None:
            raise KeyError(f"Coordinator log entry {entry_id} not found")
        return _row_to_dict(row)

    # ------------------------------------------------------------------
    # Stats
    # ------------------------------------------------------------------

    def stats(self) -> dict:
        counts = {}
        for table in (
            "findings",
            "hypotheses",
            "tasks",
            "coordinator_log",
            "runs",
            "run_events",
            "run_stages",
        ):
            counts[table] = self._conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]

        novelty = {}
        for row in self._conn.execute(
            "SELECT novelty, COUNT(*) FROM findings GROUP BY novelty"
        ).fetchall():
            novelty[str(row[0])] = row[1]

        task_status = {}
        for row in self._conn.execute(
            "SELECT status, COUNT(*) FROM tasks GROUP BY status"
        ).fetchall():
            task_status[row[0]] = row[1]

        hyp_status = {}
        for row in self._conn.execute(
            "SELECT status, COUNT(*) FROM hypotheses GROUP BY status"
        ).fetchall():
            hyp_status[row[0]] = row[1]

        run_status = {
            row[0]: row[1]
            for row in self._conn.execute(
                "SELECT status, COUNT(*) FROM runs GROUP BY status"
            ).fetchall()
        }
        stage_status = {
            row[0]: row[1]
            for row in self._conn.execute(
                "SELECT status, COUNT(*) FROM run_stages GROUP BY status"
            ).fetchall()
        }

        return {
            "counts": counts,
            "findings_by_novelty": novelty,
            "tasks_by_status": task_status,
            "hypotheses_by_status": hyp_status,
            "runs_by_status": run_status,
            "stages_by_status": stage_status,
        }

    def close(self):
        self.ledger.close()
        self._conn.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def __repr__(self):
        return f"OpsStore(db='{self.db_path}', agent_id='{self.agent_id}')"
