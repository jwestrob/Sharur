"""
Sharur Ops Server — Coordination layer for multi-agent metagenomic discovery.

Run:  uvicorn sharur.ops.server:app --host 0.0.0.0 --port 8811
Or:   python -m sharur.ops.server

Agents POST findings/hypotheses/tasks. Orchestrator queries them.
All writes serialized through a single SQLite connection in WAL mode.

Important boundary: this server stores coordination records only. Dataset-local
scientific records remain the canonical source of truth in
`survey/findings.jsonl`, `exploration/findings.jsonl`, and
`exploration/hypotheses.json`.
"""

import asyncio
import contextlib
import json
import os
import sqlite3
import threading
import time
import uuid
from pathlib import Path
from typing import Literal

from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import StreamingResponse
from pydantic import BaseModel, Field

from sharur import __version__
from sharur.ops.ledger import RunLedger
from sharur.ops.schema import DEFAULT_LEASE_SECONDS, ensure_ops_schema
from sharur.ops.store import OpsStore


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

# DB path is configurable so containers can point it at a persistent volume.
# A bare relative default resolves against the process CWD, which loses
# coordination state across container restarts — hence the env override.
DB_PATH = Path(os.environ.get("SHARUR_OPS_DB_PATH", "sharur_ops.db"))
app = FastAPI(title="Sharur Ops", version=__version__)

# ---------------------------------------------------------------------------
# Database setup — single connection, WAL mode, serialized writes via lock
# ---------------------------------------------------------------------------

_db_lock = threading.Lock()


def _get_db() -> sqlite3.Connection:
    DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(DB_PATH), check_same_thread=False)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA busy_timeout=5000")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.row_factory = sqlite3.Row
    return conn


_conn: sqlite3.Connection = _get_db()


def _init_db():
    with _db_lock:
        ensure_ops_schema(_conn)


# ---------------------------------------------------------------------------
# SSE event bus — orchestrator subscribes, agents publish
# ---------------------------------------------------------------------------


class EventBus:
    """Simple in-process pub/sub for SSE streaming."""

    def __init__(self):
        self._subscribers: list[asyncio.Queue] = []
        self._lock = threading.Lock()

    def subscribe(self) -> asyncio.Queue:
        q: asyncio.Queue = asyncio.Queue()
        with self._lock:
            self._subscribers.append(q)
        return q

    def unsubscribe(self, q: asyncio.Queue):
        with self._lock:
            self._subscribers = [s for s in self._subscribers if s is not q]

    def publish(self, event_type: str, data: dict):
        payload = json.dumps({"type": event_type, **data})
        with self._lock:
            for q in self._subscribers:
                with contextlib.suppress(asyncio.QueueFull):
                    q.put_nowait(payload)


_bus = EventBus()

# ---------------------------------------------------------------------------
# Pydantic models for coordination records only
# ---------------------------------------------------------------------------


class FindingIn(BaseModel):
    agent_id: str
    finding_type: str = Field(min_length=1)
    domain: str = Field(min_length=1)
    summary: str
    evidence: dict = Field(default_factory=dict)
    confidence: float = Field(default=0.5, ge=0.0, le=1.0)
    novelty: int = Field(default=0, ge=0, le=3)
    parent_finding_id: str | None = None
    reasoning: str = ""


class HypothesisIn(BaseModel):
    source_agent_id: str
    source_finding_ids: list[str] = Field(default_factory=list)
    hypothesis: str
    domains_relevant: list[str] = Field(default_factory=list)


class HypothesisUpdate(BaseModel):
    status: Literal["proposed", "investigating", "supported", "refuted", "inconclusive"] | None = (
        None
    )
    assigned_agent_id: str | None = None
    evidence_for: list[str] | None = None
    evidence_against: list[str] | None = None
    resolution_summary: str | None = None


class TaskIn(BaseModel):
    created_by: str = "coordinator"
    task_type: str  # investigate, validate_hypothesis, cross_domain_search, follow_up, survey
    description: str
    params: dict = Field(default_factory=dict)
    priority: int = Field(default=1, ge=0, le=3)
    domain_hint: str | None = None
    assigned_to: str | None = None
    run_id: str | None = None
    idempotency_key: str | None = None
    depends_on: list[str] = Field(default_factory=list)
    max_attempts: int = Field(default=3, ge=1)
    lease_seconds: int = Field(default=DEFAULT_LEASE_SECONDS, ge=1)


class TaskUpdate(BaseModel):
    status: Literal[
        "pending",
        "claimed",
        "in_progress",
        "retry_wait",
        "complete",
        "failed",
    ]
    assigned_to: str | None = None
    agent_id: str | None = None
    result_finding_ids: list[str] | None = None
    lease_seconds: int | None = Field(default=None, ge=1)
    error: str | None = None
    retryable: bool = False
    retry_delay_seconds: int = Field(default=0, ge=0)


class RunIn(BaseModel):
    created_by: str = "operator"
    run_type: str = Field(min_length=1)
    dataset_path: str = Field(min_length=1)
    config: dict = Field(default_factory=dict)
    idempotency_key: str | None = None
    parent_run_id: str | None = None


class RunUpdate(BaseModel):
    status: Literal["complete", "failed"]
    result: dict = Field(default_factory=dict)
    error: str | None = None


class CoordinatorLogIn(BaseModel):
    action_type: str  # synthesis, task_assignment, escalation, hypothesis_generated, checkpoint
    reasoning: str
    referenced_finding_ids: list[str] = Field(default_factory=list)
    referenced_hypothesis_ids: list[str] = Field(default_factory=list)
    decisions_made: dict = Field(default_factory=dict)


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------


def _row_to_dict(row: sqlite3.Row) -> dict:
    d = dict(row)
    # Parse JSON fields back to objects for the response
    for key in (
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
        "config",
        "result",
        "payload",
    ):
        if key in d and isinstance(d[key], str):
            with contextlib.suppress(json.JSONDecodeError, TypeError):
                d[key] = json.loads(d[key])
    return d


def _rows_to_list(cursor) -> list[dict]:
    return [_row_to_dict(r) for r in cursor.fetchall()]


# ---------------------------------------------------------------------------
# FINDINGS endpoints
# ---------------------------------------------------------------------------


@app.post("/findings", status_code=201)
def create_finding(f: FindingIn):
    fid = str(uuid.uuid4())
    now = time.time()
    with _db_lock:
        _conn.execute(
            """INSERT INTO findings (id, agent_id, ts, finding_type, domain, summary,
               evidence, confidence, novelty, parent_finding_id, reasoning)
               VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            (
                fid,
                f.agent_id,
                now,
                f.finding_type,
                f.domain,
                f.summary,
                json.dumps(f.evidence),
                f.confidence,
                f.novelty,
                f.parent_finding_id,
                f.reasoning,
            ),
        )
        _conn.commit()
    _bus.publish(
        "finding", {"id": fid, "agent_id": f.agent_id, "novelty": f.novelty, "summary": f.summary}
    )
    return {"id": fid, "ts": now}


@app.get("/findings")
def list_findings(
    since: float = 0,
    min_novelty: int = 0,
    finding_type: str | None = None,
    domain: str | None = None,
    agent_id: str | None = None,
    limit: int = Query(default=50, le=500),
):
    clauses = ["ts > ?"]
    params: list = [since]
    if min_novelty > 0:
        clauses.append("novelty >= ?")
        params.append(min_novelty)
    if finding_type:
        clauses.append("finding_type = ?")
        params.append(finding_type)
    if domain:
        clauses.append("domain = ?")
        params.append(domain)
    if agent_id:
        clauses.append("agent_id = ?")
        params.append(agent_id)
    where = " AND ".join(clauses)
    params.append(limit)
    rows = _conn.execute(f"SELECT * FROM findings WHERE {where} ORDER BY ts DESC LIMIT ?", params)
    return _rows_to_list(rows)


@app.get("/findings/{finding_id}")
def get_finding(finding_id: str):
    row = _conn.execute("SELECT * FROM findings WHERE id = ?", (finding_id,)).fetchone()
    if not row:
        raise HTTPException(404, "Finding not found")
    return _row_to_dict(row)


@app.get("/findings/search/{text}")
def search_findings(text: str, limit: int = Query(default=20, le=100)):
    """Full-text search across summary, reasoning, and evidence."""
    pattern = f"%{text}%"
    rows = _conn.execute(
        """SELECT * FROM findings
           WHERE summary LIKE ? OR reasoning LIKE ? OR evidence LIKE ?
           ORDER BY novelty DESC, ts DESC LIMIT ?""",
        (pattern, pattern, pattern, limit),
    )
    return _rows_to_list(rows)


# ---------------------------------------------------------------------------
# HYPOTHESES endpoints
# ---------------------------------------------------------------------------


@app.post("/hypotheses", status_code=201)
def create_hypothesis(h: HypothesisIn):
    hid = str(uuid.uuid4())
    now = time.time()
    with _db_lock:
        _conn.execute(
            """INSERT INTO hypotheses (id, source_agent_id, source_finding_ids, ts,
               hypothesis, status, domains_relevant)
               VALUES (?, ?, ?, ?, ?, 'proposed', ?)""",
            (
                hid,
                h.source_agent_id,
                json.dumps(h.source_finding_ids),
                now,
                h.hypothesis,
                json.dumps(h.domains_relevant),
            ),
        )
        _conn.commit()
    _bus.publish("hypothesis", {"id": hid, "hypothesis": h.hypothesis})
    return {"id": hid, "ts": now}


@app.get("/hypotheses")
def list_hypotheses(
    status: str | None = None,
    unassigned: bool = False,
    limit: int = Query(default=50, le=200),
):
    clauses = []
    params: list = []
    if status:
        clauses.append("status = ?")
        params.append(status)
    if unassigned:
        clauses.append("assigned_agent_id IS NULL")
    where = (" WHERE " + " AND ".join(clauses)) if clauses else ""
    params.append(limit)
    rows = _conn.execute(f"SELECT * FROM hypotheses{where} ORDER BY ts DESC LIMIT ?", params)
    return _rows_to_list(rows)


@app.patch("/hypotheses/{hyp_id}")
def update_hypothesis(hyp_id: str, u: HypothesisUpdate):
    with _db_lock:
        existing_row = _conn.execute(
            "SELECT * FROM hypotheses WHERE id = ?",
            (hyp_id,),
        ).fetchone()
        if existing_row is None:
            raise HTTPException(404, "Hypothesis not found")
        existing = _row_to_dict(existing_row)

        sets = []
        params: list = []
        if u.status is not None:
            sets.append("status = ?")
            params.append(u.status)
        if u.assigned_agent_id is not None:
            sets.append("assigned_agent_id = ?")
            params.append(u.assigned_agent_id)
        if u.evidence_for is not None:
            evidence_for = list(
                dict.fromkeys([*(existing.get("evidence_for") or []), *u.evidence_for])
            )
            sets.append("evidence_for = ?")
            params.append(json.dumps(evidence_for))
        if u.evidence_against is not None:
            evidence_against = list(
                dict.fromkeys([*(existing.get("evidence_against") or []), *u.evidence_against])
            )
            sets.append("evidence_against = ?")
            params.append(json.dumps(evidence_against))
        if u.resolution_summary is not None:
            sets.append("resolution_summary = ?")
            params.append(u.resolution_summary)
        if not sets:
            raise HTTPException(400, "No hypothesis updates supplied")

        params.append(hyp_id)
        cur = _conn.execute(f"UPDATE hypotheses SET {', '.join(sets)} WHERE id = ?", params)
        _conn.commit()
    if cur.rowcount == 0:
        raise HTTPException(404, "Hypothesis not found")
    return {"updated": hyp_id}


# ---------------------------------------------------------------------------
# TASKS endpoints
# ---------------------------------------------------------------------------


@app.post("/tasks", status_code=201)
def create_task(t: TaskIn):
    store = OpsStore(DB_PATH, agent_id=t.created_by)
    try:
        tid = store.create_task(
            t.task_type,
            t.description,
            params=t.params,
            priority=t.priority,
            domain_hint=t.domain_hint,
            assigned_to=t.assigned_to,
            run_id=t.run_id,
            idempotency_key=t.idempotency_key,
            depends_on=t.depends_on,
            max_attempts=t.max_attempts,
            lease_seconds=t.lease_seconds,
        )
        task = store._get_task(tid)
    except ValueError as exc:
        raise HTTPException(409, str(exc)) from exc
    finally:
        store.close()
    _bus.publish("task", {"id": tid, "description": t.description, "priority": t.priority})
    return {"id": tid, "ts": task["ts"], "task": task}


@app.get("/tasks")
def list_tasks(
    status: str | None = None,
    assigned_to: str | None = None,
    unassigned: bool = False,
    limit: int = Query(default=50, le=200),
):
    if unassigned:
        store = OpsStore(DB_PATH, agent_id="server")
        try:
            return store.available_tasks()[:limit]
        finally:
            store.close()

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
    rows = _conn.execute(
        f"SELECT * FROM tasks{where} ORDER BY priority DESC, ts ASC LIMIT ?", params
    )
    return _rows_to_list(rows)


@app.patch("/tasks/{task_id}")
def update_task(task_id: str, u: TaskUpdate):
    if not u.agent_id:
        raise HTTPException(400, "agent_id is required for task updates")
    store = OpsStore(DB_PATH, agent_id=u.agent_id)
    try:
        if u.status == "complete":
            return store.complete_task(task_id, u.result_finding_ids)
        if u.status == "failed":
            return store.fail_task(
                task_id,
                error=u.error,
                retryable=u.retryable,
                retry_delay_seconds=u.retry_delay_seconds,
            )
        if u.status == "in_progress":
            return store.heartbeat_task(
                task_id,
                lease_seconds=u.lease_seconds,
                in_progress=True,
            )
        raise HTTPException(
            400,
            "Use the claim endpoint for ownership; only in_progress, complete, "
            "and failed are valid worker updates.",
        )
    except (KeyError, ValueError) as exc:
        raise HTTPException(409, str(exc)) from exc
    finally:
        store.close()


@app.post("/tasks/{task_id}/claim")
def claim_task(
    task_id: str,
    agent_id: str,
    lease_seconds: int | None = Query(default=None, ge=1),
):
    """Atomically claim a ready task and start a finite worker lease."""
    store = OpsStore(DB_PATH, agent_id=agent_id)
    try:
        return store.claim_task(task_id, lease_seconds=lease_seconds)
    except ValueError as exc:
        raise HTTPException(409, str(exc)) from exc
    finally:
        store.close()


@app.post("/tasks/{task_id}/heartbeat")
def heartbeat_task(
    task_id: str,
    agent_id: str,
    lease_seconds: int | None = Query(default=None, ge=1),
):
    """Extend an owned task lease and mark the task in progress."""
    store = OpsStore(DB_PATH, agent_id=agent_id)
    try:
        return store.heartbeat_task(task_id, lease_seconds=lease_seconds)
    except ValueError as exc:
        raise HTTPException(409, str(exc)) from exc
    finally:
        store.close()


@app.post("/tasks/recover")
def recover_tasks(now: float | None = None):
    """Recover expired worker leases; primarily for the coordinator."""
    store = OpsStore(DB_PATH, agent_id="coordinator")
    try:
        return store.recover_expired_tasks(now=now)
    finally:
        store.close()


# ---------------------------------------------------------------------------
# RUN LEDGER endpoints
# ---------------------------------------------------------------------------


@app.post("/runs", status_code=201)
def create_run(run: RunIn):
    ledger = RunLedger(DB_PATH)
    try:
        run_id = ledger.create_run(
            run.run_type,
            run.dataset_path,
            created_by=run.created_by,
            config=run.config,
            idempotency_key=run.idempotency_key,
            parent_run_id=run.parent_run_id,
        )
        return ledger.get_run(run_id)
    except ValueError as exc:
        raise HTTPException(409, str(exc)) from exc
    finally:
        ledger.close()


@app.get("/runs")
def list_runs(
    dataset_path: str | None = None,
    run_type: str | None = None,
    status: str | None = None,
    limit: int = Query(default=50, le=200),
):
    ledger = RunLedger(DB_PATH)
    try:
        return ledger.list_runs(
            dataset_path=dataset_path,
            run_type=run_type,
            status=status,
            limit=limit,
        )
    finally:
        ledger.close()


@app.get("/runs/{run_id}")
def get_run(run_id: str):
    ledger = RunLedger(DB_PATH)
    try:
        return ledger.get_run(run_id)
    except KeyError as exc:
        raise HTTPException(404, str(exc)) from exc
    finally:
        ledger.close()


@app.post("/runs/{run_id}/start")
def start_run(run_id: str):
    ledger = RunLedger(DB_PATH)
    try:
        return ledger.start_run(run_id)
    except ValueError as exc:
        raise HTTPException(409, str(exc)) from exc
    finally:
        ledger.close()


@app.post("/runs/{run_id}/submit")
def submit_run(run_id: str):
    ledger = RunLedger(DB_PATH)
    try:
        return ledger.submit_run(run_id)
    except ValueError as exc:
        raise HTTPException(409, str(exc)) from exc
    finally:
        ledger.close()


@app.post("/runs/{run_id}/heartbeat")
def heartbeat_run(run_id: str):
    ledger = RunLedger(DB_PATH)
    try:
        return ledger.heartbeat_run(run_id)
    except ValueError as exc:
        raise HTTPException(409, str(exc)) from exc
    finally:
        ledger.close()


@app.patch("/runs/{run_id}")
def update_run(run_id: str, update: RunUpdate):
    ledger = RunLedger(DB_PATH)
    try:
        if update.status == "complete":
            return ledger.complete_run(run_id, result=update.result)
        if not update.error:
            raise HTTPException(400, "error is required when failing a run")
        return ledger.fail_run(run_id, update.error)
    except ValueError as exc:
        raise HTTPException(409, str(exc)) from exc
    finally:
        ledger.close()


@app.get("/runs/{run_id}/events")
def get_run_events(run_id: str):
    ledger = RunLedger(DB_PATH)
    try:
        ledger.get_run(run_id)
        return ledger.events(run_id)
    except KeyError as exc:
        raise HTTPException(404, str(exc)) from exc
    finally:
        ledger.close()


@app.get("/runs/{run_id}/stages")
def get_run_stages(run_id: str, stage_id: str | None = None):
    ledger = RunLedger(DB_PATH)
    try:
        ledger.get_run(run_id)
        return ledger.list_stages(run_id, stage_id=stage_id)
    except KeyError as exc:
        raise HTTPException(404, str(exc)) from exc
    finally:
        ledger.close()


@app.post("/runs/recover")
def recover_runs(
    stale_after_seconds: int = Query(default=300, ge=1),
    now: float | None = None,
):
    ledger = RunLedger(DB_PATH)
    try:
        return ledger.recover_stale_runs(
            stale_after_seconds=stale_after_seconds,
            now=now,
        )
    finally:
        ledger.close()


# ---------------------------------------------------------------------------
# COORDINATOR LOG endpoints
# ---------------------------------------------------------------------------


@app.post("/coordinator/log", status_code=201)
def create_log_entry(entry: CoordinatorLogIn):
    eid = str(uuid.uuid4())
    now = time.time()
    with _db_lock:
        _conn.execute(
            """INSERT INTO coordinator_log (id, ts, action_type, reasoning,
               referenced_finding_ids, referenced_hypothesis_ids, decisions_made)
               VALUES (?, ?, ?, ?, ?, ?, ?)""",
            (
                eid,
                now,
                entry.action_type,
                entry.reasoning,
                json.dumps(entry.referenced_finding_ids),
                json.dumps(entry.referenced_hypothesis_ids),
                json.dumps(entry.decisions_made),
            ),
        )
        _conn.commit()
    return {"id": eid, "ts": now}


@app.get("/coordinator/log")
def get_log(limit: int = Query(default=20, le=100), since: float = 0):
    rows = _conn.execute(
        "SELECT * FROM coordinator_log WHERE ts > ? ORDER BY ts DESC LIMIT ?",
        (since, limit),
    )
    return _rows_to_list(rows)


# ---------------------------------------------------------------------------
# SSE stream — orchestrator subscribes to real-time events
# ---------------------------------------------------------------------------


@app.get("/stream")
async def event_stream():
    """SSE endpoint. Orchestrator connects here to get push notifications
    for new findings, hypotheses, and tasks."""
    q = _bus.subscribe()

    async def generate():
        try:
            while True:
                try:
                    data = await asyncio.wait_for(q.get(), timeout=30)
                    yield f"data: {data}\n\n"
                except asyncio.TimeoutError:
                    yield ": keepalive\n\n"  # Prevents connection timeout
        except asyncio.CancelledError:
            pass
        finally:
            _bus.unsubscribe(q)

    return StreamingResponse(generate(), media_type="text/event-stream")


# ---------------------------------------------------------------------------
# Stats / health
# ---------------------------------------------------------------------------


@app.get("/stats")
def stats():
    """Quick overview for the orchestrator or for you."""
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
        row = _conn.execute(f"SELECT COUNT(*) as n FROM {table}").fetchone()
        counts[table] = row["n"]

    # Breakdown
    novelty_dist = _conn.execute(
        "SELECT novelty, COUNT(*) as n FROM findings GROUP BY novelty ORDER BY novelty"
    ).fetchall()
    task_status = _conn.execute(
        "SELECT status, COUNT(*) as n FROM tasks GROUP BY status"
    ).fetchall()
    hyp_status = _conn.execute(
        "SELECT status, COUNT(*) as n FROM hypotheses GROUP BY status"
    ).fetchall()
    run_status = _conn.execute("SELECT status, COUNT(*) as n FROM runs GROUP BY status").fetchall()
    stage_status = _conn.execute(
        "SELECT status, COUNT(*) as n FROM run_stages GROUP BY status"
    ).fetchall()

    return {
        "counts": counts,
        "findings_by_novelty": {r["novelty"]: r["n"] for r in novelty_dist},
        "tasks_by_status": {r["status"]: r["n"] for r in task_status},
        "hypotheses_by_status": {r["status"]: r["n"] for r in hyp_status},
        "runs_by_status": {r["status"]: r["n"] for r in run_status},
        "stages_by_status": {r["status"]: r["n"] for r in stage_status},
    }


@app.get("/health")
def health():
    return {"status": "ok", "db": str(DB_PATH), "ts": time.time()}


# ---------------------------------------------------------------------------
# Run directly
# ---------------------------------------------------------------------------

_init_db()

if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=8811)
