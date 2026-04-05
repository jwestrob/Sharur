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

import sqlite3
import json
import time
import uuid
import threading
import asyncio
from pathlib import Path
from contextlib import contextmanager
from typing import Optional
from enum import Enum

from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import StreamingResponse
from pydantic import BaseModel, Field

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

DB_PATH = Path("sharur_ops.db")
app = FastAPI(title="Sharur Ops", version="0.1.0")

# ---------------------------------------------------------------------------
# Database setup — single connection, WAL mode, serialized writes via lock
# ---------------------------------------------------------------------------

_db_lock = threading.Lock()


def _get_db() -> sqlite3.Connection:
    conn = sqlite3.connect(str(DB_PATH), check_same_thread=False)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA busy_timeout=5000")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.row_factory = sqlite3.Row
    return conn


_conn: sqlite3.Connection = _get_db()


def _init_db():
    with _db_lock:
        _conn.executescript("""
        CREATE TABLE IF NOT EXISTS findings (
            id TEXT PRIMARY KEY,
            agent_id TEXT NOT NULL,
            ts REAL NOT NULL,
            finding_type TEXT NOT NULL,
            domain TEXT NOT NULL,
            summary TEXT NOT NULL,
            evidence TEXT NOT NULL DEFAULT '{}',
            confidence REAL NOT NULL DEFAULT 0.5,
            novelty INTEGER NOT NULL DEFAULT 0,
            parent_finding_id TEXT,
            reasoning TEXT NOT NULL DEFAULT ''
        );
        CREATE INDEX IF NOT EXISTS idx_findings_ts ON findings(ts);
        CREATE INDEX IF NOT EXISTS idx_findings_novelty ON findings(novelty);
        CREATE INDEX IF NOT EXISTS idx_findings_agent ON findings(agent_id);
        CREATE INDEX IF NOT EXISTS idx_findings_type ON findings(finding_type);
        CREATE INDEX IF NOT EXISTS idx_findings_domain ON findings(domain);

        CREATE TABLE IF NOT EXISTS hypotheses (
            id TEXT PRIMARY KEY,
            source_agent_id TEXT NOT NULL,
            source_finding_ids TEXT NOT NULL DEFAULT '[]',
            ts REAL NOT NULL,
            hypothesis TEXT NOT NULL,
            status TEXT NOT NULL DEFAULT 'proposed',
            assigned_agent_id TEXT,
            domains_relevant TEXT NOT NULL DEFAULT '[]',
            evidence_for TEXT NOT NULL DEFAULT '[]',
            evidence_against TEXT NOT NULL DEFAULT '[]',
            resolution_summary TEXT
        );
        CREATE INDEX IF NOT EXISTS idx_hyp_status ON hypotheses(status);

        CREATE TABLE IF NOT EXISTS tasks (
            id TEXT PRIMARY KEY,
            created_by TEXT NOT NULL,
            assigned_to TEXT,
            ts REAL NOT NULL,
            claimed_ts REAL,
            completed_ts REAL,
            status TEXT NOT NULL DEFAULT 'pending',
            priority INTEGER NOT NULL DEFAULT 1,
            task_type TEXT NOT NULL,
            description TEXT NOT NULL,
            params TEXT NOT NULL DEFAULT '{}',
            domain_hint TEXT,
            result_finding_ids TEXT NOT NULL DEFAULT '[]'
        );
        CREATE INDEX IF NOT EXISTS idx_tasks_status ON tasks(status);
        CREATE INDEX IF NOT EXISTS idx_tasks_assigned ON tasks(assigned_to);
        CREATE INDEX IF NOT EXISTS idx_tasks_priority ON tasks(priority);

        CREATE TABLE IF NOT EXISTS coordinator_log (
            id TEXT PRIMARY KEY,
            ts REAL NOT NULL,
            action_type TEXT NOT NULL,
            reasoning TEXT NOT NULL,
            referenced_finding_ids TEXT NOT NULL DEFAULT '[]',
            referenced_hypothesis_ids TEXT NOT NULL DEFAULT '[]',
            decisions_made TEXT NOT NULL DEFAULT '{}'
        );
        CREATE INDEX IF NOT EXISTS idx_coordlog_ts ON coordinator_log(ts);
        """)
        _conn.commit()


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
                try:
                    q.put_nowait(payload)
                except asyncio.QueueFull:
                    pass  # Slow consumer drops events


_bus = EventBus()

# ---------------------------------------------------------------------------
# Pydantic models for coordination records only
# ---------------------------------------------------------------------------

class FindingType(str, Enum):
    gene = "gene"
    neighborhood = "neighborhood"
    cassette = "cassette"
    domain_architecture = "domain_architecture"
    phylogenetic = "phylogenetic"
    observation = "observation"


class DomainName(str, Enum):
    omnitrophota = "omnitrophota"
    dpann = "dpann"
    bathyarchaeia = "bathyarchaeia"
    hinthialibacterota = "hinthialibacterota"
    cross_domain = "cross_domain"


class FindingIn(BaseModel):
    agent_id: str
    finding_type: FindingType
    domain: DomainName
    summary: str
    evidence: dict = Field(default_factory=dict)
    confidence: float = Field(default=0.5, ge=0.0, le=1.0)
    novelty: int = Field(default=0, ge=0, le=3)
    parent_finding_id: Optional[str] = None
    reasoning: str = ""


class HypothesisIn(BaseModel):
    source_agent_id: str
    source_finding_ids: list[str] = Field(default_factory=list)
    hypothesis: str
    domains_relevant: list[str] = Field(default_factory=list)


class HypothesisUpdate(BaseModel):
    status: str  # proposed, investigating, supported, refuted, inconclusive
    assigned_agent_id: Optional[str] = None
    evidence_for: Optional[list[str]] = None
    evidence_against: Optional[list[str]] = None
    resolution_summary: Optional[str] = None


class TaskIn(BaseModel):
    created_by: str = "coordinator"
    task_type: str  # investigate, validate_hypothesis, cross_domain_search, follow_up, survey
    description: str
    params: dict = Field(default_factory=dict)
    priority: int = Field(default=1, ge=0, le=3)
    domain_hint: Optional[str] = None
    assigned_to: Optional[str] = None


class TaskUpdate(BaseModel):
    status: str  # pending, claimed, in_progress, complete, failed
    assigned_to: Optional[str] = None
    result_finding_ids: Optional[list[str]] = None


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
    for key in ("evidence", "source_finding_ids", "domains_relevant",
                "evidence_for", "evidence_against", "params",
                "result_finding_ids", "referenced_finding_ids",
                "referenced_hypothesis_ids", "decisions_made"):
        if key in d and isinstance(d[key], str):
            try:
                d[key] = json.loads(d[key])
            except (json.JSONDecodeError, TypeError):
                pass
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
            (fid, f.agent_id, now, f.finding_type.value, f.domain.value, f.summary,
             json.dumps(f.evidence), f.confidence, f.novelty, f.parent_finding_id, f.reasoning),
        )
        _conn.commit()
    _bus.publish("finding", {"id": fid, "agent_id": f.agent_id,
                             "novelty": f.novelty, "summary": f.summary})
    return {"id": fid, "ts": now}


@app.get("/findings")
def list_findings(
    since: float = 0,
    min_novelty: int = 0,
    finding_type: Optional[str] = None,
    domain: Optional[str] = None,
    agent_id: Optional[str] = None,
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
    rows = _conn.execute(
        f"SELECT * FROM findings WHERE {where} ORDER BY ts DESC LIMIT ?", params
    )
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
            (hid, h.source_agent_id, json.dumps(h.source_finding_ids), now,
             h.hypothesis, json.dumps(h.domains_relevant)),
        )
        _conn.commit()
    _bus.publish("hypothesis", {"id": hid, "hypothesis": h.hypothesis})
    return {"id": hid, "ts": now}


@app.get("/hypotheses")
def list_hypotheses(
    status: Optional[str] = None,
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
    rows = _conn.execute(
        f"SELECT * FROM hypotheses{where} ORDER BY ts DESC LIMIT ?", params
    )
    return _rows_to_list(rows)


@app.patch("/hypotheses/{hyp_id}")
def update_hypothesis(hyp_id: str, u: HypothesisUpdate):
    sets = ["status = ?"]
    params: list = [u.status]
    if u.assigned_agent_id is not None:
        sets.append("assigned_agent_id = ?")
        params.append(u.assigned_agent_id)
    if u.evidence_for is not None:
        sets.append("evidence_for = ?")
        params.append(json.dumps(u.evidence_for))
    if u.evidence_against is not None:
        sets.append("evidence_against = ?")
        params.append(json.dumps(u.evidence_against))
    if u.resolution_summary is not None:
        sets.append("resolution_summary = ?")
        params.append(u.resolution_summary)
    params.append(hyp_id)
    with _db_lock:
        cur = _conn.execute(
            f"UPDATE hypotheses SET {', '.join(sets)} WHERE id = ?", params
        )
        _conn.commit()
    if cur.rowcount == 0:
        raise HTTPException(404, "Hypothesis not found")
    return {"updated": hyp_id}


# ---------------------------------------------------------------------------
# TASKS endpoints
# ---------------------------------------------------------------------------

@app.post("/tasks", status_code=201)
def create_task(t: TaskIn):
    tid = str(uuid.uuid4())
    now = time.time()
    with _db_lock:
        _conn.execute(
            """INSERT INTO tasks (id, created_by, assigned_to, ts, status, priority,
               task_type, description, params, domain_hint)
               VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            (tid, t.created_by, t.assigned_to, now,
             "claimed" if t.assigned_to else "pending",
             t.priority, t.task_type, t.description,
             json.dumps(t.params), t.domain_hint),
        )
        _conn.commit()
    _bus.publish("task", {"id": tid, "description": t.description, "priority": t.priority})
    return {"id": tid, "ts": now}


@app.get("/tasks")
def list_tasks(
    status: Optional[str] = None,
    assigned_to: Optional[str] = None,
    unassigned: bool = False,
    limit: int = Query(default=50, le=200),
):
    clauses = []
    params: list = []
    if status:
        clauses.append("status = ?")
        params.append(status)
    if assigned_to:
        clauses.append("assigned_to = ?")
        params.append(assigned_to)
    if unassigned:
        clauses.append("assigned_to IS NULL AND status = 'pending'")
    where = (" WHERE " + " AND ".join(clauses)) if clauses else ""
    params.append(limit)
    rows = _conn.execute(
        f"SELECT * FROM tasks{where} ORDER BY priority DESC, ts ASC LIMIT ?", params
    )
    return _rows_to_list(rows)


@app.patch("/tasks/{task_id}")
def update_task(task_id: str, u: TaskUpdate):
    sets = ["status = ?"]
    params: list = [u.status]
    now = time.time()
    if u.status == "claimed" and u.assigned_to:
        sets.append("assigned_to = ?")
        params.append(u.assigned_to)
        sets.append("claimed_ts = ?")
        params.append(now)
    if u.status in ("complete", "failed"):
        sets.append("completed_ts = ?")
        params.append(now)
    if u.result_finding_ids is not None:
        sets.append("result_finding_ids = ?")
        params.append(json.dumps(u.result_finding_ids))
    params.append(task_id)
    with _db_lock:
        cur = _conn.execute(
            f"UPDATE tasks SET {', '.join(sets)} WHERE id = ?", params
        )
        _conn.commit()
    if cur.rowcount == 0:
        raise HTTPException(404, "Task not found")
    return {"updated": task_id}


@app.post("/tasks/{task_id}/claim")
def claim_task(task_id: str, agent_id: str):
    """Atomic claim: only succeeds if task is still pending/unassigned."""
    now = time.time()
    with _db_lock:
        cur = _conn.execute(
            """UPDATE tasks SET status = 'claimed', assigned_to = ?, claimed_ts = ?
               WHERE id = ? AND status = 'pending'""",
            (agent_id, now, task_id),
        )
        _conn.commit()
    if cur.rowcount == 0:
        raise HTTPException(409, "Task already claimed or not found")
    return {"claimed": task_id, "agent_id": agent_id}


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
            (eid, now, entry.action_type, entry.reasoning,
             json.dumps(entry.referenced_finding_ids),
             json.dumps(entry.referenced_hypothesis_ids),
             json.dumps(entry.decisions_made)),
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
                    yield f": keepalive\n\n"  # Prevents connection timeout
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
    for table in ("findings", "hypotheses", "tasks", "coordinator_log"):
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

    return {
        "counts": counts,
        "findings_by_novelty": {r["novelty"]: r["n"] for r in novelty_dist},
        "tasks_by_status": {r["status"]: r["n"] for r in task_status},
        "hypotheses_by_status": {r["status"]: r["n"] for r in hyp_status},
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
