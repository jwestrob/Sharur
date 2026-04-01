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

import json
import sqlite3
import threading
import time
import uuid
from pathlib import Path
from typing import Optional


_SCHEMA = """
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
"""


def _row_to_dict(row: sqlite3.Row) -> dict:
    d = dict(row)
    for k in ("evidence", "source_finding_ids", "domains_relevant",
              "evidence_for", "evidence_against", "params",
              "result_finding_ids", "referenced_finding_ids",
              "referenced_hypothesis_ids", "decisions_made"):
        if k in d and isinstance(d[k], str):
            try:
                d[k] = json.loads(d[k])
            except (json.JSONDecodeError, TypeError):
                pass
    return d


class OpsStore:
    """Direct SQLite ops store. Same interface as the HTTP client."""

    def __init__(self, db_path: str = "sharur_ops.db", agent_id: str = "unnamed"):
        self.db_path = Path(db_path)
        self.agent_id = agent_id
        self._lock = threading.Lock()
        self._conn = sqlite3.connect(str(self.db_path), check_same_thread=False)
        self._conn.execute("PRAGMA journal_mode=WAL")
        self._conn.execute("PRAGMA busy_timeout=5000")
        self._conn.execute("PRAGMA synchronous=NORMAL")
        self._conn.row_factory = sqlite3.Row
        with self._lock:
            self._conn.executescript(_SCHEMA)
            self._conn.commit()

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
                (fid, self.agent_id, ts, finding_type, domain, summary,
                 json.dumps(evidence or {}), confidence, novelty,
                 parent_finding_id, reasoning),
            )
            self._conn.commit()
        return fid

    def recent_findings(
        self,
        since: float = 0,
        min_novelty: int = 0,
        finding_type: str = None,
        domain: str = None,
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
        sql += " ORDER BY ts DESC LIMIT ?"
        params.append(limit)
        rows = self._conn.execute(sql, params).fetchall()
        return [_row_to_dict(r) for r in rows]

    def get_finding(self, finding_id: str) -> dict:
        row = self._conn.execute(
            "SELECT * FROM findings WHERE id = ?", (finding_id,)
        ).fetchone()
        if not row:
            raise KeyError(f"Finding {finding_id} not found")
        return _row_to_dict(row)

    def search_findings(self, text: str, limit: int = 20) -> list[dict]:
        pattern = f"%{text}%"
        rows = self._conn.execute(
            """SELECT * FROM findings
               WHERE summary LIKE ? OR reasoning LIKE ? OR evidence LIKE ?
               ORDER BY ts DESC LIMIT ?""",
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
                (hid, self.agent_id, json.dumps(source_finding_ids or []),
                 ts, hypothesis, json.dumps(domains_relevant or [])),
            )
            self._conn.commit()
        return hid

    def open_hypotheses(self, unassigned: bool = True) -> list[dict]:
        if unassigned:
            rows = self._conn.execute(
                "SELECT * FROM hypotheses WHERE status = 'proposed' AND assigned_agent_id IS NULL"
            ).fetchall()
        else:
            rows = self._conn.execute(
                "SELECT * FROM hypotheses WHERE status = 'proposed'"
            ).fetchall()
        return [_row_to_dict(r) for r in rows]

    def update_hypothesis(self, hyp_id: str, **kwargs) -> dict:
        sets = []
        params = []
        for k, v in kwargs.items():
            if k in ("evidence_for", "evidence_against"):
                # Append to existing list
                existing = self._conn.execute(
                    f"SELECT {k} FROM hypotheses WHERE id = ?", (hyp_id,)
                ).fetchone()
                if existing:
                    current = json.loads(existing[0])
                    current.extend(v)
                    v = current
                sets.append(f"{k} = ?")
                params.append(json.dumps(v))
            elif isinstance(v, (list, dict)):
                sets.append(f"{k} = ?")
                params.append(json.dumps(v))
            else:
                sets.append(f"{k} = ?")
                params.append(v)
        params.append(hyp_id)
        with self._lock:
            self._conn.execute(
                f"UPDATE hypotheses SET {', '.join(sets)} WHERE id = ?", params
            )
            self._conn.commit()
        return _row_to_dict(
            self._conn.execute("SELECT * FROM hypotheses WHERE id = ?", (hyp_id,)).fetchone()
        )

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
    ) -> str:
        tid = str(uuid.uuid4())
        ts = time.time()
        with self._lock:
            self._conn.execute(
                """INSERT INTO tasks
                   (id, created_by, assigned_to, ts, status, priority,
                    task_type, description, params, domain_hint)
                   VALUES (?, ?, ?, ?, 'pending', ?, ?, ?, ?, ?)""",
                (tid, self.agent_id, assigned_to, ts, priority,
                 task_type, description, json.dumps(params or {}), domain_hint),
            )
            self._conn.commit()
        return tid

    def my_tasks(self, status: str = None) -> list[dict]:
        sql = "SELECT * FROM tasks WHERE assigned_to = ?"
        params: list = [self.agent_id]
        if status:
            sql += " AND status = ?"
            params.append(status)
        return [_row_to_dict(r) for r in self._conn.execute(sql, params).fetchall()]

    def available_tasks(self) -> list[dict]:
        rows = self._conn.execute(
            "SELECT * FROM tasks WHERE status = 'pending' AND assigned_to IS NULL ORDER BY priority DESC"
        ).fetchall()
        return [_row_to_dict(r) for r in rows]

    def claim_task(self, task_id: str) -> dict:
        with self._lock:
            cur = self._conn.execute(
                "UPDATE tasks SET status = 'claimed', assigned_to = ?, claimed_ts = ? WHERE id = ? AND status = 'pending'",
                (self.agent_id, time.time(), task_id),
            )
            self._conn.commit()
            if cur.rowcount == 0:
                raise ValueError(f"Task {task_id} already claimed or not found")
        return _row_to_dict(
            self._conn.execute("SELECT * FROM tasks WHERE id = ?", (task_id,)).fetchone()
        )

    def complete_task(self, task_id: str, result_finding_ids: list[str] = None) -> dict:
        with self._lock:
            self._conn.execute(
                "UPDATE tasks SET status = 'complete', completed_ts = ?, result_finding_ids = ? WHERE id = ?",
                (time.time(), json.dumps(result_finding_ids or []), task_id),
            )
            self._conn.commit()
        return _row_to_dict(
            self._conn.execute("SELECT * FROM tasks WHERE id = ?", (task_id,)).fetchone()
        )

    def fail_task(self, task_id: str) -> dict:
        with self._lock:
            self._conn.execute(
                "UPDATE tasks SET status = 'failed', completed_ts = ? WHERE id = ?",
                (time.time(), task_id),
            )
            self._conn.commit()
        return _row_to_dict(
            self._conn.execute("SELECT * FROM tasks WHERE id = ?", (task_id,)).fetchone()
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
                (lid, ts, action_type, reasoning,
                 json.dumps(referenced_finding_ids or []),
                 json.dumps(referenced_hypothesis_ids or []),
                 json.dumps(decisions_made or {})),
            )
            self._conn.commit()
        return lid

    def recent_log(self, limit: int = 20, since: float = 0) -> list[dict]:
        rows = self._conn.execute(
            "SELECT * FROM coordinator_log WHERE ts > ? ORDER BY ts DESC LIMIT ?",
            (since, limit),
        ).fetchall()
        return [_row_to_dict(r) for r in rows]

    # ------------------------------------------------------------------
    # Stats
    # ------------------------------------------------------------------

    def stats(self) -> dict:
        counts = {}
        for table in ("findings", "hypotheses", "tasks", "coordinator_log"):
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

        return {
            "counts": counts,
            "findings_by_novelty": novelty,
            "tasks_by_status": task_status,
            "hypotheses_by_status": hyp_status,
        }

    def close(self):
        self._conn.close()

    def __repr__(self):
        return f"OpsStore(db='{self.db_path}', agent_id='{self.agent_id}')"
