"""Shared, migration-safe SQLite schema for Sharur operational state."""

from __future__ import annotations

from typing import TYPE_CHECKING


if TYPE_CHECKING:
    import sqlite3


OPS_SCHEMA_VERSION = 2
DEFAULT_LEASE_SECONDS = 900

OPS_SCHEMA = """
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
    result_finding_ids TEXT NOT NULL DEFAULT '[]',
    run_id TEXT,
    idempotency_key TEXT,
    dependency_ids TEXT NOT NULL DEFAULT '[]',
    attempt_count INTEGER NOT NULL DEFAULT 0,
    max_attempts INTEGER NOT NULL DEFAULT 3,
    lease_seconds INTEGER NOT NULL DEFAULT 900,
    lease_expires_ts REAL,
    heartbeat_ts REAL,
    retry_after_ts REAL,
    last_error TEXT
);
CREATE INDEX IF NOT EXISTS idx_tasks_status ON tasks(status);
CREATE INDEX IF NOT EXISTS idx_tasks_assigned ON tasks(assigned_to);
CREATE INDEX IF NOT EXISTS idx_tasks_priority ON tasks(priority);
CREATE INDEX IF NOT EXISTS idx_tasks_run ON tasks(run_id);
CREATE INDEX IF NOT EXISTS idx_tasks_lease ON tasks(status, lease_expires_ts);
CREATE UNIQUE INDEX IF NOT EXISTS idx_tasks_idempotency
    ON tasks(idempotency_key) WHERE idempotency_key IS NOT NULL;

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

CREATE TABLE IF NOT EXISTS runs (
    id TEXT PRIMARY KEY,
    idempotency_key TEXT,
    run_type TEXT NOT NULL,
    dataset_path TEXT NOT NULL,
    created_by TEXT NOT NULL,
    status TEXT NOT NULL DEFAULT 'pending',
    created_ts REAL NOT NULL,
    started_ts REAL,
    heartbeat_ts REAL,
    completed_ts REAL,
    parent_run_id TEXT,
    current_stage TEXT,
    config TEXT NOT NULL DEFAULT '{}',
    result TEXT NOT NULL DEFAULT '{}',
    error TEXT
);
CREATE UNIQUE INDEX IF NOT EXISTS idx_runs_idempotency
    ON runs(idempotency_key) WHERE idempotency_key IS NOT NULL;
CREATE INDEX IF NOT EXISTS idx_runs_dataset ON runs(dataset_path, run_type, created_ts);
CREATE INDEX IF NOT EXISTS idx_runs_status ON runs(status);

CREATE TABLE IF NOT EXISTS run_events (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    run_id TEXT NOT NULL,
    ts REAL NOT NULL,
    event_type TEXT NOT NULL,
    stage_id TEXT,
    attempt INTEGER,
    payload TEXT NOT NULL DEFAULT '{}'
);
CREATE INDEX IF NOT EXISTS idx_run_events_run ON run_events(run_id, id);

CREATE TABLE IF NOT EXISTS run_stages (
    run_id TEXT NOT NULL,
    stage_id TEXT NOT NULL,
    attempt INTEGER NOT NULL,
    signature TEXT NOT NULL,
    status TEXT NOT NULL,
    started_ts REAL,
    heartbeat_ts REAL,
    completed_ts REAL,
    command TEXT NOT NULL DEFAULT '[]',
    inputs TEXT NOT NULL DEFAULT '{}',
    outputs TEXT NOT NULL DEFAULT '{}',
    resource_profile TEXT NOT NULL DEFAULT '{}',
    reused_from_run_id TEXT,
    reused_from_attempt INTEGER,
    error TEXT,
    PRIMARY KEY (run_id, stage_id, attempt)
);
CREATE INDEX IF NOT EXISTS idx_run_stages_signature
    ON run_stages(stage_id, signature, status, completed_ts);

CREATE TABLE IF NOT EXISTS ops_schema_meta (
    version INTEGER PRIMARY KEY,
    applied_ts REAL NOT NULL
);
"""


_TASK_MIGRATION_COLUMNS = {
    "run_id": "TEXT",
    "idempotency_key": "TEXT",
    "dependency_ids": "TEXT NOT NULL DEFAULT '[]'",
    "attempt_count": "INTEGER NOT NULL DEFAULT 0",
    "max_attempts": "INTEGER NOT NULL DEFAULT 3",
    "lease_seconds": f"INTEGER NOT NULL DEFAULT {DEFAULT_LEASE_SECONDS}",
    "lease_expires_ts": "REAL",
    "heartbeat_ts": "REAL",
    "retry_after_ts": "REAL",
    "last_error": "TEXT",
}


def _existing_columns(conn: sqlite3.Connection, table: str) -> set[str]:
    return {row[1] for row in conn.execute(f"PRAGMA table_info({table})")}


def ensure_ops_schema(conn: sqlite3.Connection) -> None:
    """Create current tables and migrate legacy task queues in place."""
    # Legacy `tasks` must gain its columns before indexes referencing those
    # columns can be created, so execute the table definitions in two passes.
    statements = [statement.strip() for statement in OPS_SCHEMA.split(";") if statement.strip()]
    deferred = [
        statement
        for statement in statements
        if "idx_tasks_run" in statement
        or "idx_tasks_lease" in statement
        or "idx_tasks_idempotency" in statement
    ]
    for statement in statements:
        if statement not in deferred:
            conn.execute(statement)

    columns = _existing_columns(conn, "tasks")
    for name, definition in _TASK_MIGRATION_COLUMNS.items():
        if name not in columns:
            conn.execute(f"ALTER TABLE tasks ADD COLUMN {name} {definition}")
    for statement in deferred:
        conn.execute(statement)

    # Legacy active tasks gain finite leases. On the next queue read, genuinely
    # abandoned work is recoverable instead of remaining claimed forever.
    conn.execute(
        """
        UPDATE tasks
        SET attempt_count = CASE WHEN attempt_count < 1 THEN 1 ELSE attempt_count END,
            heartbeat_ts = COALESCE(heartbeat_ts, claimed_ts, ts),
            lease_expires_ts = COALESCE(
                lease_expires_ts,
                COALESCE(claimed_ts, ts) + lease_seconds
            )
        WHERE status IN ('claimed', 'in_progress')
        """
    )
    conn.execute(
        "INSERT OR REPLACE INTO ops_schema_meta(version, applied_ts) "
        "VALUES (?, unixepoch('subsec'))",
        (OPS_SCHEMA_VERSION,),
    )
    conn.commit()


__all__ = [
    "DEFAULT_LEASE_SECONDS",
    "OPS_SCHEMA",
    "OPS_SCHEMA_VERSION",
    "ensure_ops_schema",
]
