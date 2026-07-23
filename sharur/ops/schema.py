"""Shared, migration-safe SQLite schema for Sharur operational state.

Version 3 turns the original single-process coordination database into a
durable HTTP control-plane store.  The migration is deliberately additive:
legacy JSON columns remain readable while normalized relationship tables,
attempt fencing, campaigns, durable events, agent identities, and artifact
metadata become authoritative for new writes.
"""

from __future__ import annotations

import contextlib
import json
import time
from typing import TYPE_CHECKING


if TYPE_CHECKING:
    import sqlite3


OPS_SCHEMA_VERSION = 3
DEFAULT_LEASE_SECONDS = 900


_TABLES_SQL = """
CREATE TABLE IF NOT EXISTS campaigns (
    id TEXT PRIMARY KEY,
    name TEXT NOT NULL,
    description TEXT NOT NULL DEFAULT '',
    dataset_path TEXT,
    created_by TEXT NOT NULL,
    status TEXT NOT NULL DEFAULT 'active'
        CHECK(status IN ('active', 'paused', 'complete', 'failed', 'archived')),
    created_ts REAL NOT NULL,
    completed_ts REAL,
    metadata TEXT NOT NULL DEFAULT '{}' CHECK(json_valid(metadata)),
    idempotency_key TEXT
);

CREATE TABLE IF NOT EXISTS agents (
    id TEXT PRIMARY KEY,
    role TEXT NOT NULL
        CHECK(role IN ('reader', 'worker', 'coordinator', 'operator')),
    token_hash TEXT,
    token_hint TEXT,
    capabilities TEXT NOT NULL DEFAULT '[]' CHECK(json_valid(capabilities)),
    max_concurrent_tasks INTEGER NOT NULL DEFAULT 1
        CHECK(max_concurrent_tasks >= 0),
    capacity_cpu_slots INTEGER NOT NULL DEFAULT 1 CHECK(capacity_cpu_slots >= 0),
    capacity_memory_mb INTEGER NOT NULL DEFAULT 0 CHECK(capacity_memory_mb >= 0),
    capacity_accelerator_slots INTEGER NOT NULL DEFAULT 0
        CHECK(capacity_accelerator_slots >= 0),
    host TEXT,
    status TEXT NOT NULL DEFAULT 'active'
        CHECK(status IN ('active', 'draining', 'offline', 'revoked')),
    registered_ts REAL NOT NULL,
    heartbeat_ts REAL,
    metadata TEXT NOT NULL DEFAULT '{}' CHECK(json_valid(metadata))
);

CREATE TABLE IF NOT EXISTS findings (
    id TEXT PRIMARY KEY,
    agent_id TEXT NOT NULL,
    ts REAL NOT NULL,
    finding_type TEXT NOT NULL,
    domain TEXT NOT NULL,
    summary TEXT NOT NULL,
    evidence TEXT NOT NULL DEFAULT '{}' CHECK(json_valid(evidence)),
    confidence REAL NOT NULL DEFAULT 0.5 CHECK(confidence BETWEEN 0.0 AND 1.0),
    novelty INTEGER NOT NULL DEFAULT 0 CHECK(novelty BETWEEN 0 AND 3),
    parent_finding_id TEXT,
    reasoning TEXT NOT NULL DEFAULT '',
    campaign_id TEXT,
    task_id TEXT,
    idempotency_key TEXT,
    schema_version INTEGER NOT NULL DEFAULT 1,
    validation_status TEXT NOT NULL DEFAULT 'unreviewed',
    created_event_id INTEGER
);

CREATE TABLE IF NOT EXISTS hypotheses (
    id TEXT PRIMARY KEY,
    source_agent_id TEXT NOT NULL,
    source_finding_ids TEXT NOT NULL DEFAULT '[]' CHECK(json_valid(source_finding_ids)),
    ts REAL NOT NULL,
    hypothesis TEXT NOT NULL,
    status TEXT NOT NULL DEFAULT 'proposed'
        CHECK(status IN ('proposed', 'investigating', 'supported', 'refuted', 'inconclusive')),
    assigned_agent_id TEXT,
    domains_relevant TEXT NOT NULL DEFAULT '[]' CHECK(json_valid(domains_relevant)),
    evidence_for TEXT NOT NULL DEFAULT '[]' CHECK(json_valid(evidence_for)),
    evidence_against TEXT NOT NULL DEFAULT '[]' CHECK(json_valid(evidence_against)),
    resolution_summary TEXT,
    campaign_id TEXT,
    idempotency_key TEXT,
    schema_version INTEGER NOT NULL DEFAULT 1,
    created_event_id INTEGER
);

CREATE TABLE IF NOT EXISTS tasks (
    id TEXT PRIMARY KEY,
    created_by TEXT NOT NULL,
    assigned_to TEXT,
    ts REAL NOT NULL,
    claimed_ts REAL,
    completed_ts REAL,
    status TEXT NOT NULL DEFAULT 'pending'
        CHECK(status IN ('pending', 'claimed', 'in_progress', 'retry_wait', 'complete', 'failed', 'cancelled')),
    priority INTEGER NOT NULL DEFAULT 1 CHECK(priority BETWEEN 0 AND 3),
    task_type TEXT NOT NULL,
    description TEXT NOT NULL,
    params TEXT NOT NULL DEFAULT '{}' CHECK(json_valid(params)),
    domain_hint TEXT,
    result_finding_ids TEXT NOT NULL DEFAULT '[]' CHECK(json_valid(result_finding_ids)),
    run_id TEXT,
    idempotency_key TEXT,
    dependency_ids TEXT NOT NULL DEFAULT '[]' CHECK(json_valid(dependency_ids)),
    attempt_count INTEGER NOT NULL DEFAULT 0 CHECK(attempt_count >= 0),
    max_attempts INTEGER NOT NULL DEFAULT 3 CHECK(max_attempts >= 1),
    lease_seconds INTEGER NOT NULL DEFAULT 900 CHECK(lease_seconds >= 1),
    lease_expires_ts REAL,
    heartbeat_ts REAL,
    retry_after_ts REAL,
    last_error TEXT,
    lease_token_hash TEXT,
    reserved_for TEXT,
    campaign_id TEXT,
    required_capabilities TEXT NOT NULL DEFAULT '[]'
        CHECK(json_valid(required_capabilities)),
    resource_request TEXT NOT NULL DEFAULT '{}' CHECK(json_valid(resource_request)),
    schema_version INTEGER NOT NULL DEFAULT 1,
    created_event_id INTEGER
);

CREATE TABLE IF NOT EXISTS coordinator_log (
    id TEXT PRIMARY KEY,
    ts REAL NOT NULL,
    action_type TEXT NOT NULL,
    reasoning TEXT NOT NULL,
    referenced_finding_ids TEXT NOT NULL DEFAULT '[]'
        CHECK(json_valid(referenced_finding_ids)),
    referenced_hypothesis_ids TEXT NOT NULL DEFAULT '[]'
        CHECK(json_valid(referenced_hypothesis_ids)),
    decisions_made TEXT NOT NULL DEFAULT '{}' CHECK(json_valid(decisions_made)),
    actor_agent_id TEXT,
    campaign_id TEXT,
    idempotency_key TEXT,
    created_event_id INTEGER
);

CREATE TABLE IF NOT EXISTS runs (
    id TEXT PRIMARY KEY,
    idempotency_key TEXT,
    run_type TEXT NOT NULL,
    dataset_path TEXT NOT NULL,
    created_by TEXT NOT NULL,
    status TEXT NOT NULL DEFAULT 'pending'
        CHECK(status IN ('pending', 'submitted', 'running', 'complete', 'failed', 'cancelled')),
    created_ts REAL NOT NULL,
    started_ts REAL,
    heartbeat_ts REAL,
    completed_ts REAL,
    parent_run_id TEXT,
    current_stage TEXT,
    config TEXT NOT NULL DEFAULT '{}' CHECK(json_valid(config)),
    result TEXT NOT NULL DEFAULT '{}' CHECK(json_valid(result)),
    error TEXT,
    campaign_id TEXT,
    created_event_id INTEGER
);

CREATE TABLE IF NOT EXISTS run_events (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    run_id TEXT NOT NULL,
    ts REAL NOT NULL,
    event_type TEXT NOT NULL,
    stage_id TEXT,
    attempt INTEGER,
    payload TEXT NOT NULL DEFAULT '{}' CHECK(json_valid(payload)),
    FOREIGN KEY(run_id) REFERENCES runs(id) ON DELETE CASCADE
);

CREATE TABLE IF NOT EXISTS run_stages (
    run_id TEXT NOT NULL,
    stage_id TEXT NOT NULL,
    attempt INTEGER NOT NULL,
    signature TEXT NOT NULL,
    status TEXT NOT NULL
        CHECK(status IN ('running', 'complete', 'failed', 'reused', 'cancelled')),
    started_ts REAL,
    heartbeat_ts REAL,
    completed_ts REAL,
    command TEXT NOT NULL DEFAULT '[]' CHECK(json_valid(command)),
    inputs TEXT NOT NULL DEFAULT '{}' CHECK(json_valid(inputs)),
    outputs TEXT NOT NULL DEFAULT '{}' CHECK(json_valid(outputs)),
    resource_profile TEXT NOT NULL DEFAULT '{}' CHECK(json_valid(resource_profile)),
    reused_from_run_id TEXT,
    reused_from_attempt INTEGER,
    error TEXT,
    owner_agent_id TEXT,
    PRIMARY KEY (run_id, stage_id, attempt),
    FOREIGN KEY(run_id) REFERENCES runs(id) ON DELETE CASCADE
);

CREATE TABLE IF NOT EXISTS task_dependencies (
    task_id TEXT NOT NULL,
    depends_on_task_id TEXT NOT NULL,
    PRIMARY KEY(task_id, depends_on_task_id),
    CHECK(task_id <> depends_on_task_id),
    FOREIGN KEY(task_id) REFERENCES tasks(id) ON DELETE CASCADE,
    FOREIGN KEY(depends_on_task_id) REFERENCES tasks(id) ON DELETE RESTRICT
);

CREATE TABLE IF NOT EXISTS task_result_findings (
    task_id TEXT NOT NULL,
    finding_id TEXT NOT NULL,
    PRIMARY KEY(task_id, finding_id),
    FOREIGN KEY(task_id) REFERENCES tasks(id) ON DELETE CASCADE,
    FOREIGN KEY(finding_id) REFERENCES findings(id) ON DELETE RESTRICT
);

CREATE TABLE IF NOT EXISTS hypothesis_findings (
    hypothesis_id TEXT NOT NULL,
    finding_id TEXT NOT NULL,
    relation TEXT NOT NULL CHECK(relation IN ('source', 'for', 'against')),
    PRIMARY KEY(hypothesis_id, finding_id, relation),
    FOREIGN KEY(hypothesis_id) REFERENCES hypotheses(id) ON DELETE CASCADE,
    FOREIGN KEY(finding_id) REFERENCES findings(id) ON DELETE RESTRICT
);

CREATE TABLE IF NOT EXISTS finding_links (
    finding_id TEXT NOT NULL,
    related_finding_id TEXT NOT NULL,
    relation TEXT NOT NULL,
    created_ts REAL NOT NULL,
    created_by TEXT NOT NULL,
    PRIMARY KEY(finding_id, related_finding_id, relation),
    CHECK(finding_id <> related_finding_id),
    FOREIGN KEY(finding_id) REFERENCES findings(id) ON DELETE CASCADE,
    FOREIGN KEY(related_finding_id) REFERENCES findings(id) ON DELETE CASCADE
);

CREATE TABLE IF NOT EXISTS coordinator_log_findings (
    log_id TEXT NOT NULL,
    finding_id TEXT NOT NULL,
    PRIMARY KEY(log_id, finding_id),
    FOREIGN KEY(log_id) REFERENCES coordinator_log(id) ON DELETE CASCADE,
    FOREIGN KEY(finding_id) REFERENCES findings(id) ON DELETE RESTRICT
);

CREATE TABLE IF NOT EXISTS coordinator_log_hypotheses (
    log_id TEXT NOT NULL,
    hypothesis_id TEXT NOT NULL,
    PRIMARY KEY(log_id, hypothesis_id),
    FOREIGN KEY(log_id) REFERENCES coordinator_log(id) ON DELETE CASCADE,
    FOREIGN KEY(hypothesis_id) REFERENCES hypotheses(id) ON DELETE RESTRICT
);

CREATE TABLE IF NOT EXISTS coordination_events (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    ts REAL NOT NULL,
    event_type TEXT NOT NULL,
    actor_agent_id TEXT,
    campaign_id TEXT,
    entity_type TEXT NOT NULL,
    entity_id TEXT NOT NULL,
    run_id TEXT,
    task_id TEXT,
    payload TEXT NOT NULL DEFAULT '{}' CHECK(json_valid(payload)),
    schema_version INTEGER NOT NULL DEFAULT 1
);

CREATE TABLE IF NOT EXISTS artifacts (
    id TEXT PRIMARY KEY,
    content_hash TEXT NOT NULL UNIQUE,
    uri TEXT NOT NULL,
    size_bytes INTEGER NOT NULL CHECK(size_bytes >= 0),
    media_type TEXT NOT NULL DEFAULT 'application/octet-stream',
    created_by TEXT NOT NULL,
    created_ts REAL NOT NULL,
    campaign_id TEXT,
    task_id TEXT,
    run_id TEXT,
    metadata TEXT NOT NULL DEFAULT '{}' CHECK(json_valid(metadata))
);

CREATE TABLE IF NOT EXISTS finding_artifacts (
    finding_id TEXT NOT NULL,
    artifact_id TEXT NOT NULL,
    relation TEXT NOT NULL DEFAULT 'evidence',
    PRIMARY KEY(finding_id, artifact_id, relation),
    FOREIGN KEY(finding_id) REFERENCES findings(id) ON DELETE CASCADE,
    FOREIGN KEY(artifact_id) REFERENCES artifacts(id) ON DELETE RESTRICT
);

CREATE TABLE IF NOT EXISTS ops_schema_meta (
    version INTEGER PRIMARY KEY,
    applied_ts REAL NOT NULL
);
"""


_INDEXES_SQL = """
CREATE INDEX IF NOT EXISTS idx_campaigns_status
    ON campaigns(status, created_ts);
CREATE UNIQUE INDEX IF NOT EXISTS idx_campaigns_idempotency
    ON campaigns(idempotency_key) WHERE idempotency_key IS NOT NULL;
CREATE UNIQUE INDEX IF NOT EXISTS idx_agents_token_hash
    ON agents(token_hash) WHERE token_hash IS NOT NULL;
CREATE INDEX IF NOT EXISTS idx_agents_status
    ON agents(status, heartbeat_ts);

CREATE INDEX IF NOT EXISTS idx_findings_ts ON findings(ts);
CREATE INDEX IF NOT EXISTS idx_findings_novelty ON findings(novelty);
CREATE INDEX IF NOT EXISTS idx_findings_agent ON findings(agent_id);
CREATE INDEX IF NOT EXISTS idx_findings_type ON findings(finding_type);
CREATE INDEX IF NOT EXISTS idx_findings_domain ON findings(domain);
CREATE INDEX IF NOT EXISTS idx_findings_campaign
    ON findings(campaign_id, ts DESC, id);
CREATE INDEX IF NOT EXISTS idx_findings_task ON findings(task_id);
CREATE UNIQUE INDEX IF NOT EXISTS idx_findings_idempotency
    ON findings(agent_id, idempotency_key) WHERE idempotency_key IS NOT NULL;

CREATE INDEX IF NOT EXISTS idx_hyp_status ON hypotheses(status);
CREATE INDEX IF NOT EXISTS idx_hyp_campaign
    ON hypotheses(campaign_id, status, ts DESC, id);
CREATE UNIQUE INDEX IF NOT EXISTS idx_hyp_idempotency
    ON hypotheses(source_agent_id, idempotency_key)
    WHERE idempotency_key IS NOT NULL;

CREATE INDEX IF NOT EXISTS idx_tasks_status ON tasks(status);
CREATE INDEX IF NOT EXISTS idx_tasks_assigned ON tasks(assigned_to);
CREATE INDEX IF NOT EXISTS idx_tasks_priority ON tasks(priority);
CREATE INDEX IF NOT EXISTS idx_tasks_run ON tasks(run_id);
CREATE INDEX IF NOT EXISTS idx_tasks_lease ON tasks(status, lease_expires_ts);
CREATE INDEX IF NOT EXISTS idx_tasks_queue
    ON tasks(status, priority DESC, ts, id);
CREATE INDEX IF NOT EXISTS idx_tasks_reserved_queue
    ON tasks(reserved_for, status, priority DESC, ts, id);
CREATE INDEX IF NOT EXISTS idx_tasks_campaign
    ON tasks(campaign_id, status, priority DESC, ts, id);
CREATE UNIQUE INDEX IF NOT EXISTS idx_tasks_idempotency
    ON tasks(idempotency_key) WHERE idempotency_key IS NOT NULL;
CREATE INDEX IF NOT EXISTS idx_task_dependencies_parent
    ON task_dependencies(depends_on_task_id, task_id);

CREATE INDEX IF NOT EXISTS idx_coordlog_ts ON coordinator_log(ts);
CREATE INDEX IF NOT EXISTS idx_coordlog_campaign
    ON coordinator_log(campaign_id, ts DESC, id);
CREATE UNIQUE INDEX IF NOT EXISTS idx_coordlog_idempotency
    ON coordinator_log(actor_agent_id, idempotency_key)
    WHERE idempotency_key IS NOT NULL;

CREATE UNIQUE INDEX IF NOT EXISTS idx_runs_idempotency
    ON runs(idempotency_key) WHERE idempotency_key IS NOT NULL;
CREATE INDEX IF NOT EXISTS idx_runs_dataset
    ON runs(dataset_path, run_type, created_ts);
CREATE INDEX IF NOT EXISTS idx_runs_status ON runs(status);
CREATE INDEX IF NOT EXISTS idx_runs_campaign
    ON runs(campaign_id, status, created_ts DESC, id);
CREATE INDEX IF NOT EXISTS idx_run_events_run ON run_events(run_id, id);
CREATE INDEX IF NOT EXISTS idx_run_stages_signature
    ON run_stages(stage_id, signature, status, completed_ts);
CREATE UNIQUE INDEX IF NOT EXISTS idx_one_running_stage_attempt
    ON run_stages(run_id, stage_id) WHERE status = 'running';

CREATE INDEX IF NOT EXISTS idx_coordination_events_campaign
    ON coordination_events(campaign_id, id);
CREATE INDEX IF NOT EXISTS idx_coordination_events_entity
    ON coordination_events(entity_type, entity_id, id);
CREATE INDEX IF NOT EXISTS idx_coordination_events_run
    ON coordination_events(run_id, id);
CREATE INDEX IF NOT EXISTS idx_coordination_events_task
    ON coordination_events(task_id, id);
CREATE INDEX IF NOT EXISTS idx_artifacts_campaign
    ON artifacts(campaign_id, created_ts DESC, id);
"""


_VALIDATION_TRIGGERS = (
    """
    CREATE TRIGGER IF NOT EXISTS trg_tasks_v3_insert
    BEFORE INSERT ON tasks
    WHEN NEW.status NOT IN (
            'pending', 'claimed', 'in_progress', 'retry_wait',
            'complete', 'failed', 'cancelled'
         )
         OR NEW.priority < 0 OR NEW.priority > 3
         OR NEW.attempt_count < 0 OR NEW.max_attempts < 1
         OR NEW.lease_seconds < 1
         OR NOT json_valid(NEW.params)
         OR NOT json_valid(NEW.result_finding_ids)
         OR NOT json_valid(NEW.dependency_ids)
         OR NOT json_valid(NEW.required_capabilities)
         OR NOT json_valid(NEW.resource_request)
    BEGIN
        SELECT RAISE(ABORT, 'invalid Sharur Ops v3 task invariant');
    END
    """,
    """
    CREATE TRIGGER IF NOT EXISTS trg_tasks_v3_update
    BEFORE UPDATE ON tasks
    WHEN NEW.status NOT IN (
            'pending', 'claimed', 'in_progress', 'retry_wait',
            'complete', 'failed', 'cancelled'
         )
         OR NEW.priority < 0 OR NEW.priority > 3
         OR NEW.attempt_count < 0 OR NEW.max_attempts < 1
         OR NEW.lease_seconds < 1
         OR NOT json_valid(NEW.params)
         OR NOT json_valid(NEW.result_finding_ids)
         OR NOT json_valid(NEW.dependency_ids)
         OR NOT json_valid(NEW.required_capabilities)
         OR NOT json_valid(NEW.resource_request)
    BEGIN
        SELECT RAISE(ABORT, 'invalid Sharur Ops v3 task invariant');
    END
    """,
    """
    CREATE TRIGGER IF NOT EXISTS trg_findings_v3_insert
    BEFORE INSERT ON findings
    WHEN NEW.confidence < 0.0 OR NEW.confidence > 1.0
         OR NEW.novelty < 0 OR NEW.novelty > 3
         OR NOT json_valid(NEW.evidence)
    BEGIN
        SELECT RAISE(ABORT, 'invalid Sharur Ops v3 finding invariant');
    END
    """,
    """
    CREATE TRIGGER IF NOT EXISTS trg_hypotheses_v3_insert
    BEFORE INSERT ON hypotheses
    WHEN NEW.status NOT IN (
            'proposed', 'investigating', 'supported', 'refuted', 'inconclusive'
         )
         OR NOT json_valid(NEW.source_finding_ids)
         OR NOT json_valid(NEW.domains_relevant)
         OR NOT json_valid(NEW.evidence_for)
         OR NOT json_valid(NEW.evidence_against)
    BEGIN
        SELECT RAISE(ABORT, 'invalid Sharur Ops v3 hypothesis invariant');
    END
    """,
    """
    CREATE TRIGGER IF NOT EXISTS trg_runs_v3_update
    BEFORE UPDATE ON runs
    WHEN NEW.status NOT IN (
            'pending', 'submitted', 'running', 'complete', 'failed', 'cancelled'
         )
         OR NOT json_valid(NEW.config)
         OR NOT json_valid(NEW.result)
    BEGIN
        SELECT RAISE(ABORT, 'invalid Sharur Ops v3 run invariant');
    END
    """,
)


# Public compatibility alias used by callers and older tests.
OPS_SCHEMA = _TABLES_SQL + _INDEXES_SQL


_MIGRATION_COLUMNS: dict[str, dict[str, str]] = {
    "findings": {
        "campaign_id": "TEXT",
        "task_id": "TEXT",
        "idempotency_key": "TEXT",
        "schema_version": "INTEGER NOT NULL DEFAULT 1",
        "validation_status": "TEXT NOT NULL DEFAULT 'unreviewed'",
        "created_event_id": "INTEGER",
    },
    "hypotheses": {
        "campaign_id": "TEXT",
        "idempotency_key": "TEXT",
        "schema_version": "INTEGER NOT NULL DEFAULT 1",
        "created_event_id": "INTEGER",
    },
    "tasks": {
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
        "lease_token_hash": "TEXT",
        "reserved_for": "TEXT",
        "campaign_id": "TEXT",
        "required_capabilities": "TEXT NOT NULL DEFAULT '[]'",
        "resource_request": "TEXT NOT NULL DEFAULT '{}'",
        "schema_version": "INTEGER NOT NULL DEFAULT 1",
        "created_event_id": "INTEGER",
    },
    "coordinator_log": {
        "actor_agent_id": "TEXT",
        "campaign_id": "TEXT",
        "idempotency_key": "TEXT",
        "created_event_id": "INTEGER",
    },
    "runs": {
        "campaign_id": "TEXT",
        "created_event_id": "INTEGER",
    },
    "run_stages": {
        "owner_agent_id": "TEXT",
    },
}


def _statements(sql: str) -> list[str]:
    return [statement.strip() for statement in sql.split(";") if statement.strip()]


def _table_exists(conn: sqlite3.Connection, table: str) -> bool:
    return (
        conn.execute(
            "SELECT 1 FROM sqlite_master WHERE type IN ('table', 'view') AND name = ?",
            (table,),
        ).fetchone()
        is not None
    )


def _existing_columns(conn: sqlite3.Connection, table: str) -> set[str]:
    if not _table_exists(conn, table):
        return set()
    return {row[1] for row in conn.execute(f"PRAGMA table_info({table})")}


def _schema_version(conn: sqlite3.Connection) -> int:
    if not _table_exists(conn, "ops_schema_meta"):
        return 0
    row = conn.execute("SELECT COALESCE(MAX(version), 0) FROM ops_schema_meta").fetchone()
    return int(row[0]) if row else 0


def _backfill_normalized_links(conn: sqlite3.Connection) -> None:
    """Best-effort normalization of valid legacy references.

    Older coordination records occasionally used external identifiers in the
    JSON reference fields.  Those values remain in the compatibility columns;
    only references resolving to local primary keys enter the FK-backed tables.
    """

    task_ids = {row[0] for row in conn.execute("SELECT id FROM tasks")}
    finding_ids = {row[0] for row in conn.execute("SELECT id FROM findings")}
    hypothesis_ids = {row[0] for row in conn.execute("SELECT id FROM hypotheses")}

    for task_id, raw_dependencies, raw_results in conn.execute(
        "SELECT id, dependency_ids, result_finding_ids FROM tasks"
    ):
        with contextlib.suppress(json.JSONDecodeError, TypeError):
            for dependency_id in json.loads(raw_dependencies or "[]"):
                if dependency_id in task_ids and dependency_id != task_id:
                    conn.execute(
                        "INSERT OR IGNORE INTO task_dependencies(task_id, depends_on_task_id) "
                        "VALUES (?, ?)",
                        (task_id, dependency_id),
                    )
        with contextlib.suppress(json.JSONDecodeError, TypeError):
            for finding_id in json.loads(raw_results or "[]"):
                if finding_id in finding_ids:
                    conn.execute(
                        "INSERT OR IGNORE INTO task_result_findings(task_id, finding_id) "
                        "VALUES (?, ?)",
                        (task_id, finding_id),
                    )

    for hypothesis_id, sources, supporting, opposing in conn.execute(
        "SELECT id, source_finding_ids, evidence_for, evidence_against FROM hypotheses"
    ):
        for relation, raw_ids in (
            ("source", sources),
            ("for", supporting),
            ("against", opposing),
        ):
            with contextlib.suppress(json.JSONDecodeError, TypeError):
                for finding_id in json.loads(raw_ids or "[]"):
                    if finding_id in finding_ids:
                        conn.execute(
                            "INSERT OR IGNORE INTO hypothesis_findings"
                            "(hypothesis_id, finding_id, relation) VALUES (?, ?, ?)",
                            (hypothesis_id, finding_id, relation),
                        )

    for log_id, raw_findings, raw_hypotheses in conn.execute(
        "SELECT id, referenced_finding_ids, referenced_hypothesis_ids "
        "FROM coordinator_log"
    ):
        with contextlib.suppress(json.JSONDecodeError, TypeError):
            for finding_id in json.loads(raw_findings or "[]"):
                if finding_id in finding_ids:
                    conn.execute(
                        "INSERT OR IGNORE INTO coordinator_log_findings(log_id, finding_id) "
                        "VALUES (?, ?)",
                        (log_id, finding_id),
                    )
        with contextlib.suppress(json.JSONDecodeError, TypeError):
            for hypothesis_id in json.loads(raw_hypotheses or "[]"):
                if hypothesis_id in hypothesis_ids:
                    conn.execute(
                        "INSERT OR IGNORE INTO coordinator_log_hypotheses"
                        "(log_id, hypothesis_id) VALUES (?, ?)",
                        (log_id, hypothesis_id),
                    )


def _ensure_fts(conn: sqlite3.Connection) -> None:
    """Create and backfill the optional FTS5 index when SQLite provides it."""

    try:
        conn.execute(
            """
            CREATE VIRTUAL TABLE IF NOT EXISTS findings_fts USING fts5(
                finding_id UNINDEXED,
                summary,
                reasoning,
                evidence,
                tokenize='unicode61'
            )
            """
        )
        conn.execute(
            """
            INSERT INTO findings_fts(finding_id, summary, reasoning, evidence)
            SELECT f.id, f.summary, f.reasoning, f.evidence
            FROM findings AS f
            WHERE NOT EXISTS (
                SELECT 1 FROM findings_fts AS x WHERE x.finding_id = f.id
            )
            """
        )
    except Exception as exc:
        # FTS5 is optional in SQLite builds. Only suppress the feature-missing
        # error; surface malformed schemas and disk failures.
        if "fts5" not in str(exc).lower() and "no such module" not in str(exc).lower():
            raise


def ensure_ops_schema(conn: sqlite3.Connection) -> bool:
    """Create or migrate the ops schema.

    Returns ``True`` when a migration was applied and ``False`` when the
    connection already targeted the current schema.  A current-schema open
    performs PRAGMA configuration and reads only; request-scoped connections
    therefore avoid the former schema-write penalty.
    """

    conn.execute("PRAGMA foreign_keys=ON")
    current = _schema_version(conn)
    if current > OPS_SCHEMA_VERSION:
        raise RuntimeError(
            f"Sharur Ops schema {current} is newer than supported version "
            f"{OPS_SCHEMA_VERSION}"
        )
    if current == OPS_SCHEMA_VERSION:
        return False

    conn.execute("BEGIN IMMEDIATE")
    try:
        # Tables first, then additive columns needed by current indexes.
        for statement in _statements(_TABLES_SQL):
            conn.execute(statement)
        for table, columns in _MIGRATION_COLUMNS.items():
            existing = _existing_columns(conn, table)
            for name, definition in columns.items():
                if name not in existing:
                    conn.execute(f"ALTER TABLE {table} ADD COLUMN {name} {definition}")

        # Legacy active tasks gain finite leases. They can be recovered safely;
        # the first subsequent claim mints an attempt-specific fencing token.
        conn.execute(
            """
            UPDATE tasks
            SET attempt_count = CASE WHEN attempt_count < 1 THEN 1 ELSE attempt_count END,
                heartbeat_ts = COALESCE(heartbeat_ts, claimed_ts, ts),
                reserved_for = COALESCE(reserved_for, assigned_to),
                lease_expires_ts = COALESCE(
                    lease_expires_ts,
                    COALESCE(claimed_ts, ts) + lease_seconds
                )
            WHERE status IN ('claimed', 'in_progress')
            """
        )

        for statement in _statements(_INDEXES_SQL):
            conn.execute(statement)
        for trigger in _VALIDATION_TRIGGERS:
            conn.execute(trigger)
        _backfill_normalized_links(conn)
        _ensure_fts(conn)

        now = time.time()
        conn.execute(
            "INSERT INTO ops_schema_meta(version, applied_ts) VALUES (?, ?)",
            (OPS_SCHEMA_VERSION, now),
        )
        conn.execute(
            """
            INSERT INTO coordination_events(
                ts, event_type, actor_agent_id, entity_type, entity_id, payload
            ) VALUES (?, 'schema_migrated', 'system', 'schema', ?, ?)
            """,
            (
                now,
                str(OPS_SCHEMA_VERSION),
                json.dumps({"from_version": current, "to_version": OPS_SCHEMA_VERSION}),
            ),
        )
        conn.commit()
        return True
    except Exception:
        conn.rollback()
        raise


__all__ = [
    "DEFAULT_LEASE_SECONDS",
    "OPS_SCHEMA",
    "OPS_SCHEMA_VERSION",
    "ensure_ops_schema",
]
