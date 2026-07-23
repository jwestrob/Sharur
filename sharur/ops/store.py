"""Durable SQLite control-plane store for coordinated Sharur agents.

Direct access remains useful for one local process. Distributed agents should
use the HTTP service so exactly one server process owns SQLite, especially
when the database resides on a network filesystem.
"""

from __future__ import annotations

import contextlib
import hashlib
import json
import re
import secrets
import sqlite3
import threading
import time
import uuid
from pathlib import Path
from typing import TYPE_CHECKING, Any

from sharur.ops.db import SQLiteDirectAccessLock, open_ops_connection
from sharur.ops.ledger import RunLedger
from sharur.ops.schema import (
    DEFAULT_LEASE_SECONDS,
    OPS_SCHEMA,
    ensure_ops_schema,
)


if TYPE_CHECKING:
    from collections.abc import Callable, Iterator


_SCHEMA = OPS_SCHEMA
MAX_INLINE_JSON_BYTES = 256 * 1024
MAX_INLINE_TEXT_BYTES = 256 * 1024
DEFAULT_LIST_LIMIT = 50
MAX_LIST_LIMIT = 1_000
AGENT_ROLES = frozenset({"reader", "worker", "coordinator", "operator"})
AGENT_STATUSES = frozenset({"active", "draining", "offline", "revoked"})
RESOURCE_KEYS = frozenset({"cpu_slots", "memory_mb", "accelerator_slots"})

_JSON_FIELDS = {
    "evidence",
    "source_finding_ids",
    "domains_relevant",
    "evidence_for",
    "evidence_against",
    "params",
    "result_finding_ids",
    "dependency_ids",
    "required_capabilities",
    "resource_request",
    "referenced_finding_ids",
    "referenced_hypothesis_ids",
    "decisions_made",
    "capabilities",
    "metadata",
    "payload",
    "config",
    "result",
    "command",
    "inputs",
    "outputs",
    "resource_profile",
}
_SENSITIVE_FIELDS = {"token_hash", "lease_token_hash"}


class LeaseFenceError(ValueError):
    """A worker presented a missing or stale task-attempt fencing token."""


class IdempotencyConflictError(ValueError):
    """An idempotency key was reused with a different immutable payload."""


def _canonical_json(value: Any, *, field: str = "payload") -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)
    if len(encoded.encode("utf-8")) > MAX_INLINE_JSON_BYTES:
        raise ValueError(
            f"{field} exceeds the {MAX_INLINE_JSON_BYTES}-byte inline JSON limit; "
            "register a content-addressed artifact and reference its ID"
        )
    return encoded


def _bounded_text(value: str, *, field: str) -> str:
    if len(value.encode("utf-8")) > MAX_INLINE_TEXT_BYTES:
        raise ValueError(f"{field} exceeds the {MAX_INLINE_TEXT_BYTES}-byte inline limit")
    return value


def _hash_secret(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _validate_idempotency_key(value: str | None) -> None:
    if value is not None and not value.strip():
        raise ValueError("idempotency_key must be non-empty when provided")


def _row_to_dict(row: sqlite3.Row | None) -> dict[str, Any] | None:
    if row is None:
        return None
    result = dict(row)
    for field in _SENSITIVE_FIELDS:
        result.pop(field, None)
    for field in _JSON_FIELDS:
        value = result.get(field)
        if isinstance(value, str):
            with contextlib.suppress(json.JSONDecodeError, TypeError):
                result[field] = json.loads(value)
    return result


def _validate_limit(limit: int) -> int:
    if not 1 <= int(limit) <= MAX_LIST_LIMIT:
        raise ValueError(f"limit must be between 1 and {MAX_LIST_LIMIT}")
    return int(limit)


def _normalize_capabilities(values: list[str] | None) -> list[str]:
    capabilities: list[str] = []
    for value in values or []:
        normalized = str(value).strip()
        if not normalized:
            raise ValueError("required capabilities must be non-empty strings")
        if normalized not in capabilities:
            capabilities.append(normalized)
    return capabilities


def _normalize_resources(value: dict[str, Any] | None) -> dict[str, int]:
    request = dict(value or {})
    unknown = sorted(set(request) - RESOURCE_KEYS)
    if unknown:
        raise ValueError(f"Unsupported resource fields: {', '.join(unknown)}")
    normalized: dict[str, int] = {}
    for key in RESOURCE_KEYS:
        raw = request.get(key, 0)
        if isinstance(raw, bool) or not isinstance(raw, int) or raw < 0:
            raise ValueError(f"{key} must be a non-negative integer")
        if raw:
            normalized[key] = raw
    return normalized


class OpsStore:
    """SQLite operational store with transactional coordination semantics."""

    def __init__(
        self,
        db_path: str | Path = "sharur_ops.db",
        agent_id: str = "unnamed",
        *,
        connection: sqlite3.Connection | None = None,
        initialize: bool = True,
        close_callback: Callable[[sqlite3.Connection], None] | None = None,
        transaction_wait_observer: Callable[[float], None] | None = None,
    ):
        self.db_path = Path(db_path).expanduser().resolve()
        self.agent_id = agent_id
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._lock = threading.RLock()
        self._lease_tokens: dict[str, tuple[str, int]] = {}
        self._close_callback = close_callback
        self._transaction_wait_observer = transaction_wait_observer
        self._owns_connection = connection is None
        self._closed = False
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
        self.ledger = RunLedger(
            self.db_path,
            connection=self._conn,
            lock=self._lock,
            initialize=False,
            agent_id=self.agent_id,
            transaction_wait_observer=transaction_wait_observer,
        )

    @contextlib.contextmanager
    def _transaction(self, *, immediate: bool = True) -> Iterator[None]:
        started = time.perf_counter()
        self._conn.execute("BEGIN IMMEDIATE" if immediate else "BEGIN")
        if immediate and self._transaction_wait_observer is not None:
            self._transaction_wait_observer(time.perf_counter() - started)
        try:
            yield
            self._conn.commit()
        except Exception:
            self._conn.rollback()
            raise

    def _event_locked(
        self,
        event_type: str,
        entity_type: str,
        entity_id: str,
        *,
        campaign_id: str | None = None,
        run_id: str | None = None,
        task_id: str | None = None,
        payload: dict[str, Any] | None = None,
        actor_agent_id: str | None = None,
    ) -> int:
        cursor = self._conn.execute(
            """
            INSERT INTO coordination_events(
                ts, event_type, actor_agent_id, campaign_id, entity_type,
                entity_id, run_id, task_id, payload
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                time.time(),
                event_type,
                actor_agent_id or self.agent_id,
                campaign_id,
                entity_type,
                entity_id,
                run_id,
                task_id,
                _canonical_json(payload or {}, field="event payload"),
            ),
        )
        return int(cursor.lastrowid)

    def _validate_campaign_locked(self, campaign_id: str | None) -> None:
        if campaign_id is None:
            return
        if (
            self._conn.execute(
                "SELECT 1 FROM campaigns WHERE id = ?",
                (campaign_id,),
            ).fetchone()
            is None
        ):
            raise ValueError(f"Campaign {campaign_id} not found")

    def _validate_ids_locked(
        self,
        table: str,
        ids: list[str] | None,
        *,
        label: str,
    ) -> list[str]:
        unique = list(dict.fromkeys(ids or []))
        if not unique:
            return unique
        placeholders = ",".join("?" for _ in unique)
        found = {
            row[0]
            for row in self._conn.execute(
                f"SELECT id FROM {table} WHERE id IN ({placeholders})",
                unique,
            )
        }
        missing = [value for value in unique if value not in found]
        if missing:
            raise ValueError(f"Unknown {label}: {', '.join(missing)}")
        return unique

    # ------------------------------------------------------------------
    # Campaigns
    # ------------------------------------------------------------------

    def create_campaign(
        self,
        name: str,
        *,
        description: str = "",
        dataset_path: str | Path | None = None,
        metadata: dict[str, Any] | None = None,
        idempotency_key: str | None = None,
    ) -> str:
        name = name.strip()
        if not name:
            raise ValueError("campaign name must be non-empty")
        if idempotency_key is not None and not idempotency_key.strip():
            raise ValueError("idempotency_key must be non-empty when provided")
        normalized_path = (
            str(Path(dataset_path).expanduser().resolve())
            if dataset_path is not None
            else None
        )
        metadata_json = _canonical_json(metadata or {}, field="campaign metadata")
        with self._lock, self._transaction():
            if idempotency_key:
                existing = self._conn.execute(
                    "SELECT * FROM campaigns WHERE idempotency_key = ?",
                    (idempotency_key,),
                ).fetchone()
                if existing is not None:
                    immutable = (name, description, normalized_path, metadata_json)
                    current = (
                        existing["name"],
                        existing["description"],
                        existing["dataset_path"],
                        existing["metadata"],
                    )
                    if immutable != current:
                        raise IdempotencyConflictError(
                            "Campaign idempotency key already exists with a different payload"
                        )
                    return str(existing["id"])
            campaign_id = str(uuid.uuid4())
            now = time.time()
            self._conn.execute(
                """
                INSERT INTO campaigns(
                    id, name, description, dataset_path, created_by, created_ts,
                    metadata, idempotency_key
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    campaign_id,
                    name,
                    _bounded_text(description, field="campaign description"),
                    normalized_path,
                    self.agent_id,
                    now,
                    metadata_json,
                    idempotency_key,
                ),
            )
            self._event_locked(
                "campaign_created",
                "campaign",
                campaign_id,
                campaign_id=campaign_id,
                payload={"name": name, "dataset_path": normalized_path},
            )
            return campaign_id

    def get_campaign(self, campaign_id: str) -> dict[str, Any]:
        with self._lock:
            row = self._conn.execute(
                "SELECT * FROM campaigns WHERE id = ?",
                (campaign_id,),
            ).fetchone()
        result = _row_to_dict(row)
        if result is None:
            raise KeyError(f"Campaign {campaign_id} not found")
        return result

    def list_campaigns(
        self,
        *,
        status: str | None = None,
        before_ts: float | None = None,
        limit: int = DEFAULT_LIST_LIMIT,
    ) -> list[dict[str, Any]]:
        clauses: list[str] = []
        params: list[Any] = []
        if status is not None:
            clauses.append("status = ?")
            params.append(status)
        if before_ts is not None:
            clauses.append("created_ts < ?")
            params.append(before_ts)
        where = f"WHERE {' AND '.join(clauses)}" if clauses else ""
        params.append(_validate_limit(limit))
        with self._lock:
            rows = self._conn.execute(
                f"SELECT * FROM campaigns {where} "
                "ORDER BY created_ts DESC, id DESC LIMIT ?",
                params,
            ).fetchall()
        return [_row_to_dict(row) for row in rows]  # type: ignore[misc]

    def update_campaign(
        self,
        campaign_id: str,
        *,
        status: str,
    ) -> dict[str, Any]:
        allowed = {"active", "paused", "complete", "failed", "archived"}
        if status not in allowed:
            raise ValueError(f"Unsupported campaign status: {status}")
        now = time.time()
        with self._lock, self._transaction():
            cursor = self._conn.execute(
                """
                UPDATE campaigns
                SET status = ?,
                    completed_ts = CASE
                        WHEN ? IN ('complete', 'failed', 'archived') THEN ?
                        ELSE NULL
                    END
                WHERE id = ?
                """,
                (status, status, now, campaign_id),
            )
            if cursor.rowcount == 0:
                raise KeyError(f"Campaign {campaign_id} not found")
            self._event_locked(
                "campaign_status_changed",
                "campaign",
                campaign_id,
                campaign_id=campaign_id,
                payload={"status": status},
            )
        return self.get_campaign(campaign_id)

    # ------------------------------------------------------------------
    # Agent identity, capabilities, and capacity
    # ------------------------------------------------------------------

    def register_agent(
        self,
        agent_id: str,
        *,
        role: str = "worker",
        capabilities: list[str] | None = None,
        max_concurrent_tasks: int = 1,
        capacity_cpu_slots: int = 1,
        capacity_memory_mb: int = 0,
        capacity_accelerator_slots: int = 0,
        host: str | None = None,
        metadata: dict[str, Any] | None = None,
        rotate_token: bool = False,
    ) -> dict[str, Any]:
        agent_id = agent_id.strip()
        if not agent_id:
            raise ValueError("agent_id must be non-empty")
        if role not in AGENT_ROLES:
            raise ValueError(f"Unsupported agent role: {role}")
        for field, value in (
            ("max_concurrent_tasks", max_concurrent_tasks),
            ("capacity_cpu_slots", capacity_cpu_slots),
            ("capacity_memory_mb", capacity_memory_mb),
            ("capacity_accelerator_slots", capacity_accelerator_slots),
        ):
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ValueError(f"{field} must be a non-negative integer")
        capabilities_json = _canonical_json(
            _normalize_capabilities(capabilities),
            field="agent capabilities",
        )
        metadata_json = _canonical_json(metadata or {}, field="agent metadata")
        now = time.time()
        with self._lock, self._transaction():
            existing = self._conn.execute(
                "SELECT * FROM agents WHERE id = ?",
                (agent_id,),
            ).fetchone()
            issue_token = existing is None or rotate_token or existing["token_hash"] is None
            raw_token = secrets.token_urlsafe(32) if issue_token else None
            token_hash = _hash_secret(raw_token) if raw_token is not None else existing["token_hash"]
            token_hint = raw_token[-6:] if raw_token is not None else existing["token_hint"]
            if existing is None:
                self._conn.execute(
                    """
                    INSERT INTO agents(
                        id, role, token_hash, token_hint, capabilities,
                        max_concurrent_tasks, capacity_cpu_slots,
                        capacity_memory_mb, capacity_accelerator_slots, host,
                        status, registered_ts, heartbeat_ts, metadata
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 'active', ?, ?, ?)
                    """,
                    (
                        agent_id,
                        role,
                        token_hash,
                        token_hint,
                        capabilities_json,
                        max_concurrent_tasks,
                        capacity_cpu_slots,
                        capacity_memory_mb,
                        capacity_accelerator_slots,
                        host,
                        now,
                        now,
                        metadata_json,
                    ),
                )
                event_type = "agent_registered"
            else:
                self._conn.execute(
                    """
                    UPDATE agents
                    SET role = ?, token_hash = ?, token_hint = ?, capabilities = ?,
                        max_concurrent_tasks = ?, capacity_cpu_slots = ?,
                        capacity_memory_mb = ?, capacity_accelerator_slots = ?,
                        host = ?, status = 'active', heartbeat_ts = ?, metadata = ?
                    WHERE id = ?
                    """,
                    (
                        role,
                        token_hash,
                        token_hint,
                        capabilities_json,
                        max_concurrent_tasks,
                        capacity_cpu_slots,
                        capacity_memory_mb,
                        capacity_accelerator_slots,
                        host,
                        now,
                        metadata_json,
                        agent_id,
                    ),
                )
                event_type = "agent_token_rotated" if rotate_token else "agent_updated"
            self._event_locked(
                event_type,
                "agent",
                agent_id,
                payload={
                    "role": role,
                    "capabilities": json.loads(capabilities_json),
                    "max_concurrent_tasks": max_concurrent_tasks,
                },
            )
            row = self._conn.execute(
                "SELECT * FROM agents WHERE id = ?",
                (agent_id,),
            ).fetchone()
        result = _row_to_dict(row)
        assert result is not None
        if raw_token is not None:
            result["token"] = raw_token
        return result

    def authenticate_agent(self, token: str) -> dict[str, Any] | None:
        token_hash = _hash_secret(token)
        with self._lock:
            row = self._conn.execute(
                """
                SELECT * FROM agents
                WHERE token_hash = ? AND status IN ('active', 'draining')
                """,
                (token_hash,),
            ).fetchone()
        return _row_to_dict(row)

    def get_agent(self, agent_id: str) -> dict[str, Any]:
        with self._lock:
            row = self._conn.execute(
                "SELECT * FROM agents WHERE id = ?",
                (agent_id,),
            ).fetchone()
        result = _row_to_dict(row)
        if result is None:
            raise KeyError(f"Agent {agent_id} not found")
        return result

    def list_agents(
        self,
        *,
        status: str | None = None,
        limit: int = DEFAULT_LIST_LIMIT,
    ) -> list[dict[str, Any]]:
        params: list[Any] = []
        where = ""
        if status is not None:
            where = "WHERE status = ?"
            params.append(status)
        params.append(_validate_limit(limit))
        with self._lock:
            rows = self._conn.execute(
                f"SELECT * FROM agents {where} ORDER BY registered_ts, id LIMIT ?",
                params,
            ).fetchall()
        return [_row_to_dict(row) for row in rows]  # type: ignore[misc]

    def heartbeat_agent(
        self,
        *,
        status: str = "active",
        host: str | None = None,
    ) -> dict[str, Any]:
        if status not in AGENT_STATUSES - {"revoked"}:
            raise ValueError(f"Unsupported heartbeat status: {status}")
        with self._lock, self._transaction():
            cursor = self._conn.execute(
                """
                UPDATE agents
                SET heartbeat_ts = ?, status = ?, host = COALESCE(?, host)
                WHERE id = ? AND status <> 'revoked'
                """,
                (time.time(), status, host, self.agent_id),
            )
            if cursor.rowcount == 0:
                raise KeyError(f"Agent {self.agent_id} not found or revoked")
            self._event_locked(
                "agent_heartbeat",
                "agent",
                self.agent_id,
                payload={"status": status, "host": host},
            )
        return self.get_agent(self.agent_id)

    def set_agent_status(self, agent_id: str, status: str) -> dict[str, Any]:
        if status not in AGENT_STATUSES:
            raise ValueError(f"Unsupported agent status: {status}")
        with self._lock, self._transaction():
            cursor = self._conn.execute(
                "UPDATE agents SET status = ? WHERE id = ?",
                (status, agent_id),
            )
            if cursor.rowcount == 0:
                raise KeyError(f"Agent {agent_id} not found")
            self._event_locked(
                "agent_status_changed",
                "agent",
                agent_id,
                payload={"status": status},
            )
        return self.get_agent(agent_id)

    # ------------------------------------------------------------------
    # Findings and artifacts
    # ------------------------------------------------------------------

    def finding(
        self,
        finding_type: str,
        domain: str,
        summary: str,
        evidence: dict | None = None,
        confidence: float = 0.5,
        novelty: int = 0,
        parent_finding_id: str | None = None,
        reasoning: str = "",
        *,
        campaign_id: str | None = None,
        task_id: str | None = None,
        idempotency_key: str | None = None,
        schema_version: int = 1,
        validation_status: str = "unreviewed",
    ) -> str:
        if not finding_type.strip() or not domain.strip():
            raise ValueError("finding_type and domain must be non-empty")
        _validate_idempotency_key(idempotency_key)
        if not 0.0 <= confidence <= 1.0:
            raise ValueError("confidence must be between 0 and 1")
        if not 0 <= novelty <= 3:
            raise ValueError("novelty must be between 0 and 3")
        if schema_version < 1:
            raise ValueError("schema_version must be positive")
        evidence_json = _canonical_json(evidence or {}, field="finding evidence")
        summary = _bounded_text(summary, field="finding summary")
        reasoning = _bounded_text(reasoning, field="finding reasoning")
        immutable = {
            "finding_type": finding_type,
            "domain": domain,
            "summary": summary,
            "evidence": evidence_json,
            "confidence": confidence,
            "novelty": novelty,
            "parent_finding_id": parent_finding_id,
            "reasoning": reasoning,
            "campaign_id": campaign_id,
            "task_id": task_id,
            "schema_version": schema_version,
            "validation_status": validation_status,
        }
        with self._lock, self._transaction():
            self._validate_campaign_locked(campaign_id)
            if task_id is not None:
                self._validate_ids_locked("tasks", [task_id], label="tasks")
            if parent_finding_id is not None:
                self._validate_ids_locked(
                    "findings",
                    [parent_finding_id],
                    label="parent findings",
                )
            if idempotency_key:
                existing = self._conn.execute(
                    """
                    SELECT * FROM findings
                    WHERE agent_id = ? AND idempotency_key = ?
                    """,
                    (self.agent_id, idempotency_key),
                ).fetchone()
                if existing is not None:
                    conflicts = [
                        key for key, value in immutable.items() if existing[key] != value
                    ]
                    if conflicts:
                        raise IdempotencyConflictError(
                            "Finding idempotency key already exists with a different "
                            f"payload: {', '.join(conflicts)}"
                        )
                    return str(existing["id"])

            finding_id = str(uuid.uuid4())
            now = time.time()
            self._conn.execute(
                """
                INSERT INTO findings(
                    id, agent_id, ts, finding_type, domain, summary, evidence,
                    confidence, novelty, parent_finding_id, reasoning,
                    campaign_id, task_id, idempotency_key, schema_version,
                    validation_status
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    finding_id,
                    self.agent_id,
                    now,
                    finding_type,
                    domain,
                    summary,
                    evidence_json,
                    confidence,
                    novelty,
                    parent_finding_id,
                    reasoning,
                    campaign_id,
                    task_id,
                    idempotency_key,
                    schema_version,
                    validation_status,
                ),
            )
            event_id = self._event_locked(
                "finding_created",
                "finding",
                finding_id,
                campaign_id=campaign_id,
                task_id=task_id,
                payload={
                    "finding_type": finding_type,
                    "domain": domain,
                    "novelty": novelty,
                },
            )
            self._conn.execute(
                "UPDATE findings SET created_event_id = ? WHERE id = ?",
                (event_id, finding_id),
            )
            if self.fts_available:
                self._conn.execute(
                    """
                    INSERT INTO findings_fts(finding_id, summary, reasoning, evidence)
                    VALUES (?, ?, ?, ?)
                    """,
                    (finding_id, summary, reasoning, evidence_json),
                )
            return finding_id

    @property
    def fts_available(self) -> bool:
        return (
            self._conn.execute(
                "SELECT 1 FROM sqlite_master "
                "WHERE type = 'table' AND name = 'findings_fts'"
            ).fetchone()
            is not None
        )

    def recent_findings(
        self,
        since: float = 0,
        min_novelty: int = 0,
        finding_type: str | None = None,
        domain: str | None = None,
        agent_id: str | None = None,
        campaign_id: str | None = None,
        before_ts: float | None = None,
        limit: int = DEFAULT_LIST_LIMIT,
    ) -> list[dict[str, Any]]:
        sql = "SELECT * FROM findings WHERE ts > ? AND novelty >= ?"
        params: list[Any] = [since, min_novelty]
        for column, value in (
            ("finding_type", finding_type),
            ("domain", domain),
            ("agent_id", agent_id),
            ("campaign_id", campaign_id),
        ):
            if value is not None:
                sql += f" AND {column} = ?"
                params.append(value)
        if before_ts is not None:
            sql += " AND ts < ?"
            params.append(before_ts)
        sql += " ORDER BY ts DESC, id DESC LIMIT ?"
        params.append(_validate_limit(limit))
        with self._lock:
            rows = self._conn.execute(sql, params).fetchall()
        return [_row_to_dict(row) for row in rows]  # type: ignore[misc]

    def get_finding(self, finding_id: str) -> dict[str, Any]:
        with self._lock:
            row = self._conn.execute(
                "SELECT * FROM findings WHERE id = ?",
                (finding_id,),
            ).fetchone()
        result = _row_to_dict(row)
        if result is None:
            raise KeyError(f"Finding {finding_id} not found")
        return result

    def search_findings(
        self,
        text: str,
        limit: int = 20,
        *,
        campaign_id: str | None = None,
    ) -> list[dict[str, Any]]:
        limit = _validate_limit(limit)
        tokens = re.findall(r"[\w.-]+", text, flags=re.UNICODE)
        with self._lock:
            if self.fts_available and tokens:
                query = " AND ".join(f'"{token.replace(chr(34), chr(34) * 2)}"' for token in tokens)
                clauses = ["findings_fts MATCH ?"]
                params: list[Any] = [query]
                if campaign_id is not None:
                    clauses.append("f.campaign_id = ?")
                    params.append(campaign_id)
                params.append(limit)
                rows = self._conn.execute(
                    f"""
                    SELECT f.*
                    FROM findings_fts
                    JOIN findings AS f ON f.id = findings_fts.finding_id
                    WHERE {' AND '.join(clauses)}
                    ORDER BY bm25(findings_fts), f.novelty DESC, f.ts DESC
                    LIMIT ?
                    """,
                    params,
                ).fetchall()
            else:
                pattern = f"%{text}%"
                campaign_clause = " AND campaign_id = ?" if campaign_id is not None else ""
                params = [pattern, pattern, pattern]
                if campaign_id is not None:
                    params.append(campaign_id)
                params.append(limit)
                rows = self._conn.execute(
                    f"""
                    SELECT * FROM findings
                    WHERE (summary LIKE ? OR reasoning LIKE ? OR evidence LIKE ?)
                    {campaign_clause}
                    ORDER BY novelty DESC, ts DESC LIMIT ?
                    """,
                    params,
                ).fetchall()
        return [_row_to_dict(row) for row in rows]  # type: ignore[misc]

    def link_findings(
        self,
        finding_id: str,
        related_finding_id: str,
        *,
        relation: str,
    ) -> None:
        relation = relation.strip()
        if not relation:
            raise ValueError("relation must be non-empty")
        with self._lock, self._transaction():
            self._validate_ids_locked(
                "findings",
                [finding_id, related_finding_id],
                label="findings",
            )
            cursor = self._conn.execute(
                """
                INSERT OR IGNORE INTO finding_links(
                    finding_id, related_finding_id, relation, created_ts, created_by
                ) VALUES (?, ?, ?, ?, ?)
                """,
                (finding_id, related_finding_id, relation, time.time(), self.agent_id),
            )
            if cursor.rowcount:
                finding = self._conn.execute(
                    "SELECT campaign_id FROM findings WHERE id = ?",
                    (finding_id,),
                ).fetchone()
                self._event_locked(
                    "findings_linked",
                    "finding",
                    finding_id,
                    campaign_id=finding["campaign_id"],
                    payload={
                        "related_finding_id": related_finding_id,
                        "relation": relation,
                    },
                )

    def register_artifact(
        self,
        content_hash: str,
        uri: str,
        size_bytes: int,
        *,
        media_type: str = "application/octet-stream",
        campaign_id: str | None = None,
        task_id: str | None = None,
        run_id: str | None = None,
        metadata: dict[str, Any] | None = None,
    ) -> str:
        content_hash = content_hash.strip().lower()
        if re.fullmatch(r"[0-9a-f]{64}", content_hash):
            content_hash = f"sha256:{content_hash}"
        if re.fullmatch(r"[a-z0-9_-]+:[0-9a-f]{16,}", content_hash) is None:
            raise ValueError("content_hash must be algorithm:hex (or a 64-character SHA-256)")
        if not uri.strip():
            raise ValueError("artifact uri must be non-empty")
        if isinstance(size_bytes, bool) or not isinstance(size_bytes, int) or size_bytes < 0:
            raise ValueError("size_bytes must be a non-negative integer")
        metadata_json = _canonical_json(metadata or {}, field="artifact metadata")
        with self._lock, self._transaction():
            self._validate_campaign_locked(campaign_id)
            if task_id:
                self._validate_ids_locked("tasks", [task_id], label="tasks")
            if run_id:
                self._validate_ids_locked("runs", [run_id], label="runs")
            existing = self._conn.execute(
                "SELECT * FROM artifacts WHERE content_hash = ?",
                (content_hash,),
            ).fetchone()
            if existing is not None:
                if int(existing["size_bytes"]) != size_bytes:
                    raise IdempotencyConflictError(
                        "Artifact content hash already exists with a different size"
                    )
                return str(existing["id"])
            artifact_id = str(uuid.uuid4())
            self._conn.execute(
                """
                INSERT INTO artifacts(
                    id, content_hash, uri, size_bytes, media_type, created_by,
                    created_ts, campaign_id, task_id, run_id, metadata
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    artifact_id,
                    content_hash,
                    uri,
                    size_bytes,
                    media_type,
                    self.agent_id,
                    time.time(),
                    campaign_id,
                    task_id,
                    run_id,
                    metadata_json,
                ),
            )
            self._event_locked(
                "artifact_registered",
                "artifact",
                artifact_id,
                campaign_id=campaign_id,
                task_id=task_id,
                run_id=run_id,
                payload={
                    "content_hash": content_hash,
                    "size_bytes": size_bytes,
                    "media_type": media_type,
                },
            )
            return artifact_id

    def attach_artifact(
        self,
        finding_id: str,
        artifact_id: str,
        *,
        relation: str = "evidence",
    ) -> None:
        if not relation.strip():
            raise ValueError("relation must be non-empty")
        with self._lock, self._transaction():
            self._validate_ids_locked("findings", [finding_id], label="findings")
            self._validate_ids_locked("artifacts", [artifact_id], label="artifacts")
            cursor = self._conn.execute(
                """
                INSERT OR IGNORE INTO finding_artifacts(finding_id, artifact_id, relation)
                VALUES (?, ?, ?)
                """,
                (finding_id, artifact_id, relation),
            )
            if cursor.rowcount:
                finding = self._conn.execute(
                    "SELECT campaign_id, task_id FROM findings WHERE id = ?",
                    (finding_id,),
                ).fetchone()
                self._event_locked(
                    "artifact_attached",
                    "finding",
                    finding_id,
                    campaign_id=finding["campaign_id"],
                    task_id=finding["task_id"],
                    payload={"artifact_id": artifact_id, "relation": relation},
                )

    def get_artifact(self, artifact_id: str) -> dict[str, Any]:
        with self._lock:
            row = self._conn.execute(
                "SELECT * FROM artifacts WHERE id = ?",
                (artifact_id,),
            ).fetchone()
        result = _row_to_dict(row)
        if result is None:
            raise KeyError(f"Artifact {artifact_id} not found")
        return result

    # ------------------------------------------------------------------
    # Hypotheses
    # ------------------------------------------------------------------

    def hypothesis(
        self,
        hypothesis: str,
        source_finding_ids: list[str] | None = None,
        domains_relevant: list[str] | None = None,
        *,
        campaign_id: str | None = None,
        idempotency_key: str | None = None,
        schema_version: int = 1,
    ) -> str:
        hypothesis = _bounded_text(hypothesis, field="hypothesis")
        if not hypothesis.strip():
            raise ValueError("hypothesis must be non-empty")
        _validate_idempotency_key(idempotency_key)
        domains = list(dict.fromkeys(domains_relevant or []))
        sources_json = _canonical_json(
            list(dict.fromkeys(source_finding_ids or [])),
            field="hypothesis source findings",
        )
        domains_json = _canonical_json(domains, field="hypothesis domains")
        with self._lock, self._transaction():
            self._validate_campaign_locked(campaign_id)
            sources = self._validate_ids_locked(
                "findings",
                source_finding_ids,
                label="source findings",
            )
            if idempotency_key:
                existing = self._conn.execute(
                    """
                    SELECT * FROM hypotheses
                    WHERE source_agent_id = ? AND idempotency_key = ?
                    """,
                    (self.agent_id, idempotency_key),
                ).fetchone()
                if existing is not None:
                    immutable = (
                        existing["hypothesis"],
                        existing["source_finding_ids"],
                        existing["domains_relevant"],
                        existing["campaign_id"],
                        existing["schema_version"],
                    )
                    requested = (
                        hypothesis,
                        sources_json,
                        domains_json,
                        campaign_id,
                        schema_version,
                    )
                    if immutable != requested:
                        raise IdempotencyConflictError(
                            "Hypothesis idempotency key already exists with a different payload"
                        )
                    return str(existing["id"])
            hypothesis_id = str(uuid.uuid4())
            self._conn.execute(
                """
                INSERT INTO hypotheses(
                    id, source_agent_id, source_finding_ids, ts, hypothesis,
                    domains_relevant, campaign_id, idempotency_key, schema_version
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    hypothesis_id,
                    self.agent_id,
                    sources_json,
                    time.time(),
                    hypothesis,
                    domains_json,
                    campaign_id,
                    idempotency_key,
                    schema_version,
                ),
            )
            for finding_id in sources:
                self._conn.execute(
                    """
                    INSERT INTO hypothesis_findings(hypothesis_id, finding_id, relation)
                    VALUES (?, ?, 'source')
                    """,
                    (hypothesis_id, finding_id),
                )
            event_id = self._event_locked(
                "hypothesis_created",
                "hypothesis",
                hypothesis_id,
                campaign_id=campaign_id,
                payload={"source_finding_count": len(sources)},
            )
            self._conn.execute(
                "UPDATE hypotheses SET created_event_id = ? WHERE id = ?",
                (event_id, hypothesis_id),
            )
            return hypothesis_id

    def list_hypotheses(
        self,
        *,
        status: str | None = None,
        unassigned: bool = False,
        campaign_id: str | None = None,
        before_ts: float | None = None,
        limit: int = DEFAULT_LIST_LIMIT,
    ) -> list[dict[str, Any]]:
        clauses: list[str] = []
        params: list[Any] = []
        if status is not None:
            clauses.append("status = ?")
            params.append(status)
        if unassigned:
            clauses.append("assigned_agent_id IS NULL")
        if campaign_id is not None:
            clauses.append("campaign_id = ?")
            params.append(campaign_id)
        if before_ts is not None:
            clauses.append("ts < ?")
            params.append(before_ts)
        where = f"WHERE {' AND '.join(clauses)}" if clauses else ""
        params.append(_validate_limit(limit))
        with self._lock:
            rows = self._conn.execute(
                f"SELECT * FROM hypotheses {where} "
                "ORDER BY ts DESC, id DESC LIMIT ?",
                params,
            ).fetchall()
        return [_row_to_dict(row) for row in rows]  # type: ignore[misc]

    def open_hypotheses(self, unassigned: bool = True) -> list[dict[str, Any]]:
        return self.list_hypotheses(status="proposed", unassigned=unassigned)

    def get_hypothesis(self, hypothesis_id: str) -> dict[str, Any]:
        with self._lock:
            row = self._conn.execute(
                "SELECT * FROM hypotheses WHERE id = ?",
                (hypothesis_id,),
            ).fetchone()
        result = _row_to_dict(row)
        if result is None:
            raise KeyError(f"Hypothesis {hypothesis_id} not found")
        return result

    def update_hypothesis(self, hyp_id: str, **kwargs: Any) -> dict[str, Any]:
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
        with self._lock, self._transaction():
            row = self._conn.execute(
                "SELECT * FROM hypotheses WHERE id = ?",
                (hyp_id,),
            ).fetchone()
            if row is None:
                raise KeyError(f"Hypothesis {hyp_id} not found")
            existing = _row_to_dict(row)
            assert existing is not None
            sets: list[str] = []
            params: list[Any] = []
            for key, value in kwargs.items():
                updated = value
                if key in {"evidence_for", "evidence_against"}:
                    relation = "for" if key == "evidence_for" else "against"
                    finding_ids = self._validate_ids_locked(
                        "findings",
                        value or [],
                        label=f"{relation} evidence findings",
                    )
                    combined = list(existing.get(key) or [])
                    combined.extend(finding_ids)
                    updated = list(dict.fromkeys(combined))
                    for finding_id in finding_ids:
                        self._conn.execute(
                            """
                            INSERT OR IGNORE INTO hypothesis_findings(
                                hypothesis_id, finding_id, relation
                            ) VALUES (?, ?, ?)
                            """,
                            (hyp_id, finding_id, relation),
                        )
                if isinstance(updated, (list, dict)):
                    updated = _canonical_json(updated, field=f"hypothesis {key}")
                if isinstance(updated, str):
                    updated = _bounded_text(updated, field=f"hypothesis {key}")
                sets.append(f"{key} = ?")
                params.append(updated)
            params.append(hyp_id)
            self._conn.execute(
                f"UPDATE hypotheses SET {', '.join(sets)} WHERE id = ?",
                params,
            )
            self._event_locked(
                "hypothesis_updated",
                "hypothesis",
                hyp_id,
                campaign_id=row["campaign_id"],
                payload={"fields": sorted(kwargs)},
            )
            updated_row = self._conn.execute(
                "SELECT * FROM hypotheses WHERE id = ?",
                (hyp_id,),
            ).fetchone()
        result = _row_to_dict(updated_row)
        assert result is not None
        return result

    # ------------------------------------------------------------------
    # Tasks
    # ------------------------------------------------------------------

    def create_task(
        self,
        task_type: str,
        description: str,
        params: dict | None = None,
        priority: int = 1,
        domain_hint: str | None = None,
        assigned_to: str | None = None,
        *,
        run_id: str | None = None,
        campaign_id: str | None = None,
        idempotency_key: str | None = None,
        depends_on: list[str] | None = None,
        required_capabilities: list[str] | None = None,
        resource_request: dict[str, Any] | None = None,
        max_attempts: int = 3,
        lease_seconds: int = DEFAULT_LEASE_SECONDS,
    ) -> str:
        if not task_type.strip() or not description.strip():
            raise ValueError("task_type and description must be non-empty")
        if not 0 <= priority <= 3:
            raise ValueError("priority must be between 0 and 3")
        if max_attempts < 1:
            raise ValueError("max_attempts must be positive")
        if lease_seconds < 1:
            raise ValueError("lease_seconds must be positive")
        if idempotency_key is not None and not idempotency_key.strip():
            raise ValueError("idempotency_key must be non-empty when provided")
        dependencies = list(dict.fromkeys(depends_on or []))
        capabilities = _normalize_capabilities(required_capabilities)
        resources = _normalize_resources(resource_request)
        params_json = _canonical_json(params or {}, field="task params")
        dependencies_json = _canonical_json(dependencies, field="task dependencies")
        capabilities_json = _canonical_json(capabilities, field="task capabilities")
        resources_json = _canonical_json(resources, field="task resources")
        description = _bounded_text(description, field="task description")
        immutable = {
            "task_type": task_type,
            "description": description,
            "params": params_json,
            "priority": priority,
            "domain_hint": domain_hint,
            "run_id": run_id,
            "campaign_id": campaign_id,
            "dependency_ids": dependencies_json,
            "required_capabilities": capabilities_json,
            "resource_request": resources_json,
            "reserved_for": assigned_to,
            "max_attempts": max_attempts,
            "lease_seconds": lease_seconds,
        }
        with self._lock, self._transaction():
            self._validate_campaign_locked(campaign_id)
            if run_id is not None:
                self._validate_ids_locked("runs", [run_id], label="runs")
            self._validate_dependencies_locked(dependencies)
            if idempotency_key:
                existing = self._conn.execute(
                    "SELECT * FROM tasks WHERE idempotency_key = ?",
                    (idempotency_key,),
                ).fetchone()
                if existing is not None:
                    conflicts = [
                        key for key, value in immutable.items() if existing[key] != value
                    ]
                    if conflicts:
                        raise IdempotencyConflictError(
                            "Task idempotency key already exists with a different "
                            f"payload: {', '.join(conflicts)}"
                        )
                    return str(existing["id"])

            task_id = str(uuid.uuid4())
            now = time.time()
            # Preassignment is a reservation. The worker explicitly claims it
            # and receives the attempt-specific lease token.
            self._conn.execute(
                """
                INSERT INTO tasks(
                    id, created_by, assigned_to, reserved_for, ts, status,
                    priority, task_type, description, params, domain_hint,
                    run_id, campaign_id, idempotency_key, dependency_ids,
                    required_capabilities, resource_request, attempt_count,
                    max_attempts, lease_seconds
                ) VALUES (?, ?, ?, ?, ?, 'pending', ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 0, ?, ?)
                """,
                (
                    task_id,
                    self.agent_id,
                    assigned_to,
                    assigned_to,
                    now,
                    priority,
                    task_type,
                    description,
                    params_json,
                    domain_hint,
                    run_id,
                    campaign_id,
                    idempotency_key,
                    dependencies_json,
                    capabilities_json,
                    resources_json,
                    max_attempts,
                    lease_seconds,
                ),
            )
            for dependency_id in dependencies:
                self._conn.execute(
                    """
                    INSERT INTO task_dependencies(task_id, depends_on_task_id)
                    VALUES (?, ?)
                    """,
                    (task_id, dependency_id),
                )
            event_id = self._event_locked(
                "task_created",
                "task",
                task_id,
                campaign_id=campaign_id,
                run_id=run_id,
                task_id=task_id,
                payload={
                    "task_type": task_type,
                    "priority": priority,
                    "reserved_for": assigned_to,
                    "required_capabilities": capabilities,
                    "resource_request": resources,
                },
            )
            self._conn.execute(
                "UPDATE tasks SET created_event_id = ? WHERE id = ?",
                (event_id, task_id),
            )
            return task_id

    def my_tasks(
        self,
        status: str | None = None,
        *,
        limit: int = DEFAULT_LIST_LIMIT,
    ) -> list[dict[str, Any]]:
        return self.list_tasks(status=status, assigned_to=self.agent_id, limit=limit)

    def list_tasks(
        self,
        *,
        status: str | None = None,
        assigned_to: str | None = None,
        unassigned: bool = False,
        campaign_id: str | None = None,
        before_ts: float | None = None,
        limit: int = DEFAULT_LIST_LIMIT,
    ) -> list[dict[str, Any]]:
        if unassigned:
            return self.available_tasks(
                campaign_id=campaign_id,
                before_ts=before_ts,
                limit=limit,
            )
        clauses: list[str] = []
        params: list[Any] = []
        if status is not None:
            clauses.append("status = ?")
            params.append(status)
        if assigned_to is not None:
            clauses.append("assigned_to = ?")
            params.append(assigned_to)
        if campaign_id is not None:
            clauses.append("campaign_id = ?")
            params.append(campaign_id)
        if before_ts is not None:
            clauses.append("ts < ?")
            params.append(before_ts)
        where = f"WHERE {' AND '.join(clauses)}" if clauses else ""
        params.append(_validate_limit(limit))
        with self._lock:
            rows = self._conn.execute(
                f"SELECT * FROM tasks {where} "
                "ORDER BY priority DESC, ts ASC, id ASC LIMIT ?",
                params,
            ).fetchall()
        return [_row_to_dict(row) for row in rows]  # type: ignore[misc]

    def available_tasks(
        self,
        *,
        campaign_id: str | None = None,
        before_ts: float | None = None,
        limit: int = DEFAULT_LIST_LIMIT,
    ) -> list[dict[str, Any]]:
        """Read a bounded queue view; recovery and claim occur separately."""

        clauses = [
            "t.status = 'pending'",
            "t.assigned_to IS NULL",
            "t.attempt_count < t.max_attempts",
            """
            NOT EXISTS (
                SELECT 1
                FROM task_dependencies AS d
                JOIN tasks AS parent ON parent.id = d.depends_on_task_id
                WHERE d.task_id = t.id AND parent.status <> 'complete'
            )
            """,
        ]
        params: list[Any] = []
        if campaign_id is not None:
            clauses.append("t.campaign_id = ?")
            params.append(campaign_id)
        if before_ts is not None:
            clauses.append("t.ts < ?")
            params.append(before_ts)
        params.append(_validate_limit(limit))
        with self._lock:
            rows = self._conn.execute(
                f"""
                SELECT t.* FROM tasks AS t
                WHERE {' AND '.join(clauses)}
                ORDER BY t.priority DESC, t.ts ASC, t.id ASC
                LIMIT ?
                """,
                params,
            ).fetchall()
        return [_row_to_dict(row) for row in rows]  # type: ignore[misc]

    def _agent_dispatch_profile_locked(
        self,
        capabilities: list[str] | None,
        capacity: dict[str, int] | None,
    ) -> tuple[set[str], int | None, dict[str, int]]:
        row = self._conn.execute(
            "SELECT * FROM agents WHERE id = ?",
            (self.agent_id,),
        ).fetchone()
        if row is None:
            return (
                set(_normalize_capabilities(capabilities)),
                None,
                _normalize_resources(capacity),
            )
        if row["status"] != "active":
            raise ValueError(f"Agent {self.agent_id} is {row['status']} and cannot claim work")
        registered = set(json.loads(row["capabilities"] or "[]"))
        if capabilities is not None and not set(capabilities) <= registered:
            raise ValueError("Claim capabilities exceed the registered agent capabilities")
        limits = {
            "cpu_slots": int(row["capacity_cpu_slots"]),
            "memory_mb": int(row["capacity_memory_mb"]),
            "accelerator_slots": int(row["capacity_accelerator_slots"]),
        }
        return registered, int(row["max_concurrent_tasks"]), limits

    def _active_resource_usage_locked(self) -> tuple[int, dict[str, int]]:
        rows = self._conn.execute(
            """
            SELECT resource_request
            FROM tasks
            WHERE assigned_to = ? AND status IN ('claimed', 'in_progress')
            """,
            (self.agent_id,),
        ).fetchall()
        totals = dict.fromkeys(RESOURCE_KEYS, 0)
        for row in rows:
            request = _normalize_resources(json.loads(row["resource_request"] or "{}"))
            for key, amount in request.items():
                totals[key] += amount
        return len(rows), totals

    @staticmethod
    def _fits_dispatch_profile(
        row: sqlite3.Row,
        *,
        capabilities: set[str],
        max_concurrent: int | None,
        limits: dict[str, int],
        active_count: int,
        usage: dict[str, int],
    ) -> bool:
        required = set(json.loads(row["required_capabilities"] or "[]"))
        if not required <= capabilities:
            return False
        if max_concurrent is not None and active_count >= max_concurrent:
            return False
        request = _normalize_resources(json.loads(row["resource_request"] or "{}"))
        for key, amount in request.items():
            if key in limits and usage.get(key, 0) + amount > limits[key]:
                return False
        return True

    def _claim_row_locked(
        self,
        row: sqlite3.Row,
        *,
        lease_seconds: int | None,
    ) -> dict[str, Any]:
        now = time.time()
        duration = int(lease_seconds or row["lease_seconds"])
        if duration < 1:
            raise ValueError("lease_seconds must be positive")
        raw_token = secrets.token_urlsafe(32)
        token_hash = _hash_secret(raw_token)
        cursor = self._conn.execute(
            """
            UPDATE tasks
            SET status = 'claimed', assigned_to = ?, claimed_ts = ?,
                heartbeat_ts = ?, lease_expires_ts = ?,
                lease_token_hash = ?, attempt_count = attempt_count + 1,
                retry_after_ts = NULL, completed_ts = NULL
            WHERE id = ? AND status = 'pending'
              AND (reserved_for IS NULL OR reserved_for = ?)
              AND (assigned_to IS NULL OR assigned_to = ?)
              AND attempt_count < max_attempts
            """,
            (
                self.agent_id,
                now,
                now,
                now + duration,
                token_hash,
                row["id"],
                self.agent_id,
                self.agent_id,
            ),
        )
        if cursor.rowcount == 0:
            raise ValueError(f"Task {row['id']} already claimed, exhausted, or reserved")
        claimed = self._conn.execute(
            "SELECT * FROM tasks WHERE id = ?",
            (row["id"],),
        ).fetchone()
        assert claimed is not None
        attempt = int(claimed["attempt_count"])
        self._event_locked(
            "task_claimed",
            "task",
            row["id"],
            campaign_id=claimed["campaign_id"],
            run_id=claimed["run_id"],
            task_id=row["id"],
            payload={"attempt": attempt, "lease_expires_ts": now + duration},
        )
        self._lease_tokens[str(row["id"])] = (raw_token, attempt)
        result = _row_to_dict(claimed)
        assert result is not None
        result["lease_token"] = raw_token
        result["lease_attempt"] = attempt
        return result

    def claim_task(
        self,
        task_id: str,
        *,
        lease_seconds: int | None = None,
        capabilities: list[str] | None = None,
        capacity: dict[str, int] | None = None,
    ) -> dict[str, Any]:
        with self._lock, self._transaction():
            now = time.time()
            self._recover_expired_locked(now)
            self._promote_retries_locked(now)
            row = self._conn.execute(
                "SELECT * FROM tasks WHERE id = ?",
                (task_id,),
            ).fetchone()
            if row is None:
                raise ValueError(f"Task {task_id} not found")
            if not self._dependencies_complete_locked(task_id):
                raise ValueError(f"Task {task_id} has incomplete dependencies")
            profile = self._agent_dispatch_profile_locked(capabilities, capacity)
            active_count, usage = self._active_resource_usage_locked()
            if not self._fits_dispatch_profile(
                row,
                capabilities=profile[0],
                max_concurrent=profile[1],
                limits=profile[2],
                active_count=active_count,
                usage=usage,
            ):
                raise ValueError(
                    f"Task {task_id} exceeds agent capabilities or available capacity"
                )
            return self._claim_row_locked(row, lease_seconds=lease_seconds)

    def claim_next_task(
        self,
        *,
        campaign_id: str | None = None,
        task_types: list[str] | None = None,
        lease_seconds: int | None = None,
        capabilities: list[str] | None = None,
        capacity: dict[str, int] | None = None,
        candidate_limit: int = 100,
    ) -> dict[str, Any] | None:
        candidate_limit = _validate_limit(candidate_limit)
        with self._lock, self._transaction():
            now = time.time()
            self._recover_expired_locked(now)
            self._promote_retries_locked(now)
            profile = self._agent_dispatch_profile_locked(capabilities, capacity)
            active_count, usage = self._active_resource_usage_locked()
            if profile[1] is not None and active_count >= profile[1]:
                return None
            clauses = [
                "t.status = 'pending'",
                "t.attempt_count < t.max_attempts",
                "(t.reserved_for IS NULL OR t.reserved_for = ?)",
                "(t.assigned_to IS NULL OR t.assigned_to = ?)",
                """
                NOT EXISTS (
                    SELECT 1
                    FROM task_dependencies AS d
                    JOIN tasks AS parent ON parent.id = d.depends_on_task_id
                    WHERE d.task_id = t.id AND parent.status <> 'complete'
                )
                """,
            ]
            params: list[Any] = [self.agent_id, self.agent_id]
            available_capabilities = sorted(profile[0])
            if available_capabilities:
                placeholders = ",".join("?" for _ in available_capabilities)
                clauses.append(
                    "NOT EXISTS ("
                    "SELECT 1 FROM json_each(t.required_capabilities) AS required "
                    f"WHERE required.value NOT IN ({placeholders})"
                    ")"
                )
                params.extend(available_capabilities)
            else:
                clauses.append("json_array_length(t.required_capabilities) = 0")
            for resource in sorted(profile[2]):
                available = max(0, profile[2][resource] - usage.get(resource, 0))
                clauses.append(
                    "COALESCE(CAST(json_extract("
                    f"t.resource_request, '$.{resource}') AS INTEGER), 0) <= ?"
                )
                params.append(available)
            if campaign_id is not None:
                clauses.append("t.campaign_id = ?")
                params.append(campaign_id)
            if task_types:
                unique_types = list(dict.fromkeys(task_types))
                clauses.append(
                    f"t.task_type IN ({','.join('?' for _ in unique_types)})"
                )
                params.extend(unique_types)
            params.append(candidate_limit)
            rows = self._conn.execute(
                f"""
                SELECT t.* FROM tasks AS t
                WHERE {' AND '.join(clauses)}
                ORDER BY t.priority DESC, t.ts ASC, t.id ASC
                LIMIT ?
                """,
                params,
            ).fetchall()
            for row in rows:
                if self._fits_dispatch_profile(
                    row,
                    capabilities=profile[0],
                    max_concurrent=profile[1],
                    limits=profile[2],
                    active_count=active_count,
                    usage=usage,
                ):
                    return self._claim_row_locked(row, lease_seconds=lease_seconds)
            return None

    def _resolve_lease(
        self,
        task_id: str,
        lease_token: str | None,
        attempt: int | None,
    ) -> tuple[str, int]:
        cached = self._lease_tokens.get(task_id)
        token = lease_token or (cached[0] if cached else None)
        resolved_attempt = attempt if attempt is not None else (cached[1] if cached else None)
        if token is None or resolved_attempt is None:
            raise LeaseFenceError(
                f"Task {task_id} requires the lease_token and lease_attempt "
                "returned by its claim"
            )
        return token, int(resolved_attempt)

    def heartbeat_task(
        self,
        task_id: str,
        *,
        lease_seconds: int | None = None,
        in_progress: bool = True,
        lease_token: str | None = None,
        attempt: int | None = None,
    ) -> dict[str, Any]:
        token, attempt_number = self._resolve_lease(task_id, lease_token, attempt)
        with self._lock, self._transaction():
            now = time.time()
            row = self._conn.execute(
                "SELECT * FROM tasks WHERE id = ?",
                (task_id,),
            ).fetchone()
            if row is None:
                raise ValueError(f"Task {task_id} not found")
            duration = int(lease_seconds or row["lease_seconds"])
            if duration < 1:
                raise ValueError("lease_seconds must be positive")
            status = "in_progress" if in_progress else "claimed"
            cursor = self._conn.execute(
                """
                UPDATE tasks
                SET status = ?, heartbeat_ts = ?, lease_expires_ts = ?
                WHERE id = ? AND assigned_to = ?
                  AND attempt_count = ? AND lease_token_hash = ?
                  AND status IN ('claimed', 'in_progress')
                  AND lease_expires_ts > ?
                """,
                (
                    status,
                    now,
                    now + duration,
                    task_id,
                    self.agent_id,
                    attempt_number,
                    _hash_secret(token),
                    now,
                ),
            )
            if cursor.rowcount == 0:
                raise LeaseFenceError(
                    f"Task {task_id} has no live attempt {attempt_number} lease "
                    f"owned by {self.agent_id}"
                )
            updated = self._conn.execute(
                "SELECT * FROM tasks WHERE id = ?",
                (task_id,),
            ).fetchone()
            self._event_locked(
                "task_heartbeat",
                "task",
                task_id,
                campaign_id=updated["campaign_id"],
                run_id=updated["run_id"],
                task_id=task_id,
                payload={
                    "attempt": attempt_number,
                    "status": status,
                    "lease_expires_ts": now + duration,
                },
            )
        result = _row_to_dict(updated)
        assert result is not None
        return result

    def put_task_checkpoint(
        self,
        task_id: str,
        checkpoint_key: str,
        *,
        cursor: str | None = None,
        payload: dict[str, Any] | None = None,
        lease_token: str | None = None,
        attempt: int | None = None,
    ) -> dict[str, Any]:
        """Persist retry-visible progress under an active attempt fence."""
        normalized_key = checkpoint_key.strip()
        if not normalized_key or len(normalized_key) > 256:
            raise ValueError("checkpoint_key must contain 1-256 characters")
        if cursor is not None:
            if len(cursor) > 4_096:
                raise ValueError("checkpoint cursor exceeds 4,096 characters")
            cursor = _bounded_text(cursor, field="task checkpoint cursor")
        payload_json = _canonical_json(
            payload or {},
            field="task checkpoint payload",
        )
        token, attempt_number = self._resolve_lease(
            task_id,
            lease_token,
            attempt,
        )
        with self._lock, self._transaction():
            now = time.time()
            row = self._conn.execute(
                """
                SELECT campaign_id, run_id
                FROM tasks
                WHERE id = ? AND assigned_to = ?
                  AND attempt_count = ? AND lease_token_hash = ?
                  AND status IN ('claimed', 'in_progress')
                  AND lease_expires_ts > ?
                """,
                (
                    task_id,
                    self.agent_id,
                    attempt_number,
                    _hash_secret(token),
                    now,
                ),
            ).fetchone()
            if row is None:
                raise LeaseFenceError(
                    f"Task {task_id} has no live attempt {attempt_number} lease "
                    f"owned by {self.agent_id}"
                )
            self._conn.execute(
                """
                INSERT INTO task_checkpoints(
                    task_id, checkpoint_key, attempt, agent_id,
                    cursor, payload, updated_ts
                ) VALUES (?, ?, ?, ?, ?, ?, ?)
                ON CONFLICT(task_id, checkpoint_key) DO UPDATE SET
                    attempt = excluded.attempt,
                    agent_id = excluded.agent_id,
                    cursor = excluded.cursor,
                    payload = excluded.payload,
                    updated_ts = excluded.updated_ts
                """,
                (
                    task_id,
                    normalized_key,
                    attempt_number,
                    self.agent_id,
                    cursor,
                    payload_json,
                    now,
                ),
            )
            checkpoint = self._conn.execute(
                """
                SELECT *
                FROM task_checkpoints
                WHERE task_id = ? AND checkpoint_key = ?
                """,
                (task_id, normalized_key),
            ).fetchone()
        result = _row_to_dict(checkpoint)
        assert result is not None
        return result

    def get_task_checkpoint(
        self,
        task_id: str,
        checkpoint_key: str,
    ) -> dict[str, Any] | None:
        """Read the latest checkpoint, including one written by a prior attempt."""
        normalized_key = checkpoint_key.strip()
        if not normalized_key or len(normalized_key) > 256:
            raise ValueError("checkpoint_key must contain 1-256 characters")
        with self._lock:
            checkpoint = self._conn.execute(
                """
                SELECT *
                FROM task_checkpoints
                WHERE task_id = ? AND checkpoint_key = ?
                """,
                (task_id, normalized_key),
            ).fetchone()
        return _row_to_dict(checkpoint)

    def list_task_checkpoints(
        self,
        task_id: str,
        *,
        limit: int = DEFAULT_LIST_LIMIT,
    ) -> list[dict[str, Any]]:
        """List bounded retry-visible checkpoints for one task."""
        with self._lock:
            if (
                self._conn.execute(
                    "SELECT 1 FROM tasks WHERE id = ?",
                    (task_id,),
                ).fetchone()
                is None
            ):
                raise KeyError(f"Task {task_id} not found")
            rows = self._conn.execute(
                """
                SELECT *
                FROM task_checkpoints
                WHERE task_id = ?
                ORDER BY updated_ts DESC, checkpoint_key
                LIMIT ?
                """,
                (task_id, _validate_limit(limit)),
            ).fetchall()
        return [_row_to_dict(row) for row in rows]  # type: ignore[misc]

    def complete_task(
        self,
        task_id: str,
        result_finding_ids: list[str] | None = None,
        *,
        lease_token: str | None = None,
        attempt: int | None = None,
    ) -> dict[str, Any]:
        token, attempt_number = self._resolve_lease(task_id, lease_token, attempt)
        with self._lock, self._transaction():
            finding_ids = self._validate_ids_locked(
                "findings",
                result_finding_ids,
                label="result findings",
            )
            now = time.time()
            cursor = self._conn.execute(
                """
                UPDATE tasks
                SET status = 'complete', completed_ts = ?,
                    result_finding_ids = ?, heartbeat_ts = ?,
                    lease_expires_ts = NULL, lease_token_hash = NULL
                WHERE id = ? AND assigned_to = ?
                  AND attempt_count = ? AND lease_token_hash = ?
                  AND status IN ('claimed', 'in_progress')
                  AND lease_expires_ts > ?
                """,
                (
                    now,
                    _canonical_json(finding_ids, field="task result findings"),
                    now,
                    task_id,
                    self.agent_id,
                    attempt_number,
                    _hash_secret(token),
                    now,
                ),
            )
            if cursor.rowcount == 0:
                raise LeaseFenceError(
                    f"Task {task_id} is not active under attempt {attempt_number} "
                    f"for {self.agent_id}"
                )
            for finding_id in finding_ids:
                self._conn.execute(
                    """
                    INSERT OR IGNORE INTO task_result_findings(task_id, finding_id)
                    VALUES (?, ?)
                    """,
                    (task_id, finding_id),
                )
            row = self._conn.execute(
                "SELECT * FROM tasks WHERE id = ?",
                (task_id,),
            ).fetchone()
            self._event_locked(
                "task_completed",
                "task",
                task_id,
                campaign_id=row["campaign_id"],
                run_id=row["run_id"],
                task_id=task_id,
                payload={
                    "attempt": attempt_number,
                    "result_finding_ids": finding_ids,
                },
            )
        self._lease_tokens.pop(task_id, None)
        result = _row_to_dict(row)
        assert result is not None
        return result

    def fail_task(
        self,
        task_id: str,
        *,
        error: str | None = None,
        retryable: bool = False,
        retry_delay_seconds: int = 0,
        lease_token: str | None = None,
        attempt: int | None = None,
    ) -> dict[str, Any]:
        token, attempt_number = self._resolve_lease(task_id, lease_token, attempt)
        if retry_delay_seconds < 0:
            raise ValueError("retry_delay_seconds must be non-negative")
        with self._lock, self._transaction():
            now = time.time()
            row = self._conn.execute(
                "SELECT * FROM tasks WHERE id = ?",
                (task_id,),
            ).fetchone()
            if (
                row is None
                or row["assigned_to"] != self.agent_id
                or row["attempt_count"] != attempt_number
                or row["lease_token_hash"] != _hash_secret(token)
                or row["status"] not in {"claimed", "in_progress"}
                or row["lease_expires_ts"] is None
                or row["lease_expires_ts"] <= now
            ):
                raise LeaseFenceError(
                    f"Task {task_id} is not active under attempt {attempt_number} "
                    f"for {self.agent_id}"
                )
            can_retry = retryable and row["attempt_count"] < row["max_attempts"]
            status = "retry_wait" if can_retry else "failed"
            completed_ts = None if can_retry else now
            retry_after = now + retry_delay_seconds if can_retry else None
            assigned_to = row["reserved_for"] if can_retry else row["assigned_to"]
            self._conn.execute(
                """
                UPDATE tasks
                SET status = ?, assigned_to = ?, completed_ts = ?,
                    heartbeat_ts = ?, lease_expires_ts = NULL,
                    lease_token_hash = NULL, retry_after_ts = ?, last_error = ?
                WHERE id = ?
                """,
                (
                    status,
                    assigned_to,
                    completed_ts,
                    now,
                    retry_after,
                    error,
                    task_id,
                ),
            )
            self._event_locked(
                "task_retry_scheduled" if can_retry else "task_failed",
                "task",
                task_id,
                campaign_id=row["campaign_id"],
                run_id=row["run_id"],
                task_id=task_id,
                payload={
                    "attempt": attempt_number,
                    "error": error,
                    "retry_after_ts": retry_after,
                },
            )
            updated = self._conn.execute(
                "SELECT * FROM tasks WHERE id = ?",
                (task_id,),
            ).fetchone()
        self._lease_tokens.pop(task_id, None)
        result = _row_to_dict(updated)
        assert result is not None
        return result

    def recover_expired_tasks(self, *, now: float | None = None) -> dict[str, list[str]]:
        with self._lock, self._transaction():
            recovered, exhausted = self._recover_expired_locked(
                time.time() if now is None else now
            )
        return {"recovered": recovered, "exhausted": exhausted}

    def _recover_expired_locked(self, now: float) -> tuple[list[str], list[str]]:
        rows = self._conn.execute(
            """
            SELECT * FROM tasks
            WHERE status IN ('claimed', 'in_progress')
              AND lease_expires_ts IS NOT NULL
              AND lease_expires_ts <= ?
            ORDER BY lease_expires_ts, id
            """,
            (now,),
        ).fetchall()
        recovered: list[str] = []
        exhausted: list[str] = []
        for row in rows:
            if row["attempt_count"] < row["max_attempts"]:
                status = "pending"
                assigned_to = row["reserved_for"]
                completed_ts = None
                recovered.append(row["id"])
                event_type = "task_lease_recovered"
            else:
                status = "failed"
                assigned_to = row["assigned_to"]
                completed_ts = now
                exhausted.append(row["id"])
                event_type = "task_attempts_exhausted"
            self._conn.execute(
                """
                UPDATE tasks
                SET status = ?, assigned_to = ?, completed_ts = ?,
                    lease_expires_ts = NULL, lease_token_hash = NULL,
                    retry_after_ts = NULL,
                    last_error = COALESCE(last_error, 'worker lease expired')
                WHERE id = ?
                """,
                (status, assigned_to, completed_ts, row["id"]),
            )
            self._event_locked(
                event_type,
                "task",
                row["id"],
                campaign_id=row["campaign_id"],
                run_id=row["run_id"],
                task_id=row["id"],
                payload={"attempt": row["attempt_count"]},
                actor_agent_id="system",
            )
        return recovered, exhausted

    def _promote_retries_locked(self, now: float) -> None:
        rows = self._conn.execute(
            """
            SELECT * FROM tasks
            WHERE status = 'retry_wait'
              AND (retry_after_ts IS NULL OR retry_after_ts <= ?)
              AND attempt_count < max_attempts
            """,
            (now,),
        ).fetchall()
        for row in rows:
            self._conn.execute(
                """
                UPDATE tasks
                SET status = 'pending', assigned_to = reserved_for,
                    retry_after_ts = NULL
                WHERE id = ?
                """,
                (row["id"],),
            )
            self._event_locked(
                "task_retry_ready",
                "task",
                row["id"],
                campaign_id=row["campaign_id"],
                run_id=row["run_id"],
                task_id=row["id"],
                payload={"attempt_count": row["attempt_count"]},
                actor_agent_id="system",
            )

    def _validate_dependencies_locked(self, dependency_ids: list[str]) -> None:
        self._validate_ids_locked("tasks", dependency_ids, label="task dependencies")

    def _dependencies_complete_locked(self, task_id: str) -> bool:
        return (
            self._conn.execute(
                """
                SELECT 1
                FROM task_dependencies AS d
                JOIN tasks AS parent ON parent.id = d.depends_on_task_id
                WHERE d.task_id = ? AND parent.status <> 'complete'
                LIMIT 1
                """,
                (task_id,),
            ).fetchone()
            is None
        )

    def _get_task(self, task_id: str) -> dict[str, Any]:
        row = self._conn.execute(
            "SELECT * FROM tasks WHERE id = ?",
            (task_id,),
        ).fetchone()
        result = _row_to_dict(row)
        if result is None:
            raise KeyError(f"Task {task_id} not found")
        return result

    def get_task(self, task_id: str) -> dict[str, Any]:
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
        config: dict | None = None,
        idempotency_key: str | None = None,
        parent_run_id: str | None = None,
        campaign_id: str | None = None,
    ) -> str:
        with self._lock:
            self._validate_campaign_locked(campaign_id)
        return self.ledger.create_run(
            run_type,
            dataset_path,
            created_by=self.agent_id,
            config=config,
            idempotency_key=idempotency_key,
            parent_run_id=parent_run_id,
            campaign_id=campaign_id,
        )

    def start_run(self, run_id: str) -> dict[str, Any]:
        return self.ledger.start_run(run_id)

    def submit_run(self, run_id: str) -> dict[str, Any]:
        return self.ledger.submit_run(run_id)

    def wait_for_scheduler(self, run_id: str) -> dict[str, Any]:
        return self.ledger.wait_for_scheduler(run_id)

    def heartbeat_run(self, run_id: str) -> dict[str, Any]:
        return self.ledger.heartbeat_run(run_id)

    def complete_run(self, run_id: str, result: dict | None = None) -> dict[str, Any]:
        _canonical_json(result or {}, field="run result")
        return self.ledger.complete_run(run_id, result=result)

    def fail_run(self, run_id: str, error: str) -> dict[str, Any]:
        return self.ledger.fail_run(run_id, _bounded_text(error, field="run error"))

    def get_run(self, run_id: str) -> dict[str, Any]:
        return self.ledger.get_run(run_id)

    def list_runs(self, **filters: Any) -> list[dict[str, Any]]:
        return self.ledger.list_runs(**filters)

    def list_run_stages(
        self,
        run_id: str,
        *,
        stage_id: str | None = None,
    ) -> list[dict[str, Any]]:
        return self.ledger.list_stages(run_id, stage_id=stage_id)

    def run_events(
        self,
        run_id: str,
        *,
        after_id: int = 0,
        limit: int = MAX_LIST_LIMIT,
    ) -> list[dict[str, Any]]:
        self.ledger.get_run(run_id)
        return self.ledger.events(run_id, after_id=after_id, limit=_validate_limit(limit))

    def recover_stale_runs(
        self,
        *,
        stale_after_seconds: int = 300,
        now: float | None = None,
    ) -> dict[str, list[Any]]:
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
        referenced_finding_ids: list[str] | None = None,
        referenced_hypothesis_ids: list[str] | None = None,
        decisions_made: dict | None = None,
        *,
        campaign_id: str | None = None,
        idempotency_key: str | None = None,
    ) -> str:
        if not action_type.strip():
            raise ValueError("action_type must be non-empty")
        _validate_idempotency_key(idempotency_key)
        reasoning = _bounded_text(reasoning, field="coordinator reasoning")
        decisions_json = _canonical_json(decisions_made or {}, field="coordinator decisions")
        with self._lock, self._transaction():
            self._validate_campaign_locked(campaign_id)
            finding_ids = self._validate_ids_locked(
                "findings",
                referenced_finding_ids,
                label="referenced findings",
            )
            hypothesis_ids = self._validate_ids_locked(
                "hypotheses",
                referenced_hypothesis_ids,
                label="referenced hypotheses",
            )
            findings_json = _canonical_json(finding_ids, field="log finding references")
            hypotheses_json = _canonical_json(
                hypothesis_ids,
                field="log hypothesis references",
            )
            if idempotency_key:
                existing = self._conn.execute(
                    """
                    SELECT * FROM coordinator_log
                    WHERE actor_agent_id = ? AND idempotency_key = ?
                    """,
                    (self.agent_id, idempotency_key),
                ).fetchone()
                if existing is not None:
                    requested = (
                        action_type,
                        reasoning,
                        findings_json,
                        hypotheses_json,
                        decisions_json,
                        campaign_id,
                    )
                    current = (
                        existing["action_type"],
                        existing["reasoning"],
                        existing["referenced_finding_ids"],
                        existing["referenced_hypothesis_ids"],
                        existing["decisions_made"],
                        existing["campaign_id"],
                    )
                    if requested != current:
                        raise IdempotencyConflictError(
                            "Log idempotency key already exists with a different payload"
                        )
                    return str(existing["id"])
            log_id = str(uuid.uuid4())
            self._conn.execute(
                """
                INSERT INTO coordinator_log(
                    id, ts, action_type, reasoning, referenced_finding_ids,
                    referenced_hypothesis_ids, decisions_made, actor_agent_id,
                    campaign_id, idempotency_key
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    log_id,
                    time.time(),
                    action_type,
                    reasoning,
                    findings_json,
                    hypotheses_json,
                    decisions_json,
                    self.agent_id,
                    campaign_id,
                    idempotency_key,
                ),
            )
            for finding_id in finding_ids:
                self._conn.execute(
                    """
                    INSERT INTO coordinator_log_findings(log_id, finding_id)
                    VALUES (?, ?)
                    """,
                    (log_id, finding_id),
                )
            for hypothesis_id in hypothesis_ids:
                self._conn.execute(
                    """
                    INSERT INTO coordinator_log_hypotheses(log_id, hypothesis_id)
                    VALUES (?, ?)
                    """,
                    (log_id, hypothesis_id),
                )
            event_id = self._event_locked(
                "coordinator_log_created",
                "coordinator_log",
                log_id,
                campaign_id=campaign_id,
                payload={"action_type": action_type},
            )
            self._conn.execute(
                "UPDATE coordinator_log SET created_event_id = ? WHERE id = ?",
                (event_id, log_id),
            )
            return log_id

    def recent_log(
        self,
        limit: int = 20,
        since: float = 0,
        *,
        campaign_id: str | None = None,
        before_ts: float | None = None,
    ) -> list[dict[str, Any]]:
        clauses = ["ts > ?"]
        params: list[Any] = [since]
        if campaign_id is not None:
            clauses.append("campaign_id = ?")
            params.append(campaign_id)
        if before_ts is not None:
            clauses.append("ts < ?")
            params.append(before_ts)
        params.append(_validate_limit(limit))
        with self._lock:
            rows = self._conn.execute(
                f"""
                SELECT * FROM coordinator_log
                WHERE {' AND '.join(clauses)}
                ORDER BY ts DESC, id DESC LIMIT ?
                """,
                params,
            ).fetchall()
        return [_row_to_dict(row) for row in rows]  # type: ignore[misc]

    def get_log_entry(self, entry_id: str) -> dict[str, Any]:
        with self._lock:
            row = self._conn.execute(
                "SELECT * FROM coordinator_log WHERE id = ?",
                (entry_id,),
            ).fetchone()
        result = _row_to_dict(row)
        if result is None:
            raise KeyError(f"Coordinator log entry {entry_id} not found")
        return result

    # ------------------------------------------------------------------
    # Durable events, diagnostics, and backups
    # ------------------------------------------------------------------

    def events(
        self,
        *,
        after_id: int = 0,
        campaign_id: str | None = None,
        entity_type: str | None = None,
        limit: int = MAX_LIST_LIMIT,
    ) -> list[dict[str, Any]]:
        clauses = ["id > ?"]
        params: list[Any] = [max(0, int(after_id))]
        if campaign_id is not None:
            clauses.append("campaign_id = ?")
            params.append(campaign_id)
        if entity_type is not None:
            clauses.append("entity_type = ?")
            params.append(entity_type)
        params.append(_validate_limit(limit))
        with self._lock:
            rows = self._conn.execute(
                f"""
                SELECT * FROM coordination_events
                WHERE {' AND '.join(clauses)}
                ORDER BY id LIMIT ?
                """,
                params,
            ).fetchall()
        return [_row_to_dict(row) for row in rows]  # type: ignore[misc]

    def latest_event_id(self) -> int:
        with self._lock:
            return int(
                self._conn.execute(
                    "SELECT COALESCE(MAX(id), 0) FROM coordination_events"
                ).fetchone()[0]
            )

    def integrity_check(self) -> dict[str, Any]:
        with self._lock:
            rows = [
                row[0] for row in self._conn.execute("PRAGMA integrity_check").fetchall()
            ]
            foreign_keys = [
                dict(row)
                for row in self._conn.execute("PRAGMA foreign_key_check").fetchall()
            ]
        return {
            "ok": rows == ["ok"] and not foreign_keys,
            "integrity": rows,
            "foreign_key_violations": foreign_keys,
        }

    def backup(self, destination: str | Path) -> dict[str, Any]:
        destination_path = Path(destination).expanduser().resolve()
        if destination_path == self.db_path:
            raise ValueError("Backup destination must differ from the live database")
        destination_path.parent.mkdir(parents=True, exist_ok=True)
        target = sqlite3.connect(str(destination_path))
        try:
            with self._lock:
                self._conn.backup(target)
            target.execute("PRAGMA wal_checkpoint(TRUNCATE)")
            integrity = target.execute("PRAGMA integrity_check").fetchone()[0]
        finally:
            target.close()
        if integrity != "ok":
            raise RuntimeError(f"Backup integrity check failed: {integrity}")
        return {
            "path": str(destination_path),
            "size_bytes": destination_path.stat().st_size,
            "integrity": integrity,
            "created_ts": time.time(),
        }

    def stats(self) -> dict[str, Any]:
        with self._lock:
            tables = (
                "campaigns",
                "agents",
                "findings",
                "hypotheses",
                "tasks",
                "coordinator_log",
                "runs",
                "run_events",
                "run_stages",
                "coordination_events",
                "artifacts",
            )
            counts = {
                table: int(
                    self._conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
                )
                for table in tables
            }
            grouped = {}
            for key, table, column in (
                ("findings_by_novelty", "findings", "novelty"),
                ("tasks_by_status", "tasks", "status"),
                ("hypotheses_by_status", "hypotheses", "status"),
                ("runs_by_status", "runs", "status"),
                ("stages_by_status", "run_stages", "status"),
                ("agents_by_status", "agents", "status"),
            ):
                grouped[key] = {
                    str(row[0]): int(row[1])
                    for row in self._conn.execute(
                        f"SELECT {column}, COUNT(*) FROM {table} GROUP BY {column}"
                    ).fetchall()
                }
            now = time.time()
            oldest = self._conn.execute(
                """
                SELECT MIN(ts) FROM tasks
                WHERE status IN ('pending', 'retry_wait')
                """
            ).fetchone()[0]
            active_leases = int(
                self._conn.execute(
                    """
                    SELECT COUNT(*) FROM tasks
                    WHERE status IN ('claimed', 'in_progress')
                      AND lease_expires_ts > ?
                    """,
                    (now,),
                ).fetchone()[0]
            )
            expired_leases = int(
                self._conn.execute(
                    """
                    SELECT COUNT(*) FROM tasks
                    WHERE status IN ('claimed', 'in_progress')
                      AND lease_expires_ts <= ?
                    """,
                    (now,),
                ).fetchone()[0]
            )
            failed_row = self._conn.execute(
                """
                SELECT COUNT(*) AS total,
                       COALESCE(SUM(attempt_count >= max_attempts), 0) AS exhausted,
                       MIN(completed_ts) AS oldest_completed_ts
                FROM tasks
                WHERE status = 'failed'
                """
            ).fetchone()
            page_count = int(self._conn.execute("PRAGMA page_count").fetchone()[0])
            page_size = int(self._conn.execute("PRAGMA page_size").fetchone()[0])
            schema_version = int(
                self._conn.execute(
                    "SELECT COALESCE(MAX(version), 0) FROM ops_schema_meta"
                ).fetchone()[0]
            )
            active_resources: dict[str, dict[str, int]] = {}
            for row in self._conn.execute(
                """
                SELECT assigned_to, resource_request
                FROM tasks
                WHERE status IN ('claimed', 'in_progress')
                  AND assigned_to IS NOT NULL
                """
            ).fetchall():
                usage = active_resources.setdefault(
                    str(row["assigned_to"]),
                    {
                        "active_tasks": 0,
                        "cpu_slots": 0,
                        "memory_mb": 0,
                        "accelerator_slots": 0,
                    },
                )
                usage["active_tasks"] += 1
                request = _normalize_resources(
                    json.loads(row["resource_request"] or "{}")
                )
                for key, amount in request.items():
                    usage[key] += amount
        wal_path = Path(f"{self.db_path}-wal")
        shm_path = Path(f"{self.db_path}-shm")
        database = {
            "path": str(self.db_path),
            "schema_version": schema_version,
            "page_bytes": page_count * page_size,
            "file_bytes": self.db_path.stat().st_size if self.db_path.exists() else 0,
            "wal_bytes": wal_path.stat().st_size if wal_path.exists() else 0,
            "shm_bytes": shm_path.stat().st_size if shm_path.exists() else 0,
            "fts5": self.fts_available,
        }
        return {
            "counts": counts,
            **grouped,
            "queue": {
                "oldest_age_seconds": max(0.0, now - oldest) if oldest is not None else 0.0,
                "active_leases": active_leases,
                "expired_leases": expired_leases,
            },
            "dead_letters": {
                "count": int(failed_row["total"]),
                "attempts_exhausted": int(failed_row["exhausted"]),
                "oldest_age_seconds": (
                    max(0.0, now - float(failed_row["oldest_completed_ts"]))
                    if failed_row["oldest_completed_ts"] is not None
                    else 0.0
                ),
            },
            "active_resources_by_agent": active_resources,
            "database": database,
        }

    def close(self) -> None:
        if self._closed:
            return
        self._closed = True
        self.ledger.close()
        if self._close_callback is not None:
            self._close_callback(self._conn)
            self._close_callback = None
        elif self._owns_connection:
            self._conn.close()
        if self._direct_lock is not None:
            self._direct_lock.release()
            self._direct_lock = None

    def __enter__(self) -> OpsStore:
        return self

    def __exit__(self, *_exc_info: object) -> None:
        self.close()

    def __repr__(self) -> str:
        return f"OpsStore(db='{self.db_path}', agent_id='{self.agent_id}')"


__all__ = [
    "AGENT_ROLES",
    "DEFAULT_LEASE_SECONDS",
    "MAX_INLINE_JSON_BYTES",
    "IdempotencyConflictError",
    "LeaseFenceError",
    "OpsStore",
]
