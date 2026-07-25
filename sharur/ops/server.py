"""HTTP control plane for coordinated Sharur agents.

One server process owns the SQLite database. Remote workers communicate only
through this API; this ownership rule is especially important when the
database path is on NFS or another network filesystem.
"""

from __future__ import annotations

import argparse
import asyncio
import contextlib
import json
import logging
import os
import secrets
import sqlite3
import threading
import time
from collections import deque
from contextlib import asynccontextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import Annotated, Any, Literal

from fastapi import APIRouter, Depends, FastAPI, HTTPException, Query, Request
from fastapi.responses import JSONResponse, PlainTextResponse, StreamingResponse
from fastapi.security import HTTPAuthorizationCredentials, HTTPBearer
from pydantic import BaseModel, ConfigDict, Field

from sharur import __version__
from sharur.ops.db import SQLiteConnectionPool, SQLiteServerLock
from sharur.ops.schema import DEFAULT_LEASE_SECONDS
from sharur.ops.store import (
    IdempotencyConflictError,
    LeaseFenceError,
    OpsStore,
)
from sharur.review.controller import ReviewController
from sharur.review.metrics import review_campaign_metrics
from sharur.review.models import ReviewPolicy, load_review_policy
from sharur.review.reducer import ExactCandidateReducer


DEFAULT_HOST = "127.0.0.1"
DEFAULT_PORT = 8811
DEFAULT_SSE_QUEUE_SIZE = 1_000
DEFAULT_POOL_SIZE = 4
DEFAULT_MAX_REQUEST_BYTES = 1024 * 1024
DEFAULT_MAINTENANCE_INTERVAL_SECONDS = 15
DEFAULT_RUN_STALE_AFTER_SECONDS = 300


def _env_flag(name: str, *, default: bool = False) -> bool:
    value = os.environ.get(name)
    if value is None:
        return default
    return value.strip().lower() in {"1", "true", "yes", "on"}


def _is_loopback(host: str | None) -> bool:
    if host is None:
        return False
    return host.strip().lower() in {
        "127.0.0.1",
        "::1",
        "localhost",
        "testclient",
    }


@dataclass(frozen=True)
class Principal:
    agent_id: str
    role: str
    bootstrap: bool = False


@dataclass(frozen=True)
class _Subscriber:
    queue: asyncio.Queue[str]
    loop: asyncio.AbstractEventLoop


class EventBus:
    """Bounded wake-up channel; durable events remain authoritative."""

    def __init__(self, *, queue_size: int = DEFAULT_SSE_QUEUE_SIZE):
        if queue_size < 1:
            raise ValueError("SSE queue size must be positive")
        self._queue_size = queue_size
        self._subscribers: list[_Subscriber] = []
        self._lock = threading.Lock()
        self._dropped_notifications = 0

    def subscribe(self) -> _Subscriber:
        subscriber = _Subscriber(
            queue=asyncio.Queue(maxsize=self._queue_size),
            loop=asyncio.get_running_loop(),
        )
        with self._lock:
            self._subscribers.append(subscriber)
        return subscriber

    def unsubscribe(self, subscriber: _Subscriber) -> None:
        with self._lock:
            self._subscribers = [
                current
                for current in self._subscribers
                if current is not subscriber
            ]

    def _offer(self, queue: asyncio.Queue[str], payload: str) -> None:
        if queue.full():
            try:
                queue.get_nowait()
            except asyncio.QueueEmpty:
                pass
            else:
                with self._lock:
                    self._dropped_notifications += 1
        with contextlib.suppress(asyncio.QueueFull):
            queue.put_nowait(payload)

    def publish(self, event_type: str = "changed", data: dict | None = None) -> None:
        payload = json.dumps(
            {"type": event_type, **(data or {})},
            separators=(",", ":"),
        )
        with self._lock:
            subscribers = list(self._subscribers)
        for subscriber in subscribers:
            with contextlib.suppress(RuntimeError):
                subscriber.loop.call_soon_threadsafe(
                    self._offer,
                    subscriber.queue,
                    payload,
                )

    def stats(self) -> dict[str, int]:
        with self._lock:
            return {
                "subscribers": len(self._subscribers),
                "dropped_wakeups": self._dropped_notifications,
            }


@dataclass
class OpsRuntime:
    db_path: Path
    api_token: str | None
    allow_insecure_remote: bool
    event_bus: EventBus
    pool_size: int
    max_request_bytes: int
    maintenance_interval_seconds: int
    run_stale_after_seconds: int
    backup_dir: Path | None
    backup_interval_seconds: int
    pool: SQLiteConnectionPool | None = None
    owner_lock: SQLiteServerLock | None = None
    maintenance_task: asyncio.Task[None] | None = None
    last_backup_ts: float = 0.0
    state_lock: threading.Lock = field(default_factory=threading.Lock)
    request_count: int = 0
    request_error_count: int = 0
    request_duration_seconds: float = 0.0
    request_durations: deque[float] = field(
        default_factory=lambda: deque(maxlen=4_096)
    )
    sqlite_busy_errors: int = 0
    backup_failure_count: int = 0
    last_backup_error: str | None = None
    request_duration_buckets: dict[float, int] = field(
        default_factory=lambda: {
            0.01: 0,
            0.025: 0,
            0.05: 0,
            0.1: 0,
            0.25: 0,
            0.5: 0,
            1.0: 0,
            2.5: 0,
            5.0: 0,
        }
    )


class StrictModel(BaseModel):
    model_config = ConfigDict(extra="forbid")


class CampaignIn(StrictModel):
    name: str = Field(min_length=1, max_length=256)
    description: str = Field(default="", max_length=65_536)
    dataset_path: str | None = Field(default=None, max_length=4_096)
    metadata: dict = Field(default_factory=dict)
    idempotency_key: str | None = Field(default=None, min_length=1, max_length=512)


class CampaignUpdate(StrictModel):
    status: Literal["active", "paused", "complete", "failed", "archived"]


class AgentIn(StrictModel):
    agent_id: str = Field(min_length=1, max_length=256)
    role: Literal["reader", "worker", "coordinator", "operator"] = "worker"
    capabilities: list[str] = Field(default_factory=list, max_length=256)
    max_concurrent_tasks: int = Field(default=1, ge=0, le=100_000)
    capacity_cpu_slots: int = Field(default=1, ge=0)
    capacity_memory_mb: int = Field(default=0, ge=0)
    capacity_accelerator_slots: int = Field(default=0, ge=0)
    host: str | None = Field(default=None, max_length=1_024)
    metadata: dict = Field(default_factory=dict)
    rotate_token: bool = False


class AgentHeartbeatIn(StrictModel):
    status: Literal["active", "draining", "offline"] = "active"
    host: str | None = Field(default=None, max_length=1_024)


class AgentStatusIn(StrictModel):
    status: Literal["active", "draining", "offline", "revoked"]


class FindingIn(StrictModel):
    agent_id: str | None = Field(default=None, min_length=1, max_length=256)
    finding_type: str = Field(min_length=1, max_length=256)
    domain: str = Field(min_length=1, max_length=256)
    summary: str = Field(max_length=262_144)
    evidence: dict = Field(default_factory=dict)
    confidence: float = Field(default=0.5, ge=0.0, le=1.0)
    novelty: int = Field(default=0, ge=0, le=3)
    parent_finding_id: str | None = None
    reasoning: str = Field(default="", max_length=262_144)
    campaign_id: str | None = None
    task_id: str | None = None
    idempotency_key: str | None = Field(default=None, min_length=1, max_length=512)
    schema_version: int = Field(default=1, ge=1)
    validation_status: str = Field(default="unreviewed", max_length=128)


class FindingLinkIn(StrictModel):
    related_finding_id: str = Field(min_length=1, max_length=256)
    relation: str = Field(min_length=1, max_length=256)


class UnitDispositionIn(StrictModel):
    agent_id: str | None = Field(default=None, min_length=1, max_length=256)
    campaign_id: str
    unit_id: str = Field(min_length=1, max_length=4_096)
    dataset_id: str = Field(min_length=1, max_length=256)
    genome_id: str = Field(min_length=1, max_length=4_096)
    coverage_hash: str = Field(min_length=1, max_length=256)
    candidate_count: int = Field(ge=0)
    disposition: Literal["clear", "candidate", "incomplete", "failed"]
    evidence_bundle_hash: str = Field(min_length=1, max_length=256)
    task_id: str | None = None
    reason_codes: list[str] = Field(default_factory=list, max_length=1_000)
    strata: dict = Field(default_factory=dict)
    provenance: dict = Field(default_factory=dict)
    supersedes_disposition_id: str | None = None
    idempotency_key: str = Field(min_length=1, max_length=512)
    schema_version: int = Field(default=1, ge=1)


class CandidateOccurrenceIn(StrictModel):
    agent_id: str | None = Field(default=None, min_length=1, max_length=256)
    campaign_id: str
    dataset_id: str = Field(min_length=1, max_length=256)
    unit_id: str = Field(min_length=1, max_length=4_096)
    genome_id: str = Field(min_length=1, max_length=4_096)
    candidate_type: str = Field(min_length=1, max_length=256)
    signature_schema: str = Field(min_length=1, max_length=256)
    signature: dict
    evidence: dict
    verification: list[dict] = Field(default_factory=list, max_length=10_000)
    subject_refs: dict
    task_id: str | None = None
    reason_codes: list[str] = Field(default_factory=list, max_length=1_000)
    uncertainty: dict = Field(default_factory=dict)
    reduction_features: dict = Field(default_factory=dict)
    provenance: dict = Field(default_factory=dict)
    evidence_bundle_hash: str | None = Field(
        default=None, min_length=1, max_length=256
    )
    idempotency_key: str = Field(min_length=1, max_length=512)
    schema_version: int = Field(default=1, ge=1)


class CandidateClusterIn(StrictModel):
    campaign_id: str
    dataset_id: str = Field(min_length=1, max_length=256)
    candidate_type: str = Field(min_length=1, max_length=256)
    signature_schema: str = Field(min_length=1, max_length=256)
    member_ids: list[str] = Field(min_length=1, max_length=100_000)
    reducer_name: str = Field(min_length=1, max_length=256)
    reducer_version: str = Field(min_length=1, max_length=256)
    reducer_config_hash: str = Field(min_length=1, max_length=256)
    summary: dict
    counts: dict
    roles: dict[str, str] = Field(default_factory=dict)
    logical_cluster_id: str | None = Field(default=None, max_length=512)
    version: int = Field(default=1, ge=1)
    idempotency_key: str = Field(min_length=1, max_length=512)
    schema_version: int = Field(default=1, ge=1)


class ClusterLineageIn(StrictModel):
    child_cluster_id: str
    relation: Literal["supersedes", "split_from", "merged_from", "refines"]


class ClusterFindingIn(StrictModel):
    finding_id: str
    relation: Literal["materializes", "supports", "counterexample"] = "materializes"


class FindingReviewIn(StrictModel):
    reviewer_agent_id: str | None = Field(
        default=None, min_length=1, max_length=256
    )
    campaign_id: str
    dataset_id: str = Field(min_length=1, max_length=256)
    review_tier: str = Field(min_length=1, max_length=256)
    execution_profile: str = Field(min_length=1, max_length=256)
    provider: str = Field(min_length=1, max_length=256)
    model: str = Field(min_length=1, max_length=256)
    reasoning_effort: str = Field(min_length=1, max_length=128)
    prompt_hash: str = Field(min_length=1, max_length=256)
    rubric_version: str = Field(min_length=1, max_length=256)
    input_bundle_hash: str = Field(min_length=1, max_length=256)
    verdict: Literal[
        "promote", "hold", "needs_data", "reject", "duplicate", "split"
    ]
    confidence: float = Field(ge=0.0, le=1.0)
    finding_id: str | None = None
    cluster_id: str | None = None
    unit_disposition_id: str | None = None
    task_id: str | None = None
    run_id: str | None = None
    reconstructed_observations: dict = Field(default_factory=dict)
    claim_assessment: dict = Field(default_factory=dict)
    verification_summary: dict = Field(default_factory=dict)
    discrepancies: list[dict] = Field(default_factory=list, max_length=10_000)
    proposed_tasks: list[dict] = Field(default_factory=list, max_length=100)
    blind_to_prior_scores: bool = True
    blind_to_other_reviews: bool = True
    parent_review_id: str | None = None
    idempotency_key: str = Field(min_length=1, max_length=512)
    schema_version: int = Field(default=1, ge=1)


class ReviewVerificationIn(StrictModel):
    agent_id: str | None = Field(default=None, min_length=1, max_length=256)
    claim_key: str = Field(min_length=1, max_length=1_024)
    engine: str = Field(min_length=1, max_length=256)
    specification: dict
    dataset_id: str = Field(min_length=1, max_length=256)
    expected: Any
    status: Literal["pending", "pass", "fail", "error", "skipped"]
    actual: Any = None
    executed_ts: float | None = None
    code_commit: str | None = Field(default=None, max_length=256)
    artifact_id: str | None = None
    error: str | None = Field(default=None, max_length=16_384)
    supersedes_verification_id: str | None = None
    idempotency_key: str = Field(min_length=1, max_length=512)


class PromotionDecisionIn(StrictModel):
    campaign_id: str
    decision: Literal[
        "promote",
        "hold",
        "needs_data",
        "reject",
        "duplicate",
        "split",
        "merge",
        "publish",
        "reopen",
    ]
    source_tier: str = Field(min_length=1, max_length=256)
    target_tier: str | None = Field(default=None, max_length=256)
    policy_name: str = Field(min_length=1, max_length=256)
    policy_version: str = Field(min_length=1, max_length=256)
    policy_hash: str = Field(min_length=1, max_length=256)
    rationale: str = Field(min_length=1, max_length=262_144)
    finding_id: str | None = None
    cluster_id: str | None = None
    review_ids: list[str] = Field(default_factory=list, max_length=1_000)
    created_task_ids: list[str] = Field(default_factory=list, max_length=1_000)
    audit_sample: bool = False
    audit_stratum: dict = Field(default_factory=dict)
    idempotency_key: str = Field(min_length=1, max_length=512)
    schema_version: int = Field(default=1, ge=1)


class CanonicalPublicationIn(StrictModel):
    campaign_id: str
    finding_id: str
    decision_id: str
    dataset_id: str = Field(min_length=1, max_length=256)
    canonical_uri: str = Field(min_length=1, max_length=16_384)
    canonical_record_id: str = Field(min_length=1, max_length=4_096)
    canonical_record_hash: str = Field(min_length=1, max_length=256)
    metadata: dict = Field(default_factory=dict)
    idempotency_key: str = Field(min_length=1, max_length=512)


class ReviewReduceIn(StrictModel):
    campaign_id: str
    dataset_id: str | None = None
    candidate_type: str | None = None
    batch_size: int = Field(default=1_000, ge=1, le=10_000)


class ReviewControllerTickIn(StrictModel):
    campaign_id: str
    policy: dict | None = None


class ArtifactIn(StrictModel):
    content_hash: str = Field(min_length=16, max_length=256)
    uri: str = Field(min_length=1, max_length=16_384)
    size_bytes: int = Field(ge=0)
    media_type: str = Field(default="application/octet-stream", max_length=512)
    campaign_id: str | None = None
    task_id: str | None = None
    run_id: str | None = None
    metadata: dict = Field(default_factory=dict)


class ArtifactAttachIn(StrictModel):
    artifact_id: str = Field(min_length=1)
    relation: str = Field(default="evidence", min_length=1, max_length=256)


class HypothesisIn(StrictModel):
    source_agent_id: str | None = Field(default=None, min_length=1, max_length=256)
    source_finding_ids: list[str] = Field(default_factory=list)
    hypothesis: str = Field(min_length=1, max_length=262_144)
    domains_relevant: list[str] = Field(default_factory=list)
    campaign_id: str | None = None
    idempotency_key: str | None = Field(default=None, min_length=1, max_length=512)
    schema_version: int = Field(default=1, ge=1)


class HypothesisUpdate(StrictModel):
    status: Literal[
        "proposed",
        "investigating",
        "supported",
        "refuted",
        "inconclusive",
    ] | None = None
    assigned_agent_id: str | None = None
    domains_relevant: list[str] | None = None
    evidence_for: list[str] | None = None
    evidence_against: list[str] | None = None
    resolution_summary: str | None = Field(default=None, max_length=262_144)


class TaskIn(StrictModel):
    created_by: str | None = Field(default=None, min_length=1, max_length=256)
    task_type: str = Field(min_length=1, max_length=256)
    description: str = Field(min_length=1, max_length=262_144)
    params: dict = Field(default_factory=dict)
    priority: int = Field(default=1, ge=0, le=3)
    domain_hint: str | None = Field(default=None, max_length=256)
    assigned_to: str | None = Field(default=None, max_length=256)
    run_id: str | None = None
    campaign_id: str | None = None
    idempotency_key: str | None = Field(default=None, min_length=1, max_length=512)
    depends_on: list[str] = Field(default_factory=list)
    required_capabilities: list[str] = Field(default_factory=list, max_length=256)
    resource_request: dict[str, int] = Field(default_factory=dict)
    max_attempts: int = Field(default=3, ge=1)
    lease_seconds: int = Field(default=DEFAULT_LEASE_SECONDS, ge=1)


class ClaimIn(StrictModel):
    agent_id: str | None = Field(default=None, min_length=1, max_length=256)
    lease_seconds: int | None = Field(default=None, ge=1)
    campaign_id: str | None = None
    task_types: list[str] | None = None


class LeaseHeartbeatIn(StrictModel):
    agent_id: str | None = Field(default=None, min_length=1, max_length=256)
    lease_token: str = Field(min_length=16, max_length=512)
    lease_attempt: int = Field(ge=1)
    lease_seconds: int | None = Field(default=None, ge=1)


class TaskUpdate(StrictModel):
    status: Literal["in_progress", "complete", "failed"]
    agent_id: str | None = Field(default=None, min_length=1, max_length=256)
    lease_token: str = Field(min_length=16, max_length=512)
    lease_attempt: int = Field(ge=1)
    result_finding_ids: list[str] | None = None
    lease_seconds: int | None = Field(default=None, ge=1)
    error: str | None = Field(default=None, max_length=262_144)
    retryable: bool = False
    retry_delay_seconds: int = Field(default=0, ge=0)


class TaskCheckpointIn(StrictModel):
    agent_id: str | None = Field(default=None, min_length=1, max_length=256)
    checkpoint_key: str = Field(min_length=1, max_length=256)
    cursor: str | None = Field(default=None, max_length=4_096)
    payload: dict = Field(default_factory=dict)
    lease_token: str = Field(min_length=16, max_length=512)
    lease_attempt: int = Field(ge=1)


class RunIn(StrictModel):
    created_by: str | None = Field(default=None, min_length=1, max_length=256)
    run_type: str = Field(min_length=1, max_length=256)
    dataset_path: str = Field(min_length=1, max_length=4_096)
    config: dict = Field(default_factory=dict)
    idempotency_key: str | None = Field(default=None, min_length=1, max_length=512)
    parent_run_id: str | None = None
    campaign_id: str | None = None


class RunUpdate(StrictModel):
    status: Literal["complete", "failed"]
    result: dict = Field(default_factory=dict)
    error: str | None = Field(default=None, max_length=262_144)


class CoordinatorLogIn(StrictModel):
    action_type: str = Field(min_length=1, max_length=256)
    reasoning: str = Field(max_length=262_144)
    referenced_finding_ids: list[str] = Field(default_factory=list)
    referenced_hypothesis_ids: list[str] = Field(default_factory=list)
    decisions_made: dict = Field(default_factory=dict)
    campaign_id: str | None = None
    idempotency_key: str | None = Field(default=None, min_length=1, max_length=512)


class RequestSizeLimitMiddleware:
    """Reject declared and streamed request bodies above the configured cap."""

    def __init__(self, app: Any, *, max_bytes: int):
        self.app = app
        self.max_bytes = max_bytes

    async def __call__(self, scope: dict, receive: Any, send: Any) -> None:
        if scope["type"] == "http":
            headers = dict(scope.get("headers") or [])
            raw_length = headers.get(b"content-length")
            if raw_length is not None:
                with contextlib.suppress(ValueError):
                    if int(raw_length) > self.max_bytes:
                        response = JSONResponse(
                            {"detail": f"Request body exceeds {self.max_bytes} bytes"},
                            status_code=413,
                        )
                        await response(scope, receive, send)
                        return

            # Buffer the bounded body before routing. This gives chunked
            # requests the same deterministic 413 response as requests with a
            # declared Content-Length.
            messages: list[dict] = []
            received = 0
            while True:
                message = await receive()
                messages.append(message)
                if message["type"] != "http.request":
                    break
                received += len(message.get("body", b""))
                if received > self.max_bytes:
                    response = JSONResponse(
                        {"detail": f"Request body exceeds {self.max_bytes} bytes"},
                        status_code=413,
                    )
                    await response(scope, receive, send)
                    return
                if not message.get("more_body", False):
                    break

            async def replay_receive():
                if messages:
                    return messages.pop(0)
                return await receive()

            await self.app(scope, replay_receive, send)
            return
        await self.app(scope, receive, send)


class RequestMetricsMiddleware:
    """Low-cardinality request totals and aggregate service time."""

    def __init__(self, app: Any, *, runtime: OpsRuntime):
        self.app = app
        self.runtime = runtime

    async def __call__(self, scope: dict, receive: Any, send: Any) -> None:
        if scope["type"] != "http":
            await self.app(scope, receive, send)
            return
        started = time.perf_counter()
        status = 500
        sqlite_busy = False

        async def observe_send(message: dict) -> None:
            nonlocal status
            if message["type"] == "http.response.start":
                status = int(message["status"])
            await send(message)

        try:
            await self.app(scope, receive, observe_send)
        except sqlite3.OperationalError as exc:
            sqlite_busy = any(
                marker in str(exc).lower()
                for marker in ("database is locked", "database is busy")
            )
            raise
        finally:
            duration = time.perf_counter() - started
            with self.runtime.state_lock:
                self.runtime.request_count += 1
                self.runtime.request_duration_seconds += duration
                self.runtime.request_durations.append(duration)
                if status >= 400:
                    self.runtime.request_error_count += 1
                if sqlite_busy:
                    self.runtime.sqlite_busy_errors += 1
                for upper_bound in self.runtime.request_duration_buckets:
                    if duration <= upper_bound:
                        self.runtime.request_duration_buckets[upper_bound] += 1


def _runtime(request: Request) -> OpsRuntime:
    return request.app.state.ops_runtime


def _store(request: Request, *, agent_id: str = "server") -> OpsStore:
    runtime = _runtime(request)
    if runtime.pool is None:
        raise RuntimeError("Sharur Ops database pool is unavailable")
    connection = runtime.pool.acquire()
    return OpsStore(
        runtime.db_path,
        agent_id=agent_id,
        connection=connection,
        initialize=False,
        close_callback=runtime.pool.release,
        transaction_wait_observer=runtime.pool.observe_sqlite_write_wait,
    )


_bearer = HTTPBearer(auto_error=False)


def _authorize(
    request: Request,
    credentials: Annotated[HTTPAuthorizationCredentials | None, Depends(_bearer)],
) -> Principal:
    runtime = _runtime(request)
    supplied = (
        credentials.credentials
        if credentials is not None and credentials.scheme.lower() == "bearer"
        else None
    )
    if (
        runtime.api_token is not None
        and supplied is not None
        and secrets.compare_digest(supplied, runtime.api_token)
    ):
        principal = Principal("bootstrap-operator", "operator", bootstrap=True)
        request.state.ops_principal = principal
        return principal

    if supplied is not None:
        with _store(request, agent_id="auth") as store:
            agent = store.authenticate_agent(supplied)
        if agent is not None:
            principal = Principal(str(agent["id"]), str(agent["role"]))
            request.state.ops_principal = principal
            return principal

    client_host = request.client.host if request.client is not None else None
    if runtime.api_token is None and (
        runtime.allow_insecure_remote or _is_loopback(client_host)
    ):
        principal = Principal("local-operator", "operator", bootstrap=True)
        request.state.ops_principal = principal
        return principal

    if runtime.api_token is None:
        raise HTTPException(
            status_code=403,
            detail=(
                "Unauthenticated Sharur Ops access is restricted to loopback; "
                "set SHARUR_OPS_TOKEN for remote access"
            ),
        )
    raise HTTPException(
        status_code=401,
        detail="Missing or invalid Sharur Ops bearer token",
        headers={"WWW-Authenticate": "Bearer"},
    )


def _principal(request: Request) -> Principal:
    return request.state.ops_principal


def _require(request: Request, *roles: str) -> Principal:
    principal = _principal(request)
    if principal.role not in roles:
        raise HTTPException(
            status_code=403,
            detail=f"Role {principal.role} lacks permission for this operation",
        )
    return principal


def _actor_id(request: Request, requested: str | None = None) -> str:
    principal = _principal(request)
    if requested is None or requested == principal.agent_id:
        return requested or principal.agent_id
    if principal.role == "operator":
        return requested
    raise HTTPException(
        status_code=403,
        detail=(
            f"Credential identity {principal.agent_id} cannot act as {requested}"
        ),
    )


def _notify(request: Request, event_type: str, entity_id: str | None = None) -> None:
    _runtime(request).event_bus.publish(event_type, {"id": entity_id})


def _conflict(exc: Exception) -> HTTPException:
    return HTTPException(status_code=409, detail=str(exc))


router = APIRouter(dependencies=[Depends(_authorize)])


@router.post("/campaigns", status_code=201)
def create_campaign(campaign: CampaignIn, request: Request):
    principal = _require(request, "coordinator", "operator")
    with _store(request, agent_id=principal.agent_id) as store:
        try:
            campaign_id = store.create_campaign(
                campaign.name,
                description=campaign.description,
                dataset_path=campaign.dataset_path,
                metadata=campaign.metadata,
                idempotency_key=campaign.idempotency_key,
            )
            result = store.get_campaign(campaign_id)
        except (ValueError, IdempotencyConflictError) as exc:
            raise _conflict(exc) from exc
    _notify(request, "campaign", campaign_id)
    return result


@router.get("/campaigns")
def list_campaigns(
    request: Request,
    status: str | None = None,
    before_ts: float | None = None,
    limit: int = Query(default=50, ge=1, le=500),
):
    with _store(request) as store:
        return store.list_campaigns(status=status, before_ts=before_ts, limit=limit)


@router.get("/campaigns/{campaign_id}")
def get_campaign(campaign_id: str, request: Request):
    with _store(request) as store:
        try:
            return store.get_campaign(campaign_id)
        except KeyError as exc:
            raise HTTPException(404, str(exc)) from exc


@router.patch("/campaigns/{campaign_id}")
def update_campaign(campaign_id: str, update: CampaignUpdate, request: Request):
    principal = _require(request, "coordinator", "operator")
    with _store(request, agent_id=principal.agent_id) as store:
        try:
            result = store.update_campaign(campaign_id, status=update.status)
        except KeyError as exc:
            raise HTTPException(404, str(exc)) from exc
    _notify(request, "campaign", campaign_id)
    return result


@router.post("/agents", status_code=201)
def register_agent(agent: AgentIn, request: Request):
    principal = _require(request, "operator")
    with _store(request, agent_id=principal.agent_id) as store:
        try:
            result = store.register_agent(
                agent.agent_id,
                role=agent.role,
                capabilities=agent.capabilities,
                max_concurrent_tasks=agent.max_concurrent_tasks,
                capacity_cpu_slots=agent.capacity_cpu_slots,
                capacity_memory_mb=agent.capacity_memory_mb,
                capacity_accelerator_slots=agent.capacity_accelerator_slots,
                host=agent.host,
                metadata=agent.metadata,
                rotate_token=agent.rotate_token,
            )
        except ValueError as exc:
            raise HTTPException(400, str(exc)) from exc
    _notify(request, "agent", agent.agent_id)
    return result


@router.get("/agents")
def list_agents(
    request: Request,
    status: str | None = None,
    limit: int = Query(default=50, ge=1, le=500),
):
    _require(request, "coordinator", "operator")
    with _store(request) as store:
        return store.list_agents(status=status, limit=limit)


@router.post("/agents/me/heartbeat")
def heartbeat_agent(heartbeat: AgentHeartbeatIn, request: Request):
    principal = _require(request, "worker", "coordinator", "operator")
    if principal.bootstrap:
        raise HTTPException(400, "Bootstrap credentials do not represent a registered agent")
    with _store(request, agent_id=principal.agent_id) as store:
        try:
            result = store.heartbeat_agent(status=heartbeat.status, host=heartbeat.host)
        except KeyError as exc:
            raise HTTPException(404, str(exc)) from exc
    _notify(request, "agent", principal.agent_id)
    return result


@router.patch("/agents/{agent_id}")
def update_agent_status(agent_id: str, update: AgentStatusIn, request: Request):
    principal = _require(request, "operator")
    with _store(request, agent_id=principal.agent_id) as store:
        try:
            result = store.set_agent_status(agent_id, update.status)
        except KeyError as exc:
            raise HTTPException(404, str(exc)) from exc
    _notify(request, "agent", agent_id)
    return result


@router.post("/findings", status_code=201)
def create_finding(finding: FindingIn, request: Request):
    _require(request, "worker", "coordinator", "operator")
    actor = _actor_id(request, finding.agent_id)
    with _store(request, agent_id=actor) as store:
        try:
            finding_id = store.finding(
                finding.finding_type,
                finding.domain,
                finding.summary,
                evidence=finding.evidence,
                confidence=finding.confidence,
                novelty=finding.novelty,
                parent_finding_id=finding.parent_finding_id,
                reasoning=finding.reasoning,
                campaign_id=finding.campaign_id,
                task_id=finding.task_id,
                idempotency_key=finding.idempotency_key,
                schema_version=finding.schema_version,
                validation_status=finding.validation_status,
            )
            created = store.get_finding(finding_id)
        except (ValueError, IdempotencyConflictError) as exc:
            raise _conflict(exc) from exc
    _notify(request, "finding", finding_id)
    return {"id": finding_id, "ts": created["ts"], "finding": created}


@router.get("/findings")
def list_findings(
    request: Request,
    since: float = 0,
    min_novelty: int = 0,
    finding_type: str | None = None,
    domain: str | None = None,
    agent_id: str | None = None,
    campaign_id: str | None = None,
    before_ts: float | None = None,
    limit: int = Query(default=50, ge=1, le=500),
):
    with _store(request) as store:
        return store.recent_findings(
            since=since,
            min_novelty=min_novelty,
            finding_type=finding_type,
            domain=domain,
            agent_id=agent_id,
            campaign_id=campaign_id,
            before_ts=before_ts,
            limit=limit,
        )


@router.get("/findings/search/{text}")
def search_findings(
    text: str,
    request: Request,
    campaign_id: str | None = None,
    limit: int = Query(default=20, ge=1, le=100),
):
    with _store(request) as store:
        return store.search_findings(text, limit=limit, campaign_id=campaign_id)


@router.get("/findings/{finding_id}")
def get_finding(finding_id: str, request: Request):
    with _store(request) as store:
        try:
            return store.get_finding(finding_id)
        except KeyError as exc:
            raise HTTPException(404, str(exc)) from exc


@router.post("/findings/{finding_id}/links", status_code=201)
def link_findings(finding_id: str, link: FindingLinkIn, request: Request):
    principal = _require(request, "worker", "coordinator", "operator")
    with _store(request, agent_id=principal.agent_id) as store:
        try:
            store.link_findings(
                finding_id,
                link.related_finding_id,
                relation=link.relation,
            )
        except ValueError as exc:
            raise _conflict(exc) from exc
    _notify(request, "finding", finding_id)
    return {"finding_id": finding_id, **link.model_dump()}


@router.post("/review/unit-dispositions", status_code=201)
def create_unit_disposition(record: UnitDispositionIn, request: Request):
    _require(request, "worker", "coordinator", "operator")
    actor = _actor_id(request, record.agent_id)
    with _store(request, agent_id=actor) as store:
        try:
            record_id = store.record_unit_disposition(
                campaign_id=record.campaign_id,
                unit_id=record.unit_id,
                dataset_id=record.dataset_id,
                genome_id=record.genome_id,
                coverage_hash=record.coverage_hash,
                candidate_count=record.candidate_count,
                disposition=record.disposition,
                evidence_bundle_hash=record.evidence_bundle_hash,
                task_id=record.task_id,
                reason_codes=record.reason_codes,
                strata=record.strata,
                provenance=record.provenance,
                supersedes_disposition_id=record.supersedes_disposition_id,
                idempotency_key=record.idempotency_key,
                schema_version=record.schema_version,
            )
            result = store.get_unit_disposition(record_id)
        except ValueError as exc:
            raise _conflict(exc) from exc
    _notify(request, "unit_disposition", record_id)
    return result


@router.get("/review/unit-dispositions")
def list_unit_dispositions(
    request: Request,
    campaign_id: str,
    disposition: str | None = None,
    active_only: bool = True,
    limit: int = Query(default=500, ge=1, le=1_000),
):
    with _store(request) as store:
        return store.list_unit_dispositions(
            campaign_id=campaign_id,
            disposition=disposition,
            active_only=active_only,
            limit=limit,
        )


@router.get("/review/unit-dispositions/{disposition_id}")
def get_unit_disposition(disposition_id: str, request: Request):
    with _store(request) as store:
        try:
            return store.get_unit_disposition(disposition_id)
        except KeyError as exc:
            raise HTTPException(404, str(exc)) from exc


@router.post("/review/candidates", status_code=201)
def create_candidate_occurrence(candidate: CandidateOccurrenceIn, request: Request):
    _require(request, "worker", "coordinator", "operator")
    actor = _actor_id(request, candidate.agent_id)
    with _store(request, agent_id=actor) as store:
        try:
            candidate_id = store.create_candidate_occurrence(
                campaign_id=candidate.campaign_id,
                dataset_id=candidate.dataset_id,
                unit_id=candidate.unit_id,
                genome_id=candidate.genome_id,
                candidate_type=candidate.candidate_type,
                signature_schema=candidate.signature_schema,
                signature=candidate.signature,
                evidence=candidate.evidence,
                verification=candidate.verification,
                subject_refs=candidate.subject_refs,
                task_id=candidate.task_id,
                reason_codes=candidate.reason_codes,
                uncertainty=candidate.uncertainty,
                reduction_features=candidate.reduction_features,
                provenance=candidate.provenance,
                evidence_bundle_hash=candidate.evidence_bundle_hash,
                idempotency_key=candidate.idempotency_key,
                schema_version=candidate.schema_version,
            )
            result = store.get_candidate_occurrence(candidate_id)
        except ValueError as exc:
            raise _conflict(exc) from exc
    _notify(request, "candidate_occurrence", candidate_id)
    return result


@router.get("/review/candidates")
def list_candidate_occurrences(
    request: Request,
    campaign_id: str,
    candidate_type: str | None = None,
    unclustered_only: bool = False,
    limit: int = Query(default=500, ge=1, le=1_000),
):
    with _store(request) as store:
        return store.list_candidate_occurrences(
            campaign_id=campaign_id,
            candidate_type=candidate_type,
            unclustered_only=unclustered_only,
            limit=limit,
        )


@router.get("/review/candidates/{candidate_id}")
def get_candidate_occurrence(candidate_id: str, request: Request):
    with _store(request) as store:
        try:
            return store.get_candidate_occurrence(candidate_id)
        except KeyError as exc:
            raise HTTPException(404, str(exc)) from exc


@router.post("/review/clusters", status_code=201)
def create_candidate_cluster(cluster: CandidateClusterIn, request: Request):
    principal = _require(request, "coordinator", "operator")
    with _store(request, agent_id=principal.agent_id) as store:
        try:
            cluster_id = store.create_candidate_cluster(**cluster.model_dump())
            result = store.get_candidate_cluster(cluster_id)
        except ValueError as exc:
            raise _conflict(exc) from exc
    _notify(request, "candidate_cluster", cluster_id)
    return result


@router.get("/review/clusters")
def list_candidate_clusters(
    request: Request,
    campaign_id: str,
    candidate_type: str | None = None,
    status: str | None = "active",
    limit: int = Query(default=500, ge=1, le=1_000),
):
    with _store(request) as store:
        return store.list_candidate_clusters(
            campaign_id=campaign_id,
            candidate_type=candidate_type,
            status=status,
            limit=limit,
        )


@router.get("/review/clusters/{cluster_id}")
def get_candidate_cluster(
    cluster_id: str,
    request: Request,
    member_limit: int = Query(default=100, ge=1, le=1_000),
):
    with _store(request) as store:
        try:
            return store.get_candidate_cluster(
                cluster_id,
                member_limit=member_limit,
            )
        except KeyError as exc:
            raise HTTPException(404, str(exc)) from exc


@router.get("/review/clusters/{cluster_id}/members")
def list_candidate_cluster_members(
    cluster_id: str,
    request: Request,
    after_candidate_id: str | None = None,
    limit: int = Query(default=500, ge=1, le=1_000),
):
    with _store(request) as store:
        try:
            return store.list_candidate_cluster_members(
                cluster_id,
                after_candidate_id=after_candidate_id,
                limit=limit,
            )
        except KeyError as exc:
            raise HTTPException(404, str(exc)) from exc


@router.post("/review/clusters/{cluster_id}/lineage", status_code=201)
def link_candidate_cluster(
    cluster_id: str,
    link: ClusterLineageIn,
    request: Request,
):
    principal = _require(request, "coordinator", "operator")
    with _store(request, agent_id=principal.agent_id) as store:
        try:
            if link.relation == "supersedes":
                store.supersede_candidate_cluster(
                    cluster_id,
                    link.child_cluster_id,
                )
            else:
                store.link_candidate_clusters(
                    cluster_id,
                    link.child_cluster_id,
                    relation=link.relation,
                )
        except ValueError as exc:
            raise _conflict(exc) from exc
    _notify(request, "candidate_cluster", link.child_cluster_id)
    return {
        "parent_cluster_id": cluster_id,
        "child_cluster_id": link.child_cluster_id,
        "relation": link.relation,
    }


@router.post("/review/clusters/{cluster_id}/findings", status_code=201)
def link_cluster_finding(
    cluster_id: str,
    link: ClusterFindingIn,
    request: Request,
):
    principal = _require(request, "worker", "coordinator", "operator")
    with _store(request, agent_id=principal.agent_id) as store:
        try:
            store.link_cluster_finding(
                cluster_id,
                link.finding_id,
                relation=link.relation,
            )
        except ValueError as exc:
            raise _conflict(exc) from exc
    _notify(request, "finding", link.finding_id)
    return {"cluster_id": cluster_id, **link.model_dump()}


@router.get("/review/cluster-findings")
def list_cluster_findings(
    request: Request,
    cluster_id: str | None = None,
    finding_id: str | None = None,
):
    with _store(request) as store:
        try:
            return store.list_cluster_findings(
                cluster_id=cluster_id,
                finding_id=finding_id,
            )
        except ValueError as exc:
            raise HTTPException(400, str(exc)) from exc


@router.post("/review/reviews", status_code=201)
def create_finding_review(review: FindingReviewIn, request: Request):
    _require(request, "worker", "coordinator", "operator")
    actor = _actor_id(request, review.reviewer_agent_id)
    payload = review.model_dump(exclude={"reviewer_agent_id"})
    with _store(request, agent_id=actor) as store:
        try:
            review_id = store.create_finding_review(**payload)
            result = store.get_finding_review(review_id)
        except ValueError as exc:
            raise _conflict(exc) from exc
    _notify(request, "finding_review", review_id)
    return result


@router.get("/review/reviews")
def list_finding_reviews(
    request: Request,
    campaign_id: str,
    finding_id: str | None = None,
    cluster_id: str | None = None,
    unit_disposition_id: str | None = None,
    review_tier: str | None = None,
    verdict: str | None = None,
    limit: int = Query(default=500, ge=1, le=1_000),
):
    with _store(request) as store:
        return store.list_finding_reviews(
            campaign_id=campaign_id,
            finding_id=finding_id,
            cluster_id=cluster_id,
            unit_disposition_id=unit_disposition_id,
            review_tier=review_tier,
            verdict=verdict,
            limit=limit,
        )


@router.get("/review/reviews/{review_id}")
def get_finding_review(review_id: str, request: Request):
    with _store(request) as store:
        try:
            return store.get_finding_review(review_id)
        except KeyError as exc:
            raise HTTPException(404, str(exc)) from exc


@router.post("/review/reviews/{review_id}/verifications", status_code=201)
def create_review_verification(
    review_id: str,
    verification: ReviewVerificationIn,
    request: Request,
):
    _require(request, "worker", "coordinator", "operator")
    actor = _actor_id(request, verification.agent_id)
    payload = verification.model_dump(exclude={"agent_id"})
    with _store(request, agent_id=actor) as store:
        try:
            verification_id = store.record_review_verification(
                review_id=review_id,
                **payload,
            )
        except ValueError as exc:
            raise _conflict(exc) from exc
    _notify(request, "review_verification", verification_id)
    return {"id": verification_id}


@router.get("/review/reviews/{review_id}/verifications")
def list_review_verifications(review_id: str, request: Request):
    with _store(request) as store:
        return store.list_review_verifications(review_id)


@router.post("/review/decisions", status_code=201)
def create_promotion_decision(decision: PromotionDecisionIn, request: Request):
    principal = _require(request, "coordinator", "operator")
    with _store(request, agent_id=principal.agent_id) as store:
        try:
            decision_id = store.create_promotion_decision(**decision.model_dump())
            result = store.get_promotion_decision(decision_id)
        except ValueError as exc:
            raise _conflict(exc) from exc
    _notify(request, "promotion_decision", decision_id)
    return result


@router.get("/review/decisions")
def list_promotion_decisions(
    request: Request,
    campaign_id: str,
    finding_id: str | None = None,
    cluster_id: str | None = None,
    decision: str | None = None,
    limit: int = Query(default=500, ge=1, le=1_000),
):
    with _store(request) as store:
        return store.list_promotion_decisions(
            campaign_id=campaign_id,
            finding_id=finding_id,
            cluster_id=cluster_id,
            decision=decision,
            limit=limit,
        )


@router.get("/review/decisions/{decision_id}")
def get_promotion_decision(decision_id: str, request: Request):
    with _store(request) as store:
        try:
            return store.get_promotion_decision(decision_id)
        except KeyError as exc:
            raise HTTPException(404, str(exc)) from exc


@router.post("/review/publications", status_code=201)
def create_canonical_publication(
    publication: CanonicalPublicationIn,
    request: Request,
):
    principal = _require(request, "coordinator", "operator")
    with _store(request, agent_id=principal.agent_id) as store:
        try:
            publication_id = store.record_canonical_publication(
                **publication.model_dump()
            )
        except ValueError as exc:
            raise _conflict(exc) from exc
    _notify(request, "canonical_publication", publication_id)
    return {"id": publication_id}


@router.get("/review/publications")
def list_canonical_publications(
    request: Request,
    campaign_id: str,
    finding_id: str | None = None,
    limit: int = Query(default=500, ge=1, le=1_000),
):
    with _store(request) as store:
        return store.list_canonical_publications(
            campaign_id=campaign_id,
            finding_id=finding_id,
            limit=limit,
        )


@router.post("/review/reduce")
def reduce_candidates(reduction: ReviewReduceIn, request: Request):
    principal = _require(request, "coordinator", "operator")
    with _store(request, agent_id=principal.agent_id) as store:
        try:
            result = ExactCandidateReducer().reduce_campaign(
                store,
                reduction.campaign_id,
                dataset_id=reduction.dataset_id,
                candidate_type=reduction.candidate_type,
                batch_size=reduction.batch_size,
            )
        except (KeyError, ValueError) as exc:
            raise _conflict(exc) from exc
    _notify(request, "candidate_reduction")
    return result.to_dict()


@router.post("/review/controller/tick")
def tick_review_controller(tick: ReviewControllerTickIn, request: Request):
    principal = _require(request, "coordinator", "operator")
    try:
        policy = (
            ReviewPolicy.model_validate(tick.policy)
            if tick.policy is not None
            else load_review_policy()
        )
    except ValueError as exc:
        raise HTTPException(400, str(exc)) from exc
    with _store(request, agent_id=principal.agent_id) as store:
        try:
            result = ReviewController(store, policy).tick(tick.campaign_id)
        except (KeyError, ValueError) as exc:
            raise _conflict(exc) from exc
    _notify(request, "review_controller")
    return result.to_dict()


@router.get("/review/status")
def get_review_status(campaign_id: str, request: Request):
    with _store(request) as store:
        try:
            return review_campaign_metrics(store, campaign_id)
        except KeyError as exc:
            raise HTTPException(404, str(exc)) from exc


@router.post("/artifacts", status_code=201)
def register_artifact(artifact: ArtifactIn, request: Request):
    principal = _require(request, "worker", "coordinator", "operator")
    with _store(request, agent_id=principal.agent_id) as store:
        try:
            artifact_id = store.register_artifact(
                artifact.content_hash,
                artifact.uri,
                artifact.size_bytes,
                media_type=artifact.media_type,
                campaign_id=artifact.campaign_id,
                task_id=artifact.task_id,
                run_id=artifact.run_id,
                metadata=artifact.metadata,
            )
            result = store.get_artifact(artifact_id)
        except (ValueError, IdempotencyConflictError) as exc:
            raise _conflict(exc) from exc
    _notify(request, "artifact", artifact_id)
    return result


@router.get("/artifacts/{artifact_id}")
def get_artifact(artifact_id: str, request: Request):
    with _store(request) as store:
        try:
            return store.get_artifact(artifact_id)
        except KeyError as exc:
            raise HTTPException(404, str(exc)) from exc


@router.post("/findings/{finding_id}/artifacts", status_code=201)
def attach_artifact(finding_id: str, link: ArtifactAttachIn, request: Request):
    principal = _require(request, "worker", "coordinator", "operator")
    with _store(request, agent_id=principal.agent_id) as store:
        try:
            store.attach_artifact(
                finding_id,
                link.artifact_id,
                relation=link.relation,
            )
        except ValueError as exc:
            raise _conflict(exc) from exc
    _notify(request, "finding", finding_id)
    return {"finding_id": finding_id, **link.model_dump()}


@router.post("/hypotheses", status_code=201)
def create_hypothesis(hypothesis: HypothesisIn, request: Request):
    _require(request, "worker", "coordinator", "operator")
    actor = _actor_id(request, hypothesis.source_agent_id)
    with _store(request, agent_id=actor) as store:
        try:
            hypothesis_id = store.hypothesis(
                hypothesis.hypothesis,
                source_finding_ids=hypothesis.source_finding_ids,
                domains_relevant=hypothesis.domains_relevant,
                campaign_id=hypothesis.campaign_id,
                idempotency_key=hypothesis.idempotency_key,
                schema_version=hypothesis.schema_version,
            )
            created = store.get_hypothesis(hypothesis_id)
        except (ValueError, IdempotencyConflictError) as exc:
            raise _conflict(exc) from exc
    _notify(request, "hypothesis", hypothesis_id)
    return {"id": hypothesis_id, "ts": created["ts"], "hypothesis": created}


@router.get("/hypotheses")
def list_hypotheses(
    request: Request,
    status: str | None = None,
    unassigned: bool = False,
    campaign_id: str | None = None,
    before_ts: float | None = None,
    limit: int = Query(default=50, ge=1, le=200),
):
    with _store(request) as store:
        return store.list_hypotheses(
            status=status,
            unassigned=unassigned,
            campaign_id=campaign_id,
            before_ts=before_ts,
            limit=limit,
        )


@router.patch("/hypotheses/{hypothesis_id}")
def update_hypothesis(
    hypothesis_id: str,
    update: HypothesisUpdate,
    request: Request,
):
    principal = _require(request, "worker", "coordinator", "operator")
    values = update.model_dump(exclude_none=True)
    if not values:
        raise HTTPException(400, "No hypothesis updates supplied")
    with _store(request, agent_id=principal.agent_id) as store:
        try:
            result = store.update_hypothesis(hypothesis_id, **values)
        except KeyError as exc:
            raise HTTPException(404, str(exc)) from exc
        except ValueError as exc:
            raise _conflict(exc) from exc
    _notify(request, "hypothesis", hypothesis_id)
    return result


@router.post("/tasks", status_code=201)
def create_task(task: TaskIn, request: Request):
    _require(request, "coordinator", "operator")
    actor = _actor_id(request, task.created_by)
    with _store(request, agent_id=actor) as store:
        try:
            task_id = store.create_task(
                task.task_type,
                task.description,
                params=task.params,
                priority=task.priority,
                domain_hint=task.domain_hint,
                assigned_to=task.assigned_to,
                run_id=task.run_id,
                campaign_id=task.campaign_id,
                idempotency_key=task.idempotency_key,
                depends_on=task.depends_on,
                required_capabilities=task.required_capabilities,
                resource_request=task.resource_request,
                max_attempts=task.max_attempts,
                lease_seconds=task.lease_seconds,
            )
            created = store.get_task(task_id)
        except (ValueError, IdempotencyConflictError) as exc:
            raise _conflict(exc) from exc
    _notify(request, "task", task_id)
    return {"id": task_id, "ts": created["ts"], "task": created}


@router.get("/tasks")
def list_tasks(
    request: Request,
    status: str | None = None,
    assigned_to: str | None = None,
    unassigned: bool = False,
    campaign_id: str | None = None,
    before_ts: float | None = None,
    limit: int = Query(default=50, ge=1, le=200),
):
    with _store(request) as store:
        return store.list_tasks(
            status=status,
            assigned_to=assigned_to,
            unassigned=unassigned,
            campaign_id=campaign_id,
            before_ts=before_ts,
            limit=limit,
        )


@router.post("/tasks/claim-next")
def claim_next_task(claim: ClaimIn, request: Request):
    _require(request, "worker", "coordinator", "operator")
    actor = _actor_id(request, claim.agent_id)
    with _store(request, agent_id=actor) as store:
        try:
            task = store.claim_next_task(
                campaign_id=claim.campaign_id,
                task_types=claim.task_types,
                lease_seconds=claim.lease_seconds,
            )
        except ValueError as exc:
            raise _conflict(exc) from exc
    if task is not None:
        _notify(request, "task", task["id"])
    return task


@router.post("/tasks/recover")
def recover_tasks(request: Request):
    principal = _require(request, "coordinator", "operator")
    with _store(request, agent_id=principal.agent_id) as store:
        result = store.recover_expired_tasks()
    _notify(request, "task_recovery")
    return result


@router.post("/tasks/reset-failed")
def reset_failed_tasks(payload: ResetFailedIn, request: Request):
    principal = _require(request, "coordinator", "operator")
    with _store(request, agent_id=principal.agent_id) as store:
        result = store.reset_failed_tasks(
            campaign_id=payload.campaign_id,
            task_ids=payload.task_ids,
            only_transient=payload.only_transient,
            extra_attempts=payload.extra_attempts,
        )
    _notify(request, "task_recovery")
    return result


@router.patch("/tasks/{task_id}")
def update_task(task_id: str, update: TaskUpdate, request: Request):
    _require(request, "worker", "coordinator", "operator")
    actor = _actor_id(request, update.agent_id)
    with _store(request, agent_id=actor) as store:
        try:
            if update.status == "complete":
                result = store.complete_task(
                    task_id,
                    update.result_finding_ids,
                    lease_token=update.lease_token,
                    attempt=update.lease_attempt,
                )
            elif update.status == "failed":
                result = store.fail_task(
                    task_id,
                    error=update.error,
                    retryable=update.retryable,
                    retry_delay_seconds=update.retry_delay_seconds,
                    lease_token=update.lease_token,
                    attempt=update.lease_attempt,
                )
            else:
                result = store.heartbeat_task(
                    task_id,
                    lease_seconds=update.lease_seconds,
                    in_progress=True,
                    lease_token=update.lease_token,
                    attempt=update.lease_attempt,
                )
        except (KeyError, ValueError, LeaseFenceError) as exc:
            raise _conflict(exc) from exc
    _notify(request, "task", task_id)
    return result


@router.post("/tasks/{task_id}/claim")
def claim_task(
    task_id: str,
    request: Request,
    claim: ClaimIn | None = None,
    agent_id: str | None = None,
    lease_seconds: int | None = Query(default=None, ge=1),
):
    _require(request, "worker", "coordinator", "operator")
    requested = (claim.agent_id if claim else None) or agent_id
    actor = _actor_id(request, requested)
    duration = (claim.lease_seconds if claim else None) or lease_seconds
    with _store(request, agent_id=actor) as store:
        try:
            result = store.claim_task(task_id, lease_seconds=duration)
        except ValueError as exc:
            raise _conflict(exc) from exc
    _notify(request, "task", task_id)
    return result


@router.post("/tasks/{task_id}/heartbeat")
def heartbeat_task(
    task_id: str,
    heartbeat: LeaseHeartbeatIn,
    request: Request,
):
    _require(request, "worker", "coordinator", "operator")
    actor = _actor_id(request, heartbeat.agent_id)
    with _store(request, agent_id=actor) as store:
        try:
            result = store.heartbeat_task(
                task_id,
                lease_seconds=heartbeat.lease_seconds,
                lease_token=heartbeat.lease_token,
                attempt=heartbeat.lease_attempt,
            )
        except (ValueError, LeaseFenceError) as exc:
            raise _conflict(exc) from exc
    _notify(request, "task", task_id)
    return result


@router.put("/tasks/{task_id}/checkpoint")
def put_task_checkpoint(
    task_id: str,
    checkpoint: TaskCheckpointIn,
    request: Request,
):
    _require(request, "worker", "coordinator", "operator")
    actor = _actor_id(request, checkpoint.agent_id)
    with _store(request, agent_id=actor) as store:
        try:
            result = store.put_task_checkpoint(
                task_id,
                checkpoint.checkpoint_key,
                cursor=checkpoint.cursor,
                payload=checkpoint.payload,
                lease_token=checkpoint.lease_token,
                attempt=checkpoint.lease_attempt,
            )
        except (KeyError, ValueError, LeaseFenceError) as exc:
            raise _conflict(exc) from exc
    _notify(request, "task_checkpoint", task_id)
    return result


@router.get("/tasks/{task_id}/checkpoints")
def get_task_checkpoints(
    task_id: str,
    request: Request,
    checkpoint_key: str | None = Query(default=None, min_length=1, max_length=256),
    limit: int = Query(default=50, ge=1, le=1_000),
):
    _require(request, "worker", "coordinator", "operator")
    with _store(request) as store:
        try:
            if checkpoint_key is not None:
                result = store.get_task_checkpoint(task_id, checkpoint_key)
                if result is None:
                    raise HTTPException(404, "Task checkpoint was not found")
                return result
            return store.list_task_checkpoints(task_id, limit=limit)
        except KeyError as exc:
            raise HTTPException(404, str(exc)) from exc
        except ValueError as exc:
            raise _conflict(exc) from exc


@router.post("/runs", status_code=201)
def create_run(run: RunIn, request: Request):
    _require(request, "coordinator", "operator")
    actor = _actor_id(request, run.created_by)
    with _store(request, agent_id=actor) as store:
        try:
            run_id = store.create_run(
                run.run_type,
                run.dataset_path,
                config=run.config,
                idempotency_key=run.idempotency_key,
                parent_run_id=run.parent_run_id,
                campaign_id=run.campaign_id,
            )
            result = store.get_run(run_id)
        except ValueError as exc:
            raise _conflict(exc) from exc
    _notify(request, "run", run_id)
    return result


@router.get("/runs")
def list_runs(
    request: Request,
    dataset_path: str | None = None,
    run_type: str | None = None,
    status: str | None = None,
    campaign_id: str | None = None,
    before_ts: float | None = None,
    limit: int = Query(default=50, ge=1, le=200),
):
    with _store(request) as store:
        return store.list_runs(
            dataset_path=dataset_path,
            run_type=run_type,
            status=status,
            campaign_id=campaign_id,
            before_ts=before_ts,
            limit=limit,
        )


@router.post("/runs/recover")
def recover_runs(
    request: Request,
    stale_after_seconds: int = Query(default=300, ge=1),
):
    principal = _require(request, "coordinator", "operator")
    with _store(request, agent_id=principal.agent_id) as store:
        result = store.recover_stale_runs(
            stale_after_seconds=stale_after_seconds,
        )
    _notify(request, "run_recovery")
    return result


@router.get("/runs/{run_id}")
def get_run(run_id: str, request: Request):
    with _store(request) as store:
        try:
            return store.get_run(run_id)
        except KeyError as exc:
            raise HTTPException(404, str(exc)) from exc


@router.post("/runs/{run_id}/start")
def start_run(run_id: str, request: Request):
    principal = _require(request, "coordinator", "operator")
    with _store(request, agent_id=principal.agent_id) as store:
        try:
            result = store.start_run(run_id)
        except KeyError as exc:
            raise HTTPException(404, str(exc)) from exc
        except ValueError as exc:
            raise _conflict(exc) from exc
    _notify(request, "run", run_id)
    return result


@router.post("/runs/{run_id}/submit")
def submit_run(run_id: str, request: Request):
    principal = _require(request, "coordinator", "operator")
    with _store(request, agent_id=principal.agent_id) as store:
        try:
            result = store.submit_run(run_id)
        except KeyError as exc:
            raise HTTPException(404, str(exc)) from exc
        except ValueError as exc:
            raise _conflict(exc) from exc
    _notify(request, "run", run_id)
    return result


@router.post("/runs/{run_id}/heartbeat")
def heartbeat_run(run_id: str, request: Request):
    principal = _require(request, "worker", "coordinator", "operator")
    with _store(request, agent_id=principal.agent_id) as store:
        try:
            result = store.heartbeat_run(run_id)
        except KeyError as exc:
            raise HTTPException(404, str(exc)) from exc
        except ValueError as exc:
            raise _conflict(exc) from exc
    _notify(request, "run", run_id)
    return result


@router.patch("/runs/{run_id}")
def update_run(run_id: str, update: RunUpdate, request: Request):
    principal = _require(request, "coordinator", "operator")
    with _store(request, agent_id=principal.agent_id) as store:
        try:
            if update.status == "complete":
                result = store.complete_run(run_id, result=update.result)
            else:
                if not update.error:
                    raise HTTPException(400, "error is required when failing a run")
                result = store.fail_run(run_id, update.error)
        except KeyError as exc:
            raise HTTPException(404, str(exc)) from exc
        except ValueError as exc:
            raise _conflict(exc) from exc
    _notify(request, "run", run_id)
    return result


@router.get("/runs/{run_id}/events")
def get_run_events(
    run_id: str,
    request: Request,
    after_id: int = Query(default=0, ge=0),
    limit: int = Query(default=1_000, ge=1, le=1_000),
):
    with _store(request) as store:
        try:
            return store.run_events(run_id, after_id=after_id, limit=limit)
        except KeyError as exc:
            raise HTTPException(404, str(exc)) from exc


@router.get("/runs/{run_id}/stages")
def get_run_stages(
    run_id: str,
    request: Request,
    stage_id: str | None = None,
):
    with _store(request) as store:
        try:
            store.get_run(run_id)
            return store.list_run_stages(run_id, stage_id=stage_id)
        except KeyError as exc:
            raise HTTPException(404, str(exc)) from exc


@router.post("/coordinator/log", status_code=201)
def create_log_entry(entry: CoordinatorLogIn, request: Request):
    principal = _require(request, "coordinator", "operator")
    with _store(request, agent_id=principal.agent_id) as store:
        try:
            entry_id = store.log(
                entry.action_type,
                entry.reasoning,
                referenced_finding_ids=entry.referenced_finding_ids,
                referenced_hypothesis_ids=entry.referenced_hypothesis_ids,
                decisions_made=entry.decisions_made,
                campaign_id=entry.campaign_id,
                idempotency_key=entry.idempotency_key,
            )
            created = store.get_log_entry(entry_id)
        except (ValueError, IdempotencyConflictError) as exc:
            raise _conflict(exc) from exc
    _notify(request, "coordinator_log", entry_id)
    return {"id": entry_id, "ts": created["ts"], "entry": created}


@router.get("/coordinator/log")
def get_log(
    request: Request,
    limit: int = Query(default=20, ge=1, le=100),
    since: float = 0,
    campaign_id: str | None = None,
    before_ts: float | None = None,
):
    with _store(request) as store:
        return store.recent_log(
            limit=limit,
            since=since,
            campaign_id=campaign_id,
            before_ts=before_ts,
        )


@router.get("/events")
def get_events(
    request: Request,
    after_id: int = Query(default=0, ge=0),
    campaign_id: str | None = None,
    entity_type: str | None = None,
    limit: int = Query(default=1_000, ge=1, le=1_000),
):
    with _store(request) as store:
        return store.events(
            after_id=after_id,
            campaign_id=campaign_id,
            entity_type=entity_type,
            limit=limit,
        )


@router.get("/stream")
async def event_stream(
    request: Request,
    after_id: int | None = Query(default=None, ge=0),
    campaign_id: str | None = None,
):
    """Replay durable events and then tail them using the bus as a wake-up."""

    header_cursor = request.headers.get("last-event-id")
    if after_id is not None:
        cursor = after_id
    elif header_cursor is not None:
        try:
            cursor = max(0, int(header_cursor))
        except ValueError as exc:
            raise HTTPException(400, "Last-Event-ID must be an integer") from exc
    else:
        cursor = 0
    bus = _runtime(request).event_bus
    subscriber = bus.subscribe()

    def read_events(current: int) -> list[dict[str, Any]]:
        with _store(request) as store:
            return store.events(
                after_id=current,
                campaign_id=campaign_id,
                limit=1_000,
            )

    async def generate():
        nonlocal cursor
        try:
            while True:
                events = await asyncio.to_thread(read_events, cursor)
                if events:
                    for event in events:
                        cursor = int(event["id"])
                        event_type = str(event["event_type"]).replace("\n", "_")
                        data = json.dumps(event, separators=(",", ":"), default=str)
                        yield f"id: {cursor}\nevent: {event_type}\ndata: {data}\n\n"
                    continue
                try:
                    await asyncio.wait_for(subscriber.queue.get(), timeout=15)
                except asyncio.TimeoutError:
                    yield ": keepalive\n\n"
        except asyncio.CancelledError:
            return
        finally:
            bus.unsubscribe(subscriber)

    return StreamingResponse(
        generate(),
        media_type="text/event-stream",
        headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"},
    )


def _backup_destination(runtime: OpsRuntime, directory: Path) -> Path:
    now_ns = time.time_ns()
    timestamp = time.strftime(
        "%Y%m%dT%H%M%SZ",
        time.gmtime(now_ns / 1_000_000_000),
    )
    return directory / (
        f"{runtime.db_path.stem}-{timestamp}-{now_ns % 1_000_000_000:09d}.db"
    )


@router.post("/admin/backup")
def backup(request: Request):
    principal = _require(request, "operator")
    runtime = _runtime(request)
    backup_dir = runtime.backup_dir or runtime.db_path.parent / "ops_backups"
    try:
        with _store(request, agent_id=principal.agent_id) as store:
            result = store.backup(_backup_destination(runtime, backup_dir))
    except Exception as exc:
        with runtime.state_lock:
            runtime.backup_failure_count += 1
            runtime.last_backup_error = str(exc)
        raise
    with runtime.state_lock:
        runtime.last_backup_ts = float(result["created_ts"])
        runtime.last_backup_error = None
    return result


@router.get("/admin/integrity")
def integrity(request: Request):
    _require(request, "operator")
    with _store(request) as store:
        result = store.integrity_check()
    if not result["ok"]:
        raise HTTPException(503, detail=result)
    return result


def _runtime_observability(runtime: OpsRuntime) -> dict[str, Any]:
    now = time.time()
    with runtime.state_lock:
        durations = sorted(runtime.request_durations)
        if durations:
            p95_index = max(0, (95 * len(durations) + 99) // 100 - 1)
            p95_duration = durations[p95_index]
        else:
            p95_duration = 0.0
        last_backup_ts = runtime.last_backup_ts or None
        return {
            "http": {
                "requests": runtime.request_count,
                "errors": runtime.request_error_count,
                "duration_seconds": runtime.request_duration_seconds,
                "duration_buckets": dict(runtime.request_duration_buckets),
                "p95_duration_seconds": p95_duration,
                "latency_window_requests": len(durations),
                "sqlite_busy_errors": runtime.sqlite_busy_errors,
            },
            "backup": {
                "enabled": (
                    runtime.backup_dir is not None
                    and runtime.backup_interval_seconds > 0
                ),
                "interval_seconds": runtime.backup_interval_seconds,
                "last_success_ts": last_backup_ts,
                "age_seconds": (
                    max(0.0, now - last_backup_ts)
                    if last_backup_ts is not None
                    else None
                ),
                "failure_count": runtime.backup_failure_count,
                "last_error": runtime.last_backup_error,
            },
        }


@router.get("/stats")
def stats(request: Request):
    runtime = _runtime(request)
    with _store(request) as store:
        result = store.stats()
    result["pool"] = runtime.pool.stats() if runtime.pool is not None else None
    result.update(_runtime_observability(runtime))
    result["event_bus"] = runtime.event_bus.stats()
    return result


def _prometheus(stats_payload: dict[str, Any]) -> str:
    lines = [
        "# HELP sharur_ops_records Durable coordination records.",
        "# TYPE sharur_ops_records gauge",
    ]
    for table, count in sorted(stats_payload["counts"].items()):
        lines.append(f'sharur_ops_records{{table="{table}"}} {count}')
    lines.extend(
        [
            "# HELP sharur_ops_tasks Task records by lifecycle status.",
            "# TYPE sharur_ops_tasks gauge",
        ]
    )
    for status, count in sorted(stats_payload["tasks_by_status"].items()):
        lines.append(f'sharur_ops_tasks{{status="{status}"}} {count}')
    queue = stats_payload["queue"]
    lines.extend(
        [
            "# TYPE sharur_ops_queue_oldest_age_seconds gauge",
            f"sharur_ops_queue_oldest_age_seconds {queue['oldest_age_seconds']}",
            "# TYPE sharur_ops_active_leases gauge",
            f"sharur_ops_active_leases {queue['active_leases']}",
            "# TYPE sharur_ops_expired_leases gauge",
            f"sharur_ops_expired_leases {queue['expired_leases']}",
        ]
    )
    dead_letters = stats_payload["dead_letters"]
    lines.extend(
        [
            "# TYPE sharur_ops_dead_letter_tasks gauge",
            f"sharur_ops_dead_letter_tasks {dead_letters['count']}",
            "# TYPE sharur_ops_attempts_exhausted_tasks gauge",
            "sharur_ops_attempts_exhausted_tasks "
            f"{dead_letters['attempts_exhausted']}",
            "# TYPE sharur_ops_dead_letter_oldest_age_seconds gauge",
            "sharur_ops_dead_letter_oldest_age_seconds "
            f"{dead_letters['oldest_age_seconds']}",
        ]
    )
    pool = stats_payload.get("pool") or {}
    for key in (
        "checked_out",
        "high_watermark",
        "wait_count",
        "wait_seconds_total",
        "wait_seconds_max",
        "sqlite_write_transactions",
        "sqlite_write_wait_seconds",
        "sqlite_write_wait_max_seconds",
    ):
        if key in pool:
            lines.append(f"sharur_ops_pool_{key} {pool[key]}")
    event_bus = stats_payload.get("event_bus") or {}
    if event_bus:
        lines.extend(
            [
                "# TYPE sharur_ops_sse_subscribers gauge",
                f"sharur_ops_sse_subscribers {event_bus['subscribers']}",
                "# TYPE sharur_ops_sse_dropped_wakeups counter",
                "sharur_ops_sse_dropped_wakeups "
                f"{event_bus['dropped_wakeups']}",
            ]
        )
    backup = stats_payload.get("backup") or {}
    if backup:
        lines.extend(
            [
                "# TYPE sharur_ops_backup_enabled gauge",
                f"sharur_ops_backup_enabled {int(backup['enabled'])}",
                "# TYPE sharur_ops_backup_failures_total counter",
                f"sharur_ops_backup_failures_total {backup['failure_count']}",
            ]
        )
        if backup["age_seconds"] is not None:
            lines.extend(
                [
                    "# TYPE sharur_ops_backup_age_seconds gauge",
                    f"sharur_ops_backup_age_seconds {backup['age_seconds']}",
                ]
            )
    http = stats_payload.get("http") or {}
    if http:
        lines.extend(
            [
                "# TYPE sharur_ops_http_requests_total counter",
                f"sharur_ops_http_requests_total {http['requests']}",
                "# TYPE sharur_ops_http_errors_total counter",
                f"sharur_ops_http_errors_total {http['errors']}",
                "# TYPE sharur_ops_sqlite_busy_errors_total counter",
                "sharur_ops_sqlite_busy_errors_total "
                f"{http['sqlite_busy_errors']}",
                "# TYPE sharur_ops_http_p95_duration_seconds gauge",
                "sharur_ops_http_p95_duration_seconds "
                f"{http['p95_duration_seconds']}",
                "# TYPE sharur_ops_http_duration_seconds histogram",
            ]
        )
        for upper_bound, count in sorted(http["duration_buckets"].items()):
            lines.append(
                "sharur_ops_http_duration_seconds_bucket"
                f'{{le="{upper_bound}"}} {count}'
            )
        lines.append(
            "sharur_ops_http_duration_seconds_bucket"
            f'{{le="+Inf"}} {http["requests"]}'
        )
        lines.append(
            f"sharur_ops_http_duration_seconds_sum {http['duration_seconds']}"
        )
        lines.append(
            f"sharur_ops_http_duration_seconds_count {http['requests']}"
        )
    return "\n".join(lines) + "\n"


@router.get("/metrics")
def metrics(request: Request):
    runtime = _runtime(request)
    with _store(request) as store:
        result = store.stats()
    result["pool"] = runtime.pool.stats() if runtime.pool is not None else None
    result.update(_runtime_observability(runtime))
    result["event_bus"] = runtime.event_bus.stats()
    return PlainTextResponse(_prometheus(result), media_type="text/plain; version=0.0.4")


@router.get("/auth/whoami")
def whoami(request: Request):
    """Expose the authenticated principal for peer-service token introspection."""

    principal = _principal(request)
    return {
        "agent_id": principal.agent_id,
        "role": principal.role,
        "bootstrap": principal.bootstrap,
    }


@router.get("/health")
def health(request: Request):
    runtime = _runtime(request)
    with _store(request) as store:
        database = store.stats()["database"]
    return {
        "status": "ok",
        "db": str(runtime.db_path),
        "schema_version": database["schema_version"],
        "auth_required": runtime.api_token is not None,
        "sqlite_owner": "single-http-server",
        "pool": runtime.pool.stats() if runtime.pool is not None else None,
        "ts": time.time(),
    }


def _maintenance_once(runtime: OpsRuntime) -> None:
    if runtime.pool is None:
        return
    connection = runtime.pool.acquire()
    store: OpsStore | None = None
    try:
        store = OpsStore(
            runtime.db_path,
            agent_id="system-maintenance",
            connection=connection,
            initialize=False,
            close_callback=runtime.pool.release,
            transaction_wait_observer=runtime.pool.observe_sqlite_write_wait,
        )
        with store:
            store.recover_expired_tasks()
            store.recover_stale_runs(
                stale_after_seconds=runtime.run_stale_after_seconds,
            )
            if runtime.backup_dir is not None and runtime.backup_interval_seconds > 0:
                now = time.time()
                with runtime.state_lock:
                    due = now - runtime.last_backup_ts >= runtime.backup_interval_seconds
                if due:
                    try:
                        result = store.backup(
                            _backup_destination(runtime, runtime.backup_dir)
                        )
                    except Exception as exc:
                        with runtime.state_lock:
                            runtime.backup_failure_count += 1
                            runtime.last_backup_error = str(exc)
                        raise
                    with runtime.state_lock:
                        runtime.last_backup_ts = float(result["created_ts"])
                        runtime.last_backup_error = None
    except Exception:
        if store is None and runtime.pool is not None:
            runtime.pool.release(connection)
        raise


async def _maintenance_loop(runtime: OpsRuntime) -> None:
    logger = logging.getLogger("sharur.ops")
    while True:
        await asyncio.sleep(runtime.maintenance_interval_seconds)
        try:
            await asyncio.to_thread(_maintenance_once, runtime)
            runtime.event_bus.publish("maintenance")
        except asyncio.CancelledError:
            raise
        except Exception:
            logger.exception("Sharur Ops maintenance cycle failed")


def create_app(
    *,
    db_path: str | Path | None = None,
    api_token: str | None = None,
    allow_insecure_remote: bool | None = None,
    sse_queue_size: int = DEFAULT_SSE_QUEUE_SIZE,
    pool_size: int = DEFAULT_POOL_SIZE,
    max_request_bytes: int = DEFAULT_MAX_REQUEST_BYTES,
    maintenance_interval_seconds: int = DEFAULT_MAINTENANCE_INTERVAL_SECONDS,
    run_stale_after_seconds: int = DEFAULT_RUN_STALE_AFTER_SECONDS,
    backup_dir: str | Path | None = None,
    backup_interval_seconds: int = 0,
) -> FastAPI:
    """Build an isolated app; SQLite opens only during application startup."""

    if pool_size < 1:
        raise ValueError("pool_size must be positive")
    if max_request_bytes < 1:
        raise ValueError("max_request_bytes must be positive")
    if maintenance_interval_seconds < 1:
        raise ValueError("maintenance_interval_seconds must be positive")
    configured_path = (
        Path(db_path)
        if db_path is not None
        else Path(os.environ.get("SHARUR_OPS_DB_PATH", "sharur_ops.db"))
    )
    configured_token = (
        api_token if api_token is not None else os.environ.get("SHARUR_OPS_TOKEN")
    )
    normalized_token = configured_token.strip() if configured_token else None
    insecure_remote = (
        _env_flag("SHARUR_OPS_ALLOW_INSECURE_REMOTE")
        if allow_insecure_remote is None
        else allow_insecure_remote
    )
    configured_backup_dir = (
        Path(backup_dir).expanduser().resolve()
        if backup_dir is not None
        else None
    )
    runtime = OpsRuntime(
        db_path=configured_path.expanduser().resolve(),
        api_token=normalized_token,
        allow_insecure_remote=insecure_remote,
        event_bus=EventBus(queue_size=sse_queue_size),
        pool_size=pool_size,
        max_request_bytes=max_request_bytes,
        maintenance_interval_seconds=maintenance_interval_seconds,
        run_stale_after_seconds=run_stale_after_seconds,
        backup_dir=configured_backup_dir,
        backup_interval_seconds=max(0, int(backup_interval_seconds)),
    )

    @asynccontextmanager
    async def lifespan(application: FastAPI):
        application.state.ops_runtime = runtime
        runtime.owner_lock = SQLiteServerLock(runtime.db_path)
        runtime.owner_lock.acquire()
        try:
            if runtime.backup_dir is not None:
                runtime.backup_dir.mkdir(parents=True, exist_ok=True)
                prior_backups = runtime.backup_dir.glob(
                    f"{runtime.db_path.stem}-*.db"
                )
                runtime.last_backup_ts = max(
                    (path.stat().st_mtime for path in prior_backups),
                    default=0.0,
                )
            runtime.pool = SQLiteConnectionPool(
                runtime.db_path,
                size=runtime.pool_size,
            )
            runtime.maintenance_task = asyncio.create_task(
                _maintenance_loop(runtime),
                name="sharur-ops-maintenance",
            )
            yield
        finally:
            if runtime.maintenance_task is not None:
                runtime.maintenance_task.cancel()
                with contextlib.suppress(asyncio.CancelledError):
                    await runtime.maintenance_task
                runtime.maintenance_task = None
            if runtime.pool is not None:
                runtime.pool.close()
                runtime.pool = None
            if runtime.owner_lock is not None:
                runtime.owner_lock.release()
                runtime.owner_lock = None

    application = FastAPI(
        title="Sharur Ops",
        version=__version__,
        lifespan=lifespan,
    )
    application.add_middleware(
        RequestSizeLimitMiddleware,
        max_bytes=max_request_bytes,
    )
    application.add_middleware(
        RequestMetricsMiddleware,
        runtime=runtime,
    )
    application.state.ops_runtime = runtime
    application.include_router(router)
    return application


app = create_app()


def main() -> None:
    """Run one coordination owner with guarded network defaults."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--host",
        default=os.environ.get("SHARUR_OPS_HOST", DEFAULT_HOST),
        help="Bind host (default: 127.0.0.1)",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=int(os.environ.get("SHARUR_OPS_PORT", DEFAULT_PORT)),
    )
    parser.add_argument(
        "--db",
        type=Path,
        default=Path(os.environ.get("SHARUR_OPS_DB_PATH", "sharur_ops.db")),
    )
    parser.add_argument(
        "--pool-size",
        type=int,
        default=int(os.environ.get("SHARUR_OPS_POOL_SIZE", DEFAULT_POOL_SIZE)),
    )
    parser.add_argument(
        "--backup-dir",
        type=Path,
        default=(
            Path(os.environ["SHARUR_OPS_BACKUP_DIR"])
            if os.environ.get("SHARUR_OPS_BACKUP_DIR")
            else None
        ),
    )
    parser.add_argument(
        "--backup-interval",
        type=int,
        default=int(os.environ.get("SHARUR_OPS_BACKUP_INTERVAL_SECONDS", "0")),
        help="Seconds between online SQLite backups; 0 disables scheduled backups",
    )
    parser.add_argument(
        "--allow-insecure-remote",
        action="store_true",
        default=_env_flag("SHARUR_OPS_ALLOW_INSECURE_REMOTE"),
        help="Allow remote clients without SHARUR_OPS_TOKEN (unsafe)",
    )
    args = parser.parse_args()
    token = os.environ.get("SHARUR_OPS_TOKEN")
    if not _is_loopback(args.host) and not token and not args.allow_insecure_remote:
        parser.error(
            "remote binding requires SHARUR_OPS_TOKEN; "
            "use --allow-insecure-remote only on a trusted network"
        )

    import uvicorn  # noqa: PLC0415

    server_app = create_app(
        db_path=args.db,
        api_token=token,
        allow_insecure_remote=args.allow_insecure_remote,
        pool_size=args.pool_size,
        backup_dir=args.backup_dir,
        backup_interval_seconds=args.backup_interval,
    )
    uvicorn.run(server_app, host=args.host, port=args.port, workers=1)


__all__ = [
    "EventBus",
    "OpsRuntime",
    "Principal",
    "app",
    "create_app",
    "main",
]


if __name__ == "__main__":
    main()
