"""Optional HTTP coordination service for distributed Sharur agents.

The service is a transport around :class:`sharur.ops.store.OpsStore`. It keeps
coordination state in ``sharur_ops.db``; canonical scientific findings remain
in the dataset's JSONL records and DuckDB.

Use ``sharur-ops`` (or ``python -m sharur.ops.server``) for the guarded server
launcher. It binds to loopback by default and refuses unauthenticated remote
clients unless an explicit insecure override is enabled.
"""

from __future__ import annotations

import argparse
import asyncio
import contextlib
import json
import os
import secrets
import threading
import time
from contextlib import asynccontextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated, Literal

from fastapi import APIRouter, Depends, FastAPI, HTTPException, Query, Request
from fastapi.responses import StreamingResponse
from fastapi.security import HTTPAuthorizationCredentials, HTTPBearer
from pydantic import BaseModel, Field

from sharur import __version__
from sharur.ops.schema import DEFAULT_LEASE_SECONDS
from sharur.ops.store import OpsStore


DEFAULT_HOST = "127.0.0.1"
DEFAULT_PORT = 8811
DEFAULT_SSE_QUEUE_SIZE = 1_000


def _env_flag(name: str, *, default: bool = False) -> bool:
    value = os.environ.get(name)
    if value is None:
        return default
    return value.strip().lower() in {"1", "true", "yes", "on"}


def _is_loopback(host: str | None) -> bool:
    if host is None:
        return False
    normalized = host.strip().lower()
    return normalized in {"127.0.0.1", "::1", "localhost", "testclient"}


@dataclass(frozen=True)
class _Subscriber:
    queue: asyncio.Queue[str]
    loop: asyncio.AbstractEventLoop


class EventBus:
    """Bounded, thread-safe in-process pub/sub for the SSE endpoint."""

    def __init__(self, *, queue_size: int = DEFAULT_SSE_QUEUE_SIZE):
        if queue_size < 1:
            raise ValueError("SSE queue size must be positive")
        self._queue_size = queue_size
        self._subscribers: list[_Subscriber] = []
        self._lock = threading.Lock()

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
                current for current in self._subscribers if current is not subscriber
            ]

    @staticmethod
    def _offer(queue: asyncio.Queue[str], payload: str) -> None:
        if queue.full():
            with contextlib.suppress(asyncio.QueueEmpty):
                queue.get_nowait()
        with contextlib.suppress(asyncio.QueueFull):
            queue.put_nowait(payload)

    def publish(self, event_type: str, data: dict) -> None:
        payload = json.dumps({"type": event_type, **data}, separators=(",", ":"))
        with self._lock:
            subscribers = list(self._subscribers)
        for subscriber in subscribers:
            with contextlib.suppress(RuntimeError):
                subscriber.loop.call_soon_threadsafe(
                    self._offer,
                    subscriber.queue,
                    payload,
                )


@dataclass(frozen=True)
class OpsRuntime:
    """Application-scoped configuration with no open database connection."""

    db_path: Path
    api_token: str | None
    allow_insecure_remote: bool
    event_bus: EventBus


class FindingIn(BaseModel):
    agent_id: str = Field(min_length=1)
    finding_type: str = Field(min_length=1)
    domain: str = Field(min_length=1)
    summary: str
    evidence: dict = Field(default_factory=dict)
    confidence: float = Field(default=0.5, ge=0.0, le=1.0)
    novelty: int = Field(default=0, ge=0, le=3)
    parent_finding_id: str | None = None
    reasoning: str = ""


class HypothesisIn(BaseModel):
    source_agent_id: str = Field(min_length=1)
    source_finding_ids: list[str] = Field(default_factory=list)
    hypothesis: str = Field(min_length=1)
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
    task_type: str = Field(min_length=1)
    description: str = Field(min_length=1)
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
    action_type: str = Field(min_length=1)
    reasoning: str
    referenced_finding_ids: list[str] = Field(default_factory=list)
    referenced_hypothesis_ids: list[str] = Field(default_factory=list)
    decisions_made: dict = Field(default_factory=dict)


def _runtime(request: Request) -> OpsRuntime:
    return request.app.state.ops_runtime


def _store(request: Request, *, agent_id: str = "server") -> OpsStore:
    return OpsStore(_runtime(request).db_path, agent_id=agent_id)


_bearer = HTTPBearer(auto_error=False)


def _authorize(
    request: Request,
    credentials: Annotated[HTTPAuthorizationCredentials | None, Depends(_bearer)],
) -> None:
    runtime = _runtime(request)
    if runtime.api_token is not None:
        valid = (
            credentials is not None
            and credentials.scheme.lower() == "bearer"
            and secrets.compare_digest(credentials.credentials, runtime.api_token)
        )
        if not valid:
            raise HTTPException(
                status_code=401,
                detail="Missing or invalid Sharur Ops bearer token",
                headers={"WWW-Authenticate": "Bearer"},
            )
        return

    client_host = request.client.host if request.client is not None else None
    if not runtime.allow_insecure_remote and not _is_loopback(client_host):
        raise HTTPException(
            status_code=403,
            detail=(
                "Unauthenticated Sharur Ops access is restricted to loopback; "
                "set SHARUR_OPS_TOKEN for remote access"
            ),
        )


router = APIRouter(dependencies=[Depends(_authorize)])


@router.post("/findings", status_code=201)
def create_finding(finding: FindingIn, request: Request):
    with _store(request, agent_id=finding.agent_id) as store:
        finding_id = store.finding(
            finding.finding_type,
            finding.domain,
            finding.summary,
            evidence=finding.evidence,
            confidence=finding.confidence,
            novelty=finding.novelty,
            parent_finding_id=finding.parent_finding_id,
            reasoning=finding.reasoning,
        )
        created = store.get_finding(finding_id)
    _runtime(request).event_bus.publish(
        "finding",
        {
            "id": finding_id,
            "agent_id": finding.agent_id,
            "novelty": finding.novelty,
            "summary": finding.summary,
        },
    )
    return {"id": finding_id, "ts": created["ts"]}


@router.get("/findings")
def list_findings(
    request: Request,
    since: float = 0,
    min_novelty: int = 0,
    finding_type: str | None = None,
    domain: str | None = None,
    agent_id: str | None = None,
    limit: int = Query(default=50, ge=1, le=500),
):
    with _store(request) as store:
        return store.recent_findings(
            since=since,
            min_novelty=min_novelty,
            finding_type=finding_type,
            domain=domain,
            agent_id=agent_id,
            limit=limit,
        )


@router.get("/findings/search/{text}")
def search_findings(
    text: str,
    request: Request,
    limit: int = Query(default=20, ge=1, le=100),
):
    with _store(request) as store:
        return store.search_findings(text, limit=limit)


@router.get("/findings/{finding_id}")
def get_finding(finding_id: str, request: Request):
    with _store(request) as store:
        try:
            return store.get_finding(finding_id)
        except KeyError as exc:
            raise HTTPException(404, str(exc)) from exc


@router.post("/hypotheses", status_code=201)
def create_hypothesis(hypothesis: HypothesisIn, request: Request):
    with _store(request, agent_id=hypothesis.source_agent_id) as store:
        hypothesis_id = store.hypothesis(
            hypothesis.hypothesis,
            source_finding_ids=hypothesis.source_finding_ids,
            domains_relevant=hypothesis.domains_relevant,
        )
        created = store.get_hypothesis(hypothesis_id)
    _runtime(request).event_bus.publish(
        "hypothesis",
        {"id": hypothesis_id, "hypothesis": hypothesis.hypothesis},
    )
    return {"id": hypothesis_id, "ts": created["ts"]}


@router.get("/hypotheses")
def list_hypotheses(
    request: Request,
    status: str | None = None,
    unassigned: bool = False,
    limit: int = Query(default=50, ge=1, le=200),
):
    with _store(request) as store:
        return store.list_hypotheses(
            status=status,
            unassigned=unassigned,
            limit=limit,
        )


@router.patch("/hypotheses/{hypothesis_id}")
def update_hypothesis(
    hypothesis_id: str,
    update: HypothesisUpdate,
    request: Request,
):
    values = update.model_dump(exclude_none=True)
    if not values:
        raise HTTPException(400, "No hypothesis updates supplied")
    with _store(request) as store:
        try:
            return store.update_hypothesis(hypothesis_id, **values)
        except KeyError as exc:
            raise HTTPException(404, str(exc)) from exc
        except ValueError as exc:
            raise HTTPException(400, str(exc)) from exc


@router.post("/tasks", status_code=201)
def create_task(task: TaskIn, request: Request):
    with _store(request, agent_id=task.created_by) as store:
        try:
            task_id = store.create_task(
                task.task_type,
                task.description,
                params=task.params,
                priority=task.priority,
                domain_hint=task.domain_hint,
                assigned_to=task.assigned_to,
                run_id=task.run_id,
                idempotency_key=task.idempotency_key,
                depends_on=task.depends_on,
                max_attempts=task.max_attempts,
                lease_seconds=task.lease_seconds,
            )
            created = store.get_task(task_id)
        except ValueError as exc:
            raise HTTPException(409, str(exc)) from exc
    _runtime(request).event_bus.publish(
        "task",
        {"id": task_id, "description": task.description, "priority": task.priority},
    )
    return {"id": task_id, "ts": created["ts"], "task": created}


@router.get("/tasks")
def list_tasks(
    request: Request,
    status: str | None = None,
    assigned_to: str | None = None,
    unassigned: bool = False,
    limit: int = Query(default=50, ge=1, le=200),
):
    with _store(request) as store:
        return store.list_tasks(
            status=status,
            assigned_to=assigned_to,
            unassigned=unassigned,
            limit=limit,
        )


@router.post("/tasks/recover")
def recover_tasks(request: Request, now: float | None = None):
    with _store(request, agent_id="coordinator") as store:
        return store.recover_expired_tasks(now=now)


@router.patch("/tasks/{task_id}")
def update_task(task_id: str, update: TaskUpdate, request: Request):
    if not update.agent_id:
        raise HTTPException(400, "agent_id is required for task updates")
    with _store(request, agent_id=update.agent_id) as store:
        try:
            if update.status == "complete":
                return store.complete_task(task_id, update.result_finding_ids)
            if update.status == "failed":
                return store.fail_task(
                    task_id,
                    error=update.error,
                    retryable=update.retryable,
                    retry_delay_seconds=update.retry_delay_seconds,
                )
            if update.status == "in_progress":
                return store.heartbeat_task(
                    task_id,
                    lease_seconds=update.lease_seconds,
                    in_progress=True,
                )
            raise HTTPException(
                400,
                "Use the claim endpoint for ownership; only in_progress, complete, "
                "and failed are valid worker updates.",
            )
        except (KeyError, ValueError) as exc:
            raise HTTPException(409, str(exc)) from exc


@router.post("/tasks/{task_id}/claim")
def claim_task(
    task_id: str,
    request: Request,
    agent_id: str,
    lease_seconds: int | None = Query(default=None, ge=1),
):
    with _store(request, agent_id=agent_id) as store:
        try:
            return store.claim_task(task_id, lease_seconds=lease_seconds)
        except ValueError as exc:
            raise HTTPException(409, str(exc)) from exc


@router.post("/tasks/{task_id}/heartbeat")
def heartbeat_task(
    task_id: str,
    request: Request,
    agent_id: str,
    lease_seconds: int | None = Query(default=None, ge=1),
):
    with _store(request, agent_id=agent_id) as store:
        try:
            return store.heartbeat_task(task_id, lease_seconds=lease_seconds)
        except ValueError as exc:
            raise HTTPException(409, str(exc)) from exc


@router.post("/runs", status_code=201)
def create_run(run: RunIn, request: Request):
    with _store(request, agent_id=run.created_by) as store:
        try:
            run_id = store.create_run(
                run.run_type,
                run.dataset_path,
                config=run.config,
                idempotency_key=run.idempotency_key,
                parent_run_id=run.parent_run_id,
            )
            return store.get_run(run_id)
        except ValueError as exc:
            raise HTTPException(409, str(exc)) from exc


@router.get("/runs")
def list_runs(
    request: Request,
    dataset_path: str | None = None,
    run_type: str | None = None,
    status: str | None = None,
    limit: int = Query(default=50, ge=1, le=200),
):
    with _store(request) as store:
        return store.list_runs(
            dataset_path=dataset_path,
            run_type=run_type,
            status=status,
            limit=limit,
        )


@router.post("/runs/recover")
def recover_runs(
    request: Request,
    stale_after_seconds: int = Query(default=300, ge=1),
    now: float | None = None,
):
    with _store(request) as store:
        return store.recover_stale_runs(
            stale_after_seconds=stale_after_seconds,
            now=now,
        )


@router.get("/runs/{run_id}")
def get_run(run_id: str, request: Request):
    with _store(request) as store:
        try:
            return store.get_run(run_id)
        except KeyError as exc:
            raise HTTPException(404, str(exc)) from exc


@router.post("/runs/{run_id}/start")
def start_run(run_id: str, request: Request):
    with _store(request) as store:
        try:
            return store.start_run(run_id)
        except KeyError as exc:
            raise HTTPException(404, str(exc)) from exc
        except ValueError as exc:
            raise HTTPException(409, str(exc)) from exc


@router.post("/runs/{run_id}/submit")
def submit_run(run_id: str, request: Request):
    with _store(request) as store:
        try:
            return store.submit_run(run_id)
        except KeyError as exc:
            raise HTTPException(404, str(exc)) from exc
        except ValueError as exc:
            raise HTTPException(409, str(exc)) from exc


@router.post("/runs/{run_id}/heartbeat")
def heartbeat_run(run_id: str, request: Request):
    with _store(request) as store:
        try:
            return store.heartbeat_run(run_id)
        except KeyError as exc:
            raise HTTPException(404, str(exc)) from exc
        except ValueError as exc:
            raise HTTPException(409, str(exc)) from exc


@router.patch("/runs/{run_id}")
def update_run(run_id: str, update: RunUpdate, request: Request):
    with _store(request) as store:
        try:
            if update.status == "complete":
                return store.complete_run(run_id, result=update.result)
            if not update.error:
                raise HTTPException(400, "error is required when failing a run")
            return store.fail_run(run_id, update.error)
        except KeyError as exc:
            raise HTTPException(404, str(exc)) from exc
        except ValueError as exc:
            raise HTTPException(409, str(exc)) from exc


@router.get("/runs/{run_id}/events")
def get_run_events(run_id: str, request: Request):
    with _store(request) as store:
        try:
            return store.run_events(run_id)
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
    with _store(request, agent_id="coordinator") as store:
        entry_id = store.log(
            entry.action_type,
            entry.reasoning,
            referenced_finding_ids=entry.referenced_finding_ids,
            referenced_hypothesis_ids=entry.referenced_hypothesis_ids,
            decisions_made=entry.decisions_made,
        )
        created = store.get_log_entry(entry_id)
    return {"id": entry_id, "ts": created["ts"]}


@router.get("/coordinator/log")
def get_log(
    request: Request,
    limit: int = Query(default=20, ge=1, le=100),
    since: float = 0,
):
    with _store(request) as store:
        return store.recent_log(limit=limit, since=since)


@router.get("/stream")
async def event_stream(request: Request):
    """Stream bounded in-process notifications for new coordination records."""
    bus = _runtime(request).event_bus
    subscriber = bus.subscribe()

    async def generate():
        try:
            while True:
                try:
                    data = await asyncio.wait_for(subscriber.queue.get(), timeout=30)
                    yield f"data: {data}\n\n"
                except asyncio.TimeoutError:
                    yield ": keepalive\n\n"
        except asyncio.CancelledError:
            return
        finally:
            bus.unsubscribe(subscriber)

    return StreamingResponse(generate(), media_type="text/event-stream")


@router.get("/stats")
def stats(request: Request):
    with _store(request) as store:
        return store.stats()


@router.get("/health")
def health(request: Request):
    runtime = _runtime(request)
    with _store(request):
        pass
    return {
        "status": "ok",
        "db": str(runtime.db_path),
        "auth_required": runtime.api_token is not None,
        "ts": time.time(),
    }


def create_app(
    *,
    db_path: str | Path | None = None,
    api_token: str | None = None,
    allow_insecure_remote: bool | None = None,
    sse_queue_size: int = DEFAULT_SSE_QUEUE_SIZE,
) -> FastAPI:
    """Build an isolated app without opening SQLite until application startup."""
    configured_path = (
        Path(db_path)
        if db_path is not None
        else Path(os.environ.get("SHARUR_OPS_DB_PATH", "sharur_ops.db"))
    )
    configured_token = api_token if api_token is not None else os.environ.get("SHARUR_OPS_TOKEN")
    normalized_token = configured_token.strip() if configured_token else None
    insecure_remote = (
        _env_flag("SHARUR_OPS_ALLOW_INSECURE_REMOTE")
        if allow_insecure_remote is None
        else allow_insecure_remote
    )
    runtime = OpsRuntime(
        db_path=configured_path.expanduser().resolve(),
        api_token=normalized_token,
        allow_insecure_remote=insecure_remote,
        event_bus=EventBus(queue_size=sse_queue_size),
    )

    @asynccontextmanager
    async def lifespan(application: FastAPI):
        application.state.ops_runtime = runtime
        with OpsStore(runtime.db_path, agent_id="server-bootstrap"):
            pass
        yield

    application = FastAPI(
        title="Sharur Ops",
        version=__version__,
        lifespan=lifespan,
    )
    application.state.ops_runtime = runtime
    application.include_router(router)
    return application


app = create_app()


def main() -> None:
    """Run the single-worker coordination server with safe network defaults."""
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

    import uvicorn  # noqa: PLC0415 - optional HTTP dependency

    server_app = create_app(
        db_path=args.db,
        api_token=token,
        allow_insecure_remote=args.allow_insecure_remote,
    )
    uvicorn.run(server_app, host=args.host, port=args.port, workers=1)


__all__ = ["EventBus", "OpsRuntime", "app", "create_app", "main"]


if __name__ == "__main__":
    main()
