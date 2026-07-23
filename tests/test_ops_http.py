"""Contract tests for the optional Sharur Ops HTTP transport."""

from __future__ import annotations

import asyncio
import subprocess
import sys
from typing import TYPE_CHECKING

import pytest
import requests
from fastapi.testclient import TestClient

from sharur.ops.client import DEFAULT_TIMEOUT, SharurOps
from sharur.ops.server import EventBus, create_app
from sharur.ops.store import OpsStore


if TYPE_CHECKING:
    from pathlib import Path

TOKEN = "test-secret"


def _headers() -> dict[str, str]:
    return {"Authorization": f"Bearer {TOKEN}"}


def test_import_has_no_default_database_side_effect(tmp_path: Path) -> None:
    process = subprocess.run(
        [sys.executable, "-c", "import sharur.ops.server"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )
    assert process.returncode == 0, process.stderr
    assert not (tmp_path / "sharur_ops.db").exists()


def test_app_lifecycle_is_lazy_isolated_and_authenticated(tmp_path: Path) -> None:
    first_db = tmp_path / "first" / "ops.db"
    second_db = tmp_path / "second" / "ops.db"
    first_app = create_app(db_path=first_db, api_token=TOKEN)
    second_app = create_app(db_path=second_db, api_token=TOKEN)

    assert not first_db.exists()
    assert not second_db.exists()

    with TestClient(first_app) as first:
        assert first_db.is_file()
        assert first.get("/health").status_code == 401
        assert first.get("/health", headers={"Authorization": "Bearer wrong"}).status_code == 401
        health = first.get("/health", headers=_headers())
        assert health.status_code == 200
        assert health.json()["db"] == str(first_db.resolve())
        assert health.json()["auth_required"] is True
        principal = first.get("/auth/whoami", headers=_headers())
        assert principal.json() == {
            "agent_id": "bootstrap-operator",
            "role": "operator",
            "bootstrap": True,
        }

        created = first.post(
            "/findings",
            headers=_headers(),
            json={
                "agent_id": "agent-a",
                "finding_type": "observation",
                "domain": "fixture",
                "summary": "bounded transport works",
                "novelty": 1,
            },
        )
        assert created.status_code == 201

    assert not second_db.exists()
    with TestClient(second_app) as second:
        assert second.get("/findings", headers=_headers()).json() == []
    assert second_db.is_file()


def test_unauthenticated_remote_clients_fail_closed(tmp_path: Path) -> None:
    remote = ("198.51.100.10", 50000)
    guarded = create_app(db_path=tmp_path / "guarded.db")
    with TestClient(guarded, client=remote) as client:
        response = client.get("/health")
        assert response.status_code == 403
        assert "restricted to loopback" in response.json()["detail"]

    explicitly_insecure = create_app(
        db_path=tmp_path / "insecure.db",
        allow_insecure_remote=True,
    )
    with TestClient(explicitly_insecure, client=remote) as client:
        assert client.get("/health").status_code == 200


def test_http_task_and_run_lifecycles_match_direct_store_contract(tmp_path: Path) -> None:
    app = create_app(db_path=tmp_path / "ops.db", api_token=TOKEN)
    headers = _headers()
    with TestClient(app) as client:
        finding_id = client.post(
            "/findings",
            headers=headers,
            json={
                "agent_id": "worker",
                "finding_type": "observation",
                "domain": "fixture",
                "summary": "HTTP finding",
                "evidence": {"source": "test"},
                "confidence": 0.8,
                "novelty": 2,
            },
        ).json()["id"]
        filtered = client.get(
            "/findings",
            headers=headers,
            params={"agent_id": "worker", "min_novelty": 2},
        ).json()
        assert [row["id"] for row in filtered] == [finding_id]

        hypothesis_id = client.post(
            "/hypotheses",
            headers=headers,
            json={
                "source_agent_id": "worker",
                "source_finding_ids": [finding_id],
                "hypothesis": "transport parity",
            },
        ).json()["id"]
        updated_hypothesis = client.patch(
            f"/hypotheses/{hypothesis_id}",
            headers=headers,
            json={"status": "supported", "evidence_for": [finding_id]},
        )
        assert updated_hypothesis.status_code == 200
        assert updated_hypothesis.json()["evidence_for"] == [finding_id]

        task_id = client.post(
            "/tasks",
            headers=headers,
            json={
                "created_by": "coordinator",
                "task_type": "survey",
                "description": "exercise HTTP lease lifecycle",
                "idempotency_key": "http-task",
                "lease_seconds": 60,
            },
        ).json()["id"]
        claimed = client.post(
            f"/tasks/{task_id}/claim",
            headers=headers,
            params={"agent_id": "worker"},
        )
        assert claimed.status_code == 200
        assert claimed.json()["status"] == "claimed"
        lease_token = claimed.json()["lease_token"]
        lease_attempt = claimed.json()["lease_attempt"]
        heartbeat = client.post(
            f"/tasks/{task_id}/heartbeat",
            headers=headers,
            json={
                "agent_id": "worker",
                "lease_token": lease_token,
                "lease_attempt": lease_attempt,
            },
        )
        assert heartbeat.json()["status"] == "in_progress"
        completed = client.patch(
            f"/tasks/{task_id}",
            headers=headers,
            json={
                "status": "complete",
                "agent_id": "worker",
                "lease_token": lease_token,
                "lease_attempt": lease_attempt,
                "result_finding_ids": [finding_id],
            },
        )
        assert completed.json()["status"] == "complete"

        run_id = client.post(
            "/runs",
            headers=headers,
            json={
                "created_by": "operator",
                "run_type": "fixture",
                "dataset_path": str(tmp_path / "dataset"),
                "idempotency_key": "http-run",
            },
        ).json()["id"]
        assert client.post(f"/runs/{run_id}/start", headers=headers).json()["status"] == "running"
        assert client.post(f"/runs/{run_id}/heartbeat", headers=headers).status_code == 200
        finished = client.patch(
            f"/runs/{run_id}",
            headers=headers,
            json={"status": "complete", "result": {"ok": True}},
        )
        assert finished.json()["status"] == "complete"
        events = client.get(f"/runs/{run_id}/events", headers=headers).json()
        assert events[-1]["event_type"] == "run_completed"

        stats = client.get("/stats", headers=headers).json()
        assert stats["counts"]["findings"] == 1
        assert stats["counts"]["hypotheses"] == 1
        assert stats["counts"]["tasks"] == 1
        assert stats["counts"]["runs"] == 1


class _Response:
    def __init__(self, payload: dict):
        self._payload = payload

    def raise_for_status(self) -> None:
        return None

    def json(self) -> dict:
        return self._payload


class _RecordingSession:
    def __init__(self):
        self.headers: dict[str, str] = {}
        self.calls: list[tuple[str, str, dict]] = []

    def request(self, method: str, url: str, **kwargs):
        self.calls.append((method, url, kwargs))
        return _Response({"status": "ok", "counts": {}})


class _FlakyMutationSession(_RecordingSession):
    def request(self, method: str, url: str, **kwargs):
        self.calls.append((method, url, kwargs))
        if len(self.calls) == 1:
            return _Response({"status": "ok"})
        if len(self.calls) == 2:
            raise requests.Timeout("response lost after commit")
        return _Response({"id": "stable-finding-id"})


def test_http_client_centralizes_auth_and_bounded_timeouts() -> None:
    session = _RecordingSession()
    client = SharurOps(
        "http://ops.example",
        agent_id="worker",
        api_token=TOKEN,
        session=session,
    )
    client.stats()

    assert session.headers["Authorization"] == f"Bearer {TOKEN}"
    assert [(method, url) for method, url, _ in session.calls] == [
        ("GET", "http://ops.example/health"),
        ("GET", "http://ops.example/stats"),
    ]
    assert all(call[2]["timeout"] == DEFAULT_TIMEOUT for call in session.calls)


def test_default_client_retries_only_idempotent_methods() -> None:
    client = SharurOps("http://127.0.0.1:9", verify_connection=False)
    try:
        retry = client._session.get_adapter("http://").max_retries
        assert retry.allowed_methods == frozenset({"GET", "HEAD", "OPTIONS"})
        assert "POST" not in retry.allowed_methods
    finally:
        client.close()


def test_client_retries_keyed_mutation_after_lost_response() -> None:
    session = _FlakyMutationSession()
    client = SharurOps(
        "http://ops.example",
        agent_id="worker",
        session=session,
        retries=1,
        backoff_factor=0,
    )

    finding_id = client.finding(
        "observation",
        "fixture",
        "safe replay",
        idempotency_key="finding:stable",
    )

    assert finding_id == "stable-finding-id"
    post_calls = [call for call in session.calls if call[0] == "POST"]
    assert len(post_calls) == 2
    assert all(
        call[2]["json"]["idempotency_key"] == "finding:stable"
        for call in post_calls
    )


def test_registered_worker_identity_roles_and_server_clock_recovery(tmp_path: Path) -> None:
    app = create_app(
        db_path=tmp_path / "ops.db",
        api_token=TOKEN,
        maintenance_interval_seconds=60,
    )
    with TestClient(app) as client:
        worker = client.post(
            "/agents",
            headers=_headers(),
            json={
                "agent_id": "worker",
                "role": "worker",
                "capabilities": ["duckdb"],
                "capacity_cpu_slots": 2,
            },
        ).json()
        reader = client.post(
            "/agents",
            headers=_headers(),
            json={"agent_id": "reader", "role": "reader"},
        ).json()
        worker_headers = {"Authorization": f"Bearer {worker['token']}"}
        reader_headers = {"Authorization": f"Bearer {reader['token']}"}

        spoof = client.post(
            "/findings",
            headers=worker_headers,
            json={
                "agent_id": "another-agent",
                "finding_type": "observation",
                "domain": "fixture",
                "summary": "spoof attempt",
            },
        )
        assert spoof.status_code == 403

        created = client.post(
            "/findings",
            headers=worker_headers,
            json={
                "finding_type": "observation",
                "domain": "fixture",
                "summary": "credential-derived identity",
            },
        )
        assert created.status_code == 201
        assert created.json()["finding"]["agent_id"] == "worker"
        assert client.post(
            "/tasks",
            headers=worker_headers,
            json={"task_type": "survey", "description": "forbidden"},
        ).status_code == 403
        assert client.get("/findings", headers=reader_headers).status_code == 200
        assert client.post(
            "/findings",
            headers=reader_headers,
            json={
                "finding_type": "observation",
                "domain": "fixture",
                "summary": "reader write",
            },
        ).status_code == 403

        task_id = client.post(
            "/tasks",
            headers=_headers(),
            json={
                "task_type": "query",
                "description": "live lease",
                "required_capabilities": ["duckdb"],
                "resource_request": {"cpu_slots": 1},
                "lease_seconds": 60,
            },
        ).json()["id"]
        claim = client.post(
            "/tasks/claim-next",
            headers=worker_headers,
            json={},
        )
        assert claim.json()["id"] == task_id
        assert client.post(
            "/tasks/recover",
            headers=worker_headers,
        ).status_code == 403

        # Caller-supplied time is absent from the contract and has no effect.
        recovered = client.post(
            "/tasks/recover",
            headers=_headers(),
            params={"now": 10**15},
        ).json()
        assert recovered == {"recovered": [], "exhausted": []}
        active = client.get(
            "/tasks",
            headers=_headers(),
            params={"assigned_to": "worker"},
        ).json()[0]
        assert active["status"] == "claimed"


def test_durable_event_cursor_and_request_size_limit(tmp_path: Path) -> None:
    app = create_app(
        db_path=tmp_path / "ops.db",
        api_token=TOKEN,
        max_request_bytes=512,
        maintenance_interval_seconds=60,
    )
    with TestClient(app) as client:
        baseline = client.get("/events", headers=_headers()).json()
        cursor = baseline[-1]["id"]
        created = client.post(
            "/findings",
            headers=_headers(),
            json={
                "agent_id": "worker",
                "finding_type": "observation",
                "domain": "fixture",
                "summary": "durable event",
            },
        )
        assert created.status_code == 201
        replay = client.get(
            "/events",
            headers=_headers(),
            params={"after_id": cursor},
        ).json()
        assert [event["event_type"] for event in replay] == ["finding_created"]
        assert replay[0]["entity_id"] == created.json()["id"]

        oversized = client.post(
            "/findings",
            headers=_headers(),
            json={
                "agent_id": "worker",
                "finding_type": "observation",
                "domain": "fixture",
                "summary": "x" * 2_000,
            },
        )
        assert oversized.status_code == 413

        chunked = client.post(
            "/findings",
            headers={**_headers(), "Content-Type": "application/json"},
            content=iter([b"{" + b"x" * 600, b"x" * 600 + b"}"]),
        )
        assert chunked.status_code == 413


def test_event_replay_survives_server_restart(tmp_path: Path) -> None:
    db_path = tmp_path / "ops.db"
    first_app = create_app(db_path=db_path, api_token=TOKEN)
    with TestClient(first_app) as client:
        cursor = client.get("/events", headers=_headers()).json()[-1]["id"]
        finding_id = client.post(
            "/findings",
            headers=_headers(),
            json={
                "agent_id": "worker",
                "finding_type": "observation",
                "domain": "fixture",
                "summary": "restart-safe replay",
            },
        ).json()["id"]

    second_app = create_app(db_path=db_path, api_token=TOKEN)
    with TestClient(second_app) as client:
        replay = client.get(
            "/events",
            headers=_headers(),
            params={"after_id": cursor},
        ).json()
    assert [event["entity_id"] for event in replay] == [finding_id]


def test_sse_wakeup_overflow_is_visible_and_durable_feed_remains_authoritative() -> None:
    async def exercise() -> tuple[dict[str, int], str]:
        bus = EventBus(queue_size=1)
        subscriber = bus.subscribe()
        bus.publish("first")
        bus.publish("second")
        await asyncio.sleep(0)
        return bus.stats(), subscriber.queue.get_nowait()

    stats, queued = asyncio.run(exercise())
    assert stats == {"subscribers": 1, "dropped_wakeups": 1}
    assert '"type":"second"' in queued


def test_metrics_and_backups_expose_operational_health(tmp_path: Path) -> None:
    backup_dir = tmp_path / "backups"
    app = create_app(
        db_path=tmp_path / "ops.db",
        api_token=TOKEN,
        backup_dir=backup_dir,
        backup_interval_seconds=3600,
        maintenance_interval_seconds=60,
    )
    with TestClient(app) as client:
        assert client.post(
            "/findings",
            headers=_headers(),
            json={
                "agent_id": "worker",
                "finding_type": "observation",
                "domain": "fixture",
                "summary": "instrumented write",
            },
        ).status_code == 201
        first_backup = client.post("/admin/backup", headers=_headers()).json()
        second_backup = client.post("/admin/backup", headers=_headers()).json()
        assert first_backup["path"] != second_backup["path"]

        stats = client.get("/stats", headers=_headers()).json()
        assert stats["backup"]["last_success_ts"] is not None
        assert stats["backup"]["age_seconds"] >= 0
        assert stats["dead_letters"]["count"] == 0
        assert stats["pool"]["sqlite_write_transactions"] >= 1
        assert stats["http"]["latency_window_requests"] >= 2

        metrics = client.get("/metrics", headers=_headers()).text
        assert "sharur_ops_http_p95_duration_seconds" in metrics
        assert "sharur_ops_backup_age_seconds" in metrics
        assert "sharur_ops_dead_letter_tasks 0" in metrics
        assert "sharur_ops_pool_sqlite_write_wait_seconds" in metrics


def test_only_one_http_server_can_own_a_database(tmp_path: Path) -> None:
    db_path = tmp_path / "ops.db"
    first = create_app(db_path=db_path, api_token=TOKEN)
    second = create_app(db_path=db_path, api_token=TOKEN)

    with (
        TestClient(first),
        pytest.raises(RuntimeError, match="already has an HTTP owner"),
        TestClient(second),
    ):
        pass

    with TestClient(first), pytest.raises(
        RuntimeError,
        match="Direct SQLite access is unavailable",
    ):
        OpsStore(db_path)
