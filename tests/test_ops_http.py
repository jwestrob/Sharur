"""Contract tests for the optional Sharur Ops HTTP transport."""

from __future__ import annotations

import subprocess
import sys
from typing import TYPE_CHECKING

from fastapi.testclient import TestClient

from sharur.ops.client import DEFAULT_TIMEOUT, SharurOps
from sharur.ops.server import create_app


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
        heartbeat = client.post(
            f"/tasks/{task_id}/heartbeat",
            headers=headers,
            params={"agent_id": "worker"},
        )
        assert heartbeat.json()["status"] == "in_progress"
        completed = client.patch(
            f"/tasks/{task_id}",
            headers=headers,
            json={
                "status": "complete",
                "agent_id": "worker",
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
