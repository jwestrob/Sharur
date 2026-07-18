# ruff: noqa: RUF013
"""
Sharur Ops Client — drop this into any agent's environment.

Usage:
    from sharur.ops.client import SharurOps
    ops = SharurOps("http://localhost:8811", agent_id="dpann_surveyor")

    # Post a finding
    fid = ops.finding(
        finding_type="cassette",
        domain="dpann",
        summary="Novel CHAT protease defense cassette: GHG + ELV-DEIG + PTW",
        evidence={"genes": [...], "conservation_count": 44, "phyla_observed": [...]},
        confidence=0.85,
        novelty=3,
        reasoning="Found via ELSA neighborhood search at tau=0.85..."
    )

    # Check for tasks
    tasks = ops.my_tasks()

    # Claim an unassigned task
    ops.claim_task(task_id)

    # Post a hypothesis for cross-domain validation
    ops.hypothesis(
        hypothesis="CHAT protease cassette is present in Bathyarchaeia",
        source_finding_ids=[fid],
        domains_relevant=["bathyarchaeia", "omnitrophota"]
    )

    # Check what other agents have found
    new_stuff = ops.recent_findings(min_novelty=2)
"""

from typing import Any

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


DEFAULT_TIMEOUT = (3.05, 30.0)


class SharurOps:
    def __init__(
        self,
        base_url: str = "http://localhost:8811",
        agent_id: str = "unnamed",
        *,
        api_token: str = None,
        timeout: float | tuple[float, float] = DEFAULT_TIMEOUT,
        retries: int = 2,
        backoff_factor: float = 0.25,
        session: requests.Session = None,
        verify_connection: bool = True,
    ):
        if retries < 0:
            raise ValueError("retries must be non-negative")
        if isinstance(timeout, tuple):
            if len(timeout) != 2 or any(value <= 0 for value in timeout):
                raise ValueError("timeout tuple must contain positive connect/read values")
        elif timeout <= 0:
            raise ValueError("timeout must be positive")
        self.base = base_url.rstrip("/")
        self.agent_id = agent_id
        self.timeout = timeout
        self._owns_session = session is None
        self._session = session or requests.Session()
        if api_token:
            self._session.headers.update({"Authorization": f"Bearer {api_token}"})
        if self._owns_session:
            retry_policy = Retry(
                total=retries,
                connect=retries,
                read=retries,
                status=retries,
                allowed_methods=frozenset({"GET", "HEAD", "OPTIONS"}),
                status_forcelist=(429, 502, 503, 504),
                backoff_factor=backoff_factor,
                respect_retry_after_header=True,
            )
            adapter = HTTPAdapter(max_retries=retry_policy)
            self._session.mount("http://", adapter)
            self._session.mount("https://", adapter)
        if verify_connection:
            self._request("GET", "/health")

    def _request(self, method: str, path: str, **kwargs: Any) -> requests.Response:
        """Send one bounded request; only idempotent reads are retried."""
        url = path if path.startswith(("http://", "https://")) else f"{self.base}{path}"
        kwargs.setdefault("timeout", self.timeout)
        response = self._session.request(method, url, **kwargs)
        response.raise_for_status()
        return response

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
        """Post a finding. Returns the finding ID."""
        r = self._request(
            "POST",
            "/findings",
            json={
                "agent_id": self.agent_id,
                "finding_type": finding_type,
                "domain": domain,
                "summary": summary,
                "evidence": evidence or {},
                "confidence": confidence,
                "novelty": novelty,
                "parent_finding_id": parent_finding_id,
                "reasoning": reasoning,
            },
        )
        return r.json()["id"]

    def recent_findings(
        self,
        since: float = 0,
        min_novelty: int = 0,
        finding_type: str = None,
        domain: str = None,
        agent_id: str = None,
        limit: int = 50,
    ) -> list[dict]:
        """Query recent findings, optionally filtered."""
        params = {"since": since, "min_novelty": min_novelty, "limit": limit}
        if finding_type:
            params["finding_type"] = finding_type
        if domain:
            params["domain"] = domain
        if agent_id:
            params["agent_id"] = agent_id
        r = self._request("GET", "/findings", params=params)
        return r.json()

    def get_finding(self, finding_id: str) -> dict:
        r = self._request("GET", f"/findings/{finding_id}")
        return r.json()

    def search_findings(self, text: str, limit: int = 20) -> list[dict]:
        """Full-text search across findings."""
        r = self._request("GET", f"/findings/search/{text}", params={"limit": limit})
        return r.json()

    # ------------------------------------------------------------------
    # Hypotheses
    # ------------------------------------------------------------------

    def hypothesis(
        self,
        hypothesis: str,
        source_finding_ids: list[str] = None,
        domains_relevant: list[str] = None,
    ) -> str:
        """Propose a hypothesis. Returns the hypothesis ID."""
        r = self._request(
            "POST",
            "/hypotheses",
            json={
                "source_agent_id": self.agent_id,
                "source_finding_ids": source_finding_ids or [],
                "hypothesis": hypothesis,
                "domains_relevant": domains_relevant or [],
            },
        )
        return r.json()["id"]

    def open_hypotheses(self, unassigned: bool = True) -> list[dict]:
        """Get hypotheses needing investigation."""
        r = self._request(
            "GET",
            "/hypotheses",
            params={"status": "proposed", "unassigned": unassigned},
        )
        return r.json()

    def update_hypothesis(self, hyp_id: str, **kwargs) -> dict:
        """Update hypothesis status/evidence."""
        r = self._request("PATCH", f"/hypotheses/{hyp_id}", json=kwargs)
        return r.json()

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
        *,
        run_id: str = None,
        idempotency_key: str = None,
        depends_on: list[str] = None,
        max_attempts: int = 3,
        lease_seconds: int = 900,
    ) -> str:
        """Create a task (usually called by orchestrator). Returns task ID."""
        r = self._request(
            "POST",
            "/tasks",
            json={
                "created_by": self.agent_id,
                "task_type": task_type,
                "description": description,
                "params": params or {},
                "priority": priority,
                "domain_hint": domain_hint,
                "assigned_to": assigned_to,
                "run_id": run_id,
                "idempotency_key": idempotency_key,
                "depends_on": depends_on or [],
                "max_attempts": max_attempts,
                "lease_seconds": lease_seconds,
            },
        )
        return r.json()["id"]

    def my_tasks(self, status: str = None) -> list[dict]:
        """Get tasks assigned to this agent."""
        params = {"assigned_to": self.agent_id}
        if status:
            params["status"] = status
        r = self._request("GET", "/tasks", params=params)
        return r.json()

    def available_tasks(self) -> list[dict]:
        """Get unassigned pending tasks."""
        r = self._request("GET", "/tasks", params={"unassigned": True})
        return r.json()

    def claim_task(self, task_id: str, lease_seconds: int = None) -> dict:
        """Atomically claim a pending task."""
        params = {"agent_id": self.agent_id}
        if lease_seconds is not None:
            params["lease_seconds"] = lease_seconds
        r = self._request("POST", f"/tasks/{task_id}/claim", params=params)
        return r.json()

    def heartbeat_task(self, task_id: str, lease_seconds: int = None) -> dict:
        """Extend this agent's task lease and mark the task in progress."""
        params = {"agent_id": self.agent_id}
        if lease_seconds is not None:
            params["lease_seconds"] = lease_seconds
        r = self._request(
            "POST",
            f"/tasks/{task_id}/heartbeat",
            params=params,
        )
        return r.json()

    def complete_task(self, task_id: str, result_finding_ids: list[str] = None) -> dict:
        """Mark a task as complete with optional result findings."""
        r = self._request(
            "PATCH",
            f"/tasks/{task_id}",
            json={
                "status": "complete",
                "agent_id": self.agent_id,
                "result_finding_ids": result_finding_ids or [],
            },
        )
        return r.json()

    def fail_task(
        self,
        task_id: str,
        *,
        error: str = None,
        retryable: bool = False,
        retry_delay_seconds: int = 0,
    ) -> dict:
        r = self._request(
            "PATCH",
            f"/tasks/{task_id}",
            json={
                "status": "failed",
                "agent_id": self.agent_id,
                "error": error,
                "retryable": retryable,
                "retry_delay_seconds": retry_delay_seconds,
            },
        )
        return r.json()

    def recover_expired_tasks(self, now: float = None) -> dict:
        params = {}
        if now is not None:
            params["now"] = now
        r = self._request("POST", "/tasks/recover", params=params)
        return r.json()

    # ------------------------------------------------------------------
    # Runs
    # ------------------------------------------------------------------

    def create_run(
        self,
        run_type: str,
        dataset_path: str,
        *,
        config: dict = None,
        idempotency_key: str = None,
        parent_run_id: str = None,
    ) -> str:
        r = self._request(
            "POST",
            "/runs",
            json={
                "created_by": self.agent_id,
                "run_type": run_type,
                "dataset_path": dataset_path,
                "config": config or {},
                "idempotency_key": idempotency_key,
                "parent_run_id": parent_run_id,
            },
        )
        return r.json()["id"]

    def start_run(self, run_id: str) -> dict:
        r = self._request("POST", f"/runs/{run_id}/start")
        return r.json()

    def submit_run(self, run_id: str) -> dict:
        r = self._request("POST", f"/runs/{run_id}/submit")
        return r.json()

    def heartbeat_run(self, run_id: str) -> dict:
        r = self._request("POST", f"/runs/{run_id}/heartbeat")
        return r.json()

    def complete_run(self, run_id: str, result: dict = None) -> dict:
        r = self._request(
            "PATCH",
            f"/runs/{run_id}",
            json={"status": "complete", "result": result or {}},
        )
        return r.json()

    def fail_run(self, run_id: str, error: str) -> dict:
        r = self._request(
            "PATCH",
            f"/runs/{run_id}",
            json={"status": "failed", "error": error},
        )
        return r.json()

    def get_run(self, run_id: str) -> dict:
        r = self._request("GET", f"/runs/{run_id}")
        return r.json()

    def list_runs(self, **filters) -> list[dict]:
        r = self._request("GET", "/runs", params=filters)
        return r.json()

    def run_events(self, run_id: str) -> list[dict]:
        r = self._request("GET", f"/runs/{run_id}/events")
        return r.json()

    def run_stages(self, run_id: str, stage_id: str = None) -> list[dict]:
        params = {"stage_id": stage_id} if stage_id is not None else None
        r = self._request("GET", f"/runs/{run_id}/stages", params=params)
        return r.json()

    def recover_stale_runs(
        self,
        *,
        stale_after_seconds: int = 300,
        now: float = None,
    ) -> dict:
        params = {"stale_after_seconds": stale_after_seconds}
        if now is not None:
            params["now"] = now
        r = self._request("POST", "/runs/recover", params=params)
        return r.json()

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
        """Write a coordinator log entry (primarily for orchestrator use)."""
        r = self._request(
            "POST",
            "/coordinator/log",
            json={
                "action_type": action_type,
                "reasoning": reasoning,
                "referenced_finding_ids": referenced_finding_ids or [],
                "referenced_hypothesis_ids": referenced_hypothesis_ids or [],
                "decisions_made": decisions_made or {},
            },
        )
        return r.json()["id"]

    def recent_log(self, limit: int = 20, since: float = 0) -> list[dict]:
        r = self._request(
            "GET",
            "/coordinator/log",
            params={"limit": limit, "since": since},
        )
        return r.json()

    # ------------------------------------------------------------------
    # Stats
    # ------------------------------------------------------------------

    def stats(self) -> dict:
        r = self._request("GET", "/stats")
        return r.json()

    def close(self) -> None:
        if self._owns_session:
            self._session.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def __repr__(self):
        return f"SharurOps(base='{self.base}', agent_id='{self.agent_id}')"
