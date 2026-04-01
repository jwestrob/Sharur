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

import requests
import json
from typing import Optional


class SharurOps:
    def __init__(self, base_url: str = "http://localhost:8811", agent_id: str = "unnamed"):
        self.base = base_url.rstrip("/")
        self.agent_id = agent_id
        # Verify connection
        r = requests.get(f"{self.base}/health", timeout=5)
        r.raise_for_status()

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
        r = requests.post(f"{self.base}/findings", json={
            "agent_id": self.agent_id,
            "finding_type": finding_type,
            "domain": domain,
            "summary": summary,
            "evidence": evidence or {},
            "confidence": confidence,
            "novelty": novelty,
            "parent_finding_id": parent_finding_id,
            "reasoning": reasoning,
        })
        r.raise_for_status()
        return r.json()["id"]

    def recent_findings(
        self,
        since: float = 0,
        min_novelty: int = 0,
        finding_type: str = None,
        domain: str = None,
        limit: int = 50,
    ) -> list[dict]:
        """Query recent findings, optionally filtered."""
        params = {"since": since, "min_novelty": min_novelty, "limit": limit}
        if finding_type:
            params["finding_type"] = finding_type
        if domain:
            params["domain"] = domain
        r = requests.get(f"{self.base}/findings", params=params)
        r.raise_for_status()
        return r.json()

    def get_finding(self, finding_id: str) -> dict:
        r = requests.get(f"{self.base}/findings/{finding_id}")
        r.raise_for_status()
        return r.json()

    def search_findings(self, text: str, limit: int = 20) -> list[dict]:
        """Full-text search across findings."""
        r = requests.get(f"{self.base}/findings/search/{text}", params={"limit": limit})
        r.raise_for_status()
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
        r = requests.post(f"{self.base}/hypotheses", json={
            "source_agent_id": self.agent_id,
            "source_finding_ids": source_finding_ids or [],
            "hypothesis": hypothesis,
            "domains_relevant": domains_relevant or [],
        })
        r.raise_for_status()
        return r.json()["id"]

    def open_hypotheses(self, unassigned: bool = True) -> list[dict]:
        """Get hypotheses needing investigation."""
        r = requests.get(f"{self.base}/hypotheses",
                         params={"status": "proposed", "unassigned": unassigned})
        r.raise_for_status()
        return r.json()

    def update_hypothesis(self, hyp_id: str, **kwargs) -> dict:
        """Update hypothesis status/evidence."""
        r = requests.patch(f"{self.base}/hypotheses/{hyp_id}", json=kwargs)
        r.raise_for_status()
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
    ) -> str:
        """Create a task (usually called by orchestrator). Returns task ID."""
        r = requests.post(f"{self.base}/tasks", json={
            "created_by": self.agent_id,
            "task_type": task_type,
            "description": description,
            "params": params or {},
            "priority": priority,
            "domain_hint": domain_hint,
            "assigned_to": assigned_to,
        })
        r.raise_for_status()
        return r.json()["id"]

    def my_tasks(self, status: str = None) -> list[dict]:
        """Get tasks assigned to this agent."""
        params = {"assigned_to": self.agent_id}
        if status:
            params["status"] = status
        r = requests.get(f"{self.base}/tasks", params=params)
        r.raise_for_status()
        return r.json()

    def available_tasks(self) -> list[dict]:
        """Get unassigned pending tasks."""
        r = requests.get(f"{self.base}/tasks", params={"unassigned": True})
        r.raise_for_status()
        return r.json()

    def claim_task(self, task_id: str) -> dict:
        """Atomically claim a pending task."""
        r = requests.post(f"{self.base}/tasks/{task_id}/claim",
                          params={"agent_id": self.agent_id})
        r.raise_for_status()
        return r.json()

    def complete_task(self, task_id: str, result_finding_ids: list[str] = None) -> dict:
        """Mark a task as complete with optional result findings."""
        r = requests.patch(f"{self.base}/tasks/{task_id}", json={
            "status": "complete",
            "result_finding_ids": result_finding_ids or [],
        })
        r.raise_for_status()
        return r.json()

    def fail_task(self, task_id: str) -> dict:
        r = requests.patch(f"{self.base}/tasks/{task_id}", json={"status": "failed"})
        r.raise_for_status()
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
        r = requests.post(f"{self.base}/coordinator/log", json={
            "action_type": action_type,
            "reasoning": reasoning,
            "referenced_finding_ids": referenced_finding_ids or [],
            "referenced_hypothesis_ids": referenced_hypothesis_ids or [],
            "decisions_made": decisions_made or {},
        })
        r.raise_for_status()
        return r.json()["id"]

    def recent_log(self, limit: int = 20, since: float = 0) -> list[dict]:
        r = requests.get(f"{self.base}/coordinator/log",
                         params={"limit": limit, "since": since})
        r.raise_for_status()
        return r.json()

    # ------------------------------------------------------------------
    # Stats
    # ------------------------------------------------------------------

    def stats(self) -> dict:
        r = requests.get(f"{self.base}/stats")
        r.raise_for_status()
        return r.json()

    def __repr__(self):
        return f"SharurOps(base='{self.base}', agent_id='{self.agent_id}')"
