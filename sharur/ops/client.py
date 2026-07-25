"""Bounded HTTP client for the Sharur Ops control plane."""

from __future__ import annotations

import time
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
        api_token: str | None = None,
        timeout: float | tuple[float, float] = DEFAULT_TIMEOUT,
        retries: int = 2,
        backoff_factor: float = 0.25,
        session: requests.Session | None = None,
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
        self._mutation_retries = retries
        self._backoff_factor = backoff_factor
        self._lease_tokens: dict[str, tuple[str, int]] = {}
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

    def _request(
        self,
        method: str,
        path: str,
        *,
        retry_mutation: bool = False,
        **kwargs: Any,
    ) -> requests.Response:
        """Send a bounded request.

        GET/HEAD retries live in the HTTP adapter. A mutation retries only when
        its caller supplied an idempotency key (or another natural immutable
        identity, such as an artifact content hash).
        """

        url = path if path.startswith(("http://", "https://")) else f"{self.base}{path}"
        kwargs.setdefault("timeout", self.timeout)
        attempts = self._mutation_retries + 1 if retry_mutation else 1
        for attempt in range(attempts):
            try:
                response = self._session.request(method, url, **kwargs)
                response.raise_for_status()
                return response
            except (requests.ConnectionError, requests.Timeout):
                if attempt + 1 >= attempts:
                    raise
                time.sleep(self._backoff_factor * (2**attempt))
        raise AssertionError("unreachable")

    def _remember_lease(self, task: dict[str, Any]) -> dict[str, Any]:
        token = task.get("lease_token")
        attempt = task.get("lease_attempt")
        if token is not None and attempt is not None:
            self._lease_tokens[str(task["id"])] = (str(token), int(attempt))
        return task

    def _lease(
        self,
        task_id: str,
        lease_token: str | None,
        lease_attempt: int | None,
    ) -> tuple[str, int]:
        cached = self._lease_tokens.get(task_id)
        token = lease_token or (cached[0] if cached else None)
        attempt = lease_attempt if lease_attempt is not None else (cached[1] if cached else None)
        if token is None or attempt is None:
            raise ValueError(
                f"Task {task_id} requires the lease token returned by claim_task "
                "or claim_next_task"
            )
        return token, int(attempt)

    # ------------------------------------------------------------------
    # Campaigns and agents
    # ------------------------------------------------------------------

    def create_campaign(
        self,
        name: str,
        *,
        description: str = "",
        dataset_path: str | None = None,
        metadata: dict | None = None,
        idempotency_key: str | None = None,
    ) -> str:
        response = self._request(
            "POST",
            "/campaigns",
            json={
                "name": name,
                "description": description,
                "dataset_path": dataset_path,
                "metadata": metadata or {},
                "idempotency_key": idempotency_key,
            },
            retry_mutation=idempotency_key is not None,
        )
        return response.json()["id"]

    def list_campaigns(self, **filters: Any) -> list[dict[str, Any]]:
        return self._request("GET", "/campaigns", params=filters).json()

    def get_campaign(self, campaign_id: str) -> dict[str, Any]:
        return self._request("GET", f"/campaigns/{campaign_id}").json()

    def update_campaign(self, campaign_id: str, status: str) -> dict[str, Any]:
        return self._request(
            "PATCH",
            f"/campaigns/{campaign_id}",
            json={"status": status},
        ).json()

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
        metadata: dict | None = None,
        rotate_token: bool = False,
    ) -> dict[str, Any]:
        return self._request(
            "POST",
            "/agents",
            json={
                "agent_id": agent_id,
                "role": role,
                "capabilities": capabilities or [],
                "max_concurrent_tasks": max_concurrent_tasks,
                "capacity_cpu_slots": capacity_cpu_slots,
                "capacity_memory_mb": capacity_memory_mb,
                "capacity_accelerator_slots": capacity_accelerator_slots,
                "host": host,
                "metadata": metadata or {},
                "rotate_token": rotate_token,
            },
        ).json()

    def list_agents(self, **filters: Any) -> list[dict[str, Any]]:
        return self._request("GET", "/agents", params=filters).json()

    def heartbeat_agent(
        self,
        *,
        status: str = "active",
        host: str | None = None,
    ) -> dict[str, Any]:
        return self._request(
            "POST",
            "/agents/me/heartbeat",
            json={"status": status, "host": host},
        ).json()

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
        response = self._request(
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
                "campaign_id": campaign_id,
                "task_id": task_id,
                "idempotency_key": idempotency_key,
                "schema_version": schema_version,
                "validation_status": validation_status,
            },
            retry_mutation=idempotency_key is not None,
        )
        return response.json()["id"]

    def recent_findings(
        self,
        since: float = 0,
        min_novelty: int = 0,
        finding_type: str | None = None,
        domain: str | None = None,
        agent_id: str | None = None,
        campaign_id: str | None = None,
        before_ts: float | None = None,
        limit: int = 50,
    ) -> list[dict[str, Any]]:
        params: dict[str, Any] = {
            "since": since,
            "min_novelty": min_novelty,
            "limit": limit,
        }
        for key, value in (
            ("finding_type", finding_type),
            ("domain", domain),
            ("agent_id", agent_id),
            ("campaign_id", campaign_id),
            ("before_ts", before_ts),
        ):
            if value is not None:
                params[key] = value
        return self._request("GET", "/findings", params=params).json()

    def get_finding(self, finding_id: str) -> dict[str, Any]:
        return self._request("GET", f"/findings/{finding_id}").json()

    def search_findings(
        self,
        text: str,
        limit: int = 20,
        *,
        campaign_id: str | None = None,
    ) -> list[dict[str, Any]]:
        params: dict[str, Any] = {"limit": limit}
        if campaign_id is not None:
            params["campaign_id"] = campaign_id
        return self._request(
            "GET",
            f"/findings/search/{text}",
            params=params,
        ).json()

    def link_findings(
        self,
        finding_id: str,
        related_finding_id: str,
        *,
        relation: str,
    ) -> dict[str, Any]:
        return self._request(
            "POST",
            f"/findings/{finding_id}/links",
            json={
                "related_finding_id": related_finding_id,
                "relation": relation,
            },
        ).json()

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
        metadata: dict | None = None,
    ) -> str:
        response = self._request(
            "POST",
            "/artifacts",
            json={
                "content_hash": content_hash,
                "uri": uri,
                "size_bytes": size_bytes,
                "media_type": media_type,
                "campaign_id": campaign_id,
                "task_id": task_id,
                "run_id": run_id,
                "metadata": metadata or {},
            },
            retry_mutation=True,
        )
        return response.json()["id"]

    def attach_artifact(
        self,
        finding_id: str,
        artifact_id: str,
        *,
        relation: str = "evidence",
    ) -> dict[str, Any]:
        return self._request(
            "POST",
            f"/findings/{finding_id}/artifacts",
            json={"artifact_id": artifact_id, "relation": relation},
            retry_mutation=True,
        ).json()

    # ------------------------------------------------------------------
    # Typed candidates and scientific review DAG
    # ------------------------------------------------------------------

    def record_unit_disposition(
        self,
        *,
        campaign_id: str,
        unit_id: str,
        dataset_id: str,
        genome_id: str,
        coverage_hash: str,
        candidate_count: int,
        disposition: str,
        evidence_bundle_hash: str,
        task_id: str | None = None,
        reason_codes: list[str] | None = None,
        strata: dict | None = None,
        provenance: dict | None = None,
        supersedes_disposition_id: str | None = None,
        idempotency_key: str,
        schema_version: int = 1,
    ) -> str:
        response = self._request(
            "POST",
            "/review/unit-dispositions",
            json={
                "agent_id": self.agent_id,
                "campaign_id": campaign_id,
                "unit_id": unit_id,
                "dataset_id": dataset_id,
                "genome_id": genome_id,
                "coverage_hash": coverage_hash,
                "candidate_count": candidate_count,
                "disposition": disposition,
                "evidence_bundle_hash": evidence_bundle_hash,
                "task_id": task_id,
                "reason_codes": reason_codes or [],
                "strata": strata or {},
                "provenance": provenance or {},
                "supersedes_disposition_id": supersedes_disposition_id,
                "idempotency_key": idempotency_key,
                "schema_version": schema_version,
            },
            retry_mutation=True,
        )
        return str(response.json()["id"])

    def list_unit_dispositions(self, **filters: Any) -> list[dict[str, Any]]:
        return self._request(
            "GET", "/review/unit-dispositions", params=filters
        ).json()

    def get_unit_disposition(self, disposition_id: str) -> dict[str, Any]:
        return self._request(
            "GET", f"/review/unit-dispositions/{disposition_id}"
        ).json()

    def create_candidate_occurrence(
        self,
        *,
        campaign_id: str,
        dataset_id: str,
        unit_id: str,
        genome_id: str,
        candidate_type: str,
        signature_schema: str,
        signature: dict,
        evidence: dict,
        verification: list[dict],
        subject_refs: dict,
        task_id: str | None = None,
        reason_codes: list[str] | None = None,
        uncertainty: dict | None = None,
        reduction_features: dict | None = None,
        provenance: dict | None = None,
        evidence_bundle_hash: str | None = None,
        idempotency_key: str,
        schema_version: int = 1,
    ) -> str:
        response = self._request(
            "POST",
            "/review/candidates",
            json={
                "agent_id": self.agent_id,
                "campaign_id": campaign_id,
                "dataset_id": dataset_id,
                "unit_id": unit_id,
                "genome_id": genome_id,
                "candidate_type": candidate_type,
                "signature_schema": signature_schema,
                "signature": signature,
                "evidence": evidence,
                "verification": verification,
                "subject_refs": subject_refs,
                "task_id": task_id,
                "reason_codes": reason_codes or [],
                "uncertainty": uncertainty or {},
                "reduction_features": reduction_features or {},
                "provenance": provenance or {},
                "evidence_bundle_hash": evidence_bundle_hash,
                "idempotency_key": idempotency_key,
                "schema_version": schema_version,
            },
            retry_mutation=True,
        )
        return str(response.json()["id"])

    def list_candidate_occurrences(self, **filters: Any) -> list[dict[str, Any]]:
        return self._request("GET", "/review/candidates", params=filters).json()

    def get_candidate_occurrence(self, candidate_id: str) -> dict[str, Any]:
        return self._request(
            "GET", f"/review/candidates/{candidate_id}"
        ).json()

    def reduce_candidates(
        self,
        campaign_id: str,
        *,
        dataset_id: str | None = None,
        candidate_type: str | None = None,
        batch_size: int = 1_000,
    ) -> dict[str, Any]:
        return self._request(
            "POST",
            "/review/reduce",
            json={
                "campaign_id": campaign_id,
                "dataset_id": dataset_id,
                "candidate_type": candidate_type,
                "batch_size": batch_size,
            },
        ).json()

    def list_candidate_clusters(self, **filters: Any) -> list[dict[str, Any]]:
        return self._request("GET", "/review/clusters", params=filters).json()

    def get_candidate_cluster(
        self,
        cluster_id: str,
        *,
        member_limit: int = 100,
    ) -> dict[str, Any]:
        return self._request(
            "GET",
            f"/review/clusters/{cluster_id}",
            params={"member_limit": member_limit},
        ).json()

    def list_candidate_cluster_members(
        self,
        cluster_id: str,
        *,
        after_candidate_id: str | None = None,
        limit: int = 500,
    ) -> dict[str, Any]:
        params: dict[str, Any] = {"limit": limit}
        if after_candidate_id is not None:
            params["after_candidate_id"] = after_candidate_id
        return self._request(
            "GET",
            f"/review/clusters/{cluster_id}/members",
            params=params,
        ).json()

    def link_cluster_finding(
        self,
        cluster_id: str,
        finding_id: str,
        *,
        relation: str = "materializes",
    ) -> dict[str, Any]:
        return self._request(
            "POST",
            f"/review/clusters/{cluster_id}/findings",
            json={"finding_id": finding_id, "relation": relation},
        ).json()

    def list_cluster_findings(self, **filters: Any) -> list[dict[str, Any]]:
        return self._request(
            "GET", "/review/cluster-findings", params=filters
        ).json()

    def create_finding_review(
        self,
        *,
        campaign_id: str,
        dataset_id: str,
        review_tier: str,
        execution_profile: str,
        provider: str,
        model: str,
        reasoning_effort: str,
        prompt_hash: str,
        rubric_version: str,
        input_bundle_hash: str,
        verdict: str,
        confidence: float,
        finding_id: str | None = None,
        cluster_id: str | None = None,
        unit_disposition_id: str | None = None,
        task_id: str | None = None,
        run_id: str | None = None,
        reconstructed_observations: dict | None = None,
        claim_assessment: dict | None = None,
        verification_summary: dict | None = None,
        discrepancies: list[dict] | None = None,
        proposed_tasks: list[dict] | None = None,
        blind_to_prior_scores: bool = True,
        blind_to_other_reviews: bool = True,
        parent_review_id: str | None = None,
        idempotency_key: str,
        schema_version: int = 1,
    ) -> str:
        response = self._request(
            "POST",
            "/review/reviews",
            json={
                "reviewer_agent_id": self.agent_id,
                "campaign_id": campaign_id,
                "dataset_id": dataset_id,
                "review_tier": review_tier,
                "execution_profile": execution_profile,
                "provider": provider,
                "model": model,
                "reasoning_effort": reasoning_effort,
                "prompt_hash": prompt_hash,
                "rubric_version": rubric_version,
                "input_bundle_hash": input_bundle_hash,
                "verdict": verdict,
                "confidence": confidence,
                "finding_id": finding_id,
                "cluster_id": cluster_id,
                "unit_disposition_id": unit_disposition_id,
                "task_id": task_id,
                "run_id": run_id,
                "reconstructed_observations": reconstructed_observations or {},
                "claim_assessment": claim_assessment or {},
                "verification_summary": verification_summary or {},
                "discrepancies": discrepancies or [],
                "proposed_tasks": proposed_tasks or [],
                "blind_to_prior_scores": blind_to_prior_scores,
                "blind_to_other_reviews": blind_to_other_reviews,
                "parent_review_id": parent_review_id,
                "idempotency_key": idempotency_key,
                "schema_version": schema_version,
            },
            retry_mutation=True,
        )
        return str(response.json()["id"])

    def list_finding_reviews(self, **filters: Any) -> list[dict[str, Any]]:
        return self._request("GET", "/review/reviews", params=filters).json()

    def get_finding_review(self, review_id: str) -> dict[str, Any]:
        return self._request("GET", f"/review/reviews/{review_id}").json()

    def record_review_verification(
        self,
        review_id: str,
        *,
        claim_key: str,
        engine: str,
        specification: dict,
        dataset_id: str,
        expected: Any,
        status: str,
        actual: Any = None,
        executed_ts: float | None = None,
        code_commit: str | None = None,
        artifact_id: str | None = None,
        error: str | None = None,
        supersedes_verification_id: str | None = None,
        idempotency_key: str,
    ) -> str:
        response = self._request(
            "POST",
            f"/review/reviews/{review_id}/verifications",
            json={
                "agent_id": self.agent_id,
                "claim_key": claim_key,
                "engine": engine,
                "specification": specification,
                "dataset_id": dataset_id,
                "expected": expected,
                "status": status,
                "actual": actual,
                "executed_ts": executed_ts,
                "code_commit": code_commit,
                "artifact_id": artifact_id,
                "error": error,
                "supersedes_verification_id": supersedes_verification_id,
                "idempotency_key": idempotency_key,
            },
            retry_mutation=True,
        )
        return str(response.json()["id"])

    def list_review_verifications(
        self, review_id: str
    ) -> list[dict[str, Any]]:
        return self._request(
            "GET", f"/review/reviews/{review_id}/verifications"
        ).json()

    def create_promotion_decision(self, **payload: Any) -> dict[str, Any]:
        return self._request(
            "POST",
            "/review/decisions",
            json=payload,
            retry_mutation=bool(payload.get("idempotency_key")),
        ).json()

    def list_promotion_decisions(self, **filters: Any) -> list[dict[str, Any]]:
        return self._request("GET", "/review/decisions", params=filters).json()

    def get_promotion_decision(self, decision_id: str) -> dict[str, Any]:
        return self._request("GET", f"/review/decisions/{decision_id}").json()

    def record_canonical_publication(self, **payload: Any) -> str:
        response = self._request(
            "POST",
            "/review/publications",
            json=payload,
            retry_mutation=bool(payload.get("idempotency_key")),
        )
        return str(response.json()["id"])

    def list_canonical_publications(self, **filters: Any) -> list[dict[str, Any]]:
        return self._request(
            "GET", "/review/publications", params=filters
        ).json()

    def tick_review_controller(
        self,
        campaign_id: str,
        *,
        policy: dict | None = None,
    ) -> dict[str, Any]:
        return self._request(
            "POST",
            "/review/controller/tick",
            json={"campaign_id": campaign_id, "policy": policy},
        ).json()

    def review_status(self, campaign_id: str) -> dict[str, Any]:
        return self._request(
            "GET",
            "/review/status",
            params={"campaign_id": campaign_id},
        ).json()

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
        response = self._request(
            "POST",
            "/hypotheses",
            json={
                "source_agent_id": self.agent_id,
                "source_finding_ids": source_finding_ids or [],
                "hypothesis": hypothesis,
                "domains_relevant": domains_relevant or [],
                "campaign_id": campaign_id,
                "idempotency_key": idempotency_key,
                "schema_version": schema_version,
            },
            retry_mutation=idempotency_key is not None,
        )
        return response.json()["id"]

    def open_hypotheses(
        self,
        unassigned: bool = True,
        *,
        campaign_id: str | None = None,
    ) -> list[dict[str, Any]]:
        params: dict[str, Any] = {
            "status": "proposed",
            "unassigned": unassigned,
        }
        if campaign_id is not None:
            params["campaign_id"] = campaign_id
        return self._request("GET", "/hypotheses", params=params).json()

    def update_hypothesis(self, hyp_id: str, **kwargs: Any) -> dict[str, Any]:
        return self._request(
            "PATCH",
            f"/hypotheses/{hyp_id}",
            json=kwargs,
        ).json()

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
        resource_request: dict[str, int] | None = None,
        max_attempts: int = 3,
        lease_seconds: int = 900,
    ) -> str:
        response = self._request(
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
                "campaign_id": campaign_id,
                "idempotency_key": idempotency_key,
                "depends_on": depends_on or [],
                "required_capabilities": required_capabilities or [],
                "resource_request": resource_request or {},
                "max_attempts": max_attempts,
                "lease_seconds": lease_seconds,
            },
            retry_mutation=idempotency_key is not None,
        )
        return response.json()["id"]

    def my_tasks(
        self,
        status: str | None = None,
        *,
        campaign_id: str | None = None,
        limit: int = 50,
    ) -> list[dict[str, Any]]:
        params: dict[str, Any] = {
            "assigned_to": self.agent_id,
            "limit": limit,
        }
        if status is not None:
            params["status"] = status
        if campaign_id is not None:
            params["campaign_id"] = campaign_id
        return self._request("GET", "/tasks", params=params).json()

    def available_tasks(
        self,
        *,
        campaign_id: str | None = None,
        limit: int = 50,
    ) -> list[dict[str, Any]]:
        params: dict[str, Any] = {"unassigned": True, "limit": limit}
        if campaign_id is not None:
            params["campaign_id"] = campaign_id
        return self._request("GET", "/tasks", params=params).json()

    def claim_task(
        self,
        task_id: str,
        lease_seconds: int | None = None,
    ) -> dict[str, Any]:
        task = self._request(
            "POST",
            f"/tasks/{task_id}/claim",
            json={
                "agent_id": self.agent_id,
                "lease_seconds": lease_seconds,
            },
        ).json()
        return self._remember_lease(task)

    def claim_next_task(
        self,
        *,
        campaign_id: str | None = None,
        task_types: list[str] | None = None,
        lease_seconds: int | None = None,
    ) -> dict[str, Any] | None:
        task = self._request(
            "POST",
            "/tasks/claim-next",
            json={
                "agent_id": self.agent_id,
                "campaign_id": campaign_id,
                "task_types": task_types,
                "lease_seconds": lease_seconds,
            },
        ).json()
        return self._remember_lease(task) if task is not None else None

    def heartbeat_task(
        self,
        task_id: str,
        lease_seconds: int | None = None,
        *,
        lease_token: str | None = None,
        lease_attempt: int | None = None,
    ) -> dict[str, Any]:
        token, attempt = self._lease(task_id, lease_token, lease_attempt)
        return self._request(
            "POST",
            f"/tasks/{task_id}/heartbeat",
            json={
                "agent_id": self.agent_id,
                "lease_token": token,
                "lease_attempt": attempt,
                "lease_seconds": lease_seconds,
            },
        ).json()

    def put_task_checkpoint(
        self,
        task_id: str,
        checkpoint_key: str,
        *,
        cursor: str | None = None,
        payload: dict[str, Any] | None = None,
        lease_token: str | None = None,
        lease_attempt: int | None = None,
    ) -> dict[str, Any]:
        """Persist retry-visible progress under the active task lease."""
        token, attempt = self._lease(task_id, lease_token, lease_attempt)
        return self._request(
            "PUT",
            f"/tasks/{task_id}/checkpoint",
            json={
                "agent_id": self.agent_id,
                "checkpoint_key": checkpoint_key,
                "cursor": cursor,
                "payload": payload or {},
                "lease_token": token,
                "lease_attempt": attempt,
            },
        ).json()

    def get_task_checkpoint(
        self,
        task_id: str,
        checkpoint_key: str,
    ) -> dict[str, Any]:
        """Read the latest checkpoint, including progress from a prior attempt."""
        return self._request(
            "GET",
            f"/tasks/{task_id}/checkpoints",
            params={"checkpoint_key": checkpoint_key},
        ).json()

    def list_task_checkpoints(
        self,
        task_id: str,
        *,
        limit: int = 50,
    ) -> list[dict[str, Any]]:
        return self._request(
            "GET",
            f"/tasks/{task_id}/checkpoints",
            params={"limit": limit},
        ).json()

    def complete_task(
        self,
        task_id: str,
        result_finding_ids: list[str] | None = None,
        *,
        lease_token: str | None = None,
        lease_attempt: int | None = None,
    ) -> dict[str, Any]:
        token, attempt = self._lease(task_id, lease_token, lease_attempt)
        result = self._request(
            "PATCH",
            f"/tasks/{task_id}",
            json={
                "status": "complete",
                "agent_id": self.agent_id,
                "lease_token": token,
                "lease_attempt": attempt,
                "result_finding_ids": result_finding_ids or [],
            },
        ).json()
        self._lease_tokens.pop(task_id, None)
        return result

    def fail_task(
        self,
        task_id: str,
        *,
        error: str | None = None,
        retryable: bool = False,
        retry_delay_seconds: int = 0,
        lease_token: str | None = None,
        lease_attempt: int | None = None,
    ) -> dict[str, Any]:
        token, attempt = self._lease(task_id, lease_token, lease_attempt)
        result = self._request(
            "PATCH",
            f"/tasks/{task_id}",
            json={
                "status": "failed",
                "agent_id": self.agent_id,
                "lease_token": token,
                "lease_attempt": attempt,
                "error": error,
                "retryable": retryable,
                "retry_delay_seconds": retry_delay_seconds,
            },
        ).json()
        self._lease_tokens.pop(task_id, None)
        return result

    def reset_failed_tasks(
        self,
        *,
        campaign_id: str | None = None,
        task_ids: list[str] | None = None,
        only_transient: bool = False,
        extra_attempts: int = 5,
    ) -> dict[str, Any]:
        """Return attempt-exhausted tasks to the queue."""
        body: dict[str, Any] = {
            "only_transient": only_transient,
            "extra_attempts": extra_attempts,
        }
        if campaign_id is not None:
            body["campaign_id"] = campaign_id
        if task_ids:
            body["task_ids"] = task_ids
        return self._request("POST", "/tasks/reset-failed", json=body).json()

    def recover_expired_tasks(self, now: float | None = None) -> dict[str, Any]:
        if now is not None:
            raise ValueError("HTTP recovery uses the server clock")
        return self._request("POST", "/tasks/recover").json()

    # ------------------------------------------------------------------
    # Runs
    # ------------------------------------------------------------------

    def create_run(
        self,
        run_type: str,
        dataset_path: str,
        *,
        config: dict | None = None,
        idempotency_key: str | None = None,
        parent_run_id: str | None = None,
        campaign_id: str | None = None,
    ) -> str:
        response = self._request(
            "POST",
            "/runs",
            json={
                "created_by": self.agent_id,
                "run_type": run_type,
                "dataset_path": dataset_path,
                "config": config or {},
                "idempotency_key": idempotency_key,
                "parent_run_id": parent_run_id,
                "campaign_id": campaign_id,
            },
            retry_mutation=idempotency_key is not None,
        )
        return response.json()["id"]

    def start_run(self, run_id: str) -> dict[str, Any]:
        return self._request("POST", f"/runs/{run_id}/start").json()

    def submit_run(self, run_id: str) -> dict[str, Any]:
        return self._request("POST", f"/runs/{run_id}/submit").json()

    def heartbeat_run(self, run_id: str) -> dict[str, Any]:
        return self._request("POST", f"/runs/{run_id}/heartbeat").json()

    def complete_run(self, run_id: str, result: dict | None = None) -> dict[str, Any]:
        return self._request(
            "PATCH",
            f"/runs/{run_id}",
            json={"status": "complete", "result": result or {}},
        ).json()

    def fail_run(self, run_id: str, error: str) -> dict[str, Any]:
        return self._request(
            "PATCH",
            f"/runs/{run_id}",
            json={"status": "failed", "error": error},
        ).json()

    def get_run(self, run_id: str) -> dict[str, Any]:
        return self._request("GET", f"/runs/{run_id}").json()

    def list_runs(self, **filters: Any) -> list[dict[str, Any]]:
        return self._request("GET", "/runs", params=filters).json()

    def run_events(
        self,
        run_id: str,
        *,
        after_id: int = 0,
        limit: int = 1_000,
    ) -> list[dict[str, Any]]:
        return self._request(
            "GET",
            f"/runs/{run_id}/events",
            params={"after_id": after_id, "limit": limit},
        ).json()

    def run_stages(
        self,
        run_id: str,
        stage_id: str | None = None,
    ) -> list[dict[str, Any]]:
        params = {"stage_id": stage_id} if stage_id is not None else None
        return self._request(
            "GET",
            f"/runs/{run_id}/stages",
            params=params,
        ).json()

    def recover_stale_runs(
        self,
        *,
        stale_after_seconds: int = 300,
        now: float | None = None,
    ) -> dict[str, Any]:
        if now is not None:
            raise ValueError("HTTP recovery uses the server clock")
        return self._request(
            "POST",
            "/runs/recover",
            params={"stale_after_seconds": stale_after_seconds},
        ).json()

    # ------------------------------------------------------------------
    # Coordinator log, events, and diagnostics
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
        response = self._request(
            "POST",
            "/coordinator/log",
            json={
                "action_type": action_type,
                "reasoning": reasoning,
                "referenced_finding_ids": referenced_finding_ids or [],
                "referenced_hypothesis_ids": referenced_hypothesis_ids or [],
                "decisions_made": decisions_made or {},
                "campaign_id": campaign_id,
                "idempotency_key": idempotency_key,
            },
            retry_mutation=idempotency_key is not None,
        )
        return response.json()["id"]

    def recent_log(
        self,
        limit: int = 20,
        since: float = 0,
        *,
        campaign_id: str | None = None,
        before_ts: float | None = None,
    ) -> list[dict[str, Any]]:
        params: dict[str, Any] = {"limit": limit, "since": since}
        if campaign_id is not None:
            params["campaign_id"] = campaign_id
        if before_ts is not None:
            params["before_ts"] = before_ts
        return self._request("GET", "/coordinator/log", params=params).json()

    def events(
        self,
        *,
        after_id: int = 0,
        campaign_id: str | None = None,
        entity_type: str | None = None,
        limit: int = 1_000,
    ) -> list[dict[str, Any]]:
        params: dict[str, Any] = {"after_id": after_id, "limit": limit}
        if campaign_id is not None:
            params["campaign_id"] = campaign_id
        if entity_type is not None:
            params["entity_type"] = entity_type
        return self._request("GET", "/events", params=params).json()

    def backup(self) -> dict[str, Any]:
        return self._request("POST", "/admin/backup").json()

    def integrity_check(self) -> dict[str, Any]:
        return self._request("GET", "/admin/integrity").json()

    def stats(self) -> dict[str, Any]:
        return self._request("GET", "/stats").json()

    def close(self) -> None:
        if self._owns_session:
            self._session.close()

    def __enter__(self) -> SharurOps:
        return self

    def __exit__(self, *_exc_info: object) -> None:
        self.close()

    def __repr__(self) -> str:
        return f"SharurOps(base='{self.base}', agent_id='{self.agent_id}')"


__all__ = ["DEFAULT_TIMEOUT", "SharurOps"]
