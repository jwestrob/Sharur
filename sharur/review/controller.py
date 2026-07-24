"""Event-driven policy controller for the append-only scientific review DAG."""

from __future__ import annotations

import hashlib
import json
import time
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from sharur.ops.review_store import content_hash


if TYPE_CHECKING:
    from sharur.ops.store import OpsStore
    from sharur.review.models import ReviewPolicy


REVIEW_TASK_CONTRACT = "sharur-review-task/1.0"
MATERIALIZATION_TASK_CONTRACT = "sharur-finding-materialization/1.0"


class ControllerBackpressureError(RuntimeError):
    """A configured execution-profile queue has reached its pending ceiling."""


@dataclass
class ControllerTickResult:
    campaign_id: str
    policy_hash: str
    start_cursor: int
    end_cursor: int
    events_processed: int = 0
    tasks_created: int = 0
    decisions_created: int = 0
    events_ignored: int = 0
    tasks_cancelled: int = 0
    backpressured_profile: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "campaign_id": self.campaign_id,
            "policy_hash": self.policy_hash,
            "start_cursor": self.start_cursor,
            "end_cursor": self.end_cursor,
            "events_processed": self.events_processed,
            "tasks_created": self.tasks_created,
            "decisions_created": self.decisions_created,
            "events_ignored": self.events_ignored,
            "tasks_cancelled": self.tasks_cancelled,
            "backpressured_profile": self.backpressured_profile,
        }


def deterministic_sample(
    rate: float,
    *,
    seed: str,
    campaign_id: str,
    subject_id: str,
    stratum: dict[str, Any] | None = None,
) -> tuple[bool, float]:
    """Return a reproducible Bernoulli draw and its [0, 1) value."""

    if not 0.0 <= rate <= 1.0:
        raise ValueError("Sampling rate must be between 0 and 1")
    digest = hashlib.sha256(
        json.dumps(
            {
                "seed": seed,
                "campaign_id": campaign_id,
                "subject_id": subject_id,
                "stratum": stratum or {},
            },
            sort_keys=True,
            separators=(",", ":"),
            default=str,
        ).encode()
    ).digest()
    value = int.from_bytes(digest[:8], "big") / 2**64
    return value < rate, value


class ReviewController:
    """Translate durable Ops events into idempotent, capability-routed work."""

    def __init__(self, store: OpsStore, policy: ReviewPolicy):
        self.store = store
        self.policy = policy

    @staticmethod
    def _target(review: dict[str, Any]) -> tuple[str, str]:
        for kind, field in (
            ("finding", "finding_id"),
            ("candidate_cluster", "cluster_id"),
            ("unit_disposition", "unit_disposition_id"),
        ):
            value = review.get(field)
            if value is not None:
                return kind, str(value)
        raise ValueError(f"Review {review['id']} has no target")

    def _target_record(
        self,
        target_kind: str,
        target_id: str,
    ) -> dict[str, Any]:
        if target_kind == "candidate_cluster":
            return self.store.get_candidate_cluster(target_id)
        if target_kind == "finding":
            return self.store.get_finding(target_id)
        if target_kind == "unit_disposition":
            return self.store.get_unit_disposition(target_id)
        raise ValueError(f"Unsupported review target kind {target_kind!r}")

    @staticmethod
    def _dataset_id(
        target_kind: str,
        target: dict[str, Any],
        *,
        fallback: str | None = None,
    ) -> str:
        value = target.get("dataset_id") or fallback
        if not value:
            raise ValueError(
                f"Review target {target_kind} lacks an exact dataset_id"
            )
        return str(value)

    def _task_by_key(self, idempotency_key: str) -> str | None:
        with self.store._lock:
            row = self.store._conn.execute(
                "SELECT id FROM tasks WHERE idempotency_key = ?",
                (idempotency_key,),
            ).fetchone()
        return str(row["id"]) if row is not None else None

    def _pending_for_profile(self, profile_name: str) -> int:
        with self.store._lock:
            return int(
                self.store._conn.execute(
                    """
                    SELECT COUNT(*) FROM tasks
                    WHERE status IN ('pending', 'claimed', 'in_progress')
                      AND execution_profile = ?
                    """,
                    (profile_name,),
                ).fetchone()[0]
            )

    def _cancel_superseded_cluster_work(self, cluster_id: str) -> int:
        """Fence active work whose immutable target became stale."""

        with self.store._lock, self.store._transaction():
            cluster = self.store._conn.execute(
                "SELECT campaign_id FROM candidate_clusters WHERE id = ?",
                (cluster_id,),
            ).fetchone()
            if cluster is None:
                return 0
            rows = self.store._conn.execute(
                """
                SELECT id, run_id
                FROM tasks
                WHERE campaign_id = ?
                  AND status IN (
                      'pending', 'claimed', 'in_progress', 'retry_wait'
                  )
                  AND (
                      (
                          json_extract(params, '$.target.kind')
                              = 'candidate_cluster'
                          AND json_extract(params, '$.target.id') = ?
                      )
                      OR (
                          json_extract(params, '$.target.kind') = 'finding'
                          AND EXISTS (
                              SELECT 1
                              FROM cluster_findings AS cf
                              WHERE cf.cluster_id = ?
                                AND cf.relation = 'materializes'
                                AND cf.finding_id = json_extract(
                                    tasks.params, '$.target.id'
                                )
                          )
                      )
                  )
                ORDER BY id
                """,
                (cluster["campaign_id"], cluster_id, cluster_id),
            ).fetchall()
            now = time.time()
            for row in rows:
                self.store._conn.execute(
                    """
                    UPDATE tasks
                    SET status = 'cancelled', completed_ts = ?,
                        heartbeat_ts = ?, lease_expires_ts = NULL,
                        lease_token_hash = NULL,
                        last_error = 'target candidate cluster superseded'
                    WHERE id = ?
                      AND status IN (
                          'pending', 'claimed', 'in_progress', 'retry_wait'
                      )
                    """,
                    (now, now, row["id"]),
                )
                self.store._event_locked(
                    "task_cancelled",
                    "task",
                    str(row["id"]),
                    campaign_id=str(cluster["campaign_id"]),
                    run_id=row["run_id"],
                    task_id=str(row["id"]),
                    payload={
                        "reason": "target_candidate_cluster_superseded",
                        "cluster_id": cluster_id,
                    },
                    actor_agent_id=self.store.agent_id,
                )
            return len(rows)

    def _finding_materialization_is_current(self, finding_id: str) -> bool:
        with self.store._lock:
            rows = self.store._conn.execute(
                """
                SELECT c.status
                FROM cluster_findings AS cf
                JOIN candidate_clusters AS c ON c.id = cf.cluster_id
                WHERE cf.finding_id = ? AND cf.relation = 'materializes'
                """,
                (finding_id,),
            ).fetchall()
        return bool(rows) and all(row["status"] == "active" for row in rows)

    def _route_materialized_finding(
        self,
        finding_id: str,
        result: ControllerTickResult,
    ) -> None:
        finding = self.store.get_finding(finding_id)
        task_id = finding.get("task_id")
        if task_id is not None:
            task = self.store.get_task(str(task_id))
            if (
                task["params"].get("materialization_contract") is not None
                and task["status"] != "complete"
            ):
                return
        links = [
            link
            for link in self.store.list_cluster_findings(finding_id=finding_id)
            if link["relation"] == "materializes"
        ]
        current_clusters = [
            self.store.get_candidate_cluster(str(link["cluster_id"]))
            for link in links
        ]
        if not current_clusters or any(
            cluster["status"] != "active" for cluster in current_clusters
        ):
            return
        datasets = {str(cluster["dataset_id"]) for cluster in current_clusters}
        campaigns = {str(cluster["campaign_id"]) for cluster in current_clusters}
        if len(datasets) != 1 or len(campaigns) != 1:
            raise ValueError(
                "Materialized finding crosses a campaign or dataset boundary"
            )
        dataset_id = next(iter(datasets))
        campaign_id = next(iter(campaigns))
        tier_name = self.policy.routing.canonical_tier
        for profile_name in self.policy.tier(tier_name).profiles:
            _, created = self._create_review_task(
                campaign_id=campaign_id,
                target_kind="finding",
                target_id=finding_id,
                dataset_id=dataset_id,
                tier_name=tier_name,
                profile_name=profile_name,
                priority=3,
            )
            result.tasks_created += int(created)

    def _atlas_campaign_may_be_ready(self, campaign_id: str) -> bool:
        campaign = self.store.get_campaign(campaign_id)
        metadata = dict(campaign.get("metadata") or {})
        if metadata.get("unit_count") is None:
            return False
        with self.store._lock:
            unfinished = self.store._conn.execute(
                """
                SELECT 1
                FROM tasks
                WHERE campaign_id = ?
                  AND task_type = 'atlas_genome_read'
                  AND status IN (
                      'pending', 'claimed', 'in_progress', 'retry_wait',
                      'failed', 'cancelled'
                  )
                LIMIT 1
                """,
                (campaign_id,),
            ).fetchone()
        return unfinished is None

    def _revisit_canonical_reviews(
        self,
        campaign_id: str,
        result: ControllerTickResult,
    ) -> None:
        if not self.store.review_campaign_scan_status(campaign_id)["ready"]:
            return
        with self.store._lock:
            rows = self.store._conn.execute(
                """
                SELECT id
                FROM finding_reviews
                WHERE campaign_id = ? AND review_tier = ?
                ORDER BY ts, id
                """,
                (campaign_id, self.policy.routing.canonical_tier),
            ).fetchall()
        for row in rows:
            self._process_review(str(row["id"]), result)

    def _cancel_superseded_disposition_work(self, disposition_id: str) -> int:
        """Fence audits targeting a historical unit-disposition version."""

        with self.store._lock, self.store._transaction():
            disposition = self.store._conn.execute(
                "SELECT campaign_id FROM unit_dispositions WHERE id = ?",
                (disposition_id,),
            ).fetchone()
            if disposition is None:
                return 0
            rows = self.store._conn.execute(
                """
                SELECT id, run_id
                FROM tasks
                WHERE campaign_id = ?
                  AND status IN (
                      'pending', 'claimed', 'in_progress', 'retry_wait'
                  )
                  AND json_extract(params, '$.target.kind')
                      = 'unit_disposition'
                  AND json_extract(params, '$.target.id') = ?
                ORDER BY id
                """,
                (disposition["campaign_id"], disposition_id),
            ).fetchall()
            now = time.time()
            for row in rows:
                self.store._conn.execute(
                    """
                    UPDATE tasks
                    SET status = 'cancelled', completed_ts = ?,
                        heartbeat_ts = ?, lease_expires_ts = NULL,
                        lease_token_hash = NULL,
                        last_error = 'target unit disposition superseded'
                    WHERE id = ?
                      AND status IN (
                          'pending', 'claimed', 'in_progress', 'retry_wait'
                      )
                    """,
                    (now, now, row["id"]),
                )
                self.store._event_locked(
                    "task_cancelled",
                    "task",
                    str(row["id"]),
                    campaign_id=str(disposition["campaign_id"]),
                    run_id=row["run_id"],
                    task_id=str(row["id"]),
                    payload={
                        "reason": "target_unit_disposition_superseded",
                        "unit_disposition_id": disposition_id,
                    },
                    actor_agent_id=self.store.agent_id,
                )
            return len(rows)

    def _create_review_task(
        self,
        *,
        campaign_id: str,
        target_kind: str,
        target_id: str,
        dataset_id: str,
        tier_name: str,
        profile_name: str,
        source_review_ids: list[str] | None = None,
        audit: dict[str, Any] | None = None,
        priority: int = 2,
    ) -> tuple[str, bool]:
        tier = self.policy.tier(tier_name) if tier_name in self.policy.tiers else None
        profile = self.policy.profile(profile_name)
        sources = sorted(set(source_review_ids or []))
        source_hash = content_hash(sources)
        key = (
            f"review:{self.policy.policy_hash}:{campaign_id}:{target_kind}:"
            f"{target_id}:{tier_name}:{profile_name}:{source_hash}"
        )
        existing = self._task_by_key(key)
        if existing is not None:
            return existing, False
        if self._pending_for_profile(profile_name) >= profile.max_pending_tasks:
            raise ControllerBackpressureError(profile_name)
        target = self._target_record(target_kind, target_id)
        immutable_input = {
            key: target.get(key)
            for key in (
                "id",
                "version",
                "member_manifest_hash",
                "coverage_hash",
                "evidence_bundle_hash",
                "created_event_id",
            )
            if target.get(key) is not None
        }
        params = {
            "review_contract": REVIEW_TASK_CONTRACT,
            "campaign_id": campaign_id,
            "dataset_id": dataset_id,
            "target": {"kind": target_kind, "id": target_id},
            "target_input": immutable_input,
            "review_tier": tier_name,
            "execution_profile": profile_name,
            "resolved_execution": {
                "provider": profile.provider,
                "model": profile.model,
                "reasoning_effort": profile.reasoning_effort,
            },
            "policy": {
                "name": self.policy.name,
                "version": self.policy.version,
                "hash": self.policy.policy_hash,
            },
            "rubric_version": self.policy.rubric_version,
            "blind_to_prior_scores": (
                tier.blind_to_prior_scores if tier is not None else True
            ),
            "blind_to_other_reviews": (
                tier.blind_to_other_reviews if tier is not None else True
            ),
            "source_review_ids": sources,
            "source_review_manifest_hash": source_hash,
            "audit": audit,
            "scientific_contract": {
                "reconstruct_observations": True,
                "discover_live_curated_callers_before_named_claims": True,
                "executable_verification_records": True,
                "raw_sequences_model_visible": False,
            },
        }
        task_id = self.store.create_task(
            "scientific_review",
            (
                f"Review {target_kind} {target_id} at tier {tier_name} under "
                f"profile {profile_name}"
            ),
            params=params,
            priority=priority,
            campaign_id=campaign_id,
            idempotency_key=key,
            required_capabilities=[profile.capability],
            resource_request=profile.resource_request,
            max_attempts=3,
            lease_seconds=1_800,
        )
        return task_id, True

    def _create_materialization_task(
        self,
        *,
        campaign_id: str,
        cluster_id: str,
        dataset_id: str,
        source_review_ids: list[str],
    ) -> tuple[str, bool]:
        profile_name = self.policy.routing.materialization_profile
        profile = self.policy.profile(profile_name)
        source_ids = sorted(set(source_review_ids))
        source_hash = content_hash(source_ids)
        key = (
            f"materialize:{self.policy.policy_hash}:{campaign_id}:"
            f"{cluster_id}:{source_hash}"
        )
        existing = self._task_by_key(key)
        if existing is not None:
            return existing, False
        if self._pending_for_profile(profile_name) >= profile.max_pending_tasks:
            raise ControllerBackpressureError(profile_name)
        cluster = self.store.get_candidate_cluster(cluster_id)
        params = {
            "materialization_contract": MATERIALIZATION_TASK_CONTRACT,
            "campaign_id": campaign_id,
            "dataset_id": dataset_id,
            "target": {"kind": "candidate_cluster", "id": cluster_id},
            "target_input": {
                "version": cluster["version"],
                "member_manifest_hash": cluster["member_manifest_hash"],
            },
            "source_review_ids": source_ids,
            "source_review_manifest_hash": source_hash,
            "execution_profile": profile_name,
            "resolved_execution": {
                "provider": profile.provider,
                "model": profile.model,
                "reasoning_effort": profile.reasoning_effort,
            },
            "policy": {
                "name": self.policy.name,
                "version": self.policy.version,
                "hash": self.policy.policy_hash,
            },
            "output_contract": {
                "create_finding": True,
                "link_cluster_relation": "materializes",
                "validation_status": "review_dag_pending",
                "canonical_publication": False,
            },
        }
        task_id = self.store.create_task(
            "materialize_finding",
            f"Materialize reviewed candidate cluster {cluster_id} as a finding",
            params=params,
            priority=3,
            campaign_id=campaign_id,
            idempotency_key=key,
            required_capabilities=[profile.capability],
            resource_request=profile.resource_request,
            max_attempts=3,
            lease_seconds=1_800,
        )
        return task_id, True

    def _create_evidence_tasks(
        self,
        review: dict[str, Any],
    ) -> tuple[list[str], int]:
        target_kind, target_id = self._target(review)
        proposals = list(review.get("proposed_tasks") or [])
        if not proposals:
            proposals = [
                {
                    "task_type": "collect_review_evidence",
                    "description": (
                        f"Collect evidence requested by review {review['id']}"
                    ),
                    "params": {},
                    "required_capabilities": ["review:evidence"],
                }
            ]
        proposals = proposals[
            : self.policy.limits.max_evidence_tasks_per_review
        ]
        task_ids: list[str] = []
        created = 0
        for index, proposal in enumerate(proposals):
            if not isinstance(proposal, dict):
                raise ValueError("proposed_tasks entries must be objects")
            task_type = str(
                proposal.get("task_type") or "collect_review_evidence"
            ).strip()
            description = str(
                proposal.get("description")
                or f"Collect evidence requested by review {review['id']}"
            ).strip()
            if not task_type or not description:
                raise ValueError("Proposed evidence task requires type and description")
            params = dict(proposal.get("params") or {})
            params.update(
                {
                    "source_review_id": review["id"],
                    "target": {"kind": target_kind, "id": target_id},
                    "dataset_id": review["dataset_id"],
                    "policy_hash": self.policy.policy_hash,
                }
            )
            key = (
                f"review-evidence:{self.policy.policy_hash}:{review['id']}:"
                f"{index}:{content_hash(proposal)}"
            )
            existing = self._task_by_key(key)
            if existing is not None:
                task_ids.append(existing)
                continue
            task_id = self.store.create_task(
                task_type,
                description[:262_144],
                params=params,
                priority=min(3, max(0, int(proposal.get("priority", 2)))),
                campaign_id=review["campaign_id"],
                idempotency_key=key,
                required_capabilities=[
                    str(value)
                    for value in proposal.get(
                        "required_capabilities", ["review:evidence"]
                    )
                ],
                resource_request=dict(proposal.get("resource_request") or {}),
                max_attempts=3,
                lease_seconds=1_800,
            )
            task_ids.append(task_id)
            created += 1
        return task_ids, created

    def _create_decision(
        self,
        *,
        review: dict[str, Any],
        decision: str,
        source_tier: str,
        target_tier: str | None,
        review_ids: list[str],
        task_ids: list[str],
        rationale: str,
        audit_sample: bool = False,
        audit_stratum: dict[str, Any] | None = None,
    ) -> tuple[str, bool]:
        target_kind, target_id = self._target(review)
        if target_kind == "unit_disposition":
            raise ValueError("Unit audit outcomes remain review records")
        review_ids = sorted(set(review_ids))
        task_ids = sorted(set(task_ids))
        key = (
            f"decision:{self.policy.policy_hash}:{target_kind}:{target_id}:"
            f"{source_tier}:{target_tier}:{decision}:"
            f"{content_hash({'reviews': review_ids, 'tasks': task_ids})}"
        )
        with self.store._lock:
            existing = self.store._conn.execute(
                """
                SELECT id FROM promotion_decisions
                WHERE actor_agent_id = ? AND idempotency_key = ?
                """,
                (self.store.agent_id, key),
            ).fetchone()
        if existing is not None:
            return str(existing["id"]), False
        target_argument = (
            {"finding_id": target_id}
            if target_kind == "finding"
            else {"cluster_id": target_id}
        )
        decision_id = self.store.create_promotion_decision(
            campaign_id=review["campaign_id"],
            decision=decision,
            source_tier=source_tier,
            target_tier=target_tier,
            policy_name=self.policy.name,
            policy_version=self.policy.version,
            policy_hash=self.policy.policy_hash,
            rationale=rationale,
            review_ids=review_ids,
            created_task_ids=task_ids,
            audit_sample=audit_sample,
            audit_stratum=audit_stratum,
            idempotency_key=key,
            **target_argument,
        )
        return decision_id, True

    def _verification_state(self, review_id: str) -> str:
        with self.store._lock:
            rows = self.store._conn.execute(
                """
                WITH latest AS (
                    SELECT
                        status,
                        ROW_NUMBER() OVER (
                            PARTITION BY claim_key, specification_hash
                            ORDER BY attempt DESC, created_ts DESC, id DESC
                        ) AS rank
                    FROM review_verifications
                    WHERE review_id = ?
                )
                SELECT status FROM latest WHERE rank = 1
                """,
                (review_id,),
            ).fetchall()
        if not rows:
            return "missing"
        statuses = {str(row["status"]) for row in rows}
        if statuses == {"pass"}:
            return "passed"
        if statuses & {"fail", "error"}:
            return "failed"
        return "pending"

    def _selected_reviews(
        self,
        review: dict[str, Any],
        tier_name: str,
    ) -> list[dict[str, Any]]:
        target_kind, target_id = self._target(review)
        filters: dict[str, Any] = {
            "campaign_id": review["campaign_id"],
            "review_tier": tier_name,
            f"{'cluster' if target_kind == 'candidate_cluster' else target_kind}_id": (
                target_id
            ),
        }
        rows = self.store.list_finding_reviews(**filters)
        latest_by_profile: dict[str, dict[str, Any]] = {}
        for row in rows:
            profile = str(row["execution_profile"])
            current = latest_by_profile.get(profile)
            if current is None or (row["ts"], row["id"]) > (
                current["ts"],
                current["id"],
            ):
                latest_by_profile[profile] = row
        tier = self.policy.tier(tier_name)
        return [
            latest_by_profile[profile]
            for profile in tier.profiles
            if profile in latest_by_profile
        ]

    def _route_initial_review(
        self,
        review: dict[str, Any],
        result: ControllerTickResult,
    ) -> None:
        tier_name = self.policy.routing.initial_tier
        if review["review_tier"] != tier_name:
            return
        tier = self.policy.tier(tier_name)
        verdict = str(review["verdict"])
        verification_state = self._verification_state(str(review["id"]))
        target_kind, target_id = self._target(review)
        if target_kind != "candidate_cluster":
            return
        if (
            verdict == "promote"
            and float(review["confidence"]) >= tier.minimum_confidence
            and verification_state == "passed"
        ):
            task_ids: list[str] = []
            for profile_name in self.policy.tier(
                self.policy.routing.independent_tier
            ).profiles:
                task_id, created = self._create_review_task(
                    campaign_id=review["campaign_id"],
                    target_kind=target_kind,
                    target_id=target_id,
                    dataset_id=review["dataset_id"],
                    tier_name=self.policy.routing.independent_tier,
                    profile_name=profile_name,
                    priority=2,
                )
                task_ids.append(task_id)
                result.tasks_created += int(created)
            _, created = self._create_decision(
                review=review,
                decision="promote",
                source_tier=tier_name,
                target_tier=self.policy.routing.independent_tier,
                review_ids=[review["id"]],
                task_ids=task_ids,
                rationale="Initial review met confidence and verification gates",
            )
            result.decisions_created += int(created)
            return
        if verdict == "needs_data":
            tasks, created_tasks = self._create_evidence_tasks(review)
            result.tasks_created += created_tasks
            _, created = self._create_decision(
                review=review,
                decision="needs_data",
                source_tier=tier_name,
                target_tier=tier_name,
                review_ids=[review["id"]],
                task_ids=tasks,
                rationale="Reviewer requested additional executable evidence",
            )
            result.decisions_created += int(created)
            return
        if verdict == "promote":
            if verification_state in {"missing", "pending"}:
                return
            decision = "hold"
            rationale = (
                "Promotion held because confidence or executable verification "
                "did not meet policy"
            )
        else:
            decision = verdict
            rationale = f"Initial review returned {verdict}"
        audit_task_ids: list[str] = []
        sampled = False
        stratum = {"candidate_type": self.store.get_candidate_cluster(target_id)[
            "candidate_type"
        ]}
        if decision in {"reject", "duplicate"}:
            sampled, sampling_value = deterministic_sample(
                self.policy.audit.rejected_cluster_rate,
                seed=self.policy.audit.sampling_seed,
                campaign_id=review["campaign_id"],
                subject_id=target_id,
                stratum=stratum,
            )
            stratum.update(
                {
                    "sampling_value": sampling_value,
                    "sampling_rate": self.policy.audit.rejected_cluster_rate,
                }
            )
            if sampled:
                task_id, created = self._create_review_task(
                    campaign_id=review["campaign_id"],
                    target_kind=target_kind,
                    target_id=target_id,
                    dataset_id=review["dataset_id"],
                    tier_name="audit_reject",
                    profile_name=self.policy.audit.rejected_profile,
                    audit=stratum,
                    priority=1,
                )
                audit_task_ids.append(task_id)
                result.tasks_created += int(created)
        _, created = self._create_decision(
            review=review,
            decision=decision,
            source_tier=tier_name,
            target_tier=None,
            review_ids=[review["id"]],
            task_ids=audit_task_ids,
            rationale=rationale,
            audit_sample=sampled,
            audit_stratum=stratum if sampled else {},
        )
        result.decisions_created += int(created)

    def _route_independent_reviews(
        self,
        review: dict[str, Any],
        result: ControllerTickResult,
    ) -> None:
        tier_name = self.policy.routing.independent_tier
        if review["review_tier"] != tier_name:
            return
        selected = self._selected_reviews(review, tier_name)
        tier = self.policy.tier(tier_name)
        if len(selected) < tier.minimum_reviews:
            return
        providers = {str(row["provider"]) for row in selected}
        if tier.require_distinct_providers and len(providers) < tier.minimum_reviews:
            return
        reviewers = {str(row["reviewer_agent_id"]) for row in selected}
        if tier.require_distinct_reviewers and len(reviewers) < tier.minimum_reviews:
            return
        if any(
            self._verification_state(str(row["id"])) in {"missing", "pending"}
            for row in selected
        ):
            return
        promoted = [
            row
            for row in selected
            if row["verdict"] == "promote"
            and float(row["confidence"]) >= tier.minimum_confidence
            and self._verification_state(str(row["id"])) == "passed"
        ]
        target_kind, target_id = self._target(review)
        adjudication_tier = self.policy.routing.adjudication_tier
        profile_name = self.policy.tier(adjudication_tier).profiles[0]
        source_ids = [str(row["id"]) for row in selected]
        task_id, created_task = self._create_review_task(
            campaign_id=review["campaign_id"],
            target_kind=target_kind,
            target_id=target_id,
            dataset_id=review["dataset_id"],
            tier_name=adjudication_tier,
            profile_name=profile_name,
            source_review_ids=source_ids,
            priority=3,
        )
        result.tasks_created += int(created_task)
        quorum = len(promoted) >= tier.minimum_promote_reviews
        _, created = self._create_decision(
            review=review,
            decision="promote" if quorum else "hold",
            source_tier=tier_name,
            target_tier=adjudication_tier,
            review_ids=source_ids,
            task_ids=[task_id],
            rationale=(
                "Independent cross-provider promotion quorum met"
                if quorum
                else "Independent reviews require adjudication"
            ),
        )
        result.decisions_created += int(created)

    def _route_adjudication(
        self,
        review: dict[str, Any],
        result: ControllerTickResult,
    ) -> None:
        tier_name = self.policy.routing.adjudication_tier
        if review["review_tier"] != tier_name:
            return
        target_kind, target_id = self._target(review)
        if target_kind != "candidate_cluster":
            return
        tier = self.policy.tier(tier_name)
        verdict = str(review["verdict"])
        verification_state = self._verification_state(str(review["id"]))
        if (
            verdict == "promote"
            and float(review["confidence"]) >= tier.minimum_confidence
            and verification_state == "passed"
        ):
            task_id, created_task = self._create_materialization_task(
                campaign_id=review["campaign_id"],
                cluster_id=target_id,
                dataset_id=review["dataset_id"],
                source_review_ids=[review["id"]],
            )
            result.tasks_created += int(created_task)
            _, created = self._create_decision(
                review=review,
                decision="promote",
                source_tier=tier_name,
                target_tier=self.policy.routing.canonical_tier,
                review_ids=[review["id"]],
                task_ids=[task_id],
                rationale="Adjudication met materialization gates",
            )
            result.decisions_created += int(created)
        elif verdict == "needs_data":
            tasks, created_tasks = self._create_evidence_tasks(review)
            result.tasks_created += created_tasks
            _, created = self._create_decision(
                review=review,
                decision="needs_data",
                source_tier=tier_name,
                target_tier=tier_name,
                review_ids=[review["id"]],
                task_ids=tasks,
                rationale="Adjudication requested additional evidence",
            )
            result.decisions_created += int(created)
        elif verdict == "promote" and verification_state in {"missing", "pending"}:
            return
        else:
            _, created = self._create_decision(
                review=review,
                decision="hold" if verdict == "promote" else verdict,
                source_tier=tier_name,
                target_tier=None,
                review_ids=[review["id"]],
                task_ids=[],
                rationale="Adjudication did not meet materialization gates",
            )
            result.decisions_created += int(created)

    def _route_canonical_review(
        self,
        review: dict[str, Any],
        result: ControllerTickResult,
    ) -> None:
        tier_name = self.policy.routing.canonical_tier
        if review["review_tier"] != tier_name:
            return
        target_kind, _ = self._target(review)
        if target_kind != "finding":
            return
        tier = self.policy.tier(tier_name)
        verdict = str(review["verdict"])
        verification_state = self._verification_state(str(review["id"]))
        if (
            verdict == "promote"
            and float(review["confidence"]) >= tier.minimum_confidence
            and verification_state == "passed"
        ):
            if not self.store.review_campaign_scan_status(
                str(review["campaign_id"])
            )["ready"]:
                return
            if not self._finding_materialization_is_current(
                str(review["finding_id"])
            ):
                _, created = self._create_decision(
                    review=review,
                    decision="hold",
                    source_tier=tier_name,
                    target_tier=None,
                    review_ids=[review["id"]],
                    task_ids=[],
                    rationale=(
                        "Canonical review source cluster changed after "
                        "finding materialization"
                    ),
                )
                result.decisions_created += int(created)
                return
            _, created = self._create_decision(
                review=review,
                decision="publish",
                source_tier=tier_name,
                target_tier="canonical_publication",
                review_ids=[review["id"]],
                task_ids=[],
                rationale=(
                    "Canonical review met confidence and executable verification gates"
                ),
            )
            result.decisions_created += int(created)
        elif verdict == "needs_data":
            tasks, created_tasks = self._create_evidence_tasks(review)
            result.tasks_created += created_tasks
            _, created = self._create_decision(
                review=review,
                decision="needs_data",
                source_tier=tier_name,
                target_tier=tier_name,
                review_ids=[review["id"]],
                task_ids=tasks,
                rationale="Canonical review requested additional evidence",
            )
            result.decisions_created += int(created)
        elif verdict == "promote" and verification_state in {"missing", "pending"}:
            return
        else:
            _, created = self._create_decision(
                review=review,
                decision="hold" if verdict == "promote" else verdict,
                source_tier=tier_name,
                target_tier=None,
                review_ids=[review["id"]],
                task_ids=[],
                rationale="Canonical publication gates were not met",
            )
            result.decisions_created += int(created)

    def _route_audit_review(
        self,
        review: dict[str, Any],
        result: ControllerTickResult,
    ) -> None:
        tier_name = str(review["review_tier"])
        if tier_name == "audit_clear":
            if review["verdict"] == "needs_data":
                _, created = self._create_evidence_tasks(review)
                result.tasks_created += created
            return
        if tier_name != "audit_reject":
            return
        target_kind, target_id = self._target(review)
        if target_kind != "candidate_cluster":
            return
        if (
            review["verdict"] == "promote"
            and self._verification_state(str(review["id"])) == "passed"
        ):
            task_ids: list[str] = []
            for profile_name in self.policy.tier(
                self.policy.routing.independent_tier
            ).profiles:
                task_id, created = self._create_review_task(
                    campaign_id=review["campaign_id"],
                    target_kind=target_kind,
                    target_id=target_id,
                    dataset_id=review["dataset_id"],
                    tier_name=self.policy.routing.independent_tier,
                    profile_name=profile_name,
                    source_review_ids=[review["id"]],
                    priority=3,
                )
                task_ids.append(task_id)
                result.tasks_created += int(created)
            _, created = self._create_decision(
                review=review,
                decision="reopen",
                source_tier=tier_name,
                target_tier=self.policy.routing.independent_tier,
                review_ids=[review["id"]],
                task_ids=task_ids,
                rationale="Rejected-case audit supplied verified contrary evidence",
                audit_sample=True,
            )
            result.decisions_created += int(created)

    def _process_review(
        self,
        review_id: str,
        result: ControllerTickResult,
    ) -> None:
        review = self.store.get_finding_review(review_id)
        if review.get("task_id") is not None:
            task = self.store.get_task(str(review["task_id"]))
            if task["status"] != "complete":
                return
        if review.get("cluster_id") is not None:
            cluster = self.store.get_candidate_cluster(str(review["cluster_id"]))
            if cluster["status"] != "active":
                return
        self._route_initial_review(review, result)
        self._route_independent_reviews(review, result)
        self._route_adjudication(review, result)
        self._route_canonical_review(review, result)
        self._route_audit_review(review, result)

    def _process_event(
        self,
        event: dict[str, Any],
        result: ControllerTickResult,
    ) -> bool:
        event_type = str(event["event_type"])
        if event_type == "candidate_cluster_created":
            cluster = self.store.get_candidate_cluster(str(event["entity_id"]))
            if cluster["status"] != "active":
                return True
            tier_name = self.policy.routing.initial_tier
            for profile_name in self.policy.tier(tier_name).profiles:
                _, created = self._create_review_task(
                    campaign_id=cluster["campaign_id"],
                    target_kind="candidate_cluster",
                    target_id=cluster["id"],
                    dataset_id=cluster["dataset_id"],
                    tier_name=tier_name,
                    profile_name=profile_name,
                    priority=2,
                )
                result.tasks_created += int(created)
            return True
        if event_type == "candidate_cluster_superseded":
            payload = dict(event.get("payload") or {})
            parent_cluster_id = payload.get("parent_cluster_id")
            if parent_cluster_id is not None:
                result.tasks_cancelled += self._cancel_superseded_cluster_work(
                    str(parent_cluster_id)
                )
            return True
        if event_type == "unit_disposition_created":
            disposition = self.store.get_unit_disposition(str(event["entity_id"]))
            superseded_id = disposition.get("supersedes_disposition_id")
            if superseded_id is not None:
                result.tasks_cancelled += (
                    self._cancel_superseded_disposition_work(str(superseded_id))
                )
            if self._atlas_campaign_may_be_ready(
                str(disposition["campaign_id"])
            ):
                self._revisit_canonical_reviews(
                    str(disposition["campaign_id"]),
                    result,
                )
            if (
                disposition["record_status"] != "active"
                or disposition["disposition"] != "clear"
            ):
                return True
            strata = dict(disposition.get("strata") or {})
            selected_strata = {
                key: strata.get(key)
                for key in self.policy.audit.stratification_fields
                if strata.get(key) is not None
            }
            sampled, value = deterministic_sample(
                self.policy.audit.clear_unit_rate,
                seed=self.policy.audit.sampling_seed,
                campaign_id=disposition["campaign_id"],
                subject_id=disposition["id"],
                stratum=selected_strata,
            )
            if sampled:
                audit = {
                    "sampling_rate": self.policy.audit.clear_unit_rate,
                    "sampling_value": value,
                    "stratum": selected_strata,
                    "audit_class": "clear_unit",
                }
                _, created = self._create_review_task(
                    campaign_id=disposition["campaign_id"],
                    target_kind="unit_disposition",
                    target_id=disposition["id"],
                    dataset_id=disposition["dataset_id"],
                    tier_name="audit_clear",
                    profile_name=self.policy.audit.clear_profile,
                    audit=audit,
                    priority=0,
                )
                result.tasks_created += int(created)
            return True
        if event_type == "cluster_finding_linked":
            finding_id = str(event["entity_id"])
            payload = dict(event.get("payload") or {})
            if payload.get("relation") != "materializes":
                return True
            self._route_materialized_finding(finding_id, result)
            return True
        if event_type == "finding_review_created":
            self._process_review(str(event["entity_id"]), result)
            return True
        if event_type == "review_verification_created":
            with self.store._lock:
                row = self.store._conn.execute(
                    "SELECT review_id FROM review_verifications WHERE id = ?",
                    (event["entity_id"],),
                ).fetchone()
            if row is not None:
                self._process_review(str(row["review_id"]), result)
            return True
        if event_type == "task_completed":
            task_id = str(event["entity_id"])
            task = self.store.get_task(task_id)
            handled = False
            if task["params"].get("review_contract") is not None:
                with self.store._lock:
                    row = self.store._conn.execute(
                        "SELECT id FROM finding_reviews WHERE task_id = ?",
                        (task_id,),
                    ).fetchone()
                if row is not None:
                    self._process_review(str(row["id"]), result)
                handled = True
            if task["params"].get("materialization_contract") is not None:
                with self.store._lock:
                    rows = self.store._conn.execute(
                        "SELECT id FROM findings WHERE task_id = ? ORDER BY id",
                        (task_id,),
                    ).fetchall()
                for row in rows:
                    self._route_materialized_finding(str(row["id"]), result)
                handled = True
            if (
                task["task_type"] == "atlas_genome_read"
                and task.get("campaign_id") is not None
                and self._atlas_campaign_may_be_ready(str(task["campaign_id"]))
            ):
                self._revisit_canonical_reviews(
                    str(task["campaign_id"]),
                    result,
                )
                handled = True
            return handled
        return False

    def tick(self, campaign_id: str) -> ControllerTickResult:
        """Consume one bounded event batch and persist its durable cursor."""

        self.store.get_campaign(campaign_id)
        state = self.store.get_review_controller_state(
            self.policy.controller_id,
            campaign_id,
        )
        if state is not None and state["policy_hash"] != self.policy.policy_hash:
            raise ValueError(
                "Review policy hash changed for this controller and campaign; "
                "use a new controller_id for a versioned replay"
            )
        cursor = int(state["last_event_id"]) if state is not None else 0
        result = ControllerTickResult(
            campaign_id=campaign_id,
            policy_hash=self.policy.policy_hash,
            start_cursor=cursor,
            end_cursor=cursor,
        )
        events = self.store.events(
            after_id=cursor,
            campaign_id=campaign_id,
            limit=self.policy.limits.event_batch_size,
        )
        for event in events:
            try:
                handled = self._process_event(event, result)
            except ControllerBackpressureError as exc:
                result.backpressured_profile = str(exc)
                break
            result.events_processed += 1
            result.events_ignored += int(not handled)
            cursor = int(event["id"])
            result.end_cursor = cursor
        self.store.set_review_controller_state(
            self.policy.controller_id,
            campaign_id,
            policy_hash=self.policy.policy_hash,
            last_event_id=cursor,
        )
        return result


__all__ = [
    "MATERIALIZATION_TASK_CONTRACT",
    "REVIEW_TASK_CONTRACT",
    "ControllerBackpressureError",
    "ControllerTickResult",
    "ReviewController",
    "deterministic_sample",
]
