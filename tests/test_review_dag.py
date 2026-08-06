"""Scientific review DAG, reducer, verifier, and routing contracts."""

from __future__ import annotations

import sqlite3
import time
from collections import Counter

import duckdb
import pytest
from fastapi.testclient import TestClient

from sharur.ops.schema import OPS_SCHEMA_VERSION
from sharur.ops.server import create_app
from sharur.ops.store import OpsStore
from sharur.review.controller import ReviewController, deterministic_sample
from sharur.review.models import load_review_policy
from sharur.review.reducer import ExactCandidateReducer
from sharur.review.verification import run_duckdb_verification


TOKEN = "review-test-token"


def _campaign(ops: OpsStore) -> str:
    return ops.create_campaign("review fixture", idempotency_key="campaign")


def _candidate(
    ops: OpsStore,
    campaign_id: str,
    *,
    unit_id: str,
    signature: dict | None = None,
    schema: str = "architecture/v1",
    candidate_type: str = "architecture",
    role_hint: str | None = None,
    task_id: str | None = None,
) -> str:
    features = {"strata": {"phylum": "P"}}
    if role_hint is not None:
        features["role_hint"] = role_hint
    return ops.create_candidate_occurrence(
        campaign_id=campaign_id,
        dataset_id="dataset:sealed",
        unit_id=unit_id,
        genome_id=unit_id,
        candidate_type=candidate_type,
        signature_schema=schema,
        signature=signature or {"domains": ["A", "B"]},
        evidence={"protein_id": f"{unit_id}_p1"},
        verification=[
            {"claim": "one occurrence", "query": "SELECT 1", "expected": 1}
        ],
        subject_refs={"genome_id": unit_id, "protein_id": f"{unit_id}_p1"},
        task_id=task_id,
        reduction_features=features,
        idempotency_key=f"candidate:{unit_id}:{schema}:{candidate_type}",
    )


def _cluster_fixture(ops: OpsStore, campaign_id: str) -> str:
    _candidate(ops, campaign_id, unit_id="g1")
    return ops.create_candidate_cluster(
        campaign_id=campaign_id,
        dataset_id="dataset:sealed",
        candidate_type="architecture",
        signature_schema="architecture/v1",
        member_ids=[
            row["id"]
            for row in ops.list_candidate_occurrences(campaign_id=campaign_id)
        ],
        reducer_name="fixture",
        reducer_version="1",
        reducer_config_hash="config",
        summary={"signature": {"domains": ["A", "B"]}},
        counts={"occurrences": 1, "genomes": 1},
        idempotency_key="cluster",
    )


def _review(
    ops: OpsStore,
    campaign_id: str,
    *,
    finding_id: str | None = None,
    cluster_id: str | None = None,
    unit_disposition_id: str | None = None,
    tier: str = "canonical",
    profile: str = "canonical_review",
    provider: str = "anthropic",
    model: str = "fixture-model",
    verdict: str = "promote",
    task_id: str | None = None,
    key: str = "review",
) -> str:
    return ops.create_finding_review(
        campaign_id=campaign_id,
        dataset_id="dataset:sealed",
        review_tier=tier,
        execution_profile=profile,
        provider=provider,
        model=model,
        reasoning_effort="high",
        prompt_hash=f"prompt:{key}",
        rubric_version="scientific-review/1.0",
        input_bundle_hash=f"bundle:{key}",
        verdict=verdict,
        confidence=0.95,
        finding_id=finding_id,
        cluster_id=cluster_id,
        unit_disposition_id=unit_disposition_id,
        task_id=task_id,
        reconstructed_observations={"observed": ["A", "B"]},
        idempotency_key=key,
    )


def _pass_verification(
    ops: OpsStore,
    review_id: str,
    *,
    key: str,
    supersedes: str | None = None,
) -> str:
    return ops.record_review_verification(
        review_id=review_id,
        claim_key="fixture-count",
        engine="duckdb",
        specification={"sql": "SELECT COUNT(*) FROM proteins"},
        dataset_id="dataset:sealed",
        expected=1,
        actual=1,
        status="pass",
        executed_ts=time.time(),
        supersedes_verification_id=supersedes,
        idempotency_key=key,
    )


def test_default_policy_uses_sol_medium_for_deepening_and_clear_audits() -> None:
    policy = load_review_policy()
    for profile_name in ("finding_deepen", "audit_clear"):
        profile = policy.profile(profile_name)
        assert profile.model == "gpt-5.6-sol"
        assert profile.reasoning_effort == "medium"


def test_v4_database_migrates_additively_to_review_schema(tmp_path) -> None:
    path = tmp_path / "ops.db"
    ops = OpsStore(path, agent_id="seed")
    campaign_id = _campaign(ops)
    ops.close()

    connection = sqlite3.connect(path)
    connection.execute("PRAGMA foreign_keys=OFF")
    for table in (
        "review_controller_state",
        "canonical_publications",
        "promotion_decision_tasks",
        "promotion_decision_reviews",
        "promotion_decisions",
        "review_verifications",
        "finding_reviews",
        "cluster_findings",
        "candidate_cluster_lineage",
        "candidate_cluster_members",
        "candidate_clusters",
        "candidate_occurrences",
        "unit_dispositions",
    ):
        connection.execute(f"DROP TABLE {table}")
    connection.execute("DELETE FROM ops_schema_meta WHERE version = ?", (5,))
    connection.commit()
    connection.close()

    migrated = OpsStore(path, agent_id="reader")
    assert migrated.get_campaign(campaign_id)["name"] == "review fixture"
    assert (
        migrated._conn.execute(
            "SELECT MAX(version) FROM ops_schema_meta"
        ).fetchone()[0]
        == OPS_SCHEMA_VERSION
        == 5
    )
    assert migrated._conn.execute(
        "SELECT 1 FROM sqlite_master WHERE name = 'finding_reviews'"
    ).fetchone()
    assert migrated.integrity_check()["ok"] is True


def test_unit_dispositions_are_reconciled_versioned_and_idempotent(tmp_path) -> None:
    ops = OpsStore(tmp_path / "ops.db", agent_id="scanner")
    campaign_id = _campaign(ops)
    _candidate(ops, campaign_id, unit_id="g1")
    first = ops.record_unit_disposition(
        campaign_id=campaign_id,
        unit_id="g1",
        dataset_id="dataset:sealed",
        genome_id="g1",
        coverage_hash="coverage-1",
        candidate_count=1,
        disposition="candidate",
        evidence_bundle_hash="bundle-1",
        idempotency_key="disposition-1",
    )
    assert (
        ops.record_unit_disposition(
            campaign_id=campaign_id,
            unit_id="g1",
            dataset_id="dataset:sealed",
            genome_id="g1",
            coverage_hash="coverage-1",
            candidate_count=1,
            disposition="candidate",
            evidence_bundle_hash="bundle-1",
            idempotency_key="disposition-1",
        )
        == first
    )
    with pytest.raises(ValueError, match="active disposition"):
        ops.record_unit_disposition(
            campaign_id=campaign_id,
            unit_id="g1",
            dataset_id="dataset:sealed",
            genome_id="g1",
            coverage_hash="coverage-2",
            candidate_count=1,
            disposition="candidate",
            evidence_bundle_hash="bundle-2",
            idempotency_key="disposition-2",
        )
    ops.create_candidate_occurrence(
        campaign_id=campaign_id,
        dataset_id="dataset:sealed",
        unit_id="g1",
        genome_id="g1",
        candidate_type="architecture",
        signature_schema="architecture/v1",
        signature={"domains": ["C"]},
        evidence={"protein_id": "g1_p2"},
        verification=[],
        subject_refs={"protein_id": "g1_p2"},
        idempotency_key="candidate:g1:second",
    )
    second = ops.record_unit_disposition(
        campaign_id=campaign_id,
        unit_id="g1",
        dataset_id="dataset:sealed",
        genome_id="g1",
        coverage_hash="coverage-2",
        candidate_count=2,
        disposition="candidate",
        evidence_bundle_hash="bundle-2",
        supersedes_disposition_id=first,
        idempotency_key="disposition-2",
    )
    assert ops.get_unit_disposition(first)["record_status"] == "superseded"
    assert ops.get_unit_disposition(second)["version"] == 2
    assert [row["id"] for row in ops.list_unit_dispositions(campaign_id=campaign_id)] == [
        second
    ]


def test_disposition_supersession_cancels_stale_clear_audit(tmp_path) -> None:
    ops = OpsStore(tmp_path / "ops.db", agent_id="controller")
    campaign_id = _campaign(ops)
    clear_id = ops.record_unit_disposition(
        campaign_id=campaign_id,
        unit_id="g1",
        dataset_id="dataset:sealed",
        genome_id="g1",
        coverage_hash="clear-coverage",
        candidate_count=0,
        disposition="clear",
        evidence_bundle_hash="clear-bundle",
        idempotency_key="clear-disposition",
    )
    policy = load_review_policy()
    policy.audit.clear_unit_rate = 1.0
    controller = ReviewController(ops, policy)
    assert controller.tick(campaign_id).tasks_created == 1
    audit_task = next(
        task
        for task in ops.list_tasks(campaign_id=campaign_id, limit=100)
        if task["params"].get("review_tier") == "audit_clear"
    )
    _candidate(ops, campaign_id, unit_id="g1")
    ops.record_unit_disposition(
        campaign_id=campaign_id,
        unit_id="g1",
        dataset_id="dataset:sealed",
        genome_id="g1",
        coverage_hash="candidate-coverage",
        candidate_count=1,
        disposition="candidate",
        evidence_bundle_hash="candidate-bundle",
        supersedes_disposition_id=clear_id,
        idempotency_key="candidate-disposition",
    )
    routed = controller.tick(campaign_id)
    assert routed.tasks_cancelled == 1
    assert ops.get_task(audit_task["id"])["status"] == "cancelled"


def test_exact_reducer_preserves_occurrences_and_versions_clusters(tmp_path) -> None:
    ops = OpsStore(tmp_path / "ops.db", agent_id="reducer")
    campaign_id = _campaign(ops)
    _candidate(ops, campaign_id, unit_id="g1")
    _candidate(ops, campaign_id, unit_id="g2")
    _candidate(
        ops,
        campaign_id,
        unit_id="g3",
        signature={"domains": ["C"]},
    )
    _candidate(
        ops,
        campaign_id,
        unit_id="g4",
        schema="architecture/v2",
    )
    reducer = ExactCandidateReducer()
    first = reducer.reduce_campaign(ops, campaign_id)
    assert first.candidate_groups_seen == 3
    assert first.clusters_created == 3
    clusters = ops.list_candidate_clusters(campaign_id=campaign_id)
    two_member = next(
        cluster for cluster in clusters if cluster["counts"]["occurrences"] == 2
    )
    assert two_member["counts"]["genomes"] == 2
    assert two_member["counts"]["strata"]["phylum"] == {"P": 2}

    _candidate(
        ops,
        campaign_id,
        unit_id="g5",
        role_hint="counterexample",
    )
    second = reducer.reduce_campaign(ops, campaign_id)
    assert second.clusters_created == 1
    assert second.clusters_versioned == 1
    versions = [
        cluster
        for cluster in ops.list_candidate_clusters(
            campaign_id=campaign_id, status=None
        )
        if cluster["logical_cluster_id"] == two_member["logical_cluster_id"]
    ]
    assert [(row["version"], row["status"]) for row in versions] == [
        (1, "superseded"),
        (2, "active"),
    ]
    current = ops.get_candidate_cluster(versions[-1]["id"])
    assert current["member_count"] == 3
    assert current["members_truncated"] is False
    assert len(current["members"]) == 3
    assert Counter(member["role"] for member in current["members"]) == {
        "medoid": 1,
        "member": 1,
        "counterexample": 1,
    }
    bounded = ops.get_candidate_cluster(versions[-1]["id"], member_limit=2)
    assert bounded["member_count"] == 3
    assert len(bounded["members"]) == 2
    assert bounded["members_truncated"] is True
    first_page = ops.list_candidate_cluster_members(
        versions[-1]["id"],
        limit=2,
    )
    second_page = ops.list_candidate_cluster_members(
        versions[-1]["id"],
        after_candidate_id=first_page["next_after_candidate_id"],
        limit=2,
    )
    assert first_page["has_more"] is True
    assert second_page["has_more"] is False
    assert sorted(
        member["candidate_id"]
        for member in first_page["members"] + second_page["members"]
    ) == sorted(member["candidate_id"] for member in current["members"])
    assert reducer.reduce_campaign(ops, campaign_id).clusters_unchanged == 3


def test_cluster_supersession_cancels_stale_work_and_blocks_publication(
    tmp_path,
) -> None:
    ops = OpsStore(tmp_path / "ops.db", agent_id="controller")
    campaign_id = _campaign(ops)
    _candidate(ops, campaign_id, unit_id="g1")
    reducer = ExactCandidateReducer()
    reducer.reduce_campaign(ops, campaign_id)
    first_cluster = ops.list_candidate_clusters(campaign_id=campaign_id)[0]
    controller = ReviewController(ops, load_review_policy())
    assert controller.tick(campaign_id).tasks_created == 1
    stale_deep_task = next(
        task
        for task in ops.list_tasks(campaign_id=campaign_id, limit=100)
        if task["params"].get("target", {}).get("id") == first_cluster["id"]
    )
    finding_id = ops.finding(
        "reviewed_candidate",
        "test",
        "Materialized before cluster revision",
        campaign_id=campaign_id,
        idempotency_key="stale-finding",
    )
    ops.link_cluster_finding(first_cluster["id"], finding_id)
    assert controller.tick(campaign_id).tasks_created == 1
    _candidate(ops, campaign_id, unit_id="g2")
    reducer.reduce_campaign(ops, campaign_id)

    routed = controller.tick(campaign_id)
    assert routed.tasks_created == 1
    assert routed.tasks_cancelled == 2
    assert ops.get_task(stale_deep_task["id"])["status"] == "cancelled"
    canonical_tasks = [
        task
        for task in ops.list_tasks(campaign_id=campaign_id, limit=100)
        if task["params"].get("review_tier") == "canonical"
    ]
    assert len(canonical_tasks) == 1
    assert canonical_tasks[0]["status"] == "cancelled"

    review_id = _review(
        ops,
        campaign_id,
        finding_id=finding_id,
        key="stale-canonical-review",
    )
    _pass_verification(ops, review_id, key="stale-canonical-verification")
    with pytest.raises(ValueError, match="superseded candidate clusters"):
        ops.create_promotion_decision(
            campaign_id=campaign_id,
            decision="publish",
            source_tier="canonical",
            policy_name="fixture",
            policy_version="1",
            policy_hash="policy",
            rationale="stale source",
            finding_id=finding_id,
            review_ids=[review_id],
            idempotency_key="stale-publication-decision",
        )


def test_append_only_verification_revisions_gate_publication(tmp_path) -> None:
    ops = OpsStore(tmp_path / "ops.db", agent_id="reviewer")
    campaign_id = _campaign(ops)
    cluster_id = _cluster_fixture(ops, campaign_id)
    finding_id = ops.finding(
        "observation",
        "fixture",
        "Observed typed architecture",
        campaign_id=campaign_id,
        idempotency_key="finding",
    )
    ops.link_cluster_finding(cluster_id, finding_id)
    review_id = _review(ops, campaign_id, finding_id=finding_id)
    pending = ops.record_review_verification(
        review_id=review_id,
        claim_key="fixture-count",
        engine="duckdb",
        specification={"sql": "SELECT COUNT(*) FROM proteins"},
        dataset_id="dataset:sealed",
        expected=1,
        status="pending",
        idempotency_key="verification-pending",
    )
    with pytest.raises(ValueError, match="passed executable"):
        ops.create_promotion_decision(
            campaign_id=campaign_id,
            decision="publish",
            source_tier="canonical",
            policy_name="fixture",
            policy_version="1",
            policy_hash="policy",
            rationale="pending verification",
            finding_id=finding_id,
            review_ids=[review_id],
            idempotency_key="decision-pending",
        )
    passed = _pass_verification(
        ops,
        review_id,
        key="verification-pass",
        supersedes=pending,
    )
    decision_id = ops.create_promotion_decision(
        campaign_id=campaign_id,
        decision="publish",
        source_tier="canonical",
        policy_name="fixture",
        policy_version="1",
        policy_hash="policy",
        rationale="all gates passed",
        finding_id=finding_id,
        review_ids=[review_id],
        idempotency_key="decision-publish",
    )
    failed = ops.record_review_verification(
        review_id=review_id,
        claim_key="fixture-count",
        engine="duckdb",
        specification={"sql": "SELECT COUNT(*) FROM proteins"},
        dataset_id="dataset:sealed",
        expected=1,
        actual=2,
        status="fail",
        executed_ts=time.time(),
        supersedes_verification_id=passed,
        idempotency_key="verification-fail",
    )
    with pytest.raises(ValueError, match="unresolved verification"):
        ops.record_canonical_publication(
            campaign_id=campaign_id,
            finding_id=finding_id,
            decision_id=decision_id,
            dataset_id="dataset:sealed",
            canonical_uri="file:///findings.jsonl",
            canonical_record_id="F1",
            canonical_record_hash="record-hash",
            idempotency_key="publication",
        )
    _pass_verification(
        ops,
        review_id,
        key="verification-recovery",
        supersedes=failed,
    )
    publication_id = ops.record_canonical_publication(
        campaign_id=campaign_id,
        finding_id=finding_id,
        decision_id=decision_id,
        dataset_id="dataset:sealed",
        canonical_uri="file:///findings.jsonl",
        canonical_record_id="F1",
        canonical_record_hash="record-hash",
        idempotency_key="publication",
    )
    assert publication_id
    assert [row["attempt"] for row in ops.list_review_verifications(review_id)] == [
        1,
        2,
        3,
        4,
    ]


def test_atlas_closure_gate_requires_every_planned_unit(tmp_path) -> None:
    ops = OpsStore(tmp_path / "ops.db", agent_id="scanner")
    campaign_id = ops.create_campaign(
        "Atlas closure fixture",
        metadata={"unit_count": 2, "dataset_id": "dataset:sealed"},
        idempotency_key="atlas-campaign",
    )
    task_ids = [
        ops.create_task(
            "atlas_genome_read",
            f"Read g{index}",
            campaign_id=campaign_id,
            idempotency_key=f"atlas-task-{index}",
        )
        for index in (1, 2)
    ]
    ops.claim_task(task_ids[0])
    ops.record_unit_disposition(
        campaign_id=campaign_id,
        unit_id="g1",
        dataset_id="dataset:sealed",
        genome_id="g1",
        coverage_hash="coverage-1",
        candidate_count=0,
        disposition="clear",
        evidence_bundle_hash="bundle-1",
        task_id=task_ids[0],
        idempotency_key="disposition-1",
    )
    ops.complete_task(task_ids[0])
    finding_id = ops.finding(
        "observation",
        "fixture",
        "Coverage-gated finding",
        campaign_id=campaign_id,
        idempotency_key="coverage-finding",
    )
    review_id = _review(
        ops,
        campaign_id,
        finding_id=finding_id,
        key="coverage-review",
    )
    _pass_verification(ops, review_id, key="coverage-verification")
    with pytest.raises(ValueError, match="1/2 scan tasks complete"):
        ops.create_promotion_decision(
            campaign_id=campaign_id,
            decision="publish",
            source_tier="canonical",
            policy_name="fixture",
            policy_version="1",
            policy_hash="policy",
            rationale="coverage incomplete",
            finding_id=finding_id,
            review_ids=[review_id],
            idempotency_key="coverage-decision",
        )

    ops.claim_task(task_ids[1])
    ops.record_unit_disposition(
        campaign_id=campaign_id,
        unit_id="g2",
        dataset_id="dataset:sealed",
        genome_id="g2",
        coverage_hash="coverage-2",
        candidate_count=0,
        disposition="clear",
        evidence_bundle_hash="bundle-2",
        task_id=task_ids[1],
        idempotency_key="disposition-2",
    )
    ops.complete_task(task_ids[1])
    closure = ops.review_campaign_scan_status(campaign_id)
    assert closure == {
        "required": True,
        "ready": True,
        "expected_units": 2,
        "atlas_tasks": 2,
        "completed_atlas_tasks": 2,
        "active_unit_dispositions": 2,
        "ready_unit_dispositions": 2,
    }
    assert ops.create_promotion_decision(
        campaign_id=campaign_id,
        decision="publish",
        source_tier="canonical",
        policy_name="fixture",
        policy_version="1",
        policy_hash="policy",
        rationale="coverage complete",
        finding_id=finding_id,
        review_ids=[review_id],
        idempotency_key="coverage-decision",
    )


def test_controller_routes_blind_cross_provider_quorum(tmp_path) -> None:
    ops = OpsStore(tmp_path / "ops.db", agent_id="controller")
    campaign_id = _campaign(ops)
    _candidate(ops, campaign_id, unit_id="g1")
    cluster_id = ExactCandidateReducer().reduce_campaign(
        ops, campaign_id
    )
    assert cluster_id.clusters_created == 1
    cluster_id = ops.list_candidate_clusters(campaign_id=campaign_id)[0]["id"]
    policy = load_review_policy()
    controller = ReviewController(ops, policy)
    assert controller.tick(campaign_id).tasks_created == 1
    deep_task = next(
        task
        for task in ops.list_tasks(campaign_id=campaign_id, limit=100)
        if task["params"].get("review_tier") == "deepen"
    )
    ops.agent_id = "deep-reviewer"
    ops.claim_task(
        deep_task["id"],
        capabilities=[f"profile:{deep_task['params']['execution_profile']}"],
    )
    params = deep_task["params"]
    with pytest.raises(ValueError, match="execution identity"):
        _review(
            ops,
            campaign_id,
            cluster_id=cluster_id,
            tier="deepen",
            profile=params["execution_profile"],
            provider="wrong-provider",
            model=params["resolved_execution"]["model"],
            task_id=deep_task["id"],
            key="bad-deep-review",
        )
    deep_review = ops.create_finding_review(
        campaign_id=campaign_id,
        dataset_id="dataset:sealed",
        review_tier="deepen",
        execution_profile=params["execution_profile"],
        provider=params["resolved_execution"]["provider"],
        model=params["resolved_execution"]["model"],
        reasoning_effort=params["resolved_execution"]["reasoning_effort"],
        prompt_hash="deep-prompt",
        rubric_version=params["rubric_version"],
        input_bundle_hash="deep-input",
        verdict="promote",
        confidence=0.95,
        cluster_id=cluster_id,
        task_id=deep_task["id"],
        idempotency_key="deep-review",
    )
    _pass_verification(ops, deep_review, key="deep-verification")
    ops.agent_id = "controller"
    assert controller.tick(campaign_id).tasks_created == 0
    assert not [
        task
        for task in ops.list_tasks(campaign_id=campaign_id, limit=100)
        if task["params"].get("review_tier") == "independent"
    ]
    ops.agent_id = "deep-reviewer"
    assert ops.complete_task(deep_task["id"])["status"] == "complete"
    ops.agent_id = "controller"
    routed = controller.tick(campaign_id)
    assert routed.tasks_created == 2
    independent_tasks = [
        task
        for task in ops.list_tasks(campaign_id=campaign_id, limit=100)
        if task["params"].get("review_tier") == "independent"
    ]
    assert {
        task["params"]["resolved_execution"]["provider"]
        for task in independent_tasks
    } == {"openai", "anthropic"}
    assert all(task["params"]["blind_to_other_reviews"] for task in independent_tasks)
    for index, task in enumerate(independent_tasks):
        ops.agent_id = f"independent-{index}"
        ops.claim_task(
            task["id"],
            capabilities=[f"profile:{task['params']['execution_profile']}"],
        )
        params = task["params"]
        review_id = ops.create_finding_review(
            campaign_id=campaign_id,
            dataset_id="dataset:sealed",
            review_tier="independent",
            execution_profile=params["execution_profile"],
            provider=params["resolved_execution"]["provider"],
            model=params["resolved_execution"]["model"],
            reasoning_effort=params["resolved_execution"]["reasoning_effort"],
            prompt_hash=f"independent-prompt-{index}",
            rubric_version=params["rubric_version"],
            input_bundle_hash=f"independent-input-{index}",
            verdict="promote",
            confidence=0.95,
            cluster_id=cluster_id,
            task_id=task["id"],
            idempotency_key=f"independent-review-{index}",
        )
        _pass_verification(
            ops,
            review_id,
            key=f"independent-verification-{index}",
        )
        assert ops.complete_task(task["id"])["status"] == "complete"
    ops.agent_id = "controller"
    adjudicated = controller.tick(campaign_id)
    assert adjudicated.tasks_created == 1
    adjudication_task = next(
        task
        for task in ops.list_tasks(campaign_id=campaign_id, limit=100)
        if task["params"].get("review_tier") == "adjudication"
    )
    assert len(adjudication_task["params"]["source_review_ids"]) == 2
    decisions = ops.list_promotion_decisions(campaign_id=campaign_id)
    assert {(row["source_tier"], row["target_tier"]) for row in decisions} == {
        ("deepen", "independent"),
        ("independent", "adjudication"),
    }

    ops.agent_id = "adjudicator"
    ops.claim_task(
        adjudication_task["id"],
        capabilities=[
            f"profile:{adjudication_task['params']['execution_profile']}"
        ],
    )
    params = adjudication_task["params"]
    adjudication_review = ops.create_finding_review(
        campaign_id=campaign_id,
        dataset_id="dataset:sealed",
        review_tier="adjudication",
        execution_profile=params["execution_profile"],
        provider=params["resolved_execution"]["provider"],
        model=params["resolved_execution"]["model"],
        reasoning_effort=params["resolved_execution"]["reasoning_effort"],
        prompt_hash="adjudication-prompt",
        rubric_version=params["rubric_version"],
        input_bundle_hash="adjudication-input",
        verdict="promote",
        confidence=0.95,
        cluster_id=cluster_id,
        task_id=adjudication_task["id"],
        blind_to_prior_scores=params["blind_to_prior_scores"],
        blind_to_other_reviews=params["blind_to_other_reviews"],
        idempotency_key="adjudication-review",
    )
    _pass_verification(
        ops,
        adjudication_review,
        key="adjudication-verification",
    )
    ops.complete_task(adjudication_task["id"])

    ops.agent_id = "controller"
    assert controller.tick(campaign_id).tasks_created == 1
    materialization_task = next(
        task
        for task in ops.list_tasks(campaign_id=campaign_id, limit=100)
        if task["task_type"] == "materialize_finding"
    )
    ops.agent_id = "synthesizer"
    ops.claim_task(
        materialization_task["id"],
        capabilities=[
            f"profile:{materialization_task['params']['execution_profile']}"
        ],
    )
    finding_id = ops.finding(
        "reviewed_candidate",
        "architecture",
        "Canonical synthesis fixture",
        campaign_id=campaign_id,
        task_id=materialization_task["id"],
        idempotency_key="canonical-synthesis",
        validation_status="review_dag_pending",
    )
    ops.link_cluster_finding(cluster_id, finding_id, relation="materializes")
    ops.agent_id = "controller"
    assert controller.tick(campaign_id).tasks_created == 0
    ops.agent_id = "synthesizer"
    ops.complete_task(materialization_task["id"], [finding_id])

    ops.agent_id = "controller"
    assert controller.tick(campaign_id).tasks_created == 1
    canonical_task = next(
        task
        for task in ops.list_tasks(campaign_id=campaign_id, limit=100)
        if task["params"].get("review_tier") == "canonical"
    )
    ops.agent_id = "canonical-reviewer"
    ops.claim_task(
        canonical_task["id"],
        capabilities=[f"profile:{canonical_task['params']['execution_profile']}"],
    )
    params = canonical_task["params"]
    canonical_review = ops.create_finding_review(
        campaign_id=campaign_id,
        dataset_id="dataset:sealed",
        review_tier="canonical",
        execution_profile=params["execution_profile"],
        provider=params["resolved_execution"]["provider"],
        model=params["resolved_execution"]["model"],
        reasoning_effort=params["resolved_execution"]["reasoning_effort"],
        prompt_hash="canonical-prompt",
        rubric_version=params["rubric_version"],
        input_bundle_hash="canonical-input",
        verdict="promote",
        confidence=0.95,
        finding_id=finding_id,
        task_id=canonical_task["id"],
        blind_to_prior_scores=params["blind_to_prior_scores"],
        blind_to_other_reviews=params["blind_to_other_reviews"],
        idempotency_key="canonical-review",
    )
    _pass_verification(
        ops,
        canonical_review,
        key="canonical-verification",
    )
    ops.agent_id = "controller"
    assert controller.tick(campaign_id).decisions_created == 0
    ops.agent_id = "canonical-reviewer"
    ops.complete_task(canonical_task["id"])

    ops.agent_id = "controller"
    controller.tick(campaign_id)
    publish_decision = next(
        row
        for row in ops.list_promotion_decisions(
            campaign_id=campaign_id,
            finding_id=finding_id,
        )
        if row["decision"] == "publish"
    )
    publication_id = ops.record_canonical_publication(
        campaign_id=campaign_id,
        finding_id=finding_id,
        decision_id=publish_decision["id"],
        dataset_id="dataset:sealed",
        canonical_uri="file:///canonical/findings.jsonl",
        canonical_record_id="F-CANONICAL",
        canonical_record_hash="canonical-hash",
        idempotency_key="canonical-publication",
    )
    assert publication_id


def test_atlas_review_output_contract_blocks_incomplete_task(tmp_path) -> None:
    ops = OpsStore(tmp_path / "ops.db", agent_id="worker")
    campaign_id = _campaign(ops)
    task_id = ops.create_task(
        "atlas_genome_read",
        "Read genome g1",
        campaign_id=campaign_id,
        params={
            "unit_id": "g1",
            "dataset_id": "dataset:sealed",
            "genome_id": "g1",
            "review_output_contract": {
                "schema_version": "atlas-review-output/1.0"
            },
        },
    )
    ops.claim_task(task_id)
    with pytest.raises(ValueError, match="active unit disposition"):
        ops.complete_task(task_id)
    with pytest.raises(ValueError, match="unit_id, genome_id"):
        _candidate(ops, campaign_id, unit_id="g2", task_id=task_id)
    with pytest.raises(ValueError, match="dataset_id"):
        ops.record_unit_disposition(
            campaign_id=campaign_id,
            unit_id="g1",
            dataset_id="dataset:other",
            genome_id="g1",
            coverage_hash="wrong-target",
            candidate_count=0,
            disposition="clear",
            evidence_bundle_hash="wrong-target",
            task_id=task_id,
            idempotency_key="wrong-target",
        )
    _candidate(ops, campaign_id, unit_id="g1", task_id=task_id)
    ops.record_unit_disposition(
        campaign_id=campaign_id,
        unit_id="g1",
        dataset_id="dataset:sealed",
        genome_id="g1",
        coverage_hash="coverage",
        candidate_count=1,
        disposition="candidate",
        evidence_bundle_hash="bundle",
        task_id=task_id,
        idempotency_key="atlas-disposition",
    )
    assert ops.complete_task(task_id)["status"] == "complete"


def test_task_scoped_review_outputs_are_idempotent_across_retry_agents(
    tmp_path,
) -> None:
    path = tmp_path / "ops.db"
    first_worker = OpsStore(path, agent_id="worker-a")
    campaign_id = _campaign(first_worker)
    task_id = first_worker.create_task(
        "atlas_genome_read",
        "Read genome g1",
        campaign_id=campaign_id,
    )
    first_worker.claim_task(task_id)
    candidate_id = _candidate(
        first_worker,
        campaign_id,
        unit_id="g1",
        task_id=task_id,
    )
    disposition_id = first_worker.record_unit_disposition(
        campaign_id=campaign_id,
        unit_id="g1",
        dataset_id="dataset:sealed",
        genome_id="g1",
        coverage_hash="coverage",
        candidate_count=1,
        disposition="candidate",
        evidence_bundle_hash="bundle",
        task_id=task_id,
        idempotency_key="disposition",
    )
    first_worker.fail_task(task_id, retryable=True, error="transport failure")
    first_worker.close()

    second_worker = OpsStore(path, agent_id="worker-b")
    second_worker.claim_task(task_id)
    assert (
        _candidate(
            second_worker,
            campaign_id,
            unit_id="g1",
            task_id=task_id,
        )
        == candidate_id
    )
    assert (
        second_worker.record_unit_disposition(
            campaign_id=campaign_id,
            unit_id="g1",
            dataset_id="dataset:sealed",
            genome_id="g1",
            coverage_hash="coverage",
            candidate_count=1,
            disposition="candidate",
            evidence_bundle_hash="bundle",
            task_id=task_id,
            idempotency_key="disposition",
        )
        == disposition_id
    )
    assert second_worker.complete_task(task_id)["status"] == "complete"


def test_review_and_materialization_tasks_fail_closed_on_missing_outputs(
    tmp_path,
) -> None:
    ops = OpsStore(tmp_path / "ops.db", agent_id="worker")
    campaign_id = _campaign(ops)
    cluster_id = _cluster_fixture(ops, campaign_id)
    review_task_id = ops.create_task(
        "scientific_review",
        "Review cluster",
        campaign_id=campaign_id,
        params={
            "review_contract": "sharur-review-task/1.0",
            "target": {"kind": "candidate_cluster", "id": cluster_id},
            "review_tier": "canonical",
            "execution_profile": "canonical_review",
            "resolved_execution": {
                "provider": "anthropic",
                "model": "fixture-model",
                "reasoning_effort": "high",
            },
            "blind_to_prior_scores": True,
            "blind_to_other_reviews": True,
        },
    )
    ops.claim_task(review_task_id)
    with pytest.raises(ValueError, match="exactly one task-owned review"):
        ops.complete_task(review_task_id)
    review_id = _review(
        ops,
        campaign_id,
        cluster_id=cluster_id,
        task_id=review_task_id,
        key="task-review",
    )
    assert ops.complete_task(review_task_id)["status"] == "complete"
    assert ops.get_finding_review(review_id)["task_id"] == review_task_id

    materialization_task_id = ops.create_task(
        "materialize_finding",
        "Materialize cluster",
        campaign_id=campaign_id,
        params={
            "materialization_contract": "sharur-finding-materialization/1.0",
            "target": {"kind": "candidate_cluster", "id": cluster_id},
            "output_contract": {
                "link_cluster_relation": "materializes",
                "validation_status": "review_dag_pending",
            },
        },
    )
    ops.claim_task(materialization_task_id)
    finding_id = ops.finding(
        "reviewed_candidate",
        "test",
        "Materialized fixture",
        campaign_id=campaign_id,
        task_id=materialization_task_id,
        idempotency_key="materialized-finding",
        validation_status="review_dag_pending",
    )
    with pytest.raises(ValueError, match="source-cluster link"):
        ops.complete_task(materialization_task_id, [finding_id])
    ops.link_cluster_finding(cluster_id, finding_id)
    assert (
        ops.complete_task(materialization_task_id, [finding_id])["status"]
        == "complete"
    )


def test_verifier_is_bounded_read_only_and_sequence_safe(tmp_path) -> None:
    path = tmp_path / "fixture.duckdb"
    connection = duckdb.connect(str(path))
    connection.execute("CREATE TABLE values_fixture(value INTEGER)")
    connection.execute("INSERT INTO values_fixture VALUES (1), (2)")
    connection.close()
    passed = run_duckdb_verification(
        path,
        {"sql": "SELECT COUNT(*) FROM values_fixture"},
        2,
        dataset_id="fixture",
        verify_seal=False,
    )
    assert passed.status == "pass"
    assert passed.actual == 2
    blocked = run_duckdb_verification(
        path,
        {"sql": "SELECT repeat('A', 40) AS payload"},
        "irrelevant",
        dataset_id="fixture",
        verify_seal=False,
    )
    assert blocked.status == "error"
    assert "raw biological sequence" in blocked.error
    with pytest.raises(ValueError, match="forbidden operation"):
        run_duckdb_verification(
            path,
            {
                "sql": "SELECT * FROM read_text('/etc/hosts')",
                "result_shape": "rows",
            },
            [],
            dataset_id="fixture",
            verify_seal=False,
        )


def test_model_visible_review_records_reject_raw_sequences(tmp_path) -> None:
    ops = OpsStore(tmp_path / "ops.db", agent_id="worker")
    campaign_id = _campaign(ops)
    raw_sequence = "ACGT" * 20
    with pytest.raises(ValueError, match="raw biological sequence"):
        ops.finding(
            "unsafe",
            "test",
            raw_sequence,
            campaign_id=campaign_id,
            idempotency_key="unsafe-finding",
        )
    with pytest.raises(ValueError, match="raw biological sequence"):
        ops.create_candidate_occurrence(
            campaign_id=campaign_id,
            dataset_id="dataset:sealed",
            unit_id="g1",
            genome_id="g1",
            candidate_type="unsafe",
            signature_schema="unsafe/v1",
            signature={"domains": ["A"]},
            evidence={"sequence": raw_sequence},
            verification=[],
            subject_refs={"protein_id": "p1"},
            idempotency_key="unsafe-candidate",
        )
    cluster_id = _cluster_fixture(ops, campaign_id)
    with pytest.raises(ValueError, match="raw biological sequence"):
        ops.create_finding_review(
            campaign_id=campaign_id,
            dataset_id="dataset:sealed",
            review_tier="canonical",
            execution_profile="canonical_review",
            provider="anthropic",
            model="fixture-model",
            reasoning_effort="high",
            prompt_hash="unsafe-prompt",
            rubric_version="scientific-review/1.0",
            input_bundle_hash="unsafe-input",
            verdict="hold",
            confidence=0.5,
            cluster_id=cluster_id,
            reconstructed_observations={"sequence": raw_sequence},
            idempotency_key="unsafe-review",
        )
    review_id = _review(
        ops,
        campaign_id,
        cluster_id=cluster_id,
        key="safe-review",
    )
    with pytest.raises(ValueError, match="raw biological sequence"):
        ops.record_review_verification(
            review_id=review_id,
            claim_key="unsafe-result",
            engine="duckdb",
            specification={"sql": "SELECT 1"},
            dataset_id="dataset:sealed",
            expected=raw_sequence,
            actual=raw_sequence,
            status="pass",
            executed_ts=time.time(),
            idempotency_key="unsafe-verification",
        )


def test_deterministic_audit_sampling_and_http_review_surface(tmp_path) -> None:
    first = deterministic_sample(
        0.5,
        seed="seed",
        campaign_id="campaign",
        subject_id="unit",
        stratum={"class": "A"},
    )
    assert first == deterministic_sample(
        0.5,
        seed="seed",
        campaign_id="campaign",
        subject_id="unit",
        stratum={"class": "A"},
    )
    app = create_app(db_path=tmp_path / "http.db", api_token=TOKEN)
    headers = {"Authorization": f"Bearer {TOKEN}"}
    with TestClient(app) as client:
        campaign_id = client.post(
            "/campaigns",
            headers=headers,
            json={"name": "HTTP review", "idempotency_key": "campaign"},
        ).json()["id"]
        candidate = client.post(
            "/review/candidates",
            headers=headers,
            json={
                "agent_id": "scanner",
                "campaign_id": campaign_id,
                "dataset_id": "dataset:sealed",
                "unit_id": "g1",
                "genome_id": "g1",
                "candidate_type": "architecture",
                "signature_schema": "architecture/v1",
                "signature": {"domains": ["A"]},
                "evidence": {"protein_id": "p1"},
                "verification": [],
                "subject_refs": {"protein_id": "p1"},
                "idempotency_key": "candidate",
            },
        )
        assert candidate.status_code == 201
        disposition = client.post(
            "/review/unit-dispositions",
            headers=headers,
            json={
                "agent_id": "scanner",
                "campaign_id": campaign_id,
                "unit_id": "g1",
                "dataset_id": "dataset:sealed",
                "genome_id": "g1",
                "coverage_hash": "coverage",
                "candidate_count": 1,
                "disposition": "candidate",
                "evidence_bundle_hash": "bundle",
                "idempotency_key": "disposition",
            },
        )
        assert disposition.status_code == 201
        reduced = client.post(
            "/review/reduce",
            headers=headers,
            json={"campaign_id": campaign_id},
        )
        assert reduced.json()["clusters_created"] == 1
        routed = client.post(
            "/review/controller/tick",
            headers=headers,
            json={"campaign_id": campaign_id},
        )
        assert routed.json()["tasks_created"] == 1
        status = client.get(
            "/review/status",
            headers=headers,
            params={"campaign_id": campaign_id},
        ).json()
        assert status["reduction"] == {
            "candidate_occurrences": 1,
            "active_clusters": 1,
            "occurrences_per_active_cluster": 1.0,
        }
        assert status["coverage"]["active_unit_dispositions"] == 1
        clusters = client.get(
            "/review/clusters",
            headers=headers,
            params={"campaign_id": campaign_id},
        ).json()
        assert len(clusters) == 1
        cluster_detail = client.get(
            f"/review/clusters/{clusters[0]['id']}",
            headers=headers,
            params={"member_limit": 1},
        ).json()
        assert cluster_detail["member_count"] == 1
        assert cluster_detail["members_truncated"] is False
        member_page = client.get(
            f"/review/clusters/{clusters[0]['id']}/members",
            headers=headers,
            params={"limit": 1},
        ).json()
        assert len(member_page["members"]) == 1
        assert member_page["has_more"] is False
        finding_ids = [
            client.post(
                "/findings",
                headers=headers,
                json={
                    "finding_type": "review_fixture",
                    "domain": "test",
                    "summary": f"Finding {index}",
                    "campaign_id": campaign_id,
                    "idempotency_key": f"finding-{index}",
                },
            ).json()["id"]
            for index in range(2)
        ]
        finding_link = client.post(
            f"/findings/{finding_ids[0]}/links",
            headers=headers,
            json={
                "related_finding_id": finding_ids[1],
                "relation": "supports",
            },
        )
        assert finding_link.status_code == 201
        assert finding_link.json() == {
            "finding_id": finding_ids[0],
            "related_finding_id": finding_ids[1],
            "relation": "supports",
        }


def test_batch_candidates_share_one_write_transaction(tmp_path):
    """A genome's candidates must cost ONE write lock acquisition, not N.

    Submitting ~19 candidates per genome individually meant ~19 BEGIN IMMEDIATE
    acquisitions of SQLite's single global write lock, which under an 8-worker
    fleet drove write-lock waits to 122 seconds and starved the connection pool
    until requests failed with 500s. `transaction_wait_observer` fires once per
    immediate transaction, so it counts exactly what we are trying to reduce.
    """
    immediate: list[float] = []
    ops = OpsStore(
        tmp_path / "ops" / "sharur_ops.db",
        agent_id="worker",
        transaction_wait_observer=immediate.append,
    )
    campaign_id = _campaign(ops)

    def payload(n: int) -> list[dict]:
        return [
            {
                "campaign_id": campaign_id,
                "dataset_id": "dataset:sealed",
                "unit_id": f"g{n}",
                "genome_id": f"g{n}",
                "candidate_type": "architecture",
                "signature_schema": "architecture/v1",
                "signature": {"domains": ["A", f"B{i}"]},
                "evidence": {"protein_id": f"g{n}_p{i}"},
                "verification": [],
                "subject_refs": {"genome_id": f"g{n}", "protein_id": f"g{n}_p{i}"},
                "reduction_features": {"strata": {"phylum": "P"}},
                "idempotency_key": f"candidate:g{n}:{i}",
            }
            for i in range(19)
        ]

    # one-at-a-time: one write transaction each
    immediate.clear()
    for item in payload(1):
        ops.create_candidate_occurrence(**item)
    individual = len(immediate)

    # batched: one write transaction for all of them
    immediate.clear()
    ids = ops.create_candidate_occurrences(payload(2))
    batched = len(immediate)

    assert len(ids) == 19
    assert individual == 19, f"expected one txn per candidate, got {individual}"
    assert batched == 1, f"batch must take the write lock once, took {batched}"

    # identical rows land either way
    assert len(ops.list_candidate_occurrences(campaign_id=campaign_id)) == 38

    # and the batch stays idempotent on replay
    immediate.clear()
    again = ops.create_candidate_occurrences(payload(2))
    assert again == ids, "replayed batch must return the original ids"
    assert len(ops.list_candidate_occurrences(campaign_id=campaign_id)) == 38
