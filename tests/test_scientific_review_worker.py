"""Contract tests for the packaged scientific-review executor."""

from __future__ import annotations

import contextlib
import json
from typing import Any

import pytest

from sharur.review.models import load_review_policy
from sharur.review.verification import VerificationResult
from sharur.workers.model_cli import ModelError, ModelRun
from sharur.workers.scientific_review import (
    REVIEW_OUTPUT_SCHEMA,
    SCIENTIFIC_CONTRACT,
    ScientificReviewWorker,
    _content_hash,
    _effective_verdict,
    _prepare_review,
)


def _model_record(*, verdict: str = "promote") -> dict[str, Any]:
    return {
        "reconstructed_observations_json": json.dumps(
            {"observed": ["PF00001 on protein p1"]}
        ),
        "claim_assessment_json": json.dumps(
            {"claim-1": {"status": "supported", "basis": "observed domain"}}
        ),
        "discrepancies_json": "[]",
        "proposed_tasks_json": "[]",
        "verdict": verdict,
        "confidence": 0.9,
        "verification_requests": [
            {
                "claim_key": "claim-1",
                "sql": "SELECT COUNT(*) FROM annotations WHERE protein_id = ?",
                "parameters_json": '["p1"]',
                "expected_json": "1",
                "result_shape": "scalar",
                "comparison": "exact",
                "absolute_tolerance": 0.0,
                "relative_tolerance": 0.0,
                "max_rows": 1,
            }
        ],
    }


def _model_run(*, verdict: str = "promote") -> ModelRun:
    return ModelRun(
        provider="openai",
        model="gpt-5.6-sol",
        reasoning_effort="medium",
        record=_model_record(verdict=verdict),
        prompt_sha256="a" * 64,
        result_sha256="b" * 64,
        exit_status=0,
        stderr="",
        usage={"input_tokens": 100},
        raw_stdout_bytes=500,
    )


def _verification(status: str = "pass") -> VerificationResult:
    return VerificationResult(
        status=status,
        actual=1 if status in {"pass", "fail"} else None,
        executed_ts=1.0,
        specification_hash="c" * 64,
        query_hash="d" * 64,
        row_count=1 if status in {"pass", "fail"} else None,
        error="query failed" if status == "error" else None,
    )


def _task_params(*, target_kind: str = "candidate_cluster") -> dict[str, Any]:
    policy = load_review_policy()
    profile = policy.profile("finding_deepen")
    return {
        "review_contract": "sharur-review-task/1.0",
        "campaign_id": "campaign-1",
        "dataset_id": "dataset:sealed",
        "target": {"kind": target_kind, "id": "cluster-1"},
        "target_input": {
            "id": "cluster-1",
            "version": 1,
            "member_manifest_hash": "manifest-1",
        },
        "review_tier": "deepen",
        "execution_profile": "finding_deepen",
        "resolved_execution": {
            "provider": profile.provider,
            "model": profile.model,
            "reasoning_effort": profile.reasoning_effort,
        },
        "policy": {
            "name": policy.name,
            "version": policy.version,
            "hash": policy.policy_hash,
        },
        "rubric_version": policy.rubric_version,
        "blind_to_prior_scores": True,
        "blind_to_other_reviews": True,
        "source_review_ids": [],
        "source_review_manifest_hash": (
            "4f53cda18c2baa0c0354bb5f9a3ecbe5ed12ab4d8e11ba873c2f11161202b945"
        ),
        "audit": None,
        "scientific_contract": SCIENTIFIC_CONTRACT,
    }


def _task() -> dict[str, Any]:
    return {
        "id": "task-1",
        "task_type": "scientific_review",
        "campaign_id": "campaign-1",
        "params": _task_params(),
    }


def _bare_worker(ops: Any) -> ScientificReviewWorker:
    worker = object.__new__(ScientificReviewWorker)
    worker.agent_id = "reviewer-1"
    worker.profile_name = "finding_deepen"
    worker.policy = load_review_policy()
    worker.profile = worker.policy.profile("finding_deepen")
    worker.ops = ops
    worker.db_path = None
    worker.seal_path = None
    worker.lease_seconds = 1_800
    worker.model_timeout = 900
    worker.max_members = 24
    worker.max_input_bytes = 512 * 1024
    worker.verification_threads = 1
    worker.transient_retries = 0
    worker._dataset_context_cache = {
        "capabilities": [
            {
                "capability_id": "structured_callers",
                "state": "available",
                "evidence": {
                    "resources": [
                        {
                            "table": "defense_systems",
                            "rows": 2,
                            "emitted_types": [{"value": "caller-label", "rows": 2}],
                        }
                    ]
                },
            }
        ],
        "schema": [
            {
                "table": "annotations",
                "columns": [{"name": "protein_id", "type": "VARCHAR"}],
            }
        ],
    }
    worker._code_commit = "fixture"
    return worker


class _InputOps:
    def get_candidate_cluster(self, cluster_id, *, member_limit):
        assert cluster_id == "cluster-1"
        assert member_limit == 24
        return {
            "id": cluster_id,
            "campaign_id": "campaign-1",
            "dataset_id": "dataset:sealed",
            "version": 1,
            "member_manifest_hash": "manifest-1",
            "member_count": 1,
            "members": [{"candidate_id": "candidate-1", "role": "medoid"}],
        }

    def get_candidate_occurrence(self, candidate_id):
        assert candidate_id == "candidate-1"
        return {
            "id": candidate_id,
            "confidence": 0.98,
            "signature": {"accessions": ["PF00001"]},
            "evidence": {"protein_id": "p1"},
        }

    def get_finding_review(self, review_id):
        raise AssertionError(review_id)

    def list_review_verifications(self, review_id):
        raise AssertionError(review_id)


class _RunOps:
    def __init__(
        self,
        *,
        existing_review: dict[str, Any] | None = None,
        existing_verifications: list[dict[str, Any]] | None = None,
    ):
        self.existing_review = existing_review
        self.existing_verifications = existing_verifications or []
        self.created_review: dict[str, Any] | None = None
        self.recorded: list[dict[str, Any]] = []
        self.completed: list[str] = []

    def list_finding_reviews(self, **filters):
        assert filters["campaign_id"] == "campaign-1"
        return [self.existing_review] if self.existing_review else []

    def create_finding_review(self, **payload):
        self.created_review = payload
        return "review-1"

    def list_review_verifications(self, review_id):
        assert review_id == "review-1"
        return self.existing_verifications

    def record_review_verification(self, review_id, **payload):
        assert review_id == "review-1"
        self.recorded.append(payload)
        return "verification-1"

    def complete_task(self, task_id):
        self.completed.append(task_id)
        return {"id": task_id, "status": "complete"}


class TestScientificReviewSchema:
    def test_schema_is_openai_strict_compliant(self):
        def walk(node, path="root"):
            if not isinstance(node, dict):
                return
            if node.get("type") == "object":
                assert node.get("additionalProperties") is False, path
                assert set(node.get("required", [])) == set(node.get("properties", {})), path
                for key, child in node.get("properties", {}).items():
                    walk(child, f"{path}.{key}")
            if node.get("type") == "array":
                walk(node.get("items", {}), f"{path}[]")

        walk(REVIEW_OUTPUT_SCHEMA)

    def test_parser_builds_bounded_verification_spec(self):
        prepared = _prepare_review(_model_record(), max_checks=8)
        assert prepared.verdict == "promote"
        assert prepared.checks[0].specification.result_shape == "scalar"
        assert prepared.checks[0].expected == 1

    def test_parser_rejects_mutating_sql(self):
        record = _model_record()
        record["verification_requests"][0]["sql"] = "DELETE FROM annotations"
        with pytest.raises(ModelError, match="failed validation"):
            _prepare_review(record, max_checks=8)

    def test_parser_enforces_policy_check_limit(self):
        with pytest.raises(ModelError, match="policy permits 0"):
            _prepare_review(_model_record(), max_checks=0)


class TestReviewInput:
    def test_blind_input_keeps_member_observations_and_hides_prior_score(self):
        worker = _bare_worker(_InputOps())
        payload = worker._build_input(_task(), _task_params())

        members = payload["target_evidence"]["cluster_members"]
        assert len(members) == 1
        assert members[0]["occurrence"]["signature"] == {"accessions": ["PF00001"]}
        assert "confidence" not in members[0]["occurrence"]
        capabilities = payload["dataset_context"]["capabilities"]
        assert capabilities[0]["capability_id"] == "structured_callers"

    def test_contract_rejects_execution_drift(self):
        worker = _bare_worker(_InputOps())
        task = _task()
        task["params"]["resolved_execution"]["model"] = "different-model"
        with pytest.raises(ValueError, match="execution identity"):
            worker._validate_contract(task)


class TestReviewExecution:
    @pytest.mark.parametrize(
        ("status", "expected_verdict"),
        [("pass", "promote"), ("fail", "hold"), ("error", "hold")],
    )
    def test_verification_result_gates_decisive_verdict(
        self,
        monkeypatch,
        status,
        expected_verdict,
    ):
        ops = _RunOps()
        worker = _bare_worker(ops)
        monkeypatch.setattr(worker, "_lease_keepalive", lambda task_id: contextlib.nullcontext())
        monkeypatch.setattr(worker, "_verify_dataset", lambda dataset_id: None)
        monkeypatch.setattr(worker, "_build_input", lambda task, params: {"safe": True})
        monkeypatch.setattr(worker, "_call_model", lambda payload, task_id: _model_run())
        monkeypatch.setattr(
            worker,
            "_run_checks",
            lambda checks, dataset_id: [_verification(status)],
        )

        worker.run_task(_task())

        assert ops.created_review is not None
        assert ops.created_review["verdict"] == expected_verdict
        assert ops.created_review["verification_summary"]["requested_verdict"] == "promote"
        assert ops.recorded[0]["status"] == status
        assert ops.completed == ["task-1"]

    def test_resume_uses_committed_review_and_skips_model(self, monkeypatch):
        prepared = _prepare_review(_model_record(), max_checks=8)
        summary = {
            "checks": [
                {
                    "request": prepared.checks[0].to_request_dict(),
                    "result": {"status": "pass"},
                }
            ]
        }
        existing = {
            "id": "review-1",
            "task_id": "task-1",
            "verification_summary": summary,
        }
        ops = _RunOps(existing_review=existing)
        worker = _bare_worker(ops)
        monkeypatch.setattr(worker, "_lease_keepalive", lambda task_id: contextlib.nullcontext())
        monkeypatch.setattr(worker, "_verify_dataset", lambda dataset_id: None)
        monkeypatch.setattr(
            worker,
            "_call_model",
            lambda payload, task_id: pytest.fail("resume must skip the model"),
        )
        monkeypatch.setattr(
            worker,
            "_run_checks",
            lambda checks, dataset_id: [_verification("pass")],
        )

        worker.run_task(_task())

        assert ops.created_review is None
        assert len(ops.recorded) == 1
        assert ops.completed == ["task-1"]

    def test_resume_skips_checks_already_recorded(self, monkeypatch):
        prepared = _prepare_review(_model_record(), max_checks=8)
        check = prepared.checks[0]
        summary = {
            "checks": [
                {
                    "request": check.to_request_dict(),
                    "result": {"status": "pass"},
                }
            ]
        }
        existing = {
            "id": "review-1",
            "task_id": "task-1",
            "verification_summary": summary,
        }
        ops = _RunOps(
            existing_review=existing,
            existing_verifications=[
                {
                    "claim_key": check.claim_key,
                    "specification_hash": _content_hash(
                        check.specification.model_dump(mode="json")
                    ),
                }
            ],
        )
        worker = _bare_worker(ops)
        monkeypatch.setattr(worker, "_lease_keepalive", lambda task_id: contextlib.nullcontext())
        monkeypatch.setattr(worker, "_verify_dataset", lambda dataset_id: None)
        monkeypatch.setattr(
            worker,
            "_run_checks",
            lambda checks, dataset_id: pytest.fail("recorded checks must be skipped"),
        )

        worker.run_task(_task())

        assert ops.recorded == []
        assert ops.completed == ["task-1"]


def test_effective_verdict_preserves_needs_data_on_check_error():
    assert _effective_verdict("needs_data", [_verification("error")]) == "needs_data"
