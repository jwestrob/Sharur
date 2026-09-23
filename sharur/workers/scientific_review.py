"""Packaged executor for task-owned scientific reviews.

The review controller already creates frozen ``scientific_review`` tasks and
the Ops store already enforces their output contract.  This worker closes the
execution gap between those two surfaces:

1. claim a task under its exact execution-profile capability;
2. verify the frozen policy, target, blindness, and sealed-dataset contract;
3. assemble a bounded, sequence-free record from Ops plus the live caller and
   table schema discovered from the dataset;
4. ask the resolved model to reconstruct observations, assess claims, and
   specify executable checks;
5. run those checks through the bounded read-only DuckDB verifier;
6. append one task-owned review and its exact verification results;
7. complete the task under the same lease fence.

Review records are immutable.  The worker therefore stores each requested
check in ``verification_summary`` before appending the corresponding result.
A replacement attempt can recover a committed review, rerun only missing
checks against the same sealed dataset, and complete the task without another
model call.
"""

from __future__ import annotations

import contextlib
import hashlib
import json
import logging
import random
import signal
import subprocess
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import duckdb
import requests
from pydantic import BaseModel, ConfigDict, Field, ValidationError, model_validator

from sharur.capabilities import build_capability_brief
from sharur.dataset_seal import DEFAULT_SEAL_NAME, verify_dataset_seal
from sharur.ops.client import SharurOps
from sharur.ops.review_store import assert_sequence_free
from sharur.review.models import ReviewPolicy, load_review_policy
from sharur.review.verification import (
    RAW_SEQUENCE_COLUMN,
    VerificationResult,
    VerificationSpec,
    run_duckdb_verification,
)
from sharur.workers.model_cli import (
    ModelError,
    ModelQuotaExhausted,
    ModelRateLimited,
    ModelRun,
    ModelTransient,
    run_profile,
)


LOGGER = logging.getLogger("sharur.workers.scientific_review")

TASK_TYPE = "scientific_review"
REVIEW_TASK_CONTRACT = "sharur-review-task/1.0"
INPUT_SCHEMA = "sharur-scientific-review-input/1.0"
DEFAULT_MAX_INPUT_BYTES = 512 * 1024
MAX_SQL_BYTES = 16 * 1024
TARGET_KINDS = frozenset({"finding", "candidate_cluster", "unit_disposition"})
SCIENTIFIC_CONTRACT = {
    "reconstruct_observations": True,
    "discover_live_curated_callers_before_named_claims": True,
    "executable_verification_records": True,
    "raw_sequences_model_visible": False,
}


class _StrictModel(BaseModel):
    model_config = ConfigDict(extra="forbid")


class _ReviewCheckCarrier(_StrictModel):
    claim_key: str = Field(min_length=1, max_length=1_024)
    sql: str = Field(min_length=1, max_length=MAX_SQL_BYTES)
    parameters_json: str = Field(min_length=1, max_length=65_536)
    expected_json: str = Field(min_length=1, max_length=65_536)
    result_shape: str
    comparison: str
    absolute_tolerance: float = Field(ge=0.0)
    relative_tolerance: float = Field(ge=0.0)
    max_rows: int = Field(ge=1, le=100)


class _ReviewOutputCarrier(_StrictModel):
    reconstructed_observations_json: str = Field(min_length=2, max_length=262_144)
    claim_assessment_json: str = Field(min_length=2, max_length=262_144)
    discrepancies_json: str = Field(min_length=2, max_length=262_144)
    proposed_tasks_json: str = Field(min_length=2, max_length=262_144)
    verdict: str
    confidence: float = Field(ge=0.0, le=1.0)
    verification_requests: list[_ReviewCheckCarrier] = Field(min_length=1, max_length=8)

    @model_validator(mode="after")
    def validate_enums(self) -> _ReviewOutputCarrier:
        verdicts = {"promote", "hold", "needs_data", "reject", "duplicate", "split"}
        if self.verdict not in verdicts:
            raise ValueError(f"Unsupported verdict {self.verdict!r}")
        return self


REVIEW_OUTPUT_SCHEMA: dict[str, Any] = {
    "type": "object",
    "additionalProperties": False,
    "required": [
        "reconstructed_observations_json",
        "claim_assessment_json",
        "discrepancies_json",
        "proposed_tasks_json",
        "verdict",
        "confidence",
        "verification_requests",
    ],
    "properties": {
        "reconstructed_observations_json": {
            "type": "string",
            "description": "JSON object containing reconstructed observations",
        },
        "claim_assessment_json": {
            "type": "string",
            "description": "JSON object mapping claims to support and uncertainty",
        },
        "discrepancies_json": {
            "type": "string",
            "description": "JSON array of discrepancies or counterevidence",
        },
        "proposed_tasks_json": {
            "type": "string",
            "description": "JSON array of bounded follow-up tasks",
        },
        "verdict": {
            "type": "string",
            "enum": ["promote", "hold", "needs_data", "reject", "duplicate", "split"],
        },
        "confidence": {"type": "number", "minimum": 0.0, "maximum": 1.0},
        "verification_requests": {
            "type": "array",
            "minItems": 1,
            "maxItems": 8,
            "items": {
                "type": "object",
                "additionalProperties": False,
                "required": [
                    "claim_key",
                    "sql",
                    "parameters_json",
                    "expected_json",
                    "result_shape",
                    "comparison",
                    "absolute_tolerance",
                    "relative_tolerance",
                    "max_rows",
                ],
                "properties": {
                    "claim_key": {"type": "string"},
                    "sql": {"type": "string", "maxLength": MAX_SQL_BYTES},
                    "parameters_json": {
                        "type": "string",
                        "description": "JSON array for positional SQL parameters",
                    },
                    "expected_json": {
                        "type": "string",
                        "description": "JSON value expected from the query",
                    },
                    "result_shape": {
                        "type": "string",
                        "enum": ["scalar", "row", "rows"],
                    },
                    "comparison": {
                        "type": "string",
                        "enum": ["exact", "approx", "contains", "set_equal"],
                    },
                    "absolute_tolerance": {"type": "number", "minimum": 0.0},
                    "relative_tolerance": {"type": "number", "minimum": 0.0},
                    "max_rows": {"type": "integer", "minimum": 1, "maximum": 100},
                },
            },
        },
    },
}


REVIEW_SYSTEM_PROMPT = """You are the scientific reviewer for one frozen Sharur review task.

Reconstruct the observations represented by the supplied target evidence, assess each claim at
the resolution supported by those observations, and return one structured review. The input is
sequence-free and bounded. Work from the supplied records and live dataset schema.

Keep OBSERVED and NAMED claims separate. OBSERVED claims describe exact domain hits, caller
rows, predicates, identifiers, and local arrangements. A NAMED biological system, family,
pathway, or mechanism requires an exact emitted call from a populated purpose-built caller
resource shown in dataset_context. When that caller evidence is absent, state the observations,
describe the compatible interpretation, and mark the interpretation UNVERIFIED. Inspect
whatever structured caller resources the live capability record provides before assigning a
name.

Supply executable DuckDB checks for the material claims in your assessment. Every concrete
count and every identifier-dependent assertion needs a corresponding check. Each check uses one
deterministic SELECT or WITH statement, exact live table and column names, positional parameters,
and a bounded result. Use COUNT(DISTINCT protein_id) when domain multiplicity could inflate a
protein count. Use joins for context queries. Verification SQL receives read-only access to the
sealed dataset with external access disabled.

The verdict describes the target after reconstructing its observations. `promote` means the
evidence supports advancement to the next review tier. `hold` preserves the target while a
material discrepancy is resolved. `needs_data` specifies a concrete missing evidence task.
`reject`, `duplicate`, and `split` apply when the supplied evidence directly supports those
dispositions. Sharur executes every requested check before persisting the review and converts a
decisive verdict to `hold` if any check fails or errors.

Return JSON matching the supplied schema. The four fields ending in `_json` must themselves
contain valid JSON of the stated container type. Biological sequences stay outside every field.
"""


@dataclass(frozen=True)
class PreparedCheck:
    claim_key: str
    specification: VerificationSpec
    expected: Any

    def to_request_dict(self) -> dict[str, Any]:
        return {
            "claim_key": self.claim_key,
            "specification": self.specification.model_dump(mode="json"),
            "expected": self.expected,
        }


@dataclass(frozen=True)
class PreparedReview:
    reconstructed_observations: dict[str, Any]
    claim_assessment: dict[str, Any]
    discrepancies: list[dict[str, Any]]
    proposed_tasks: list[dict[str, Any]]
    verdict: str
    confidence: float
    checks: tuple[PreparedCheck, ...]


def _canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)


def _content_hash(value: Any) -> str:
    return hashlib.sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def _json_bytes(value: Any) -> int:
    return len(_canonical_json(value).encode("utf-8"))


def _decode_json_field(raw: str, expected_type: type, field: str) -> Any:
    try:
        value = json.loads(raw)
    except json.JSONDecodeError as exc:
        raise ModelError(f"{field} was not valid JSON: {exc}") from exc
    if not isinstance(value, expected_type):
        raise ModelError(f"{field} must encode a {expected_type.__name__}")
    assert_sequence_free(value, field=field)
    return value


def _prepare_review(record: dict[str, Any], *, max_checks: int) -> PreparedReview:
    try:
        carrier = _ReviewOutputCarrier.model_validate(record)
    except ValidationError as exc:
        raise ModelError(f"scientific-review output failed validation: {exc}") from exc
    if len(carrier.verification_requests) > max_checks:
        raise ModelError(
            f"review requested {len(carrier.verification_requests)} checks; "
            f"policy permits {max_checks}"
        )

    observations = _decode_json_field(
        carrier.reconstructed_observations_json,
        dict,
        "reconstructed_observations_json",
    )
    assessment = _decode_json_field(
        carrier.claim_assessment_json,
        dict,
        "claim_assessment_json",
    )
    discrepancies = _decode_json_field(
        carrier.discrepancies_json,
        list,
        "discrepancies_json",
    )
    proposed_tasks = _decode_json_field(
        carrier.proposed_tasks_json,
        list,
        "proposed_tasks_json",
    )
    if any(not isinstance(item, dict) for item in discrepancies):
        raise ModelError("discrepancies_json must encode an array of objects")
    if any(not isinstance(item, dict) for item in proposed_tasks):
        raise ModelError("proposed_tasks_json must encode an array of objects")

    checks: list[PreparedCheck] = []
    seen_claims: set[str] = set()
    for request in carrier.verification_requests:
        if request.claim_key in seen_claims:
            raise ModelError(f"duplicate verification claim_key {request.claim_key!r}")
        seen_claims.add(request.claim_key)
        parameters = _decode_json_field(
            request.parameters_json,
            list,
            f"verification {request.claim_key} parameters_json",
        )
        try:
            expected = json.loads(request.expected_json)
        except json.JSONDecodeError as exc:
            raise ModelError(
                f"verification {request.claim_key} expected_json was invalid: {exc}"
            ) from exc
        assert_sequence_free(
            expected,
            field=f"verification {request.claim_key} expected_json",
        )
        try:
            specification = VerificationSpec.model_validate(
                {
                    "sql": request.sql,
                    "parameters": parameters,
                    "result_shape": request.result_shape,
                    "comparison": request.comparison,
                    "absolute_tolerance": request.absolute_tolerance,
                    "relative_tolerance": request.relative_tolerance,
                    "max_rows": request.max_rows,
                }
            )
        except ValidationError as exc:
            raise ModelError(
                f"verification {request.claim_key} failed validation: {exc}"
            ) from exc
        checks.append(
            PreparedCheck(
                claim_key=request.claim_key,
                specification=specification,
                expected=expected,
            )
        )

    prepared = PreparedReview(
        reconstructed_observations=observations,
        claim_assessment=assessment,
        discrepancies=discrepancies,
        proposed_tasks=proposed_tasks,
        verdict=carrier.verdict,
        confidence=carrier.confidence,
        checks=tuple(checks),
    )
    assert_sequence_free(
        {
            "reconstructed_observations": prepared.reconstructed_observations,
            "claim_assessment": prepared.claim_assessment,
            "discrepancies": prepared.discrepancies,
            "proposed_tasks": prepared.proposed_tasks,
        },
        field="Prepared scientific review",
    )
    return prepared


def _strip_prior_assessments(value: Any) -> Any:
    """Remove prior evaluative fields while retaining primary observations."""

    hidden = {
        "claim_assessment",
        "confidence",
        "novelty",
        "score",
        "scores",
        "uncertainty",
        "validation_status",
        "verdict",
        "verification_summary",
    }
    if isinstance(value, dict):
        return {
            key: _strip_prior_assessments(item)
            for key, item in value.items()
            if str(key).lower() not in hidden
        }
    if isinstance(value, list):
        return [_strip_prior_assessments(item) for item in value]
    return value


def _schema_catalog(db_path: Path) -> list[dict[str, Any]]:
    connection = duckdb.connect(
        str(db_path),
        read_only=True,
        config={"enable_external_access": "false"},
    )
    try:
        rows = connection.execute(
            """
            SELECT table_name, column_name, data_type
            FROM information_schema.columns
            WHERE table_schema = 'main'
            ORDER BY table_name, ordinal_position
            """
        ).fetchall()
    finally:
        connection.close()
    tables: dict[str, list[dict[str, str]]] = {}
    for table_name, column_name, data_type in rows:
        if RAW_SEQUENCE_COLUMN.fullmatch(str(column_name)):
            continue
        tables.setdefault(str(table_name), []).append(
            {"name": str(column_name), "type": str(data_type)}
        )
    return [
        {"table": table_name, "columns": columns}
        for table_name, columns in sorted(tables.items())
    ]


def _scientific_capabilities(db_path: Path) -> list[dict[str, Any]]:
    brief = build_capability_brief(
        db_path,
        include_tools=False,
        include_execution=False,
    )
    selected = {
        "annotation_sources",
        "structured_callers",
        "semantic_v2",
        "predicate_compatibility",
    }
    return [
        capability.to_dict()
        for capability in brief.capabilities
        if capability.capability_id in selected
    ]


def _git_commit() -> str | None:
    root = Path(__file__).resolve().parents[2]
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=root,
            check=True,
            capture_output=True,
            text=True,
            timeout=5,
        ).stdout.strip()
        dirty = subprocess.run(
            ["git", "diff", "--quiet", "--"],
            cwd=root,
            check=False,
            timeout=5,
        ).returncode
    except (OSError, subprocess.SubprocessError):
        return None
    return f"{commit}+dirty" if dirty else commit


def _effective_verdict(requested: str, results: list[VerificationResult]) -> str:
    if all(result.status == "pass" for result in results):
        return requested
    return "needs_data" if requested == "needs_data" else "hold"


class ScientificReviewWorker:
    def __init__(
        self,
        *,
        ops_url: str,
        db_path: str | Path,
        agent_id: str,
        profile: str,
        policy_path: str | None = None,
        campaign_id: str | None = None,
        ops_token: str | None = None,
        seal_path: str | Path | None = None,
        lease_seconds: int = 1_800,
        model_timeout: int = 900,
        max_members: int = 24,
        max_input_bytes: int = DEFAULT_MAX_INPUT_BYTES,
        verification_threads: int = 1,
        transient_retries: int = 4,
    ) -> None:
        if max_members < 1:
            raise ValueError("max_members must be positive")
        if max_input_bytes < 16_384:
            raise ValueError("max_input_bytes must be at least 16384")
        if verification_threads < 1:
            raise ValueError("verification_threads must be positive")

        self.agent_id = agent_id
        self.profile_name = profile
        self.campaign_id = campaign_id
        self.db_path = Path(db_path).expanduser().resolve()
        self.seal_path = (
            Path(seal_path).expanduser().resolve()
            if seal_path is not None
            else self.db_path.parent / DEFAULT_SEAL_NAME
        )
        self.lease_seconds = lease_seconds
        self.model_timeout = model_timeout
        self.max_members = max_members
        self.max_input_bytes = max_input_bytes
        self.verification_threads = verification_threads
        self.transient_retries = transient_retries
        self.policy: ReviewPolicy = load_review_policy(
            Path(policy_path) if policy_path else None
        )
        if profile not in self.policy.profiles:
            raise SystemExit(
                f"profile {profile!r} not in policy {self.policy.name!r}; "
                f"available: {sorted(self.policy.profiles)}"
            )
        self.profile = self.policy.profile(profile)
        self.ops = SharurOps(
            ops_url,
            agent_id=agent_id,
            api_token=ops_token,
            timeout=(5.0, 180.0),
        )
        self._stop = False
        self._dataset_context_cache: dict[str, Any] | None = None
        self._code_commit = _git_commit()

    # ---------------------------------------------------------------- lifecycle

    def install_signal_handlers(self) -> None:
        def _handle(signum, _frame):
            LOGGER.info("signal %s received; finishing current task then stopping", signum)
            self._stop = True

        signal.signal(signal.SIGTERM, _handle)
        signal.signal(signal.SIGINT, _handle)

    def register(self) -> None:
        capabilities = ["scientific_reviewer", self.profile.capability]
        self.ops.register_agent(
            self.agent_id,
            role="worker",
            capabilities=capabilities,
            max_concurrent_tasks=1,
            capacity_cpu_slots=max(1, self.verification_threads),
        )
        LOGGER.info(
            "registered %s for %s -> %s/%s effort=%s",
            self.agent_id,
            self.profile_name,
            self.profile.provider,
            self.profile.model,
            self.profile.reasoning_effort,
        )

    def run_forever(
        self,
        *,
        idle_sleep: float = 5.0,
        max_tasks: int | None = None,
    ) -> int:
        self.register()
        completed = 0
        rate_limit_backoff = 0.0
        while not self._stop:
            if max_tasks is not None and completed >= max_tasks:
                break
            try:
                task = self.ops.claim_next_task(
                    campaign_id=self.campaign_id,
                    task_types=[TASK_TYPE],
                    lease_seconds=self.lease_seconds,
                )
            except requests.RequestException as exc:
                LOGGER.warning("review-task claim failed: %s", exc)
                time.sleep(idle_sleep)
                continue
            if task is None:
                time.sleep(idle_sleep)
                continue

            try:
                self.run_task(task)
                completed += 1
                rate_limit_backoff = 0.0
            except ModelQuotaExhausted as exc:
                LOGGER.error(
                    "usage limit reached%s; releasing %s and stopping: %s",
                    f" (resets {exc.reset_at})" if exc.reset_at else "",
                    task.get("id"),
                    str(exc)[:300],
                )
                self._release(task, f"usage limit: {exc}", retry_delay=300)
                self._stop = True
            except ModelRateLimited as exc:
                rate_limit_backoff = min(max(rate_limit_backoff * 2, 60.0), 1_800.0)
                delay = rate_limit_backoff + random.uniform(0, rate_limit_backoff * 0.1)
                LOGGER.warning(
                    "rate limited; releasing task %s and sleeping %.0fs",
                    task.get("id"),
                    delay,
                )
                self._release(task, str(exc), retry_delay=int(rate_limit_backoff))
                time.sleep(delay)
            except ModelTransient as exc:
                LOGGER.warning("transport remained unavailable for task %s", task.get("id"))
                self._release(task, f"transient transport: {exc}", retry_delay=120)
                time.sleep(120)
            except Exception as exc:
                LOGGER.exception("scientific-review task %s failed", task.get("id"))
                self._release(
                    task,
                    f"{type(exc).__name__}: {exc}",
                    retry_delay=30,
                )
        return completed

    def _release(self, task: dict[str, Any], error: str, *, retry_delay: int) -> None:
        try:
            self.ops.fail_task(
                str(task["id"]),
                error=error[:2_000],
                retryable=True,
                retry_delay_seconds=retry_delay,
            )
        except Exception:
            LOGGER.warning("could not release task %s", task.get("id"), exc_info=True)

    # -------------------------------------------------------------- task input

    def _validate_contract(self, task: dict[str, Any]) -> dict[str, Any]:
        if task.get("task_type") not in {None, TASK_TYPE}:
            raise ValueError(f"worker received task type {task.get('task_type')!r}")
        params = dict(task.get("params") or {})
        if params.get("review_contract") != REVIEW_TASK_CONTRACT:
            raise ValueError("scientific-review task uses an unsupported contract")
        campaign_id = task.get("campaign_id") or self.campaign_id
        if not campaign_id or params.get("campaign_id") != campaign_id:
            raise ValueError("review campaign differs from its frozen task parameters")
        if params.get("execution_profile") != self.profile_name:
            raise ValueError("review task execution profile differs from this worker")
        resolved = params.get("resolved_execution") or {}
        expected_execution = {
            "provider": self.profile.provider,
            "model": self.profile.model,
            "reasoning_effort": self.profile.reasoning_effort,
        }
        if resolved != expected_execution:
            raise ValueError("review task execution identity differs from the loaded policy")
        expected_policy = {
            "name": self.policy.name,
            "version": self.policy.version,
            "hash": self.policy.policy_hash,
        }
        if params.get("policy") != expected_policy:
            raise ValueError("review task policy differs from the loaded policy")
        if params.get("rubric_version") != self.policy.rubric_version:
            raise ValueError("review task rubric differs from the loaded policy")
        if params.get("scientific_contract") != SCIENTIFIC_CONTRACT:
            raise ValueError("review task scientific contract is incomplete")
        target = params.get("target") or {}
        if target.get("kind") not in TARGET_KINDS or not target.get("id"):
            raise ValueError("review task has an unsupported target")
        target_input = params.get("target_input")
        if not isinstance(target_input, dict) or target_input.get("id") != target.get("id"):
            raise ValueError("review task lacks a frozen target input")
        if not isinstance(params.get("review_tier"), str) or not params["review_tier"]:
            raise ValueError("review task lacks a review tier")
        for flag in ("blind_to_prior_scores", "blind_to_other_reviews"):
            if not isinstance(params.get(flag), bool):
                raise ValueError(f"review task {flag} must be boolean")
        source_ids = sorted(set(params.get("source_review_ids") or []))
        if params.get("source_review_manifest_hash") != _content_hash(source_ids):
            raise ValueError("review source manifest differs from its frozen hash")
        dataset_id = params.get("dataset_id")
        if not isinstance(dataset_id, str) or not dataset_id:
            raise ValueError("review task lacks a sealed dataset identity")
        return params

    def _verify_dataset(self, dataset_id: str) -> None:
        result = verify_dataset_seal(self.seal_path, db_path=self.db_path)
        if not result.valid:
            changed = ", ".join(result.changed_sections) or "dataset identity"
            raise ValueError(f"dataset seal verification failed; changed: {changed}")
        if result.expected_dataset_id != dataset_id:
            raise ValueError("review task dataset differs from the verified dataset seal")

    def _dataset_context(self) -> dict[str, Any]:
        if self._dataset_context_cache is None:
            self._dataset_context_cache = {
                "capabilities": _scientific_capabilities(self.db_path),
                "schema": _schema_catalog(self.db_path),
            }
        return self._dataset_context_cache

    @staticmethod
    def _validate_target_snapshot(
        record: dict[str, Any],
        target_input: dict[str, Any],
    ) -> None:
        drift = [
            key
            for key, expected in target_input.items()
            if record.get(key) != expected
        ]
        if drift:
            raise ValueError(
                "review target differs from its frozen input: " + ", ".join(drift)
            )

    @staticmethod
    def _validate_target_scope(
        record: dict[str, Any],
        *,
        campaign_id: str,
        dataset_id: str,
    ) -> None:
        if record.get("campaign_id") not in {None, campaign_id}:
            raise ValueError("review target crosses the frozen campaign")
        if record.get("dataset_id") not in {None, dataset_id}:
            raise ValueError("review target crosses the sealed dataset")

    def _cluster_summary(
        self,
        cluster_id: str,
        *,
        member_limit: int,
    ) -> tuple[dict[str, Any], list[dict[str, Any]]]:
        cluster = dict(
            self.ops.get_candidate_cluster(
                cluster_id,
                member_limit=max(1, member_limit),
            )
        )
        member_refs = list(cluster.pop("members", [])) if member_limit else []
        return cluster, member_refs[:member_limit]

    def _source_reviews(self, params: dict[str, Any]) -> list[dict[str, Any]]:
        if bool(params.get("blind_to_other_reviews")):
            return []
        result: list[dict[str, Any]] = []
        for review_id in sorted(set(params.get("source_review_ids") or [])):
            review = dict(self.ops.get_finding_review(review_id))
            review["verifications"] = self.ops.list_review_verifications(review_id)
            result.append(review)
        return result

    def _build_input(self, task: dict[str, Any], params: dict[str, Any]) -> dict[str, Any]:
        target = params["target"]
        kind = str(target["kind"])
        target_id = str(target["id"])
        campaign_id = str(task.get("campaign_id") or params["campaign_id"])
        member_slots: list[tuple[str, dict[str, Any]]] = []

        if kind == "candidate_cluster":
            cluster, refs = self._cluster_summary(
                target_id,
                member_limit=self.max_members,
            )
            self._validate_target_snapshot(cluster, params.get("target_input") or {})
            self._validate_target_scope(
                cluster,
                campaign_id=campaign_id,
                dataset_id=params["dataset_id"],
            )
            target_evidence: dict[str, Any] = {
                "record": cluster,
                "cluster_members": [],
                "member_sample": {
                    "available": int(cluster.get("member_count") or len(refs)),
                    "included": 0,
                },
            }
            for ref in refs:
                member_slots.append((target_id, ref))
        elif kind == "finding":
            finding = dict(self.ops.get_finding(target_id))
            self._validate_target_snapshot(finding, params.get("target_input") or {})
            self._validate_target_scope(
                finding,
                campaign_id=campaign_id,
                dataset_id=params["dataset_id"],
            )
            links = self.ops.list_cluster_findings(finding_id=target_id)
            target_evidence = {
                "record": finding,
                "cluster_links": links,
                "linked_clusters": [],
            }
            remaining = self.max_members
            for link in links:
                cluster, refs = self._cluster_summary(
                    str(link["cluster_id"]),
                    member_limit=remaining,
                )
                if cluster.get("campaign_id") != campaign_id:
                    raise ValueError("linked source cluster crosses the review campaign")
                if cluster.get("dataset_id") != params["dataset_id"]:
                    raise ValueError("linked source cluster crosses the review dataset")
                bundle = {
                    "link": link,
                    "record": cluster,
                    "cluster_members": [],
                    "member_sample": {
                        "available": int(cluster.get("member_count") or len(refs)),
                        "included": 0,
                    },
                }
                target_evidence["linked_clusters"].append(bundle)
                for ref in refs:
                    member_slots.append((str(link["cluster_id"]), ref))
                remaining = max(0, remaining - len(refs))
        else:
            disposition = dict(self.ops.get_unit_disposition(target_id))
            self._validate_target_snapshot(
                disposition,
                params.get("target_input") or {},
            )
            self._validate_target_scope(
                disposition,
                campaign_id=campaign_id,
                dataset_id=params["dataset_id"],
            )
            target_evidence = {"record": disposition}

        source_reviews = self._source_reviews(params)
        if bool(params.get("blind_to_prior_scores")):
            target_evidence = _strip_prior_assessments(target_evidence)
            source_reviews = _strip_prior_assessments(source_reviews)

        containers_by_cluster: dict[str, dict[str, Any]] = {}
        if kind == "candidate_cluster":
            containers_by_cluster[target_id] = target_evidence
        elif kind == "finding":
            containers_by_cluster = {
                str(bundle["record"]["id"]): bundle
                for bundle in target_evidence["linked_clusters"]
            }

        payload: dict[str, Any] = {
            "schema_version": INPUT_SCHEMA,
            "dataset_id": params["dataset_id"],
            "task_contract": {
                "review_tier": params["review_tier"],
                "rubric_version": params["rubric_version"],
                "blind_to_prior_scores": bool(params["blind_to_prior_scores"]),
                "blind_to_other_reviews": bool(params["blind_to_other_reviews"]),
                "scientific_contract": params["scientific_contract"],
                "audit": params.get("audit"),
            },
            "target": {"kind": kind, "id": target_id},
            "target_evidence": target_evidence,
            "source_reviews": source_reviews,
            "dataset_context": self._dataset_context(),
        }
        assert_sequence_free(payload, field="Scientific-review model input")
        if _json_bytes(payload) > self.max_input_bytes:
            raise ValueError(
                "scientific-review base input exceeds the configured byte bound"
            )

        for cluster_id, ref in member_slots:
            container = containers_by_cluster[cluster_id]
            occurrence = dict(
                self.ops.get_candidate_occurrence(str(ref["candidate_id"]))
            )
            if bool(params.get("blind_to_prior_scores")):
                occurrence = _strip_prior_assessments(occurrence)
            candidate = {"role": ref.get("role"), "occurrence": occurrence}
            members = container["cluster_members"]
            members.append(candidate)
            if _json_bytes(payload) > self.max_input_bytes:
                members.pop()
                continue
            container["member_sample"]["included"] += 1

        assert_sequence_free(payload, field="Scientific-review model input")
        return payload

    # -------------------------------------------------------------- model/checks

    def _call_model(self, payload: dict[str, Any], task_id: str) -> ModelRun:
        last: Exception | None = None
        payload_text = _canonical_json(payload)
        for attempt in range(self.transient_retries + 1):
            try:
                return run_profile(
                    provider=self.profile.provider,
                    model=self.profile.model,
                    reasoning_effort=self.profile.reasoning_effort,
                    system_prompt=REVIEW_SYSTEM_PROMPT,
                    payload_text=payload_text,
                    output_schema=REVIEW_OUTPUT_SCHEMA,
                    payload_label="SCIENTIFIC REVIEW INPUT (JSON)",
                    timeout=self.model_timeout,
                )
            except ModelTransient as exc:
                last = exc
                if attempt >= self.transient_retries:
                    break
                delay = min(2**attempt, 30)
                LOGGER.warning(
                    "transient model failure on task %s; retrying in %ss",
                    task_id,
                    delay,
                )
                time.sleep(delay)
        assert last is not None
        raise last

    def _run_checks(
        self,
        checks: tuple[PreparedCheck, ...],
        *,
        dataset_id: str,
    ) -> list[VerificationResult]:
        return [
            run_duckdb_verification(
                self.db_path,
                check.specification,
                check.expected,
                dataset_id=dataset_id,
                seal_path=self.seal_path,
                verify_seal=False,
                threads=self.verification_threads,
            )
            for check in checks
        ]

    @staticmethod
    def _verification_summary(
        prepared: PreparedReview,
        results: list[VerificationResult],
        run: ModelRun,
    ) -> dict[str, Any]:
        effective = _effective_verdict(prepared.verdict, results)
        return {
            "requested_verdict": prepared.verdict,
            "effective_verdict": effective,
            "all_checks_passed": all(result.status == "pass" for result in results),
            "checks": [
                {
                    "request": check.to_request_dict(),
                    "result": {
                        "status": result.status,
                        "specification_hash": result.specification_hash,
                        "query_hash": result.query_hash,
                        "row_count": result.row_count,
                        "error": result.error,
                    },
                }
                for check, result in zip(prepared.checks, results, strict=True)
            ],
            "model_execution": {
                "provider": run.provider,
                "model": run.model,
                "reasoning_effort": run.reasoning_effort,
                "prompt_sha256": run.prompt_sha256,
                "result_sha256": run.result_sha256,
                "exit_status": run.exit_status,
                "stderr_sha256": hashlib.sha256(run.stderr.encode("utf-8")).hexdigest(),
                "usage": run.usage,
                "raw_stdout_bytes": run.raw_stdout_bytes,
            },
        }

    def _target_kwargs(self, params: dict[str, Any]) -> dict[str, str]:
        kind = params["target"]["kind"]
        return {
            {
                "finding": "finding_id",
                "candidate_cluster": "cluster_id",
                "unit_disposition": "unit_disposition_id",
            }[kind]: str(params["target"]["id"])
        }

    def _existing_review(
        self,
        task_id: str,
        params: dict[str, Any],
    ) -> dict[str, Any] | None:
        rows = self.ops.list_finding_reviews(
            campaign_id=params["campaign_id"],
            review_tier=params["review_tier"],
            **self._target_kwargs(params),
        )
        matches = [row for row in rows if row.get("task_id") == task_id]
        if len(matches) > 1:
            raise ValueError("scientific-review task owns multiple review records")
        return matches[0] if matches else None

    @staticmethod
    def _checks_from_summary(review: dict[str, Any]) -> tuple[PreparedCheck, ...]:
        summary = review.get("verification_summary") or {}
        checks: list[PreparedCheck] = []
        for item in summary.get("checks") or []:
            request = item.get("request") or {}
            checks.append(
                PreparedCheck(
                    claim_key=str(request["claim_key"]),
                    specification=VerificationSpec.model_validate(
                        request["specification"]
                    ),
                    expected=request.get("expected"),
                )
            )
        if not checks:
            raise ValueError("task-owned review lacks executable verification requests")
        return tuple(checks)

    def _record_checks(
        self,
        review_id: str,
        task_id: str,
        checks: tuple[PreparedCheck, ...],
        results: list[VerificationResult],
        *,
        dataset_id: str,
    ) -> None:
        existing = self.ops.list_review_verifications(review_id)
        existing_keys = {
            (row.get("claim_key"), row.get("specification_hash"))
            for row in existing
        }
        for check, result in zip(checks, results, strict=True):
            specification = check.specification.model_dump(mode="json")
            spec_hash = _content_hash(specification)
            if (check.claim_key, spec_hash) in existing_keys:
                continue
            self.ops.record_review_verification(
                review_id,
                claim_key=check.claim_key,
                engine="duckdb",
                specification=specification,
                dataset_id=dataset_id,
                expected=check.expected,
                actual=result.actual,
                status=result.status,
                executed_ts=result.executed_ts,
                code_commit=self._code_commit,
                error=result.error,
                idempotency_key=(
                    "scientific-review-verification:"
                    + _content_hash(
                        {
                            "task_id": task_id,
                            "claim_key": check.claim_key,
                            "specification": specification,
                            "expected": check.expected,
                        }
                    )
                ),
            )

    def _missing_checks(
        self,
        review_id: str,
        checks: tuple[PreparedCheck, ...],
    ) -> tuple[PreparedCheck, ...]:
        existing = self.ops.list_review_verifications(review_id)
        existing_keys = {
            (row.get("claim_key"), row.get("specification_hash"))
            for row in existing
        }
        return tuple(
            check
            for check in checks
            if (
                check.claim_key,
                _content_hash(check.specification.model_dump(mode="json")),
            )
            not in existing_keys
        )

    # ---------------------------------------------------------------- execution

    def run_task(self, task: dict[str, Any]) -> None:
        task_id = str(task["id"])
        params = self._validate_contract(task)
        dataset_id = str(params["dataset_id"])

        with self._lease_keepalive(task_id):
            existing = self._existing_review(task_id, params)
            if existing is not None:
                checks = self._checks_from_summary(existing)
                missing = self._missing_checks(str(existing["id"]), checks)
                if missing:
                    self._verify_dataset(dataset_id)
                    results = self._run_checks(missing, dataset_id=dataset_id)
                else:
                    results = []
                self._record_checks(
                    str(existing["id"]),
                    task_id,
                    missing,
                    results,
                    dataset_id=dataset_id,
                )
                self.ops.complete_task(task_id)
                LOGGER.info("resumed and completed scientific-review task %s", task_id)
                return

            self._verify_dataset(dataset_id)
            payload = self._build_input(task, params)
            run = self._call_model(payload, task_id)
            prepared = _prepare_review(
                run.record,
                max_checks=self.policy.limits.max_evidence_tasks_per_review,
            )
            self._verify_dataset(dataset_id)
            results = self._run_checks(prepared.checks, dataset_id=dataset_id)
            verification_summary = self._verification_summary(prepared, results, run)
            verdict = str(verification_summary["effective_verdict"])
            review_id = self.ops.create_finding_review(
                campaign_id=params["campaign_id"],
                dataset_id=dataset_id,
                review_tier=params["review_tier"],
                execution_profile=params["execution_profile"],
                provider=self.profile.provider,
                model=self.profile.model,
                reasoning_effort=self.profile.reasoning_effort,
                prompt_hash=run.prompt_sha256,
                rubric_version=params["rubric_version"],
                input_bundle_hash=_content_hash(payload),
                verdict=verdict,
                confidence=prepared.confidence,
                task_id=task_id,
                reconstructed_observations=prepared.reconstructed_observations,
                claim_assessment=prepared.claim_assessment,
                verification_summary=verification_summary,
                discrepancies=prepared.discrepancies,
                proposed_tasks=prepared.proposed_tasks,
                blind_to_prior_scores=bool(params["blind_to_prior_scores"]),
                blind_to_other_reviews=bool(params["blind_to_other_reviews"]),
                idempotency_key=f"scientific-review:{task_id}",
                **self._target_kwargs(params),
            )
            self._record_checks(
                review_id,
                task_id,
                prepared.checks,
                results,
                dataset_id=dataset_id,
            )
            self.ops.complete_task(task_id)
            LOGGER.info(
                "scientific-review task %s complete: verdict=%s, checks=%s",
                task_id,
                verdict,
                len(prepared.checks),
            )

    @contextlib.contextmanager
    def _lease_keepalive(self, task_id: str):
        stop = threading.Event()

        def beat() -> None:
            while not stop.wait(self.lease_seconds / 3.0):
                try:
                    self.ops.heartbeat_task(
                        task_id,
                        lease_seconds=self.lease_seconds,
                    )
                except Exception:
                    LOGGER.debug("review keepalive failed for %s", task_id, exc_info=True)
                    return

        thread = threading.Thread(
            target=beat,
            name=f"review-keepalive-{task_id[:8]}",
            daemon=True,
        )
        thread.start()
        try:
            yield
        finally:
            stop.set()
            thread.join(timeout=5.0)


__all__ = [
    "DEFAULT_MAX_INPUT_BYTES",
    "INPUT_SCHEMA",
    "REVIEW_OUTPUT_SCHEMA",
    "REVIEW_SYSTEM_PROMPT",
    "REVIEW_TASK_CONTRACT",
    "SCIENTIFIC_CONTRACT",
    "PreparedCheck",
    "PreparedReview",
    "ScientificReviewWorker",
]
