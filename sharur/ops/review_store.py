"""Append-only scientific review-DAG operations for :mod:`sharur.ops.store`."""

from __future__ import annotations

import contextlib
import hashlib
import json
import time
import uuid
from typing import Any


MAX_REVIEW_JSON_BYTES = 256 * 1024
REVIEW_VERDICTS = frozenset(
    {"promote", "hold", "needs_data", "reject", "duplicate", "split"}
)
PROMOTION_DECISIONS = frozenset(
    {
        "promote",
        "hold",
        "needs_data",
        "reject",
        "duplicate",
        "split",
        "merge",
        "publish",
        "reopen",
    }
)
_JSON_COLUMNS = frozenset(
    {
        "reason_codes",
        "strata",
        "provenance",
        "signature",
        "evidence",
        "verification",
        "uncertainty",
        "subject_refs",
        "reduction_features",
        "summary",
        "counts",
        "reconstructed_observations",
        "claim_assessment",
        "verification_summary",
        "discrepancies",
        "proposed_tasks",
        "specification",
        "expected",
        "actual",
        "audit_stratum",
        "metadata",
    }
)
_RAW_SEQUENCE_KEYS = frozenset(
    {
        "sequence",
        "seq",
        "aa_sequence",
        "nt_sequence",
        "protein_sequence",
        "amino_acid_sequence",
        "nucleotide_sequence",
    }
)


def contains_raw_biological_sequence(value: Any) -> bool:
    """Detect model-visible sequence fields and sequence-like strings."""

    if isinstance(value, str):
        compact = "".join(value.split())
        if len(compact) < 40:
            return False
        letters = set(compact.upper())
        nucleotide = set("ACGTUNRYKMSWBDHV-")
        amino_acid = set("ABCDEFGHIKLMNPQRSTVWXYZ*-")
        return bool(compact) and (
            letters <= nucleotide or letters <= amino_acid
        )
    if isinstance(value, (list, tuple)):
        return any(contains_raw_biological_sequence(item) for item in value)
    if isinstance(value, dict):
        if any(str(key).strip().lower() in _RAW_SEQUENCE_KEYS for key in value):
            return True
        return any(contains_raw_biological_sequence(item) for item in value.values())
    return False


def assert_sequence_free(value: Any, *, field: str) -> None:
    if contains_raw_biological_sequence(value):
        raise ValueError(f"{field} contains model-visible raw biological sequence")


def canonical_json(value: Any, *, field: str) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)
    if len(encoded.encode("utf-8")) > MAX_REVIEW_JSON_BYTES:
        raise ValueError(
            f"{field} exceeds the {MAX_REVIEW_JSON_BYTES}-byte inline JSON limit"
        )
    return encoded


def content_hash(value: Any) -> str:
    return hashlib.sha256(canonical_json(value, field="hash input").encode()).hexdigest()


def decode_review_row(row: Any) -> dict[str, Any] | None:
    if row is None:
        return None
    result = dict(row)
    for field in _JSON_COLUMNS:
        value = result.get(field)
        if isinstance(value, str):
            with contextlib.suppress(json.JSONDecodeError, TypeError):
                result[field] = json.loads(value)
    for field in ("blind_to_prior_scores", "blind_to_other_reviews", "audit_sample"):
        if field in result:
            result[field] = bool(result[field])
    return result


def _nonempty(value: str, field: str) -> str:
    normalized = str(value).strip()
    if not normalized:
        raise ValueError(f"{field} must be non-empty")
    return normalized


def _string_list(values: list[str] | None, field: str) -> list[str]:
    result: list[str] = []
    for value in values or []:
        normalized = _nonempty(value, field)
        if normalized not in result:
            result.append(normalized)
    return result


def _idempotent_id(
    existing: Any,
    immutable: dict[str, Any],
    *,
    record_type: str,
) -> str | None:
    if existing is None:
        return None
    conflicts = [
        field for field, expected in immutable.items() if existing[field] != expected
    ]
    if conflicts:
        raise ValueError(
            f"{record_type} idempotency key already exists with a different "
            f"payload: {', '.join(conflicts)}"
        )
    return str(existing["id"])


class ReviewStoreMixin:
    """Review methods mixed into :class:`~sharur.ops.store.OpsStore`."""

    def _review_output_task_locked(
        self,
        task_id: str,
        *,
        campaign_id: str,
        dataset_id: str,
        unit_id: str,
        genome_id: str,
    ) -> Any:
        task = self._conn.execute(
            """
            SELECT campaign_id, assigned_to, status, lease_expires_ts, params
            FROM tasks WHERE id = ?
            """,
            (task_id,),
        ).fetchone()
        if task is None:
            raise ValueError(f"Task {task_id} not found")
        if task["campaign_id"] != campaign_id:
            raise ValueError("Review-output task belongs to a different campaign")
        params = json.loads(task["params"])
        if params.get("review_output_contract") is not None:
            expected = {
                "dataset_id": dataset_id,
                "unit_id": unit_id,
                "genome_id": genome_id,
            }
            mismatches = [
                field for field, value in expected.items() if params.get(field) != value
            ]
            if mismatches:
                raise ValueError(
                    "Review output differs from its task target: "
                    + ", ".join(mismatches)
                )
        return task

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
        strata: dict[str, Any] | None = None,
        provenance: dict[str, Any] | None = None,
        supersedes_disposition_id: str | None = None,
        idempotency_key: str,
        schema_version: int = 1,
    ) -> str:
        if disposition not in {"clear", "candidate", "incomplete", "failed"}:
            raise ValueError(f"Unsupported disposition: {disposition}")
        if isinstance(candidate_count, bool) or candidate_count < 0:
            raise ValueError("candidate_count must be a non-negative integer")
        unit_id = _nonempty(unit_id, "unit_id")
        dataset_id = _nonempty(dataset_id, "dataset_id")
        genome_id = _nonempty(genome_id, "genome_id")
        coverage_hash = _nonempty(coverage_hash, "coverage_hash")
        evidence_bundle_hash = _nonempty(
            evidence_bundle_hash, "evidence_bundle_hash"
        )
        values = {
            "reason_codes": canonical_json(
                _string_list(reason_codes, "reason code"), field="reason_codes"
            ),
            "strata": canonical_json(strata or {}, field="strata"),
            "provenance": canonical_json(provenance or {}, field="provenance"),
        }
        assert_sequence_free(
            {
                "reason_codes": reason_codes or [],
                "strata": strata or {},
                "provenance": provenance or {},
            },
            field="Unit disposition",
        )
        with self._lock, self._transaction():
            self._validate_campaign_locked(campaign_id)
            task = (
                self._review_output_task_locked(
                    task_id,
                    campaign_id=campaign_id,
                    dataset_id=dataset_id,
                    unit_id=unit_id,
                    genome_id=genome_id,
                )
                if task_id is not None
                else None
            )
            if task_id is not None:
                existing = self._conn.execute(
                    """
                    SELECT * FROM unit_dispositions
                    WHERE task_id = ? AND idempotency_key = ?
                    """,
                    (task_id, idempotency_key),
                ).fetchone()
            else:
                existing = self._conn.execute(
                    """
                    SELECT * FROM unit_dispositions
                    WHERE task_id IS NULL
                      AND agent_id = ? AND idempotency_key = ?
                    """,
                    (self.agent_id, idempotency_key),
                ).fetchone()
            existing_id = _idempotent_id(
                existing,
                {
                    "campaign_id": campaign_id,
                    "task_id": task_id,
                    "unit_id": unit_id,
                    "dataset_id": dataset_id,
                    "genome_id": genome_id,
                    "coverage_hash": coverage_hash,
                    "candidate_count": candidate_count,
                    "disposition": disposition,
                    "reason_codes": values["reason_codes"],
                    "strata": values["strata"],
                    "evidence_bundle_hash": evidence_bundle_hash,
                    "provenance": values["provenance"],
                    "supersedes_disposition_id": supersedes_disposition_id,
                    "schema_version": schema_version,
                },
                record_type="Unit-disposition",
            )
            if existing_id is not None:
                return existing_id
            version = 1
            if supersedes_disposition_id is not None:
                previous = self._conn.execute(
                    "SELECT * FROM unit_dispositions WHERE id = ?",
                    (supersedes_disposition_id,),
                ).fetchone()
                if previous is None:
                    raise ValueError(
                        f"Unit disposition {supersedes_disposition_id} not found"
                    )
                if (
                    previous["campaign_id"] != campaign_id
                    or previous["unit_id"] != unit_id
                    or previous["dataset_id"] != dataset_id
                    or previous["genome_id"] != genome_id
                    or previous["record_status"] != "active"
                ):
                    raise ValueError(
                        "Superseded disposition is not the active version of "
                        "the same unit"
                    )
                version = int(previous["version"]) + 1
            else:
                active = self._conn.execute(
                    """
                    SELECT id FROM unit_dispositions
                    WHERE campaign_id = ? AND unit_id = ? AND record_status = 'active'
                    """,
                    (campaign_id, unit_id),
                ).fetchone()
                if active is not None:
                    raise ValueError(
                        "Unit already has an active disposition; supply "
                        "supersedes_disposition_id"
                    )
            if task is not None and (
                task["assigned_to"] != self.agent_id
                or task["status"] not in {"claimed", "in_progress"}
                or float(task["lease_expires_ts"] or 0) <= time.time()
            ):
                raise ValueError(
                    "Disposition task has no active lease owned by the "
                    "submitting agent"
                )
            candidate_sql = """
                SELECT COUNT(*) FROM candidate_occurrences
                WHERE campaign_id = ? AND dataset_id = ?
                  AND unit_id = ? AND genome_id = ?
            """
            candidate_params: list[Any] = [
                campaign_id,
                dataset_id,
                unit_id,
                genome_id,
            ]
            if task_id is not None:
                candidate_sql += " AND task_id = ?"
                candidate_params.append(task_id)
            observed = self._conn.execute(
                candidate_sql,
                candidate_params,
            ).fetchone()[0]
            if int(observed) != int(candidate_count):
                raise ValueError(
                    f"candidate_count {candidate_count} differs from recorded "
                    f"unit candidates {observed}"
                )
            record_id = str(uuid.uuid4())
            self._conn.execute(
                """
                INSERT INTO unit_dispositions(
                    id, campaign_id, task_id, unit_id, dataset_id, genome_id,
                    agent_id, ts, coverage_hash, candidate_count, disposition,
                    reason_codes, strata, evidence_bundle_hash, provenance,
                    idempotency_key, version, supersedes_disposition_id,
                    record_status, schema_version
                ) VALUES (
                    ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?
                )
                """,
                (
                    record_id,
                    campaign_id,
                    task_id,
                    unit_id,
                    dataset_id,
                    genome_id,
                    self.agent_id,
                    time.time(),
                    coverage_hash,
                    candidate_count,
                    disposition,
                    values["reason_codes"],
                    values["strata"],
                    evidence_bundle_hash,
                    values["provenance"],
                    idempotency_key,
                    version,
                    supersedes_disposition_id,
                    "active",
                    schema_version,
                ),
            )
            if supersedes_disposition_id is not None:
                self._conn.execute(
                    """
                    UPDATE unit_dispositions SET record_status = 'superseded'
                    WHERE id = ?
                    """,
                    (supersedes_disposition_id,),
                )
            event_id = self._event_locked(
                "unit_disposition_created",
                "unit_disposition",
                record_id,
                campaign_id=campaign_id,
                task_id=task_id,
                payload={
                    "unit_id": unit_id,
                    "disposition": disposition,
                    "candidate_count": candidate_count,
                    "version": version,
                    "supersedes_disposition_id": supersedes_disposition_id,
                },
            )
            self._conn.execute(
                "UPDATE unit_dispositions SET created_event_id = ? WHERE id = ?",
                (event_id, record_id),
            )
            return record_id

    def list_unit_dispositions(
        self,
        *,
        campaign_id: str,
        disposition: str | None = None,
        active_only: bool = True,
        limit: int = 1_000,
    ) -> list[dict[str, Any]]:
        sql = "SELECT * FROM unit_dispositions WHERE campaign_id = ?"
        params: list[Any] = [campaign_id]
        if disposition is not None:
            sql += " AND disposition = ?"
            params.append(disposition)
        if active_only:
            sql += " AND record_status = 'active'"
        sql += " ORDER BY ts, id LIMIT ?"
        params.append(limit)
        with self._lock:
            rows = self._conn.execute(sql, params).fetchall()
        return [decode_review_row(row) for row in rows]

    def get_unit_disposition(self, disposition_id: str) -> dict[str, Any]:
        with self._lock:
            row = self._conn.execute(
                "SELECT * FROM unit_dispositions WHERE id = ?",
                (disposition_id,),
            ).fetchone()
        result = decode_review_row(row)
        if result is None:
            raise KeyError(f"Unit disposition {disposition_id} not found")
        return result

    def create_candidate_occurrence(
        self,
        *,
        campaign_id: str,
        dataset_id: str,
        unit_id: str,
        genome_id: str,
        candidate_type: str,
        signature_schema: str,
        signature: dict[str, Any],
        evidence: dict[str, Any],
        verification: list[dict[str, Any]],
        subject_refs: dict[str, Any],
        task_id: str | None = None,
        reason_codes: list[str] | None = None,
        uncertainty: dict[str, Any] | None = None,
        reduction_features: dict[str, Any] | None = None,
        provenance: dict[str, Any] | None = None,
        evidence_bundle_hash: str | None = None,
        idempotency_key: str,
        schema_version: int = 1,
    ) -> str:
        candidate_type = _nonempty(candidate_type, "candidate_type")
        signature_schema = _nonempty(signature_schema, "signature_schema")
        dataset_id = _nonempty(dataset_id, "dataset_id")
        unit_id = _nonempty(unit_id, "unit_id")
        genome_id = _nonempty(genome_id, "genome_id")
        assert_sequence_free(
            {
                "signature": signature,
                "evidence": evidence,
                "verification": verification,
                "subject_refs": subject_refs,
                "reason_codes": reason_codes or [],
                "uncertainty": uncertainty or {},
                "reduction_features": reduction_features or {},
                "provenance": provenance or {},
            },
            field="Candidate occurrence",
        )
        signature_json = canonical_json(signature, field="signature")
        signature_hash = hashlib.sha256(signature_json.encode()).hexdigest()
        fields = {
            "evidence": canonical_json(evidence, field="evidence"),
            "verification": canonical_json(verification, field="verification"),
            "subject_refs": canonical_json(subject_refs, field="subject_refs"),
            "reason_codes": canonical_json(
                _string_list(reason_codes, "reason code"), field="reason_codes"
            ),
            "uncertainty": canonical_json(uncertainty or {}, field="uncertainty"),
            "reduction_features": canonical_json(
                reduction_features or {}, field="reduction_features"
            ),
            "provenance": canonical_json(provenance or {}, field="provenance"),
        }
        bundle_hash = evidence_bundle_hash or content_hash(
            {
                "dataset_id": dataset_id,
                "unit_id": unit_id,
                "genome_id": genome_id,
                "candidate_type": candidate_type,
                "signature_schema": signature_schema,
                "signature": signature,
                "evidence": evidence,
                "verification": verification,
                "subject_refs": subject_refs,
                "reason_codes": reason_codes or [],
                "uncertainty": uncertainty or {},
                "reduction_features": reduction_features or {},
                "provenance": provenance or {},
            }
        )
        with self._lock, self._transaction():
            self._validate_campaign_locked(campaign_id)
            task = (
                self._review_output_task_locked(
                    task_id,
                    campaign_id=campaign_id,
                    dataset_id=dataset_id,
                    unit_id=unit_id,
                    genome_id=genome_id,
                )
                if task_id is not None
                else None
            )
            if task_id is not None:
                existing = self._conn.execute(
                    """
                    SELECT * FROM candidate_occurrences
                    WHERE task_id = ? AND idempotency_key = ?
                    """,
                    (task_id, idempotency_key),
                ).fetchone()
            else:
                existing = self._conn.execute(
                    """
                    SELECT * FROM candidate_occurrences
                    WHERE task_id IS NULL
                      AND agent_id = ? AND idempotency_key = ?
                    """,
                    (self.agent_id, idempotency_key),
                ).fetchone()
            existing_id = _idempotent_id(
                existing,
                {
                    "campaign_id": campaign_id,
                    "task_id": task_id,
                    "dataset_id": dataset_id,
                    "unit_id": unit_id,
                    "genome_id": genome_id,
                    "candidate_type": candidate_type,
                    "signature_schema": signature_schema,
                    "signature": signature_json,
                    "signature_hash": signature_hash,
                    "evidence": fields["evidence"],
                    "verification": fields["verification"],
                    "reason_codes": fields["reason_codes"],
                    "uncertainty": fields["uncertainty"],
                    "subject_refs": fields["subject_refs"],
                    "reduction_features": fields["reduction_features"],
                    "evidence_bundle_hash": bundle_hash,
                    "provenance": fields["provenance"],
                    "schema_version": schema_version,
                },
                record_type="Candidate-occurrence",
            )
            if existing_id is not None:
                return existing_id
            if task is not None and (
                task["assigned_to"] != self.agent_id
                or task["status"] not in {"claimed", "in_progress"}
                or float(task["lease_expires_ts"] or 0) <= time.time()
            ):
                raise ValueError(
                    "Candidate task has no active lease owned by the "
                    "submitting agent"
                )
            candidate_id = str(uuid.uuid4())
            self._conn.execute(
                """
                INSERT INTO candidate_occurrences(
                    id, campaign_id, task_id, agent_id, ts, dataset_id, unit_id,
                    genome_id,
                    candidate_type, signature_schema, signature, signature_hash,
                    evidence, verification, reason_codes, uncertainty, subject_refs,
                    reduction_features, evidence_bundle_hash, provenance,
                    idempotency_key, schema_version
                ) VALUES (
                    ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?
                )
                """,
                (
                    candidate_id,
                    campaign_id,
                    task_id,
                    self.agent_id,
                    time.time(),
                    dataset_id,
                    unit_id,
                    genome_id,
                    candidate_type,
                    signature_schema,
                    signature_json,
                    signature_hash,
                    fields["evidence"],
                    fields["verification"],
                    fields["reason_codes"],
                    fields["uncertainty"],
                    fields["subject_refs"],
                    fields["reduction_features"],
                    bundle_hash,
                    fields["provenance"],
                    idempotency_key,
                    schema_version,
                ),
            )
            event_id = self._event_locked(
                "candidate_occurrence_created",
                "candidate_occurrence",
                candidate_id,
                campaign_id=campaign_id,
                task_id=task_id,
                payload={
                    "candidate_type": candidate_type,
                    "signature_schema": signature_schema,
                    "signature_hash": signature_hash,
                },
            )
            self._conn.execute(
                "UPDATE candidate_occurrences SET created_event_id = ? WHERE id = ?",
                (event_id, candidate_id),
            )
            return candidate_id

    def get_candidate_occurrence(self, candidate_id: str) -> dict[str, Any]:
        with self._lock:
            row = self._conn.execute(
                "SELECT * FROM candidate_occurrences WHERE id = ?", (candidate_id,)
            ).fetchone()
        result = decode_review_row(row)
        if result is None:
            raise KeyError(f"Candidate {candidate_id} not found")
        return result

    def list_candidate_occurrences(
        self,
        *,
        campaign_id: str,
        candidate_type: str | None = None,
        task_id: str | None = None,
        unclustered_only: bool = False,
        limit: int = 1_000,
    ) -> list[dict[str, Any]]:
        sql = "SELECT c.* FROM candidate_occurrences AS c"
        params: list[Any] = []
        if unclustered_only:
            sql += " LEFT JOIN candidate_cluster_members AS m ON m.candidate_id = c.id"
        sql += " WHERE c.campaign_id = ?"
        params.append(campaign_id)
        if candidate_type is not None:
            sql += " AND c.candidate_type = ?"
            params.append(candidate_type)
        if task_id is not None:
            # Lets a worker count exactly the rows it has persisted for its own
            # task, which is what makes finalize resume-safe: the in-memory
            # candidate list resets to empty on resume, the store does not.
            sql += " AND c.task_id = ?"
            params.append(task_id)
        if unclustered_only:
            sql += " AND m.candidate_id IS NULL"
        sql += " ORDER BY c.ts, c.id LIMIT ?"
        params.append(limit)
        with self._lock:
            rows = self._conn.execute(sql, params).fetchall()
        return [decode_review_row(row) for row in rows]

    def create_candidate_cluster(
        self,
        *,
        campaign_id: str,
        dataset_id: str,
        candidate_type: str,
        signature_schema: str,
        member_ids: list[str],
        reducer_name: str,
        reducer_version: str,
        reducer_config_hash: str,
        summary: dict[str, Any],
        counts: dict[str, Any],
        roles: dict[str, str] | None = None,
        logical_cluster_id: str | None = None,
        version: int = 1,
        idempotency_key: str,
        schema_version: int = 1,
    ) -> str:
        members = sorted(set(member_ids))
        if not members:
            raise ValueError("Candidate cluster requires at least one member")
        dataset_id = _nonempty(dataset_id, "dataset_id")
        candidate_type = _nonempty(candidate_type, "candidate_type")
        signature_schema = _nonempty(signature_schema, "signature_schema")
        reducer_name = _nonempty(reducer_name, "reducer_name")
        reducer_version = _nonempty(reducer_version, "reducer_version")
        reducer_config_hash = _nonempty(
            reducer_config_hash, "reducer_config_hash"
        )
        if version < 1:
            raise ValueError("Candidate-cluster version must be positive")
        manifest_hash = content_hash(members)
        logical_id = logical_cluster_id or (
            "candidate-cluster-"
            + content_hash(
                {
                    "campaign_id": campaign_id,
                    "candidate_type": candidate_type,
                    "signature_schema": signature_schema,
                    "members": members,
                }
            )[:24]
        )
        with self._lock, self._transaction():
            self._validate_campaign_locked(campaign_id)
            self._validate_ids_locked(
                "candidate_occurrences", members, label="candidate occurrences"
            )
            placeholders = ",".join("?" for _ in members)
            rows = self._conn.execute(
                f"""
                SELECT campaign_id, dataset_id, candidate_type, signature_schema
                FROM candidate_occurrences WHERE id IN ({placeholders})
                """,
                members,
            ).fetchall()
            observed = {
                (
                    row["campaign_id"],
                    row["dataset_id"],
                    row["candidate_type"],
                    row["signature_schema"],
                )
                for row in rows
            }
            if observed != {
                (campaign_id, dataset_id, candidate_type, signature_schema)
            }:
                raise ValueError("Cluster members cross a typed evidence boundary")
            existing = self._conn.execute(
                """
                SELECT * FROM candidate_clusters
                WHERE created_by = ? AND idempotency_key = ?
                """,
                (self.agent_id, idempotency_key),
            ).fetchone()
            summary_json = canonical_json(summary, field="cluster summary")
            counts_json = canonical_json(counts, field="cluster counts")
            assert_sequence_free(
                {
                    "summary": summary,
                    "counts": counts,
                    "roles": roles or {},
                },
                field="Candidate cluster",
            )
            existing_id = _idempotent_id(
                existing,
                {
                    "logical_cluster_id": logical_id,
                    "version": version,
                    "campaign_id": campaign_id,
                    "dataset_id": dataset_id,
                    "candidate_type": candidate_type,
                    "signature_schema": signature_schema,
                    "reducer_name": reducer_name,
                    "reducer_version": reducer_version,
                    "reducer_config_hash": reducer_config_hash,
                    "member_manifest_hash": manifest_hash,
                    "summary": summary_json,
                    "counts": counts_json,
                    "schema_version": schema_version,
                },
                record_type="Candidate-cluster",
            )
            if existing_id is not None:
                return existing_id
            parent_cluster_id: str | None = None
            if version > 1:
                parent = self._conn.execute(
                    """
                    SELECT id, campaign_id, status
                    FROM candidate_clusters
                    WHERE logical_cluster_id = ? AND version = ?
                    """,
                    (logical_id, version - 1),
                ).fetchone()
                if (
                    parent is None
                    or parent["campaign_id"] != campaign_id
                    or parent["status"] != "active"
                ):
                    raise ValueError(
                        "Candidate-cluster revision requires its active prior version"
                    )
                parent_cluster_id = str(parent["id"])
            cluster_id = str(uuid.uuid4())
            self._conn.execute(
                """
                INSERT INTO candidate_clusters(
                    id, logical_cluster_id, version, campaign_id, dataset_id,
                    candidate_type, signature_schema, reducer_name, reducer_version,
                    reducer_config_hash, member_manifest_hash, summary, counts,
                    created_by, created_ts, idempotency_key, schema_version
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    cluster_id,
                    logical_id,
                    version,
                    campaign_id,
                    dataset_id,
                    candidate_type,
                    signature_schema,
                    reducer_name,
                    reducer_version,
                    reducer_config_hash,
                    manifest_hash,
                    summary_json,
                    counts_json,
                    self.agent_id,
                    time.time(),
                    idempotency_key,
                    schema_version,
                ),
            )
            normalized_roles = dict(roles or {})
            if not normalized_roles:
                normalized_roles[members[0]] = "medoid"
            for candidate_id in members:
                role = normalized_roles.get(candidate_id, "member")
                if role not in {"member", "medoid", "outlier", "counterexample"}:
                    raise ValueError(f"Unsupported candidate role: {role}")
                self._conn.execute(
                    """
                    INSERT INTO candidate_cluster_members(cluster_id, candidate_id, role)
                    VALUES (?, ?, ?)
                    """,
                    (cluster_id, candidate_id, role),
                )
            event_id = self._event_locked(
                "candidate_cluster_created",
                "candidate_cluster",
                cluster_id,
                campaign_id=campaign_id,
                payload={
                    "candidate_type": candidate_type,
                    "member_count": len(members),
                    "member_manifest_hash": manifest_hash,
                },
            )
            self._conn.execute(
                "UPDATE candidate_clusters SET created_event_id = ? WHERE id = ?",
                (event_id, cluster_id),
            )
            if parent_cluster_id is not None:
                self._conn.execute(
                    """
                    UPDATE candidate_clusters SET status = 'superseded'
                    WHERE id = ?
                    """,
                    (parent_cluster_id,),
                )
                self._conn.execute(
                    """
                    INSERT INTO candidate_cluster_lineage(
                        parent_cluster_id, child_cluster_id, relation,
                        created_ts, created_by
                    ) VALUES (?, ?, 'supersedes', ?, ?)
                    """,
                    (
                        parent_cluster_id,
                        cluster_id,
                        time.time(),
                        self.agent_id,
                    ),
                )
                self._event_locked(
                    "candidate_cluster_superseded",
                    "candidate_cluster",
                    cluster_id,
                    campaign_id=campaign_id,
                    payload={"parent_cluster_id": parent_cluster_id},
                )
            return cluster_id

    def get_candidate_cluster(
        self,
        cluster_id: str,
        *,
        member_limit: int = 100,
    ) -> dict[str, Any]:
        if isinstance(member_limit, bool) or not 1 <= member_limit <= 1_000:
            raise ValueError("member_limit must be between 1 and 1000")
        with self._lock:
            row = self._conn.execute(
                "SELECT * FROM candidate_clusters WHERE id = ?", (cluster_id,)
            ).fetchone()
            if row is None:
                raise KeyError(f"Candidate cluster {cluster_id} not found")
            role_rows = self._conn.execute(
                """
                SELECT role, COUNT(*) AS count
                FROM candidate_cluster_members
                WHERE cluster_id = ?
                GROUP BY role ORDER BY role
                """,
                (cluster_id,),
            ).fetchall()
            role_counts = {
                str(role_row["role"]): int(role_row["count"])
                for role_row in role_rows
            }
            member_count = sum(role_counts.values())
            members = self._conn.execute(
                """
                SELECT candidate_id, role FROM candidate_cluster_members
                WHERE cluster_id = ? ORDER BY candidate_id
                LIMIT ?
                """,
                (cluster_id, member_limit),
            ).fetchall()
        result = decode_review_row(row)
        assert result is not None
        result["members"] = [dict(member) for member in members]
        result["member_count"] = member_count
        result["member_role_counts"] = role_counts
        result["members_truncated"] = member_count > len(members)
        return result

    def list_candidate_cluster_members(
        self,
        cluster_id: str,
        *,
        after_candidate_id: str | None = None,
        limit: int = 500,
    ) -> dict[str, Any]:
        if isinstance(limit, bool) or not 1 <= limit <= 1_000:
            raise ValueError("limit must be between 1 and 1000")
        with self._lock:
            if (
                self._conn.execute(
                    "SELECT 1 FROM candidate_clusters WHERE id = ?",
                    (cluster_id,),
                ).fetchone()
                is None
            ):
                raise KeyError(f"Candidate cluster {cluster_id} not found")
            sql = """
                SELECT candidate_id, role FROM candidate_cluster_members
                WHERE cluster_id = ?
            """
            params: list[Any] = [cluster_id]
            if after_candidate_id is not None:
                sql += " AND candidate_id > ?"
                params.append(after_candidate_id)
            sql += " ORDER BY candidate_id LIMIT ?"
            params.append(limit + 1)
            rows = self._conn.execute(sql, params).fetchall()
        has_more = len(rows) > limit
        page = rows[:limit]
        return {
            "cluster_id": cluster_id,
            "members": [dict(row) for row in page],
            "next_after_candidate_id": (
                str(page[-1]["candidate_id"]) if has_more and page else None
            ),
            "has_more": has_more,
        }

    def list_candidate_clusters(
        self,
        *,
        campaign_id: str,
        candidate_type: str | None = None,
        status: str | None = "active",
        limit: int = 1_000,
    ) -> list[dict[str, Any]]:
        sql = "SELECT * FROM candidate_clusters WHERE campaign_id = ?"
        params: list[Any] = [campaign_id]
        if candidate_type is not None:
            sql += " AND candidate_type = ?"
            params.append(candidate_type)
        if status is not None:
            sql += " AND status = ?"
            params.append(status)
        sql += " ORDER BY created_ts, id LIMIT ?"
        params.append(limit)
        with self._lock:
            rows = self._conn.execute(sql, params).fetchall()
        return [decode_review_row(row) for row in rows]

    def link_candidate_clusters(
        self, parent_cluster_id: str, child_cluster_id: str, *, relation: str
    ) -> None:
        if relation not in {"supersedes", "split_from", "merged_from", "refines"}:
            raise ValueError(f"Unsupported cluster-lineage relation: {relation}")
        with self._lock, self._transaction():
            self._validate_ids_locked(
                "candidate_clusters",
                [parent_cluster_id, child_cluster_id],
                label="candidate clusters",
            )
            cursor = self._conn.execute(
                """
                INSERT OR IGNORE INTO candidate_cluster_lineage(
                    parent_cluster_id, child_cluster_id, relation, created_ts, created_by
                ) VALUES (?, ?, ?, ?, ?)
                """,
                (
                    parent_cluster_id,
                    child_cluster_id,
                    relation,
                    time.time(),
                    self.agent_id,
                ),
            )
            if cursor.rowcount:
                self._event_locked(
                    "candidate_cluster_lineage_created",
                    "candidate_cluster",
                    child_cluster_id,
                    campaign_id=self._conn.execute(
                        "SELECT campaign_id FROM candidate_clusters WHERE id = ?",
                        (child_cluster_id,),
                    ).fetchone()["campaign_id"],
                    payload={
                        "parent_cluster_id": parent_cluster_id,
                        "relation": relation,
                    },
                )

    def supersede_candidate_cluster(
        self,
        parent_cluster_id: str,
        child_cluster_id: str,
    ) -> None:
        """Make a new immutable cluster version active and retire its parent."""

        with self._lock, self._transaction():
            self._validate_ids_locked(
                "candidate_clusters",
                [parent_cluster_id, child_cluster_id],
                label="candidate clusters",
            )
            rows = self._conn.execute(
                """
                SELECT id, logical_cluster_id, version, campaign_id
                FROM candidate_clusters WHERE id IN (?, ?)
                """,
                (parent_cluster_id, child_cluster_id),
            ).fetchall()
            by_id = {row["id"]: row for row in rows}
            parent = by_id[parent_cluster_id]
            child = by_id[child_cluster_id]
            if (
                parent["logical_cluster_id"] != child["logical_cluster_id"]
                or int(child["version"]) <= int(parent["version"])
                or parent["campaign_id"] != child["campaign_id"]
            ):
                raise ValueError("Invalid candidate-cluster supersession")
            self._conn.execute(
                """
                UPDATE candidate_clusters SET status = 'superseded'
                WHERE id = ? AND status = 'active'
                """,
                (parent_cluster_id,),
            )
            cursor = self._conn.execute(
                """
                INSERT OR IGNORE INTO candidate_cluster_lineage(
                    parent_cluster_id, child_cluster_id, relation, created_ts, created_by
                ) VALUES (?, ?, 'supersedes', ?, ?)
                """,
                (parent_cluster_id, child_cluster_id, time.time(), self.agent_id),
            )
            if cursor.rowcount:
                self._event_locked(
                    "candidate_cluster_superseded",
                    "candidate_cluster",
                    child_cluster_id,
                    campaign_id=child["campaign_id"],
                    payload={"parent_cluster_id": parent_cluster_id},
                )

    def link_cluster_finding(
        self,
        cluster_id: str,
        finding_id: str,
        *,
        relation: str = "materializes",
    ) -> None:
        if relation not in {"materializes", "supports", "counterexample"}:
            raise ValueError(f"Unsupported cluster-finding relation: {relation}")
        with self._lock, self._transaction():
            self._validate_ids_locked(
                "candidate_clusters", [cluster_id], label="candidate clusters"
            )
            self._validate_ids_locked("findings", [finding_id], label="findings")
            cluster = self._conn.execute(
                "SELECT campaign_id, status FROM candidate_clusters WHERE id = ?",
                (cluster_id,),
            ).fetchone()
            finding = self._conn.execute(
                "SELECT campaign_id FROM findings WHERE id = ?", (finding_id,)
            ).fetchone()
            if cluster["campaign_id"] != finding["campaign_id"]:
                raise ValueError("Cluster and finding belong to different campaigns")
            if relation == "materializes" and cluster["status"] != "active":
                raise ValueError(
                    "Only an active candidate-cluster version can materialize a finding"
                )
            cursor = self._conn.execute(
                """
                INSERT OR IGNORE INTO cluster_findings(
                    cluster_id, finding_id, relation, created_ts, created_by
                ) VALUES (?, ?, ?, ?, ?)
                """,
                (cluster_id, finding_id, relation, time.time(), self.agent_id),
            )
            if cursor.rowcount:
                self._event_locked(
                    "cluster_finding_linked",
                    "finding",
                    finding_id,
                    campaign_id=cluster["campaign_id"],
                    payload={"cluster_id": cluster_id, "relation": relation},
                )

    def list_cluster_findings(
        self,
        *,
        cluster_id: str | None = None,
        finding_id: str | None = None,
    ) -> list[dict[str, Any]]:
        if cluster_id is None and finding_id is None:
            raise ValueError("cluster_id or finding_id is required")
        clauses: list[str] = []
        params: list[Any] = []
        if cluster_id is not None:
            clauses.append("cluster_id = ?")
            params.append(cluster_id)
        if finding_id is not None:
            clauses.append("finding_id = ?")
            params.append(finding_id)
        with self._lock:
            rows = self._conn.execute(
                f"""
                SELECT * FROM cluster_findings
                WHERE {' AND '.join(clauses)}
                ORDER BY created_ts, cluster_id, finding_id, relation
                """,
                params,
            ).fetchall()
        return [dict(row) for row in rows]

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
        reconstructed_observations: dict[str, Any] | None = None,
        claim_assessment: dict[str, Any] | None = None,
        verification_summary: dict[str, Any] | None = None,
        discrepancies: list[dict[str, Any]] | None = None,
        proposed_tasks: list[dict[str, Any]] | None = None,
        blind_to_prior_scores: bool = True,
        blind_to_other_reviews: bool = True,
        parent_review_id: str | None = None,
        idempotency_key: str,
        schema_version: int = 1,
    ) -> str:
        if (
            sum(
                target is not None
                for target in (finding_id, cluster_id, unit_disposition_id)
            )
            != 1
        ):
            raise ValueError(
                "Review requires exactly one finding_id, cluster_id, or "
                "unit_disposition_id"
            )
        if verdict not in REVIEW_VERDICTS:
            raise ValueError(f"Unsupported review verdict: {verdict}")
        if not 0.0 <= confidence <= 1.0:
            raise ValueError("Review confidence must be between 0 and 1")
        dataset_id = _nonempty(dataset_id, "dataset_id")
        review_tier = _nonempty(review_tier, "review_tier")
        execution_profile = _nonempty(execution_profile, "execution_profile")
        provider = _nonempty(provider, "provider")
        model = _nonempty(model, "model")
        reasoning_effort = _nonempty(reasoning_effort, "reasoning_effort")
        prompt_hash = _nonempty(prompt_hash, "prompt_hash")
        rubric_version = _nonempty(rubric_version, "rubric_version")
        input_bundle_hash = _nonempty(input_bundle_hash, "input_bundle_hash")
        assert_sequence_free(
            {
                "reconstructed_observations": reconstructed_observations or {},
                "claim_assessment": claim_assessment or {},
                "verification_summary": verification_summary or {},
                "discrepancies": discrepancies or [],
                "proposed_tasks": proposed_tasks or [],
            },
            field="Finding review",
        )
        json_fields = (
            canonical_json(
                reconstructed_observations or {}, field="reconstructed observations"
            ),
            canonical_json(claim_assessment or {}, field="claim assessment"),
            canonical_json(verification_summary or {}, field="verification summary"),
            canonical_json(discrepancies or [], field="discrepancies"),
            canonical_json(proposed_tasks or [], field="proposed tasks"),
        )
        with self._lock, self._transaction():
            self._validate_campaign_locked(campaign_id)
            task = None
            if finding_id is not None:
                self._validate_ids_locked("findings", [finding_id], label="findings")
                target = self._conn.execute(
                    "SELECT campaign_id FROM findings WHERE id = ?", (finding_id,)
                ).fetchone()
                if target["campaign_id"] != campaign_id:
                    raise ValueError("Review finding belongs to a different campaign")
            if cluster_id is not None:
                self._validate_ids_locked(
                    "candidate_clusters", [cluster_id], label="candidate clusters"
                )
                target = self._conn.execute(
                    """
                    SELECT campaign_id, dataset_id FROM candidate_clusters
                    WHERE id = ?
                    """,
                    (cluster_id,),
                ).fetchone()
                if (
                    target["campaign_id"] != campaign_id
                    or target["dataset_id"] != dataset_id
                ):
                    raise ValueError(
                        "Review cluster belongs to a different campaign or dataset"
                    )
            if unit_disposition_id is not None:
                self._validate_ids_locked(
                    "unit_dispositions",
                    [unit_disposition_id],
                    label="unit dispositions",
                )
                target = self._conn.execute(
                    """
                    SELECT campaign_id, dataset_id FROM unit_dispositions
                    WHERE id = ?
                    """,
                    (unit_disposition_id,),
                ).fetchone()
                if (
                    target["campaign_id"] != campaign_id
                    or target["dataset_id"] != dataset_id
                ):
                    raise ValueError(
                        "Review unit belongs to a different campaign or dataset"
                    )
            if task_id is not None:
                self._validate_ids_locked("tasks", [task_id], label="tasks")
                task = self._conn.execute(
                    """
                    SELECT
                        campaign_id, params, assigned_to, status,
                        lease_expires_ts
                    FROM tasks WHERE id = ?
                    """,
                    (task_id,),
                ).fetchone()
                if task["campaign_id"] != campaign_id:
                    raise ValueError("Review task belongs to a different campaign")
                if (
                    task["assigned_to"] != self.agent_id
                    or task["status"]
                    not in {"claimed", "in_progress", "complete"}
                ):
                    raise ValueError(
                        "Review task is not owned by the submitting reviewer"
                    )
                task_params = json.loads(task["params"])
                contract = task_params.get("review_contract")
                if contract is not None:
                    target = task_params.get("target", {})
                    target_kind = target.get("kind")
                    target_id = target.get("id")
                    expected_target = {
                        "finding": finding_id,
                        "candidate_cluster": cluster_id,
                        "unit_disposition": unit_disposition_id,
                    }
                    if expected_target.get(target_kind) != target_id:
                        raise ValueError("Review task targets a different subject")
                    if task_params.get("review_tier") != review_tier:
                        raise ValueError("Review tier differs from its task contract")
                    if task_params.get("execution_profile") != execution_profile:
                        raise ValueError(
                            "Execution profile differs from its task contract"
                        )
                    resolved = task_params.get("resolved_execution", {})
                    if (
                        resolved.get("provider") != provider
                        or resolved.get("model") != model
                        or resolved.get("reasoning_effort") != reasoning_effort
                    ):
                        raise ValueError(
                            "Review execution identity differs from its task contract"
                        )
                    if (
                        bool(task_params.get("blind_to_prior_scores"))
                        != blind_to_prior_scores
                        or bool(task_params.get("blind_to_other_reviews"))
                        != blind_to_other_reviews
                    ):
                        raise ValueError(
                            "Review blindness flags differ from its task contract"
                        )
            if run_id is not None:
                self._validate_ids_locked("runs", [run_id], label="runs")
            if parent_review_id is not None:
                self._validate_ids_locked(
                    "finding_reviews", [parent_review_id], label="finding reviews"
                )
                parent = self._conn.execute(
                    """
                    SELECT campaign_id, finding_id, cluster_id, unit_disposition_id
                    FROM finding_reviews WHERE id = ?
                    """,
                    (parent_review_id,),
                ).fetchone()
                if (
                    parent["campaign_id"] != campaign_id
                    or parent["finding_id"] != finding_id
                    or parent["cluster_id"] != cluster_id
                    or parent["unit_disposition_id"] != unit_disposition_id
                ):
                    raise ValueError("Parent review targets a different subject")
            if task_id is not None:
                existing = self._conn.execute(
                    "SELECT * FROM finding_reviews WHERE task_id = ?",
                    (task_id,),
                ).fetchone()
            else:
                existing = self._conn.execute(
                    """
                    SELECT * FROM finding_reviews
                    WHERE task_id IS NULL
                      AND reviewer_agent_id = ? AND idempotency_key = ?
                    """,
                    (self.agent_id, idempotency_key),
                ).fetchone()
            existing_id = _idempotent_id(
                existing,
                {
                    "campaign_id": campaign_id,
                    "finding_id": finding_id,
                    "cluster_id": cluster_id,
                    "unit_disposition_id": unit_disposition_id,
                    "task_id": task_id,
                    "run_id": run_id,
                    "review_tier": review_tier,
                    "execution_profile": execution_profile,
                    "provider": provider,
                    "model": model,
                    "reasoning_effort": reasoning_effort,
                    "prompt_hash": prompt_hash,
                    "rubric_version": rubric_version,
                    "input_bundle_hash": input_bundle_hash,
                    "dataset_id": dataset_id,
                    "verdict": verdict,
                    "reconstructed_observations": json_fields[0],
                    "claim_assessment": json_fields[1],
                    "verification_summary": json_fields[2],
                    "discrepancies": json_fields[3],
                    "proposed_tasks": json_fields[4],
                    "confidence": confidence,
                    "blind_to_prior_scores": int(blind_to_prior_scores),
                    "blind_to_other_reviews": int(blind_to_other_reviews),
                    "parent_review_id": parent_review_id,
                    "schema_version": schema_version,
                },
                record_type="Finding-review",
            )
            if existing_id is not None:
                return existing_id
            if task is not None and (
                task["status"] == "complete"
                or float(task["lease_expires_ts"] or 0) <= time.time()
            ):
                raise ValueError(
                    "Review task has no active lease owned by the submitting reviewer"
                )
            review_id = str(uuid.uuid4())
            self._conn.execute(
                """
                INSERT INTO finding_reviews(
                    id, campaign_id, finding_id, cluster_id, unit_disposition_id,
                    reviewer_agent_id, task_id, run_id, ts, review_tier,
                    execution_profile, provider,
                    model, reasoning_effort, prompt_hash, rubric_version,
                    input_bundle_hash, dataset_id, verdict,
                    reconstructed_observations, claim_assessment,
                    verification_summary, discrepancies, proposed_tasks, confidence,
                    blind_to_prior_scores, blind_to_other_reviews, parent_review_id,
                    idempotency_key, schema_version
                ) VALUES (
                    ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,
                    ?, ?, ?, ?, ?, ?, ?, ?, ?, ?
                )
                """,
                (
                    review_id,
                    campaign_id,
                    finding_id,
                    cluster_id,
                    unit_disposition_id,
                    self.agent_id,
                    task_id,
                    run_id,
                    time.time(),
                    review_tier,
                    execution_profile,
                    provider,
                    model,
                    reasoning_effort,
                    prompt_hash,
                    rubric_version,
                    input_bundle_hash,
                    dataset_id,
                    verdict,
                    *json_fields,
                    confidence,
                    int(blind_to_prior_scores),
                    int(blind_to_other_reviews),
                    parent_review_id,
                    idempotency_key,
                    schema_version,
                ),
            )
            event_id = self._event_locked(
                "finding_review_created",
                "finding_review",
                review_id,
                campaign_id=campaign_id,
                task_id=task_id,
                run_id=run_id,
                payload={
                    "finding_id": finding_id,
                    "cluster_id": cluster_id,
                    "unit_disposition_id": unit_disposition_id,
                    "review_tier": review_tier,
                    "verdict": verdict,
                    "execution_profile": execution_profile,
                },
            )
            self._conn.execute(
                "UPDATE finding_reviews SET created_event_id = ? WHERE id = ?",
                (event_id, review_id),
            )
            return review_id

    def get_finding_review(self, review_id: str) -> dict[str, Any]:
        with self._lock:
            row = self._conn.execute(
                "SELECT * FROM finding_reviews WHERE id = ?", (review_id,)
            ).fetchone()
        result = decode_review_row(row)
        if result is None:
            raise KeyError(f"Finding review {review_id} not found")
        return result

    def list_finding_reviews(
        self,
        *,
        campaign_id: str,
        finding_id: str | None = None,
        cluster_id: str | None = None,
        unit_disposition_id: str | None = None,
        review_tier: str | None = None,
        verdict: str | None = None,
        limit: int = 1_000,
    ) -> list[dict[str, Any]]:
        sql = "SELECT * FROM finding_reviews WHERE campaign_id = ?"
        params: list[Any] = [campaign_id]
        for column, value in (
            ("finding_id", finding_id),
            ("cluster_id", cluster_id),
            ("unit_disposition_id", unit_disposition_id),
            ("review_tier", review_tier),
            ("verdict", verdict),
        ):
            if value is not None:
                sql += f" AND {column} = ?"
                params.append(value)
        sql += " ORDER BY ts, id LIMIT ?"
        params.append(limit)
        with self._lock:
            rows = self._conn.execute(sql, params).fetchall()
        return [decode_review_row(row) for row in rows]

    def record_review_verification(
        self,
        *,
        review_id: str,
        claim_key: str,
        engine: str,
        specification: dict[str, Any],
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
        if status not in {"pending", "pass", "fail", "error", "skipped"}:
            raise ValueError(f"Unsupported verification status: {status}")
        if status in {"pass", "fail"} and actual is None:
            raise ValueError(f"Verification status {status} requires an actual result")
        if status in {"pass", "fail", "error", "skipped"} and executed_ts is None:
            raise ValueError(f"Verification status {status} requires executed_ts")
        if status == "pending" and executed_ts is not None:
            raise ValueError("Pending verification cannot have executed_ts")
        dataset_id = _nonempty(dataset_id, "dataset_id")
        claim_key = _nonempty(claim_key, "claim_key")
        engine = _nonempty(engine, "engine")
        specification_json = canonical_json(
            specification, field="verification specification"
        )
        specification_hash = hashlib.sha256(specification_json.encode()).hexdigest()
        expected_json = canonical_json(expected, field="verification expected")
        actual_json = (
            None
            if actual is None
            else canonical_json(actual, field="verification actual")
        )
        assert_sequence_free(
            {
                "specification": specification,
                "expected": expected,
                "actual": actual,
            },
            field="Review verification",
        )
        with self._lock, self._transaction():
            self._validate_ids_locked(
                "finding_reviews", [review_id], label="finding reviews"
            )
            review = self._conn.execute(
                """
                SELECT campaign_id, task_id, run_id, dataset_id
                FROM finding_reviews WHERE id = ?
                """,
                (review_id,),
            ).fetchone()
            if review["dataset_id"] != dataset_id:
                raise ValueError("Verification dataset differs from review dataset")
            if artifact_id is not None:
                self._validate_ids_locked("artifacts", [artifact_id], label="artifacts")
            attempt = 1
            previous = None
            if supersedes_verification_id is not None:
                previous = self._conn.execute(
                    """
                    SELECT * FROM review_verifications
                    WHERE id = ?
                    """,
                    (supersedes_verification_id,),
                ).fetchone()
                if previous is None:
                    raise ValueError(
                        f"Review verification {supersedes_verification_id} not found"
                    )
                if (
                    previous["review_id"] != review_id
                    or previous["claim_key"] != claim_key
                    or previous["specification_hash"] != specification_hash
                ):
                    raise ValueError(
                        "Superseded verification belongs to a different claim "
                        "or specification"
                    )
                attempt = int(previous["attempt"]) + 1
            immutable = {
                "review_id": review_id,
                "claim_key": claim_key,
                "engine": engine,
                "specification": specification_json,
                "dataset_id": dataset_id,
                "specification_hash": specification_hash,
                "expected": expected_json,
                "actual": actual_json,
                "status": status,
                "attempt": attempt,
                "supersedes_verification_id": supersedes_verification_id,
                "executed_ts": executed_ts,
                "code_commit": code_commit,
                "artifact_id": artifact_id,
                "error": error,
            }
            existing_by_key = self._conn.execute(
                """
                SELECT * FROM review_verifications
                WHERE created_by = ? AND idempotency_key = ?
                """,
                (self.agent_id, idempotency_key),
            ).fetchone()
            existing_id = _idempotent_id(
                existing_by_key,
                immutable,
                record_type="Review-verification",
            )
            if existing_id is not None:
                return existing_id
            existing_attempt = self._conn.execute(
                """
                SELECT * FROM review_verifications
                WHERE review_id = ? AND claim_key = ?
                  AND specification_hash = ? AND attempt = ?
                """,
                (review_id, claim_key, specification_hash, attempt),
            ).fetchone()
            existing_id = _idempotent_id(
                existing_attempt,
                immutable,
                record_type="Review-verification attempt",
            )
            if existing_id is not None:
                return existing_id
            if previous is not None:
                latest = self._conn.execute(
                    """
                    SELECT id FROM review_verifications
                    WHERE review_id = ? AND claim_key = ?
                      AND specification_hash = ?
                    ORDER BY attempt DESC, created_ts DESC, id DESC
                    LIMIT 1
                    """,
                    (review_id, claim_key, specification_hash),
                ).fetchone()
                if latest is None or latest["id"] != previous["id"]:
                    raise ValueError(
                        "Verification revision must supersede the latest attempt"
                    )
            verification_id = str(uuid.uuid4())
            self._conn.execute(
                """
                INSERT INTO review_verifications(
                    id, review_id, claim_key, engine, specification, dataset_id,
                    specification_hash, expected, actual, status, attempt,
                    supersedes_verification_id, executed_ts, code_commit, artifact_id,
                    error, created_by, created_ts, idempotency_key
                ) VALUES (
                    ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?
                )
                """,
                (
                    verification_id,
                    review_id,
                    claim_key,
                    engine,
                    specification_json,
                    dataset_id,
                    specification_hash,
                    expected_json,
                    actual_json,
                    status,
                    attempt,
                    supersedes_verification_id,
                    executed_ts,
                    code_commit,
                    artifact_id,
                    error,
                    self.agent_id,
                    time.time(),
                    idempotency_key,
                ),
            )
            event_id = self._event_locked(
                "review_verification_created",
                "review_verification",
                verification_id,
                campaign_id=review["campaign_id"],
                task_id=review["task_id"],
                run_id=review["run_id"],
                payload={
                    "review_id": review_id,
                    "claim_key": claim_key,
                    "status": status,
                    "attempt": attempt,
                },
            )
            self._conn.execute(
                "UPDATE review_verifications SET created_event_id = ? WHERE id = ?",
                (event_id, verification_id),
            )
            return verification_id

    def list_review_verifications(self, review_id: str) -> list[dict[str, Any]]:
        with self._lock:
            rows = self._conn.execute(
                """
                SELECT * FROM review_verifications
                WHERE review_id = ?
                ORDER BY claim_key, specification_hash, attempt, created_ts, id
                """,
                (review_id,),
            ).fetchall()
        return [decode_review_row(row) for row in rows]

    def _assert_materialized_sources_current_locked(self, finding_id: str) -> None:
        rows = self._conn.execute(
            """
            SELECT c.id, c.status
            FROM cluster_findings AS cf
            JOIN candidate_clusters AS c ON c.id = cf.cluster_id
            WHERE cf.finding_id = ? AND cf.relation = 'materializes'
            """,
            (finding_id,),
        ).fetchall()
        stale = [str(row["id"]) for row in rows if row["status"] != "active"]
        if stale:
            raise ValueError(
                "Finding materialization references superseded candidate "
                f"clusters: {', '.join(stale)}"
            )

    def _review_campaign_scan_status_locked(
        self,
        campaign_id: str,
    ) -> dict[str, Any]:
        campaign = self._conn.execute(
            "SELECT metadata FROM campaigns WHERE id = ?",
            (campaign_id,),
        ).fetchone()
        if campaign is None:
            raise KeyError(f"Campaign {campaign_id} not found")
        metadata = json.loads(campaign["metadata"] or "{}")
        expected = metadata.get("unit_count")
        if expected is None:
            return {
                "required": False,
                "ready": True,
                "expected_units": None,
                "atlas_tasks": None,
                "completed_atlas_tasks": None,
                "active_unit_dispositions": None,
                "ready_unit_dispositions": None,
            }
        if isinstance(expected, bool) or not isinstance(expected, int) or expected < 0:
            raise ValueError("Campaign metadata unit_count must be a non-negative integer")
        task_counts = self._conn.execute(
            """
            SELECT
                COUNT(*) AS total,
                SUM(CASE WHEN status = 'complete' THEN 1 ELSE 0 END) AS complete
            FROM tasks
            WHERE campaign_id = ? AND task_type = 'atlas_genome_read'
            """,
            (campaign_id,),
        ).fetchone()
        dataset_id = metadata.get("dataset_id")
        disposition_sql = """
            SELECT
                COUNT(DISTINCT unit_id) AS active,
                COUNT(
                    DISTINCT CASE
                        WHEN disposition IN ('clear', 'candidate') THEN unit_id
                    END
                ) AS ready
            FROM unit_dispositions
            WHERE campaign_id = ? AND record_status = 'active'
        """
        disposition_params: list[Any] = [campaign_id]
        if dataset_id is not None:
            disposition_sql += " AND dataset_id = ?"
            disposition_params.append(str(dataset_id))
        disposition_counts = self._conn.execute(
            disposition_sql,
            disposition_params,
        ).fetchone()
        total_tasks = int(task_counts["total"] or 0)
        complete_tasks = int(task_counts["complete"] or 0)
        active_units = int(disposition_counts["active"] or 0)
        ready_units = int(disposition_counts["ready"] or 0)
        return {
            "required": True,
            "ready": (
                total_tasks == expected
                and complete_tasks == expected
                and active_units == expected
                and ready_units == expected
            ),
            "expected_units": expected,
            "atlas_tasks": total_tasks,
            "completed_atlas_tasks": complete_tasks,
            "active_unit_dispositions": active_units,
            "ready_unit_dispositions": ready_units,
        }

    def review_campaign_scan_status(self, campaign_id: str) -> dict[str, Any]:
        """Return the exact Atlas closure gate for canonical publication."""

        with self._lock:
            return self._review_campaign_scan_status_locked(campaign_id)

    def _assert_review_campaign_scan_ready_locked(self, campaign_id: str) -> None:
        status = self._review_campaign_scan_status_locked(campaign_id)
        if not status["ready"]:
            raise ValueError(
                "Canonical publication requires complete Atlas coverage: "
                f"{status['completed_atlas_tasks']}/{status['expected_units']} "
                "scan tasks complete and "
                f"{status['ready_unit_dispositions']}/"
                f"{status['expected_units']} unit dispositions ready"
            )

    def create_promotion_decision(
        self,
        *,
        campaign_id: str,
        decision: str,
        source_tier: str,
        policy_name: str,
        policy_version: str,
        policy_hash: str,
        rationale: str,
        finding_id: str | None = None,
        cluster_id: str | None = None,
        target_tier: str | None = None,
        review_ids: list[str] | None = None,
        created_task_ids: list[str] | None = None,
        audit_sample: bool = False,
        audit_stratum: dict[str, Any] | None = None,
        idempotency_key: str,
        schema_version: int = 1,
    ) -> str:
        if (finding_id is None) == (cluster_id is None):
            raise ValueError("Promotion requires exactly one finding_id or cluster_id")
        if decision not in PROMOTION_DECISIONS:
            raise ValueError(f"Unsupported promotion decision: {decision}")
        source_tier = _nonempty(source_tier, "source_tier")
        policy_name = _nonempty(policy_name, "policy_name")
        policy_version = _nonempty(policy_version, "policy_version")
        policy_hash = _nonempty(policy_hash, "policy_hash")
        rationale = _nonempty(rationale, "rationale")
        reviews = sorted(set(review_ids or []))
        tasks = sorted(set(created_task_ids or []))
        assert_sequence_free(
            {
                "rationale": rationale,
                "audit_stratum": audit_stratum or {},
            },
            field="Promotion decision",
        )
        if decision == "publish" and (finding_id is None or not reviews):
            raise ValueError(
                "Publish decisions require a finding target and supporting reviews"
            )
        with self._lock, self._transaction():
            self._validate_campaign_locked(campaign_id)
            if finding_id is not None:
                self._validate_ids_locked("findings", [finding_id], label="findings")
                target = self._conn.execute(
                    "SELECT campaign_id FROM findings WHERE id = ?", (finding_id,)
                ).fetchone()
                if target["campaign_id"] != campaign_id:
                    raise ValueError("Decision finding belongs to a different campaign")
            if cluster_id is not None:
                self._validate_ids_locked(
                    "candidate_clusters", [cluster_id], label="candidate clusters"
                )
                target = self._conn.execute(
                    """
                    SELECT campaign_id, status FROM candidate_clusters
                    WHERE id = ?
                    """,
                    (cluster_id,),
                ).fetchone()
                if target["campaign_id"] != campaign_id:
                    raise ValueError("Decision cluster belongs to a different campaign")
                if (
                    decision in {"promote", "reopen"}
                    and target["status"] != "active"
                ):
                    raise ValueError(
                        "Promotion requires the active candidate-cluster version"
                    )
            self._validate_ids_locked("finding_reviews", reviews, label="finding reviews")
            self._validate_ids_locked("tasks", tasks, label="tasks")
            for review_id in reviews:
                review = self._conn.execute(
                    "SELECT finding_id, cluster_id FROM finding_reviews WHERE id = ?",
                    (review_id,),
                ).fetchone()
                if review["finding_id"] != finding_id or review["cluster_id"] != cluster_id:
                    raise ValueError("Promotion reviews target a different subject")
            if decision == "publish":
                assert finding_id is not None
                self._assert_review_campaign_scan_ready_locked(campaign_id)
                self._assert_materialized_sources_current_locked(finding_id)
                for review_id in reviews:
                    review = self._conn.execute(
                        "SELECT verdict FROM finding_reviews WHERE id = ?",
                        (review_id,),
                    ).fetchone()
                    if review["verdict"] != "promote":
                        raise ValueError("Publish decision includes a non-promote review")
                    verification_rows = self._conn.execute(
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
                    if not verification_rows or any(
                        row["status"] != "pass" for row in verification_rows
                    ):
                        raise ValueError(
                            "Publish decision requires passed executable "
                            "verification records for every supporting review"
                        )
            existing = self._conn.execute(
                """
                SELECT * FROM promotion_decisions
                WHERE actor_agent_id = ? AND idempotency_key = ?
                """,
                (self.agent_id, idempotency_key),
            ).fetchone()
            audit_stratum_json = canonical_json(
                audit_stratum or {}, field="audit stratum"
            )
            existing_id = _idempotent_id(
                existing,
                {
                    "campaign_id": campaign_id,
                    "finding_id": finding_id,
                    "cluster_id": cluster_id,
                    "decision": decision,
                    "source_tier": source_tier,
                    "target_tier": target_tier,
                    "policy_name": policy_name,
                    "policy_version": policy_version,
                    "policy_hash": policy_hash,
                    "rationale": rationale,
                    "audit_sample": int(audit_sample),
                    "audit_stratum": audit_stratum_json,
                    "schema_version": schema_version,
                },
                record_type="Promotion-decision",
            )
            if existing_id is not None:
                current_reviews = [
                    row["review_id"]
                    for row in self._conn.execute(
                        """
                        SELECT review_id FROM promotion_decision_reviews
                        WHERE decision_id = ? ORDER BY review_id
                        """,
                        (existing_id,),
                    ).fetchall()
                ]
                current_tasks = [
                    row["task_id"]
                    for row in self._conn.execute(
                        """
                        SELECT task_id FROM promotion_decision_tasks
                        WHERE decision_id = ? ORDER BY task_id
                        """,
                        (existing_id,),
                    ).fetchall()
                ]
                if current_reviews != reviews or current_tasks != tasks:
                    raise ValueError(
                        "Promotion-decision idempotency key has different "
                        "review or task links"
                    )
                return existing_id
            decision_id = str(uuid.uuid4())
            self._conn.execute(
                """
                INSERT INTO promotion_decisions(
                    id, campaign_id, finding_id, cluster_id, actor_agent_id, ts,
                    decision, source_tier, target_tier, policy_name, policy_version,
                    policy_hash, rationale, audit_sample, audit_stratum,
                    idempotency_key, schema_version
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    decision_id,
                    campaign_id,
                    finding_id,
                    cluster_id,
                    self.agent_id,
                    time.time(),
                    decision,
                    source_tier,
                    target_tier,
                    policy_name,
                    policy_version,
                    policy_hash,
                    rationale,
                    int(audit_sample),
                    audit_stratum_json,
                    idempotency_key,
                    schema_version,
                ),
            )
            for review_id in reviews:
                self._conn.execute(
                    """
                    INSERT INTO promotion_decision_reviews(decision_id, review_id)
                    VALUES (?, ?)
                    """,
                    (decision_id, review_id),
                )
            for task_id in tasks:
                self._conn.execute(
                    """
                    INSERT INTO promotion_decision_tasks(decision_id, task_id, relation)
                    VALUES (?, ?, 'created')
                    """,
                    (decision_id, task_id),
                )
            event_id = self._event_locked(
                "promotion_decision_created",
                "promotion_decision",
                decision_id,
                campaign_id=campaign_id,
                payload={
                    "finding_id": finding_id,
                    "cluster_id": cluster_id,
                    "decision": decision,
                    "source_tier": source_tier,
                    "target_tier": target_tier,
                    "audit_sample": audit_sample,
                },
            )
            self._conn.execute(
                "UPDATE promotion_decisions SET created_event_id = ? WHERE id = ?",
                (event_id, decision_id),
            )
            return decision_id

    def get_promotion_decision(self, decision_id: str) -> dict[str, Any]:
        with self._lock:
            row = self._conn.execute(
                "SELECT * FROM promotion_decisions WHERE id = ?", (decision_id,)
            ).fetchone()
            reviews = self._conn.execute(
                """
                SELECT review_id FROM promotion_decision_reviews
                WHERE decision_id = ? ORDER BY review_id
                """,
                (decision_id,),
            ).fetchall()
            tasks = self._conn.execute(
                """
                SELECT task_id, relation FROM promotion_decision_tasks
                WHERE decision_id = ? ORDER BY task_id, relation
                """,
                (decision_id,),
            ).fetchall()
        result = decode_review_row(row)
        if result is None:
            raise KeyError(f"Promotion decision {decision_id} not found")
        result["review_ids"] = [item["review_id"] for item in reviews]
        result["tasks"] = [dict(item) for item in tasks]
        return result

    def list_promotion_decisions(
        self,
        *,
        campaign_id: str,
        finding_id: str | None = None,
        cluster_id: str | None = None,
        decision: str | None = None,
        limit: int = 1_000,
    ) -> list[dict[str, Any]]:
        sql = "SELECT * FROM promotion_decisions WHERE campaign_id = ?"
        params: list[Any] = [campaign_id]
        for column, value in (
            ("finding_id", finding_id),
            ("cluster_id", cluster_id),
            ("decision", decision),
        ):
            if value is not None:
                sql += f" AND {column} = ?"
                params.append(value)
        sql += " ORDER BY ts, id LIMIT ?"
        params.append(limit)
        with self._lock:
            rows = self._conn.execute(sql, params).fetchall()
        return [decode_review_row(row) for row in rows]

    def record_canonical_publication(
        self,
        *,
        campaign_id: str,
        finding_id: str,
        decision_id: str,
        dataset_id: str,
        canonical_uri: str,
        canonical_record_id: str,
        canonical_record_hash: str,
        metadata: dict[str, Any] | None = None,
        idempotency_key: str,
    ) -> str:
        dataset_id = _nonempty(dataset_id, "dataset_id")
        canonical_uri = _nonempty(canonical_uri, "canonical_uri")
        canonical_record_id = _nonempty(
            canonical_record_id, "canonical_record_id"
        )
        canonical_record_hash = _nonempty(
            canonical_record_hash, "canonical_record_hash"
        )
        assert_sequence_free(
            metadata or {},
            field="Canonical publication metadata",
        )
        with self._lock, self._transaction():
            self._validate_campaign_locked(campaign_id)
            self._validate_ids_locked("findings", [finding_id], label="findings")
            decision = self._conn.execute(
                "SELECT * FROM promotion_decisions WHERE id = ?", (decision_id,)
            ).fetchone()
            if (
                decision is None
                or decision["decision"] != "publish"
                or decision["finding_id"] != finding_id
                or decision["campaign_id"] != campaign_id
            ):
                raise ValueError("Canonical publication requires a matching publish decision")
            self._assert_review_campaign_scan_ready_locked(campaign_id)
            self._assert_materialized_sources_current_locked(finding_id)
            review_datasets = {
                row["dataset_id"]
                for row in self._conn.execute(
                    """
                    SELECT r.dataset_id
                    FROM promotion_decision_reviews AS d
                    JOIN finding_reviews AS r ON r.id = d.review_id
                    WHERE d.decision_id = ?
                    """,
                    (decision_id,),
                ).fetchall()
            }
            if review_datasets != {dataset_id}:
                raise ValueError(
                    "Canonical publication dataset differs from supporting reviews"
                )
            verification_gate = self._conn.execute(
                """
                WITH linked_reviews AS (
                    SELECT r.review_id
                    FROM promotion_decision_reviews AS r
                    WHERE r.decision_id = ?
                ),
                latest AS (
                    SELECT
                        v.review_id,
                        v.status,
                        ROW_NUMBER() OVER (
                            PARTITION BY
                                v.review_id, v.claim_key, v.specification_hash
                            ORDER BY v.attempt DESC, v.created_ts DESC, v.id DESC
                        ) AS rank
                    FROM review_verifications AS v
                    JOIN linked_reviews AS r ON r.review_id = v.review_id
                ),
                per_review AS (
                    SELECT
                        r.review_id,
                        COUNT(l.status) AS verification_count,
                        SUM(CASE WHEN l.status = 'pass' THEN 0 ELSE 1 END) AS blocked
                    FROM linked_reviews AS r
                    LEFT JOIN latest AS l
                        ON l.review_id = r.review_id AND l.rank = 1
                    GROUP BY r.review_id
                )
                SELECT
                    COUNT(*) AS review_count,
                    SUM(
                        CASE
                            WHEN verification_count = 0 OR blocked > 0 THEN 1
                            ELSE 0
                        END
                    ) AS blocked_reviews
                FROM per_review
                """,
                (decision_id,),
            ).fetchone()
            if (
                not verification_gate
                or int(verification_gate["review_count"] or 0) == 0
                or int(verification_gate["blocked_reviews"] or 0) > 0
            ):
                raise ValueError("Canonical publication has unresolved verification records")
            existing = self._conn.execute(
                """
                SELECT * FROM canonical_publications
                WHERE published_by = ? AND idempotency_key = ?
                """,
                (self.agent_id, idempotency_key),
            ).fetchone()
            metadata_json = canonical_json(
                metadata or {}, field="publication metadata"
            )
            existing_id = _idempotent_id(
                existing,
                {
                    "campaign_id": campaign_id,
                    "finding_id": finding_id,
                    "decision_id": decision_id,
                    "dataset_id": dataset_id,
                    "canonical_uri": canonical_uri,
                    "canonical_record_id": canonical_record_id,
                    "canonical_record_hash": canonical_record_hash,
                    "metadata": metadata_json,
                },
                record_type="Canonical-publication",
            )
            if existing_id is not None:
                return existing_id
            publication_id = str(uuid.uuid4())
            self._conn.execute(
                """
                INSERT INTO canonical_publications(
                    id, campaign_id, finding_id, decision_id, dataset_id,
                    canonical_uri, canonical_record_id, canonical_record_hash,
                    published_by, published_ts, metadata, idempotency_key
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    publication_id,
                    campaign_id,
                    finding_id,
                    decision_id,
                    dataset_id,
                    canonical_uri,
                    canonical_record_id,
                    canonical_record_hash,
                    self.agent_id,
                    time.time(),
                    metadata_json,
                    idempotency_key,
                ),
            )
            event_id = self._event_locked(
                "canonical_publication_created",
                "canonical_publication",
                publication_id,
                campaign_id=campaign_id,
                payload={
                    "finding_id": finding_id,
                    "decision_id": decision_id,
                    "canonical_record_id": canonical_record_id,
                },
            )
            self._conn.execute(
                "UPDATE canonical_publications SET created_event_id = ? WHERE id = ?",
                (event_id, publication_id),
            )
            return publication_id

    def list_canonical_publications(
        self,
        *,
        campaign_id: str,
        finding_id: str | None = None,
        limit: int = 1_000,
    ) -> list[dict[str, Any]]:
        sql = "SELECT * FROM canonical_publications WHERE campaign_id = ?"
        params: list[Any] = [campaign_id]
        if finding_id is not None:
            sql += " AND finding_id = ?"
            params.append(finding_id)
        sql += " ORDER BY published_ts, id LIMIT ?"
        params.append(limit)
        with self._lock:
            rows = self._conn.execute(sql, params).fetchall()
        return [decode_review_row(row) for row in rows]

    def get_review_controller_state(
        self, controller_id: str, campaign_id: str
    ) -> dict[str, Any] | None:
        with self._lock:
            row = self._conn.execute(
                """
                SELECT * FROM review_controller_state
                WHERE controller_id = ? AND campaign_id = ?
                """,
                (controller_id, campaign_id),
            ).fetchone()
        return decode_review_row(row)

    def set_review_controller_state(
        self,
        controller_id: str,
        campaign_id: str,
        *,
        policy_hash: str,
        last_event_id: int,
    ) -> dict[str, Any]:
        if last_event_id < 0:
            raise ValueError("last_event_id must be non-negative")
        with self._lock, self._transaction():
            self._validate_campaign_locked(campaign_id)
            current = self._conn.execute(
                """
                SELECT last_event_id FROM review_controller_state
                WHERE controller_id = ? AND campaign_id = ?
                """,
                (controller_id, campaign_id),
            ).fetchone()
            if current is not None and int(current["last_event_id"]) > last_event_id:
                raise ValueError("Review controller cursor cannot move backwards")
            self._conn.execute(
                """
                INSERT INTO review_controller_state(
                    controller_id, campaign_id, policy_hash, last_event_id, updated_ts
                ) VALUES (?, ?, ?, ?, ?)
                ON CONFLICT(controller_id, campaign_id) DO UPDATE SET
                    policy_hash = excluded.policy_hash,
                    last_event_id = excluded.last_event_id,
                    updated_ts = excluded.updated_ts
                """,
                (
                    controller_id,
                    campaign_id,
                    policy_hash,
                    last_event_id,
                    time.time(),
                ),
            )
        return self.get_review_controller_state(controller_id, campaign_id)
