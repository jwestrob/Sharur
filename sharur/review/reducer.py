"""Deterministic, lossless reduction of typed candidate occurrences."""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Protocol

from sharur.ops.review_store import content_hash, decode_review_row


if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence

    from sharur.ops.store import OpsStore


EXACT_REDUCER_NAME = "typed-exact-signature"
EXACT_REDUCER_VERSION = "1.0"
CLUSTER_ROLE_HINTS = frozenset({"member", "medoid", "outlier", "counterexample"})


class CandidateReductionAdapter(Protocol):
    """Optional typed logic layered over exact signature equivalence."""

    name: str
    version: str

    def validate(self, candidate: dict[str, Any]) -> None:
        """Reject malformed typed candidates before they enter a cluster."""

    def role_hint(self, candidate: dict[str, Any]) -> str | None:
        """Return a typed role hint while retaining the full occurrence."""

    def strata(self, candidate: dict[str, Any]) -> dict[str, str]:
        """Return bounded categorical strata used in exact cluster counts."""


class DefaultReductionAdapter:
    """Generic adapter that reads explicit, scanner-authored reduction fields."""

    name = "default"
    version = "1.0"

    def validate(self, candidate: dict[str, Any]) -> None:
        if not isinstance(candidate.get("signature"), dict):
            raise ValueError("Candidate signature must be a JSON object")
        if not isinstance(candidate.get("subject_refs"), dict):
            raise ValueError("Candidate subject_refs must be a JSON object")
        if not isinstance(candidate.get("reduction_features"), dict):
            raise ValueError("Candidate reduction_features must be a JSON object")

    def role_hint(self, candidate: dict[str, Any]) -> str | None:
        value = candidate["reduction_features"].get("role_hint")
        if value is None:
            return None
        role = str(value)
        if role not in CLUSTER_ROLE_HINTS:
            raise ValueError(f"Unsupported candidate role_hint {role!r}")
        return role

    def strata(self, candidate: dict[str, Any]) -> dict[str, str]:
        features = candidate["reduction_features"]
        explicit = features.get("strata", {})
        if explicit is None:
            explicit = {}
        if not isinstance(explicit, dict):
            raise ValueError("reduction_features.strata must be an object")
        result = {
            str(key): str(value)
            for key, value in explicit.items()
            if value is not None
        }
        refs = candidate["subject_refs"]
        for key in ("phylum", "class", "order", "family", "genus"):
            value = refs.get(key)
            if value is not None and key not in result:
                result[key] = str(value)
        return result


@dataclass(frozen=True)
class ReductionResult:
    campaign_id: str
    reducer_name: str
    reducer_version: str
    candidate_groups_seen: int
    clusters_created: int
    clusters_unchanged: int
    clusters_versioned: int
    occurrences_represented: int

    def to_dict(self) -> dict[str, Any]:
        return {
            "campaign_id": self.campaign_id,
            "reducer_name": self.reducer_name,
            "reducer_version": self.reducer_version,
            "candidate_groups_seen": self.candidate_groups_seen,
            "clusters_created": self.clusters_created,
            "clusters_unchanged": self.clusters_unchanged,
            "clusters_versioned": self.clusters_versioned,
            "occurrences_represented": self.occurrences_represented,
        }


def _quality(candidate: dict[str, Any]) -> tuple[float, int, str]:
    features = candidate["reduction_features"]
    raw_score = features.get("quality_score", 0.0)
    score = (
        float(raw_score)
        if isinstance(raw_score, (int, float)) and not isinstance(raw_score, bool)
        else 0.0
    )
    verification_count = len(candidate.get("verification") or [])
    # min() over the negated quality tuple selects the best quality and then
    # the lexicographically stable candidate ID.
    return (-score, -verification_count, str(candidate["id"]))


def _roles(
    candidates: Sequence[dict[str, Any]],
    adapter: CandidateReductionAdapter,
) -> tuple[dict[str, str], str]:
    hints: dict[str, str] = {}
    ordinary: list[dict[str, Any]] = []
    explicit_medoids: list[dict[str, Any]] = []
    for candidate in candidates:
        hint = adapter.role_hint(candidate)
        if hint is not None:
            hints[str(candidate["id"])] = hint
        if hint == "medoid":
            explicit_medoids.append(candidate)
        elif hint not in {"outlier", "counterexample"}:
            ordinary.append(candidate)
    pool = explicit_medoids or ordinary or list(candidates)
    medoid = min(pool, key=_quality)
    medoid_id = str(medoid["id"])
    roles = {
        str(candidate["id"]): hints.get(str(candidate["id"]), "member")
        for candidate in candidates
    }
    for candidate in explicit_medoids:
        roles[str(candidate["id"])] = "member"
    roles[medoid_id] = "medoid"
    return roles, medoid_id


def _categorical_counts(
    candidates: Sequence[dict[str, Any]],
    adapter: CandidateReductionAdapter,
) -> dict[str, dict[str, int]]:
    counters: dict[str, Counter[str]] = {}
    for candidate in candidates:
        for field, value in adapter.strata(candidate).items():
            counters.setdefault(field, Counter())[value] += 1
    return {
        field: dict(sorted(counter.items()))
        for field, counter in sorted(counters.items())
    }


class ExactCandidateReducer:
    """Cluster candidates only across an exact typed-signature boundary.

    The reducer performs no prose similarity or model inference. Every
    occurrence remains a member row, and cluster revisions form an explicit
    supersession lineage when later scans add occurrences.
    """

    def __init__(
        self,
        *,
        adapters: dict[str, CandidateReductionAdapter] | None = None,
    ):
        self._default_adapter = DefaultReductionAdapter()
        self._adapters = dict(adapters or {})

    def adapter_for(self, candidate_type: str) -> CandidateReductionAdapter:
        return self._adapters.get(candidate_type, self._default_adapter)

    def _config_hash(
        self,
        candidate_type: str,
        adapter: CandidateReductionAdapter,
    ) -> str:
        return content_hash(
            {
                "reducer_name": EXACT_REDUCER_NAME,
                "reducer_version": EXACT_REDUCER_VERSION,
                "candidate_type": candidate_type,
                "adapter_name": adapter.name,
                "adapter_version": adapter.version,
                "equivalence": (
                    "dataset_id,candidate_type,signature_schema,signature_hash"
                ),
                "representative": (
                    "explicit-medoid-or-highest-quality-then-verification-count-then-id"
                ),
            }
        )

    @staticmethod
    def _logical_cluster_id(
        *,
        campaign_id: str,
        dataset_id: str,
        candidate_type: str,
        signature_schema: str,
        signature_hash: str,
    ) -> str:
        identity = {
            "campaign_id": campaign_id,
            "dataset_id": dataset_id,
            "candidate_type": candidate_type,
            "signature_schema": signature_schema,
            "signature_hash": signature_hash,
        }
        digest = hashlib.sha256(
            json.dumps(identity, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest()
        return f"candidate-cluster-{digest[:32]}"

    @staticmethod
    def _signature_groups(
        store: OpsStore,
        campaign_id: str,
        *,
        dataset_id: str | None,
        candidate_type: str | None,
        batch_size: int,
    ) -> Iterable[tuple[str, str, str, str]]:
        cursor: tuple[str, str, str, str] | None = None
        while True:
            clauses = ["campaign_id = ?"]
            params: list[Any] = [campaign_id]
            if dataset_id is not None:
                clauses.append("dataset_id = ?")
                params.append(dataset_id)
            if candidate_type is not None:
                clauses.append("candidate_type = ?")
                params.append(candidate_type)
            if cursor is not None:
                clauses.append(
                    "(dataset_id, candidate_type, signature_schema, signature_hash) "
                    "> (?, ?, ?, ?)"
                )
                params.extend(cursor)
            params.append(batch_size)
            with store._lock:
                rows = store._conn.execute(
                    f"""
                    SELECT DISTINCT
                        dataset_id, candidate_type, signature_schema, signature_hash
                    FROM candidate_occurrences
                    WHERE {' AND '.join(clauses)}
                    ORDER BY
                        dataset_id, candidate_type, signature_schema, signature_hash
                    LIMIT ?
                    """,
                    params,
                ).fetchall()
            if not rows:
                return
            for row in rows:
                cursor = (
                    str(row["dataset_id"]),
                    str(row["candidate_type"]),
                    str(row["signature_schema"]),
                    str(row["signature_hash"]),
                )
                yield cursor

    @staticmethod
    def _members(
        store: OpsStore,
        campaign_id: str,
        key: tuple[str, str, str, str],
    ) -> list[dict[str, Any]]:
        dataset_id, candidate_type, signature_schema, signature_hash = key
        with store._lock:
            rows = store._conn.execute(
                """
                SELECT * FROM candidate_occurrences
                WHERE campaign_id = ?
                  AND dataset_id = ?
                  AND candidate_type = ?
                  AND signature_schema = ?
                  AND signature_hash = ?
                ORDER BY id
                """,
                (
                    campaign_id,
                    dataset_id,
                    candidate_type,
                    signature_schema,
                    signature_hash,
                ),
            ).fetchall()
        return [decode_review_row(row) for row in rows]

    def reduce_campaign(
        self,
        store: OpsStore,
        campaign_id: str,
        *,
        dataset_id: str | None = None,
        candidate_type: str | None = None,
        batch_size: int = 1_000,
    ) -> ReductionResult:
        if not 1 <= batch_size <= 10_000:
            raise ValueError("batch_size must be between 1 and 10000")
        store.get_campaign(campaign_id)
        groups_seen = 0
        created = 0
        unchanged = 0
        versioned = 0
        represented = 0
        for key in self._signature_groups(
            store,
            campaign_id,
            dataset_id=dataset_id,
            candidate_type=candidate_type,
            batch_size=batch_size,
        ):
            group_dataset, group_type, signature_schema, signature_hash = key
            candidates = self._members(store, campaign_id, key)
            if not candidates:
                continue
            groups_seen += 1
            represented += len(candidates)
            adapter = self.adapter_for(group_type)
            for candidate in candidates:
                adapter.validate(candidate)
            signatures = {
                json.dumps(
                    candidate["signature"],
                    sort_keys=True,
                    separators=(",", ":"),
                )
                for candidate in candidates
            }
            if len(signatures) != 1:
                raise RuntimeError(
                    "Signature-hash collision crossed distinct candidate signatures"
                )
            roles, medoid_id = _roles(candidates, adapter)
            role_counts = dict(sorted(Counter(roles.values()).items()))
            config_hash = self._config_hash(group_type, adapter)
            member_ids = [str(candidate["id"]) for candidate in candidates]
            manifest_hash = content_hash(sorted(member_ids))
            logical_id = self._logical_cluster_id(
                campaign_id=campaign_id,
                dataset_id=group_dataset,
                candidate_type=group_type,
                signature_schema=signature_schema,
                signature_hash=signature_hash,
            )
            with store._lock:
                existing = store._conn.execute(
                    """
                    SELECT * FROM candidate_clusters
                    WHERE logical_cluster_id = ?
                    ORDER BY version DESC LIMIT 1
                    """,
                    (logical_id,),
                ).fetchone()
            if (
                existing is not None
                and existing["member_manifest_hash"] == manifest_hash
                and existing["reducer_config_hash"] == config_hash
                and existing["status"] == "active"
            ):
                unchanged += 1
                continue
            next_version = int(existing["version"]) + 1 if existing else 1
            medoid = next(
                candidate
                for candidate in candidates
                if str(candidate["id"]) == medoid_id
            )
            counts = {
                "occurrences": len(candidates),
                "genomes": len(
                    {str(candidate["genome_id"]) for candidate in candidates}
                ),
                "units": len({str(candidate["unit_id"]) for candidate in candidates}),
                "tasks": len(
                    {
                        str(candidate["task_id"])
                        for candidate in candidates
                        if candidate.get("task_id") is not None
                    }
                ),
                "roles": role_counts,
                "strata": _categorical_counts(candidates, adapter),
            }
            summary = {
                "signature": medoid["signature"],
                "signature_hash": signature_hash,
                "representative_candidate_id": medoid_id,
                "representative_evidence_bundle_hash": medoid[
                    "evidence_bundle_hash"
                ],
                "equivalence": "exact_typed_signature",
            }
            cluster_id = store.create_candidate_cluster(
                campaign_id=campaign_id,
                dataset_id=group_dataset,
                candidate_type=group_type,
                signature_schema=signature_schema,
                member_ids=member_ids,
                reducer_name=EXACT_REDUCER_NAME,
                reducer_version=EXACT_REDUCER_VERSION,
                reducer_config_hash=config_hash,
                summary=summary,
                counts=counts,
                roles=roles,
                logical_cluster_id=logical_id,
                version=next_version,
                idempotency_key=(
                    f"reduce:{logical_id}:{next_version}:{manifest_hash}:{config_hash}"
                ),
            )
            created += 1
            if existing is not None:
                store.supersede_candidate_cluster(str(existing["id"]), cluster_id)
                versioned += 1
        return ReductionResult(
            campaign_id=campaign_id,
            reducer_name=EXACT_REDUCER_NAME,
            reducer_version=EXACT_REDUCER_VERSION,
            candidate_groups_seen=groups_seen,
            clusters_created=created,
            clusters_unchanged=unchanged,
            clusters_versioned=versioned,
            occurrences_represented=represented,
        )


__all__ = [
    "EXACT_REDUCER_NAME",
    "EXACT_REDUCER_VERSION",
    "CandidateReductionAdapter",
    "DefaultReductionAdapter",
    "ExactCandidateReducer",
    "ReductionResult",
]
