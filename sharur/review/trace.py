"""Bounded provenance reconstruction for review-DAG subjects."""

from __future__ import annotations

from collections import Counter
from typing import Any, Protocol


class ReviewAccess(Protocol):
    def get_candidate_cluster(self, cluster_id: str) -> dict[str, Any]: ...

    def get_finding(self, finding_id: str) -> dict[str, Any]: ...

    def get_unit_disposition(self, disposition_id: str) -> dict[str, Any]: ...

    def list_cluster_findings(self, **filters: Any) -> list[dict[str, Any]]: ...

    def list_finding_reviews(self, **filters: Any) -> list[dict[str, Any]]: ...

    def list_review_verifications(
        self, review_id: str
    ) -> list[dict[str, Any]]: ...

    def list_promotion_decisions(self, **filters: Any) -> list[dict[str, Any]]: ...

    def list_canonical_publications(
        self, **filters: Any
    ) -> list[dict[str, Any]]: ...


def _review_trace(
    access: ReviewAccess,
    reviews: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    return [
        {
            **review,
            "verifications": access.list_review_verifications(str(review["id"])),
        }
        for review in reviews
    ]


def _cluster_summary(cluster: dict[str, Any]) -> dict[str, Any]:
    members = list(cluster.get("members") or [])
    member_count = int(cluster.get("member_count", len(members)))
    exact_role_counts = cluster.get("member_role_counts") or (
        cluster.get("counts") or {}
    ).get("roles")
    role_counts = (
        {str(key): int(value) for key, value in exact_role_counts.items()}
        if isinstance(exact_role_counts, dict)
        else dict(sorted(Counter(str(member["role"]) for member in members).items()))
    )
    return {
        **{
            key: value
            for key, value in cluster.items()
            if key not in {"members", "member_role_counts", "members_truncated"}
        },
        "member_count": member_count,
        "member_role_counts": role_counts,
        "member_sample": members[:20],
        "member_sample_truncated": member_count > min(len(members), 20),
    }


def trace_review_subject(
    access: ReviewAccess,
    *,
    campaign_id: str,
    subject_kind: str,
    subject_id: str,
) -> dict[str, Any]:
    """Reconstruct a bounded, machine-readable review lineage."""

    if subject_kind == "candidate_cluster":
        cluster = access.get_candidate_cluster(subject_id)
        reviews = access.list_finding_reviews(
            campaign_id=campaign_id,
            cluster_id=subject_id,
            limit=1_000,
        )
        links = access.list_cluster_findings(cluster_id=subject_id)
        return {
            "subject_kind": subject_kind,
            "subject_id": subject_id,
            "cluster": _cluster_summary(cluster),
            "reviews": _review_trace(access, reviews),
            "decisions": access.list_promotion_decisions(
                campaign_id=campaign_id,
                cluster_id=subject_id,
                limit=1_000,
            ),
            "finding_links": links,
            "materialized_findings": [
                access.get_finding(str(link["finding_id"])) for link in links
            ],
        }
    if subject_kind == "finding":
        finding = access.get_finding(subject_id)
        reviews = access.list_finding_reviews(
            campaign_id=campaign_id,
            finding_id=subject_id,
            limit=1_000,
        )
        links = access.list_cluster_findings(finding_id=subject_id)
        return {
            "subject_kind": subject_kind,
            "subject_id": subject_id,
            "finding": finding,
            "reviews": _review_trace(access, reviews),
            "decisions": access.list_promotion_decisions(
                campaign_id=campaign_id,
                finding_id=subject_id,
                limit=1_000,
            ),
            "cluster_links": links,
            "source_clusters": [
                _cluster_summary(
                    access.get_candidate_cluster(str(link["cluster_id"]))
                )
                for link in links
            ],
            "publications": access.list_canonical_publications(
                campaign_id=campaign_id,
                finding_id=subject_id,
                limit=1_000,
            ),
        }
    if subject_kind == "unit_disposition":
        disposition = access.get_unit_disposition(subject_id)
        reviews = access.list_finding_reviews(
            campaign_id=campaign_id,
            unit_disposition_id=subject_id,
            limit=1_000,
        )
        return {
            "subject_kind": subject_kind,
            "subject_id": subject_id,
            "unit_disposition": disposition,
            "reviews": _review_trace(access, reviews),
        }
    raise ValueError(
        "subject_kind must be candidate_cluster, finding, or unit_disposition"
    )


__all__ = ["ReviewAccess", "trace_review_subject"]
