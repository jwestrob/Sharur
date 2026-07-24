"""Exact campaign-level funnel, audit, and queue metrics for review operations."""

from __future__ import annotations

import time
from typing import TYPE_CHECKING, Any


if TYPE_CHECKING:
    from sharur.ops.store import OpsStore


def _grouped(
    store: OpsStore,
    sql: str,
    params: tuple[Any, ...],
) -> dict[str, int]:
    return {
        str(row[0]): int(row[1])
        for row in store._conn.execute(sql, params).fetchall()
    }


def review_campaign_metrics(
    store: OpsStore,
    campaign_id: str,
) -> dict[str, Any]:
    """Compute exact review-funnel metrics from normalized Ops tables."""

    campaign = store.get_campaign(campaign_id)
    scan_status = store.review_campaign_scan_status(campaign_id)
    with store._lock:
        dispositions = _grouped(
            store,
            """
            SELECT disposition, COUNT(*)
            FROM unit_dispositions
            WHERE campaign_id = ? AND record_status = 'active'
            GROUP BY disposition
            """,
            (campaign_id,),
        )
        reviews_by_tier = _grouped(
            store,
            """
            SELECT review_tier, COUNT(*)
            FROM finding_reviews
            WHERE campaign_id = ?
            GROUP BY review_tier
            """,
            (campaign_id,),
        )
        reviews_by_verdict = _grouped(
            store,
            """
            SELECT verdict, COUNT(*)
            FROM finding_reviews
            WHERE campaign_id = ?
            GROUP BY verdict
            """,
            (campaign_id,),
        )
        reviews_by_provider = _grouped(
            store,
            """
            SELECT provider, COUNT(*)
            FROM finding_reviews
            WHERE campaign_id = ?
            GROUP BY provider
            """,
            (campaign_id,),
        )
        decisions = _grouped(
            store,
            """
            SELECT decision, COUNT(*)
            FROM promotion_decisions
            WHERE campaign_id = ?
            GROUP BY decision
            """,
            (campaign_id,),
        )
        queue_rows = store._conn.execute(
            """
            SELECT
                COALESCE(execution_profile, 'unprofiled') AS profile,
                status,
                COUNT(*) AS count
            FROM tasks
            WHERE campaign_id = ?
            GROUP BY profile, status
            ORDER BY profile, status
            """,
            (campaign_id,),
        ).fetchall()
        queue: dict[str, dict[str, int]] = {}
        for row in queue_rows:
            queue.setdefault(str(row["profile"]), {})[str(row["status"])] = int(
                row["count"]
            )
        totals = store._conn.execute(
            """
            SELECT
                (SELECT COUNT(*) FROM unit_dispositions
                 WHERE campaign_id = ? AND record_status = 'active') AS units,
                (SELECT COUNT(*) FROM candidate_occurrences
                 WHERE campaign_id = ?) AS candidates,
                (SELECT COUNT(*) FROM candidate_clusters
                 WHERE campaign_id = ? AND status = 'active') AS active_clusters,
                (SELECT COUNT(*) FROM cluster_findings AS cf
                 JOIN candidate_clusters AS c ON c.id = cf.cluster_id
                 WHERE c.campaign_id = ? AND cf.relation = 'materializes')
                    AS materialized_findings,
                (SELECT COUNT(*) FROM canonical_publications
                 WHERE campaign_id = ?) AS publications
            """,
            (campaign_id,) * 5,
        ).fetchone()
        verification_rows = store._conn.execute(
            """
            WITH latest AS (
                SELECT
                    v.status,
                    ROW_NUMBER() OVER (
                        PARTITION BY
                            v.review_id, v.claim_key, v.specification_hash
                        ORDER BY v.attempt DESC, v.created_ts DESC, v.id DESC
                    ) AS rank
                FROM review_verifications AS v
                JOIN finding_reviews AS r ON r.id = v.review_id
                WHERE r.campaign_id = ?
            )
            SELECT status, COUNT(*) FROM latest
            WHERE rank = 1
            GROUP BY status
            """,
            (campaign_id,),
        ).fetchall()
        latest_verifications = {
            str(row[0]): int(row[1]) for row in verification_rows
        }
        audit = store._conn.execute(
            """
            SELECT
                (
                    SELECT COUNT(*) FROM tasks
                    WHERE campaign_id = ?
                      AND json_extract(params, '$.review_tier') = 'audit_clear'
                ) AS clear_tasks,
                (
                    SELECT COUNT(*) FROM finding_reviews
                    WHERE campaign_id = ? AND review_tier = 'audit_clear'
                ) AS clear_reviews,
                (
                    SELECT COUNT(*) FROM finding_reviews
                    WHERE campaign_id = ? AND review_tier = 'audit_clear'
                      AND verdict <> 'promote'
                ) AS clear_contrary,
                (
                    SELECT COUNT(*) FROM tasks
                    WHERE campaign_id = ?
                      AND json_extract(params, '$.review_tier') = 'audit_reject'
                ) AS rejected_tasks,
                (
                    SELECT COUNT(*) FROM finding_reviews
                    WHERE campaign_id = ? AND review_tier = 'audit_reject'
                ) AS rejected_reviews,
                (
                    SELECT COUNT(*) FROM promotion_decisions
                    WHERE campaign_id = ? AND source_tier = 'audit_reject'
                      AND decision = 'reopen'
                ) AS rejected_reopened
            """,
            (campaign_id,) * 6,
        ).fetchone()
        disagreement_count = int(
            store._conn.execute(
                """
                SELECT COUNT(*) FROM (
                    SELECT cluster_id
                    FROM finding_reviews
                    WHERE campaign_id = ?
                      AND review_tier = 'independent'
                      AND cluster_id IS NOT NULL
                    GROUP BY cluster_id
                    HAVING COUNT(DISTINCT verdict) > 1
                )
                """,
                (campaign_id,),
            ).fetchone()[0]
        )
    unit_count = int(totals["units"])
    candidate_count = int(totals["candidates"])
    cluster_count = int(totals["active_clusters"])
    clear_units = int(dispositions.get("clear", 0))
    clear_tasks = int(audit["clear_tasks"])
    return {
        "campaign_id": campaign_id,
        "campaign_status": campaign["status"],
        "generated_ts": time.time(),
        "coverage": {
            "active_unit_dispositions": unit_count,
            "by_disposition": dispositions,
            "atlas_closure": scan_status,
        },
        "reduction": {
            "candidate_occurrences": candidate_count,
            "active_clusters": cluster_count,
            "occurrences_per_active_cluster": (
                candidate_count / cluster_count if cluster_count else None
            ),
        },
        "review": {
            "by_tier": reviews_by_tier,
            "by_verdict": reviews_by_verdict,
            "by_provider": reviews_by_provider,
            "independent_disagreement_subjects": disagreement_count,
            "latest_verifications_by_status": latest_verifications,
        },
        "audit": {
            "clear_units": clear_units,
            "clear_tasks": clear_tasks,
            "clear_realized_rate": (
                clear_tasks / clear_units if clear_units else None
            ),
            "clear_reviews": int(audit["clear_reviews"]),
            "clear_contrary_reviews": int(audit["clear_contrary"]),
            "rejected_tasks": int(audit["rejected_tasks"]),
            "rejected_reviews": int(audit["rejected_reviews"]),
            "rejected_reopened": int(audit["rejected_reopened"]),
        },
        "promotion": {
            "by_decision": decisions,
            "materialized_findings": int(totals["materialized_findings"]),
            "canonical_publications": int(totals["publications"]),
        },
        "queue": queue,
    }


__all__ = ["review_campaign_metrics"]
