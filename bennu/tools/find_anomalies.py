"""find_anomalies tool (Part 4.6)."""

from __future__ import annotations

from time import perf_counter
from typing import List, Optional

from pydantic import BaseModel, Field

from bennu.core.session import ExplorationSession
from bennu.tools.registry import ToolResult


class FindAnomaliesParams(BaseModel):
    scope: str = Field(default="genome", description="genome|bin|set")
    scope_id: Optional[str] = None
    anomaly_types: List[str] = Field(default_factory=list)
    signals: List[str] = Field(default_factory=list)
    min_score: float = 0.7
    limit: int = 50


class FindAnomaliesTool:
    name = "find_anomalies"
    description = "Surface statistically unusual proteins using precomputed feature_store."
    tier: int = 2

    def execute(self, params: FindAnomaliesParams, session: ExplorationSession) -> ToolResult:
        start = perf_counter()
        db = session.db
        warnings: List[str] = []

        # Determine scope filter
        where_clauses = []
        sql_params: list = []
        if params.scope == "bin" and params.scope_id:
            where_clauses.append("p.bin_id = ?")
            sql_params.append(params.scope_id)
        elif params.scope == "set" and params.scope_id:
            # scope_id for set is WorkingSet set_id hex or name; try both
            ws = session.get_set(params.scope_id) or session.get_set(params.scope_id.lower()) if hasattr(session, "get_set") else None
            if ws:
                placeholders = ",".join(["?"] * len(ws.protein_ids))
                where_clauses.append(f"p.protein_id IN ({placeholders})")
                sql_params.extend(ws.protein_ids)
            else:
                warnings.append("Working set not found; returning empty results.")
        # anomaly types / signals map to metric_name patterns
        metric_clause = ""
        if params.signals:
            placeholders = ",".join(["?"] * len(params.signals))
            metric_clause = f"AND f.metric_name IN ({placeholders})"
            sql_params.extend(params.signals)

        where_sql = "WHERE " + " AND ".join(where_clauses) if where_clauses else ""

        query = f"""
            SELECT p.protein_id, p.bin_id, p.contig_id, f.metric_name, f.metric_value
            FROM feature_store f
            JOIN proteins p ON p.protein_id = f.protein_id
            {where_sql}
            {metric_clause}
            ORDER BY f.metric_value DESC
            LIMIT {int(params.limit)}
        """
        try:
            rows = db.execute(query, sql_params if sql_params else None)
        except Exception as exc:
            return ToolResult(
                success=False,
                data=None,
                summary=f"find_anomalies failed: {exc}",
                count=0,
                warnings=[str(exc)],
                suggestions=[],
            )

        results = [
            {"protein_id": r[0], "bin_id": r[1], "contig_id": r[2], "signal": r[3], "score": r[4]}
            for r in rows
            if r[4] is None or r[4] >= params.min_score
        ]

        summary = f"Found {len(results)} anomalous proteins"
        duration_ms = int((perf_counter() - start) * 1000)
        session.log_query(
            query="find anomalies",
            tool_calls=[{"tool": self.name, "params": params.model_dump()}],
            results_summary=summary,
            duration_ms=duration_ms,
            error=None,
        )

        return ToolResult(
            success=True,
            data=results,
            summary=summary,
            count=len(results),
            truncated=len(results) >= params.limit,
            warnings=warnings,
            suggestions=[],
        )


__all__ = ["FindAnomaliesTool", "FindAnomaliesParams"]
