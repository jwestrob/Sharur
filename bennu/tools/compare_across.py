"""compare_across tool (Part 4.7)."""

from __future__ import annotations

from time import perf_counter
from typing import Optional

from pydantic import BaseModel, Field

from bennu.core.session import ExplorationSession
from bennu.tools.registry import ToolResult


class CompareAcrossParams(BaseModel):
    feature: str = Field(default="domain", description="domain|function|pathway|locus_type")
    feature_id: str = ""
    group_by: str = Field(default="bin", description="bin|taxonomy|sample")
    taxonomy_level: int = 2
    metric: str = Field(default="count", description="presence|count|abundance")


class CompareAcrossTool:
    name = "compare_across"
    description = "Compare feature presence across bins or taxonomy groups."
    tier: int = 2

    def execute(self, params: CompareAcrossParams, session: ExplorationSession) -> ToolResult:
        start = perf_counter()
        db = session.db

        if params.feature == "domain":
            if params.feature_id:
                where = "WHERE a.source = 'pfam' AND a.accession = ?"
                sql_params = [params.feature_id]
            else:
                where = "WHERE a.source = 'pfam'"
                sql_params = []
            if params.group_by == "bin":
                query = f"""
                    SELECT p.bin_id as group_key, COUNT(DISTINCT a.protein_id) as cnt
                    FROM annotations a
                    JOIN proteins p ON p.protein_id = a.protein_id
                    {where}
                    GROUP BY 1
                    ORDER BY cnt DESC
                """
            else:  # taxonomy
                query = f"""
                    SELECT substr(b.taxonomy,1,100) as group_key, COUNT(DISTINCT a.protein_id) as cnt
                    FROM annotations a
                    JOIN proteins p ON p.protein_id = a.protein_id
                    LEFT JOIN bins b ON p.bin_id = b.bin_id
                    {where}
                    GROUP BY 1
                    ORDER BY cnt DESC
                """
        elif params.feature == "locus_type":
            query = """
                SELECT locus_type as group_key, COUNT(*) as cnt
                FROM loci
                GROUP BY locus_type
                ORDER BY cnt DESC
            """
            sql_params = []
        else:
            # fallback
            query = """
                SELECT 'unknown' as group_key, 0 as cnt
            """
            sql_params = []

        rows = db.execute(query, sql_params if sql_params else None)
        data = [{"group": r[0] or "unknown", "count": r[1]} for r in rows]
        summary = f"Compared {params.feature} across {len(data)} groups"
        duration_ms = int((perf_counter() - start) * 1000)
        session.log_query(
            query=f"compare {params.feature}",
            tool_calls=[{"tool": self.name, "params": params.model_dump()}],
            results_summary=summary,
            duration_ms=duration_ms,
        )

        return ToolResult(
            success=True,
            data=data,
            summary=summary,
            count=len(data),
            truncated=False,
            warnings=[],
            suggestions=[],
        )


__all__ = ["CompareAcrossTool", "CompareAcrossParams"]
