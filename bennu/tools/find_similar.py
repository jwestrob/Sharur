"""find_similar tool (Part 4.5)."""

from __future__ import annotations

from time import perf_counter
from typing import List, Optional

from pydantic import BaseModel, Field

from bennu.core.session import ExplorationSession
from bennu.tools.registry import ToolResult


class FindSimilarParams(BaseModel):
    query_id: Optional[str] = None
    level: str = Field(default="protein", description="protein|locus")
    n: int = 20
    threshold: float = 0.0
    exclude_same_bin: bool = False
    taxonomy_filter: Optional[str] = None


class FindSimilarTool:
    name = "find_similar"
    description = "Find proteins or loci similar by embeddings."
    tier: int = 2

    def execute(self, params: FindSimilarParams, session: ExplorationSession) -> ToolResult:
        start = perf_counter()
        db = session.db
        warnings: List[str] = []

        pid = params.query_id or session.resolve_reference("it")
        if not pid:
            return ToolResult(
                success=False,
                data=None,
                summary="No anchor protein specified.",
                count=0,
                warnings=["query_id missing"],
                suggestions=["Provide a protein id to compare against."],
            )

        vs = session.vector_store
        if vs is None:
            warnings.append("No vector store configured; returning empty results.")
            return ToolResult(
                success=True, data=[], summary="Vector search unavailable.", count=0, warnings=warnings, suggestions=[]
            )

        neighbors = vs.query(pid, k=params.n, threshold=params.threshold or None, include_distances=True)

        if params.exclude_same_bin or params.taxonomy_filter:
            filtered = []
            for nid, score in neighbors:
                row = db.execute(
                    "SELECT bin_id FROM proteins WHERE protein_id = ?", [nid]
                )
                bin_id = row[0][0] if row else None
                if params.exclude_same_bin and bin_id:
                    anchor_bin = db.execute("SELECT bin_id FROM proteins WHERE protein_id = ?", [pid])
                    if anchor_bin and anchor_bin[0][0] == bin_id:
                        continue
                filtered.append((nid, score))
            neighbors = filtered

        protein_rows = db.get_proteins([nid for nid, _ in neighbors]) if neighbors else []
        proteins_by_id = {p.protein_id: p for p in protein_rows}
        data = [
            {"protein": proteins_by_id.get(nid), "similarity": float(score)}
            for nid, score in neighbors
            if nid in proteins_by_id
        ]

        summary = f"Found {len(data)} similar proteins to {pid}"
        duration_ms = int((perf_counter() - start) * 1000)
        session.log_query(
            query=f"find similar to {pid}",
            tool_calls=[{"tool": self.name, "params": params.model_dump()}],
            results_summary=summary,
            duration_ms=duration_ms,
        )

        return ToolResult(
            success=True,
            data=data,
            summary=summary,
            count=len(data),
            truncated=len(data) >= params.n,
            warnings=warnings,
            suggestions=[],
        )


__all__ = ["FindSimilarTool", "FindSimilarParams"]
