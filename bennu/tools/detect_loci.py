"""detect_loci tool (Part 4.4)."""

from __future__ import annotations

from time import perf_counter
from typing import List, Optional

from pydantic import BaseModel, Field

from bennu.core.session import ExplorationSession
from bennu.core.types import LocusType
from bennu.tools.registry import ToolResult


class DetectLociParams(BaseModel):
    locus_type: str = Field(default="bgc")
    scope: str = Field(default="genome", description="genome|bin|contig")
    scope_id: Optional[str] = None
    required_domains: List[str] = Field(default_factory=list)
    optional_domains: List[str] = Field(default_factory=list)


class DetectLociTool:
    name = "detect_loci"
    description = "Retrieve precomputed loci of specified type."
    tier: int = 2

    def execute(self, params: DetectLociParams, session: ExplorationSession) -> ToolResult:
        start = perf_counter()
        db = session.db
        where = ["locus_type = ?"]
        sql_params: List = [params.locus_type]

        if params.scope == "bin" and params.scope_id:
            where.append("c.bin_id = ?")
            sql_params.append(params.scope_id)
        elif params.scope == "contig" and params.scope_id:
            where.append("l.contig_id = ?")
            sql_params.append(params.scope_id)

        where_sql = "WHERE " + " AND ".join(where)
        query = f"""
            SELECT l.locus_id, l.locus_type, l.contig_id, l.start, l.end_coord, l.confidence, l.metadata,
                   lp.protein_id, lp.position
            FROM loci l
            LEFT JOIN contigs c ON l.contig_id = c.contig_id
            LEFT JOIN locus_proteins lp ON l.locus_id = lp.locus_id
            {where_sql}
            ORDER BY l.locus_id, lp.position
        """
        rows = db.execute(query, sql_params)

        loci: dict[str, dict] = {}
        for row in rows:
            locus_id, locus_type, contig_id, start, end, conf, metadata, pid, pos = row
            if locus_id not in loci:
                loci[locus_id] = {
                    "locus_id": locus_id,
                    "locus_type": locus_type,
                    "contig_id": contig_id,
                    "start": start,
                    "end": end,
                    "confidence": conf,
                    "metadata": metadata or {},
                    "proteins": [],
                }
            if pid:
                loci[locus_id]["proteins"].append(pid)

        result_list = list(loci.values())
        summary = f"Found {len(result_list)} {params.locus_type} loci"
        duration_ms = int((perf_counter() - start) * 1000)
        session.log_query(
            query=f"detect loci {params.locus_type}",
            tool_calls=[{"tool": self.name, "params": params.model_dump()}],
            results_summary=summary,
            duration_ms=duration_ms,
        )

        return ToolResult(
            success=True,
            data=result_list,
            summary=summary,
            count=len(result_list),
            truncated=False,
            warnings=[],
            suggestions=[],
        )


__all__ = ["DetectLociTool", "DetectLociParams"]
