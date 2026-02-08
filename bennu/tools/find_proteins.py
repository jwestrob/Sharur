"""find_proteins tool (Part 4.2)."""

from __future__ import annotations

from time import perf_counter
import sys
from typing import List, Optional, Dict, Set

from pydantic import BaseModel, Field

from bennu.core.session import ExplorationSession
from bennu.core.types import Protein
from bennu.tools.registry import Tool, ToolResult

_PROTEIN_COLS = "protein_id, contig_id, bin_id, start, end_coord, strand, sequence"


def _rows_to_proteins(rows: list[tuple]) -> list[Protein]:
    return [
        Protein(
            protein_id=r[0], contig_id=r[1], bin_id=r[2],
            start=r[3], end=r[4], strand=r[5], sequence=r[6],
        )
        for r in rows
    ]


class FindProteinsParams(BaseModel):
    domains: List[str] = Field(default_factory=list, description="Pfam accessions")
    functions: List[str] = Field(default_factory=list, description="Function keywords to search in annotations")
    kegg: List[str] = Field(default_factory=list, description="KEGG ortholog IDs")
    taxonomy: Optional[str] = None
    in_set: Optional[str] = None
    exclude_set: Optional[str] = None
    limit: Optional[int] = Field(default=None, description="Optional max results; None means no cap")


class FindProteinsTool:
    name = "find_proteins"
    description = "Search proteins by domain, function keyword, KEGG, and taxonomy."
    tier: int = 2

    def execute(self, params: FindProteinsParams, session: ExplorationSession) -> ToolResult:
        start_time = perf_counter()
        store = session.db

        applied_filters: list[str] = []
        exclude_ids: Set[str] = set()

        # Resolve set constraints
        seed_ids: Optional[list[str]] = None
        if params.in_set:
            ws = session.get_set(params.in_set)
            if ws:
                seed_ids = list(ws.protein_ids)
                applied_filters.append(f"in_set={params.in_set}")
        if params.exclude_set:
            ws_ex = session.get_set(params.exclude_set)
            if ws_ex:
                exclude_ids = set(ws_ex.protein_ids)
                applied_filters.append(f"exclude_set={params.exclude_set}")

        # Build WHERE clauses for hard constraints
        where_parts: list[str] = []
        query_params: list = []

        if seed_ids is not None:
            placeholders = ",".join(["?"] * len(seed_ids))
            where_parts.append(f"p.protein_id IN ({placeholders})")
            query_params.extend(seed_ids)

        if params.taxonomy:
            where_parts.append("EXISTS (SELECT 1 FROM contigs c WHERE c.contig_id = p.contig_id AND c.taxonomy LIKE ?)")
            query_params.append(f"%{params.taxonomy}%")
            applied_filters.append(f"taxonomy~{params.taxonomy}")

        # Collect annotation-based modality subqueries (OR across modalities)
        modality_subqueries: list[tuple[str, str]] = []

        if params.domains:
            domain_placeholders = ",".join(["?"] * len(params.domains))
            modality_subqueries.append((
                "pfam",
                f"SELECT DISTINCT protein_id FROM annotations WHERE source = 'pfam' AND accession IN ({domain_placeholders})",
            ))
            query_params.extend(params.domains)
            applied_filters.append(f"pfam={params.domains}")

        if params.kegg:
            kegg_placeholders = ",".join(["?"] * len(params.kegg))
            modality_subqueries.append((
                "kegg",
                f"SELECT DISTINCT protein_id FROM annotations WHERE source = 'kegg' AND accession IN ({kegg_placeholders})",
            ))
            query_params.extend(params.kegg)
            applied_filters.append(f"kegg={params.kegg}")

        if params.functions:
            terms = []
            for fn in params.functions:
                terms.append(fn)
                if fn.endswith("ses") or fn.endswith("ases"):
                    terms.append(fn[:-1])
                elif fn.endswith("es"):
                    terms.append(fn[:-2])
                elif fn.endswith("s"):
                    terms.append(fn[:-1])
            terms = list(dict.fromkeys(terms))

            like_clauses = []
            for term in terms:
                like_clauses.append("(LOWER(a.name) LIKE LOWER(?) OR LOWER(a.description) LIKE LOWER(?))")
                query_params.extend([f"%{term}%", f"%{term}%"])

            modality_subqueries.append((
                "functions",
                f"SELECT DISTINCT a.protein_id FROM annotations a WHERE {' OR '.join(like_clauses)}",
            ))
            applied_filters.append(f"functions={params.functions}")

        # Combine: if modalities exist, union them as an IN constraint
        if modality_subqueries:
            union_sql = " UNION ".join(sq for _, sq in modality_subqueries)
            where_parts.append(f"p.protein_id IN ({union_sql})")

        where_clause = " AND ".join(where_parts) if where_parts else "1=1"
        sql = f"SELECT {_PROTEIN_COLS} FROM proteins p WHERE {where_clause}"

        rows = store.execute(sql, query_params) if query_params else store.execute(sql)
        proteins = _rows_to_proteins(rows)

        pre_limit_count = len(proteins)
        if exclude_ids:
            proteins = [p for p in proteins if p.protein_id not in exclude_ids]
        post_exclude_count = len(proteins)
        truncated = False
        if params.limit is not None and len(proteins) > params.limit:
            proteins = proteins[: params.limit]
            truncated = True

        duration_ms = int((perf_counter() - start_time) * 1000)
        summary = f"Found {len(proteins)} proteins"
        if params.taxonomy:
            summary += f" in {params.taxonomy}"
        if params.domains:
            summary += f" with domains {params.domains}"
        debug_line = (
            f"[debug] find_proteins filters={applied_filters or 'none'} "
            f"count_before_limit={pre_limit_count} after_exclude={post_exclude_count} "
            f"returned={len(proteins)} truncated={truncated}"
        )
        print(debug_line, file=sys.stderr)

        # Update session focus/provenance
        if proteins:
            session.push_focus("protein", proteins[0].protein_id)
            session.push_focus("set", f"last_results_{len(session.get_provenance())+1}")
        session.log_query(
            query=summary,
            tool_calls=[{"tool": self.name, "params": params.model_dump()}],
            results_summary=summary,
            duration_ms=duration_ms,
        )

        return ToolResult(
            success=True,
            data=proteins,
            summary=summary,
            count=len(proteins),
            truncated=truncated,
            warnings=[],
            suggestions=[],
        )


__all__ = ["FindProteinsTool", "FindProteinsParams"]
