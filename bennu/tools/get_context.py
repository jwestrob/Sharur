"""get_context tool (Part 4.3)."""

from __future__ import annotations

from time import perf_counter
from typing import List, Optional

from pydantic import BaseModel, Field

from bennu.core.session import ExplorationSession
from bennu.core.types import Protein
from bennu.tools.registry import Tool, ToolResult


class GetContextParams(BaseModel):
    protein_id: Optional[str] = None
    window_genes: int = 10
    window_bp: int = 0
    include_annotations: bool = True
    expand_operon: bool = False


class GetContextTool:
    name = "get_context"
    description = "Retrieve genomic neighborhood around a protein."
    tier: int = 2

    def execute(self, params: GetContextParams, session: ExplorationSession) -> ToolResult:
        start = perf_counter()

        pid = params.protein_id or session.resolve_reference("it")
        if not pid:
            return ToolResult(
                success=False,
                data=None,
                summary="No protein specified or resolvable from context.",
                count=0,
                warnings=["protein_id missing"],
                suggestions=["Specify a protein id."],
            )

        store = session.db
        proteins = self._get_neighborhood(store, pid, params)
        proteins = sorted(proteins, key=lambda p: (p.contig_id, p.start))

        summary = f"Context around {pid}: {len(proteins)} genes"
        duration_ms = int((perf_counter() - start) * 1000)

        # Update session focus and provenance
        session.push_focus("protein", pid)
        session.push_focus("set", f"context_{pid}")
        session.log_query(
            query=f"get context for {pid}",
            tool_calls=[{"tool": self.name, "params": params.model_dump()}],
            results_summary=summary,
            duration_ms=duration_ms,
        )

        return ToolResult(
            success=True,
            data=proteins,
            summary=summary,
            count=len(proteins),
            truncated=False,
            warnings=[],
            suggestions=[],
        )

    @staticmethod
    def _get_neighborhood(store, pid: str, params: GetContextParams) -> List[Protein]:
        """Retrieve proteins in neighborhood using direct store queries."""
        # Get anchor protein's contig
        rows = store.execute(
            "SELECT protein_id, contig_id, bin_id, start, end_coord, strand, sequence "
            "FROM proteins WHERE protein_id = ?",
            [pid],
        )
        if not rows:
            return []
        anchor = Protein(
            protein_id=rows[0][0], contig_id=rows[0][1], bin_id=rows[0][2],
            start=rows[0][3], end=rows[0][4], strand=rows[0][5], sequence=rows[0][6],
        )

        if params.window_bp > 0:
            # Coordinate-based window
            w_start = max(0, anchor.start - params.window_bp)
            w_end = anchor.end + params.window_bp
            rows = store.execute(
                "SELECT protein_id, contig_id, bin_id, start, end_coord, strand, sequence "
                "FROM proteins WHERE contig_id = ? AND end_coord >= ? AND start <= ? "
                "ORDER BY start",
                [anchor.contig_id, w_start, w_end],
            )
        else:
            # Gene-count window
            contig_rows = store.execute(
                "SELECT protein_id, contig_id, bin_id, start, end_coord, strand, sequence "
                "FROM proteins WHERE contig_id = ? ORDER BY start",
                [anchor.contig_id],
            )
            idx = next((i for i, r in enumerate(contig_rows) if r[0] == pid), None)
            if idx is None:
                return [anchor]
            start_idx = max(0, idx - params.window_genes)
            end_idx = min(len(contig_rows), idx + params.window_genes + 1)
            rows = contig_rows[start_idx:end_idx]

        proteins = [
            Protein(
                protein_id=r[0], contig_id=r[1], bin_id=r[2],
                start=r[3], end=r[4], strand=r[5], sequence=r[6],
            )
            for r in rows
        ]

        if params.expand_operon:
            proteins = GetContextTool._expand_operon(store, proteins, anchor)

        return proteins

    @staticmethod
    def _expand_operon(store, proteins: List[Protein], anchor: Protein, gap_bp: int = 300) -> List[Protein]:
        """Expand to include co-transcribed genes on same strand."""
        contig_rows = store.execute(
            "SELECT protein_id, contig_id, bin_id, start, end_coord, strand, sequence "
            "FROM proteins WHERE contig_id = ? ORDER BY start",
            [anchor.contig_id],
        )
        all_proteins = [
            Protein(
                protein_id=r[0], contig_id=r[1], bin_id=r[2],
                start=r[3], end=r[4], strand=r[5], sequence=r[6],
            )
            for r in contig_rows
        ]
        idx = next((i for i, p in enumerate(all_proteins) if p.protein_id == anchor.protein_id), None)
        if idx is None:
            return proteins

        operon_ids = {p.protein_id for p in proteins}

        # Extend left
        i = idx - 1
        while i >= 0:
            p = all_proteins[i]
            if p.strand != anchor.strand:
                break
            if all_proteins[i + 1].start - p.end > gap_bp:
                break
            operon_ids.add(p.protein_id)
            i -= 1

        # Extend right
        i = idx + 1
        while i < len(all_proteins):
            p = all_proteins[i]
            if p.strand != anchor.strand:
                break
            if p.start - all_proteins[i - 1].end > gap_bp:
                break
            operon_ids.add(p.protein_id)
            i += 1

        return [p for p in all_proteins if p.protein_id in operon_ids]


__all__ = ["GetContextTool", "GetContextParams"]
