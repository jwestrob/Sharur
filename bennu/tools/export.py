"""export tool (Part 4.9)."""

from __future__ import annotations

from io import StringIO
from typing import Optional

from pydantic import BaseModel, Field

from bennu.core.session import ExplorationSession
from bennu.tools.registry import ToolResult


class ExportParams(BaseModel):
    format: str = Field(default="tsv", description="fasta|gff|tsv|json|notebook")
    source: str = Field(default="set", description="set|query|loci")
    source_id: Optional[str] = None
    include_sequences: bool = True
    include_annotations: bool = True


class ExportTool:
    name = "export"
    description = "Export working sets or last query results."
    tier: int = 1

    def execute(self, params: ExportParams, session: ExplorationSession) -> ToolResult:
        # Identify protein IDs
        if params.source == "set":
            ws = session.get_set(params.source_id or "")
            if not ws:
                return ToolResult(
                    success=False,
                    data=None,
                    summary="Set not found for export.",
                    count=0,
                    warnings=["Invalid set name"],
                    suggestions=[],
                )
            protein_ids = ws.protein_ids
        else:
            # fallback: export all proteins limited to 100 for safety
            rows = session.db.execute("SELECT protein_id FROM proteins LIMIT 100")
            protein_ids = [r[0] for r in rows]

        proteins = session.db.get_proteins(protein_ids)
        annotations = {}
        if params.include_annotations:
            for pid in protein_ids:
                annotations[pid] = session.db.get_annotations(pid)

        fmt = params.format.lower()
        if fmt == "json":
            import json

            data = [
                {
                    "protein": p.model_dump(),
                    "annotations": [a.model_dump() for a in annotations.get(p.protein_id, [])],
                }
                for p in proteins
            ]
            payload = json.dumps(data, indent=2)
        elif fmt == "tsv":
            buf = StringIO()
            header = ["protein_id", "contig_id", "bin_id", "start", "end", "strand"]
            buf.write("\t".join(header) + "\n")
            for p in proteins:
                buf.write(
                    "\t".join(
                        [
                            p.protein_id,
                            p.contig_id,
                            p.bin_id or "",
                            str(p.start),
                            str(p.end),
                            p.strand.value if hasattr(p.strand, "value") else str(p.strand),
                        ]
                    )
                    + "\n"
                )
            payload = buf.getvalue()
        elif fmt == "fasta":
            buf = StringIO()
            for p in proteins:
                if not p.sequence:
                    continue
                buf.write(f">{p.protein_id}\n{p.sequence}\n")
            payload = buf.getvalue()
        elif fmt == "notebook":
            payload = session.to_notebook()
        else:
            return ToolResult(
                success=False,
                data=None,
                summary=f"Unsupported export format '{params.format}'",
                count=0,
                warnings=[f"Unsupported format {params.format}"],
                suggestions=["Use one of: fasta, gff, tsv, json, notebook"],
            )

        summary = f"Exported {len(proteins)} proteins to {fmt}"
        return ToolResult(
            success=True,
            data=payload,
            summary=summary,
            count=len(proteins),
            truncated=False,
            warnings=[],
            suggestions=[],
        )


__all__ = ["ExportTool", "ExportParams"]
