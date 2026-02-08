"""
Export operators for extracting sequences and data from Bennu.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Optional

from bennu.operators.base import BennuResult, OperatorContext

if TYPE_CHECKING:
    from bennu.storage.duckdb_store import DuckDBStore


def export_fasta(
    store: "DuckDBStore",
    protein_ids: list[str],
    output_path: Optional[str] = None,
    include_annotations: bool = False,
) -> BennuResult:
    """
    Export proteins as FASTA format.

    Args:
        store: DuckDB store
        protein_ids: List of protein IDs to export
        output_path: Optional file path to write FASTA
        include_annotations: Include top annotation in header

    Returns:
        BennuResult with FASTA content
    """
    params = {
        "protein_ids": protein_ids,
        "output_path": output_path,
        "include_annotations": include_annotations,
    }

    with OperatorContext("export_fasta", params) as ctx:
        if not protein_ids:
            return ctx.make_result(
                data="No protein IDs provided",
                rows=0,
            )

        # Get sequences
        placeholders = ",".join(["?"] * len(protein_ids))
        query = f"""
            SELECT protein_id, sequence
            FROM proteins
            WHERE protein_id IN ({placeholders})
        """
        rows = store.execute(query, protein_ids)

        if not rows:
            return ctx.make_result(
                data="No sequences found for provided protein IDs",
                rows=0,
            )

        # Get annotations if requested
        annotations = {}
        if include_annotations:
            ann_query = f"""
                SELECT protein_id, name, accession
                FROM annotations
                WHERE protein_id IN ({placeholders})
                ORDER BY evalue NULLS LAST
            """
            ann_rows = store.execute(ann_query, protein_ids)
            for pid, name, acc in ann_rows:
                if pid not in annotations:
                    annotations[pid] = f"{name} ({acc})" if name else acc

        # Build FASTA
        fasta_lines = []
        for protein_id, sequence in rows:
            if sequence:
                header = f">{protein_id}"
                if include_annotations and protein_id in annotations:
                    header += f" {annotations[protein_id]}"
                fasta_lines.append(header)
                # Wrap sequence at 80 characters
                for i in range(0, len(sequence), 80):
                    fasta_lines.append(sequence[i:i+80])

        fasta_content = "\n".join(fasta_lines)

        # Write to file if path provided
        if output_path:
            Path(output_path).write_text(fasta_content)
            summary = f"Exported {len(rows)} sequences to {output_path}"
        else:
            summary = f"# FASTA Export ({len(rows)} sequences)\n\n{fasta_content}"

        return ctx.make_result(
            data=summary,
            rows=len(rows),
            raw={"fasta": fasta_content, "count": len(rows)},
        )


def export_neighborhood_fasta(
    store: "DuckDBStore",
    protein_id: str,
    window: int = 10,
    output_path: Optional[str] = None,
) -> BennuResult:
    """
    Export genomic neighborhood as multi-FASTA.

    Args:
        store: DuckDB store
        protein_id: Center protein ID
        window: Number of genes on each side
        output_path: Optional file path to write FASTA

    Returns:
        BennuResult with FASTA content
    """
    params = {
        "protein_id": protein_id,
        "window": window,
        "output_path": output_path,
    }

    with OperatorContext("export_neighborhood_fasta", params) as ctx:
        # Get anchor protein info
        anchor = store.execute(
            "SELECT contig_id, gene_index FROM proteins WHERE protein_id = ?",
            [protein_id],
        )
        if not anchor:
            return ctx.make_result(
                data=f"Protein {protein_id} not found",
                rows=0,
            )

        contig_id, gene_index = anchor[0]

        # Get neighborhood proteins
        query = """
            SELECT p.protein_id, p.sequence, p.gene_index, p.strand,
                   (SELECT a.name || ' (' || a.accession || ')'
                    FROM annotations a
                    WHERE a.protein_id = p.protein_id
                    ORDER BY a.evalue NULLS LAST
                    LIMIT 1) as annotation
            FROM proteins p
            WHERE p.contig_id = ?
              AND p.gene_index BETWEEN ? AND ?
            ORDER BY p.gene_index
        """
        rows = store.execute(
            query,
            [contig_id, gene_index - window, gene_index + window],
        )

        if not rows:
            return ctx.make_result(
                data=f"No proteins found in neighborhood of {protein_id}",
                rows=0,
            )

        # Build FASTA with position info
        fasta_lines = []
        for pid, sequence, order, strand, annotation in rows:
            if sequence:
                position = order - gene_index
                pos_str = f"+{position}" if position >= 0 else str(position)
                marker = " [ANCHOR]" if pid == protein_id else ""
                header = f">{pid} pos={pos_str} strand={strand}{marker}"
                if annotation:
                    header += f" {annotation}"
                fasta_lines.append(header)
                for i in range(0, len(sequence), 80):
                    fasta_lines.append(sequence[i:i+80])

        fasta_content = "\n".join(fasta_lines)

        if output_path:
            Path(output_path).write_text(fasta_content)
            summary = f"Exported {len(rows)} sequences from neighborhood to {output_path}"
        else:
            summary = f"# Neighborhood FASTA ({len(rows)} sequences, window={window})\n\n{fasta_content}"

        return ctx.make_result(
            data=summary,
            rows=len(rows),
            raw={"fasta": fasta_content, "count": len(rows), "contig_id": contig_id},
        )


def get_sequence(
    store: "DuckDBStore",
    protein_id: str,
) -> BennuResult:
    """
    Get raw sequence for a single protein.

    Args:
        store: DuckDB store
        protein_id: Protein ID

    Returns:
        BennuResult with sequence
    """
    with OperatorContext("get_sequence", {"protein_id": protein_id}) as ctx:
        rows = store.execute(
            "SELECT sequence, sequence_length FROM proteins WHERE protein_id = ?",
            [protein_id],
        )

        if not rows or not rows[0][0]:
            return ctx.make_result(
                data=f"No sequence found for {protein_id}",
                rows=0,
            )

        sequence, length = rows[0]

        return ctx.make_result(
            data=sequence,
            rows=1,
            raw={"protein_id": protein_id, "sequence": sequence, "length": length},
        )


__all__ = ["export_fasta", "export_neighborhood_fasta", "get_sequence"]
