"""
Visualization operators for gene diagrams and pathway context.

Generates publication-ready gene arrow diagrams using dna_features_viewer.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Optional

from sharur.operators.base import SharurResult, OperatorContext

if TYPE_CHECKING:
    from sharur.storage.duckdb_store import DuckDBStore


# Color scheme for different annotation types
ANNOTATION_COLORS = {
    "unannotated": "#CCCCCC",       # Gray
    "hypothetical": "#AAAAAA",       # Dark gray
    "transporter": "#4CAF50",        # Green
    "oxidoreductase": "#2196F3",     # Blue
    "transferase": "#9C27B0",        # Purple
    "hydrolase": "#FF9800",          # Orange
    "kinase": "#E91E63",             # Pink
    "phage": "#F44336",              # Red
    "transposase": "#795548",        # Brown
    "crispr": "#00BCD4",             # Cyan
    "metal_binding": "#FFC107",      # Amber
    "membrane": "#8BC34A",           # Light green
    "default": "#607D8B",            # Blue gray
}


def _get_color_for_protein(annotation: Optional[str], predicates: Optional[list] = None) -> str:
    """Determine color based on annotation and predicates."""
    if not annotation or annotation == "NO HITS":
        return ANNOTATION_COLORS["unannotated"]

    ann_lower = annotation.lower()

    if "hypothetical" in ann_lower or "uncharacterized" in ann_lower:
        return ANNOTATION_COLORS["hypothetical"]
    if "transporter" in ann_lower or "permease" in ann_lower:
        return ANNOTATION_COLORS["transporter"]
    if "oxidoreductase" in ann_lower or "dehydrogenase" in ann_lower:
        return ANNOTATION_COLORS["oxidoreductase"]
    if "transferase" in ann_lower:
        return ANNOTATION_COLORS["transferase"]
    if "hydrolase" in ann_lower or "peptidase" in ann_lower:
        return ANNOTATION_COLORS["hydrolase"]
    if "kinase" in ann_lower:
        return ANNOTATION_COLORS["kinase"]
    if "phage" in ann_lower or "prophage" in ann_lower:
        return ANNOTATION_COLORS["phage"]
    if "transposase" in ann_lower or "integrase" in ann_lower:
        return ANNOTATION_COLORS["transposase"]
    if "crispr" in ann_lower or "cas" in ann_lower:
        return ANNOTATION_COLORS["crispr"]

    return ANNOTATION_COLORS["default"]


def _format_label(name: str, max_width: int = 15) -> str:
    """Format a label with word wrapping."""
    if not name or len(name) <= max_width:
        return name or "?"

    # Wrap at ~max_width chars, preserve whole words
    words = name.replace("-", "- ").replace("_", "_ ").split()
    lines = []
    current = ""
    for word in words:
        test = f"{current} {word}".strip() if current else word
        if len(test) <= max_width:
            current = test
        else:
            if current:
                lines.append(current)
            current = word
    if current:
        lines.append(current)
    return "\n".join(lines[:3])  # Max 3 lines


def visualize_neighborhood(
    store: "DuckDBStore",
    protein_id: str,
    window: int = 10,
    output_path: Optional[str] = None,
    figure_width: int = 14,
    show_labels: bool = True,
    show_crispr_arrays: bool = True,
) -> SharurResult:
    """
    Generate gene arrow diagram for a genomic neighborhood.

    Args:
        store: DuckDB store
        protein_id: Center protein ID
        window: Number of genes on each side
        output_path: Path to save PNG/SVG (None for temp file)
        figure_width: Width of figure in inches
        show_labels: Show gene labels
        show_crispr_arrays: Show CRISPR arrays from loci table

    Returns:
        SharurResult with path to generated image
    """
    params = {
        "protein_id": protein_id,
        "window": window,
        "output_path": output_path,
        "figure_width": figure_width,
    }

    with OperatorContext("visualize_neighborhood", params) as ctx:
        try:
            from dna_features_viewer import GraphicFeature, GraphicRecord
            import matplotlib
            matplotlib.use('Agg')  # Non-interactive backend
            import matplotlib.pyplot as plt
        except ImportError:
            return ctx.make_result(
                data="dna_features_viewer not installed. Run: pip install dna_features_viewer",
                rows=0,
            )

        # Get anchor protein info including its annotation
        anchor = store.execute(
            """SELECT p.contig_id, p.gene_index, p.start, p.end_coord,
                      (SELECT a.name FROM annotations a
                       WHERE a.protein_id = p.protein_id
                       ORDER BY a.evalue NULLS LAST LIMIT 1) as annotation
               FROM proteins p WHERE p.protein_id = ?""",
            [protein_id],
        )
        if not anchor:
            return ctx.make_result(
                data=f"Protein {protein_id} not found",
                rows=0,
            )

        contig_id, gene_index, anchor_start, anchor_end, anchor_annotation = anchor[0]

        # Get neighborhood proteins
        query = """
            SELECT
                p.protein_id,
                p.start,
                p.end_coord,
                p.strand,
                p.gene_index,
                (SELECT a.name
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
                data=f"No proteins found in neighborhood",
                rows=0,
            )

        # Calculate coordinate offset (normalize to start at 0)
        min_coord = min(r[1] for r in rows)
        max_coord = max(r[2] for r in rows)

        # Get CRISPR arrays in this region (if table has data)
        # Use expanded range to catch arrays in intergenic regions near the proteins
        crispr_arrays = []
        if show_crispr_arrays:
            try:
                # Expand search range by 1kb on each side to catch nearby arrays
                expanded_min = max(0, min_coord - 1000)
                expanded_max = max_coord + 1000
                crispr_arrays = store.execute(
                    """SELECT locus_id, start, end_coord, metadata
                       FROM loci
                       WHERE contig_id = ?
                         AND locus_type = 'crispr_array'
                         AND start <= ? AND end_coord >= ?""",
                    [contig_id, expanded_max, expanded_min],
                )
                # Extend view to include any found arrays
                for _, arr_start, arr_end, _ in crispr_arrays:
                    min_coord = min(min_coord, arr_start)
                    max_coord = max(max_coord, arr_end)
            except Exception:
                pass  # Table may not exist or be empty

        # Build features
        features = []

        # Add CRISPR arrays first (as background boxes)
        for locus_id, arr_start, arr_end, metadata in crispr_arrays:
            features.append(
                GraphicFeature(
                    start=max(arr_start - min_coord, 0),
                    end=min(arr_end - min_coord, max_coord - min_coord),
                    strand=0,  # No strand for arrays
                    color="#E3F2FD",  # Light blue background
                    label="CRISPR\narray",
                    linewidth=2,
                    linecolor="#1976D2",
                )
            )

        # Add protein features
        for pid, start, end, strand, order, annotation in rows:
            strand_int = 1 if strand == "+" else -1
            color = _get_color_for_protein(annotation)

            # Format label - anchor gets highlighted with its real annotation
            if pid == protein_id:
                # Use actual annotation for anchor, with marker
                if anchor_annotation and anchor_annotation != "NO HITS":
                    ann_name = anchor_annotation.split("(")[0].strip()
                    label = f">>> {_format_label(ann_name, 12)} <<<"
                else:
                    label = ">>> QUERY <<<"
                color = "#FF0000"  # Red for anchor
            elif annotation and annotation != "NO HITS":
                name = annotation.split("(")[0].strip()
                label = _format_label(name, 15)
            else:
                label = "?"

            features.append(
                GraphicFeature(
                    start=start - min_coord,
                    end=end - min_coord,
                    strand=strand_int,
                    color=color,
                    label=label if show_labels else None,
                )
            )

        # Create record and plot
        record = GraphicRecord(
            sequence_length=max_coord - min_coord,
            features=features,
        )

        # Generate figure (extra height for wrapped labels - up to 3 lines)
        fig, ax = plt.subplots(1, 1, figsize=(figure_width, 5))
        record.plot(ax=ax, with_ruler=True, strand_in_label_threshold=4)

        # Relabel x-axis ticks to show absolute genome coordinates
        ticks = ax.get_xticks()
        ax.set_xticklabels([f"{int(t + min_coord):,}" for t in ticks])

        # Build title with anchor annotation
        anchor_desc = anchor_annotation if anchor_annotation and anchor_annotation != "NO HITS" else "unannotated"
        if len(anchor_desc) > 40:
            anchor_desc = anchor_desc[:37] + "..."
        ax.set_title(
            f"Neighborhood of {protein_id}\n{anchor_desc} ({len(rows)} genes, {window}-gene window)",
            fontsize=10
        )

        # Save or return path
        if output_path is None:
            output_path = f"/tmp/bennu_neighborhood_{protein_id[:30]}.png"

        plt.tight_layout()
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()

        return ctx.make_result(
            data=f"Gene diagram saved to: {output_path}\n\nView with: open {output_path}",
            rows=len(rows),
            raw={"output_path": output_path, "gene_count": len(rows), "crispr_arrays": len(crispr_arrays)},
        )


def visualize_domain_architecture(
    store: "DuckDBStore",
    protein_id: str,
    output_path: Optional[str] = None,
) -> SharurResult:
    """
    Generate domain architecture diagram for a protein.

    Args:
        store: DuckDB store
        protein_id: Protein ID
        output_path: Path to save image

    Returns:
        SharurResult with path to generated image
    """
    params = {"protein_id": protein_id, "output_path": output_path}

    with OperatorContext("visualize_domain_architecture", params) as ctx:
        try:
            from dna_features_viewer import GraphicFeature, GraphicRecord
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
        except ImportError:
            return ctx.make_result(
                data="dna_features_viewer not installed",
                rows=0,
            )

        # Get protein length
        protein = store.execute(
            "SELECT sequence_length FROM proteins WHERE protein_id = ?",
            [protein_id],
        )
        if not protein or not protein[0][0]:
            return ctx.make_result(
                data=f"Protein {protein_id} not found or has no length",
                rows=0,
            )

        length = protein[0][0]

        # Get PFAM domain annotations with coordinates
        # Note: This assumes annotations have domain coordinates stored
        # If not available, we'll create a simple representation
        annotations = store.execute(
            """
            SELECT name, accession, description
            FROM annotations
            WHERE protein_id = ? AND source = 'pfam'
            ORDER BY evalue NULLS LAST
            """,
            [protein_id],
        )

        if not annotations:
            return ctx.make_result(
                data=f"No PFAM domains found for {protein_id}",
                rows=0,
            )

        # Create features (evenly spaced if no coordinates)
        features = []
        n_domains = len(annotations)
        segment_length = length // (n_domains + 1)

        colors = ['#2196F3', '#4CAF50', '#FF9800', '#9C27B0', '#E91E63', '#00BCD4']

        for i, (name, accession, description) in enumerate(annotations):
            start = (i + 1) * segment_length - segment_length // 2
            end = start + segment_length // 2
            features.append(
                GraphicFeature(
                    start=start,
                    end=min(end, length),
                    strand=0,
                    color=colors[i % len(colors)],
                    label=f"{name}\n{accession}",
                )
            )

        record = GraphicRecord(
            sequence_length=length,
            features=features,
        )

        fig, ax = plt.subplots(1, 1, figsize=(12, 2))
        record.plot(ax=ax, with_ruler=True)
        ax.set_title(f"Domain Architecture: {protein_id[:50]}... ({length} aa)", fontsize=10)

        if output_path is None:
            output_path = f"/tmp/bennu_domains_{protein_id[:30]}.png"

        plt.tight_layout()
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()

        return ctx.make_result(
            data=f"Domain diagram saved to: {output_path}\n\nView with: open {output_path}",
            rows=len(annotations),
            raw={"output_path": output_path, "domain_count": len(annotations)},
        )


def get_kegg_pathway_context(
    store: "DuckDBStore",
    protein_id: str,
) -> SharurResult:
    """
    Get KEGG pathway context for a protein.

    Args:
        store: DuckDB store
        protein_id: Protein ID

    Returns:
        SharurResult with KEGG pathway information
    """
    params = {"protein_id": protein_id}

    with OperatorContext("get_kegg_pathway_context", params) as ctx:
        # Get KEGG annotations for protein
        kegg_anns = store.execute(
            """
            SELECT accession, name, description
            FROM annotations
            WHERE protein_id = ? AND source = 'kegg'
            """,
            [protein_id],
        )

        if not kegg_anns:
            return ctx.make_result(
                data=f"No KEGG annotations found for {protein_id}",
                rows=0,
            )

        lines = [
            f"# KEGG Pathway Context for {protein_id[:40]}...",
            "",
            "## KEGG Orthologs",
            "",
        ]

        # KEGG pathway mappings (common pathways)
        pathway_hints = {
            "dehydrogenase": "Energy metabolism",
            "kinase": "Signal transduction",
            "transferase": "Biosynthesis",
            "synthase": "Biosynthesis",
            "reductase": "Redox reactions",
            "oxidase": "Oxidative metabolism",
            "transporter": "Membrane transport",
            "permease": "Membrane transport",
            "peptidase": "Protein processing",
            "lyase": "Carbon-carbon bond cleavage",
        }

        pathways = []
        for accession, name, description in kegg_anns:
            lines.append(f"- **{accession}**: {name or description or 'Unknown'}")

            # Infer pathway from description
            desc_lower = (description or "").lower()
            for keyword, pathway in pathway_hints.items():
                if keyword in desc_lower:
                    pathways.append(pathway)

            # Add KEGG link
            lines.append(f"  - Link: https://www.kegg.jp/entry/{accession}")

        if pathways:
            lines.append("")
            lines.append("## Likely Pathways")
            for p in set(pathways):
                lines.append(f"- {p}")

        lines.append("")
        lines.append("## External Resources")
        lines.append(f"- [KEGG Orthology Search](https://www.kegg.jp/kegg/ko.html)")
        lines.append(f"- [KEGG Pathway Maps](https://www.kegg.jp/kegg/pathway.html)")

        return ctx.make_result(
            data="\n".join(lines),
            rows=len(kegg_anns),
            raw={"kegg_orthologs": [{"accession": a, "name": n, "description": d} for a, n, d in kegg_anns]},
        )


__all__ = [
    "visualize_neighborhood",
    "visualize_domain_architecture",
    "get_kegg_pathway_context",
]
