"""
Navigation operators for browsing genomes, proteins, and neighborhoods.

These operators enable exploration of the dataset through listing,
filtering, and contextual views.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Optional

from bennu.operators.base import BennuResult, OperatorContext
from bennu.predicates.evaluator import compute_predicates_for_protein

if TYPE_CHECKING:
    from bennu.storage.duckdb_store import DuckDBStore


def list_genomes(
    store: "DuckDBStore",
    taxonomy_filter: Optional[str] = None,
    min_completeness: Optional[float] = None,
    max_contamination: Optional[float] = None,
    limit: int = 20,
    offset: int = 0,
) -> BennuResult:
    """
    List genomes (MAGs) with optional filtering.

    Args:
        store: DuckDB store
        taxonomy_filter: Filter by taxonomy substring (e.g., "Archaea")
        min_completeness: Minimum completeness percentage
        max_contamination: Maximum contamination percentage
        limit: Maximum results to return
        offset: Pagination offset

    Returns:
        BennuResult with formatted genome list
    """
    params = {
        "taxonomy_filter": taxonomy_filter,
        "min_completeness": min_completeness,
        "max_contamination": max_contamination,
        "limit": limit,
        "offset": offset,
    }

    with OperatorContext("list_genomes", params) as ctx:
        # Build query
        clauses = []
        query_params = []

        if taxonomy_filter:
            clauses.append("LOWER(taxonomy) LIKE LOWER(?)")
            query_params.append(f"%{taxonomy_filter}%")

        if min_completeness is not None:
            clauses.append("completeness >= ?")
            query_params.append(min_completeness)

        if max_contamination is not None:
            clauses.append("contamination <= ?")
            query_params.append(max_contamination)

        where = f"WHERE {' AND '.join(clauses)}" if clauses else ""

        # Get total count
        count_query = f"SELECT COUNT(*) FROM bins {where}"
        total_count = store.execute(count_query, query_params or None)[0][0]

        # Get results
        query = f"""
            SELECT bin_id, completeness, contamination, taxonomy, n_contigs, total_length
            FROM bins
            {where}
            ORDER BY completeness DESC NULLS LAST
            LIMIT ? OFFSET ?
        """
        query_params.extend([limit, offset])
        rows = store.execute(query, query_params)

        # Format output
        lines = [
            "# Genomes",
            f"Showing {len(rows)} of {total_count:,} genomes",
            "",
            "| Bin ID | Comp% | Cont% | Contigs | Length | Taxonomy |",
            "|--------|-------|-------|---------|--------|----------|",
        ]

        genomes = []
        for row in rows:
            bin_id, comp, cont, tax, n_contigs, total_len = row
            # Truncate taxonomy for display
            tax_short = _truncate_taxonomy(tax) if tax else "Unknown"
            comp_str = f"{comp:.1f}" if comp is not None else "-"
            cont_str = f"{cont:.1f}" if cont is not None else "-"
            len_str = _format_bp(total_len) if total_len else "-"

            lines.append(
                f"| {bin_id[:20]} | {comp_str:>5} | {cont_str:>5} | {n_contigs or 0:>7} | {len_str:>6} | {tax_short[:30]} |"
            )
            genomes.append({
                "bin_id": bin_id,
                "completeness": comp,
                "contamination": cont,
                "taxonomy": tax,
                "n_contigs": n_contigs,
                "total_length": total_len,
            })

        data = "\n".join(lines)
        truncated = len(rows) < total_count

        return ctx.make_result(
            data=data,
            rows=len(rows),
            total_rows=total_count,
            truncated=truncated,
            raw=genomes,
        )


def list_proteins(
    store: "DuckDBStore",
    genome_id: Optional[str] = None,
    contig_id: Optional[str] = None,
    min_length: Optional[int] = None,
    max_length: Optional[int] = None,
    has_annotation: Optional[bool] = None,
    limit: int = 50,
    offset: int = 0,
) -> BennuResult:
    """
    List proteins with optional filtering.

    Args:
        store: DuckDB store
        genome_id: Filter to specific genome (bin_id)
        contig_id: Filter to specific contig
        min_length: Minimum protein length (aa)
        max_length: Maximum protein length (aa)
        has_annotation: Filter by annotation status
        limit: Maximum results to return
        offset: Pagination offset

    Returns:
        BennuResult with formatted protein list
    """
    params = {
        "genome_id": genome_id,
        "contig_id": contig_id,
        "min_length": min_length,
        "max_length": max_length,
        "has_annotation": has_annotation,
        "limit": limit,
        "offset": offset,
    }

    with OperatorContext("list_proteins", params) as ctx:
        # Build query
        clauses = []
        query_params = []

        if genome_id:
            clauses.append("p.bin_id = ?")
            query_params.append(genome_id)

        if contig_id:
            clauses.append("p.contig_id = ?")
            query_params.append(contig_id)

        if min_length is not None:
            clauses.append("(p.sequence_length >= ? OR (p.end_coord - p.start) / 3 >= ?)")
            query_params.extend([min_length, min_length])

        if max_length is not None:
            clauses.append("(p.sequence_length <= ? OR (p.end_coord - p.start) / 3 <= ?)")
            query_params.extend([max_length, max_length])

        where = f"WHERE {' AND '.join(clauses)}" if clauses else ""

        # Handle annotation filter with subquery
        if has_annotation is True:
            if where:
                where += " AND p.protein_id IN (SELECT DISTINCT protein_id FROM annotations)"
            else:
                where = "WHERE p.protein_id IN (SELECT DISTINCT protein_id FROM annotations)"
        elif has_annotation is False:
            if where:
                where += " AND p.protein_id NOT IN (SELECT DISTINCT protein_id FROM annotations)"
            else:
                where = "WHERE p.protein_id NOT IN (SELECT DISTINCT protein_id FROM annotations)"

        # Get total count
        count_query = f"SELECT COUNT(*) FROM proteins p {where}"
        total_count = store.execute(count_query, query_params or None)[0][0]

        # Get results with best annotation
        query = f"""
            SELECT
                p.protein_id,
                p.contig_id,
                p.bin_id,
                p.start,
                p.end_coord,
                p.strand,
                COALESCE(p.sequence_length, (p.end_coord - p.start) / 3) as length_aa,
                (
                    SELECT a.name || ' (' || a.accession || ')'
                    FROM annotations a
                    WHERE a.protein_id = p.protein_id
                    ORDER BY a.evalue NULLS LAST
                    LIMIT 1
                ) as best_annotation
            FROM proteins p
            {where}
            ORDER BY p.contig_id, p.start
            LIMIT ? OFFSET ?
        """
        query_params.extend([limit, offset])
        rows = store.execute(query, query_params)

        # Format output
        lines = [
            "# Proteins",
            f"Showing {len(rows)} of {total_count:,} proteins",
            "",
            "| Protein ID | Length | Strand | Annotation |",
            "|------------|--------|--------|------------|",
        ]

        proteins = []
        for row in rows:
            protein_id, contig_id, bin_id, start, end, strand, length_aa, annotation = row
            ann_str = (annotation or "NO HITS")[:40]
            id_short = protein_id[:25] if len(protein_id) > 25 else protein_id

            lines.append(
                f"| {id_short} | {length_aa:>6}aa | {strand:>6} | {ann_str} |"
            )
            proteins.append({
                "protein_id": protein_id,
                "contig_id": contig_id,
                "bin_id": bin_id,
                "start": start,
                "end": end,
                "strand": strand,
                "length_aa": length_aa,
                "annotation": annotation,
            })

        data = "\n".join(lines)
        truncated = len(rows) < total_count

        return ctx.make_result(
            data=data,
            rows=len(rows),
            total_rows=total_count,
            truncated=truncated,
            raw=proteins,
        )


def get_genome(store: "DuckDBStore", genome_id: str, verbosity: int = 1) -> BennuResult:
    """
    Get detailed information about a specific genome.

    Args:
        store: DuckDB store
        genome_id: Bin ID to retrieve
        verbosity: 0=minimal, 1=standard, 2=detailed

    Returns:
        BennuResult with genome details
    """
    params = {"genome_id": genome_id, "verbosity": verbosity}

    with OperatorContext("get_genome", params) as ctx:
        # Get genome info
        row = store.execute(
            """
            SELECT bin_id, completeness, contamination, taxonomy, n_contigs, total_length
            FROM bins
            WHERE bin_id = ?
            """,
            [genome_id],
        )

        if not row:
            return ctx.make_result(
                data=f"Genome not found: {genome_id}",
                rows=0,
                total_rows=0,
            )

        bin_id, comp, cont, tax, n_contigs, total_len = row[0]

        # Get protein count
        protein_count = store.execute(
            "SELECT COUNT(*) FROM proteins WHERE bin_id = ?",
            [genome_id],
        )[0][0]

        # Get annotation stats
        annotation_stats = store.execute(
            """
            SELECT source, COUNT(*) as count
            FROM annotations
            WHERE protein_id IN (SELECT protein_id FROM proteins WHERE bin_id = ?)
            GROUP BY source
            ORDER BY count DESC
            """,
            [genome_id],
        )

        lines = [
            f"# Genome: {bin_id}",
            "",
            "## Quality",
            f"- Completeness: {comp:.1f}%" if comp is not None else "- Completeness: Unknown",
            f"- Contamination: {cont:.1f}%" if cont is not None else "- Contamination: Unknown",
            "",
            "## Size",
            f"- Contigs: {n_contigs:,}" if n_contigs else "- Contigs: Unknown",
            f"- Total length: {_format_bp(total_len)}" if total_len else "- Total length: Unknown",
            f"- Proteins: {protein_count:,}",
            "",
            "## Taxonomy",
            f"{tax or 'Unknown'}",
            "",
        ]

        if verbosity >= 1 and annotation_stats:
            lines.extend([
                "## Annotations by Source",
            ])
            for source, count in annotation_stats:
                lines.append(f"- {source}: {count:,}")
            lines.append("")

        data = "\n".join(lines)

        return ctx.make_result(
            data=data,
            rows=1,
            total_rows=1,
            raw={
                "bin_id": bin_id,
                "completeness": comp,
                "contamination": cont,
                "taxonomy": tax,
                "n_contigs": n_contigs,
                "total_length": total_len,
                "protein_count": protein_count,
            },
        )


def get_protein(store: "DuckDBStore", protein_id: str, verbosity: int = 1) -> BennuResult:
    """
    Get detailed information about a specific protein.

    Args:
        store: DuckDB store
        protein_id: Protein ID to retrieve
        verbosity: 0=minimal, 1=standard, 2=detailed

    Returns:
        BennuResult with protein details
    """
    params = {"protein_id": protein_id, "verbosity": verbosity}

    with OperatorContext("get_protein", params) as ctx:
        # Get protein info
        row = store.execute(
            """
            SELECT protein_id, contig_id, bin_id, start, end_coord, strand,
                   sequence, sequence_length, gc_content
            FROM proteins
            WHERE protein_id = ?
            """,
            [protein_id],
        )

        if not row:
            return ctx.make_result(
                data=f"Protein not found: {protein_id}",
                rows=0,
                total_rows=0,
            )

        pid, contig_id, bin_id, start, end, strand, seq, seq_len, gc = row[0]
        length_aa = seq_len or ((end - start) // 3)

        # Get annotations
        annotations = store.execute(
            """
            SELECT source, accession, name, description, evalue, score
            FROM annotations
            WHERE protein_id = ?
            ORDER BY evalue NULLS LAST
            """,
            [protein_id],
        )

        # Compute predicates
        predicates = compute_predicates_for_protein(protein_id, store)

        lines = [
            f"# Protein: {pid}",
            "",
            "## Location",
            f"- Contig: {contig_id}",
            f"- Genome: {bin_id or 'Unknown'}",
            f"- Coordinates: {start:,}-{end:,} ({strand})",
            f"- Length: {length_aa:,} aa",
        ]

        if gc is not None:
            lines.append(f"- GC content: {gc:.1%}")

        lines.append("")

        if predicates:
            lines.extend([
                "## Predicates",
                ", ".join(predicates),
                "",
            ])

        if annotations:
            lines.extend([
                "## Annotations",
                "| Source | Accession | Name | E-value |",
                "|--------|-----------|------|---------|",
            ])
            for ann in annotations:
                source, acc, name, desc, evalue, score = ann
                ev_str = f"{evalue:.1e}" if evalue is not None else "-"
                name_str = (name or desc or "-")[:30]
                lines.append(f"| {source} | {acc} | {name_str} | {ev_str} |")
            lines.append("")
        else:
            lines.extend([
                "## Annotations",
                "No annotations found.",
                "",
            ])

        if verbosity >= 2 and seq:
            lines.extend([
                "## Sequence",
                f"```",
                _format_sequence(seq),
                "```",
                "",
            ])

        data = "\n".join(lines)

        return ctx.make_result(
            data=data,
            rows=1,
            total_rows=1,
            raw={
                "protein_id": pid,
                "contig_id": contig_id,
                "bin_id": bin_id,
                "start": start,
                "end": end,
                "strand": strand,
                "length_aa": length_aa,
                "gc_content": gc,
                "predicates": predicates,
                "annotations": [
                    {"source": a[0], "accession": a[1], "name": a[2], "description": a[3], "evalue": a[4], "score": a[5]}
                    for a in annotations
                ],
            },
        )


def get_neighborhood(
    store: "DuckDBStore",
    entity_id: str,
    window: int = 10,
    verbosity: int = 1,
    all_annotations: bool = False,
) -> BennuResult:
    """
    Get genomic neighborhood around a protein.

    Args:
        store: DuckDB store
        entity_id: Protein ID as anchor
        window: Number of genes on each side
        verbosity: 0=minimal, 1=standard, 2=detailed
        all_annotations: If True, return all annotation sources per gene
            (PFAM, KEGG, VOGdb, DefenseFinder, HydDB, CAZy) instead of
            just the best hit. Essential for context-based functional
            interpretation — the domain tells you the fold, the full
            annotation context tells you the function.

    Returns:
        BennuResult with neighborhood as ASCII table
    """
    params = {
        "entity_id": entity_id,
        "window": window,
        "verbosity": verbosity,
        "all_annotations": all_annotations,
    }

    with OperatorContext("get_neighborhood", params) as ctx:
        # Get anchor protein
        anchor = store.execute(
            """
            SELECT protein_id, contig_id, bin_id, start, end_coord, strand, gene_index
            FROM proteins
            WHERE protein_id = ?
            """,
            [entity_id],
        )

        if not anchor:
            return ctx.make_result(
                data=f"Protein not found: {entity_id}",
                rows=0,
                total_rows=0,
            )

        anchor_pid, contig_id, bin_id, anchor_start, anchor_end, anchor_strand, anchor_idx = anchor[0]

        # Get neighboring proteins by position
        neighbors = store.execute(
            """
            SELECT
                p.protein_id,
                p.start,
                p.end_coord,
                p.strand,
                COALESCE(p.sequence_length, (p.end_coord - p.start) / 3) as length_aa,
                p.gene_index,
                (
                    SELECT a.name || COALESCE(' (' || a.accession || ')', '')
                    FROM annotations a
                    WHERE a.protein_id = p.protein_id
                    ORDER BY a.evalue NULLS LAST
                    LIMIT 1
                ) as best_annotation
            FROM proteins p
            WHERE p.contig_id = ?
            ORDER BY p.start
            """,
            [contig_id],
        )

        if not neighbors:
            return ctx.make_result(
                data="No proteins found on contig",
                rows=0,
                total_rows=0,
            )

        # Find anchor index and slice window
        anchor_list_idx = None
        for i, row in enumerate(neighbors):
            if row[0] == entity_id:
                anchor_list_idx = i
                break

        if anchor_list_idx is None:
            return ctx.make_result(
                data="Anchor protein not found in contig",
                rows=0,
                total_rows=0,
            )

        start_idx = max(0, anchor_list_idx - window)
        end_idx = min(len(neighbors), anchor_list_idx + window + 1)
        window_proteins = neighbors[start_idx:end_idx]

        # Calculate region bounds
        region_start = window_proteins[0][1]
        region_end = window_proteins[-1][2]

        # If all_annotations requested, fetch full annotation sets
        annotations_by_protein: dict[str, list[dict[str, Any]]] = {}
        if all_annotations:
            protein_ids = [p[0] for p in window_proteins]
            # Batch fetch all annotations for proteins in the window
            placeholders = ", ".join(["?"] * len(protein_ids))
            all_annots = store.execute(
                f"""
                SELECT
                    a.protein_id,
                    a.source,
                    a.accession,
                    a.name,
                    a.description,
                    a.evalue,
                    a.score
                FROM annotations a
                WHERE a.protein_id IN ({placeholders})
                ORDER BY a.protein_id, a.source, a.evalue NULLS LAST
                """,
                protein_ids,
            )
            for row in all_annots:
                pid = row[0]
                if pid not in annotations_by_protein:
                    annotations_by_protein[pid] = []
                annotations_by_protein[pid].append({
                    "source": row[1],
                    "accession": row[2],
                    "name": row[3],
                    "description": row[4],
                    "evalue": row[5],
                    "score": row[6],
                })

        # Format header
        lines = [
            f"# Neighborhood: {entity_id}",
            f"**Contig:** {contig_id}",
            f"**Region:** {region_start:,}-{region_end:,} bp | {len(window_proteins)} genes",
            "",
        ]

        # Format as ASCII table
        if all_annotations:
            lines.extend(_format_neighborhood_table_full(
                window_proteins,
                entity_id,
                start_idx,
                annotations_by_protein,
            ))
        else:
            lines.extend(_format_neighborhood_table(
                window_proteins,
                entity_id,
                start_idx,
                store,
                verbosity,
            ))

        data = "\n".join(lines)

        # Build protein dicts for raw output
        protein_dicts = []
        for p in window_proteins:
            d: dict[str, Any] = {
                "protein_id": p[0],
                "start": p[1],
                "end": p[2],
                "strand": p[3],
                "length_aa": p[4],
                "gene_index": p[5],
                "annotation": p[6],
                "is_anchor": p[0] == entity_id,
            }
            if all_annotations:
                d["annotations"] = _group_annotations_by_source(
                    annotations_by_protein.get(p[0], [])
                )
            protein_dicts.append(d)

        return ctx.make_result(
            data=data,
            rows=len(window_proteins),
            total_rows=len(neighbors),
            raw={
                "anchor_protein_id": entity_id,
                "contig_id": contig_id,
                "bin_id": bin_id,
                "region_start": region_start,
                "region_end": region_end,
                "proteins": protein_dicts,
            },
        )


def _group_annotations_by_source(
    annotations: list[dict[str, Any]],
) -> dict[str, list[dict[str, Any]]]:
    """Group a flat list of annotations by source (pfam, kegg, etc.)."""
    by_source: dict[str, list[dict[str, Any]]] = {}
    for ann in annotations:
        source = ann["source"]
        if source not in by_source:
            by_source[source] = []
        by_source[source].append({
            "accession": ann["accession"],
            "name": ann["name"],
            "description": ann.get("description"),
            "evalue": ann.get("evalue"),
        })
    return by_source


def _format_neighborhood_table_full(
    proteins: list,
    anchor_id: str,
    start_offset: int,
    annotations_by_protein: dict[str, list[dict[str, Any]]],
) -> list[str]:
    """Format neighborhood with all annotation sources per gene.

    Shows each gene with its annotations grouped by source, giving
    full context for functional interpretation. This is the view that
    answers "what IS this gene" rather than just "what's the best hit."
    """
    lines: list[str] = []

    for i, row in enumerate(proteins):
        protein_id, start, end, strand, length_aa, gene_idx, best_annotation = row
        is_anchor = protein_id == anchor_id

        marker = ">>>" if is_anchor else f"g{gene_idx}"
        lines.append(
            f"{'=' * 70}" if is_anchor else f"{'─' * 70}"
        )
        lines.append(
            f" {marker}  {protein_id}  |  {length_aa}aa  {strand}  "
            f"{start:,}-{end:,} bp"
        )

        annots = annotations_by_protein.get(protein_id, [])
        if not annots:
            lines.append("     NO ANNOTATIONS")
        else:
            # Group by source and display
            by_source: dict[str, list] = {}
            for ann in annots:
                src = ann["source"]
                if src not in by_source:
                    by_source[src] = []
                by_source[src].append(ann)

            # Display order: pfam, kegg, hyddb, defensefinder, vogdb, cazy, other
            source_order = [
                "pfam", "kegg", "hyddb", "defensefinder",
                "vogdb", "cazy",
            ]
            seen_sources: set[str] = set()
            for src in source_order:
                if src in by_source:
                    seen_sources.add(src)
                    _format_source_annotations(lines, src, by_source[src])
            # Any remaining sources not in the ordered list
            for src in sorted(by_source.keys()):
                if src not in seen_sources:
                    _format_source_annotations(lines, src, by_source[src])

    lines.append("=" * 70 if proteins and proteins[-1][0] == anchor_id else "─" * 70)
    return lines


def _format_source_annotations(
    lines: list[str],
    source: str,
    annotations: list[dict[str, Any]],
) -> None:
    """Format annotations from a single source for the full neighborhood view."""
    source_label = source.upper()
    for j, ann in enumerate(annotations):
        acc = ann.get("accession", "?")
        name = ann.get("name", "?")
        evalue = ann.get("evalue")
        ev_str = f"  e={evalue:.1e}" if evalue is not None else ""
        prefix = f"     {source_label}:" if j == 0 else f"     {'':>{len(source_label)}} "
        lines.append(f"{prefix} {acc} {name}{ev_str}")


def _format_neighborhood_table(
    proteins: list,
    anchor_id: str,
    start_offset: int,
    store: "DuckDBStore",
    verbosity: int,
) -> list[str]:
    """Format neighborhood as ASCII table."""
    lines = [
        "```",
        " #   Start      Len   Str  Annotation              Predicates",
        "─" * 70,
    ]

    for i, row in enumerate(proteins):
        protein_id, start, end, strand, length_aa, gene_idx, annotation = row
        is_anchor = protein_id == anchor_id

        # Get predicates if verbose
        if verbosity >= 1:
            try:
                predicates = compute_predicates_for_protein(protein_id, store)
                pred_str = ", ".join(predicates[:2]) if predicates else ""
            except Exception:
                pred_str = ""
        else:
            pred_str = ""

        # Format row
        marker = "→" if is_anchor else " "
        idx_str = f"{start_offset + i + 1:>2}"
        ann_str = (annotation or "NO HITS")[:22].ljust(22)
        pred_str = pred_str[:15]

        line = f"{marker}{idx_str}  {start:>8}  {length_aa:>4}aa  {strand:>2}  {ann_str}  {pred_str}"
        lines.append(line)

    lines.append("```")
    return lines


def _truncate_taxonomy(taxonomy: str, max_len: int = 30) -> str:
    """Truncate taxonomy string, keeping most specific level."""
    if not taxonomy:
        return "Unknown"
    # GTDB format: d__Bacteria;p__Proteobacteria;c__...
    parts = taxonomy.split(";")
    # Return last non-empty level
    for part in reversed(parts):
        if part and not part.endswith("__"):
            # Remove prefix (e.g., "s__")
            clean = part.split("__")[-1] if "__" in part else part
            if len(clean) > max_len:
                return clean[:max_len-3] + "..."
            return clean
    return "Unknown"


def _format_bp(bp: int) -> str:
    """Format base pairs with appropriate unit."""
    if bp >= 1_000_000:
        return f"{bp / 1_000_000:.1f}Mb"
    elif bp >= 1_000:
        return f"{bp / 1_000:.1f}kb"
    return f"{bp}bp"


def _format_sequence(seq: str, width: int = 60) -> str:
    """Format sequence with line breaks."""
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))


__all__ = [
    "list_genomes",
    "list_proteins",
    "get_genome",
    "get_protein",
    "get_neighborhood",
]
