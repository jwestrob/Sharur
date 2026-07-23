"""
Navigation operators for browsing genomes, proteins, and neighborhoods.

These operators enable exploration of the dataset through listing,
filtering, and contextual views.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from sharur.operators.base import OperatorContext, SharurResult
from sharur.operators.cases import structured_projection_sources
from sharur.operators.predicates_v2 import explain
from sharur.operators.semantics import (
    get_active_predicates,
    get_active_predicates_for_protein,
)


if TYPE_CHECKING:
    from sharur.storage.duckdb_store import DuckDBStore


def _load_annotation_rows(
    store: DuckDBStore,
    protein_ids: list[str],
    *,
    excluded_sources: list[str] | None = None,
) -> dict[str, list[dict[str, Any]]]:
    """Batch observed annotation rows for a bounded protein cohort."""
    if not protein_ids:
        return {}
    placeholders = ", ".join(["?"] * len(protein_ids))
    excluded = excluded_sources or []
    exclusion_sql = ""
    params: list[Any] = [*protein_ids]
    if excluded:
        source_placeholders = ", ".join(["?"] * len(excluded))
        exclusion_sql = f"AND a.source NOT IN ({source_placeholders})"
        params.extend(excluded)
    rows = store.execute(
        f"""
        SELECT
            a.protein_id,
            a.source,
            a.accession,
            a.name,
            a.description,
            a.evalue,
            a.score
        FROM annotations AS a
        WHERE a.protein_id IN ({placeholders})
          {exclusion_sql}
        ORDER BY
            a.protein_id,
            a.evalue NULLS LAST,
            a.score DESC NULLS LAST,
            a.annotation_id
        """,
        params,
    )
    by_protein: dict[str, list[dict[str, Any]]] = {}
    for protein_id, source, accession, name, description, evalue, score in rows:
        by_protein.setdefault(str(protein_id), []).append(
            {
                "source": source,
                "accession": accession,
                "name": name,
                "description": description,
                "evalue": evalue,
                "score": score,
            }
        )
    return by_protein


def _annotation_label(annotation: dict[str, Any] | None) -> str | None:
    if annotation is None:
        return None
    name = annotation.get("name") or annotation.get("description")
    accession = annotation.get("accession")
    if name and accession:
        return f"{name} ({accession})"
    return str(name or accession) if name or accession else None


def list_genomes(
    store: DuckDBStore,
    taxonomy_filter: str | None = None,
    min_completeness: float | None = None,
    max_contamination: float | None = None,
    limit: int = 20,
    offset: int = 0,
) -> SharurResult:
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
        SharurResult with formatted genome list
    """
    params = {
        "taxonomy_filter": taxonomy_filter,
        "min_completeness": min_completeness,
        "max_contamination": max_contamination,
        "limit": limit,
        "offset": offset,
    }

    with OperatorContext("list_genomes", params, store=store) as ctx:
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
            ORDER BY completeness DESC NULLS LAST, bin_id
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
    store: DuckDBStore,
    genome_id: str | None = None,
    contig_id: str | None = None,
    min_length: int | None = None,
    max_length: int | None = None,
    has_annotation: bool | None = None,
    limit: int = 50,
    offset: int = 0,
) -> SharurResult:
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
        SharurResult with formatted protein list
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

    with OperatorContext("list_proteins", params, store=store) as ctx:
        clauses: list[str] = []
        filter_params: list[Any] = []

        if genome_id:
            clauses.append("p.bin_id = ?")
            filter_params.append(genome_id)

        if contig_id:
            clauses.append("p.contig_id = ?")
            filter_params.append(contig_id)

        if min_length is not None:
            clauses.append("(p.sequence_length >= ? OR (p.end_coord - p.start) / 3 >= ?)")
            filter_params.extend([min_length, min_length])

        if max_length is not None:
            clauses.append("(p.sequence_length <= ? OR (p.end_coord - p.start) / 3 <= ?)")
            filter_params.extend([max_length, max_length])

        where = f"WHERE {' AND '.join(clauses)}" if clauses else ""
        excluded_sources = structured_projection_sources(store)
        annotation_relation = "SELECT DISTINCT protein_id FROM annotations"
        annotation_params: list[Any] = []
        if excluded_sources:
            placeholders = ", ".join(["?"] * len(excluded_sources))
            annotation_relation += f" WHERE source NOT IN ({placeholders})"
            annotation_params.extend(excluded_sources)
        annotation_join = ""
        if has_annotation is True:
            annotation_join = (
                f"SEMI JOIN ({annotation_relation}) AS annotated "
                "ON annotated.protein_id = p.protein_id"
            )
        elif has_annotation is False:
            annotation_join = (
                f"ANTI JOIN ({annotation_relation}) AS annotated "
                "ON annotated.protein_id = p.protein_id"
            )

        query_params = [
            *(annotation_params if annotation_join else []),
            *filter_params,
        ]
        total_count = store.execute(
            f"""
            SELECT COUNT(*)
            FROM proteins AS p
            {annotation_join}
            {where}
            """,
            query_params or None,
        )[0][0]

        query = f"""
            SELECT
                p.protein_id,
                p.contig_id,
                p.bin_id,
                p.start,
                p.end_coord,
                p.strand,
                CAST(
                    COALESCE(p.sequence_length, (p.end_coord - p.start) / 3)
                    AS BIGINT
                ) AS length_aa
            FROM proteins AS p
            {annotation_join}
            {where}
            ORDER BY p.contig_id, p.start, p.protein_id
            LIMIT ? OFFSET ?
        """
        rows = store.execute(query, [*query_params, limit, offset])
        annotations = _load_annotation_rows(
            store,
            [str(row[0]) for row in rows],
            excluded_sources=excluded_sources,
        )

        lines = [
            "# Proteins",
            f"Showing {len(rows)} of {total_count:,} proteins",
            "",
            "| Protein ID | Length | Strand | Annotation |",
            "|------------|--------|--------|------------|",
        ]

        proteins = []
        for row in rows:
            protein_id, contig_id, bin_id, start, end, strand, length_aa = row
            protein_annotations = annotations.get(str(protein_id), [])
            annotation = _annotation_label(
                protein_annotations[0] if protein_annotations else None
            )
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


def get_genome(store: DuckDBStore, genome_id: str, verbosity: int = 1) -> SharurResult:
    """
    Get detailed information about a specific genome.

    Args:
        store: DuckDB store
        genome_id: Bin ID to retrieve
        verbosity: 0=minimal, 1=standard, 2=detailed

    Returns:
        SharurResult with genome details
    """
    params = {"genome_id": genome_id, "verbosity": verbosity}

    with OperatorContext("get_genome", params, store=store) as ctx:
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

        excluded_sources = structured_projection_sources(store)
        exclusion_sql = ""
        annotation_params: list[Any] = [genome_id]
        if excluded_sources:
            placeholders = ", ".join(["?"] * len(excluded_sources))
            exclusion_sql = f"AND a.source NOT IN ({placeholders})"
            annotation_params.extend(excluded_sources)
        annotation_stats = store.execute(
            f"""
            SELECT
                a.source,
                COUNT(*) AS domain_hits,
                COUNT(DISTINCT a.protein_id) AS proteins
            FROM proteins AS p
            JOIN annotations AS a ON a.protein_id = p.protein_id
            WHERE p.bin_id = ?
              {exclusion_sql}
            GROUP BY a.source
            ORDER BY proteins DESC, domain_hits DESC, a.source
            """,
            annotation_params,
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
                "| Source | Proteins | Domain hits |",
                "|--------|---------:|------------:|",
            ])
            for source, domain_hits, annotated_proteins in annotation_stats:
                lines.append(
                    f"| {source} | {annotated_proteins:,} | {domain_hits:,} |"
                )
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
                "annotation_stats": [
                    {
                        "source": source,
                        "domain_hits": domain_hits,
                        "proteins": annotated_proteins,
                    }
                    for source, domain_hits, annotated_proteins in annotation_stats
                ],
            },
        )


def get_protein(store: DuckDBStore, protein_id: str, verbosity: int = 1) -> SharurResult:
    """
    Get detailed information about a specific protein.

    Args:
        store: DuckDB store
        protein_id: Protein ID to retrieve
        verbosity: 0=minimal, 1=standard, 2=detailed

    Returns:
        SharurResult with protein details
    """
    params = {"protein_id": protein_id, "verbosity": verbosity}

    with OperatorContext("get_protein", params, store=store) as ctx:
        # Get protein info
        row = store.execute(
            """
            SELECT protein_id, contig_id, bin_id, start, end_coord, strand,
                   sequence_length, gc_content
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

        pid, contig_id, bin_id, start, end, strand, seq_len, gc = row[0]
        length_aa = seq_len or ((end - start) // 3)

        annotations = _load_annotation_rows(
            store,
            [protein_id],
            excluded_sources=structured_projection_sources(store),
        ).get(protein_id, [])

        # Read the same persisted V2-compatible predicate view used by search.
        predicates = get_active_predicates_for_protein(store, protein_id)

        try:
            semantic_explanation = explain(store, protein_id)
        except Exception as exc:
            semantic_explanation = {
                "status": "unavailable",
                "error": f"{type(exc).__name__}: {exc}",
            }

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
                "## Active Semantic Predicates",
                ", ".join(predicates),
                "",
            ])

        validated_systems = semantic_explanation.get("validated_systems", [])
        if validated_systems:
            lines.append("## Validated System Membership")
            for system in validated_systems:
                details = [
                    str(value)
                    for value in (
                        system.get("profile_name"),
                        f"position {system.get('position')}"
                        if system.get("position") is not None
                        else None,
                    )
                    if value
                ]
                suffix = f" ({', '.join(details)})" if details else ""
                lines.append(
                    f"- {system.get('system_source')}: "
                    f"{system.get('system_id')}{suffix}"
                )
            lines.append("")

        if semantic_explanation.get("status") == "unavailable":
            lines.extend([
                "## Semantic Status",
                f"Unavailable: {semantic_explanation['error']}",
                "",
            ])

        if annotations:
            lines.extend([
                "## Annotations",
                "| Source | Accession | Name | E-value |",
                "|--------|-----------|------|---------|",
            ])
            for ann in annotations:
                source = ann["source"]
                acc = ann["accession"]
                name = ann["name"]
                desc = ann["description"]
                evalue = ann["evalue"]
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

        if verbosity >= 2:
            lines.extend([
                "## Sequence Access",
                "Raw sequence is compute-only. Use `get_sequence()` or an export "
                "operator with a local output path.",
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
                "semantic": semantic_explanation,
                "annotations": annotations,
            },
            warnings=(
                ["Raw sequence is intentionally omitted from model-visible output."]
                if verbosity >= 2
                else None
            ),
        )


def get_neighborhood(
    store: DuckDBStore,
    entity_id: str,
    window: int = 10,
    verbosity: int = 1,
    all_annotations: bool = False,
) -> SharurResult:
    """
    Get genomic neighborhood around a protein.

    Args:
        store: DuckDB store
        entity_id: Protein ID as anchor
        window: Number of genes on each side
        verbosity: 0=minimal, 1=standard, 2=detailed
        all_annotations: If True, return all live observed annotation sources
            per gene instead of only the best hit. Structured caller output is
            kept separate from these per-domain observations.

    Returns:
        SharurResult with neighborhood as ASCII table
    """
    params = {
        "entity_id": entity_id,
        "window": window,
        "verbosity": verbosity,
        "all_annotations": all_annotations,
    }

    with OperatorContext("get_neighborhood", params, store=store) as ctx:
        if window < 0:
            raise ValueError("window must be non-negative")
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

        _anchor_pid, contig_id, bin_id, *_anchor_fields = anchor[0]
        neighbor_rows = store.execute(
            """
            WITH ranked AS (
                SELECT
                    p.protein_id,
                    p.start,
                    p.end_coord,
                    p.strand,
                    CAST(
                        COALESCE(
                            p.sequence_length,
                            (p.end_coord - p.start) / 3
                        )
                        AS BIGINT
                    ) AS length_aa,
                    p.gene_index,
                    ROW_NUMBER() OVER (
                        ORDER BY p.start, p.protein_id
                    ) AS row_number,
                    COUNT(*) OVER () AS total_rows
                FROM proteins AS p
                WHERE p.contig_id = ?
                  AND p.bin_id IS NOT DISTINCT FROM ?
            ),
            anchor_position AS (
                SELECT row_number
                FROM ranked
                WHERE protein_id = ?
            )
            SELECT
                ranked.protein_id,
                ranked.start,
                ranked.end_coord,
                ranked.strand,
                ranked.length_aa,
                ranked.gene_index,
                ranked.row_number,
                ranked.total_rows
            FROM ranked
            CROSS JOIN anchor_position
            WHERE ranked.row_number BETWEEN
                  anchor_position.row_number - ?
              AND anchor_position.row_number + ?
            ORDER BY ranked.row_number
            """,
            [contig_id, bin_id, entity_id, window, window],
        )

        if not neighbor_rows:
            return ctx.make_result(
                data="No proteins found on contig",
                rows=0,
                total_rows=0,
            )

        excluded_sources = structured_projection_sources(store)
        annotations_by_protein = _load_annotation_rows(
            store,
            [str(row[0]) for row in neighbor_rows],
            excluded_sources=excluded_sources,
        )
        window_proteins = []
        for row in neighbor_rows:
            protein_annotations = annotations_by_protein.get(str(row[0]), [])
            best_annotation = _annotation_label(
                protein_annotations[0] if protein_annotations else None
            )
            window_proteins.append((*row[:6], best_annotation))
        start_idx = int(neighbor_rows[0][6]) - 1
        total_proteins = int(neighbor_rows[0][7])
        predicates_by_protein = (
            get_active_predicates(store, [row[0] for row in window_proteins])
            if verbosity >= 1
            else {}
        )

        # Calculate region bounds
        region_start = window_proteins[0][1]
        region_end = window_proteins[-1][2]

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
                verbosity,
                predicates_by_protein,
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
                "predicates": predicates_by_protein.get(p[0], []),
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
            total_rows=total_proteins,
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

    for row in proteins:
        protein_id, start, end, strand, length_aa, gene_idx, _best_annotation = row
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
    verbosity: int,
    predicates_by_protein: dict[str, list[str]],
) -> list[str]:
    """Format neighborhood as ASCII table."""
    lines = [
        "```",
        " #   Start      Len   Str  Annotation              Predicates",
        "─" * 70,
    ]

    for i, row in enumerate(proteins):
        protein_id, start, _end, strand, length_aa, _gene_idx, annotation = row
        is_anchor = protein_id == anchor_id

        if verbosity >= 1:
            predicates = predicates_by_protein.get(protein_id, [])
            pred_str = ", ".join(predicates[:2]) if predicates else ""
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
    if bp >= 1_000:
        return f"{bp / 1_000:.1f}kb"
    return f"{bp}bp"


__all__ = [
    "get_genome",
    "get_neighborhood",
    "get_protein",
    "list_genomes",
    "list_proteins",
]
