"""
Search operators for finding proteins by various criteria.

Includes predicate-based search using DuckDB array functions.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Optional

from sharur.operators.base import SharurResult, OperatorContext
from sharur.predicates.registry import get_registry
from sharur.predicates.evaluator import evaluate_predicate
from sharur.predicates.vocabulary import PREDICATE_BY_ID

if TYPE_CHECKING:
    from sharur.storage.duckdb_store import DuckDBStore


def search_by_predicates(
    store: "DuckDBStore",
    has: Optional[list[str]] = None,
    lacks: Optional[list[str]] = None,
    limit: int = 50,
    offset: int = 0,
) -> SharurResult:
    """
    Search proteins by predicate membership.

    Uses set intersection/difference logic:
    - has: Protein must have ALL of these predicates (AND)
    - lacks: Protein must have NONE of these predicates (AND NOT)

    Args:
        store: DuckDB store
        has: List of predicates that must be true
        lacks: List of predicates that must be false
        limit: Maximum results
        offset: Pagination offset

    Returns:
        SharurResult with matching proteins
    """
    has = has or []
    lacks = lacks or []

    params = {"has": has, "lacks": lacks, "limit": limit, "offset": offset}

    with OperatorContext("search_by_predicates", params) as ctx:
        registry = get_registry()

        # Validate predicates - accept vocabulary, direct-access, or legacy registry
        def is_valid_predicate(pred_id: str) -> bool:
            # Check vocabulary
            if pred_id in PREDICATE_BY_ID:
                return True
            # Check direct-access predicates (pfam:*, kegg:*, cazy:*)
            if ":" in pred_id:
                return True
            # Check legacy registry
            if registry.exists(pred_id):
                return True
            return False

        invalid_preds = [p for p in has + lacks if not is_valid_predicate(p)]
        if invalid_preds:
            # List some example predicates from vocabulary
            example_preds = list(PREDICATE_BY_ID.keys())[:20]
            return ctx.make_result(
                data=f"Unknown predicate(s): {', '.join(invalid_preds)}\n\nExample predicates: {', '.join(example_preds)}...\n\nUse direct access: pfam:PF00001, kegg:K00001, cazy:GH1",
                rows=0,
                total_rows=0,
            )

        # Check if protein_predicates table has data
        has_precomputed = False
        try:
            count = store.execute("SELECT COUNT(*) FROM protein_predicates")[0][0]
            has_precomputed = count > 0
        except Exception:
            pass

        if has_precomputed:
            # Use pre-computed predicates with DuckDB array functions
            return _search_precomputed(store, has, lacks, limit, offset, ctx)
        else:
            # Evaluate predicates on-the-fly
            return _search_dynamic(store, has, lacks, limit, offset, ctx)


def _search_precomputed(
    store: "DuckDBStore",
    has: list[str],
    lacks: list[str],
    limit: int,
    offset: int,
    ctx: OperatorContext,
) -> SharurResult:
    """Search using pre-computed protein_predicates table."""

    # Build query using DuckDB array functions
    # list_has_all checks all elements are present
    # list_contains checks single element presence
    clauses = []
    params = []

    if has:
        # All predicates in 'has' must be present
        clauses.append("list_has_all(pp.predicates, ?)")
        params.append(has)

    if lacks:
        # None of the predicates in 'lacks' should be present
        for pred in lacks:
            clauses.append("NOT list_contains(pp.predicates, ?)")
            params.append(pred)

    where = f"WHERE {' AND '.join(clauses)}" if clauses else ""

    # Get total count
    count_query = f"""
        SELECT COUNT(*)
        FROM protein_predicates pp
        {where}
    """
    total_count = store.execute(count_query, params or None)[0][0]

    # Get results with protein details
    query = f"""
        SELECT
            p.protein_id,
            p.contig_id,
            p.bin_id,
            COALESCE(p.sequence_length, (p.end_coord - p.start) / 3) as length_aa,
            p.strand,
            pp.predicates,
            (
                SELECT a.name || ' (' || a.accession || ')'
                FROM annotations a
                WHERE a.protein_id = p.protein_id
                ORDER BY a.evalue NULLS LAST
                LIMIT 1
            ) as best_annotation
        FROM protein_predicates pp
        JOIN proteins p ON pp.protein_id = p.protein_id
        {where}
        ORDER BY length_aa DESC
        LIMIT ? OFFSET ?
    """
    params.extend([limit, offset])
    rows = store.execute(query, params)

    return _format_search_results(rows, total_count, has, lacks, ctx)


def _search_dynamic(
    store: "DuckDBStore",
    has: list[str],
    lacks: list[str],
    limit: int,
    offset: int,
    ctx: OperatorContext,
) -> SharurResult:
    """Search by evaluating predicates dynamically."""

    # Start with all proteins
    all_proteins = set(row[0] for row in store.execute("SELECT protein_id FROM proteins"))

    # Apply 'has' predicates (intersection)
    matching = all_proteins
    for pred_id in has:
        pred_matches = evaluate_predicate(pred_id, store)
        matching = matching.intersection(pred_matches)

    # Apply 'lacks' predicates (difference)
    for pred_id in lacks:
        pred_matches = evaluate_predicate(pred_id, store)
        matching = matching.difference(pred_matches)

    total_count = len(matching)
    matching_list = sorted(matching)[offset : offset + limit]

    if not matching_list:
        return ctx.make_result(
            data=f"No proteins found matching predicates.\nHas: {has}\nLacks: {lacks}",
            rows=0,
            total_rows=0,
            raw=[],
        )

    # Get protein details
    placeholders = ",".join(["?"] * len(matching_list))
    query = f"""
        SELECT
            p.protein_id,
            p.contig_id,
            p.bin_id,
            COALESCE(p.sequence_length, (p.end_coord - p.start) / 3) as length_aa,
            p.strand,
            (
                SELECT a.name || ' (' || a.accession || ')'
                FROM annotations a
                WHERE a.protein_id = p.protein_id
                ORDER BY a.evalue NULLS LAST
                LIMIT 1
            ) as best_annotation
        FROM proteins p
        WHERE p.protein_id IN ({placeholders})
        ORDER BY length_aa DESC
    """
    rows = store.execute(query, matching_list)

    # Add predicate info (we know they match 'has' and don't match 'lacks')
    rows_with_preds = [
        (r[0], r[1], r[2], r[3], r[4], has, r[5])
        for r in rows
    ]

    return _format_search_results(rows_with_preds, total_count, has, lacks, ctx)


def _format_search_results(
    rows: list,
    total_count: int,
    has: list[str],
    lacks: list[str],
    ctx: OperatorContext,
) -> SharurResult:
    """Format search results into SharurResult."""

    # Build header
    filter_desc = []
    if has:
        filter_desc.append(f"has: {', '.join(has)}")
    if lacks:
        filter_desc.append(f"lacks: {', '.join(lacks)}")

    lines = [
        "# Predicate Search Results",
        f"**Filters:** {' | '.join(filter_desc) or 'none'}",
        f"**Found:** {total_count:,} proteins",
        "",
        "| Protein ID | Length | Strand | Predicates | Annotation |",
        "|------------|--------|--------|------------|------------|",
    ]

    proteins = []
    for row in rows:
        protein_id, contig_id, bin_id, length_aa, strand, predicates, annotation = row
        id_short = protein_id[:22] if len(protein_id) > 22 else protein_id

        # Handle predicates as list or already formatted
        if isinstance(predicates, list):
            pred_str = ", ".join(predicates[:2])
            if len(predicates) > 2:
                pred_str += "..."
        else:
            pred_str = str(predicates)[:15] if predicates else ""

        ann_str = (annotation or "NO HITS")[:25]

        lines.append(
            f"| {id_short} | {length_aa:>6}aa | {strand:>6} | {pred_str[:15]} | {ann_str} |"
        )

        proteins.append({
            "protein_id": protein_id,
            "contig_id": contig_id,
            "bin_id": bin_id,
            "length_aa": length_aa,
            "strand": strand,
            "predicates": list(predicates) if isinstance(predicates, (list, tuple)) else [],
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


def search_proteins(
    store: "DuckDBStore",
    annotation_pattern: Optional[str] = None,
    accession: Optional[str] = None,
    taxonomy_filter: Optional[str] = None,
    min_length: Optional[int] = None,
    max_length: Optional[int] = None,
    limit: int = 50,
    offset: int = 0,
) -> SharurResult:
    """
    Search proteins by annotation, accession, or taxonomy.

    Args:
        store: DuckDB store
        annotation_pattern: Pattern to match in annotation name/description
        accession: Exact accession to match (e.g., "PF00142")
        taxonomy_filter: Filter by genome taxonomy
        min_length: Minimum protein length (aa)
        max_length: Maximum protein length (aa)
        limit: Maximum results
        offset: Pagination offset

    Returns:
        SharurResult with matching proteins
    """
    params = {
        "annotation_pattern": annotation_pattern,
        "accession": accession,
        "taxonomy_filter": taxonomy_filter,
        "min_length": min_length,
        "max_length": max_length,
        "limit": limit,
        "offset": offset,
    }

    with OperatorContext("search_proteins", params) as ctx:
        clauses = []
        query_params = []

        # Build base query depending on filters
        if annotation_pattern or accession:
            # Join with annotations
            base = """
                SELECT DISTINCT
                    p.protein_id,
                    p.contig_id,
                    p.bin_id,
                    COALESCE(p.sequence_length, (p.end_coord - p.start) / 3) as length_aa,
                    p.strand,
                    a.name || ' (' || a.accession || ')' as annotation
                FROM proteins p
                JOIN annotations a ON p.protein_id = a.protein_id
            """

            if annotation_pattern:
                clauses.append(
                    "(LOWER(a.name) LIKE LOWER(?) OR LOWER(a.description) LIKE LOWER(?))"
                )
                pattern = f"%{annotation_pattern}%"
                query_params.extend([pattern, pattern])

            if accession:
                clauses.append("a.accession = ?")
                query_params.append(accession)
        else:
            # No annotation filter, simpler query
            base = """
                SELECT
                    p.protein_id,
                    p.contig_id,
                    p.bin_id,
                    COALESCE(p.sequence_length, (p.end_coord - p.start) / 3) as length_aa,
                    p.strand,
                    (
                        SELECT a.name || ' (' || a.accession || ')'
                        FROM annotations a
                        WHERE a.protein_id = p.protein_id
                        ORDER BY a.evalue NULLS LAST
                        LIMIT 1
                    ) as annotation
                FROM proteins p
            """

        # Add taxonomy filter (requires join with bins)
        if taxonomy_filter:
            if "JOIN bins" not in base:
                base = base.replace("FROM proteins p", "FROM proteins p JOIN bins b ON p.bin_id = b.bin_id")
            clauses.append("LOWER(b.taxonomy) LIKE LOWER(?)")
            query_params.append(f"%{taxonomy_filter}%")

        # Length filters
        if min_length is not None:
            clauses.append("(p.sequence_length >= ? OR (p.end_coord - p.start) / 3 >= ?)")
            query_params.extend([min_length, min_length])

        if max_length is not None:
            clauses.append("(p.sequence_length <= ? OR (p.end_coord - p.start) / 3 <= ?)")
            query_params.extend([max_length, max_length])

        where = f"WHERE {' AND '.join(clauses)}" if clauses else ""

        # Get total count
        count_query = f"""
            SELECT COUNT(*) FROM (
                {base} {where}
            ) sub
        """
        total_count = store.execute(count_query, query_params or None)[0][0]

        # Get results
        query = f"""
            {base}
            {where}
            ORDER BY length_aa DESC
            LIMIT ? OFFSET ?
        """
        query_params.extend([limit, offset])
        rows = store.execute(query, query_params)

        # Format results
        lines = [
            "# Protein Search Results",
            f"**Found:** {total_count:,} proteins",
            "",
            "| Protein ID | Length | Strand | Annotation |",
            "|------------|--------|--------|------------|",
        ]

        proteins = []
        for row in rows:
            protein_id, contig_id, bin_id, length_aa, strand, annotation = row
            id_short = protein_id[:25] if len(protein_id) > 25 else protein_id
            ann_str = (annotation or "NO HITS")[:35]

            lines.append(f"| {id_short} | {length_aa:>6}aa | {strand:>6} | {ann_str} |")

            proteins.append({
                "protein_id": protein_id,
                "contig_id": contig_id,
                "bin_id": bin_id,
                "length_aa": length_aa,
                "strand": strand,
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


__all__ = ["search_by_predicates", "search_proteins"]
