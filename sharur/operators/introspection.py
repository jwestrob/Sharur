"""
Introspection operators for dataset overview and schema description.

These operators provide high-level summaries of the dataset without
requiring specific entity queries.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from sharur.operators.base import SharurResult, OperatorContext
from sharur.predicates.registry import get_registry
from sharur.predicates.vocabulary import (
    PREDICATE_BY_ID,
    PREDICATES_BY_CATEGORY,
)

if TYPE_CHECKING:
    from sharur.storage.duckdb_store import DuckDBStore


def overview(store: "DuckDBStore") -> SharurResult:
    """
    Generate a dataset overview (~400-600 tokens).

    Returns summary statistics including:
    - Genome and protein counts
    - Annotation coverage
    - Taxonomy distribution (top phyla)
    - Predicate summary
    """
    with OperatorContext("overview", {}, store=store) as ctx:
        stats = _gather_stats(store)
        taxonomy = _gather_taxonomy(store)
        annotation_summary = _gather_annotation_summary(store)
        predicate_summary = _gather_predicate_summary(store)

        lines = [
            "# Dataset Overview",
            "",
            "## Summary",
            f"- **Genomes (MAGs):** {stats['genome_count']:,}",
            f"- **Contigs:** {stats['contig_count']:,}",
            f"- **Proteins:** {stats['protein_count']:,}",
            f"- **Annotations:** {stats['annotation_count']:,}",
            f"- **Annotation rate:** {stats['annotation_rate']:.1%}",
            "",
        ]

        # Taxonomy distribution
        if taxonomy:
            lines.extend([
                "## Taxonomy Distribution (Top Phyla)",
            ])
            for phylum, count in taxonomy[:5]:
                pct = count / stats["genome_count"] * 100 if stats["genome_count"] else 0
                lines.append(f"- {phylum}: {count:,} ({pct:.1f}%)")
            lines.append("")

        # Annotation summary
        if annotation_summary:
            lines.extend([
                "## Annotation Sources",
            ])
            for source, count in annotation_summary[:5]:
                lines.append(f"- {source}: {count:,}")
            lines.append("")

        # Predicate summary
        if predicate_summary:
            lines.extend([
                "## Predicate Summary",
            ])
            for category, preds in predicate_summary.items():
                lines.append(f"**{category.title()}:**")
                for pred_id, count in preds[:3]:
                    lines.append(f"  - {pred_id}: {count:,}")
            lines.append("")

        # Quality metrics
        if stats.get("completeness_mean") is not None:
            lines.extend([
                "## Genome Quality",
                f"- Mean completeness: {stats['completeness_mean']:.1f}%",
                f"- Mean contamination: {stats['contamination_mean']:.1f}%",
                f"- High-quality (>90% complete, <5% contam): {stats['hq_count']:,}",
                "",
            ])

        data = "\n".join(lines)

        return ctx.make_result(
            data=data,
            rows=1,
            total_rows=1,
            raw=stats,
        )


def describe_schema(store: "DuckDBStore") -> SharurResult:
    """
    Describe the database schema and available tables.

    Useful for understanding what data is available.
    """
    with OperatorContext("describe_schema", {}, store=store) as ctx:
        # Get table info
        tables = store.execute(
            """
            SELECT table_name, column_count
            FROM (
                SELECT table_name, COUNT(*) as column_count
                FROM information_schema.columns
                WHERE table_schema = 'main'
                GROUP BY table_name
            )
            ORDER BY table_name
            """
        )

        lines = [
            "# Database Schema",
            "",
        ]

        table_summaries = []
        for table_name, col_count in tables:
            # Get row count
            try:
                count = store.execute(f"SELECT COUNT(*) FROM {table_name}")[0][0]
                count_status = "ok"
                count_error = None
            except Exception as exc:
                count = None
                count_status = "unavailable"
                count_error = f"{type(exc).__name__}: {exc}"

            lines.append(f"## {table_name}")
            lines.append(f"- Columns: {col_count}")
            if count is None:
                lines.append(f"- Rows: unavailable ({count_error})")
            else:
                lines.append(f"- Rows: {count:,}")
            lines.append("")
            table_summaries.append(
                {
                    "table_name": table_name,
                    "column_count": col_count,
                    "row_count": count,
                    "status": count_status,
                    "error": count_error,
                }
            )

        lines.extend([
            "# Available Predicate Vocabulary",
            "",
            f"- {len(PREDICATE_BY_ID):,} declared V2-compatible predicates",
            "- Runtime direct-access terms such as `pfam:*`, `kegg:*`, and "
            "`cazy:*` are discovered from the live data and are not a closed list.",
            "",
        ])
        predicate_categories: dict[str, list[str]] = {}
        for category, entries in sorted(PREDICATES_BY_CATEGORY.items()):
            predicate_ids = [entry.predicate_id for entry in entries]
            predicate_categories[category] = predicate_ids
            preview = ", ".join(predicate_ids[:12])
            if len(predicate_ids) > 12:
                preview += ", …"
            lines.append(f"- **{category} ({len(predicate_ids)}):** {preview}")

        lines.append("")
        data = "\n".join(lines)

        return ctx.make_result(
            data=data,
            rows=len(tables) + len(PREDICATE_BY_ID),
            total_rows=len(tables) + len(PREDICATE_BY_ID),
            raw={
                "tables": table_summaries,
                "predicates": sorted(PREDICATE_BY_ID),
                "predicate_categories": predicate_categories,
            },
        )


def _gather_stats(store: "DuckDBStore") -> dict[str, Any]:
    """Gather basic dataset statistics."""
    stats: dict[str, Any] = {}

    # Counts
    try:
        stats["genome_count"] = store.execute("SELECT COUNT(*) FROM bins")[0][0]
    except Exception:
        stats["genome_count"] = 0

    try:
        stats["contig_count"] = store.execute("SELECT COUNT(*) FROM contigs")[0][0]
    except Exception:
        stats["contig_count"] = 0

    try:
        stats["protein_count"] = store.execute("SELECT COUNT(*) FROM proteins")[0][0]
    except Exception:
        stats["protein_count"] = 0

    try:
        stats["annotation_count"] = store.execute("SELECT COUNT(*) FROM annotations")[0][0]
    except Exception:
        stats["annotation_count"] = 0

    # Annotation rate
    try:
        annotated = store.execute(
            """
            SELECT COUNT(DISTINCT protein_id) FROM annotations
            """
        )[0][0]
        stats["annotation_rate"] = (
            annotated / stats["protein_count"] if stats["protein_count"] else 0
        )
    except Exception:
        stats["annotation_rate"] = 0

    # Quality metrics
    try:
        quality = store.execute(
            """
            SELECT
                AVG(completeness) as mean_comp,
                AVG(contamination) as mean_contam,
                SUM(CASE WHEN completeness > 90 AND contamination < 5 THEN 1 ELSE 0 END) as hq
            FROM bins
            WHERE completeness IS NOT NULL
            """
        )[0]
        stats["completeness_mean"] = quality[0]
        stats["contamination_mean"] = quality[1]
        stats["hq_count"] = quality[2] or 0
    except Exception:
        stats["completeness_mean"] = None
        stats["contamination_mean"] = None
        stats["hq_count"] = 0

    return stats


def _gather_taxonomy(store: "DuckDBStore") -> list[tuple[str, int]]:
    """Gather taxonomy distribution by phylum."""
    try:
        # Extract phylum from GTDB taxonomy string (e.g., "d__Bacteria;p__Proteobacteria;...")
        rows = store.execute(
            """
            SELECT
                CASE
                    WHEN taxonomy LIKE '%p__%' THEN
                        REGEXP_EXTRACT(taxonomy, 'p__([^;]+)', 1)
                    ELSE 'Unknown'
                END as phylum,
                COUNT(*) as count
            FROM bins
            WHERE taxonomy IS NOT NULL
            GROUP BY phylum
            ORDER BY count DESC
            LIMIT 10
            """
        )
        return [(row[0] or "Unknown", row[1]) for row in rows]
    except Exception:
        return []


def _gather_annotation_summary(store: "DuckDBStore") -> list[tuple[str, int]]:
    """Gather annotation counts by source."""
    try:
        rows = store.execute(
            """
            SELECT source, COUNT(*) as count
            FROM annotations
            GROUP BY source
            ORDER BY count DESC
            """
        )
        return list(rows)
    except Exception:
        return []


def _gather_predicate_summary(store: "DuckDBStore") -> dict[str, list[tuple[str, int]]]:
    """Gather a compact summary from the active persisted predicate surface."""
    summary: dict[str, list[tuple[str, int]]] = {}

    summary_predicates = [
        "giant",
        "massive",
        "unannotated",
        "well_annotated",
        "multi_domain",
        "defense_system_validated",
        "secretion_system_validated",
    ]

    try:
        expressions = ", ".join(
            f"SUM(CASE WHEN list_contains(predicates, '{predicate}') "
            f"THEN 1 ELSE 0 END) AS {predicate}"
            for predicate in summary_predicates
        )
        row = store.execute(f"SELECT {expressions} FROM protein_predicates")[0]
        for predicate_id, count in zip(summary_predicates, row):
            definition = PREDICATE_BY_ID.get(predicate_id)
            category = definition.category if definition else "runtime"
            summary.setdefault(category, []).append(
                (predicate_id, int(count or 0))
            )
    except Exception:
        registry = get_registry()
        for pred in registry.list_predicates():
            try:
                if pred.eval_query:
                    count = store.execute(
                        f"SELECT COUNT(*) FROM ({pred.eval_query.strip()}) sub"
                    )[0][0]
                else:
                    count = 0
                summary.setdefault(pred.category, []).append(
                    (pred.predicate_id, count)
                )
            except Exception:
                continue

    for category in summary:
        summary[category].sort(key=lambda item: -item[1])

    return summary


__all__ = ["overview", "describe_schema"]
