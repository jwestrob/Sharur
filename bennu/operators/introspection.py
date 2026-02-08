"""
Introspection operators for dataset overview and schema description.

These operators provide high-level summaries of the dataset without
requiring specific entity queries.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from bennu.operators.base import BennuResult, OperatorContext
from bennu.predicates.registry import get_registry

if TYPE_CHECKING:
    from bennu.storage.duckdb_store import DuckDBStore


def overview(store: "DuckDBStore") -> BennuResult:
    """
    Generate a dataset overview (~400-600 tokens).

    Returns summary statistics including:
    - Genome and protein counts
    - Annotation coverage
    - Taxonomy distribution (top phyla)
    - Predicate summary
    """
    with OperatorContext("overview", {}) as ctx:
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


def describe_schema(store: "DuckDBStore") -> BennuResult:
    """
    Describe the database schema and available tables.

    Useful for understanding what data is available.
    """
    with OperatorContext("describe_schema", {}) as ctx:
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

        for table_name, col_count in tables:
            # Get row count
            try:
                count = store.execute(f"SELECT COUNT(*) FROM {table_name}")[0][0]
            except Exception:
                count = 0

            lines.append(f"## {table_name}")
            lines.append(f"- Columns: {col_count}")
            lines.append(f"- Rows: {count:,}")
            lines.append("")

        # List available predicates
        registry = get_registry()
        predicates = registry.list_predicates()

        lines.extend([
            "# Available Predicates",
            "",
        ])
        for pred in predicates:
            lines.append(f"- **{pred.predicate_id}** ({pred.category}): {pred.description}")

        lines.append("")
        data = "\n".join(lines)

        return ctx.make_result(
            data=data,
            rows=len(tables) + len(predicates),
            total_rows=len(tables) + len(predicates),
            raw={"tables": tables, "predicates": [p.predicate_id for p in predicates]},
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
    """Gather predicate counts, grouped by category."""
    registry = get_registry()
    summary: dict[str, list[tuple[str, int]]] = {}

    for pred in registry.list_predicates():
        try:
            # Quick count via eval_query
            if pred.eval_query:
                count = store.execute(
                    f"SELECT COUNT(*) FROM ({pred.eval_query.strip()}) sub"
                )[0][0]
            else:
                count = 0

            if pred.category not in summary:
                summary[pred.category] = []
            summary[pred.category].append((pred.predicate_id, count))
        except Exception:
            # Skip predicates that fail (e.g., missing tables)
            pass

    # Sort each category by count descending
    for cat in summary:
        summary[cat].sort(key=lambda x: -x[1])

    return summary


__all__ = ["overview", "describe_schema"]
