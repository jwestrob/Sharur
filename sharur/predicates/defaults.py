"""
Default predicate definitions for Sharur.

These predicates cover common filtering use cases:
- Size-based: giant, massive, tiny
- Annotation-based: unannotated, hypothetical, confident_hit
- Composition-based: gc_outlier
"""

from sharur.predicates.registry import PredicateDefinition, PredicateRegistry


def register_defaults(registry: PredicateRegistry) -> None:
    """Register default predicates with the registry."""

    # Size predicates
    registry.register(
        PredicateDefinition(
            predicate_id="giant",
            name="Giant protein",
            description="Protein with length > 2000 amino acids",
            category="size",
            eval_query="""
                SELECT protein_id FROM proteins
                WHERE sequence_length > 2000
                   OR (end_coord - start) / 3 > 2000
            """,
        )
    )

    registry.register(
        PredicateDefinition(
            predicate_id="massive",
            name="Massive protein",
            description="Protein with length > 5000 amino acids (extremely large)",
            category="size",
            eval_query="""
                SELECT protein_id FROM proteins
                WHERE sequence_length > 5000
                   OR (end_coord - start) / 3 > 5000
            """,
        )
    )

    registry.register(
        PredicateDefinition(
            predicate_id="tiny",
            name="Tiny protein",
            description="Protein with length < 50 amino acids (small ORF)",
            category="size",
            eval_query="""
                SELECT protein_id FROM proteins
                WHERE sequence_length < 50
                   OR (sequence_length IS NULL AND (end_coord - start) / 3 < 50)
            """,
        )
    )

    # Annotation predicates
    registry.register(
        PredicateDefinition(
            predicate_id="unannotated",
            name="Unannotated",
            description="Protein with no functional annotations",
            category="annotation",
            eval_query="""
                SELECT p.protein_id FROM proteins p
                LEFT JOIN annotations a ON p.protein_id = a.protein_id
                WHERE a.protein_id IS NULL
            """,
        )
    )

    registry.register(
        PredicateDefinition(
            predicate_id="hypothetical",
            name="Hypothetical",
            description="Best annotation contains 'hypothetical' or similar",
            category="annotation",
            eval_query="""
                SELECT DISTINCT a.protein_id FROM annotations a
                WHERE LOWER(a.name) LIKE '%hypothetical%'
                   OR LOWER(a.description) LIKE '%hypothetical%'
                   OR LOWER(a.name) LIKE '%uncharacterized%'
                   OR LOWER(a.description) LIKE '%uncharacterized%'
                   OR LOWER(a.name) LIKE '%unknown function%'
                   OR LOWER(a.description) LIKE '%unknown function%'
            """,
        )
    )

    registry.register(
        PredicateDefinition(
            predicate_id="confident_hit",
            name="Confident annotation",
            description="Has at least one annotation with e-value < 1e-10",
            category="annotation",
            eval_query="""
                SELECT DISTINCT protein_id FROM annotations
                WHERE evalue IS NOT NULL AND evalue < 1e-10
            """,
        )
    )

    registry.register(
        PredicateDefinition(
            predicate_id="multi_domain",
            name="Multi-domain",
            description="Protein with 3+ distinct domain annotations",
            category="annotation",
            eval_query="""
                SELECT protein_id FROM (
                    SELECT protein_id, COUNT(DISTINCT accession) as domain_count
                    FROM annotations
                    WHERE source = 'pfam'
                    GROUP BY protein_id
                ) sub
                WHERE domain_count >= 3
            """,
        )
    )

    # Composition predicates
    registry.register(
        PredicateDefinition(
            predicate_id="gc_outlier",
            name="GC outlier",
            description="Protein with GC content deviating >2 std from genome mean",
            category="composition",
            eval_query="""
                WITH genome_stats AS (
                    SELECT bin_id, AVG(gc_content) as mean_gc, STDDEV(gc_content) as std_gc
                    FROM proteins
                    WHERE gc_content IS NOT NULL AND bin_id IS NOT NULL
                    GROUP BY bin_id
                )
                SELECT p.protein_id FROM proteins p
                JOIN genome_stats gs ON p.bin_id = gs.bin_id
                WHERE gs.std_gc > 0
                  AND ABS(p.gc_content - gs.mean_gc) > 2 * gs.std_gc
            """,
        )
    )


__all__ = ["register_defaults"]
