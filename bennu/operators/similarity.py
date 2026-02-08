"""
Similarity search operators using vector embeddings.

Uses LanceDB with ESM2 protein embeddings for finding structurally/functionally
similar proteins across the database.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Optional

from bennu.operators.base import BennuResult, OperatorContext

if TYPE_CHECKING:
    from bennu.core.session import ExplorationSession
    from bennu.storage.duckdb_store import DuckDBStore


def find_similar(
    store: "DuckDBStore",
    vector_store,
    protein_id: str,
    k: int = 10,
    threshold: Optional[float] = None,
    include_self_genome: bool = True,
) -> BennuResult:
    """
    Find proteins with similar ESM2 embeddings.

    Uses LanceDB kNN search to find proteins with similar sequence embeddings,
    which often correlates with structural and functional similarity.

    Args:
        store: DuckDB store for protein metadata
        vector_store: LanceDB vector store with embeddings
        protein_id: Query protein ID
        k: Number of similar proteins to return
        threshold: Minimum similarity score (0-1, cosine similarity)
        include_self_genome: Include hits from the same genome

    Returns:
        BennuResult with similar proteins ranked by similarity
    """
    params = {
        "protein_id": protein_id,
        "k": k,
        "threshold": threshold,
        "include_self_genome": include_self_genome,
    }

    with OperatorContext("find_similar", params) as ctx:
        if vector_store is None:
            return ctx.make_result(
                data="Vector store not available. Embeddings may not be computed for this database.",
                rows=0,
            )

        # Get query protein's genome for filtering
        query_genome = None
        if not include_self_genome:
            rows = store.execute(
                "SELECT bin_id FROM proteins WHERE protein_id = ?",
                [protein_id],
            )
            if rows:
                query_genome = rows[0][0]

        # Find similar proteins
        try:
            # Request more than k to allow filtering
            results = vector_store.query(
                protein_id,
                k=k * 3 if not include_self_genome else k + 1,
                threshold=threshold,
                include_distances=True,
            )
        except Exception as e:
            return ctx.make_result(
                data=f"Error querying vector store: {e}",
                rows=0,
            )

        if not results:
            return ctx.make_result(
                data=f"No similar proteins found for {protein_id}",
                rows=0,
            )

        # Get protein IDs from results
        similar_ids = [r[0] if isinstance(r, tuple) else r for r in results]
        scores = {r[0]: r[1] for r in results} if results and isinstance(results[0], tuple) else {}

        # Get metadata for similar proteins
        placeholders = ",".join(["?"] * len(similar_ids))
        query = f"""
            SELECT
                p.protein_id,
                p.bin_id,
                p.contig_id,
                COALESCE(p.sequence_length, (p.end_coord - p.start) / 3) as length_aa,
                (
                    SELECT a.name || ' (' || a.accession || ')'
                    FROM annotations a
                    WHERE a.protein_id = p.protein_id
                    ORDER BY a.evalue NULLS LAST
                    LIMIT 1
                ) as annotation
            FROM proteins p
            WHERE p.protein_id IN ({placeholders})
        """
        rows = store.execute(query, similar_ids)

        # Build result with scores
        protein_data = {row[0]: row for row in rows}

        lines = [
            f"# Similar Proteins to {protein_id[:40]}...",
            f"**Method:** ESM2 embedding similarity (cosine)",
            f"**Found:** {len(results)} similar proteins",
            "",
            "| Rank | Similarity | Length | Genome | Annotation |",
            "|------|------------|--------|--------|------------|",
        ]

        output_proteins = []
        rank = 0
        for pid, score in results:
            if pid == protein_id:
                continue  # Skip self
            if pid not in protein_data:
                continue

            _, bin_id, contig_id, length, annotation = protein_data[pid]

            # Filter by genome if requested
            if not include_self_genome and bin_id == query_genome:
                continue

            rank += 1
            if rank > k:
                break

            genome_short = (bin_id or "")[:20]
            ann_short = (annotation or "NO HITS")[:30]

            lines.append(
                f"| {rank} | {score:.3f} | {length:>6}aa | {genome_short} | {ann_short} |"
            )

            output_proteins.append({
                "rank": rank,
                "protein_id": pid,
                "similarity": score,
                "bin_id": bin_id,
                "contig_id": contig_id,
                "length_aa": length,
                "annotation": annotation,
            })

        if not output_proteins:
            return ctx.make_result(
                data=f"No similar proteins found for {protein_id} with given filters",
                rows=0,
            )

        return ctx.make_result(
            data="\n".join(lines),
            rows=len(output_proteins),
            raw=output_proteins,
        )


def find_similar_to_set(
    store: "DuckDBStore",
    vector_store,
    protein_ids: list[str],
    k: int = 10,
    threshold: Optional[float] = 0.7,
) -> BennuResult:
    """
    Find proteins similar to ANY in a set (union of similar proteins).

    Useful for expanding a working set or finding related proteins.

    Args:
        store: DuckDB store
        vector_store: LanceDB vector store
        protein_ids: Set of query protein IDs
        k: Number of similar proteins per query
        threshold: Minimum similarity score

    Returns:
        BennuResult with proteins similar to the set
    """
    params = {
        "protein_ids": protein_ids,
        "k": k,
        "threshold": threshold,
    }

    with OperatorContext("find_similar_to_set", params) as ctx:
        if vector_store is None:
            return ctx.make_result(
                data="Vector store not available.",
                rows=0,
            )

        # Collect all similar proteins
        all_similar: dict[str, float] = {}  # pid -> max similarity
        query_set = set(protein_ids)

        for pid in protein_ids:
            try:
                results = vector_store.query(
                    pid,
                    k=k * 2,
                    threshold=threshold,
                    include_distances=True,
                )
                for similar_id, score in results:
                    if similar_id not in query_set:
                        all_similar[similar_id] = max(
                            all_similar.get(similar_id, 0),
                            score,
                        )
            except Exception:
                continue

        if not all_similar:
            return ctx.make_result(
                data="No similar proteins found outside the query set",
                rows=0,
            )

        # Sort by score and take top k
        sorted_similar = sorted(
            all_similar.items(),
            key=lambda x: x[1],
            reverse=True,
        )[:k]

        # Get metadata
        similar_ids = [pid for pid, _ in sorted_similar]
        placeholders = ",".join(["?"] * len(similar_ids))
        query = f"""
            SELECT
                p.protein_id,
                p.bin_id,
                COALESCE(p.sequence_length, (p.end_coord - p.start) / 3) as length_aa,
                (
                    SELECT a.name
                    FROM annotations a
                    WHERE a.protein_id = p.protein_id
                    ORDER BY a.evalue NULLS LAST
                    LIMIT 1
                ) as annotation
            FROM proteins p
            WHERE p.protein_id IN ({placeholders})
        """
        rows = store.execute(query, similar_ids)
        protein_data = {row[0]: row for row in rows}

        lines = [
            f"# Proteins Similar to Set ({len(protein_ids)} queries)",
            f"**Threshold:** {threshold}",
            f"**Found:** {len(sorted_similar)} unique similar proteins",
            "",
            "| Rank | Similarity | Length | Annotation |",
            "|------|------------|--------|------------|",
        ]

        for rank, (pid, score) in enumerate(sorted_similar, 1):
            if pid in protein_data:
                _, bin_id, length, annotation = protein_data[pid]
                ann_short = (annotation or "NO HITS")[:35]
                lines.append(f"| {rank} | {score:.3f} | {length:>6}aa | {ann_short} |")

        return ctx.make_result(
            data="\n".join(lines),
            rows=len(sorted_similar),
            raw=[{"protein_id": pid, "similarity": score} for pid, score in sorted_similar],
        )


__all__ = ["find_similar", "find_similar_to_set"]
