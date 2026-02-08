"""
Predicate evaluator for computing predicates on proteins.

Supports both SQL-based evaluation (for bulk computation) and
Python callable evaluation (for custom predicates).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Optional

from bennu.predicates.registry import PredicateDefinition, get_registry

if TYPE_CHECKING:
    from bennu.storage.duckdb_store import DuckDBStore


def evaluate_predicate(
    predicate_id: str,
    store: "DuckDBStore",
    protein_ids: Optional[list[str]] = None,
) -> set[str]:
    """
    Evaluate a predicate and return protein IDs where it's true.

    Args:
        predicate_id: ID of predicate to evaluate
        store: DuckDB store to query
        protein_ids: Optional subset of proteins to check (None = all)

    Returns:
        Set of protein IDs where predicate is true
    """
    registry = get_registry()
    pred = registry.get(predicate_id)
    if pred is None:
        raise ValueError(f"Unknown predicate: {predicate_id}")

    if pred.eval_query:
        # SQL-based evaluation
        query = pred.eval_query.strip()

        # If filtering to specific proteins, wrap query
        if protein_ids:
            placeholders = ",".join(["?"] * len(protein_ids))
            query = f"""
                SELECT protein_id FROM ({query}) sub
                WHERE protein_id IN ({placeholders})
            """
            rows = store.execute(query, protein_ids)
        else:
            rows = store.execute(query)

        return {row[0] for row in rows}

    elif pred.eval_func:
        # Python callable evaluation
        if protein_ids is None:
            # Get all protein IDs
            rows = store.execute("SELECT protein_id FROM proteins")
            protein_ids = [row[0] for row in rows]

        result = set()
        for pid in protein_ids:
            # Fetch protein data for evaluation
            protein_data = _get_protein_dict(store, pid)
            if protein_data and pred.eval_func(protein_data):
                result.add(pid)
        return result

    raise ValueError(f"Predicate {predicate_id} has no evaluation method")


def compute_predicates_for_protein(
    protein_id: str,
    store: "DuckDBStore",
    predicate_ids: Optional[list[str]] = None,
) -> list[str]:
    """
    Compute all predicates for a single protein.

    Args:
        protein_id: Protein to evaluate
        store: DuckDB store
        predicate_ids: Optional subset of predicates (None = all)

    Returns:
        List of predicate IDs that are true for this protein
    """
    registry = get_registry()

    if predicate_ids is None:
        predicate_ids = [p.predicate_id for p in registry.list_predicates()]

    result = []
    for pred_id in predicate_ids:
        matching = evaluate_predicate(pred_id, store, [protein_id])
        if protein_id in matching:
            result.append(pred_id)

    return result


def compute_all_predicates(
    store: "DuckDBStore",
    predicate_ids: Optional[list[str]] = None,
    batch_size: int = 1000,
) -> dict[str, list[str]]:
    """
    Compute all predicates for all proteins.

    Args:
        store: DuckDB store
        predicate_ids: Optional subset of predicates (None = all)
        batch_size: Unused, kept for API compatibility

    Returns:
        Dict mapping protein_id -> list of true predicate IDs
    """
    registry = get_registry()

    if predicate_ids is None:
        predicate_ids = [p.predicate_id for p in registry.list_predicates()]

    # Build up predicate membership per protein
    protein_predicates: dict[str, list[str]] = {}

    for pred_id in predicate_ids:
        matching = evaluate_predicate(pred_id, store)
        for protein_id in matching:
            if protein_id not in protein_predicates:
                protein_predicates[protein_id] = []
            protein_predicates[protein_id].append(pred_id)

    return protein_predicates


def persist_predicates(
    store: "DuckDBStore",
    protein_predicates: Optional[dict[str, list[str]]] = None,
) -> int:
    """
    Persist computed predicates to the protein_predicates table.

    Args:
        store: DuckDB store
        protein_predicates: Pre-computed predicates, or None to compute fresh

    Returns:
        Number of proteins updated
    """
    if protein_predicates is None:
        protein_predicates = compute_all_predicates(store)

    # Clear existing and insert new
    store.execute("DELETE FROM protein_predicates")

    count = 0
    for protein_id, predicates in protein_predicates.items():
        store.execute(
            """
            INSERT INTO protein_predicates (protein_id, predicates, updated_at)
            VALUES (?, ?, CURRENT_TIMESTAMP)
            """,
            [protein_id, predicates],
        )
        count += 1

    return count


def _get_protein_dict(store: "DuckDBStore", protein_id: str) -> Optional[dict]:
    """Fetch protein data as dict for Python callable evaluation."""
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
        return None

    r = row[0]
    return {
        "protein_id": r[0],
        "contig_id": r[1],
        "bin_id": r[2],
        "start": r[3],
        "end": r[4],
        "strand": r[5],
        "sequence": r[6],
        "sequence_length": r[7],
        "gc_content": r[8],
    }


__all__ = [
    "evaluate_predicate",
    "compute_predicates_for_protein",
    "compute_all_predicates",
    "persist_predicates",
]
