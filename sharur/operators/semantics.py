"""Authoritative semantic reads shared by navigation operators.

The persisted V2 compatibility cache is the normal read path. Databases that
have V2 state but no cache row are resolved from V2 in memory; legacy dynamic
evaluation is used only when no V2 state is available for that protein.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from sharur.operators.predicates_v2 import get_atoms, get_semantic_state
from sharur.predicates.evaluator import compute_predicates_for_protein
from sharur.predicates_v2.compat import semantic_state_to_predicates


if TYPE_CHECKING:
    from collections.abc import Iterable

    from sharur.storage.duckdb_store import DuckDBStore


def get_active_predicates(
    store: DuckDBStore,
    protein_ids: Iterable[str],
) -> dict[str, list[str]]:
    """Return the active predicate view for each requested protein."""
    ordered_ids = list(dict.fromkeys(str(protein_id) for protein_id in protein_ids))
    if not ordered_ids:
        return {}

    table_names = _table_names(store)
    predicates_by_protein: dict[str, list[str]] = {}

    if "protein_predicates" in table_names:
        placeholders = ",".join("?" for _ in ordered_ids)
        rows = store.execute(
            f"""
            SELECT protein_id, predicates
            FROM protein_predicates
            WHERE protein_id IN ({placeholders})
            """,
            ordered_ids,
        )
        predicates_by_protein.update(
            {str(protein_id): sorted(predicates or []) for protein_id, predicates in rows}
        )

    missing_ids = [
        protein_id for protein_id in ordered_ids if protein_id not in predicates_by_protein
    ]

    if missing_ids and "semantic_state" in table_names:
        still_missing: list[str] = []
        for protein_id in missing_ids:
            state = get_semantic_state(store, protein_id)
            if state is None:
                still_missing.append(protein_id)
                continue
            atoms = get_atoms(store, protein_id) if "semantic_atoms" in table_names else []
            predicates_by_protein[protein_id] = semantic_state_to_predicates(
                state,
                atoms=atoms,
            )
        missing_ids = still_missing

    for protein_id in missing_ids:
        predicates_by_protein[protein_id] = sorted(
            compute_predicates_for_protein(protein_id, store)
        )

    return {protein_id: predicates_by_protein.get(protein_id, []) for protein_id in ordered_ids}


def get_active_predicates_for_protein(
    store: DuckDBStore,
    protein_id: str,
) -> list[str]:
    """Return the authoritative flat compatibility view for one protein."""
    return get_active_predicates(store, [protein_id]).get(protein_id, [])


def _table_names(store: DuckDBStore) -> set[str]:
    cached = getattr(store, "_sharur_semantic_table_names", None)
    if cached is not None:
        return cached
    try:
        table_names = {str(row[0]) for row in store.execute("SHOW TABLES")}
    except Exception:
        return set()
    if getattr(store, "read_only", False):
        store._sharur_semantic_table_names = table_names
    return table_names


__all__ = ["get_active_predicates", "get_active_predicates_for_protein"]
