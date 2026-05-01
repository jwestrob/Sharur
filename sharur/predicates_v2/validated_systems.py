"""Validated defense/secretion system evidence for V2 atom generation."""

from __future__ import annotations

import re
from collections import defaultdict
from typing import TYPE_CHECKING

from sharur.predicates.generator import AnnotationRecord


if TYPE_CHECKING:
    from sharur.storage.duckdb_store import DuckDBStore


def fetch_validated_system_annotations(
    store: DuckDBStore,
    protein_ids: set[str] | None = None,
) -> dict[str, list[AnnotationRecord]]:
    """Return synthetic annotation records from validated system tables.

    The ingest validation path writes system summaries into ``defense_systems``
    and ``secretion_systems`` and may also write per-gene annotation rows.
    V2 should not depend on those annotation rows being present: the validated
    tables are the durable authority for system calls.
    """
    annotations: dict[str, list[AnnotationRecord]] = defaultdict(list)
    available_tables = _table_names(store)

    if "system_proteins" in available_tables and _system_proteins_has_rows(store):
        if "defense_systems" in available_tables:
            _collect_normalized_system_annotations(
                store=store,
                table_name="defense_systems",
                source="defensefinder_system",
                annotations=annotations,
                protein_ids=protein_ids,
            )

        if "secretion_systems" in available_tables:
            _collect_normalized_system_annotations(
                store=store,
                table_name="secretion_systems",
                source="txsscan_system",
                annotations=annotations,
                protein_ids=protein_ids,
            )

        return dict(annotations)

    if "defense_systems" in available_tables:
        _collect_system_annotations(
            store=store,
            table_name="defense_systems",
            source="defensefinder_system",
            annotations=annotations,
            protein_ids=protein_ids,
        )

    if "secretion_systems" in available_tables:
        _collect_system_annotations(
            store=store,
            table_name="secretion_systems",
            source="txsscan_system",
            annotations=annotations,
            protein_ids=protein_ids,
        )

    return dict(annotations)


def materialize_system_proteins(store: DuckDBStore) -> int:
    """Populate normalized member rows from validated system summary tables."""
    available_tables = _table_names(store)
    if "defense_systems" not in available_tables and "secretion_systems" not in available_tables:
        return 0

    _ensure_system_proteins_table(store)
    rows: list[tuple[str, str, str, int, str | None, float | None]] = []

    if "defense_systems" in available_tables:
        store.execute(
            "DELETE FROM system_proteins WHERE system_source = 'defensefinder_system'"
        )
        rows.extend(
            _system_member_rows(
                store=store,
                table_name="defense_systems",
                source="defensefinder_system",
            )
        )

    if "secretion_systems" in available_tables:
        store.execute("DELETE FROM system_proteins WHERE system_source = 'txsscan_system'")
        rows.extend(
            _system_member_rows(
                store=store,
                table_name="secretion_systems",
                source="txsscan_system",
            )
        )

    if rows:
        store.conn.executemany(
            """
            INSERT INTO system_proteins (
                system_id, protein_id, system_source, position, profile_name, score
            )
            VALUES (?, ?, ?, ?, ?, ?)
            """,
            rows,
        )

    return len(rows)


def _table_names(store: DuckDBStore) -> set[str]:
    """Return existing DuckDB table names."""
    return {row[0] for row in store.execute("SHOW TABLES")}


def _ensure_system_proteins_table(store: DuckDBStore) -> None:
    """Create normalized validated-system membership table if needed."""
    store.execute("""
        CREATE TABLE IF NOT EXISTS system_proteins (
            system_id VARCHAR NOT NULL,
            protein_id VARCHAR NOT NULL,
            system_source VARCHAR NOT NULL,
            position INTEGER,
            profile_name VARCHAR,
            score DOUBLE,
            PRIMARY KEY (system_id, protein_id, system_source),
            FOREIGN KEY (protein_id) REFERENCES proteins(protein_id)
        );
    """)
    store.execute(
        "CREATE INDEX IF NOT EXISTS idx_system_proteins_protein "
        "ON system_proteins(protein_id);"
    )
    store.execute(
        "CREATE INDEX IF NOT EXISTS idx_system_proteins_system "
        "ON system_proteins(system_id);"
    )
    store.execute(
        "CREATE INDEX IF NOT EXISTS idx_system_proteins_source "
        "ON system_proteins(system_source);"
    )


def _system_proteins_has_rows(store: DuckDBStore) -> bool:
    """Return whether normalized system membership has any rows."""
    return bool(store.execute("SELECT 1 FROM system_proteins LIMIT 1"))


def _collect_normalized_system_annotations(
    *,
    store: DuckDBStore,
    table_name: str,
    source: str,
    annotations: dict[str, list[AnnotationRecord]],
    protein_ids: set[str] | None,
) -> None:
    """Load synthetic annotations from normalized validated-system members."""
    protein_filter = ""
    params: list[object] = [source]
    if protein_ids is not None:
        if not protein_ids:
            return
        ids = sorted(protein_ids)
        placeholders = ",".join(["?"] * len(ids))
        protein_filter = f"AND sp.protein_id IN ({placeholders})"
        params.extend(ids)

    rows = store.execute(
        f"""
        SELECT
            sp.protein_id,
            sp.system_id,
            st.system_type,
            st.system_subtype,
            st.genes_count,
            sp.profile_name,
            sp.score
        FROM system_proteins sp
        JOIN {table_name} st ON st.system_id = sp.system_id
        WHERE sp.system_source = ?
        {protein_filter}
        """,
        params,
    )

    for protein_id, system_id, system_type, system_subtype, genes_count, profile_name, score in rows:
        name = _system_name(system_type, system_subtype)
        description = f"Validated system: {system_id}"
        if profile_name:
            description = f"{description}; member profile: {profile_name}"

        annotations[protein_id].append(
            AnnotationRecord(
                source=source,
                accession=str(system_id or name),
                name=name,
                description=description,
                evalue=None,
                score=score if score is not None else (
                    float(genes_count) if genes_count is not None else None
                ),
            )
        )


def _system_member_rows(
    *,
    store: DuckDBStore,
    table_name: str,
    source: str,
) -> list[tuple[str, str, str, int, str | None, float | None]]:
    """Parse summary table membership strings into normalized rows."""
    rows = store.execute(
        f"""
        SELECT system_id, system_type, system_subtype, genes_count,
               protein_ids, profile_names
        FROM {table_name}
        WHERE protein_ids IS NOT NULL AND protein_ids != ''
        """
    )

    normalized: list[tuple[str, str, str, int, str | None, float | None]] = []
    seen: set[tuple[str, str, str]] = set()
    for system_id, system_type, system_subtype, genes_count, protein_text, profile_text in rows:
        system_name = str(system_id or _system_name(system_type, system_subtype))
        score = float(genes_count) if genes_count is not None else None
        profiles = _split_protein_ids(profile_text)

        for position, protein_id in enumerate(_split_protein_ids(protein_text), start=1):
            key = (system_name, protein_id, source)
            if key in seen:
                continue
            seen.add(key)
            profile_name = profiles[position - 1] if position <= len(profiles) else None
            normalized.append(
                (system_name, protein_id, source, position, profile_name, score)
            )

    return normalized


def _collect_system_annotations(
    *,
    store: DuckDBStore,
    table_name: str,
    source: str,
    annotations: dict[str, list[AnnotationRecord]],
    protein_ids: set[str] | None,
) -> None:
    """Load one validated system table into synthetic annotations."""
    rows = store.execute(
        f"""
        SELECT system_id, system_type, system_subtype, genes_count, protein_ids
        FROM {table_name}
        WHERE protein_ids IS NOT NULL AND protein_ids != ''
        """
    )

    for system_id, system_type, system_subtype, genes_count, protein_text in rows:
        name = _system_name(system_type, system_subtype)
        description = f"Validated system: {system_id}"
        score = float(genes_count) if genes_count is not None else None

        for pid in _split_protein_ids(protein_text):
            if protein_ids is not None and pid not in protein_ids:
                continue

            annotations[pid].append(
                AnnotationRecord(
                    source=source,
                    accession=str(system_id or name),
                    name=name,
                    description=description,
                    evalue=None,
                    score=score,
                )
            )


def _system_name(system_type: str | None, system_subtype: str | None) -> str:
    """Format system type/subtype like existing system annotation rows."""
    system_type = (system_type or "").strip()
    system_subtype = (system_subtype or "").strip()
    if system_type and system_subtype and system_subtype != system_type:
        return f"{system_type}/{system_subtype}"
    return system_type or system_subtype


def _split_protein_ids(protein_text: object) -> list[str]:
    """Parse comma/semicolon/whitespace-delimited protein ID lists."""
    if protein_text is None:
        return []
    if isinstance(protein_text, list | tuple):
        raw_items = [str(item) for item in protein_text]
    else:
        raw = str(protein_text).strip()
        if not raw:
            return []
        raw_items = re.split(r"[,;\s]+", raw.strip("[](){}"))

    protein_ids: list[str] = []
    for item in raw_items:
        pid = item.strip().strip("'\"")
        if pid:
            protein_ids.append(pid)
    return protein_ids


__all__ = ["fetch_validated_system_annotations", "materialize_system_proteins"]
