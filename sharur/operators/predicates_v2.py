"""Facade-facing Predicate V2 operations."""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path
from typing import TYPE_CHECKING

from sharur.predicates_v2.model import (
    ClaimRelation,
    SemanticAtom,
    SemanticFacet,
    SemanticState,
)

if TYPE_CHECKING:
    from sharur.storage.duckdb_store import DuckDBStore


def generate_v2(
    store: DuckDBStore,
    protein_ids: list[str] | None = None,
    output_review_queue: str | None = None,
    chunk_size: int = 25_000,
    workers: int = 1,
    worker_batch_size: int | None = None,
    pipeline_depth: int = 2,
    resume: bool = False,
    update_legacy_predicates: bool = False,
    return_states: bool | None = None,
    predict_topology: bool = False,
) -> dict:
    """Generate V2 semantic atoms + states and persist to DuckDB."""
    from sharur.predicates_v2.persistence import generate_and_persist_v2

    if return_states is None:
        return_states = protein_ids is not None

    return generate_and_persist_v2(
        store,
        protein_ids=protein_ids,
        output_review_queue=output_review_queue,
        chunk_size=chunk_size,
        workers=workers,
        worker_batch_size=worker_batch_size,
        pipeline_depth=pipeline_depth,
        resume=resume,
        update_legacy_predicates=update_legacy_predicates,
        return_states=return_states,
        predict_topology=predict_topology,
    )


def get_semantic_state(
    store: DuckDBStore,
    protein_id: str,
) -> SemanticState | None:
    """Get the V2 semantic state for a protein."""
    rows = store.execute(
        "SELECT activities, roles, architecture, localization, "
        "topology, size_class, quality_flags, composite_predicates, "
        "unresolved_count FROM semantic_state WHERE protein_id = ?",
        [protein_id],
    )
    if not rows:
        return None

    row = rows[0]
    topology = json.loads(row[4]) if isinstance(row[4], str) else (row[4] or {})
    return SemanticState(
        protein_id=protein_id,
        activities=list(row[0] or []),
        roles=list(row[1] or []),
        architecture=list(row[2] or []),
        localization=list(row[3] or []),
        topology=topology,
        size_class=row[5] or "",
        quality_flags=list(row[6] or []),
        composite_predicates=list(row[7] or []),
    )


def get_atoms(store: DuckDBStore, protein_id: str) -> list[SemanticAtom]:
    """Get all V2 semantic atoms for a protein."""
    rows = store.execute(
        "SELECT atom_id, facet, relation, source_accession, "
        "source_db, evidence_evalue, evidence_score "
        "FROM semantic_atoms WHERE protein_id = ?",
        [protein_id],
    )
    return [
        SemanticAtom(
            protein_id=protein_id,
            atom_id=row[0],
            facet=SemanticFacet(row[1]),
            relation=ClaimRelation(row[2]),
            source_accession=row[3],
            source_db=row[4],
            evidence_evalue=row[5],
            evidence_score=row[6],
        )
        for row in rows
    ]


def explain(store: DuckDBStore, protein_id: str) -> dict:
    """Explain one protein through the active predicate tables."""
    from sharur.predicates_v2.composites import explain_composites

    table_names = _table_names(store)
    state = (
        get_semantic_state(store, protein_id)
        if "semantic_state" in table_names
        else None
    )
    atoms = get_atoms(store, protein_id) if "semantic_atoms" in table_names else []
    terms = _fetch_semantic_terms(store, protein_id, table_names)

    direct_access_terms = sorted(
        {
            term["term_id"]
            for term in terms
            if term["term_kind"] == "direct_access"
        }
    )
    composite_terms = sorted(
        {
            term["term_id"]
            for term in terms
            if term["term_kind"] == "composite"
        }
    )
    if not composite_terms and state is not None:
        composite_terms = sorted(state.composite_predicates)

    return {
        "protein_id": protein_id,
        "protein": _fetch_protein_summary(store, protein_id, table_names),
        "semantic_state": state.to_dict() if state is not None else None,
        "atoms": [atom.to_dict() for atom in atoms],
        "terms": terms,
        "direct_access_terms": direct_access_terms,
        "composite_terms": composite_terms,
        "composite_explanations": explain_composites(
            atoms,
            topology=state.topology if state is not None else {},
            only=composite_terms,
        ) if atoms and composite_terms else {},
        "unresolved": _unresolved_atoms(atoms),
        "conflicts": _conflicting_atoms(atoms),
        "validated_systems": _fetch_validated_system_membership(
            store,
            protein_id,
            table_names,
        ),
        "legacy_predicates": _fetch_legacy_predicates(store, protein_id, table_names),
    }


def search_by_facet(
    store: DuckDBStore,
    facet: str,
    atom_ids: list[str] | None = None,
    relation: str | None = None,
    source_db: str | None = None,
    source_accession: str | None = None,
    limit: int = 50,
) -> list[tuple]:
    """Search proteins by V2 semantic facet."""
    atom_column = "atom_id"
    table_name = "semantic_atoms"
    params: list[object] = [facet]
    clauses = ["facet = ?"]

    if atom_ids:
        placeholders = ",".join(["?"] * len(atom_ids))
        clauses.append(f"{atom_column} IN ({placeholders})")
        params.extend(atom_ids)

    if relation:
        clauses.append("relation = ?")
        params.append(relation)

    if source_db:
        clauses.append("source_db = ?")
        params.append(source_db)

    if source_accession:
        clauses.append("source_accession = ?")
        params.append(source_accession)

    where = " AND ".join(clauses)
    rows = store.execute(
        f"SELECT DISTINCT protein_id, {atom_column}, relation "
        f"FROM {table_name} WHERE {where} "
        "ORDER BY protein_id LIMIT ?",
        [*params, limit],
    )
    return [tuple(row) for row in rows]


def search_atoms(
    store: DuckDBStore,
    atom_id: str | list[str] | None = None,
    facet: str | None = None,
    relation: str | None = None,
    source_db: str | None = None,
    source_accession: str | None = None,
    limit: int = 50,
) -> list[dict]:
    """Search raw V2 atom evidence with optional source/relation filters."""
    table_name = "semantic_atoms"
    atom_column = "atom_id"
    clauses: list[str] = []
    params: list[object] = []

    if atom_id:
        atom_ids = [atom_id] if isinstance(atom_id, str) else list(atom_id)
        placeholders = ",".join(["?"] * len(atom_ids))
        clauses.append(f"{atom_column} IN ({placeholders})")
        params.extend(atom_ids)

    if facet:
        clauses.append("facet = ?")
        params.append(facet)

    if relation:
        clauses.append("relation = ?")
        params.append(relation)

    if source_db:
        clauses.append("source_db = ?")
        params.append(source_db)

    if source_accession:
        clauses.append("source_accession = ?")
        params.append(source_accession)

    where = f"WHERE {' AND '.join(clauses)}" if clauses else ""
    rows = store.execute(
        f"""
        SELECT DISTINCT protein_id, {atom_column}, facet, relation,
               source_db, source_accession
        FROM {table_name}
        {where}
        ORDER BY protein_id, {atom_column}, source_db, source_accession
        LIMIT ?
        """,
        [*params, limit],
    )
    return [
        {
            "protein_id": row[0],
            "atom_id": row[1],
            "facet": row[2],
            "relation": row[3],
            "source_db": row[4],
            "source_accession": row[5] or None,
        }
        for row in rows
    ]


def search_by_atoms(
    store: DuckDBStore,
    has: list[str] | None = None,
    lacks: list[str] | None = None,
    limit: int = 50,
) -> list[str]:
    """Search proteins by V2 atom or composite predicate presence/absence."""
    has = has or []
    lacks = lacks or []
    if not has and not lacks:
        return []

    if _semantic_terms_has_rows(store):
        return _search_by_atoms_from_semantic_terms(store, has, lacks, limit)

    return _search_by_atoms_from_legacy_union(store, has, lacks, limit)


def list_composites() -> list:
    """List available composite predicate definitions."""
    from sharur.predicates_v2.composites import load_composites

    return load_composites()


def v2_review_queue(
    store: DuckDBStore,
    limit: int = 50,
    source: str | list[str] | None = None,
    min_proteins: int = 1,
    exclude_raw_system_profiles: bool = False,
    output_tsv: str | Path | None = None,
) -> list[dict]:
    """Get unmapped accession review queue."""
    from sharur.predicates_v2.review_queue import (
        format_review_queue_tsv,
        suggest_facet,
    )

    if min_proteins < 1:
        raise ValueError("min_proteins must be >= 1")

    where_clauses = [
        "s.relation = 'unresolved'",
        "s.source_accession IS NOT NULL",
        "s.source_db IS NOT NULL",
        "substr(s.source_accession, 1, 1) != '_'",
    ]
    params: list[object] = []

    if source:
        sources = [source] if isinstance(source, str) else list(source)
        if sources:
            placeholders = ",".join(["?"] * len(sources))
            where_clauses.append(f"s.source_db IN ({placeholders})")
            params.extend(sources)

    if exclude_raw_system_profiles:
        where_clauses.append("s.source_db NOT IN ('defensefinder', 'txsscan')")

    where_sql = " AND ".join(where_clauses)
    params.extend([min_proteins, limit])

    sql = f"""
        WITH unresolved AS (
            SELECT
                s.source_accession AS accession,
                s.source_db,
                s.protein_id,
                p.bin_id
            FROM semantic_atoms s
            LEFT JOIN proteins p ON p.protein_id = s.protein_id
            WHERE {where_sql}
        ),
        distinct_hits AS (
            SELECT DISTINCT accession, source_db, protein_id, bin_id
            FROM unresolved
        ),
        agg AS (
            SELECT
                accession,
                source_db,
                COUNT(DISTINCT protein_id) AS n_proteins,
                COUNT(DISTINCT bin_id) AS n_genomes
            FROM distinct_hits
            GROUP BY accession, source_db
            HAVING COUNT(DISTINCT protein_id) >= ?
        ),
        ranked_examples AS (
            SELECT
                accession,
                source_db,
                protein_id,
                row_number() OVER (
                    PARTITION BY accession, source_db
                    ORDER BY protein_id
                ) AS rn
            FROM (
                SELECT DISTINCT accession, source_db, protein_id
                FROM unresolved
            )
        ),
        examples AS (
            SELECT
                accession,
                source_db,
                string_agg(protein_id, ';' ORDER BY protein_id)
                    AS example_protein_ids
            FROM ranked_examples
            WHERE rn <= 5
            GROUP BY accession, source_db
        ),
        total AS (
            SELECT COUNT(*) AS total_proteins FROM proteins
        )
        SELECT
            agg.accession,
            agg.source_db,
            agg.n_proteins,
            agg.n_genomes,
            ROUND(
                agg.n_proteins * 100.0 / NULLIF(total.total_proteins, 0),
                2
            ) AS pct_proteome,
            COALESCE(examples.example_protein_ids, '')
                AS example_protein_ids
        FROM agg
        LEFT JOIN examples
            ON examples.accession = agg.accession
            AND examples.source_db = agg.source_db
        CROSS JOIN total
        ORDER BY
            agg.n_proteins * GREATEST(agg.n_genomes, 1) DESC,
            agg.n_proteins DESC,
            agg.source_db,
            agg.accession
        LIMIT ?
    """

    rows = store.execute(sql, params)
    queue = [
        {
            "accession": row[0],
            "source_db": row[1],
            "n_proteins": row[2],
            "n_genomes": row[3],
            "pct_proteome": row[4] or 0.0,
            "example_protein_ids": row[5],
            "suggested_facet": suggest_facet(row[0]),
        }
        for row in rows
    ]

    if output_tsv:
        Path(output_tsv).parent.mkdir(parents=True, exist_ok=True)
        Path(output_tsv).write_text(format_review_queue_tsv(queue))

    return queue


def run_shadow_diff(
    store: DuckDBStore,
    protein_ids: list[str] | None = None,
    output_report: str | None = None,
    output_jsonl: str | None = None,
) -> dict:
    """Run V1 vs V2 shadow comparison."""
    from sharur.predicates_v2.shadow import run_shadow_diff_on_store

    return run_shadow_diff_on_store(
        store,
        protein_ids=protein_ids,
        output_report=output_report,
        output_jsonl=output_jsonl,
    )


def _semantic_terms_has_rows(store: DuckDBStore) -> bool:
    """Return whether semantic_terms exists and has materialized rows."""
    try:
        return bool(store.execute("SELECT 1 FROM semantic_terms LIMIT 1"))
    except Exception:
        return False


def _table_exists(store: DuckDBStore, table_name: str) -> bool:
    """Return whether a table exists in the current DuckDB connection."""
    return table_name in _table_names(store)


def _table_names(store: DuckDBStore) -> set[str]:
    """Return existing DuckDB table names."""
    try:
        return {row[0] for row in store.execute("SHOW TABLES")}
    except Exception:
        return set()


def _fetch_protein_summary(
    store: DuckDBStore,
    protein_id: str,
    table_names: set[str],
) -> dict | None:
    """Fetch stable, compact protein metadata for explain()."""
    if "proteins" not in table_names:
        return None

    rows = store.execute(
        """
        SELECT protein_id, contig_id, bin_id, start, end_coord, strand,
               gene_index, sequence_length, gc_content
        FROM proteins
        WHERE protein_id = ?
        """,
        [protein_id],
    )
    if not rows:
        return None

    row = rows[0]
    return {
        "protein_id": row[0],
        "contig_id": row[1],
        "bin_id": row[2],
        "start": row[3],
        "end_coord": row[4],
        "strand": row[5],
        "gene_index": row[6],
        "sequence_length": row[7],
        "gc_content": row[8],
    }


def _fetch_semantic_terms(
    store: DuckDBStore,
    protein_id: str,
    table_names: set[str],
) -> list[dict]:
    """Fetch materialized atom, accession, and composite terms."""
    if "semantic_terms" not in table_names:
        return []

    rows = store.execute(
        """
        SELECT term_id, term_kind, facet, relation, source_db, source_accession
        FROM semantic_terms
        WHERE protein_id = ?
        ORDER BY term_kind, term_id, source_db, source_accession
        """,
        [protein_id],
    )
    return [
        {
            "term_id": row[0],
            "term_kind": row[1],
            "facet": row[2],
            "relation": row[3],
            "source_db": row[4] or None,
            "source_accession": row[5] or None,
        }
        for row in rows
    ]


def _fetch_validated_system_membership(
    store: DuckDBStore,
    protein_id: str,
    table_names: set[str],
) -> list[dict]:
    """Fetch normalized validated-system rows for a protein."""
    if "system_proteins" not in table_names:
        return []

    rows = store.execute(
        """
        SELECT system_id, system_source, position, profile_name, score
        FROM system_proteins
        WHERE protein_id = ?
        ORDER BY system_source, system_id, position
        """,
        [protein_id],
    )
    return [
        {
            "system_id": row[0],
            "system_source": row[1],
            "position": row[2],
            "profile_name": row[3],
            "score": row[4],
        }
        for row in rows
    ]


def _fetch_legacy_predicates(
    store: DuckDBStore,
    protein_id: str,
    table_names: set[str],
) -> list[str]:
    """Fetch the flat compatibility predicates when present."""
    if "protein_predicates" not in table_names:
        return []

    rows = store.execute(
        "SELECT predicates FROM protein_predicates WHERE protein_id = ?",
        [protein_id],
    )
    if not rows:
        return []
    return sorted(rows[0][0] or [])


def _unresolved_atoms(atoms: list[SemanticAtom]) -> list[dict]:
    """Return unresolved source accessions for explain()."""
    return [
        {
            "atom_id": atom.atom_id,
            "source_db": atom.source_db,
            "source_accession": atom.source_accession,
        }
        for atom in atoms
        if atom.relation == ClaimRelation.unresolved
    ]


def _conflicting_atoms(atoms: list[SemanticAtom]) -> list[dict]:
    """Summarize atoms with both positive and excluding evidence."""
    grouped: dict[str, list[SemanticAtom]] = defaultdict(list)
    for atom in atoms:
        grouped[atom.atom_id].append(atom)

    conflicts = []
    for atom_id, atom_group in sorted(grouped.items()):
        relations = {atom.relation.value for atom in atom_group}
        if "implies" not in relations or "excludes" not in relations:
            continue
        conflicts.append(
            {
                "atom_id": atom_id,
                "relations": sorted(relations),
                "sources": sorted(_format_atom_source(atom) for atom in atom_group),
            }
        )
    return conflicts


def _format_atom_source(atom: SemanticAtom) -> str:
    """Return a compact source key for an atom."""
    if atom.source_accession:
        return f"{atom.source_db}:{atom.source_accession}"
    return atom.source_db


def _search_by_atoms_from_semantic_terms(
    store: DuckDBStore,
    has: list[str],
    lacks: list[str],
    limit: int,
) -> list[str]:
    """Search using the materialized V2 term table."""
    active_filter = "(term_kind != 'atom' OR relation != 'excludes')"
    if len(has) == 1 and not lacks:
        rows = store.execute(
            "SELECT DISTINCT protein_id FROM semantic_terms "
            f"WHERE term_id = ? AND {active_filter} LIMIT ?",
            [has[0], limit],
        )
        return [row[0] for row in rows]

    with_active = (
        "active_terms AS ("
        "SELECT protein_id, term_id FROM semantic_terms "
        f"WHERE {active_filter}"
        ")"
    )
    params: list[object] = []

    if has:
        placeholders = ",".join(["?"] * len(has))
        base_sql = (
            "SELECT protein_id FROM active_terms "
            f"WHERE term_id IN ({placeholders}) "
            "GROUP BY protein_id "
            "HAVING COUNT(DISTINCT term_id) = ?"
        )
        params.extend(has)
        params.append(len(has))
    else:
        base_sql = "SELECT DISTINCT protein_id FROM active_terms"

    if lacks:
        placeholders = ",".join(["?"] * len(lacks))
        query = (
            f"WITH {with_active}, base AS ({base_sql}) "
            "SELECT protein_id FROM base "
            "WHERE NOT EXISTS ("
            "  SELECT 1 FROM active_terms t "
            "  WHERE t.protein_id = base.protein_id "
            f"  AND t.term_id IN ({placeholders})"
            ") LIMIT ?"
        )
        params.extend(lacks)
        params.append(limit)
    else:
        query = f"WITH {with_active} {base_sql} LIMIT ?"
        params.append(limit)

    rows = store.execute(query, params)
    return [row[0] for row in rows]


def _search_by_atoms_from_legacy_union(
    store: DuckDBStore,
    has: list[str],
    lacks: list[str],
    limit: int,
) -> list[str]:
    """Fallback search for databases that have not materialized semantic_terms."""
    params: list[object] = []
    unified = (
        "SELECT protein_id, atom_id FROM semantic_atoms "
        "UNION ALL "
        "SELECT protein_id, "
        "  CASE source_db "
        "    WHEN 'kofam' THEN 'kegg' "
        "    WHEN 'vogdb' THEN 'vog' "
        "    WHEN 'defensefinder_system' THEN 'defensefinder' "
        "    WHEN 'txsscan_system' THEN 'txsscan' "
        "    ELSE source_db "
        "  END || ':' || source_accession AS atom_id "
        "FROM semantic_atoms "
        "WHERE source_db IN ("
        "  'pfam', 'kegg', 'kofam', 'cazy', 'vog', 'vogdb', 'hyddb', "
        "  'hyddb_subgroup', 'defensefinder', 'defensefinder_system', "
        "  'txsscan', 'txsscan_system', 'cant_hyd'"
        ") "
        "AND source_accession IS NOT NULL "
        "AND LEFT(source_accession, 1) != '_' "
        "UNION ALL "
        "SELECT protein_id, UNNEST(composite_predicates) as atom_id "
        "FROM semantic_state WHERE LEN(composite_predicates) > 0"
    )

    if has:
        placeholders = ",".join(["?"] * len(has))
        base = (
            f"SELECT protein_id FROM ({unified}) u "
            f"WHERE atom_id IN ({placeholders}) "
            "GROUP BY protein_id "
            "HAVING COUNT(DISTINCT atom_id) = ?"
        )
        params.extend(has)
        params.append(len(has))
    else:
        base = f"SELECT DISTINCT protein_id FROM ({unified}) u"

    if lacks:
        placeholders = ",".join(["?"] * len(lacks))
        query = (
            f"WITH base AS ({base}) "
            "SELECT protein_id FROM base "
            "WHERE protein_id NOT IN ("
            f"  SELECT DISTINCT protein_id FROM ({unified}) u2 "
            f"  WHERE atom_id IN ({placeholders})"
            ") LIMIT ?"
        )
        params.extend(lacks)
        params.append(limit)
    else:
        query = f"{base} LIMIT ?"
        params.append(limit)

    rows = store.execute(query, params)
    return [row[0] for row in rows]


__all__ = [
    "explain",
    "generate_v2",
    "get_atoms",
    "get_semantic_state",
    "list_composites",
    "run_shadow_diff",
    "search_atoms",
    "search_by_atoms",
    "search_by_facet",
    "v2_review_queue",
]
