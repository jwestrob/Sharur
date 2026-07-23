"""First-class biological case inspection and genomic-context comparison."""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal

from sharur.assembly_evidence import (
    AssemblyEvidenceStore,
    discover_assembly_evidence,
)
from sharur.core.case_models import (
    DEFAULT_CONTEXT_ORFS,
    AnnotationEvidence,
    CaseEntityType,
    CaseRecord,
    ContextComparison,
    ContextFeature,
    ContextGroupSummary,
    ContextStatistic,
    ContextWindow,
    EntityReference,
    NamedCallEvidence,
    ProteinContextRecord,
    SyntenyMembershipEvidence,
)
from sharur.operators.semantics import get_active_predicates
from sharur.synteny import (
    SyntenyDatasetMismatchError,
    SyntenyStore,
    discover_synteny_sidecar,
)


if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence

    from sharur.core.claim_compiler import FindingDraft
    from sharur.storage.duckdb_store import DuckDBStore


_TAXONOMY_PREFIX = {
    "domain": "d__",
    "phylum": "p__",
    "class": "c__",
    "order": "o__",
    "family": "f__",
    "genus": "g__",
    "species": "s__",
}


@dataclass(frozen=True)
class _ResolvedEntity:
    entity_id: str
    entity_type: CaseEntityType
    source_table: str
    row: dict[str, Any]
    name: str | None = None
    subtype: str | None = None


@dataclass(frozen=True)
class _Anchor:
    entity_id: str
    entity_name: str | None
    bin_id: str
    contig_id: str
    min_gene_index: int
    max_gene_index: int
    orientation: Literal["+", "-", "mixed", "unknown"]
    component_count: int
    taxonomy: str | None
    source_table: str


@dataclass(frozen=True, slots=True)
class _ContextAnnotation:
    """Minimal annotation projection used by cohort comparisons."""

    source: str
    accession: str
    name: str | None
    description: str | None


def _quote_identifier(value: str) -> str:
    return '"' + value.replace('"', '""') + '"'


def _table_catalog(store: DuckDBStore) -> dict[str, set[str]]:
    """Return stored table schemas, excluding convenience/projection views."""
    cached = getattr(store, "_sharur_table_catalog", None)
    if cached is not None:
        return cached
    rows = store.execute(
        """
        SELECT columns.table_name, columns.column_name
        FROM information_schema.columns AS columns
        JOIN information_schema.tables AS tables
          ON tables.table_catalog = columns.table_catalog
         AND tables.table_schema = columns.table_schema
         AND tables.table_name = columns.table_name
        WHERE columns.table_schema = 'main'
          AND tables.table_type = 'BASE TABLE'
        ORDER BY columns.table_name, columns.ordinal_position
        """
    )
    catalog: dict[str, set[str]] = {}
    for table_name, column_name in rows:
        catalog.setdefault(str(table_name), set()).add(str(column_name))
    if getattr(store, "read_only", False):
        store._sharur_table_catalog = catalog
    return catalog


def _fetch_dicts(
    store: DuckDBStore,
    query: str,
    params: Sequence[Any] | None = None,
) -> list[dict[str, Any]]:
    cursor = store.conn.execute(query, list(params)) if params else store.conn.execute(query)
    columns = [str(item[0]) for item in cursor.description]
    return [dict(zip(columns, row, strict=True)) for row in cursor.fetchall()]


def _fetch_one_dict(
    store: DuckDBStore,
    query: str,
    params: Sequence[Any] | None = None,
) -> dict[str, Any] | None:
    rows = _fetch_dicts(store, query, params)
    return rows[0] if rows else None


def _system_tables(catalog: Mapping[str, set[str]]) -> list[str]:
    return sorted(
        table_name
        for table_name, columns in catalog.items()
        if {"system_id", "system_type"} <= columns and table_name != "system_proteins"
    )


def _locus_tables(catalog: Mapping[str, set[str]]) -> list[str]:
    return sorted(
        table_name
        for table_name, columns in catalog.items()
        if {"locus_id", "locus_type"} <= columns
    )


def _structured_projection_sources(
    store: DuckDBStore,
    catalog: Mapping[str, set[str]],
) -> list[str]:
    """Discover annotation sources used for structured system projections."""
    cached = getattr(store, "_sharur_projection_sources", None)
    if cached is not None:
        return cached
    if "system_proteins" not in catalog:
        if getattr(store, "read_only", False):
            store._sharur_projection_sources = []
        return []
    sources = [
        str(row[0])
        for row in store.execute(
            """
            SELECT DISTINCT system_source
            FROM system_proteins
            WHERE system_source IS NOT NULL AND system_source != ''
            ORDER BY system_source
            """
        )
    ]
    if getattr(store, "read_only", False):
        store._sharur_projection_sources = sources
    return sources


def _resolve_windows(
    window: int,
    upstream_orfs: int | None,
    downstream_orfs: int | None,
) -> tuple[int, int, int]:
    if window < 0:
        raise ValueError("window must be non-negative")
    upstream = window if upstream_orfs is None else upstream_orfs
    downstream = window if downstream_orfs is None else downstream_orfs
    if upstream < 0 or downstream < 0:
        raise ValueError("upstream_orfs and downstream_orfs must be non-negative")
    return window, upstream, downstream


def _resolve_entity(
    store: DuckDBStore,
    entity_id: str,
    *,
    entity_type: CaseEntityType | str | None = None,
    bin_id: str | None = None,
    source_table: str | None = None,
) -> _ResolvedEntity:
    catalog = _table_catalog(store)
    requested_type = CaseEntityType(entity_type) if entity_type is not None else None
    candidates: list[_ResolvedEntity] = []

    if requested_type in {None, CaseEntityType.protein} and "proteins" in catalog:
        row = _fetch_one_dict(
            store,
            "SELECT * FROM proteins WHERE protein_id = ?",
            [entity_id],
        )
        if row is not None:
            candidates.append(
                _ResolvedEntity(
                    entity_id,
                    CaseEntityType.protein,
                    "proteins",
                    row,
                )
            )

    if requested_type in {None, CaseEntityType.system}:
        for table_name in _system_tables(catalog):
            if source_table is not None and table_name != source_table:
                continue
            row = _fetch_one_dict(
                store,
                f"SELECT * FROM {_quote_identifier(table_name)} WHERE system_id = ?",
                [entity_id],
            )
            if row is not None:
                candidates.append(
                    _ResolvedEntity(
                        entity_id,
                        CaseEntityType.system,
                        table_name,
                        row,
                        name=str(row.get("system_type") or "") or None,
                        subtype=str(row.get("system_subtype") or "") or None,
                    )
                )

    if requested_type in {None, CaseEntityType.locus}:
        for table_name in _locus_tables(catalog):
            if source_table is not None and table_name != source_table:
                continue
            row = _fetch_one_dict(
                store,
                f"SELECT * FROM {_quote_identifier(table_name)} WHERE locus_id = ?",
                [entity_id],
            )
            if row is not None:
                candidates.append(
                    _ResolvedEntity(
                        entity_id,
                        CaseEntityType.locus,
                        table_name,
                        row,
                        name=str(row.get("locus_type") or "") or None,
                    )
                )

    if requested_type in {None, CaseEntityType.bin} and "bins" in catalog:
        row = _fetch_one_dict(
            store,
            "SELECT * FROM bins WHERE bin_id = ?",
            [entity_id],
        )
        if row is not None:
            candidates.append(_ResolvedEntity(entity_id, CaseEntityType.bin, "bins", row))

    if requested_type in {None, CaseEntityType.contig} and "proteins" in catalog:
        rows = _fetch_dicts(
            store,
            """
            SELECT DISTINCT contig_id, bin_id
            FROM proteins
            WHERE contig_id = ?
              AND (? IS NULL OR bin_id = ?)
            ORDER BY bin_id
            """,
            [entity_id, bin_id, bin_id],
        )
        if len(rows) > 1 and bin_id is None:
            bins = ", ".join(str(row["bin_id"]) for row in rows)
            raise ValueError(
                f"Contig ID {entity_id!r} occurs in multiple bins ({bins}); pass bin_id explicitly"
            )
        if rows:
            row = rows[0]
            candidates.append(
                _ResolvedEntity(
                    entity_id,
                    CaseEntityType.contig,
                    "contigs",
                    row,
                )
            )

    if not candidates:
        qualifier = f" as {requested_type.value}" if requested_type else ""
        raise KeyError(f"Entity {entity_id!r} was not found{qualifier}")
    if len(candidates) > 1:
        descriptions = ", ".join(
            f"{candidate.entity_type.value}:{candidate.source_table}" for candidate in candidates
        )
        raise ValueError(
            f"Entity ID {entity_id!r} is ambiguous across {descriptions}; "
            "pass entity_type and, for caller resources, source_table"
        )
    return candidates[0]


def _parse_serialized_protein_ids(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return [str(item) for item in value if str(item).strip()]
    return [token for token in re.split(r"[,;\s]+", str(value).strip()) if token]


def _component_ids(
    store: DuckDBStore,
    entity: _ResolvedEntity,
    catalog: Mapping[str, set[str]],
) -> list[str]:
    if entity.entity_type == CaseEntityType.protein:
        return [entity.entity_id]
    if entity.entity_type == CaseEntityType.system:
        call_rows = _call_type_map(store, [entity.entity_id], catalog).get(
            entity.entity_id,
            [],
        )
        unambiguous_membership = (
            len(call_rows) == 1 and call_rows[0].get("source_table") == entity.source_table
        )
        if "system_proteins" in catalog and unambiguous_membership:
            rows = store.execute(
                """
                SELECT protein_id
                FROM system_proteins
                WHERE system_id = ?
                ORDER BY position NULLS LAST, protein_id
                """,
                [entity.entity_id],
            )
            if rows:
                return [str(row[0]) for row in rows]
        return _parse_serialized_protein_ids(entity.row.get("protein_ids"))
    if entity.entity_type == CaseEntityType.locus and "locus_proteins" in catalog:
        rows = store.execute(
            """
            SELECT protein_id
            FROM locus_proteins
            WHERE locus_id = ?
            ORDER BY position, protein_id
            """,
            [entity.entity_id],
        )
        return [str(row[0]) for row in rows]
    return []


def _load_proteins(
    store: DuckDBStore,
    protein_ids: Sequence[str],
) -> list[dict[str, Any]]:
    if not protein_ids:
        return []
    placeholders = ", ".join(["?"] * len(protein_ids))
    return _fetch_dicts(
        store,
        f"""
        SELECT protein_id, contig_id, bin_id, start, end_coord, strand,
               gene_index, sequence, sequence_length, gc_content
        FROM proteins
        WHERE protein_id IN ({placeholders})
        ORDER BY bin_id, contig_id, gene_index NULLS LAST, start
        """,
        protein_ids,
    )


def _normalize_strand(value: Any) -> Literal["+", "-"] | None:
    """Normalize the strand encodings present in legacy and current datasets."""
    normalized = str(value).strip().lower()
    if normalized in {"+", "+1", "1", "forward", "f"}:
        return "+"
    if normalized in {"-", "-1", "reverse", "r"}:
        return "-"
    return None


def _orientation(rows: Sequence[Mapping[str, Any]]) -> Literal["+", "-", "mixed", "unknown"]:
    strands = {
        strand for row in rows if (strand := _normalize_strand(row.get("strand"))) is not None
    }
    if len(strands) == 1:
        return next(iter(strands))  # type: ignore[return-value]
    if len(strands) > 1:
        return "mixed"
    return "unknown"


def _window_bounds(
    anchor_min: int,
    anchor_max: int,
    orientation: str,
    upstream: int,
    downstream: int,
) -> tuple[int, int]:
    if orientation == "-":
        return anchor_min - downstream, anchor_max + upstream
    return anchor_min - upstream, anchor_max + downstream


def _relative_orf(
    gene_index: int,
    anchor_min: int,
    anchor_max: int,
    orientation: str,
) -> int:
    if anchor_min <= gene_index <= anchor_max:
        return 0
    if orientation == "-":
        if gene_index > anchor_max:
            return -(gene_index - anchor_max)
        return anchor_min - gene_index
    if gene_index < anchor_min:
        return gene_index - anchor_min
    return gene_index - anchor_max


def _context_completeness(
    relative_orfs: Sequence[int],
    upstream_orfs: int,
    downstream_orfs: int,
) -> tuple[int, int, bool, bool]:
    """Count distinct flanking indices and require every requested position."""
    positions = set(relative_orfs)
    realized_upstream = len({position for position in positions if position < 0})
    realized_downstream = len({position for position in positions if position > 0})
    complete_upstream = all(-offset in positions for offset in range(1, upstream_orfs + 1))
    complete_downstream = all(offset in positions for offset in range(1, downstream_orfs + 1))
    return (
        realized_upstream,
        realized_downstream,
        complete_upstream,
        complete_downstream,
    )


def _load_annotation_map(
    store: DuckDBStore,
    protein_ids: Sequence[str],
    catalog: Mapping[str, set[str]],
    *,
    batch_size: int = 5_000,
) -> dict[str, list[AnnotationEvidence]]:
    result: dict[str, list[AnnotationEvidence]] = {}
    excluded_sources = _structured_projection_sources(store, catalog)
    for start in range(0, len(protein_ids), batch_size):
        batch = list(protein_ids[start : start + batch_size])
        if not batch:
            continue
        placeholders = ", ".join(["?"] * len(batch))
        exclusion_sql = ""
        if excluded_sources:
            exclusion_placeholders = ", ".join(["?"] * len(excluded_sources))
            exclusion_sql = f"AND source NOT IN ({exclusion_placeholders})"
        rows = _fetch_dicts(
            store,
            f"""
            SELECT protein_id, source, accession, name, description, evalue,
                   score, start_aa, end_aa
            FROM annotations
            WHERE protein_id IN ({placeholders})
              {exclusion_sql}
            ORDER BY protein_id, source, evalue NULLS LAST, score DESC NULLS LAST
            """,
            [*batch, *excluded_sources],
        )
        for row in rows:
            evidence = AnnotationEvidence(**row)
            result.setdefault(evidence.protein_id, []).append(evidence)
    return result


def _load_feature_annotation_map(
    store: DuckDBStore,
    protein_ids: Sequence[str],
    features: Sequence[ContextFeature],
    catalog: Mapping[str, set[str]],
    *,
    batch_size: int = 50_000,
) -> dict[str, list[_ContextAnnotation]]:
    """Load only annotation rows capable of satisfying a requested feature."""
    annotation_features = [
        feature for feature in features if feature.kind in {"annotation", "annotation_name"}
    ]
    if not annotation_features or not protein_ids:
        return {}

    excluded_sources = _structured_projection_sources(store, catalog)
    conditions: list[str] = []
    condition_params: list[Any] = []
    for feature in annotation_features:
        if feature.kind == "annotation":
            clause = "accession = ?"
            condition_params.append(feature.value)
            if feature.source is not None:
                clause += " AND LOWER(source) = LOWER(?)"
                condition_params.append(feature.source)
            conditions.append(f"({clause})")
        else:
            conditions.append(
                "(LOWER(COALESCE(name, '')) LIKE ? OR LOWER(COALESCE(description, '')) LIKE ?)"
            )
            needle = f"%{(feature.value or '').lower()}%"
            condition_params.extend([needle, needle])

    result: dict[str, list[_ContextAnnotation]] = {}
    for start in range(0, len(protein_ids), batch_size):
        batch = list(protein_ids[start : start + batch_size])
        placeholders = ", ".join(["?"] * len(batch))
        exclusion_sql = ""
        if excluded_sources:
            exclusion_placeholders = ", ".join(["?"] * len(excluded_sources))
            exclusion_sql = f"AND source NOT IN ({exclusion_placeholders})"
        rows = _fetch_dicts(
            store,
            f"""
            SELECT protein_id, source, accession, name, description
            FROM annotations
            WHERE protein_id IN ({placeholders})
              {exclusion_sql}
              AND ({" OR ".join(conditions)})
            ORDER BY protein_id, source, accession
            """,
            [*batch, *excluded_sources, *condition_params],
        )
        for row in rows:
            result.setdefault(str(row["protein_id"]), []).append(
                _ContextAnnotation(
                    source=str(row["source"]),
                    accession=str(row["accession"]),
                    name=str(row["name"]) if row.get("name") is not None else None,
                    description=(
                        str(row["description"]) if row.get("description") is not None else None
                    ),
                )
            )
    return result


def _call_type_map(
    store: DuckDBStore,
    call_ids: Sequence[str],
    catalog: Mapping[str, set[str]],
) -> dict[str, list[dict[str, Any]]]:
    result: dict[str, list[dict[str, Any]]] = {}
    if not call_ids:
        return result
    unique_ids = sorted(set(call_ids))
    for table_name in _system_tables(catalog):
        for start in range(0, len(unique_ids), 5_000):
            batch = unique_ids[start : start + 5_000]
            placeholders = ", ".join(["?"] * len(batch))
            columns = catalog[table_name]
            subtype_select = (
                "system_subtype" if "system_subtype" in columns else "NULL AS system_subtype"
            )
            rows = _fetch_dicts(
                store,
                f"""
                SELECT system_id, system_type, {subtype_select}
                FROM {_quote_identifier(table_name)}
                WHERE system_id IN ({placeholders})
                """,
                batch,
            )
            for row in rows:
                row["source_table"] = table_name
                result.setdefault(str(row["system_id"]), []).append(row)
    return result


def _load_named_call_map(
    store: DuckDBStore,
    protein_ids: Sequence[str],
    catalog: Mapping[str, set[str]],
    *,
    batch_size: int = 5_000,
) -> dict[str, list[NamedCallEvidence]]:
    result: dict[str, list[NamedCallEvidence]] = {}
    if "system_proteins" not in catalog or not protein_ids:
        return result

    membership_rows: list[dict[str, Any]] = []
    for start in range(0, len(protein_ids), batch_size):
        batch = list(protein_ids[start : start + batch_size])
        placeholders = ", ".join(["?"] * len(batch))
        membership_rows.extend(
            _fetch_dicts(
                store,
                f"""
                SELECT system_id, protein_id, system_source, position,
                       profile_name, score
                FROM system_proteins
                WHERE protein_id IN ({placeholders})
                ORDER BY protein_id, system_id
                """,
                batch,
            )
        )
    type_map = _call_type_map(
        store,
        [str(row["system_id"]) for row in membership_rows],
        catalog,
    )
    for row in membership_rows:
        call_rows = type_map.get(str(row["system_id"]), [])
        if len(call_rows) != 1:
            # A membership source cannot be mapped generically onto two
            # structured tables. Failing closed avoids assigning the wrong
            # caller-emitted name when system IDs collide across resources.
            continue
        for call_row in call_rows:
            evidence = NamedCallEvidence(
                call_id=str(row["system_id"]),
                call_type=str(call_row["system_type"]),
                call_subtype=(
                    str(call_row["system_subtype"]) if call_row.get("system_subtype") else None
                ),
                source_table=str(call_row["source_table"]),
                system_source=(str(row["system_source"]) if row.get("system_source") else None),
                protein_id=str(row["protein_id"]),
                profile_name=(str(row["profile_name"]) if row.get("profile_name") else None),
                position=int(row["position"]) if row.get("position") is not None else None,
                score=float(row["score"]) if row.get("score") is not None else None,
            )
            result.setdefault(evidence.protein_id or "", []).append(evidence)
    return result


def _load_locus_memberships(
    store: DuckDBStore,
    protein_ids: Sequence[str],
    catalog: Mapping[str, set[str]],
    *,
    batch_size: int = 5_000,
) -> dict[str, set[str]]:
    result: dict[str, set[str]] = {}
    if "locus_proteins" not in catalog or "loci" not in catalog:
        return result
    for start in range(0, len(protein_ids), batch_size):
        batch = list(protein_ids[start : start + batch_size])
        placeholders = ", ".join(["?"] * len(batch))
        rows = store.execute(
            f"""
            SELECT lp.protein_id, l.locus_type
            FROM locus_proteins lp
            JOIN loci l ON l.locus_id = lp.locus_id
            WHERE lp.protein_id IN ({placeholders})
            """,
            batch,
        )
        for protein_id, locus_type in rows:
            result.setdefault(str(protein_id), set()).add(str(locus_type))
    return result


def structured_projection_sources(store: DuckDBStore) -> list[str]:
    """Return raw-annotation sources that project structured caller output.

    These sources belong on the caller-named side of the provenance boundary,
    so navigation and packet operators can keep them separate from observed
    per-domain annotations.
    """
    catalog = _table_catalog(store)
    return _structured_projection_sources(store, catalog)


def load_context_evidence(
    store: DuckDBStore,
    protein_ids: Sequence[str],
) -> dict[str, Any]:
    """Batch observed annotations and exact named memberships by protein.

    System names are returned only when ``system_proteins`` membership maps
    unambiguously to one live structured caller table. Locus memberships come
    directly from the normalized ``loci``/``locus_proteins`` resources.
    """
    unique_ids = sorted(set(protein_ids))
    result: dict[str, Any] = {
        "projection_sources": [],
        "observed_annotations": {},
        "named_calls": {},
        "loci": {},
    }
    if not unique_ids:
        return result

    catalog = _table_catalog(store)
    result["projection_sources"] = _structured_projection_sources(store, catalog)
    annotation_map = _load_annotation_map(store, unique_ids, catalog)
    result["observed_annotations"] = {
        protein_id: [
            annotation.model_dump(mode="json")
            for annotation in annotations
        ]
        for protein_id, annotations in annotation_map.items()
    }
    named_calls = _load_named_call_map(store, unique_ids, catalog)
    result["named_calls"] = {
        protein_id: [
            call.model_dump(mode="json")
            for call in calls
        ]
        for protein_id, calls in named_calls.items()
    }

    if "locus_proteins" not in catalog or "loci" not in catalog:
        return result
    loci_by_protein: dict[str, list[dict[str, Any]]] = {}
    for start in range(0, len(unique_ids), 5_000):
        batch = unique_ids[start : start + 5_000]
        placeholders = ", ".join(["?"] * len(batch))
        rows = store.execute(
            f"""
            SELECT
                lp.protein_id,
                l.locus_id,
                l.locus_type,
                lp.position,
                l.confidence
            FROM locus_proteins AS lp
            JOIN loci AS l ON l.locus_id = lp.locus_id
            WHERE lp.protein_id IN ({placeholders})
            ORDER BY lp.protein_id, l.locus_type, l.locus_id, lp.position
            """,
            batch,
        )
        for protein_id, locus_id, locus_type, position, confidence in rows:
            loci_by_protein.setdefault(str(protein_id), []).append(
                {
                    "locus_id": str(locus_id),
                    "locus_type": str(locus_type),
                    "source_table": "loci",
                    "position": int(position),
                    "confidence": (
                        float(confidence) if confidence is not None else None
                    ),
                    "evidence_level": "caller_named",
                }
            )
    result["loci"] = loci_by_protein
    return result


def load_structured_context_evidence(
    store: DuckDBStore,
    protein_ids: Sequence[str],
) -> dict[str, Any]:
    """Compatibility alias for the complete packet-oriented evidence loader."""
    return load_context_evidence(store, protein_ids)


def _bin_record(store: DuckDBStore, bin_id: str | None) -> dict[str, Any] | None:
    if bin_id is None:
        return None
    return _fetch_one_dict(store, "SELECT * FROM bins WHERE bin_id = ?", [bin_id])


def _contig_record(
    store: DuckDBStore,
    contig_id: str | None,
    bin_id: str | None,
) -> dict[str, Any] | None:
    if contig_id is None:
        return None
    try:
        rows = _fetch_dicts(
            store,
            """
            SELECT *
            FROM contigs
            WHERE contig_id = ?
              AND (? IS NULL OR bin_id = ?)
            """,
            [contig_id, bin_id, bin_id],
        )
    except Exception:
        return None
    return rows[0] if len(rows) == 1 else None


def inspect_case(
    store: DuckDBStore,
    entity_id: str,
    *,
    entity_type: CaseEntityType | str | None = None,
    bin_id: str | None = None,
    source_table: str | None = None,
    window: int = DEFAULT_CONTEXT_ORFS,
    upstream_orfs: int | None = None,
    downstream_orfs: int | None = None,
    include_sequences: bool = False,
    assembly_evidence_path: str | Path | None = None,
    synteny_path: str | Path | None = None,
    allow_stale_synteny: bool = False,
) -> BiologicalCase:
    """Resolve an entity into a typed, provenance-separated biological case."""
    default_window, upstream, downstream = _resolve_windows(
        window,
        upstream_orfs,
        downstream_orfs,
    )
    catalog = _table_catalog(store)
    entity = _resolve_entity(
        store,
        entity_id,
        entity_type=entity_type,
        bin_id=bin_id,
        source_table=source_table,
    )
    limitations: list[str] = []
    component_ids = _component_ids(store, entity, catalog)
    component_rows = _load_proteins(store, component_ids)
    missing_component_ids = sorted(
        set(component_ids) - {str(row["protein_id"]) for row in component_rows}
    )
    if missing_component_ids:
        preview = ", ".join(missing_component_ids[:5])
        limitations.append(
            f"{len(missing_component_ids)} caller component IDs do not resolve "
            f"to proteins, including {preview}."
        )

    resolved_bin_id = bin_id
    resolved_contig_id: str | None = None
    if component_rows:
        replicons = {(str(row.get("bin_id")), str(row.get("contig_id"))) for row in component_rows}
        if len(replicons) > 1:
            limitations.append(
                "Anchor components span multiple (bin_id, contig_id) replicons; "
                "neighborhood construction failed closed."
            )
        else:
            resolved_bin_id, resolved_contig_id = next(iter(replicons))
    elif entity.entity_type == CaseEntityType.system:
        resolved_bin_id = (
            str(entity.row.get("genome_id"))
            if entity.row.get("genome_id") is not None
            else resolved_bin_id
        )
        resolved_contig_id = (
            str(entity.row.get("contig_id")) if entity.row.get("contig_id") is not None else None
        )
    elif entity.entity_type == CaseEntityType.locus:
        resolved_contig_id = (
            str(entity.row.get("contig_id")) if entity.row.get("contig_id") is not None else None
        )
        if resolved_contig_id and resolved_bin_id is None:
            bins = store.execute(
                "SELECT DISTINCT bin_id FROM proteins WHERE contig_id = ?",
                [resolved_contig_id],
            )
            if len(bins) == 1:
                resolved_bin_id = str(bins[0][0])
            elif len(bins) > 1:
                limitations.append(
                    "Locus contig label maps to multiple bins; pass bin_id explicitly."
                )
    elif entity.entity_type == CaseEntityType.contig:
        resolved_bin_id = str(entity.row.get("bin_id") or bin_id)
        resolved_contig_id = entity.entity_id
    elif entity.entity_type == CaseEntityType.bin:
        resolved_bin_id = entity.entity_id

    context_window: ContextWindow | None = None
    context_rows: list[dict[str, Any]] = []
    orientation: Literal["+", "-", "mixed", "unknown"] = _orientation(component_rows)
    valid_component_rows = [row for row in component_rows if row.get("gene_index") is not None]
    if (
        resolved_bin_id is not None
        and resolved_contig_id is not None
        and valid_component_rows
        and len({(row.get("bin_id"), row.get("contig_id")) for row in valid_component_rows}) == 1
    ):
        anchor_min = min(int(row["gene_index"]) for row in valid_component_rows)
        anchor_max = max(int(row["gene_index"]) for row in valid_component_rows)
        requested_min, requested_max = _window_bounds(
            anchor_min,
            anchor_max,
            orientation,
            upstream,
            downstream,
        )
        context_rows = _fetch_dicts(
            store,
            """
            SELECT protein_id, contig_id, bin_id, start, end_coord, strand,
                   gene_index, sequence, sequence_length, gc_content
            FROM proteins
            WHERE bin_id = ? AND contig_id = ?
              AND gene_index BETWEEN ? AND ?
            ORDER BY gene_index, start
            """,
            [
                resolved_bin_id,
                resolved_contig_id,
                requested_min,
                requested_max,
            ],
        )
        relative_values = [
            _relative_orf(
                int(row["gene_index"]),
                anchor_min,
                anchor_max,
                orientation,
            )
            for row in context_rows
            if row.get("gene_index") is not None
        ]
        (
            realized_upstream,
            realized_downstream,
            complete_upstream,
            complete_downstream,
        ) = _context_completeness(relative_values, upstream, downstream)
        context_window = ContextWindow(
            default_orfs=default_window,
            upstream_orfs=upstream,
            downstream_orfs=downstream,
            orientation=orientation,
            orientation_basis=(
                "co-oriented anchor component strands"
                if orientation in {"+", "-"}
                else "coordinate-order fallback for mixed/unknown anchor orientation"
            ),
            requested_min_gene_index=requested_min,
            requested_max_gene_index=requested_max,
            realized_upstream_orfs=realized_upstream,
            realized_downstream_orfs=realized_downstream,
            complete_upstream=complete_upstream,
            complete_downstream=complete_downstream,
        )
        if orientation in {"mixed", "unknown"}:
            limitations.append(
                "Anchor orientation is mixed or unknown; upstream/downstream use "
                "ascending coordinate order."
            )
        if not context_window.complete_upstream or not context_window.complete_downstream:
            limitations.append(
                "Requested ORF context is censored by a contig edge or missing gene indices."
            )
    elif component_ids and not valid_component_rows:
        limitations.append(
            "Anchor proteins lack gene_index values; ORF neighborhood is unavailable."
        )

    context_ids = [str(row["protein_id"]) for row in context_rows]
    annotation_map = _load_annotation_map(store, context_ids, catalog)
    named_call_map = _load_named_call_map(store, context_ids, catalog)
    predicate_map = get_active_predicates(store, context_ids) if context_ids else {}
    component_set = set(component_ids)
    if valid_component_rows:
        anchor_min = min(int(row["gene_index"]) for row in valid_component_rows)
        anchor_max = max(int(row["gene_index"]) for row in valid_component_rows)
    else:
        anchor_min = anchor_max = 0

    proteins: list[ProteinContextRecord] = []
    for row in context_rows:
        relative = (
            _relative_orf(
                int(row["gene_index"]),
                anchor_min,
                anchor_max,
                orientation,
            )
            if row.get("gene_index") is not None
            else None
        )
        role: Literal["upstream", "anchor", "downstream", "context"]
        if relative is None:
            role = "context"
        elif relative < 0:
            role = "upstream"
        elif relative > 0:
            role = "downstream"
        else:
            role = "anchor"
        proteins.append(
            ProteinContextRecord(
                protein_id=str(row["protein_id"]),
                bin_id=str(row["bin_id"]) if row.get("bin_id") is not None else None,
                contig_id=str(row["contig_id"]),
                start=int(row["start"]),
                end=int(row["end_coord"]),
                strand=_normalize_strand(row.get("strand")) or "unknown",
                gene_index=(int(row["gene_index"]) if row.get("gene_index") is not None else None),
                sequence_length=(
                    int(row["sequence_length"]) if row.get("sequence_length") is not None else None
                ),
                sequence=(
                    str(row["sequence"])
                    if include_sequences and row.get("sequence") is not None
                    else None
                ),
                relative_orf=relative,
                region_role=role,
                is_component=str(row["protein_id"]) in component_set,
                predicates=predicate_map.get(str(row["protein_id"]), []),
                annotations=annotation_map.get(str(row["protein_id"]), []),
                named_calls=named_call_map.get(str(row["protein_id"]), []),
            )
        )

    named_calls = [
        call
        for protein_calls in named_call_map.values()
        for call in protein_calls
        if call.call_id == entity.entity_id
    ]
    if entity.entity_type == CaseEntityType.system and not named_calls:
        named_calls = [
            NamedCallEvidence(
                call_id=entity.entity_id,
                call_type=entity.name or "unknown",
                call_subtype=entity.subtype,
                source_table=entity.source_table,
            )
        ]
    if entity.entity_type == CaseEntityType.locus:
        named_calls.append(
            NamedCallEvidence(
                call_id=entity.entity_id,
                call_type=entity.name or "unknown",
                source_table=entity.source_table,
            )
        )

    assembly_record = None
    assembly_state: Literal["available", "unavailable", "failed"] = "unavailable"
    db_path = getattr(store, "db_path", None)
    discovered_sidecar = None
    if db_path is not None:
        discovered_sidecar = discover_assembly_evidence(
            db_path,
            explicit_path=assembly_evidence_path,
        )
    if (
        discovered_sidecar is not None
        and resolved_bin_id is not None
        and resolved_contig_id is not None
    ):
        try:
            with AssemblyEvidenceStore(discovered_sidecar, read_only=True) as evidence_store:
                assembly_record = evidence_store.get(
                    resolved_bin_id,
                    resolved_contig_id,
                )
            assembly_state = "available" if assembly_record is not None else "unavailable"
            if assembly_record is None:
                limitations.append(
                    "Assembly-evidence sidecar exists but has no row for this replicon."
                )
        except Exception as exc:
            assembly_state = "failed"
            limitations.append(
                f"Assembly-evidence sidecar could not be read: {type(exc).__name__}: {exc}"
            )
    elif assembly_evidence_path is not None and discovered_sidecar is None:
        assembly_state = "failed"
        limitations.append(
            "Explicit assembly-evidence sidecar does not exist: "
            f"{Path(assembly_evidence_path).expanduser().resolve()}"
        )
    elif resolved_contig_id is not None:
        limitations.append(
            "Optional assembly/host-assignment evidence is not available for this dataset."
        )

    synteny_memberships: list[SyntenyMembershipEvidence] = []
    synteny_state: Literal["available", "unavailable", "stale", "failed"] = (
        "unavailable"
    )
    discovered_synteny = None
    if db_path is not None:
        discovered_synteny = discover_synteny_sidecar(
            db_path,
            explicit_path=synteny_path,
        )
        if discovered_synteny is not None:
            try:
                if component_ids:
                    with SyntenyStore(
                        discovered_synteny,
                        core_db_path=db_path,
                        allow_stale=allow_stale_synteny,
                    ) as synteny_store:
                        membership_rows, membership_total, _run_id = (
                            synteny_store.protein_memberships(
                                component_ids,
                                limit=25,
                            )
                        )
                    synteny_memberships = [
                        SyntenyMembershipEvidence(
                            run_id=str(row["run_id"]),
                            protein_id=str(row["protein_id"]),
                            member_role=str(row["member_role"]),
                            cluster_key=str(row["cluster_key"]),
                            source_cluster_id=(
                                int(row["source_cluster_id"])
                                if row.get("source_cluster_id") is not None
                                else None
                            ),
                            cluster_kind=str(row["cluster_kind"]),
                            block_count=int(row["block_count"]),
                            genome_support=int(row["genome_support"]),
                            locus_key=str(row["locus_key"]),
                            locus_genome_id=str(row["genome_id"]),
                            locus_contig_id=str(row["contig_id"]),
                            locus_start_position_index=int(
                                row["start_position_index"]
                            ),
                            locus_end_position_index=int(
                                row["end_position_index"]
                            ),
                        )
                        for row in membership_rows
                    ]
                    if membership_total > len(membership_rows):
                        limitations.append(
                            "ELSA membership summary is paginated: "
                            f"returned {len(membership_rows)} of {membership_total}; "
                            "use synteny_for_proteins() for the full result."
                        )
                synteny_state = "available"
            except SyntenyDatasetMismatchError as exc:
                synteny_state = "stale"
                limitations.append(str(exc))
            except Exception as exc:
                synteny_state = "failed"
                limitations.append(
                    "ELSA sidecar could not be read: "
                    f"{type(exc).__name__}: {exc}"
                )
        elif synteny_path is not None:
            synteny_state = "failed"
            limitations.append(
                "Explicit ELSA sidecar does not exist: "
                f"{Path(synteny_path).expanduser().resolve()}"
            )

    trace: dict[str, Any] = {
        "operator": "inspect",
        "parameters": {
            "entity_id": entity_id,
            "entity_type": entity.entity_type.value,
            "source_table": entity.source_table,
            "bin_id": bin_id,
            "window": default_window,
            "upstream_orfs": upstream,
            "downstream_orfs": downstream,
            "include_sequences": include_sequences,
            "allow_stale_synteny": allow_stale_synteny,
        },
        "database": str(db_path) if db_path is not None else "memory",
        "assembly_evidence_sidecar": (
            str(discovered_sidecar) if discovered_sidecar is not None else None
        ),
        "synteny_sidecar": (
            str(discovered_synteny) if discovered_synteny is not None else None
        ),
    }
    try:
        trace["schema_version"] = store.execute("SELECT MAX(version) FROM schema_version")[0][0]
    except Exception:
        trace["schema_version"] = None

    record = CaseRecord(
        entity=EntityReference(
            entity_id=entity.entity_id,
            entity_type=entity.entity_type,
            source_table=entity.source_table,
            name=entity.name,
            subtype=entity.subtype,
        ),
        bin=_bin_record(store, resolved_bin_id),
        contig=_contig_record(store, resolved_contig_id, resolved_bin_id),
        components=component_ids,
        context_window=context_window,
        proteins=proteins,
        named_calls=named_calls,
        assembly_evidence=assembly_record,
        assembly_evidence_state=assembly_state,
        synteny_memberships=synteny_memberships,
        synteny_state=synteny_state,
        limitations=list(dict.fromkeys(limitations)),
        trace=trace,
    )
    return BiologicalCase(
        store=store,
        record=record,
        assembly_evidence_path=discovered_sidecar,
    )


def _taxonomy_value(taxonomy: str | None, rank: str) -> str | None:
    prefix = _TAXONOMY_PREFIX.get(rank.lower())
    if prefix is None:
        raise ValueError(
            f"Unknown taxonomy rank {rank!r}; choose one of " + ", ".join(_TAXONOMY_PREFIX)
        )
    for part in (taxonomy or "").split(";"):
        if part.startswith(prefix) and part != prefix:
            return part
    return None


def _load_system_anchors(
    store: DuckDBStore,
    source_table: str,
    catalog: Mapping[str, set[str]],
) -> tuple[list[_Anchor], list[str]]:
    limitations: list[str] = []
    table = _quote_identifier(source_table)
    system_rows = _fetch_dicts(
        store,
        f"SELECT system_id, system_type FROM {table} ORDER BY system_id",
    )
    call_locations = _call_type_map(
        store,
        [str(row["system_id"]) for row in system_rows],
        catalog,
    )
    if "system_proteins" not in catalog:
        raise ValueError("Context comparison requires normalized system_proteins")
    memberships = _fetch_dicts(
        store,
        f"""
        SELECT s.system_id, s.system_type, p.protein_id, p.bin_id, p.contig_id,
               p.gene_index, p.strand
        FROM {table} s
        JOIN system_proteins sp ON sp.system_id = s.system_id
        JOIN proteins p ON p.protein_id = sp.protein_id
        WHERE p.gene_index IS NOT NULL
        ORDER BY s.system_id, sp.position NULLS LAST, p.gene_index
        """,
    )
    by_id: dict[str, list[dict[str, Any]]] = {}
    for row in memberships:
        by_id.setdefault(str(row["system_id"]), []).append(row)
    taxonomy_rows = {
        str(bin_id): taxonomy
        for bin_id, taxonomy in store.execute("SELECT bin_id, taxonomy FROM bins")
    }
    anchors: list[_Anchor] = []
    for system_row in system_rows:
        system_id = str(system_row["system_id"])
        locations = call_locations.get(system_id, [])
        if len(locations) != 1 or locations[0].get("source_table") != source_table:
            limitations.append(
                f"Skipped {system_id}: system ID is ambiguous across structured caller resources."
            )
            continue
        rows = by_id.get(system_id, [])
        replicons = {(str(row["bin_id"]), str(row["contig_id"])) for row in rows}
        if len(replicons) != 1:
            limitations.append(
                f"Skipped {system_id}: expected one replicon, found {len(replicons)}."
            )
            continue
        bin_id, contig_id = next(iter(replicons))
        anchors.append(
            _Anchor(
                entity_id=system_id,
                entity_name=str(system_row["system_type"]),
                bin_id=bin_id,
                contig_id=contig_id,
                min_gene_index=min(int(row["gene_index"]) for row in rows),
                max_gene_index=max(int(row["gene_index"]) for row in rows),
                orientation=_orientation(rows),
                component_count=len({str(row["protein_id"]) for row in rows}),
                taxonomy=(
                    str(taxonomy_rows[bin_id]) if taxonomy_rows.get(bin_id) is not None else None
                ),
                source_table=source_table,
            )
        )
    return anchors, limitations


def _load_locus_anchors(
    store: DuckDBStore,
    source_table: str,
    catalog: Mapping[str, set[str]],
) -> tuple[list[_Anchor], list[str]]:
    if source_table != "loci":
        raise ValueError("Context comparison currently requires locus_proteins-backed `loci` rows")
    if "locus_proteins" not in catalog:
        raise ValueError("Context comparison requires locus_proteins")
    rows = _fetch_dicts(
        store,
        """
        SELECT l.locus_id, l.locus_type, p.protein_id, p.bin_id, p.contig_id,
               p.gene_index, p.strand
        FROM loci l
        JOIN locus_proteins lp ON lp.locus_id = l.locus_id
        JOIN proteins p ON p.protein_id = lp.protein_id
        WHERE p.gene_index IS NOT NULL
        ORDER BY l.locus_id, lp.position, p.gene_index
        """,
    )
    by_id: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        by_id.setdefault(str(row["locus_id"]), []).append(row)
    taxonomy_rows = {
        str(bin_id): taxonomy
        for bin_id, taxonomy in store.execute("SELECT bin_id, taxonomy FROM bins")
    }
    anchors: list[_Anchor] = []
    limitations: list[str] = []
    for locus_id, locus_rows in by_id.items():
        replicons = {(str(row["bin_id"]), str(row["contig_id"])) for row in locus_rows}
        if len(replicons) != 1:
            limitations.append(
                f"Skipped {locus_id}: expected one replicon, found {len(replicons)}."
            )
            continue
        bin_id, contig_id = next(iter(replicons))
        anchors.append(
            _Anchor(
                entity_id=locus_id,
                entity_name=str(locus_rows[0]["locus_type"]),
                bin_id=bin_id,
                contig_id=contig_id,
                min_gene_index=min(int(row["gene_index"]) for row in locus_rows),
                max_gene_index=max(int(row["gene_index"]) for row in locus_rows),
                orientation=_orientation(locus_rows),
                component_count=len({str(row["protein_id"]) for row in locus_rows}),
                taxonomy=(
                    str(taxonomy_rows[bin_id]) if taxonomy_rows.get(bin_id) is not None else None
                ),
                source_table=source_table,
            )
        )
    return anchors, limitations


def _load_protein_anchors(
    store: DuckDBStore,
    entity_ids: Sequence[str],
) -> tuple[list[_Anchor], list[str]]:
    rows = _load_proteins(store, entity_ids)
    taxonomy_rows = {
        str(bin_id): taxonomy
        for bin_id, taxonomy in store.execute("SELECT bin_id, taxonomy FROM bins")
    }
    anchors = [
        _Anchor(
            entity_id=str(row["protein_id"]),
            entity_name=None,
            bin_id=str(row["bin_id"]),
            contig_id=str(row["contig_id"]),
            min_gene_index=int(row["gene_index"]),
            max_gene_index=int(row["gene_index"]),
            orientation=_orientation([row]),
            component_count=1,
            taxonomy=(
                str(taxonomy_rows[str(row["bin_id"])])
                if taxonomy_rows.get(str(row["bin_id"])) is not None
                else None
            ),
            source_table="proteins",
        )
        for row in rows
        if row.get("gene_index") is not None and row.get("bin_id") is not None
    ]
    return anchors, []


def _load_anchor_contexts(
    store: DuckDBStore,
    anchors: Sequence[_Anchor],
    *,
    upstream_orfs: int,
    downstream_orfs: int,
    batch_size: int = 5_000,
) -> dict[str, list[dict[str, Any]]]:
    """Fetch only requested windows, never complete replicons, for each anchor."""
    result: dict[str, list[dict[str, Any]]] = {}
    for start in range(0, len(anchors), batch_size):
        batch = list(anchors[start : start + batch_size])
        values = ", ".join(["(?, ?, ?, ?, ?)"] * len(batch))
        params: list[Any] = []
        for anchor in batch:
            requested_min, requested_max = _window_bounds(
                anchor.min_gene_index,
                anchor.max_gene_index,
                anchor.orientation,
                upstream_orfs,
                downstream_orfs,
            )
            params.extend(
                [
                    anchor.entity_id,
                    anchor.bin_id,
                    anchor.contig_id,
                    requested_min,
                    requested_max,
                ]
            )
        rows = _fetch_dicts(
            store,
            f"""
            WITH targets(
                entity_id, bin_id, contig_id, requested_min, requested_max
            ) AS (VALUES {values})
            SELECT t.entity_id, p.protein_id, p.bin_id, p.contig_id, p.gene_index,
                   p.start, p.end_coord, p.strand
            FROM targets t
            JOIN proteins p
              ON t.bin_id = p.bin_id AND t.contig_id = p.contig_id
             AND p.gene_index BETWEEN t.requested_min AND t.requested_max
            WHERE p.gene_index IS NOT NULL
            ORDER BY t.entity_id, p.gene_index, p.start
            """,
            params,
        )
        for row in rows:
            result.setdefault(str(row["entity_id"]), []).append(row)
    return result


def _feature_matches(
    feature: ContextFeature,
    protein_id: str,
    relative_orf: int,
    anchor: _Anchor,
    annotation_map: Mapping[
        str,
        Sequence[AnnotationEvidence | _ContextAnnotation],
    ],
    predicate_map: Mapping[str, list[str]],
    named_call_map: Mapping[str, list[NamedCallEvidence]],
    locus_map: Mapping[str, set[str]],
) -> bool:
    if feature.side == "upstream" and relative_orf >= 0:
        return False
    if feature.side == "downstream" and relative_orf <= 0:
        return False
    if feature.side == "anchor" and relative_orf != 0:
        return False
    if feature.max_orfs is not None and abs(relative_orf) > feature.max_orfs:
        return False

    if feature.kind == "annotation":
        return any(
            annotation.accession == feature.value
            and (feature.source is None or annotation.source.lower() == feature.source.lower())
            for annotation in annotation_map.get(protein_id, [])
        )
    if feature.kind == "annotation_name":
        needle = (feature.value or "").lower()
        return any(
            needle in (annotation.name or "").lower()
            or needle in (annotation.description or "").lower()
            for annotation in annotation_map.get(protein_id, [])
        )
    if feature.kind == "predicate":
        return feature.value in predicate_map.get(protein_id, [])
    if feature.kind in {"called_system", "other_called_system"}:
        for call in named_call_map.get(protein_id, []):
            if not feature.include_anchor_call and call.call_id == anchor.entity_id:
                continue
            if feature.kind == "other_called_system":
                return True
            if call.call_type == feature.value:
                return True
        return False
    if feature.kind == "locus_type":
        return feature.value in locus_map.get(protein_id, set())
    return False


def _log_choose(n: int, k: int) -> float:
    if k < 0 or k > n:
        return -math.inf
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def _hypergeom_probability(x: int, row_one: int, col_one: int, total: int) -> float:
    return math.exp(_log_hypergeom_probability(x, row_one, col_one, total))


def _log_hypergeom_probability(
    x: int,
    row_one: int,
    col_one: int,
    total: int,
) -> float:
    return (
        _log_choose(col_one, x)
        + _log_choose(total - col_one, row_one - x)
        - _log_choose(total, row_one)
    )


def _fisher_exact(
    a: int,
    b: int,
    c: int,
    d: int,
    alternative: Literal["greater", "less", "two-sided"],
) -> float:
    row_one = a + b
    col_one = a + c
    total = a + b + c + d
    minimum = max(0, row_one - (total - col_one))
    maximum = min(row_one, col_one)
    observed_log_probability = _log_hypergeom_probability(
        a,
        row_one,
        col_one,
        total,
    )
    if alternative == "greater":
        values = range(a, maximum + 1)
    elif alternative == "less":
        values = range(minimum, a + 1)
    else:
        values = [
            value
            for value in range(minimum, maximum + 1)
            if _log_hypergeom_probability(value, row_one, col_one, total)
            <= observed_log_probability + math.log1p(1e-12)
        ]
    log_probabilities = [
        _log_hypergeom_probability(value, row_one, col_one, total) for value in values
    ]
    maximum_log_probability = max(log_probabilities)
    log_sum = maximum_log_probability + math.log(
        sum(
            math.exp(log_probability - maximum_log_probability)
            for log_probability in log_probabilities
        )
    )
    return min(1.0, math.exp(log_sum))


def _odds_ratio_ci(
    a: int,
    b: int,
    c: int,
    d: int,
) -> tuple[float, tuple[float, float]]:
    aa, bb, cc, dd = (float(value) for value in (a, b, c, d))
    if 0 in {aa, bb, cc, dd}:
        aa += 0.5
        bb += 0.5
        cc += 0.5
        dd += 0.5
    odds_ratio = (aa * dd) / (bb * cc)
    standard_error = math.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    log_or = math.log(odds_ratio)
    return odds_ratio, (
        math.exp(log_or - 1.96 * standard_error),
        math.exp(log_or + 1.96 * standard_error),
    )


def _benjamini_hochberg(p_values: Sequence[float]) -> list[float]:
    if not p_values:
        return []
    order = sorted(range(len(p_values)), key=lambda index: p_values[index])
    adjusted = [1.0] * len(p_values)
    running = 1.0
    total = len(p_values)
    for reverse_rank, index in enumerate(reversed(order), start=1):
        rank = total - reverse_rank + 1
        running = min(running, p_values[index] * total / rank)
        adjusted[index] = min(1.0, running)
    return adjusted


def _statistic(
    feature_id: str,
    foreground_positive: int,
    foreground_total: int,
    background_positive: int,
    background_total: int,
    alternative: Literal["greater", "less", "two-sided"],
) -> ContextStatistic:
    a = foreground_positive
    b = foreground_total - foreground_positive
    c = background_positive
    d = background_total - background_positive
    odds_ratio, interval = _odds_ratio_ci(a, b, c, d)
    return ContextStatistic(
        feature_id=feature_id,
        foreground_positive=foreground_positive,
        foreground_total=foreground_total,
        background_positive=background_positive,
        background_total=background_total,
        odds_ratio=odds_ratio,
        odds_ratio_ci95=interval,
        fisher_p=_fisher_exact(a, b, c, d, alternative),
        alternative=alternative,
    )


def _comparison_context(
    case: BiologicalCase,
    *,
    features: Sequence[ContextFeature | str | dict[str, Any]],
    window: int,
    upstream_orfs: int | None,
    downstream_orfs: int | None,
    foreground_ids: Sequence[str] | None,
    background_ids: Sequence[str] | None,
    combine: Literal["all", "any"],
    min_components: int,
    require_full_context: bool,
    deduplicate_by: Literal["entity", "replicon", "bin"],
    exclude_foreground_units: bool,
    taxonomy_filter: str | None,
    same_taxonomy_rank: str | None,
    alternative: Literal["greater", "less", "two-sided"],
) -> ContextComparison:
    if not features:
        raise ValueError("At least one context feature is required")
    if combine not in {"all", "any"}:
        raise ValueError("combine must be 'all' or 'any'")
    if deduplicate_by not in {"entity", "replicon", "bin"}:
        raise ValueError("deduplicate_by must be 'entity', 'replicon', or 'bin'")
    if alternative not in {"greater", "less", "two-sided"}:
        raise ValueError("alternative must be 'greater', 'less', or 'two-sided'")
    parsed_features = [ContextFeature.parse(feature) for feature in features]
    feature_ids = [feature.feature_id for feature in parsed_features]
    if len(set(feature_ids)) != len(feature_ids):
        raise ValueError(
            "Context feature IDs must be unique; set explicit feature_id values "
            "when testing multiple variants of one feature"
        )
    explicit_overlap = set(foreground_ids or []) & set(background_ids or [])
    if explicit_overlap:
        preview = ", ".join(sorted(explicit_overlap)[:5])
        raise ValueError(f"Explicit foreground_ids and background_ids overlap, including {preview}")
    default_window, upstream, downstream = _resolve_windows(
        window,
        upstream_orfs,
        downstream_orfs,
    )
    catalog = _table_catalog(case.store)
    entity_type = case.record.entity.entity_type
    limitations: list[str] = []

    if entity_type == CaseEntityType.system:
        anchors, anchor_limitations = _load_system_anchors(
            case.store,
            case.record.entity.source_table,
            catalog,
        )
    elif entity_type == CaseEntityType.locus:
        anchors, anchor_limitations = _load_locus_anchors(
            case.store,
            case.record.entity.source_table,
            catalog,
        )
    elif entity_type == CaseEntityType.protein:
        if foreground_ids is None or background_ids is None:
            raise ValueError(
                "Protein context comparisons require explicit foreground_ids "
                "and background_ids; a single protein has no defensible automatic cohort."
            )
        anchors, anchor_limitations = _load_protein_anchors(
            case.store,
            list(foreground_ids) + list(background_ids),
        )
    else:
        raise ValueError("compare_context is supported for protein, system, and locus cases")
    limitations.extend(anchor_limitations)
    requested_ids = set(foreground_ids or []) | set(background_ids or [])
    if requested_ids:
        missing_ids = sorted(requested_ids - {anchor.entity_id for anchor in anchors})
        if missing_ids:
            preview = ", ".join(missing_ids[:5])
            raise ValueError(
                f"{len(missing_ids)} explicit cohort IDs are unavailable after "
                f"anchor validation, including {preview}"
            )

    focus_name = case.record.entity.name
    foreground_id_set = set(foreground_ids or [])
    background_id_set = set(background_ids or [])
    focus_taxonomy = (
        str(case.record.bin.get("taxonomy"))
        if case.record.bin and case.record.bin.get("taxonomy")
        else None
    )
    rank_value = _taxonomy_value(focus_taxonomy, same_taxonomy_rank) if same_taxonomy_rank else None
    if same_taxonomy_rank and rank_value is None:
        raise ValueError(
            f"The inspected case has no populated {same_taxonomy_rank} rank "
            "for same_taxonomy_rank matching"
        )

    selected: list[tuple[str, _Anchor]] = []
    for anchor in anchors:
        if anchor.component_count < min_components:
            continue
        if taxonomy_filter and taxonomy_filter.lower() not in (anchor.taxonomy or "").lower():
            continue
        if rank_value and _taxonomy_value(anchor.taxonomy, same_taxonomy_rank or "") != rank_value:
            continue
        if foreground_ids is not None:
            group = "foreground" if anchor.entity_id in foreground_id_set else None
        elif entity_type in {CaseEntityType.system, CaseEntityType.locus}:
            group = "foreground" if anchor.entity_name == focus_name else None
        else:
            group = None
        if group is None:
            if background_ids is not None:
                group = "background" if anchor.entity_id in background_id_set else None
            elif entity_type in {CaseEntityType.system, CaseEntityType.locus}:
                group = "background" if anchor.entity_name != focus_name else None
        if group is not None:
            selected.append((group, anchor))

    context_map = _load_anchor_contexts(
        case.store,
        [anchor for _, anchor in selected],
        upstream_orfs=upstream,
        downstream_orfs=downstream,
    )
    eligible: list[tuple[str, _Anchor, list[dict[str, Any]]]] = []
    edge_censored = {"foreground": 0, "background": 0}
    for group, anchor in selected:
        proteins = context_map.get(anchor.entity_id, [])
        if not proteins:
            continue
        relative_orfs = [
            _relative_orf(
                int(row["gene_index"]),
                anchor.min_gene_index,
                anchor.max_gene_index,
                anchor.orientation,
            )
            for row in proteins
        ]
        _, _, complete_upstream, complete_downstream = _context_completeness(
            relative_orfs,
            upstream,
            downstream,
        )
        complete = complete_upstream and complete_downstream
        if require_full_context and not complete:
            edge_censored[group] += 1
            continue
        eligible.append((group, anchor, proteins))

    if not any(group == "foreground" for group, _, _ in eligible):
        raise ValueError("No foreground entities remain after filters/context-edge requirements")
    if not any(group == "background" for group, _, _ in eligible):
        raise ValueError("No background entities remain after filters/context-edge requirements")

    protein_ids = sorted(
        {str(row["protein_id"]) for _, _, proteins in eligible for row in proteins}
    )
    annotation_map = (
        _load_feature_annotation_map(
            case.store,
            protein_ids,
            parsed_features,
            catalog,
        )
        if any(feature.kind in {"annotation", "annotation_name"} for feature in parsed_features)
        else {}
    )
    predicate_map: dict[str, list[str]] = {}
    if any(feature.kind == "predicate" for feature in parsed_features):
        for start in range(0, len(protein_ids), 5_000):
            predicate_map.update(
                get_active_predicates(
                    case.store,
                    protein_ids[start : start + 5_000],
                )
            )
    named_call_map = (
        _load_named_call_map(case.store, protein_ids, catalog)
        if any(
            feature.kind in {"called_system", "other_called_system"} for feature in parsed_features
        )
        else {}
    )
    locus_map = (
        _load_locus_memberships(case.store, protein_ids, catalog)
        if any(feature.kind == "locus_type" for feature in parsed_features)
        else {}
    )

    raw_matrix: list[dict[str, Any]] = []
    for group, anchor, proteins in eligible:
        feature_values: dict[str, bool] = {}
        for feature in parsed_features:
            feature_values[feature.feature_id or "feature"] = any(
                _feature_matches(
                    feature,
                    str(row["protein_id"]),
                    _relative_orf(
                        int(row["gene_index"]),
                        anchor.min_gene_index,
                        anchor.max_gene_index,
                        anchor.orientation,
                    ),
                    anchor,
                    annotation_map,
                    predicate_map,
                    named_call_map,
                    locus_map,
                )
                for row in proteins
            )
        composite = (
            all(feature_values.values()) if combine == "all" else any(feature_values.values())
        )
        if deduplicate_by == "entity":
            unit_key = anchor.entity_id
        elif deduplicate_by == "bin":
            unit_key = anchor.bin_id
        else:
            unit_key = f"{anchor.bin_id}|{anchor.contig_id}"
        raw_matrix.append(
            {
                "group": group,
                "unit_key": unit_key,
                "entity_ids": [anchor.entity_id],
                "bin_id": anchor.bin_id,
                "contig_id": anchor.contig_id,
                "features": feature_values,
                "composite": composite,
            }
        )

    grouped_matrix: dict[tuple[str, str], dict[str, Any]] = {}
    for row in raw_matrix:
        key = (str(row["group"]), str(row["unit_key"]))
        if key not in grouped_matrix:
            grouped_matrix[key] = row
            continue
        existing = grouped_matrix[key]
        existing["entity_ids"] = sorted(set(existing["entity_ids"]) | set(row["entity_ids"]))
        for feature_id, value in row["features"].items():
            existing["features"][feature_id] = existing["features"].get(feature_id, False) or value
        existing["composite"] = (
            all(existing["features"].values())
            if combine == "all"
            else any(existing["features"].values())
        )

    foreground_keys = {unit_key for group, unit_key in grouped_matrix if group == "foreground"}
    background_keys = {unit_key for group, unit_key in grouped_matrix if group == "background"}
    overlapping_units = foreground_keys & background_keys
    overlap_removed = 0
    if exclude_foreground_units:
        for key in list(grouped_matrix):
            group, unit_key = key
            if group == "background" and unit_key in foreground_keys:
                grouped_matrix.pop(key)
                overlap_removed += 1

    matrix = sorted(
        grouped_matrix.values(),
        key=lambda row: (str(row["group"]), str(row["unit_key"])),
    )
    foreground_rows = [row for row in matrix if row["group"] == "foreground"]
    background_rows = [row for row in matrix if row["group"] == "background"]
    if not background_rows:
        raise ValueError("No background units remain after foreground-overlap exclusion")

    def summary(
        group: Literal["foreground", "background"], rows: list[dict[str, Any]]
    ) -> ContextGroupSummary:
        return ContextGroupSummary(
            group=group,
            positive=sum(bool(row["composite"]) for row in rows),
            total=len(rows),
            entity_ids=sorted({entity_id for row in rows for entity_id in row["entity_ids"]}),
            unit_keys=[str(row["unit_key"]) for row in rows],
        )

    foreground_summary = summary("foreground", foreground_rows)
    background_summary = summary("background", background_rows)
    composite_statistic = _statistic(
        "composite",
        foreground_summary.positive,
        foreground_summary.total,
        background_summary.positive,
        background_summary.total,
        alternative,
    )
    per_feature: list[ContextStatistic] = []
    for feature in parsed_features:
        feature_id = feature.feature_id or "feature"
        per_feature.append(
            _statistic(
                feature_id,
                sum(bool(row["features"][feature_id]) for row in foreground_rows),
                len(foreground_rows),
                sum(bool(row["features"][feature_id]) for row in background_rows),
                len(background_rows),
                alternative,
            )
        )
    all_statistics = [composite_statistic, *per_feature]
    adjusted = _benjamini_hochberg([statistic.fisher_p for statistic in all_statistics])
    for statistic, q_value in zip(all_statistics, adjusted, strict=True):
        statistic.fisher_q = q_value

    if edge_censored["foreground"] or edge_censored["background"]:
        limitations.append(
            "Edge-censored entities excluded: "
            f"{edge_censored['foreground']} foreground and "
            f"{edge_censored['background']} background."
        )
    if overlap_removed:
        limitations.append(
            f"Removed {overlap_removed} background units that also contained a foreground anchor."
        )
    elif overlapping_units:
        limitations.append(
            f"Foreground and background share {len(overlapping_units)} counting "
            "units because overlap exclusion was disabled."
        )
    limitations.append(
        "Units are locus/replicon observations and are not guaranteed to be "
        "phylogenetically independent; exact p-values are not prevalence estimates."
    )

    recipe = {
        "method": "Sharur.inspect(...).compare_context",
        "case_entity_id": case.record.entity.entity_id,
        "case_entity_type": case.record.entity.entity_type.value,
        "source_table": case.record.entity.source_table,
        "features": [feature.model_dump(mode="json") for feature in parsed_features],
        "window": default_window,
        "upstream_orfs": upstream,
        "downstream_orfs": downstream,
        "foreground_ids": list(foreground_ids) if foreground_ids is not None else None,
        "background_ids": list(background_ids) if background_ids is not None else None,
        "combine": combine,
        "min_components": min_components,
        "require_full_context": require_full_context,
        "deduplicate_by": deduplicate_by,
        "exclude_foreground_units": exclude_foreground_units,
        "taxonomy_filter": taxonomy_filter,
        "same_taxonomy_rank": same_taxonomy_rank,
        "alternative": alternative,
    }
    context_window = ContextWindow(
        default_orfs=default_window,
        upstream_orfs=upstream,
        downstream_orfs=downstream,
        orientation=(
            case.record.context_window.orientation
            if case.record.context_window is not None
            else "unknown"
        ),
        orientation_basis=(
            "per-anchor co-oriented strands; coordinate-order fallback for mixed anchors"
        ),
        complete_upstream=require_full_context,
        complete_downstream=require_full_context,
    )
    return ContextComparison(
        case_entity_id=case.record.entity.entity_id,
        case_entity_type=case.record.entity.entity_type,
        foreground_definition={
            "mode": (
                "explicit_ids" if foreground_ids is not None else "same caller-emitted type as case"
            ),
            "entity_name": focus_name,
            "taxonomy_filter": taxonomy_filter,
            "same_taxonomy_rank": same_taxonomy_rank,
            "min_components": min_components,
        },
        background_definition={
            "mode": (
                "explicit_ids"
                if background_ids is not None
                else "other caller-emitted types in the same structured resource"
            ),
            "source_table": case.record.entity.source_table,
            "exclude_foreground_units": exclude_foreground_units,
        },
        features=parsed_features,
        combine=combine,
        context_window=context_window,
        require_full_context=require_full_context,
        deduplicate_by=deduplicate_by,
        foreground=foreground_summary,
        background=background_summary,
        composite=composite_statistic,
        per_feature=per_feature,
        entity_matrix=matrix,
        recipe=recipe,
        limitations=list(dict.fromkeys(limitations)),
    )


class BiologicalCase:
    """Connected case object returned by :meth:`Sharur.inspect`.

    The serializable evidence is available through :attr:`record`.  Methods
    that need the live dataset, such as comparison and plotting, retain the
    originating read-only store.
    """

    def __init__(
        self,
        *,
        store: DuckDBStore,
        record: CaseRecord,
        assembly_evidence_path: str | Path | None = None,
    ):
        self.store = store
        self.record = record
        self.assembly_evidence_path = (
            Path(assembly_evidence_path) if assembly_evidence_path is not None else None
        )

    @property
    def entity_id(self) -> str:
        return self.record.entity.entity_id

    def evidence(self) -> CaseRecord:
        """Return the typed, serializable evidence record."""
        return self.record

    def to_dict(self) -> dict[str, Any]:
        return self.record.model_dump(mode="json")

    def to_json(self, *, indent: int | None = 2) -> str:
        return self.record.to_json(indent=indent)

    def to_markdown(self) -> str:
        entity = self.record.entity

        def escape(value: Any) -> str:
            return str(value).replace("|", "\\|").replace("\n", " ")

        lines = [
            f"# Case: {entity.entity_id}",
            "",
            f"- Entity type: `{entity.entity_type.value}`",
            f"- Structured source: `{entity.source_table}`",
        ]
        if entity.name:
            lines.append(f"- Caller-emitted name: `{entity.name}`")
        if self.record.bin:
            lines.append(f"- Bin: `{self.record.bin.get('bin_id')}`")
            if self.record.bin.get("taxonomy"):
                lines.append(f"- Taxonomy: {self.record.bin['taxonomy']}")
        if self.record.contig:
            lines.append(f"- Contig: `{self.record.contig.get('contig_id')}`")
        if self.record.context_window:
            context = self.record.context_window
            lines.extend(
                [
                    "",
                    "## Context",
                    "",
                    f"- Requested: {context.upstream_orfs} upstream ORFs / "
                    f"{context.downstream_orfs} downstream ORFs",
                    f"- Realized: {context.realized_upstream_orfs} upstream / "
                    f"{context.realized_downstream_orfs} downstream",
                    f"- Orientation: `{context.orientation}` ({context.orientation_basis})",
                    f"- Proteins returned: {len(self.record.proteins)}",
                ]
            )
        if self.record.proteins:
            lines.extend(
                [
                    "",
                    "## Resolved proteins",
                    "",
                    "| Relative ORF | Component | Protein | Strand | Evidence |",
                    "|---:|:---:|---|:---:|---|",
                ]
            )
            for protein in self.record.proteins:
                evidence: list[str] = []
                seen_calls: set[tuple[str, str]] = set()
                for call in protein.named_calls:
                    key = (call.call_id, call.source_table)
                    if key in seen_calls:
                        continue
                    seen_calls.add(key)
                    evidence.append(
                        f"NAMED `{escape(call.call_type)}` (`{escape(call.source_table)}`)"
                    )
                evidence.extend(
                    f"OBSERVED `{escape(annotation.source)}:{escape(annotation.accession)}`"
                    for annotation in protein.annotations[:3]
                )
                if len(protein.annotations) > 3:
                    evidence.append(f"+{len(protein.annotations) - 3} observations")
                relative = "NA" if protein.relative_orf is None else str(protein.relative_orf)
                lines.append(
                    f"| {relative} | {'yes' if protein.is_component else ''} | "
                    f"`{escape(protein.protein_id)}` | {escape(protein.strand)} | "
                    f"{'; '.join(evidence) if evidence else 'none'} |"
                )
        if self.record.named_calls:
            lines.extend(["", "## Caller-emitted names", ""])
            seen = set()
            for call in self.record.named_calls:
                key = (call.call_id, call.call_type, call.source_table)
                if key in seen:
                    continue
                seen.add(key)
                lines.append(f"- `{call.call_type}` from `{call.source_table}` (`{call.call_id}`)")
        lines.extend(
            [
                "",
                "## Assembly evidence",
                "",
                f"- State: `{self.record.assembly_evidence_state}`",
            ]
        )
        if self.record.assembly_evidence is not None:
            payload = self.record.assembly_evidence.model_dump()
            fields = self.record.assembly_evidence.available_fields
            lines.append(f"- Populated fields: {', '.join(fields) if fields else 'none'}")
            for field_name in fields:
                lines.append(f"- `{field_name}`: {escape(payload[field_name])}")
        lines.extend(
            [
                "",
                "## ELSA synteny",
                "",
                f"- State: `{self.record.synteny_state}`",
            ]
        )
        if self.record.synteny_memberships:
            lines.extend(
                [
                    "",
                    "| Protein | Role | Run | Cluster | Blocks | Genomes |",
                    "|---|---|---|---|---:|---:|",
                ]
            )
            for membership in self.record.synteny_memberships:
                source = (
                    f" / source {membership.source_cluster_id}"
                    if membership.source_cluster_id is not None
                    else ""
                )
                lines.append(
                    f"| `{escape(membership.protein_id)}` | "
                    f"{membership.member_role} | `{membership.run_id}` | "
                    f"`{membership.cluster_key}{source}` | "
                    f"{membership.block_count} | {membership.genome_support} |"
                )
        if self.record.limitations:
            lines.extend(["", "## Limitations", ""])
            lines.extend(f"- {limitation}" for limitation in self.record.limitations)
        return "\n".join(lines)

    def __str__(self) -> str:
        return self.to_markdown()

    def compare_context(
        self,
        *,
        features: Sequence[ContextFeature | str | dict[str, Any]],
        window: int = DEFAULT_CONTEXT_ORFS,
        upstream_orfs: int | None = None,
        downstream_orfs: int | None = None,
        foreground_ids: Sequence[str] | None = None,
        background_ids: Sequence[str] | None = None,
        combine: Literal["all", "any"] = "all",
        min_components: int = 1,
        require_full_context: bool = True,
        deduplicate_by: Literal["entity", "replicon", "bin"] = "replicon",
        exclude_foreground_units: bool = True,
        taxonomy_filter: str | None = None,
        same_taxonomy_rank: str | None = None,
        alternative: Literal["greater", "less", "two-sided"] = "greater",
    ) -> ContextComparison:
        """Compare feature incidence in foreground versus background contexts.

        ``window`` is the symmetric default.  ``upstream_orfs`` and
        ``downstream_orfs`` independently override it.  For co-oriented
        anchors, biological strand defines upstream/downstream; mixed-strand
        anchors use ascending coordinate order and report that limitation.
        """
        if min_components < 1:
            raise ValueError("min_components must be at least 1")
        return _comparison_context(
            self,
            features=features,
            window=window,
            upstream_orfs=upstream_orfs,
            downstream_orfs=downstream_orfs,
            foreground_ids=foreground_ids,
            background_ids=background_ids,
            combine=combine,
            min_components=min_components,
            require_full_context=require_full_context,
            deduplicate_by=deduplicate_by,
            exclude_foreground_units=exclude_foreground_units,
            taxonomy_filter=taxonomy_filter,
            same_taxonomy_rank=same_taxonomy_rank,
            alternative=alternative,
        )

    def plot(
        self,
        output_path: str | Path | None = None,
        *,
        title: str | None = None,
        figure_width: float = 16,
    ) -> Path:
        """Render the already-resolved case context without re-querying anchors."""
        from sharur.operators.visualization import visualize_case  # noqa: PLC0415

        return visualize_case(
            self.record,
            output_path=output_path,
            title=title,
            figure_width=figure_width,
        )

    def propose_finding(
        self,
        *,
        title: str,
        category: str,
        description: str,
        claims: Sequence[Mapping[str, Any]],
        novelty: int = 0,
        falsification: str | None = None,
        literature_status: str = "unchecked",
        literature_references: Sequence[str] | None = None,
        comparison: ContextComparison | None = None,
    ) -> FindingDraft:
        """Create a claim-checked finding draft attached to this case."""
        from sharur.core.claim_compiler import FindingDraft  # noqa: PLC0415

        return FindingDraft.from_case(
            self,
            title=title,
            category=category,
            description=description,
            claims=claims,
            novelty=novelty,
            falsification=falsification,
            literature_status=literature_status,
            literature_references=literature_references,
            comparison=comparison,
        )

    def export(
        self,
        output_dir: str | Path,
        *,
        comparison: ContextComparison | None = None,
        finding: FindingDraft | None = None,
        include_sequences: bool = True,
        include_plot: bool = False,
        overwrite: bool = False,
    ) -> Path:
        """Write an atomic, compact evidence bundle."""
        from sharur.evidence_bundle import export_evidence_bundle  # noqa: PLC0415

        return export_evidence_bundle(
            self,
            output_dir,
            comparison=comparison,
            finding=finding,
            include_sequences=include_sequences,
            include_plot=include_plot,
            overwrite=overwrite,
        )


__all__ = [
    "BiologicalCase",
    "inspect_case",
    "load_context_evidence",
    "load_structured_context_evidence",
    "structured_projection_sources",
]
