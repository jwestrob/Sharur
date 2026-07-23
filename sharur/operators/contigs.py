"""Contig-native navigation and bounded genome-reading packets."""

from __future__ import annotations

import base64
import json
from typing import TYPE_CHECKING, Any

from sharur.operators.base import OperatorContext, SharurResult
from sharur.operators.cases import (
    load_context_evidence,
    structured_projection_sources,
)
from sharur.operators.semantics import get_active_predicates


if TYPE_CHECKING:
    from sharur.storage.duckdb_store import DuckDBStore


PACKET_VERSION = "1.0"
_NULL_GENE_INDEX = 9_223_372_036_854_775_807


def _dataset_id(store: DuckDBStore) -> str | None:
    value = getattr(store, "dataset_id", None)
    return str(value) if value else None


def _require_contig_navigation_index(store: DuckDBStore) -> None:
    """Fail before a genome-filtered full scan on legacy large databases."""
    cached = getattr(store, "_sharur_has_contig_navigation_index", None)
    if cached is True:
        return
    rows = store.execute(
        """
        SELECT 1
        FROM duckdb_indexes()
        WHERE index_name = 'idx_contigs_bin'
        LIMIT 1
        """
    )
    if not rows:
        raise RuntimeError(
            "Contig pagination requires idx_contigs_bin. Run `sharur migrate "
            "--db PATH`, then rebuild the dataset seal and query replica."
        )
    if getattr(store, "read_only", False):
        store._sharur_has_contig_navigation_index = True


def _encode_cursor(payload: dict[str, Any]) -> str:
    encoded = base64.urlsafe_b64encode(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).decode("ascii")
    return encoded.rstrip("=")


def _decode_cursor(
    token: str,
    *,
    kind: str,
    store: DuckDBStore,
    genome_id: str,
    contig_id: str | None = None,
) -> dict[str, Any]:
    if len(token) > 4_096:
        raise ValueError("Cursor exceeds 4,096 characters")
    try:
        padding = "=" * (-len(token) % 4)
        payload = json.loads(
            base64.urlsafe_b64decode(token + padding).decode("utf-8")
        )
    except (ValueError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError("Cursor is malformed") from exc
    if not isinstance(payload, dict) or payload.get("v") != 1:
        raise ValueError("Cursor version is unsupported")
    expected = {
        "kind": kind,
        "dataset_id": _dataset_id(store),
        "genome_id": genome_id,
    }
    if contig_id is not None:
        expected["contig_id"] = contig_id
    for key, value in expected.items():
        if payload.get(key) != value:
            raise ValueError(f"Cursor {key} does not match this request")
    return payload


def list_contigs(
    store: DuckDBStore,
    genome_id: str,
    *,
    limit: int = 100,
    cursor: str | None = None,
) -> SharurResult:
    """List one genome's contigs with deterministic keyset pagination."""
    if limit < 1 or limit > 1_000:
        raise ValueError("limit must be between 1 and 1,000")
    params = {
        "genome_id": genome_id,
        "limit": limit,
        "cursor": cursor,
    }
    with OperatorContext("list_contigs", params, store=store) as ctx:
        _require_contig_navigation_index(store)
        if not store.execute("SELECT 1 FROM bins WHERE bin_id = ?", [genome_id]):
            return ctx.make_result(
                data=f"Genome not found: {genome_id}",
                rows=0,
                total_rows=0,
                raw=[],
            )
        after_contig = ""
        if cursor:
            payload = _decode_cursor(
                cursor,
                kind="contigs",
                store=store,
                genome_id=genome_id,
            )
            after_contig = payload.get("last_contig_id")
            if not isinstance(after_contig, str):
                raise ValueError("Contig cursor lacks a valid position")

        total_rows = int(
            store.execute(
                "SELECT COUNT(*) FROM contigs WHERE bin_id = ?",
                [genome_id],
            )[0][0]
        )
        rows = store.execute(
            """
            WITH selected_contigs AS (
                SELECT
                    contig_id,
                    bin_id,
                    length,
                    gc_content,
                    is_circular,
                    taxonomy
                FROM contigs
                WHERE bin_id = ?
                  AND contig_id > ?
                ORDER BY contig_id
                LIMIT ?
            )
            SELECT
                c.contig_id,
                c.length,
                c.gc_content,
                c.is_circular,
                c.taxonomy,
                COUNT(p.protein_id) AS protein_count
            FROM selected_contigs AS c
            LEFT JOIN proteins AS p
              ON p.bin_id = c.bin_id AND p.contig_id = c.contig_id
            GROUP BY
                c.contig_id,
                c.length,
                c.gc_content,
                c.is_circular,
                c.taxonomy
            ORDER BY c.contig_id
            """,
            [genome_id, after_contig, limit + 1],
        )
        has_more = len(rows) > limit
        page = rows[:limit]
        next_cursor = None
        if has_more and page:
            next_cursor = _encode_cursor(
                {
                    "v": 1,
                    "kind": "contigs",
                    "dataset_id": _dataset_id(store),
                    "genome_id": genome_id,
                    "last_contig_id": str(page[-1][0]),
                }
            )
        records = [
            {
                "genome_id": genome_id,
                "contig_id": str(contig_id),
                "length": int(length),
                "gc_content": float(gc_content) if gc_content is not None else None,
                "is_circular": bool(is_circular),
                "taxonomy": taxonomy,
                "protein_count": int(protein_count),
            }
            for (
                contig_id,
                length,
                gc_content,
                is_circular,
                taxonomy,
                protein_count,
            ) in page
        ]
        lines = [
            f"# Contigs: {genome_id}",
            f"Page contains {len(records)} of {total_rows:,} contigs",
            "",
            "| Contig ID | Length | Proteins | Circular |",
            "|-----------|-------:|---------:|:--------:|",
        ]
        for record in records:
            lines.append(
                f"| {record['contig_id']} | {record['length']:,} | "
                f"{record['protein_count']:,} | "
                f"{'yes' if record['is_circular'] else 'no'} |"
            )
        return ctx.make_result(
            data="\n".join(lines),
            rows=len(records),
            total_rows=total_rows,
            truncated=has_more,
            ref=next_cursor,
            raw=records,
        )


def get_contig(
    store: DuckDBStore,
    genome_id: str,
    contig_id: str,
    *,
    verbosity: int = 1,
) -> SharurResult:
    """Return contig metadata and observed-annotation coverage."""
    if verbosity < 0 or verbosity > 2:
        raise ValueError("verbosity must be between 0 and 2")
    params = {
        "genome_id": genome_id,
        "contig_id": contig_id,
        "verbosity": verbosity,
    }
    with OperatorContext("get_contig", params, store=store) as ctx:
        rows = store.execute(
            """
            SELECT contig_id, bin_id, length, gc_content, is_circular, taxonomy
            FROM contigs
            WHERE bin_id = ? AND contig_id = ?
            """,
            [genome_id, contig_id],
        )
        if not rows:
            return ctx.make_result(
                data=f"Contig not found in genome {genome_id}: {contig_id}",
                rows=0,
                total_rows=0,
            )
        _, bin_id, length, gc_content, is_circular, taxonomy = rows[0]
        protein_count = int(
            store.execute(
                """
                SELECT COUNT(*)
                FROM proteins
                WHERE bin_id = ? AND contig_id = ?
                """,
                [genome_id, contig_id],
            )[0][0]
        )
        excluded_sources = structured_projection_sources(store)
        exclusion_sql = ""
        query_params: list[Any] = [genome_id, contig_id]
        if excluded_sources:
            placeholders = ", ".join(["?"] * len(excluded_sources))
            exclusion_sql = f"AND a.source NOT IN ({placeholders})"
            query_params.extend(excluded_sources)
        stats = store.execute(
            f"""
            SELECT
                a.source,
                COUNT(*) AS domain_hits,
                COUNT(DISTINCT a.protein_id) AS proteins
            FROM proteins AS p
            JOIN annotations AS a ON a.protein_id = p.protein_id
            WHERE p.bin_id = ? AND p.contig_id = ?
              {exclusion_sql}
            GROUP BY a.source
            ORDER BY proteins DESC, domain_hits DESC, a.source
            """,
            query_params,
        )
        annotation_stats = [
            {
                "source": str(source),
                "domain_hits": int(domain_hits),
                "proteins": int(annotated_proteins),
            }
            for source, domain_hits, annotated_proteins in stats
        ]
        lines = [
            f"# Contig: {contig_id}",
            f"- Genome: {bin_id}",
            f"- Length: {int(length):,} bp",
            f"- Proteins: {protein_count:,}",
            f"- Circular: {'yes' if is_circular else 'no'}",
        ]
        if gc_content is not None:
            lines.append(f"- GC content: {float(gc_content):.1%}")
        if taxonomy:
            lines.append(f"- Taxonomy: {taxonomy}")
        if verbosity >= 1 and annotation_stats:
            lines.extend(
                [
                    "",
                    "## Observed annotations",
                    "| Source | Proteins | Domain hits |",
                    "|--------|---------:|------------:|",
                ]
            )
            for stat in annotation_stats:
                lines.append(
                    f"| {stat['source']} | {stat['proteins']:,} | "
                    f"{stat['domain_hits']:,} |"
                )
        raw = {
            "genome_id": str(bin_id),
            "contig_id": contig_id,
            "length": int(length),
            "gc_content": float(gc_content) if gc_content is not None else None,
            "is_circular": bool(is_circular),
            "taxonomy": taxonomy,
            "protein_count": protein_count,
            "annotation_stats": annotation_stats,
        }
        return ctx.make_result(
            data="\n".join(lines),
            rows=1,
            total_rows=1,
            raw=raw,
        )


def get_contig_packet(
    store: DuckDBStore,
    genome_id: str,
    contig_id: str,
    *,
    cursor: str | None = None,
    limit: int = 100,
    all_annotations: bool = True,
) -> SharurResult:
    """Return one bounded, sequence-free packet for exhaustive genome reading."""
    if limit < 1 or limit > 250:
        raise ValueError("limit must be between 1 and 250")
    params = {
        "genome_id": genome_id,
        "contig_id": contig_id,
        "cursor": cursor,
        "limit": limit,
        "all_annotations": all_annotations,
        "packet_version": PACKET_VERSION,
    }
    with OperatorContext("get_contig_packet", params, store=store) as ctx:
        contig_rows = store.execute(
            """
            SELECT c.length, COUNT(p.protein_id)
            FROM contigs AS c
            LEFT JOIN proteins AS p
              ON p.contig_id = c.contig_id AND p.bin_id = c.bin_id
            WHERE c.bin_id = ? AND c.contig_id = ?
            GROUP BY c.length
            """,
            [genome_id, contig_id],
        )
        if not contig_rows:
            return ctx.make_result(
                data=f"Contig not found in genome {genome_id}: {contig_id}",
                rows=0,
                total_rows=0,
                raw={
                    "packet_version": PACKET_VERSION,
                    "dataset_id": _dataset_id(store),
                    "genome_id": genome_id,
                    "contig_id": contig_id,
                    "cursor": cursor,
                    "next_cursor": None,
                    "complete": True,
                    "proteins": [],
                },
            )
        contig_length, total_rows = contig_rows[0]
        last_position = (0, -1, -1, "")
        if cursor:
            payload = _decode_cursor(
                cursor,
                kind="contig_packet",
                store=store,
                genome_id=genome_id,
                contig_id=contig_id,
            )
            position = payload.get("position")
            if (
                not isinstance(position, list)
                or len(position) != 4
                or not all(isinstance(value, int) for value in position[:3])
                or not isinstance(position[3], str)
            ):
                raise ValueError("Packet cursor lacks a valid position")
            last_position = (
                int(position[0]),
                int(position[1]),
                int(position[2]),
                position[3],
            )

        rows = store.execute(
            """
            SELECT
                p.protein_id,
                p.gene_index,
                p.start,
                p.end_coord,
                p.strand,
                CAST(
                    COALESCE(
                        p.sequence_length,
                        (p.end_coord - p.start) / 3
                    )
                    AS BIGINT
                ) AS length_aa,
                p.gc_content,
                CAST(p.gene_index IS NULL AS INTEGER) AS missing_gene_index,
                COALESCE(CAST(p.gene_index AS BIGINT), ?) AS gene_order
            FROM proteins AS p
            WHERE p.bin_id = ? AND p.contig_id = ?
              AND (
                    CAST(p.gene_index IS NULL AS INTEGER),
                    COALESCE(CAST(p.gene_index AS BIGINT), ?),
                    p.start,
                    p.protein_id
                  ) > (?, ?, ?, ?)
            ORDER BY
                p.gene_index NULLS LAST,
                p.start,
                p.protein_id
            LIMIT ?
            """,
            [
                _NULL_GENE_INDEX,
                genome_id,
                contig_id,
                _NULL_GENE_INDEX,
                *last_position,
                limit + 1,
            ],
        )
        has_more = len(rows) > limit
        page = rows[:limit]
        protein_ids = [str(row[0]) for row in page]
        context = load_context_evidence(store, protein_ids)
        predicates = get_active_predicates(store, protein_ids)
        observed = context["observed_annotations"]
        named_calls = context["named_calls"]
        loci = context["loci"]
        proteins = []
        for (
            protein_id,
            gene_index,
            start,
            end_coord,
            strand,
            length_aa,
            gc_content,
            _missing_gene_index,
            _gene_order,
        ) in page:
            protein_annotations = observed.get(str(protein_id), [])
            if not all_annotations:
                protein_annotations = protein_annotations[:1]
            proteins.append(
                {
                    "protein_id": str(protein_id),
                    "gene_index": (
                        int(gene_index) if gene_index is not None else None
                    ),
                    "start": int(start),
                    "end": int(end_coord),
                    "strand": str(strand),
                    "length_aa": int(length_aa),
                    "gc_content": (
                        float(gc_content) if gc_content is not None else None
                    ),
                    "observed_annotations": protein_annotations,
                    "predicates": predicates.get(str(protein_id), []),
                    "named_calls": named_calls.get(str(protein_id), []),
                    "loci": loci.get(str(protein_id), []),
                }
            )

        next_cursor = None
        if has_more and page:
            final = page[-1]
            next_cursor = _encode_cursor(
                {
                    "v": 1,
                    "kind": "contig_packet",
                    "dataset_id": _dataset_id(store),
                    "genome_id": genome_id,
                    "contig_id": contig_id,
                    "position": [
                        int(final[7]),
                        int(final[8]),
                        int(final[2]),
                        str(final[0]),
                    ],
                }
            )
        complete = not has_more
        raw = {
            "packet_version": PACKET_VERSION,
            "dataset_id": _dataset_id(store),
            "genome_id": genome_id,
            "contig_id": contig_id,
            "contig_length": int(contig_length),
            "cursor": cursor,
            "next_cursor": next_cursor,
            "complete": complete,
            "proteins": proteins,
        }
        lines = [
            f"# Contig packet: {genome_id}/{contig_id}",
            f"- Packet version: {PACKET_VERSION}",
            f"- Proteins in packet: {len(proteins):,}",
            f"- Proteins on contig: {int(total_rows):,}",
            f"- Complete: {'yes' if complete else 'no'}",
            "- Sequence payload: compute-only",
        ]
        return ctx.make_result(
            data="\n".join(lines),
            rows=len(proteins),
            total_rows=int(total_rows),
            truncated=has_more,
            ref=next_cursor,
            raw=raw,
        )


__all__ = [
    "PACKET_VERSION",
    "get_contig",
    "get_contig_packet",
    "list_contigs",
]
