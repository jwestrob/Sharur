"""Contig-native navigation and bounded genome-reading packets."""

from __future__ import annotations

import base64
import hashlib
import json
import math
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
GENOME_PACKET_VERSION = "genome-packet/2.0"
GENOME_MODEL_PAYLOAD_SCHEMA = "atlas-model-packet/1.0"
GENOME_PACKET_PACKING_SCHEMA = "genome-packet-packing/1.1"
DEFAULT_GENOME_PACKET_BYTES = 512 * 1024
MAX_GENOME_PACKET_BYTES = 1_500_000
_NULL_GENE_INDEX = 9_223_372_036_854_775_807


def _dataset_id(store: DuckDBStore) -> str | None:
    value = getattr(store, "dataset_id", None)
    return str(value) if value else None


def _canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    ).encode("utf-8")


def _sha256(value: Any) -> str:
    return hashlib.sha256(_canonical_bytes(value)).hexdigest()


def _token_scenarios(payload_bytes: int) -> dict[str, int]:
    """Return transparent byte-to-token scenarios, not tokenizer claims."""
    return {
        "at_2_bytes_per_token": math.ceil(payload_bytes / 2),
        "at_3_bytes_per_token": math.ceil(payload_bytes / 3),
        "at_4_bytes_per_token": math.ceil(payload_bytes / 4),
    }


def minimum_genome_packet_record_bytes() -> dict[str, int]:
    """Return canonical schema floors for proportional request guardrails."""
    protein = {
        "protein_id": "",
        "gene_index": 0,
        "start": 0,
        "end": 0,
        "strand": "",
        "length_aa": 0,
        "gc_content": 0,
        "observed_annotations": [],
        "predicates": [],
        "named_calls": [],
        "loci": [],
    }
    contig = {
        "bin_id": "",
        "contig_id": "",
        "length": 0,
        "gc_content": 0,
        "is_circular": False,
        "taxonomy": "",
        "total_protein_count": 0,
        "protein_offset_start": 0,
        "protein_offset_end": 0,
        "segment_index": 0,
        "complete": True,
        "proteins": [],
    }
    return {
        "protein": len(_canonical_bytes(protein)),
        "contig": len(_canonical_bytes(contig)),
    }


_GENOME_PACKET_RECORD_FLOORS = minimum_genome_packet_record_bytes()
DEFAULT_GENOME_PACKET_CONTIGS = (
    DEFAULT_GENOME_PACKET_BYTES // _GENOME_PACKET_RECORD_FLOORS["contig"]
)
DEFAULT_GENOME_PACKET_PROTEINS = (
    DEFAULT_GENOME_PACKET_BYTES // _GENOME_PACKET_RECORD_FLOORS["protein"]
)
MAX_GENOME_PACKET_CONTIGS = (
    MAX_GENOME_PACKET_BYTES // _GENOME_PACKET_RECORD_FLOORS["contig"]
)
MAX_GENOME_PACKET_PROTEINS = (
    MAX_GENOME_PACKET_BYTES // _GENOME_PACKET_RECORD_FLOORS["protein"]
)


def genome_packet_packing_contract(
    *,
    max_contigs: int | None = None,
    max_proteins: int | None = None,
    max_model_payload_bytes: int = DEFAULT_GENOME_PACKET_BYTES,
    all_annotations: bool = True,
) -> dict[str, Any]:
    """Build the immutable contract used to pack one bin's model packets."""
    if max_model_payload_bytes < 1_024:
        raise ValueError("max_model_payload_bytes must be at least 1,024")
    if max_model_payload_bytes > MAX_GENOME_PACKET_BYTES:
        raise ValueError(
            "max_model_payload_bytes must be at most "
            f"{MAX_GENOME_PACKET_BYTES:,}"
        )
    record_floors = dict(_GENOME_PACKET_RECORD_FLOORS)
    proportional_contig_limit = max(
        1,
        max_model_payload_bytes // record_floors["contig"],
    )
    proportional_protein_limit = max(
        1,
        max_model_payload_bytes // record_floors["protein"],
    )
    resolved_contigs = (
        proportional_contig_limit if max_contigs is None else max_contigs
    )
    resolved_proteins = (
        proportional_protein_limit if max_proteins is None else max_proteins
    )
    if resolved_contigs < 1 or resolved_contigs > MAX_GENOME_PACKET_CONTIGS:
        raise ValueError(
            f"max_contigs must be between 1 and {MAX_GENOME_PACKET_CONTIGS:,}"
        )
    if resolved_proteins < 1 or resolved_proteins > MAX_GENOME_PACKET_PROTEINS:
        raise ValueError(
            f"max_proteins must be between 1 and {MAX_GENOME_PACKET_PROTEINS:,}"
        )
    if resolved_contigs > proportional_contig_limit:
        raise ValueError(
            "max_contigs exceeds the canonical payload-proportional safety "
            f"limit of {proportional_contig_limit:,}"
        )
    if resolved_proteins > proportional_protein_limit:
        raise ValueError(
            "max_proteins exceeds the canonical payload-proportional safety "
            f"limit of {proportional_protein_limit:,}"
        )
    return {
        "schema_version": GENOME_PACKET_PACKING_SCHEMA,
        "packet_version": GENOME_PACKET_VERSION,
        "strategy": "ordered-whole-contig-first-split-oversize/1.0",
        "bin_isolation": "exact bins.bin_id",
        "contig_order": "contig_id",
        "protein_order": "gene_index NULLS LAST, start, protein_id",
        "max_contigs": resolved_contigs,
        "max_proteins": resolved_proteins,
        "max_model_payload_bytes": max_model_payload_bytes,
        "all_annotations": all_annotations,
        "serialized_record_floor_bytes": record_floors,
        "count_limit_rule": (
            "floor(max_model_payload_bytes / canonical serialized record floor)"
        ),
    }


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


PROTEIN_COLUMNS = ["protein_id", "idx", "start", "end", "strand", "len", "predicates", "annotations"]


def _compact_protein_row(protein: dict[str, Any]) -> list[Any]:
    """One protein as a positional row.

    A dict per protein repeats every key name once per record. At ~480 proteins
    a frame that is 16% of the payload, paid for information the column header
    already carries. Redundancy inside each record costs more again: the
    annotation objects repeat the parent `protein_id`, duplicate `accession`
    into `name`, carry empty `description`, and serialise scores at float64
    precision that no consumer reads.

    Everything the scan prompt refers to is preserved. `evidence_level` is
    emitted only when it is not the `observed` default, which keeps the
    purpose-built-caller boundary visible while dropping it from the ~98% of
    annotations where it says nothing.
    """
    annotations: list[str] = []
    for ann in protein.get("observed_annotations") or []:
        text = str(ann.get("accession") or "")
        if ann.get("start_aa") is not None:
            text += f"@{ann['start_aa']}-{ann['end_aa']}"
        level = ann.get("evidence_level")
        if level and level != "observed":
            text += f"!{level}"
        annotations.append(text)
    for call in protein.get("named_calls") or []:
        annotations.append(f"NAMED:{call.get('call_type')}")
    for locus in protein.get("loci") or []:
        annotations.append(f"LOCUS:{locus.get('locus_type')}")
    # `unannotated` restates an empty annotation list, which is already visible.
    predicates = [p for p in (protein.get("predicates") or []) if p != "unannotated"]
    return [
        protein.get("protein_id"),
        protein.get("gene_index"),
        protein.get("start"),
        protein.get("end"),
        protein.get("strand"),
        protein.get("length_aa"),
        ",".join(predicates),
        ";".join(annotations),
    ]


def _genome_model_payload(
    *,
    dataset_id: str | None,
    genome: dict[str, Any],
    frame_index: int,
    contigs: list[dict[str, Any]],
) -> dict[str, Any]:
    compact_contigs = [
        {
            # bin_id is retained per contig: the packet census uses it as a
            # mixed-bin guard, and at one field per contig it costs nothing.
            "bin_id": contig.get("bin_id"),
            "contig_id": contig.get("contig_id"),
            "proteins": [_compact_protein_row(p) for p in contig.get("proteins") or []],
        }
        for contig in contigs
    ]
    return {
        "schema_version": GENOME_MODEL_PAYLOAD_SCHEMA,
        "packet_version": GENOME_PACKET_VERSION,
        "dataset_id": dataset_id,
        "bin_id": genome["bin_id"],
        "frame_index": frame_index,
        "genome": genome,
        "protein_columns": PROTEIN_COLUMNS,
        "contigs": compact_contigs,
    }


def _packet_position(value: Any, *, field: str) -> tuple[int, int, int, str]:
    if (
        not isinstance(value, list)
        or len(value) != 4
        or not all(isinstance(item, int) for item in value[:3])
        or not isinstance(value[3], str)
    ):
        raise ValueError(f"Genome packet cursor lacks a valid {field}")
    return int(value[0]), int(value[1]), int(value[2]), value[3]


def _genome_packet_metadata(
    store: DuckDBStore,
    genome_id: str,
) -> tuple[Any, ...] | None:
    cache = getattr(store, "_sharur_genome_packet_metadata", None)
    if cache is not None and genome_id in cache:
        return cache[genome_id]
    rows = store.execute(
        """
        SELECT
            b.completeness,
            b.contamination,
            b.taxonomy,
            (SELECT COUNT(*) FROM contigs c WHERE c.bin_id = b.bin_id),
            (SELECT COUNT(*) FROM proteins p WHERE p.bin_id = b.bin_id)
        FROM bins AS b
        WHERE b.bin_id = ?
        """,
        [genome_id],
    )
    value = rows[0] if rows else None
    if getattr(store, "read_only", False):
        if cache is None:
            cache = {}
            store._sharur_genome_packet_metadata = cache
        cache[genome_id] = value
    return value


def get_genome_packet(
    store: DuckDBStore,
    genome_id: str,
    *,
    cursor: str | None = None,
    max_contigs: int | None = None,
    max_proteins: int | None = None,
    max_model_payload_bytes: int = DEFAULT_GENOME_PACKET_BYTES,
    all_annotations: bool = True,
) -> SharurResult:
    """Pack consecutive sequence-free contig records from exactly one bin.

    Whole contigs are retained when they fit a fresh packet. A contig is split
    only when its remaining protein/evidence payload exceeds the packing
    contract. The opaque continuation cursor is bound to the sealed dataset,
    exact ``bins.bin_id``, and exact packing-contract hash.
    """
    packing = genome_packet_packing_contract(
        max_contigs=max_contigs,
        max_proteins=max_proteins,
        max_model_payload_bytes=max_model_payload_bytes,
        all_annotations=all_annotations,
    )
    max_contigs = int(packing["max_contigs"])
    max_proteins = int(packing["max_proteins"])
    packing_hash = _sha256(packing)
    params = {
        "genome_id": genome_id,
        "cursor": cursor,
        **packing,
        "packing_contract_hash": packing_hash,
    }
    with OperatorContext("get_genome_packet", params, store=store) as ctx:
        _require_contig_navigation_index(store)
        bin_row = _genome_packet_metadata(store, genome_id)
        if bin_row is None:
            genome = {
                "bin_id": genome_id,
                "completeness": None,
                "contamination": None,
                "taxonomy": None,
                "total_contig_count": 0,
                "total_protein_count": 0,
            }
            model_payload = _genome_model_payload(
                dataset_id=_dataset_id(store),
                genome=genome,
                frame_index=0,
                contigs=[],
            )
            payload_bytes = len(_canonical_bytes(model_payload))
            raw = {
                "packet_version": GENOME_PACKET_VERSION,
                "dataset_id": _dataset_id(store),
                "genome_id": genome_id,
                "bin_id": genome_id,
                "frame_id": None,
                "frame_index": 0,
                "cursor": cursor,
                "next_cursor": None,
                "complete": True,
                "packing_contract": packing,
                "packing_contract_hash": packing_hash,
                "model_payload": model_payload,
                "model_payload_bytes": payload_bytes,
                "model_payload_sha256": _sha256(model_payload),
                "token_scenarios": _token_scenarios(payload_bytes),
                "target_exceeded": False,
                "coverage_receipts": [],
            }
            return ctx.make_result(
                data=f"Genome not found: {genome_id}",
                rows=0,
                total_rows=0,
                raw=raw,
            )

        (
            completeness,
            contamination,
            taxonomy,
            total_contigs_raw,
            total_proteins_raw,
        ) = bin_row
        total_contigs = int(total_contigs_raw)
        total_proteins = int(total_proteins_raw)
        genome = {
            "bin_id": genome_id,
            "completeness": (
                float(completeness) if completeness is not None else None
            ),
            "contamination": (
                float(contamination) if contamination is not None else None
            ),
            "taxonomy": taxonomy,
            "total_contig_count": total_contigs,
            "total_protein_count": total_proteins,
        }

        frame_index = 0
        after_contig_id = ""
        current_contig_id: str | None = None
        current_position: tuple[int, int, int, str] | None = None
        consumed_in_current = 0
        current_segment_index = 0
        completed_contigs = 0
        completed_proteins = 0
        if cursor:
            payload = _decode_cursor(
                cursor,
                kind="genome_packet",
                store=store,
                genome_id=genome_id,
            )
            if payload.get("packing_contract_hash") != packing_hash:
                raise ValueError(
                    "Cursor packing_contract_hash does not match this request"
                )
            integer_fields = {
                "frame_index": payload.get("frame_index"),
                "consumed_in_current": payload.get("consumed_in_current"),
                "current_segment_index": payload.get("current_segment_index"),
                "completed_contigs": payload.get("completed_contigs"),
                "completed_proteins": payload.get("completed_proteins"),
            }
            if any(
                isinstance(value, bool) or not isinstance(value, int) or value < 0
                for value in integer_fields.values()
            ):
                raise ValueError("Genome packet cursor has invalid counters")
            frame_index = int(integer_fields["frame_index"])
            consumed_in_current = int(integer_fields["consumed_in_current"])
            current_segment_index = int(integer_fields["current_segment_index"])
            completed_contigs = int(integer_fields["completed_contigs"])
            completed_proteins = int(integer_fields["completed_proteins"])
            after_value = payload.get("after_contig_id")
            if not isinstance(after_value, str):
                raise ValueError("Genome packet cursor lacks after_contig_id")
            after_contig_id = after_value
            current_value = payload.get("current_contig_id")
            if current_value is not None and not isinstance(current_value, str):
                raise ValueError("Genome packet cursor has invalid current_contig_id")
            current_contig_id = current_value
            if current_contig_id is None:
                if (
                    payload.get("position") is not None
                    or consumed_in_current != 0
                    or current_segment_index != 0
                ):
                    raise ValueError(
                        "Genome packet cursor has partial-contig state without "
                        "current_contig_id"
                    )
            else:
                current_position = _packet_position(
                    payload.get("position"),
                    field="position",
                )
                if consumed_in_current < 1 or current_segment_index < 1:
                    raise ValueError("Genome packet cursor has incomplete split state")
                if current_contig_id <= after_contig_id:
                    raise ValueError(
                        "Genome packet cursor current contig is outside its ordered suffix"
                    )
            if completed_contigs > total_contigs or completed_proteins > total_proteins:
                raise ValueError("Genome packet cursor counters exceed live bin totals")

        if current_contig_id is not None:
            contig_rows = store.execute(
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
                    WHERE bin_id = ? AND contig_id >= ?
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
                [genome_id, current_contig_id, max_contigs],
            )
        else:
            contig_rows = store.execute(
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
                    WHERE bin_id = ? AND contig_id > ?
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
                [genome_id, after_contig_id, max_contigs],
            )
        if current_contig_id is not None and (
            not contig_rows or str(contig_rows[0][0]) != current_contig_id
        ):
            raise ValueError("Genome packet cursor current contig is absent from this bin")

        selected_contig_ids = [str(row[0]) for row in contig_rows]
        protein_rows: list[tuple[Any, ...]] = []
        if selected_contig_ids:
            placeholders = ", ".join(["?"] * len(selected_contig_ids))
            position = current_position or (0, -1, -1, "")
            protein_rows = store.execute(
                f"""
                SELECT
                    p.bin_id,
                    p.contig_id,
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
                WHERE p.bin_id = ?
                  AND p.contig_id IN ({placeholders})
                  AND (
                        ? IS NULL
                        OR p.contig_id != ?
                        OR (
                            CAST(p.gene_index IS NULL AS INTEGER),
                            COALESCE(CAST(p.gene_index AS BIGINT), ?),
                            p.start,
                            p.protein_id
                        ) > (?, ?, ?, ?)
                      )
                ORDER BY
                    p.contig_id,
                    p.gene_index NULLS LAST,
                    p.start,
                    p.protein_id
                LIMIT ?
                """,
                [
                    _NULL_GENE_INDEX,
                    genome_id,
                    *selected_contig_ids,
                    current_contig_id,
                    current_contig_id,
                    _NULL_GENE_INDEX,
                    *position,
                    max_proteins + 1,
                ],
            )

        protein_ids = [str(row[2]) for row in protein_rows]
        context = load_context_evidence(store, protein_ids)
        predicates = get_active_predicates(store, protein_ids)
        observed = context["observed_annotations"]
        named_calls = context["named_calls"]
        loci = context["loci"]
        proteins_by_contig: dict[
            str,
            list[tuple[dict[str, Any], tuple[int, int, int, str]]],
        ] = {}
        for row in protein_rows:
            (
                protein_bin_id,
                protein_contig_id,
                protein_id,
                gene_index,
                start,
                end_coord,
                strand,
                length_aa,
                gc_content,
                missing_gene_index,
                gene_order,
            ) = row
            if str(protein_bin_id) != genome_id:
                raise RuntimeError("Genome packet query returned a cross-bin protein")
            contig_key = str(protein_contig_id)
            if contig_key not in selected_contig_ids:
                raise RuntimeError("Genome packet query returned an unselected contig")
            protein_key = str(protein_id)
            protein_annotations = observed.get(protein_key, [])
            if not all_annotations:
                protein_annotations = protein_annotations[:1]
            record = {
                "protein_id": protein_key,
                "gene_index": int(gene_index) if gene_index is not None else None,
                "start": int(start),
                "end": int(end_coord),
                "strand": str(strand),
                "length_aa": int(length_aa),
                "gc_content": (
                    float(gc_content) if gc_content is not None else None
                ),
                "observed_annotations": protein_annotations,
                "predicates": predicates.get(protein_key, []),
                "named_calls": named_calls.get(protein_key, []),
                "loci": loci.get(protein_key, []),
            }
            packet_position = (
                int(missing_gene_index),
                int(gene_order),
                int(start),
                protein_key,
            )
            proteins_by_contig.setdefault(contig_key, []).append(
                (record, packet_position)
            )

        frame_contigs: list[dict[str, Any]] = []
        frame_proteins = 0
        state_after_contig = after_contig_id
        state_current_contig = current_contig_id
        state_position = current_position
        state_consumed = consumed_in_current
        state_segment_index = current_segment_index
        state_completed_contigs = completed_contigs
        state_completed_proteins = completed_proteins

        def segment_record(
            *,
            contig_row: tuple[Any, ...],
            records: list[tuple[dict[str, Any], tuple[int, int, int, str]]],
            prefix: int,
            start_offset: int,
            segment_index: int,
        ) -> dict[str, Any]:
            (
                contig_id,
                length,
                contig_gc_content,
                is_circular,
                contig_taxonomy,
                protein_count,
            ) = contig_row
            end_offset = start_offset + prefix
            return {
                "bin_id": genome_id,
                "contig_id": str(contig_id),
                "length": int(length),
                "gc_content": (
                    float(contig_gc_content)
                    if contig_gc_content is not None
                    else None
                ),
                "is_circular": bool(is_circular),
                "taxonomy": contig_taxonomy,
                "total_protein_count": int(protein_count),
                "protein_offset_start": start_offset,
                "protein_offset_end": end_offset,
                "segment_index": segment_index,
                "complete": end_offset == int(protein_count),
                "proteins": [record for record, _position in records[:prefix]],
            }

        # Frame size is computed incrementally rather than by re-serialising the
        # whole accumulated frame once per contig considered.
        #
        # Canonical JSON (sort_keys, separators=(",", ":")) is compositional, so
        #   len(payload with contigs=[c1..cn])
        #     == len(payload with contigs=[]) + sum(len(ci)) + (n - 1)
        # where the (n - 1) accounts for the separating commas. This is
        # byte-identical to the previous full re-serialisation — verified in
        # tests/test_operators/test_contigs.py — so it changes no packing
        # decision and no packing-contract hash.
        #
        # The previous form was O(n^2) in contigs per frame and GIL-bound, which
        # dominated the offline packet census: ~1.79M contig-serialisations for a
        # 1890-contig frame versus 1890 here.
        _frame_envelope_bytes = len(
            _canonical_bytes(
                _genome_model_payload(
                    dataset_id=_dataset_id(store),
                    genome=genome,
                    frame_index=frame_index,
                    contigs=[],
                )
            )
        )
        _accumulated_contig_bytes = [0]

        def candidate_bytes_with(segment: dict[str, Any]) -> int:
            """Serialised size of the frame if `segment` were appended."""
            count = len(frame_contigs) + 1
            return (
                _frame_envelope_bytes
                + _accumulated_contig_bytes[0]
                + len(_canonical_bytes(segment))
                + (count - 1)
            )

        def append_segment(segment: dict[str, Any]) -> None:
            """Append to the frame and keep the running byte total in step.

            Appending only through here is what guarantees the accumulator can
            never drift from `frame_contigs`.
            """
            frame_contigs.append(segment)
            _accumulated_contig_bytes[0] += len(_canonical_bytes(segment))

        for contig_row in contig_rows:
            contig_id = str(contig_row[0])
            total_on_contig = int(contig_row[5])
            is_current = state_current_contig == contig_id
            start_offset = state_consumed if is_current else 0
            segment_index = state_segment_index if is_current else 0
            if start_offset > total_on_contig:
                raise ValueError("Genome packet cursor exceeds current contig total")
            remaining_on_contig = total_on_contig - start_offset
            available = proteins_by_contig.get(contig_id, [])
            whole_loaded = len(available) == remaining_on_contig

            if remaining_on_contig == 0 and total_on_contig > 0:
                raise ValueError("Genome packet cursor points past a completed contig")
            if remaining_on_contig == 0:
                full_segment = segment_record(
                    contig_row=contig_row,
                    records=[],
                    prefix=0,
                    start_offset=0,
                    segment_index=0,
                )
                if (
                    candidate_bytes_with(full_segment) > max_model_payload_bytes
                    and frame_contigs
                ):
                    break
                append_segment(full_segment)
                state_after_contig = contig_id
                state_current_contig = None
                state_position = None
                state_consumed = 0
                state_segment_index = 0
                state_completed_contigs += 1
                continue

            full_segment = (
                segment_record(
                    contig_row=contig_row,
                    records=available,
                    prefix=remaining_on_contig,
                    start_offset=start_offset,
                    segment_index=segment_index,
                )
                if whole_loaded
                else None
            )
            full_fits = (
                full_segment is not None
                and frame_proteins + remaining_on_contig <= max_proteins
                and candidate_bytes_with(full_segment) <= max_model_payload_bytes
            )
            if full_fits:
                append_segment(full_segment)
                frame_proteins += remaining_on_contig
                state_after_contig = contig_id
                state_current_contig = None
                state_position = None
                state_consumed = 0
                state_segment_index = 0
                state_completed_contigs += 1
                state_completed_proteins += remaining_on_contig
                continue

            standalone_fits = False
            if full_segment is not None and remaining_on_contig <= max_proteins:
                standalone_payload = _genome_model_payload(
                    dataset_id=_dataset_id(store),
                    genome=genome,
                    frame_index=frame_index + 1,
                    contigs=[full_segment],
                )
                standalone_fits = (
                    len(_canonical_bytes(standalone_payload))
                    <= max_model_payload_bytes
                )
            if frame_contigs and standalone_fits:
                break
            if (
                frame_contigs
                and not whole_loaded
                and remaining_on_contig <= max_proteins
            ):
                break

            available_budget = max_proteins - frame_proteins
            max_prefix = min(
                available_budget,
                len(available),
                remaining_on_contig,
            )
            if max_prefix < 1:
                break
            low = 1
            high = max_prefix
            best = 0
            while low <= high:
                middle = (low + high) // 2
                segment = segment_record(
                    contig_row=contig_row,
                    records=available,
                    prefix=middle,
                    start_offset=start_offset,
                    segment_index=segment_index,
                )
                if candidate_bytes_with(segment) <= max_model_payload_bytes:
                    best = middle
                    low = middle + 1
                else:
                    high = middle - 1
            if best == 0:
                if frame_contigs:
                    break
                best = 1
            segment = segment_record(
                contig_row=contig_row,
                records=available,
                prefix=best,
                start_offset=start_offset,
                segment_index=segment_index,
            )
            append_segment(segment)
            frame_proteins += best
            state_completed_proteins += best
            state_consumed = start_offset + best
            state_position = available[best - 1][1]
            if segment["complete"]:
                state_after_contig = contig_id
                state_current_contig = None
                state_position = None
                state_consumed = 0
                state_segment_index = 0
                state_completed_contigs += 1
                continue
            state_current_contig = contig_id
            state_segment_index = segment_index + 1
            break

        complete = (
            state_completed_contigs == total_contigs
            and state_completed_proteins == total_proteins
        )
        if (
            state_completed_contigs == total_contigs
        ) != (state_completed_proteins == total_proteins):
            raise RuntimeError(
                "Genome packet traversal reached inconsistent contig/protein totals"
            )
        if not frame_contigs and not complete:
            raise RuntimeError(
                "Genome packet packing made no progress under the requested contract"
            )
        model_payload = _genome_model_payload(
            dataset_id=_dataset_id(store),
            genome=genome,
            frame_index=frame_index,
            contigs=frame_contigs,
        )
        payload_bytes = len(_canonical_bytes(model_payload))
        payload_sha256 = _sha256(model_payload)
        target_exceeded = payload_bytes > max_model_payload_bytes
        next_cursor = None
        if not complete:
            next_cursor = _encode_cursor(
                {
                    "v": 1,
                    "kind": "genome_packet",
                    "dataset_id": _dataset_id(store),
                    "genome_id": genome_id,
                    "packing_contract_hash": packing_hash,
                    "frame_index": frame_index + 1,
                    "after_contig_id": state_after_contig,
                    "current_contig_id": state_current_contig,
                    "position": (
                        list(state_position) if state_position is not None else None
                    ),
                    "consumed_in_current": state_consumed,
                    "current_segment_index": state_segment_index,
                    "completed_contigs": state_completed_contigs,
                    "completed_proteins": state_completed_proteins,
                }
            )
        frame_identity = {
            "packet_version": GENOME_PACKET_VERSION,
            "dataset_id": _dataset_id(store),
            "bin_id": genome_id,
            "packing_contract_hash": packing_hash,
            "frame_index": frame_index,
            "model_payload_sha256": payload_sha256,
        }
        coverage_receipts = [
            {
                "bin_id": genome_id,
                "contig_id": contig["contig_id"],
                "total_protein_count": contig["total_protein_count"],
                "protein_offset_start": contig["protein_offset_start"],
                "protein_offset_end": contig["protein_offset_end"],
                "segment_index": contig["segment_index"],
                "complete": contig["complete"],
            }
            for contig in frame_contigs
        ]
        raw = {
            "packet_version": GENOME_PACKET_VERSION,
            "dataset_id": _dataset_id(store),
            "genome_id": genome_id,
            "bin_id": genome_id,
            "frame_id": f"genome-frame-{_sha256(frame_identity)[:24]}",
            "frame_index": frame_index,
            "cursor": cursor,
            "next_cursor": next_cursor,
            "complete": complete,
            "packing_contract": packing,
            "packing_contract_hash": packing_hash,
            "model_payload": model_payload,
            "model_payload_bytes": payload_bytes,
            "model_payload_sha256": payload_sha256,
            "token_scenarios": _token_scenarios(payload_bytes),
            "target_exceeded": target_exceeded,
            "coverage_receipts": coverage_receipts,
        }
        lines = [
            f"# Genome packet: {genome_id}",
            f"- Packet version: {GENOME_PACKET_VERSION}",
            f"- Frame: {frame_index:,}",
            f"- Contig segments: {len(frame_contigs):,}",
            f"- Proteins: {frame_proteins:,}",
            f"- Model payload: {payload_bytes:,} bytes",
            f"- Complete bin: {'yes' if complete else 'no'}",
            "- Bin isolation: exact",
            "- Sequence payload: compute-only",
        ]
        warnings = []
        if target_exceeded:
            warnings.append(
                "One indivisible protein record exceeds max_model_payload_bytes; "
                "the Atlas census will block launch until the contract is repaired."
            )
        return ctx.make_result(
            data="\n".join(lines),
            rows=frame_proteins,
            total_rows=total_proteins,
            truncated=not complete,
            ref=next_cursor,
            raw=raw,
            index_used="idx_contigs_bin",
            status="ok",
            warnings=warnings,
        )


__all__ = [
    "DEFAULT_GENOME_PACKET_BYTES",
    "DEFAULT_GENOME_PACKET_CONTIGS",
    "DEFAULT_GENOME_PACKET_PROTEINS",
    "GENOME_MODEL_PAYLOAD_SCHEMA",
    "GENOME_PACKET_PACKING_SCHEMA",
    "GENOME_PACKET_VERSION",
    "MAX_GENOME_PACKET_BYTES",
    "MAX_GENOME_PACKET_CONTIGS",
    "MAX_GENOME_PACKET_PROTEINS",
    "PACKET_VERSION",
    "genome_packet_packing_contract",
    "get_contig",
    "get_contig_packet",
    "get_genome_packet",
    "list_contigs",
    "minimum_genome_packet_record_bytes",
]
