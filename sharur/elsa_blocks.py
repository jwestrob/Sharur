"""Orientation-aware parsing of ELSA synteny block anchors.

ELSA serializes the query and target anchor arrays independently in ascending
genomic-coordinate order.  Consequently, equal list positions are homologous
anchor pairs only for forward blocks.  For an inverted block, query position
``i`` pairs with target position ``n - 1 - i``.

Use :func:`elsa_anchor_pairs_from_block` whenever pair identity matters.  Raw
anchor arrays remain suitable for side-specific membership tests, but must not
be zipped directly.
"""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from typing import Any


class ELSABlockFormatError(ValueError):
    """Raised when an ELSA block cannot be paired without ambiguity."""


_GENE_ID_KEY_PAIRS = (
    ("anchor_query_gene_ids", "anchor_target_gene_ids"),
    ("query_anchor_genes", "target_anchor_genes"),
)

_FORWARD_ORIENTATIONS = {"1", "+1", "same", "forward", "fwd"}
_INVERTED_ORIENTATIONS = {"-1", "inverted", "inverse", "reverse", "rev"}


def parse_elsa_anchor_gene_ids(
    value: str | Sequence[str],
    *,
    field_name: str = "anchor_gene_ids",
) -> tuple[str, ...]:
    """Parse one ELSA JSON/list field into a validated tuple of protein IDs."""

    parsed: Any = value
    if isinstance(value, str):
        try:
            parsed = json.loads(value)
        except json.JSONDecodeError as exc:
            raise ELSABlockFormatError(
                f"{field_name} is not valid JSON: {exc.msg}"
            ) from exc

    if isinstance(parsed, (str, bytes)) or not isinstance(parsed, Sequence):
        raise ELSABlockFormatError(
            f"{field_name} must be a JSON array or sequence of protein IDs"
        )

    gene_ids: list[str] = []
    for index, gene_id in enumerate(parsed):
        if not isinstance(gene_id, str) or not gene_id:
            raise ELSABlockFormatError(
                f"{field_name}[{index}] must be a non-empty string"
            )
        gene_ids.append(gene_id)
    return tuple(gene_ids)


def normalize_elsa_orientation(value: object) -> int:
    """Normalize ELSA integer or textual orientation to ``1`` or ``-1``."""

    if isinstance(value, bool):
        raise ELSABlockFormatError("ELSA orientation cannot be boolean")
    if value == 1:
        return 1
    if value == -1:
        return -1
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in _FORWARD_ORIENTATIONS:
            return 1
        if normalized in _INVERTED_ORIENTATIONS:
            return -1
    raise ELSABlockFormatError(
        "ELSA orientation must be 1/-1 or same/inverted; "
        f"received {value!r}"
    )


def pair_elsa_anchor_gene_ids(
    query_gene_ids: str | Sequence[str],
    target_gene_ids: str | Sequence[str],
    orientation: object,
) -> tuple[tuple[str, str], ...]:
    """Return query-ascending homologous anchor pairs for one ELSA block.

    Both serialized arrays are assumed to follow the ELSA output contract:
    independently sorted in ascending genomic-coordinate order.  Target IDs
    are therefore reversed for inverted blocks before pairing.
    """

    query = parse_elsa_anchor_gene_ids(
        query_gene_ids,
        field_name="query anchor gene IDs",
    )
    target = parse_elsa_anchor_gene_ids(
        target_gene_ids,
        field_name="target anchor gene IDs",
    )
    if len(query) != len(target):
        raise ELSABlockFormatError(
            "ELSA query and target anchor arrays must have equal length; "
            f"received {len(query)} and {len(target)}"
        )

    if normalize_elsa_orientation(orientation) == -1:
        target = tuple(reversed(target))
    return tuple(zip(query, target, strict=True))


def elsa_anchor_pairs_from_block(
    block: Mapping[str, Any],
    *,
    query_key: str | None = None,
    target_key: str | None = None,
    orientation_key: str = "orientation",
) -> tuple[tuple[str, str], ...]:
    """Parse and pair anchors from an ELSA CSV/parquet block mapping.

    When keys are omitted, both supported ELSA schemas are recognized:
    ``anchor_query_gene_ids``/``anchor_target_gene_ids`` and
    ``query_anchor_genes``/``target_anchor_genes``.
    """

    if (query_key is None) != (target_key is None):
        raise ELSABlockFormatError(
            "query_key and target_key must be supplied together"
        )
    if query_key is None:
        for candidate_query, candidate_target in _GENE_ID_KEY_PAIRS:
            if candidate_query in block and candidate_target in block:
                query_key = candidate_query
                target_key = candidate_target
                break
    if query_key is None or target_key is None:
        expected = " or ".join(f"{query}/{target}" for query, target in _GENE_ID_KEY_PAIRS)
        raise ELSABlockFormatError(
            f"ELSA block is missing an anchor field pair; expected {expected}"
        )
    if orientation_key not in block:
        raise ELSABlockFormatError(
            f"ELSA block is missing orientation field {orientation_key!r}"
        )

    pairs = pair_elsa_anchor_gene_ids(
        block[query_key],
        block[target_key],
        block[orientation_key],
    )
    if "n_anchors" in block:
        try:
            expected_count = int(block["n_anchors"])
        except (TypeError, ValueError) as exc:
            raise ELSABlockFormatError(
                f"ELSA n_anchors must be an integer; received {block['n_anchors']!r}"
            ) from exc
        if expected_count != len(pairs):
            raise ELSABlockFormatError(
                "ELSA n_anchors disagrees with serialized anchor arrays; "
                f"received {expected_count} and {len(pairs)}"
            )
    return pairs


__all__ = [
    "ELSABlockFormatError",
    "elsa_anchor_pairs_from_block",
    "normalize_elsa_orientation",
    "pair_elsa_anchor_gene_ids",
    "parse_elsa_anchor_gene_ids",
]
