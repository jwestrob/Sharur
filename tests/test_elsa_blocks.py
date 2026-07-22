"""Regression tests for orientation-aware ELSA anchor pairing."""

from __future__ import annotations

import pytest

from sharur.elsa_blocks import (
    ELSABlockFormatError,
    elsa_anchor_pairs_from_block,
    pair_elsa_anchor_gene_ids,
)


def test_forward_block_pairs_equal_positions() -> None:
    pairs = pair_elsa_anchor_gene_ids(
        '["query_1", "query_3"]',
        '["target_2", "target_4"]',
        1,
    )

    assert pairs == (("query_1", "target_2"), ("query_3", "target_4"))


def test_inverted_block_reverses_independently_sorted_target_array() -> None:
    """Regression for DPANN block 1065691 and its formerly crossed pairs."""

    block = {
        "anchor_query_gene_ids": (
            '["JAGVFW010000009.1_28", "JAGVFW010000009.1_30"]'
        ),
        "anchor_target_gene_ids": '["MEW5955379.1", "MEW5955381.1"]',
        "orientation": -1,
    }

    assert elsa_anchor_pairs_from_block(block) == (
        ("JAGVFW010000009.1_28", "MEW5955381.1"),
        ("JAGVFW010000009.1_30", "MEW5955379.1"),
    )


def test_csv_schema_and_text_orientation_are_supported() -> None:
    block = {
        "query_anchor_genes": ["query_low", "query_high"],
        "target_anchor_genes": ["target_low", "target_high"],
        "orientation": "inverted",
    }

    assert elsa_anchor_pairs_from_block(block) == (
        ("query_low", "target_high"),
        ("query_high", "target_low"),
    )


def test_mismatched_anchor_counts_fail_loudly() -> None:
    with pytest.raises(ELSABlockFormatError, match="equal length"):
        pair_elsa_anchor_gene_ids(
            ["query_1", "query_2"],
            ["target_1"],
            "same",
        )


def test_declared_anchor_count_mismatch_fails_loudly() -> None:
    block = {
        "query_anchor_genes": ["query_1", "query_2"],
        "target_anchor_genes": ["target_1", "target_2"],
        "orientation": "same",
        "n_anchors": "3",
    }

    with pytest.raises(ELSABlockFormatError, match="n_anchors disagrees"):
        elsa_anchor_pairs_from_block(block)


@pytest.mark.parametrize("orientation", [0, "unknown", True, None])
def test_unknown_orientation_fails_loudly(orientation: object) -> None:
    with pytest.raises(ELSABlockFormatError, match="orientation"):
        pair_elsa_anchor_gene_ids(["query"], ["target"], orientation)


def test_malformed_anchor_json_fails_loudly() -> None:
    with pytest.raises(ELSABlockFormatError, match="valid JSON"):
        pair_elsa_anchor_gene_ids('["query"', '["target"]', 1)
