from __future__ import annotations

import pytest

from sharur.operators.export import export_neighborhood_fasta
from sharur.operators.navigation import get_neighborhood
from sharur.operators.validation import validate_context


@pytest.fixture
def colliding_store(store):
    """Add a legacy-style contig label reused by a different genome."""
    store.execute(
        """
        INSERT INTO proteins (
            protein_id,
            contig_id,
            bin_id,
            start,
            end_coord,
            strand,
            gene_index,
            sequence,
            sequence_length,
            gc_content
        )
        VALUES (
            'prot_foreign',
            'contig_001',
            'bin_002',
            5200,
            5450,
            '+',
            4,
            'MMMM',
            4,
            0.5
        )
        """
    )
    store.execute(
        """
        INSERT INTO annotations (
            annotation_id,
            protein_id,
            source,
            accession,
            name,
            description,
            evalue,
            score
        )
        VALUES (
            1000,
            'prot_foreign',
            'pfam',
            'PF99999',
            'foreign_marker',
            'Protein from a different bin with a reused contig label',
            1e-20,
            80
        )
        """
    )
    return store


def test_navigation_neighborhood_is_scoped_by_bin(colliding_store) -> None:
    result = get_neighborhood(colliding_store, "prot_003", window=10)
    protein_ids = {row["protein_id"] for row in result._raw["proteins"]}

    assert "prot_foreign" not in protein_ids
    assert protein_ids == {
        "prot_001",
        "prot_002",
        "prot_003",
        "prot_004",
        "prot_005",
    }


def test_neighborhood_fasta_is_scoped_by_bin(colliding_store) -> None:
    result = export_neighborhood_fasta(
        colliding_store,
        "prot_003",
        window=10,
    )

    assert "prot_foreign" not in result._raw["fasta"]
    assert result._raw["bin_id"] == "bin_001"


def test_context_validation_is_scoped_by_bin(colliding_store) -> None:
    result = validate_context(
        colliding_store,
        "prot_003",
        expected_neighbors=["foreign_marker"],
        window=10,
    )

    assert result._raw["neighbors_found"] == []
    assert result._raw["neighbors_missing"] == ["foreign_marker"]


def test_coordinate_window_fails_closed_for_ambiguous_contig(
    colliding_store,
) -> None:
    with pytest.raises(ValueError, match="multiple bins"):
        colliding_store.get_proteins_in_window(
            "contig_001",
            0,
            20_000,
        )

    proteins = colliding_store.get_proteins_in_window(
        "contig_001",
        0,
        20_000,
        bin_id="bin_001",
    )
    assert {protein.protein_id for protein in proteins} == {
        "prot_001",
        "prot_002",
        "prot_003",
        "prot_004",
        "prot_005",
    }
