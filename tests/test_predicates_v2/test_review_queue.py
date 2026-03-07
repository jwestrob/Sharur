"""Tests for V2 review queue."""

import pytest

from sharur.predicates_v2.model import ClaimRelation, SemanticAtom, SemanticFacet
from sharur.predicates_v2.review_queue import (
    build_review_queue,
    format_review_queue_tsv,
)


def _unresolved_atom(
    accession: str,
    protein_id: str = "prot_001",
    source_db: str = "pfam",
) -> SemanticAtom:
    """Helper to create unresolved atoms."""
    return SemanticAtom(
        protein_id=protein_id,
        atom_id=f"unresolved:{accession}",
        facet=SemanticFacet.quality_flag,
        relation=ClaimRelation.unresolved,
        source_accession=accession,
        source_db=source_db,
    )


def _resolved_atom(
    atom_id: str = "hydrogenase",
    protein_id: str = "prot_001",
) -> SemanticAtom:
    """Helper to create resolved atoms."""
    return SemanticAtom(
        protein_id=protein_id,
        atom_id=atom_id,
        facet=SemanticFacet.activity,
        relation=ClaimRelation.implies,
        source_accession="K00532",
        source_db="kegg",
    )


class TestBuildReviewQueue:
    """Tests for building the review queue."""

    def test_empty_atoms(self):
        """No atoms should produce empty queue."""
        queue = build_review_queue([])
        assert queue == []

    def test_no_unresolved(self):
        """Only resolved atoms should produce empty queue."""
        atoms = [_resolved_atom(), _resolved_atom(protein_id="prot_002")]
        queue = build_review_queue(atoms)
        assert queue == []

    def test_single_unresolved(self):
        """Should surface a single unresolved accession."""
        atoms = [_unresolved_atom("PF12345")]
        queue = build_review_queue(atoms, total_proteins=100)
        assert len(queue) == 1
        assert queue[0]["accession"] == "PF12345"
        assert queue[0]["source_db"] == "pfam"
        assert queue[0]["n_proteins"] == 1
        assert queue[0]["pct_proteome"] == 1.0

    def test_multiple_proteins_same_accession(self):
        """Should aggregate by accession across proteins."""
        atoms = [
            _unresolved_atom("PF12345", protein_id="prot_001"),
            _unresolved_atom("PF12345", protein_id="prot_002"),
            _unresolved_atom("PF12345", protein_id="prot_003"),
        ]
        queue = build_review_queue(atoms, total_proteins=100)
        assert len(queue) == 1
        assert queue[0]["n_proteins"] == 3
        assert queue[0]["pct_proteome"] == 3.0

    def test_priority_ordering(self):
        """Should sort by n_proteins * n_genomes descending."""
        atoms = [
            _unresolved_atom("PF_LOW", protein_id="prot_001"),
            _unresolved_atom("PF_HIGH", protein_id="prot_001"),
            _unresolved_atom("PF_HIGH", protein_id="prot_002"),
            _unresolved_atom("PF_HIGH", protein_id="prot_003"),
        ]
        queue = build_review_queue(atoms, total_proteins=100)
        assert len(queue) == 2
        assert queue[0]["accession"] == "PF_HIGH"
        assert queue[1]["accession"] == "PF_LOW"

    def test_genome_counting(self):
        """Should count genomes when mapping provided."""
        atoms = [
            _unresolved_atom("PF12345", protein_id="prot_001"),
            _unresolved_atom("PF12345", protein_id="prot_002"),
        ]
        protein_genomes = {
            "prot_001": "genome_A",
            "prot_002": "genome_B",
        }
        queue = build_review_queue(
            atoms, protein_genomes=protein_genomes, total_proteins=100,
        )
        assert queue[0]["n_genomes"] == 2

    def test_example_protein_ids(self):
        """Should include up to 5 example protein IDs."""
        atoms = [
            _unresolved_atom("PF12345", protein_id=f"prot_{i:03d}")
            for i in range(10)
        ]
        queue = build_review_queue(atoms, total_proteins=100)
        examples = queue[0]["example_protein_ids"].split(";")
        assert len(examples) == 5

    def test_suggested_facet(self):
        """Should suggest facet based on accession keywords."""
        atoms = [
            _unresolved_atom("some_kinase_domain"),
            _unresolved_atom("membrane_protein"),
            _unresolved_atom("PF12345"),  # No keyword
        ]
        queue = build_review_queue(atoms, total_proteins=100)
        # Build dict for easier lookup
        by_acc = {r["accession"]: r for r in queue}
        assert by_acc["some_kinase_domain"]["suggested_facet"] == "activity"
        assert by_acc["membrane_protein"]["suggested_facet"] == "localization"
        assert by_acc["PF12345"]["suggested_facet"] == ""


class TestFormatReviewQueueTsv:
    """Tests for TSV formatting."""

    def test_empty_queue(self):
        """Should produce header-only TSV for empty queue."""
        tsv = format_review_queue_tsv([])
        lines = tsv.strip().split("\n")
        assert len(lines) == 1  # Header only
        assert "accession" in lines[0]

    def test_format_with_data(self):
        """Should produce properly formatted TSV."""
        rows = [
            {
                "accession": "PF12345",
                "source_db": "pfam",
                "n_proteins": 42,
                "n_genomes": 7,
                "pct_proteome": 1.5,
                "example_protein_ids": "prot_001;prot_002",
                "suggested_facet": "activity",
            }
        ]
        tsv = format_review_queue_tsv(rows)
        lines = tsv.strip().split("\n")
        assert len(lines) == 2  # Header + 1 row
        fields = lines[1].split("\t")
        assert fields[0] == "PF12345"
        assert fields[2] == "42"
