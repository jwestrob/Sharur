"""Tests for navigation operators (list_*, get_*)."""

import pytest

from sharur.operators.navigation import (
    list_genomes,
    list_proteins,
    get_genome,
    get_protein,
    get_neighborhood,
)
from sharur.operators.base import SharurResult


class TestListGenomes:
    """Tests for list_genomes() operator."""

    def test_returns_bennu_result(self, store):
        """list_genomes() should return a SharurResult."""
        result = list_genomes(store)
        assert isinstance(result, SharurResult)

    def test_returns_all_genomes(self, store):
        """list_genomes() should return all genomes by default."""
        result = list_genomes(store)
        assert result.meta.rows == 3
        assert result.meta.total_rows == 3

    def test_taxonomy_filter(self, store):
        """list_genomes() should filter by taxonomy."""
        result = list_genomes(store, taxonomy_filter="Archaea")
        assert result.meta.rows == 1
        assert "bin_001" in result.data

    def test_completeness_filter(self, store):
        """list_genomes() should filter by completeness."""
        result = list_genomes(store, min_completeness=90)
        assert result.meta.rows == 1
        assert "bin_001" in result.data

    def test_contamination_filter(self, store):
        """list_genomes() should filter by contamination."""
        result = list_genomes(store, max_contamination=3)
        assert result.meta.rows == 1
        assert "bin_001" in result.data

    def test_limit(self, store):
        """list_genomes() should respect limit."""
        result = list_genomes(store, limit=2)
        assert result.meta.rows == 2
        assert result.meta.total_rows == 3
        assert result.meta.truncated is True

    def test_raw_data(self, store):
        """list_genomes() should include raw data."""
        result = list_genomes(store)
        assert result._raw is not None
        assert len(result._raw) == 3
        assert all("bin_id" in g for g in result._raw)


class TestListProteins:
    """Tests for list_proteins() operator."""

    def test_returns_bennu_result(self, store):
        """list_proteins() should return a SharurResult."""
        result = list_proteins(store)
        assert isinstance(result, SharurResult)

    def test_returns_all_proteins(self, store):
        """list_proteins() should return all proteins by default."""
        result = list_proteins(store)
        assert result.meta.total_rows == 10

    def test_genome_filter(self, store):
        """list_proteins() should filter by genome."""
        result = list_proteins(store, genome_id="bin_001")
        assert result.meta.total_rows == 5

    def test_contig_filter(self, store):
        """list_proteins() should filter by contig."""
        result = list_proteins(store, contig_id="contig_003")
        assert result.meta.total_rows == 3

    def test_length_filter(self, store):
        """list_proteins() should filter by length."""
        # min_length >= 2000 should match prot_004 (2333) and prot_007 (5333)
        result = list_proteins(store, min_length=2000)
        assert result.meta.total_rows == 2

    def test_has_annotation_filter(self, store):
        """list_proteins() should filter by annotation status."""
        # Annotated proteins: prot_001, 002, 003, 006, 008, 009
        result = list_proteins(store, has_annotation=True)
        assert result.meta.total_rows == 6

        # Unannotated: prot_004, 005, 007, 010
        result = list_proteins(store, has_annotation=False)
        assert result.meta.total_rows == 4


class TestGetGenome:
    """Tests for get_genome() operator."""

    def test_returns_bennu_result(self, store):
        """get_genome() should return a SharurResult."""
        result = get_genome(store, "bin_001")
        assert isinstance(result, SharurResult)

    def test_genome_not_found(self, store):
        """get_genome() should handle missing genome."""
        result = get_genome(store, "nonexistent")
        assert "not found" in result.data.lower()
        assert result.meta.rows == 0

    def test_contains_quality_info(self, store):
        """get_genome() should include quality metrics."""
        result = get_genome(store, "bin_001")
        assert "95.5" in result.data or "Completeness" in result.data
        assert "2.1" in result.data or "Contamination" in result.data

    def test_contains_taxonomy(self, store):
        """get_genome() should include taxonomy."""
        result = get_genome(store, "bin_001")
        assert "Archaea" in result.data or "Euryarchaeota" in result.data


class TestGetProtein:
    """Tests for get_protein() operator."""

    def test_returns_bennu_result(self, store):
        """get_protein() should return a SharurResult."""
        result = get_protein(store, "prot_001")
        assert isinstance(result, SharurResult)

    def test_protein_not_found(self, store):
        """get_protein() should handle missing protein."""
        result = get_protein(store, "nonexistent")
        assert "not found" in result.data.lower()
        assert result.meta.rows == 0

    def test_contains_location(self, store):
        """get_protein() should include location info."""
        result = get_protein(store, "prot_001")
        assert "contig_001" in result.data
        assert "1,000" in result.data or "1000" in result.data

    def test_contains_annotations(self, store):
        """get_protein() should include annotations."""
        result = get_protein(store, "prot_001")
        assert "PF00142" in result.data or "NiFe-hydrogenase" in result.data


class TestGetNeighborhood:
    """Tests for get_neighborhood() operator."""

    def test_returns_bennu_result(self, store):
        """get_neighborhood() should return a SharurResult."""
        result = get_neighborhood(store, "prot_003")  # Middle of contig_001
        assert isinstance(result, SharurResult)

    def test_protein_not_found(self, store):
        """get_neighborhood() should handle missing protein."""
        result = get_neighborhood(store, "nonexistent")
        assert "not found" in result.data.lower()
        assert result.meta.rows == 0

    def test_returns_neighbors(self, store):
        """get_neighborhood() should return neighboring proteins."""
        # prot_003 is gene_index 3, should have neighbors
        result = get_neighborhood(store, "prot_003", window=2)
        # Should include prot_001, prot_002, prot_003, prot_004, prot_005
        assert result.meta.rows >= 3

    def test_marks_anchor(self, store):
        """get_neighborhood() should mark anchor protein."""
        result = get_neighborhood(store, "prot_003")
        # Arrow marker for anchor
        assert "→" in result.data or "prot_003" in result.data

    def test_raw_includes_proteins(self, store):
        """get_neighborhood() raw data should include proteins."""
        result = get_neighborhood(store, "prot_003", window=2)
        assert result._raw is not None
        assert "proteins" in result._raw
        assert any(p["is_anchor"] for p in result._raw["proteins"])
