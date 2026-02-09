"""Tests for Sharur facade class."""

import pytest

from sharur.operators import Sharur, SharurResult


class TestSharurFacade:
    """Tests for Sharur facade class."""

    def test_init_with_path(self, tmp_path):
        """Sharur should accept db_path."""
        db_file = tmp_path / "test.duckdb"
        b = Sharur(db_file)
        assert b._db_path == db_file

    def test_init_without_path(self):
        """Sharur should work without db_path (in-memory)."""
        b = Sharur()
        assert b._db_path is None

    def test_overview(self, sharur):
        """Sharur.overview() should return SharurResult."""
        result = sharur.overview()
        assert isinstance(result, SharurResult)
        assert result.meta.rows > 0

    def test_list_genomes(self, sharur):
        """Sharur.list_genomes() should return SharurResult."""
        result = sharur.list_genomes()
        assert isinstance(result, SharurResult)
        assert result.meta.rows == 3

    def test_list_genomes_with_filter(self, sharur):
        """Sharur.list_genomes() should accept filters."""
        result = sharur.list_genomes(taxonomy_filter="Archaea")
        assert result.meta.rows == 1

    def test_list_proteins(self, sharur):
        """Sharur.list_proteins() should return SharurResult."""
        result = sharur.list_proteins()
        assert isinstance(result, SharurResult)
        assert result.meta.total_rows == 10

    def test_list_proteins_with_filter(self, sharur):
        """Sharur.list_proteins() should accept filters."""
        result = sharur.list_proteins(genome_id="bin_001")
        assert result.meta.total_rows == 5

    def test_get_genome(self, sharur):
        """Sharur.get_genome() should return SharurResult."""
        result = sharur.get_genome("bin_001")
        assert isinstance(result, SharurResult)
        assert "bin_001" in result.data

    def test_get_protein(self, sharur):
        """Sharur.get_protein() should return SharurResult."""
        result = sharur.get_protein("prot_001")
        assert isinstance(result, SharurResult)
        assert "prot_001" in result.data

    def test_get_neighborhood(self, sharur):
        """Sharur.get_neighborhood() should return SharurResult."""
        result = sharur.get_neighborhood("prot_003", window=2)
        assert isinstance(result, SharurResult)
        assert result.meta.rows > 0

    def test_search_by_predicates(self, sharur):
        """Sharur.search_by_predicates() should return SharurResult."""
        result = sharur.search_by_predicates(has=["giant"])
        assert isinstance(result, SharurResult)
        assert result.meta.total_rows >= 1

    def test_search_proteins(self, sharur):
        """Sharur.search_proteins() should return SharurResult."""
        result = sharur.search_proteins(annotation_pattern="hydrogenase")
        assert isinstance(result, SharurResult)
        assert result.meta.total_rows >= 1


class TestSharurResultStr:
    """Tests for SharurResult string representation."""

    def test_str_returns_data(self, sharur):
        """str(SharurResult) should return data."""
        result = sharur.overview()
        assert str(result) == result.data

    def test_repr(self, sharur):
        """repr(SharurResult) should be informative."""
        result = sharur.list_genomes()
        repr_str = repr(result)
        assert "SharurResult" in repr_str
        assert "rows=" in repr_str
