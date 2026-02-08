"""Tests for Bennu facade class."""

import pytest

from bennu.operators import Bennu, BennuResult


class TestBennuFacade:
    """Tests for Bennu facade class."""

    def test_init_with_path(self, tmp_path):
        """Bennu should accept db_path."""
        db_file = tmp_path / "test.duckdb"
        b = Bennu(db_file)
        assert b._db_path == db_file

    def test_init_without_path(self):
        """Bennu should work without db_path (in-memory)."""
        b = Bennu()
        assert b._db_path is None

    def test_overview(self, bennu):
        """Bennu.overview() should return BennuResult."""
        result = bennu.overview()
        assert isinstance(result, BennuResult)
        assert result.meta.rows > 0

    def test_list_genomes(self, bennu):
        """Bennu.list_genomes() should return BennuResult."""
        result = bennu.list_genomes()
        assert isinstance(result, BennuResult)
        assert result.meta.rows == 3

    def test_list_genomes_with_filter(self, bennu):
        """Bennu.list_genomes() should accept filters."""
        result = bennu.list_genomes(taxonomy_filter="Archaea")
        assert result.meta.rows == 1

    def test_list_proteins(self, bennu):
        """Bennu.list_proteins() should return BennuResult."""
        result = bennu.list_proteins()
        assert isinstance(result, BennuResult)
        assert result.meta.total_rows == 10

    def test_list_proteins_with_filter(self, bennu):
        """Bennu.list_proteins() should accept filters."""
        result = bennu.list_proteins(genome_id="bin_001")
        assert result.meta.total_rows == 5

    def test_get_genome(self, bennu):
        """Bennu.get_genome() should return BennuResult."""
        result = bennu.get_genome("bin_001")
        assert isinstance(result, BennuResult)
        assert "bin_001" in result.data

    def test_get_protein(self, bennu):
        """Bennu.get_protein() should return BennuResult."""
        result = bennu.get_protein("prot_001")
        assert isinstance(result, BennuResult)
        assert "prot_001" in result.data

    def test_get_neighborhood(self, bennu):
        """Bennu.get_neighborhood() should return BennuResult."""
        result = bennu.get_neighborhood("prot_003", window=2)
        assert isinstance(result, BennuResult)
        assert result.meta.rows > 0

    def test_search_by_predicates(self, bennu):
        """Bennu.search_by_predicates() should return BennuResult."""
        result = bennu.search_by_predicates(has=["giant"])
        assert isinstance(result, BennuResult)
        assert result.meta.total_rows >= 1

    def test_search_proteins(self, bennu):
        """Bennu.search_proteins() should return BennuResult."""
        result = bennu.search_proteins(annotation_pattern="hydrogenase")
        assert isinstance(result, BennuResult)
        assert result.meta.total_rows >= 1


class TestBennuResultStr:
    """Tests for BennuResult string representation."""

    def test_str_returns_data(self, bennu):
        """str(BennuResult) should return data."""
        result = bennu.overview()
        assert str(result) == result.data

    def test_repr(self, bennu):
        """repr(BennuResult) should be informative."""
        result = bennu.list_genomes()
        repr_str = repr(result)
        assert "BennuResult" in repr_str
        assert "rows=" in repr_str
