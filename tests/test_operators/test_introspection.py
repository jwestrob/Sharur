"""Tests for introspection operators (overview, describe_schema)."""

import pytest

from bennu.operators.introspection import overview, describe_schema
from bennu.operators.base import BennuResult


class TestOverview:
    """Tests for overview() operator."""

    def test_returns_bennu_result(self, store):
        """overview() should return a BennuResult."""
        result = overview(store)
        assert isinstance(result, BennuResult)

    def test_contains_genome_count(self, store):
        """overview() should include genome count."""
        result = overview(store)
        assert "3" in result.data  # We have 3 bins
        assert "Genomes" in result.data or "genome" in result.data.lower()

    def test_contains_protein_count(self, store):
        """overview() should include protein count."""
        result = overview(store)
        assert "10" in result.data  # We have 10 proteins
        assert "Protein" in result.data or "protein" in result.data.lower()

    def test_contains_taxonomy_distribution(self, store):
        """overview() should include taxonomy info."""
        result = overview(store)
        # Should mention Archaea or Bacteria
        assert "Archaea" in result.data or "Bacteria" in result.data or "rchaeota" in result.data

    def test_meta_has_timing(self, store):
        """Result meta should include timing info."""
        result = overview(store)
        assert result.meta.time_ms >= 0

    def test_raw_data_contains_stats(self, store):
        """Raw data should contain stats dict."""
        result = overview(store)
        assert result._raw is not None
        assert "genome_count" in result._raw
        assert result._raw["genome_count"] == 3


class TestDescribeSchema:
    """Tests for describe_schema() operator."""

    def test_returns_bennu_result(self, store):
        """describe_schema() should return a BennuResult."""
        result = describe_schema(store)
        assert isinstance(result, BennuResult)

    def test_lists_tables(self, store):
        """describe_schema() should list database tables."""
        result = describe_schema(store)
        assert "proteins" in result.data.lower()
        assert "annotations" in result.data.lower()
        assert "bins" in result.data.lower()

    def test_lists_predicates(self, store):
        """describe_schema() should list available predicates."""
        result = describe_schema(store)
        assert "giant" in result.data.lower()
        assert "unannotated" in result.data.lower()
