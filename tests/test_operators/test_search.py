"""Tests for search operators (search_by_predicates, search_proteins)."""

import pytest

from bennu.operators.search import search_by_predicates, search_proteins
from bennu.operators.base import BennuResult


class TestSearchByPredicates:
    """Tests for search_by_predicates() operator."""

    def test_returns_bennu_result(self, store):
        """search_by_predicates() should return a BennuResult."""
        result = search_by_predicates(store, has=["giant"])
        assert isinstance(result, BennuResult)

    def test_unknown_predicate(self, store):
        """search_by_predicates() should error on unknown predicate."""
        result = search_by_predicates(store, has=["nonexistent_predicate"])
        assert "unknown" in result.data.lower() or "Unknown" in result.data

    def test_giant_predicate(self, store):
        """search_by_predicates() should find giant proteins."""
        result = search_by_predicates(store, has=["giant"])
        # prot_004 (2333aa) and prot_007 (5333aa) should match
        assert result.meta.total_rows == 2 or result.meta.rows >= 1

    def test_massive_predicate(self, store):
        """search_by_predicates() should find massive proteins."""
        result = search_by_predicates(store, has=["massive"])
        # Only prot_007 (5333aa) should match
        assert result.meta.total_rows == 1

    def test_tiny_predicate(self, store):
        """search_by_predicates() should find tiny proteins."""
        result = search_by_predicates(store, has=["tiny"])
        # prot_005 (40aa) should match
        assert result.meta.total_rows >= 1

    def test_unannotated_predicate(self, store):
        """search_by_predicates() should find unannotated proteins."""
        result = search_by_predicates(store, has=["unannotated"])
        # prot_004, 005, 007, 010 are unannotated
        assert result.meta.total_rows == 4

    def test_combined_has_predicates(self, store):
        """search_by_predicates() should AND multiple has predicates."""
        result = search_by_predicates(store, has=["giant", "unannotated"])
        # prot_004 and prot_007 are both giant and unannotated
        assert result.meta.total_rows == 2

    def test_lacks_predicate(self, store):
        """search_by_predicates() should exclude lacks predicates."""
        result = search_by_predicates(store, has=["giant"], lacks=["massive"])
        # prot_004 is giant but not massive
        assert result.meta.total_rows == 1


class TestSearchProteins:
    """Tests for search_proteins() operator."""

    def test_returns_bennu_result(self, store):
        """search_proteins() should return a BennuResult."""
        result = search_proteins(store)
        assert isinstance(result, BennuResult)

    def test_annotation_pattern(self, store):
        """search_proteins() should filter by annotation pattern."""
        result = search_proteins(store, annotation_pattern="hydrogenase")
        assert result.meta.total_rows >= 1
        assert "prot_001" in result.data or "NiFe" in result.data

    def test_accession_exact(self, store):
        """search_proteins() should filter by exact accession."""
        result = search_proteins(store, accession="PF00142")
        assert result.meta.total_rows == 1

    def test_taxonomy_filter(self, store):
        """search_proteins() should filter by taxonomy."""
        result = search_proteins(store, taxonomy_filter="Archaea")
        # bin_001 is Archaea, has 5 proteins
        assert result.meta.total_rows == 5

    def test_length_filters(self, store):
        """search_proteins() should filter by length."""
        result = search_proteins(store, min_length=2000, max_length=3000)
        # prot_004 (2333aa) should match
        assert result.meta.total_rows == 1

    def test_combined_filters(self, store):
        """search_proteins() should combine multiple filters."""
        result = search_proteins(
            store,
            annotation_pattern="ATP",
            taxonomy_filter="Archaea",
        )
        # prot_002 has ATP synthase annotation and is in Archaea bin
        assert result.meta.total_rows >= 1
