"""Tests for search operators (search_by_predicates, search_proteins)."""

import pytest

from sharur.operators.search import search_by_predicates, search_proteins
from sharur.operators.base import SharurResult
from sharur.predicates_v2.persistence import generate_and_persist_v2


class TestSearchByPredicates:
    """Tests for search_by_predicates() operator."""

    def test_returns_sharur_result(self, store):
        """search_by_predicates() should return a SharurResult."""
        result = search_by_predicates(store, has=["giant"])
        assert isinstance(result, SharurResult)

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

    def test_repairs_partial_legacy_cache_from_complete_v2_state(self, store):
        """Partial protein_predicates should be rebuilt from complete V2 state."""
        generate_and_persist_v2(
            store,
            chunk_size=5,
            update_legacy_predicates=True,
            return_states=False,
        )
        store.execute("DELETE FROM protein_predicates WHERE protein_id = 'prot_001'")

        result = search_by_predicates(store, has=["giant"])

        assert result.meta.total_rows == 2
        protein_count = store.execute("SELECT COUNT(*) FROM proteins")[0][0]
        legacy_count = store.execute("SELECT COUNT(*) FROM protein_predicates")[0][0]
        assert legacy_count == protein_count

    def test_repairs_stale_legacy_cache_from_complete_v2_state(self, store):
        """Complete but older protein_predicates should be rebuilt from V2."""
        generate_and_persist_v2(
            store,
            chunk_size=5,
            update_legacy_predicates=True,
            return_states=False,
        )
        store.execute("""
            UPDATE protein_predicates
            SET updated_at = updated_at - INTERVAL 1 HOUR
            WHERE protein_id = 'prot_001'
        """)

        result = search_by_predicates(store, has=["giant"])

        assert result.meta.total_rows == 2
        stale = store.execute("""
            SELECT COUNT(*)
            FROM semantic_state ss
            JOIN protein_predicates pp ON pp.protein_id = ss.protein_id
            WHERE pp.updated_at < ss.updated_at
        """)[0][0]
        assert stale == 0

    def test_read_only_search_does_not_repair_complete_stale_cache(self, store):
        """A read query must not mutate a complete compatibility snapshot."""
        generate_and_persist_v2(
            store,
            chunk_size=5,
            update_legacy_predicates=True,
            return_states=False,
        )
        store.execute("""
            UPDATE protein_predicates
            SET updated_at = updated_at - INTERVAL 1 HOUR
            WHERE protein_id = 'prot_001'
        """)
        store.read_only = True

        result = search_by_predicates(store, has=["giant"])

        assert result.meta.total_rows == 2
        stale = store.execute("""
            SELECT COUNT(*)
            FROM semantic_state ss
            JOIN protein_predicates pp ON pp.protein_id = ss.protein_id
            WHERE pp.updated_at < ss.updated_at
        """)[0][0]
        assert stale == 1

    def test_read_only_partial_cache_fails_without_attempting_repair(self, store):
        """An incomplete read-only snapshot gets an actionable failure."""
        generate_and_persist_v2(
            store,
            chunk_size=5,
            update_legacy_predicates=True,
            return_states=False,
        )
        store.execute("DELETE FROM protein_predicates WHERE protein_id = 'prot_001'")
        store.read_only = True

        with pytest.raises(RuntimeError, match="cannot be repaired in a read-only session"):
            search_by_predicates(store, has=["giant"])

    def test_partial_legacy_cache_without_complete_v2_state_fails(self, store):
        """Partial precomputed predicates should not silently truncate search."""
        store.execute("""
            INSERT INTO protein_predicates (protein_id, predicates, updated_at)
            VALUES ('prot_004', ['giant'], CURRENT_TIMESTAMP)
        """)

        with pytest.raises(RuntimeError, match="compatibility cache is incomplete"):
            search_by_predicates(store, has=["giant"])


class TestSearchProteins:
    """Tests for search_proteins() operator."""

    def test_returns_sharur_result(self, store):
        """search_proteins() should return a SharurResult."""
        result = search_proteins(store)
        assert isinstance(result, SharurResult)

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
