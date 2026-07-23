"""Tests for navigation operators (list_*, get_*)."""


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

    def test_returns_sharur_result(self, store):
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

    def test_returns_sharur_result(self, store):
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

    def test_returns_sharur_result(self, store):
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

    def test_annotation_stats_separate_proteins_from_repeat_domain_hits(self, store):
        result = get_genome(store, "bin_003")
        pfam = next(
            stat
            for stat in result.raw["annotation_stats"]
            if stat["source"] == "pfam"
        )

        assert pfam == {"source": "pfam", "domain_hits": 3, "proteins": 1}


class TestGetProtein:
    """Tests for get_protein() operator."""

    def test_returns_sharur_result(self, store):
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

    def test_detailed_view_keeps_sequence_compute_only(self, store):
        store.execute(
            """
            UPDATE proteins
            SET sequence = 'SENSITIVE_SEQUENCE_PAYLOAD'
            WHERE protein_id = 'prot_001'
            """
        )

        result = get_protein(store, "prot_001", verbosity=2)

        assert "SENSITIVE_SEQUENCE_PAYLOAD" not in result.data
        assert "sequence" not in result.raw
        assert "compute-only" in result.data

    def test_uses_persisted_v2_compatibility_predicates(self, store):
        """Protein views must use the same persisted predicate surface as search."""
        store.execute(
            """
            INSERT INTO protein_predicates (protein_id, predicates)
            VALUES ('prot_001', ['v2_only_predicate'])
            """
        )

        result = get_protein(store, "prot_001")

        assert result._raw["predicates"] == ["v2_only_predicate"]
        assert "v2_only_predicate" in result.data

    def test_derives_predicates_from_v2_state_when_cache_row_is_missing(self, store):
        """A missing compatibility row must not force the legacy evaluator."""
        store.execute(
            """
            INSERT INTO semantic_state (
                protein_id, activities, roles, architecture, localization,
                topology, size_class, quality_flags, composite_predicates,
                unresolved_count
            )
            VALUES (
                'prot_001', [], ['v2_state_role'], [], [], '{}', 'standard',
                [], ['v2_composite'], 0
            )
            """
        )

        result = get_protein(store, "prot_001")

        assert "v2_state_role" in result._raw["predicates"]
        assert "v2_composite" in result._raw["predicates"]


class TestGetNeighborhood:
    """Tests for get_neighborhood() operator."""

    def test_returns_sharur_result(self, store):
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

    def test_neighborhood_uses_persisted_predicates(self, store):
        """Neighborhood rows must agree with the persisted V2 search view."""
        store.execute(
            """
            INSERT INTO protein_predicates (protein_id, predicates)
            VALUES ('prot_003', ['v2_neighborhood_predicate'])
            """
        )

        result = get_neighborhood(store, "prot_003", window=1)
        anchor = next(p for p in result._raw["proteins"] if p["is_anchor"])

        assert anchor["predicates"] == ["v2_neighborhood_predicate"]
        assert "v2_neighborhood" in result.data
