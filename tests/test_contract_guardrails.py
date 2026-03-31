"""
Guardrail tests to prevent schema/API drift.

These tests verify that documented interfaces match the actual codebase,
catching drift before it reaches users or agent prompts.
"""

import inspect

import pytest


class TestFacadeMethodsExist:
    """Verify that all documented Sharur facade methods actually exist."""

    EXPECTED_METHODS = [
        # Navigation
        "overview",
        "describe_schema",
        "list_genomes",
        "list_proteins",
        "get_genome",
        "get_protein",
        "get_neighborhood",
        # Search
        "search_by_predicates",
        "search_proteins",
        # Export
        "export_fasta",
        "export_neighborhood_fasta",
        "get_sequence",
        # Similarity
        "find_similar",
        "find_similar_to_set",
        # Visualization
        "visualize_neighborhood",
        "visualize_domains",
        "get_pathway_context",
        # Structure
        "predict_structure",
        "batch_predict_structures",
        "search_foldseek",
        "search_foldseek_for_protein",
        "list_foldseek_databases",
        # Validation
        "validate_annotation",
        "validate_context",
        "analyze_crispr_systems",
        "detect_annotation_errors",
        # Hypothesis
        "propose_hypothesis",
        "add_evidence",
        "list_hypotheses",
        "hypothesis_summary",
        # Provenance
        "render_provenance",
        "log_provenance",
    ]

    def test_all_documented_methods_exist(self):
        from sharur.operators import Sharur

        missing = []
        for method_name in self.EXPECTED_METHODS:
            if not hasattr(Sharur, method_name):
                missing.append(method_name)
        assert not missing, f"Facade missing documented methods: {missing}"

    def test_no_stale_methods_documented(self):
        """Methods that should NOT exist (catch stale docs)."""
        from sharur.operators import Sharur

        MUST_NOT_EXIST = ["search", "total_proteins", "list_predicates", "regenerate_predicates", "ask"]
        for method_name in MUST_NOT_EXIST:
            assert not hasattr(Sharur, method_name), (
                f"Sharur.{method_name} exists but should not — "
                f"if this was intentionally added, update this test"
            )


class TestAnnotationSourceVocabulary:
    """Verify annotation source enums cover real database values."""

    CANONICAL_SOURCES = {
        "pfam", "kegg", "cazy", "hyddb", "defensefinder", "vogdb",
        "cant_hyd", "txsscan", "defensefinder_system", "txsscan_system",
        "pfam_relaxed", "foldseek", "crispr", "custom",
    }

    def test_models_enum_covers_canonical(self):
        from sharur.core.models import AnnotationSource

        enum_values = {e.value for e in AnnotationSource}
        missing = self.CANONICAL_SOURCES - enum_values
        assert not missing, f"AnnotationSource enum missing: {missing}"

    def test_types_annotation_accepts_all_sources(self):
        """types.Annotation.source should accept any string (not a restrictive Literal)."""
        from sharur.core.types import Annotation

        for source in self.CANONICAL_SOURCES:
            a = Annotation(
                protein_id="test",
                source=source,
                accession="TEST001",
            )
            assert a.source == source


class TestLocusTypeVocabulary:
    """Verify locus type enums cover real database values."""

    CANONICAL_LOCUS_TYPES = {
        "crispr", "bgc", "prophage", "island", "defense_island",
        "viral_contig", "operon", "metabolic", "transposon", "custom",
    }

    def test_types_enum_covers_canonical(self):
        from sharur.core.types import LocusType

        enum_values = {e.value for e in LocusType}
        missing = self.CANONICAL_LOCUS_TYPES - enum_values
        assert not missing, f"types.LocusType enum missing: {missing}"

    def test_models_enum_covers_canonical(self):
        from sharur.core.models import LocusType

        enum_values = {e.value for e in LocusType}
        missing = self.CANONICAL_LOCUS_TYPES - enum_values
        assert not missing, f"models.LocusType enum missing: {missing}"


class TestDuckDBStoreContract:
    """Verify DuckDBStore.execute() returns list[tuple], not a cursor."""

    def test_execute_returns_list(self):
        from sharur.storage.duckdb_store import DuckDBStore

        store = DuckDBStore()  # in-memory
        result = store.execute("SELECT 1 AS x, 2 AS y")
        assert isinstance(result, list), f"execute() returned {type(result)}, expected list"
        assert len(result) == 1
        assert isinstance(result[0], tuple)
        assert result[0] == (1, 2)

    def test_execute_empty_result(self):
        from sharur.storage.duckdb_store import DuckDBStore

        store = DuckDBStore()
        result = store.execute("SELECT 1 WHERE 1 = 0")
        assert isinstance(result, list)
        assert len(result) == 0


class TestFAISSH5Path:
    """Verify session loads FAISS from canonical H5 path."""

    def test_session_uses_canonical_h5(self):
        """session.py should reference protein_embeddings.h5."""
        import sharur.core.session as session_mod

        source = inspect.getsource(session_mod)
        assert "protein_embeddings.h5" in source


class TestSchemaColumnNames:
    """Verify schema uses the correct column names."""

    def test_proteins_uses_end_coord(self):
        from sharur.storage.schema import SCHEMA

        assert "end_coord" in SCHEMA
        # 'end' as a standalone column name should not appear (it's a reserved word)

    def test_proteins_uses_sequence_length(self):
        from sharur.storage.schema import SCHEMA

        assert "sequence_length" in SCHEMA

    def test_annotations_uses_start_aa_end_aa(self):
        from sharur.storage.schema import SCHEMA

        assert "start_aa" in SCHEMA
        assert "end_aa" in SCHEMA
