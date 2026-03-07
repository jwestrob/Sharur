"""Tests for V2 rule loader."""

import pytest

from sharur.predicates_v2.model import ClaimRelation, SemanticFacet
from sharur.predicates_v2.rules import clear_caches, get_facet, get_relation


@pytest.fixture(autouse=True)
def _reset_caches():
    """Clear rule caches before each test."""
    clear_caches()
    yield
    clear_caches()


class TestGetFacet:
    """Tests for facet assignment lookup."""

    def test_activity_predicates(self):
        """Activity predicates should map to activity facet."""
        assert get_facet("hydrogenase") == SemanticFacet.activity
        assert get_facet("kinase") == SemanticFacet.activity
        assert get_facet("protease") == SemanticFacet.activity
        assert get_facet("oxidoreductase") == SemanticFacet.activity

    def test_role_predicates(self):
        """Role predicates should map to role facet."""
        assert get_facet("transporter") == SemanticFacet.role
        assert get_facet("defense_system") == SemanticFacet.role
        assert get_facet("regulator") == SemanticFacet.role
        assert get_facet("chaperone") == SemanticFacet.role

    def test_localization_overrides(self):
        """Transport-category localization predicates should be overridden."""
        assert get_facet("membrane") == SemanticFacet.localization
        assert get_facet("periplasmic") == SemanticFacet.localization
        assert get_facet("outer_membrane") == SemanticFacet.localization
        assert get_facet("inner_membrane") == SemanticFacet.localization
        assert get_facet("cytoplasmic") == SemanticFacet.localization
        assert get_facet("secreted") == SemanticFacet.localization
        assert get_facet("cell_surface") == SemanticFacet.localization

    def test_architecture_predicates(self):
        """Architecture predicates should map to architecture facet."""
        assert get_facet("multi_domain") == SemanticFacet.architecture
        assert get_facet("beta_barrel") == SemanticFacet.architecture
        assert get_facet("coiled_coil") == SemanticFacet.architecture

    def test_architecture_overrides(self):
        """Binding-category structural features should be architecture."""
        assert get_facet("zinc_finger") == SemanticFacet.architecture
        assert get_facet("helix_turn_helix") == SemanticFacet.architecture
        assert get_facet("p_loop") == SemanticFacet.architecture

    def test_size_predicates(self):
        """Size predicates should map to size_class facet."""
        assert get_facet("tiny") == SemanticFacet.size_class
        assert get_facet("giant") == SemanticFacet.size_class
        assert get_facet("massive") == SemanticFacet.size_class

    def test_quality_flag_predicates(self):
        """Annotation status predicates should map to quality_flag."""
        assert get_facet("unannotated") == SemanticFacet.quality_flag
        assert get_facet("well_annotated") == SemanticFacet.quality_flag
        assert get_facet("gc_outlier") == SemanticFacet.quality_flag

    def test_topology_predicates(self):
        """Topology predicates should map to topology facet."""
        assert get_facet("transmembrane_predicted") == SemanticFacet.topology
        assert get_facet("soluble_predicted") == SemanticFacet.topology
        assert get_facet("single_pass_membrane") == SemanticFacet.topology

    def test_unknown_predicate_defaults_to_role(self):
        """Unknown predicates should default to role."""
        assert get_facet("some_unknown_predicate_xyz") == SemanticFacet.role

    def test_envelope_predicates(self):
        """Envelope predicates should map to localization."""
        assert get_facet("cell_wall") == SemanticFacet.localization
        assert get_facet("peptidoglycan") == SemanticFacet.localization
        assert get_facet("s_layer") == SemanticFacet.localization


class TestGetRelation:
    """Tests for relation override lookup."""

    def test_source_defaults(self):
        """Should return correct source-level defaults."""
        assert get_relation("kegg") == ClaimRelation.implies
        assert get_relation("pfam") == ClaimRelation.supports
        assert get_relation("cazy") == ClaimRelation.supports
        assert get_relation("vogdb") == ClaimRelation.flags
        assert get_relation("hyddb") == ClaimRelation.implies
        assert get_relation("defensefinder_system") == ClaimRelation.implies
        assert get_relation("defensefinder") == ClaimRelation.flags

    def test_accession_overrides(self):
        """Accession overrides should take priority over source defaults."""
        # PF00374 is promoted to implies even though pfam default is supports
        assert get_relation("pfam", "PF00374") == ClaimRelation.implies
        # PF02906 is promoted to implies
        assert get_relation("pfam", "PF02906") == ClaimRelation.implies
        # PF04055 (Radical_SAM) is demoted to flags
        assert get_relation("pfam", "PF04055") == ClaimRelation.flags

    def test_accession_without_override_uses_source_default(self):
        """Accession without override should use source default."""
        assert get_relation("pfam", "PF99999") == ClaimRelation.supports
        assert get_relation("kegg", "K99999") == ClaimRelation.implies

    def test_unknown_source_defaults_to_supports(self):
        """Unknown sources should default to supports."""
        assert get_relation("some_unknown_db") == ClaimRelation.supports

    def test_property_source(self):
        """Internal property source should be implies."""
        assert get_relation("_property") == ClaimRelation.implies
        assert get_relation("_computed") == ClaimRelation.implies

    def test_source_name_case_insensitive(self):
        """Source lookup should be case-insensitive."""
        assert get_relation("PFAM") == ClaimRelation.supports
        assert get_relation("Kegg") == ClaimRelation.implies

    def test_vog_alias(self):
        """'vog' should be treated same as 'vogdb'."""
        assert get_relation("vog") == ClaimRelation.flags
