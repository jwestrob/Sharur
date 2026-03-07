"""Tests for V2 compatibility layer."""

import pytest

from sharur.predicates_v2.compat import semantic_state_to_predicates
from sharur.predicates_v2.model import SemanticState


class TestSemanticStateToPredicates:
    """Tests for V2 -> V1 predicate conversion."""

    def test_empty_state(self):
        """Empty state should produce empty predicate list."""
        state = SemanticState(protein_id="prot_001")
        preds = semantic_state_to_predicates(state)
        assert preds == []

    def test_activities(self):
        """Activities should appear in flat list."""
        state = SemanticState(
            protein_id="prot_001",
            activities=["hydrogenase", "oxidoreductase"],
        )
        preds = semantic_state_to_predicates(state)
        assert "hydrogenase" in preds
        assert "oxidoreductase" in preds

    def test_roles(self):
        """Roles should appear in flat list."""
        state = SemanticState(
            protein_id="prot_001",
            roles=["defense_system", "transporter"],
        )
        preds = semantic_state_to_predicates(state)
        assert "defense_system" in preds
        assert "transporter" in preds

    def test_architecture(self):
        """Architecture atoms should appear in flat list."""
        state = SemanticState(
            protein_id="prot_001",
            architecture=["multi_domain", "beta_barrel"],
        )
        preds = semantic_state_to_predicates(state)
        assert "multi_domain" in preds
        assert "beta_barrel" in preds

    def test_localization(self):
        """Localization atoms should appear in flat list."""
        state = SemanticState(
            protein_id="prot_001",
            localization=["membrane", "outer_membrane"],
        )
        preds = semantic_state_to_predicates(state)
        assert "membrane" in preds
        assert "outer_membrane" in preds

    def test_size_class(self):
        """Size class should appear in flat list."""
        state = SemanticState(
            protein_id="prot_001",
            size_class="giant",
        )
        preds = semantic_state_to_predicates(state)
        assert "giant" in preds

    def test_empty_size_class(self):
        """Empty size class should not add anything."""
        state = SemanticState(
            protein_id="prot_001",
            size_class="",
        )
        preds = semantic_state_to_predicates(state)
        assert "" not in preds

    def test_quality_flags(self):
        """Quality flags should appear in flat list."""
        state = SemanticState(
            protein_id="prot_001",
            quality_flags=["well_annotated", "multi_source"],
        )
        preds = semantic_state_to_predicates(state)
        assert "well_annotated" in preds
        assert "multi_source" in preds

    def test_composite_predicates(self):
        """Composite predicates should appear in flat list."""
        state = SemanticState(
            protein_id="prot_001",
            composite_predicates=["giant_unannotated", "novel_membrane_protein"],
        )
        preds = semantic_state_to_predicates(state)
        assert "giant_unannotated" in preds
        assert "novel_membrane_protein" in preds

    def test_topology_transmembrane(self):
        """Transmembrane topology should add transmembrane_predicted."""
        state = SemanticState(
            protein_id="prot_001",
            topology={"tm_count": 7, "signal_peptide": False},
        )
        preds = semantic_state_to_predicates(state)
        assert "transmembrane_predicted" in preds

    def test_topology_no_tm(self):
        """Non-transmembrane topology should not add transmembrane_predicted."""
        state = SemanticState(
            protein_id="prot_001",
            topology={"tm_count": 0, "signal_peptide": False},
        )
        preds = semantic_state_to_predicates(state)
        assert "transmembrane_predicted" not in preds

    def test_topology_signal_peptide(self):
        """Signal peptide should add signal_peptide predicate."""
        state = SemanticState(
            protein_id="prot_001",
            topology={"tm_count": 0, "signal_peptide": True},
        )
        preds = semantic_state_to_predicates(state)
        assert "signal_peptide" in preds

    def test_sorted_output(self):
        """Output should be sorted."""
        state = SemanticState(
            protein_id="prot_001",
            activities=["kinase", "atpase"],
            roles=["regulator"],
            size_class="large",
        )
        preds = semantic_state_to_predicates(state)
        assert preds == sorted(preds)

    def test_no_duplicates(self):
        """Output should not contain duplicates."""
        state = SemanticState(
            protein_id="prot_001",
            activities=["kinase"],
            roles=["kinase"],  # Same predicate in different facets
        )
        preds = semantic_state_to_predicates(state)
        assert preds.count("kinase") == 1

    def test_full_state(self):
        """Should convert a fully populated state."""
        state = SemanticState(
            protein_id="prot_001",
            activities=["hydrogenase", "oxidoreductase"],
            roles=["defense_system"],
            architecture=["multi_domain"],
            localization=["membrane"],
            topology={"tm_count": 3, "signal_peptide": False},
            size_class="large",
            quality_flags=["well_annotated"],
            composite_predicates=["nife_hydrogenase_validated"],
        )
        preds = semantic_state_to_predicates(state)
        assert "hydrogenase" in preds
        assert "defense_system" in preds
        assert "multi_domain" in preds
        assert "membrane" in preds
        assert "transmembrane_predicted" in preds
        assert "large" in preds
        assert "well_annotated" in preds
        assert "nife_hydrogenase_validated" in preds

    def test_unresolved_not_in_output(self):
        """Unresolved accessions should NOT appear in flat predicates."""
        state = SemanticState(
            protein_id="prot_001",
            unresolved_accessions=["PF12345", "PF67890"],
        )
        preds = semantic_state_to_predicates(state)
        assert "PF12345" not in preds
        assert "PF67890" not in preds
