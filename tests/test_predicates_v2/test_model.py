"""Tests for predicate V2 data model."""

import pytest

from sharur.predicates_v2.model import (
    ClaimRelation,
    SemanticAtom,
    SemanticFacet,
    SemanticState,
)


class TestSemanticFacet:
    """Tests for SemanticFacet enum."""

    def test_all_facets_exist(self):
        """Should have all 7 facet types."""
        expected = {
            "activity",
            "role",
            "architecture",
            "localization",
            "topology",
            "size_class",
            "quality_flag",
        }
        actual = {f.value for f in SemanticFacet}
        assert actual == expected

    def test_facet_is_str_enum(self):
        """Facet values should be usable as strings."""
        assert SemanticFacet.activity == "activity"
        assert str(SemanticFacet.role) == "SemanticFacet.role"
        assert SemanticFacet.activity.value == "activity"

    def test_facet_from_string(self):
        """Should construct facet from string value."""
        f = SemanticFacet("activity")
        assert f is SemanticFacet.activity

    def test_invalid_facet_raises(self):
        """Should raise ValueError for invalid facet."""
        with pytest.raises(ValueError):
            SemanticFacet("nonexistent")


class TestClaimRelation:
    """Tests for ClaimRelation enum."""

    def test_all_relations_exist(self):
        """Should have all 5 relation types."""
        expected = {"implies", "supports", "flags", "excludes", "unresolved"}
        actual = {r.value for r in ClaimRelation}
        assert actual == expected

    def test_relation_is_str_enum(self):
        """Relation values should be usable as strings."""
        assert ClaimRelation.implies == "implies"
        assert ClaimRelation.implies.value == "implies"

    def test_relation_from_string(self):
        """Should construct relation from string value."""
        r = ClaimRelation("excludes")
        assert r is ClaimRelation.excludes


class TestSemanticAtom:
    """Tests for SemanticAtom dataclass."""

    def test_basic_construction(self):
        """Should construct an atom with all fields."""
        atom = SemanticAtom(
            protein_id="prot_001",
            atom_id="hydrogenase",
            facet=SemanticFacet.activity,
            relation=ClaimRelation.implies,
            source_accession="K00532",
            source_db="kegg",
            evidence_evalue=1e-50,
            evidence_score=200.5,
        )
        assert atom.protein_id == "prot_001"
        assert atom.atom_id == "hydrogenase"
        assert atom.facet == SemanticFacet.activity
        assert atom.relation == ClaimRelation.implies
        assert atom.source_accession == "K00532"
        assert atom.source_db == "kegg"
        assert atom.evidence_evalue == 1e-50
        assert atom.evidence_score == 200.5

    def test_optional_evidence_fields(self):
        """Evidence fields should default to None."""
        atom = SemanticAtom(
            protein_id="prot_001",
            atom_id="giant",
            facet=SemanticFacet.size_class,
            relation=ClaimRelation.implies,
            source_accession="_computed",
            source_db="_property",
        )
        assert atom.evidence_evalue is None
        assert atom.evidence_score is None

    def test_to_dict(self):
        """Should serialize to dictionary."""
        atom = SemanticAtom(
            protein_id="prot_001",
            atom_id="hydrogenase",
            facet=SemanticFacet.activity,
            relation=ClaimRelation.implies,
            source_accession="K00532",
            source_db="kegg",
            evidence_evalue=1e-50,
            evidence_score=200.5,
        )
        d = atom.to_dict()
        assert d["protein_id"] == "prot_001"
        assert d["atom_id"] == "hydrogenase"
        assert d["facet"] == "activity"
        assert d["relation"] == "implies"
        assert d["source_accession"] == "K00532"
        assert d["source_db"] == "kegg"
        assert d["evidence_evalue"] == 1e-50
        assert d["evidence_score"] == 200.5

    def test_from_dict(self):
        """Should deserialize from dictionary."""
        d = {
            "protein_id": "prot_001",
            "atom_id": "hydrogenase",
            "facet": "activity",
            "relation": "implies",
            "source_accession": "K00532",
            "source_db": "kegg",
            "evidence_evalue": 1e-50,
            "evidence_score": 200.5,
        }
        atom = SemanticAtom.from_dict(d)
        assert atom.protein_id == "prot_001"
        assert atom.facet == SemanticFacet.activity
        assert atom.relation == ClaimRelation.implies

    def test_roundtrip(self):
        """to_dict -> from_dict should produce equal atom."""
        original = SemanticAtom(
            protein_id="prot_002",
            atom_id="membrane",
            facet=SemanticFacet.localization,
            relation=ClaimRelation.supports,
            source_accession="PF07690",
            source_db="pfam",
            evidence_evalue=1e-20,
            evidence_score=100.0,
        )
        reconstructed = SemanticAtom.from_dict(original.to_dict())
        assert original == reconstructed


class TestSemanticState:
    """Tests for SemanticState dataclass."""

    def test_empty_state(self):
        """Should create empty state with defaults."""
        state = SemanticState(protein_id="prot_001")
        assert state.protein_id == "prot_001"
        assert state.activities == []
        assert state.roles == []
        assert state.architecture == []
        assert state.localization == []
        assert state.topology == {}
        assert state.size_class == ""
        assert state.quality_flags == []
        assert state.composite_predicates == []
        assert state.unresolved_accessions == []

    def test_populated_state(self):
        """Should create state with populated fields."""
        state = SemanticState(
            protein_id="prot_001",
            activities=["hydrogenase", "oxidoreductase"],
            roles=["defense_system"],
            architecture=["multi_domain"],
            localization=["membrane"],
            topology={"tm_count": 7, "signal_peptide": False},
            size_class="large",
            quality_flags=["well_annotated"],
            composite_predicates=["nife_hydrogenase_validated"],
            unresolved_accessions=["PF12345"],
        )
        assert len(state.activities) == 2
        assert state.topology["tm_count"] == 7
        assert state.size_class == "large"

    def test_to_dict(self):
        """Should serialize to dictionary."""
        state = SemanticState(
            protein_id="prot_001",
            activities=["hydrogenase"],
            size_class="large",
        )
        d = state.to_dict()
        assert d["protein_id"] == "prot_001"
        assert d["activities"] == ["hydrogenase"]
        assert d["size_class"] == "large"
        assert d["topology"] == {}

    def test_from_dict(self):
        """Should deserialize from dictionary."""
        d = {
            "protein_id": "prot_001",
            "activities": ["hydrogenase"],
            "roles": [],
            "architecture": [],
            "localization": ["membrane"],
            "topology": {"tm_count": 3},
            "size_class": "medium",
            "quality_flags": [],
            "composite_predicates": [],
            "unresolved_accessions": [],
        }
        state = SemanticState.from_dict(d)
        assert state.protein_id == "prot_001"
        assert state.activities == ["hydrogenase"]
        assert state.localization == ["membrane"]
        assert state.topology["tm_count"] == 3

    def test_roundtrip(self):
        """to_dict -> from_dict should produce equal state."""
        original = SemanticState(
            protein_id="prot_002",
            activities=["kinase"],
            roles=["regulator"],
            size_class="medium",
            quality_flags=["confident_hit"],
        )
        reconstructed = SemanticState.from_dict(original.to_dict())
        assert original == reconstructed

    def test_from_dict_missing_fields(self):
        """Should handle missing optional fields gracefully."""
        d = {"protein_id": "prot_003"}
        state = SemanticState.from_dict(d)
        assert state.protein_id == "prot_003"
        assert state.activities == []
        assert state.size_class == ""
