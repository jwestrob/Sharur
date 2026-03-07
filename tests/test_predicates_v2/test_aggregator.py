"""Tests for V2 aggregator."""

import pytest

from sharur.predicates_v2.aggregator import aggregate_atoms
from sharur.predicates_v2.model import (
    ClaimRelation,
    SemanticAtom,
    SemanticFacet,
    SemanticState,
)


def _atom(
    atom_id: str,
    facet: SemanticFacet,
    relation: ClaimRelation = ClaimRelation.supports,
    source_accession: str = "TEST",
    source_db: str = "pfam",
    protein_id: str = "prot_001",
    evalue: float | None = None,
    score: float | None = None,
) -> SemanticAtom:
    """Helper to create test atoms."""
    return SemanticAtom(
        protein_id=protein_id,
        atom_id=atom_id,
        facet=facet,
        relation=relation,
        source_accession=source_accession,
        source_db=source_db,
        evidence_evalue=evalue,
        evidence_score=score,
    )


class TestAggregateAtoms:
    """Tests for aggregate_atoms function."""

    def test_empty_atoms(self):
        """Should produce empty state from no atoms."""
        state = aggregate_atoms("prot_001", [])
        assert state.protein_id == "prot_001"
        assert state.activities == []
        assert state.roles == []
        assert state.size_class == ""

    def test_single_activity_atom(self):
        """Should aggregate a single activity atom."""
        atoms = [
            _atom("hydrogenase", SemanticFacet.activity, ClaimRelation.implies),
        ]
        state = aggregate_atoms("prot_001", atoms)
        assert "hydrogenase" in state.activities

    def test_multiple_facets(self):
        """Should aggregate atoms across different facets."""
        atoms = [
            _atom("hydrogenase", SemanticFacet.activity, ClaimRelation.implies),
            _atom("membrane", SemanticFacet.localization, ClaimRelation.supports),
            _atom("giant", SemanticFacet.size_class, ClaimRelation.implies),
        ]
        state = aggregate_atoms("prot_001", atoms)
        assert "hydrogenase" in state.activities
        assert "membrane" in state.localization
        assert state.size_class == "giant"

    def test_excludes_beats_supports(self):
        """Excludes should remove an atom that only has supports."""
        atoms = [
            _atom("hydrogenase", SemanticFacet.activity, ClaimRelation.supports,
                  source_accession="PF00374"),
            _atom("hydrogenase", SemanticFacet.activity, ClaimRelation.excludes,
                  source_accession="PF00346"),
        ]
        state = aggregate_atoms("prot_001", atoms)
        assert "hydrogenase" not in state.activities

    def test_excludes_beats_flags(self):
        """Excludes should remove an atom that only has flags."""
        atoms = [
            _atom("defense_system", SemanticFacet.role, ClaimRelation.flags),
            _atom("defense_system", SemanticFacet.role, ClaimRelation.excludes,
                  source_accession="OVERRIDE"),
        ]
        state = aggregate_atoms("prot_001", atoms)
        assert "defense_system" not in state.roles

    def test_implies_plus_excludes_is_conflict(self):
        """Implies + excludes = conflict, surfaced as unresolved."""
        atoms = [
            _atom("hydrogenase", SemanticFacet.activity, ClaimRelation.implies,
                  source_accession="HYDDB_NiFe"),
            _atom("hydrogenase", SemanticFacet.activity, ClaimRelation.excludes,
                  source_accession="PF00346"),
        ]
        state = aggregate_atoms("prot_001", atoms)
        # Conflict means the atom is NOT resolved
        assert "hydrogenase" not in state.activities
        # But it IS surfaced in unresolved
        assert "conflict:hydrogenase" in state.unresolved_accessions

    def test_multiple_supports_convergent(self):
        """Multiple supports for the same atom should keep it."""
        atoms = [
            _atom("kinase", SemanticFacet.activity, ClaimRelation.supports,
                  source_accession="PF00069"),
            _atom("kinase", SemanticFacet.activity, ClaimRelation.supports,
                  source_accession="K00001"),
        ]
        state = aggregate_atoms("prot_001", atoms)
        assert "kinase" in state.activities

    def test_deduplication(self):
        """Same atom_id should appear only once in resolved state."""
        atoms = [
            _atom("transporter", SemanticFacet.role, ClaimRelation.supports,
                  source_accession="PF00005"),
            _atom("transporter", SemanticFacet.role, ClaimRelation.supports,
                  source_accession="K02500"),
        ]
        state = aggregate_atoms("prot_001", atoms)
        assert state.roles.count("transporter") == 1

    def test_unresolved_atoms_not_in_state(self):
        """Unresolved atoms should not appear in resolved facets."""
        atoms = [
            _atom("hydrogenase", SemanticFacet.activity, ClaimRelation.implies),
            SemanticAtom(
                protein_id="prot_001",
                atom_id="unresolved:PF12345",
                facet=SemanticFacet.quality_flag,
                relation=ClaimRelation.unresolved,
                source_accession="PF12345",
                source_db="pfam",
            ),
        ]
        state = aggregate_atoms("prot_001", atoms)
        assert "hydrogenase" in state.activities
        assert "PF12345" in state.unresolved_accessions

    def test_size_priority(self):
        """Should pick the largest size class when multiple are present."""
        atoms = [
            _atom("giant", SemanticFacet.size_class, ClaimRelation.implies),
            _atom("massive", SemanticFacet.size_class, ClaimRelation.implies),
        ]
        state = aggregate_atoms("prot_001", atoms)
        assert state.size_class == "massive"

    def test_topology_transmembrane(self):
        """Should build topology dict for transmembrane proteins."""
        atoms = [
            _atom("transmembrane_predicted", SemanticFacet.topology),
            _atom("multi_pass_membrane", SemanticFacet.topology),
            _atom("n_in_topology", SemanticFacet.topology),
        ]
        state = aggregate_atoms("prot_001", atoms)
        assert state.topology["tm_count"] >= 2
        assert state.topology["n_terminus"] == "inside"

    def test_topology_soluble(self):
        """Should build topology dict for soluble proteins."""
        atoms = [
            _atom("soluble_predicted", SemanticFacet.topology),
        ]
        state = aggregate_atoms("prot_001", atoms)
        assert state.topology["tm_count"] == 0

    def test_topology_polytopic(self):
        """Should recognize polytopic membrane proteins."""
        atoms = [
            _atom("transmembrane_predicted", SemanticFacet.topology),
            _atom("multi_pass_membrane", SemanticFacet.topology),
            _atom("polytopic_membrane", SemanticFacet.topology),
        ]
        state = aggregate_atoms("prot_001", atoms)
        assert state.topology["tm_count"] >= 4

    def test_sorted_output(self):
        """Resolved lists should be sorted."""
        atoms = [
            _atom("kinase", SemanticFacet.activity),
            _atom("hydrolase", SemanticFacet.activity),
            _atom("atpase", SemanticFacet.activity),
        ]
        state = aggregate_atoms("prot_001", atoms)
        assert state.activities == sorted(state.activities)

    def test_quality_flags(self):
        """Should aggregate quality flag atoms."""
        atoms = [
            _atom("well_annotated", SemanticFacet.quality_flag, ClaimRelation.implies),
            _atom("multi_source", SemanticFacet.quality_flag, ClaimRelation.implies),
        ]
        state = aggregate_atoms("prot_001", atoms)
        assert "well_annotated" in state.quality_flags
        assert "multi_source" in state.quality_flags

    def test_multiple_unresolved(self):
        """Should collect all unique unresolved accessions."""
        atoms = [
            SemanticAtom(
                protein_id="prot_001",
                atom_id="unresolved:PF12345",
                facet=SemanticFacet.quality_flag,
                relation=ClaimRelation.unresolved,
                source_accession="PF12345",
                source_db="pfam",
            ),
            SemanticAtom(
                protein_id="prot_001",
                atom_id="unresolved:PF67890",
                facet=SemanticFacet.quality_flag,
                relation=ClaimRelation.unresolved,
                source_accession="PF67890",
                source_db="pfam",
            ),
        ]
        state = aggregate_atoms("prot_001", atoms)
        assert "PF12345" in state.unresolved_accessions
        assert "PF67890" in state.unresolved_accessions

    def test_architecture_atoms(self):
        """Should aggregate architecture atoms."""
        atoms = [
            _atom("multi_domain", SemanticFacet.architecture),
            _atom("beta_barrel", SemanticFacet.architecture),
        ]
        state = aggregate_atoms("prot_001", atoms)
        assert "multi_domain" in state.architecture
        assert "beta_barrel" in state.architecture
