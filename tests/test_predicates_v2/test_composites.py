"""Tests for V2 composite predicate evaluator."""

import pytest

from sharur.predicates_v2.composites import (
    CompositeDefinition,
    clear_composites_cache,
    evaluate_composites,
    load_composites,
)
from sharur.predicates_v2.model import ClaimRelation, SemanticAtom, SemanticFacet


@pytest.fixture(autouse=True)
def _reset_caches():
    """Clear composite cache before each test."""
    clear_composites_cache()
    yield
    clear_composites_cache()


def _atom(
    atom_id: str,
    facet: SemanticFacet = SemanticFacet.activity,
    relation: ClaimRelation = ClaimRelation.supports,
    source_db: str = "pfam",
    source_accession: str = "TEST",
) -> SemanticAtom:
    """Helper to create test atoms."""
    return SemanticAtom(
        protein_id="prot_001",
        atom_id=atom_id,
        facet=facet,
        relation=relation,
        source_accession=source_accession,
        source_db=source_db,
    )


class TestLoadComposites:
    """Tests for loading composite definitions from YAML."""

    def test_loads_from_default_config(self):
        """Should load composites from default YAML."""
        composites = load_composites()
        assert len(composites) > 0

    def test_composite_has_required_fields(self):
        """Each composite should have name, description, facet, requires."""
        composites = load_composites()
        for comp in composites:
            assert comp.name, f"Composite missing name"
            assert comp.facet, f"Composite {comp.name} missing facet"
            assert comp.requires, f"Composite {comp.name} missing requires"


class TestHasAtomOperator:
    """Tests for the has_atom DSL operator."""

    def test_basic_has_atom(self):
        """Should match when atom exists."""
        comp = CompositeDefinition(
            name="test",
            description="test",
            facet="activity",
            requires={"all_of": [{"has_atom": "hydrogenase"}]},
        )
        atoms = [_atom("hydrogenase")]
        result = evaluate_composites(atoms, [comp])
        assert "test" in result

    def test_missing_atom(self):
        """Should not match when atom is missing."""
        comp = CompositeDefinition(
            name="test",
            description="test",
            facet="activity",
            requires={"all_of": [{"has_atom": "hydrogenase"}]},
        )
        atoms = [_atom("kinase")]
        result = evaluate_composites(atoms, [comp])
        assert "test" not in result

    def test_has_atom_with_relation(self):
        """Should filter by exact relation."""
        comp = CompositeDefinition(
            name="test",
            description="test",
            facet="activity",
            requires={
                "all_of": [{"has_atom": "defense_system", "relation": "implies"}],
            },
        )
        # Only supports, not implies
        atoms = [_atom("defense_system", relation=ClaimRelation.supports)]
        result = evaluate_composites(atoms, [comp])
        assert "test" not in result

        # With implies
        atoms = [_atom("defense_system", relation=ClaimRelation.implies)]
        result = evaluate_composites(atoms, [comp])
        assert "test" in result

    def test_has_atom_with_relation_in(self):
        """Should filter by relation list."""
        comp = CompositeDefinition(
            name="test",
            description="test",
            facet="activity",
            requires={
                "all_of": [
                    {"has_atom": "hydrogenase", "relation_in": ["implies", "supports"]},
                ],
            },
        )
        # flags should not match
        atoms = [_atom("hydrogenase", relation=ClaimRelation.flags)]
        result = evaluate_composites(atoms, [comp])
        assert "test" not in result

        # supports should match
        atoms = [_atom("hydrogenase", relation=ClaimRelation.supports)]
        result = evaluate_composites(atoms, [comp])
        assert "test" in result

    def test_has_atom_with_source_db(self):
        """Should filter by source database."""
        comp = CompositeDefinition(
            name="test",
            description="test",
            facet="role",
            requires={
                "all_of": [
                    {"has_atom": "defense_system", "source_db": "defensefinder_system"},
                ],
            },
        )
        # Wrong source
        atoms = [_atom("defense_system", source_db="defensefinder")]
        result = evaluate_composites(atoms, [comp])
        assert "test" not in result

        # Correct source
        atoms = [_atom("defense_system", source_db="defensefinder_system")]
        result = evaluate_composites(atoms, [comp])
        assert "test" in result


class TestAllOfOperator:
    """Tests for the all_of DSL operator."""

    def test_all_conditions_met(self):
        """Should match when all conditions are true."""
        comp = CompositeDefinition(
            name="test",
            description="test",
            facet="quality_flag",
            requires={
                "all_of": [
                    {"has_atom": "giant"},
                    {"has_atom": "unannotated"},
                ],
            },
        )
        atoms = [
            _atom("giant", facet=SemanticFacet.size_class),
            _atom("unannotated", facet=SemanticFacet.quality_flag),
        ]
        result = evaluate_composites(atoms, [comp])
        assert "test" in result

    def test_not_all_conditions_met(self):
        """Should not match when any condition is false."""
        comp = CompositeDefinition(
            name="test",
            description="test",
            facet="quality_flag",
            requires={
                "all_of": [
                    {"has_atom": "giant"},
                    {"has_atom": "unannotated"},
                ],
            },
        )
        atoms = [_atom("giant", facet=SemanticFacet.size_class)]
        result = evaluate_composites(atoms, [comp])
        assert "test" not in result


class TestAnyOfOperator:
    """Tests for the any_of DSL operator."""

    def test_one_condition_met(self):
        """Should match when at least one condition is true."""
        comp = CompositeDefinition(
            name="test",
            description="test",
            facet="activity",
            requires={
                "any_of": [
                    {"has_atom": "nife_group4"},
                    {"has_atom": "mbh_hydrogenase"},
                    {"has_atom": "ech_hydrogenase"},
                ],
            },
        )
        atoms = [_atom("mbh_hydrogenase")]
        result = evaluate_composites(atoms, [comp])
        assert "test" in result

    def test_no_condition_met(self):
        """Should not match when no condition is true."""
        comp = CompositeDefinition(
            name="test",
            description="test",
            facet="activity",
            requires={
                "any_of": [
                    {"has_atom": "nife_group4"},
                    {"has_atom": "mbh_hydrogenase"},
                ],
            },
        )
        atoms = [_atom("kinase")]
        result = evaluate_composites(atoms, [comp])
        assert "test" not in result


class TestNoneOfOperator:
    """Tests for the none_of DSL operator."""

    def test_none_matched(self):
        """Should match when no exclusion condition is true."""
        comp = CompositeDefinition(
            name="test",
            description="test",
            facet="activity",
            requires={
                "all_of": [
                    {"has_atom": "hydrogenase"},
                ],
                "none_of": [
                    {"has_atom": "complex_I_subunit", "relation": "implies"},
                ],
            },
        )
        atoms = [_atom("hydrogenase")]
        result = evaluate_composites(atoms, [comp])
        assert "test" in result

    def test_exclusion_triggered(self):
        """Should not match when an exclusion condition is true."""
        comp = CompositeDefinition(
            name="test",
            description="test",
            facet="activity",
            requires={
                "all_of": [
                    {"has_atom": "hydrogenase"},
                ],
                "none_of": [
                    {"has_atom": "complex_I_subunit", "relation": "implies"},
                ],
            },
        )
        atoms = [
            _atom("hydrogenase"),
            _atom("complex_I_subunit", relation=ClaimRelation.implies),
        ]
        result = evaluate_composites(atoms, [comp])
        assert "test" not in result

    def test_exclusion_wrong_relation(self):
        """Exclusion should not trigger if relation doesn't match."""
        comp = CompositeDefinition(
            name="test",
            description="test",
            facet="activity",
            requires={
                "all_of": [
                    {"has_atom": "hydrogenase"},
                ],
                "none_of": [
                    {"has_atom": "complex_I_subunit", "relation": "implies"},
                ],
            },
        )
        atoms = [
            _atom("hydrogenase"),
            _atom("complex_I_subunit", relation=ClaimRelation.flags),
        ]
        result = evaluate_composites(atoms, [comp])
        assert "test" in result  # flags != implies, so exclusion doesn't trigger


class TestPropertyOperators:
    """Tests for property_equals, property_gte, property_lte operators."""

    def test_property_equals(self):
        """Should match when topology property equals value."""
        comp = CompositeDefinition(
            name="test",
            description="test",
            facet="topology",
            requires={"property_equals": {"tm_count": 7}},
        )
        atoms = []
        result = evaluate_composites(atoms, [comp], topology={"tm_count": 7})
        assert "test" in result

        result = evaluate_composites(atoms, [comp], topology={"tm_count": 3})
        assert "test" not in result

    def test_property_gte(self):
        """Should match when topology property >= value."""
        comp = CompositeDefinition(
            name="test",
            description="test",
            facet="topology",
            requires={"property_gte": {"tm_count": 4}},
        )
        atoms = []
        result = evaluate_composites(atoms, [comp], topology={"tm_count": 7})
        assert "test" in result

        result = evaluate_composites(atoms, [comp], topology={"tm_count": 4})
        assert "test" in result

        result = evaluate_composites(atoms, [comp], topology={"tm_count": 2})
        assert "test" not in result

    def test_property_lte(self):
        """Should match when topology property <= value."""
        comp = CompositeDefinition(
            name="test",
            description="test",
            facet="topology",
            requires={"property_lte": {"tm_count": 2}},
        )
        atoms = []
        result = evaluate_composites(atoms, [comp], topology={"tm_count": 1})
        assert "test" in result

        result = evaluate_composites(atoms, [comp], topology={"tm_count": 5})
        assert "test" not in result

    def test_property_missing_key(self):
        """Missing property key should not match gte/lte."""
        comp = CompositeDefinition(
            name="test",
            description="test",
            facet="topology",
            requires={"property_gte": {"tm_count": 1}},
        )
        atoms = []
        result = evaluate_composites(atoms, [comp], topology={})
        assert "test" not in result


class TestCompositeEvaluation:
    """Integration tests for composite evaluation."""

    def test_giant_unannotated_from_yaml(self):
        """Should evaluate giant_unannotated from YAML config."""
        composites = load_composites()
        atoms = [
            _atom("giant", facet=SemanticFacet.size_class, relation=ClaimRelation.implies),
            _atom("unannotated", facet=SemanticFacet.quality_flag, relation=ClaimRelation.implies),
        ]
        result = evaluate_composites(atoms, composites)
        assert "giant_unannotated" in result

    def test_energy_conserving_hydrogenase_from_yaml(self):
        """Should evaluate energy_conserving_hydrogenase from YAML config."""
        composites = load_composites()
        atoms = [
            _atom("nife_group4", relation=ClaimRelation.implies),
        ]
        result = evaluate_composites(atoms, composites)
        assert "energy_conserving_hydrogenase" in result

    def test_novel_membrane_protein_from_yaml(self):
        """Should evaluate novel_membrane_protein from YAML config."""
        composites = load_composites()
        atoms = [
            _atom("unannotated", facet=SemanticFacet.quality_flag),
            _atom("transmembrane_predicted", facet=SemanticFacet.topology),
        ]
        result = evaluate_composites(atoms, composites)
        assert "novel_membrane_protein" in result

    def test_defense_system_validated_from_yaml(self):
        """Should evaluate defense_system_validated from YAML config."""
        composites = load_composites()
        atoms = [
            _atom("defense_system", facet=SemanticFacet.role,
                  relation=ClaimRelation.implies,
                  source_db="defensefinder_system"),
        ]
        result = evaluate_composites(atoms, composites)
        assert "defense_system_validated" in result

    def test_multiple_composites_can_match(self):
        """Multiple composites should be able to match for one protein."""
        composites = load_composites()
        atoms = [
            _atom("giant", facet=SemanticFacet.size_class, relation=ClaimRelation.implies),
            _atom("unannotated", facet=SemanticFacet.quality_flag, relation=ClaimRelation.implies),
            _atom("transmembrane_predicted", facet=SemanticFacet.topology),
        ]
        result = evaluate_composites(atoms, composites)
        assert "giant_unannotated" in result
        assert "novel_membrane_protein" in result

    def test_sorted_output(self):
        """Output should be sorted."""
        composites = load_composites()
        atoms = [
            _atom("giant", facet=SemanticFacet.size_class, relation=ClaimRelation.implies),
            _atom("unannotated", facet=SemanticFacet.quality_flag, relation=ClaimRelation.implies),
            _atom("transmembrane_predicted", facet=SemanticFacet.topology),
        ]
        result = evaluate_composites(atoms, composites)
        assert result == sorted(result)
