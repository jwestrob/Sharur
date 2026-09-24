"""Component-level evidence never becomes a system-level claim."""

import pytest

from sharur.predicates.generator import AnnotationRecord, PredicateGenerator, ProteinRecord
from sharur.predicates.mappings.kegg_evidence import load_brite_rules, load_module_rules
from sharur.predicates.mappings.pfam_evidence import resolve
from sharur.predicates.mappings.pfam_map import PFAM_TO_PREDICATES
from sharur.predicates.vocabulary import (
    ALL_PREDICATES,
    COMPONENT,
    COMPONENT_EQUIVALENT,
    PREDICATE_BY_ID,
    SYSTEM,
    SYSTEM_PREDICATES,
    component_level,
    get_hierarchy,
)
from sharur.predicates_v2.generator import AtomGenerator


def test_every_predicate_has_a_level():
    assert {p.level for p in ALL_PREDICATES} == {COMPONENT, SYSTEM}


def test_no_component_predicate_expands_into_a_system_claim():
    offenders = [(p.predicate_id, ancestor) for p in ALL_PREDICATES if p.level == COMPONENT
                 for ancestor in get_hierarchy(p.predicate_id)[1:]
                 if PREDICATE_BY_ID[ancestor].level == SYSTEM]
    assert not offenders


def test_component_equivalents_are_component_level():
    for system, component in COMPONENT_EQUIVALENT.items():
        assert PREDICATE_BY_ID[system].level == SYSTEM
        assert PREDICATE_BY_ID[component].level == COMPONENT, system


def test_part_of_points_at_systems_and_never_expands():
    for p in ALL_PREDICATES:
        if p.part_of:
            assert PREDICATE_BY_ID[p.part_of].level == SYSTEM, p.predicate_id
            assert p.part_of not in get_hierarchy(p.predicate_id), p.predicate_id


def test_named_system_descriptions_are_system_level():
    # Descriptions that assert a confirmed or paired system belong to system-level predicates.
    for p in ALL_PREDICATES:
        if any(word in p.description.lower() for word in ("confirmed", "locus", "paired")):
            assert p.level == SYSTEM, p.predicate_id


def test_maps_and_rules_carry_component_predicates_only():
    assert not {p for preds in PFAM_TO_PREDICATES.values() for p in preds} & SYSTEM_PREDICATES
    assert not {p for rule in load_brite_rules() for p in rule.predicates} & SYSTEM_PREDICATES
    assert not {p for preds in load_module_rules().values() for p in preds} & SYSTEM_PREDICATES


@pytest.mark.parametrize(("proposed", "resolved"), [
    ("toxin_antitoxin", "defense_component"),
    ("abortive_infection", "abi_domain"),
    ("rm_type_ii", "restriction_modification"),
    ("type_iv_secretion", "t4ss_component"),
    ("kinase", "kinase"),
])
def test_system_proposals_resolve_to_components(proposed, resolved):
    assert resolve(proposed, PREDICATE_BY_ID) == resolved
    assert component_level(proposed) == resolved


def _preds(source, accession, name=""):
    gen = PredicateGenerator(include_direct_access=False)
    ann = AnnotationRecord(source=source, accession=accession, name=name or accession, evalue=1e-30)
    return set(gen.generate_for_protein(ProteinRecord(protein_id="p", sequence_length=300), [ann]))


def test_raw_defense_hmm_hits_stay_component_level():
    for accession in ("MazF_toxin", "AbiEii", "RM_Type_I_HsdR", "CBASS_Cap2", "Cas3"):
        preds = _preds("defensefinder", accession)
        assert not preds & SYSTEM_PREDICATES, (accession, preds & SYSTEM_PREDICATES)
        assert "defense_component" in preds, accession


def test_validated_system_calls_keep_system_predicates():
    preds = _preds("defensefinder_system", "RM_Type_I_HsdR", "RM/RM_Type_I")
    assert {"defense_system", "restriction_modification"} <= preds
    assert "secretion_system" in _preds("txsscan_system", "T4SS_virB4", "T4SS/T4SS_typeT")
    assert "type_iv_secretion" in _preds("txsscan_system", "T4SS_virB4", "T4SS/T4SS_typeT")


def test_v2_atoms_follow_the_same_rule():
    gen = AtomGenerator()
    protein = ProteinRecord(protein_id="p", sequence_length=300)
    raw = gen.generate_atoms(protein, [AnnotationRecord(source="defensefinder", accession="MazF_toxin",
                                                         name="MazF_toxin", evalue=1e-30)])
    assert not {a.atom_id for a in raw} & SYSTEM_PREDICATES


@pytest.mark.parametrize("complex_member", ["ech_hydrogenase", "mbh_hydrogenase"])
def test_named_complex_membership_implies_no_catalytic_group(complex_member):
    # Membrane subunits of Ech/Mbh belong to the complex; HydDB groups come from catalytic subunits.
    assert not {"nife_group4", "nife_hydrogenase"} & set(get_hierarchy(complex_member))
