"""Pfam domains match by accession or profile name.

Stage 07 stores the profile name in the accession column when its reference
map lacks the profile, so supporting-domain checks must recognize both forms.
"""

import pytest

from sharur.predicates.generator import AnnotationRecord, PredicateGenerator, ProteinRecord
from sharur.predicates.mappings.pfam_map import get_predicates_for_pfam
from sharur.predicates.pfam_identity import (
    COMPLEX1,
    NIFESE_HASES,
    has_pfam_domain,
    normalize_pfam_accession,
    pfam_keys,
)
from sharur.predicates_v2.generator import AtomGenerator
from sharur.predicates_v2.model import ClaimRelation
from sharur.predicates_v2.rules import clear_caches, get_relation


@pytest.fixture(autouse=True)
def _reset_caches():
    clear_caches()
    yield
    clear_caches()


def _pfam(accession, name):
    return AnnotationRecord(
        source="pfam", accession=accession, name=name, evalue=1e-40, score=150.0,
    )


def _hyddb_nife():
    return AnnotationRecord(source="hyddb", accession="NiFe", evalue=1e-50, score=300.0)


def _excluded_ids(atoms):
    return {a.atom_id for a in atoms if a.relation == ClaimRelation.excludes}


class TestPfamIdentity:
    def test_normalize_strips_version(self):
        assert normalize_pfam_accession("PF00374.26") == "PF00374"
        assert normalize_pfam_accession(" NiFeSe_Hases ") == "NiFeSe_Hases"
        assert normalize_pfam_accession(None) == ""

    def test_keys_cover_accession_and_name(self):
        assert pfam_keys("PF00374.26", "NiFeSe_Hases") == ("PF00374", "NiFeSe_Hases")
        assert pfam_keys("NiFeSe_Hases", "NiFeSe_Hases") == ("NiFeSe_Hases",)

    @pytest.mark.parametrize("identifier", ["PF00374", "PF00374.26", "NiFeSe_Hases"])
    def test_domain_matches_every_stored_form(self, identifier):
        assert has_pfam_domain([identifier], NIFESE_HASES)
        assert not has_pfam_domain([identifier], COMPLEX1)


def _review_atoms(atoms):
    return [a for a in atoms if a.atom_id == "hydrogenase_complex1_review"]


class TestV1HydrogenaseValidation:
    @pytest.mark.parametrize("complex1", ["PF00346", "Complex1_49kDa", "PF00329", "Complex1_30kDa"])
    def test_complex1_without_nifese_flags_for_review(self, complex1):
        gen = PredicateGenerator()
        preds = gen._validate_hydrogenase_calls(
            {"hyddb:NiFe", "nife_hydrogenase", "hydrogenase", f"pfam:{complex1}"}
        )
        assert "hydrogenase_complex1_review" in preds
        assert "nife_hydrogenase" in preds

    @pytest.mark.parametrize("nifese", ["PF00374", "NiFeSe_Hases"])
    def test_nifese_by_either_form_clears_the_domain_check(self, nifese):
        gen = PredicateGenerator()
        preds = gen._validate_hydrogenase_calls(
            {"hyddb:NiFe", "nife_hydrogenase", "hydrogenase",
             "pfam:Complex1_49kDa", f"pfam:{nifese}"}
        )
        assert "hydrogenase_complex1_review" not in preds


class TestV2HydrogenaseValidation:
    def test_name_form_nifese_clears_review(self):
        atoms = AtomGenerator().generate_atoms(
            ProteinRecord(protein_id="p", sequence_length=550),
            [_hyddb_nife(), _pfam("NiFeSe_Hases", "NiFeSe_Hases"),
             _pfam("Complex1_49kDa", "Complex1_49kDa")],
        )
        assert not _review_atoms(atoms)
        assert not _excluded_ids(atoms)

    def test_name_form_complex1_flags_review_without_exclusion(self):
        atoms = AtomGenerator().generate_atoms(
            ProteinRecord(protein_id="p", sequence_length=550),
            [_hyddb_nife(), _pfam("Complex1_49kDa", "Complex1_49kDa")],
        )
        review = _review_atoms(atoms)
        assert [a.relation for a in review] == [ClaimRelation.flags]
        assert review[0].source_accession == "Complex1_49kDa"
        assert not _excluded_ids(atoms)
        assert "nife_hydrogenase" in {a.atom_id for a in atoms}

    def test_equivalent_identifier_forms_give_identical_decisions(self):
        def decision(nifese, complex1):
            atoms = AtomGenerator().generate_atoms(
                ProteinRecord(protein_id="p", sequence_length=550),
                [_hyddb_nife(), _pfam(*nifese), _pfam(*complex1)],
            )
            return bool(_review_atoms(atoms))

        forms = [
            (("PF00374", "NiFeSe_Hases"), ("PF00346", "Complex1_49kDa")),
            (("PF00374.26", "NiFeSe_Hases"), ("PF00346.29", "Complex1_49kDa")),
            (("NiFeSe_Hases", "NiFeSe_Hases"), ("Complex1_49kDa", "Complex1_49kDa")),
            ((" NiFeSe_Hases ", None), (" Complex1_49kDa", None)),
        ]
        assert {decision(n, c) for n, c in forms} == {False}

    def test_kegg_only_nife_is_not_reviewed(self):
        atoms = AtomGenerator().generate_atoms(
            ProteinRecord(protein_id="p", sequence_length=550),
            [AnnotationRecord(source="kegg", accession="K00437", evalue=1e-50, score=300.0),
             _pfam("PF00346", "Complex1_49kDa")],
        )
        assert not _review_atoms(atoms)


class TestRelationOverrides:
    @pytest.mark.parametrize(
        ("accession", "name"),
        [("PF00374", None), ("PF00374.26", None), ("NiFeSe_Hases", "NiFeSe_Hases")],
    )
    def test_nifese_override_applies_to_every_form(self, accession, name):
        assert get_relation("pfam", accession, name) == ClaimRelation.implies

    def test_unlisted_profile_uses_source_default(self):
        assert get_relation("pfam", "Unlisted_domain", "Unlisted_domain") == ClaimRelation.supports


class TestDirectMapping:
    def test_versioned_accession_uses_accession_mapping(self):
        assert "nife_hydrogenase" in get_predicates_for_pfam("PF00374.26", "NiFeSe_Hases")

    def test_name_keyed_mapping_applies_to_name_form_rows(self, monkeypatch):
        import sharur.predicates.mappings.pfam_map as pfam_map

        monkeypatch.setitem(pfam_map.PFAM_TO_PREDICATES, "Example_domain", ["transporter"])
        assert "transporter" in get_predicates_for_pfam("Example_domain", "Example_domain")
        assert "transporter" not in get_predicates_for_pfam("PF99999", "Example_domain")


def test_both_complex1_domains_cite_the_49kda_subunit():
    atoms = AtomGenerator().generate_atoms(
        ProteinRecord(protein_id="p", sequence_length=550),
        [_hyddb_nife(), _pfam("PF00329", "Complex1_30kDa"), _pfam("PF00346", "Complex1_49kDa")],
    )
    assert {a.source_accession for a in _review_atoms(atoms)} == {"PF00346"}
