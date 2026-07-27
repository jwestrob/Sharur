"""Tests for the deterministic candidate digest.

The digest exists so a model is never asked to do arithmetic. That only holds if
the arithmetic is right, so the reductions that change counts -- census removal,
nested folding, relabelling -- are pinned here.
"""

from __future__ import annotations

import json

from sharur.atlas_triage import (
    CENSUS_FRACTION,
    build_digest,
    flag_confounds,
    relabel,
    render,
)


def occ(system, accessions, genome, subtype="", hypothesis=""):
    return {
        "system": system,
        "genome": genome,
        "candidate_type": "co-located-pathway",
        "accessions": tuple(accessions),
        "subtype": subtype,
        "evidence": json.dumps({"hypothesis": hypothesis}) if hypothesis else "",
    }


class TestConfoundFlagging:
    def test_rubisco_without_prk_is_flagged(self):
        assert "RuBisCO-without-PRK" in flag_confounds({"RuBisCO_large", "RuBisCO_small"})

    def test_rubisco_with_prk_is_clean(self):
        assert flag_confounds({"RuBisCO_large", "RuBisCO_small", "PRK"}) == []

    def test_group4_hydrogenase_without_nifese_is_flagged(self):
        assert "hydrogenase-vs-ComplexI" in flag_confounds({"nife_group4", "Complex1_49kDa"})

    def test_methyl_branch_is_not_wood_ljungdahl(self):
        assert "WoodLjungdahl-methyl-branch" in flag_confounds({"FTHFS", "MTHFR"})
        assert flag_confounds({"FTHFS", "MTHFR", "AcsB"}) == []


class TestRelabelling:
    def test_pilus_components_are_not_counted_as_competence(self):
        """A scanner's label is evidence of what it thought, not of what is there."""
        assert relabel("natural-competence", {"Tad", "TadE"}) == "pilus-assembly(relabelled)"

    def test_genuine_competence_keeps_its_label(self):
        assert relabel("natural-competence", {"Competence", "K02237"}) == "natural-competence"

    def test_mixed_evidence_keeps_the_original_label(self):
        """With a real competence marker present, the pilus genes are context."""
        assert relabel("natural-competence", {"Tad", "Competence"}) == "natural-competence"


class TestDigest:
    def test_census_systems_are_identified_and_deranked(self):
        genomes = [f"g{i}" for i in range(10)]
        occs = [occ("ubiquitous", ["X", "Y"], g) for g in genomes]
        occs += [occ("rare", ["P", "Q"], g) for g in genomes[:3]]
        d = build_digest(occs)
        assert "ubiquitous" in d["census"]
        assert "rare" not in d["census"]
        assert d["entries"][0]["system"] == "rare", "census must not outrank a real finding"

    def test_singletons_are_dropped_by_min_genomes(self):
        occs = [occ("s", ["A", "B"], "g1")] + [occ("t", ["C", "D"], f"g{i}") for i in (2, 3)]
        systems = {e["system"] for e in build_digest(occs)["entries"]}
        assert systems == {"t"}

    def test_prevalence_and_counts_are_computed_not_guessed(self):
        occs = [occ("s", ["A"], f"g{i}") for i in range(4)]
        occs += [occ("other", ["Z"], f"g{i}") for i in range(10)]
        entry = next(e for e in build_digest(occs)["entries"] if e["system"] == "s")
        assert entry["genomes"] == 4
        assert entry["prevalence"] == 4 / 10

    def test_empty_input_is_not_an_error(self):
        d = build_digest([])
        assert d["total"] == 0
        assert render(d) == "no candidates"

    def test_render_reports_the_census_fraction(self):
        genomes = [f"g{i}" for i in range(10)]
        occs = [occ("ubiquitous", ["X"], g) for g in genomes]
        occs += [occ("rare", ["P", "Q"], g) for g in genomes[:2]]
        text = render(build_digest(occs), top=5)
        assert "inventory, not findings" in text
        assert f"{CENSUS_FRACTION:.0%}" in text
