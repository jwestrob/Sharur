"""Every shipped Pfam -> predicate pair is supported by Pfam's own information.

``pfam_predicates.tsv`` records one evidence item per pair; this suite
re-verifies each item against ``pfam_evidence_snapshot.tsv`` (the family's Pfam
name/description, InterPro GO anchors, and ENZYME names at build time) using the
current evidence rules. Swiss-Prot consensus pairs are checked against the
consensus thresholds here; set ``SHARUR_SWISSPROT`` (uniprot_sprot.dat.gz) and
``SHARUR_GO_OBO`` (go-basic.obo) to recount them from the reviewed proteins. Rebuild with ``scripts/build_pfam_predicate_map.py``
after editing proposals or evidence definitions.
"""

import os
import re
from pathlib import Path

import pytest

from sharur.predicates.mappings import pfam_map
from sharur.predicates.mappings.kegg_map import get_predicates_for_ec
from sharur.predicates.mappings.pfam_evidence import (
    EVIDENCE,
    SWISSPROT_COVERAGE,
    SWISSPROT_EXCLUDED,
    SWISSPROT_KO_MAJORITY,
    SWISSPROT_MIN_LOWER_BOUND,
    SWISSPROT_SINGLE,
    judge,
    text_evidence,
    wilson_lower_bound,
)
from sharur.predicates.mappings.pfam_evidence_spec import HOMONYMS, SUBSET_NAMES, TARGET_NAMES
from sharur.predicates.vocabulary import PREDICATE_BY_ID


DATA = Path(pfam_map.__file__).with_name("data")


def _rows(name):
    return [line.rstrip("\n").split("\t") for line in (DATA / name).read_text().splitlines()
            if line and not line.startswith("#")]


MAP = {acc: (name, preds.split(","), dict(i.split("=", 1) for i in ev.split(";")))
       for acc, name, preds, ev in _rows("pfam_predicates.tsv")}
SNAPSHOT = {acc: {"name": name, "desc": desc, "go": set(filter(None, go.split(","))),
                  "enzymes": {p: ecs.split("|") for p, ecs in (x.split("=", 1) for x in enz.split(";") if x)}}
            for acc, name, desc, go, enz in _rows("pfam_evidence_snapshot.tsv")}
PAIRS = [(acc, pred, ev) for acc, (_, _, evidence) in MAP.items() for pred, ev in evidence.items()]
SWISSPROT = re.compile(r"swissprot:(\d+)/(\d+)(?: KOs (\d+)/(\d+))? (?:single-domain (\d+)/(\d+)|co-domains excluded)")


def _meets_consensus(pred, detail):
    m = SWISSPROT.fullmatch(f"swissprot:{detail}")
    if not m or pred in SWISSPROT_EXCLUDED:
        return False
    k, n = int(m.group(1)), int(m.group(2))
    covered = n >= SWISSPROT_COVERAGE[0] and k >= SWISSPROT_COVERAGE[1] * n \
        and wilson_lower_bound(k, n) >= SWISSPROT_MIN_LOWER_BOUND
    if m.group(3) is not None:
        covered = covered and int(m.group(4)) >= 2 and int(m.group(3)) > SWISSPROT_KO_MAJORITY * int(m.group(4))
    if m.group(5) is None:
        return covered
    ks, ns = int(m.group(5)), int(m.group(6))
    return covered and ns >= SWISSPROT_SINGLE[0] and ks >= SWISSPROT_SINGLE[1] * ns


def test_map_and_snapshot_cover_the_same_families():
    assert set(MAP) == set(SNAPSHOT)


def test_every_predicate_is_in_the_vocabulary():
    assert not {pred for _, pred, _ in PAIRS} - set(PREDICATE_BY_ID)


def test_every_pair_reverifies_against_the_snapshot():
    failures = []
    for acc, pred, ev in PAIRS:
        snap = SNAPSHOT[acc]
        kind, _, detail = ev.partition(":")
        if kind == "go":
            ok = detail in EVIDENCE[pred]["go"] and detail in snap["go"]
        elif kind == "text":
            ok = text_evidence(pred, snap["name"], snap["desc"]) is not None
        elif kind == "enzyme":
            m = re.fullmatch(r"(\S+) \((.+)\)", detail)
            ok = bool(m) and m.group(1) in snap["enzymes"].get(m.group(2), ()) \
                and pred in get_predicates_for_ec(m.group(1))
        elif kind == "swissprot":
            ok = _meets_consensus(pred, detail)
        else:
            ok = False
        subset = kind in ("text", "enzyme") and (
            (snap["name"], pred) in SUBSET_NAMES or (snap["name"], pred) in TARGET_NAMES)
        if not ok or subset or (snap["name"], pred) in HOMONYMS:
            failures.append((acc, snap["name"], pred, ev))
    assert not failures, failures[:20]


def test_swissprot_pairs_record_their_release():
    header = (DATA / "pfam_predicates.tsv").read_text().split("# accession", 1)[0]
    if any(ev.startswith("swissprot:") for _, _, ev in PAIRS):
        assert "# Swiss-Prot: " in header


@pytest.mark.skipif(not (os.environ.get("SHARUR_SWISSPROT") and os.environ.get("SHARUR_GO_OBO")),
                    reason="set SHARUR_SWISSPROT and SHARUR_GO_OBO (and SHARUR_SWISSPROT_KEGG when the "
                           "map was built with KEGG gene links) to recount Swiss-Prot consensus")
def test_swissprot_counts_recount_from_reviewed_proteins():
    from sharur.predicates.mappings.kegg_map import KEGG_EVIDENCE
    from sharur.predicates.mappings.swissprot_evidence import (
        go_ancestors,
        ortholog_support,
        pfam_view,
        read_gene_ko_links,
        read_swissprot,
    )

    reviewed = read_swissprot(Path(os.environ["SHARUR_SWISSPROT"]), go_ancestors(Path(os.environ["SHARUR_GO_OBO"])))
    kegg_links = os.environ.get("SHARUR_SWISSPROT_KEGG")
    gene_ko = read_gene_ko_links(Path(kegg_links)) if kegg_links else None
    ko_predicates = {ko: {p for p, e in ev.items() if not e.startswith("swissprot:")}
                     for ko, ev in KEGG_EVIDENCE.items()}
    proteins = pfam_view(reviewed, gene_ko, ko_predicates)
    mismatched = []
    for acc, pred, ev in PAIRS:
        m = SWISSPROT.fullmatch(ev)
        if not m:
            continue
        carriers = [(fams, preds, kos) for fams, preds, kos in proteins if acc in fams]
        counted = [sum(pred in p for _, p, _ in carriers), len(carriers)]
        k_kos, n_kos = ortholog_support([(p, kos) for _, p, kos in carriers], pred)
        if n_kos >= 2:
            counted += [k_kos, n_kos]
        if m.group(5) is not None:
            single = [p for fams, p, _ in carriers if len(fams) == 1]
            counted += [sum(pred in p for p in single), len(single)]
        if [int(g) for g in m.groups() if g is not None] != counted:
            mismatched.append((acc, pred, ev, counted))
    assert not mismatched, mismatched[:20]


def test_runtime_uses_only_the_generated_map():
    assert not hasattr(pfam_map, "PFAM_PATTERNS")
    assert not hasattr(pfam_map, "EXTERNAL_PFAM_MAPPING_CANDIDATES")
    assert pfam_map.get_predicates_for_pfam("PF99999", "x", "mercury transporter") == []


def test_profile_name_resolves_to_the_installed_release_family():
    assert pfam_map.get_predicates_for_pfam("HTH_65") == pfam_map.get_predicates_for_pfam("PF21320")
    assert pfam_map.get_predicates_for_pfam("NiFeSe_Hases") == pfam_map.get_predicates_for_pfam("PF00374.26")


def _families_with(pred):
    return {MAP[acc][0] for acc, p, _ in PAIRS if p == pred}


@pytest.mark.parametrize(("family", "absent"), [
    ("Rav1p_C", {"hydrogenase", "nife_hydrogenase", "nife_group4"}),   # yeast V-ATPase assembly
    ("PFO_beta_C", {"toxin", "hydrolase"}),                            # pyruvate:ferredoxin oxidoreductase
    ("MMPL", {"mercury_resistance", "heavy_metal_resistance"}),        # lipid transporter
    ("Cass2", {"crispr_associated"}),                                  # integron effector-binding protein
    ("TRAM_LAG1_CLN8", {"t3ss_component"}),
    ("TrbI", {"t6ss_component"}),
    ("T6SS_VasE", {"t4ss_component"}),
    ("Thioredoxin", {"arsenic_resistance"}),
    ("2-Hacid_dh", {"2fe2s", "fad_binding"}),
    ("Cas1_AcylT", {"defense_system", "crispr_associated"}),            # fungal Cas1p
    ("NIF", {"nitrogen_fixation"}),                                     # NLI-interacting-factor phosphatase
    ("Oxidored_q2", {"cofactor_biosynthesis"}),                         # uses ubiquinone
    ("GreA_GreB", {"translation_factor"}),                              # transcription elongation
    ("RasGAP", {"gtpase"}),
    ("Transpeptidase", {"protease"}),
    ("PIN", {"esterase"}),                                              # nuclease with a partial EC 3.1
    ("TACC_C", {"wd40_repeat", "lrr_repeat"}),                          # GO protein binding is function-agnostic
    ("ATP-cone", {"kinase"}),                                           # mostly ribonucleotide reductases
    ("HATPase_c", {"sensor_kinase"}),                                   # also gyrase, Hsp90, MutL
    ("GATase", {"lyase"}),
    ("Fer4", {"ferredoxin"}),                                           # 4Fe-4S binding domain, many enzymes
    ("HNOB", {"cytochrome"}),                                           # heme binding without electron carriage
    ("Lipid_desat", {"monooxygenase"}),                                 # EC 1.14.19 desaturase
    ("TauD", {"monooxygenase"}),                                        # EC 1.14.11 2OG-dependent dioxygenase
    ("MCRA", {"mcr_complex", "methanogenesis"}),                        # myosin cross-reactive antigen
    # Families named after some of their members (SUBSET_NAMES) and name tokens
    # that state another subunit's cofactor.
    ("NiFeSe_Hases", {"nifese_hydrogenase", "selenium_metabolism"}),    # all [NiFe] large subunits
    ("SelR", {"selenium_metabolism", "trna_synthetase"}),               # MsrB
    ("SepSecS", {"selenium_metabolism"}),                               # also SepCysS
    ("MbhD", {"nife_hydrogenase", "nife_group4", "mbh_hydrogenase"}),   # Mrp/Mnh antiporter subunit A
    ("MbhE", {"nife_hydrogenase", "nife_group4", "mbh_hydrogenase"}),
    ("CbiA", {"cobalamin_biosynthesis", "chromosome_partitioning", "cell_division"}),  # CobQ/CobB/MinD/ParA
    ("cobW", {"cobalamin_biosynthesis", "nickel_binding"}),             # CobW/HypB/UreG/zinc chaperones
    ("CobW_C", {"cobalamin_biosynthesis", "cobalamin_binding", "cobalt_binding"}),  # COG0523 zinc chaperones
    ("CobW_C_DAAF9", {"cobalamin_biosynthesis"}),                       # dynein assembly factor DAAF9
    ("CbiQ", {"cobalamin_biosynthesis", "cobalt_binding"}),             # ECF transporter T component
    ("CobN-Mg_chel", {"cobalamin_biosynthesis"}),                       # also Mg-chelatase H
    ("TraB_PrgY_gumN", {"conjugation"}),                                # TIKI/YbaP metalloprotease-like
    ("DegT_DnrJ_EryC1", {"secondary_metabolism"}),                      # mostly LPS/glycan aminotransferases
    ("NiFe_hyd_SSU_C", {"nickel_binding", "oxygen_binding"}),           # [NiFe] small subunit
    ("NiFe_hyd_3_EhaA", {"nickel_binding"}),                            # Eha membrane subunits
    ("Prd_NiFe_hyd_3_EhaI", {"nickel_binding"}),
    ("Prd_NiFe_hyd_3_EhaK", {"nickel_binding"}),
    ("NiFe-hyd_HybE", {"nickel_binding"}),                              # large-subunit chaperone
    ("DsrC", {"sulfate_reduction"}),                                    # also TusE/YccK and RspA
    # Families named after the protein they act on (TARGET_NAMES).
    ("DRAT", {"nitrogenase", "oxidoreductase", "reductase"}),           # DraT ADP-ribosylates NifH
])
def test_known_false_claims_stay_out(family, absent):
    accs = [acc for acc, (name, _, _) in MAP.items() if name == family]
    for acc in accs:
        assert not absent & set(MAP[acc][1]), (family, MAP[acc][1])


def test_subset_names_take_no_wording_evidence():
    desc = "Nickel-dependent hydrogenase"
    # The wording matches; judge() refuses it for families named after some members.
    assert text_evidence("selenium_metabolism", "NiFeSe_Hases", desc) == "NiFeSe"
    assert judge("selenium_metabolism", "NiFeSe_Hases", desc, set(), {}) is None
    assert judge("nifese_hydrogenase", "NiFeSe_Hases", desc, set(), {}) is None
    # Family-wide properties keep their wording and GO evidence.
    assert judge("hydrogen_metabolism", "NiFeSe_Hases", desc, set(), {}) == "text:NiFe"
    assert judge("nickel_binding", "NiFeSe_Hases", desc, {"GO:0016151"}, {}) == "go:GO:0016151"
    assert judge("selenium_metabolism", "SelA", "L-seryl-tRNA selenium transferase", set(), {}) == "text:SelA"
    assert {pred for _, pred in SUBSET_NAMES} <= set(PREDICATE_BY_ID)
    worded = {(MAP[acc][0], pred) for acc, pred, ev in PAIRS if ev.startswith(("text:", "enzyme:"))}
    assert not worded & set(SUBSET_NAMES)


def test_target_names_take_no_wording_evidence():
    # DRAT's description names its target, the nitrogenase Fe protein; the family is a transferase.
    desc = "Dinitrogenase reductase ADP-ribosyltransferase (DRAT)"
    assert text_evidence("nitrogenase", "DRAT", desc) is not None
    assert judge("nitrogenase", "DRAT", desc, set(), {}) is None
    assert judge("reductase", "DRAT", desc, set(), {}) is None
    # Process-level GO evidence still counts.
    assert judge("nitrogen_fixation", "DRAT", desc, {"GO:0009399"}, {}) == "go:GO:0009399"
    assert {pred for _, pred in TARGET_NAMES} <= set(PREDICATE_BY_ID)
    worded = {(MAP[acc][0], pred) for acc, pred, ev in PAIRS if ev.startswith(("text:", "enzyme:"))}
    assert not worded & set(TARGET_NAMES)


def test_dsrc_family_keeps_sulfur_relay_claims_only():
    # PF04358 holds DsrC and the tRNA-thiolation relay TusE; only family-wide claims ship.
    (acc,) = [acc for acc, (name, _, _) in MAP.items() if name == "DsrC"]
    assert "sulfur_metabolism" in MAP[acc][1]
    assert "sulfate_reduction" not in MAP[acc][1]


def test_subset_names_leave_kegg_wording_alone():
    # KEGG passes KO symbols as the name; K02234's symbol 'cobW' equals a Pfam family name.
    assert text_evidence("cobalamin_biosynthesis", "cobW", "cobalamin biosynthesis protein CobW") is not None


def test_nifese_hases_carries_only_family_wide_hydrogenase_claims():
    # PF00374 covers every [NiFe] large subunit; selenocysteine enzymes are a few members.
    (acc,) = [acc for acc, (name, _, _) in MAP.items() if name == "NiFeSe_Hases"]
    preds = set(MAP[acc][1])
    assert {"hydrogenase", "nife_hydrogenase", "nickel_binding"} <= preds
    assert not {"selenium_metabolism", "nifese_hydrogenase"} & preds


def test_substring_collisions_do_not_label_families():
    # 'mer[ABCDE]' once matched isomerase/polymerase; 'AraC' matched 'characterized'.
    assert len(_families_with("mercury_resistance")) < 10
    assert all("Mer" in name for name in _families_with("mercury_resistance"))
    assert not any("haracteri" in SNAPSHOT[acc]["desc"] for acc, p, _ in PAIRS if p == "arac_family")
