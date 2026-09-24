"""Every shipped Pfam -> predicate pair is supported by Pfam's own information.

``pfam_predicates.tsv`` records one evidence item per pair; this suite
re-verifies each item against ``pfam_evidence_snapshot.tsv`` (the family's Pfam
name/description, InterPro GO anchors, and ENZYME names at build time) using the
current evidence rules. Rebuild with ``scripts/build_pfam_predicate_map.py``
after editing proposals or evidence definitions.
"""

import re
from pathlib import Path

import pytest

from sharur.predicates.mappings import pfam_map
from sharur.predicates.mappings.kegg_map import get_predicates_for_ec
from sharur.predicates.mappings.pfam_evidence import EVIDENCE, text_evidence
from sharur.predicates.mappings.pfam_evidence_spec import HOMONYMS
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
        else:
            ok = False
        if not ok or (snap["name"], pred) in HOMONYMS:
            failures.append((acc, snap["name"], pred, ev))
    assert not failures, failures[:20]


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
])
def test_known_false_claims_stay_out(family, absent):
    accs = [acc for acc, (name, _, _) in MAP.items() if name == family]
    for acc in accs:
        assert not absent & set(MAP[acc][1]), (family, MAP[acc][1])


def test_substring_collisions_do_not_label_families():
    # 'mer[ABCDE]' once matched isomerase/polymerase; 'AraC' matched 'characterized'.
    assert len(_families_with("mercury_resistance")) < 10
    assert all("Mer" in name for name in _families_with("mercury_resistance"))
    assert not any("haracteri" in SNAPSHOT[acc]["desc"] for acc, p, _ in PAIRS if p == "arac_family")
