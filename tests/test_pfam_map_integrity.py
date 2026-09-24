"""Every shipped Pfam -> predicate pair is supported by Pfam's own information.

``pfam_predicates.tsv`` records one evidence item per pair; this suite
re-verifies each item against ``pfam_evidence_snapshot.tsv`` (the family's Pfam
name/description, InterPro GO anchors, and ENZYME names at build time) using the
current evidence rules. Swiss-Prot consensus pairs are checked against the
consensus thresholds here; set ``SHARUR_SWISSPROT`` (uniprot_sprot.dat.gz) and
``SHARUR_GO_OBO`` (go-basic.obo) to recount them from the reviewed proteins. Rebuild with ``scripts/build_pfam_predicate_map.py``
after editing proposals or evidence definitions.
"""

import importlib.util
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
    SWISSPROT_MIN_LOWER_BOUND,
    SWISSPROT_SINGLE,
    text_evidence,
    wilson_lower_bound,
)
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
SWISSPROT = re.compile(r"swissprot:(\d+)/(\d+) (?:single-domain (\d+)/(\d+)|co-domains excluded)")


def _meets_consensus(pred, detail):
    m = SWISSPROT.fullmatch(f"swissprot:{detail}")
    if not m or pred in SWISSPROT_EXCLUDED:
        return False
    k, n = int(m.group(1)), int(m.group(2))
    covered = n >= SWISSPROT_COVERAGE[0] and k >= SWISSPROT_COVERAGE[1] * n \
        and wilson_lower_bound(k, n) >= SWISSPROT_MIN_LOWER_BOUND
    if m.group(3) is None:
        return covered
    ks, ns = int(m.group(3)), int(m.group(4))
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
        if not ok or (snap["name"], pred) in HOMONYMS:
            failures.append((acc, snap["name"], pred, ev))
    assert not failures, failures[:20]


def test_swissprot_pairs_record_their_release():
    header = (DATA / "pfam_predicates.tsv").read_text().split("# accession", 1)[0]
    if any(ev.startswith("swissprot:") for _, _, ev in PAIRS):
        assert "# Swiss-Prot: " in header


@pytest.mark.skipif(not (os.environ.get("SHARUR_SWISSPROT") and os.environ.get("SHARUR_GO_OBO")),
                    reason="set SHARUR_SWISSPROT and SHARUR_GO_OBO to recount Swiss-Prot consensus")
def test_swissprot_counts_recount_from_reviewed_proteins():
    script = Path(__file__).resolve().parents[1] / "scripts/build_pfam_predicate_map.py"
    spec = importlib.util.spec_from_file_location("build_pfam_predicate_map", script)
    build = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(build)
    obo = Path(os.environ["SHARUR_GO_OBO"])
    build.read_go(Path(os.devnull), obo)
    proteins = build.read_swissprot(Path(os.environ["SHARUR_SWISSPROT"]), build.read_go.ancestors)
    mismatched = []
    for acc, pred, ev in PAIRS:
        m = SWISSPROT.fullmatch(ev)
        if not m:
            continue
        carriers = [(fams, preds) for fams, preds in proteins if acc in fams]
        counted = [sum(pred in p for _, p in carriers), len(carriers)]
        if m.group(3) is not None:
            single = [p for fams, p in carriers if len(fams) == 1]
            counted += [sum(pred in p for p in single), len(single)]
        if list(map(int, filter(None, m.groups()))) != counted:
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
