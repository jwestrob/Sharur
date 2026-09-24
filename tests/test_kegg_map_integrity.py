"""Every pair of the local KEGG build is supported by KEGG's own information.

``sharur setup-kegg`` builds ``kegg_predicates.tsv`` with one evidence item per
pair; this suite re-verifies each item against the build's
``kegg_evidence_snapshot.tsv`` (KEGG's symbols, name, BRITE placements and
module memberships), the shipped ``kegg_hyddb_snapshot.tsv`` and
``kegg_swissprot_consensus.tsv``, and the current rules. It skips when no local
build exists. Set ``SHARUR_SWISSPROT``, ``SHARUR_SWISSPROT_KEGG`` and
``SHARUR_GO_OBO`` to recount the Swiss-Prot consensus from reviewed proteins.
"""

import os
from pathlib import Path

import pytest

from sharur.predicates.mappings import kegg_map
from sharur.predicates.mappings.kegg_build import load_swissprot_consensus
from sharur.predicates.mappings.kegg_evidence import (
    PATH_SEP,
    brite_evidence,
    hyddb_evidence,
    ko_definition,
    ko_text_evidence,
    load_brite_rules,
    load_hyddb_snapshot,
    load_module_rules,
    module_evidence,
    names_hydrogenase,
)
from sharur.predicates.mappings.kegg_map import (
    KEGG_EVIDENCE,
    KEGG_TO_PREDICATES,
    get_predicates_for_ec,
    parse_ec_numbers,
)
from sharur.predicates.mappings.pfam_evidence import (
    SWISSPROT_COVERAGE,
    SWISSPROT_EXCLUDED,
    SWISSPROT_MIN_LOWER_BOUND,
    wilson_lower_bound,
)
from sharur.predicates.mappings.swissprot_evidence import (
    go_ancestors,
    protein_kos,
    read_gene_ko_links,
    read_swissprot,
)


DATA = Path(kegg_map.__file__).with_name("data")
LOCAL = kegg_map.kegg_dir()

pytestmark = pytest.mark.skipif(LOCAL is None, reason="no local KEGG build; run `sharur setup-kegg`")


def _rows(path):
    if not path.exists():
        return []
    return [line.rstrip("\n").split("\t") for line in path.read_text().splitlines()
            if line and not line.startswith("#")]


def _placements(brite):
    out = []
    for item in filter(None, brite.split("|")):
        hierarchy, _, path = item.partition(":")
        out.append((hierarchy, tuple(path.split(PATH_SEP)) if path else ()))
    return out


def _memberships(modules):
    return [(m.rstrip("+"), m.endswith("+")) for m in filter(None, modules.split(","))]


SNAPSHOT = {ko: {"symbols": symbols, "name": name, "brite": _placements(brite),
                 "modules": _memberships(modules)}
            for ko, symbols, name, brite, modules in _rows((LOCAL or DATA) / "kegg_evidence_snapshot.tsv")}
PAIRS = [(ko, pred, ev) for ko, evidence in KEGG_EVIDENCE.items() for pred, ev in evidence.items()]
BRITE_RULES = load_brite_rules()
MODULE_RULES = load_module_rules()
HYDDB = load_hyddb_snapshot()
CONSENSUS = load_swissprot_consensus()


def test_map_and_snapshot_cover_the_same_kos():
    assert set(KEGG_TO_PREDICATES) == set(SNAPSHOT)


def test_every_pair_reverifies_against_the_snapshot():
    failures = []
    for ko, pred, ev in PAIRS:
        snap = SNAPSHOT[ko]
        kind, _, detail = ev.partition(":")
        if kind == "ec":
            ok = detail in parse_ec_numbers(snap["name"]) and pred in get_predicates_for_ec(detail)
        elif kind == "brite":
            found = brite_evidence(BRITE_RULES, snap["brite"], ko_definition(snap["symbols"], snap["name"]))
            ok = pred in found and any(ev == f"brite:{h} {PATH_SEP.join(labels)}" for h, labels in snap["brite"])
        elif kind == "module":
            module = detail
            memberships = [m for m in snap["modules"] if m[0] == module]
            ok = bool(memberships) and pred in module_evidence(MODULE_RULES, memberships)
        elif kind == "hyddb":
            ok = names_hydrogenase(snap["name"]) and hyddb_evidence(HYDDB.get(ko, {})).get(pred) == ev
        elif kind == "text":
            ok = ko_text_evidence(pred, snap["symbols"], snap["name"]) is not None
        elif kind == "swissprot":
            k, n = map(int, detail.split("/"))
            ok = pred not in SWISSPROT_EXCLUDED and n >= SWISSPROT_COVERAGE[0] \
                and k >= SWISSPROT_COVERAGE[1] * n and wilson_lower_bound(k, n) >= SWISSPROT_MIN_LOWER_BOUND \
                and CONSENSUS.get(ko, {}).get(pred) == ev
        else:
            ok = False
        if not ok:
            failures.append((ko, snap["symbols"], pred, ev))
    assert not failures, failures[:20]


RECOUNT = ("SHARUR_SWISSPROT", "SHARUR_SWISSPROT_KEGG", "SHARUR_GO_OBO")


@pytest.mark.skipif(not all(os.environ.get(v) for v in RECOUNT),
                    reason="set " + ", ".join(RECOUNT) + " to recount Swiss-Prot consensus")
def test_swissprot_counts_recount_from_reviewed_proteins():
    proteins = read_swissprot(Path(os.environ["SHARUR_SWISSPROT"]), go_ancestors(Path(os.environ["SHARUR_GO_OBO"])))
    gene_ko = read_gene_ko_links(Path(os.environ["SHARUR_SWISSPROT_KEGG"]))
    members: dict[str, list] = {}
    for protein in proteins:
        for ko in protein_kos(protein, gene_ko):
            members.setdefault(ko, []).append(protein.predicates)
    mismatched = []
    for ko, pred, ev in PAIRS:
        if ev.startswith("swissprot:"):
            group = members.get(ko, [])
            counted = f"swissprot:{sum(pred in p for p in group)}/{len(group)}"
            if counted != ev:
                mismatched.append((ko, pred, ev, counted))
    assert not mismatched, mismatched[:20]


def test_every_rule_is_used():
    brite_used = {(r.hierarchy, r.levels, r.definition and r.definition.pattern) for r in BRITE_RULES
                  if any(r.matches(h, labels, ko_definition(s["symbols"], s["name"]))
                         for s in SNAPSHOT.values() for h, labels in s["brite"])}
    unused = [(r.hierarchy, r.levels) for r in BRITE_RULES
              if (r.hierarchy, r.levels, r.definition and r.definition.pattern) not in brite_used]
    assert not unused
    modules_seen = {m for s in SNAPSHOT.values() for m, _ in s["modules"]}
    assert not set(MODULE_RULES) - modules_seen


def test_proposal_comments_start_with_the_kegg_symbol():
    mismatched = {}
    for ko, _, comment in _rows(DATA / "kegg_predicate_proposals.tsv"):
        if ko in SNAPSHOT and SNAPSHOT[ko]["symbols"]:
            symbols = {s.strip().lower() for s in SNAPSHOT[ko]["symbols"].split(",")}
            if comment.split(";")[0].strip().lower() not in symbols:
                mismatched[ko] = comment
    assert not mismatched


@pytest.mark.parametrize(("ko", "present", "absent"), [
    ("K03529", {"chromosome_partitioning"}, {"crispr_associated", "defense_system"}),  # smc
    ("K03595", {"gtp_binding"}, {"chromosome_partitioning"}),                           # era
    ("K01356", {"repressor"}, {"integrase", "mobile_element"}),                         # lexA
    ("K15830", {"formate_coupled", "nife_group4"}, {"ech_hydrogenase"}),                # hycE: HydDB Group 4a
    ("K15828", set(), {"nife_hydrogenase", "formate_coupled"}),                         # hycC: membrane subunit
    ("K14087", {"ech_hydrogenase"}, set()),                                             # echB
    ("K14090", {"ech_hydrogenase", "nife_group4", "h2_evolving"}, {"ferredoxin_coupled"}),  # echE: 4e 130/170
    ("K14106", {"nife_group4"}, {"ech_hydrogenase"}),                                   # ehaO: HydDB Group 4h
    ("K14126", {"nife_group3", "heterodisulfide_reductase_linked"}, {"mbh_hydrogenase", "nife_group4"}),  # mvhA
    ("K18016", {"mbh_hydrogenase", "nife_group4", "ferredoxin_coupled"}, set()),        # mbhL: 4d 33/39
    ("K00436", {"nad_coupled", "nife_group3"}, {"uptake_hydrogenase"}),                 # hoxH
    ("K06281", {"nife_group1", "uptake_hydrogenase"}, {"nife_group2"}),                 # hyaB, hybC
    ("K23549", {"nife_group2"}, {"uptake_hydrogenase", "h2_sensor"}),                   # hupV: 2a/2b/2c/2e
    ("K02588", {"nitrogenase", "nitrogen_fixation"}, set()),                            # nifH
    ("K02040", {"phosphate_transporter"}, {"sulfate_transporter"}),                     # pstS
    ("K03407", {"sensor_kinase", "chemotaxis"}, {"response_regulator"}),                # cheA
    ("K07659", {"response_regulator", "ompr_family"}, {"sensor_kinase"}),               # ompR
    ("K07642", {"sensor_kinase"}, {"ompr_family", "response_regulator"}),               # baeS
    ("K03194", {"t4ss_component", "conjugation"}, {"t6ss_component"}),                  # virB1
    ("K17836", {"beta_lactamase"}, {"aminoglycoside_resistance"}),                      # penP
    ("K00330", {"complex_subunit", "electron_transport"}, set()),                       # nuoA
    ("K00558", {"dna_methylase"}, set()),                                               # DNMT1 [EC:2.1.1.37]
    ("K19159", {"antitoxin"}, {"toxin"}),                                               # yefM
    ("K19158", {"toxin"}, {"antitoxin"}),                                               # yoeB
])
def test_corrected_entries(ko, present, absent):
    preds = set(KEGG_TO_PREDICATES[ko])
    assert present <= preds, (ko, present - preds)
    assert not absent & preds, (ko, absent & preds)


@pytest.mark.parametrize("ko", ["K15826", "K14129", "K14130"])
def test_non_hydrogenase_kos_carry_no_hydrogenase_claims(ko):
    assert not {"hydrogenase", "nife_hydrogenase", "ech_hydrogenase"} & set(KEGG_TO_PREDICATES.get(ko, ()))
