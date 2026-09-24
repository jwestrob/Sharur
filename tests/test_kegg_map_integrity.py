"""Every shipped KO -> predicate pair is supported by KEGG's own information.

``kegg_predicates.tsv`` records one evidence item per pair; this suite
re-verifies each item against ``kegg_evidence_snapshot.tsv`` (KEGG's symbols,
name, BRITE placements and module memberships at build time),
``kegg_hyddb_snapshot.tsv`` and the current rules. Swiss-Prot consensus pairs are
checked against the consensus thresholds; set ``SHARUR_SWISSPROT``,
``SHARUR_SWISSPROT_KEGG`` and ``SHARUR_GO_OBO`` to recount them. Rebuild with
``scripts/build_kegg_predicate_map.py`` after editing proposals or rules.
"""

import os
from pathlib import Path

import pytest

from sharur.predicates.mappings import kegg_map
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
    parse_module_definition,
)
from sharur.predicates.mappings.kegg_map import (
    EC_TO_PREDICATES,
    KEGG_EVIDENCE,
    KEGG_TO_PREDICATES,
    get_predicates_for_ec,
    get_predicates_for_kegg,
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
from sharur.predicates.vocabulary import ALL_PREDICATES, PREDICATE_BY_ID


DATA = Path(kegg_map.__file__).with_name("data")


def _rows(name):
    return [line.rstrip("\n").split("\t") for line in (DATA / name).read_text().splitlines()
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
            for ko, symbols, name, brite, modules in _rows("kegg_evidence_snapshot.tsv")}
PAIRS = [(ko, pred, ev) for ko, evidence in KEGG_EVIDENCE.items() for pred, ev in evidence.items()]
BRITE_RULES = load_brite_rules()
MODULE_RULES = load_module_rules()
HYDDB = load_hyddb_snapshot()


def test_map_and_snapshot_cover_the_same_kos():
    assert set(KEGG_TO_PREDICATES) == set(SNAPSHOT)


def test_vocabulary_ids_are_unique():
    ids = [p.predicate_id for p in ALL_PREDICATES]
    assert len(ids) == len(set(ids))


@pytest.mark.parametrize("table", [KEGG_TO_PREDICATES, EC_TO_PREDICATES])
def test_every_mapped_predicate_is_in_the_vocabulary(table):
    assert not {p for preds in table.values() for p in preds if p not in PREDICATE_BY_ID}


def test_rules_name_vocabulary_predicates():
    named = {p for r in BRITE_RULES for p in r.predicates} | {p for _, ps in MODULE_RULES.values() for p in ps}
    assert not named - set(PREDICATE_BY_ID)


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
                and k >= SWISSPROT_COVERAGE[1] * n and wilson_lower_bound(k, n) >= SWISSPROT_MIN_LOWER_BOUND
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
    for ko, _, comment in _rows("kegg_predicate_proposals.tsv"):
        if ko in SNAPSHOT and SNAPSHOT[ko]["symbols"]:
            symbols = {s.strip().lower() for s in SNAPSHOT[ko]["symbols"].split(",")}
            if comment.split(";")[0].strip().lower() not in symbols:
                mismatched[ko] = comment
    assert not mismatched


def test_runtime_uses_only_the_generated_map():
    assert not hasattr(kegg_map, "KEGG_PATTERNS")
    assert get_predicates_for_kegg("K99999", "hydrogenase membrane protein") == []
    assert set(get_predicates_for_kegg("K99999", "x [EC:2.7.1.1]")) == set(get_predicates_for_ec("2.7.1.1"))


def test_module_definition_parser():
    parsed = parse_module_definition("K00330+(K00331+K00332,K13380)+K00334 (K00844,K12407) K01803-(K1,K00002)")
    assert all(parsed[k] for k in ("K13380", "K00334", "K00002"))
    assert not any(parsed[k] for k in ("K00844", "K12407"))


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


@pytest.mark.parametrize(("ec", "present", "absent"), [
    ("2.1.2.1", {"transferase"}, {"methyltransferase"}),
    ("2.1.1.37", {"methyltransferase", "dna_methylase"}, set()),
    ("2.1.1.45", {"methyltransferase"}, {"dna_methylase"}),
    ("3.6.5.2", {"gtp_binding"}, {"atpase"}),
    ("3.6.1.1", {"hydrolase"}, {"atpase", "atp_binding"}),
    ("3.1.21.4", {"nuclease", "endonuclease"}, {"esterase"}),
    ("3.1.-.-", {"hydrolase"}, {"esterase"}),
    ("4.1.3.1", {"lyase"}, {"decarboxylase"}),
    ("7.6.2.1", {"atp_binding"}, {"gtp_binding"}),
    ("2.7.8.5", {"transferase"}, {"kinase"}),
    ("1.4.1.3", {"oxidoreductase"}, {"aminotransferase"}),
    ("2.6.1.1", {"aminotransferase", "plp_binding"}, set()),
    ("1.2.1.12", {"nad_binding"}, set()),                          # ENZYME 1.2.1: NAD(+)/NADP(+) acceptor
    ("1.3.8.7", {"flavin_binding"}, {"nad_binding"}),              # ENZYME 1.3.8: flavin acceptor
    ("1.14.12.10", {"dioxygenase", "nad_binding"}, {"monooxygenase"}),
    ("1.14.13.1", {"monooxygenase", "nad_binding"}, {"dioxygenase"}),
    ("1.14.11.2", {"dioxygenase"}, {"monooxygenase"}),
    ("1.14.19.1", {"oxygenase"}, {"monooxygenase", "dioxygenase"}),
])
def test_ec_classes(ec, present, absent):
    preds = set(get_predicates_for_ec(ec))
    assert present <= preds
    assert not absent & preds
