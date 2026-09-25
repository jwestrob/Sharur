"""Every shipped CAZy family -> predicate pair is supported by its recorded evidence.

``class:`` evidence is re-derived from CAZy's class definitions; ``swissprot:``
evidence is checked against the consensus thresholds, and recounted from the
reviewed proteins when ``SHARUR_SWISSPROT`` and ``SHARUR_GO_OBO`` are set.
"""

import os
import re
from pathlib import Path

import pytest

from sharur.predicates.mappings.cazy_map import (
    CAZY_CLASS_PREDICATES,
    CAZY_EVIDENCE,
    cazy_class,
    cazy_family,
    get_predicates_for_cazy,
)
from sharur.predicates.mappings.pfam_evidence import (
    SWISSPROT_COVERAGE,
    SWISSPROT_EXCLUDED,
    SWISSPROT_MIN_LOWER_BOUND,
    SWISSPROT_SINGLE,
    wilson_lower_bound,
)
from sharur.predicates.mappings.swissprot_evidence import go_ancestors, ortholog_support, read_swissprot
from sharur.predicates.vocabulary import PREDICATE_BY_ID, SYSTEM_PREDICATES


PAIRS = [(family, pred, ev) for family, evidence in CAZY_EVIDENCE.items() for pred, ev in evidence.items()]
SWISSPROT = re.compile(r"swissprot:(\d+)/(\d+)(?: KOs (\d+)/(\d+))? (?:single-domain (\d+)/(\d+)|co-domains excluded)")


def _consensus_ok(pred: str, ev: str) -> bool:
    m = SWISSPROT.fullmatch(ev)
    if not m or pred in SWISSPROT_EXCLUDED:
        return False
    k, n = int(m.group(1)), int(m.group(2))
    ok = n >= SWISSPROT_COVERAGE[0] and k >= SWISSPROT_COVERAGE[1] * n and wilson_lower_bound(k, n) >= SWISSPROT_MIN_LOWER_BOUND
    if m.group(5) is not None:
        ok = ok and int(m.group(6)) >= SWISSPROT_SINGLE[0] and int(m.group(5)) >= SWISSPROT_SINGLE[1] * int(m.group(6))
    return ok


def test_every_predicate_is_a_component_vocabulary_predicate():
    preds = {pred for _, pred, _ in PAIRS}
    assert preds <= set(PREDICATE_BY_ID)
    assert not preds & SYSTEM_PREDICATES


def test_every_pair_reverifies():
    failures = []
    for family, pred, ev in PAIRS:
        kind, _, detail = ev.partition(":")
        if kind == "class":
            ok = detail == cazy_class(family) and pred in CAZY_CLASS_PREDICATES[detail]
        elif kind == "swissprot":
            ok = cazy_class(family) != "CBM" and _consensus_ok(pred, ev)
        else:
            ok = False
        if not ok:
            failures.append((family, pred, ev))
    assert not failures, failures[:20]


def test_every_family_carries_its_class_predicates():
    for family in CAZY_EVIDENCE:
        assert set(CAZY_CLASS_PREDICATES[cazy_class(family)]) <= set(get_predicates_for_cazy(family)), family


@pytest.mark.parametrize(("family", "absent"), [
    ("GH5", {"cellulase"}),                    # polyspecific family
    ("GH13", {"amylase"}),                     # alpha-amylase superfamily includes many other activities
    ("GH2", {"mannanase"}),
    ("GH20", {"chitinase"}),                   # beta-hexosaminidases
    ("CBM20", {"amylase"}),                    # non-catalytic module
    ("CBM5", {"chitinase"}),
    ("CE11", {"esterase"}),                    # LpxC, a de-N-acetylase
])
def test_known_false_claims_stay_out(family, absent):
    assert not absent & set(get_predicates_for_cazy(family))


@pytest.mark.parametrize(("family", "present"), [
    ("GH7", {"cellulase"}), ("GH11", {"xylanase"}), ("GH19", {"chitinase"}),
    ("GH22", {"lysozyme"}), ("AA9", {"lytic_polysaccharide_monooxygenase"}),
])
def test_supported_substrate_claims(family, present):
    assert present <= set(get_predicates_for_cazy(family))


def test_subfamilies_and_unknown_families_resolve():
    assert cazy_family("GH5_12") == "GH5"
    assert get_predicates_for_cazy("GH5_12") == get_predicates_for_cazy("GH5")
    assert set(get_predicates_for_cazy("GH999")) == set(CAZY_CLASS_PREDICATES["GH"])
    assert get_predicates_for_cazy("not_a_family") == []


@pytest.mark.skipif(not (os.environ.get("SHARUR_SWISSPROT") and os.environ.get("SHARUR_GO_OBO")),
                    reason="set SHARUR_SWISSPROT and SHARUR_GO_OBO to recount CAZy consensus")
def test_swissprot_counts_recount_from_reviewed_proteins():
    reviewed = read_swissprot(Path(os.environ["SHARUR_SWISSPROT"]), go_ancestors(Path(os.environ["SHARUR_GO_OBO"])))
    proteins = []
    for p in reviewed:
        catalytic = frozenset(cazy_family(f) for f in p.cazy if cazy_class(f) not in (None, "CBM"))
        if catalytic:
            proteins.append((catalytic, p.predicates))
    mismatched = []
    for family, pred, ev in PAIRS:
        m = SWISSPROT.fullmatch(ev)
        if not m:
            continue
        carriers = [(fams, preds) for fams, preds in proteins if family in fams]
        counted = [sum(pred in p for _, p in carriers), len(carriers)]
        k_kos, n_kos = ortholog_support([(p, frozenset()) for _, p in carriers], pred)
        if n_kos >= 2:
            counted += [k_kos, n_kos]
        if m.group(5) is not None:
            single = [p for fams, p in carriers if len(fams) == 1]
            counted += [sum(pred in p for p in single), len(single)]
        if [int(g) for g in m.groups() if g is not None] != counted:
            mismatched.append((family, pred, ev, counted))
    assert not mismatched, mismatched[:20]
