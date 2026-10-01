"""Search results: scoring, grouping, KO and EC lookups, near misses, highlighting."""

import threading
from types import SimpleNamespace

import pytest

from sharur.browser import routes_matrix as mx
from sharur.browser import search_page as sp
from sharur.browser.catalog import Catalog, Genome


class _Names:
    names = {"K00001": ("adh", "alcohol dehydrogenase [EC:1.1.1.1]"),
             "K00172": ("porC", "pyruvate ferredoxin oxidoreductase gamma subunit [EC:1.2.7.1]"),
             "K09999": ("xyz", "not in this dataset")}

    def get(self, ko):
        return self.names.get(ko)


class _Store:
    def execute(self, sql, params=()):
        return []


def _ctx():
    catalog = Catalog()
    for i, clade in enumerate(["Alpha", "Alpha", "Beta"]):
        g = Genome(i, f"g{i}", {"phylum": clade, "order": f"{clade}ales"}, 90.0, 1.0, 1, 1000, 1000, 10, 9, 300)
        catalog.genomes.append(g)
        catalog.by_bin[g.bin_id] = g
    catalog.domains = {"PF00005": {"accession": "PF00005", "name": "ABC_tran", "description": "ABC transporter",
                                   "proteins": 4, "genomes": 3, "hits": 4}}
    catalog.ready.set()
    catalog.notable = {"giants": []}
    ids, bins, feats, values = mx._sparse(catalog, [("g0", "K00001", 1.0), ("g1", "K00001", 1.0),
                                                    ("g2", "K00172", 1.0)])
    fs = mx.FeatureSet("ko", ids, ids, ids, bins, feats, values)
    return SimpleNamespace(catalog=catalog, lock=threading.Lock(), store=_Store(), ko_names=_Names(),
                           feature_sets=SimpleNamespace(get=lambda kind: fs))


def test_groups_best_match_and_tie_break():
    out = sp.run(_ctx(), "Alpha")
    assert [g["kind"] for g in out["groups"]] == ["taxon"]
    assert out["best"]["entry"].label == "Alpha"                 # the exact name beats "Alphaales"
    assert out["best"]["share"] == pytest.approx(2 / 3)


def test_ko_id_ec_number_and_missing_ko():
    ctx = _ctx()
    ko = sp.run(ctx, "K00001")
    assert ko["best"]["entry"].id == "K00001" and "1.1.1.1" in str(ko["best"]["text"])
    assert [e.id for e in sp.run(ctx, "EC 1.2.7.1")["groups"][0]["rows"]] == ["K00172"]
    assert [e.id for e in sp.run(ctx, "1.2.7.-")["groups"][0]["rows"]] == ["K00172"]
    missing = sp.run(ctx, "K09999")
    assert missing["missing_ko"]["symbols"] == "xyz" and not missing["groups"]


def test_near_misses_and_highlight():
    assert "abc_tran" in sp.run(_ctx(), "ABC_trn")["near"]
    assert str(sp.highlight("ABC <b> transporter", "abc")) == "<mark>ABC</mark> &lt;b&gt; transporter"


def test_rank_word_does_not_match_every_taxon_of_that_rank():
    assert not any(g["kind"] == "taxon" for g in sp.run(_ctx(), "order")["groups"])
