"""Presence/absence matrix: statistics, the gene-call guard, feature picks and ordering."""

import threading
from types import SimpleNamespace

import numpy as np
import pytest

from sharur.browser import routes_matrix as mx
from sharur.browser.catalog import Catalog, Genome


def test_poisson_binomial_matches_binomial_for_equal_probabilities():
    from scipy.stats import binom

    cdf = mx.poisson_binomial_cdf(np.full(12, 0.8))
    assert cdf == pytest.approx(binom.cdf(np.arange(13), 12, 0.8))
    assert cdf[-1] == pytest.approx(1.0)


def test_poisson_binomial_mixed_probabilities():
    # two trials, p = 0.5 and 0.9: P(0) = 0.05, P(<=1) = 0.05 + 0.5*0.9 + 0.5*0.1 = 0.55
    assert mx.poisson_binomial_cdf(np.array([0.5, 0.9])) == pytest.approx([0.05, 0.55, 1.0])


def test_bh_qvalues_treat_untested_as_p_one():
    q = mx.bh_qvalues(np.array([0.01, 0.04, 0.03]), 3)
    assert q == pytest.approx([0.03, 0.04, 0.04])
    # the same p-values among 10 tests (7 untested) give larger q
    assert mx.bh_qvalues(np.array([0.01]), 10) == pytest.approx([0.1])


def _genome(i, clade, completeness, proteins=1000, sub="x"):
    return Genome(i, f"g{i}", {"phylum": clade, "class": f"{clade}{sub}"}, completeness, 1.0, 1, 1_000_000,
                  1000, proteins, proteins // 2, 300)


def test_protein_deficits_prefer_assembly_density():
    genomes = [_genome(i, "A", 90.0) for i in range(6)] + [_genome(6, "A", 90.0, proteins=20)]
    assert mx.protein_deficits(genomes, genomes) == {"g6"}             # completeness proxy
    assert mx.protein_deficits(genomes, genomes, {"g6": 30.0}) == set()  # small assembly: 20 proteins on 30 kb
    assert mx.protein_deficits(genomes, genomes, {"g0": 4000.0}) == {"g0", "g6"}


def _ctx():
    """10 genomes: A (6, complete) and B (4); KO K1 in all of A, K2 only in B, K3 in half of A."""
    catalog = Catalog()
    clades = ["A"] * 6 + ["B"] * 4
    for i, clade in enumerate(clades):
        g = _genome(i, clade, 95.0, sub="1" if i % 2 else "2")
        catalog.genomes.append(g)
        catalog.by_bin[g.bin_id] = g
    catalog.ready.set()
    pairs = [(f"g{i}", "K1", 1.0) for i in range(6)] + [(f"g{i}", "K2", 2.0) for i in range(6, 10)] + \
        [(f"g{i}", "K3", 1.0) for i in (0, 1, 2)]
    ids, bins, feats, values = mx._sparse(catalog, pairs)
    fs = mx.FeatureSet("ko", ids, ids, [f"{k} name" for k in ids], bins, feats, values)
    ctx = SimpleNamespace(catalog=catalog, lock=threading.Lock(), store=None)
    return ctx, fs


def test_two_groups_rank_differences_with_fisher():
    ctx, fs = _ctx()
    a = mx.Group("phylum:A", "A", "phylum", ctx.catalog.clade("phylum", "A"))
    b = mx.Group("phylum:B", "B", "phylum", ctx.catalog.clade("phylum", "B"))
    out = mx.build_matrix(ctx, fs, [a, b], n=3)
    assert out["pick"] == "differential"
    ids = [f["id"] for f in out["features"]]
    assert ids[0] == "K1" and ids[-1] == "K2"      # A-enriched first, B-enriched last
    k1 = out["features"][0]
    assert k1["delta"] == pytest.approx(1.0) and k1["fisher_p"] == pytest.approx(1 / 210)
    assert out["values"].shape == (10, 3)
    assert [g["group"] for g in out["genomes"]] == [0] * 6 + [1] * 4


def test_single_group_variable_pick_and_completeness_test():
    ctx, fs = _ctx()
    a = mx.Group("phylum:A", "A", "phylum", ctx.catalog.clade("phylum", "A"))
    out = mx.build_matrix(ctx, fs, [a], n=5)
    assert out["pick"] == "variable" and [f["id"] for f in out["features"]] == ["K3"]
    k3 = out["features"][0]["groups"][0]
    # 3 of 6 genomes at 95% completeness: far below the 5.7 expected
    assert k3["expected"] == pytest.approx(5.7) and k3["p"] < 0.01 and k3["q"] < 0.05


def test_absences_pick_needs_completeness():
    ctx, fs = _ctx()
    for g in ctx.catalog.genomes:
        g.completeness = None
    a = mx.Group("phylum:A", "A", "phylum", ctx.catalog.clade("phylum", "A"))
    assert "import-quality" in mx.build_matrix(ctx, fs, [a], pick="absences")["empty"]


def test_explicit_features_and_completeness_filter():
    ctx, fs = _ctx()
    ctx.catalog.genomes[0].completeness = 40.0
    a = mx.Group("all", "All", None, list(ctx.catalog.genomes))
    out = mx.build_matrix(ctx, fs, [a], features=["K2", "K9"], min_completeness=50)
    assert [f["id"] for f in out["features"]] == ["K2"]
    assert any("K9" in note for note in out["notes"]) and any("1 genomes below" in note for note in out["notes"])
    assert len(out["genomes"]) == 9
    tsv = mx._tsv(out)
    assert tsv.splitlines()[0].split("\t")[-1] == "K2" and len(tsv.splitlines()) == 10
