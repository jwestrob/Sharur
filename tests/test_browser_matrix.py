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


def test_cluster_order_groups_similar_rows_and_orders_leaves_optimally():
    import numpy as np
    from scipy.spatial.distance import pdist, squareform

    from sharur.browser.routes_matrix import _jaccard_order

    rng = np.random.default_rng(0)
    blocks = [rng.random(30) < 0.5 for _ in range(4)]
    rows = []
    for i in range(24):                       # four profiles, interleaved, with one flipped bit each
        row = blocks[i % 4].copy()
        row[rng.integers(30)] ^= True
        rows.append(row)
    presence = np.array(rows)
    order = _jaccard_order(presence)
    assert sorted(order) == list(range(24))
    labels = [i % 4 for i in order]
    assert sum(a != b for a, b in zip(labels, labels[1:])) == 3     # each profile is one contiguous run
    d = squareform(pdist(presence, "jaccard"))
    adjacent = sum(d[a, b] for a, b in zip(order, order[1:]))
    from scipy.cluster.hierarchy import leaves_list, linkage
    plain = leaves_list(linkage(pdist(presence, "jaccard"), "average"))
    assert adjacent <= sum(d[a, b] for a, b in zip(plain, plain[1:])) + 1e-12


def _mixed_sets():
    """KO set from _ctx plus a curated-system set and a module-completeness set."""
    ctx, ko = _ctx()
    sys_pairs = [(f"g{i}", "secretion:Archaeal-T4P", 1.0) for i in range(5)]
    ids, bins, feats, values = mx._sparse(ctx.catalog, sys_pairs)
    system = mx.FeatureSet("system", ids, ["Archaeal-T4P"], ["Archaeal-T4P (secretion system)"], bins, feats, values)
    mod_pairs = [(f"g{i}", "M00052", v) for i, v in enumerate([1.0, 1.0, 0.8, 0.5, 0.5, 0.0, 1.0, 1.0, 1.0, 1.0])]
    ids, bins, feats, values = mx._sparse(ctx.catalog, mod_pairs)
    module = mx.FeatureSet("module", ids, ids, ["M00052 name"], bins, feats, values, threshold=mx.MODULE_PRESENT)
    by_kind = {"ko": ko, "system": system, "module": module}

    class Sets:
        def get(self, kind):
            if kind not in by_kind:
                raise LookupError(kind)
            return by_kind[kind]
    return ctx, Sets()


def test_features_resolve_across_kinds_and_loose_spelling():
    _, sets = _mixed_sets()
    assert mx.resolve_feature(sets, "system", "M00052")[0].kind == "module"         # other kind
    assert mx.resolve_feature(sets, "system", "secretion:Archaeal_T4P")[1] == "secretion:Archaeal-T4P"
    assert mx.resolve_feature(sets, "ko", "system:secretion:Archaeal-T4P")[0].kind == "system"
    assert mx.resolve_feature(sets, "ko", "module:K1") is None                      # explicit kind is binding
    assert mx.resolve_feature(sets, "ko", "nothing") is None


def test_mixed_columns_keep_each_kinds_presence_threshold():
    ctx, sets = _mixed_sets()
    fs = mx.merge_columns([mx.resolve_feature(sets, "system", t) for t in ("secretion:Archaeal-T4P", "M00052", "K2")])
    assert fs.kind == "mixed" and fs.kinds == ["system", "module", "ko"]
    a = mx.Group("all", "All", None, list(ctx.catalog.genomes))
    out = mx.build_matrix(ctx, fs, [a], features=list(fs.ids))
    t4p, m52, k2 = out["features"]
    assert [t4p["kind"], m52["kind"], k2["kind"]] == ["system", "module", "ko"]
    assert t4p["groups"][0]["carriers"] == 5
    assert m52["groups"][0]["carriers"] == 7 and m52["threshold"] == pytest.approx(mx.MODULE_PRESENT)  # 0.5 is absent
    assert k2["groups"][0]["carriers"] == 4
    assert out["tested"]                                                    # a KO column is present
    assert "p" not in t4p["groups"][0] and "p" not in m52["groups"][0]     # multi-gene kinds skip the test


def test_categories_resolve_by_id_or_name_and_use_curated_mappings():
    assert mx.resolve_category("ribosomal_protein") == mx.resolve_category("Ribosomal protein") == "ribosomal_protein"
    assert mx.resolve_category("no such label") is None
    pfams = mx.category_members("pfam", "ribosomal_protein")
    assert {"PF00163", "Ribosomal_S4"} <= pfams              # by accession, and by the name some annotations store
    assert mx.category_members("function", "translation") >= {"translation", "ribosomal_protein"}
    assert mx.category_members("system", "ribosomal_protein") == set()


def test_ribosomal_pfam_panel_holds_ribosomal_protein_families():
    panel = mx.resolve_panel("Ribosomal proteins (Pfam)")
    assert panel is mx.resolve_panel("ribosomal_proteins_pfam") and panel.kind == "pfam"
    assert mx.resolve_panel("ribosomal_protein") is None                 # the function label stays a category
    rows = [line.split("\t") for line in (mx.PANEL_DIR / "ribosomal_proteins_pfam.tsv").read_text().splitlines()
            if line and not line.startswith("#")]
    accessions = [r[0] for r in rows if r[0][:2] == "PF"]
    assert len({r[0] for r in rows}) == len(rows) and all(len(r) in (4, 5) for r in rows)
    assert {r[0] for r in rows} - set(accessions) == {"ko:K02986"} and panel.extra == (("ko", "K02986"),)
    assert panel.proteins["PF03947"] == panel.proteins["Ribosomal_L2_C"] == panel.proteins["Ribosomal_L2"] == "uL2"
    # universal names carry the traditional ones, whose numbers differ
    assert panel.describe("PF01247") == "eL33 (L35Ae)" and panel.describe("PF00252") == "uL16 (L10e)"
    assert panel.describe("K02986") == "uS4 (S4p)" and panel.describe("PF29996") == "aL49"
    # the panel refines the curated Pfam category: its families sit inside it, plus RS4NT (the eS4
    # N-terminal domain, named "RS4NT (NUC023) domain"), and the category's enzymes, factors and
    # shared domains sit outside the panel
    assert set(accessions) - mx.category_members("pfam", "ribosomal_protein") == {"PF08071"}
    outside = {"PF00467": "KOW", "PF00575": "S1", "PF06325": "PrmA", "PF00398": "RrnaAD", "PF07091": "FmrO",
               "PF02482": "Ribosomal_S30AE", "PF02410": "RsfS", "PF01926": "MMR_HSR1", "PF14404": "Strep_pep"}
    assert not (set(outside) | set(outside.values())) & set(panel.proteins)


def test_keyword_and_category_columns_are_ordered_most_widespread_first():
    ctx, fs = _ctx()
    assert mx.keyword_members(fs, "K3 NAME") == {"K3"} and mx.keyword_members(fs, "name") == {"K1", "K2", "K3"}
    words = mx.FeatureSet("ko", ["Ka", "Kb"], ["Ka", "Kb"], ["nonribosomal peptide synthetase", "ribosomal-protein-serine acetyltransferase"],
                          np.zeros(0, np.int32), np.zeros(0, np.int32), np.zeros(0, np.float32))
    assert mx.keyword_members(words, "ribosomal") == {"Kb"}            # word starts only
    a = mx.Group("all", "All", None, list(ctx.catalog.genomes))
    out = mx.build_matrix(ctx, fs, [a], candidates={"K3", "K1", "K9"}, pick="category", selection="category X")
    assert [f["id"] for f in out["features"]] == ["K1", "K3"]          # 6 then 3 carriers; K9 occurs nowhere
    assert out["pick"] == "category" and out["selection"] == "category X"
