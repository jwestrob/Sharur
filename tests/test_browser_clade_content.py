"""Clade gene content: spectrum, core, signature and absent families."""

import threading
from types import SimpleNamespace

from sharur.browser import routes_matrix as mx
from sharur.browser.catalog import Catalog, Genome
from sharur.browser.clade_content import clade_content


def _setup():
    catalog = Catalog()
    clades = ["A"] * 5 + ["B"] * 10
    for i, clade in enumerate(clades):
        g = Genome(i, f"g{i}", {"phylum": clade}, 95.0, 1.0, 1, 10000, 1000, 1000, 900, 300)
        catalog.genomes.append(g)
        catalog.by_bin[g.bin_id] = g
    pairs = [(f"g{i}", "K_SIG", 1.0) for i in range(5)]                 # every A genome, no B genome
    pairs += [(f"g{i}", "K_ALL", 1.0) for i in range(15)]                # everywhere
    pairs += [(f"g{i}", "K_LOST", 1.0) for i in range(5, 15)]            # every B genome, no A genome
    pairs += [("g0", "K_RARE", 1.0)]
    ids, bins, feats, values = mx._sparse(catalog, pairs)
    fs = mx.FeatureSet("ko", ids, ids, [f"{k} name" for k in ids], bins, feats, values)
    ctx = SimpleNamespace(catalog=catalog, lock=threading.Lock(), _assembly_kb={})
    return ctx, SimpleNamespace(get=lambda kind: fs), catalog.clade("phylum", "A")


def test_clade_content_spectrum_core_and_lists():
    ctx, sets, genomes = _setup()
    cc = clade_content(ctx, sets, genomes)
    assert cc["genomes"] == 5 and cc["rest"] == 10 and cc["observed"] == 3
    assert cc["every"] == 2 and cc["rare"] == 0
    assert cc["spectrum"][-1]["count"] == 2 and cc["spectrum"][1]["count"] == 1   # K_RARE in 1 of 5 = 20%
    assert [f["id"] for f in cc["signature"]] == ["K_SIG"]
    assert [f["id"] for f in cc["absent"]] == ["K_LOST"]
    assert cc["consistent"] == 2                                                     # K_SIG and K_ALL


def test_small_clades_are_skipped():
    ctx, sets, genomes = _setup()
    assert clade_content(ctx, sets, genomes[:2]) is None
