"""Clade-level KEGG module completeness aggregation."""

import numpy as np
import pytest

from sharur.browser.catalog import Catalog, Genome


class _Module:
    def __init__(self, module, name):
        self.module, self.name, self.module_class = module, name, "Pathway modules; Energy; Test"


def _catalog():
    catalog = Catalog()
    for i, clade in enumerate(["A", "A", "A", "B", "B"]):
        g = Genome(i, f"g{i}", {"phylum": clade}, None, None, 1, 1000, 1000, 10, 5, 100)
        catalog.genomes.append(g)
        catalog.by_bin[g.bin_id] = g
    catalog.module_ids = ["M1", "M2", "M3"]
    catalog.modules = {m: _Module(m, f"module {m}") for m in catalog.module_ids}
    # rows: genomes; columns: M1 (complete in A), M2 (partial everywhere), M3 (only in B)
    catalog.module_completeness = np.array([
        [1.0, 0.5, 0.0],
        [0.8, 0.5, 0.0],
        [1.0, 0.25, 0.0],
        [0.0, 0.5, 1.0],
        [0.2, 0.5, 1.0],
    ], dtype=np.float32)
    return catalog


def test_median_mode_keeps_half_complete_modules():
    catalog = _catalog()
    rows = {r["module"]: r for r in catalog.clade_modules(catalog.clade("phylum", "A"))}
    assert set(rows) == {"M1", "M2"}
    assert rows["M1"]["median"] == pytest.approx(1.0)
    assert rows["M1"]["prevalence"] == pytest.approx(1.0)       # all three >= 0.75
    assert rows["M1"]["rest_prevalence"] == pytest.approx(0.0)
    assert rows["M2"]["median"] == pytest.approx(0.5) and rows["M2"]["prevalence"] == 0
    assert rows["M1"]["class"] == "Test"


def test_present_mode_and_distinctive():
    catalog = _catalog()
    clade_b = catalog.clade("phylum", "B")
    present = {r["module"] for r in catalog.clade_modules(clade_b, mode="present")}
    assert present == {"M1", "M2", "M3"}  # M1 has a step in 1 of 2 genomes
    distinctive = catalog.distinctive_modules(clade_b)
    assert [r["module"] for r in distinctive] == ["M3"]
    assert catalog.distinctive_modules(catalog.genomes) == []  # whole dataset has no "rest"


def test_missing_matrix_returns_nothing():
    catalog = _catalog()
    catalog.module_completeness = None
    assert catalog.clade_modules(catalog.genomes) == [] and catalog.distinctive_modules(catalog.genomes[:2]) == []
