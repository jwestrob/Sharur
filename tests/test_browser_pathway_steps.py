"""Pathway step views: per-genome step table, clade gaps, one-step-short genomes, the diagram."""

import xml.etree.ElementTree as ET

from sharur.browser import pathway_steps as ps
from sharur.browser.catalog import Catalog, Genome
from sharur.modules import ModuleDefinition

# step 1: K00001 or K00002; step 2: K00003 + K00004 (K00005 optional); step 3: K00006; step 4: no KO (left out)
DEF = ModuleDefinition("M90001", "Test pathway", "Pathway modules; Test", "(K00001,K00002) K00003+K00004-K00005 K00006 --")
REF = ModuleDefinition("M90002", "Uses a module", "Pathway modules; Test", "K00007 M90001")


def _catalog():
    catalog = Catalog()
    kos = {
        "g0": {"K00001", "K00003", "K00004", "K00006"},          # complete
        "g1": {"K00002", "K00003", "K00004", "K00006", "K00005"},    # complete (optional K00005 too)
        "g2": {"K00001", "K00003", "K00006"},                # one step short: step 2 (needs K00004)
        "g3": {"K00001", "K00003", "K00004"},                # one step short: step 3 (needs K00006)
        "g4": {"K00001", "K00003", "K00004"},                # one step short: step 3
        "g5": set(),                             # nothing
    }
    for i, (b, present) in enumerate(kos.items()):
        g = Genome(i, b, {"order": "A" if i < 3 else "B"}, 90.0 + i, 1.0, 1, 1000, 1000, 100, 90, 300)
        catalog.genomes.append(g)
        catalog.by_bin[b] = g
        if present:
            catalog.ko_sets[b] = present
    catalog.modules = {DEF.module: DEF, REF.module: REF}
    return catalog


def test_step_table_and_shares():
    catalog = _catalog()
    table = ps.step_table(catalog, DEF)
    assert table.steps == [1, 2, 3]                       # the KO-less step is left out
    assert table.complete.sum(axis=0).tolist() == [5, 4, 3]
    assert table.relevant.tolist() == [True] * 5 + [False]
    assert ps.step_shares(table, catalog.genomes[:3]) == [1.0, 2 / 3, 1.0]
    assert ps.ko_shares(catalog, {"K00001", "K00005"}, catalog.genomes) == {"K00001": 4 / 6, "K00005": 1 / 6}


def test_clade_rows_mark_the_usual_gap():
    catalog = _catalog()
    table = ps.step_table(catalog, DEF)
    rows = {r["clade"]: r for r in ps.clade_rows(catalog, table, catalog.genomes, "order")}
    assert rows["A"]["shares"] == [1.0, 2 / 3, 1.0] and rows["A"]["gap"] is None
    assert rows["B"]["shares"] == [2 / 3, 2 / 3, 0.0] and rows["B"]["gap"] == 2   # step 3 is B's gap
    assert rows["B"]["complete"] == 0.0


def test_one_step_short_groups():
    catalog = _catalog()
    table = ps.step_table(catalog, DEF)
    groups = {g["step"]: g for g in ps.nearly_complete(catalog, table, catalog.genomes)}
    assert groups[3]["count"] == 2 and groups[3]["missing"] == [("K00006", 2)]
    assert groups[2]["count"] == 1 and groups[2]["missing"] == [("K00004", 1)]
    assert [g.bin_id for g in groups[3]["genomes"]] == ["g4", "g3"]      # most complete first


def test_nested_module_counts_and_diagram_svg():
    catalog = _catalog()
    assert ps._all_kos(REF, catalog.modules) == {"K00001", "K00002", "K00003", "K00004", "K00005", "K00006", "K00007"}
    table = ps.step_table(catalog, REF)
    svg = ps.diagram(DEF, shares={"K00001": 1.0}, step_complete=[1.0, 0.5, 0.25], steps=[1, 2, 3],
                     names=lambda ko: ("sym" + ko, "a name") if ko == "K00001" else None, ko_url=lambda ko: "/search?q=" + ko)
    root = ET.fromstring(str(svg))
    assert root.get("class") == "pw-diagram"
    text = str(svg)
    assert "symK00001" in text and 'class="pw-alt"' in text and "pw-optional" in text and "Step 4" in text
    heat = ps.heatmap(ps.clade_rows(catalog, ps.step_table(catalog, DEF), catalog.genomes, "order"), [1, 2, 3],
                      "order", lambda *p: "/" + "/".join(p))
    ET.fromstring(str(heat))
    assert table.n_steps == 2
