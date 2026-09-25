"""KEGG module completeness on synthetic module definitions (no KEGG data)."""

import pytest

from sharur.modules import (
    ModuleDefinition,
    evaluate,
    genome_modules,
    locus_modules,
    modules_markdown,
    parse,
)
from sharur.storage.duckdb_store import DuckDBStore


DEFINITIONS = {
    # three steps; step 2 has alternatives; step 3 is a nested sequence alternative
    "M90001": "K90001 (K90002,K90003) ((K90004 K90005),K90006)",
    # one complex step with an optional component
    "M90002": "K90011+K90012+K90013-K90014",
    # a step referring to another module and an undefined step
    "M90003": "K90021 M90002 --",
}


def _modules():
    return {m: ModuleDefinition(m, f"test module {m}", "Pathway modules; Test", d) for m, d in DEFINITIONS.items()}


@pytest.fixture
def kegg_dir(tmp_path):
    lines = ["# module\tname\tclass\tdefinition"] + [f"{m}\ttest module {m}\tPathway modules; Test\t{d}"
                                                    for m, d in DEFINITIONS.items()]
    (tmp_path / "kegg_modules.tsv").write_text("\n".join(lines) + "\n")
    return tmp_path


def test_grammar():
    tree = parse(DEFINITIONS["M90001"])
    assert [c.kind for c in tree.children] == ["ko", "alt", "alt"]
    complex_ = parse(DEFINITIONS["M90002"]).children[0]
    assert complex_.kind == "complex"
    assert complex_.optional == (False, False, False, True)


@pytest.mark.parametrize(("present", "complete", "score"), [
    ({"K90001", "K90002", "K90006"}, 1.0, 1.0),
    ({"K90001", "K90003", "K90004", "K90005"}, 1.0, 1.0),
    ({"K90001", "K90004"}, 1 / 3, (1 + 0 + 0.5) / 3),
    (set(), 0.0, 0.0),
])
def test_stepwise_completeness(present, complete, score):
    result = evaluate(_modules()["M90001"], present, _modules())
    assert result.completeness == pytest.approx(complete)
    assert result.score == pytest.approx(score)


def test_complex_optional_components_and_missing():
    modules = _modules()
    assert evaluate(modules["M90002"], {"K90011", "K90012", "K90013"}, modules).completeness == 1.0
    partial = evaluate(modules["M90002"], {"K90011", "K90012"}, modules)
    assert partial.score == pytest.approx(2 / 3)
    assert partial.steps[0].missing == ["K90013"]


def test_module_references_and_undefined_steps():
    modules = _modules()
    result = evaluate(modules["M90003"], {"K90021", "K90011", "K90012", "K90013"}, modules)
    assert (result.steps_complete, result.steps_total) == (2, 2)  # '--' is left out


@pytest.fixture
def store():
    store = DuckDBStore()
    store.conn.execute("""
        INSERT INTO bins (bin_id, completeness, contamination) VALUES ('b1', 62.0, 1.0), ('b2', 95.0, 0.5);
        INSERT INTO contigs (contig_id, bin_id, length) VALUES ('c1', 'b1', 50000), ('c2', 'b2', 50000);
        INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, gene_index, sequence_length)
        VALUES ('p1', 'c1', 'b1', 1, 900, '+', 0, 300), ('p2', 'c1', 'b1', 1000, 1900, '+', 1, 300),
               ('p3', 'c1', 'b1', 2000, 2900, '+', 2, 300), ('p9', 'c1', 'b1', 90000, 90900, '+', 40, 300),
               ('q1', 'c2', 'b2', 1, 900, '+', 0, 300);
        INSERT INTO annotations (annotation_id, protein_id, source, accession, name, evalue)
        VALUES (1, 'p1', 'kofam', 'K90011', 'aA', 1e-50), (2, 'p2', 'kofam', 'K90012', 'aB', 1e-50),
               (3, 'p3', 'kofam', 'K90013', 'aC', 1e-50), (4, 'p9', 'kofam', 'K90001', 'x', 1e-50),
               (5, 'q1', 'kegg', 'K90001', 'x', 1e-50), (6, 'p1', 'pfam', 'PF00001', 'y', 1e-50);
    """)
    return store


def test_genome_modules_report_proteins_and_bin_quality(store, kegg_dir):
    rows = genome_modules(store, directory=kegg_dir)
    complex_row = next(r for r in rows if r["bin_id"] == "b1" and r["module"] == "M90002")
    assert complex_row["completeness"] == 1.0
    assert complex_row["steps"][0]["proteins"] == {"K90011": ["p1"], "K90012": ["p2"], "K90013": ["p3"]}
    assert complex_row["bin_completeness"] == 62.0
    assert {r["module"] for r in rows if r["bin_id"] == "b2"} == {"M90001"}  # zero-score modules are omitted
    assert "62.0% complete" in modules_markdown(rows, "t")


def test_min_completeness_and_module_filter(store, kegg_dir):
    rows = genome_modules(store, directory=kegg_dir, modules=["M90002"], min_completeness=1.0)
    assert [(r["bin_id"], r["module"]) for r in rows] == [("b1", "M90002")]


def test_locus_modules_stay_within_the_window(store, kegg_dir):
    near = locus_modules(store, "p2", window=2, directory=kegg_dir)
    assert near[0]["module"] == "M90002"
    assert near[0]["completeness"] == 1.0
    assert "M90001" not in {r["module"] for r in near}  # its KO sits 38 genes away
    assert near[0]["steps"][0]["proteins"]["K90011"] == ["p1 (-1)"]


def test_missing_definitions_raise(store, tmp_path):
    with pytest.raises(FileNotFoundError, match="setup-kegg"):
        genome_modules(store, directory=tmp_path)
