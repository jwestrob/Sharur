"""Domain-architecture patterns and search."""

import pytest

from sharur.architecture import (
    Domain,
    PatternError,
    architecture_markdown,
    compact,
    compile_pattern,
    encode,
    resolve,
    search_architecture,
)
from sharur.operators.cards import card, card_markdown
from sharur.storage.duckdb_store import DuckDBStore


def _domains(names):
    return [Domain(n, f"PF{90000 + i}", 1 + 100 * i, 90 + 100 * i, 1e-10, "pfam") for i, n in enumerate(names)]


ARCH = ["SP", "Big_2", "Big_2", "Big_2", "Big_2", "Big_2", "Big_3", "VWA", "CHU_C"]


@pytest.mark.parametrize(("pattern", "matches"), [
    ("Big_2{5,} . VWA", True),
    ("Big_2{6,}", False),
    ("^ SP Big_2+", True),
    ("^ Big_2", False),
    ("VWA CHU_C $", True),
    ("Big_3 $", False),
    ("Big_* +", True),
    ("( Big_3 | Cadherin ) VWA", True),
    ("( Cadherin | Fn3 )", False),
    ("Big_2 ? VWA", True),
    ("PF90007", True),            # accession of the VWA domain
    ("^ . {9} $", True),
    ("^ . {10,} $", False),
])
def test_pattern_semantics(pattern, matches):
    text, _ = encode(_domains(ARCH))
    assert bool(compile_pattern(pattern).regex.search(text)) is matches


def test_prefilter_clauses_and_minimum_length():
    compiled = compile_pattern("^ SP Big_2{5,} ( VWA | VWA_2 ) Fn3 ? Big_* +")
    assert set(compiled.clauses) == {frozenset({"SP"}), frozenset({"Big_2"}), frozenset({"VWA", "VWA_2"}),
                                     frozenset({"Big_*"})}
    assert compiled.min_domains == 1 + 5 + 1 + 0 + 1
    assert compile_pattern("Big_2 {0,3} VWA").clauses == (frozenset({"VWA"}),)


@pytest.mark.parametrize("bad", ["", "( VWA", "VWA )", "A | B", "+ VWA", "VWA ^", "( | VWA )"])
def test_malformed_patterns_raise(bad):
    with pytest.raises(PatternError):
        compile_pattern(bad)


def test_overlaps_resolve_by_evalue():
    hits = [Domain("TPR_1", "PF00515", 10, 44, 1e-8, "pfam"),
            Domain("TPR_2", "PF07719", 12, 45, 1e-5, "pfam"),     # overlaps TPR_1, weaker: dropped
            Domain("TPR_1", "PF00515", 50, 84, 1e-6, "pfam"),
            Domain("CAZy", "GT2", None, None, None, "cazy")]       # no coordinates: left out
    assert [(d.name, d.start) for d in resolve(hits)] == [("TPR_1", 10), ("TPR_1", 50)]
    assert len(resolve(hits, max_overlap=1.0)) == 3
    assert compact(["A", "B", "B", "B", "A"]) == "A - B x3 - A"


@pytest.fixture
def store():
    store = DuckDBStore()
    rows, n = [], 0
    for pid, names in {"giant": ARCH, "small": ["VWA"], "other": ["Big_2", "Big_2", "VWA"]}.items():
        for i, name in enumerate(names):
            n += 1
            rows.append(f"({n}, '{pid}', 'pfam', 'PF{90000 + ARCH.index(name) if name in ARCH else 99999}', "
                        f"'{name}', 1e-10, {1 + 100 * i}, {90 + 100 * i})")
    store.conn.execute(f"""
        INSERT INTO bins (bin_id) VALUES ('b1'), ('b2');
        INSERT INTO contigs (contig_id, bin_id, length) VALUES ('c1', 'b1', 9000), ('c2', 'b2', 9000);
        INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, gene_index, sequence_length)
        VALUES ('giant', 'c1', 'b1', 1, 3000, '+', 0, 1000), ('small', 'c1', 'b1', 4000, 4600, '+', 1, 200),
               ('other', 'c2', 'b2', 1, 1200, '+', 0, 400);
        INSERT INTO annotations (annotation_id, protein_id, source, accession, name, evalue, start_aa, end_aa)
        VALUES {", ".join(rows)};
    """)
    return store


def test_search_reports_matches_with_coordinates(store):
    result = search_architecture(store, "Big_2{2,} . ? VWA")
    assert result["total"] == 2
    giant = next(r for r in result["records"] if r["protein_id"] == "giant")
    assert giant["architecture"] == "SP - Big_2 x5 - Big_3 - VWA - CHU_C"
    assert (giant["match"]["start_aa"], giant["match"]["end_aa"], giant["match"]["n_domains"]) == (101, 790, 7)
    assert giant["bin_id"] == "b1"
    assert "2 proteins match" in architecture_markdown(result)


def test_search_limit_bins_and_globs(store):
    assert search_architecture(store, "VWA", limit=1)["total"] == 3
    assert len(search_architecture(store, "VWA", limit=1)["records"]) == 1
    assert [r["protein_id"] for r in search_architecture(store, "VWA", bins=["b2"])["records"]] == ["other"]
    assert search_architecture(store, "^ Big_* {2} VWA $")["total"] == 1


def test_card_shows_architecture(store):
    assert "Pfam architecture: SP - Big_2 x5 - Big_3 - VWA - CHU_C" in card_markdown(card(store, "giant"))
