"""Compact semantic-term index: provider contract, generation integrity, and every consumer route.

The SQL provider over ``semantic_terms`` and the compact provider must return identical
results; both are also checked against a brute-force reading of the stored rows with SQL
NULL semantics for the active-membership expression.
"""

from __future__ import annotations

import json
import os
import random
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

import duckdb
import pytest


pytest.importorskip("pyroaring")

from sharur.semantic_index import (
    CompactSemanticProvider,
    GenerationError,
    SearchRequest,
    SqlSemanticProvider,
    StaleGenerationError,
    assemble,
    attach,
    resolve,
    select,
)
from sharur.semantic_index.build import (
    build_forward,
    build_genome_scope,
    build_membership,
)


UNKNOWN = "__unknown_term__"


class _Store:
    """Minimal read-only store over a fixture DuckDB."""

    def __init__(self, path):
        self.con = duckdb.connect(str(path), read_only=True)

    def execute(self, query, params=None):
        return self.con.execute(query, params or []).fetchall()

    def close(self):
        self.con.close()


def _contract_rows():
    rows = []
    for i in range(3000):
        rows.append((f"p{i:04d}", "dense", "atom", "activity", "implies", "pfam", "PF1"))
    rows += [
        ("p0001", "sparse:a", "atom", "activity", "implies", "pfam", "PF2"),
        ("p0100", "sparse:a", "atom", "activity", "implies", "kofam", "K1"),
        ("é2", "sparse:a", "atom", "activity", "supports", "", ""),
        ("é2", "sparse:a", "atom", "activity", "supports", "", ""),          # exact duplicate row
        ("é2", "sparse:a", "atom", "activity", "supports", None, None),      # NULL beside blank
        ("ζ-1", "kind_null_implies", None, "role", "implies", "pfam", "PF3"),
        ("ζ-1", "kind_null_excludes", None, "role", "excludes", "pfam", "PF3"),
        ("Z9", "atom_null_rel", "atom", None, None, "pfam", "PF4"),
        ("Z9", "composite_null_rel", "composite", None, None, "_composite", ""),
        ("a_1", "atom_excl", "atom", "activity", "excludes", "pfam", "PF5"),
        ("p0002", "atom_excl", "atom", "activity", "implies", "pfam", "PF5"),
        ("p0002", "atom_excl", "atom", "activity", "excludes", "hyddb", "X"),
        ("a_1", "τ:unicode", "direct_access", "architecture", "implies", "pfam", "PF6"),
        ("nullbin", "sparse:a", "atom", "activity", "implies", "pfam", "PF2"),
        ("emptybin", "sparse:a", "atom", "activity", "implies", "pfam", "PF2"),
        ("outside", "τ:unicode", "direct_access", "architecture", "implies", "pfam", "PF6"),
        ("p2999", "τ:unicode", "composite", "architecture", "flags", "_composite", ""),
    ]
    return rows


PROTEINS = ([(f"p{i:04d}", "g1") for i in range(3000)]
            + [("é2", "g2"), ("ζ-1", "g2"), ("Z9", "g2"), ("a_1", "g2"), ("zero_hit", "g2"),
               ("nullbin", None), ("emptybin", ""), ("outside", "g_outside")])
BINS = ["g1", "g2", "", "g_empty"]


def _make_contract_db(path):
    con = duckdb.connect(str(path))
    con.execute("CREATE TABLE bins(bin_id VARCHAR)")
    con.execute("CREATE TABLE proteins(protein_id VARCHAR, bin_id VARCHAR)")
    con.execute("CREATE TABLE semantic_terms(protein_id VARCHAR, term_id VARCHAR, term_kind VARCHAR, "
                "facet VARCHAR, relation VARCHAR, source_db VARCHAR, source_accession VARCHAR)")
    con.executemany("INSERT INTO bins VALUES (?)", [(b,) for b in BINS])
    con.executemany("INSERT INTO proteins VALUES (?, ?)", PROTEINS)
    rows = _contract_rows()
    random.Random(3).shuffle(rows)           # physical order must not matter
    con.executemany("INSERT INTO semantic_terms VALUES (?, ?, ?, ?, ?, ?, ?)", rows)
    con.close()
    return path


def _active(kind, relation):
    # SQL: (term_kind != 'atom' OR relation != 'excludes'); NULL operands are unknown, not true
    left = None if kind is None else kind != "atom"
    right = None if relation is None else relation != "excludes"
    return left is True or right is True


def _membership():
    out = defaultdict(set)
    for pid, term, kind, _, relation, _, _ in _contract_rows():
        if _active(kind, relation):
            out[term].add(pid)
    return out


def _build_generation(tmp, db, index, *, seal=None, name="build"):
    work = tmp / name
    build_membership(db, work / "membership", seal=seal)
    build_forward(db, work / "membership", work / "forward")
    build_genome_scope(db, work / "membership", work / "genome_scope")
    gen = assemble(index, db, membership=work / "membership", forward=work / "forward",
                   genome_scope=work / "genome_scope", seal_path=seal)
    select(index, gen.generation_id)
    return gen


@pytest.fixture(scope="module")
def contract(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("contract")
    db = _make_contract_db(tmp / "contract.duckdb")
    index = tmp / "index"
    _build_generation(tmp, db, index)
    store = _Store(db)
    compact = CompactSemanticProvider.open(index, db, store=store)
    yield {"db": db, "index": index, "tmp": tmp, "sql": SqlSemanticProvider(store), "compact": compact}
    compact.close()
    store.close()


def _brute(request: SearchRequest):
    members = _membership()
    owner = dict(PROTEINS)
    if request.empty:
        return []
    universe = set().union(*members.values())
    result = set(universe)
    for term in request.has:
        result &= members.get(term, set())
    if request.any_of:
        result &= set().union(*(members.get(t, set()) for t in request.any_of))
    for term in request.lacks:
        result -= members.get(term, set())
    if request.genomes is not None:
        result = {p for p in result if owner.get(p) is not None and owner[p] in request.genomes}
    return sorted(result, key=lambda s: s.encode())


def test_both_codecs_and_counts_are_exercised(contract):
    manifest = contract["compact"].membership.manifest
    assert manifest["list_terms"] > 0 and manifest["roaring_terms"] > 0
    scope = contract["compact"].scope.manifest
    assert scope["roaring_bins"] > 0 and scope["list_bins"] > 0 and scope["empty_bins"] == 1
    assert scope["null_ownership_count"] == 1 and scope["owners_outside_declared_bin_catalog"] == 1
    # NULL-kind/excludes, atom/NULL-relation and atom/excludes-only terms are outside membership
    catalog = dict(contract["compact"].term_catalog())
    assert catalog == dict(contract["sql"].term_catalog())
    assert {"kind_null_excludes", "atom_null_rel"}.isdisjoint(catalog)
    assert catalog["kind_null_implies"] == 1 and catalog["composite_null_rel"] == 1
    assert catalog["atom_excl"] == 1 and catalog["sparse:a"] == 5


@pytest.mark.parametrize("has,lacks,expected", [
    ([], [], []),
    (["dense"], [], "dense"),
    (["sparse:a", "dense"], [], ["p0001", "p0100"]),
    (["sparse:a", "sparse:a"], [], []),                     # repeated positive: unsatisfiable (live quirk)
    ([UNKNOWN], [], []),
    (["sparse:a", UNKNOWN], [], []),
    (["sparse:a"], [UNKNOWN], "sparse:a"),
    (["sparse:a"], ["sparse:a"], []),
    (["sparse:a"], ["dense", "dense"], ["emptybin", "nullbin", "é2"]),
    ([], ["dense"], "universe-dense"),
    ([], [UNKNOWN], "universe"),
    (["atom_excl"], [], ["p0002"]),
    (["kind_null_excludes"], [], []),
])
def test_search_compatible_matches_the_live_contract(contract, has, lacks, expected):
    members = _membership()
    universe = set().union(*members.values())
    if expected == "universe":
        expected = universe
    elif expected == "universe-dense":
        expected = universe - members["dense"]
    elif isinstance(expected, str):
        expected = members[expected]
    for name in ("sql", "compact"):
        full = contract[name].search_compatible(has, lacks, 100_000)
        assert len(full) == len(set(full)) and set(full) == set(expected), name
        limited = contract[name].search_compatible(has, lacks, 7)
        assert len(limited) == min(7, len(expected)) and len(set(limited)) == len(limited)
        assert set(limited) <= set(expected)


def test_operator_routes_through_an_attached_provider(contract):
    from sharur.operators.predicates_v2 import (  # noqa: PLC0415
        _fetch_semantic_terms,
        search_by_atoms,
    )

    sql_store = _Store(contract["db"])
    compact_store = _Store(contract["db"])
    attach(compact_store, contract["compact"])
    try:
        for has, lacks in ((["sparse:a"], []), (["sparse:a", "dense"], ["τ:unicode"]), ([], ["dense"]),
                           (["sparse:a", "sparse:a"], []), ([], [])):
            assert (set(search_by_atoms(sql_store, has, lacks, limit=10_000))
                    == set(search_by_atoms(compact_store, has, lacks, limit=10_000)))
        before = contract["compact"].counters.snapshot().get("search_compatible", 0)
        search_by_atoms(compact_store, ["dense"], [], limit=3)
        assert contract["compact"].counters.snapshot()["search_compatible"] == before + 1
        for pid in ("é2", "ζ-1", "Z9", "a_1", "p0002", "zero_hit", "missing"):
            assert (_fetch_semantic_terms(sql_store, pid, {"semantic_terms"})
                    == _fetch_semantic_terms(compact_store, pid, set()))
    finally:
        sql_store.close()
        compact_store.close()


def test_search_extension_matches_brute_force(contract):
    terms = sorted(_membership()) + [UNKNOWN]
    genomes = BINS + ["g_outside", "unknown_bin"]
    rng = random.Random(11)
    for _ in range(160):
        request = SearchRequest.build(
            has=rng.sample(terms, rng.randint(0, 2)), any_of=rng.sample(terms, rng.randint(0, 2)),
            lacks=rng.sample(terms, rng.randint(0, 2)),
            genomes=None if rng.random() < 0.4 else rng.choices(genomes, k=rng.randint(0, 3)))
        expected = _brute(request)
        offset = rng.choice([0, 0, 1, 5, 2999])
        for name in ("sql", "compact"):
            page = contract[name].search(request, offset=offset, limit=50)
            assert page.total == len(expected), (name, request)
            assert page.protein_ids == expected[offset:offset + 50], (name, request)


def test_genome_scope_semantics(contract):
    for name in ("sql", "compact"):
        p = contract[name]
        unscoped = p.search(SearchRequest.build(has=["sparse:a"])).protein_ids
        assert "nullbin" in unscoped and "emptybin" in unscoped
        assert p.search(SearchRequest.build(has=["sparse:a"], genomes=[])).total == 0
        assert p.search(SearchRequest.build(has=["sparse:a"], genomes=["unknown_bin"])).total == 0
        assert p.search(SearchRequest.build(has=["sparse:a"], genomes=[""])).protein_ids == ["emptybin"]
        assert (p.search(SearchRequest.build(has=["dense"], genomes=["g1", "g1"])).total
                == p.search(SearchRequest.build(has=["dense"], genomes=["g1"])).total == 3000)
        assert p.search(SearchRequest.build(has=["τ:unicode"], genomes=["g_outside"])).protein_ids == ["outside"]
        assert p.search(SearchRequest.build(has=["dense"], genomes=["g_empty"])).total == 0
        assert p.search(SearchRequest.build()).total == 0
    with pytest.raises(ValueError):
        SearchRequest.build(has=["dense"], genomes=["g1", None])
    with pytest.raises(ValueError):
        SearchRequest.build(has=["dense"], genomes="g1")


def test_rich_rows_are_exact_and_identically_ordered(contract):
    expected = defaultdict(Counter)
    for row in _contract_rows():
        expected[row[0]][row[1:]] += 1
    for pid, _ in PROTEINS + [("missing", None)]:
        sql_rows = contract["sql"].rich_rows(pid)
        compact_rows = contract["compact"].rich_rows(pid)
        assert compact_rows == sql_rows, pid
        assert Counter(compact_rows) == expected.get(pid, Counter())
        assert contract["compact"].rich_row_count(pid) == contract["sql"].rich_row_count(pid) == len(sql_rows)
    rows = contract["compact"].rich_rows("é2")
    assert rows.count(("sparse:a", "atom", "activity", "supports", "", "")) == 2
    assert ("sparse:a", "atom", "activity", "supports", None, None) in rows
    assert contract["compact"].rich_rows_many(["Z9", "missing"]) == contract["sql"].rich_rows_many(["Z9", "missing"])


def test_term_counts_and_suggestions_agree(contract):
    terms = ["dense", "sparse:a", UNKNOWN, "τ:unicode"]
    assert contract["sql"].term_counts(terms) == contract["compact"].term_counts(terms)
    assert contract["sql"].suggest_terms("a", 10) == contract["compact"].suggest_terms("a", 10)
    assert contract["compact"].suggest_terms("sparse", 3)[0] == ("sparse:a", 5)


def _copy_index(src: Path, dst: Path) -> Path:
    shutil.copytree(src, dst)
    return dst


def test_generation_rejects_corrupt_and_stale_inputs(contract, tmp_path):
    db, index = contract["db"], contract["index"]
    gen = resolve(index)
    store = _Store(db)
    try:
        corrupt = _copy_index(index, tmp_path / "corrupt")
        payload = resolve(corrupt).component("membership") / "postings.bin"
        data = bytearray(payload.read_bytes())
        data[len(data) // 2] ^= 0xFF
        payload.write_bytes(bytes(data))
        with pytest.raises(GenerationError, match="Payload differs"):
            CompactSemanticProvider.open(corrupt, db)

        edited = _copy_index(index, tmp_path / "edited")
        manifest = resolve(edited).component("forward") / "manifest.json"
        manifest.write_text(manifest.read_text().replace('"rows"', '"rows" ', 1))
        with pytest.raises(GenerationError, match="manifest differs"):
            CompactSemanticProvider.open(edited, db)

        dangling = _copy_index(index, tmp_path / "dangling")
        (dangling / "CURRENT").write_text(json.dumps({"generation": "g1-missing"}))
        with pytest.raises(GenerationError, match="missing generation"):
            CompactSemanticProvider.open(dangling, db)
        with pytest.raises(GenerationError, match="No such generation"):
            select(dangling, "g1-missing")

        copy = tmp_path / "copy.duckdb"
        shutil.copyfile(db, copy)                       # same bytes, new inode and mtime
        with pytest.raises(StaleGenerationError):
            CompactSemanticProvider.open(index, copy)
        CompactSemanticProvider.open(index, copy, content_match=True).close()

        changed = tmp_path / "changed.duckdb"
        shutil.copyfile(db, changed)
        con = duckdb.connect(str(changed))
        con.execute("INSERT INTO proteins VALUES ('late', 'g2')")
        con.close()
        with pytest.raises(StaleGenerationError):
            CompactSemanticProvider.open(index, changed, content_match=True)

        wrong_store = _Store(changed)
        try:
            with pytest.raises(GenerationError, match="proteins"):
                CompactSemanticProvider.open(index, db, store=wrong_store)
        finally:
            wrong_store.close()
        assert resolve(index).generation_id == gen.generation_id
    finally:
        store.close()


def test_seal_binding(tmp_path):
    db = _make_contract_db(tmp_path / "sealed.duckdb")
    seal = tmp_path / "dataset.seal.json"
    seal.write_text(json.dumps({"dataset_id": "fixture"}))
    _build_generation(tmp_path, db, tmp_path / "index", seal=seal)
    CompactSemanticProvider.open(tmp_path / "index", db).close()
    seal.write_text(json.dumps({"dataset_id": "other"}))
    with pytest.raises(StaleGenerationError, match="seal"):
        CompactSemanticProvider.open(tmp_path / "index", db)


def test_generation_switch_is_atomic_and_open_providers_stay_pinned(tmp_path):
    db = _make_contract_db(tmp_path / "own.duckdb")      # no other connection holds this file
    index = tmp_path / "index"
    first = _build_generation(tmp_path, db, index).generation_id
    pinned = CompactSemanticProvider.open(index, db)
    second = _build_generation(tmp_path, db, index, name="rebuild")
    assert second.generation_id != first and resolve(index).generation_id == second.generation_id
    assert pinned.identity()["generation_id"] == first
    reopened = CompactSemanticProvider.open(index, db)
    request = SearchRequest.build(has=["sparse:a"], lacks=["dense"])
    assert pinned.search(request).protein_ids == reopened.search(request).protein_ids
    # a rebuild of the same source reproduces every payload byte
    for name in ("membership", "forward", "genome_scope"):
        a, b = resolve(index).component(name), pinned.generation.component(name)
        assert ({f: e["sha256"] for f, e in json.loads((a / "manifest.json").read_text())["files"].items()}
                == {f: e["sha256"] for f, e in json.loads((b / "manifest.json").read_text())["files"].items()})
    pinned.close()
    reopened.close()


def test_close_releases_readers_and_refuses_use(contract):
    provider = CompactSemanticProvider.open(contract["index"], contract["db"])
    assert provider.search_compatible(["dense"], [], 1)
    provider.close()
    provider.close()
    with pytest.raises(GenerationError, match="closed"):
        provider.search_compatible(["dense"], [], 1)
    with pytest.raises(GenerationError, match="closed"):
        provider.rich_rows("é2")


def test_package_imports_and_sql_paths_work_without_pyroaring():
    code = """
import sys
sys.modules['pyroaring'] = None
import sharur.semantic_index, sharur.operators.predicates_v2, sharur.browser.app, sharur.query.server
from sharur.semantic_index import CompactSemanticProvider
try:
    CompactSemanticProvider.open('/nonexistent', '/nonexistent')
except ModuleNotFoundError as exc:
    print('needs', exc.name)
"""
    out = subprocess.run([sys.executable, "-c", code], check=False, capture_output=True, text=True, cwd=Path(__file__).parents[1])
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip() == "needs pyroaring"


# --------------------------------------------------------------------------- #
# Browser, query service and readiness over the real schema
# --------------------------------------------------------------------------- #

SEQ = "MKTAYIAKQRQISFVKSHFSRQ" * 10


@pytest.fixture(scope="module")
def dataset(tmp_path_factory):
    from sharur.predicates_v2.persistence import generate_and_persist_v2  # noqa: PLC0415
    from sharur.storage.duckdb_store import DuckDBStore  # noqa: PLC0415

    tmp = tmp_path_factory.mktemp("dataset")
    path = tmp / "sharur.duckdb"
    store = DuckDBStore(str(path))
    store.conn.execute(f"""
        INSERT INTO bins (bin_id, completeness, contamination, taxonomy) VALUES
            ('bin|1', 91.5, 1.2, 'd__Archaea;p__Alpha'), ('bin|2', 80.0, 2.0, 'd__Archaea;p__Beta');
        INSERT INTO contigs (contig_id, bin_id, length, length_source) VALUES
            ('bin|1_c1', 'bin|1', 9000, 'assembly'), ('bin|2_c1', 'bin|2', 9000, 'assembly');
        INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, gene_index,
                              sequence, sequence_length, partial)
        VALUES ('bin|1_c1_1', 'bin|1_c1', 'bin|1', 1, 600, '+', 0, '{SEQ}', 220, '00'),
               ('bin|1_c1_2', 'bin|1_c1', 'bin|1', 700, 1400, '-', 1, '{SEQ}', 220, '00'),
               ('bin|2_c1_1', 'bin|2_c1', 'bin|2', 1, 600, '+', 0, '{SEQ}', 220, '00'),
               ('bin|2_c1_2', 'bin|2_c1', 'bin|2', 700, 1400, '+', 1, '{SEQ}', 220, '00');
        INSERT INTO annotations (annotation_id, protein_id, source, accession, name, description, evalue, score,
                                 start_aa, end_aa)
        VALUES (1, 'bin|1_c1_2', 'pfam', 'PF00005', 'ABC_tran', 'ABC transporter', 1e-40, 150.0, 20, 160),
               (2, 'bin|2_c1_2', 'pfam', 'PF00005', 'ABC_tran', 'ABC transporter', 1e-30, 120.0, 20, 160);
    """)
    generate_and_persist_v2(store, chunk_size=1, return_states=False, update_legacy_predicates=True)
    # rows the generator does not emit: NULL facet/relation, duplicate, inactive atom
    store.conn.execute("""
        INSERT INTO semantic_terms VALUES
            ('bin|2_c1_1', 'pfam:PF00005', 'atom', NULL, NULL, '', ''),
            ('bin|2_c1_2', 'fixture_dup', 'composite', 'role', 'implies', '_composite', ''),
            ('bin|2_c1_2', 'fixture_dup', 'composite', 'role', 'implies', '_composite', ''),
            ('bin|1_c1_1', 'fixture_excluded', 'atom', 'activity', 'excludes', 'pfam', 'PF00005')
    """)
    store.close()
    index = tmp / "index"
    _build_generation(tmp, path, index)
    return {"db": path, "index": index, "tmp": tmp}


@pytest.fixture(scope="module")
def clients(dataset):
    from fastapi.testclient import TestClient  # noqa: PLC0415

    from sharur.browser import create_app  # noqa: PLC0415

    sql_app = create_app(dataset["db"], background=False)
    compact_app = create_app(dataset["db"], background=False, semantic_index=dataset["index"])
    with TestClient(sql_app) as sql, TestClient(compact_app) as compact:
        yield {"sql": sql, "compact": compact, "compact_app": compact_app}


def _strip(body):
    return {k: v for k, v in body.items() if k not in ("backend", "timing_ms")}


TERM_QUERIES = [
    {"has": "pfam:PF00005"},
    {"has": "pfam:PF00005", "scope": "bin|2"},
    {"has": "pfam:PF00005", "scope": "Alpha"},
    {"any": "pfam:PF00005 fixture_dup"},
    {"lacks": "pfam:PF00005"},
    {"has": "pfam:PF00005", "lacks": "fixture_dup"},
    {"has": "pfam:PF00005", "limit": 1, "offset": 1},
    {"has": "fixture_excluded"},
    {"has": UNKNOWN},
    {"has": "pfam:PF00005", "genomes": ""},
    {"has": "pfam:PF00005", "genomes": "bin|1,bin|1"},
]


def test_browser_term_search_is_identical_across_backends(clients):
    for params in TERM_QUERIES:
        a = clients["sql"].get("/api/v1/terms", params=params)
        b = clients["compact"].get("/api/v1/terms", params=params)
        assert a.status_code == b.status_code == 200, params
        assert a.headers["X-Sharur-Semantic-Backend"] == "sql"
        assert b.headers["X-Sharur-Semantic-Backend"].startswith("compact g1-")
        assert "search;dur=" in b.headers["Server-Timing"]
        assert _strip(a.json()) == _strip(b.json()), params
    body = clients["compact"].get("/api/v1/terms", params={"has": "pfam:PF00005", "scope": "bin|2"}).json()
    # bin|2_c1_1 holds only an atom row with a NULL relation: outside membership under SQL NULL semantics
    assert [r["protein_id"] for r in body["rows"]] == ["bin|2_c1_2"]
    assert body["request"]["scope"] == {"kind": "genome", "label": "bin|2", "genomes": 1}
    explicit_empty = clients["compact"].get("/api/v1/terms", params={"has": "pfam:PF00005", "genomes": ""}).json()
    assert explicit_empty["total"] == 0 and explicit_empty["request"]["genomes"] == []
    assert clients["compact"].get("/api/v1/terms", params={"has": "x", "scope": "nowhere"}).status_code == 400
    assert clients["compact"].get("/api/v1/terms", params={"has": "x", "scope": "bin|2", "genomes": "bin|2"}
                                  ).status_code == 400


def test_browser_stored_rows_and_pages_are_identical_across_backends(clients):
    for pid in ("bin|1_c1_1", "bin|1_c1_2", "bin|2_c1_1", "bin|2_c1_2"):
        a = clients["sql"].get(f"/api/v1/protein/{pid}/terms").json()
        b = clients["compact"].get(f"/api/v1/protein/{pid}/terms").json()
        assert _strip(a) == _strip(b) and a["count"] > 0, pid
        page = clients["compact"].get(f"/protein/{pid}/terms")
        assert page.status_code == 200 and page.headers["X-Sharur-Semantic-Backend"].startswith("compact")
    dup = clients["compact"].get("/api/v1/protein/bin|2_c1_2/terms").json()["rows"]
    assert dup.count(["fixture_dup", "composite", "role", "implies", "_composite", ""]) == 2
    null_row = clients["compact"].get("/api/v1/protein/bin|2_c1_1/terms").json()["rows"]
    assert ["pfam:PF00005", "atom", None, None, "", ""] in null_row
    excluded = clients["compact"].get("/protein/bin|1_c1_1/terms").text
    assert "outside search membership" in excluded and "fixture_excluded" in excluded
    assert clients["compact"].get("/api/v1/protein/missing/terms").status_code == 404
    assert clients["compact"].get("/protein/missing/terms").status_code == 404
    assert clients["sql"].get("/api/v1/terms/suggest", params={"q": "pfam"}).json() == \
        clients["compact"].get("/api/v1/terms/suggest", params={"q": "pfam"}).json()


def test_browser_pages_render_and_name_the_backend(clients):
    page = clients["compact"].get("/terms", params={"has": "pfam:PF00005", "scope": "Beta"})
    assert page.status_code == 200 and "bin|2_c1_2" in page.text and "compact" in page.text
    assert "within Beta" in page.text and "ABC_tran" in page.text
    assert "No genome or clade matches" in clients["compact"].get("/terms", params={"has": "a", "scope": "zz"}).text
    assert clients["sql"].get("/terms").status_code == 200
    protein = clients["compact"].get("/protein/bin%7C2_c1_2")
    assert protein.status_code == 200 and "/protein/bin%7C2_c1_2/terms" in protein.text
    assert "Stored V2 term rows" in clients["sql"].get("/protein/bin%7C2_c1_2").text
    # the main search box lists V2 terms with the compact catalog; an exact term opens its search
    found = clients["compact"].get("/search", params={"q": "fixture"})
    assert "V2 semantic terms" in found.text and "fixture_dup" in found.text
    assert "V2 semantic terms" not in clients["sql"].get("/search", params={"q": "fixture"}).text
    exact = clients["compact"].get("/search", params={"q": "fixture_dup"}, follow_redirects=False)
    assert exact.status_code == 303 and exact.headers["location"] == "/terms?has=fixture_dup"
    scoped = clients["compact"].get("/search", params={"q": "pfam:PF00005 in Beta"})
    assert "is also a V2 term" in scoped.text
    status = clients["compact"].get("/api/v1/semantic").json()
    assert status["backend"]["backend"] == "compact" and status["calls"]["search"] > 0
    assert status["counts"]["proteins"] == 4
    assert clients["sql"].get("/api/v1/semantic").json()["backend"] == {"backend": "sql", "table": "semantic_terms"}


def test_browser_refuses_a_stale_index(dataset, tmp_path):
    from sharur.browser import create_app  # noqa: PLC0415

    copy = tmp_path / "sharur.duckdb"
    shutil.copyfile(dataset["db"], copy)
    with pytest.raises(StaleGenerationError):
        create_app(copy, background=False, semantic_index=dataset["index"])


def test_browser_shutdown_closes_the_provider(dataset):
    from fastapi.testclient import TestClient  # noqa: PLC0415

    from sharur.browser import create_app  # noqa: PLC0415

    app = create_app(dataset["db"], background=False, semantic_index=dataset["index"])
    with TestClient(app) as client:
        assert client.get("/api/v1/terms", params={"has": "pfam:PF00005"}).json()["total"] == 2
    with pytest.raises(GenerationError, match="closed"):
        app.state.semantic.rich_rows("bin|1_c1_1")


def test_startup_summaries_never_read_v2_terms(dataset, tmp_path):
    from sharur.browser import startup_cache  # noqa: PLC0415

    out = tmp_path / "cache.pkl"
    startup_cache.build(dataset["db"], out)
    import pickle  # noqa: PLC0415

    rows = pickle.loads(out.read_bytes())["rows"]
    assert rows and not any("semantic_terms" in query for query, _ in rows)


def test_query_service_atom_search_uses_the_compact_index(dataset, tmp_path):
    from fastapi.testclient import TestClient  # noqa: PLC0415

    from sharur.query.server import create_app  # noqa: PLC0415

    replica = tmp_path / "replica" / "sharur.duckdb"      # a staged byte copy: new inode, bound by content
    replica.parent.mkdir()
    shutil.copyfile(dataset["db"], replica)

    def app(**kw):
        return create_app(db_path=replica, temp_directory=tmp_path / "spill", threads=2,
                          memory_limit="256MB", max_temp_directory_size="512MB", allow_unsealed=True, **kw)

    requests = [{"has": ["pfam:PF00005"], "limit": 50}, {"has": ["pfam:PF00005"], "lacks": ["fixture_dup"]},
                {"lacks": ["pfam:PF00005"], "limit": 50}, {"has": ["pfam:PF00005", "pfam:PF00005"]}]
    results = {}
    for name, kw in (("sql", {}), ("compact", {"semantic_index": dataset["index"]})):
        with TestClient(app(**kw)) as client:
            health = client.get("/health").json()
            assert health["semantic_backend"]["backend"] == name
            results[name] = [sorted(client.post("/v1/atoms/proteins", json=r).json()["raw"]) for r in requests]
    assert results["sql"] == results["compact"]
    assert results["compact"][0] == ["bin|1_c1_2", "bin|2_c1_2"] and results["compact"][3] == []


def test_preflight_reports_the_semantic_index(dataset, tmp_path):
    from sharur.capabilities import build_capability_brief  # noqa: PLC0415

    brief = {c.capability_id: c for c in build_capability_brief(dataset["db"], include_execution=False,
                                                                semantic_index=dataset["index"]).capabilities}
    assert brief["semantic_index"].state.value == "available"
    assert brief["semantic_v2"].evidence["terms_backend"] == "compact"
    copy = tmp_path / "sharur.duckdb"
    shutil.copyfile(dataset["db"], copy)
    stale = {c.capability_id: c for c in build_capability_brief(copy, include_execution=False,
                                                                semantic_index=dataset["index"]).capabilities}
    assert stale["semantic_index"].state.value == "stale"
    assert "terms_backend" not in stale["semantic_v2"].evidence     # default briefs are unchanged


def test_browse_cli_refuses_a_stale_index(dataset, tmp_path):
    copy = tmp_path / "sharur.duckdb"
    shutil.copyfile(dataset["db"], copy)
    out = subprocess.run([sys.executable, "-m", "sharur.cli", "browse", "--db", str(copy), "--no-summary-cache",
                          "--semantic-index", str(dataset["index"]), "--port", "0"],
                         check=False, capture_output=True, text=True, cwd=Path(__file__).parents[1],
                         env={**os.environ, "PYTHONWARNINGS": "ignore"}, timeout=120)
    assert out.returncode == 1 and "Semantic index refused" in out.stderr


def test_facade_opts_in_to_the_compact_index(dataset, clients):
    from sharur.operators import Sharur  # noqa: PLC0415

    sql = Sharur(dataset["db"], read_only=True)
    compact = Sharur(dataset["db"], read_only=True, semantic_index_path=dataset["index"])
    assert compact.store.semantic_provider.backend == "compact"
    assert getattr(sql.store, "semantic_provider", None) is None
    for has, lacks in ((["pfam:PF00005"], []), ([], ["pfam:PF00005"]), (["pfam:PF00005", "fixture_dup"], [])):
        assert (sorted(sql.search_by_atoms(has=has, lacks=lacks, limit=100))
                == sorted(compact.search_by_atoms(has=has, lacks=lacks, limit=100)))
    for pid in ("bin|1_c1_1", "bin|2_c1_2"):
        assert sql.explain(pid)["terms"] == compact.explain(pid)["terms"]


def test_builder_connection_closes_before_its_spill_directory_is_removed(tmp_path):
    """DuckDB removes its own spill files on close; the directory must outlive the connection."""
    from sharur.semantic_index.build import _session  # noqa: PLC0415

    db = _make_contract_db(tmp_path / "spill.duckdb")
    with _session(db, 1, "256MB", tmp_path, "spill-") as con:
        spill = Path(con.execute("SELECT current_setting('temp_directory')").fetchone()[0])
        assert spill.is_dir() and spill.parent == tmp_path
    assert not spill.exists()
    with pytest.raises(duckdb.ConnectionException):
        con.execute("SELECT 1")
