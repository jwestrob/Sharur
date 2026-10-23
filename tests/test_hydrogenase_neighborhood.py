"""The hydrogenase neighborhood rule: support needs a marker that is not a Complex I homolog."""

import pytest

from sharur.hydrogenase import neighborhood as nb
from sharur.hydrogenase.ko_association import named_kos
from sharur.hydrogenase.neighborhood import neighborhood_contexts
from sharur.predicates.mappings import kegg_map
from sharur.predicates.mappings.kegg_evidence import names_hydrogenase
from sharur.storage.duckdb_store import DuckDBStore

# Each case: a Group 4 call at gene 1 of its own contig, neighbors at genes 2, 3, ... carrying (source, accession, name).
CASES = {
    "hyf_only": [("kegg", "K12137", "K12137"), ("kegg", "K12141", "K12141"), ("kegg", "K12144", "K12144")],
    "hyp": [("kegg", "K12137", "K12137"), ("kegg", "K04651", "K04651")],
    "fpo": [("kegg", "K22167", "K22167"), ("kegg", "K22168", "K22168")],
    "catalytic": [("pfam", "PF00374.24", "NiFeSe_Hases")],
    "hyf_nuo": [("kegg", "K12137", "K12137"), ("kegg", "K00341", "K00341")],
    "hyp_nuo": [("kegg", "K04652", "K04652"), ("kegg", "K00333", "K00333")],
    "eha_pfam": [("pfam", "EhaB", "EhaB")],
    "hyci_pfam": [("pfam", "HycI", "HycI")],
    # enoyl-CoA hydratase, CO dehydrogenase accessory protein, Cas4, a sigma-54 activator, a HydDB-derived label,
    # a Mrp antiporter subunit family and a [NiFe] small-subunit family shared with NuoB
    "unrelated": [("pfam", "ECH_1", "ECH_1"), ("pfam", "CooT", "CooT"), ("kegg", "K07464", "K07464"),
                  ("kegg", "K15836", "K15836"), ("hyddb_subgroup", "ech_hydrogenase", "ech_hydrogenase"),
                  ("pfam", "MbhD", "MbhD"), ("pfam", "Oxidored_q6", "Oxidored_q6")],
}
EXPECTED = {
    "hyf_only": "complex_i_homolog_only",
    "hyp": "supported",
    "fpo": "complex_i_context",
    "catalytic": "supported",
    "hyf_nuo": "complex_i_context",
    "hyp_nuo": "mixed",
    "eha_pfam": "complex_i_homolog_only",
    "hyci_pfam": "supported",
    "unrelated": "no_markers",
    "loose": "no_position",
}


@pytest.fixture
def store(tmp_path):
    path = tmp_path / "sharur.duckdb"
    s = DuckDBStore(str(path))
    c = s.conn
    c.execute("INSERT INTO bins (bin_id) VALUES ('g1')")
    n = 0
    for case, neighbors in CASES.items():
        c.execute("INSERT INTO contigs (contig_id, bin_id, length) VALUES (?, 'g1', 20000)", [case])
        for i in range(len(neighbors) + 2):
            c.execute("INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, gene_index, "
                      "sequence_length) VALUES (?, ?, 'g1', ?, ?, '+', ?, 300)",
                      [f"{case}_{i}", case, 1000 * i + 1, 1000 * i + 900, i])
        for i, (source, acc, name) in enumerate(neighbors, start=2):
            n += 1
            c.execute("INSERT INTO annotations (annotation_id, protein_id, source, accession, name) VALUES (?, ?, ?, ?, ?)",
                      [n, f"{case}_{i}", source, acc, name])
    c.execute("INSERT INTO contigs (contig_id, bin_id, length) VALUES ('loose', 'g1', 900)")
    c.execute("INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand) "
              "VALUES ('loose', 'loose', 'g1', 0, 900, '+')")
    s.close()
    ro = DuckDBStore(str(path), read_only=True)
    yield ro
    ro.close()


@pytest.fixture
def contexts(store):
    return neighborhood_contexts(store, [f"{case}_1" for case in CASES] + ["loose"])


def test_verdicts(contexts):
    assert {pid.removesuffix("_1"): ctx.verdict for pid, ctx in contexts.items()} == EXPECTED


def test_hyf_only_neighbors_are_complex_i_like(contexts):
    ctx = contexts["hyf_only_1"]
    assert not ctx.hydrogenase and not ctx.complex_i
    assert ctx.complex_i_like == {"K12137 hyfB", "K12141 hyfF", "K12144 hyfI"}
    assert ctx.complex_i_like_genes == {"hyf_only_2", "hyf_only_3", "hyf_only_4"}


def test_hyp_neighbor_supports_beside_complex_i_like_subunits(contexts):
    ctx = contexts["hyp_1"]
    assert ctx.hydrogenase == {"K04651 hypA"} and ctx.hydrogenase_genes == {"hyp_3"}
    assert ctx.complex_i_like == {"K12137 hyfB"}


def test_fpo_neighbors_are_complex_i(contexts):
    assert contexts["fpo_1"].complex_i == {"K22167 fpoL", "K22168 fpoM"}


def test_catalytic_domain_neighbor_supports(contexts):
    assert contexts["catalytic_1"].hydrogenase == {"PF00374.24 NiFeSe_Hases"}


def test_marker_tables_are_disjoint():
    tables = [set(nb.COMPLEX_I_KOS), set(nb.SPECIFIC_KOS), set(nb.COMPLEX_I_LIKE_KOS), set(nb.UNUSED_KOS)]
    assert sum(len(t) for t in tables) == len(set().union(*tables))


def test_every_kegg_hydrogenase_ko_is_classified():
    """A KO added to kegg_hyddb_snapshot.tsv must be placed on a side (or set aside) before it counts."""
    hydrogenase_side = set(nb.SPECIFIC_KOS) | set(nb.COMPLEX_I_LIKE_KOS) | set(nb.UNUSED_KOS)
    assert named_kos() == hydrogenase_side
    assert not named_kos() & set(nb.COMPLEX_I_KOS)


@pytest.mark.parametrize(("name", "named"), [
    ("cas4; CRISPR-associated exonuclease Cas4 [EC:3.1.12.1]", False),
    ("hyaB; hydrogenase large subunit [EC:1.12.99.6]", True),
    ("mvhA; F420-non-reducing hydrogenase large subunit [EC:1.12.99.- 1.8.98.5]", True),
    ("fooX; bifunctional enzyme [EC:2.7.1.1 1.12.7.2]", True),
    ("fooY; formate dehydrogenase [EC:1.17.1.9]", False),
])
def test_hydrogenase_ec_pattern_reads_whole_ec_numbers(name, named):
    assert names_hydrogenase(name) is named


def _local_snapshot():
    local = kegg_map.kegg_dir()
    path = local / "kegg_evidence_snapshot.tsv" if local else None
    if path is None or not path.exists():
        return None
    return {f[0]: f for f in (line.rstrip("\n").split("\t") for line in path.read_text().splitlines()
                              if line and not line.startswith("#"))}


@pytest.mark.skipif(_local_snapshot() is None, reason="no local KEGG build; run `sharur setup-kegg`")
def test_marker_kos_match_the_local_kegg_snapshot():
    snapshot = _local_snapshot()
    mismatched = []
    for table, hydrogenase_named in ((nb.COMPLEX_I_KOS, False), (nb.SPECIFIC_KOS, True), (nb.COMPLEX_I_LIKE_KOS, True)):
        for ko, symbol in table.items():
            row = snapshot.get(ko)
            if row is None or symbol not in [s.strip() for s in row[1].split(",")]:
                mismatched.append((ko, symbol, row[1] if row else "absent"))
            elif names_hydrogenase(row[2]) is not hydrogenase_named or "antiporter" in row[2]:
                mismatched.append((ko, symbol, "side"))
    assert not mismatched, mismatched
