"""Nitrogenase structural-locus calls and the composite that reads them."""

import json

import pytest

from sharur.predicates_v2.composites import (
    CompositeDefinition,
    clear_composites_cache,
    evaluate_composites,
    load_composites,
)
from sharur.predicates_v2.model import ClaimRelation, SemanticAtom, SemanticFacet
from sharur.predicates_v2.nitrogenase_loci import (
    CANDIDATE_TYPE,
    CONFIRMED_TYPE,
    EDGE_ATOM,
    SOURCE_DB,
    SYSTEM_ATOM,
    Contig,
    Gene,
    call_loci,
    locus_atoms,
    refresh_nitrogenase_loci,
)
from sharur.predicates_v2.persistence import (
    generate_and_persist_v2,
    refresh_composite_predicates_from_v2,
)
from sharur.storage.duckdb_store import DuckDBStore


KO = {"H": "K02588", "vH": "K22899", "D": "K02586", "K": "K02591", "vD": "K22896", "vK": "K22897",
      "vG": "K22898", "N": "K02592", "E": "K02587", "I1": "K02589", "I2": "K02590", "B": "K02585",
      "Cc": "K21610", "Cd": "K21611", "W": "K02595"}
CATALYTIC = {"D", "K", "vD", "vK", "N", "E", "Cd"}
LENGTH = {"H": 290, "vH": 290, "D": 480, "K": 510, "vD": 470, "vK": 470, "vG": 110, "N": 450,
          "E": 470, "I1": 110, "I2": 120, "B": 330, "Cc": 260, "Cd": 380, "W": 110}


def _gene(contig_id: str, index: int, token: str, strand: str) -> Gene:
    """One synthetic gene. Tokens: KO codes above, ``oNNN`` (KO-less Oxidored_nitro of NNN aa),
    ``.`` (unrelated gene); prefix ``-`` flips the strand, suffix ``!`` fails the KOfam
    threshold, suffix ``*`` marks the gene partial."""
    passes = not token.endswith("!")
    partial = token.endswith("*")
    token = token.rstrip("!*")
    if token.startswith("-"):
        strand = "-" if strand == "+" else "+"
        token = token[1:]
    ko = KO.get(token)
    pfams: set[str] = set()
    if token in CATALYTIC or token.startswith("o"):
        pfams.add("Oxidored_nitro")
    if token in ("H", "vH", "Cc"):
        pfams.add("Fer4_NifH")
    length = int(token[1:]) if token.startswith("o") else LENGTH.get(token, 200)
    start = index * 1000 + 1
    return Gene(f"{contig_id}_{index}", contig_id, index, strand, start, start + 3 * length,
                length, partial, ko, passes and ko is not None, frozenset(pfams))


def _contig(spec: list[str], contig_id: str = "c1", strand: str = "+", circular: bool = False):
    """Genes in operon order; a minus-strand operon is laid out right to left."""
    n = len(spec)
    genes = [_gene(contig_id, i if strand == "+" else n - 1 - i, token, strand) for i, token in enumerate(spec)]
    return {contig_id: Contig(contig_id, n, circular)}, genes


def _call(spec, **kwargs):
    contigs, genes = _contig(spec, **kwargs)
    return call_loci(contigs, genes)


def _only(loci):
    assert len(loci) == 1, [(locus.call, locus.locus_id) for locus in loci]
    return loci[0]


def _roles(locus):
    return {gene.protein_id: role for gene, role in locus.members}


PAD = ["."] * 8


# --------------------------------------------------------------------------- #
# Positive layouts
# --------------------------------------------------------------------------- #


def test_bacterial_nifhdk_is_confirmed():
    locus = _only(_call([".", ".", "H", "D", "K", "E", "N", *PAD]))
    assert locus.call == "confirmed"
    assert locus.locus_type == CONFIRMED_TYPE
    assert _roles(locus) == {"c1_2": "fe_protein", "c1_3": "alpha", "c1_4": "beta"}
    assert locus.metadata()["call"] == "confirmed"


def test_methanogen_type_operon_takes_the_nifn_group_beta():
    locus = _only(_call([".", "H", "I1", "I2", "D", "N", *PAD]))
    assert locus.call == "confirmed"
    assert locus.evidence["group_ii_beta"] is True
    assert _roles(locus) == {"c1_1": "fe_protein", "c1_2": "nifI", "c1_3": "nifI", "c1_4": "alpha", "c1_5": "beta"}


def test_nifn_group_beta_needs_the_nifi_genes():
    locus = _only(_call([".", "H", "D", "N", *PAD]))
    assert locus.call == "beta_slot_mismatch"
    assert locus.locus_type == CANDIDATE_TYPE


def test_minus_strand_operon_reads_right_to_left():
    locus = _only(_call([".", "H", "D", "K", *PAD], strand="-"))
    assert locus.call == "confirmed"
    assert locus.strand == "-"


def test_vanadium_operon_with_delta_subunit():
    locus = _only(_call([".", "vH", "vD", "vG", "vK", *PAD]))
    assert locus.call == "confirmed"
    assert set(_roles(locus).values()) == {"fe_protein", "alpha", "delta", "beta"}


def test_alpha_and_beta_up_to_the_window_edge():
    # nifH, four unrelated genes, nifD at +5 and nifK at +6: inside WINDOW = 6.
    assert _only(_call(["H", ".", ".", ".", ".", "D", "K", *PAD])).call == "confirmed"
    # nifK at +7 lies outside the window and forms its own gene set.
    loci = _call(["H", ".", ".", ".", ".", ".", "D", "K", *PAD])
    assert sorted(locus.call for locus in loci) == ["beta_slot_mismatch", "no_fe_protein"]


def test_circular_contig_wraps_around_its_origin():
    spec = ["K", *PAD, "H", "D"]
    assert _only(_call(spec, circular=True)).call == "confirmed"
    # Linear: nifH-nifD at one end and nifK at the other are two gene sets cut by the ends.
    assert [locus.call for locus in _call(spec, circular=False)] == ["contig_edge", "contig_edge"]


# --------------------------------------------------------------------------- #
# Paralog and partial layouts never confirm
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(("spec", "call"), [
    ([".", "H", "o348", "o377", *PAD], "short_partners"),
    ([".", "H", "o514", "o476", "B", *PAD], "unassigned_partners"),
    ([".", "H", "Cd", "o400", *PAD], "paralog_ko"),
    ([".", "H", "E", "N", *PAD], "cofactor_ko_partners"),
    ([".", "H", "o380", *PAD], "single_partner"),
    ([".", ".", ".", "o474", "o453", "H", *PAD], "upstream_partners"),
    ([".", "H", *PAD], "no_partner"),
    ([".", "H", "D", "K!", *PAD], "ko_below_threshold"),
    ([".", "H", "D", "E", "N", *PAD], "beta_slot_mismatch"),
])
def test_uncalled_layouts(spec, call):
    locus = _only(_call(spec))
    assert locus.call == call
    assert locus.locus_type == CANDIDATE_TYPE
    assert "fe_protein" not in _roles(locus).values()


def test_partners_on_the_other_strand_do_not_fill_slots():
    loci = _call([".", "H", "-D", "-K", *PAD])
    assert sorted(locus.call for locus in loci) == ["no_fe_protein", "no_partner"]


def test_cfbc_cfbd_pair_seeds_no_locus():
    assert _call([".", "Cd", "Cc", *PAD]) == []


@pytest.mark.parametrize("spec", [
    [".", ".", "H", "I1", "I2"],                 # window runs off the contig end
    [".", "H", "I1", "I2", "D", "o229*"],        # beta slot cut by the contig end
    ["H*", "I1", "I2", "D", "N", *PAD],          # Fe protein cut by the contig start
    [".", "H", "D", "K*"],                       # all slots called, the beta gene partial
])
def test_gene_sets_cut_by_a_contig_end(spec):
    locus = _only(_call(spec))
    assert locus.call == "contig_edge"
    assert locus.locus_type == CANDIDATE_TYPE


def test_two_partners_before_a_contig_end_are_classified():
    # Both slots are present, so the cut window does not hide the layout.
    assert _only(_call([".", "H", "o348", "o377"])).call == "short_partners"


def test_alpha_without_fe_protein():
    edge = _only(_call(["I2", "D", "N", *PAD]))
    assert edge.call == "contig_edge"
    assert _roles(edge) == {"c1_1": "partner", "c1_2": "partner"}  # nifD and its NifN-group neighbor
    interior = _only(_call([*PAD, "D", "K", *PAD]))
    assert interior.call == "no_fe_protein"


def test_each_fe_protein_gets_its_own_gene_set():
    loci = _call([".", "H", ".", "o488", "o436", ".", ".", "H", "o514", "o476", "B", *PAD])
    assert sorted(locus.call for locus in loci) == ["unassigned_partners", "unassigned_partners"]
    assert len({locus.locus_id for locus in loci}) == 2


# --------------------------------------------------------------------------- #
# Atoms and the composite
# --------------------------------------------------------------------------- #


def test_confirmed_structural_members_carry_the_system_atom():
    locus = _only(_call([".", "H", "I1", "I2", "D", "N", *PAD]))
    metadata = json.dumps(locus.metadata())
    by_protein = {gene.protein_id: locus_atoms(gene.protein_id, locus.locus_id, locus.locus_type, metadata)
                  for gene, _ in locus.members}
    assert {pid for pid, atoms in by_protein.items() if atoms} == {"c1_1", "c1_4", "c1_5"}
    atom = by_protein["c1_1"][0]
    assert (atom.atom_id, atom.relation, atom.source_db, atom.source_accession) == (
        SYSTEM_ATOM, ClaimRelation.implies, SOURCE_DB, locus.locus_id)


def test_contig_edge_members_carry_the_flag_only():
    locus = _only(_call([".", "H", "I1", "I2", "D", "o229*"]))
    atoms = [a for gene, _ in locus.members
             for a in locus_atoms(gene.protein_id, locus.locus_id, locus.locus_type, locus.metadata())]
    assert {a.atom_id for a in atoms} == {EDGE_ATOM}
    assert {a.protein_id for a in atoms} == {"c1_1", "c1_4", "c1_5"}
    assert all(a.relation == ClaimRelation.flags for a in atoms)


def test_other_candidates_carry_no_atom():
    locus = _only(_call([".", "H", "o348", "o377", *PAD]))
    assert not [a for gene, _ in locus.members
                for a in locus_atoms(gene.protein_id, locus.locus_id, locus.locus_type, locus.metadata())]


def _atom(atom_id, relation=ClaimRelation.implies, source_db="kegg", facet=SemanticFacet.activity):
    return SemanticAtom(protein_id="p", atom_id=atom_id, facet=facet, relation=relation,
                        source_accession="x", source_db=source_db)


def test_composite_ignores_single_gene_ko_atoms():
    clear_composites_cache()
    composites = load_composites()
    kegg_only = [_atom("nitrogenase"), _atom("nitrogen_fixation")]
    assert "nitrogen_fixation_confirmed" not in evaluate_composites(kegg_only, composites)
    assert "nitrogen_fixation_confirmed" in evaluate_composites(
        [*kegg_only, _atom(SYSTEM_ATOM, source_db=SOURCE_DB)], composites)
    flagged = [_atom(EDGE_ATOM, relation=ClaimRelation.flags, source_db=SOURCE_DB, facet=SemanticFacet.quality_flag)]
    assert "nitrogen_fixation_confirmed" not in evaluate_composites([*kegg_only, *flagged], composites)


def test_composite_requires_the_locus_caller_source():
    clear_composites_cache()
    rule = next(c for c in load_composites() if c.name == "nitrogen_fixation_confirmed")
    assert isinstance(rule, CompositeDefinition)
    assert not evaluate_composites([_atom(SYSTEM_ATOM, source_db="kegg")], [rule])


# --------------------------------------------------------------------------- #
# End to end: annotations -> loci -> atoms -> composite
# --------------------------------------------------------------------------- #


KO_LIST = (
    "knum\tthreshold\tscore_type\tprofile_type\tF-measure\tnseq\tnseq_used\talen\tmlen\teff_nseq\tre/pos\tdefinition\n"
    "K02588\t362.67\tfull\tall\t0.9\t1\t1\t1\t1\t1\t0.6\tnitrogenase iron protein NifH\n"
    "K02586\t541.60\tfull\tall\t0.9\t1\t1\t1\t1\t1\t0.6\tnitrogenase molybdenum-iron protein alpha chain\n"
    "K02591\t403.20\tdomain\tall\t0.9\t1\t1\t1\t1\t1\t0.6\tnitrogenase molybdenum-iron protein beta chain\n"
    "K02592\t471.20\tdomain\tall\t0.9\t1\t1\t1\t1\t1\t0.6\tnitrogenase molybdenum-iron protein NifN\n"
    "K02587\t501.00\tdomain\tall\t0.9\t1\t1\t1\t1\t1\t0.6\tnitrogenase molybdenum-cofactor synthesis protein NifE\n"
    "K02589\t119.40\tfull\tall\t0.9\t1\t1\t1\t1\t1\t0.6\tnitrogen regulatory protein PII 1\n"
    "K02590\t122.13\tdomain\tall\t0.9\t1\t1\t1\t1\t1\t0.6\tnitrogen regulatory protein PII 2\n"
    "K21610\t405.23\tfull\tall\t0.9\t1\t1\t1\t1\t1\t0.6\tCfbC\n"
    "K21611\t433.03\tfull\tall\t0.9\t1\t1\t1\t1\t1\t0.6\tCfbD\n"
)

# contig -> operon tokens (all plus strand); '.' genes carry no annotation.
LAYOUTS = {
    "bact_nif": [".", "H", "D", "K", "E", "N", ".", ".", ".", ".", ".", "."],
    "arch_nif": [".", "H", "I1", "I2", "D", "N", ".", ".", ".", ".", ".", "."],
    "paralog": [".", "H", "o348", "o377", ".", ".", ".", ".", ".", ".", "."],
    "edge": [".", ".", ".", "H", "I1", "I2"],
    "cfb": [".", "Cd", "Cc", ".", ".", ".", ".", ".", "."],
    "weak": [".", "H", "D", "K!", ".", ".", ".", ".", ".", ".", "."],
}


def _seed_store(tmp_path, monkeypatch) -> DuckDBStore:
    ko_list = tmp_path / "ko_list"
    ko_list.write_text(KO_LIST)
    monkeypatch.setenv("SHARUR_KOFAM_KO_LIST", str(ko_list))
    store = DuckDBStore()
    conn = store.conn
    conn.execute("INSERT INTO bins (bin_id, completeness, contamination, taxonomy) VALUES ('bin1', 90, 1, 'd__Archaea')")
    annotation_id = 0
    for contig_id, spec in LAYOUTS.items():
        conn.execute("INSERT INTO contigs (contig_id, bin_id, length, gc_content) VALUES (?, 'bin1', ?, 0.5)",
                     [contig_id, len(spec) * 1000 + 500])
        for index, token in enumerate(spec):
            gene = _gene(contig_id, index, token, "+")
            conn.execute(
                """INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, gene_index,
                                         sequence_length, gc_content, partial)
                   VALUES (?, ?, 'bin1', ?, ?, '+', ?, ?, 0.5, ?)""",
                [gene.protein_id, contig_id, gene.start, gene.end, index, gene.length,
                 "01" if gene.partial else "00"],
            )
            hits = []
            if gene.ko:
                score = 100.0 if token.endswith("!") else 900.0
                hits.append(("kegg", gene.ko, gene.ko, score, 1e-100))
            hits.extend(("pfam", pfam, pfam, 200.0, 1e-50) for pfam in sorted(gene.pfams))
            for source, accession, name, score, evalue in hits:
                annotation_id += 1
                conn.execute(
                    """INSERT INTO annotations (annotation_id, protein_id, source, accession, name, description,
                                                evalue, score)
                       VALUES (?, ?, ?, ?, ?, '', ?, ?)""",
                    [annotation_id, gene.protein_id, source, accession, name, evalue, score],
                )
    return store


def _composite_hits(store) -> set[str]:
    return {row[0] for row in store.conn.execute(
        "SELECT protein_id FROM semantic_state WHERE list_contains(composite_predicates, 'nitrogen_fixation_confirmed')"
    ).fetchall()}


def test_end_to_end_full_generation(tmp_path, monkeypatch):
    clear_composites_cache()
    store = _seed_store(tmp_path, monkeypatch)
    generate_and_persist_v2(store, chunk_size=7, update_legacy_predicates=True, return_states=False)

    loci = {row[0]: (row[1], json.loads(row[2])["call"]) for row in store.conn.execute(
        "SELECT contig_id, locus_type, metadata FROM loci WHERE locus_type LIKE 'nitrogenase%'").fetchall()}
    assert loci == {
        "bact_nif": (CONFIRMED_TYPE, "confirmed"),
        "arch_nif": (CONFIRMED_TYPE, "confirmed"),
        "paralog": (CANDIDATE_TYPE, "short_partners"),
        "edge": (CANDIDATE_TYPE, "contig_edge"),
        "weak": (CANDIDATE_TYPE, "ko_below_threshold"),
    }
    # Fe protein, alpha and beta of the two confirmed gene sets only: never NifE, the
    # unattached NifN, nifI, the paralog pair, the CfbC/CfbD pair or the truncated set.
    assert _composite_hits(store) == {"bact_nif_1", "bact_nif_2", "bact_nif_3",
                                      "arch_nif_1", "arch_nif_4", "arch_nif_5"}

    flags = {row[0] for row in store.conn.execute(
        "SELECT protein_id FROM semantic_state WHERE list_contains(quality_flags, ?)", [EDGE_ATOM]).fetchall()}
    assert flags == {"edge_3"}
    atom = store.conn.execute(
        "SELECT relation, source_db, source_accession FROM semantic_atoms WHERE protein_id = 'arch_nif_5' AND atom_id = ?",
        [SYSTEM_ATOM]).fetchone()
    assert atom[:2] == ("implies", SOURCE_DB)
    assert atom[2].startswith("nitrogenase:arch_nif:")
    legacy = store.conn.execute("SELECT predicates FROM protein_predicates WHERE protein_id = 'bact_nif_2'").fetchone()[0]
    assert "nitrogen_fixation_confirmed" in legacy

    # A composite-only refresh reads the persisted locus atoms.
    store.conn.execute("UPDATE semantic_state SET composite_predicates = []")
    refresh_composite_predicates_from_v2(store)
    assert len(_composite_hits(store)) == 6


def test_end_to_end_subset_generation_and_rerun(tmp_path, monkeypatch):
    clear_composites_cache()
    store = _seed_store(tmp_path, monkeypatch)
    states = generate_and_persist_v2(store, protein_ids=["arch_nif_5", "arch_nif_2", "paralog_1"],
                                     chunk_size=2, return_states=True)
    assert "nitrogen_fixation_confirmed" in states["arch_nif_5"].composite_predicates
    assert "nitrogen_fixation_confirmed" not in states["arch_nif_2"].composite_predicates  # nifI1
    assert "nitrogen_fixation_confirmed" not in states["paralog_1"].composite_predicates
    first = store.conn.execute("SELECT locus_id, metadata FROM loci ORDER BY 1").fetchall()
    generate_and_persist_v2(store, protein_ids=["arch_nif_5"], return_states=False)
    assert store.conn.execute("SELECT locus_id, metadata FROM loci ORDER BY 1").fetchall() == first
    members = store.conn.execute("SELECT COUNT(*) FROM locus_proteins").fetchone()[0]
    assert members == sum(len(json.loads(m)["members"]) for _, m in first)


def test_unchanged_calls_leave_the_loci_tables_untouched(tmp_path, monkeypatch):
    # The generation's prefetch reader reads `loci` while the writer commits chunks, so
    # generation checkpoints after a rewrite and skips the rewrite when nothing changed.
    store = _seed_store(tmp_path, monkeypatch)
    first = refresh_nitrogenase_loci(store)
    assert first["written"] is True
    assert first["by_call"] == {"confirmed": 2, "contig_edge": 1, "ko_below_threshold": 1, "short_partners": 1}
    again = refresh_nitrogenase_loci(store)
    assert again["written"] is False
    assert again["by_call"] == first["by_call"]
