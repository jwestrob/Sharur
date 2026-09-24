"""Hydrogenase subgroup interpretation, reconciliation, and staged refresh."""

import hashlib
from pathlib import Path

import duckdb
import pytest

from sharur.hydrogenase import classifier as hyd
from sharur.hydrogenase.classifier import (
    ASSIGNED,
    CLASS_CONFLICT,
    CLEARED,
    MISSING_SEQUENCE,
    NEEDS_CURATION,
    NO_REFERENCE_HIT,
    UNPARSED_REFERENCE_LABEL,
    HydrogenaseSearchError,
    ReferenceHit,
    classify,
    find_reference,
    parse_reference_id,
)
from sharur.hydrogenase.ko_association import (
    COMPATIBLE,
    CONFLICT,
    GROUP,
    NONE,
    SUBGROUP,
    associations,
    ko_support,
    load_associations,
)
from sharur.hydrogenase.refresh import refresh_hydrogenases
from sharur.hydrogenase.subgroups import (
    CONTEXT_DEPENDENT_SUBGROUPS,
    PUTATIVE,
    RETIRED_TERMS,
    SUBGROUPS,
    UNRESOLVED,
    UNVERIFIED,
    lookup,
)
from sharur.predicates.generator import AnnotationRecord, PredicateGenerator, ProteinRecord
from sharur.predicates.vocabulary import PREDICATE_BY_ID, get_hierarchy
from sharur.predicates_v2.model import ClaimRelation
from sharur.predicates_v2.persistence import generate_and_persist_v2
from sharur.predicates_v2.rules import clear_caches, get_facet, get_relation
from sharur.storage.duckdb_store import DuckDBStore


@pytest.fixture(autouse=True)
def _reset_caches():
    clear_caches()
    yield
    clear_caches()


# --------------------------------------------------------------------------- #
# Subgroup interpretations and ontology
# --------------------------------------------------------------------------- #


def _preds(hyd_type, subgroup):
    return set(lookup(hyd_type, subgroup).predicates)


class TestSubgroupInterpretations:
    def test_2a_is_respiratory_uptake_not_a_sensor(self):
        assert "uptake_hydrogenase" in _preds("NiFe", "Group_2a")
        assert "h2_sensor" not in _preds("NiFe", "Group_2a")
        assert "h2_sensor" in _preds("NiFe", "Group_2b")
        for sub in ("Group_2c", "Group_2d", "Group_2e"):
            assert lookup("NiFe", sub).functional == ()

    def test_3b_nadp_and_3d_nad_are_distinct(self):
        assert "nadp_coupled" in _preds("NiFe", "Group_3b")
        assert "nad_coupled" not in _preds("NiFe", "Group_3b")
        assert "nad_coupled" in _preds("NiFe", "Group_3d")
        assert "nadp_coupled" not in _preds("NiFe", "Group_3d")

    def test_3c_is_heterodisulfide_linked(self):
        preds = _preds("NiFe", "Group_3c")
        assert "heterodisulfide_reductase_linked" in preds
        assert "methyl_viologen_reducing" not in preds

    def test_group4_ferredoxin_and_ech_assignments(self):
        assert "co_coupled" not in _preds("NiFe", "Group_4d")
        assert "ferredoxin_coupled" in _preds("NiFe", "Group_4d")
        assert "ech_hydrogenase" in _preds("NiFe", "Group_4e")
        assert lookup("NiFe", "Group_4g").functional == ()
        assert "ech_hydrogenase" not in _preds("NiFe", "Group_4g")
        assert "co_coupled" in _preds("NiFe", "Group_4c")

    def test_fefe_group_a_subtypes_carry_structure_only(self):
        for key in CONTEXT_DEPENDENT_SUBGROUPS:
            assert SUBGROUPS[key].functional == ()
            assert SUBGROUPS[key].structural == ("fefe_groupA",)
        assert lookup("FeFe", "Group_A3").reference_name == "Bifurcating"

    def test_fefe_b_and_c_carry_no_mechanism(self):
        for sub in ("Group_B", "Group_C1", "Group_C2", "Group_C3"):
            assert lookup("FeFe", sub).functional == ()
            assert lookup("FeFe", sub).status == PUTATIVE

    def test_uncharacterized_subgroups_emit_no_traits(self):
        for sub in SUBGROUPS.values():
            if sub.status in (PUTATIVE, UNRESOLVED, UNVERIFIED):
                assert sub.functional == (), sub.label

    def test_1l_is_retained_as_unverified(self):
        sub = lookup("NiFe", "Group_1l")
        assert sub.status == UNVERIFIED
        assert sub.predicates == ("nife_group1",)

    def test_every_emitted_predicate_is_in_the_vocabulary_and_faceted(self):
        emitted = {p for sub in SUBGROUPS.values() for p in sub.predicates}
        emitted |= {"hyddb_needs_curation", "hyddb_class_conflict", "hydrogenase_complex1_review"}
        assert emitted <= set(PREDICATE_BY_ID)
        assert not emitted & RETIRED_TERMS
        for pred in emitted:
            assert get_facet(pred) is not None

    def test_retired_terms_left_the_vocabulary(self):
        assert not RETIRED_TERMS & set(PREDICATE_BY_ID)

    def test_installed_reference_labels_are_all_interpreted(self):
        faa = Path(__file__).resolve().parents[1] / "data/reference/hyddb/HydDB_all_hydrogenases.faa"
        if not faa.exists():
            pytest.skip("HydDB reference not installed")
        labels = {line.rstrip().split("|")[-1] for line in faa.open() if line.startswith(">")}
        for label in labels:
            _, _, _, hyd_type, subgroup = parse_reference_id(f"x|y|{label}")
            assert lookup(hyd_type, subgroup) is not None, label


class TestHierarchy:
    def test_functional_traits_do_not_imply_groups(self):
        assert "nife_group1" not in get_hierarchy("uptake_hydrogenase")
        assert "nife_group2" not in get_hierarchy("h2_sensor")
        assert "fefe_groupB" not in get_hierarchy("bifurcating_hydrogenase")
        for trait in ("nad_coupled", "nadp_coupled", "f420_reducing", "ferredoxin_coupled",
                      "heterodisulfide_reductase_linked", "formate_coupled", "co_coupled",
                      "energy_conserving_hydrogenase"):
            assert not any(a.startswith(("nife_group", "fefe_group")) for a in get_hierarchy(trait)), trait

    def _expanded(self, labels):
        gen = PredicateGenerator()
        anns = [AnnotationRecord(source="hyddb_subgroup", accession=x, name=x, evalue=1e-80, score=500.0)
                for x in labels]
        return set(gen.generate_for_protein(ProteinRecord(protein_id="p", sequence_length=600), anns))

    def test_2a_uptake_expansion_stays_in_group_2(self):
        preds = self._expanded(lookup("NiFe", "Group_2a").predicates)
        assert {"nife_group2", "uptake_hydrogenase"} <= preds
        assert "nife_group1" not in preds

    def test_3c_bifurcation_expansion_stays_in_nife(self):
        preds = self._expanded(lookup("NiFe", "Group_3c").predicates)
        assert "bifurcating_hydrogenase" in preds
        assert "fefe_groupB" not in preds and "fefe_hydrogenase" not in preds

    def test_a3_expansion_stays_in_group_a(self):
        preds = self._expanded(lookup("FeFe", "Group_A3").predicates)
        assert "fefe_groupA" in preds
        assert not {"fefe_groupB", "bifurcating_hydrogenase"} & preds


class TestProvisionalRelations:
    def test_subgroup_matches_support_rather_than_imply(self):
        assert get_relation("hyddb_subgroup", "nife_group1", "nife_group1") == ClaimRelation.supports

    def test_review_flags_are_flags(self):
        for flag in ("hyddb_needs_curation", "hyddb_class_conflict", "hyddb_ko_supported", "hyddb_ko_conflict"):
            assert get_relation("hyddb_subgroup", flag, flag) == ClaimRelation.flags


# --------------------------------------------------------------------------- #
# Database fixtures
# --------------------------------------------------------------------------- #

SEQ = "M" + "ACDEFGHIKLMNPQRSTVWY" * 20

PROTEINS = [
    "p_1h", "p_2a", "p_3d", "p_4g_kegg", "p_a3", "p_conflict",
    "p_nohit", "p_noseq", "p_unparsed", "p_unrelated",
]

# A KEGG ortholog that independently supports ech_hydrogenase. Patched into the
# KEGG map so the provenance test exercises the mechanism, not a table entry.
INDEPENDENT_KO = "K99901"


@pytest.fixture(autouse=True)
def _independent_kegg_support(monkeypatch):
    from sharur.predicates.mappings import kegg_map

    monkeypatch.setitem(kegg_map.KEGG_TO_PREDICATES, INDEPENDENT_KO,
                        ["hydrogenase", "nife_hydrogenase", "ech_hydrogenase", "hydrogen_metabolism"])


ANNOTATIONS = [
    # protein, source, accession, name, score
    ("p_1h", "hyddb", "NiFe", "NiFe", 300.0),
    ("p_1h", "hyddb", "NiFe", "NiFe", 120.0),          # second HMM row, same protein
    ("p_1h", "pfam", "NiFeSe_Hases", "NiFeSe_Hases", 400.0),  # name stored as accession
    ("p_2a", "hyddb", "NiFe", "NiFe", 280.0),
    ("p_2a", "pfam", "PF00374", "NiFeSe_Hases", 380.0),
    ("p_3d", "hyddb", "NiFe", "NiFe", 250.0),
    ("p_3d", "pfam", "PF00374.26", "NiFeSe_Hases", 350.0),
    ("p_4g_kegg", "hyddb", "NiFe", "NiFe", 200.0),
    ("p_4g_kegg", "pfam", "Complex1_49kDa", "Complex1_49kDa", 250.0),
    ("p_4g_kegg", "kegg", INDEPENDENT_KO, "echA", 180.0),  # independent Ech support
    ("p_a3", "hyddb", "FeFe", "FeFe", 310.0),
    ("p_a3", "pfam", "PF02906", "Fe_hyd_lg_C", 330.0),
    ("p_conflict", "hyddb", "FeFe", "FeFe", 150.0),
    ("p_nohit", "hyddb", "NiFe", "NiFe", 90.0),
    ("p_noseq", "hyddb", "NiFe", "NiFe", 95.0),
    ("p_unparsed", "hyddb", "NiFe", "NiFe", 99.0),
    ("p_unrelated", "pfam", "PF00005", "ABC_tran", 110.0),
    # KOfam hits whose HydDB associations support or question the assignment
    ("p_3d", "kofam", "K00436", "hoxH", 600.0),      # captures Group_3d 174/174 -> subgroup
    ("p_1h", "kofam", "K15830", "hycE", 500.0),      # captures Group_4a 47/47 -> conflict
    ("p_2a", "kofam", "K23549", "hupV", 450.0),      # Group_2a/2b/2c/2e -> compatible
]

# Labels the pre-audit classifier wrote (seeded as stale state)
LEGACY_LABELS = {
    "p_3d": ["nife_group3", "bidirectional_hydrogenase", "nadp_reducing", "hyddb_needs_curation"],
    "p_4g_kegg": ["nife_group4", "h2_evolving", "energy_conserving_hydrogenase", "ech_hydrogenase",
                  "hyddb_needs_curation"],
    "p_a3": ["fefe_groupA", "monomeric_fefe", "fermentative_hydrogenase", "hyddb_needs_curation"],
}

HITS = {
    "p_1h": ("WP_1|Org|[NiFe]_Group_1h", 72.0),
    "p_2a": ("WP_2|Org|[NiFe]_Group_2a", 65.0),
    "p_3d": ("WP_3|Org|[NiFe]_Group_3d", 58.0),
    "p_4g_kegg": ("WP_4|Org|[NiFe]_Group_4g", 44.0),
    "p_a3": ("WP_5|Org|[FeFe]_Group_A3", 51.0),
    "p_conflict": ("WP_6|Org|[NiFe]_Group_1a", 38.0),
    "p_unparsed": ("WP_7|Org|mystery", 30.0),
    "p_noseq": ("WP_8|Org|[NiFe]_Group_1a", 99.0),
}


def stub_search(sequences, reference, threads):
    hits = {}
    for pid, (sid, pident) in HITS.items():
        if pid in sequences:
            acc, org, label, t, s = parse_reference_id(sid)
            hits[pid] = ReferenceHit(acc, org, label, t, s, pident, 1e-50, 400.0)
    return hits


def failing_search(sequences, reference, threads):
    raise HydrogenaseSearchError("DIAMOND exited 1: simulated")


@pytest.fixture
def reference_dir(tmp_path):
    ref = tmp_path / "hyddb"
    ref.mkdir()
    (ref / "HydDB_all.dmnd").write_bytes(b"stub reference")
    (ref / "NiFe-HydDB_MM2022.hmm").write_text("")
    return ref


@pytest.fixture
def db_path(tmp_path):
    path = tmp_path / "data" / "sharur.duckdb"
    path.parent.mkdir()
    store = DuckDBStore(path)
    store.execute("INSERT INTO bins (bin_id) VALUES ('bin_1')")
    store.execute("INSERT INTO contigs (contig_id, bin_id, length) VALUES ('c1', 'bin_1', 100000)")
    for i, pid in enumerate(PROTEINS):
        store.execute(
            """INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand,
               gene_index, sequence, sequence_length) VALUES (?, 'c1', 'bin_1', ?, ?, '+', ?, ?, ?)""",
            [pid, i * 2000 + 1, i * 2000 + 1200, i, None if pid == "p_noseq" else SEQ, len(SEQ)],
        )
    rows = [(i + 1, p, s, a, n, "", 1e-40, sc) for i, (p, s, a, n, sc) in enumerate(ANNOTATIONS)]
    next_id = len(rows)
    for pid, labels in LEGACY_LABELS.items():
        for label in labels:
            next_id += 1
            rows.append((next_id, pid, "hyddb_subgroup", label, label, "HydDB DIAMOND: legacy", 1e-50, 400.0))
    store.conn.executemany(
        "INSERT INTO annotations (annotation_id, protein_id, source, accession, name, description, evalue, score) "
        "VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
        rows,
    )
    store.execute("CREATE TABLE rna_expression (protein_id VARCHAR, sample VARCHAR, tpm DOUBLE)")
    store.execute("INSERT INTO rna_expression VALUES ('p_1h', 'PLM2_5cm', 225.73), ('p_unrelated', 'PLM1', 3.0)")
    generate_and_persist_v2(store, update_legacy_predicates=True, return_states=False)
    store.close()
    return path


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _query(path, sql, params=()):
    conn = duckdb.connect(str(path), read_only=True)
    try:
        return conn.execute(sql, list(params)).fetchall()
    finally:
        conn.close()


def _atoms(path, pid):
    return {r[0] for r in _query(path, "SELECT atom_id FROM semantic_atoms WHERE protein_id = ?", [pid])}


def _legacy(path, pid):
    return set(_query(path, "SELECT predicates FROM protein_predicates WHERE protein_id = ?", [pid])[0][0])


# --------------------------------------------------------------------------- #
# Classification and reconciliation
# --------------------------------------------------------------------------- #


class TestClassification:
    @pytest.fixture
    def results(self, db_path, reference_dir):
        conn = duckdb.connect(str(db_path))
        try:
            items = classify(conn, find_reference(reference_dir), search=stub_search)
        finally:
            conn.close()
        return {item.protein_id: item for item in items}

    def test_one_record_per_protein_keeps_all_discovery_rows(self, results):
        assert len(results) == 9
        assert results["p_1h"].discovery_hit_count == 2
        assert results["p_1h"].discovery_best_score == 300.0
        assert results["p_1h"].discovery_classes == ("NiFe",)

    def test_distinguishable_outcomes(self, results):
        assert results["p_1h"].outcome == ASSIGNED
        assert results["p_conflict"].outcome == CLASS_CONFLICT
        assert results["p_nohit"].outcome == NO_REFERENCE_HIT
        assert results["p_noseq"].outcome == MISSING_SEQUENCE
        assert results["p_unparsed"].outcome == UNPARSED_REFERENCE_LABEL

    def test_class_conflict_keeps_both_classes_and_emits_no_subgroup(self, results):
        item = results["p_conflict"]
        assert item.hit.hyd_type == "NiFe" and item.discovery_classes == ("FeFe",)
        assert item.derived_labels == ("hyddb_class_conflict",)
        assert "FeFe" in item.curation_reason and "NiFe" in item.curation_reason

    def test_equivalent_identifiers_clear_the_same_check(self, results):
        # name-as-accession, canonical accession, versioned accession
        for pid in ("p_1h", "p_2a", "p_3d"):
            assert results[pid].curation_status == CLEARED, pid
        assert results["p_a3"].curation_status == CLEARED

    def test_complex1_only_is_reviewed(self, results):
        item = results["p_4g_kegg"]
        assert item.curation_status == NEEDS_CURATION
        assert "Complex1" in item.curation_reason
        assert "hyddb_needs_curation" in item.derived_labels

    def test_provisional_match_retains_evidence(self, results):
        item = results["p_a3"]
        assert item.hit.reference_accession == "WP_5"
        assert item.hit.label == "[FeFe]_Group_A3"
        assert item.subgroup.reference_name == "Bifurcating"
        assert "bifurcating_hydrogenase" not in item.derived_labels

    def test_ko_support_is_recorded_beside_the_assignment(self, results):
        agree, conflict, mixed = results["p_3d"], results["p_1h"], results["p_2a"]
        assert (agree.ko_support.status, agree.ko_support.detail) == (SUBGROUP, "K00436 [NiFe]_Group_3d 174/174")
        assert "hyddb_ko_supported" in agree.derived_labels
        assert conflict.ko_support.status == CONFLICT
        assert "hyddb_ko_conflict" in conflict.derived_labels
        assert mixed.ko_support.status == COMPATIBLE
        assert not {"hyddb_ko_supported", "hyddb_ko_conflict"} & set(mixed.derived_labels)
        assert results["p_a3"].ko_support.status == NONE

    def test_ko_support_leaves_assignment_and_curation_unchanged(self, results):
        item = results["p_1h"]
        assert (item.outcome, item.hit.subgroup, item.curation_status) == (ASSIGNED, "Group_1h", CLEARED)
        assert {"nife_group1", "uptake_hydrogenase"} <= set(item.derived_labels)

    def test_raw_hmm_rows_are_read_only(self, db_path, reference_dir):
        before = _query(db_path, "SELECT * FROM annotations WHERE source = 'hyddb' ORDER BY annotation_id")
        hyd.classify_database(db_path, reference_dir=reference_dir, search=stub_search)
        after = _query(db_path, "SELECT * FROM annotations WHERE source = 'hyddb' ORDER BY annotation_id")
        assert before == after


class TestDiamondFailure:
    def test_nonzero_exit_raises(self, monkeypatch, reference_dir):
        class Result:
            returncode, stdout, stderr = 1, "", "boom"

        monkeypatch.setattr(hyd.shutil, "which", lambda _: "/usr/bin/diamond")
        monkeypatch.setattr(hyd.subprocess, "run", lambda *a, **k: Result())
        with pytest.raises(HydrogenaseSearchError, match="boom"):
            hyd.run_diamond({"p": SEQ}, find_reference(reference_dir), 1)

    def test_success_with_no_hits_is_empty(self, monkeypatch, reference_dir):
        class Result:
            returncode, stdout, stderr = 0, "", ""

        monkeypatch.setattr(hyd.shutil, "which", lambda _: "/usr/bin/diamond")
        monkeypatch.setattr(hyd.subprocess, "run", lambda *a, **k: Result())
        assert hyd.run_diamond({"p": SEQ}, find_reference(reference_dir), 1) == {}

    def test_best_hsp_wins(self, monkeypatch, reference_dir):
        class Result:
            returncode, stderr = 0, ""
            stdout = ("p\tWP_1|O|[NiFe]_Group_1h\t70\t1e-90\t500\n"
                      "p\tWP_1|O|[NiFe]_Group_1h\t40\t1e-5\t60\n")

        monkeypatch.setattr(hyd.shutil, "which", lambda _: "/usr/bin/diamond")
        monkeypatch.setattr(hyd.subprocess, "run", lambda *a, **k: Result())
        assert hyd.run_diamond({"p": SEQ}, find_reference(reference_dir), 1)["p"].bitscore == 500


# --------------------------------------------------------------------------- #
# Staged refresh
# --------------------------------------------------------------------------- #


class TestRefresh:
    def test_dry_run_reports_and_leaves_database_unchanged(self, db_path, reference_dir, tmp_path):
        digest = _sha(db_path)
        tsv = tmp_path / "transitions.tsv"
        report = refresh_hydrogenases(db_path, dry_run=True, reference_dir=reference_dir,
                                      search=stub_search, transitions_path=tsv)
        assert _sha(db_path) == digest
        assert not list(db_path.parent.glob("*.staging*"))
        assert all(v.startswith("PASS") for v in report.checks.values())
        assert report.proteins_classified == 9
        by_pid = {t.protein_id: t for t in report.transitions}
        assert by_pid["p_3d"].curation_before and not by_pid["p_3d"].curation_after
        assert "nadp_reducing" in set(by_pid["p_3d"].before) - set(by_pid["p_3d"].after)
        assert tsv.read_text().count("\n") == len(report.transitions) + 1

    def test_publish_replaces_stale_claims_and_keeps_independent_ones(self, db_path, reference_dir):
        assert "nadp_reducing" in _atoms(db_path, "p_3d")
        assert "monomeric_fefe" in _atoms(db_path, "p_a3")
        unrelated_before = _query(db_path, "SELECT * FROM semantic_atoms WHERE protein_id = 'p_unrelated' ORDER BY ALL")
        rna_before = _query(db_path, "SELECT * FROM rna_expression ORDER BY ALL")
        raw_before = _query(db_path, "SELECT * FROM annotations WHERE source <> 'hyddb_subgroup' ORDER BY ALL")

        report = refresh_hydrogenases(db_path, dry_run=False, reference_dir=reference_dir, search=stub_search)

        assert report.backup is not None and report.backup.exists()
        assert "nadp_reducing" not in _atoms(db_path, "p_3d")
        assert "nad_coupled" in _atoms(db_path, "p_3d")
        assert "nadp_reducing" not in _legacy(db_path, "p_3d")
        assert not {"monomeric_fefe", "fermentative_hydrogenase", "bifurcating_hydrogenase",
                    "fefe_groupB"} & _atoms(db_path, "p_a3")
        # HydDB no longer names p_4g_kegg as Ech, but KEGG echA still does
        hyd_ech = _query(db_path, "SELECT COUNT(*) FROM semantic_atoms WHERE protein_id = 'p_4g_kegg' "
                                  "AND atom_id = 'ech_hydrogenase' AND source_db = 'hyddb_subgroup'")[0][0]
        assert hyd_ech == 0
        assert "ech_hydrogenase" in _atoms(db_path, "p_4g_kegg")
        assert "ech_hydrogenase" in _legacy(db_path, "p_4g_kegg")
        # 2a uptake without group 1
        assert "uptake_hydrogenase" in _legacy(db_path, "p_2a")
        assert "nife_group1" not in _legacy(db_path, "p_2a")
        # Preserved
        assert _query(db_path, "SELECT * FROM semantic_atoms WHERE protein_id = 'p_unrelated' ORDER BY ALL") == unrelated_before
        assert _query(db_path, "SELECT * FROM rna_expression ORDER BY ALL") == rna_before
        assert _query(db_path, "SELECT * FROM annotations WHERE source <> 'hyddb_subgroup' ORDER BY ALL") == raw_before
        rows = _query(db_path, "SELECT outcome, reference_release, reference_sha256, classifier_version "
                               "FROM hydrogenase_classifications WHERE protein_id = 'p_1h'")
        assert rows == [(ASSIGNED, "MM2022", hashlib.sha256(b"stub reference").hexdigest(),
                         hyd.CLASSIFIER_VERSION)]
        # KOfam support: recorded per protein and emitted as flag labels only
        support = dict(_query(db_path, "SELECT protein_id, ko_support FROM hydrogenase_classifications"))
        assert (support["p_3d"], support["p_1h"], support["p_2a"]) == (SUBGROUP, CONFLICT, COMPATIBLE)
        assert "hyddb_ko_supported" in _legacy(db_path, "p_3d")
        relations = set(_query(db_path, "SELECT relation FROM semantic_atoms WHERE protein_id = 'p_1h' "
                                        "AND atom_id = 'hyddb_ko_conflict'"))
        assert relations == {("flags",)}

    def test_publish_replaces_a_previous_classification_schema(self, db_path, reference_dir):
        conn = duckdb.connect(str(db_path))
        conn.execute("CREATE TABLE hydrogenase_classifications (protein_id VARCHAR PRIMARY KEY, outcome VARCHAR)")
        conn.execute("INSERT INTO hydrogenase_classifications VALUES ('stale', 'assigned')")
        conn.close()
        refresh_hydrogenases(db_path, dry_run=False, reference_dir=reference_dir, search=stub_search)
        columns = {r[0] for r in _query(db_path, "SELECT column_name FROM information_schema.columns "
                                                 "WHERE table_name = 'hydrogenase_classifications'")}
        assert {"ko_support", "ko_support_detail"} <= columns
        assert _query(db_path, "SELECT COUNT(*) FROM hydrogenase_classifications WHERE protein_id = 'stale'") == [(0,)]

    def test_refresh_is_idempotent(self, db_path, reference_dir):
        refresh_hydrogenases(db_path, dry_run=False, reference_dir=reference_dir, search=stub_search)
        state = _query(db_path, "SELECT * FROM semantic_atoms ORDER BY ALL")
        second = refresh_hydrogenases(db_path, dry_run=False, reference_dir=reference_dir, search=stub_search)
        assert second.changed == []
        assert _query(db_path, "SELECT * FROM semantic_atoms ORDER BY ALL") == state

    def test_search_failure_leaves_production_intact(self, db_path, reference_dir):
        digest = _sha(db_path)
        with pytest.raises(HydrogenaseSearchError):
            refresh_hydrogenases(db_path, dry_run=False, reference_dir=reference_dir, search=failing_search)
        assert _sha(db_path) == digest
        assert not list(db_path.parent.glob("*.staging*"))
        assert not list(db_path.parent.glob("*.pre-hydrogenase-refresh-*"))

    def test_failure_between_staging_and_publication_rolls_back(self, db_path, reference_dir):
        digest = _sha(db_path)

        def explode(staging):
            assert staging.exists()
            raise RuntimeError("injected publication failure")

        with pytest.raises(RuntimeError, match="injected"):
            refresh_hydrogenases(db_path, dry_run=False, reference_dir=reference_dir,
                                 search=stub_search, _before_publish=explode)
        assert _sha(db_path) == digest
        assert not list(db_path.parent.glob("*.staging*"))

    def test_existing_seal_is_refreshed(self, db_path, reference_dir):
        from sharur.dataset_seal import seal_dataset, verify_dataset_seal

        seal_path, _ = seal_dataset(db_path)
        report = refresh_hydrogenases(db_path, dry_run=False, reference_dir=reference_dir, search=stub_search)
        assert report.resealed
        assert verify_dataset_seal(seal_path).valid

    def test_unmigrated_database_is_refused_untouched(self, db_path, reference_dir):
        from sharur.hydrogenase.refresh import RefreshValidationError

        conn = duckdb.connect(str(db_path))
        conn.execute("DELETE FROM schema_version WHERE version = (SELECT MAX(version) FROM schema_version)")
        conn.close()
        digest = _sha(db_path)
        with pytest.raises(RefreshValidationError, match="sharur migrate"):
            refresh_hydrogenases(db_path, dry_run=True, reference_dir=reference_dir, search=stub_search)
        assert _sha(db_path) == digest


# --------------------------------------------------------------------------- #
# KO -> HydDB subgroup associations
# --------------------------------------------------------------------------- #


class TestKOAssociation:
    def test_associations_come_from_the_snapshot(self):
        rows = associations("K14090")  # echE
        assert [(a.label, a.count, a.total) for a in rows] == [
            ("[NiFe]_Group_4e", 130, 170), ("[NiFe]_Group_4c", 27, 170), ("[NiFe]_Group_4g", 13, 170)]
        assert associations("K00001") == ()

    @pytest.mark.parametrize(("kos", "subgroup", "status"), [
        (["K15830"], "Group_4a", SUBGROUP),
        (["K14090"], "Group_4e", COMPATIBLE),      # 130/170 falls short of the agreement threshold
        (["K14090"], "Group_4d", GROUP),
        (["K15830"], "Group_1h", CONFLICT),
        (["K15830", "K06281"], "Group_1h", COMPATIBLE),  # conflict needs every associated KO to conflict
        (["K00001"], "Group_1a", NONE),
        ([], "Group_1a", NONE),
    ])
    def test_support_grades(self, kos, subgroup, status):
        assert ko_support(kos, "NiFe", subgroup).status == status

    def test_every_associated_label_is_interpreted(self):
        labels = {(a.hyd_type, a.subgroup) for rows in load_associations().values() for a in rows}
        assert labels <= set(SUBGROUPS)
