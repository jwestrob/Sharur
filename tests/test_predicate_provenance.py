"""Predicate generations record which maps produced them."""

import pytest

import sharur.predicates_v2.persistence as persistence
from sharur.capabilities import CapabilityState, _predicate_map_check
from sharur.predicates import provenance
from sharur.predicates.mappings import kegg_map
from sharur.predicates.provenance import current_maps, latest_stamps, map_status
from sharur.predicates_v2.persistence import generate_and_persist_v2
from sharur.storage.duckdb_store import DuckDBStore
from sharur.storage.schema import SCHEMA_VERSION


def _store() -> DuckDBStore:
    store = DuckDBStore()
    store.conn.execute("""
        INSERT INTO bins (bin_id, completeness, contamination, taxonomy) VALUES ('bin1', 90.0, 2.0, 'd__Bacteria');
        INSERT INTO contigs (contig_id, bin_id, length, gc_content) VALUES ('contig1', 'bin1', 10000, 0.5);
        INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, gene_index,
                              sequence_length, gc_content)
        VALUES ('p1', 'contig1', 'bin1', 1, 300, '+', 1, 100, 0.5),
               ('p2', 'contig1', 'bin1', 301, 900, '+', 2, 200, 0.5);
        INSERT INTO annotations (annotation_id, protein_id, source, accession, name, description, evalue, score)
        VALUES (1, 'p1', 'pfam', 'PF00005', 'ABC_tran', 'ABC transporter', 1e-30, 100.0);
    """)
    return store


def test_new_databases_carry_the_provenance_table():
    store = DuckDBStore()
    assert store.execute("SELECT MAX(version) FROM schema_version")[0][0] == SCHEMA_VERSION == 7
    assert latest_stamps(store) == {"full": None, "subsets": []}


def test_generations_are_stamped_with_the_installed_maps():
    store = _store()
    assert map_status(store).state == "unstamped"

    generate_and_persist_v2(store, chunk_size=1, return_states=False)
    stamps = latest_stamps(store)
    full = stamps["full"]
    installed = current_maps()
    assert (full["scope"], full["protein_count"]) == ("full", 2)
    for field in provenance.COMPARED:
        assert full[field] == installed[field], field
    assert full["pfam_sources"].startswith("Pfam:")
    assert full["semantic_fingerprint"]
    assert map_status(store).state == "current"

    generate_and_persist_v2(store, protein_ids=["p1"], chunk_size=1, return_states=False)
    stamps = latest_stamps(store)
    assert [s["scope"] for s in stamps["subsets"]] == ["subset"]
    assert stamps["subsets"][0]["protein_count"] == 1


def test_a_changed_map_marks_predicates_stale():
    store = _store()
    generate_and_persist_v2(store, chunk_size=1, return_states=False)
    store.execute("UPDATE predicate_provenance SET pfam_map_sha256 = 'older-map'")
    status = map_status(store)
    assert (status.state, status.changed) == ("stale", ("pfam_map_sha256",))
    check = _predicate_map_check(store)
    assert check.state == CapabilityState.stale
    assert "Pfam map" in check.summary
    assert check.remediation


def test_subsets_on_other_maps_are_counted():
    store = _store()
    generate_and_persist_v2(store, chunk_size=1, return_states=False)
    generate_and_persist_v2(store, protein_ids=["p1"], chunk_size=1, return_states=False)
    store.execute("UPDATE predicate_provenance SET vocabulary_sha256 = 'other' WHERE scope = 'subset'")
    status = map_status(store)
    assert (status.state, status.mixed_subsets) == ("current", 1)


def test_unstamped_databases_are_reported_as_stale():
    check = _predicate_map_check(_store())
    assert check.state == CapabilityState.stale
    assert "no map provenance stamp" in check.summary


def test_kegg_map_changes_the_semantic_fingerprint(tmp_path, monkeypatch):
    fingerprint = lambda: persistence._semantic_generation_fingerprint(  # noqa: E731
        predict_topology=False, update_legacy_predicates=False)
    first = tmp_path / "a.tsv"
    first.write_text("K90001\tkinase\tkinase=ec:2.7.1.1\n")
    monkeypatch.setattr(kegg_map, "KEGG_MAPPING_FILE", first)
    a = fingerprint()
    first.write_text("K90001\tkinase,transferase\tkinase=ec:2.7.1.1;transferase=ec:2.7.1.1\n")
    assert fingerprint() != a
    monkeypatch.setattr(kegg_map, "KEGG_MAPPING_FILE", None)
    assert fingerprint() != a


@pytest.mark.parametrize("has_kegg", [True, False])
def test_stamp_records_kegg_build_presence(tmp_path, monkeypatch, has_kegg):
    if has_kegg:
        built = tmp_path / "kegg_predicates.tsv"
        built.write_text("K90001\tkinase\tkinase=ec:2.7.1.1\n")
        (tmp_path / "provenance.json").write_text('{"kegg_release": "ko 2026/01/01"}')
        monkeypatch.setattr(kegg_map, "KEGG_MAPPING_FILE", built)
    else:
        monkeypatch.setattr(kegg_map, "KEGG_MAPPING_FILE", None)
    maps = current_maps()
    assert (maps["kegg_map_sha256"] is not None) == has_kegg
    assert maps["kegg_release"] == ("ko 2026/01/01" if has_kegg else None)
