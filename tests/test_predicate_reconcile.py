"""A KEGG-map-only change can be reconciled by regenerating exactly the proteins whose KOs changed."""

import json

import pytest

from sharur.capabilities import CapabilityState, _predicate_map_check
from sharur.predicates import provenance
from sharur.predicates.mappings import kegg_map
from sharur.predicates.provenance import ReconcileError, kegg_map_changes, latest_stamps, map_status, reconcile
from sharur.predicates_v2.persistence import generate_and_persist_v2
from sharur.storage.duckdb_store import DuckDBStore

HEADER = "# ko\tpredicates\tevidence\n"


def _map(path, rows):
    path.write_text(HEADER + "".join(f"{ko}\t{','.join(p)}\t{';'.join(f'{x}=ec:1' for x in p)}\n"
                                     for ko, p in rows.items()))
    return path


def _store() -> DuckDBStore:
    store = DuckDBStore()
    store.conn.execute("""
        INSERT INTO bins (bin_id, completeness, contamination, taxonomy) VALUES ('bin1', 90.0, 2.0, 'd__Bacteria');
        INSERT INTO contigs (contig_id, bin_id, length, gc_content) VALUES ('contig1', 'bin1', 10000, 0.5);
        INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, gene_index,
                              sequence_length, gc_content)
        VALUES ('p1', 'contig1', 'bin1', 1, 300, '+', 1, 100, 0.5),
               ('p2', 'contig1', 'bin1', 301, 900, '+', 2, 200, 0.5),
               ('p3', 'contig1', 'bin1', 901, 1200, '+', 3, 100, 0.5);
        INSERT INTO annotations (annotation_id, protein_id, source, accession, name, description, evalue, score)
        VALUES (1, 'p1', 'kegg', 'K05986', 'nuclease S1', 'nuclease S1', 1e-30, 100.0),
               (2, 'p2', 'kegg', 'K00001', 'adh', 'alcohol dehydrogenase', 1e-30, 100.0),
               (3, 'p3', 'pfam', 'PF00005', 'ABC_tran', 'ABC transporter', 1e-30, 100.0);
    """)
    return store


@pytest.fixture()
def maps(tmp_path, monkeypatch):
    old = _map(tmp_path / "old.tsv", {"K05986": ["cas_nuclease", "endonuclease", "nuclease"],
                                      "K00001": ["dehydrogenase", "oxidoreductase"]})
    new = _map(tmp_path / "new.tsv", {"K05986": ["endonuclease", "nuclease"],
                                      "K00001": ["dehydrogenase", "oxidoreductase"]})
    monkeypatch.setattr(kegg_map, "KEGG_MAPPING_FILE", old)
    return old, new, monkeypatch


def test_kegg_map_changes_lists_only_changed_kos(maps):
    old, new, _ = maps
    assert kegg_map_changes(old, new) == ["K05986"]


def test_reconcile_regenerates_only_carriers_and_marks_current(maps):
    old, new, monkeypatch = maps
    store = _store()
    generate_and_persist_v2(store, chunk_size=1, return_states=False)  # stamped with the old map
    monkeypatch.setattr(kegg_map, "KEGG_MAPPING_FILE", new)  # the rejudged map is installed
    assert (map_status(store).state, map_status(store).changed) == ("stale", ("kegg_map_sha256",))
    check = _predicate_map_check(store)
    assert "reconcile-predicates" in check.remediation

    seen = []
    result = reconcile(store, old, lambda ids: seen.append(list(ids)))
    assert seen == [["p1"]]
    assert (result["state"], result["accessions"], result["proteins"]) == ("reconciled", ["K05986"], 1)
    assert map_status(store).state == "current"
    base = latest_stamps(store)["full"]
    assert base["scope"] == "reconciled"
    assert json.loads(base["reconciled_accessions"]) == ["K05986"]
    assert _predicate_map_check(store).state == CapabilityState.available


def test_reconcile_refuses_when_another_map_changed(maps):
    old, new, monkeypatch = maps
    store = _store()
    generate_and_persist_v2(store, chunk_size=1, return_states=False)
    monkeypatch.setattr(kegg_map, "KEGG_MAPPING_FILE", new)
    store.execute("UPDATE predicate_provenance SET pfam_map_sha256 = 'older-map'")
    with pytest.raises(ReconcileError, match="pfam_map_sha256"):
        reconcile(store, old, lambda ids: None)
    assert map_status(store).state == "stale"


def test_reconcile_refuses_a_map_that_is_not_the_stamped_one(maps, tmp_path):
    old, new, monkeypatch = maps
    store = _store()
    generate_and_persist_v2(store, chunk_size=1, return_states=False)
    monkeypatch.setattr(kegg_map, "KEGG_MAPPING_FILE", new)
    wrong = _map(tmp_path / "wrong.tsv", {"K05986": ["nuclease"]})
    with pytest.raises(ReconcileError, match="not the KEGG map"):
        reconcile(store, wrong, lambda ids: None)


def test_reconcile_is_a_no_op_when_current(maps):
    old, _, _ = maps
    store = _store()
    generate_and_persist_v2(store, chunk_size=1, return_states=False)
    result = reconcile(store, old, lambda ids: pytest.fail("nothing to regenerate"))
    assert result["state"] == "current"


def test_reconcile_with_real_regeneration_updates_predicates(maps):
    old, new, monkeypatch = maps
    store = _store()
    generate_and_persist_v2(store, chunk_size=1, return_states=False)
    monkeypatch.setattr(kegg_map, "KEGG_MAPPING_FILE", new)
    monkeypatch.setattr(kegg_map, "KEGG_TO_PREDICATES", kegg_map._load_kegg_mapping_file(new)[0])
    monkeypatch.setattr(kegg_map, "KEGG_EVIDENCE", kegg_map._load_kegg_mapping_file(new)[1])
    result = reconcile(store, old, lambda ids: generate_and_persist_v2(
        store, protein_ids=ids, chunk_size=1, update_legacy_predicates=True, return_states=False))
    assert result["proteins"] == 1
    assert map_status(store).state == "current"
    scopes = [r[0] for r in store.execute("SELECT scope FROM predicate_provenance ORDER BY generation_id")]
    assert scopes == ["full", "subset", "reconciled"]


def test_old_provenance_tables_gain_the_reconcile_columns():
    store = DuckDBStore()
    store.execute("ALTER TABLE predicate_provenance DROP COLUMN reconciled_accessions")
    assert latest_stamps(store) == {"full": None, "subsets": []}  # reads a table without the new columns
    provenance.ensure_table(store)
    cols = {r[0] for r in store.execute(
        "SELECT column_name FROM information_schema.columns WHERE table_name = 'predicate_provenance'")}
    assert {"reconciled_from", "reconciled_accessions"} <= cols
