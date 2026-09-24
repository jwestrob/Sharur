"""Protein cards and predicate explanations."""

import json

import pytest

from sharur.operators.cards import card, card_markdown, mapping_evidence, why, why_markdown
from sharur.predicates_v2.persistence import generate_and_persist_v2
from sharur.storage.duckdb_store import DuckDBStore


SEQ = "MKTAYIAKQRQISFVKSHFSRQ" * 10


@pytest.fixture
def store():
    store = DuckDBStore()
    store.conn.execute(f"""
        INSERT INTO bins (bin_id, completeness, contamination, taxonomy)
        VALUES ('bin1', 91.5, 1.2, 'd__Archaea;p__Micrarchaeota');
        INSERT INTO contigs (contig_id, bin_id, length, gc_content) VALUES ('c1', 'bin1', 20000, 0.4);
        INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, gene_index,
                              sequence, sequence_length, gc_content)
        VALUES ('p0', 'c1', 'bin1', 1, 600, '+', 0, '{SEQ}', 220, 0.4),
               ('p1', 'c1', 'bin1', 700, 1400, '+', 1, '{SEQ}', 220, 0.4),
               ('p2', 'c1', 'bin1', 1500, 2100, '-', 2, '{SEQ}', 220, 0.4);
        INSERT INTO annotations (annotation_id, protein_id, source, accession, name, description, evalue, score)
        VALUES (1, 'p1', 'pfam', 'PF00005', 'ABC_tran', 'ABC transporter', 1e-40, 150.0),
               (2, 'p0', 'pfam', 'PF00106', 'adh_short', 'short chain dehydrogenase', 1e-30, 120.0);
    """)
    generate_and_persist_v2(store, chunk_size=1, return_states=False)
    return store


def test_direct_mapping_cites_map_evidence(store):
    result = why(store, "p1", "abc_transporter")
    assert result["present"]
    (path,) = [p for p in result["paths"] if p["source_db"] == "pfam"]
    assert (path["accession"], path["annotation"]) == ("PF00005", "ABC_tran")
    assert path["mapping"]["kind"] == "direct"
    assert path["mapping"]["evidence"].startswith(("text:", "go:", "swissprot:", "enzyme:"))


def test_expanded_predicates_show_the_mapped_child_and_chain(store):
    result = why(store, "p1", "nucleotide_binding")
    mapping = next(p["mapping"] for p in result["paths"] if p["source_db"] == "pfam")
    assert mapping["kind"] == "expanded"
    assert mapping["chain"][0] == mapping["via"]
    assert mapping["chain"][-1] == "nucleotide_binding"
    assert "via" in why_markdown(result)


def test_absent_predicates_are_reported(store):
    result = why(store, "p1", "hydrogenase")
    assert not result["present"]
    assert result["paths"] == []
    assert "Not present" in why_markdown(result)


def test_atoms_from_an_earlier_map_are_marked_unsupported():
    assert mapping_evidence("pfam", "PF00005", "mercury_resistance")["kind"] == "unsupported"
    assert mapping_evidence("pfam", "PF99999", "kinase")["kind"] == "unsupported"
    assert mapping_evidence("_property", "_computed", "small")["kind"] == "computed"
    assert mapping_evidence("defensefinder_system", "RM_Type_I", "defense_system")["kind"] == "system_call"


def test_card_summarizes_context_evidence_and_provenance(store):
    c = card(store, "p1", window=1)
    assert c["found"]
    assert c["genome"]["taxonomy"].startswith("d__Archaea")
    assert c["annotations"]["pfam"][0]["name"] == "ABC_tran"
    roles_and_activities = {**c["predicates"]["roles"], **c["predicates"]["activities"]}
    assert "pfam:PF00005" in roles_and_activities["abc_transporter"]
    assert [n["offset"] for n in c["neighborhood"]] == [-1, 0, 1]
    assert c["neighborhood"][0]["top_annotation"] == "pfam:adh_short"
    assert c["map_status"]["state"] == "current"
    markdown = card_markdown(c)
    assert "## Neighborhood" in markdown
    assert ">> +0" in markdown


def test_cards_and_explanations_carry_no_sequence(store):
    assert SEQ[:30] not in json.dumps(card(store, "p1"), default=str)
    assert SEQ[:30] not in json.dumps(why(store, "p1", "abc_transporter"), default=str)


def test_missing_protein(store):
    assert card(store, "nope") == {"protein_id": "nope", "found": False}
    assert "not found" in card_markdown(card(store, "nope"))
