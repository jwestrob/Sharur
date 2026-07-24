"""Tests for contig-native navigation and exhaustive-reading packets."""

import json

import pytest

from sharur.operators.contigs import (
    GENOME_PACKET_VERSION,
    PACKET_VERSION,
    get_contig,
    get_contig_packet,
    get_genome_packet,
    list_contigs,
)


def test_list_contigs_uses_stable_keyset_pagination(store):
    first = list_contigs(store, "bin_001", limit=1)
    second = list_contigs(store, "bin_001", limit=1, cursor=first.ref)

    assert first.meta.total_rows == 2
    assert first.meta.truncated is True
    assert first.records[0]["contig_id"] == "contig_001"
    assert second.records[0]["contig_id"] == "contig_002"
    assert second.ref is None
    assert second.meta.truncated is False


def test_get_contig_reports_distinct_proteins_and_domain_hits(store):
    result = get_contig(store, "bin_001", "contig_001")
    pfam = next(
        stat for stat in result.raw["annotation_stats"] if stat["source"] == "pfam"
    )

    assert result.raw["protein_count"] == 5
    assert pfam == {"source": "pfam", "domain_hits": 3, "proteins": 3}


def test_contig_packet_is_complete_paginated_and_sequence_free(store):
    store.execute(
        """
        UPDATE proteins
        SET sequence = 'SENSITIVE_SEQUENCE_PAYLOAD'
        WHERE protein_id = 'prot_001'
        """
    )
    cursor = None
    protein_ids = []
    packet_count = 0
    while True:
        result = get_contig_packet(
            store,
            "bin_001",
            "contig_001",
            cursor=cursor,
            limit=2,
        )
        packet_count += 1
        protein_ids.extend(
            protein["protein_id"] for protein in result.raw["proteins"]
        )
        assert "sequence" not in json.dumps(result.raw).lower()
        if result.raw["complete"]:
            break
        cursor = result.ref

    assert packet_count == 3
    assert protein_ids == [
        "prot_001",
        "prot_002",
        "prot_003",
        "prot_004",
        "prot_005",
    ]
    assert len(protein_ids) == len(set(protein_ids))
    assert result.raw["packet_version"] == PACKET_VERSION


def test_contig_packet_separates_observed_and_caller_named_evidence(store):
    store.execute(
        """
        INSERT INTO defense_systems(
            system_id, genome_id, contig_id, system_type, system_subtype,
            genes_count, protein_ids
        )
        VALUES (
            'call_001', 'bin_001', 'contig_001', 'CallerType',
            'CallerSubtype', 1, 'prot_002'
        )
        """
    )
    store.execute(
        """
        INSERT INTO system_proteins(
            system_id, protein_id, system_source, position, profile_name, score
        )
        VALUES (
            'call_001', 'prot_002', 'caller_projection', 1,
            'component_profile', 101
        )
        """
    )
    store.execute(
        """
        INSERT INTO annotations(
            annotation_id, protein_id, source, accession, name, evalue, score
        )
        VALUES (
            100, 'prot_002', 'caller_projection', 'PROJECTED',
            'Projected name', 1e-30, 101
        )
        """
    )
    store.execute(
        """
        INSERT INTO loci(
            locus_id, locus_type, contig_id, start, end_coord, confidence
        )
        VALUES ('locus_001', 'CallerLocus', 'contig_001', 3000, 4200, 0.9)
        """
    )
    store.execute(
        """
        INSERT INTO locus_proteins(locus_id, protein_id, position)
        VALUES ('locus_001', 'prot_002', 1)
        """
    )

    result = get_contig_packet(
        store,
        "bin_001",
        "contig_001",
        limit=5,
    )
    protein = next(
        item for item in result.raw["proteins"] if item["protein_id"] == "prot_002"
    )

    assert {row["source"] for row in protein["observed_annotations"]} == {
        "kegg",
        "pfam",
    }
    assert protein["named_calls"][0]["call_id"] == "call_001"
    assert protein["named_calls"][0]["call_type"] == "CallerType"
    assert protein["named_calls"][0]["evidence_level"] == "caller_named"
    assert protein["loci"][0]["locus_id"] == "locus_001"


def test_packet_cursor_is_scoped_to_exact_genome_and_contig(store):
    first = get_contig_packet(
        store,
        "bin_001",
        "contig_001",
        limit=1,
    )

    with pytest.raises(ValueError, match="contig_id"):
        get_contig_packet(
            store,
            "bin_001",
            "contig_002",
            cursor=first.ref,
            limit=1,
        )


def test_genome_packet_combines_contigs_and_never_crosses_bins(store):
    store.execute(
        """
        UPDATE proteins
        SET sequence = 'SENSITIVE_SEQUENCE_PAYLOAD'
        WHERE protein_id = 'prot_001'
        """
    )

    result = get_genome_packet(
        store,
        "bin_001",
        max_contigs=10,
        max_proteins=10,
        max_model_payload_bytes=100_000,
    )

    payload = result.raw["model_payload"]
    assert result.raw["packet_version"] == GENOME_PACKET_VERSION
    assert result.raw["complete"] is True
    assert payload["bin_id"] == "bin_001"
    assert [contig["contig_id"] for contig in payload["contigs"]] == [
        "contig_001",
        "contig_002",
    ]
    assert {contig["bin_id"] for contig in payload["contigs"]} == {"bin_001"}
    protein_ids = [
        protein["protein_id"]
        for contig in payload["contigs"]
        for protein in contig["proteins"]
    ]
    assert protein_ids == [
        "prot_001",
        "prot_002",
        "prot_003",
        "prot_004",
        "prot_005",
    ]
    assert "sequence" not in json.dumps(result.raw).lower()
    assert result.raw["model_payload_bytes"] == len(
        json.dumps(
            payload,
            sort_keys=True,
            separators=(",", ":"),
            default=str,
        ).encode()
    )


def test_genome_packet_splits_only_oversized_contig_with_contiguous_receipts(store):
    cursor = None
    receipts = []
    protein_ids = []
    frame_count = 0
    while True:
        result = get_genome_packet(
            store,
            "bin_001",
            cursor=cursor,
            max_contigs=10,
            max_proteins=2,
            max_model_payload_bytes=100_000,
        )
        frame_count += 1
        receipts.extend(result.raw["coverage_receipts"])
        protein_ids.extend(
            protein["protein_id"]
            for contig in result.raw["model_payload"]["contigs"]
            for protein in contig["proteins"]
        )
        if result.raw["complete"]:
            break
        cursor = result.raw["next_cursor"]

    assert frame_count == 3
    assert protein_ids == [
        "prot_001",
        "prot_002",
        "prot_003",
        "prot_004",
        "prot_005",
    ]
    contig_001 = [
        receipt for receipt in receipts if receipt["contig_id"] == "contig_001"
    ]
    assert [
        (
            receipt["protein_offset_start"],
            receipt["protein_offset_end"],
            receipt["segment_index"],
            receipt["complete"],
        )
        for receipt in contig_001
    ] == [
        (0, 2, 0, False),
        (2, 4, 1, False),
        (4, 5, 2, True),
    ]
    assert receipts[-1] == {
        "bin_id": "bin_001",
        "contig_id": "contig_002",
        "total_protein_count": 0,
        "protein_offset_start": 0,
        "protein_offset_end": 0,
        "segment_index": 0,
        "complete": True,
    }


def test_genome_packet_cursor_is_bound_to_bin_and_packing_contract(store):
    first = get_genome_packet(
        store,
        "bin_001",
        max_contigs=10,
        max_proteins=1,
        max_model_payload_bytes=100_000,
    )

    with pytest.raises(ValueError, match="genome_id"):
        get_genome_packet(
            store,
            "bin_002",
            cursor=first.ref,
            max_contigs=10,
            max_proteins=1,
            max_model_payload_bytes=100_000,
        )
    with pytest.raises(ValueError, match="packing_contract_hash"):
        get_genome_packet(
            store,
            "bin_001",
            cursor=first.ref,
            max_contigs=10,
            max_proteins=2,
            max_model_payload_bytes=100_000,
        )


def test_missing_genome_or_contig_returns_typed_empty_result(store):
    missing_genome = list_contigs(store, "missing")
    missing_contig = get_contig_packet(store, "bin_001", "missing")
    missing_genome_packet = get_genome_packet(store, "missing")

    assert missing_genome.status == "empty"
    assert missing_contig.status == "empty"
    assert missing_contig.raw["complete"] is True
    assert missing_genome_packet.status == "empty"
    assert missing_genome_packet.raw["complete"] is True
