"""Tests for contig-native navigation and exhaustive-reading packets."""

import json

import pytest

from sharur.operators import contigs as contigs_module
from sharur.operators.contigs import (
    GENOME_PACKET_VERSION,
    PACKET_VERSION,
    genome_packet_packing_contract,
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
    # Proteins are positional rows under the compact encoding; resolve the
    # column by name so this stays correct if the column order changes.
    pid_col = payload["protein_columns"].index("protein_id")
    protein_ids = [
        protein[pid_col]
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
        pid_col = result.raw["model_payload"]["protein_columns"].index("protein_id")
        protein_ids.extend(
            protein[pid_col]
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


def test_genome_packet_count_guards_scale_with_payload_budget():
    automatic = genome_packet_packing_contract(
        max_model_payload_bytes=1_024,
    )
    floors = automatic["serialized_record_floor_bytes"]
    assert automatic["max_contigs"] == 1_024 // floors["contig"]
    assert automatic["max_proteins"] == 1_024 // floors["protein"]

    contract = genome_packet_packing_contract(
        max_contigs=4,
        max_proteins=6,
        max_model_payload_bytes=1_024,
    )

    assert contract["serialized_record_floor_bytes"]["contig"] > 0
    assert contract["serialized_record_floor_bytes"]["protein"] > 0
    with pytest.raises(ValueError, match="payload-proportional"):
        genome_packet_packing_contract(
            max_contigs=100,
            max_proteins=6,
            max_model_payload_bytes=1_024,
        )
    with pytest.raises(ValueError, match="payload-proportional"):
        genome_packet_packing_contract(
            max_contigs=4,
            max_proteins=100,
            max_model_payload_bytes=1_024,
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


class TestIncrementalFrameByteAccounting:
    """Frame size is a running sum, not a full re-serialisation per contig.

    The packing loop used to re-serialise the entire accumulated frame once per
    contig considered, which is O(n^2) in contigs per frame and GIL-bound (the
    C JSON encoder holds the GIL). That dominated the offline packet census.

    The replacement relies on canonical JSON being compositional. If that
    identity ever breaks, frames pack differently, the packing-contract hash
    changes, and every existing plan and census is silently invalidated — so it
    is asserted directly here rather than only end-to-end.
    """

    @staticmethod
    def _payload(contigs):
        return {
            "schema_version": "atlas-model-packet/1.0",
            "packet_version": "genome-packet/1.1",
            "dataset_id": "d" * 64,
            "bin_id": "genome_a",
            "frame_index": 3,
            "genome": {"bin_id": "genome_a", "taxonomy": "d__Bacteria;p__Chloroflexota"},
            "contigs": contigs,
        }

    @staticmethod
    def _contig(index, n_proteins):
        return {
            "bin_id": "genome_a",
            "contig_id": f"contig_{index}",
            "length": 40000 + index,
            "gc_content": 0.5512,
            "is_circular": False,
            "taxonomy": None,
            "total_protein_count": n_proteins,
            "protein_offset_start": 0,
            "protein_offset_end": n_proteins,
            "segment_index": 0,
            "complete": True,
            "proteins": [
                {
                    "protein_id": f"contig_{index}_p{j}",
                    "gene_index": j,
                    "start": j * 900,
                    "end": j * 900 + 720,
                    "strand": "+" if j % 2 else "-",
                    "length_aa": 240,
                    "gc_content": 0.53,
                    "observed_annotations": [
                        {"accession": "PF00466", "name": "Ribosomal_L10", "evalue": 1e-30}
                    ],
                    "predicates": ["ribosomal"],
                    "named_calls": [],
                    "loci": [],
                }
                for j in range(n_proteins)
            ],
        }

    def _incremental(self, contigs):
        envelope = len(contigs_module._canonical_bytes(self._payload([])))
        return (
            envelope
            + sum(len(contigs_module._canonical_bytes(c)) for c in contigs)
            + (len(contigs) - 1 if contigs else 0)
        )

    @pytest.mark.parametrize("n_contigs", [1, 2, 3, 17, 64, 250])
    def test_running_sum_equals_full_serialisation(self, n_contigs):
        contigs = [self._contig(i, (i % 7) + 1) for i in range(n_contigs)]
        full = len(contigs_module._canonical_bytes(self._payload(contigs)))
        assert self._incremental(contigs) == full

    def test_holds_for_an_empty_frame(self):
        empty = len(contigs_module._canonical_bytes(self._payload([])))
        assert self._incremental([]) == empty

    def test_holds_with_unicode_and_nulls(self):
        """default=str and non-ASCII must not break the composition."""
        contigs = [self._contig(0, 2), self._contig(1, 2)]
        contigs[0]["taxonomy"] = "d__Bacteria;p__Chloroflexota;s__Ca. Dormibacter µ"
        contigs[1]["gc_content"] = None
        full = len(contigs_module._canonical_bytes(self._payload(contigs)))
        assert self._incremental(contigs) == full