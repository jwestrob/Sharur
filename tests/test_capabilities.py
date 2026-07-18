"""Typed dataset capability/preflight checks."""

from types import SimpleNamespace

import h5py
import numpy as np

from sharur import capabilities
from sharur.capabilities import CapabilityState, build_capability_brief
from sharur.storage.duckdb_store import DuckDBStore


def _build_ready_database(path):
    store = DuckDBStore(path)
    store.execute("INSERT INTO bins(bin_id, n_contigs, total_length) VALUES ('bin1', 1, 300)")
    store.execute("INSERT INTO contigs(contig_id, bin_id, length) VALUES ('contig1', 'bin1', 300)")
    store.execute(
        """
        INSERT INTO proteins(
            protein_id, contig_id, bin_id, start, end_coord, strand, sequence,
            sequence_length
        ) VALUES ('p1', 'contig1', 'bin1', 1, 300, '+', 'MA', 2)
        """
    )
    store.execute(
        """
        INSERT INTO annotations(
            annotation_id, protein_id, source, accession, name
        ) VALUES (1, 'p1', 'test_domains', 'T0001', 'Observed test domain')
        """
    )
    store.execute(
        """
        INSERT INTO semantic_state(
            protein_id, activities, roles, architecture, localization, topology,
            size_class, quality_flags, composite_predicates, unresolved_count
        ) VALUES ('p1', [], [], [], [], '{}', 'standard', [], [], 0)
        """
    )
    store.execute(
        """
        INSERT INTO semantic_terms(
            protein_id, term_id, term_kind, facet, relation, source_db, source_accession
        ) VALUES ('p1', 'standard', 'state', 'size', 'observed', '', '')
        """
    )
    store.execute("INSERT INTO protein_predicates(protein_id, predicates) VALUES ('p1', [])")
    return store


def test_preflight_uses_only_typed_states_and_reports_live_sources(tmp_path):
    db_path = tmp_path / "sharur.duckdb"
    store = _build_ready_database(db_path)
    store.conn.close()

    brief = build_capability_brief(db_path)

    assert brief.overall_state == CapabilityState.available
    assert {capability.state for capability in brief.capabilities} <= set(CapabilityState)
    annotation_capability = brief.get("annotation_sources")
    assert annotation_capability.evidence["claim_level"] == "mixed"
    sources = annotation_capability.evidence["sources"]
    assert sources == [{"source": "test_domains", "rows": 1, "proteins": 1}]
    assert brief.get("similarity_index").state == CapabilityState.unavailable


def test_preflight_marks_partial_semantic_materialization_stale(tmp_path):
    db_path = tmp_path / "sharur.duckdb"
    store = _build_ready_database(db_path)
    store.execute(
        """
        INSERT INTO proteins(
            protein_id, contig_id, bin_id, start, end_coord, strand, sequence,
            sequence_length
        ) VALUES ('p2', 'contig1', 'bin1', 301, 600, '+', 'MV', 2)
        """
    )
    store.conn.close()

    brief = build_capability_brief(db_path)

    assert brief.get("semantic_v2").state == CapabilityState.stale
    assert brief.get("predicate_compatibility").state == CapabilityState.stale
    assert brief.overall_state == CapabilityState.stale


def test_preflight_marks_embedding_row_count_drift_stale(tmp_path):
    db_path = tmp_path / "sharur.duckdb"
    store = _build_ready_database(db_path)
    store.close()
    embeddings_dir = tmp_path / "embeddings"
    embeddings_dir.mkdir()
    with h5py.File(embeddings_dir / "protein_embeddings.h5", "w") as handle:
        handle.create_dataset("protein_ids", data=np.asarray(["p1", "extra"], dtype="S"))
        handle.create_dataset(
            "embeddings",
            data=np.asarray([[1.0, 0.0], [0.0, 1.0]], dtype=np.float32),
        )

    brief = build_capability_brief(db_path)

    embeddings = brief.get("embeddings")
    assert embeddings.state == CapabilityState.stale
    assert embeddings.evidence["count"] == 2
    assert embeddings.evidence["sequence_protein_count"] == 1
    assert embeddings.evidence["row_count_ratio"] == 2.0
    assert brief.get("similarity_index").remediation == (
        "Regenerate current embeddings before rebuilding the vector index."
    )


def test_preflight_marks_empty_protein_rows_as_dataset_stale(tmp_path):
    db_path = tmp_path / "sharur.duckdb"
    store = _build_ready_database(db_path)
    store.execute(
        """
        INSERT INTO proteins(
            protein_id, contig_id, bin_id, start, end_coord, strand
        ) VALUES ('empty', 'contig1', 'bin1', 301, 600, '+')
        """
    )
    store.close()

    brief = build_capability_brief(db_path)

    dataset = brief.get("dataset")
    assert dataset.state == CapabilityState.stale
    assert dataset.evidence["protein_count"] == 2
    assert dataset.evidence["sequence_protein_count"] == 1
    assert dataset.evidence["invalid_sequence_rows"] == 1
    assert brief.overall_state == CapabilityState.stale


def test_preflight_discovers_structured_caller_tables_by_live_schema(tmp_path):
    db_path = tmp_path / "sharur.duckdb"
    store = _build_ready_database(db_path)
    store.execute(
        """
        CREATE TABLE custom_caller_output (
            system_id VARCHAR,
            system_type VARCHAR
        )
        """
    )
    store.execute("INSERT INTO custom_caller_output VALUES ('call-1', 'caller-emitted-name')")
    store.conn.close()

    brief = build_capability_brief(db_path)
    resources = brief.get("structured_callers").evidence["resources"]
    custom = next(item for item in resources if item["table"] == "custom_caller_output")

    assert custom["rows"] == 1
    assert custom["emitted_types"] == [{"value": "caller-emitted-name", "rows": 1}]


def test_preflight_missing_database_is_machine_readable(tmp_path):
    brief = build_capability_brief(tmp_path / "missing.duckdb")

    assert brief.overall_state == CapabilityState.unavailable
    payload = brief.to_dict()
    assert payload["capabilities"][0]["state"] == "unavailable"


def test_mps_probe_reports_unavailable_as_a_typed_absence(monkeypatch):
    monkeypatch.setattr(capabilities.sys, "platform", "darwin")
    monkeypatch.setattr(
        capabilities.subprocess,
        "run",
        lambda *args, **kwargs: SimpleNamespace(stdout='{"available": false, "version": "test"}'),
    )

    check = next(
        item for item in capabilities._execution_checks() if item.capability_id == "execution_mps"
    )

    assert check.state == CapabilityState.unavailable
