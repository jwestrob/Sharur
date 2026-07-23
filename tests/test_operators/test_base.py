"""Tests for the shared operator result and trace contract."""

import json

from sharur.operators.base import OperatorContext
from sharur.storage.duckdb_store import DuckDBStore


def _result_for_store(store):
    with OperatorContext("same_operator", {"limit": 5}, store=store) as context:
        store.execute("SELECT COUNT(*) FROM proteins")
        return context.make_result(data="ok", rows=1)


def test_trace_uses_live_schema_version(tmp_path):
    dataset_dir = tmp_path / "dataset"
    dataset_dir.mkdir()
    store = DuckDBStore(dataset_dir / "sharur.duckdb")

    result = _result_for_store(store)

    assert result.trace.schema_version == "6"
    assert "sharur.duckdb@size=" in result.trace.dataset_version


def test_trace_hash_is_scoped_to_dataset_instance(tmp_path):
    left_dir = tmp_path / "left"
    right_dir = tmp_path / "right"
    left_dir.mkdir()
    right_dir.mkdir()
    left = DuckDBStore(left_dir / "sharur.duckdb")
    right = DuckDBStore(right_dir / "sharur.duckdb")

    left_result = _result_for_store(left)
    right_result = _result_for_store(right)

    assert left_result.trace.trace_hash != right_result.trace.trace_hash


def test_non_dataset_operator_does_not_claim_fake_versions():
    with OperatorContext("standalone", {}) as context:
        result = context.make_result(data="ok", rows=1)

    assert result.trace.dataset_version == "not_applicable"
    assert result.trace.schema_version == "not_applicable"


def test_result_exposes_public_programmatic_payload():
    with OperatorContext("records", {}) as context:
        result = context.make_result(
            data="formatted",
            rows=1,
            raw=[{"protein_id": "p1"}],
        )

    assert result.raw == [{"protein_id": "p1"}]
    assert result.records == [{"protein_id": "p1"}]
    assert result.to_dataframe().iloc[0]["protein_id"] == "p1"
    payload = json.loads(result.to_json())
    assert payload["status"] == "ok"
    assert payload["raw"][0]["protein_id"] == "p1"


def test_zero_row_result_has_explicit_empty_status():
    with OperatorContext("empty", {}) as context:
        result = context.make_result(data="No matches", rows=0, raw=[])

    assert result.status == "empty"
