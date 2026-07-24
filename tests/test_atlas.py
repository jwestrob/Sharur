"""Deterministic Atlas planning, dispatch, and coverage contracts."""

from __future__ import annotations

from typing import TYPE_CHECKING

import duckdb
import pytest

from sharur.atlas import (
    build_atlas_plan,
    enqueue_atlas_plan,
    load_atlas_plan,
    verify_atlas_coverage,
    write_genome_coverage_manifest,
)
from sharur.dataset_seal import build_dataset_seal, write_dataset_seal
from sharur.storage.duckdb_store import DuckDBStore


if TYPE_CHECKING:
    from pathlib import Path


def _atlas_database(root: Path) -> Path:
    database = root / "sharur.duckdb"
    store = DuckDBStore(database)
    store.execute(
        """
        INSERT INTO bins(
            bin_id, completeness, contamination, taxonomy, n_contigs, total_length
        ) VALUES
            ('genome_a', 90, 2, 'd__Bacteria;p__A', 5, 3000),
            ('genome_b', 80, 3, 'd__Bacteria;p__B', 1, 1500)
        """
    )
    store.execute(
        """
        INSERT INTO contigs(contig_id, bin_id, length)
        VALUES
            ('a_1', 'genome_a', 1000),
            ('a_2', 'genome_a', 2000),
            ('b_1', 'genome_b', 1500)
        """
    )
    store.execute(
        """
        INSERT INTO proteins(
            protein_id, contig_id, bin_id, start, end_coord, strand,
            gene_index, sequence_length
        ) VALUES
            ('a_p1', 'a_1', 'genome_a', 1, 300, '+', 1, 100),
            ('a_p2', 'a_2', 'genome_a', 1, 600, '+', 1, 200),
            ('b_p1', 'b_1', 'genome_b', 1, 450, '-', 1, 150)
        """
    )
    store.close()
    seal = build_dataset_seal(database, max_hash_bytes=0)
    write_dataset_seal(seal, root / "dataset.seal.json")
    return database


def test_atlas_plan_is_deterministic_and_uses_actual_live_counts(tmp_path):
    database = _atlas_database(tmp_path)
    first = build_atlas_plan(database, tmp_path / "first", threads=1)
    second = build_atlas_plan(database, tmp_path / "second", threads=1)
    first_manifest, first_units = load_atlas_plan(tmp_path / "first")
    _second_manifest, second_units = load_atlas_plan(tmp_path / "second")

    assert first["plan_id"] == second["plan_id"]
    assert first["units_sha256"] == second["units_sha256"]
    assert first_manifest["unit_count"] == 2
    assert first_manifest["total_contigs"] == 3
    assert first_manifest["total_proteins"] == 3
    assert first_manifest["declared_contig_count_mismatches"] == 1
    assert [unit["genome_id"] for unit in first_units] == [
        "genome_a",
        "genome_b",
    ]
    assert [unit["unit_id"] for unit in first_units] == [
        unit["unit_id"] for unit in second_units
    ]
    assert first_units[0]["work_weight"] == 2 + 32 * 2


def test_atlas_plan_fails_closed_when_unit_stream_changes(tmp_path):
    database = _atlas_database(tmp_path)
    build_atlas_plan(database, tmp_path / "plan", threads=1)
    units_path = tmp_path / "plan" / "units.jsonl"
    units_path.write_bytes(units_path.read_bytes() + b"\n")

    with pytest.raises(ValueError, match="digest"):
        load_atlas_plan(tmp_path / "plan")


def test_atlas_plan_requires_contig_navigation_index(tmp_path):
    database = _atlas_database(tmp_path)
    connection = duckdb.connect(str(database))
    connection.execute("DROP INDEX idx_contigs_bin")
    connection.close()
    seal = build_dataset_seal(database, max_hash_bytes=0)
    write_dataset_seal(seal, tmp_path / "dataset.seal.json", force=True)

    with pytest.raises(RuntimeError, match="sharur migrate"):
        build_atlas_plan(database, tmp_path / "plan", threads=1)


def test_atlas_coverage_verifies_every_owned_genome_and_exact_counts(tmp_path):
    database = _atlas_database(tmp_path)
    build_atlas_plan(database, tmp_path / "plan", threads=1)
    _manifest, units = load_atlas_plan(tmp_path / "plan")
    coverage = tmp_path / "plan" / "coverage"

    records = {
        "genome_a": [
            {
                "contig_id": "a_1",
                "protein_count": 1,
                "packet_count": 1,
                "complete": True,
            },
            {
                "contig_id": "a_2",
                "protein_count": 1,
                "packet_count": 1,
                "complete": True,
            },
        ],
        "genome_b": [
            {
                "contig_id": "b_1",
                "protein_count": 1,
                "packet_count": 1,
                "complete": True,
            }
        ],
    }
    for unit in units:
        write_genome_coverage_manifest(
            unit,
            records[unit["genome_id"]],
            coverage / f"{unit['unit_id']}.json",
        )

    result = verify_atlas_coverage(tmp_path / "plan")

    assert result["status"] == "complete"
    assert result["complete_units"] == 2
    assert result["observed_contigs"] == result["expected_contigs"] == 3
    assert result["observed_proteins"] == result["expected_proteins"] == 3


def test_atlas_coverage_reports_missing_and_incomplete_units(tmp_path):
    database = _atlas_database(tmp_path)
    build_atlas_plan(database, tmp_path / "plan", threads=1)
    _manifest, units = load_atlas_plan(tmp_path / "plan")
    first = units[0]
    write_genome_coverage_manifest(
        first,
        [
            {
                "contig_id": "a_1",
                "protein_count": 1,
                "packet_count": 1,
                "complete": True,
            }
        ],
        tmp_path / "plan" / "coverage" / f"{first['unit_id']}.json",
    )

    result = verify_atlas_coverage(tmp_path / "plan")

    assert result["status"] == "incomplete"
    assert first["unit_id"] in result["invalid_units"]
    assert units[1]["unit_id"] in result["missing_units"]


class _FakeOps:
    def __init__(self):
        self.campaigns = []
        self.tasks = []

    def create_campaign(self, name, **kwargs):
        self.campaigns.append((name, kwargs))
        return "campaign-001"

    def create_task(self, task_type, description, **kwargs):
        self.tasks.append((task_type, description, kwargs))
        return f"task-{len(self.tasks):03d}"


def test_atlas_enqueue_creates_one_idempotent_genome_task(tmp_path):
    database = _atlas_database(tmp_path)
    plan = build_atlas_plan(database, tmp_path / "plan", threads=1)
    ops = _FakeOps()

    result = enqueue_atlas_plan(
        tmp_path / "plan",
        query_url="http://query:8812",
        ops=ops,
    )

    assert result["campaign_id"] == "campaign-001"
    assert result["task_count"] == plan["unit_count"] == 2
    assert len(ops.tasks) == 2
    assert all(task[0] == "atlas_genome_read" for task in ops.tasks)
    assert all(
        task[2]["required_capabilities"]
        == ["atlas_reader", "profile:atlas_scan"]
        for task in ops.tasks
    )
    assert all(
        task[2]["params"]["review_output_contract"]["schema_version"]
        == "atlas-review-output/1.0"
        for task in ops.tasks
    )
    assert all(
        task[2]["params"]["execution_profile"] == "atlas_scan"
        for task in ops.tasks
    )
    assert all(
        task[2]["idempotency_key"].startswith(
            f"atlas-task:{plan['plan_id']}:"
        )
        for task in ops.tasks
    )
