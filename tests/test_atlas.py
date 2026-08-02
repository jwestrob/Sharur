"""Deterministic Atlas planning, dispatch, and coverage contracts."""

from __future__ import annotations

from typing import TYPE_CHECKING

import duckdb
import pytest

from sharur.atlas import (
    build_atlas_packet_census,
    build_atlas_plan,
    build_genome_coverage_manifest,
    enqueue_atlas_plan,
    load_atlas_plan,
    verify_atlas_coverage,
    verify_atlas_packet_census,
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


def _coverage_frames(unit, contigs):
    return [
        {
            "frame_id": f"frame-{unit['unit_id']}-0",
            "frame_index": 0,
            "bin_id": unit["genome_id"],
            "packet_packing_contract_hash": unit[
                "packet_packing_contract_hash"
            ],
            "model_payload_sha256": "a" * 64,
            "model_payload_bytes": 1_024,
            "contig_ids": [record["contig_id"] for record in contigs],
            "segments": [
                {
                    "bin_id": unit["genome_id"],
                    "contig_id": record["contig_id"],
                    "total_protein_count": record["protein_count"],
                    "protein_offset_start": 0,
                    "protein_offset_end": record["protein_count"],
                    "segment_index": 0,
                    "complete": True,
                }
                for record in contigs
            ],
            "protein_count": sum(record["protein_count"] for record in contigs),
            "target_exceeded": False,
        }
    ]


def test_atlas_plan_is_deterministic_and_uses_actual_live_counts(tmp_path):
    database = _atlas_database(tmp_path)
    first = build_atlas_plan(
        database,
        tmp_path / "first",
        packet_model_payload_bytes=524_288,
        threads=1,
    )
    second = build_atlas_plan(
        database,
        tmp_path / "second",
        packet_model_payload_bytes=524_288,
        threads=1,
    )
    first_manifest, first_units = load_atlas_plan(tmp_path / "first")
    _second_manifest, second_units = load_atlas_plan(tmp_path / "second")

    assert first["plan_id"] == second["plan_id"]
    assert first["units_sha256"] == second["units_sha256"]
    assert first_manifest["unit_count"] == 2
    assert first_manifest["total_contigs"] == 3
    assert first_manifest["total_proteins"] == 3
    assert first_manifest["declared_contig_count_mismatches"] == 1
    calibration = first_manifest["packet_calibration"]
    assert calibration["model_calls_made"] == 0
    assert calibration["sampled_genomes"] == 2
    assert calibration["sample_genome_count_source"] == "ceil(sqrt(live bin count))"
    assert calibration["sampled_proteins"] == 3
    assert calibration["sampled_model_payload_bytes"] > 0
    assert calibration["mean_model_payload_bytes_per_protein"] > 0
    assert first_manifest["packet_limit_sources"] == {
        "contigs": "payload-proportional schema ceiling",
        "proteins": "payload-proportional schema ceiling",
        "model_payload_bytes": "explicit model-facing target",
    }
    assert (
        first_manifest["packet_protein_limit"]
        == calibration["payload_proportional_protein_limit"]
    )
    assert (
        first_manifest["packet_contig_limit"]
        == calibration["payload_proportional_contig_limit"]
    )
    assert calibration["payload_proportional_protein_limit"] == (
        first_manifest["packet_model_payload_bytes"]
        // calibration["canonical_protein_record_floor_bytes"]
    )
    assert calibration["payload_proportional_contig_limit"] == (
        first_manifest["packet_model_payload_bytes"]
        // calibration["canonical_contig_segment_floor_bytes"]
    )
    assert (
        first_manifest["model_invocation_projection"][
            "at_calibrated_mean_density"
        ]
        > 0
    )
    assert [unit["genome_id"] for unit in first_units] == [
        "genome_a",
        "genome_b",
    ]
    assert [unit["unit_id"] for unit in first_units] == [
        unit["unit_id"] for unit in second_units
    ]
    assert first_units[0]["work_weight"] == 2 + 32 * 2


def test_atlas_plan_preserves_explicit_count_overrides(tmp_path):
    database = _atlas_database(tmp_path)
    plan = build_atlas_plan(
        database,
        tmp_path / "plan",
        packet_model_payload_bytes=524_288,
        packet_contig_limit=37,
        packet_protein_limit=41,
        threads=1,
    )

    assert plan["packet_contig_limit"] == 37
    assert plan["packet_protein_limit"] == 41
    assert plan["packet_limit_sources"]["contigs"] == "explicit"
    assert plan["packet_limit_sources"]["proteins"] == "explicit"
    assert plan["packet_calibration"]["model_calls_made"] == 0


def test_atlas_plan_fails_closed_when_unit_stream_changes(tmp_path):
    database = _atlas_database(tmp_path)
    build_atlas_plan(
        database,
        tmp_path / "plan",
        packet_model_payload_bytes=524_288,
        threads=1,
    )
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
        build_atlas_plan(
            database,
            tmp_path / "plan",
            packet_model_payload_bytes=524_288,
            threads=1,
        )


def test_atlas_packet_census_is_exact_bin_scoped_and_zero_model_call(tmp_path):
    database = _atlas_database(tmp_path)
    build_atlas_plan(
        database,
        tmp_path / "plan",
        packet_model_payload_bytes=524_288,
        packet_protein_limit=1,
        threads=1,
    )
    progress = []

    first = build_atlas_packet_census(
        tmp_path / "plan",
        workers=2,
        threads=2,
        progress_callback=progress.append,
    )
    second = build_atlas_packet_census(
        tmp_path / "plan",
        workers=2,
        threads=2,
    )
    verified = verify_atlas_packet_census(tmp_path / "plan", deep=True)

    assert first["status"] == second["status"] == "complete"
    assert first["census_sha256"] == second["census_sha256"]
    assert first["model_calls_made"] == 0
    assert first["census_resources"]["workers"] == 2
    assert progress[-1]["completed_units"] == 2
    assert progress[-1]["cumulative_frame_count"] == 3
    assert first["model_invocation_count"] == 3
    assert first["maximum_bins_per_frame"] == 1
    assert first["mixed_bin_frames"] == 0
    assert first["query_result_budget_exceeded_frames"] == 0
    assert first["observed_contigs"] == first["expected_contigs"] == 3
    assert first["observed_proteins"] == first["expected_proteins"] == 3
    assert first["frame_protein_distribution"]["maximum"] == 1
    assert verified["verification"] == "passed"
    assert verified["deep_verified"] is True


def test_atlas_packet_census_blocks_target_exceeding_singleton(tmp_path):
    database = _atlas_database(tmp_path)
    store = DuckDBStore(database)
    store.execute(
        """
        INSERT INTO annotations(
            annotation_id, protein_id, source, accession, name, description,
            evalue, score
        ) VALUES (1, 'a_p1', 'pfam', ?, 'large_record', '', 1e-20, 80)
        """,
        # The oversize must sit in a field the model payload actually carries.
        # `description` is emitted by no source in practice (empty in 4.79M of
        # 4.88M rows) and the compact encoding drops it, so planting the bulk
        # there would exercise nothing.
        ["PF-LARGE-" + "9" * 5_000],  # non-residue filler: this test is about size, not sequence detection
    )
    store.close()
    seal = build_dataset_seal(database, max_hash_bytes=0)
    write_dataset_seal(seal, tmp_path / "dataset.seal.json", force=True)
    build_atlas_plan(
        database,
        tmp_path / "plan",
        packet_model_payload_bytes=1_024,
        threads=1,
    )

    census = build_atlas_packet_census(tmp_path / "plan", threads=1)

    assert census["status"] == "blocked"
    assert census["target_exceeded_frames"] == 1
    assert census["model_calls_made"] == 0
    with pytest.raises(RuntimeError, match="blocks enqueue"):
        verify_atlas_packet_census(tmp_path / "plan")


def test_atlas_packet_census_blocks_query_result_overflow(tmp_path):
    database = _atlas_database(tmp_path)
    build_atlas_plan(
        database,
        tmp_path / "plan",
        packet_model_payload_bytes=524_288,
        query_result_limit_bytes=4_097,
        threads=1,
    )

    census = build_atlas_packet_census(tmp_path / "plan", threads=1)

    assert census["status"] == "blocked"
    assert census["query_result_budget_exceeded_frames"] == 2
    with pytest.raises(RuntimeError, match="query-result overflows"):
        verify_atlas_packet_census(tmp_path / "plan")


def test_atlas_enqueue_requires_completed_packet_census(tmp_path):
    database = _atlas_database(tmp_path)
    build_atlas_plan(
        database,
        tmp_path / "plan",
        packet_model_payload_bytes=524_288,
        threads=1,
    )

    with pytest.raises(ValueError, match="packet_census"):
        enqueue_atlas_plan(
            tmp_path / "plan",
            query_url="http://query:8812",
            ops=_FakeOps(),
        )


def test_atlas_coverage_verifies_every_owned_genome_and_exact_counts(tmp_path):
    database = _atlas_database(tmp_path)
    build_atlas_plan(
        database,
        tmp_path / "plan",
        packet_model_payload_bytes=524_288,
        threads=1,
    )
    _manifest, units = load_atlas_plan(tmp_path / "plan")
    coverage = tmp_path / "plan" / "coverage"

    records = {
        "genome_a": [
            {
                "contig_id": "a_1",
                "protein_count": 1,
                "segment_count": 1,
                "complete": True,
            },
            {
                "contig_id": "a_2",
                "protein_count": 1,
                "segment_count": 1,
                "complete": True,
            },
        ],
        "genome_b": [
            {
                "contig_id": "b_1",
                "protein_count": 1,
                "segment_count": 1,
                "complete": True,
            }
        ],
    }
    for unit in units:
        write_genome_coverage_manifest(
            unit,
            records[unit["genome_id"]],
            _coverage_frames(unit, records[unit["genome_id"]]),
            coverage / f"{unit['unit_id']}.json",
        )

    result = verify_atlas_coverage(tmp_path / "plan")

    assert result["status"] == "complete"
    assert result["complete_units"] == 2
    assert result["observed_contigs"] == result["expected_contigs"] == 3
    assert result["observed_proteins"] == result["expected_proteins"] == 3


def test_atlas_coverage_reports_missing_and_incomplete_units(tmp_path):
    database = _atlas_database(tmp_path)
    build_atlas_plan(
        database,
        tmp_path / "plan",
        packet_model_payload_bytes=524_288,
        threads=1,
    )
    _manifest, units = load_atlas_plan(tmp_path / "plan")
    first = units[0]
    write_genome_coverage_manifest(
        first,
        [
            {
                "contig_id": "a_1",
                "protein_count": 1,
                "segment_count": 1,
                "complete": True,
            }
        ],
        _coverage_frames(
            first,
            [
                {
                    "contig_id": "a_1",
                    "protein_count": 1,
                    "segment_count": 1,
                    "complete": True,
                }
            ],
        ),
        tmp_path / "plan" / "coverage" / f"{first['unit_id']}.json",
    )

    result = verify_atlas_coverage(tmp_path / "plan")

    assert result["status"] == "incomplete"
    assert first["unit_id"] in result["invalid_units"]
    assert units[1]["unit_id"] in result["missing_units"]


def test_atlas_coverage_rejects_cross_bin_frame_receipt(tmp_path):
    database = _atlas_database(tmp_path)
    build_atlas_plan(
        database,
        tmp_path / "plan",
        packet_model_payload_bytes=524_288,
        threads=1,
    )
    _manifest, units = load_atlas_plan(tmp_path / "plan")
    unit = units[1]
    contigs = [
        {
            "contig_id": "b_1",
            "protein_count": 1,
            "segment_count": 1,
            "complete": True,
        }
    ]
    frames = _coverage_frames(unit, contigs)
    frames[0]["bin_id"] = "genome_a"
    frames[0]["segments"][0]["bin_id"] = "genome_a"

    coverage = build_genome_coverage_manifest(unit, contigs, frames)

    assert coverage["coverage_status"] == "incomplete"
    assert any("belongs to bin genome_a" in error for error in coverage["errors"])


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
    plan = build_atlas_plan(
        database,
        tmp_path / "plan",
        packet_model_payload_bytes=524_288,
        threads=1,
    )
    census = build_atlas_packet_census(tmp_path / "plan", threads=1)
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
        task[2]["params"]["packet_census_sha256"] == census["census_sha256"]
        for task in ops.tasks
    )
    assert all(
        task[2]["idempotency_key"].startswith(
            f"atlas-task:{plan['plan_id']}:"
        )
        for task in ops.tasks
    )
