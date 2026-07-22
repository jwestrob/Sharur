"""Exact, run-scoped ELSA sidecar access from Sharur."""

from __future__ import annotations

import json

import duckdb
import pytest

from sharur.capabilities import CapabilityState, build_capability_brief
from sharur.operators import Sharur
from sharur.storage.duckdb_store import DuckDBStore
from sharur.synteny import (
    SyntenyDatasetMismatchError,
    SyntenyStore,
    inspect_synteny_sidecar,
)


def _build_core_database(path) -> None:
    store = DuckDBStore(path)
    store.execute(
        "INSERT INTO bins(bin_id, n_contigs, total_length) VALUES ('genome_a', 1, 900)"
    )
    store.execute(
        """
        INSERT INTO contigs(contig_id, bin_id, length)
        VALUES ('contig_a', 'genome_a', 900)
        """
    )
    store.conn.executemany(
        """
        INSERT INTO proteins(
            protein_id, contig_id, bin_id, start, end_coord, strand,
            gene_index, sequence, sequence_length
        ) VALUES (?, 'contig_a', 'genome_a', ?, ?, '+', ?, 'MA', 2)
        """,
        [
            ("p1", 1, 270, 0),
            ("p10", 301, 570, 1),
        ],
    )
    store.conn.executemany(
        "INSERT INTO protein_predicates(protein_id, predicates) VALUES (?, [])",
        [("p1",), ("p10",)],
    )
    store.close()


def _build_synteny_sidecar(path, *, dataset_id: str = "dataset-1") -> str:
    run_id = "elsa-run-a"
    connection = duckdb.connect(str(path))
    connection.execute(
        """
        CREATE TABLE elsa_schema_version(version INTEGER, installed_at TIMESTAMPTZ);
        INSERT INTO elsa_schema_version VALUES (1, current_timestamp);

        CREATE TABLE elsa_runs(
            run_id VARCHAR,
            run_label VARCHAR,
            created_at TIMESTAMPTZ,
            status VARCHAR,
            is_active BOOLEAN,
            dataset_id VARCHAR,
            elsa_version VARCHAR,
            elsa_commit VARCHAR,
            mapping_version VARCHAR,
            gene_count BIGINT,
            block_count BIGINT,
            cluster_count BIGINT,
            singleton_count BIGINT,
            anchor_pair_count BIGINT,
            locus_count BIGINT,
            member_count BIGINT,
            validation_json VARCHAR
        );
        CREATE TABLE elsa_dataset_state(singleton BOOLEAN, active_run_id VARCHAR);
        CREATE TABLE elsa_genes(
            run_id VARCHAR,
            protein_id VARCHAR,
            genome_id VARCHAR,
            contig_id VARCHAR,
            position_index BIGINT,
            start_bp BIGINT,
            end_bp BIGINT,
            strand BIGINT
        );
        CREATE TABLE elsa_blocks(
            run_id VARCHAR,
            block_id BIGINT,
            cluster_key VARCHAR,
            source_cluster_id BIGINT,
            query_genome VARCHAR,
            target_genome VARCHAR,
            query_contig VARCHAR,
            target_contig VARCHAR,
            query_start BIGINT,
            query_end BIGINT,
            target_start BIGINT,
            target_end BIGINT,
            query_start_bp BIGINT,
            query_end_bp BIGINT,
            target_start_bp BIGINT,
            target_end_bp BIGINT,
            n_anchors BIGINT,
            chain_score DOUBLE,
            orientation SMALLINT
        );
        CREATE TABLE elsa_anchor_pairs(
            run_id VARCHAR,
            block_id BIGINT,
            cluster_key VARCHAR,
            pair_order BIGINT,
            query_protein_id VARCHAR,
            target_protein_id VARCHAR,
            query_position_index BIGINT,
            target_position_index BIGINT
        );
        CREATE TABLE elsa_clusters(
            run_id VARCHAR,
            cluster_key VARCHAR,
            source_cluster_id BIGINT,
            cluster_kind VARCHAR,
            size BIGINT,
            genome_support BIGINT,
            mean_anchor_count DOUBLE,
            mean_chain_score DOUBLE,
            locus_count BIGINT,
            member_count BIGINT,
            anchor_member_count BIGINT,
            source_metadata_present BOOLEAN,
            reported_size BIGINT,
            reported_genome_support BIGINT,
            reported_mean_chain_length DOUBLE
        );
        CREATE TABLE elsa_cluster_loci(
            run_id VARCHAR,
            cluster_key VARCHAR,
            locus_key VARCHAR,
            locus_index BIGINT,
            genome_id VARCHAR,
            contig_id VARCHAR,
            start_position_index BIGINT,
            end_position_index BIGINT,
            start_bp BIGINT,
            end_bp BIGINT,
            start_protein_id VARCHAR,
            end_protein_id VARCHAR,
            n_genes BIGINT,
            block_support BIGINT
        );
        CREATE TABLE elsa_cluster_members(
            run_id VARCHAR,
            cluster_key VARCHAR,
            locus_key VARCHAR,
            protein_id VARCHAR,
            genome_id VARCHAR,
            contig_id VARCHAR,
            position_index BIGINT,
            locus_order BIGINT,
            member_role VARCHAR,
            start_bp BIGINT,
            end_bp BIGINT,
            strand BIGINT
        );
        """
    )
    connection.execute(
        """
        INSERT INTO elsa_runs VALUES (
            ?, 'test run', current_timestamp, 'ready', TRUE, ?, '0.2.0',
            'commit-a', 'coordinate-sorted-position-index-v1',
            6, 2, 2, 1, 3, 3, 6, '{"exact_anchor_resolution":true}'
        )
        """,
        [run_id, dataset_id],
    )
    connection.execute(
        "INSERT INTO elsa_dataset_state VALUES (TRUE, ?)",
        [run_id],
    )
    connection.execute(
        """
        INSERT INTO elsa_genes VALUES
            ('elsa-run-a', 'p1', 'genome_a', 'contig_a', 0, 1, 270, 1),
            ('elsa-run-a', 'p10', 'genome_a', 'contig_a', 1, 301, 570, 1),
            ('elsa-run-a', 'q2', 'genome_b', 'contig_b', 2, 601, 870, -1);

        INSERT INTO elsa_blocks VALUES
            ('elsa-run-a', 100, 'cluster:7', 7, 'genome_a', 'genome_b',
             'contig_a', 'contig_b', 0, 1, 1, 2, 1, 570, 301, 870,
             2, 1.9, -1),
            ('elsa-run-a', 101, 'cluster:8', 8, 'genome_a', 'genome_c',
             'contig_a', 'contig_c', 0, 0, 0, 0, 1, 270, 1, 270,
             1, 0.9, 1);

        INSERT INTO elsa_anchor_pairs VALUES
            ('elsa-run-a', 100, 'cluster:7', 0, 'p1', 'q2', 0, 2),
            ('elsa-run-a', 100, 'cluster:7', 1, 'p10', 'q1', 1, 1),
            ('elsa-run-a', 101, 'cluster:8', 0, 'r1', 's1', 0, 0);

        INSERT INTO elsa_clusters VALUES
            ('elsa-run-a', 'cluster:7', 7, 'cluster', 2, 3, 2.0, 1.9,
             2, 4, 3, TRUE, 2, 3, 2.0),
            ('elsa-run-a', 'cluster:8', 8, 'singleton', 1, 2, 1.0, 0.9,
             1, 2, 1, FALSE, NULL, NULL, NULL);

        INSERT INTO elsa_cluster_loci VALUES
            ('elsa-run-a', 'cluster:7', 'cluster:7:locus:0', 0,
             'genome_a', 'contig_a',
             0, 1, 1, 570, 'p1', 'p10', 2, 2),
            ('elsa-run-a', 'cluster:7', 'cluster:7:locus:1', 1,
             'genome_b', 'contig_b',
             1, 2, 301, 870, 'q1', 'q2', 2, 2),
            ('elsa-run-a', 'cluster:8', 'cluster:8:locus:0', 0,
             'genome_a', 'contig_a',
             0, 1, 1, 570, 'p1', 'p10', 2, 1);

        INSERT INTO elsa_cluster_members VALUES
            ('elsa-run-a', 'cluster:7', 'cluster:7:locus:0', 'p1', 'genome_a',
             'contig_a', 0, 0, 'anchor', 1, 270, 1),
            ('elsa-run-a', 'cluster:7', 'cluster:7:locus:0', 'p10', 'genome_a',
             'contig_a', 1, 1, 'anchor', 301, 570, 1),
            ('elsa-run-a', 'cluster:7', 'cluster:7:locus:1', 'q1', 'genome_b',
             'contig_b', 1, 0, 'anchor', 301, 570, -1),
            ('elsa-run-a', 'cluster:7', 'cluster:7:locus:1', 'q2', 'genome_b',
             'contig_b', 2, 1, 'anchor', 601, 870, -1),
            ('elsa-run-a', 'cluster:8', 'cluster:8:locus:0', 'p1', 'genome_a',
             'contig_a', 0, 0, 'context', 1, 270, 1),
            ('elsa-run-a', 'cluster:8', 'cluster:8:locus:0', 'p10', 'genome_a',
             'contig_a', 1, 1, 'anchor', 301, 570, 1);
        """
    )
    connection.execute(
        """
        CREATE VIEW current_elsa_run AS
        SELECT runs.*
        FROM elsa_runs AS runs
        JOIN elsa_dataset_state AS state
          ON state.active_run_id = runs.run_id
        WHERE state.singleton = TRUE
        """
    )
    connection.close()
    return run_id


def test_exact_membership_preserves_many_to_many_and_avoids_substrings(tmp_path):
    sidecar_path = tmp_path / "synteny.duckdb"
    run_id = _build_synteny_sidecar(sidecar_path)

    with SyntenyStore(sidecar_path) as store:
        rows, total, resolved_run = store.protein_memberships(["p1"])
        anchor_rows, anchor_total, _ = store.anchor_blocks("p1")

    assert resolved_run == run_id
    assert total == 2
    assert {row["cluster_key"] for row in rows} == {"cluster:7", "cluster:8"}
    assert {row["protein_id"] for row in rows} == {"p1"}
    assert anchor_total == 1
    assert anchor_rows[0]["partner_protein_id"] == "q2"
    assert anchor_rows[0]["orientation"] == -1


def test_facade_queries_run_scoped_cluster_and_inspect_enriches_case(tmp_path):
    db_path = tmp_path / "sharur.duckdb"
    sidecar_path = tmp_path / "synteny.duckdb"
    _build_core_database(db_path)
    run_id = _build_synteny_sidecar(sidecar_path)

    browser = Sharur(db_path, read_only=True)
    memberships = browser.synteny_for_protein("p1")
    cluster = browser.get_synteny_cluster(7, run_id=run_id, member_limit=2)
    case = browser.inspect("p1", entity_type="protein", window=0)

    assert memberships.meta.total_rows == 2
    assert memberships.trace is not None
    assert f"elsa_run={run_id}" in memberships.trace.dataset_version
    assert cluster.raw["cluster_key"] == "cluster:7"
    assert cluster.meta.truncated
    assert case.record.synteny_state == "available"
    assert len(case.record.synteny_memberships) == 2
    assert {item.evidence_level.value for item in case.record.synteny_memberships} == {
        "inferred"
    }


def test_capability_reports_schema_counts_and_dataset_identity(tmp_path):
    db_path = tmp_path / "sharur.duckdb"
    sidecar_path = tmp_path / "synteny.duckdb"
    _build_core_database(db_path)
    _build_synteny_sidecar(sidecar_path)
    (tmp_path / "dataset.seal.json").write_text(
        json.dumps({"dataset_id": "dataset-1"})
    )

    inspection = inspect_synteny_sidecar(sidecar_path)
    capability = build_capability_brief(db_path).get("elsa_synteny")

    assert inspection.state == "available"
    assert inspection.cluster_count == 2
    assert capability.state == CapabilityState.available
    assert capability.evidence["dataset_identity_check"] == "match"
    assert capability.evidence["anchor_pair_count"] == 3

    (tmp_path / "dataset.seal.json").write_text(
        json.dumps({"dataset_id": "dataset-2"})
    )
    stale = build_capability_brief(db_path).get("elsa_synteny")
    assert stale.state == CapabilityState.stale
    assert stale.evidence["dataset_identity_check"] == "mismatch"


def test_dataset_mismatch_fails_closed_with_explicit_historical_opt_in(tmp_path):
    db_path = tmp_path / "sharur.duckdb"
    sidecar_path = tmp_path / "synteny.duckdb"
    _build_core_database(db_path)
    _build_synteny_sidecar(sidecar_path, dataset_id="historical-dataset")
    (tmp_path / "dataset.seal.json").write_text(
        json.dumps({"dataset_id": "current-dataset"})
    )

    browser = Sharur(db_path, read_only=True)
    with pytest.raises(SyntenyDatasetMismatchError, match="identity is mismatch"):
        browser.synteny_for_protein("p1")
    case = browser.inspect("p1", entity_type="protein", window=0)
    assert case.record.synteny_state == "stale"
    assert case.record.synteny_memberships == []

    historical = Sharur(
        db_path,
        read_only=True,
        allow_stale_synteny=True,
    )
    memberships = historical.synteny_for_protein("p1")
    historical_case = historical.inspect("p1", entity_type="protein", window=0)
    assert memberships.meta.total_rows == 2
    assert historical_case.record.synteny_state == "available"
    assert len(historical_case.record.synteny_memberships) == 2
