"""Adversarial, migration, and scale contracts for Sharur Ops v3."""

from __future__ import annotations

import sqlite3
import time
import uuid
from concurrent.futures import ThreadPoolExecutor

import pytest

from sharur.ops.db import SQLiteConnectionPool, SQLiteServerLock, open_ops_connection
from sharur.ops.schema import OPS_SCHEMA_VERSION, ensure_ops_schema
from sharur.ops.store import (
    MAX_INLINE_JSON_BYTES,
    IdempotencyConflictError,
    LeaseFenceError,
    OpsStore,
)


def test_attempt_token_fences_stale_same_agent_process(tmp_path):
    db_path = tmp_path / "ops.db"
    coordinator = OpsStore(db_path, agent_id="coordinator")
    task_id = coordinator.create_task(
        "survey",
        "Fence replacement attempts",
        lease_seconds=10,
        max_attempts=2,
    )
    stale_process = OpsStore(db_path, agent_id="same-agent")
    first = stale_process.claim_task(task_id)

    coordinator.recover_expired_tasks(now=first["lease_expires_ts"] + 1)
    replacement_process = OpsStore(db_path, agent_id="same-agent")
    replacement = replacement_process.claim_task(task_id)

    assert replacement["lease_attempt"] == first["lease_attempt"] + 1
    with pytest.raises(LeaseFenceError, match="attempt 1"):
        stale_process.complete_task(task_id)
    completed = replacement_process.complete_task(task_id)
    assert completed["status"] == "complete"


def test_concurrent_hypothesis_evidence_updates_preserve_both_writers(tmp_path):
    db_path = tmp_path / "ops.db"
    seed = OpsStore(db_path, agent_id="seed")
    finding_ids = [
        seed.finding("observation", "fixture", f"evidence {index}")
        for index in range(2)
    ]
    hypothesis_id = seed.hypothesis("Concurrent evidence remains lossless")
    writers = [
        OpsStore(db_path, agent_id="writer-a"),
        OpsStore(db_path, agent_id="writer-b"),
    ]

    def update(args):
        writer, finding_id = args
        return writer.update_hypothesis(
            hypothesis_id,
            evidence_for=[finding_id],
        )

    with ThreadPoolExecutor(max_workers=2) as executor:
        list(executor.map(update, zip(writers, finding_ids, strict=True)))

    updated = seed.get_hypothesis(hypothesis_id)
    assert set(updated["evidence_for"]) == set(finding_ids)
    normalized = seed._conn.execute(
        """
        SELECT finding_id FROM hypothesis_findings
        WHERE hypothesis_id = ? AND relation = 'for'
        """,
        (hypothesis_id,),
    ).fetchall()
    assert {row[0] for row in normalized} == set(finding_ids)


def test_idempotency_covers_findings_hypotheses_and_log(tmp_path):
    ops = OpsStore(tmp_path / "ops.db", agent_id="coordinator")
    first = ops.finding(
        "observation",
        "fixture",
        "idempotent finding",
        idempotency_key="finding:1",
    )
    assert (
        ops.finding(
            "observation",
            "fixture",
            "idempotent finding",
            idempotency_key="finding:1",
        )
        == first
    )
    with pytest.raises(IdempotencyConflictError):
        ops.finding(
            "observation",
            "fixture",
            "changed payload",
            idempotency_key="finding:1",
        )

    hypothesis = ops.hypothesis(
        "idempotent hypothesis",
        [first],
        idempotency_key="hypothesis:1",
    )
    assert (
        ops.hypothesis(
            "idempotent hypothesis",
            [first],
            idempotency_key="hypothesis:1",
        )
        == hypothesis
    )
    log_id = ops.log(
        "dispatch",
        "idempotent decision",
        referenced_finding_ids=[first],
        referenced_hypothesis_ids=[hypothesis],
        idempotency_key="log:1",
    )
    assert (
        ops.log(
            "dispatch",
            "idempotent decision",
            referenced_finding_ids=[first],
            referenced_hypothesis_ids=[hypothesis],
            idempotency_key="log:1",
        )
        == log_id
    )


def test_normalized_task_dependencies_and_results_enforce_references(tmp_path):
    ops = OpsStore(tmp_path / "ops.db", agent_id="coordinator")
    parent = ops.create_task("prepare", "prepare")
    child = ops.create_task("analyze", "analyze", depends_on=[parent])
    relation = ops._conn.execute(
        "SELECT task_id, depends_on_task_id FROM task_dependencies"
    ).fetchone()
    assert tuple(relation) == (child, parent)

    worker = OpsStore(ops.db_path, agent_id="worker")
    worker.claim_task(parent)
    worker.complete_task(parent)
    worker.claim_task(child)
    with pytest.raises(ValueError, match="Unknown result findings"):
        worker.complete_task(child, ["missing"])

    finding_id = worker.finding("observation", "fixture", "task result")
    completed = worker.complete_task(child, [finding_id])
    assert completed["result_finding_ids"] == [finding_id]
    assert (
        ops._conn.execute(
            "SELECT finding_id FROM task_result_findings WHERE task_id = ?",
            (child,),
        ).fetchone()[0]
        == finding_id
    )


def test_claim_next_respects_capabilities_capacity_and_concurrency(tmp_path):
    db_path = tmp_path / "ops.db"
    coordinator = OpsStore(db_path, agent_id="operator")
    coordinator.register_agent(
        "worker",
        capabilities=["duckdb"],
        max_concurrent_tasks=2,
        capacity_cpu_slots=4,
        capacity_memory_mb=4096,
        capacity_accelerator_slots=0,
    )
    incompatible = coordinator.create_task(
        "gpu",
        "needs an unavailable accelerator",
        priority=3,
        required_capabilities=["accelerated"],
        resource_request={"accelerator_slots": 1},
    )
    first_id = coordinator.create_task(
        "query",
        "three slots",
        priority=2,
        required_capabilities=["duckdb"],
        resource_request={"cpu_slots": 3, "memory_mb": 1024},
    )
    second_id = coordinator.create_task(
        "query",
        "two slots",
        priority=1,
        required_capabilities=["duckdb"],
        resource_request={"cpu_slots": 2},
    )
    worker = OpsStore(db_path, agent_id="worker")

    first = worker.claim_next_task()
    assert first["id"] == first_id
    assert worker.claim_next_task() is None
    worker.complete_task(first_id)
    second = worker.claim_next_task()
    assert second["id"] == second_id
    assert coordinator.get_task(incompatible)["status"] == "pending"


def test_claim_next_filters_incompatible_work_before_bounded_candidate_scan(tmp_path):
    db_path = tmp_path / "ops.db"
    coordinator = OpsStore(db_path, agent_id="operator")
    coordinator.register_agent(
        "worker",
        capabilities=["duckdb"],
        capacity_cpu_slots=2,
    )
    for index in range(120):
        coordinator.create_task(
            "accelerated",
            f"incompatible task {index}",
            priority=3,
            required_capabilities=["accelerated"],
            resource_request={"accelerator_slots": 1},
        )
    compatible = coordinator.create_task(
        "query",
        "compatible task after the bounded prefix",
        priority=1,
        required_capabilities=["duckdb"],
        resource_request={"cpu_slots": 1},
    )

    worker = OpsStore(db_path, agent_id="worker")
    claimed = worker.claim_next_task(candidate_limit=100)
    assert claimed["id"] == compatible


def test_content_addressed_artifacts_fts_and_inline_limit(tmp_path):
    ops = OpsStore(tmp_path / "ops.db", agent_id="worker")
    finding_id = ops.finding(
        "observation",
        "fixture",
        "rare coordination motif",
        evidence={"method": "fixture"},
    )
    matches = ops.search_findings("rare coordination")
    assert [row["id"] for row in matches] == [finding_id]

    artifact_id = ops.register_artifact(
        "a" * 64,
        "file:///analysis/result.parquet",
        42,
        media_type="application/vnd.apache.parquet",
    )
    assert (
        ops.register_artifact(
            "sha256:" + "a" * 64,
            "file:///analysis/result.parquet",
            42,
        )
        == artifact_id
    )
    ops.attach_artifact(finding_id, artifact_id)
    ops.attach_artifact(finding_id, artifact_id)
    assert (
        ops._conn.execute(
            "SELECT artifact_id FROM finding_artifacts WHERE finding_id = ?",
            (finding_id,),
        ).fetchone()[0]
        == artifact_id
    )
    assert (
        ops._conn.execute(
            """
            SELECT COUNT(*) FROM coordination_events
            WHERE event_type = 'artifact_attached' AND entity_id = ?
            """,
            (finding_id,),
        ).fetchone()[0]
        == 1
    )
    with pytest.raises(ValueError, match="inline JSON limit"):
        ops.finding(
            "observation",
            "fixture",
            "oversized",
            evidence={"blob": "x" * MAX_INLINE_JSON_BYTES},
        )


def test_schema_current_open_has_zero_schema_writes(tmp_path):
    db_path = tmp_path / "ops.db"
    OpsStore(db_path).close()
    conn = open_ops_connection(db_path)
    before = conn.total_changes
    assert ensure_ops_schema(conn) is False
    assert conn.total_changes == before
    assert (
        conn.execute("SELECT MAX(version) FROM ops_schema_meta").fetchone()[0]
        == OPS_SCHEMA_VERSION
    )
    conn.close()


def test_online_backup_is_integral_and_independent(tmp_path):
    ops = OpsStore(tmp_path / "ops.db", agent_id="operator")
    finding_id = ops.finding("observation", "fixture", "back me up")
    backup_path = tmp_path / "backups" / "ops.db"
    result = ops.backup(backup_path)

    assert result["integrity"] == "ok"
    backup = OpsStore(backup_path, agent_id="reader")
    assert backup.get_finding(finding_id)["summary"] == "back me up"
    assert backup.integrity_check()["ok"] is True


def test_stats_expose_terminal_failures_as_dead_letters(tmp_path):
    ops = OpsStore(tmp_path / "ops.db", agent_id="worker")
    task_id = ops.create_task(
        "fixture",
        "exhaust the only attempt",
        max_attempts=1,
        lease_seconds=1,
    )
    claimed = ops.claim_task(task_id)
    recovered = ops.recover_expired_tasks(now=claimed["lease_expires_ts"] + 1)

    assert recovered["exhausted"] == [task_id]
    dead_letters = ops.stats()["dead_letters"]
    assert dead_letters["count"] == 1
    assert dead_letters["attempts_exhausted"] == 1
    assert dead_letters["oldest_age_seconds"] >= 0


def test_queue_read_is_bounded_at_scale_and_uses_no_writer_transaction(tmp_path):
    ops = OpsStore(tmp_path / "ops.db", agent_id="reader")
    now = time.time()
    with ops._lock, ops._transaction():
        ops._conn.executemany(
            """
            INSERT INTO tasks(
                id, created_by, ts, status, priority, task_type, description
            ) VALUES (?, 'seed', ?, 'pending', 1, 'fixture', 'fixture')
            """,
            [(str(uuid.uuid4()), now + index / 10_000) for index in range(5_000)],
        )
    before = ops._conn.total_changes
    started = time.perf_counter()
    available = ops.available_tasks(limit=25)
    elapsed = time.perf_counter() - started

    assert len(available) == 25
    assert ops._conn.total_changes == before
    assert elapsed < 1.0


def test_concurrent_claim_next_returns_unique_tasks(tmp_path):
    db_path = tmp_path / "ops.db"
    coordinator = OpsStore(db_path, agent_id="coordinator")
    task_ids = {
        coordinator.create_task("fixture", f"task {index}")
        for index in range(24)
    }
    workers = [OpsStore(db_path, agent_id=f"worker-{index}") for index in range(12)]

    def claim_two(worker):
        return [worker.claim_next_task()["id"] for _ in range(2)]

    with ThreadPoolExecutor(max_workers=12) as executor:
        claimed = [
            task_id
            for pair in executor.map(claim_two, workers)
            for task_id in pair
        ]
    assert len(claimed) == 24
    assert set(claimed) == task_ids


def test_pool_is_bounded_and_server_owner_lock_is_exclusive(tmp_path):
    db_path = tmp_path / "ops.db"
    pool = SQLiteConnectionPool(db_path, size=2)
    first = pool.acquire()
    second = pool.acquire()
    assert pool.stats()["checked_out"] == 2
    pool.release(first)
    pool.release(second)
    pool.observe_sqlite_write_wait(0.125)
    stats = pool.stats()
    assert stats["available"] == 2
    assert stats["wait_seconds_total"] >= 0
    assert stats["sqlite_write_transactions"] == 1
    assert stats["sqlite_write_wait_seconds"] == pytest.approx(0.125)
    pool.close()

    first_lock = SQLiteServerLock(db_path)
    second_lock = SQLiteServerLock(db_path)
    first_lock.acquire()
    try:
        with pytest.raises(RuntimeError, match="already has an HTTP owner"):
            second_lock.acquire()
    finally:
        first_lock.release()


def test_parallel_run_stages_keep_a_live_current_stage(tmp_path):
    ops = OpsStore(tmp_path / "ops.db", agent_id="operator")
    run_id = ops.create_run("analysis", tmp_path / "dataset")
    ops.start_run(run_id)
    first = ops.ledger.start_stage(run_id, "a", "sig-a")
    second = ops.ledger.start_stage(run_id, "b", "sig-b")

    ops.ledger.complete_stage(run_id, "b", second)
    assert ops.get_run(run_id)["current_stage"] == "a"
    ops.ledger.complete_stage(run_id, "a", first)
    assert ops.get_run(run_id)["current_stage"] is None


def test_migrated_legacy_database_enables_foreign_keys(tmp_path):
    db_path = tmp_path / "ops.db"
    conn = sqlite3.connect(db_path)
    conn.execute(
        """
        CREATE TABLE findings (
            id TEXT PRIMARY KEY, agent_id TEXT NOT NULL, ts REAL NOT NULL,
            finding_type TEXT NOT NULL, domain TEXT NOT NULL, summary TEXT NOT NULL,
            evidence TEXT NOT NULL DEFAULT '{}', confidence REAL NOT NULL DEFAULT 0.5,
            novelty INTEGER NOT NULL DEFAULT 0, parent_finding_id TEXT,
            reasoning TEXT NOT NULL DEFAULT ''
        )
        """
    )
    conn.commit()
    conn.close()

    ops = OpsStore(db_path)
    assert ops._conn.execute("PRAGMA foreign_keys").fetchone()[0] == 1
    assert ensure_ops_schema(ops._conn) is False
