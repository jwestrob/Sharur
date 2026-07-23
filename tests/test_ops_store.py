"""Crash- and concurrency-sensitive tests for the direct ops backend."""

import sqlite3
import sys
from concurrent.futures import ThreadPoolExecutor

import pytest

from sharur.ingest.resources import ResourceRequest
from sharur.ingest.stage_runner import execute_stage
from sharur.ops.ledger import RunLedger
from sharur.ops.store import LeaseFenceError, OpsStore


def test_preassigned_task_is_reserved_until_worker_claims(tmp_path):
    coordinator = OpsStore(
        tmp_path / "ops" / "sharur_ops.db",
        agent_id="coordinator",
    )
    task_id = coordinator.create_task(
        "survey",
        "Inspect a dataset",
        assigned_to="worker",
    )

    worker = OpsStore(
        coordinator.db_path,
        agent_id="worker",
    )
    task = worker.my_tasks()[0]

    assert task["id"] == task_id
    assert task["status"] == "pending"
    assert task["reserved_for"] == "worker"
    assert task["claimed_ts"] is None
    assert coordinator.available_tasks() == []
    claimed = worker.claim_task(task_id)
    assert claimed["status"] == "claimed"
    assert claimed["lease_attempt"] == 1


def test_only_one_agent_can_claim_task(tmp_path):
    db_path = tmp_path / "sharur_ops.db"
    coordinator = OpsStore(db_path, agent_id="coordinator")
    task_id = coordinator.create_task("survey", "Inspect a dataset")
    workers = [
        OpsStore(db_path, agent_id="worker_a"),
        OpsStore(db_path, agent_id="worker_b"),
    ]

    def claim(worker):
        try:
            return worker.claim_task(task_id)["assigned_to"]
        except ValueError:
            return None

    with ThreadPoolExecutor(max_workers=2) as executor:
        outcomes = list(executor.map(claim, workers))

    assert sum(outcome is not None for outcome in outcomes) == 1


def test_task_completion_requires_current_owner(tmp_path):
    db_path = tmp_path / "sharur_ops.db"
    coordinator = OpsStore(db_path, agent_id="coordinator")
    task_id = coordinator.create_task(
        "survey",
        "Inspect a dataset",
        assigned_to="worker",
    )

    worker = OpsStore(db_path, agent_id="worker")
    claimed = worker.claim_task(task_id)
    intruder = OpsStore(db_path, agent_id="intruder")
    with pytest.raises(LeaseFenceError, match="not active"):
        intruder.complete_task(
            task_id,
            lease_token=claimed["lease_token"],
            attempt=claimed["lease_attempt"],
        )

    finding_id = worker.finding("observation", "fixture", "result")
    completed = worker.complete_task(task_id, [finding_id])
    assert completed["status"] == "complete"
    assert completed["result_finding_ids"] == [finding_id]


def test_completing_missing_task_has_typed_error(tmp_path):
    ops = OpsStore(tmp_path / "sharur_ops.db", agent_id="worker")

    with pytest.raises(LeaseFenceError, match="not active"):
        ops.complete_task("missing", lease_token="x" * 32, attempt=1)


def test_hypothesis_evidence_updates_append_and_deduplicate(tmp_path):
    ops = OpsStore(tmp_path / "sharur_ops.db", agent_id="worker")
    hypothesis_id = ops.hypothesis("A testable claim")
    f1 = ops.finding("observation", "fixture", "first")
    f2 = ops.finding("observation", "fixture", "second")

    ops.update_hypothesis(hypothesis_id, evidence_for=[f1])
    updated = ops.update_hypothesis(
        hypothesis_id,
        evidence_for=[f1, f2],
    )

    assert updated["evidence_for"] == [f1, f2]


def test_updating_missing_hypothesis_raises_key_error(tmp_path):
    ops = OpsStore(tmp_path / "sharur_ops.db", agent_id="worker")

    with pytest.raises(KeyError, match="not found"):
        ops.update_hypothesis("missing", status="investigating")


def test_task_dependencies_block_claim_until_parent_completes(tmp_path):
    db_path = tmp_path / "sharur_ops.db"
    coordinator = OpsStore(db_path, agent_id="coordinator")
    parent = coordinator.create_task("prepare", "Prepare inputs")
    child = coordinator.create_task(
        "analyze",
        "Analyze prepared inputs",
        depends_on=[parent],
    )
    worker = OpsStore(db_path, agent_id="worker")

    assert [task["id"] for task in worker.available_tasks()] == [parent]
    with pytest.raises(ValueError, match="incomplete dependencies"):
        worker.claim_task(child)

    worker.claim_task(parent)
    worker.complete_task(parent)
    assert [task["id"] for task in worker.available_tasks()] == [child]


def test_expired_lease_is_recovered_then_exhausted(tmp_path):
    db_path = tmp_path / "sharur_ops.db"
    coordinator = OpsStore(db_path, agent_id="coordinator")
    task_id = coordinator.create_task(
        "survey",
        "Recover me",
        max_attempts=2,
        lease_seconds=10,
    )

    first = OpsStore(db_path, agent_id="worker_a").claim_task(task_id)
    recovered = coordinator.recover_expired_tasks(
        now=first["lease_expires_ts"] + 1,
    )
    assert recovered == {"recovered": [task_id], "exhausted": []}

    second = OpsStore(db_path, agent_id="worker_b").claim_task(task_id)
    exhausted = coordinator.recover_expired_tasks(
        now=second["lease_expires_ts"] + 1,
    )
    assert exhausted == {"recovered": [], "exhausted": [task_id]}
    assert coordinator._get_task(task_id)["status"] == "failed"


def test_heartbeat_extends_owned_lease(tmp_path):
    db_path = tmp_path / "sharur_ops.db"
    coordinator = OpsStore(db_path, agent_id="coordinator")
    task_id = coordinator.create_task("survey", "Long task", lease_seconds=10)
    worker = OpsStore(db_path, agent_id="worker")
    claimed = worker.claim_task(task_id)

    heartbeat = worker.heartbeat_task(task_id, lease_seconds=60)

    assert heartbeat["status"] == "in_progress"
    assert heartbeat["heartbeat_ts"] >= claimed["heartbeat_ts"]
    assert heartbeat["lease_expires_ts"] > claimed["lease_expires_ts"]


def test_retryable_failure_requeues_without_losing_attempt_count(tmp_path):
    db_path = tmp_path / "sharur_ops.db"
    coordinator = OpsStore(db_path, agent_id="coordinator")
    task_id = coordinator.create_task("survey", "Flaky", max_attempts=2)
    worker = OpsStore(db_path, agent_id="worker")
    worker.claim_task(task_id)

    waiting = worker.fail_task(task_id, error="transient", retryable=True)
    assert waiting["status"] == "retry_wait"
    assert waiting["attempt_count"] == 1

    claimed_again = OpsStore(db_path, agent_id="worker_2").claim_task(task_id)
    assert claimed_again["attempt_count"] == 2


def test_task_idempotency_reuses_exact_request_and_rejects_conflict(tmp_path):
    ops = OpsStore(tmp_path / "sharur_ops.db", agent_id="coordinator")
    first = ops.create_task(
        "survey",
        "Inspect",
        params={"dataset": "x"},
        idempotency_key="survey:x",
    )
    second = ops.create_task(
        "survey",
        "Inspect",
        params={"dataset": "x"},
        idempotency_key="survey:x",
    )
    assert second == first

    with pytest.raises(ValueError, match="different payload"):
        ops.create_task(
            "survey",
            "Different",
            idempotency_key="survey:x",
        )
    with pytest.raises(ValueError, match="non-empty"):
        ops.create_task("survey", "Empty key", idempotency_key="")


def test_legacy_task_schema_migrates_in_place(tmp_path):
    db_path = tmp_path / "sharur_ops.db"
    conn = sqlite3.connect(db_path)
    conn.execute(
        """
        CREATE TABLE tasks (
            id TEXT PRIMARY KEY, created_by TEXT NOT NULL, assigned_to TEXT,
            ts REAL NOT NULL, claimed_ts REAL, completed_ts REAL,
            status TEXT NOT NULL DEFAULT 'pending', priority INTEGER NOT NULL,
            task_type TEXT NOT NULL, description TEXT NOT NULL,
            params TEXT NOT NULL DEFAULT '{}', domain_hint TEXT,
            result_finding_ids TEXT NOT NULL DEFAULT '[]'
        )
        """
    )
    conn.execute(
        """
        INSERT INTO tasks(
            id, created_by, assigned_to, ts, claimed_ts, status, priority,
            task_type, description
        ) VALUES ('legacy', 'old', 'dead_worker', 1, 1, 'claimed', 1, 'x', 'old')
        """
    )
    conn.commit()
    conn.close()

    ops = OpsStore(db_path, agent_id="operator")
    columns = {row[1] for row in ops._conn.execute("PRAGMA table_info(tasks)").fetchall()}
    assert {"lease_expires_ts", "attempt_count", "dependency_ids"} <= columns
    assert ops.recover_expired_tasks(now=10_000)["recovered"] == ["legacy"]


def test_run_ledger_records_stage_attempts_and_reuse(tmp_path):
    ledger = RunLedger(tmp_path / "sharur_ops.db")
    dataset = tmp_path / "dataset"
    run_id = ledger.create_run(
        "ingest",
        dataset,
        config={"profile": "local"},
        idempotency_key="ingest:test",
    )
    assert (
        ledger.create_run(
            "ingest",
            dataset,
            config={"profile": "local"},
            idempotency_key="ingest:test",
        )
        == run_id
    )
    with pytest.raises(ValueError, match="different payload"):
        ledger.create_run(
            "ingest",
            dataset,
            config={"profile": "local"},
            idempotency_key="ingest:test",
            parent_run_id="different-parent",
        )
    with pytest.raises(ValueError, match="non-empty"):
        ledger.create_run("ingest", dataset, idempotency_key="")
    ledger.start_run(run_id)
    attempt = ledger.start_stage(
        run_id,
        "03",
        "signature-a",
        command=["stage03"],
    )
    ledger.complete_stage(run_id, "03", attempt, outputs={"path": "stage03"})
    ledger.complete_run(run_id)

    source = ledger.find_reusable_stage(dataset, "ingest", "03", "signature-a")
    assert source is not None
    next_run = ledger.create_run("ingest", dataset)
    ledger.start_run(next_run)
    ledger.reuse_stage(next_run, "03", "signature-a", source)
    ledger.reuse_stage(next_run, "03", "signature-a", source)
    events = ledger.events(next_run)

    assert [event["event_type"] for event in events] == [
        "run_created",
        "run_started",
        "stage_reused",
    ]
    assert len(ledger.list_stages(next_run, stage_id="03")) == 1


def test_run_creation_and_start_are_atomic_across_connections(tmp_path):
    db_path = tmp_path / "sharur_ops.db"
    ledgers = [RunLedger(db_path), RunLedger(db_path)]

    def create(ledger):
        return ledger.create_run(
            "ingest",
            tmp_path / "dataset",
            config={"profile": "local"},
            idempotency_key="same-run",
        )

    with ThreadPoolExecutor(max_workers=2) as executor:
        run_ids = list(executor.map(create, ledgers))
    assert len(set(run_ids)) == 1

    def start(ledger):
        try:
            ledger.start_run(run_ids[0])
            return True
        except ValueError:
            return False

    with ThreadPoolExecutor(max_workers=2) as executor:
        starts = list(executor.map(start, ledgers))
    assert starts.count(True) == 1
    assert [event["event_type"] for event in ledgers[0].events(run_ids[0])].count(
        "run_started"
    ) == 1

    for ledger in ledgers:
        ledger.close()


def test_duplicate_stage_launcher_is_suppressed_before_cleanup(tmp_path):
    ledger = RunLedger(tmp_path / "sharur_ops.db")
    run_id = ledger.create_run("ingest", tmp_path / "dataset")
    ledger.start_run(run_id)
    script = tmp_path / "stage.py"
    output = tmp_path / "output.txt"
    counter = tmp_path / "counter.txt"
    script.write_text(
        "from pathlib import Path\n"
        "import sys\n"
        "output, counter = map(Path, sys.argv[1:])\n"
        "count = int(counter.read_text()) + 1 if counter.exists() else 1\n"
        "counter.write_text(str(count))\n"
        "output.write_text(f'run-{count}')\n"
    )
    kwargs = {
        "ledger": ledger,
        "run_id": run_id,
        "stage_id": "03",
        "signature": "same-stage",
        "command": [sys.executable, str(script), str(output), str(counter)],
        "outputs": [output],
        "cleanup": [output, counter],
        "resource": ResourceRequest(1, 1, "00:10:00"),
    }

    execute_stage(**kwargs)
    execute_stage(**kwargs)

    assert output.read_text() == "run-1"
    assert counter.read_text() == "1"
    assert any(
        event["event_type"] == "duplicate_stage_launch_suppressed"
        for event in ledger.events(run_id)
    )
    output.write_text("tampered-but-nonempty")
    with pytest.raises(RuntimeError, match="outputs changed"):
        execute_stage(**kwargs)
    assert counter.read_text() == "1"
    ledger.close()


def test_submitted_run_waits_for_scheduler_stage_without_stale_recovery(tmp_path):
    ledger = RunLedger(tmp_path / "sharur_ops.db")
    run_id = ledger.create_run("ingest", tmp_path / "dataset")
    submitted = ledger.submit_run(run_id)

    assert submitted["status"] == "submitted"
    assert ledger.recover_stale_runs(stale_after_seconds=1, now=10**12) == {
        "runs": [],
        "stages": [],
    }
    attempt = ledger.start_stage(run_id, "03", "signature")
    assert attempt == 1
    assert ledger.get_run(run_id)["status"] == "running"
    ledger.complete_stage(run_id, "03", attempt)
    waiting = ledger.wait_for_scheduler(run_id)
    assert waiting["status"] == "submitted"
    assert ledger.recover_stale_runs(stale_after_seconds=1, now=10**12) == {
        "runs": [],
        "stages": [],
    }
    ledger.close()


def test_ops_stats_include_run_and_stage_lifecycle(tmp_path):
    ops = OpsStore(tmp_path / "sharur_ops.db", agent_id="operator")
    run_id = ops.create_run("ingest", tmp_path / "dataset")
    ops.start_run(run_id)
    attempt = ops.ledger.start_stage(run_id, "03", "signature")
    ops.ledger.complete_stage(run_id, "03", attempt)
    ops.complete_run(run_id)

    stats = ops.stats()

    assert stats["counts"]["run_stages"] == 1
    assert stats["runs_by_status"] == {"complete": 1}
    assert stats["stages_by_status"] == {"complete": 1}
    ops.close()


def test_run_ledger_recovers_stage_after_hard_crash_heartbeat_gap(tmp_path):
    ledger = RunLedger(tmp_path / "sharur_ops.db")
    run_id = ledger.create_run("ingest", tmp_path / "dataset")
    ledger.start_run(run_id)
    attempt = ledger.start_stage(run_id, "04", "signature", command=["astra"])
    ledger._conn.execute(
        """
        UPDATE run_stages SET heartbeat_ts = 10
        WHERE run_id = ? AND stage_id = '04' AND attempt = ?
        """,
        (run_id, attempt),
    )
    ledger._conn.execute(
        "UPDATE runs SET heartbeat_ts = 10 WHERE id = ?",
        (run_id,),
    )
    ledger._conn.commit()

    recovered = ledger.recover_stale_runs(
        stale_after_seconds=300,
        now=1_000,
    )

    assert recovered["runs"] == [run_id]
    assert recovered["stages"] == [{"run_id": run_id, "stage_id": "04", "attempt": attempt}]
    assert ledger.get_run(run_id)["status"] == "failed"
    assert (
        ledger.find_reusable_stage(
            tmp_path / "dataset",
            "ingest",
            "04",
            "signature",
        )
        is None
    )
