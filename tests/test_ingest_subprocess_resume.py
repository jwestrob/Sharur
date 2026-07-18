"""Process-boundary failure/resume coverage for the public ingest CLI."""

from __future__ import annotations

import os
import signal
import sqlite3
import subprocess
import sys
import time
from typing import TYPE_CHECKING

from sharur.ops.ledger import RunLedger


if TYPE_CHECKING:
    from pathlib import Path


def _fast_ingest_command(input_dir: Path, data_dir: Path, output: Path) -> list[str]:
    return [
        sys.executable,
        "-m",
        "sharur.ingest_cli",
        "--mode",
        "fast",
        "--profile",
        "local",
        "--input-dir",
        str(input_dir),
        "--data-dir",
        str(data_dir),
        "--output",
        str(output),
    ]


def _run_fast_ingest(input_dir: Path, data_dir: Path, output: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        _fast_ingest_command(input_dir, data_dir, output),
        cwd=data_dir.parent,
        text=True,
        capture_output=True,
        timeout=60,
        check=False,
    )


def _wait_for_running_stage(ledger_path: Path, stage_id: str, process: subprocess.Popen) -> None:
    deadline = time.monotonic() + 30
    while time.monotonic() < deadline:
        if process.poll() is not None:
            stdout, stderr = process.communicate()
            raise AssertionError(
                f"ingest exited before stage {stage_id} could be killed: "
                f"returncode={process.returncode}\nstdout={stdout}\nstderr={stderr}"
            )
        if ledger_path.is_file():
            try:
                connection = sqlite3.connect(ledger_path, timeout=1)
                try:
                    row = connection.execute(
                        "SELECT 1 FROM run_stages WHERE stage_id = ? AND status = 'running'",
                        (stage_id,),
                    ).fetchone()
                finally:
                    connection.close()
                if row is not None:
                    return
            except sqlite3.Error:
                pass
        time.sleep(0.01)
    raise AssertionError(f"timed out waiting for running stage {stage_id}")


def test_cli_hard_crash_resume_reuses_completed_process_stage(tmp_path: Path) -> None:
    input_dir = tmp_path / "input"
    data_dir = tmp_path / "dataset"
    input_dir.mkdir()
    records = "".join(f">contig{index}\n{'ATG' * 100}\n" for index in range(20))
    (input_dir / "tiny.fna").write_text(records)

    output = data_dir / "sharur.duckdb"
    crashed = subprocess.Popen(
        _fast_ingest_command(input_dir, data_dir, output),
        cwd=data_dir.parent,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        start_new_session=True,
    )
    ledger_path = data_dir / "sharur_ops.db"
    try:
        _wait_for_running_stage(ledger_path, "07", crashed)
    finally:
        if crashed.poll() is None:
            os.killpg(crashed.pid, signal.SIGKILL)
        crashed.communicate(timeout=10)
    assert crashed.returncode == -signal.SIGKILL

    # Advance the durable clock without a five-minute test sleep. The next
    # public CLI invocation must recover this real abandoned process attempt.
    connection = sqlite3.connect(ledger_path)
    try:
        connection.execute(
            "UPDATE run_stages SET heartbeat_ts = 0, started_ts = 0 "
            "WHERE status = 'running'"
        )
        connection.execute(
            "UPDATE runs SET heartbeat_ts = 0, started_ts = 0 "
            "WHERE status = 'running'"
        )
        connection.commit()
    finally:
        connection.close()

    resumed = _run_fast_ingest(input_dir, data_dir, output)
    assert resumed.returncode == 0, resumed.stderr
    assert output.is_file()

    fully_reused = _run_fast_ingest(input_dir, data_dir, output)
    assert fully_reused.returncode == 0, fully_reused.stderr

    ledger = RunLedger(data_dir / "sharur_ops.db")
    try:
        runs = list(
            reversed(
                ledger.list_runs(
                    dataset_path=data_dir,
                    run_type="ingest",
                    limit=10,
                )
            )
        )
        assert [run["status"] for run in runs] == ["failed", "complete", "complete"]
        first_stages = ledger.list_stages(runs[0]["id"])
        resumed_stages = ledger.list_stages(runs[1]["id"])
        final_stages = ledger.list_stages(runs[2]["id"])
        first_events = ledger.events(runs[0]["id"])
    finally:
        ledger.close()

    assert {stage["stage_id"]: stage["status"] for stage in first_stages} == {
        "00": "complete",
        "07": "failed",
    }
    assert any(
        event["event_type"] == "stage_recovered_after_crash"
        for event in first_events
    )
    assert {stage["stage_id"]: stage["status"] for stage in resumed_stages} == {
        "00": "reused",
        "07": "complete",
    }
    assert {stage["stage_id"]: stage["status"] for stage in final_stages} == {
        "00": "reused",
        "07": "reused",
    }
