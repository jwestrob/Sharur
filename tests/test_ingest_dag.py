"""Dependency-aware ingest resume and resource-profile tests."""

import subprocess
from pathlib import Path

import pytest

from sharur.ingest.dag import IngestDAG, StageNode
from sharur.ingest.resources import resolve_resource_profile
from sharur.ingest_cli import (
    StagePaths,
    _build_tools_dag,
    _emit_slurm_bundle,
    _execute_local_dag,
)
from sharur.ops.ledger import RunLedger


def _single_stage_dag(tmp_path: Path) -> tuple[IngestDAG, Path, Path]:
    script = tmp_path / "stage.py"
    output = tmp_path / "dataset" / "artifact.txt"
    counter = tmp_path / "counter.txt"
    script.write_text(
        "from pathlib import Path\n"
        "import sys\n"
        "output, counter = map(Path, sys.argv[1:])\n"
        "count = int(counter.read_text()) + 1 if counter.exists() else 1\n"
        "counter.write_text(str(count))\n"
        "output.parent.mkdir(parents=True, exist_ok=True)\n"
        "output.write_text(f'run-{count}')\n"
    )
    profile = resolve_resource_profile("local")
    node = StageNode(
        stage_id="00",
        label="fixture",
        command=[str(script), str(output), str(counter)],
        required_outputs=(output,),
        input_paths=(script,),
        resource=profile.request("00"),
    )
    dag = IngestDAG([node])
    dag.compute_signatures(profile)
    return dag, output, counter


def test_resume_reuses_only_matching_successful_stage_with_live_outputs(tmp_path):
    dag, output, counter = _single_stage_dag(tmp_path)
    data_dir = tmp_path / "dataset"
    profile = resolve_resource_profile("local")

    _execute_local_dag(
        dag,
        data_dir=data_dir,
        profile=profile,
        resume=True,
        force=False,
        idempotency_key=None,
    )
    _execute_local_dag(
        dag,
        data_dir=data_dir,
        profile=profile,
        resume=True,
        force=False,
        idempotency_key=None,
    )

    assert output.read_text() == "run-1"
    assert counter.read_text() == "1"
    ledger = RunLedger(data_dir / "sharur_ops.db")
    try:
        runs = ledger.list_runs(dataset_path=data_dir, run_type="ingest")
        newest_events = ledger.events(runs[0]["id"])
    finally:
        ledger.close()
    assert newest_events[-1]["event_type"] == "run_completed"
    assert any(event["event_type"] == "stage_reused" for event in newest_events)


def test_missing_output_invalidates_otherwise_matching_resume(tmp_path):
    dag, output, counter = _single_stage_dag(tmp_path)
    data_dir = tmp_path / "dataset"
    profile = resolve_resource_profile("local")
    kwargs = {
        "dag": dag,
        "data_dir": data_dir,
        "profile": profile,
        "resume": True,
        "force": False,
        "idempotency_key": None,
    }
    _execute_local_dag(**kwargs)
    output.unlink()
    _execute_local_dag(**kwargs)

    assert output.read_text() == "run-2"
    assert counter.read_text() == "2"


def test_modified_output_invalidates_otherwise_matching_resume(tmp_path):
    dag, output, counter = _single_stage_dag(tmp_path)
    data_dir = tmp_path / "dataset"
    profile = resolve_resource_profile("local")
    kwargs = {
        "dag": dag,
        "data_dir": data_dir,
        "profile": profile,
        "resume": True,
        "force": False,
        "idempotency_key": None,
    }
    _execute_local_dag(**kwargs)
    output.write_text("tampered-but-nonempty")
    _execute_local_dag(**kwargs)

    assert output.read_text() == "run-2"
    assert counter.read_text() == "2"


def test_rerun_upstream_invalidates_downstream_reuse(tmp_path):
    data_dir = tmp_path / "dataset"
    profile = resolve_resource_profile("local")
    upstream_output = data_dir / "upstream.txt"
    downstream_output = data_dir / "downstream.txt"
    upstream_counter = tmp_path / "upstream-counter.txt"
    downstream_counter = tmp_path / "downstream-counter.txt"
    script = tmp_path / "stage.py"
    script.write_text(
        "from pathlib import Path\n"
        "import sys\n"
        "output, counter = map(Path, sys.argv[1:])\n"
        "count = int(counter.read_text()) + 1 if counter.exists() else 1\n"
        "counter.write_text(str(count))\n"
        "output.parent.mkdir(parents=True, exist_ok=True)\n"
        "output.write_text(f'run-{count}')\n"
    )
    nodes = [
        StageNode(
            stage_id="00",
            label="upstream",
            command=[str(script), str(upstream_output), str(upstream_counter)],
            required_outputs=(upstream_output,),
            resource=profile.request("00"),
        ),
        StageNode(
            stage_id="03",
            label="downstream",
            command=[str(script), str(downstream_output), str(downstream_counter)],
            dependencies=("00",),
            required_outputs=(downstream_output,),
            resource=profile.request("03"),
        ),
    ]
    dag = IngestDAG(nodes)
    dag.compute_signatures(profile)
    kwargs = {
        "dag": dag,
        "data_dir": data_dir,
        "profile": profile,
        "resume": True,
        "force": False,
        "idempotency_key": None,
    }

    _execute_local_dag(**kwargs)
    upstream_output.write_text("tampered")
    _execute_local_dag(**kwargs)

    assert upstream_counter.read_text() == "2"
    assert downstream_counter.read_text() == "2"


def test_recursive_directory_snapshot_detects_nested_file_change(tmp_path):
    directory = tmp_path / "output"
    nested = directory / "genome" / "proteins.faa"
    nested.parent.mkdir(parents=True)
    nested.write_text(">p1\nMA\n")
    node = StageNode(
        stage_id="03",
        label="nested",
        command=["unused"],
        required_outputs=(directory,),
    )
    recorded = node.output_snapshot()

    nested.write_text(">p1\nVV\n")

    assert node.outputs_match(recorded) is False


def test_resource_profiles_encode_mps_exclusivity_and_slurm_limits(monkeypatch):
    monkeypatch.setattr("sharur.ingest.resources._mps_available", lambda: True)
    mps = resolve_resource_profile("mps")
    slurm = resolve_resource_profile("slurm")

    assert mps.request("06").accelerator == "mps"
    assert mps.request("06").exclusive_accelerator is True
    assert slurm.request("04").walltime == "72:00:00"
    assert slurm.request("06").gpus == 1
    assert slurm.request("06i").gpus == 0
    assert slurm.request("06i").accelerator == "cpu"
    assert slurm.request("05c").executor == "local"
    assert slurm.request("05c").cpus == 1


def test_auto_profile_falls_back_to_cpu_without_usable_mps(monkeypatch):
    monkeypatch.setattr("sharur.ingest.resources._mps_available", lambda: False)
    assert resolve_resource_profile("auto").name == "local"


def test_slurm_bundle_preserves_dependencies_and_keeps_minced_local(tmp_path):
    data_dir = tmp_path / "dataset"
    stages = StagePaths.from_root(data_dir)
    profile = resolve_resource_profile("slurm")
    dag = _build_tools_dag(
        input_dir=tmp_path / "inputs",
        data_dir=data_dir,
        output=data_dir / "sharur.duckdb",
        stages=stages,
        profile=profile,
        skip_quast=True,
        skip_dfast=True,
        skip_prodigal=False,
        skip_astra=False,
        skip_gecco=True,
        skip_dbcan=True,
        skip_crispr=False,
        skip_embeddings=False,
        enable_cazymes=False,
    )

    submit = _emit_slurm_bundle(
        dag,
        data_dir=data_dir,
        profile=profile,
        resume=False,
        force=False,
        idempotency_key=None,
        submit=False,
    )

    assert (data_dir / "slurm" / "05c.local.sh").is_file()
    astra_script = (data_dir / "slurm" / "04.sbatch").read_text()
    assert "#SBATCH --time=72:00:00" in astra_script
    assert "--gres=gpu" in (data_dir / "slurm" / "06.sbatch").read_text()
    assert "--gres=gpu" not in (data_dir / "slurm" / "06i.sbatch").read_text()
    assert "--skip-index" in (data_dir / "slurm" / "06.sbatch").read_text()
    submit_text = submit.read_text()
    assert "bash" in submit_text
    assert "05c.local.sh" in submit_text
    assert "afterok:${jid_03}" in submit_text
    assert "--scheduler-managed" in (data_dir / "slurm" / "05c.local.sh").read_text()
    assert "--scheduler-managed" in astra_script


def test_slurm_submission_idempotency_suppresses_duplicate_scheduler_launch(
    tmp_path,
    monkeypatch,
):
    data_dir = tmp_path / "dataset"
    stages = StagePaths.from_root(data_dir)
    profile = resolve_resource_profile("slurm")
    dag = _build_tools_dag(
        input_dir=tmp_path / "inputs",
        data_dir=data_dir,
        output=data_dir / "sharur.duckdb",
        stages=stages,
        profile=profile,
        skip_quast=True,
        skip_dfast=True,
        skip_prodigal=False,
        skip_astra=False,
        skip_gecco=True,
        skip_dbcan=True,
        skip_crispr=False,
        skip_embeddings=False,
        enable_cazymes=False,
    )
    submissions = []
    monkeypatch.setattr(
        "sharur.ingest_cli.subprocess.run",
        lambda command, check: submissions.append((command, check)),
    )
    kwargs = {
        "dag": dag,
        "data_dir": data_dir,
        "profile": profile,
        "resume": True,
        "force": False,
        "idempotency_key": "stable-slurm-run",
        "submit": True,
    }

    first = _emit_slurm_bundle(**kwargs)
    second = _emit_slurm_bundle(**kwargs)

    assert first == second
    assert len(submissions) == 1
    ledger = RunLedger(data_dir / "sharur_ops.db")
    try:
        runs = ledger.list_runs(dataset_path=data_dir, run_type="ingest")
    finally:
        ledger.close()
    assert len(runs) == 1
    assert runs[0]["status"] == "submitted"


def test_slurm_submission_failure_is_terminally_recorded(tmp_path, monkeypatch):
    data_dir = tmp_path / "dataset"
    stages = StagePaths.from_root(data_dir)
    profile = resolve_resource_profile("slurm")
    dag = _build_tools_dag(
        input_dir=tmp_path / "inputs",
        data_dir=data_dir,
        output=data_dir / "sharur.duckdb",
        stages=stages,
        profile=profile,
        skip_quast=True,
        skip_dfast=True,
        skip_prodigal=False,
        skip_astra=False,
        skip_gecco=True,
        skip_dbcan=True,
        skip_crispr=False,
        skip_embeddings=True,
        enable_cazymes=False,
    )

    def fail_submission(command, check):
        raise subprocess.CalledProcessError(1, command)

    monkeypatch.setattr("sharur.ingest_cli.subprocess.run", fail_submission)

    with pytest.raises(subprocess.CalledProcessError):
        _emit_slurm_bundle(
            dag,
            data_dir=data_dir,
            profile=profile,
            resume=False,
            force=False,
            idempotency_key="failed-submit",
            submit=True,
        )

    ledger = RunLedger(data_dir / "sharur_ops.db")
    try:
        run = ledger.list_runs(dataset_path=data_dir, run_type="ingest")[0]
    finally:
        ledger.close()
    assert run["status"] == "failed"
    assert "scheduler submission failed" in run["error"]
