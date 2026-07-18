"""Execute one ledger-tracked ingest stage (local or inside a SLURM job)."""

from __future__ import annotations

import argparse
import contextlib
import json
import shutil
import subprocess
from pathlib import Path

from sharur.ingest.dag import snapshot_path
from sharur.ingest.resources import ResourceRequest, accelerator_lock
from sharur.ops.ledger import RunLedger, StageAlreadyCompleteError


def _clean(path: Path) -> None:
    if path.is_dir() and not path.is_symlink():
        shutil.rmtree(path)
    else:
        path.unlink(missing_ok=True)


def _validate_outputs(outputs: list[Path]) -> None:
    missing = [str(path) for path in outputs if not path.exists()]
    empty = [str(path) for path in outputs if path.is_file() and path.stat().st_size == 0]
    if missing or empty:
        raise RuntimeError(f"Stage output validation failed; missing={missing}, empty={empty}")


def _validate_output_snapshots(
    outputs: list[Path],
    expected: dict[str, object],
) -> None:
    _validate_outputs(outputs)
    current = {str(path): snapshot_path(path) for path in outputs}
    mismatched = [path for path, snapshot in current.items() if expected.get(path) != snapshot]
    if mismatched:
        raise RuntimeError(
            "Completed stage outputs changed before a duplicate launch: " + ", ".join(mismatched)
        )


def execute_stage(
    *,
    ledger: RunLedger,
    run_id: str,
    stage_id: str,
    signature: str,
    command: list[str],
    outputs: list[Path],
    cleanup: list[Path],
    resource: ResourceRequest,
    complete_run: bool = False,
    scheduler_managed: bool = False,
) -> int:
    """Run a subprocess with periodic durable heartbeats."""
    try:
        attempt = ledger.start_stage(
            run_id,
            stage_id,
            signature,
            command=command,
            outputs={str(path): {"state": "planned"} for path in outputs},
            resource_profile=resource.to_dict(),
        )
    except StageAlreadyCompleteError as exc:
        _validate_output_snapshots(outputs, exc.stage.get("outputs") or {})
        ledger.append_event(
            run_id,
            "duplicate_stage_launch_suppressed",
            stage_id=stage_id,
            attempt=exc.stage["attempt"],
            payload={"signature": signature},
        )
        run_status = ledger.get_run(run_id)["status"]
        if complete_run and run_status != "complete":
            ledger.complete_run(
                run_id,
                result={"last_stage": stage_id, "outputs": exc.stage.get("outputs") or {}},
            )
        elif scheduler_managed and run_status in {"submitted", "running"}:
            ledger.wait_for_scheduler(run_id)
        return 0

    try:
        for path in cleanup:
            _clean(path)
        with accelerator_lock(
            resource,
            on_wait=lambda: ledger.heartbeat_stage(run_id, stage_id, attempt),
        ):
            process = subprocess.Popen(command)
            while True:
                try:
                    return_code = process.wait(timeout=30)
                    break
                except subprocess.TimeoutExpired:
                    ledger.heartbeat_stage(run_id, stage_id, attempt)
        if return_code != 0:
            raise subprocess.CalledProcessError(return_code, command)

        _validate_outputs(outputs)
        result = {str(path): snapshot_path(path) for path in outputs}
        ledger.complete_stage(
            run_id,
            stage_id,
            attempt,
            outputs=result,
        )
        if complete_run:
            ledger.complete_run(
                run_id,
                result={"last_stage": stage_id, "outputs": result},
            )
        elif scheduler_managed:
            ledger.wait_for_scheduler(run_id)
        return 0
    except Exception as exc:
        message = f"{type(exc).__name__}: {exc}"
        ledger.fail_stage(run_id, stage_id, attempt, message)
        with contextlib.suppress(ValueError):
            ledger.fail_run(run_id, message)
        raise


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", required=True, type=Path)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--stage-id", required=True)
    parser.add_argument("--signature", required=True)
    parser.add_argument("--resource-json", required=True)
    parser.add_argument("--output", action="append", default=[], type=Path)
    parser.add_argument("--cleanup", action="append", default=[], type=Path)
    parser.add_argument("--complete-run", action="store_true")
    parser.add_argument("--scheduler-managed", action="store_true")
    parser.add_argument("command", nargs=argparse.REMAINDER)
    return parser


def main() -> int:
    args = _parser().parse_args()
    command = list(args.command)
    if command and command[0] == "--":
        command = command[1:]
    if not command:
        raise SystemExit("A command is required after --")
    resource = ResourceRequest(**json.loads(args.resource_json))
    ledger = RunLedger(args.ledger)
    try:
        return execute_stage(
            ledger=ledger,
            run_id=args.run_id,
            stage_id=args.stage_id,
            signature=args.signature,
            command=command,
            outputs=args.output,
            cleanup=args.cleanup,
            resource=resource,
            complete_run=args.complete_run,
            scheduler_managed=args.scheduler_managed,
        )
    finally:
        ledger.close()


if __name__ == "__main__":
    raise SystemExit(main())
