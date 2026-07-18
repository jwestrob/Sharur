"""Ingest DAG, signatures, output validation, and SLURM bundle rendering."""

from __future__ import annotations

import hashlib
import json
import os
import shlex
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any


if TYPE_CHECKING:
    from collections.abc import Iterable

    from sharur.ingest.resources import ResourceProfile, ResourceRequest


SIGNATURE_VERSION = 1


@dataclass
class StageNode:
    stage_id: str
    label: str
    command: list[str]
    dependencies: tuple[str, ...] = ()
    required_outputs: tuple[Path, ...] = ()
    input_paths: tuple[Path, ...] = ()
    cleanup_paths: tuple[Path, ...] = ()
    resource: ResourceRequest | None = None
    signature: str = field(default="", init=False)

    def output_snapshot(self) -> dict[str, Any]:
        return {str(path): snapshot_path(path) for path in self.required_outputs}

    def outputs_match(self, recorded: dict[str, Any] | None) -> bool:
        """Require live outputs to match the successful attempt's snapshots."""
        if not self.outputs_ready() or not recorded:
            return False
        current = self.output_snapshot()
        return all(recorded.get(path) == snapshot for path, snapshot in current.items())

    def outputs_ready(self) -> bool:
        for path in self.required_outputs:
            if not path.exists():
                return False
            if path.is_file() and path.stat().st_size == 0:
                return False
        return True


class IngestDAG:
    """Validated stage graph with deterministic topological ordering."""

    def __init__(self, nodes: Iterable[StageNode]):
        node_list = list(nodes)
        self.nodes = {node.stage_id: node for node in node_list}
        if len(self.nodes) != len(node_list):
            raise ValueError("Stage IDs must be unique")
        for node in node_list:
            missing = set(node.dependencies) - self.nodes.keys()
            if missing:
                raise ValueError(
                    f"Stage {node.stage_id} has missing dependencies: {sorted(missing)}"
                )
        self._ordered = self._topological_order()

    def _topological_order(self) -> list[StageNode]:
        remaining = set(self.nodes)
        complete: set[str] = set()
        ordered: list[StageNode] = []
        while remaining:
            ready = sorted(
                stage_id
                for stage_id in remaining
                if set(self.nodes[stage_id].dependencies) <= complete
            )
            if not ready:
                raise ValueError("Ingest stage graph contains a dependency cycle")
            for stage_id in ready:
                ordered.append(self.nodes[stage_id])
                complete.add(stage_id)
                remaining.remove(stage_id)
        return ordered

    def ordered(self) -> list[StageNode]:
        return list(self._ordered)

    def compute_signatures(self, profile: ResourceProfile) -> dict[str, str]:
        signatures: dict[str, str] = {}
        for node in self._ordered:
            payload = {
                "version": SIGNATURE_VERSION,
                "stage_id": node.stage_id,
                "command": node.command,
                "dependencies": {
                    dependency: signatures[dependency] for dependency in node.dependencies
                },
                "inputs": {str(path): snapshot_path(path) for path in node.input_paths},
                "script": snapshot_path(Path(node.command[0])),
                "resource": (node.resource or profile.request(node.stage_id)).to_dict(),
            }
            encoded = json.dumps(
                payload,
                sort_keys=True,
                separators=(",", ":"),
                default=str,
            ).encode("utf-8")
            node.signature = hashlib.sha256(encoded).hexdigest()
            signatures[node.stage_id] = node.signature
        return signatures


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def snapshot_path(path: Path) -> dict[str, Any]:
    """Deterministic artifact snapshot used for safe resume decisions."""
    resolved = path.expanduser().resolve()
    if not resolved.exists():
        return {"state": "missing", "path": str(resolved)}
    stat = resolved.stat()
    if resolved.is_file():
        result: dict[str, Any] = {
            "state": "file",
            "path": str(resolved),
            "size": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
        }
        if stat.st_size <= 10 * 1024 * 1024 or resolved.suffix in {
            ".json",
            ".py",
            ".toml",
            ".yaml",
            ".yml",
        }:
            result["sha256"] = _file_sha256(resolved)
        return result

    digest = hashlib.sha256()
    entry_count = 0
    file_count = 0
    directory_count = 0
    symlink_count = 0
    total_file_size = 0

    def walk(directory: Path) -> None:
        nonlocal entry_count, file_count, directory_count, symlink_count
        nonlocal total_file_size
        with os.scandir(directory) as iterator:
            entries = sorted(iterator, key=lambda item: item.name)
        for entry in entries:
            entry_path = Path(entry.path)
            relative = str(entry_path.relative_to(resolved))
            record: dict[str, Any] = {"path": relative}
            try:
                if entry.is_symlink():
                    symlink_count += 1
                    link_stat = entry.stat(follow_symlinks=False)
                    record.update(
                        {
                            "kind": "symlink",
                            "link_target": os.readlink(entry.path),
                            "link_mtime_ns": link_stat.st_mtime_ns,
                        }
                    )
                    try:
                        target_stat = entry.stat(follow_symlinks=True)
                        record.update(
                            {
                                "target_size": target_stat.st_size,
                                "target_mtime_ns": target_stat.st_mtime_ns,
                            }
                        )
                    except OSError as exc:
                        record["target_error"] = f"{type(exc).__name__}: {exc}"
                elif entry.is_dir(follow_symlinks=False):
                    directory_count += 1
                    child_stat = entry.stat(follow_symlinks=False)
                    record.update(
                        {
                            "kind": "directory",
                            "mtime_ns": child_stat.st_mtime_ns,
                        }
                    )
                else:
                    file_count += 1
                    child_stat = entry.stat(follow_symlinks=False)
                    total_file_size += child_stat.st_size
                    record.update(
                        {
                            "kind": "file",
                            "size": child_stat.st_size,
                            "mtime_ns": child_stat.st_mtime_ns,
                        }
                    )
            except OSError as exc:
                record["error"] = f"{type(exc).__name__}: {exc}"
            digest.update(
                json.dumps(
                    record,
                    sort_keys=True,
                    separators=(",", ":"),
                ).encode("utf-8")
            )
            entry_count += 1
            if record.get("kind") == "directory":
                walk(entry_path)

    walk(resolved)
    return {
        "state": "directory",
        "path": str(resolved),
        "mtime_ns": stat.st_mtime_ns,
        "entry_count": entry_count,
        "file_count": file_count,
        "directory_count": directory_count,
        "symlink_count": symlink_count,
        "total_file_size": total_file_size,
        "tree_sha256": digest.hexdigest(),
    }


def stage_runner_command(
    node: StageNode,
    *,
    ledger_path: Path,
    run_id: str,
    complete_run: bool,
) -> list[str]:
    request = node.resource
    if request is None:
        raise ValueError(f"Stage {node.stage_id} has no resource request")
    command = [
        sys.executable,
        "-m",
        "sharur.ingest.stage_runner",
        "--ledger",
        str(ledger_path),
        "--run-id",
        run_id,
        "--stage-id",
        node.stage_id,
        "--signature",
        node.signature,
        "--resource-json",
        json.dumps(request.to_dict(), separators=(",", ":")),
    ]
    for output in node.required_outputs:
        command.extend(["--output", str(output)])
    for cleanup in node.cleanup_paths:
        command.extend(["--cleanup", str(cleanup)])
    if complete_run:
        command.append("--complete-run")
    command.append("--scheduler-managed")
    command.extend(["--", sys.executable, *node.command])
    return command


def render_slurm_script(
    node: StageNode,
    runner_command: list[str],
    *,
    log_dir: Path,
) -> str:
    request = node.resource
    if request is None:
        raise ValueError(f"Stage {node.stage_id} has no resource request")
    lines = [
        "#!/usr/bin/env bash",
        f"#SBATCH --job-name=sharur-{node.stage_id}",
        f"#SBATCH --cpus-per-task={request.cpus}",
        f"#SBATCH --mem={request.memory_gb}G",
        f"#SBATCH --time={request.walltime}",
        f"#SBATCH --output={log_dir / (node.stage_id + '-%j.log')}",
    ]
    if request.gpus:
        lines.append(f"#SBATCH --gres=gpu:{request.gpus}")
    lines.extend(["set -euo pipefail", shlex.join(runner_command), ""])
    return "\n".join(lines)


__all__ = [
    "IngestDAG",
    "StageNode",
    "render_slurm_script",
    "snapshot_path",
    "stage_runner_command",
]
