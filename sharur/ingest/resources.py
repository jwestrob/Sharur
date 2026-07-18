"""Explicit local CPU, Apple MPS, and SLURM ingest resource profiles."""

from __future__ import annotations

import contextlib
import os
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path
from typing import TYPE_CHECKING


if TYPE_CHECKING:
    from collections.abc import Callable, Iterator


@dataclass(frozen=True)
class ResourceRequest:
    cpus: int
    memory_gb: int
    walltime: str
    executor: str = "local"  # local | slurm
    accelerator: str = "cpu"  # cpu | mps | cuda
    gpus: int = 0
    exclusive_accelerator: bool = False

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass(frozen=True)
class ResourceProfile:
    name: str
    requests: dict[str, ResourceRequest]
    max_workers: int
    annotation_threads: int
    index_threads: int

    def request(self, stage_id: str) -> ResourceRequest:
        return self.requests[stage_id]

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "max_workers": self.max_workers,
            "annotation_threads": self.annotation_threads,
            "index_threads": self.index_threads,
            "requests": {stage: request.to_dict() for stage, request in self.requests.items()},
        }


_STAGES = (
    "00",
    "01",
    "02",
    "03",
    "04",
    "05a",
    "05b",
    "05c",
    "07",
    "06",
    "06i",
)


@lru_cache(maxsize=1)
def _mps_available() -> bool:
    """Probe MPS outside this process to avoid Torch/FAISS native-runtime conflicts."""
    if sys.platform != "darwin":
        return False
    try:
        process = subprocess.run(
            [
                sys.executable,
                "-c",
                (
                    "import torch; "
                    "raise SystemExit(0 if "
                    "torch.backends.mps.is_built() and "
                    "torch.backends.mps.is_available() else 1)"
                ),
            ],
            check=False,
            capture_output=True,
            timeout=15,
        )
    except (OSError, subprocess.SubprocessError):
        return False
    return process.returncode == 0


def _local_profile(*, mps: bool = False) -> ResourceProfile:
    cpus = max(1, os.cpu_count() or 1)
    workers = min(cpus, 16)
    annotation_threads = min(cpus, 16)
    common = {
        stage: ResourceRequest(
            cpus=1,
            memory_gb=4,
            walltime="04:00:00",
        )
        for stage in _STAGES
    }
    common.update(
        {
            "01": ResourceRequest(workers, 16, "12:00:00"),
            "02": ResourceRequest(workers, 24, "24:00:00"),
            "03": ResourceRequest(workers, 16, "12:00:00"),
            "04": ResourceRequest(annotation_threads, 48, "72:00:00"),
            "05a": ResourceRequest(min(workers, 8), 24, "24:00:00"),
            "05b": ResourceRequest(min(workers, 8), 24, "24:00:00"),
            # MinCED is intentionally one CPU and local.
            "05c": ResourceRequest(1, 4, "24:00:00"),
            "07": ResourceRequest(min(cpus, 8), 32, "24:00:00"),
            "06": ResourceRequest(
                cpus=min(cpus, 8),
                memory_gb=48,
                walltime="72:00:00",
                accelerator="mps" if mps else "cpu",
                gpus=1 if mps else 0,
                exclusive_accelerator=mps,
            ),
            "06i": ResourceRequest(
                cpus=min(cpus, 8),
                memory_gb=48,
                walltime="24:00:00",
            ),
        }
    )
    return ResourceProfile(
        name="mps" if mps else "local",
        requests=common,
        max_workers=workers,
        annotation_threads=annotation_threads,
        index_threads=min(cpus, 8),
    )


def _slurm_profile() -> ResourceProfile:
    requests = {
        "00": ResourceRequest(1, 4, "04:00:00", executor="local"),
        "01": ResourceRequest(8, 24, "12:00:00", executor="slurm"),
        "02": ResourceRequest(16, 48, "24:00:00", executor="slurm"),
        "03": ResourceRequest(16, 32, "24:00:00", executor="slurm"),
        # KOFAM can dominate this stage; use a deliberately generous walltime.
        "04": ResourceRequest(16, 64, "72:00:00", executor="slurm"),
        "05a": ResourceRequest(8, 32, "24:00:00", executor="slurm"),
        "05b": ResourceRequest(8, 32, "24:00:00", executor="slurm"),
        # Project policy: single-threaded MinCED runs on the login node.
        "05c": ResourceRequest(1, 8, "24:00:00", executor="local"),
        "07": ResourceRequest(16, 64, "24:00:00", executor="slurm"),
        "06": ResourceRequest(
            8,
            64,
            "72:00:00",
            executor="slurm",
            accelerator="cuda",
            gpus=1,
            exclusive_accelerator=True,
        ),
        "06i": ResourceRequest(
            8,
            64,
            "24:00:00",
            executor="slurm",
        ),
    }
    return ResourceProfile(
        name="slurm",
        requests=requests,
        max_workers=16,
        annotation_threads=16,
        index_threads=8,
    )


def resolve_resource_profile(name: str) -> ResourceProfile:
    normalized = name.lower().strip()
    if normalized == "auto":
        normalized = "mps" if _mps_available() else "local"
    if normalized == "local":
        return _local_profile(mps=False)
    if normalized == "mps":
        if not _mps_available():
            raise ValueError("The mps resource profile requires a usable PyTorch MPS backend")
        return _local_profile(mps=True)
    if normalized == "slurm":
        return _slurm_profile()
    raise ValueError("profile must be one of: auto, local, mps, slurm")


@contextlib.contextmanager
def accelerator_lock(
    request: ResourceRequest,
    lock_path: Path | None = None,
    on_wait: Callable[[], None] | None = None,
) -> Iterator[None]:
    """Enforce one local MPS workload across Sharur processes."""
    if request.accelerator != "mps" or not request.exclusive_accelerator:
        yield
        return

    import fcntl  # noqa: PLC0415 - unavailable on non-POSIX platforms

    path = lock_path or Path.home() / ".cache" / "sharur" / "mps.lock"
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "a+b") as handle:
        while True:
            try:
                fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except BlockingIOError:
                if on_wait is not None:
                    on_wait()
                time.sleep(5)
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


__all__ = [
    "ResourceProfile",
    "ResourceRequest",
    "accelerator_lock",
    "resolve_resource_profile",
]
