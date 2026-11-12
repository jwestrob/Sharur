#!/usr/bin/env python3
"""Block until Slurm jobs reach terminal states, then print their outcome.

An agent harness that notifies on background-process exit (Claude Code's
``run_in_background``) turns this into a completion watcher: start it right
after ``sbatch`` and the session hears back when the jobs end, fail, or wedge.

    python scripts/slurm_wait.py 21001 21002_[0-15] --host user@cluster.example.edu \
        --log '/path/to/logs/cazy_%j.out'

It polls ``sacct`` (over SSH when ``--host`` is given) every ``--interval``
seconds. Array jobs are tracked task by task. A pending job whose reason is
``DependencyNeverSatisfied`` counts as failed, since Slurm never runs it. The
exit status is 0 when every job COMPLETED, 1 when any ended otherwise, and 2
when the jobs could not be queried for ``--give-up`` seconds.
"""

from __future__ import annotations

import argparse
import re
import shlex
import subprocess
import sys
import time
from collections import Counter

TERMINAL_OK = {"COMPLETED"}
TERMINAL_BAD = {"FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY", "NODE_FAIL", "PREEMPTED",
                "BOOT_FAIL", "DEADLINE", "REVOKED"}


def run_remote(host: str | None, command: str) -> str:
    argv = ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=20", host, command] if host else \
        ["/bin/sh", "-c", command]
    return subprocess.run(argv, capture_output=True, text=True, check=True).stdout


def base_ids(jobs: list[str]) -> list[str]:
    return sorted({re.match(r"\d+", j).group(0) for j in jobs if re.match(r"\d+", j)})


def poll(host: str | None, jobs: list[str]) -> dict[str, tuple[str, str]]:
    """job/task id -> (state, reason) for every job and array task."""
    ids = ",".join(base_ids(jobs))
    out = run_remote(host, f"sacct -X -n -P -j {shlex.quote(ids)} --format=JobID,State,Reason")
    states: dict[str, tuple[str, str]] = {}
    for line in out.splitlines():
        parts = line.split("|")
        if len(parts) < 2 or not parts[0]:
            continue
        job_id, state = parts[0], parts[1].split()[0]
        reason = parts[2] if len(parts) > 2 else ""
        states[job_id] = (state.rstrip("+"), reason)
    # sacct reports a pending array as one row (``123_[0-15]``); squeue names the reason
    pending = [j for j, (s, _) in states.items() if s == "PENDING"]
    if pending:
        sq = run_remote(host, f"squeue -h -j {shlex.quote(ids)} -o '%i|%T|%r'")
        for line in sq.splitlines():
            parts = line.split("|")
            if len(parts) == 3:
                states[parts[0]] = (parts[1], parts[2])
    return states


def summarize(states: dict[str, tuple[str, str]]) -> str:
    counts = Counter(s for s, _ in states.values())
    return ", ".join(f"{k} {v}" for k, v in sorted(counts.items()))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("jobs", nargs="+", help="Job IDs (array jobs by their base ID)")
    ap.add_argument("--host", default=None, help="SSH target running Slurm (default: local)")
    ap.add_argument("--interval", type=int, default=120, help="Seconds between polls (default 120)")
    ap.add_argument("--give-up", type=int, default=3600,
                    help="Exit 2 after this many seconds of failed queries (default 3600)")
    ap.add_argument("--log", action="append", default=[],
                    help="Log path pattern to tail at the end; %%j is replaced by each job/task ID "
                         "and %%A/%%a by array job/task")
    ap.add_argument("--tail", type=int, default=15, help="Lines of each failed job's log to print")
    args = ap.parse_args()

    last_ok = time.time()
    last_line = ""
    while True:
        try:
            states = poll(args.host, args.jobs)
            last_ok = time.time()
        except subprocess.CalledProcessError as exc:
            if time.time() - last_ok > args.give_up:
                print(f"giving up: sacct unreachable for {args.give_up}s: {exc.stderr.strip()[:300]}")
                return 2
            time.sleep(args.interval)
            continue
        wedged = {j for j, (s, r) in states.items() if s == "PENDING" and "DependencyNeverSatisfied" in r}
        live = {j for j, (s, _) in states.items()
                if s not in TERMINAL_OK | TERMINAL_BAD and j not in wedged}
        line = summarize(states)
        if line != last_line:
            print(f"{time.strftime('%H:%M:%S')} {line}", flush=True)
            last_line = line
        if states and not live:
            break
        time.sleep(args.interval)

    bad = {j: s for j, (s, _) in states.items() if s not in TERMINAL_OK}
    print(f"FINAL {summarize(states)}" + (f"; DependencyNeverSatisfied: {sorted(wedged)}" if wedged else ""))
    def order(job_id: str) -> list[tuple[int, int | str]]:
        return [(0, int(x)) if x.isdigit() else (1, x) for x in re.split(r"[_\[\]]", job_id) if x]

    for job_id in sorted(bad, key=order):
        print(f"  {job_id}: {bad[job_id]}")
        if job_id in wedged:
            continue
        array, _, task = job_id.partition("_")
        for pattern in args.log:
            path = pattern.replace("%j", job_id).replace("%A", array).replace("%a", task)
            try:
                text = run_remote(args.host, f"tail -n {args.tail} {shlex.quote(path)}")
            except subprocess.CalledProcessError:
                continue
            print("    " + text.strip().replace("\n", "\n    "))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
