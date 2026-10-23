#!/usr/bin/env python3
"""Hydrogenase neighborhood verdict counts under the current rule and an earlier revision's rule.

The earlier rule is ``sharur/`` at ``--old-rev`` (``git archive``: its neighborhood module and the
mapping files it reads), run in a subprocess; the current rule is this checkout's. Every database opens
read-only. Writes one TSV row per dataset, call set and verdict to stdout; transitions and the markers
behind each lost ``supported`` verdict go to stderr.

Usage:
    python scripts/compare_hydrogenase_neighborhood.py --old-rev REV \\
        --db data/A/sharur.duckdb --db data/B/sharur.duckdb > IMPACT_hydrogenase.tsv
"""

from __future__ import annotations

import argparse
import io
import json
import os
import subprocess
import sys
import tarfile
import tempfile
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CALL_SETS = (
    ("all classified", "TRUE"),
    ("assigned, needs curation", "outcome = 'assigned' AND curation_status = 'needs_curation'"),
    ("assigned Group 4", "outcome = 'assigned' AND reference_class = 'NiFe' AND reference_subgroup LIKE 'Group_4%'"),
    ("assigned Group 4, needs curation", "outcome = 'assigned' AND reference_class = 'NiFe' "
                                         "AND reference_subgroup LIKE 'Group_4%' AND curation_status = 'needs_curation'"),
)


class ReadOnlyStore:
    """The ``execute`` surface ``neighborhood_contexts`` needs, on a read-only DuckDB connection."""

    def __init__(self, path: str, threads: int, memory: str):
        import duckdb  # noqa: PLC0415

        self.conn = duckdb.connect(path, read_only=True, config={"threads": threads, "memory_limit": memory})

    def execute(self, sql, params=None):
        return self.conn.execute(sql, params or []).fetchall()


def verdicts(db: str, threads: int, memory: str) -> dict[str, list]:
    """protein_id -> [verdict, hydrogenase-side markers] under the ``sharur`` first on sys.path."""
    from sharur.hydrogenase.neighborhood import neighborhood_contexts  # noqa: PLC0415

    store = ReadOnlyStore(db, threads, memory)
    ids = [pid for (pid,) in store.execute("SELECT protein_id FROM hydrogenase_classifications")]
    return {pid: [ctx.verdict, sorted(ctx.hydrogenase)] for pid, ctx in neighborhood_contexts(store, ids).items()}


def old_verdicts(rev: str, db: str, threads: int, memory: str) -> dict[str, list]:
    with tempfile.TemporaryDirectory() as tree:
        archive = subprocess.run(["git", "-C", str(ROOT), "archive", rev, "sharur"], capture_output=True, check=True)
        with tarfile.open(fileobj=io.BytesIO(archive.stdout)) as tar:
            tar.extractall(tree, filter="data")
        env = os.environ | {"SHARUR_RULE_TREE": tree}
        out = subprocess.run([sys.executable, __file__, "--verdicts-only", "--db", db, "--threads", str(threads),
                              "--memory", memory], capture_output=True, check=True, text=True, env=env, cwd=tree)
    return json.loads(out.stdout)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--db", action="append", required=True)
    parser.add_argument("--old-rev")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--memory", default="4GB")
    parser.add_argument("--verdicts-only", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()

    if args.verdicts_only:
        sys.path.insert(0, os.environ.get("SHARUR_RULE_TREE", str(ROOT)))
        json.dump(verdicts(args.db[0], args.threads, args.memory), sys.stdout)
        return 0
    if not args.old_rev:
        parser.error("--old-rev is required")

    sys.path.insert(0, str(ROOT))
    from sharur.hydrogenase.neighborhood import VERDICTS  # noqa: PLC0415

    print("dataset\tcalls\tverdict\told_rule\tnew_rule\tchange")
    for db in args.db:
        name = Path(db).parent.name
        old = old_verdicts(args.old_rev, db, args.threads, args.memory)
        new = verdicts(db, args.threads, args.memory)
        store = ReadOnlyStore(db, args.threads, args.memory)
        for label, where in CALL_SETS:
            rows = store.execute(f"SELECT protein_id, curation_status, outcome FROM hydrogenase_classifications "
                                 f"WHERE {where}")
            ids = [r[0] for r in rows]
            before, after = Counter(old[p][0] for p in ids), Counter(new[p][0] for p in ids)
            for verdict in list(VERDICTS) + sorted(set(before) - set(VERDICTS)):
                if before[verdict] or after[verdict]:
                    print(f"{name}\t{label}\t{verdict}\t{before[verdict]}\t{after[verdict]}\t"
                          f"{after[verdict] - before[verdict]:+d}")
            if label.startswith("assigned"):
                cleared = {p for p, status, _ in rows if status == "domain_check_cleared"}
                credible = [sum(p in cleared or v[p][0] == "supported" for p in ids) for v in (old, new)]
                print(f"{name}\t{label}\tcredible (domain cleared or supported)\t{credible[0]}\t{credible[1]}\t"
                      f"{credible[1] - credible[0]:+d}")
            moves = Counter((old[p][0], new[p][0]) for p in ids if old[p][0] != new[p][0])
            for (a, b), n in sorted(moves.items(), key=lambda kv: -kv[1]):
                print(f"# {name} | {label} | {a} -> {b}: {n}", file=sys.stderr)
            if label == "all classified":
                lost = [p for p in ids if old[p][0] == "supported" and new[p][0] != "supported"]
                markers = Counter(m for p in lost for m in old[p][1])
                print(f"# {name} | {len(lost)} calls lose 'supported'; their old hydrogenase-side markers: "
                      + ", ".join(f"{m} {n}" for m, n in markers.most_common(15)), file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
