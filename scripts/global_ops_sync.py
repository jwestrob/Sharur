#!/usr/bin/env python3
"""Sync findings from per-dataset OpsStores into the global mirror.

Use:
    # Sync one dataset
    python scripts/global_ops_sync.py --dataset data/coronamine_boiler_100nm

    # Sync every dataset under data/
    python scripts/global_ops_sync.py --all

    # Backfill from findings.jsonl when no per-dataset OpsStore exists
    python scripts/global_ops_sync.py --dataset data/spicy_lams --backfill-jsonl

    # Show what's in the global mirror
    python scripts/global_ops_sync.py --status
    python scripts/global_ops_sync.py --query-accession DUF6088
    python scripts/global_ops_sync.py --co-occurring DUF6088
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from sharur.ops.global_store import (
    GlobalOpsStore,
    mirror_from_local,
    mirror_from_findings_jsonl,
    DEFAULT_GLOBAL_DB,
)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--dataset", help="One dataset directory to sync (data/<name>).")
    p.add_argument("--all", action="store_true",
                   help="Sync every dataset under data/ that has a sharur_ops.db.")
    p.add_argument("--data-root", default="data",
                   help="Root containing dataset subdirs. Default: data/")
    p.add_argument("--min-novelty", type=int, default=2,
                   help="Minimum novelty to mirror. Default: 2.")
    p.add_argument("--global-db", default=str(DEFAULT_GLOBAL_DB),
                   help=f"Path to global mirror DB. Default: {DEFAULT_GLOBAL_DB}")
    p.add_argument("--backfill-jsonl", action="store_true",
                   help="If no sharur_ops.db, backfill from findings.jsonl files.")
    p.add_argument("--status", action="store_true",
                   help="Print global mirror status and exit.")
    p.add_argument("--query-accession", metavar="ACC",
                   help="Show all findings mentioning ACC across datasets.")
    p.add_argument("--co-occurring", metavar="ACC",
                   help="Show accessions that co-occur with ACC in 2+ datasets.")
    p.add_argument("--recent", action="store_true",
                   help="Show the 20 most recent high-novelty findings.")
    args = p.parse_args()

    g = GlobalOpsStore(args.global_db)

    if args.status:
        datasets = g.datasets()
        if not datasets:
            print("Global mirror is empty. Run with --dataset <path> or --all to populate.")
            return 0
        print(f"Global mirror: {args.global_db}")
        print()
        print(f"{'dataset':<40} {'findings':>10} {'last sync':>20}")
        for d in datasets:
            ts = d["last_synced_ts"] or 0
            import datetime
            sync_str = datetime.datetime.fromtimestamp(ts).strftime("%Y-%m-%d %H:%M") if ts else "(never)"
            print(f"{d['dataset_name']:<40} {d['n_findings']:>10} {sync_str:>20}")
        return 0

    if args.query_accession:
        hits = g.find_by_accession(args.query_accession)
        if not hits:
            print(f"No findings mention {args.query_accession} in the global mirror.")
            return 0
        print(f"Findings mentioning {args.query_accession}:")
        seen_datasets = set()
        for h in hits:
            ds_name = Path(h["source_dataset"]).name
            seen_datasets.add(ds_name)
            print(f"  [{h['novelty']}] {ds_name}/{h['finding_id']}: {h['summary'][:80]}")
        print(f"\n  Across {len(seen_datasets)} dataset(s): {', '.join(sorted(seen_datasets))}")
        return 0

    if args.co_occurring:
        hits = g.co_occurring_accessions(args.co_occurring)
        if not hits:
            print(f"No accessions co-occur with {args.co_occurring} in 2+ datasets.")
            return 0
        print(f"Accessions co-occurring with {args.co_occurring}:")
        print(f"  {'accession':<15} {'#datasets':>10} {'#findings':>10}")
        for h in hits[:30]:
            print(f"  {h['accession']:<15} {h['n_datasets']:>10} {h['n_findings']:>10}")
        return 0

    if args.recent:
        for f in g.recent(min_novelty=args.min_novelty, limit=20):
            ds_name = Path(f["source_dataset"]).name
            print(f"  [{f['novelty']}] {ds_name}/{f['finding_id']}: {f['summary'][:80]}")
        return 0

    # Sync paths
    targets: list[Path] = []
    if args.all:
        root = Path(args.data_root)
        if not root.exists():
            print(f"ERROR: {root} not found", file=sys.stderr)
            return 2
        for ds in sorted(root.iterdir()):
            if not ds.is_dir():
                continue
            if (ds / "sharur_ops.db").exists() or args.backfill_jsonl:
                targets.append(ds)
    elif args.dataset:
        targets.append(Path(args.dataset))
    else:
        p.print_help()
        return 2

    total = 0
    for ds in targets:
        ds = ds.resolve()
        ops_db = ds / "sharur_ops.db"
        n = 0
        if ops_db.exists():
            n = mirror_from_local(g, ops_db, min_novelty=args.min_novelty, dataset_root=ds)
            print(f"  {ds.name}: synced {n} findings from sharur_ops.db")
        elif args.backfill_jsonl:
            for phase in ("survey", "exploration"):
                jp = ds / phase / "findings.jsonl"
                if jp.exists():
                    sub = mirror_from_findings_jsonl(g, jp, dataset_root=ds, min_novelty=args.min_novelty)
                    n += sub
                    print(f"  {ds.name}/{phase}: backfilled {sub} findings from jsonl")
        else:
            print(f"  {ds.name}: no sharur_ops.db; skipped (--backfill-jsonl to use findings.jsonl)")
        total += n

    print(f"\nDone. Mirrored {total} findings into {args.global_db}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
