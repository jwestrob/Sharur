#!/usr/bin/env python3
"""Compare Sharur's co-location engine against the reference implementation.

The engine reimplements MacSyFinder so that system-level validation runs in
seconds against DuckDB instead of minutes-to-hours in a subprocess. That is only
worth having if the two agree, and agreement has to be measured rather than
assumed -- especially for model sets that use multiple loci, loner genes or
per-gene spacing overrides, where the reference behaviour is subtle.

Identity is compared on ``(genome_id, system_type, frozenset(protein_ids))``.
Counts alone are not sufficient: two implementations can produce the same number
of systems while disagreeing about which proteins are in them, and that error is
invisible to a totals check.

Both sides must be derived from the SAME hits. If the engine's input annotations
were produced by a different search than the reference run, a mismatch is
uninformative -- it measures search settings, not co-location logic.

Usage:
    # secretion systems, against MacSyFinder per-genome output
    check_colocation_concordance.py DATASET --source txsscan \
        --reference-dir DATASET/txsscan/results

    # defense systems, against the table a real defense-finder run produced
    check_colocation_concordance.py DATASET --source defensefinder \
        --reference-table defense_systems
"""

from __future__ import annotations

import argparse
import glob
import json
import logging
import sys
from collections import defaultdict
from pathlib import Path

import duckdb
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

System = tuple[str, str, frozenset]


def reference_from_macsyfinder(results_dir: Path) -> set[System]:
    """Systems as reported in per-genome ``best_solution.tsv`` files."""
    out: set[System] = set()
    for path in glob.glob(str(results_dir / "*" / "best_solution.tsv")):
        try:
            df = pd.read_csv(path, sep="\t", comment="#")
        except Exception:  # noqa: BLE001 - an empty file means "no systems"
            continue
        if df.empty or "sys_id" not in df.columns:
            continue
        genome = Path(path).parent.name
        for _sid, grp in df.groupby("sys_id"):
            model = str(grp["model_fqn"].iloc[0]).split("/")[-1]
            out.add((genome, model, frozenset(str(h) for h in grp["hit_id"])))
    return out


def reference_from_table(db_path: Path, table: str) -> set[System]:
    """Systems as stored by a previous reference run."""
    con = duckdb.connect(str(db_path), read_only=True, config={"threads": 4})
    rows = con.execute(f"SELECT genome_id, system_type, protein_ids FROM {table}").fetchall()
    con.close()
    return {
        (str(g), str(t), frozenset(str(p).split(",")))
        for g, t, p in rows
    }


def engine_systems(db_path: Path, source: str) -> set[System]:
    from sharur.colocation import validate_systems

    systems_df, genes_df = validate_systems(db_path=str(db_path), source=source, verbose=False)
    if systems_df.empty:
        return set()
    by_system = genes_df.groupby("system_id")["protein_id"].apply(
        lambda x: frozenset(map(str, x))
    )
    return {
        (str(r["genome_id"]), str(r["type"]), by_system.get(r["system_id"], frozenset()))
        for _, r in systems_df.iterrows()
    }


def report(ref: set[System], ours: set[System]) -> dict:
    matched = ref & ours
    ref_only = ref - ours
    ours_only = ours - ref

    # A system present in both with a differing protein set shows up once in each
    # difference set. Pairing them by (genome, type) separates "we disagree about
    # membership" from "one side missed it entirely" -- different bugs.
    key = lambda s: (s[0], s[1])  # noqa: E731
    ref_idx = defaultdict(list)
    for s in ref_only:
        ref_idx[key(s)].append(s)
    membership, missing, extra = [], [], []
    for s in ours_only:
        if ref_idx.get(key(s)):
            membership.append((ref_idx[key(s)].pop(0), s))
        else:
            extra.append(s)
    missing = [s for group in ref_idx.values() for s in group]

    pct = len(matched) / len(ref) * 100 if ref else 100.0
    return {
        "reference": len(ref),
        "engine": len(ours),
        "exact_matches": len(matched),
        "concordance_pct": round(pct, 4),
        "membership_disagreements": len(membership),
        "missing_from_engine": len(missing),
        "extra_in_engine": len(extra),
        "_membership": membership,
        "_missing": missing,
        "_extra": extra,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dataset", type=Path)
    ap.add_argument("--source", required=True, help="annotation source the engine validates")
    ap.add_argument("--reference-dir", type=Path, default=None, help="MacSyFinder per-genome results root")
    ap.add_argument("--reference-table", default=None, help="DuckDB table holding a prior reference run")
    ap.add_argument("--json-out", type=Path, default=None)
    ap.add_argument("--fail-under", type=float, default=None, help="Exit non-zero below this concordance %%")
    args = ap.parse_args()

    db_path = args.dataset / "sharur.duckdb"
    if args.reference_dir:
        ref = reference_from_macsyfinder(args.reference_dir)
    elif args.reference_table:
        ref = reference_from_table(db_path, args.reference_table)
    else:
        logger.error("need --reference-dir or --reference-table")
        return 2

    ours = engine_systems(db_path, args.source)
    res = report(ref, ours)

    print(f"\n  source                   {args.source}")
    print(f"  reference systems        {res['reference']}")
    print(f"  engine systems           {res['engine']}")
    print(f"  exact matches            {res['exact_matches']}  ({res['concordance_pct']}%)")
    print(f"  membership disagreements {res['membership_disagreements']}")
    print(f"  missing from engine      {res['missing_from_engine']}")
    print(f"  extra in engine          {res['extra_in_engine']}")
    for label, items in (("membership", res["_membership"]),):
        for a, b in items[:10]:
            only_ref = sorted(a[2] - b[2])
            only_our = sorted(b[2] - a[2])
            print(f"    {label}: {a[1]} {a[0]}  ref-only={only_ref}  engine-only={only_our}")
    for label, items in (("missing", res["_missing"]), ("extra", res["_extra"])):
        for s in items[:10]:
            print(f"    {label}: {s[1]:<14}{s[0]:<22}{len(s[2])} genes")

    if args.json_out:
        args.json_out.write_text(json.dumps({k: v for k, v in res.items() if not k.startswith("_")}, indent=2))

    if args.fail_under is not None and res["concordance_pct"] < args.fail_under:
        logger.error("concordance %.4f%% below threshold %.4f%%", res["concordance_pct"], args.fail_under)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
