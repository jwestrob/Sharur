#!/usr/bin/env python3
"""Run TXSScan (MacSyFinder) secretion-system calling and load it into Sharur.

Sibling of ``run_defensefinder_systems.py`` and deliberately the same shape:
Astra's stage-04 HMM scan gives per-protein profile hits and nothing more,
while the system-level call requires co-localization, mandatory/accessory
quorum, and a model. DefenseFinder supplies that for anti-phage systems; TXSScan
supplies it for the 14 secretion systems and surface appendages (T1SS-T9SS,
T4aP/T4bP, Tad pilus, Flagellum, Com pilus, MSH, archaeal TFF).

Why this matters here rather than being a nice-to-have: without a system-level
caller, a reading agent sees `T2SSE`/`T2SSF` and calls a T2SS. Those PFAMs are
shared between T2SS and the type-IV-filament superfamily and do not
discriminate, so a shared-fold hit gets promoted to a named system. TXSScan is
the model set that separates them.

Input comes from the dataset DuckDB (``proteins.sequence``) rather than a
stage03 directory, so it works for subset builds that were assembled straight
into a database. Pass --data-dir for a full ingest tree if one exists.

Usage:
    python scripts/run_txsscan_systems.py DATASET_DIR [--workers N]
        [--models TXSScan] [--shard i/N] [--faa-dir DIR] [--dry-run]
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import duckdb
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# TXSScan model set. `all` runs every model in the package.
DEFAULT_MODELS = ("TXSScan", "all")


def export_genome_faas(db_path: Path, faa_dir: Path, shard: tuple[int, int] | None) -> list[Path]:
    """Write one gene-ordered FAA per bin from the dataset DuckDB.

    Order matters: MacSyFinder's `ordered_replicon` mode uses adjacency to
    enforce co-localization, so proteins must be emitted in genomic order.
    Sorting is (contig_id, gene_index) and the contig of every protein is
    recorded alongside, because concatenating a MAG's contigs into one pseudo
    -replicon lets multi_loci models gather components from different contigs.
    Under the default `gembase` db-type each contig is its own replicon, so that
    cannot happen; the multi-contig filter in `load_results` then acts as an
    assertion that should drop zero systems.
    """
    faa_dir.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(str(db_path), read_only=True, config={"threads": 4})
    bins = [r[0] for r in con.execute("SELECT DISTINCT bin_id FROM proteins ORDER BY bin_id").fetchall()]
    if shard is not None:
        idx, total = shard
        bins = [b for i, b in enumerate(bins) if i % total == idx]
        logger.info("shard %s/%s -> %d bins", idx, total, len(bins))

    written: list[Path] = []
    for bin_id in bins:
        out = faa_dir / f"{bin_id}.faa"
        if out.exists() and out.stat().st_size > 0:
            written.append(out)
            continue
        rows = con.execute(
            """
            SELECT protein_id, contig_id, sequence
            FROM proteins
            WHERE bin_id = ? AND sequence IS NOT NULL AND length(sequence) > 0
            ORDER BY contig_id, gene_index, start
            """,
            [bin_id],
        ).fetchall()
        if not rows:
            continue
        tmp = out.with_suffix(".faa.tmp")
        with tmp.open("w") as fh:
            for protein_id, _contig, seq in rows:
                fh.write(f">{protein_id}\n{seq}\n")
        tmp.replace(out)
        written.append(out)
    con.close()
    logger.info("exported %d genome FAAs to %s", len(written), faa_dir)
    return written


def run_txsscan_on_genome(
    faa_path: Path, output_dir: Path, workers: int, models: tuple[str, ...],
    db_type: str = "gembase",
) -> pd.DataFrame:
    """Run macsyfinder on one genome; return its best_solution rows."""
    genome_name = faa_path.stem
    genome_out = output_dir / genome_name
    best = genome_out / "best_solution.tsv"

    if best.exists():  # resume: already done
        try:
            df = pd.read_csv(best, sep="\t", comment="#")
            df["genome_id"] = genome_name
            return df
        except Exception:  # noqa: BLE001 - fall through and recompute
            pass
    if genome_out.exists():
        shutil.rmtree(genome_out)

    cmd = [
        "macsyfinder",
        "--models", *models,
        "--sequence-db", str(faa_path),
        "--db-type", db_type,
        "--out-dir", str(genome_out),
        "--worker", str(workers),
        "--mute",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0 and not best.exists():
        logger.warning("  macsyfinder rc=%s for %s: %s", result.returncode, genome_name, (result.stderr or "")[-300:])
        return pd.DataFrame()
    if not best.exists() or best.stat().st_size == 0:
        return pd.DataFrame()
    try:
        df = pd.read_csv(best, sep="\t", comment="#")
    except Exception as exc:  # noqa: BLE001
        logger.warning("  unreadable best_solution for %s: %s", genome_name, exc)
        return pd.DataFrame()
    if df.empty:
        return df
    df["genome_id"] = genome_name
    return df


def load_results(db_path: Path, genes: pd.DataFrame, source: str = "txsscan") -> dict[str, int]:
    """Populate `secretion_systems` (+ membership) from MacSyFinder gene rows.

    Systems whose members span more than one contig are dropped: they are an
    artifact of concatenating a MAG's contigs into a single pseudo-replicon, not
    biology. The count is reported so the discard is never silent.
    """
    if genes.empty:
        logger.warning("no TXSScan hits to load")
        return {"systems": 0, "dropped_multicontig": 0, "members": 0}

    con = duckdb.connect(str(db_path))
    contig_of = dict(con.execute("SELECT protein_id, contig_id FROM proteins").fetchall())

    systems: list[tuple] = []
    members: list[tuple] = []
    dropped = 0
    now = datetime.now(timezone.utc)

    for (genome_id, sys_id), grp in genes.groupby(["genome_id", "sys_id"], sort=True):
        hits = [str(h) for h in grp["hit_id"].tolist()]
        contigs = {contig_of.get(h) for h in hits}
        contigs.discard(None)
        if len(contigs) != 1:
            dropped += 1
            continue
        contig_id = contigs.pop()
        model = str(grp["model_fqn"].iloc[0]) if "model_fqn" in grp else ""
        # model_fqn looks like TXSScan/bacterial/T6SSi -> type=T6SSi
        parts = [p for p in model.split("/") if p]
        system_type = parts[-1] if parts else "unknown"
        subtype = parts[-2] if len(parts) > 1 else ""
        profiles = ",".join(str(g) for g in grp.get("gene_name", pd.Series(dtype=str)).tolist())
        positions = sorted(int(p) for p in grp.get("hit_pos", pd.Series(dtype=int)).tolist()) or [0]
        systems.append(
            (
                f"{genome_id}__{sys_id}",
                genome_id,
                contig_id,
                system_type,
                subtype,
                len(hits),
                ",".join(hits),
                profiles,
                str(positions[0]),
                str(positions[-1]),
                now,
            )
        )
        for pos, hit in enumerate(hits, start=1):
            members.append((f"{genome_id}__{sys_id}", hit, source, pos, "", None))

    con.execute(
        """
        CREATE TABLE IF NOT EXISTS secretion_systems (
            system_id VARCHAR PRIMARY KEY, genome_id VARCHAR, contig_id VARCHAR,
            system_type VARCHAR, system_subtype VARCHAR, genes_count INTEGER,
            protein_ids VARCHAR, profile_names VARCHAR, sys_beg VARCHAR,
            sys_end VARCHAR, created_at TIMESTAMP
        )
        """
    )
    con.execute("DELETE FROM secretion_systems")
    if systems:
        con.executemany(
            "INSERT INTO secretion_systems VALUES (?,?,?,?,?,?,?,?,?,?,?)", systems
        )
    # system_proteins is shared with DefenseFinder; only replace our own source.
    try:
        con.execute("DELETE FROM system_proteins WHERE system_source = ?", [source])
        if members:
            con.executemany(
                "INSERT INTO system_proteins(system_id, protein_id, system_source, position, profile_name, score) "
                "VALUES (?,?,?,?,?,?)",
                members,
            )
    except Exception as exc:  # noqa: BLE001 - membership is a bonus, not the point
        logger.warning("could not write system_proteins: %s", exc)
    con.close()

    logger.info(
        "loaded %d secretion systems (%d members); dropped %d multi-contig artifacts",
        len(systems), len(members), dropped,
    )
    return {"systems": len(systems), "dropped_multicontig": dropped, "members": len(members)}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("data_dir", type=Path, help="Dataset dir containing sharur.duckdb")
    ap.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_CPUS_ON_NODE", 8)))
    ap.add_argument("--models", nargs="+", default=list(DEFAULT_MODELS))
    ap.add_argument("--faa-dir", type=Path, default=None, help="Where to stage per-genome FAAs")
    ap.add_argument("--out-dir", type=Path, default=None, help="MacSyFinder output root")
    ap.add_argument("--shard", type=str, default=None, help="i/N — process every Nth genome (SLURM arrays)")
    ap.add_argument(
        "--db-type",
        default="gembase",
        choices=["gembase", "ordered_replicon", "unordered"],
        help=(
            "gembase scopes each CONTIG as its own replicon and is correct for MAGs. "
            "ordered_replicon declares the whole bin one replicon, which lets models "
            "permitting multiple loci assemble a system from clusters on different "
            "contigs -- an adjacency that does not exist in the genome."
        ),
    )
    ap.add_argument("--load-only", action="store_true", help="Skip running; load existing results")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    db_path = args.data_dir / "sharur.duckdb"
    if not db_path.exists():
        logger.error("no database at %s", db_path)
        return 1

    faa_dir = args.faa_dir or (args.data_dir / "txsscan" / "faa")
    out_dir = args.out_dir or (args.data_dir / "txsscan" / "results")
    out_dir.mkdir(parents=True, exist_ok=True)

    shard = None
    if args.shard:
        i, n = args.shard.split("/")
        shard = (int(i), int(n))

    if not args.load_only:
        faas = export_genome_faas(db_path, faa_dir, shard)
        if args.dry_run:
            logger.info("dry run: would run TXSScan on %d genomes with %d workers", len(faas), args.workers)
            return 0
        for idx, faa in enumerate(faas, start=1):
            run_txsscan_on_genome(faa, out_dir, args.workers, tuple(args.models), args.db_type)
            if idx % 25 == 0:
                logger.info("  %d/%d genomes", idx, len(faas))

    # Collect every shard's results, then load once.
    frames: list[pd.DataFrame] = []
    for best in sorted(out_dir.glob("*/best_solution.tsv")):
        try:
            df = pd.read_csv(best, sep="\t", comment="#")
        except Exception:  # noqa: BLE001
            continue
        if df.empty:
            continue
        df["genome_id"] = best.parent.name
        frames.append(df)

    if shard is not None and not args.load_only:
        logger.info("shard complete; %d genomes with hits. Run --load-only to integrate.", len(frames))
        return 0

    genes = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    stats = load_results(db_path, genes)
    (out_dir.parent / "txsscan_manifest.json").write_text(
        json.dumps({"generated": datetime.now(timezone.utc).isoformat(), "models": args.models, **stats}, indent=2)
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
