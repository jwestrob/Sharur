#!/usr/bin/env python3
"""Call biosynthetic gene clusters with GECCO on a dataset's own gene calls.

GECCO calls genes itself by default, and the protein IDs it then reports exist
nowhere in the database. This runner gives GECCO the database's gene calls
instead. For each genome it writes GenBank records whose CDS features carry the
stored protein ID (``locus_tag``), coordinates, strand and translation, and runs
``gecco run --cds-feature CDS --locus-tag locus_tag``. Every cluster member is
therefore a database protein, and ``sharur load-bgcs`` loads the calls as
``loci`` rows with ``locus_type='bgc'``.

A shard writes all of its genomes into one GenBank file and makes a single GECCO
run, so the domain annotation runs as one multithreaded pass. GenBank caps
LOCUS names, so records get short names (``r<shard>_<n>``); ``record_map.tsv``
maps each one back to ``(bin_id, contig_id)``. Each genome contributes only its
own proteins, as stored in the database, on the contigs that carry them.

Usage:
    python scripts/run_gecco_bgcs.py DB --assemblies MAP.tsv --out-dir DIR \\
        [--shard i/N] [--jobs N]
    sharur load-bgcs --db DB --gecco-dir DIR

MAP.tsv has two tab-separated columns, bin_id and the assembly FASTA path,
with no header.
"""

from __future__ import annotations

import argparse
import gzip
import logging
import os
import shutil
import subprocess
import sys
from pathlib import Path

import duckdb

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def read_fasta(path: Path, wanted: set[str]) -> dict[str, str]:
    """Record ID (first header word) -> sequence, for the wanted IDs only."""
    opener = gzip.open if path.suffix == ".gz" else open
    seqs: dict[str, list[str]] = {}
    current = None
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                rid = line[1:].split(None, 1)[0] if line[1:].strip() else ""
                current = rid if rid in wanted else None
                if current is not None:
                    seqs[current] = []
            elif current is not None:
                seqs[current].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def write_shard_genbank(db_path: Path, assemblies: dict[str, Path], bins: list[str],
                        shard_dir: Path, shard_tag: str) -> dict[str, int]:
    """Write input.gbk and record_map.tsv for ``bins``; return counts."""
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO

    con = duckdb.connect(str(db_path), read_only=True, config={"threads": 4})
    stats = {"genomes": 0, "records": 0, "cds": 0, "genomes_without_assembly": 0,
             "genomes_without_proteins": 0, "contigs_missing_from_assembly": 0,
             "cds_out_of_bounds": 0}
    tmp_gbk = shard_dir / "input.gbk.tmp"
    tmp_map = shard_dir / "record_map.tsv.tmp"
    n_record = 0
    with tmp_gbk.open("w") as gbk, tmp_map.open("w") as rmap:
        rmap.write("record\tbin_id\tcontig_id\n")
        for bin_id in bins:
            path = assemblies.get(bin_id)
            if path is None:
                stats["genomes_without_assembly"] += 1
                logger.warning("no assembly for %s", bin_id)
                continue
            rows = con.execute(
                """SELECT protein_id, contig_id, start, end_coord, strand, sequence
                   FROM proteins
                   WHERE bin_id = ? AND start > 0 AND sequence IS NOT NULL AND length(sequence) > 0
                   ORDER BY contig_id, start, end_coord""", [bin_id]).fetchall()
            if not rows:
                stats["genomes_without_proteins"] += 1
                continue
            by_contig: dict[str, list[tuple]] = {}
            for row in rows:
                by_contig.setdefault(row[1], []).append(row)
            seqs = read_fasta(path, set(by_contig))
            stats["genomes"] += 1
            records = []
            for contig_id, genes in by_contig.items():
                seq = seqs.get(contig_id)
                if seq is None:
                    stats["contigs_missing_from_assembly"] += 1
                    continue
                name = f"r{shard_tag}_{n_record}"
                n_record += 1
                features = []
                for protein_id, _c, start, end, strand, aa in genes:
                    if start < 1 or end > len(seq) or start > end:
                        stats["cds_out_of_bounds"] += 1
                        continue
                    features.append(SeqFeature(
                        FeatureLocation(start - 1, end, strand=1 if strand == "+" else -1),
                        type="CDS",
                        qualifiers={"locus_tag": [protein_id], "translation": [aa.rstrip("*")]},
                    ))
                record = SeqRecord(Seq(seq), id=name, name=name, description=f"{bin_id} {contig_id}",
                                   features=features, annotations={"molecule_type": "DNA"})
                records.append(record)
                rmap.write(f"{name}\t{bin_id}\t{contig_id}\n")
                stats["records"] += 1
                stats["cds"] += len(features)
            SeqIO.write(records, gbk, "genbank")
    con.close()
    tmp_gbk.replace(shard_dir / "input.gbk")
    tmp_map.replace(shard_dir / "record_map.tsv")
    return stats


def run_gecco(shard_dir: Path, jobs: int) -> None:
    """One GECCO run over the shard's GenBank file, using the stored gene calls."""
    version = subprocess.run(["gecco", "--version"], capture_output=True, text=True)
    (shard_dir / "gecco_version.txt").write_text((version.stdout or version.stderr).strip() + "\n")
    out = shard_dir / "gecco_out"
    if out.exists():
        shutil.rmtree(out)
    cmd = ["gecco", "run", "-v", "--genome", str(shard_dir / "input.gbk"), "--format", "genbank",
           "--cds-feature", "CDS", "--locus-tag", "locus_tag",
           "-o", str(out), "--jobs", str(jobs)]
    logger.info("running %s", " ".join(cmd))
    with (shard_dir / "gecco.log").open("w") as log:
        result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        raise RuntimeError(f"gecco exited {result.returncode}; see {shard_dir / 'gecco.log'}")
    tables = sorted(out.glob("*.clusters.tsv"))
    if not tables:
        raise RuntimeError(f"gecco wrote no clusters table in {out}")
    shutil.copy2(tables[0], shard_dir / "input.clusters.tsv")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("db", type=Path, help="Dataset sharur.duckdb (opened read-only)")
    ap.add_argument("--assemblies", type=Path, required=True, help="bin_id<TAB>assembly FASTA path")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--shard", default=None, help="i/N: every Nth genome in bin_id order")
    ap.add_argument("--jobs", type=int, default=int(os.environ.get("SLURM_CPUS_ON_NODE", os.cpu_count() or 1)))
    args = ap.parse_args()

    assemblies = {}
    for line in args.assemblies.read_text().splitlines():
        if line.strip():
            bin_id, path = line.split("\t")[:2]
            assemblies[bin_id] = Path(path)
    con = duckdb.connect(str(args.db), read_only=True, config={"threads": 4})
    bins = [r[0] for r in con.execute("SELECT DISTINCT bin_id FROM proteins ORDER BY bin_id").fetchall()]
    con.close()
    tag = "all"
    if args.shard:
        i, n = (int(x) for x in args.shard.split("/"))
        bins = [b for k, b in enumerate(bins) if k % n == i]
        tag = f"{i:03d}"
    shard_dir = args.out_dir / f"shard_{tag}"
    shard_dir.mkdir(parents=True, exist_ok=True)
    if (shard_dir / "input.clusters.tsv").is_file():
        logger.info("shard %s already has a clusters table; nothing to do", tag)
        return 0
    stats = write_shard_genbank(args.db, assemblies, bins, shard_dir, tag)
    logger.info("shard %s: %s", tag, stats)
    (shard_dir / "export_stats.tsv").write_text("".join(f"{k}\t{v}\n" for k, v in stats.items()))
    if stats["records"] == 0:  # nothing to scan: an empty table marks the shard finished
        (shard_dir / "input.clusters.tsv").write_text("sequence_id\tbgc_id\tstart\tend\ttype\tproteins\n")
        return 0
    run_gecco(shard_dir, args.jobs)
    logger.info("shard %s done", tag)
    return 0


if __name__ == "__main__":
    sys.exit(main())
