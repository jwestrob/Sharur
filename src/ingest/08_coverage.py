#!/usr/bin/env python3
"""
Stage 8 (optional): per-sample read coverage with CoverM.

Maps each sample's reads to the dataset's prepared assemblies (concatenated
into one reference) with ``coverm contig`` and imports mean depth, covered
fraction, read count and contig length into the abundance sidecar
(``abundance.duckdb`` beside the dataset DuckDB).

The reads table is tab-separated with a header: ``sample_id``, ``read1`` and
optionally ``read2`` (paired) or ``interleaved`` (true/false); any other columns
are kept as sample metadata.

CoverM runs with every allocated CPU (``$SLURM_CPUS_ON_NODE`` when set).
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import shutil
import subprocess
from pathlib import Path

from sharur.abundance import default_path, import_coverage
from sharur.contig_context import assemblies_from_stage00
from sharur.storage.duckdb_store import DuckDBStore


def read_samples(path: Path) -> list[dict[str, str]]:
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or not {"sample_id", "read1"} <= set(rows[0]):
        raise SystemExit(f"{path}: needs a header with sample_id and read1 (read2 or interleaved optional)")
    for row in rows:
        for column in ("read1", "read2"):
            if row.get(column) and not Path(row[column]).is_file():
                raise SystemExit(f"Read file not found for {row['sample_id']}: {row[column]}")
    return rows


def write_reference(assemblies: dict[str, Path], target: Path) -> Path:
    """Concatenate the prepared assemblies (contig IDs are unique across the dataset)."""
    target.parent.mkdir(parents=True, exist_ok=True)
    with open(target, "w") as out:
        for _, path in sorted(assemblies.items()):
            opener = gzip.open if path.suffix == ".gz" else open
            with opener(path, "rt") as handle:
                shutil.copyfileobj(handle, out)
    return target


def coverm_command(reference: Path, sample: dict[str, str], output: Path, threads: int) -> list[str]:
    command = ["coverm", "contig", "--reference", str(reference), "-t", str(threads),
               "-m", "mean", "covered_fraction", "count", "length", "-o", str(output)]
    if sample.get("read2"):
        command += ["-1", sample["read1"], "-2", sample["read2"]]
    elif str(sample.get("interleaved", "")).lower() in ("1", "true", "yes"):
        command += ["--interleaved", sample["read1"]]
    else:
        command += ["--single", sample["read1"]]
    return command


def rename_sample_columns(path: Path, sample_id: str) -> None:
    """CoverM names columns after the read file; rename them to the sample ID."""
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    suffixes = (" Mean", " Covered Fraction", " Read Count", " Length")
    header = [header[0]] + [next((sample_id + s for s in suffixes if h.endswith(s)), h) for h in header[1:]]
    path.write_text("\n".join(["\t".join(header), *lines[1:]]) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description="Per-sample read coverage with CoverM.")
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--db", type=Path, required=True)
    parser.add_argument("--reads", type=Path, required=True, help="TSV: sample_id, read1[, read2 | interleaved]")
    parser.add_argument("--threads", type=int, default=None)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    if shutil.which("coverm") is None:
        raise SystemExit("coverm is not on PATH; install CoverM to compute coverage")
    threads = args.threads or int(os.environ.get("SLURM_CPUS_ON_NODE") or os.cpu_count() or 1)
    samples = read_samples(args.reads)
    assemblies = assemblies_from_stage00(args.data_dir / "stage00_prepared")
    if not assemblies:
        raise SystemExit("No prepared assemblies found in stage00_prepared")
    out_dir = args.data_dir / "stage08_coverage"
    reference = out_dir / "reference.fna"
    if args.force or not reference.is_file():
        write_reference(assemblies, reference)

    tables = []
    for sample in samples:
        table = out_dir / f"{sample['sample_id']}.coverm.tsv"
        if args.force or not table.is_file():
            command = coverm_command(reference, sample, table, threads)
            print("+", " ".join(command), flush=True)
            subprocess.run(command, check=True)
            rename_sample_columns(table, sample["sample_id"])
        tables.append(table)

    metadata = out_dir / "sample_metadata.tsv"
    extra = [c for c in samples[0] if c not in ("sample_id", "read1", "read2", "interleaved")]
    with open(metadata, "w") as handle:
        handle.write("\t".join(["sample_id", *extra]) + "\n")
        for sample in samples:
            handle.write("\t".join([sample["sample_id"], *(sample.get(c, "") for c in extra)]) + "\n")

    with DuckDBStore(args.db, read_only=True) as store:
        result = import_coverage(store, default_path(args.db), tables, fmt="coverm",
                                 sample_metadata=metadata, source="coverm")
    manifest = {"stage": "stage08_coverage", "samples": [s["sample_id"] for s in samples],
                "threads": threads, "reference": str(reference), "import": result}
    (out_dir / "processing_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"Coverage imported for {len(result['samples'])} samples ({result['rows']:,} contig rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
