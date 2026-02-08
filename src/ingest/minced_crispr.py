#!/usr/bin/env python3
"""
Stage 05c: CRISPR array detection via minced.

Runs minced on each genome FASTA and emits *_crispr_arrays.json expected by Stage07.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path
from typing import List, Dict, Any

import typer
from rich.console import Console
from rich.progress import Progress

console = Console()


def run_minced(fasta: Path, out_dir: Path) -> Dict[str, Any]:
    """Run minced on a single fasta and return parsed arrays."""
    exe = shutil.which("minced")
    if exe is None:
        raise FileNotFoundError("minced executable not found on PATH")

    out_dir.mkdir(parents=True, exist_ok=True)
    gff_out = out_dir / f"{fasta.stem}_crispr.gff"
    txt_out = out_dir / f"{fasta.stem}_crispr.txt"

    # MinCED usage: minced input.fa output.txt output.gff
    cmd = [exe, str(fasta), str(txt_out), str(gff_out)]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    if result.returncode != 0:
        raise RuntimeError(f"minced failed: {result.stderr}")

    arrays = parse_minced_gff(gff_out)
    return {"genome": fasta.stem, "arrays": arrays}


def parse_minced_gff(gff_path: Path) -> List[Dict[str, Any]]:
    """Parse minced GFF; arrays are marked as 'repeat_region' features with rpt_family=CRISPR."""
    arrays = []
    if not gff_path.exists():
        return arrays
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            seqid, _, ftype, start, end, score, strand, phase, attrs = parts
            # MinCED outputs 'repeat_region' as feature type with rpt_family=CRISPR in attributes
            if ftype != "repeat_region":
                continue
            if "rpt_family=CRISPR" not in attrs:
                continue
            attr_dict = {}
            for kv in attrs.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    attr_dict[k] = v
            arrays.append(
                {
                    "id": attr_dict.get("ID", f"{seqid}_{start}_{end}"),
                    "contig": seqid,
                    "startCoordinate": int(start),
                    "endCoordinate": int(end),
                    "strand": strand,
                    "metadata": attr_dict,
                }
            )
    return arrays


def run(
    input_dir: Path = typer.Option(..., "--input-dir", "-i"),
    output_dir: Path = typer.Option(..., "--output-dir", "-o"),
    force: bool = typer.Option(False, "--force"),
) -> None:
    if force and output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fasta_files = []
    for ext in ("*.fna", "*.fa", "*.fasta"):
        fasta_files.extend(input_dir.glob(ext))
    if not fasta_files:
        console.print("[red]No FASTA files found for CRISPR detection.[/red]")
        raise typer.Exit(1)

    combined = []
    with Progress() as progress:
        task = progress.add_task("Running minced...", total=len(fasta_files))
        for fasta in fasta_files:
            try:
                result = run_minced(fasta, output_dir)
                combined.append(result)
                progress.advance(task)
            except Exception as exc:
                console.print(f"[yellow]minced failed on {fasta.name}: {exc}[/yellow]")
                progress.advance(task)

    # Write per-genome JSONs
    for entry in combined:
        path = output_dir / f"{entry['genome']}_crispr_arrays.json"
        path.write_text(json.dumps({"arrays": entry["arrays"]}, indent=2))

    console.print(f"[green]CRISPR detection complete[/green]: wrote {len(combined)} JSONs to {output_dir}")


if __name__ == "__main__":
    typer.run(run)
