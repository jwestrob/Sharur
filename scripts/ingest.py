#!/usr/bin/env python3
"""
Thin orchestration CLI for Bennu ingest stages (00→07).

Two modes:
- fast (default): synthesize minimal stage outputs from contig FASTAs (no external tools)
- tools: call the existing stage scripts; individual stages can be skipped

Goal: build a DuckDB at the end (Stage 07) with one command, especially against
`dummy_dataset/` for local smoke tests.
"""

from __future__ import annotations

import json
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated, Iterable, Tuple

import typer
from rich.console import Console

console = Console()
app = typer.Typer(no_args_is_help=True, add_completion=False)

# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #


def _run(script_path: Path, args: list[str]) -> list[str]:
    """Run a stage script via python <path> ... and return the command list."""
    cmd = [sys.executable, str(script_path), *args]
    console.print(f"[dim]$ {' '.join(cmd)}[/dim]")
    subprocess.run(cmd, check=True)
    return cmd


def _force_clean(path: Path, enable: bool) -> None:
    """If force is set and path exists, remove it to avoid stage overwrite errors."""
    if enable and path.exists():
        import shutil

        shutil.rmtree(path)


def _parse_contigs(fasta_path: Path) -> Iterable[Tuple[str, str]]:
    """Yield (contig_id, sequence) pairs from FASTA."""
    header = None
    seq_parts: list[str] = []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
        if header is not None:
            yield header, "".join(seq_parts)


@dataclass
class StagePaths:
    stage00: Path
    stage01: Path
    stage02: Path
    stage03: Path
    stage04: Path
    stage05a: Path
    stage05b: Path
    stage05c: Path
    stage06: Path

    @classmethod
    def from_root(cls, root: Path) -> "StagePaths":
        return cls(
            stage00=root / "stage00_prepared",
            stage01=root / "stage01_quast",
            stage02=root / "stage02_dfast_qc",
            stage03=root / "stage03_prodigal",
            stage04=root / "stage04_astra",
            stage05a=root / "stage05a_gecco",
            stage05b=root / "stage05b_dbcan",
            stage05c=root / "stage05c_crispr",
            stage06=root / "stage06_embeddings",
        )


# --------------------------------------------------------------------------- #
# Fast (synthetic) pipeline
# --------------------------------------------------------------------------- #


def _synthesize_pipeline(raw_dir: Path, data_dir: Path) -> None:
    """
    Build minimal stage outputs from contig FASTAs (mirrors integration tests).

    Produces:
      - stage02 manifest
      - stage03 Prodigal-like FAA (1 gene per contig)
      - stage04 PFAM-like TSV (one annotation per bin)
      - stage05a GECCO-like BGC JSON (one cluster per bin)
      - stage05b dbCAN-like CAZy JSON (empty ok)
      - stage05c CRISPR JSON (empty)
    """
    stages = StagePaths.from_root(data_dir)
    stages.stage02.mkdir(parents=True, exist_ok=True)
    stages.stage03.mkdir(parents=True, exist_ok=True)
    stages.stage04.mkdir(parents=True, exist_ok=True)
    stages.stage05a.mkdir(parents=True, exist_ok=True)
    stages.stage05b.mkdir(parents=True, exist_ok=True)
    stages.stage05c.mkdir(parents=True, exist_ok=True)

    genomes_manifest = {"genomes": []}
    pfam_rows = []
    gecco_clusters = []

    fasta_paths = sorted(raw_dir.glob("*.fna")) + sorted(raw_dir.glob("*.fa")) + sorted(raw_dir.glob("*.fasta"))
    if not fasta_paths:
        raise FileNotFoundError(f"No FASTA files found in {raw_dir}")

    for fasta_path in fasta_paths:
        raw_bin_id = fasta_path.stem.replace(".contigs", "")
        bin_id = raw_bin_id.replace("_", "")
        contigs = list(_parse_contigs(fasta_path))

        genomes_manifest["genomes"].append(
            {
                "genome_id": bin_id,
                "taxonomy": {"name": "unknown", "completeness": 90.0, "contamination": 1.0},
                "n_contigs": len(contigs),
                "total_length": sum(len(seq) for _, seq in contigs),
            }
        )

        # Stage03 proteins (one per contig)
        bin_dir = stages.stage03 / "genomes" / bin_id
        bin_dir.mkdir(parents=True, exist_ok=True)
        faa_path = bin_dir / f"{bin_id}.faa"
        with open(faa_path, "w") as faa:
            for idx, (contig_id, seq) in enumerate(contigs):
                pid = f"{bin_id}_ctg{idx:05d}_00001"
                faa.write(f">{pid} # 1 # {min(len(seq), 300)} # 1 # ID={pid};partial=00\n")
                faa.write("M" * min(len(seq), 100) + "\n")
                if idx == 0:
                    pfam_rows.append(
                        {
                            "sequence_id": pid,
                            "hmm_name": "PF00001",
                            "e_value": 1e-5,
                            "bit_score": 42.0,
                            "ali_from": 1,
                            "ali_to": 30,
                            "name": "MockDomain",
                            "description": "Synthetic annotation",
                        }
                    )
                    gecco_clusters.append(
                        {
                            "cluster_id": f"{bin_id}_cluster1",
                            "contig": pid.rsplit("_", 1)[0],
                            "start": 1,
                            "end": 100,
                            "bgc_type": "nrps",
                            "protein_list": [pid],
                        }
                    )

    # Stage02 manifest
    (stages.stage02 / "processing_manifest.json").write_text(json.dumps(genomes_manifest))

    # Stage04 PFAM TSV (if any)
    if pfam_rows:
        header = [
            "sequence_id",
            "hmm_name",
            "e_value",
            "bit_score",
            "ali_from",
            "ali_to",
            "name",
            "description",
        ]
        lines = ["\t".join(header)]
        for row in pfam_rows:
            lines.append("\t".join(str(row[h]) for h in header))
        (stages.stage04 / "synthetic_hits_df.tsv").write_text("\n".join(lines) + "\n")

    # Stage05a GECCO clusters
    (stages.stage05a / "combined_bgc_data.json").write_text(json.dumps({"clusters": gecco_clusters}))

    # Stage05b CAZy annotations (empty allowed)
    (stages.stage05b / "synthetic_cazyme_results.json").write_text(json.dumps({"annotations": []}))

    # Stage05c CRISPR arrays (empty)
    (stages.stage05c / "synthetic_crispr_arrays.json").write_text(json.dumps({"arrays": []}))

    console.print(f"[green]Synthetic stages written under {data_dir}[/green]")


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #


@app.command()
def run(
    input_dir: Annotated[Path, typer.Option("--input-dir", "-i", help="Raw genomes (.fna/.fa/.fasta)")] = Path("dummy_dataset"),
    data_dir: Annotated[Path, typer.Option("--data-dir", "-d", help="Working directory for stage outputs")] = Path("data"),
    output: Annotated[Path, typer.Option("--output", "-o", help="Destination DuckDB path")] = Path("data/bennu.duckdb"),
    mode: Annotated[str, typer.Option("--mode", "-m", help="tools (real pipeline) or fast (synthetic smoke)")] = "tools",
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing DuckDB")] = False,
    skip_quast: Annotated[bool, typer.Option(help="Skip Stage 01 (QUAST)")] = False,
    skip_dfast: Annotated[bool, typer.Option(help="Skip Stage 02 (DFAST QC)")] = False,
    skip_prodigal: Annotated[bool, typer.Option(help="Skip Stage 03 (Prodigal)")] = False,
    skip_astra: Annotated[bool, typer.Option(help="Skip Stage 04 (Astra scan / PFAM/KOFAM)")] = False,
    skip_gecco: Annotated[bool, typer.Option(help="Skip Stage 05a (GECCO BGC)")] = False,
    skip_dbcan: Annotated[bool, typer.Option(help="Skip Stage 05b (dbCAN CAZy)")] = False,
    skip_crispr: Annotated[bool, typer.Option(help="Skip Stage 05c (CRISPR)")] = False,
    skip_embeddings: Annotated[bool, typer.Option(help="Skip Stage 06 (ESM2 embeddings)")] = False,
    dry_run: Annotated[bool, typer.Option(help="Plan commands without executing or building DB")] = False,
) -> None:
    """
    Build a Bennu DuckDB. Default is the real tool chain; use --mode fast only for test smoke.
    """
    data_dir.mkdir(parents=True, exist_ok=True)
    stages = StagePaths.from_root(data_dir)
    planned: list[list[str]] = []

    if mode == "fast":
        if dry_run:
            console.print("[yellow]Dry run: fast mode would synthesize stage outputs (no commands to run).[/yellow]")
            return
        _synthesize_pipeline(input_dir, data_dir)
    elif mode == "tools":
        repo_root = Path(__file__).resolve().parents[1]
        # Stage 00
        if force or not stages.stage00.exists():
            _force_clean(stages.stage00, force)
            cmd = [str(repo_root / "src" / "ingest" / "00_prepare_inputs.py"), "-i", str(input_dir), "-o", str(stages.stage00), "--force"]
            planned.append(cmd)
            if not dry_run:
                _run(Path(cmd[0]), cmd[1:])
        # Stage 01
        if not skip_quast:
            _force_clean(stages.stage01, force)
            cmd = [str(repo_root / "src" / "ingest" / "01_run_quast.py"), "-i", str(stages.stage00), "-o", str(stages.stage01)]
            if force:
                cmd.append("--force")
            planned.append(cmd)
            if not dry_run:
                _run(Path(cmd[0]), cmd[1:])
        # Stage 02
        if not skip_dfast:
            _force_clean(stages.stage02, force)
            cmd = [str(repo_root / "src" / "ingest" / "02_dfast_qc.py"), "-i", str(stages.stage00), "-o", str(stages.stage02)]
            if force:
                cmd.append("--force")
            planned.append(cmd)
            if not dry_run:
                _run(Path(cmd[0]), cmd[1:])
        # Stage 03
        if not skip_prodigal:
            _force_clean(stages.stage03, force)
            cmd = [str(repo_root / "src" / "ingest" / "03_prodigal.py"), "-i", str(stages.stage00), "-o", str(stages.stage03)]
            if force:
                cmd.append("--force")
            planned.append(cmd)
            if not dry_run:
                _run(Path(cmd[0]), cmd[1:])
        # Stage 04
        if not skip_astra:
            _force_clean(stages.stage04, force)
            cmd = [str(repo_root / "src" / "ingest" / "04_astra_scan.py"), "-i", str(stages.stage03), "-o", str(stages.stage04)]
            if force:
                cmd.append("--force")
            planned.append(cmd)
            if not dry_run:
                _run(Path(cmd[0]), cmd[1:])
        # Stage 05a
        if not skip_gecco:
            _force_clean(stages.stage05a, force)
            cmd = [str(repo_root / "src" / "ingest" / "gecco_bgc.py"), str(stages.stage00), str(stages.stage05a)]
            if force:
                cmd.append("--force")
            planned.append(cmd)
            if not dry_run:
                _run(Path(cmd[0]), cmd[1:])
        # Stage 05b
        if not skip_dbcan:
            _force_clean(stages.stage05b, force)
            cmd = [
                str(repo_root / "src" / "ingest" / "dbcan_cazyme.py"),
                "--input-dir",
                str(stages.stage03 / "genomes" / "all_protein_symlinks"),
                "--output-dir",
                str(stages.stage05b),
            ]
            planned.append(cmd)
            if not dry_run:
                _run(Path(cmd[0]), cmd[1:])
        # Stage 05c: no script shipped; expect external
        if not skip_crispr:
            _force_clean(stages.stage05c, force)
            cmd = [
                str(repo_root / "src" / "ingest" / "minced_crispr.py"),
                "--input-dir",
                str(stages.stage00),
                "--output-dir",
                str(stages.stage05c),
            ]
            if force:
                cmd.append("--force")
            planned.append(cmd)
            if not dry_run:
                _run(Path(cmd[0]), cmd[1:])
        # Stage 06
        if not skip_embeddings:
            _force_clean(stages.stage06, force)
            cmd = [str(repo_root / "src" / "ingest" / "06_esm2_embeddings.py"), str(stages.stage03), str(stages.stage06)]
            planned.append(cmd)
            if not dry_run:
                _run(Path(cmd[0]), cmd[1:])
    else:
        raise typer.BadParameter("mode must be 'fast' or 'tools'")

    if dry_run:
        console.print("[cyan]Dry run plan (commands in order):[/cyan]")
        for cmd in planned:
            console.print("  " + " ".join(cmd))
        return planned

    # Stage 07: build DuckDB
    # Import builder dynamically because module name starts with a number
    import importlib.util
    repo_root = Path(__file__).resolve().parents[1]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))
    builder_path = repo_root / "src" / "ingest" / "07_build_knowledge_base.py"
    spec = importlib.util.spec_from_file_location("kb_build", builder_path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Unable to load stage07 builder.")
    kb = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = kb  # type: ignore[arg-type]
    spec.loader.exec_module(kb)  # type: ignore

    outputs = kb.PipelineOutputs(
        stage00_dir=stages.stage00,
        stage01_dir=stages.stage01,
        stage02_dir=stages.stage02,
        stage03_dir=stages.stage03,
        stage04_dir=stages.stage04,
        stage05a_dir=stages.stage05a,
        stage05b_dir=stages.stage05b,
        stage05c_dir=stages.stage05c,
        stage06_dir=stages.stage06,
    )
    builder = kb.KnowledgeBaseBuilder(outputs, output, force=force)
    stats = builder.build()
    console.print(f"[green]Build complete[/green]: {stats}")


if __name__ == "__main__":
    app()
