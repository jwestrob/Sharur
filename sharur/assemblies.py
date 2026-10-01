"""Where a dataset's assembly FASTAs live, and what their sizes say about gene calls.

Shared by the browser (CRISPR repeats, the presence/absence matrix) and the
health checks, so they agree on which assembly belongs to which genome and on
when a genome's gene calls are mostly missing.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from statistics import median

ASSEMBLY_SUFFIXES = (".fna", ".fa", ".fasta", ".fas", ".fna.gz", ".fa.gz", ".fasta.gz")
ASSEMBLY_DIRS = ("stage00_prepared/genomes", "genomes_fna", "genomes_fna_new", "source", "assemblies")
# Prokaryotic genomes carry about one gene per kb; under half that, most gene calls are missing.
DEFICIT = 0.5
# A FASTA file holds about 1.3% newlines and headers beyond its bases.
_FASTA_OVERHEAD = 1013.0


def find_assemblies(dataset_dir: Path, extra: list[Path] | None = None) -> dict[str, Path]:
    """Genome ID -> assembly FASTA: the stage 00 manifest, then FASTAs named after their genome."""
    dataset_dir = Path(dataset_dir)
    paths: dict[str, Path] = {}
    manifest = dataset_dir / "stage00_prepared" / "processing_manifest.json"
    if manifest.is_file():
        try:
            for entry in json.loads(manifest.read_text()).get("genomes", []):
                if entry.get("genome_id") and entry.get("output_path") and Path(entry["output_path"]).is_file():
                    paths.setdefault(entry["genome_id"], Path(entry["output_path"]))
        except (ValueError, OSError):
            pass
    for d in [*(extra or []), *(dataset_dir / name for name in ASSEMBLY_DIRS)]:
        if not Path(d).is_dir():
            continue
        for f in Path(d).iterdir():
            for suffix in ASSEMBLY_SUFFIXES:
                if f.name.endswith(suffix):
                    paths.setdefault(f.name[: -len(suffix)], f)
                    break
    return paths


def assembly_kb(paths: dict[str, Path]) -> dict[str, float]:
    """Approximate assembly length (kb) per genome from uncompressed FASTA file sizes."""
    sizes: dict[str, float] = {}
    for bin_id, path in paths.items():
        if not str(path).endswith(".gz"):
            try:
                sizes[bin_id] = os.path.getsize(path) / _FASTA_OVERHEAD
            except OSError:
                pass
    return sizes


def gene_call_deficits(proteins: dict[str, int], completeness: dict[str, float | None],
                       kb: dict[str, float] | None = None, reference: list[str] | None = None) -> set[str]:
    """Genomes whose gene calls are mostly missing.

    With the assembly at hand: under ``DEFICIT`` proteins per kb. Otherwise:
    under ``DEFICIT`` × the proteins its completeness implies, using the median
    proteins-per-completeness of ``reference`` genomes (default: all).
    """
    kb = kb or {}
    ref = reference if reference is not None else list(proteins)
    ratios = [proteins[b] / (completeness[b] / 100.0) for b in ref
              if proteins.get(b) and completeness.get(b)]
    per_complete = median(ratios) if len(ratios) >= 5 else None
    flagged = set()
    for b, n in proteins.items():
        if kb.get(b):
            if n < DEFICIT * kb[b]:
                flagged.add(b)
        elif per_complete and completeness.get(b) and n < DEFICIT * per_complete * completeness[b] / 100.0:
            flagged.add(b)
    return flagged


__all__ = ["ASSEMBLY_DIRS", "ASSEMBLY_SUFFIXES", "DEFICIT", "assembly_kb", "find_assemblies", "gene_call_deficits"]
