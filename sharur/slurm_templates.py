"""SLURM script templates for the Sharur staged ingest pipeline.

These templates encode the CLAUDE.md SLURM rules:

  - `#SBATCH --exclusive` — never partial `--mem` or `--cpus-per-task`
  - no `#SBATCH --time` / `-t` (biotite has TimeLimit=UNLIMITED on the relevant
    partitions; setting a walltime can only kill the job early)
  - `--dependency=afterok:` for chaining (never poll)
  - Astra array of 5, `--exclusive` per task, `--force` per task per
    `04_astra_scan.py` requirements

`emit_pipeline()` writes the 5 scripts plus `submit.sh` into a target directory.
Designed to be called from `sharur-ingest --emit-slurm`, but can be invoked
directly for ad-hoc dataset bootstrapping.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


@dataclass
class PipelineSpec:
    """Inputs needed to render a SLURM pipeline."""

    dataset_name: str           # e.g. "srvp_bacteria_pb"
    sharur_root: Path           # repo root, used in `cd $SHARUR`
    input_dir: Path             # source FASTAs
    data_dir: Path              # data/<dataset_name>, possibly a symlink
    log_dir: Path               # where SLURM .out files land
    pipeline_dir: Path          # where the .sh files go

    # Per CLAUDE.md: standard for CPU, gpu for embeddings, gpu_h200 only if
    # large memory needed
    cpu_partition: str = "standard"
    gpu_partition: str = "gpu"

    # The 5 Astra DBs to scan in parallel.
    astra_dbs: tuple[str, ...] = ("PFAM", "KOFAM", "HydDB", "DefenseFinder", "dbCAN")

    # Stage skip flags. When True, the corresponding script is NOT emitted and
    # the submit.sh chain skips it.
    skip_astra: bool = False
    skip_crispr: bool = False
    skip_embeddings: bool = True   # default off because GPU; user opts in

    # Optional comment line prepended to each script.
    emitted_by: str = "sharur-ingest --emit-slurm"

    # Files that get written (computed in __post_init__).
    files: dict[str, Path] = field(default_factory=dict)


def _header(spec: PipelineSpec, *, job_name: str, partition: str, log_basename: str,
            extra_sbatch: Optional[list[str]] = None) -> str:
    """Common SBATCH header block honoring CLAUDE.md SLURM rules."""
    lines = [
        "#!/bin/bash",
        f"#SBATCH -J {job_name}",
        f"#SBATCH -p {partition}",
        "#SBATCH -N 1",
        "#SBATCH --exclusive",
        f"#SBATCH -o {spec.log_dir}/{log_basename}",
    ]
    if extra_sbatch:
        lines.extend(extra_sbatch)
    lines.append("")
    lines.append("set -euo pipefail")
    lines.append("")
    lines.append(f"SHARUR={spec.sharur_root}")
    lines.append(f"DATASET={spec.dataset_name}")
    lines.append('cd "$SHARUR"')
    return "\n".join(lines)


def _prep_script(spec: PipelineSpec) -> str:
    """Stage 00 prepare + Stage 03 Prodigal in one job."""
    head = _header(
        spec, job_name=f"{spec.dataset_name[:12]}_prep",
        partition=spec.cpu_partition, log_basename="01_prep_%j.out",
    )
    return f"""{head}

echo "[$(date)] Stage 00 + Stage 03 on $DATASET"
echo "  input:  {spec.input_dir}"
echo "  data:   data/$DATASET -> $(readlink -f data/$DATASET)"
echo "  threads: $SLURM_CPUS_ON_NODE"

# Stage 00: validate + organize bins
python src/ingest/00_prepare_inputs.py \\
    -i "{spec.input_dir}" \\
    -o "data/$DATASET/stage00_prepared"

# Stage 03: Prodigal gene calling
python src/ingest/03_prodigal.py \\
    -i "data/$DATASET/stage00_prepared" \\
    -o "data/$DATASET/stage03_prodigal" \\
    --max-workers "$SLURM_CPUS_ON_NODE"

echo "[$(date)] Done. stage00 + stage03 ready under data/$DATASET/"
"""


def _astra_script(spec: PipelineSpec) -> str:
    """Stage 04 Astra array — one task per DB, --exclusive, --force.

    Per CLAUDE.md: --force is REQUIRED for parallel array submission so
    concurrent tasks don't see stage04_astra/ exist and bail. With --force,
    only the specific per-DB subdir is wiped — not the parent — so concurrent
    tasks coexist cleanly.
    """
    dbs = list(spec.astra_dbs)
    n = len(dbs)
    dbs_bash = " ".join(dbs)
    head = _header(
        spec, job_name=f"{spec.dataset_name[:12]}_astra",
        partition=spec.cpu_partition, log_basename="02_astra_%A_%a.out",
        extra_sbatch=[f"#SBATCH --array=0-{n-1}"],
    )
    return f"""{head}

DBS=({dbs_bash})
DB="${{DBS[$SLURM_ARRAY_TASK_ID]}}"

echo "[$(date)] Stage 04 / Astra : DB=$DB (task $SLURM_ARRAY_TASK_ID)"
echo "  threads: $SLURM_CPUS_ON_NODE"

python src/ingest/04_astra_scan.py \\
    -i "data/$DATASET/stage03_prodigal" \\
    -o "data/$DATASET/stage04_astra" \\
    -d "$DB" \\
    -t "$SLURM_CPUS_ON_NODE" \\
    --force

echo "[$(date)] Done with $DB"
"""


def _crispr_script(spec: PipelineSpec) -> str:
    """Stage 05c MinCED CRISPR. Single-threaded but submitted via SLURM so
    Stage 07 can depend on it cleanly via --dependency=afterok."""
    head = _header(
        spec, job_name=f"{spec.dataset_name[:12]}_crispr",
        partition=spec.cpu_partition, log_basename="03_crispr_%j.out",
    )
    return f"""{head}

echo "[$(date)] Stage 05c / MinCED CRISPR on $DATASET"

python src/ingest/minced_crispr.py \\
    -i "data/$DATASET/stage00_prepared" \\
    -o "data/$DATASET/stage05c_crispr"

echo "[$(date)] Done."
"""


def _build_kb_script(spec: PipelineSpec) -> str:
    """Stage 07 build_knowledge_base. Depends on Astra + CRISPR completing."""
    head = _header(
        spec, job_name=f"{spec.dataset_name[:12]}_kb",
        partition=spec.cpu_partition, log_basename="04_build_kb_%j.out",
    )
    return f"""{head}

echo "[$(date)] Stage 07 / build_knowledge_base on $DATASET"

python src/ingest/07_build_knowledge_base.py \\
    --data-dir "data/$DATASET" \\
    --output "data/$DATASET/sharur.duckdb" \\
    --force

echo "[$(date)] Done. sharur.duckdb built."
"""


def _embeddings_script(spec: PipelineSpec) -> str:
    """Stage 06 ESM2 embeddings. GPU job; user opts in by submitting separately."""
    head = _header(
        spec, job_name=f"{spec.dataset_name[:12]}_esm2",
        partition=spec.gpu_partition, log_basename="05_embeddings_%j.out",
        extra_sbatch=["#SBATCH --gres=gpu:1"],
    )
    return f"""{head}

# Workaround for xet backend (predecessor_session note 7); harmless if unset
export HF_HUB_ENABLE_XET=0

echo "[$(date)] Stage 06 / ESM2 embeddings on $DATASET"

python src/ingest/06_esm2_embeddings.py \\
    "data/$DATASET/stage03_prodigal" \\
    "data/$DATASET/embeddings" \\
    --force

echo "[$(date)] Done."
"""


def _submit_script(spec: PipelineSpec) -> str:
    """The submit driver: chains jobs with --dependency=afterok per CLAUDE.md.

    Omits any stage whose skip_* flag is True. Embeddings is always emitted
    as a one-off (not chained) because most runs don't need it immediately.
    """
    pipe = spec.pipeline_dir

    lines = [
        "#!/bin/bash",
        f"# Submit the {spec.dataset_name} ingest pipeline as a dependency chain.",
        "# CLAUDE.md SLURM rules: --exclusive, no --time, --dependency=afterok",
        "# (no polling). Emitted by " + spec.emitted_by + ".",
        "",
        "set -euo pipefail",
        "",
        f"PIPE={pipe}",
        "",
        f'echo "Submitting {spec.dataset_name} ingest chain..."',
        "",
        "JOB_A=$(sbatch --parsable \"$PIPE/01_prep_genes.sh\")",
        'echo "  A. prep+prodigal       : $JOB_A"',
        "",
    ]
    after_a_jobs = []
    if not spec.skip_astra:
        lines += [
            "JOB_B=$(sbatch --parsable --dependency=afterok:$JOB_A \"$PIPE/02_astra_array.sh\")",
            'echo "  B. astra array (0-' + str(len(spec.astra_dbs) - 1) + ')   : $JOB_B  (after A)"',
            "",
        ]
        after_a_jobs.append("$JOB_B")
    if not spec.skip_crispr:
        lines += [
            "JOB_C=$(sbatch --parsable --dependency=afterok:$JOB_A \"$PIPE/03_crispr.sh\")",
            'echo "  C. minced crispr       : $JOB_C  (after A)"',
            "",
        ]
        after_a_jobs.append("$JOB_C")

    dep = ":".join(after_a_jobs) if after_a_jobs else "$JOB_A"
    lines += [
        f"JOB_D=$(sbatch --parsable --dependency=afterok:{dep} \"$PIPE/04_build_kb.sh\")",
        f'echo "  D. build_kb            : $JOB_D  (after {",".join(j.lstrip("$") for j in after_a_jobs) or "A"})"',
        "",
    ]

    # Embeddings: emit a fire-once line, not chained, in a comment.
    if not spec.skip_embeddings:
        lines += [
            "JOB_E=$(sbatch --parsable --dependency=afterok:$JOB_D \"$PIPE/05_embeddings.sh\")",
            'echo "  E. esm2 embeddings     : $JOB_E  (after D)"',
            "",
        ]
    else:
        lines += [
            "# Stage 06 (ESM2 embeddings) NOT chained by default. Fire it manually:",
            f'#   sbatch --dependency=afterok:$JOB_D "$PIPE/05_embeddings.sh"',
            "",
        ]

    lines += [
        "echo",
        'echo "All jobs queued."',
        f'echo "Logs: {spec.log_dir}/"',
    ]
    return "\n".join(lines)


def emit_pipeline(spec: PipelineSpec) -> dict[str, Path]:
    """Write all scripts. Returns mapping of logical name -> path."""
    spec.pipeline_dir.mkdir(parents=True, exist_ok=True)
    spec.log_dir.mkdir(parents=True, exist_ok=True)

    written: dict[str, Path] = {}

    def write(name: str, content: str) -> None:
        path = spec.pipeline_dir / name
        path.write_text(content if content.endswith("\n") else content + "\n")
        os.chmod(path, 0o755)
        written[name] = path

    write("01_prep_genes.sh", _prep_script(spec))
    if not spec.skip_astra:
        write("02_astra_array.sh", _astra_script(spec))
    if not spec.skip_crispr:
        write("03_crispr.sh", _crispr_script(spec))
    write("04_build_kb.sh", _build_kb_script(spec))
    write("05_embeddings.sh", _embeddings_script(spec))
    write("submit.sh", _submit_script(spec))

    spec.files = written
    return written


__all__ = ["PipelineSpec", "emit_pipeline"]
