#!/usr/bin/env python3
"""
Primary command-line interface for the Sharur staged ingest pipeline.

Use this module through the packaged `sharur-ingest` entrypoint for standard
dataset ingestion, dry-run planning, and local orchestration. The supported
documentation entry points are:

- QUICKSTART.md
- src/ingest/README.md

Two modes are available:
- tools (default): call the existing staged ingest scripts
- fast: synthesize minimal stage outputs for local smoke tests

The canonical scientific workflow still runs through the stage scripts
themselves, especially 00_prepare_inputs.py, 04_astra_scan.py,
minced_crispr.py, 07_build_knowledge_base.py, and 06_esm2_embeddings.py.
"""

from __future__ import annotations

import contextlib
import shlex
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console

from sharur.ingest.dag import (
    IngestDAG,
    StageNode,
    render_slurm_script,
    stage_runner_command,
)
from sharur.ingest.resources import ResourceProfile, resolve_resource_profile
from sharur.ingest.stage_runner import execute_stage
from sharur.ops.ledger import RunLedger


console = Console()
app = typer.Typer(
    no_args_is_help=True,
    add_completion=False,
    help=(
        "Primary CLI for the staged Sharur ingest pipeline. "
        "Use QUICKSTART.md for the default workflow and src/ingest/README.md "
        "for manual stage reference."
    ),
)


def _stage_script_dir() -> Path:
    """Resolve ingest stages from a wheel or the editable source tree."""
    packaged = Path(__file__).resolve().parent / "ingest" / "stages"
    if packaged.is_dir():
        return packaged
    source_tree = Path(__file__).resolve().parents[1] / "src" / "ingest"
    if source_tree.is_dir():
        return source_tree
    raise FileNotFoundError(
        "Sharur ingest stage scripts are missing from this installation; "
        "reinstall Sharur from a current wheel or editable checkout"
    )


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
    def from_root(cls, root: Path) -> StagePaths:
        return cls(
            stage00=root / "stage00_prepared",
            stage01=root / "stage01_quast",
            stage02=root / "stage02_dfast_qc",
            stage03=root / "stage03_prodigal",
            stage04=root / "stage04_astra",
            stage05a=root / "stage05a_gecco",
            stage05b=root / "stage05b_dbcan",
            stage05c=root / "stage05c_crispr",
            stage06=root / "embeddings",
        )


def _build_tools_dag(
    *,
    input_dir: Path,
    data_dir: Path,
    output: Path,
    stages: StagePaths,
    profile: ResourceProfile,
    skip_quast: bool,
    skip_dfast: bool,
    skip_prodigal: bool,
    skip_astra: bool,
    skip_gecco: bool,
    skip_dbcan: bool,
    skip_crispr: bool,
    skip_embeddings: bool,
    enable_cazymes: bool,
) -> IngestDAG:
    """Build the selected stage graph and its concrete resource-tuned commands."""
    stage_dir = _stage_script_dir()
    nodes: list[StageNode] = []

    def add(
        stage_id: str,
        label: str,
        command: list[str],
        *,
        dependencies: tuple[str, ...] = (),
        outputs: tuple[Path, ...] = (),
        inputs: tuple[Path, ...] = (),
        cleanup: tuple[Path, ...] = (),
    ) -> None:
        nodes.append(
            StageNode(
                stage_id=stage_id,
                label=label,
                command=command,
                dependencies=tuple(
                    dependency
                    for dependency in dependencies
                    if any(node.stage_id == dependency for node in nodes)
                ),
                required_outputs=outputs,
                input_paths=inputs,
                cleanup_paths=cleanup,
                resource=profile.request(stage_id),
            )
        )

    add(
        "00",
        "Prepare inputs",
        [
            str(stage_dir / "00_prepare_inputs.py"),
            "-i",
            str(input_dir),
            "-o",
            str(stages.stage00),
            "--force",
        ],
        outputs=(stages.stage00 / "processing_manifest.json",),
        inputs=(input_dir,),
    )
    if not skip_quast:
        add(
            "01",
            "QUAST",
            [
                str(stage_dir / "01_run_quast.py"),
                "-i",
                str(stages.stage00),
                "-o",
                str(stages.stage01),
                "--threads-per-genome",
                "1",
                "--max-workers",
                str(profile.max_workers),
                "--force",
            ],
            dependencies=("00",),
            outputs=(stages.stage01 / "processing_manifest.json",),
        )
    if not skip_dfast:
        add(
            "02",
            "DFAST QC",
            [
                str(stage_dir / "02_dfast_qc.py"),
                "-i",
                str(stages.stage00),
                "-o",
                str(stages.stage02),
                "--threads",
                "1",
                "--max-workers",
                str(profile.max_workers),
                "--force",
            ],
            dependencies=("00",),
            outputs=(stages.stage02 / "processing_manifest.json",),
        )
    if not skip_prodigal:
        add(
            "03",
            "Prodigal",
            [
                str(stage_dir / "03_prodigal.py"),
                "-i",
                str(stages.stage00),
                "-o",
                str(stages.stage03),
                "--max-workers",
                str(profile.max_workers),
                "--force",
            ],
            dependencies=("00",),
            outputs=(
                stages.stage03 / "processing_manifest.json",
                stages.stage03 / "genomes" / "all_protein_symlinks",
            ),
        )
    if not skip_astra:
        add(
            "04",
            "Astra annotation",
            [
                str(stage_dir / "04_astra_scan.py"),
                "-i",
                str(stages.stage03),
                "-o",
                str(stages.stage04),
                "--threads",
                str(profile.annotation_threads),
                "--force",
            ],
            dependencies=("03",),
            outputs=(stages.stage04 / "processing_manifest.json",),
            inputs=(() if not skip_prodigal else (stages.stage03 / "processing_manifest.json",)),
        )
    if not skip_gecco:
        add(
            "05a",
            "GECCO",
            [
                str(stage_dir / "gecco_bgc.py"),
                str(stages.stage00),
                str(stages.stage05a),
                "--threads",
                "1",
                "--max-workers",
                str(min(profile.max_workers, 8)),
                "--force",
            ],
            dependencies=("00",),
            outputs=(
                stages.stage05a / "processing_manifest.json",
                stages.stage05a / "combined_bgc_data.json",
            ),
        )
    if not skip_dbcan:
        add(
            "05b",
            "Legacy dbCAN",
            [
                str(stage_dir / "dbcan_cazyme.py"),
                "--input-dir",
                str(stages.stage03 / "genomes" / "all_protein_symlinks"),
                "--output-dir",
                str(stages.stage05b),
                "--max-workers",
                str(min(profile.max_workers, 4)),
                "--threads-per-job",
                "1",
            ],
            dependencies=("03",),
            outputs=(stages.stage05b / "processing_manifest.json",),
            inputs=(() if not skip_prodigal else (stages.stage03 / "processing_manifest.json",)),
            cleanup=(stages.stage05b,),
        )
    if not skip_crispr:
        add(
            "05c",
            "MinCED arrays",
            [
                str(stage_dir / "minced_crispr.py"),
                "--input-dir",
                str(stages.stage00),
                "--output-dir",
                str(stages.stage05c),
                "--force",
            ],
            dependencies=("00",),
            outputs=(stages.stage05c / "processing_manifest.json",),
        )

    upstream_ids = tuple(node.stage_id for node in nodes)
    builder_command = [
        str(stage_dir / "07_build_knowledge_base.py"),
        "--data-dir",
        str(data_dir),
        "--output",
        str(output),
        "--force",
    ]
    if enable_cazymes:
        builder_command.append("--enable-cazymes")
    skipped_inputs = tuple(
        path
        for enabled, path in (
            (skip_quast, stages.stage01 / "processing_manifest.json"),
            (skip_dfast, stages.stage02 / "processing_manifest.json"),
            (skip_prodigal, stages.stage03 / "processing_manifest.json"),
            (skip_astra, stages.stage04 / "processing_manifest.json"),
            (skip_gecco, stages.stage05a / "processing_manifest.json"),
            (skip_dbcan, stages.stage05b / "processing_manifest.json"),
            (skip_crispr, stages.stage05c / "processing_manifest.json"),
        )
        if enabled and path.exists()
    )
    add(
        "07",
        "Build knowledge base",
        builder_command,
        dependencies=upstream_ids,
        outputs=(output,),
        inputs=skipped_inputs,
    )
    if not skip_embeddings:
        request = profile.request("06")
        add(
            "06",
            "ESM2 embeddings",
            [
                str(stage_dir / "06_esm2_embeddings.py"),
                "--device",
                request.accelerator,
                "--skip-index",
                "--force",
                str(stages.stage03),
                str(stages.stage06),
            ],
            dependencies=("03", "07"),
            outputs=(
                stages.stage06 / "protein_embeddings.h5",
                stages.stage06 / "embedding_manifest.json",
            ),
            inputs=(() if not skip_prodigal else (stages.stage03 / "processing_manifest.json",)),
        )
        add(
            "06i",
            "Persistent FAISS index",
            [
                str(Path(__file__).resolve().parent / "ingest" / "vector_index_runner.py"),
                "--embeddings",
                str(stages.stage06 / "protein_embeddings.h5"),
                "--threads",
                str(profile.index_threads),
            ],
            dependencies=("06",),
            outputs=(stages.stage06 / "protein_embeddings.index.json",),
        )

    dag = IngestDAG(nodes)
    dag.compute_signatures(profile)
    return dag


def _execute_local_dag(
    dag: IngestDAG,
    *,
    data_dir: Path,
    profile: ResourceProfile,
    resume: bool,
    force: bool,
    idempotency_key: str | None,
) -> None:
    ledger_path = data_dir / "sharur_ops.db"
    ledger = RunLedger(ledger_path)
    ledger.recover_stale_runs(stale_after_seconds=300)
    run_id = ledger.create_run(
        "ingest",
        data_dir,
        created_by="sharur-ingest",
        config={
            "profile": profile.to_dict(),
            "resume": resume,
            "force": force,
            "stages": [node.stage_id for node in dag.ordered()],
            "stage_signatures": {node.stage_id: node.signature for node in dag.ordered()},
        },
        idempotency_key=idempotency_key,
    )
    existing_run = ledger.get_run(run_id)
    if existing_run["status"] == "complete":
        console.print(f"[green]idempotent run already complete[/green]: {run_id}")
        ledger.close()
        return
    if existing_run["status"] == "failed":
        ledger.close()
        raise RuntimeError(f"Idempotent run {run_id} already failed; use a new key to retry")
    if existing_run["status"] in {"submitted", "running"}:
        console.print(f"[green]idempotent run already {existing_run['status']}[/green]: {run_id}")
        ledger.close()
        return
    try:
        ledger.start_run(run_id)
    except ValueError:
        raced = ledger.get_run(run_id)
        if raced["status"] in {"submitted", "running", "complete"}:
            console.print(f"[green]idempotent run already {raced['status']}[/green]: {run_id}")
            ledger.close()
            return
        ledger.close()
        raise
    try:
        executed_stages: set[str] = set()
        for node in dag.ordered():
            source = None
            dependency_executed = bool(set(node.dependencies) & executed_stages)
            if resume and not force and not dependency_executed:
                source = ledger.find_reusable_stage(
                    data_dir,
                    "ingest",
                    node.stage_id,
                    node.signature,
                )
            if source is not None and node.outputs_match(source.get("outputs")):
                ledger.reuse_stage(run_id, node.stage_id, node.signature, source)
                console.print(
                    f"[green]reuse[/green] {node.stage_id} {node.label} "
                    f"(from run {source['run_id']})"
                )
                continue

            console.print(
                f"[cyan]run[/cyan] {node.stage_id} {node.label} "
                f"[dim]({profile.name}, {node.resource.cpus} CPU)[/dim]"
            )
            execute_stage(
                ledger=ledger,
                run_id=run_id,
                stage_id=node.stage_id,
                signature=node.signature,
                command=[sys.executable, *node.command],
                outputs=list(node.required_outputs),
                cleanup=list(node.cleanup_paths),
                resource=node.resource,
            )
            executed_stages.add(node.stage_id)
        ledger.complete_run(
            run_id,
            result={"stages": [node.stage_id for node in dag.ordered()]},
        )
    except Exception:
        # execute_stage records both stage and run failure.
        raise
    finally:
        ledger.close()


def _emit_slurm_bundle(
    dag: IngestDAG,
    *,
    data_dir: Path,
    profile: ResourceProfile,
    resume: bool,
    force: bool,
    idempotency_key: str | None,
    submit: bool,
) -> Path:
    bundle_dir = data_dir / "slurm"
    log_dir = bundle_dir / "logs"
    bundle_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    submit_path = bundle_dir / "submit.sh"
    ledger_path = data_dir / "sharur_ops.db"
    ledger = RunLedger(ledger_path)
    ledger.recover_stale_runs(stale_after_seconds=300)
    run_id = ledger.create_run(
        "ingest",
        data_dir,
        created_by="sharur-ingest",
        config={
            "profile": profile.to_dict(),
            "resume": resume,
            "force": force,
            "submission": "slurm",
            "stages": [node.stage_id for node in dag.ordered()],
            "stage_signatures": {node.stage_id: node.signature for node in dag.ordered()},
        },
        idempotency_key=idempotency_key,
    )
    existing_run = ledger.get_run(run_id)
    if existing_run["status"] == "complete":
        if not submit_path.is_file():
            submit_path.write_text(
                "#!/usr/bin/env bash\nset -euo pipefail\n"
                "echo 'idempotent ingest run already complete; nothing to submit'\n"
            )
            submit_path.chmod(0o755)
        console.print(f"[green]idempotent run already complete[/green]: {run_id}")
        ledger.close()
        return submit_path
    if existing_run["status"] == "failed":
        ledger.close()
        raise RuntimeError(f"Idempotent run {run_id} already failed; use a new key to retry")
    if existing_run["status"] in {"submitted", "running"}:
        if not submit_path.is_file():
            ledger.close()
            raise RuntimeError(
                f"Run {run_id} is {existing_run['status']} but its submit script is missing"
            )
        console.print(f"[green]idempotent run already {existing_run['status']}[/green]: {run_id}")
        ledger.close()
        return submit_path

    active: list[StageNode] = []
    active_ids: set[str] = set()
    try:
        for node in dag.ordered():
            existing_stage = ledger.latest_stage(run_id, node.stage_id)
            if existing_stage is not None:
                if existing_stage["signature"] != node.signature:
                    raise RuntimeError(
                        f"Run {run_id} already records stage {node.stage_id} "
                        "with a different signature"
                    )
                if existing_stage["status"] in {"complete", "reused"}:
                    if not node.outputs_match(existing_stage.get("outputs")):
                        raise RuntimeError(
                            f"Idempotent run {run_id} has a successful stage "
                            f"{node.stage_id}, but its outputs no longer match"
                        )
                    continue
                raise RuntimeError(
                    f"Run {run_id} already records stage {node.stage_id} "
                    f"as {existing_stage['status']}"
                )
            source = None
            dependency_active = bool(set(node.dependencies) & active_ids)
            if resume and not force and not dependency_active:
                source = ledger.find_reusable_stage(
                    data_dir,
                    "ingest",
                    node.stage_id,
                    node.signature,
                )
            if source is not None and node.outputs_match(source.get("outputs")):
                ledger.reuse_stage(run_id, node.stage_id, node.signature, source)
            else:
                active.append(node)
                active_ids.add(node.stage_id)

        if not active:
            ledger.complete_run(run_id, result={"stages": "all_reused"})
            submit_path.write_text(
                "#!/usr/bin/env bash\nset -euo pipefail\n"
                "echo 'all ingest stages reused; nothing to submit'\n"
            )
            submit_path.chmod(0o755)
            return submit_path

        script_paths: dict[str, Path] = {}
        for index, node in enumerate(active):
            runner = stage_runner_command(
                node,
                ledger_path=ledger_path,
                run_id=run_id,
                complete_run=index == len(active) - 1,
            )
            suffix = "sbatch" if node.resource.executor == "slurm" else "local.sh"
            script_path = bundle_dir / f"{node.stage_id}.{suffix}"
            if node.resource.executor == "slurm":
                content = render_slurm_script(node, runner, log_dir=log_dir)
            else:
                content = "#!/usr/bin/env bash\nset -euo pipefail\n" + shlex.join(runner) + "\n"
            script_path.write_text(content)
            script_path.chmod(0o755)
            script_paths[node.stage_id] = script_path

        lines = ["#!/usr/bin/env bash", "set -euo pipefail"]
        for node in active:
            path = script_paths[node.stage_id]
            if node.resource.executor == "local":
                lines.append(f"bash {shlex.quote(str(path))}")
                continue
            dependency_vars = [
                f"${{jid_{dependency}}}"
                for dependency in node.dependencies
                if dependency in active_ids and dag.nodes[dependency].resource.executor == "slurm"
            ]
            dependency_arg = (
                f"--dependency=afterok:{':'.join(dependency_vars)} " if dependency_vars else ""
            )
            lines.append(
                f"jid_{node.stage_id}=$(sbatch --parsable {dependency_arg}{shlex.quote(str(path))})"
            )
            lines.append(f"echo 'submitted {node.stage_id} as '${{jid_{node.stage_id}}}")
        submit_path.write_text("\n".join(lines) + "\n")
        submit_path.chmod(0o755)
        if submit:
            ledger.submit_run(run_id)
    finally:
        ledger.close()

    console.print(f"[green]SLURM bundle ready[/green]: {submit_path}")
    if submit:
        try:
            subprocess.run(["bash", str(submit_path)], check=True)
        except Exception as exc:
            failure_ledger = RunLedger(ledger_path)
            try:
                with contextlib.suppress(ValueError):
                    failure_ledger.fail_run(
                        run_id,
                        f"scheduler submission failed: {type(exc).__name__}: {exc}",
                    )
            finally:
                failure_ledger.close()
            raise
    return submit_path


@app.command()
def run(
    input_dir: Annotated[
        Path,
        typer.Option("--input-dir", "-i", help="Directory of genome FASTAs (.fna/.fa/.fasta)"),
    ] = Path("dummy_dataset"),
    data_dir: Annotated[
        Path,
        typer.Option("--data-dir", "-d", help="Dataset root where stageXX outputs will be written"),
    ] = Path("data"),
    output: Annotated[
        Path,
        typer.Option(
            "--output", "-o", help="Destination DuckDB path, usually DATASET/sharur.duckdb"
        ),
    ] = Path("data/sharur.duckdb"),
    mode: Annotated[
        str,
        typer.Option(
            "--mode",
            "-m",
            help="tools = run the standard staged pipeline; fast = synthetic smoke mode",
        ),
    ] = "tools",
    force: Annotated[
        bool,
        typer.Option("--force", help="Remove existing stage outputs and rebuild"),
    ] = False,
    skip_quast: Annotated[
        bool,
        typer.Option(
            "--skip-quast/--with-quast",
            help="Skip optional Stage 01 by default; use --with-quast to enable it",
        ),
    ] = True,
    skip_dfast: Annotated[
        bool,
        typer.Option(
            "--skip-dfast/--with-dfast",
            help="Skip optional Stage 02 by default; use --with-dfast to enable it",
        ),
    ] = True,
    skip_prodigal: Annotated[bool, typer.Option(help="Skip standard Stage 03 (Prodigal)")] = False,
    skip_astra: Annotated[
        bool,
        typer.Option(help="Skip standard Stage 04 (Astra via 04_astra_scan.py)"),
    ] = False,
    skip_gecco: Annotated[
        bool,
        typer.Option(
            "--skip-gecco/--with-gecco",
            help="Skip optional Stage 05a by default; use --with-gecco to enable it",
        ),
    ] = True,
    skip_dbcan: Annotated[
        bool,
        typer.Option(
            "--skip-dbcan/--with-legacy-dbcan",
            help="Skip deprecated Stage 05b by default; use --with-legacy-dbcan to enable it",
        ),
    ] = True,
    enable_cazymes: Annotated[
        bool,
        typer.Option(
            "--enable-cazymes",
            help="Run the optional Stage 07 dbCAN three-tool consensus classifier",
        ),
    ] = False,
    skip_crispr: Annotated[
        bool,
        typer.Option(help="Skip standard Stage 05c (CRISPR arrays via minced_crispr.py)"),
    ] = False,
    skip_embeddings: Annotated[
        bool,
        typer.Option(help="Skip Stage 06 post-build embeddings (06_esm2_embeddings.py)"),
    ] = False,
    profile: Annotated[
        str,
        typer.Option(
            "--profile",
            help="Execution profile: auto, local (CPU), mps, or slurm",
        ),
    ] = "auto",
    resume: Annotated[
        bool,
        typer.Option(
            "--resume/--no-resume",
            help="Reuse only ledger-verified stages with matching signatures and outputs",
        ),
    ] = True,
    submit_slurm: Annotated[
        bool,
        typer.Option(
            "--submit-slurm",
            help="Submit the generated SLURM bundle; otherwise only write it",
        ),
    ] = False,
    run_idempotency_key: Annotated[
        str | None,
        typer.Option(
            "--run-idempotency-key",
            help="Optional caller-stable key that deduplicates identical run creation",
        ),
    ] = None,
    dry_run: Annotated[
        bool,
        typer.Option(help="Print the planned stage commands without executing them"),
    ] = False,
) -> list[list[str]] | None:
    """
    Run the staged ingest pipeline from the primary CLI entrypoint.

    `mode=tools` is the default and mirrors the documented standard pipeline.
    `mode=fast` exists only for local smoke tests and synthetic fixtures.

    For the user-facing workflow, see QUICKSTART.md or src/ingest/README.md.
    """
    try:
        resource_profile = resolve_resource_profile(profile)
    except ValueError as exc:
        raise typer.BadParameter(str(exc)) from exc
    if submit_slurm and resource_profile.name != "slurm":
        raise typer.BadParameter("--submit-slurm requires --profile slurm")

    input_dir = input_dir.expanduser().resolve()
    data_dir = data_dir.expanduser().resolve()
    output = output.expanduser().resolve()
    stages = StagePaths.from_root(data_dir)
    if mode == "fast":
        if resource_profile.name == "slurm":
            raise typer.BadParameter("fast smoke mode does not support the slurm profile")
        stage_dir = _stage_script_dir()
        synthetic_command = [
            str(Path(__file__).resolve().parent / "ingest" / "synthetic_stage.py"),
            "--input-dir",
            str(input_dir),
            "--data-dir",
            str(data_dir),
        ]
        fixture_node = StageNode(
            stage_id="00",
            label="Synthesize smoke-test fixtures",
            command=synthetic_command,
            required_outputs=(
                stages.stage02 / "processing_manifest.json",
                stages.stage03 / "genomes",
                stages.stage04 / "synthetic_hits_df.tsv",
                stages.stage05a / "combined_bgc_data.json",
                stages.stage05b / "synthetic_cazyme_results.json",
                stages.stage05c / "synthetic_crispr_arrays.json",
            ),
            input_paths=(input_dir,),
            resource=resource_profile.request("00"),
        )
        builder_command = [
            str(stage_dir / "07_build_knowledge_base.py"),
            "--data-dir",
            str(data_dir),
            "--output",
            str(output),
            "--force",
        ]
        if enable_cazymes:
            builder_command.append("--enable-cazymes")
        builder_node = StageNode(
            stage_id="07",
            label="Build synthetic knowledge base",
            command=builder_command,
            dependencies=("00",),
            required_outputs=(output,),
            resource=resource_profile.request("07"),
        )
        dag = IngestDAG([fixture_node, builder_node])
        dag.compute_signatures(resource_profile)
    elif mode == "tools":
        dag = _build_tools_dag(
            input_dir=input_dir,
            data_dir=data_dir,
            output=output,
            stages=stages,
            profile=resource_profile,
            skip_quast=skip_quast,
            skip_dfast=skip_dfast,
            skip_prodigal=skip_prodigal,
            skip_astra=skip_astra,
            skip_gecco=skip_gecco,
            skip_dbcan=skip_dbcan,
            skip_crispr=skip_crispr,
            skip_embeddings=skip_embeddings,
            enable_cazymes=enable_cazymes,
        )
    else:
        raise typer.BadParameter("mode must be 'tools' or 'fast'")

    planned = [node.command for node in dag.ordered()]
    if dry_run:
        console.print(f"[cyan]Dry run DAG ({resource_profile.name} profile):[/cyan]")
        for node in dag.ordered():
            dependencies = ",".join(node.dependencies) or "none"
            request = node.resource
            console.print(
                f"  [{node.stage_id}] deps={dependencies} "
                f"resources={request.cpus}cpu/{request.memory_gb}GB/"
                f"{request.accelerator}: {' '.join(node.command)}"
            )
        return planned

    data_dir.mkdir(parents=True, exist_ok=True)
    if resource_profile.name == "slurm":
        _emit_slurm_bundle(
            dag,
            data_dir=data_dir,
            profile=resource_profile,
            resume=resume,
            force=force,
            idempotency_key=run_idempotency_key,
            submit=submit_slurm,
        )
        return planned

    _execute_local_dag(
        dag,
        data_dir=data_dir,
        profile=resource_profile,
        resume=resume,
        force=force,
        idempotency_key=run_idempotency_key,
    )
    return None


def main() -> None:
    app()


__all__ = ["app", "main", "run"]


if __name__ == "__main__":
    main()
