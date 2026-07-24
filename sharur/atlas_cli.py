"""Command-line entry point for deterministic Atlas campaigns."""

from __future__ import annotations

import json
from pathlib import Path  # noqa: TC003 - Typer resolves runtime annotations
from typing import Annotated

import typer

from sharur.atlas import (
    DEFAULT_CENSUS_MAX_TEMP_SIZE,
    DEFAULT_CENSUS_MEMORY_LIMIT,
    DEFAULT_CENSUS_WORKERS,
    DEFAULT_CHECKPOINT_INTERVAL_FRAMES,
    DEFAULT_QUERY_RESULT_LIMIT_BYTES,
    build_atlas_packet_census,
    build_atlas_plan,
    enqueue_atlas_plan,
    verify_atlas_coverage,
    verify_atlas_packet_census,
)
from sharur.operators.contigs import (
    MAX_GENOME_PACKET_CONTIGS,
    MAX_GENOME_PACKET_PROTEINS,
)
from sharur.ops.client import SharurOps


app = typer.Typer(
    no_args_is_help=True,
    add_completion=False,
    help="Plan, enqueue, and verify exhaustive genome-owned Atlas reading.",
)


@app.command("plan")
def plan_command(
    db: Annotated[Path, typer.Option("--db", exists=True, dir_okay=False)],
    output_dir: Annotated[Path, typer.Option("--output-dir")],
    packet_bytes: Annotated[
        int,
        typer.Option(
            "--packet-bytes",
            min=1_024,
            max=1_500_000,
            help=(
                "Required model-facing canonical payload budget; derive it "
                "from the executor context budget."
            ),
        ),
    ],
    seal: Annotated[
        Path | None,
        typer.Option("--seal", exists=True, dir_okay=False),
    ] = None,
    packet_contigs: Annotated[
        int | None,
        typer.Option(
            "--packet-contigs",
            min=1,
            max=MAX_GENOME_PACKET_CONTIGS,
            help=(
                "Diagnostic hard cap; omitted values use the byte-proportional "
                "schema ceiling."
            ),
        ),
    ] = None,
    packet_proteins: Annotated[
        int | None,
        typer.Option(
            "--packet-proteins",
            min=1,
            max=MAX_GENOME_PACKET_PROTEINS,
            help=(
                "Diagnostic hard cap; omitted values use the byte-proportional "
                "schema ceiling."
            ),
        ),
    ] = None,
    all_annotations: Annotated[
        bool,
        typer.Option("--all-annotations/--top-annotation-only"),
    ] = True,
    calibration_genomes: Annotated[
        int | None,
        typer.Option(
            "--calibration-genomes",
            min=1,
            help=(
                "Explicit calibration size; omitted values use "
                "ceil(sqrt(live bin count))."
            ),
        ),
    ] = None,
    checkpoint_interval_frames: Annotated[
        int,
        typer.Option("--checkpoint-interval-frames", min=1),
    ] = DEFAULT_CHECKPOINT_INTERVAL_FRAMES,
    query_result_bytes: Annotated[
        int,
        typer.Option("--query-result-bytes", min=4_097),
    ] = DEFAULT_QUERY_RESULT_LIMIT_BYTES,
    threads: Annotated[int, typer.Option("--threads", min=1)] = 4,
    verify_seal: Annotated[
        bool,
        typer.Option("--verify-seal/--skip-seal-verification"),
    ] = True,
) -> None:
    """Build stable one-genome work units from a sealed DuckDB."""
    result = build_atlas_plan(
        db,
        output_dir,
        seal_path=seal,
        packet_contig_limit=packet_contigs,
        packet_protein_limit=packet_proteins,
        packet_model_payload_bytes=packet_bytes,
        packet_all_annotations=all_annotations,
        packet_calibration_genomes=calibration_genomes,
        checkpoint_interval_frames=checkpoint_interval_frames,
        query_result_limit_bytes=query_result_bytes,
        threads=threads,
        verify_seal=verify_seal,
    )
    typer.echo(json.dumps(result, indent=2, sort_keys=True))


@app.command("packet-census")
def packet_census_command(
    plan_dir: Annotated[Path, typer.Option("--plan-dir", exists=True, file_okay=False)],
    output_dir: Annotated[Path | None, typer.Option("--output-dir")] = None,
    workers: Annotated[int, typer.Option("--workers", min=1)] = DEFAULT_CENSUS_WORKERS,
    threads: Annotated[int, typer.Option("--threads", min=1)] = 4,
    memory_limit: Annotated[
        str,
        typer.Option("--memory-limit"),
    ] = DEFAULT_CENSUS_MEMORY_LIMIT,
    temp_directory: Annotated[Path | None, typer.Option("--temp-directory")] = None,
    max_temp_size: Annotated[
        str,
        typer.Option("--max-temp-size"),
    ] = DEFAULT_CENSUS_MAX_TEMP_SIZE,
    resume: Annotated[bool, typer.Option("--resume/--recompute")] = True,
    verify_seal: Annotated[
        bool,
        typer.Option("--verify-seal/--skip-seal-verification"),
    ] = True,
) -> None:
    """Count exact bin-scoped packets and payload sizes with zero model calls."""

    def report_progress(progress: dict) -> None:
        completed = int(progress["completed_units"])
        total = int(progress["total_units"])
        if (
            completed == total
            or completed % 100 == 0
            or progress["status"] == "resumed"
        ):
            percent = 100.0 * completed / total if total else 100.0
            typer.echo(
                (
                    f"packet-census {completed:,}/{total:,} units "
                    f"({percent:.2f}%); "
                    f"{int(progress['cumulative_frame_count']):,} frames"
                ),
                err=True,
            )

    result = build_atlas_packet_census(
        plan_dir,
        output_dir=output_dir,
        workers=workers,
        threads=threads,
        memory_limit=memory_limit,
        temp_directory=temp_directory,
        max_temp_directory_size=max_temp_size,
        resume=resume,
        verify_seal=verify_seal,
        progress_callback=report_progress,
    )
    typer.echo(json.dumps(result, indent=2, sort_keys=True))
    if result["status"] != "complete":
        raise typer.Exit(code=1)


@app.command("verify-packet-census")
def verify_packet_census_command(
    plan_dir: Annotated[Path, typer.Option("--plan-dir", exists=True, file_okay=False)],
    census_dir: Annotated[
        Path | None,
        typer.Option("--census-dir", exists=True, file_okay=False),
    ] = None,
    deep: Annotated[bool, typer.Option("--deep/--summary-only")] = False,
) -> None:
    """Verify the zero-model-call launch gate and optional unit records."""
    result = verify_atlas_packet_census(
        plan_dir,
        census_dir=census_dir,
        deep=deep,
    )
    typer.echo(json.dumps(result, indent=2, sort_keys=True))


@app.command("enqueue")
def enqueue_command(
    plan_dir: Annotated[Path, typer.Option("--plan-dir", exists=True, file_okay=False)],
    query_url: Annotated[str, typer.Option("--query-url")],
    ops_url: Annotated[str, typer.Option("--ops-url")] = "http://localhost:8811",
    agent_id: Annotated[str, typer.Option("--agent-id")] = "atlas-coordinator",
    api_token: Annotated[
        str | None,
        typer.Option("--api-token", envvar="SHARUR_OPS_TOKEN"),
    ] = None,
    priority: Annotated[int, typer.Option("--priority", min=0, max=3)] = 1,
    max_attempts: Annotated[int, typer.Option("--max-attempts", min=1)] = 5,
    lease_seconds: Annotated[int, typer.Option("--lease-seconds", min=1)] = 900,
    scan_execution_profile: Annotated[
        str,
        typer.Option("--scan-execution-profile"),
    ] = "atlas_scan",
) -> None:
    """Create an idempotent Ops campaign and one task per genome."""
    with SharurOps(
        ops_url,
        agent_id=agent_id,
        api_token=api_token,
    ) as ops:
        result = enqueue_atlas_plan(
            plan_dir,
            query_url=query_url,
            ops=ops,
            priority=priority,
            max_attempts=max_attempts,
            lease_seconds=lease_seconds,
            scan_execution_profile=scan_execution_profile,
        )
    typer.echo(json.dumps(result, indent=2, sort_keys=True))


@app.command("verify-coverage")
def verify_coverage_command(
    plan_dir: Annotated[Path, typer.Option("--plan-dir", exists=True, file_okay=False)],
    coverage_dir: Annotated[
        Path | None,
        typer.Option("--coverage-dir", exists=True, file_okay=False),
    ] = None,
) -> None:
    """Validate every assigned genome and exact contig/protein totals."""
    result = verify_atlas_coverage(plan_dir, coverage_dir=coverage_dir)
    typer.echo(json.dumps(result, indent=2, sort_keys=True))
    if result["status"] != "complete":
        raise typer.Exit(code=1)


def main() -> None:
    app()


__all__ = ["app", "main"]


if __name__ == "__main__":
    main()
