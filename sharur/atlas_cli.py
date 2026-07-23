"""Command-line entry point for deterministic Atlas campaigns."""

from __future__ import annotations

import json
from pathlib import Path  # noqa: TC003 - Typer resolves runtime annotations
from typing import Annotated

import typer

from sharur.atlas import (
    DEFAULT_CHECKPOINT_INTERVAL_CONTIGS,
    DEFAULT_PACKET_PROTEINS,
    build_atlas_plan,
    enqueue_atlas_plan,
    verify_atlas_coverage,
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
    seal: Annotated[
        Path | None,
        typer.Option("--seal", exists=True, dir_okay=False),
    ] = None,
    packet_proteins: Annotated[
        int,
        typer.Option("--packet-proteins", min=1, max=250),
    ] = DEFAULT_PACKET_PROTEINS,
    checkpoint_interval_contigs: Annotated[
        int,
        typer.Option("--checkpoint-interval-contigs", min=1),
    ] = DEFAULT_CHECKPOINT_INTERVAL_CONTIGS,
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
        packet_protein_limit=packet_proteins,
        checkpoint_interval_contigs=checkpoint_interval_contigs,
        threads=threads,
        verify_seal=verify_seal,
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
