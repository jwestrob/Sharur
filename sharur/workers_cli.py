"""CLI for Sharur worker executors (`sharur-worker`)."""

from __future__ import annotations

import logging
import os

import typer

app = typer.Typer(help="Sharur model-worker executors", no_args_is_help=True)


@app.callback()
def _root() -> None:
    """Sharur worker executors.

    A callback keeps subcommands explicit even while only one exists, so
    `atlas-scan` stays a stable invocation as review-tier workers are added.
    """


@app.command("atlas-scan")
def atlas_scan(
    ops_url: str = typer.Option("http://127.0.0.1:8811", help="Sharur Ops base URL"),
    query_url: str = typer.Option("http://127.0.0.1:8812", help="Sharur Query base URL"),
    agent_id: str = typer.Option(..., help="Distinct agent identity for this worker"),
    profile: str = typer.Option("atlas_scan", help="Execution profile from the review policy"),
    campaign_id: str | None = typer.Option(None, help="Restrict claims to one campaign"),
    policy: str | None = typer.Option(None, help="Path to a review policy YAML (default: packaged)"),
    lease_seconds: int = typer.Option(900, min=60, help="Lease duration requested per claim"),
    model_timeout: int = typer.Option(1800, min=60, help="Per-frame model timeout (seconds)"),
    max_tasks: int | None = typer.Option(None, min=1, help="Exit after N completed genomes"),
    max_frames: int | None = typer.Option(None, min=1, help="Stop each genome after N frames (smoke tests)"),
    idle_sleep: float = typer.Option(5.0, help="Seconds to sleep when the queue is empty"),
    dry_run: bool = typer.Option(
        False,
        "--dry-run",
        help="Walk packets and build coverage without calling any model (zero model calls)",
    ),
    verbose: bool = typer.Option(False, "-v", "--verbose"),
) -> None:
    """Claim and execute `atlas_genome_read` tasks.

    `--dry-run` exercises the full claim -> packet -> coverage -> disposition ->
    complete path with zero model calls, which is the right smoke test before
    spending any subscription budget.
    """
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)-7s %(name)s: %(message)s",
    )

    from sharur.workers.atlas_scan import AtlasScanWorker

    worker = AtlasScanWorker(
        ops_url=ops_url,
        query_url=query_url,
        agent_id=agent_id,
        profile=profile,
        policy_path=policy,
        campaign_id=campaign_id,
        ops_token=os.environ.get("SHARUR_OPS_TOKEN"),
        query_token=os.environ.get("SHARUR_QUERY_TOKEN"),
        lease_seconds=lease_seconds,
        model_timeout=model_timeout,
        max_frames=max_frames,
        dry_run=dry_run,
    )
    worker.install_signal_handlers()
    done = worker.run_forever(idle_sleep=idle_sleep, max_tasks=max_tasks)
    typer.echo(f"completed {done} task(s)")


def main() -> None:
    app()


if __name__ == "__main__":
    main()
