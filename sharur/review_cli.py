"""Command-line interface for reduction, routing, verification, and tracing."""

from __future__ import annotations

import json
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Annotated, Any

import typer
import yaml

from sharur.ops.client import SharurOps
from sharur.ops.review_store import content_hash
from sharur.ops.store import OpsStore
from sharur.review.controller import ReviewController
from sharur.review.metrics import review_campaign_metrics
from sharur.review.models import load_review_policy
from sharur.review.reducer import ExactCandidateReducer
from sharur.review.trace import trace_review_subject
from sharur.review.verification import VerificationSpec, run_duckdb_verification


app = typer.Typer(
    no_args_is_help=True,
    add_completion=False,
    help="Operate Sharur's typed, policy-driven scientific review DAG.",
)


def _json(value: Any) -> str:
    return json.dumps(value, indent=2, sort_keys=True, default=str)


def _read_value(value: str) -> Any:
    path = Path(value).expanduser()
    if path.is_file():
        try:
            return yaml.safe_load(path.read_text())
        except (OSError, yaml.YAMLError) as exc:
            raise typer.BadParameter(f"Could not read {path}: {exc}") from exc
    try:
        return json.loads(value)
    except json.JSONDecodeError as exc:
        raise typer.BadParameter(
            "Value must be JSON or a path to a YAML/JSON file"
        ) from exc


@contextmanager
def _ops_access(
    *,
    ops_db: Path | None,
    ops_url: str | None,
    agent_id: str,
    api_token: str | None,
):
    if (ops_db is None) == (ops_url is None):
        raise typer.BadParameter("Supply exactly one of --ops-db or --ops-url")
    if ops_db is not None:
        with OpsStore(ops_db, agent_id=agent_id) as store:
            yield store
    else:
        assert ops_url is not None
        with SharurOps(
            ops_url,
            agent_id=agent_id,
            api_token=api_token,
        ) as client:
            yield client


@app.command("policy-check")
def policy_check_command(
    policy_path: Annotated[
        Path | None,
        typer.Option("--policy", exists=True, dir_okay=False),
    ] = None,
) -> None:
    """Validate a policy and print its immutable execution contract."""

    policy = load_review_policy(policy_path)
    typer.echo(
        _json(
            {
                "name": policy.name,
                "version": policy.version,
                "schema_version": policy.schema_version,
                "policy_hash": policy.policy_hash,
                "profiles": {
                    name: profile.model_dump(mode="json")
                    for name, profile in policy.profiles.items()
                },
                "tiers": {
                    name: tier.model_dump(mode="json")
                    for name, tier in policy.tiers.items()
                },
                "audit": policy.audit.model_dump(mode="json"),
            }
        )
    )


@app.command("reduce")
def reduce_command(
    campaign_id: Annotated[str, typer.Option("--campaign-id")],
    ops_db: Annotated[
        Path | None,
        typer.Option("--ops-db", exists=True, dir_okay=False),
    ] = None,
    ops_url: Annotated[str | None, typer.Option("--ops-url")] = None,
    agent_id: Annotated[str, typer.Option("--agent-id")] = "review-reducer",
    api_token: Annotated[
        str | None,
        typer.Option("--api-token", envvar="SHARUR_OPS_TOKEN"),
    ] = None,
    dataset_id: Annotated[str | None, typer.Option("--dataset-id")] = None,
    candidate_type: Annotated[str | None, typer.Option("--candidate-type")] = None,
    batch_size: Annotated[
        int, typer.Option("--batch-size", min=1, max=10_000)
    ] = 1_000,
) -> None:
    """Reduce exact typed signatures into lossless versioned clusters."""

    with _ops_access(
        ops_db=ops_db,
        ops_url=ops_url,
        agent_id=agent_id,
        api_token=api_token,
    ) as access:
        if isinstance(access, OpsStore):
            result = ExactCandidateReducer().reduce_campaign(
                access,
                campaign_id,
                dataset_id=dataset_id,
                candidate_type=candidate_type,
                batch_size=batch_size,
            ).to_dict()
        else:
            result = access.reduce_candidates(
                campaign_id,
                dataset_id=dataset_id,
                candidate_type=candidate_type,
                batch_size=batch_size,
            )
    typer.echo(_json(result))


@app.command("route")
def route_command(
    campaign_id: Annotated[str, typer.Option("--campaign-id")],
    ops_db: Annotated[
        Path | None,
        typer.Option("--ops-db", exists=True, dir_okay=False),
    ] = None,
    ops_url: Annotated[str | None, typer.Option("--ops-url")] = None,
    agent_id: Annotated[str, typer.Option("--agent-id")] = "review-controller",
    api_token: Annotated[
        str | None,
        typer.Option("--api-token", envvar="SHARUR_OPS_TOKEN"),
    ] = None,
    policy_path: Annotated[
        Path | None,
        typer.Option("--policy", exists=True, dir_okay=False),
    ] = None,
    watch: Annotated[bool, typer.Option("--watch/--once")] = False,
    interval_seconds: Annotated[
        float, typer.Option("--interval-seconds", min=0.1, max=60.0)
    ] = 2.0,
) -> None:
    """Consume durable events and create idempotent review work."""

    policy = load_review_policy(policy_path)
    policy_payload = (
        policy.model_dump(mode="json") if policy_path is not None else None
    )
    with _ops_access(
        ops_db=ops_db,
        ops_url=ops_url,
        agent_id=agent_id,
        api_token=api_token,
    ) as access:
        try:
            while True:
                if isinstance(access, OpsStore):
                    result = ReviewController(access, policy).tick(
                        campaign_id
                    ).to_dict()
                else:
                    result = access.tick_review_controller(
                        campaign_id,
                        policy=policy_payload,
                    )
                typer.echo(_json(result))
                if not watch:
                    break
                time.sleep(interval_seconds)
        except KeyboardInterrupt:
            raise typer.Exit(code=130) from None


@app.command("verify")
def verify_command(
    review_id: Annotated[str, typer.Option("--review-id")],
    claim_key: Annotated[str, typer.Option("--claim-key")],
    db: Annotated[Path, typer.Option("--db", exists=True, dir_okay=False)],
    dataset_id: Annotated[str, typer.Option("--dataset-id")],
    specification: Annotated[
        str,
        typer.Option(
            "--specification",
            help="Inline JSON or a path to a YAML/JSON verification spec.",
        ),
    ],
    expected: Annotated[
        str,
        typer.Option(
            "--expected",
            help="Inline JSON or a path to a YAML/JSON expected value.",
        ),
    ],
    ops_db: Annotated[
        Path | None,
        typer.Option("--ops-db", exists=True, dir_okay=False),
    ] = None,
    ops_url: Annotated[str | None, typer.Option("--ops-url")] = None,
    agent_id: Annotated[str, typer.Option("--agent-id")] = "review-verifier",
    api_token: Annotated[
        str | None,
        typer.Option("--api-token", envvar="SHARUR_OPS_TOKEN"),
    ] = None,
    seal: Annotated[
        Path | None,
        typer.Option("--seal", exists=True, dir_okay=False),
    ] = None,
    verify_seal: Annotated[
        bool,
        typer.Option("--verify-seal/--skip-seal-verification"),
    ] = True,
    threads: Annotated[int, typer.Option("--threads", min=1)] = 1,
    code_commit: Annotated[str | None, typer.Option("--code-commit")] = None,
    supersedes_verification_id: Annotated[
        str | None, typer.Option("--supersedes")
    ] = None,
) -> None:
    """Execute and append one bounded read-only DuckDB verification."""

    spec = VerificationSpec.model_validate(_read_value(specification))
    expected_value = _read_value(expected)
    result = run_duckdb_verification(
        db,
        spec,
        expected_value,
        dataset_id=dataset_id,
        seal_path=seal,
        verify_seal=verify_seal,
        threads=threads,
    )
    payload = {
        "claim_key": claim_key,
        "engine": "duckdb",
        "specification": spec.model_dump(mode="json"),
        "dataset_id": dataset_id,
        "expected": expected_value,
        "status": result.status,
        "actual": result.actual,
        "executed_ts": result.executed_ts,
        "code_commit": code_commit,
        "error": result.error,
        "supersedes_verification_id": supersedes_verification_id,
    }
    payload["idempotency_key"] = "verify:" + content_hash(
        {
            "review_id": review_id,
            **payload,
        }
    )
    with _ops_access(
        ops_db=ops_db,
        ops_url=ops_url,
        agent_id=agent_id,
        api_token=api_token,
    ) as access:
        if isinstance(access, OpsStore):
            verification_id = access.record_review_verification(
                review_id=review_id,
                **payload,
            )
        else:
            verification_id = access.record_review_verification(
                review_id,
                **payload,
            )
    typer.echo(
        _json(
            {
                "verification_id": verification_id,
                **result.to_dict(),
            }
        )
    )
    if result.status != "pass":
        raise typer.Exit(code=1)


@app.command("trace")
def trace_command(
    campaign_id: Annotated[str, typer.Option("--campaign-id")],
    subject_kind: Annotated[
        str,
        typer.Option(
            "--subject-kind",
            help="candidate_cluster, finding, or unit_disposition",
        ),
    ],
    subject_id: Annotated[str, typer.Option("--subject-id")],
    ops_db: Annotated[
        Path | None,
        typer.Option("--ops-db", exists=True, dir_okay=False),
    ] = None,
    ops_url: Annotated[str | None, typer.Option("--ops-url")] = None,
    agent_id: Annotated[str, typer.Option("--agent-id")] = "review-tracer",
    api_token: Annotated[
        str | None,
        typer.Option("--api-token", envvar="SHARUR_OPS_TOKEN"),
    ] = None,
) -> None:
    """Reconstruct a bounded candidate-to-publication provenance graph."""

    with _ops_access(
        ops_db=ops_db,
        ops_url=ops_url,
        agent_id=agent_id,
        api_token=api_token,
    ) as access:
        result = trace_review_subject(
            access,
            campaign_id=campaign_id,
            subject_kind=subject_kind,
            subject_id=subject_id,
        )
    typer.echo(_json(result))


@app.command("status")
def status_command(
    campaign_id: Annotated[str, typer.Option("--campaign-id")],
    ops_db: Annotated[
        Path | None,
        typer.Option("--ops-db", exists=True, dir_okay=False),
    ] = None,
    ops_url: Annotated[str | None, typer.Option("--ops-url")] = None,
    agent_id: Annotated[str, typer.Option("--agent-id")] = "review-observer",
    api_token: Annotated[
        str | None,
        typer.Option("--api-token", envvar="SHARUR_OPS_TOKEN"),
    ] = None,
) -> None:
    """Report exact funnel, audit, verification, and queue metrics."""

    with _ops_access(
        ops_db=ops_db,
        ops_url=ops_url,
        agent_id=agent_id,
        api_token=api_token,
    ) as access:
        result = (
            review_campaign_metrics(access, campaign_id)
            if isinstance(access, OpsStore)
            else access.review_status(campaign_id)
        )
    typer.echo(_json(result))


def main() -> None:
    app()


__all__ = ["app", "main"]


if __name__ == "__main__":
    main()
