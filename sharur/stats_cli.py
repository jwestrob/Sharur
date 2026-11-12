"""``sharur stats``: pre-registration records and static checks of statistical claims."""

from __future__ import annotations

import json
from pathlib import Path

import typer


stats_app = typer.Typer(help="Pre-register predictions and check statistical claims in findings.",
                        no_args_is_help=True)


def _load_spec(path: Path) -> dict:
    text = path.read_text()
    if path.suffix.lower() in {".yaml", ".yml"}:
        import yaml  # noqa: PLC0415

        return yaml.safe_load(text)
    return json.loads(text)


@stats_app.command("preregister")
def preregister_command(
    spec: Path = typer.Argument(..., help="JSON or YAML: id, title, predictions[], correction, design_inputs."),
    out_dir: Path = typer.Option(..., "--out-dir", "-o", help="Directory for the read-only record."),
    db: Path | None = typer.Option(
        None, "--db", "-d", help="Dataset DuckDB; its dataset.seal.json fills `dataset` (read only)."),
    document: Path | None = typer.Option(
        None, "--document", help="A longer pre-registration text to bind by SHA-256."),
):
    """Freeze predictions, tests and verdict rules before the analysis runs.

    Writes <id>_<UTC time>.prereg.json (read-only, never replaced) and prints its content_sha256,
    which every statistical check implementing it cites.
    """
    from sharur.preregistration import (  # noqa: PLC0415
        PreregError,
        Preregistration,
        dataset_ref_from_db,
        document_ref,
        write_preregistration,
    )

    data = _load_spec(spec)
    try:
        if db is not None:
            data["dataset"] = dataset_ref_from_db(db).model_dump(exclude_none=True)
        if document is not None:
            data["document"] = document_ref(document).model_dump(exclude_none=True)
        record = Preregistration.model_validate(data)
        path, written = write_preregistration(record, out_dir)
    except (PreregError, ValueError, KeyError, OSError) as exc:
        typer.echo(f"preregister: {exc}", err=True)
        raise typer.Exit(1) from exc
    typer.echo(json.dumps({"path": str(path), "id": written.id, "written_at": written.written_at,
                           "content_sha256": written.content_sha256,
                           "predictions": [p.id for p in written.predictions]}, indent=1))


@stats_app.command("verify-prereg")
def verify_prereg_command(
    record: Path = typer.Argument(..., help="A .prereg.json record."),
    seal: bool = typer.Option(False, "--seal", help="Also check the dataset against its seal (slow, read-only)."),
    db: Path | None = typer.Option(None, "--db", "-d", help="Current database path for --seal."),
    output_format: str = typer.Option("markdown", "--format", "-f", help="markdown or json"),
):
    """Recompute a record's hash and check its document and, with --seal, its dataset."""
    from sharur.preregistration import PreregError, verify_preregistration  # noqa: PLC0415

    try:
        result = verify_preregistration(record, check_seal=seal, db_path=db)
    except PreregError as exc:
        typer.echo(f"verify-prereg: {exc}", err=True)
        raise typer.Exit(1) from exc
    typer.echo(result.model_dump_json(indent=1) if output_format == "json" else result.to_markdown())
    if not result.valid:
        raise typer.Exit(1)


@stats_app.command("check")
def check_command(
    findings: list[Path] = typer.Argument(..., help="findings.jsonl files or dataset directories."),
    root: Path | None = typer.Option(None, "--root", help="Base for relative paths (default: repo root)."),
    output_format: str = typer.Option("markdown", "--format", "-f", help="markdown or json"),
    strict: bool = typer.Option(False, "--strict", help="Exit 1 when any statistical check has issues."),
):
    """Check every statistical claim as recorded: verdict rules, interval, controls, pre-registration,
    dataset seal and input hashes. Replaying the numbers is `sharur reverify`.
    """
    from sharur.reverify import find_root, resolve_findings  # noqa: PLC0415
    from sharur.stat_checks import check_findings  # noqa: PLC0415

    files = resolve_findings(findings)
    base = (root or find_root(files[0])).resolve()
    rows = check_findings(files, root=base)
    if output_format == "json":
        typer.echo(json.dumps(rows, indent=1))
    else:
        bad = [r for r in rows if r["issues"]]
        typer.echo(f"# Statistical checks: {len(rows)} found, {len(rows) - len(bad)} consistent, {len(bad)} with issues")
        for r in bad:
            typer.echo(f"\n## {r['finding_id']} #{r['index']} ({r['source']}:{r['line']}) {r['claim'][:120]}")
            for issue in r["issues"]:
                typer.echo(f"- {issue}")
    if strict and any(r["issues"] for r in rows):
        raise typer.Exit(1)
