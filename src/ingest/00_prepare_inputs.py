#!/usr/bin/env python3
"""Stage 00: validate and organize input genome assemblies."""

from __future__ import annotations

import hashlib
import json
import logging
import shutil
from pathlib import Path
from typing import Any

import typer
from rich.console import Console
from rich.table import Table

from sharur.ingest.input_integrity import (
    DEFAULT_FASTA_EXTENSIONS,
    IntegrityIssue,
    audit_fasta_file,
    audit_input_directory,
    calculate_checksums,
    write_json_atomic,
)
from sharur.ingest.input_integrity import (
    find_fasta_files as _find_fasta_files,
)
from sharur.ingest.input_integrity import (
    generate_genome_id as _generate_genome_id,
)


console = Console()
logger = logging.getLogger(__name__)


def validate_fasta_format(file_path: Path) -> dict[str, Any]:
    """Compatibility wrapper around the Stage-00 integrity auditor."""
    audit = audit_fasta_file(file_path)
    errors = [issue for issue in audit.issues if issue.severity == "error"]
    return {
        "is_valid": audit.valid,
        "error_message": errors[0].message if errors else None,
        "sequence_count": audit.sequence_count or 0,
        "total_length": audit.total_length or 0,
        "sequence_ids": audit.sequence_ids,
        "duplicate_ids": [
            issue.record_id
            for issue in errors
            if issue.code == "duplicate_record_id" and issue.record_id is not None
        ],
        "invalid_characters": [],
        "issues": [issue.to_dict() for issue in audit.issues],
    }


def calculate_file_checksum(file_path: Path) -> str:
    """Return the legacy MD5 checksum retained in processing manifests."""
    md5, _sha256 = calculate_checksums(file_path)
    return md5


def find_genome_files(input_dir: Path, extensions: list[str]) -> list[Path]:
    """Compatibility wrapper for callers of the original Stage-00 helper."""
    return _find_fasta_files(input_dir, extensions)


def generate_genome_id(file_path: Path) -> str:
    """Compatibility wrapper for callers of the original Stage-00 helper."""
    return _generate_genome_id(file_path)


def _remove_existing_output(output_dir: Path) -> None:
    if output_dir.is_dir() and not output_dir.is_symlink():
        shutil.rmtree(output_dir)
    else:
        output_dir.unlink(missing_ok=True)


def _source_set_sha256(files: list[dict[str, Any]]) -> str:
    records = [
        {
            "filename": entry["filename"],
            "genome_id": entry["genome_id"],
            "file_size": entry["file_size"],
            "sha256": entry["sha256"],
        }
        for entry in files
    ]
    encoded = json.dumps(
        records,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _print_integrity_summary(report_payload: dict[str, Any], report_path: Path) -> None:
    summary = report_payload["summary"]
    table = Table(title="Stage 00 input-integrity contract")
    table.add_column("Metric", style="cyan")
    table.add_column("Value", style="magenta")
    table.add_row("Status", report_payload["status"])
    table.add_row("Files scanned", str(summary["files_scanned"]))
    table.add_row("Records scanned", str(summary["records_scanned"]))
    table.add_row("Total bases", str(summary["total_bases"]))
    table.add_row("Errors", str(summary["errors"]))
    table.add_row("Warnings", str(summary["warnings"]))
    table.add_row("Report", str(report_path))
    console.print(table)

    for issue in report_payload["issues"]:
        location = issue["file"]
        if issue.get("line") is not None:
            location += f":{issue['line']}"
        related = (
            f" (first seen in {issue['related_file']})"
            if issue.get("related_file")
            else ""
        )
        style = "red" if issue["severity"] == "error" else "yellow"
        console.print(
            f"[{style}]{issue['severity'].upper()} "
            f"{issue['code']} at {location}: {issue['message']}{related}[/{style}]"
        )


def prepare_inputs(
    input_dir: Path = typer.Option(
        Path("data/raw"),
        "--input-dir",
        "-i",
        help="Directory containing input genome assemblies",
    ),
    output_dir: Path = typer.Option(
        Path("data/stage00_prepared"),
        "--output-dir",
        "-o",
        help="Output directory for validated assemblies",
    ),
    file_extensions: list[str] = typer.Option(
        list(DEFAULT_FASTA_EXTENSIONS),
        "--extensions",
        "-e",
        help="File extensions to search for",
    ),
    validate_format: bool = typer.Option(
        True,
        "--validate/--no-validate",
        help=(
            "Validate FASTA structure, residues, empty records, and global IDs. "
            "Disabling this is only for explicit legacy recovery."
        ),
    ),
    copy_files: bool = typer.Option(
        False,
        "--copy/--symlink",
        help="Copy files instead of creating symlinks",
    ),
    force: bool = typer.Option(
        False,
        "--force",
        help="Overwrite an existing Stage-00 output directory",
    ),
) -> None:
    """Audit the complete input set before exposing it to downstream stages."""
    console.print("[bold blue]Stage 00: Input Preparation[/bold blue]")
    resolved_input = input_dir.expanduser().resolve()
    resolved_output = output_dir.expanduser().resolve()

    if not resolved_input.exists():
        console.print(f"[red]Input directory does not exist: {resolved_input}[/red]")
        raise typer.Exit(1)
    if not resolved_input.is_dir():
        console.print(f"[red]Input path is not a directory: {resolved_input}[/red]")
        raise typer.Exit(1)
    if resolved_output == resolved_input or resolved_input.is_relative_to(resolved_output):
        console.print(
            "[red]Output directory cannot equal or contain the input directory.[/red]"
        )
        raise typer.Exit(1)

    if resolved_output.exists():
        if not force:
            console.print(f"[red]Output already exists: {resolved_output}[/red]")
            console.print("[yellow]Use --force to replace it.[/yellow]")
            raise typer.Exit(1)
        _remove_existing_output(resolved_output)
    resolved_output.mkdir(parents=True)

    report = audit_input_directory(
        resolved_input,
        extensions=file_extensions,
        validate_sequences=validate_format,
    )
    report_path = resolved_output / "input_integrity.json"
    report_payload = report.to_dict()
    write_json_atomic(report_payload, report_path)
    _print_integrity_summary(report_payload, report_path)

    if report.status == "failed":
        console.print(
            "[bold red]Stage 00 rejected the input set. "
            "No downstream processing manifest or assembly links were created.[/bold red]"
        )
        raise typer.Exit(1)
    if report.status == "skipped":
        console.print(
            "[yellow]FASTA content validation was explicitly disabled; "
            "the integrity contract is marked skipped.[/yellow]"
        )

    genome_entries: list[dict[str, Any]] = []
    for audit in report.files:
        output_file = resolved_output / audit.path.name
        entry = {
            "filename": audit.path.name,
            "genome_id": audit.genome_id,
            "file_path": str(audit.path),
            "file_size": audit.file_size,
            "checksum": audit.md5,
            "checksum_algorithm": "md5",
            "sha256": audit.sha256,
            "format_valid": audit.valid if validate_format else True,
            "validation_errors": [
                issue.message
                for issue in audit.issues
                if issue.severity == "error"
            ],
            "sequence_count": audit.sequence_count,
            "total_length": audit.total_length,
            "sequence_ids": audit.sequence_ids,
        }
        try:
            if copy_files:
                shutil.copy2(audit.path, output_file)
            else:
                output_file.symlink_to(audit.path)
            entry["output_path"] = str(output_file)
            entry["linked_successfully"] = True
        except OSError as exc:
            entry["linked_successfully"] = False
            entry["link_error"] = f"{type(exc).__name__}: {exc}"
            report.global_issues.append(
                IntegrityIssue(
                    severity="error",
                    code="assembly_publish_failed",
                    file=audit.path.name,
                    message=entry["link_error"],
                )
            )
        genome_entries.append(entry)

    if report.status == "failed":
        for entry in genome_entries:
            output_path = entry.get("output_path")
            if isinstance(output_path, str):
                Path(output_path).unlink(missing_ok=True)
        write_json_atomic(report.to_dict(), report_path)
        console.print(
            "[bold red]Stage 00 could not publish every validated assembly. "
            "No processing manifest was created.[/bold red]"
        )
        raise typer.Exit(1)

    manifest_files = [
        {
            "filename": entry["filename"],
            "genome_id": entry["genome_id"],
            "file_size": entry["file_size"],
            "sha256": entry["sha256"],
        }
        for entry in genome_entries
    ]
    manifest = {
        "version": "0.2.0",
        "stage": "stage00_prepare_inputs",
        "timestamp": report.generated_at,
        "input_dir": str(resolved_input),
        "output_dir": str(resolved_output),
        "file_extensions": list(report.extensions),
        "validate_format": validate_format,
        "copy_files": copy_files,
        "integrity_report": str(report_path),
        "integrity_status": report.status,
        "source_set_sha256": _source_set_sha256(manifest_files),
        "summary": report_payload["summary"],
        "genomes": genome_entries,
    }
    manifest_path = resolved_output / "processing_manifest.json"
    write_json_atomic(manifest, manifest_path)

    console.print(
        f"[bold green]Stage 00 accepted {len(genome_entries)} assemblies.[/bold green]"
    )
    console.print(f"Processing manifest: {manifest_path}")
    logger.info(
        "Stage 00 completed: %s assemblies, %s warnings",
        len(genome_entries),
        report.warning_count,
    )


def main() -> None:
    """Entry point for standalone and packaged execution."""
    logging.basicConfig(level=logging.INFO)
    typer.run(prepare_inputs)


if __name__ == "__main__":
    main()
