"""Deterministic, fail-early FASTA integrity checks for Stage 00."""

from __future__ import annotations

import hashlib
import json
import os
import re
import tempfile
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Literal


INTEGRITY_SCHEMA_VERSION = 1
DEFAULT_FASTA_EXTENSIONS = (".fasta", ".fa", ".fna")
VALID_NUCLEOTIDE_CHARACTERS = frozenset("ACGTUNRYSWKMBDHVX-")


@dataclass(frozen=True)
class IntegrityIssue:
    """One actionable problem or warning found during input inspection."""

    severity: Literal["error", "warning"]
    code: str
    file: str
    message: str
    line: int | None = None
    record_id: str | None = None
    related_file: str | None = None

    def to_dict(self) -> dict:
        return {
            key: value
            for key, value in asdict(self).items()
            if value is not None
        }


@dataclass
class FastaFileAudit:
    """Integrity result for one assembly FASTA."""

    path: Path
    genome_id: str
    file_size: int = 0
    md5: str | None = None
    sha256: str | None = None
    sequence_count: int | None = None
    total_length: int | None = None
    sequence_ids: list[str] = field(default_factory=list)
    issues: list[IntegrityIssue] = field(default_factory=list)

    @property
    def valid(self) -> bool:
        return not any(issue.severity == "error" for issue in self.issues)

    def to_dict(self, *, input_dir: Path) -> dict:
        try:
            relative_path = str(self.path.relative_to(input_dir))
        except ValueError:
            relative_path = self.path.name
        return {
            "filename": self.path.name,
            "relative_path": relative_path,
            "genome_id": self.genome_id,
            "file_size": self.file_size,
            "checksums": {
                "md5": self.md5,
                "sha256": self.sha256,
            },
            "sequence_count": self.sequence_count,
            "total_length": self.total_length,
            "sequence_ids": self.sequence_ids,
            "valid": self.valid,
            "issue_count": len(self.issues),
        }


@dataclass
class InputIntegrityReport:
    """Complete Stage-00 input contract for one directory."""

    input_dir: Path
    extensions: tuple[str, ...]
    validation_enabled: bool
    files: list[FastaFileAudit]
    global_issues: list[IntegrityIssue] = field(default_factory=list)
    generated_at: str = field(
        default_factory=lambda: datetime.now(timezone.utc).isoformat()
    )

    @property
    def issues(self) -> list[IntegrityIssue]:
        return [
            *(issue for audit in self.files for issue in audit.issues),
            *self.global_issues,
        ]

    @property
    def error_count(self) -> int:
        return sum(issue.severity == "error" for issue in self.issues)

    @property
    def warning_count(self) -> int:
        return sum(issue.severity == "warning" for issue in self.issues)

    @property
    def status(self) -> Literal["passed", "failed", "skipped"]:
        if self.error_count:
            return "failed"
        if not self.validation_enabled:
            return "skipped"
        return "passed"

    def to_dict(self) -> dict:
        return {
            "schema_version": INTEGRITY_SCHEMA_VERSION,
            "generated_at": self.generated_at,
            "status": self.status,
            "validation_enabled": self.validation_enabled,
            "input_dir": str(self.input_dir),
            "extensions": list(self.extensions),
            "summary": {
                "files_scanned": len(self.files),
                "records_scanned": sum(
                    audit.sequence_count or 0 for audit in self.files
                ),
                "total_bases": sum(
                    audit.total_length or 0 for audit in self.files
                ),
                "errors": self.error_count,
                "warnings": self.warning_count,
            },
            "files": [
                audit.to_dict(input_dir=self.input_dir) for audit in self.files
            ],
            "issues": [issue.to_dict() for issue in self.issues],
        }


def normalize_extensions(extensions: list[str] | tuple[str, ...]) -> tuple[str, ...]:
    """Normalize extensions to unique, lowercase, dot-prefixed values."""
    normalized = {
        extension.lower()
        if extension.startswith(".")
        else f".{extension.lower()}"
        for extension in extensions
    }
    return tuple(sorted(normalized))


def find_fasta_files(
    input_dir: Path,
    extensions: list[str] | tuple[str, ...] = DEFAULT_FASTA_EXTENSIONS,
) -> list[Path]:
    """Return unique top-level FASTA inputs in deterministic order."""
    normalized = normalize_extensions(extensions)
    return sorted(
        {
            path
            for path in input_dir.iterdir()
            if path.is_file() and path.suffix.lower() in normalized
        },
        key=lambda path: path.name,
    )


def generate_genome_id(file_path: Path) -> str:
    """Generate the stable Stage-00 genome ID from an input filename."""
    name = file_path.stem
    genome_id = re.sub(r"[^a-zA-Z0-9_]", "_", name)
    genome_id = re.sub(r"_+", "_", genome_id)
    return genome_id.strip("_")


def calculate_checksums(file_path: Path) -> tuple[str, str]:
    """Calculate MD5 for legacy manifests and SHA-256 for the integrity contract."""
    md5 = hashlib.md5(usedforsecurity=False)
    sha256 = hashlib.sha256()
    with file_path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            md5.update(chunk)
            sha256.update(chunk)
    return md5.hexdigest(), sha256.hexdigest()


def audit_fasta_file(
    file_path: Path,
    *,
    validate_sequences: bool = True,
) -> FastaFileAudit:
    """Audit one nucleotide FASTA, including terminal empty records."""
    audit = FastaFileAudit(path=file_path, genome_id=generate_genome_id(file_path))
    try:
        audit.file_size = file_path.stat().st_size
        audit.md5, audit.sha256 = calculate_checksums(file_path)
    except OSError as exc:
        audit.issues.append(
            IntegrityIssue(
                severity="error",
                code="file_unreadable",
                file=file_path.name,
                message=f"{type(exc).__name__}: {exc}",
            )
        )
        return audit

    if not validate_sequences:
        return audit

    audit.sequence_count = 0
    audit.total_length = 0
    current_record_id: str | None = None
    current_header_line: int | None = None
    current_length = 0
    record_open = False
    seen_ids: set[str] = set()

    def finish_record() -> None:
        nonlocal current_length, record_open
        if not record_open:
            return
        if current_length == 0:
            audit.issues.append(
                IntegrityIssue(
                    severity="error",
                    code="empty_sequence",
                    file=file_path.name,
                    line=current_header_line,
                    record_id=current_record_id,
                    message="FASTA header has no sequence residues.",
                )
            )
        audit.total_length = (audit.total_length or 0) + current_length
        current_length = 0
        record_open = False

    try:
        with file_path.open("r", encoding="utf-8", errors="strict") as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                line = raw_line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    finish_record()
                    audit.sequence_count = (audit.sequence_count or 0) + 1
                    current_header_line = line_number
                    record_open = True
                    header = line[1:].strip()
                    current_record_id = header.split()[0] if header else None
                    if current_record_id is None:
                        audit.issues.append(
                            IntegrityIssue(
                                severity="error",
                                code="empty_header",
                                file=file_path.name,
                                line=line_number,
                                message="FASTA header does not contain a record identifier.",
                            )
                        )
                        continue
                    audit.sequence_ids.append(current_record_id)
                    if current_record_id in seen_ids:
                        audit.issues.append(
                            IntegrityIssue(
                                severity="error",
                                code="duplicate_record_id",
                                file=file_path.name,
                                line=line_number,
                                record_id=current_record_id,
                                message="Record identifier is duplicated within this file.",
                            )
                        )
                    seen_ids.add(current_record_id)
                    continue

                if not record_open:
                    audit.issues.append(
                        IntegrityIssue(
                            severity="error",
                            code="sequence_before_header",
                            file=file_path.name,
                            line=line_number,
                            message="Sequence residues occur before the first FASTA header.",
                        )
                    )
                    continue

                invalid = sorted(
                    set(line.upper()) - VALID_NUCLEOTIDE_CHARACTERS
                )
                if invalid:
                    audit.issues.append(
                        IntegrityIssue(
                            severity="error",
                            code="invalid_sequence_characters",
                            file=file_path.name,
                            line=line_number,
                            record_id=current_record_id,
                            message=(
                                "Sequence contains invalid nucleotide characters: "
                                + ", ".join(invalid)
                            ),
                        )
                    )
                current_length += len(line)
        finish_record()
    except (OSError, UnicodeError) as exc:
        audit.issues.append(
            IntegrityIssue(
                severity="error",
                code="file_read_error",
                file=file_path.name,
                message=f"{type(exc).__name__}: {exc}",
            )
        )

    if audit.sequence_count == 0:
        audit.issues.append(
            IntegrityIssue(
                severity="error",
                code="no_fasta_records",
                file=file_path.name,
                message="File contains no FASTA records.",
            )
        )
    return audit


def audit_input_directory(
    input_dir: Path,
    *,
    extensions: list[str] | tuple[str, ...] = DEFAULT_FASTA_EXTENSIONS,
    validate_sequences: bool = True,
) -> InputIntegrityReport:
    """Audit all input FASTAs and enforce cross-file identifier uniqueness."""
    resolved = input_dir.expanduser().resolve()
    normalized_extensions = normalize_extensions(extensions)
    files = find_fasta_files(resolved, normalized_extensions)
    audits = [
        audit_fasta_file(path, validate_sequences=validate_sequences)
        for path in files
    ]
    report = InputIntegrityReport(
        input_dir=resolved,
        extensions=normalized_extensions,
        validation_enabled=validate_sequences,
        files=audits,
    )

    if not audits:
        report.global_issues.append(
            IntegrityIssue(
                severity="error",
                code="no_input_files",
                file=str(resolved),
                message=(
                    "No top-level FASTA files matched: "
                    + ", ".join(normalized_extensions)
                ),
            )
        )
        return report

    genome_ids: dict[str, Path] = {}
    record_ids: dict[str, Path] = {}
    content_hashes: dict[str, Path] = {}
    for audit in audits:
        if not audit.genome_id:
            report.global_issues.append(
                IntegrityIssue(
                    severity="error",
                    code="empty_genome_id",
                    file=audit.path.name,
                    message="Filename does not produce a usable genome identifier.",
                )
            )
        elif audit.genome_id in genome_ids:
            report.global_issues.append(
                IntegrityIssue(
                    severity="error",
                    code="duplicate_genome_id",
                    file=audit.path.name,
                    related_file=genome_ids[audit.genome_id].name,
                    message=(
                        f"Filename normalizes to duplicate genome ID "
                        f"{audit.genome_id!r}."
                    ),
                )
            )
        else:
            genome_ids[audit.genome_id] = audit.path

        if validate_sequences:
            for record_id in set(audit.sequence_ids):
                first_path = record_ids.get(record_id)
                if first_path is not None and first_path != audit.path:
                    report.global_issues.append(
                        IntegrityIssue(
                            severity="error",
                            code="duplicate_record_id_across_files",
                            file=audit.path.name,
                            related_file=first_path.name,
                            record_id=record_id,
                            message=(
                                "Record identifiers must be globally unique across "
                                "input assemblies."
                            ),
                        )
                    )
                else:
                    record_ids[record_id] = audit.path

        if audit.sha256 is not None:
            first_path = content_hashes.get(audit.sha256)
            if first_path is not None:
                report.global_issues.append(
                    IntegrityIssue(
                        severity="warning",
                        code="duplicate_file_content",
                        file=audit.path.name,
                        related_file=first_path.name,
                        message="Input files have identical SHA-256 content.",
                    )
                )
            else:
                content_hashes[audit.sha256] = audit.path

    return report


def write_json_atomic(payload: dict, path: Path) -> None:
    """Write JSON durably without exposing a partial target file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w",
        encoding="utf-8",
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temporary = Path(handle.name)
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    try:
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


__all__ = [
    "DEFAULT_FASTA_EXTENSIONS",
    "INTEGRITY_SCHEMA_VERSION",
    "VALID_NUCLEOTIDE_CHARACTERS",
    "FastaFileAudit",
    "InputIntegrityReport",
    "IntegrityIssue",
    "audit_fasta_file",
    "audit_input_directory",
    "calculate_checksums",
    "find_fasta_files",
    "generate_genome_id",
    "normalize_extensions",
    "write_json_atomic",
]
