"""Pre-registration records: predictions, tests and verdict rules frozen before an analysis.

A record states, before any test statistic exists, what will be tested and how the verdict
will be decided: per prediction its text, the test, the predicted direction and ordered
verdict rules (expressions over named values, see :mod:`sharur.stat_checks`), plus the
correction family, the dataset seal it refers to and, optionally, the SHA-256 of a longer
pre-registration document. ``content_sha256`` covers the whole canonical record.

:func:`write_preregistration` stamps ``written_at`` (UTC), computes the hash, writes
``<id>_<UTC timestamp>.prereg.json`` atomically, refuses to replace an existing file and
makes the file read-only. A later statistical check cites ``{id, content_sha256}``;
``sharur stats check`` and ``sharur reverify`` then confirm that the record is unchanged, was
written before the check was computed, names the same dataset and holds the same rules.

    sharur stats preregister spec.json --db data/DATASET/sharur.duckdb --document study_plan.txt \\
        --out-dir analyses/study-x/prereg
    sharur stats verify-prereg analyses/study-x/prereg/study-x_20260115T100000Z.prereg.json
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import stat
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator


SCHEMA = "sharur.preregistration/1"
SUFFIX = ".prereg.json"
_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")


class PreregError(ValueError):
    """A pre-registration record is malformed, altered or cannot be written."""


class _Model(BaseModel):
    model_config = ConfigDict(extra="forbid")


class VerdictRule(_Model):
    """``label`` applies when ``when`` (an expression over named values) is true; first match wins.

    ``label`` may name values in braces, e.g. ``"{placebos_at_or_beyond} of 10 placebos"``.
    """

    label: str = Field(min_length=1)
    when: str = Field(min_length=1)


class Prediction(_Model):
    id: str = Field(min_length=1)
    text: str = Field(min_length=1)
    test: str = Field(min_length=1, description="The statistical test, as it will be run.")
    statistic: str | None = Field(default=None, description="Name of the effect or statistic predicted.")
    direction: Literal["greater", "less", "two-sided"] | None = None
    inferential: bool = True
    verdict_rules: list[VerdictRule] = Field(min_length=1)
    extra_rules: dict[str, list[VerdictRule]] = Field(
        default_factory=dict, description="Further named rule sets (e.g. a calibration verdict).")


class Correction(_Model):
    method: Literal["BH", "Holm", "Bonferroni", "none"] = "BH"
    family: list[str] = Field(default_factory=list, description="Prediction IDs in the correction family.")
    alpha: float = Field(default=0.05, gt=0, lt=1)


class DatasetRef(_Model):
    dataset_id: str = Field(min_length=16)
    seal_path: str | None = None
    seal_strength: str | None = None
    seal_generated_at: str | None = None
    database: str | None = None


class DocumentRef(_Model):
    path: str
    sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    bytes: int | None = None


class Preregistration(_Model):
    schema_: str = Field(default=SCHEMA, alias="schema")
    id: str
    title: str = Field(min_length=1)
    written_at: str | None = None
    author: str | None = None
    dataset: DatasetRef
    document: DocumentRef | None = None
    design_inputs: list[str] = Field(default_factory=list,
                                     description="What was already in hand when this was written.")
    predictions: list[Prediction] = Field(min_length=1)
    correction: Correction = Field(default_factory=Correction)
    notes: str | None = None
    content_sha256: str | None = None

    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    @field_validator("id")
    @classmethod
    def _valid_id(cls, value: str) -> str:
        if not _ID.match(value):
            raise ValueError("id must be 1-128 characters of letters, digits, '.', '_' or '-'")
        return value

    @field_validator("predictions")
    @classmethod
    def _unique_predictions(cls, value: list[Prediction]) -> list[Prediction]:
        ids = [p.id for p in value]
        if len(ids) != len(set(ids)):
            raise ValueError("prediction ids must be unique")
        return value

    def prediction(self, prediction_id: str) -> Prediction:
        for p in self.predictions:
            if p.id == prediction_id:
                return p
        raise KeyError(prediction_id)

    def payload(self) -> dict[str, Any]:
        """The record as JSON, keys as written to disk."""
        return self.model_dump(mode="json", by_alias=True, exclude_none=True)

    def compute_hash(self) -> str:
        body = {k: v for k, v in self.payload().items() if k != "content_sha256"}
        return hashlib.sha256(canonical_bytes(body)).hexdigest()


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, default=str).encode()


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def parse_time(value: str) -> datetime:
    """ISO 8601 (``Z`` allowed) as an aware datetime; naive values are taken as UTC."""
    dt = datetime.fromisoformat(value.replace("Z", "+00:00"))
    return dt if dt.tzinfo else dt.replace(tzinfo=timezone.utc)


def dataset_ref_from_db(db_path: str | Path) -> DatasetRef:
    """The dataset seal beside a database (``dataset.seal.json``), as a reference. Reads only."""
    from sharur.dataset_seal import DEFAULT_SEAL_NAME  # noqa: PLC0415

    db = Path(db_path).expanduser().resolve()
    seal_path = db.parent / DEFAULT_SEAL_NAME
    if not seal_path.is_file():
        raise PreregError(f"No dataset seal at {seal_path}; run `sharur seal --db {db}` first")
    seal = json.loads(seal_path.read_text())
    return DatasetRef(dataset_id=str(seal["dataset_id"]), seal_path=str(seal_path),
                      seal_strength=seal.get("seal_strength"), seal_generated_at=seal.get("generated_at"),
                      database=str(db))


def document_ref(path: str | Path) -> DocumentRef:
    p = Path(path).expanduser().resolve()
    return DocumentRef(path=str(p), sha256=sha256_file(p), bytes=p.stat().st_size)


def finalize(record: Preregistration, *, written_at: str | None = None) -> Preregistration:
    """A copy with ``written_at`` (now, unless given) and ``content_sha256`` filled in."""
    data = record.model_copy(update={"written_at": written_at or record.written_at or utc_now(),
                                     "content_sha256": None})
    return data.model_copy(update={"content_sha256": data.compute_hash()})


def write_preregistration(record: Preregistration | dict[str, Any], out_dir: str | Path, *,
                          written_at: str | None = None) -> tuple[Path, Preregistration]:
    """Stamp, hash and write a record as a new read-only file; never replaces one."""
    rec = record if isinstance(record, Preregistration) else Preregistration.model_validate(record)
    if rec.content_sha256:
        raise PreregError("record already carries content_sha256; a written record is not rewritten")
    rec = finalize(rec, written_at=written_at)
    stamp = parse_time(rec.written_at).astimezone(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    directory = Path(out_dir).expanduser().resolve()
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{rec.id}_{stamp}{SUFFIX}"
    if path.exists():
        raise PreregError(f"{path} exists; pre-registrations are never replaced")
    tmp = directory / f".{path.name}.{os.getpid()}.tmp"
    tmp.write_text(json.dumps(rec.payload(), indent=1, sort_keys=False) + "\n")
    os.chmod(tmp, stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH)
    os.link(tmp, path)  # fails if the name appeared meanwhile
    tmp.unlink()
    return path, rec


def load_preregistration(path: str | Path) -> Preregistration:
    try:
        return Preregistration.model_validate(json.loads(Path(path).read_text()))
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        raise PreregError(f"{path}: {exc}") from exc


class PreregVerification(_Model):
    path: str | None = None
    id: str
    valid: bool
    content_sha256: str | None
    recomputed_sha256: str
    issues: list[str] = Field(default_factory=list)

    def to_markdown(self) -> str:
        lines = [f"# Pre-registration {self.id}: {'verified' if self.valid else 'PROBLEMS'}",
                 "", f"- content_sha256 recorded `{self.content_sha256}`, recomputed `{self.recomputed_sha256}`"]
        lines += [f"- {issue}" for issue in self.issues]
        return "\n".join(lines)


def verify_preregistration(source: str | Path | Preregistration, *, check_document: bool = True,
                           check_seal: bool = False, db_path: str | Path | None = None) -> PreregVerification:
    """Recompute the record hash; check the document hash and, on request, the dataset seal.

    ``check_seal`` compares the record's ``dataset_id`` with the seal file beside the database
    and runs ``verify-seal`` on it (slow for large datasets; read-only).
    """
    path = None
    if isinstance(source, Preregistration):
        rec = source
    else:
        path = str(Path(source).expanduser().resolve())
        rec = load_preregistration(path)
    issues: list[str] = []
    recomputed = rec.compute_hash()
    if not rec.content_sha256:
        issues.append("record has no content_sha256 (never written with write_preregistration)")
    elif rec.content_sha256 != recomputed:
        issues.append("content_sha256 does not match the record: it was edited after it was written")
    if not rec.written_at:
        issues.append("record has no written_at")
    if check_document and rec.document is not None:
        doc = Path(rec.document.path)
        if not doc.is_file():
            issues.append(f"document not found: {doc}")
        elif sha256_file(doc) != rec.document.sha256:
            issues.append(f"document {doc} changed since it was registered")
    if check_seal:
        from sharur.dataset_seal import verify_dataset_seal  # noqa: PLC0415

        seal_path = Path(rec.dataset.seal_path) if rec.dataset.seal_path else None
        if db_path is not None:
            seal_path = Path(db_path).expanduser().resolve().parent / "dataset.seal.json"
        if seal_path is None or not seal_path.is_file():
            issues.append("dataset seal file not found")
        else:
            current = json.loads(seal_path.read_text()).get("dataset_id")
            if current != rec.dataset.dataset_id:
                issues.append(f"dataset seal now names {current}, the record names {rec.dataset.dataset_id}")
            result = verify_dataset_seal(seal_path, db_path=db_path)
            if not result.valid:
                issues.append("dataset no longer matches its seal: " + ", ".join(result.changed_sections))
    return PreregVerification(path=path, id=rec.id, valid=not issues, content_sha256=rec.content_sha256,
                              recomputed_sha256=recomputed, issues=issues)


__all__ = [
    "SCHEMA",
    "Correction",
    "DatasetRef",
    "DocumentRef",
    "Prediction",
    "PreregError",
    "PreregVerification",
    "Preregistration",
    "VerdictRule",
    "dataset_ref_from_db",
    "document_ref",
    "finalize",
    "load_preregistration",
    "parse_time",
    "sha256_file",
    "utc_now",
    "verify_preregistration",
    "write_preregistration",
]
