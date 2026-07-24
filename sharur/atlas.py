"""Deterministic Atlas planning, task dispatch, and coverage verification."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING, Any

import duckdb

from sharur.dataset_seal import (
    DEFAULT_SEAL_NAME,
    DatasetSealError,
    verify_dataset_seal,
)
from sharur.ingest.input_integrity import write_json_atomic
from sharur.operators.contigs import PACKET_VERSION


if TYPE_CHECKING:
    from collections.abc import Iterable

    from sharur.ops.client import SharurOps


ATLAS_PLAN_SCHEMA = "atlas-plan/1.0"
ATLAS_UNIT_SCHEMA = "atlas-unit/1.0"
ATLAS_COVERAGE_SCHEMA = "atlas-coverage/1.0"
ATLAS_REVIEW_OUTPUT_SCHEMA = "atlas-review-output/1.0"
DEFAULT_PACKET_PROTEINS = 100
DEFAULT_CHECKPOINT_INTERVAL_CONTIGS = 25
WORK_WEIGHT_FORMULA = "n_proteins + 32 * n_contigs"


def _canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    ).encode("utf-8")


def _sha256(value: Any) -> str:
    return hashlib.sha256(_canonical_bytes(value)).hexdigest()


def _read_json_object(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"Could not read JSON object {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def _write_jsonl_atomic(
    rows: Iterable[dict[str, Any]],
    path: Path,
) -> tuple[int, str]:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None
    count = 0
    digest = hashlib.sha256()
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".partial",
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
            for row in rows:
                encoded = _canonical_bytes(row) + b"\n"
                handle.write(encoded)
                digest.update(encoded)
                count += 1
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_path, path)
        directory_fd = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    finally:
        if temporary_path is not None and temporary_path.exists():
            temporary_path.unlink()
    return count, digest.hexdigest()


def _sealed_identity(
    database: Path,
    seal_path: Path,
    *,
    verify: bool,
) -> tuple[str, str]:
    if verify:
        verification = verify_dataset_seal(seal_path, db_path=database)
        if not verification.valid:
            changed = ", ".join(verification.changed_sections) or "dataset identity"
            raise DatasetSealError(
                f"Dataset seal verification failed; changed: {changed}"
            )
    seal = _read_json_object(seal_path)
    dataset_id = seal.get("dataset_id")
    if not isinstance(dataset_id, str) or len(dataset_id) < 12:
        raise DatasetSealError("Dataset seal lacks a usable dataset_id")
    return dataset_id, str(seal.get("seal_strength", "unknown"))


def _require_contig_navigation_index(
    connection: duckdb.DuckDBPyConnection,
) -> None:
    """Require the indexed access path used by exhaustive contig traversal."""
    present = connection.execute(
        """
        SELECT EXISTS (
            SELECT 1
            FROM duckdb_indexes()
            WHERE index_name = 'idx_contigs_bin'
        )
        """
    ).fetchone()[0]
    if not present:
        raise RuntimeError(
            "Atlas planning requires idx_contigs_bin. Run `sharur migrate "
            "--db PATH`, rebuild dataset.seal.json with `sharur seal --force`, "
            "and restage the query replica."
        )


def _genome_counts(database: Path, *, threads: int) -> list[tuple[Any, ...]]:
    if threads < 1:
        raise ValueError("threads must be positive")
    connection = duckdb.connect(str(database), read_only=True)
    try:
        connection.execute("SET threads = ?", [threads])
        _require_contig_navigation_index(connection)
        return connection.execute(
            """
            WITH contig_counts AS (
                SELECT bin_id, COUNT(*) AS n_contigs
                FROM contigs
                GROUP BY bin_id
            ),
            protein_counts AS (
                SELECT bin_id, COUNT(*) AS n_proteins
                FROM proteins
                GROUP BY bin_id
            )
            SELECT
                b.bin_id,
                b.n_contigs AS declared_n_contigs,
                COALESCE(c.n_contigs, 0) AS n_contigs,
                COALESCE(p.n_proteins, 0) AS n_proteins
            FROM bins AS b
            LEFT JOIN contig_counts AS c ON c.bin_id = b.bin_id
            LEFT JOIN protein_counts AS p ON p.bin_id = b.bin_id
            ORDER BY b.bin_id
            """
        ).fetchall()
    finally:
        connection.close()


def build_atlas_plan(
    db_path: str | Path,
    output_dir: str | Path,
    *,
    seal_path: str | Path | None = None,
    packet_protein_limit: int = DEFAULT_PACKET_PROTEINS,
    checkpoint_interval_contigs: int = DEFAULT_CHECKPOINT_INTERVAL_CONTIGS,
    threads: int = 4,
    verify_seal: bool = True,
) -> dict[str, Any]:
    """Build one stable genome-owned work unit per live database bin."""
    if not 1 <= packet_protein_limit <= 250:
        raise ValueError("packet_protein_limit must be between 1 and 250")
    if checkpoint_interval_contigs < 1:
        raise ValueError("checkpoint_interval_contigs must be positive")
    database = Path(db_path).expanduser().resolve()
    destination = Path(output_dir).expanduser().resolve()
    resolved_seal = (
        Path(seal_path).expanduser().resolve()
        if seal_path is not None
        else database.parent / DEFAULT_SEAL_NAME
    )
    if not database.is_file():
        raise FileNotFoundError(f"DuckDB file does not exist: {database}")
    dataset_id, seal_strength = _sealed_identity(
        database,
        resolved_seal,
        verify=verify_seal,
    )
    counts = _genome_counts(database, threads=threads)
    plan_identity = {
        "schema_version": ATLAS_PLAN_SCHEMA,
        "dataset_id": dataset_id,
        "packet_version": PACKET_VERSION,
        "packet_protein_limit": packet_protein_limit,
        "checkpoint_interval_contigs": checkpoint_interval_contigs,
        "ownership_unit": "genome",
        "checkpoint_order": "contig_id",
        "protein_order": "gene_index NULLS LAST, start, protein_id",
    }
    plan_id = _sha256(plan_identity)
    units: list[dict[str, Any]] = []
    declared_mismatches = 0
    for bin_id, declared_n_contigs, n_contigs, n_proteins in counts:
        genome_id = str(bin_id)
        actual_contigs = int(n_contigs)
        actual_proteins = int(n_proteins)
        declared = int(declared_n_contigs) if declared_n_contigs is not None else None
        if declared is not None and declared != actual_contigs:
            declared_mismatches += 1
        unit_identity = {
            "dataset_id": dataset_id,
            "plan_id": plan_id,
            "genome_id": genome_id,
        }
        units.append(
            {
                "schema_version": ATLAS_UNIT_SCHEMA,
                "unit_id": f"atlas-{_sha256(unit_identity)[:24]}",
                "plan_id": plan_id,
                "dataset_id": dataset_id,
                "genome_id": genome_id,
                "declared_n_contigs": declared,
                "n_contigs": actual_contigs,
                "n_proteins": actual_proteins,
                "packet_version": PACKET_VERSION,
                "packet_protein_limit": packet_protein_limit,
                "checkpoint_key": "atlas_progress",
                "checkpoint_interval_contigs": checkpoint_interval_contigs,
                "work_weight": actual_proteins + 32 * actual_contigs,
            }
        )

    destination.mkdir(parents=True, exist_ok=True)
    units_path = destination / "units.jsonl"
    unit_count, units_digest = _write_jsonl_atomic(units, units_path)
    total_contigs = sum(int(unit["n_contigs"]) for unit in units)
    total_proteins = sum(int(unit["n_proteins"]) for unit in units)
    manifest = {
        **plan_identity,
        "plan_id": plan_id,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "database": str(database),
        "seal_path": str(resolved_seal),
        "seal_strength": seal_strength,
        "units_file": units_path.name,
        "units_sha256": units_digest,
        "unit_count": unit_count,
        "total_contigs": total_contigs,
        "total_proteins": total_proteins,
        "declared_contig_count_mismatches": declared_mismatches,
        "work_weight_formula": WORK_WEIGHT_FORMULA,
        "coverage_directory": "coverage",
    }
    manifest_path = destination / "plan.json"
    write_json_atomic(manifest, manifest_path)
    return {
        "plan_path": str(manifest_path),
        "units_path": str(units_path),
        **manifest,
    }


def load_atlas_plan(
    plan_dir: str | Path,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Read a plan and verify the exact unit stream digest and count."""
    root = Path(plan_dir).expanduser().resolve()
    manifest = _read_json_object(root / "plan.json")
    if manifest.get("schema_version") != ATLAS_PLAN_SCHEMA:
        raise ValueError("Atlas plan schema is unsupported")
    units_name = manifest.get("units_file")
    if not isinstance(units_name, str):
        raise ValueError("Atlas plan lacks units_file")
    units_path = root / units_name
    encoded_units = units_path.read_bytes()
    if hashlib.sha256(encoded_units).hexdigest() != manifest.get("units_sha256"):
        raise ValueError("Atlas units digest differs from plan.json")
    units: list[dict[str, Any]] = []
    for line_number, line in enumerate(encoded_units.splitlines(), start=1):
        try:
            unit = json.loads(line)
        except json.JSONDecodeError as exc:
            raise ValueError(
                f"Invalid Atlas unit JSON on line {line_number}"
            ) from exc
        if not isinstance(unit, dict):
            raise ValueError(f"Atlas unit line {line_number} is not an object")
        units.append(unit)
    if len(units) != manifest.get("unit_count"):
        raise ValueError("Atlas unit count differs from plan.json")
    return manifest, units


def enqueue_atlas_plan(
    plan_dir: str | Path,
    *,
    query_url: str,
    ops: SharurOps,
    priority: int = 1,
    max_attempts: int = 5,
    lease_seconds: int = 900,
    scan_execution_profile: str = "atlas_scan",
) -> dict[str, Any]:
    """Create an idempotent campaign and one claimable task per genome."""
    scan_execution_profile = scan_execution_profile.strip()
    if not scan_execution_profile:
        raise ValueError("scan_execution_profile must be non-empty")
    manifest, units = load_atlas_plan(plan_dir)
    root = Path(plan_dir).expanduser().resolve()
    campaign_id = ops.create_campaign(
        f"Atlas {manifest['dataset_id'][:12]}",
        description="Exhaustive genome-owned, contig-checkpointed Atlas reading",
        dataset_path=str(manifest["database"]),
        metadata={
            "plan_id": manifest["plan_id"],
            "dataset_id": manifest["dataset_id"],
            "unit_count": manifest["unit_count"],
            "query_url": query_url,
            "scan_execution_profile": scan_execution_profile,
        },
        idempotency_key=f"atlas-campaign:{manifest['plan_id']}",
    )
    task_ids: list[str] = []
    for unit in units:
        unit_id = str(unit["unit_id"])
        task_ids.append(
            ops.create_task(
                "atlas_genome_read",
                f"Read every contig in genome {unit['genome_id']}",
                params={
                    **unit,
                    "query_url": query_url,
                    "plan_path": str(root / "plan.json"),
                    "coverage_manifest": str(root / "coverage" / f"{unit_id}.json"),
                    "execution_profile": scan_execution_profile,
                    "review_output_contract": {
                        "schema_version": ATLAS_REVIEW_OUTPUT_SCHEMA,
                        "candidate_endpoint": "/review/candidates",
                        "disposition_endpoint": "/review/unit-dispositions",
                        "ordering": "candidates_before_disposition",
                        "completion_requires_active_disposition": True,
                        "candidate_reduction": "exact_typed_signature",
                    },
                },
                priority=priority,
                domain_hint="genome_reading",
                campaign_id=campaign_id,
                idempotency_key=f"atlas-task:{manifest['plan_id']}:{unit_id}",
                required_capabilities=[
                    "atlas_reader",
                    f"profile:{scan_execution_profile}",
                ],
                resource_request={"cpu_slots": 1},
                max_attempts=max_attempts,
                lease_seconds=lease_seconds,
            )
        )
    return {
        "campaign_id": campaign_id,
        "plan_id": manifest["plan_id"],
        "dataset_id": manifest["dataset_id"],
        "scan_execution_profile": scan_execution_profile,
        "task_count": len(task_ids),
        "task_ids": task_ids,
    }


def build_genome_coverage_manifest(
    unit: dict[str, Any],
    contigs: list[dict[str, Any]],
) -> dict[str, Any]:
    """Build a complete or explicitly incomplete per-genome coverage record."""
    required = {
        "unit_id",
        "plan_id",
        "dataset_id",
        "genome_id",
        "n_contigs",
        "n_proteins",
        "packet_version",
    }
    missing_fields = sorted(required - set(unit))
    if missing_fields:
        raise ValueError(f"Atlas unit lacks fields: {', '.join(missing_fields)}")
    normalized: list[dict[str, Any]] = []
    seen: set[str] = set()
    errors: list[str] = []
    for record in sorted(contigs, key=lambda value: str(value.get("contig_id", ""))):
        contig_id = record.get("contig_id")
        protein_count = record.get("protein_count")
        packet_count = record.get("packet_count")
        if not isinstance(contig_id, str) or not contig_id:
            raise ValueError("Coverage contig_id must be a non-empty string")
        if contig_id in seen:
            raise ValueError(f"Coverage repeats contig_id {contig_id!r}")
        seen.add(contig_id)
        if (
            isinstance(protein_count, bool)
            or not isinstance(protein_count, int)
            or protein_count < 0
        ):
            raise ValueError("Coverage protein_count must be a non-negative integer")
        if (
            isinstance(packet_count, bool)
            or not isinstance(packet_count, int)
            or packet_count < 0
        ):
            raise ValueError("Coverage packet_count must be a non-negative integer")
        complete = record.get("complete") is True
        if not complete:
            errors.append(f"Contig {contig_id} lacks a terminal packet")
        normalized.append(
            {
                "contig_id": contig_id,
                "protein_count": protein_count,
                "packet_count": packet_count,
                "complete": complete,
            }
        )
    observed_contigs = len(normalized)
    observed_proteins = sum(record["protein_count"] for record in normalized)
    if observed_contigs != int(unit["n_contigs"]):
        errors.append(
            f"Observed {observed_contigs} contigs; expected {int(unit['n_contigs'])}"
        )
    if observed_proteins != int(unit["n_proteins"]):
        errors.append(
            f"Observed {observed_proteins} proteins; expected {int(unit['n_proteins'])}"
        )
    coverage_identity = {
        "schema_version": ATLAS_COVERAGE_SCHEMA,
        "unit_id": unit["unit_id"],
        "plan_id": unit["plan_id"],
        "dataset_id": unit["dataset_id"],
        "genome_id": unit["genome_id"],
        "packet_version": unit["packet_version"],
        "expected_n_contigs": int(unit["n_contigs"]),
        "observed_n_contigs": observed_contigs,
        "expected_n_proteins": int(unit["n_proteins"]),
        "observed_n_proteins": observed_proteins,
        "packet_count": sum(record["packet_count"] for record in normalized),
        "contigs": normalized,
        "coverage_status": "complete" if not errors else "incomplete",
        "errors": errors,
    }
    return {
        **coverage_identity,
        "coverage_sha256": _sha256(coverage_identity),
        "generated_at": datetime.now(timezone.utc).isoformat(),
    }


def write_genome_coverage_manifest(
    unit: dict[str, Any],
    contigs: list[dict[str, Any]],
    output_path: str | Path,
) -> dict[str, Any]:
    """Atomically write one per-genome coverage manifest."""
    manifest = build_genome_coverage_manifest(unit, contigs)
    write_json_atomic(manifest, Path(output_path).expanduser().resolve())
    return manifest


def verify_atlas_coverage(
    plan_dir: str | Path,
    *,
    coverage_dir: str | Path | None = None,
) -> dict[str, Any]:
    """Verify plan-wide ownership and exact per-genome coverage manifests."""
    root = Path(plan_dir).expanduser().resolve()
    plan, units = load_atlas_plan(root)
    coverage_root = (
        Path(coverage_dir).expanduser().resolve()
        if coverage_dir is not None
        else root / str(plan["coverage_directory"])
    )
    missing: list[str] = []
    invalid: dict[str, list[str]] = {}
    observed_contigs = 0
    observed_proteins = 0
    for unit in units:
        unit_id = str(unit["unit_id"])
        path = coverage_root / f"{unit_id}.json"
        if not path.is_file():
            missing.append(unit_id)
            continue
        manifest = _read_json_object(path)
        errors: list[str] = []
        for key in ("unit_id", "plan_id", "dataset_id", "genome_id", "packet_version"):
            if manifest.get(key) != unit.get(key):
                errors.append(f"{key} differs from the assigned unit")
        identity = {
            key: value
            for key, value in manifest.items()
            if key not in {"coverage_sha256", "generated_at"}
        }
        if manifest.get("coverage_sha256") != _sha256(identity):
            errors.append("coverage_sha256 is invalid")
        if manifest.get("coverage_status") != "complete":
            errors.extend(str(value) for value in manifest.get("errors", []))
        if manifest.get("observed_n_contigs") != unit.get("n_contigs"):
            errors.append("contig count differs from the assigned unit")
        if manifest.get("observed_n_proteins") != unit.get("n_proteins"):
            errors.append("protein count differs from the assigned unit")
        if errors:
            invalid[unit_id] = sorted(set(errors))
            continue
        observed_contigs += int(manifest["observed_n_contigs"])
        observed_proteins += int(manifest["observed_n_proteins"])
    complete_units = len(units) - len(missing) - len(invalid)
    return {
        "status": "complete" if complete_units == len(units) else "incomplete",
        "plan_id": plan["plan_id"],
        "dataset_id": plan["dataset_id"],
        "expected_units": len(units),
        "complete_units": complete_units,
        "missing_units": missing,
        "invalid_units": invalid,
        "observed_contigs": observed_contigs,
        "expected_contigs": int(plan["total_contigs"]),
        "observed_proteins": observed_proteins,
        "expected_proteins": int(plan["total_proteins"]),
    }


__all__ = [
    "ATLAS_COVERAGE_SCHEMA",
    "ATLAS_PLAN_SCHEMA",
    "ATLAS_REVIEW_OUTPUT_SCHEMA",
    "ATLAS_UNIT_SCHEMA",
    "DEFAULT_CHECKPOINT_INTERVAL_CONTIGS",
    "DEFAULT_PACKET_PROTEINS",
    "WORK_WEIGHT_FORMULA",
    "build_atlas_plan",
    "build_genome_coverage_manifest",
    "enqueue_atlas_plan",
    "load_atlas_plan",
    "verify_atlas_coverage",
    "write_genome_coverage_manifest",
]
