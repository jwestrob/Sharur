"""Deterministic Atlas planning, task dispatch, and coverage verification."""

from __future__ import annotations

import hashlib
import json
import math
import os
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
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
from sharur.operators.contigs import (
    GENOME_PACKET_VERSION,
    genome_packet_packing_contract,
    get_genome_packet,
    minimum_genome_packet_record_bytes,
)
from sharur.query.runtime import ReadOnlyDuckDBRuntime
from sharur.storage.duckdb_store import DuckDBStore


if TYPE_CHECKING:
    from collections.abc import Callable, Iterable

    from sharur.ops.client import SharurOps


ATLAS_PLAN_SCHEMA = "atlas-plan/2.1"
ATLAS_UNIT_SCHEMA = "atlas-unit/2.0"
ATLAS_COVERAGE_SCHEMA = "atlas-coverage/2.0"
ATLAS_REVIEW_OUTPUT_SCHEMA = "atlas-review-output/1.0"
ATLAS_PACKET_CENSUS_SCHEMA = "atlas-packet-census/1.0"
ATLAS_PACKET_CENSUS_UNIT_SCHEMA = "atlas-packet-census-unit/1.0"
ATLAS_PACKET_CALIBRATION_SCHEMA = "atlas-packet-calibration/1.0"
DEFAULT_CHECKPOINT_INTERVAL_FRAMES = 1
DEFAULT_QUERY_RESULT_LIMIT_BYTES = 2 * 1024 * 1024
QUERY_RESULT_ENVELOPE_RESERVE_BYTES = 4 * 1024
DEFAULT_CENSUS_WORKERS = 4
DEFAULT_CENSUS_MEMORY_LIMIT = "16GB"
DEFAULT_CENSUS_MAX_TEMP_SIZE = "128GB"
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


def _select_packet_calibration_rows(
    counts: list[tuple[Any, ...]],
    *,
    sample_genomes: int,
) -> list[tuple[Any, ...]]:
    """Select deterministic genome-size quantiles for payload calibration."""
    if sample_genomes < 1:
        raise ValueError("packet_calibration_genomes must be positive")
    eligible = [
        row
        for row in counts
        if int(row[2] or 0) > 0 or int(row[3] or 0) > 0
    ]
    ordered = sorted(
        eligible,
        key=lambda row: (int(row[3] or 0), int(row[2] or 0), str(row[0])),
    )
    if len(ordered) <= sample_genomes:
        return ordered
    if sample_genomes == 1:
        return [ordered[len(ordered) // 2]]
    indices = [
        round(position * (len(ordered) - 1) / (sample_genomes - 1))
        for position in range(sample_genomes)
    ]
    return [ordered[index] for index in dict.fromkeys(indices)]


def _calibrate_packet_limits(
    database: Path,
    *,
    dataset_id: str,
    counts: list[tuple[Any, ...]],
    packet_model_payload_bytes: int,
    packet_all_annotations: bool,
    sample_genomes: int,
    sample_genome_source: str,
    threads: int,
) -> dict[str, Any]:
    """Measure real serialized frames for a byte-defined packing contract."""
    selected = _select_packet_calibration_rows(
        counts,
        sample_genomes=sample_genomes,
    )
    protein_frame_densities: list[int] = []
    contig_frame_densities: list[int] = []
    sampled_payload_bytes = 0
    protein_bearing_payload_bytes = 0
    sampled_proteins = 0
    sampled_contig_segments = 0
    target_exceeded_frames = 0
    sampled_ids: list[str] = []
    store = DuckDBStore(database, read_only=True, threads=threads)
    store.dataset_id = dataset_id
    try:
        for row in selected:
            genome_id = str(row[0])
            sampled_ids.append(genome_id)
            packet = get_genome_packet(
                store,
                genome_id,
                max_contigs=None,
                max_proteins=None,
                max_model_payload_bytes=packet_model_payload_bytes,
                all_annotations=packet_all_annotations,
            ).raw
            model_payload = packet.get("model_payload")
            if not isinstance(model_payload, dict):
                raise RuntimeError(
                    f"Packet calibration received an invalid payload in {genome_id}"
                )
            contigs = model_payload.get("contigs")
            if not isinstance(contigs, list):
                raise RuntimeError(
                    f"Packet calibration received invalid contigs in {genome_id}"
                )
            protein_count = sum(
                len(contig.get("proteins", []))
                for contig in contigs
                if isinstance(contig, dict)
            )
            payload_bytes = packet.get("model_payload_bytes")
            if (
                isinstance(payload_bytes, bool)
                or not isinstance(payload_bytes, int)
                or payload_bytes < 0
            ):
                raise RuntimeError(
                    f"Packet calibration received invalid bytes in {genome_id}"
                )
            sampled_payload_bytes += payload_bytes
            sampled_contig_segments += len(contigs)
            if contigs:
                contig_frame_densities.append(
                    math.ceil(payload_bytes / len(contigs))
                )
            if protein_count:
                protein_bearing_payload_bytes += payload_bytes
                sampled_proteins += protein_count
                protein_frame_densities.append(
                    math.ceil(payload_bytes / protein_count)
                )
            if packet.get("target_exceeded") is True:
                target_exceeded_frames += 1
    finally:
        store.close()

    if sampled_proteins:
        mean_bytes_per_protein = (
            protein_bearing_payload_bytes / sampled_proteins
        )
    else:
        mean_bytes_per_protein = None
    mean_bytes_per_contig_segment = (
        sampled_payload_bytes / sampled_contig_segments
        if sampled_contig_segments
        else None
    )
    protein_density_distribution = _distribution(protein_frame_densities)
    contig_density_distribution = _distribution(contig_frame_densities)
    record_floors = minimum_genome_packet_record_bytes()
    proportional_contract = genome_packet_packing_contract(
        max_contigs=None,
        max_proteins=None,
        max_model_payload_bytes=packet_model_payload_bytes,
        all_annotations=packet_all_annotations,
    )
    proportional_contig_limit = int(proportional_contract["max_contigs"])
    proportional_protein_limit = int(proportional_contract["max_proteins"])
    identity = {
        "schema_version": ATLAS_PACKET_CALIBRATION_SCHEMA,
        "selection": "genome protein-count quantiles",
        "sample_genome_count_source": sample_genome_source,
        "requested_genomes": sample_genomes,
        "sampled_genomes": len(sampled_ids),
        "sampled_genome_ids": sampled_ids,
        "calibration_packet_contigs": proportional_contig_limit,
        "calibration_packet_proteins": proportional_protein_limit,
        "calibration_packet_bytes": packet_model_payload_bytes,
        "sampled_nonempty_frames": len(contig_frame_densities),
        "sampled_contig_segments": sampled_contig_segments,
        "sampled_proteins": sampled_proteins,
        "sampled_model_payload_bytes": sampled_payload_bytes,
        "protein_bearing_model_payload_bytes": protein_bearing_payload_bytes,
        "mean_model_payload_bytes_per_protein": (
            round(mean_bytes_per_protein, 3)
            if mean_bytes_per_protein is not None
            else None
        ),
        "mean_model_payload_bytes_per_contig_segment": (
            round(mean_bytes_per_contig_segment, 3)
            if mean_bytes_per_contig_segment is not None
            else None
        ),
        "frame_bytes_per_protein_distribution": protein_density_distribution,
        "frame_bytes_per_contig_segment_distribution": (
            contig_density_distribution
        ),
        "target_model_payload_bytes": packet_model_payload_bytes,
        "all_annotations": packet_all_annotations,
        "payload_proportional_protein_limit": proportional_protein_limit,
        "protein_limit_derivation": (
            "target_model_payload_bytes // "
            f"{record_floors['protein']}-byte serialized schema floor"
        ),
        "payload_proportional_contig_limit": proportional_contig_limit,
        "canonical_contig_segment_floor_bytes": record_floors["contig"],
        "canonical_protein_record_floor_bytes": record_floors["protein"],
        "contig_limit_derivation": (
            "target_model_payload_bytes // "
            f"{record_floors['contig']}-byte serialized schema floor"
        ),
        "target_exceeded_frames": target_exceeded_frames,
        "model_calls_made": 0,
    }
    return {
        **identity,
        "calibration_sha256": _sha256(identity),
    }


def _project_model_invocations(
    counts: list[tuple[Any, ...]],
    *,
    packet_contig_limit: int,
    packet_protein_limit: int,
    packet_model_payload_bytes: int,
    bytes_per_protein: int | float | None,
    bytes_per_contig_segment: int | float | None,
) -> int:
    """Project calls from count limits and measured payload densities."""
    if bytes_per_protein is None or bytes_per_protein <= 0:
        byte_limited_proteins = packet_protein_limit
    else:
        byte_limited_proteins = max(
            1,
            math.floor(packet_model_payload_bytes / bytes_per_protein),
        )
    effective_proteins = max(
        1,
        min(packet_protein_limit, byte_limited_proteins),
    )
    if bytes_per_contig_segment is None or bytes_per_contig_segment <= 0:
        byte_limited_contigs = packet_contig_limit
    else:
        byte_limited_contigs = max(
            1,
            math.floor(
                packet_model_payload_bytes / bytes_per_contig_segment
            ),
        )
    effective_contigs = max(
        1,
        min(packet_contig_limit, byte_limited_contigs),
    )
    total = 0
    for _bin_id, _declared, n_contigs, n_proteins in counts:
        contig_calls = (
            math.ceil(int(n_contigs) / effective_contigs)
            if int(n_contigs)
            else 0
        )
        protein_calls = (
            math.ceil(int(n_proteins) / effective_proteins)
            if int(n_proteins)
            else 0
        )
        total += max(contig_calls, protein_calls)
    return total


def build_atlas_plan(
    db_path: str | Path,
    output_dir: str | Path,
    *,
    packet_model_payload_bytes: int,
    seal_path: str | Path | None = None,
    packet_contig_limit: int | None = None,
    packet_protein_limit: int | None = None,
    packet_all_annotations: bool = True,
    packet_calibration_genomes: int | None = None,
    checkpoint_interval_frames: int = DEFAULT_CHECKPOINT_INTERVAL_FRAMES,
    query_result_limit_bytes: int = DEFAULT_QUERY_RESULT_LIMIT_BYTES,
    threads: int = 4,
    verify_seal: bool = True,
) -> dict[str, Any]:
    """Build one stable genome-owned work unit per live database bin."""
    genome_packet_packing_contract(
        max_contigs=packet_contig_limit,
        max_proteins=packet_protein_limit,
        max_model_payload_bytes=packet_model_payload_bytes,
        all_annotations=packet_all_annotations,
    )
    if checkpoint_interval_frames < 1:
        raise ValueError("checkpoint_interval_frames must be positive")
    if query_result_limit_bytes <= QUERY_RESULT_ENVELOPE_RESERVE_BYTES:
        raise ValueError(
            "query_result_limit_bytes must exceed the HTTP envelope reserve"
        )
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
    resolved_calibration_genomes = (
        int(packet_calibration_genomes)
        if packet_calibration_genomes is not None
        else max(1, math.ceil(math.sqrt(len(counts))))
    )
    packet_calibration = _calibrate_packet_limits(
        database,
        dataset_id=dataset_id,
        counts=counts,
        packet_model_payload_bytes=packet_model_payload_bytes,
        packet_all_annotations=packet_all_annotations,
        sample_genomes=resolved_calibration_genomes,
        sample_genome_source=(
            "explicit"
            if packet_calibration_genomes is not None
            else "ceil(sqrt(live bin count))"
        ),
        threads=threads,
    )
    resolved_contig_limit = (
        int(packet_contig_limit)
        if packet_contig_limit is not None
        else int(packet_calibration["payload_proportional_contig_limit"])
    )
    resolved_protein_limit = (
        int(packet_protein_limit)
        if packet_protein_limit is not None
        else int(packet_calibration["payload_proportional_protein_limit"])
    )
    packet_contract = genome_packet_packing_contract(
        max_contigs=resolved_contig_limit,
        max_proteins=resolved_protein_limit,
        max_model_payload_bytes=packet_model_payload_bytes,
        all_annotations=packet_all_annotations,
    )
    packet_contract_hash = _sha256(packet_contract)
    protein_density = packet_calibration[
        "frame_bytes_per_protein_distribution"
    ]
    contig_density = packet_calibration[
        "frame_bytes_per_contig_segment_distribution"
    ]
    invocation_projection = {
        "method": (
            "per-bin count ceilings using calibrated serialized-payload "
            "density; whole-contig boundary effects await exact census"
        ),
        "at_calibrated_mean_density": _project_model_invocations(
            counts,
            packet_contig_limit=resolved_contig_limit,
            packet_protein_limit=resolved_protein_limit,
            packet_model_payload_bytes=packet_model_payload_bytes,
            bytes_per_protein=packet_calibration[
                "mean_model_payload_bytes_per_protein"
            ],
            bytes_per_contig_segment=packet_calibration[
                "mean_model_payload_bytes_per_contig_segment"
            ],
        ),
        "at_calibrated_p50_frame_density": _project_model_invocations(
            counts,
            packet_contig_limit=resolved_contig_limit,
            packet_protein_limit=resolved_protein_limit,
            packet_model_payload_bytes=packet_model_payload_bytes,
            bytes_per_protein=protein_density["p50"],
            bytes_per_contig_segment=contig_density["p50"],
        ),
        "at_calibrated_p95_frame_density": _project_model_invocations(
            counts,
            packet_contig_limit=resolved_contig_limit,
            packet_protein_limit=resolved_protein_limit,
            packet_model_payload_bytes=packet_model_payload_bytes,
            bytes_per_protein=protein_density["p95"],
            bytes_per_contig_segment=contig_density["p95"],
        ),
        "exact_count_source": "sharur-atlas packet-census",
    }
    plan_identity = {
        "schema_version": ATLAS_PLAN_SCHEMA,
        "dataset_id": dataset_id,
        "packet_version": GENOME_PACKET_VERSION,
        "packet_packing_contract": packet_contract,
        "packet_packing_contract_hash": packet_contract_hash,
        "packet_contig_limit": resolved_contig_limit,
        "packet_protein_limit": resolved_protein_limit,
        "packet_model_payload_bytes": packet_model_payload_bytes,
        "packet_all_annotations": packet_all_annotations,
        "packet_limit_sources": {
            "contigs": (
                "explicit"
                if packet_contig_limit is not None
                else "payload-proportional schema ceiling"
            ),
            "proteins": (
                "explicit"
                if packet_protein_limit is not None
                else "payload-proportional schema ceiling"
            ),
            "model_payload_bytes": "explicit model-facing target",
        },
        "packet_calibration": packet_calibration,
        "model_invocation_projection": invocation_projection,
        "checkpoint_interval_frames": checkpoint_interval_frames,
        "query_result_limit_bytes": query_result_limit_bytes,
        "query_result_envelope_reserve_bytes": (
            QUERY_RESULT_ENVELOPE_RESERVE_BYTES
        ),
        "ownership_unit": "genome",
        "model_call_unit": "one adaptive packet from one exact bins.bin_id",
        "checkpoint_order": "genome-packet cursor",
        "contig_order": "contig_id",
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
                "packet_version": GENOME_PACKET_VERSION,
                "packet_packing_contract_hash": packet_contract_hash,
                "packet_contig_limit": resolved_contig_limit,
                "packet_protein_limit": resolved_protein_limit,
                "packet_model_payload_bytes": packet_model_payload_bytes,
                "packet_all_annotations": packet_all_annotations,
                "checkpoint_key": "atlas_progress",
                "checkpoint_interval_frames": checkpoint_interval_frames,
                "query_result_limit_bytes": query_result_limit_bytes,
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
        "packet_census_file": "packet_census/census.json",
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


def _packet_kwargs(plan: dict[str, Any]) -> dict[str, Any]:
    contract = plan.get("packet_packing_contract")
    if not isinstance(contract, dict):
        raise ValueError("Atlas plan lacks packet_packing_contract")
    if _sha256(contract) != plan.get("packet_packing_contract_hash"):
        raise ValueError("Atlas packet packing-contract hash is invalid")
    expected = genome_packet_packing_contract(
        max_contigs=int(plan["packet_contig_limit"]),
        max_proteins=int(plan["packet_protein_limit"]),
        max_model_payload_bytes=int(plan["packet_model_payload_bytes"]),
        all_annotations=bool(plan["packet_all_annotations"]),
    )
    if contract != expected:
        raise ValueError("Atlas packet packing contract differs from plan fields")
    return {
        "max_contigs": int(contract["max_contigs"]),
        "max_proteins": int(contract["max_proteins"]),
        "max_model_payload_bytes": int(contract["max_model_payload_bytes"]),
        "all_annotations": bool(contract["all_annotations"]),
    }


def _distribution(values: list[int]) -> dict[str, int | float | None]:
    if not values:
        return {
            "count": 0,
            "minimum": None,
            "p50": None,
            "p90": None,
            "p95": None,
            "p99": None,
            "maximum": None,
            "mean": None,
        }
    ordered = sorted(values)

    def percentile(fraction: float) -> int:
        return ordered[round((len(ordered) - 1) * fraction)]

    return {
        "count": len(ordered),
        "minimum": ordered[0],
        "p50": percentile(0.50),
        "p90": percentile(0.90),
        "p95": percentile(0.95),
        "p99": percentile(0.99),
        "maximum": ordered[-1],
        "mean": round(sum(ordered) / len(ordered), 3),
    }


def _census_unit_hash(record: dict[str, Any]) -> str:
    identity = {
        key: value
        for key, value in record.items()
        if key not in {"generated_at", "unit_census_sha256"}
    }
    return _sha256(identity)


def _validate_census_unit(
    record: dict[str, Any],
    *,
    plan: dict[str, Any],
    unit: dict[str, Any],
) -> None:
    expected = {
        "schema_version": ATLAS_PACKET_CENSUS_UNIT_SCHEMA,
        "plan_id": plan["plan_id"],
        "dataset_id": plan["dataset_id"],
        "packet_version": plan["packet_version"],
        "packet_packing_contract_hash": plan["packet_packing_contract_hash"],
        "unit_id": unit["unit_id"],
        "genome_id": unit["genome_id"],
        "expected_contigs": int(unit["n_contigs"]),
        "expected_proteins": int(unit["n_proteins"]),
        "maximum_query_result_bytes": int(plan["query_result_limit_bytes"]),
    }
    for key, value in expected.items():
        if record.get(key) != value:
            raise ValueError(
                f"Packet census unit {unit['unit_id']} has mismatched {key}"
            )
    if record.get("unit_census_sha256") != _census_unit_hash(record):
        raise ValueError(f"Packet census unit {unit['unit_id']} has an invalid hash")
    if record.get("status") not in {"complete", "blocked"}:
        raise ValueError(f"Packet census unit {unit['unit_id']} has invalid status")
    if record.get("observed_contigs") != int(unit["n_contigs"]):
        raise ValueError(f"Packet census unit {unit['unit_id']} missed contigs")
    if record.get("observed_proteins") != int(unit["n_proteins"]):
        raise ValueError(f"Packet census unit {unit['unit_id']} missed proteins")
    if record.get("mixed_bin_frames") != 0:
        raise ValueError(f"Packet census unit {unit['unit_id']} mixed bins")
    if (
        record.get("status") == "complete"
        and record.get("query_result_budget_exceeded_frames") != 0
    ):
        raise ValueError(
            f"Packet census unit {unit['unit_id']} exceeded the query-result budget"
        )
    if record.get("frame_count") != record.get("model_invocation_count"):
        raise ValueError(
            f"Packet census unit {unit['unit_id']} has inconsistent invocation count"
        )


def _census_genome_packets(
    store: Any,
    *,
    plan: dict[str, Any],
    unit: dict[str, Any],
    packet_kwargs: dict[str, Any],
) -> dict[str, Any]:
    genome_id = str(unit["genome_id"])
    cursor: str | None = None
    seen_cursors: set[str] = set()
    frame_payload_bytes: list[int] = []
    frame_http_result_bytes: list[int] = []
    frame_protein_counts: list[int] = []
    frame_contig_segment_counts: list[int] = []
    token_scenarios = {
        "at_2_bytes_per_token": 0,
        "at_3_bytes_per_token": 0,
        "at_4_bytes_per_token": 0,
    }
    contig_state: dict[str, dict[str, int | bool]] = {}
    target_exceeded_frames = 0
    query_result_budget_exceeded_frames = 0
    mixed_bin_frames = 0

    while True:
        result = get_genome_packet(
            store,
            genome_id,
            cursor=cursor,
            **packet_kwargs,
        )
        packet = result.raw
        if packet.get("dataset_id") != plan["dataset_id"]:
            raise RuntimeError(f"Packet census dataset drift in {genome_id}")
        if packet.get("bin_id") != genome_id:
            mixed_bin_frames += 1
            raise RuntimeError(f"Packet census received another bin in {genome_id}")
        if (
            packet.get("packing_contract_hash")
            != plan["packet_packing_contract_hash"]
        ):
            raise RuntimeError(f"Packet census contract drift in {genome_id}")
        model_payload = packet.get("model_payload")
        if not isinstance(model_payload, dict):
            raise RuntimeError(f"Packet census received an invalid payload in {genome_id}")
        contigs = model_payload.get("contigs")
        if not isinstance(contigs, list):
            raise RuntimeError(f"Packet census received invalid contigs in {genome_id}")
        encoded_payload = _canonical_bytes(model_payload)
        if len(encoded_payload) != packet.get("model_payload_bytes"):
            raise RuntimeError(f"Packet byte count is invalid in {genome_id}")
        if hashlib.sha256(encoded_payload).hexdigest() != packet.get(
            "model_payload_sha256"
        ):
            raise RuntimeError(f"Packet payload hash is invalid in {genome_id}")

        if contigs:
            frame_payload_bytes.append(len(encoded_payload))
            operator_payload = json.loads(result.to_json())
            operator_payload_bytes = len(
                json.dumps(
                    operator_payload,
                    default=str,
                    separators=(",", ":"),
                ).encode("utf-8")
            )
            projected_http_bytes = (
                operator_payload_bytes + QUERY_RESULT_ENVELOPE_RESERVE_BYTES
            )
            frame_http_result_bytes.append(projected_http_bytes)
            if projected_http_bytes > int(plan["query_result_limit_bytes"]):
                query_result_budget_exceeded_frames += 1
            frame_protein_counts.append(
                sum(
                    len(contig.get("proteins", []))
                    for contig in contigs
                    if isinstance(contig, dict)
                )
            )
            frame_contig_segment_counts.append(len(contigs))
            scenarios = packet.get("token_scenarios")
            if not isinstance(scenarios, dict):
                raise RuntimeError(
                    f"Packet census received invalid token scenarios in {genome_id}"
                )
            for key in token_scenarios:
                value = scenarios.get(key)
                if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                    raise RuntimeError(
                        f"Packet census received invalid {key} in {genome_id}"
                    )
                token_scenarios[key] += value
            if packet.get("target_exceeded") is True:
                target_exceeded_frames += 1

        receipts = packet.get("coverage_receipts")
        if not isinstance(receipts, list) or len(receipts) != len(contigs):
            raise RuntimeError(f"Packet census receipt mismatch in {genome_id}")
        for contig, receipt in zip(contigs, receipts, strict=True):
            if not isinstance(contig, dict) or not isinstance(receipt, dict):
                raise RuntimeError(f"Packet census record is invalid in {genome_id}")
            if contig.get("bin_id") != genome_id or receipt.get("bin_id") != genome_id:
                mixed_bin_frames += 1
                raise RuntimeError(f"Packet census detected mixed bins in {genome_id}")
            contig_id = contig.get("contig_id")
            if not isinstance(contig_id, str) or receipt.get("contig_id") != contig_id:
                raise RuntimeError(f"Packet census contig identity drift in {genome_id}")
            if any(
                "sequence" in protein
                for protein in contig.get("proteins", [])
                if isinstance(protein, dict)
            ):
                raise RuntimeError(f"Packet census exposed sequence data in {genome_id}")
            total = receipt.get("total_protein_count")
            start = receipt.get("protein_offset_start")
            end = receipt.get("protein_offset_end")
            segment_index = receipt.get("segment_index")
            if any(
                isinstance(value, bool) or not isinstance(value, int) or value < 0
                for value in (total, start, end, segment_index)
            ):
                raise RuntimeError(f"Packet census has invalid offsets in {genome_id}")
            state = contig_state.setdefault(
                contig_id,
                {
                    "total": int(total),
                    "next_offset": 0,
                    "next_segment": 0,
                    "complete": False,
                },
            )
            if state["complete"] is True:
                raise RuntimeError(
                    f"Packet census repeated completed contig {genome_id}/{contig_id}"
                )
            if state["total"] != total:
                raise RuntimeError(
                    f"Packet census changed contig total {genome_id}/{contig_id}"
                )
            if state["next_offset"] != start or state["next_segment"] != segment_index:
                raise RuntimeError(
                    f"Packet census has a gap or overlap in {genome_id}/{contig_id}"
                )
            if end < start or end > total:
                raise RuntimeError(
                    f"Packet census has an invalid range in {genome_id}/{contig_id}"
                )
            if len(contig.get("proteins", [])) != end - start:
                raise RuntimeError(
                    f"Packet census protein count differs in {genome_id}/{contig_id}"
                )
            complete_contig = receipt.get("complete") is True
            if complete_contig != (end == total):
                raise RuntimeError(
                    f"Packet census terminal receipt differs in {genome_id}/{contig_id}"
                )
            state["next_offset"] = int(end)
            state["next_segment"] = int(segment_index) + 1
            state["complete"] = complete_contig

        complete = packet.get("complete") is True
        next_cursor = packet.get("next_cursor")
        if complete:
            if next_cursor is not None:
                raise RuntimeError(
                    f"Packet census terminal frame has a cursor in {genome_id}"
                )
            break
        if not isinstance(next_cursor, str) or not next_cursor:
            raise RuntimeError(f"Packet census lacks continuation in {genome_id}")
        if next_cursor in seen_cursors:
            raise RuntimeError(f"Packet census repeated a cursor in {genome_id}")
        seen_cursors.add(next_cursor)
        cursor = next_cursor

    observed_contigs = len(contig_state)
    observed_proteins = sum(int(state["next_offset"]) for state in contig_state.values())
    if observed_contigs != int(unit["n_contigs"]):
        raise RuntimeError(
            f"Packet census observed {observed_contigs} contigs in {genome_id}; "
            f"expected {int(unit['n_contigs'])}"
        )
    if observed_proteins != int(unit["n_proteins"]):
        raise RuntimeError(
            f"Packet census observed {observed_proteins} proteins in {genome_id}; "
            f"expected {int(unit['n_proteins'])}"
        )
    incomplete_contigs = [
        contig_id
        for contig_id, state in contig_state.items()
        if state["complete"] is not True
    ]
    if incomplete_contigs:
        raise RuntimeError(f"Packet census left incomplete contigs in {genome_id}")
    split_contigs = sum(
        int(state["next_segment"]) > 1 for state in contig_state.values()
    )
    status = (
        "blocked"
        if target_exceeded_frames or query_result_budget_exceeded_frames
        else "complete"
    )
    identity = {
        "schema_version": ATLAS_PACKET_CENSUS_UNIT_SCHEMA,
        "status": status,
        "plan_id": plan["plan_id"],
        "dataset_id": plan["dataset_id"],
        "packet_version": plan["packet_version"],
        "packet_packing_contract_hash": plan["packet_packing_contract_hash"],
        "unit_id": unit["unit_id"],
        "genome_id": genome_id,
        "expected_contigs": int(unit["n_contigs"]),
        "observed_contigs": observed_contigs,
        "expected_proteins": int(unit["n_proteins"]),
        "observed_proteins": observed_proteins,
        "frame_count": len(frame_payload_bytes),
        "model_invocation_count": len(frame_payload_bytes),
        "model_invocation_rule": "one invocation per nonempty bin-scoped frame",
        "total_model_payload_bytes": sum(frame_payload_bytes),
        "maximum_query_result_bytes": int(plan["query_result_limit_bytes"]),
        "query_result_envelope_reserve_bytes": (
            QUERY_RESULT_ENVELOPE_RESERVE_BYTES
        ),
        "query_result_budget_exceeded_frames": (
            query_result_budget_exceeded_frames
        ),
        "token_scenarios": token_scenarios,
        "target_exceeded_frames": target_exceeded_frames,
        "mixed_bin_frames": mixed_bin_frames,
        "split_contigs": split_contigs,
        "frame_payload_bytes": frame_payload_bytes,
        "frame_http_result_bytes": frame_http_result_bytes,
        "frame_protein_counts": frame_protein_counts,
        "frame_contig_segment_counts": frame_contig_segment_counts,
    }
    return {
        **identity,
        "unit_census_sha256": _sha256(identity),
        "generated_at": datetime.now(timezone.utc).isoformat(),
    }


def build_atlas_packet_census(
    plan_dir: str | Path,
    *,
    output_dir: str | Path | None = None,
    workers: int = DEFAULT_CENSUS_WORKERS,
    threads: int = 4,
    memory_limit: str = DEFAULT_CENSUS_MEMORY_LIMIT,
    temp_directory: str | Path | None = None,
    max_temp_directory_size: str = DEFAULT_CENSUS_MAX_TEMP_SIZE,
    resume: bool = True,
    verify_seal: bool = True,
    progress_callback: Callable[[dict[str, Any]], None] | None = None,
) -> dict[str, Any]:
    """Enumerate exact bin-scoped model packets while making zero model calls."""
    if threads < 1:
        raise ValueError("threads must be positive")
    if workers < 1:
        raise ValueError("workers must be positive")
    effective_workers = min(workers, threads)
    root = Path(plan_dir).expanduser().resolve()
    plan, units = load_atlas_plan(root)
    packet_kwargs = _packet_kwargs(plan)
    database = Path(str(plan["database"])).expanduser().resolve()
    seal_path = Path(str(plan["seal_path"])).expanduser().resolve()
    dataset_id, _seal_strength = _sealed_identity(
        database,
        seal_path,
        verify=verify_seal,
    )
    if dataset_id != plan["dataset_id"]:
        raise DatasetSealError("Atlas plan dataset_id differs from the live seal")
    census_root = (
        Path(output_dir).expanduser().resolve()
        if output_dir is not None
        else root / "packet_census"
    )
    unit_root = census_root / "units"
    records_by_unit: dict[str, dict[str, Any]] = {}
    pending_units: list[dict[str, Any]] = []
    for unit in units:
        unit_path = unit_root / f"{unit['unit_id']}.json"
        if resume and unit_path.is_file():
            candidate = _read_json_object(unit_path)
            _validate_census_unit(candidate, plan=plan, unit=unit)
            records_by_unit[str(unit["unit_id"])] = candidate
        else:
            pending_units.append(unit)
    completed_units = len(records_by_unit)
    cumulative_frame_count = sum(
        int(record["frame_count"]) for record in records_by_unit.values()
    )
    if progress_callback is not None and completed_units:
        progress_callback(
            {
                "completed_units": completed_units,
                "total_units": len(units),
                "reused_units": completed_units,
                "computed_units": 0,
                "unit_id": None,
                "genome_id": None,
                "status": "resumed",
                "unit_frame_count": None,
                "cumulative_frame_count": cumulative_frame_count,
            }
        )

    resolved_temp_directory = (
        Path(temp_directory).expanduser().resolve()
        if temp_directory is not None
        else census_root / "spill"
    )
    if pending_units:
        runtime = ReadOnlyDuckDBRuntime(
            database,
            threads=threads,
            memory_limit=memory_limit,
            temp_directory=resolved_temp_directory,
            max_temp_directory_size=max_temp_directory_size,
            dataset_id=dataset_id,
        )
        runtime.open()
        try:

            def census_one(unit: dict[str, Any]) -> dict[str, Any]:
                record = _census_genome_packets(
                    runtime.cursor_store(),
                    plan=plan,
                    unit=unit,
                    packet_kwargs=packet_kwargs,
                )
                write_json_atomic(
                    record,
                    unit_root / f"{unit['unit_id']}.json",
                )
                return record

            with ThreadPoolExecutor(
                max_workers=min(effective_workers, len(pending_units)),
                thread_name_prefix="atlas-packet-census",
            ) as executor:
                futures = {
                    executor.submit(census_one, unit): unit
                    for unit in pending_units
                }
                for computed_units, future in enumerate(
                    as_completed(futures),
                    start=1,
                ):
                    unit = futures[future]
                    record = future.result()
                    records_by_unit[str(unit["unit_id"])] = record
                    cumulative_frame_count += int(record["frame_count"])
                    if progress_callback is not None:
                        progress_callback(
                            {
                                "completed_units": completed_units
                                + computed_units,
                                "total_units": len(units),
                                "reused_units": completed_units,
                                "computed_units": computed_units,
                                "unit_id": unit["unit_id"],
                                "genome_id": unit["genome_id"],
                                "status": record["status"],
                                "unit_frame_count": record["frame_count"],
                                "cumulative_frame_count": (
                                    cumulative_frame_count
                                ),
                            }
                        )
        finally:
            runtime.close()

    unit_records = [
        records_by_unit[str(unit["unit_id"])]
        for unit in units
    ]

    frame_payload_bytes = [
        int(value)
        for record in unit_records
        for value in record["frame_payload_bytes"]
    ]
    frame_http_result_bytes = [
        int(value)
        for record in unit_records
        for value in record["frame_http_result_bytes"]
    ]
    frame_protein_counts = [
        int(value)
        for record in unit_records
        for value in record["frame_protein_counts"]
    ]
    frame_contig_segment_counts = [
        int(value)
        for record in unit_records
        for value in record["frame_contig_segment_counts"]
    ]
    blocked_units = [
        str(record["unit_id"])
        for record in unit_records
        if record["status"] != "complete"
    ]
    total_token_scenarios = {
        key: sum(int(record["token_scenarios"][key]) for record in unit_records)
        for key in (
            "at_2_bytes_per_token",
            "at_3_bytes_per_token",
            "at_4_bytes_per_token",
        )
    }
    unit_record_digest = _sha256(
        [
            {
                "unit_id": record["unit_id"],
                "unit_census_sha256": record["unit_census_sha256"],
            }
            for record in unit_records
        ]
    )
    target_exceeded_frames = sum(
        int(record["target_exceeded_frames"]) for record in unit_records
    )
    query_result_budget_exceeded_frames = sum(
        int(record["query_result_budget_exceeded_frames"])
        for record in unit_records
    )
    mixed_bin_frames = sum(
        int(record["mixed_bin_frames"]) for record in unit_records
    )
    status = (
        "complete"
        if not blocked_units
        and target_exceeded_frames == 0
        and query_result_budget_exceeded_frames == 0
        and mixed_bin_frames == 0
        else "blocked"
    )
    identity = {
        "schema_version": ATLAS_PACKET_CENSUS_SCHEMA,
        "status": status,
        "plan_id": plan["plan_id"],
        "dataset_id": plan["dataset_id"],
        "packet_version": plan["packet_version"],
        "packet_packing_contract": plan["packet_packing_contract"],
        "packet_packing_contract_hash": plan["packet_packing_contract_hash"],
        "census_resources": {
            "requested_workers": workers,
            "workers": effective_workers,
            "duckdb_threads": threads,
            "memory_limit": memory_limit,
            "temp_directory": str(resolved_temp_directory),
            "max_temp_directory_size": max_temp_directory_size,
        },
        "model_calls_made": 0,
        "model_invocation_rule": "one invocation per nonempty bin-scoped frame",
        "model_invocation_count": len(frame_payload_bytes),
        "prompt_and_output_tokens_included": False,
        "token_accounting": (
            "deterministic serialized-payload scenarios; exact provider "
            "tokenizers are executor-specific"
        ),
        "token_scenarios": total_token_scenarios,
        "total_model_payload_bytes": sum(frame_payload_bytes),
        "maximum_query_result_bytes": int(plan["query_result_limit_bytes"]),
        "query_result_envelope_reserve_bytes": (
            QUERY_RESULT_ENVELOPE_RESERVE_BYTES
        ),
        "query_result_budget_exceeded_frames": (
            query_result_budget_exceeded_frames
        ),
        "unit_count": len(unit_records),
        "expected_units": int(plan["unit_count"]),
        "observed_contigs": sum(
            int(record["observed_contigs"]) for record in unit_records
        ),
        "expected_contigs": int(plan["total_contigs"]),
        "observed_proteins": sum(
            int(record["observed_proteins"]) for record in unit_records
        ),
        "expected_proteins": int(plan["total_proteins"]),
        "mixed_bin_frames": mixed_bin_frames,
        "maximum_bins_per_frame": 1 if frame_payload_bytes else 0,
        "target_exceeded_frames": target_exceeded_frames,
        "split_contigs": sum(int(record["split_contigs"]) for record in unit_records),
        "blocked_units": blocked_units,
        "frame_payload_byte_distribution": _distribution(frame_payload_bytes),
        "frame_http_result_byte_distribution": _distribution(
            frame_http_result_bytes
        ),
        "frame_protein_distribution": _distribution(frame_protein_counts),
        "frame_contig_segment_distribution": _distribution(
            frame_contig_segment_counts
        ),
        "unit_records_directory": "units",
        "unit_records_sha256": unit_record_digest,
    }
    summary = {
        **identity,
        "census_sha256": _sha256(identity),
        "generated_at": datetime.now(timezone.utc).isoformat(),
    }
    write_json_atomic(summary, census_root / "census.json")
    return {
        "census_path": str(census_root / "census.json"),
        **summary,
    }


def verify_atlas_packet_census(
    plan_dir: str | Path,
    *,
    census_dir: str | Path | None = None,
    deep: bool = False,
) -> dict[str, Any]:
    """Verify the launch-gating zero-model-call packet census."""
    root = Path(plan_dir).expanduser().resolve()
    plan, units = load_atlas_plan(root)
    census_root = (
        Path(census_dir).expanduser().resolve()
        if census_dir is not None
        else root / "packet_census"
    )
    summary = _read_json_object(census_root / "census.json")
    expected = {
        "schema_version": ATLAS_PACKET_CENSUS_SCHEMA,
        "plan_id": plan["plan_id"],
        "dataset_id": plan["dataset_id"],
        "packet_version": plan["packet_version"],
        "packet_packing_contract_hash": plan["packet_packing_contract_hash"],
        "model_calls_made": 0,
        "unit_count": int(plan["unit_count"]),
        "expected_units": int(plan["unit_count"]),
        "observed_contigs": int(plan["total_contigs"]),
        "expected_contigs": int(plan["total_contigs"]),
        "observed_proteins": int(plan["total_proteins"]),
        "expected_proteins": int(plan["total_proteins"]),
        "mixed_bin_frames": 0,
        "maximum_query_result_bytes": int(plan["query_result_limit_bytes"]),
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            raise ValueError(f"Atlas packet census has mismatched {key}")
    identity = {
        key: value
        for key, value in summary.items()
        if key not in {"census_sha256", "generated_at"}
    }
    if summary.get("census_sha256") != _sha256(identity):
        raise ValueError("Atlas packet census hash is invalid")
    if summary.get("status") != "complete":
        raise RuntimeError(
            "Atlas packet census blocks enqueue: "
            f"{summary.get('target_exceeded_frames', 0)} target-exceeded frames, "
            f"{summary.get('query_result_budget_exceeded_frames', 0)} "
            "query-result overflows, "
            f"{len(summary.get('blocked_units', []))} blocked units"
        )
    if summary.get("target_exceeded_frames") != 0:
        raise RuntimeError("Atlas packet census contains target-exceeded frames")
    if summary.get("query_result_budget_exceeded_frames") != 0:
        raise RuntimeError("Atlas packet census exceeds the query-result budget")
    if deep:
        unit_records: list[dict[str, Any]] = []
        for unit in units:
            path = census_root / "units" / f"{unit['unit_id']}.json"
            record = _read_json_object(path)
            _validate_census_unit(record, plan=plan, unit=unit)
            unit_records.append(record)
        digest = _sha256(
            [
                {
                    "unit_id": record["unit_id"],
                    "unit_census_sha256": record["unit_census_sha256"],
                }
                for record in unit_records
            ]
        )
        if digest != summary.get("unit_records_sha256"):
            raise ValueError("Atlas packet census unit-record digest is invalid")
    return {
        **summary,
        "verification": "passed",
        "deep_verified": deep,
        "census_path": str(census_root / "census.json"),
    }


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
    packet_census = verify_atlas_packet_census(plan_dir)
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
            "packet_census_sha256": packet_census["census_sha256"],
            "model_invocation_count": packet_census["model_invocation_count"],
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
                    "packet_census_path": packet_census["census_path"],
                    "packet_census_sha256": packet_census["census_sha256"],
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


def frame_coverage_receipt(packet: dict[str, Any]) -> dict[str, Any] | None:
    """Project one genome packet into a compact, model-call coverage receipt."""
    model_payload = packet.get("model_payload")
    if not isinstance(model_payload, dict):
        raise ValueError("Genome packet lacks model_payload")
    contigs = model_payload.get("contigs")
    if not isinstance(contigs, list):
        raise ValueError("Genome packet model_payload lacks contigs")
    if not contigs:
        return None
    segments = packet.get("coverage_receipts")
    if not isinstance(segments, list) or len(segments) != len(contigs):
        raise ValueError("Genome packet coverage receipts differ from contigs")
    return {
        "frame_id": packet.get("frame_id"),
        "frame_index": packet.get("frame_index"),
        "bin_id": packet.get("bin_id"),
        "packet_packing_contract_hash": packet.get(
            "packing_contract_hash"
        ),
        "model_payload_sha256": packet.get("model_payload_sha256"),
        "model_payload_bytes": packet.get("model_payload_bytes"),
        "contig_ids": [contig.get("contig_id") for contig in contigs],
        "segments": segments,
        "protein_count": sum(
            len(contig.get("proteins", []))
            for contig in contigs
            if isinstance(contig, dict)
        ),
        "target_exceeded": packet.get("target_exceeded") is True,
    }


def build_genome_coverage_manifest(
    unit: dict[str, Any],
    contigs: list[dict[str, Any]],
    frames: list[dict[str, Any]],
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
        "packet_packing_contract_hash",
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
        segment_count = record.get("segment_count")
        if segment_count is None:
            segment_count = record.get("packet_count")
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
            isinstance(segment_count, bool)
            or not isinstance(segment_count, int)
            or segment_count < 0
        ):
            raise ValueError("Coverage segment_count must be a non-negative integer")
        complete = record.get("complete") is True
        if not complete:
            errors.append(f"Contig {contig_id} lacks a terminal packet")
        normalized.append(
            {
                "contig_id": contig_id,
                "protein_count": protein_count,
                "segment_count": segment_count,
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
    normalized_frames: list[dict[str, Any]] = []
    seen_frames: set[str] = set()
    expected_frame_index = 0
    known_contigs = {record["contig_id"] for record in normalized}
    segment_state: dict[str, dict[str, int | bool]] = {}
    for frame in sorted(frames, key=lambda value: value.get("frame_index", -1)):
        frame_id = frame.get("frame_id")
        frame_index = frame.get("frame_index")
        bin_id = frame.get("bin_id")
        payload_hash = frame.get("model_payload_sha256")
        payload_bytes = frame.get("model_payload_bytes")
        protein_count = frame.get("protein_count")
        contig_ids = frame.get("contig_ids")
        segments = frame.get("segments")
        if not isinstance(frame_id, str) or not frame_id:
            raise ValueError("Coverage frame_id must be a non-empty string")
        if frame_id in seen_frames:
            raise ValueError(f"Coverage repeats frame_id {frame_id!r}")
        seen_frames.add(frame_id)
        if (
            isinstance(frame_index, bool)
            or not isinstance(frame_index, int)
            or frame_index < 0
        ):
            raise ValueError("Coverage frame_index must be a non-negative integer")
        if frame_index != expected_frame_index:
            errors.append(
                f"Frame index {frame_index} follows {expected_frame_index - 1}"
            )
        expected_frame_index = frame_index + 1
        if bin_id != unit["genome_id"]:
            errors.append(f"Frame {frame_id} belongs to bin {bin_id}")
        if (
            frame.get("packet_packing_contract_hash")
            != unit["packet_packing_contract_hash"]
        ):
            errors.append(f"Frame {frame_id} uses another packing contract")
        if not isinstance(payload_hash, str) or len(payload_hash) != 64:
            raise ValueError("Coverage model_payload_sha256 must be a SHA-256 hex")
        if (
            isinstance(payload_bytes, bool)
            or not isinstance(payload_bytes, int)
            or payload_bytes < 0
        ):
            raise ValueError(
                "Coverage model_payload_bytes must be a non-negative integer"
            )
        if (
            isinstance(protein_count, bool)
            or not isinstance(protein_count, int)
            or protein_count < 0
        ):
            raise ValueError("Coverage frame protein_count must be non-negative")
        if (
            not isinstance(contig_ids, list)
            or not contig_ids
            or any(not isinstance(value, str) or not value for value in contig_ids)
        ):
            raise ValueError("Coverage frame contig_ids must be non-empty strings")
        if not isinstance(segments, list) or len(segments) != len(contig_ids):
            raise ValueError("Coverage frame segments must match contig_ids")
        if any(not isinstance(segment, dict) for segment in segments):
            raise ValueError("Coverage frame segment must be an object")
        if [segment.get("contig_id") for segment in segments] != contig_ids:
            raise ValueError("Coverage frame segment order differs from contig_ids")
        unknown_contigs = sorted(set(contig_ids) - known_contigs)
        if unknown_contigs:
            errors.append(
                f"Frame {frame_id} references unknown contigs: "
                f"{', '.join(unknown_contigs[:5])}"
            )
        if frame.get("target_exceeded") is True:
            errors.append(f"Frame {frame_id} exceeds the model-payload target")
        normalized_segments: list[dict[str, Any]] = []
        frame_segment_proteins = 0
        for segment in segments:
            contig_id = segment.get("contig_id")
            segment_bin_id = segment.get("bin_id")
            total = segment.get("total_protein_count")
            start = segment.get("protein_offset_start")
            end = segment.get("protein_offset_end")
            segment_index = segment.get("segment_index")
            if segment_bin_id != unit["genome_id"]:
                errors.append(
                    f"Frame {frame_id} segment {contig_id} belongs to "
                    f"bin {segment_bin_id}"
                )
            if any(
                isinstance(value, bool) or not isinstance(value, int) or value < 0
                for value in (total, start, end, segment_index)
            ):
                raise ValueError(
                    "Coverage frame segment offsets must be non-negative integers"
                )
            state = segment_state.setdefault(
                str(contig_id),
                {
                    "total": int(total),
                    "next_offset": 0,
                    "next_segment": 0,
                    "complete": False,
                },
            )
            if state["complete"] is True:
                errors.append(
                    f"Frame {frame_id} repeats completed contig {contig_id}"
                )
            if state["total"] != total:
                errors.append(f"Frame {frame_id} changes total for {contig_id}")
            if state["next_offset"] != start or state["next_segment"] != segment_index:
                errors.append(
                    f"Frame {frame_id} has a gap or overlap for {contig_id}"
                )
            if end < start or end > total:
                errors.append(f"Frame {frame_id} has an invalid range for {contig_id}")
            complete_segment = segment.get("complete") is True
            if complete_segment != (end == total):
                errors.append(
                    f"Frame {frame_id} has an invalid terminal flag for {contig_id}"
                )
            state["next_offset"] = int(end)
            state["next_segment"] = int(segment_index) + 1
            state["complete"] = complete_segment
            frame_segment_proteins += int(end) - int(start)
            normalized_segments.append(
                {
                    "bin_id": segment_bin_id,
                    "contig_id": contig_id,
                    "total_protein_count": total,
                    "protein_offset_start": start,
                    "protein_offset_end": end,
                    "segment_index": segment_index,
                    "complete": complete_segment,
                }
            )
        if frame_segment_proteins != protein_count:
            errors.append(f"Frame {frame_id} protein count differs from its segments")
        normalized_frames.append(
            {
                "frame_id": frame_id,
                "frame_index": frame_index,
                "bin_id": bin_id,
                "packet_packing_contract_hash": frame[
                    "packet_packing_contract_hash"
                ],
                "model_payload_sha256": payload_hash,
                "model_payload_bytes": payload_bytes,
                "contig_ids": contig_ids,
                "segments": normalized_segments,
                "protein_count": protein_count,
                "target_exceeded": frame.get("target_exceeded") is True,
            }
        )
    if sum(len(frame["contig_ids"]) for frame in normalized_frames) != sum(
        record["segment_count"] for record in normalized
    ):
        errors.append("Frame contig-segment total differs from contig receipts")
    if sum(frame["protein_count"] for frame in normalized_frames) != observed_proteins:
        errors.append("Frame protein total differs from contig receipts")
    for record in normalized:
        state = segment_state.get(record["contig_id"])
        if state is None:
            errors.append(f"Contig {record['contig_id']} has no frame segment")
            continue
        if state["total"] != record["protein_count"]:
            errors.append(
                f"Contig {record['contig_id']} protein total differs from frames"
            )
        if state["next_offset"] != record["protein_count"]:
            errors.append(
                f"Contig {record['contig_id']} terminal offset differs from total"
            )
        if state["next_segment"] != record["segment_count"]:
            errors.append(
                f"Contig {record['contig_id']} segment count differs from frames"
            )
        if state["complete"] is not record["complete"]:
            errors.append(
                f"Contig {record['contig_id']} completion differs from frames"
            )
    coverage_identity = {
        "schema_version": ATLAS_COVERAGE_SCHEMA,
        "unit_id": unit["unit_id"],
        "plan_id": unit["plan_id"],
        "dataset_id": unit["dataset_id"],
        "genome_id": unit["genome_id"],
        "packet_version": unit["packet_version"],
        "packet_packing_contract_hash": unit[
            "packet_packing_contract_hash"
        ],
        "expected_n_contigs": int(unit["n_contigs"]),
        "observed_n_contigs": observed_contigs,
        "expected_n_proteins": int(unit["n_proteins"]),
        "observed_n_proteins": observed_proteins,
        "model_frame_count": len(normalized_frames),
        "contig_segment_count": sum(
            record["segment_count"] for record in normalized
        ),
        "model_payload_bytes": sum(
            frame["model_payload_bytes"] for frame in normalized_frames
        ),
        "frames": normalized_frames,
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
    frames: list[dict[str, Any]],
    output_path: str | Path,
) -> dict[str, Any]:
    """Atomically write one per-genome coverage manifest."""
    manifest = build_genome_coverage_manifest(unit, contigs, frames)
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
        for key in (
            "unit_id",
            "plan_id",
            "dataset_id",
            "genome_id",
            "packet_version",
            "packet_packing_contract_hash",
        ):
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
    "ATLAS_PACKET_CALIBRATION_SCHEMA",
    "ATLAS_PACKET_CENSUS_SCHEMA",
    "ATLAS_PACKET_CENSUS_UNIT_SCHEMA",
    "ATLAS_PLAN_SCHEMA",
    "ATLAS_REVIEW_OUTPUT_SCHEMA",
    "ATLAS_UNIT_SCHEMA",
    "DEFAULT_CENSUS_MAX_TEMP_SIZE",
    "DEFAULT_CENSUS_MEMORY_LIMIT",
    "DEFAULT_CENSUS_WORKERS",
    "DEFAULT_CHECKPOINT_INTERVAL_FRAMES",
    "DEFAULT_QUERY_RESULT_LIMIT_BYTES",
    "QUERY_RESULT_ENVELOPE_RESERVE_BYTES",
    "WORK_WEIGHT_FORMULA",
    "build_atlas_packet_census",
    "build_atlas_plan",
    "build_genome_coverage_manifest",
    "enqueue_atlas_plan",
    "frame_coverage_receipt",
    "load_atlas_plan",
    "verify_atlas_coverage",
    "verify_atlas_packet_census",
    "write_genome_coverage_manifest",
]
