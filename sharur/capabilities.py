"""Typed, machine-readable Sharur dataset and runtime preflight."""

from __future__ import annotations

import json
import os
import shutil
import sqlite3
import subprocess
import sys
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import TYPE_CHECKING, Any

from sharur.storage.duckdb_store import DuckDBStore
from sharur.storage.schema import SCHEMA_VERSION
from sharur.storage.vector_store import inspect_vector_index


if TYPE_CHECKING:
    from collections.abc import Iterable


class CapabilityState(str, Enum):
    """The four states every preflight capability must use."""

    available = "available"
    unavailable = "unavailable"
    stale = "stale"
    failed = "failed"


_STATE_SEVERITY = {
    CapabilityState.available: 0,
    CapabilityState.unavailable: 1,
    CapabilityState.stale: 2,
    CapabilityState.failed: 3,
}


@dataclass(frozen=True)
class Capability:
    """One typed capability check."""

    capability_id: str
    state: CapabilityState
    summary: str
    required: bool = False
    evidence: dict[str, Any] = field(default_factory=dict)
    remediation: str | None = None

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["state"] = self.state.value
        return payload


@dataclass(frozen=True)
class CapabilityBrief:
    """Complete operator-facing preflight for one dataset."""

    dataset_path: str
    generated_at: str
    overall_state: CapabilityState
    capabilities: tuple[Capability, ...]

    @property
    def ready(self) -> bool:
        return self.overall_state == CapabilityState.available

    def get(self, capability_id: str) -> Capability:
        for capability in self.capabilities:
            if capability.capability_id == capability_id:
                return capability
        raise KeyError(capability_id)

    def to_dict(self) -> dict[str, Any]:
        return {
            "dataset_path": self.dataset_path,
            "generated_at": self.generated_at,
            "overall_state": self.overall_state.value,
            "ready": self.ready,
            "capabilities": [item.to_dict() for item in self.capabilities],
        }

    def to_json(self, *, indent: int | None = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent, default=str)

    def to_markdown(self) -> str:
        lines = [
            "# Sharur Capability Brief",
            "",
            f"- Dataset: `{self.dataset_path}`",
            f"- Overall: **{self.overall_state.value}**",
            f"- Generated: {self.generated_at}",
            "",
            "| Capability | State | Required | Summary |",
            "|---|---|---:|---|",
        ]
        for item in self.capabilities:
            summary = item.summary.replace("|", "\\|")
            lines.append(
                f"| `{item.capability_id}` | {item.state.value} | "
                f"{'yes' if item.required else 'no'} | {summary} |"
            )
        remediations = [item for item in self.capabilities if item.remediation is not None]
        if remediations:
            lines.extend(["", "## Remediation"])
            for item in remediations:
                lines.append(f"- `{item.capability_id}`: {item.remediation}")
        return "\n".join(lines)


def _overall_state(capabilities: Iterable[Capability]) -> CapabilityState:
    required = [item.state for item in capabilities if item.required]
    if not required:
        return CapabilityState.available
    return max(required, key=lambda state: _STATE_SEVERITY[state])


def _table_catalog(store: DuckDBStore) -> dict[str, set[str]]:
    rows = store.execute(
        """
        SELECT table_name, column_name
        FROM information_schema.columns
        WHERE table_schema = 'main'
        ORDER BY table_name, ordinal_position
        """
    )
    result: dict[str, set[str]] = {}
    for table_name, column_name in rows:
        result.setdefault(table_name, set()).add(column_name)
    return result


def _count(store: DuckDBStore, table_name: str) -> int:
    escaped = table_name.replace('"', '""')
    return int(store.execute(f'SELECT COUNT(*) FROM "{escaped}"')[0][0])


def _discover_embedding_h5(db_path: Path) -> Path | None:
    for directory in ("embeddings", "stage06_embeddings"):
        candidate = db_path.parent / directory / "protein_embeddings.h5"
        if candidate.is_file():
            return candidate
    return None


def _database_checks(
    db_path: Path,
) -> tuple[list[Capability], DuckDBStore | None, dict[str, set[str]], int, int, int]:
    if not db_path.is_file():
        return (
            [
                Capability(
                    "dataset",
                    CapabilityState.unavailable,
                    "DuckDB file does not exist.",
                    required=True,
                    evidence={"path": str(db_path)},
                    remediation="Build the dataset with sharur-ingest.",
                )
            ],
            None,
            {},
            0,
            0,
            0,
        )

    store = DuckDBStore(db_path, read_only=True)
    try:
        tables = _table_catalog(store)
    except Exception as exc:
        store.close()
        return (
            [
                Capability(
                    "dataset",
                    CapabilityState.failed,
                    "DuckDB could not be opened or inspected.",
                    required=True,
                    evidence={"error": f"{type(exc).__name__}: {exc}"},
                    remediation="Inspect the database file and rebuild if it is corrupt.",
                )
            ],
            None,
            {},
            0,
            0,
            0,
        )

    capabilities: list[Capability] = []
    core_tables = {"bins", "contigs", "proteins", "annotations"}
    missing = sorted(core_tables - tables.keys())
    if missing:
        capabilities.append(
            Capability(
                "dataset",
                CapabilityState.stale,
                "DuckDB is readable but the core schema is incomplete.",
                required=True,
                evidence={"missing_tables": missing, "tables": sorted(tables)},
                remediation="Run the current Stage 07 builder or rebuild the dataset.",
            )
        )
        protein_count = 0
        sequence_protein_count = 0
        invalid_sequence_count = 0
    else:
        protein_count = _count(store, "proteins")
        sequence_protein_count = int(
            store.execute(
                """
                SELECT COUNT(*)
                FROM proteins
                WHERE sequence IS NOT NULL AND sequence != ''
                """
            )[0][0]
        )
        invalid_sequence_count = protein_count - sequence_protein_count
        if protein_count == 0:
            state = CapabilityState.unavailable
            summary = "Core DuckDB tables exist but contain no proteins."
            remediation = "Ingest non-empty assemblies before analysis."
        elif invalid_sequence_count:
            state = CapabilityState.stale
            summary = (
                "Core DuckDB tables are readable, but "
                f"{invalid_sequence_count} protein rows have no sequence."
            )
            remediation = (
                "Repair or regenerate the Stage 03 protein FASTAs, then rebuild the knowledge base."
            )
        else:
            state = CapabilityState.available
            summary = "Core DuckDB tables and protein sequences are readable."
            remediation = None
        capabilities.append(
            Capability(
                "dataset",
                state,
                summary,
                required=True,
                evidence={
                    "protein_count": protein_count,
                    "sequence_protein_count": sequence_protein_count,
                    "invalid_sequence_rows": invalid_sequence_count,
                    "bin_count": _count(store, "bins"),
                    "contig_count": _count(store, "contigs"),
                },
                remediation=remediation,
            )
        )

    try:
        version_rows = store.execute("SELECT MAX(version) FROM schema_version")
        version = int(version_rows[0][0]) if version_rows and version_rows[0][0] is not None else 0
        schema_state = (
            CapabilityState.available if version == SCHEMA_VERSION else CapabilityState.stale
        )
        summary = (
            f"Schema version {version} matches code."
            if schema_state == CapabilityState.available
            else f"Schema version {version} does not match code version {SCHEMA_VERSION}."
        )
        capabilities.append(
            Capability(
                "schema",
                schema_state,
                summary,
                required=True,
                evidence={"database_version": version, "code_version": SCHEMA_VERSION},
                remediation=(
                    None
                    if schema_state == CapabilityState.available
                    else "Open the database writable once to run migrations, or rebuild it."
                ),
            )
        )
    except Exception as exc:
        capabilities.append(
            Capability(
                "schema",
                CapabilityState.failed,
                "Schema version could not be read.",
                required=True,
                evidence={"error": f"{type(exc).__name__}: {exc}"},
                remediation="Run schema migrations on a writable database copy.",
            )
        )

    return (
        capabilities,
        store,
        tables,
        protein_count,
        sequence_protein_count,
        invalid_sequence_count,
    )


def _annotation_checks(
    store: DuckDBStore,
    tables: dict[str, set[str]],
) -> list[Capability]:
    capabilities: list[Capability] = []
    try:
        sources = [
            {"source": source, "rows": int(rows), "proteins": int(proteins)}
            for source, rows, proteins in store.execute(
                """
                SELECT source, COUNT(*) AS rows, COUNT(DISTINCT protein_id) AS proteins
                FROM annotations
                GROUP BY source
                ORDER BY source
                """
            )
        ]
        capabilities.append(
            Capability(
                "annotation_sources",
                (CapabilityState.available if sources else CapabilityState.unavailable),
                (
                    f"Discovered {len(sources)} live sources in the annotations table."
                    if sources
                    else "No rows are present in the annotations table."
                ),
                evidence={
                    "table": "annotations",
                    "claim_level": "mixed",
                    "interpretation": (
                        "Sources may contain raw observations or caller-derived annotation "
                        "rows; use structured_callers for emitted named calls."
                    ),
                    "sources": sources,
                },
                remediation=(
                    None
                    if sources
                    else "Run the configured annotation stage and rebuild the DuckDB."
                ),
            )
        )
    except Exception as exc:
        capabilities.append(
            Capability(
                "annotation_sources",
                CapabilityState.failed,
                "Annotations-table sources could not be enumerated.",
                evidence={"error": f"{type(exc).__name__}: {exc}"},
            )
        )

    structured: list[dict[str, Any]] = []
    for table_name, columns in sorted(tables.items()):
        has_system_signature = {
            "system_id",
            "system_type",
        } <= columns and table_name != "system_proteins"
        has_locus_signature = {"locus_id", "locus_type"} <= columns
        looks_like_called_loci = table_name.endswith("_loci") and (
            "locus_id" in columns or "cluster_id" in columns
        )
        if not (has_system_signature or has_locus_signature or looks_like_called_loci):
            continue
        entry: dict[str, Any] = {
            "table": table_name,
            "rows": _count(store, table_name),
            "columns": sorted(columns),
        }
        if "locus_type" in columns and entry["rows"]:
            escaped = table_name.replace('"', '""')
            entry["emitted_types"] = [
                {"value": value, "rows": int(count)}
                for value, count in store.execute(
                    f'SELECT locus_type, COUNT(*) FROM "{escaped}" '
                    "GROUP BY locus_type ORDER BY locus_type"
                )
            ]
        if "system_type" in columns and entry["rows"]:
            escaped = table_name.replace('"', '""')
            entry["emitted_types"] = [
                {"value": value, "rows": int(count)}
                for value, count in store.execute(
                    f'SELECT system_type, COUNT(*) FROM "{escaped}" '
                    "GROUP BY system_type ORDER BY system_type"
                )
            ]
        structured.append(entry)

    populated_resources = [item for item in structured if item["rows"]]
    capabilities.append(
        Capability(
            "structured_callers",
            (CapabilityState.available if populated_resources else CapabilityState.unavailable),
            (
                f"Discovered {len(populated_resources)} populated structured caller resources."
                if populated_resources
                else "No emitted calls were found in live structured caller tables."
            ),
            evidence={"resources": structured},
            remediation=(
                None
                if populated_resources
                else "Run the relevant purpose-built callers; domain hits alone do not populate this capability."
            ),
        )
    )
    return capabilities


def _semantic_checks(
    store: DuckDBStore,
    tables: dict[str, set[str]],
    protein_count: int,
) -> list[Capability]:
    checks: list[Capability] = []
    required_v2 = {"semantic_atoms", "semantic_state", "semantic_terms"}
    missing_v2 = sorted(required_v2 - tables.keys())
    if missing_v2:
        checks.append(
            Capability(
                "semantic_v2",
                CapabilityState.unavailable,
                "V2 semantic tables are missing.",
                required=True,
                evidence={"missing_tables": missing_v2},
                remediation="Run current schema migrations and compute V2 predicates.",
            )
        )
    else:
        state_count = _count(store, "semantic_state")
        atom_count = _count(store, "semantic_atoms")
        term_count = _count(store, "semantic_terms")
        complete = protein_count > 0 and state_count == protein_count and term_count > 0
        checks.append(
            Capability(
                "semantic_v2",
                (CapabilityState.available if complete else CapabilityState.stale),
                (
                    "V2 semantic state covers every protein."
                    if complete
                    else "V2 semantic state/materialized terms are incomplete."
                ),
                required=True,
                evidence={
                    "proteins": protein_count,
                    "states": state_count,
                    "atoms": atom_count,
                    "terms": term_count,
                    "coverage": (state_count / protein_count if protein_count else 0.0),
                },
                remediation=(
                    None
                    if complete
                    else "Run `sharur compute-predicates --db DATASET/sharur.duckdb`."
                ),
            )
        )

    if "protein_predicates" not in tables:
        checks.append(
            Capability(
                "predicate_compatibility",
                CapabilityState.unavailable,
                "Compatibility predicate table is missing.",
                required=True,
                remediation="Run schema migrations and compute V2 predicates.",
            )
        )
    else:
        predicate_count = _count(store, "protein_predicates")
        complete = protein_count > 0 and predicate_count == protein_count
        checks.append(
            Capability(
                "predicate_compatibility",
                (CapabilityState.available if complete else CapabilityState.stale),
                (
                    "Compatibility predicates cover every protein."
                    if complete
                    else "Compatibility predicate coverage is incomplete."
                ),
                required=True,
                evidence={
                    "proteins": protein_count,
                    "predicate_rows": predicate_count,
                    "coverage": (predicate_count / protein_count if protein_count else 0.0),
                },
                remediation=(
                    None
                    if complete
                    else "Recompute V2 predicates with compatibility materialization."
                ),
            )
        )
    return checks


def _embedding_checks(
    db_path: Path,
    sequence_protein_count: int,
    invalid_sequence_count: int,
) -> list[Capability]:
    h5_path = _discover_embedding_h5(db_path)
    if h5_path is None:
        missing = Capability(
            "embeddings",
            CapabilityState.unavailable,
            "No canonical protein embedding H5 was discovered.",
            remediation="Run Stage 06 if similarity or ELSA is required.",
        )
        index = Capability(
            "similarity_index",
            CapabilityState.unavailable,
            "Similarity index cannot exist without embeddings.",
            remediation="Run Stage 06, then build the persistent vector index.",
        )
        return [missing, index]

    inspection = inspect_vector_index(h5_path)
    if inspection.state == "failed" and inspection.count is None:
        embedding_state = CapabilityState.failed
        embedding_summary = "Embedding H5 failed structural validation."
    elif inspection.count != sequence_protein_count:
        embedding_state = CapabilityState.stale
        embedding_summary = (
            f"Embedding H5 has {inspection.count} rows for "
            f"{sequence_protein_count} sequence-bearing dataset proteins."
        )
    else:
        embedding_state = CapabilityState.available
        embedding_summary = (
            "Embedding row count matches the sequence-bearing dataset protein count."
        )
    embeddings = Capability(
        "embeddings",
        embedding_state,
        embedding_summary,
        evidence={
            "path": str(h5_path),
            "count": inspection.count,
            "sequence_protein_count": sequence_protein_count,
            "row_count_ratio": (
                inspection.count / sequence_protein_count
                if inspection.count is not None and sequence_protein_count
                else 0.0
            ),
            "dimension": inspection.dimension,
            "detail": inspection.detail if embedding_state == CapabilityState.failed else None,
        },
        remediation=(
            "Regenerate the embedding H5 from the current Stage 03 proteins."
            if embedding_state in {CapabilityState.failed, CapabilityState.stale}
            else None
        ),
    )
    index_state = CapabilityState(inspection.state)
    index_summary = inspection.summary
    if embedding_state == CapabilityState.stale and index_state == CapabilityState.available:
        index_state = CapabilityState.stale
        index_summary = "Vector sidecars match an embedding H5 that is stale for this dataset."
    index = Capability(
        "similarity_index",
        index_state,
        index_summary,
        evidence=inspection.to_dict(),
        remediation=(
            None
            if index_state == CapabilityState.available
            else (
                "Regenerate current embeddings before rebuilding the vector index."
                if embedding_state == CapabilityState.stale
                else (
                    "Repair invalid protein rows before building the vector index."
                    if invalid_sequence_count
                    else f"Run `sharur build-vector-index --embeddings {h5_path}`."
                )
            )
        ),
    )
    return [embeddings, index]


def _ledger_check(db_path: Path) -> Capability:
    ledger_path = db_path.parent / "sharur_ops.db"
    if not ledger_path.is_file():
        return Capability(
            "run_ledger",
            CapabilityState.unavailable,
            "Dataset-local operational ledger has not been created.",
            evidence={"path": str(ledger_path)},
            remediation="Start the dataset through current sharur-ingest.",
        )
    try:
        uri = f"{ledger_path.resolve().as_uri()}?mode=ro"
        conn = sqlite3.connect(uri, uri=True)
        try:
            tables = {
                row[0] for row in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")
            }
            run_count = (
                int(conn.execute("SELECT COUNT(*) FROM runs").fetchone()[0])
                if "runs" in tables
                else 0
            )
        finally:
            conn.close()
        if "runs" not in tables or "run_events" not in tables:
            return Capability(
                "run_ledger",
                CapabilityState.stale,
                "Operational SQLite exists but predates the unified run ledger.",
                evidence={"path": str(ledger_path), "tables": sorted(tables)},
                remediation="Open it once with current OpsStore or rerun sharur-ingest.",
            )
        return Capability(
            "run_ledger",
            CapabilityState.available,
            "Dataset-local run ledger is readable.",
            evidence={"path": str(ledger_path), "runs": run_count},
        )
    except Exception as exc:
        return Capability(
            "run_ledger",
            CapabilityState.failed,
            "Operational ledger could not be inspected.",
            evidence={"error": f"{type(exc).__name__}: {exc}"},
            remediation="Inspect or restore the SQLite operational ledger.",
        )


def _execution_checks() -> list[Capability]:
    from sharur.ingest.resources import resolve_resource_profile  # noqa: PLC0415

    local = resolve_resource_profile("local")
    checks = [
        Capability(
            "execution_local",
            CapabilityState.available,
            "Local CPU execution profile is available.",
            evidence={
                "cpu_count": os.cpu_count(),
                "profile": local.to_dict(),
            },
        )
    ]

    if sys.platform != "darwin":
        checks.append(
            Capability(
                "execution_mps",
                CapabilityState.unavailable,
                "Apple MPS execution is unavailable on this platform.",
            )
        )
    else:
        try:
            # Keep PyTorch out of the operator process. On macOS, importing it
            # before FAISS can load conflicting native OpenMP runtimes and abort
            # a later similarity query.
            process = subprocess.run(
                [
                    sys.executable,
                    "-c",
                    (
                        "import json, torch; "
                        "print(json.dumps({'available': bool("
                        "torch.backends.mps.is_built() and "
                        "torch.backends.mps.is_available()), "
                        "'version': torch.__version__}))"
                    ),
                ],
                check=True,
                capture_output=True,
                text=True,
                timeout=15,
            )
            probe = json.loads(process.stdout)
            available = bool(probe["available"])
            evidence: dict[str, Any] = {"torch_version": probe["version"]}
            if available:
                evidence["profile"] = resolve_resource_profile("mps").to_dict()
            checks.append(
                Capability(
                    "execution_mps",
                    (CapabilityState.available if available else CapabilityState.unavailable),
                    (
                        "Apple MPS profile is available with an exclusive process lock."
                        if available
                        else "PyTorch does not report a usable Apple MPS backend."
                    ),
                    evidence=evidence,
                    remediation=(
                        None
                        if available
                        else "Use --profile local or install an MPS-capable PyTorch build."
                    ),
                )
            )
        except Exception as exc:
            checks.append(
                Capability(
                    "execution_mps",
                    CapabilityState.failed,
                    "Apple MPS availability probe failed.",
                    evidence={"error": f"{type(exc).__name__}: {exc}"},
                )
            )

    sbatch = shutil.which("sbatch")
    checks.append(
        Capability(
            "execution_slurm",
            (CapabilityState.available if sbatch else CapabilityState.unavailable),
            (
                "SLURM submission profile is available."
                if sbatch
                else "SLURM submission tools are not on PATH."
            ),
            evidence={
                "sbatch": sbatch,
                "profile": resolve_resource_profile("slurm").to_dict(),
            },
            remediation=(
                None
                if sbatch
                else "Generate the bundle locally, then submit it on a SLURM login node."
            ),
        )
    )
    return checks


def _toolchain_check() -> Capability:
    from sharur import diagnostics  # noqa: PLC0415 - version probes are opt-in

    checks = diagnostics.run_all_checks()
    core_missing = [
        item.label for item in checks if item.core and item.status == diagnostics.MISSING
    ]
    state = CapabilityState.unavailable if core_missing else CapabilityState.available
    return Capability(
        "ingest_toolchain",
        state,
        (
            f"Core ingest components missing: {', '.join(core_missing)}."
            if core_missing
            else "Core ingest tools and reference databases are available."
        ),
        evidence={
            "checks": [
                {
                    "label": item.label,
                    "status": item.status,
                    "core": item.core,
                    "detail": item.detail,
                    "purpose": item.purpose,
                }
                for item in checks
            ]
        },
        remediation=(
            "Run `sharur doctor --strict` for installation details." if core_missing else None
        ),
    )


def build_capability_brief(
    db_path: str | Path,
    *,
    include_tools: bool = False,
) -> CapabilityBrief:
    """Build one non-mutating capability/preflight brief."""
    resolved = Path(db_path).expanduser().resolve()
    (
        capabilities,
        store,
        tables,
        protein_count,
        sequence_protein_count,
        invalid_sequence_count,
    ) = _database_checks(resolved)
    try:
        if store is not None and {"annotations", "proteins"} <= tables.keys():
            capabilities.extend(_annotation_checks(store, tables))
            capabilities.extend(_semantic_checks(store, tables, protein_count))
        else:
            for capability_id in (
                "annotation_sources",
                "structured_callers",
                "semantic_v2",
                "predicate_compatibility",
            ):
                capabilities.append(
                    Capability(
                        capability_id,
                        CapabilityState.unavailable,
                        "Unavailable because the core dataset is not usable.",
                        required=capability_id in {"semantic_v2", "predicate_compatibility"},
                    )
                )
    finally:
        if store is not None:
            store.close()
    capabilities.extend(_embedding_checks(resolved, sequence_protein_count, invalid_sequence_count))
    capabilities.append(_ledger_check(resolved))
    capabilities.extend(_execution_checks())
    if include_tools:
        capabilities.append(_toolchain_check())

    return CapabilityBrief(
        dataset_path=str(resolved),
        generated_at=datetime.now(timezone.utc).isoformat(),
        overall_state=_overall_state(capabilities),
        capabilities=tuple(capabilities),
    )


__all__ = [
    "Capability",
    "CapabilityBrief",
    "CapabilityState",
    "build_capability_brief",
]
