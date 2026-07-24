"""Executable, bounded verification records for scientific review claims."""

from __future__ import annotations

import hashlib
import json
import math
import re
import time
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal

import duckdb
from pydantic import BaseModel, ConfigDict, Field, model_validator

from sharur.dataset_seal import DEFAULT_SEAL_NAME, verify_dataset_seal
from sharur.ops.review_store import (
    contains_raw_biological_sequence,
    content_hash,
)


if TYPE_CHECKING:
    from sharur.ops.store import OpsStore


FORBIDDEN_SQL = re.compile(
    r"\b("
    r"attach|detach|copy|export|import|install|load|create|insert|update|"
    r"delete|drop|alter|call|pragma|set|vacuum|truncate|merge|read_csv|"
    r"read_csv_auto|read_parquet|read_json|read_json_auto|read_text|read_blob|"
    r"glob|csv_scan|parquet_scan|json_scan|sqlite_scan|sqlite_attach|"
    r"postgres_scan|mysql_scan|iceberg_scan|delta_scan|httpfs"
    r")\b",
    flags=re.IGNORECASE,
)
RAW_SEQUENCE_COLUMN = re.compile(
    r"\b(sequence|protein_sequence|amino_acid_sequence|nucleotide_sequence)\b",
    flags=re.IGNORECASE,
)
MAX_RESULT_BYTES = 256 * 1024


class VerificationSpec(BaseModel):
    """One declarative read-only query and its comparison rule."""

    model_config = ConfigDict(extra="forbid")

    sql: str = Field(min_length=1, max_length=262_144)
    parameters: list[Any] = Field(default_factory=list, max_length=1_000)
    result_shape: Literal["scalar", "row", "rows"] = "scalar"
    comparison: Literal["exact", "approx", "contains", "set_equal"] = "exact"
    absolute_tolerance: float = Field(default=0.0, ge=0.0)
    relative_tolerance: float = Field(default=0.0, ge=0.0)
    max_rows: int = Field(default=100, ge=1, le=1_000)

    @model_validator(mode="after")
    def validate_sql(self) -> VerificationSpec:
        sql = self.sql.strip()
        statements = [part.strip() for part in sql.split(";") if part.strip()]
        if len(statements) != 1:
            raise ValueError("Verification SQL must contain exactly one statement")
        first = statements[0].split(None, 1)[0].lower()
        if first not in {"select", "with"}:
            raise ValueError("Verification SQL must be a SELECT or WITH query")
        if FORBIDDEN_SQL.search(statements[0]):
            raise ValueError("Verification SQL contains a forbidden operation")
        if RAW_SEQUENCE_COLUMN.search(statements[0]):
            raise ValueError("Verification SQL cannot select raw sequence fields")
        return self


@dataclass(frozen=True)
class VerificationResult:
    status: Literal["pass", "fail", "error"]
    actual: Any
    executed_ts: float
    specification_hash: str
    query_hash: str
    row_count: int | None
    error: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "status": self.status,
            "actual": self.actual,
            "executed_ts": self.executed_ts,
            "specification_hash": self.specification_hash,
            "query_hash": self.query_hash,
            "row_count": self.row_count,
            "error": self.error,
        }


def _json_value(value: Any) -> Any:
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    return str(value)


def _shape_result(
    rows: list[tuple[Any, ...]],
    columns: list[str],
    shape: str,
) -> Any:
    if shape == "scalar":
        if len(rows) != 1 or len(columns) != 1:
            raise ValueError(
                "scalar verification requires exactly one row and one column"
            )
        return _json_value(rows[0][0])
    dictionaries = [
        {column: _json_value(value) for column, value in zip(columns, row, strict=True)}
        for row in rows
    ]
    if shape == "row":
        if len(dictionaries) != 1:
            raise ValueError("row verification requires exactly one row")
        return dictionaries[0]
    return dictionaries


def _compare(actual: Any, expected: Any, spec: VerificationSpec) -> bool:
    if spec.comparison == "exact":
        return actual == expected
    if spec.comparison == "approx":
        if isinstance(actual, bool) or isinstance(expected, bool):
            return False
        if not isinstance(actual, (int, float)) or not isinstance(
            expected, (int, float)
        ):
            return False
        return math.isclose(
            float(actual),
            float(expected),
            rel_tol=spec.relative_tolerance,
            abs_tol=spec.absolute_tolerance,
        )
    if spec.comparison == "contains":
        if isinstance(actual, dict) and isinstance(expected, dict):
            return all(actual.get(key) == value for key, value in expected.items())
        if isinstance(actual, list):
            expected_items = expected if isinstance(expected, list) else [expected]
            return all(item in actual for item in expected_items)
        if isinstance(actual, str) and isinstance(expected, str):
            return expected in actual
        return False
    if not isinstance(actual, list) or not isinstance(expected, list):
        return False
    return {
        json.dumps(item, sort_keys=True, separators=(",", ":"), default=str)
        for item in actual
    } == {
        json.dumps(item, sort_keys=True, separators=(",", ":"), default=str)
        for item in expected
    }


def _verify_sealed_dataset(
    db_path: Path,
    dataset_id: str,
    *,
    seal_path: str | Path | None,
) -> None:
    resolved_seal = (
        Path(seal_path).expanduser().resolve()
        if seal_path is not None
        else db_path.parent / DEFAULT_SEAL_NAME
    )
    result = verify_dataset_seal(resolved_seal, db_path=db_path)
    if not result.valid:
        changed = ", ".join(result.changed_sections) or "dataset identity"
        raise ValueError(f"Dataset seal verification failed; changed: {changed}")
    try:
        seal = json.loads(resolved_seal.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"Could not read dataset seal {resolved_seal}: {exc}") from exc
    if seal.get("dataset_id") != dataset_id:
        raise ValueError("Verification dataset_id differs from the dataset seal")


def run_duckdb_verification(
    db_path: str | Path,
    specification: VerificationSpec | dict[str, Any],
    expected: Any,
    *,
    dataset_id: str,
    seal_path: str | Path | None = None,
    verify_seal: bool = True,
    threads: int = 1,
) -> VerificationResult:
    """Execute one bounded read-only verification query."""

    started = time.time()
    spec = (
        specification
        if isinstance(specification, VerificationSpec)
        else VerificationSpec.model_validate(specification)
    )
    canonical_spec = spec.model_dump(mode="json")
    specification_hash = content_hash(canonical_spec)
    query_hash = hashlib.sha256(
        json.dumps(
            {"sql": spec.sql.strip(), "parameters": spec.parameters},
            sort_keys=True,
            separators=(",", ":"),
            default=str,
        ).encode()
    ).hexdigest()
    database = Path(db_path).expanduser().resolve()
    if threads < 1:
        raise ValueError("threads must be positive")
    if not database.is_file():
        raise FileNotFoundError(f"DuckDB file does not exist: {database}")
    try:
        if verify_seal:
            _verify_sealed_dataset(database, dataset_id, seal_path=seal_path)
        connection = duckdb.connect(
            str(database),
            read_only=True,
            config={"enable_external_access": "false"},
        )
        try:
            connection.execute("SET threads = ?", [threads])
            cursor = connection.execute(spec.sql, spec.parameters)
            columns = [str(item[0]) for item in cursor.description]
            if any(RAW_SEQUENCE_COLUMN.fullmatch(column) for column in columns):
                raise ValueError("Verification result exposes a raw sequence field")
            rows = cursor.fetchmany(spec.max_rows + 1)
            if len(rows) > spec.max_rows:
                raise ValueError(
                    f"Verification result exceeds max_rows={spec.max_rows}"
                )
        finally:
            connection.close()
        actual = _shape_result(rows, columns, spec.result_shape)
        if contains_raw_biological_sequence(actual):
            raise ValueError("Verification result resembles raw biological sequence")
        encoded = json.dumps(actual, separators=(",", ":"), default=str).encode()
        if len(encoded) > MAX_RESULT_BYTES:
            raise ValueError(
                f"Verification result exceeds {MAX_RESULT_BYTES} bytes"
            )
        status: Literal["pass", "fail", "error"] = (
            "pass" if _compare(actual, expected, spec) else "fail"
        )
        return VerificationResult(
            status=status,
            actual=actual,
            executed_ts=started,
            specification_hash=specification_hash,
            query_hash=query_hash,
            row_count=len(rows),
        )
    except Exception as exc:
        return VerificationResult(
            status="error",
            actual=None,
            executed_ts=started,
            specification_hash=specification_hash,
            query_hash=query_hash,
            row_count=None,
            error=str(exc)[:16_384],
        )


def execute_and_record_verification(
    store: OpsStore,
    *,
    review_id: str,
    claim_key: str,
    db_path: str | Path,
    specification: VerificationSpec | dict[str, Any],
    expected: Any,
    dataset_id: str,
    seal_path: str | Path | None = None,
    verify_seal: bool = True,
    threads: int = 1,
    code_commit: str | None = None,
    supersedes_verification_id: str | None = None,
) -> tuple[str, VerificationResult]:
    """Execute a query and append its exact result to the Ops review graph."""

    spec = (
        specification
        if isinstance(specification, VerificationSpec)
        else VerificationSpec.model_validate(specification)
    )
    spec_hash = content_hash(spec.model_dump(mode="json"))
    if supersedes_verification_id is None:
        existing = [
            row
            for row in store.list_review_verifications(review_id)
            if row["claim_key"] == claim_key
            and row["specification_hash"] == spec_hash
        ]
        if existing:
            latest = max(
                existing,
                key=lambda row: (
                    int(row["attempt"]),
                    float(row["created_ts"]),
                    str(row["id"]),
                ),
            )
            return str(latest["id"]), VerificationResult(
                status=latest["status"],
                actual=latest["actual"],
                executed_ts=float(
                    latest["executed_ts"] or latest["created_ts"]
                ),
                specification_hash=spec_hash,
                query_hash=hashlib.sha256(
                    json.dumps(
                        {
                            "sql": spec.sql.strip(),
                            "parameters": spec.parameters,
                        },
                        sort_keys=True,
                        separators=(",", ":"),
                        default=str,
                    ).encode()
                ).hexdigest(),
                row_count=None,
                error=latest["error"],
            )
    result = run_duckdb_verification(
        db_path,
        spec,
        expected,
        dataset_id=dataset_id,
        seal_path=seal_path,
        verify_seal=verify_seal,
        threads=threads,
    )
    idempotency_key = content_hash(
        {
            "review_id": review_id,
            "claim_key": claim_key,
            "specification": spec.model_dump(mode="json"),
            "expected": expected,
            "result": result.to_dict(),
            "supersedes_verification_id": supersedes_verification_id,
        }
    )
    verification_id = store.record_review_verification(
        review_id=review_id,
        claim_key=claim_key,
        engine="duckdb",
        specification=spec.model_dump(mode="json"),
        dataset_id=dataset_id,
        expected=expected,
        actual=result.actual,
        status=result.status,
        executed_ts=result.executed_ts,
        code_commit=code_commit,
        error=result.error,
        supersedes_verification_id=supersedes_verification_id,
        idempotency_key=f"verify:{idempotency_key}",
    )
    return verification_id, result


__all__ = [
    "VerificationResult",
    "VerificationSpec",
    "execute_and_record_verification",
    "run_duckdb_verification",
]
