"""
Base classes for Sharur operators.

SharurResult is the standard return type for all operator methods,
providing structured output with metadata for token budgeting and tracing.
"""

from __future__ import annotations

import hashlib
import json
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Optional


@dataclass
class ResultMeta:
    """Metadata about an operator result."""

    rows: int
    total_rows: Optional[int]
    bytes: int
    time_ms: int
    truncated: bool
    index_used: Optional[str] = None


@dataclass
class OperatorTrace:
    """Trace information for reproducibility and debugging."""

    operator: str
    params: dict
    dataset_version: str
    schema_version: str
    trace_hash: str = field(init=False)

    def __post_init__(self):
        """Compute a trace hash that is scoped to the dataset instance."""
        content = json.dumps(
            {
                "operator": self.operator,
                "params": self.params,
                "dataset_version": self.dataset_version,
                "schema_version": self.schema_version,
            },
            sort_keys=True,
            default=str,
            separators=(",", ":"),
        )
        self.trace_hash = hashlib.sha256(content.encode()).hexdigest()[:12]


@dataclass
class SharurResult:
    """
    Standard result from a Sharur operator.

    Attributes:
        data: Formatted output string (for display)
        meta: Result metadata (row counts, timing, truncation)
        ref: Reference ID for expand() pagination
        trace: Operator trace for reproducibility
        _raw: Raw data for programmatic access
    """

    data: str
    meta: ResultMeta
    ref: Optional[str] = None
    trace: Optional[OperatorTrace] = None
    _raw: Any = None
    status: str = "ok"
    warnings: tuple[str, ...] = ()

    def __str__(self) -> str:
        """Return formatted data for display."""
        return self.data

    def __repr__(self) -> str:
        return (
            f"SharurResult(status={self.status!r}, rows={self.meta.rows}, "
            f"truncated={self.meta.truncated})"
        )

    @property
    def raw(self) -> Any:
        """Public programmatic payload (``_raw`` remains for compatibility)."""
        return self._raw

    @property
    def records(self) -> list[Any]:
        """Return the payload as a record list for tabular consumers."""
        if self._raw is None:
            return []
        if isinstance(self._raw, list):
            return self._raw
        return [self._raw]

    def to_json(self, *, indent: Optional[int] = None) -> str:
        """Serialize the public result contract as JSON."""
        payload = {
            "status": self.status,
            "warnings": list(self.warnings),
            "meta": asdict(self.meta),
            "ref": self.ref,
            "trace": asdict(self.trace) if self.trace is not None else None,
            "raw": self._raw,
        }
        return json.dumps(payload, default=str, indent=indent)

    def to_dataframe(self):
        """Return records as a pandas DataFrame."""
        import pandas as pd

        return pd.DataFrame(self.records)


class OperatorContext:
    """
    Context for operator execution.

    Provides timing utilities and result construction helpers.
    """

    def __init__(self, operator_name: str, params: dict, *, store: Any = None):
        self.operator_name = operator_name
        self.params = params
        self.store = store
        self._start_time: Optional[float] = None

    def __enter__(self) -> "OperatorContext":
        self._start_time = time.perf_counter()
        return self

    def __exit__(self, *args):
        pass

    @property
    def elapsed_ms(self) -> int:
        """Get elapsed time in milliseconds."""
        if self._start_time is None:
            return 0
        return int((time.perf_counter() - self._start_time) * 1000)

    def make_result(
        self,
        data: str,
        rows: int,
        total_rows: Optional[int] = None,
        truncated: bool = False,
        ref: Optional[str] = None,
        raw: Any = None,
        index_used: Optional[str] = None,
        dataset_version: Optional[str] = None,
        schema_version: Optional[str] = None,
        status: Optional[str] = None,
        warnings: Optional[list[str] | tuple[str, ...]] = None,
    ) -> SharurResult:
        """Construct a SharurResult with metadata and trace."""
        detected_dataset_version, detected_schema_version = self._trace_versions()
        dataset_version = dataset_version or detected_dataset_version
        schema_version = schema_version or detected_schema_version
        meta = ResultMeta(
            rows=rows,
            total_rows=total_rows,
            bytes=len(data.encode("utf-8")),
            time_ms=self.elapsed_ms,
            truncated=truncated,
            index_used=index_used,
        )
        trace = OperatorTrace(
            operator=self.operator_name,
            params=self.params,
            dataset_version=dataset_version,
            schema_version=schema_version,
        )
        return SharurResult(
            data=data,
            meta=meta,
            ref=ref,
            trace=trace,
            _raw=raw,
            status=status or ("empty" if rows == 0 else "ok"),
            warnings=tuple(warnings or ()),
        )

    def _trace_versions(self) -> tuple[str, str]:
        """Return honest dataset-instance and live schema identifiers."""
        if self.store is None:
            return "not_applicable", "not_applicable"

        db_path = getattr(self.store, "db_path", None)
        if db_path is None:
            dataset_version = f"memory:{id(self.store)}"
        else:
            path = Path(db_path).expanduser().resolve()
            try:
                stat = path.stat()
                dataset_version = (
                    f"{path}@size={stat.st_size}:mtime_ns={stat.st_mtime_ns}"
                )
            except OSError:
                dataset_version = str(path)

        try:
            rows = self.store.execute("SELECT MAX(version) FROM schema_version")
            value = rows[0][0] if rows and rows[0][0] is not None else None
            schema_version = str(value) if value is not None else "unknown"
        except Exception:
            schema_version = "unknown"

        return dataset_version, schema_version


__all__ = ["ResultMeta", "OperatorTrace", "SharurResult", "OperatorContext"]
