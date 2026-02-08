"""
Base classes for Bennu operators.

BennuResult is the standard return type for all operator methods,
providing structured output with metadata for token budgeting and tracing.
"""

from __future__ import annotations

import hashlib
import time
from dataclasses import dataclass, field
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
        """Compute trace hash from operator + params."""
        content = f"{self.operator}:{sorted(self.params.items())}"
        self.trace_hash = hashlib.sha256(content.encode()).hexdigest()[:12]


@dataclass
class BennuResult:
    """
    Standard result from a Bennu operator.

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

    def __str__(self) -> str:
        """Return formatted data for display."""
        return self.data

    def __repr__(self) -> str:
        return f"BennuResult(rows={self.meta.rows}, truncated={self.meta.truncated})"


class OperatorContext:
    """
    Context for operator execution.

    Provides timing utilities and result construction helpers.
    """

    def __init__(self, operator_name: str, params: dict):
        self.operator_name = operator_name
        self.params = params
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
        dataset_version: str = "1.0",
        schema_version: str = "1.0",
    ) -> BennuResult:
        """Construct a BennuResult with metadata and trace."""
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
        return BennuResult(
            data=data,
            meta=meta,
            ref=ref,
            trace=trace,
            _raw=raw,
        )


__all__ = ["ResultMeta", "OperatorTrace", "BennuResult", "OperatorContext"]
