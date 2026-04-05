"""
Canonical helpers for dataset-level scientific analysis records.

This module defines stable conventions for canonical findings records.
Compatibility handling for older findings shapes lives in
`sharur.core.analysis_record_compat`.
"""

from __future__ import annotations

import re
from pathlib import Path

CANONICAL_FINDING_SCHEMA_VERSION = "2.0"

FINDING_PHASE_ALIASES = {
    "explore": "exploration",
    "deepen": "exploration",
    "char": "characterization",
    "characterize": "characterization",
}

CANONICAL_FINDING_PHASES = (
    "survey",
    "exploration",
    "characterization",
    "defense",
    "metabolism",
)

CANONICAL_FINDING_FILE_PHASES = ("survey", "exploration")

_LEGACY_FINDING_ID_RE = re.compile(r"^[A-Z]+[0-9]{3,}$")


def canonicalize_finding_phase(phase: str | None) -> str | None:
    """Return a normalized phase slug for a finding."""
    if phase is None:
        return None

    normalized = str(phase).strip().lower().replace("-", "_").replace(" ", "_")
    if not normalized:
        return None
    return FINDING_PHASE_ALIASES.get(normalized, normalized)


def infer_finding_phase_from_path(path: str | Path | None) -> str | None:
    """Infer a finding phase from a path like survey/findings.jsonl."""
    if path is None:
        return None

    path_obj = Path(path)
    for part in path_obj.parts:
        phase = canonicalize_finding_phase(part)
        if phase in CANONICAL_FINDING_PHASES:
            return phase
    return None


def canonical_finding_id(phase: str | None, ordinal: int) -> str:
    """Generate a canonical `{phase}-{NNN}` finding ID."""
    phase_slug = canonicalize_finding_phase(phase) or "finding"
    return f"{phase_slug}-{ordinal:03d}"


def is_legacy_finding_id(value: str | None) -> bool:
    """Return True if the ID looks like a legacy top-level finding ID."""
    if not value:
        return False
    return bool(_LEGACY_FINDING_ID_RE.fullmatch(value))


__all__ = [
    "CANONICAL_FINDING_FILE_PHASES",
    "CANONICAL_FINDING_PHASES",
    "CANONICAL_FINDING_SCHEMA_VERSION",
    "canonical_finding_id",
    "canonicalize_finding_phase",
    "infer_finding_phase_from_path",
    "is_legacy_finding_id",
]
