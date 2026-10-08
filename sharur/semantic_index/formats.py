"""Format identifiers and helpers shared by the readers, builders and generation layer (no native deps)."""

from __future__ import annotations

import hashlib
from typing import TYPE_CHECKING


if TYPE_CHECKING:
    from pathlib import Path

ACTIVE_EXPRESSION = "(term_kind != 'atom' OR relation != 'excludes')"
MEMBERSHIP_FORMAT = "sharur-isolated-membership-v1"
FORWARD_FORMAT = "sharur-isolated-rich-forward-v1"
SCOPE_FORMAT = "sharur-isolated-genome-scope-v1"
FIELDS = ("term_id", "term_kind", "facet", "relation", "source_db", "source_accession")
# Build-time controls for the numeric DuckDB comparison; never deployed or read.
CONTROL_FILES = frozenset({"numeric_pairs.u32", "numeric_pairs.parquet"})


def sha256_file(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(16 << 20), b""):
            h.update(block)
    return h.hexdigest()
