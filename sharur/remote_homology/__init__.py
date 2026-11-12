"""Remote-homology annotation: CDD per protein and HMM-HMM searches per protein family.

Two optional layers record *observed similarity* to characterized families, read as
"resembles":

- **CDD per protein** (``runner.py cdd``): RPS-BLAST of every unique protein sequence
  against CDD; hits at E <= 0.01 are stored in ``annotations`` with ``source = 'cdd'``.
- **Family tier** (``runner.py families``): for each selected protein family, an alignment
  of up to 250 unique members (representative first), optional HHblits enrichment against
  UniRef30, and HHsearch against Pfam and PDB profile databases. Hits are stored per family
  in ``family_annotations``.

Neither layer names a system, a pathway or a function: a row says that a family or a
protein resembles a profile. The SQL surface labels these rows ``evidence = 'resembles'``.

:mod:`sharur.remote_homology.hhr` and :mod:`sharur.remote_homology.runner` use only the
standard library, so the searches run on any host with Python 3, HH-suite, MAFFT and
BLAST+; :mod:`sharur.remote_homology.store` (DuckDB) prepares inputs and loads results.
"""

from __future__ import annotations

from sharur.remote_homology.hhr import (
    FAMILY_SOURCES,
    MSA_ENRICHED,
    MSA_FAMILY_ONLY,
    MSA_KINDS,
    RULE_SOURCES,
    SOURCE_PDB,
    SOURCE_PDB_HHBLITS,
    SOURCE_PFAM,
    SOURCE_UNIREF,
    STORE_MIN_PROBABILITY,
    UNKNOWN_FUNCTION,
    UNKNOWN_FUNCTION_PATTERN,
    is_unknown_function,
)


# Annotation sources whose rows are observed similarities ("resembles"), never annotations
# at curated thresholds. Predicate generation skips them, and the SQL surface labels them.
# sharur/sql/prelude.sql repeats this list; tests/test_remote_homology.py keeps both equal.
RESEMBLES_SOURCES: frozenset[str] = frozenset({"cdd"})

# The default "resembles" rule for family-tier hits (probability in percent).
RULE_MIN_PROBABILITY = 90.0
RULE_MIN_ALIGNED_COLS = 30
CDD_MAX_EVALUE = 0.01

__all__ = [
    "CDD_MAX_EVALUE", "FAMILY_SOURCES", "MSA_ENRICHED", "MSA_FAMILY_ONLY", "MSA_KINDS", "RESEMBLES_SOURCES",
    "RULE_MIN_ALIGNED_COLS", "RULE_MIN_PROBABILITY", "RULE_SOURCES", "SOURCE_PDB", "SOURCE_PDB_HHBLITS", "SOURCE_PFAM",
    "SOURCE_UNIREF", "STORE_MIN_PROBABILITY", "UNKNOWN_FUNCTION", "UNKNOWN_FUNCTION_PATTERN",
    "is_unknown_function",
]
