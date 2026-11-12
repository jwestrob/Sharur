"""Annotation-independent protein families (MMseqs2), an optional dataset layer.

- :mod:`sharur.protein_families.mmseqs`: clustering and search commands (standard
  library only, runnable as a script on a cluster node);
- :mod:`sharur.protein_families.store`: export, staging build, summaries,
  activation into a dataset database and incremental assignment;
- :mod:`sharur.protein_families.cli`: the ``sharur families`` commands.

Guide: ``docs/guides/protein-families.md``.
"""

from sharur.protein_families.mmseqs import (
    FAMILY_MIN_SEQ_ID,
    SUBFAMILY_MIN_SEQ_ID,
    LevelParams,
    default_levels,
)

__all__ = ["FAMILY_MIN_SEQ_ID", "SUBFAMILY_MIN_SEQ_ID", "LevelParams", "default_levels"]
