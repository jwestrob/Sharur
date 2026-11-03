"""KOfam thresholds, and which bitscore a KO annotation row stores.

KOfam's ``ko_list`` gives each KO a bitscore ``threshold`` and a ``score_type``:
``full`` compares the threshold with the full-sequence bitscore, ``domain`` with
the best single-domain bitscore (KofamScan's rule). KOs without a threshold carry
``-`` in both columns and are kept at E <= 1e-15. An Astra/Aksha KOFAM hits table
holds one row per reported domain, with the full-sequence bitscore (``bitscore``)
and that domain's bitscore (``dom_bitscore``) side by side.

Loaders store, per (protein, KO), the score the KO's ``score_type`` names
(:func:`kofam_scores`) and record which score they stored in
``annotation_source_versions`` (:func:`record_ko_score_basis`). Readers that
compare stored scores with ``ko_list`` thresholds read that record
(:func:`ko_score_bases`) and know a comparison is KOfam's own when
:func:`compares_as_kofam` holds. KO rows loaded before the record existed carry
no basis: Stage 07 stored the bitscore of one domain, other loaders the
full-sequence bitscore.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING, Any


if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping

    import pandas as pd


KO_SOURCES = ("kegg", "kofam")

#: Full-sequence bitscore for score_type ``full`` (and KOs without a threshold), best-domain bitscore for ``domain``.
SCORE_BASIS_KOFAM = "kofam_score_type"
#: Full-sequence bitscore for every KO (the loader had no ko_list).
SCORE_BASIS_FULL = "full_sequence"
SCORE_RULES = {
    SCORE_BASIS_KOFAM: "score: KOfam score_type per KO (full-sequence bitscore for 'full' and threshold-less KOs, "
                       "best-domain bitscore for 'domain'); one row per (protein, KO), coordinates of the best domain",
    SCORE_BASIS_FULL: "score: full-sequence bitscore for every KO (no ko_list at load time); one row per (protein, KO), "
                      "coordinates of the best domain",
}


@dataclass(frozen=True)
class KoThreshold:
    """One ko_list entry: the bitscore threshold (None: KOfam has none) and the score it applies to."""

    threshold: float | None
    score_type: str | None  # 'full', 'domain', or None for '-'


def read_ko_list(path: str | Path, kos: Iterable[str] | None = None) -> dict[str, KoThreshold]:
    """``ko_list`` entries for ``kos`` (every KO when None)."""
    wanted = set(kos) if kos is not None else None
    entries: dict[str, KoThreshold] = {}
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        ko_idx, thr_idx, type_idx = header.index("knum"), header.index("threshold"), header.index("score_type")
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(ko_idx, thr_idx, type_idx) or (wanted is not None and parts[ko_idx] not in wanted):
                continue
            value, score_type = parts[thr_idx].strip(), parts[type_idx].strip()
            entries[parts[ko_idx]] = KoThreshold(
                None if value in ("", "-") else float(value),
                score_type if score_type in ("full", "domain") else None,
            )
    return entries


def kofam_scores(
    df: pd.DataFrame,
    score_types: Mapping[str, str | None] | None,
    *,
    full: str = "bitscore",
    domain: str = "dom_bitscore",
    protein: str = "protein_id",
    ko: str = "accession",
) -> tuple[pd.DataFrame, str]:
    """One row per (protein, KO) carrying the score KOfam compares with the KO's threshold.

    ``df`` holds one row per reported domain, with the full-sequence bitscore in
    ``full`` and the domain's bitscore in ``domain``. The row kept for each
    (protein, KO) is its best domain (its coordinates and E-values stay); its
    ``score`` becomes the full-sequence bitscore, or for a ``domain`` KO the best
    domain bitscore. ``score_types`` maps KO to score_type; None (no ko_list)
    stores the full-sequence bitscore for every KO. Returns the frame and the
    score basis to record.
    """
    if df.empty:
        out = df.copy()
        out["score"] = out.get(full)
        return out, SCORE_BASIS_KOFAM if score_types is not None else SCORE_BASIS_FULL
    best = (df.sort_values([protein, ko, domain], ascending=[True, True, False], na_position="last")
              .drop_duplicates(subset=[protein, ko], keep="first")
              .copy())
    if score_types is None:
        best["score"] = best[full]
        return best, SCORE_BASIS_FULL
    is_domain = best[ko].map(lambda k: score_types.get(k) == "domain").astype(bool)
    best["score"] = best[full].where(~is_domain, best[domain])
    return best, SCORE_BASIS_KOFAM


def compares_as_kofam(basis: str | None, score_type: str | None) -> bool:
    """Whether a score stored under ``basis`` is the one KOfam compares for a KO of ``score_type``."""
    if basis == SCORE_BASIS_KOFAM:
        return True
    if basis == SCORE_BASIS_FULL:
        return score_type != "domain"
    return False


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def record_ko_score_basis(
    conn,
    source: str,
    basis: str | None,
    *,
    ko_list: str | Path | None = None,
    hits_paths: Iterable[str | Path] = (),
    raw_rows: int | None = None,
    loaded_rows: int | None = None,
) -> None:
    """Record which score ``source``'s KO rows store; ``basis`` None removes the record (basis unknown)."""
    from sharur.defensefinder.selection import (  # noqa: PLC0415
        _SOURCE_VERSIONS_SQL,
        SOURCE_VERSIONS_TABLE,
    )

    conn.execute(_SOURCE_VERSIONS_SQL)
    if basis is None:
        conn.execute(f"DELETE FROM {SOURCE_VERSIONS_TABLE} WHERE source = ?", [source])
        return
    ko_list = Path(ko_list) if ko_list is not None else None
    parameters: dict[str, Any] = {"score_basis": basis}
    if ko_list is not None and ko_list.is_file():
        parameters.update(ko_list=str(ko_list), ko_list_sha256=_sha256(ko_list))
    paths = [str(p) for p in hits_paths]
    conn.execute(
        f"INSERT OR REPLACE INTO {SOURCE_VERSIONS_TABLE} (source, model_set, selection_rule, parameters, hits_path, "
        "raw_rows, loaded_rows, loaded_at) VALUES (?, 'KOfam', ?, ?, ?, ?, ?, ?)",
        [source, SCORE_RULES[basis], json.dumps(parameters, sort_keys=True), "; ".join(paths) or None,
         raw_rows, loaded_rows, datetime.now(timezone.utc).replace(tzinfo=None)],
    )


def ko_score_bases(conn) -> dict[str, str]:
    """Lower-cased KO source -> recorded score basis, for the KO sources that have a record."""
    from sharur.defensefinder.selection import SOURCE_VERSIONS_TABLE  # noqa: PLC0415

    exists = conn.execute(
        "SELECT COUNT(*) FROM information_schema.tables WHERE table_name = ? AND table_schema = 'main'",
        [SOURCE_VERSIONS_TABLE],
    ).fetchone()[0]
    if not exists:
        return {}
    bases: dict[str, str] = {}
    sources = ", ".join("?" * len(KO_SOURCES))
    for source, parameters in conn.execute(
        f"SELECT LOWER(source), parameters FROM {SOURCE_VERSIONS_TABLE} WHERE LOWER(source) IN ({sources})",
        list(KO_SOURCES),
    ).fetchall():
        try:
            basis = json.loads(parameters or "{}").get("score_basis")
        except (TypeError, ValueError):
            basis = None
        if basis in SCORE_RULES:
            bases[source] = basis
    return bases


__all__ = [
    "KO_SOURCES",
    "SCORE_BASIS_FULL",
    "SCORE_BASIS_KOFAM",
    "KoThreshold",
    "compares_as_kofam",
    "ko_score_bases",
    "kofam_scores",
    "read_ko_list",
    "record_ko_score_basis",
]
