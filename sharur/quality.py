"""Genome completeness and contamination from CheckM, CheckM2 or GTDB tables.

Stage 02 (DFAST_QC) fills ``bins.completeness`` and ``bins.contamination``
during ingest. Datasets ingested without it can take the values from any
quality table that names genomes in one column and reports completeness and
contamination (percent) in others:

- CheckM2 ``quality_report.tsv``: ``Name``, ``Completeness``, ``Contamination``
- CheckM ``--tab_table`` output: ``Bin Id``, ``Completeness``, ``Contamination``
- GTDB metadata (``ar53_metadata.tsv``, ``bac120_metadata.tsv``): ``accession``
  (``RS_``/``GB_`` prefixes are dropped), ``checkm2_completeness`` (preferred)
  or ``checkm_completeness``, and the matching contamination column
"""

from __future__ import annotations

import csv
import gzip
import re
from pathlib import Path
from typing import Any

ID_COLUMNS = ("bin_id", "genome_id", "genome", "Name", "Bin Id", "Bin_Id", "user_genome", "accession")
COMPLETENESS = ("completeness", "Completeness", "checkm2_completeness", "checkm_completeness")
CONTAMINATION = ("contamination", "Contamination", "checkm2_contamination", "checkm_contamination")
_GTDB_PREFIX = re.compile(r"^(RS|GB)_")
_SUFFIX = re.compile(r"\.(fna|fa|fasta|fas)(\.gz)?$")


def _open(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def _pick(header: list[str], options: tuple[str, ...], preferred: str | None = None) -> str | None:
    if preferred:
        if preferred not in header:
            raise ValueError(f"column {preferred!r} not in table (columns: {', '.join(header[:12])}...)")
        return preferred
    # GTDB tables carry CheckM2 and CheckM columns; CheckM2 is the current estimate
    for name in sorted(options, key=lambda o: "checkm2" not in o):
        if name in header:
            return name
    return None


def normalize_id(value: str) -> str:
    return _SUFFIX.sub("", _GTDB_PREFIX.sub("", value.strip()))


def read_quality_table(path: str | Path, *, id_column: str | None = None, completeness_column: str | None = None,
                       contamination_column: str | None = None) -> tuple[dict[str, tuple[float, float | None]], dict[str, str]]:
    """``{genome: (completeness, contamination)}`` and the columns used."""
    path = Path(path)
    with _open(path) as handle:
        sample = handle.readline()
        delimiter = "\t" if "\t" in sample else ","
        handle.seek(0)
        reader = csv.DictReader(handle, delimiter=delimiter)
        header = reader.fieldnames or []
        id_col = _pick(header, ID_COLUMNS, id_column)
        comp_col = _pick(header, COMPLETENESS, completeness_column)
        cont_col = _pick(header, CONTAMINATION, contamination_column)
        if id_col is None or comp_col is None:
            raise ValueError(f"{path.name}: needs a genome column ({', '.join(ID_COLUMNS)}) and a completeness "
                             f"column ({', '.join(COMPLETENESS)}); found {', '.join(header[:12])}")
        values: dict[str, tuple[float, float | None]] = {}
        for row in reader:
            try:
                completeness = float(row[comp_col])
            except (TypeError, ValueError):
                continue
            contamination = None
            if cont_col:
                try:
                    contamination = float(row[cont_col])
                except (TypeError, ValueError):
                    contamination = None
            values[normalize_id(row[id_col])] = (completeness, contamination)
    return values, {"id": id_col, "completeness": comp_col, "contamination": cont_col or ""}


def import_quality(conn, values: dict[str, tuple[float, float | None]], *, overwrite: bool = False,
                   dry_run: bool = False) -> dict[str, Any]:
    """Fill ``bins.completeness``/``contamination`` for genomes in ``values``.

    Existing values are kept unless ``overwrite``. Returns match counts.
    """
    bins = conn.execute("SELECT bin_id, completeness FROM bins").fetchall()
    by_norm = {normalize_id(b): b for b, _ in bins}
    current = dict(bins)
    updates = []
    for key, (completeness, contamination) in values.items():
        bin_id = key if key in current else by_norm.get(key)
        if bin_id is None:
            continue
        if current[bin_id] is not None and not overwrite:
            continue
        updates.append((completeness, contamination, bin_id))
    matched = {u[2] for u in updates} | {b for b in current if current[b] is not None and
                                         (b in values or normalize_id(b) in values)}
    if updates and not dry_run:
        conn.executemany("UPDATE bins SET completeness = ?, contamination = ? WHERE bin_id = ?", updates)
    filled = [u[0] for u in updates]
    return {"bins": len(bins), "matched": len(matched), "updated": len(updates),
            "unmatched_bins": len(bins) - len(matched),
            "median_completeness": sorted(filled)[len(filled) // 2] if filled else None}
