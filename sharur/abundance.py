"""Per-sample read coverage and genome abundance (optional layer).

Coverage lives in a sidecar (``abundance.duckdb`` beside the dataset DuckDB) so
the canonical database and its seal are unchanged. Imports accept CoverM
``contig`` tables (``Contig``, then ``<sample> Mean``, ``<sample> Covered
Fraction``, ``<sample> Read Count``, ``<sample> Length`` columns, any subset,
any number of samples) or a long table with ``sample_id``, ``contig_id`` and any
of ``mean_depth``, ``covered_fraction``, ``read_count``, ``contig_length``.

Derived quantities:

- genome mean depth: contig depths weighted by contig length;
- genome relative abundance: the genome's share of reads mapped to dataset
  contigs in that sample (read counts; depth x length when counts are absent);
- coverage outliers: contigs whose depth departs from their genome's median
  contig depth by more than ``min_log2`` (log2 fold) in most samples, a common
  signature of binning errors, multicopy elements or strain variation.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import statistics
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import duckdb

DEFAULT_FILENAME = "abundance.duckdb"
SCHEMA_VERSION = 1
_COVERM_FIELDS = {"mean": "mean_depth", "covered fraction": "covered_fraction", "read count": "read_count",
                  "length": "contig_length"}

SCHEMA = """
CREATE TABLE IF NOT EXISTS abundance_meta (key VARCHAR PRIMARY KEY, value VARCHAR);
CREATE TABLE IF NOT EXISTS samples (
    sample_id VARCHAR PRIMARY KEY,
    metadata JSON,
    source VARCHAR,
    imported_at TIMESTAMP
);
CREATE TABLE IF NOT EXISTS contig_coverage (
    sample_id VARCHAR NOT NULL,
    contig_id VARCHAR NOT NULL,
    bin_id VARCHAR,
    contig_length BIGINT,
    mean_depth DOUBLE,
    covered_fraction DOUBLE,
    read_count BIGINT,
    PRIMARY KEY (sample_id, contig_id)
);
CREATE TABLE IF NOT EXISTS imports (
    path VARCHAR,
    sha256 VARCHAR,
    format VARCHAR,
    samples VARCHAR[],
    rows BIGINT,
    imported_at TIMESTAMP
);
"""


def default_path(db_path: str | Path) -> Path:
    return Path(db_path).expanduser().resolve().parent / DEFAULT_FILENAME


def connect(path: str | Path, *, read_only: bool = False) -> duckdb.DuckDBPyConnection:
    path = Path(path)
    if read_only:
        if not path.is_file():
            raise FileNotFoundError(f"No abundance sidecar at {path}; import coverage first.")
        return duckdb.connect(str(path), read_only=True)
    conn = duckdb.connect(str(path))
    conn.execute(SCHEMA)
    conn.execute("INSERT OR REPLACE INTO abundance_meta VALUES ('schema_version', ?)", [str(SCHEMA_VERSION)])
    return conn


# --------------------------------------------------------------------------- #
# Parsing
# --------------------------------------------------------------------------- #


def _open(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path, newline="")


def _number(value: str | None, kind: type) -> Any:
    if value is None or value.strip() in ("", "NA", "nan"):
        return None
    return kind(float(value)) if kind is int else kind(value)


def parse_coverm(path: Path) -> list[dict[str, Any]]:
    """Rows from a CoverM ``contig`` table (one or more samples)."""
    with _open(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if not header or header[0].strip().lower() not in ("contig", "contig_id"):
            raise ValueError(f"{path}: expected a CoverM contig table starting with 'Contig'")
        columns: list[tuple[int, str, str]] = []
        for index, name in enumerate(header[1:], start=1):
            for suffix, field in _COVERM_FIELDS.items():
                if name.lower().endswith(" " + suffix):
                    columns.append((index, name[: -len(suffix) - 1], field))
                    break
        if not columns:
            raise ValueError(f"{path}: no Mean, Covered Fraction or Read Count columns")
        rows: dict[tuple[str, str], dict[str, Any]] = {}
        for record in reader:
            if not record or not record[0]:
                continue
            for index, sample, field in columns:
                row = rows.setdefault((sample, record[0]), {"sample_id": sample, "contig_id": record[0]})
                row[field] = _number(record[index], int if field in ("read_count", "contig_length") else float)
    return list(rows.values())


def parse_long(path: Path) -> list[dict[str, Any]]:
    """Rows from a long table: sample_id, contig_id, mean_depth/covered_fraction/read_count."""
    with _open(path) as handle:
        sample = handle.read(4096)
        handle.seek(0)
        delimiter = "\t" if "\t" in sample else ","
        reader = csv.DictReader(handle, delimiter=delimiter)
        missing = {"sample_id", "contig_id"} - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing columns {sorted(missing)}")
        return [{"sample_id": r["sample_id"], "contig_id": r["contig_id"],
                 "mean_depth": _number(r.get("mean_depth"), float),
                 "covered_fraction": _number(r.get("covered_fraction"), float),
                 "read_count": _number(r.get("read_count"), int),
                 "contig_length": _number(r.get("contig_length"), int)} for r in reader]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def import_coverage(store, sidecar: Path, paths: list[Path], *, fmt: str = "coverm",
                    sample_metadata: Path | None = None, source: str | None = None,
                    allow_unknown_contigs: bool = False) -> dict[str, Any]:
    """Validate coverage rows against the core contigs and write them to the sidecar.

    Rows replace earlier ones for the same (sample, contig).
    """
    import pandas as pd  # noqa: PLC0415

    parse = {"coverm": parse_coverm, "long": parse_long}[fmt]
    rows = [row for path in paths for row in parse(path)]
    if not rows:
        raise ValueError("No coverage rows found")
    contigs = {cid: (bid, length) for cid, bid, length in
               store.execute("SELECT contig_id, bin_id, length FROM contigs")}
    unknown = sorted({r["contig_id"] for r in rows if r["contig_id"] not in contigs})
    if unknown and not allow_unknown_contigs:
        raise ValueError(f"{len(unknown):,} contigs are absent from the dataset (e.g. {unknown[:3]}); "
                         "pass --allow-unknown-contigs to keep only matching rows")
    rows = [r for r in rows if r["contig_id"] in contigs]
    frame = pd.DataFrame(rows).reindex(columns=["sample_id", "contig_id", "mean_depth", "covered_fraction",
                                                "read_count", "contig_length"])
    frame["bin_id"] = frame["contig_id"].map(lambda c: contigs[c][0])
    # Lengths reported by the mapper are the reference lengths; the core value is the fallback.
    frame["contig_length"] = frame["contig_length"].fillna(frame["contig_id"].map(lambda c: contigs[c][1]))
    frame = frame.drop_duplicates(["sample_id", "contig_id"], keep="last")
    samples = sorted(frame["sample_id"].unique())
    metadata = _read_metadata(sample_metadata) if sample_metadata else {}
    now = datetime.now(timezone.utc)
    conn = connect(sidecar)
    try:
        conn.execute("BEGIN TRANSACTION")
        conn.register("_cov", frame)
        conn.execute("""DELETE FROM contig_coverage USING _cov
                        WHERE contig_coverage.sample_id = _cov.sample_id
                          AND contig_coverage.contig_id = _cov.contig_id""")
        conn.execute("""INSERT INTO contig_coverage
                        SELECT sample_id, contig_id, bin_id, contig_length, mean_depth, covered_fraction, read_count
                        FROM _cov""")
        for sample in samples:
            conn.execute(
                """INSERT INTO samples VALUES (?, ?, ?, ?) ON CONFLICT (sample_id) DO UPDATE SET
                   metadata = COALESCE(excluded.metadata, samples.metadata), source = excluded.source,
                   imported_at = excluded.imported_at""",
                [sample, json.dumps(metadata[sample]) if sample in metadata else None, source or fmt, now])
        for path in paths:
            conn.execute("INSERT INTO imports VALUES (?, ?, ?, ?, ?, ?)",
                         [str(path), _sha256(path), fmt, samples, len(frame), now])
        conn.execute("COMMIT")
    except Exception:
        conn.execute("ROLLBACK")
        raise
    finally:
        conn.close()
    return {"sidecar": str(sidecar), "samples": samples, "rows": len(frame), "unknown_contigs": len(unknown)}


def _read_metadata(path: Path) -> dict[str, dict[str, str]]:
    with _open(path) as handle:
        sample = handle.read(4096)
        handle.seek(0)
        reader = csv.DictReader(handle, delimiter="\t" if "\t" in sample else ",")
        if "sample_id" not in (reader.fieldnames or []):
            raise ValueError(f"{path}: sample metadata needs a sample_id column")
        return {r["sample_id"]: {k: v for k, v in r.items() if k != "sample_id"} for r in reader}


# --------------------------------------------------------------------------- #
# Queries
# --------------------------------------------------------------------------- #


def samples(sidecar: Path) -> list[dict[str, Any]]:
    with connect(sidecar, read_only=True) as conn:
        return [{"sample_id": s, "metadata": json.loads(m) if m else {}, "contigs": n}
                for s, m, n in conn.execute(
                    """SELECT s.sample_id, s.metadata, COUNT(c.contig_id) FROM samples s
                       LEFT JOIN contig_coverage c USING (sample_id) GROUP BY ALL ORDER BY 1""").fetchall()]


def genome_abundance(sidecar: Path, *, sample_ids: list[str] | None = None,
                     bins: list[str] | None = None) -> list[dict[str, Any]]:
    """Per genome and sample: length-weighted depth, covered fraction, reads, relative abundance."""
    where, params = [], []
    if sample_ids:
        where.append(f"sample_id IN ({','.join('?' * len(sample_ids))})")
        params += sample_ids
    clause = ("WHERE " + " AND ".join(where)) if where else ""
    with connect(sidecar, read_only=True) as conn:
        rows = conn.execute(
            f"""
            WITH per_bin AS (
                SELECT sample_id, bin_id,
                       SUM(mean_depth * contig_length) / NULLIF(SUM(contig_length) FILTER (WHERE mean_depth IS NOT NULL), 0)
                           AS mean_depth,
                       SUM(covered_fraction * contig_length)
                           / NULLIF(SUM(contig_length) FILTER (WHERE covered_fraction IS NOT NULL), 0)
                           AS covered_fraction,
                       SUM(read_count) AS reads,
                       SUM(mean_depth * contig_length) AS depth_bases,
                       COUNT(*) AS contigs
                FROM contig_coverage {clause} GROUP BY sample_id, bin_id
            )
            SELECT sample_id, bin_id, mean_depth, covered_fraction, reads, contigs,
                   COALESCE(reads / NULLIF(SUM(reads) OVER (PARTITION BY sample_id), 0),
                            depth_bases / NULLIF(SUM(depth_bases) OVER (PARTITION BY sample_id), 0))
                       AS relative_abundance
            FROM per_bin ORDER BY sample_id, relative_abundance DESC NULLS LAST, bin_id
            """, params).fetchall()
    keep = set(bins) if bins else None
    return [{"sample_id": s, "bin_id": b, "mean_depth": d, "covered_fraction": f, "reads": r,
             "contigs": n, "relative_abundance": a}
            for s, b, d, f, r, n, a in rows if keep is None or b in keep]


def feature_abundance(store, sidecar: Path, *, predicate: str | None = None,
                      annotation: str | None = None, sample_ids: list[str] | None = None) -> dict[str, Any]:
    """Share of each sample's mapped reads in genomes that carry a predicate or annotation."""
    if (predicate is None) == (annotation is None):
        raise ValueError("Give exactly one of predicate or annotation")
    if predicate is not None:
        carriers = dict(store.execute(
            """SELECT p.bin_id, COUNT(DISTINCT pp.protein_id) FROM protein_predicates pp
               JOIN proteins p USING (protein_id) WHERE list_contains(pp.predicates, ?) GROUP BY 1""",
            [predicate]))
        feature = {"predicate": predicate}
    else:
        carriers = dict(store.execute(
            """SELECT p.bin_id, COUNT(DISTINCT a.protein_id) FROM annotations a JOIN proteins p USING (protein_id)
               WHERE a.accession = ? OR a.name = ? GROUP BY 1""", [annotation, annotation]))
        feature = {"annotation": annotation}
    genomes = genome_abundance(sidecar, sample_ids=sample_ids)
    by_sample: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in genomes:
        by_sample[row["sample_id"]].append(row)
    result = []
    for sample, rows in sorted(by_sample.items()):
        carrying = [r for r in rows if r["bin_id"] in carriers]
        result.append({
            "sample_id": sample,
            "share_in_carrier_genomes": sum(r["relative_abundance"] or 0 for r in carrying),
            "carrier_genomes_detected": sum(1 for r in carrying if (r["mean_depth"] or 0) > 0),
            "top_carriers": [{"bin_id": r["bin_id"], "relative_abundance": r["relative_abundance"],
                              "copies": carriers[r["bin_id"]]} for r in carrying[:5]],
        })
    return {**feature, "carrier_genomes": len(carriers), "samples": result}


def coverage_outliers(sidecar: Path, bin_id: str, *, min_log2: float = 1.0,
                      min_fraction: float = 0.5, min_bin_depth: float = 1.0) -> dict[str, Any]:
    """Contigs whose depth departs from their genome's median contig depth.

    Uses samples where the genome's median contig depth is at least
    ``min_bin_depth``; a contig is an outlier when |log2(depth / median)| is at
    least ``min_log2`` in at least ``min_fraction`` of those samples.
    """
    with connect(sidecar, read_only=True) as conn:
        rows = conn.execute(
            "SELECT sample_id, contig_id, contig_length, mean_depth FROM contig_coverage "
            "WHERE bin_id = ? AND mean_depth IS NOT NULL", [bin_id]).fetchall()
    per_sample: dict[str, dict[str, float]] = defaultdict(dict)
    lengths: dict[str, int] = {}
    for sample, contig, length, depth in rows:
        per_sample[sample][contig] = depth
        lengths[contig] = length
    used = {s: statistics.median(d.values()) for s, d in per_sample.items()
            if d and statistics.median(d.values()) >= min_bin_depth}
    ratios: dict[str, list[float]] = defaultdict(list)
    for sample, median in used.items():
        for contig, depth in per_sample[sample].items():
            ratios[contig].append(math.log2(max(depth, 1e-3) / median))
    outliers = []
    for contig, values in ratios.items():
        high = sum(v >= min_log2 for v in values)
        low = sum(v <= -min_log2 for v in values)
        if max(high, low) >= min_fraction * len(values):
            outliers.append({"contig_id": contig, "contig_length": lengths[contig],
                             "direction": "high" if high >= low else "low",
                             "median_log2_ratio": statistics.median(values), "samples": len(values)})
    outliers.sort(key=lambda o: -abs(o["median_log2_ratio"]))
    return {"bin_id": bin_id, "samples_used": len(used), "contigs": len(ratios),
            "min_log2": min_log2, "outliers": outliers}


# --------------------------------------------------------------------------- #
# Markdown
# --------------------------------------------------------------------------- #


def genome_abundance_markdown(rows: list[dict[str, Any]], top: int = 10) -> str:
    lines = ["# Genome abundance"]
    by_sample: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        by_sample[row["sample_id"]].append(row)
    for sample, items in by_sample.items():
        lines.append(f"\n## {sample} ({len(items)} genomes with coverage)")
        for r in items[:top]:
            share = f"{r['relative_abundance']:.2%}" if r["relative_abundance"] is not None else "n/a"
            depth = f"{r['mean_depth']:.1f}x" if r["mean_depth"] is not None else "n/a"
            cov = f", {r['covered_fraction']:.0%} covered" if r["covered_fraction"] is not None else ""
            lines.append(f"- {r['bin_id']}: {share} of mapped reads, {depth}{cov}")
    return "\n".join(lines)


def feature_abundance_markdown(result: dict[str, Any]) -> str:
    feature = result.get("predicate") or result.get("annotation")
    lines = [f"# Abundance of genomes carrying {feature}",
             f"{result['carrier_genomes']:,} genomes carry it."]
    for s in result["samples"]:
        top = ", ".join(f"{c['bin_id']} ({(c['relative_abundance'] or 0):.1%})" for c in s["top_carriers"][:3])
        lines.append(f"- {s['sample_id']}: {s['share_in_carrier_genomes']:.1%} of mapped reads in "
                     f"{s['carrier_genomes_detected']} detected carriers" + (f"; top: {top}" if top else ""))
    return "\n".join(lines)


def coverage_outliers_markdown(result: dict[str, Any]) -> str:
    lines = [f"# Coverage outliers: {result['bin_id']}",
             f"{len(result['outliers'])} of {result['contigs']} contigs depart from the genome's median depth by "
             f">= {result['min_log2']:g} log2 in most of {result['samples_used']} samples."]
    for o in result["outliers"][:25]:
        lines.append(f"- {o['contig_id']} ({o['contig_length']:,} bp): {o['direction']}, "
                     f"median log2 ratio {o['median_log2_ratio']:+.2f}")
    return "\n".join(lines)
