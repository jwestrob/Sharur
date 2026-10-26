"""Per-sample gene expression from metatranscriptome read counts (optional layer).

Counts live in a sidecar (``expression.duckdb`` beside the dataset DuckDB) so the
canonical database and its seal are unchanged. Imports accept

- a long table (TSV/CSV/Parquet) with ``sample_id``, ``protein_id`` and either
  ``sense`` (plus optional ``antisense``) or ``count_s1`` and ``count_s2``
  (featureCounts ``-s 1`` / ``-s 2`` counts; sense is the orientation that most
  genes with >= 10 fragments favor, or the one holding more fragments when fewer
  than 20 genes qualify, unless the sample metadata gives
  ``sense_orientation`` = ``s1`` or ``s2``). Strand specificity is the fraction
  of those genes with more sense than antisense fragments (about 0.5 for an
  unstranded library);
- featureCounts tables (``Geneid ... Length <sample>...``), all with one
  orientation given at import (``sense`` or ``antisense``).

Counts are fragments assigned to a protein's CDS. Derived quantities:

- genome transcript share: a genome's share of the sample's sense fragments on
  dataset genes;
- TPM: fragments per kb of CDS, scaled to 10^6 over all dataset genes in the
  sample (``tpm``) or over one genome's genes (``genome_tpm``);
- within-genome percentile: a gene's rank by fragments per kb among all CDS of
  its genome (genes without reads included), 1.0 = the most expressed; genes
  without reads get 0.

Gene flags (``import_gene_flags``) mark CDS whose counts are not mRNA, such as
ORFs called inside rRNA, tRNA, tmRNA, RNase P, SRP or 6S RNA genes; reads from
those abundant RNAs land on them. Queries leave flagged genes out of every
numerator, denominator and ranking unless ``include_flagged=True``.
"""

from __future__ import annotations

import contextlib
import csv
import gzip
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import duckdb


DEFAULT_FILENAME = "expression.duckdb"
SCHEMA_VERSION = 1
# Orientation vote and strand specificity use genes with at least this many fragments ...
VOTE_MIN_FRAGMENTS = 10
# ... when at least this many such genes exist; otherwise fragment totals decide.
VOTE_MIN_GENES = 20

SCHEMA = """
CREATE TABLE IF NOT EXISTS expression_meta (key VARCHAR PRIMARY KEY, value VARCHAR);
CREATE TABLE IF NOT EXISTS rna_samples (
    sample_id VARCHAR PRIMARY KEY,
    metadata JSON,
    sense_orientation VARCHAR,
    strand_specificity DOUBLE,
    sense_fragments BIGINT,
    antisense_fragments BIGINT,
    source VARCHAR,
    imported_at TIMESTAMP
);
CREATE TABLE IF NOT EXISTS gene_counts (
    sample_id VARCHAR NOT NULL,
    protein_id VARCHAR NOT NULL,
    bin_id VARCHAR,
    gene_length INTEGER,
    sense BIGINT,
    antisense BIGINT,
    PRIMARY KEY (sample_id, protein_id)
);
CREATE TABLE IF NOT EXISTS gene_flags (
    protein_id VARCHAR NOT NULL,
    flag VARCHAR NOT NULL,
    detail VARCHAR,
    PRIMARY KEY (protein_id, flag)
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
            raise FileNotFoundError(f"No expression sidecar at {path}; import counts first.")
        return duckdb.connect(str(path), read_only=True)
    conn = duckdb.connect(str(path))
    conn.execute(SCHEMA)
    conn.execute("INSERT OR REPLACE INTO expression_meta VALUES ('schema_version', ?)", [str(SCHEMA_VERSION)])
    return conn


# --------------------------------------------------------------------------- #
# Parsing
# --------------------------------------------------------------------------- #


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _open(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path, newline="")


def _read_metadata(path: Path) -> dict[str, dict[str, str]]:
    with _open(path) as handle:
        head = handle.read(4096)
        handle.seek(0)
        reader = csv.DictReader(handle, delimiter="\t" if "\t" in head else ",")
        if "sample_id" not in (reader.fieldnames or []):
            raise ValueError(f"{path}: sample metadata needs a sample_id column")
        return {r["sample_id"]: {k: v for k, v in r.items() if k != "sample_id"} for r in reader}


def _relation(conn: duckdb.DuckDBPyConnection, path: Path) -> str:
    """SQL source for a long table file."""
    name = path.name.lower()
    if name.endswith(".parquet"):
        return f"read_parquet('{path}')"
    with _open(path) as handle:
        head = handle.read(4096)
    return f"read_csv('{path}', delim='{chr(9) if chr(9) in head else ','}', header=true)"


def _long_frame(conn: duckdb.DuckDBPyConnection, paths: list[Path]):
    parts = []
    for path in paths:
        cols = {c.lower() for c in conn.execute(f"DESCRIBE SELECT * FROM {_relation(conn, path)}").df()["column_name"]}
        missing = {"sample_id", "protein_id"} - cols
        if missing:
            raise ValueError(f"{path}: missing columns {sorted(missing)}")
        if {"count_s1", "count_s2"} <= cols:
            parts.append(f"SELECT sample_id, protein_id, count_s1::BIGINT AS s1, count_s2::BIGINT AS s2, "
                         f"NULL::BIGINT AS sense, NULL::BIGINT AS antisense FROM {_relation(conn, path)}")
        elif "sense" in cols:
            anti = "antisense::BIGINT" if "antisense" in cols else "NULL::BIGINT"
            parts.append(f"SELECT sample_id, protein_id, NULL::BIGINT AS s1, NULL::BIGINT AS s2, sense::BIGINT AS sense, "
                         f"{anti} AS antisense FROM {_relation(conn, path)}")
        else:
            raise ValueError(f"{path}: needs sense[/antisense] or count_s1/count_s2 columns")
    return " UNION ALL ".join(parts)


def _featurecounts_frame(conn: duckdb.DuckDBPyConnection, paths: list[Path], orientation: str):
    """featureCounts tables -> long rows with the given orientation."""
    import pandas as pd  # noqa: PLC0415

    frames = []
    for path in paths:
        with _open(path) as handle:
            lines = [ln for ln in handle if not ln.startswith("#")]
        reader = csv.reader(lines, delimiter="\t")
        header = next(reader)
        if header[:6] != ["Geneid", "Chr", "Start", "End", "Strand", "Length"]:
            raise ValueError(f"{path}: not a featureCounts table (expected Geneid ... Length header)")
        samples = [Path(h).name.removesuffix(".bam") for h in header[6:]]
        rows = [(s, rec[0], int(rec[6 + i])) for rec in reader for i, s in enumerate(samples) if int(rec[6 + i])]
        frames.append(pd.DataFrame(rows, columns=["sample_id", "protein_id", "n"]))
    frame = pd.concat(frames, ignore_index=True)
    frame["sense"] = frame["n"] if orientation == "sense" else None
    frame["antisense"] = frame["n"] if orientation == "antisense" else None
    frame["s1"] = None
    frame["s2"] = None
    conn.register("_fc", frame[["sample_id", "protein_id", "s1", "s2", "sense", "antisense"]])
    return "SELECT sample_id, protein_id, s1::BIGINT s1, s2::BIGINT s2, sense::BIGINT sense, antisense::BIGINT antisense FROM _fc"


# --------------------------------------------------------------------------- #
# Import
# --------------------------------------------------------------------------- #


def import_counts(store, sidecar: Path, paths: list[Path], *, fmt: str = "long", orientation: str = "sense",
                  sample_metadata: Path | None = None, source: str | None = None,
                  allow_unknown_proteins: bool = False) -> dict[str, Any]:
    """Validate gene counts against the core proteins and write them to the sidecar.

    Rows replace earlier ones for the same (sample, protein). Rows with no reads are dropped.
    """
    import pandas as pd  # noqa: PLC0415

    if fmt not in ("long", "featurecounts"):
        raise ValueError("fmt must be long or featurecounts")
    if orientation not in ("sense", "antisense"):
        raise ValueError("orientation must be sense or antisense")
    metadata = _read_metadata(sample_metadata) if sample_metadata else {}
    genes = pd.DataFrame(store.execute(
        "SELECT protein_id, bin_id, (end_coord - start + 1)::INTEGER AS gene_length FROM proteins"),
        columns=["protein_id", "bin_id", "gene_length"])
    now = datetime.now(timezone.utc)
    conn = connect(sidecar)
    try:
        conn.register("_genes", genes)
        body = _long_frame(conn, paths) if fmt == "long" else _featurecounts_frame(conn, paths, orientation)
        conn.execute(f"CREATE TEMP TABLE _raw AS SELECT sample_id, protein_id, SUM(s1) s1, SUM(s2) s2, "
                     f"SUM(sense) sense, SUM(antisense) antisense FROM ({body}) GROUP BY 1, 2")
        unknown = conn.execute("SELECT COUNT(DISTINCT protein_id) FROM _raw r ANTI JOIN _genes g USING (protein_id)"
                               ).fetchone()[0]
        if unknown and not allow_unknown_proteins:
            example = [r[0] for r in conn.execute(
                "SELECT DISTINCT protein_id FROM _raw r ANTI JOIN _genes g USING (protein_id) LIMIT 3").fetchall()]
            raise ValueError(f"{unknown:,} proteins are absent from the dataset (e.g. {example}); "
                             "pass --allow-unknown-proteins to keep only matching rows")
        # orientation per sample for s1/s2 input: metadata, else a per-gene majority vote (genes with
        # >= VOTE_MIN_FRAGMENTS), else the orientation holding more fragments. The vote resists a few
        # loci with huge one-strand counts (ORFs inside abundant non-coding RNAs) that can outweigh
        # every mRNA in fragment totals.
        totals = conn.execute(f"""SELECT sample_id, SUM(s1), SUM(s2),
                                         COUNT(*) FILTER (WHERE s1 + s2 >= {VOTE_MIN_FRAGMENTS} AND s1 > s2),
                                         COUNT(*) FILTER (WHERE s1 + s2 >= {VOTE_MIN_FRAGMENTS} AND s2 > s1)
                                  FROM _raw GROUP BY 1 ORDER BY 1""").fetchall()
        sense_of: dict[str, str | None] = {}
        for sample, t1, t2, v1, v2 in totals:
            given = (metadata.get(sample) or {}).get("sense_orientation", "").strip()
            if t1 is None and t2 is None:
                sense_of[sample] = None
            elif given in ("s1", "s2"):
                sense_of[sample] = given
            elif v1 + v2 >= VOTE_MIN_GENES and v1 != v2:
                sense_of[sample] = "s1" if v1 > v2 else "s2"
            else:
                sense_of[sample] = "s1" if (t1 or 0) >= (t2 or 0) else "s2"
        conn.register("_orient", pd.DataFrame(list(sense_of.items()), columns=["sample_id", "o"]))
        conn.execute("""CREATE TEMP TABLE _new AS
            SELECT r.sample_id, r.protein_id, g.bin_id, g.gene_length,
                   CASE o.o WHEN 's1' THEN r.s1 WHEN 's2' THEN r.s2 ELSE r.sense END AS sense,
                   CASE o.o WHEN 's1' THEN r.s2 WHEN 's2' THEN r.s1 ELSE r.antisense END AS antisense
            FROM _raw r JOIN _genes g USING (protein_id) JOIN _orient o USING (sample_id)
            WHERE COALESCE(CASE o.o WHEN 's1' THEN r.s1 WHEN 's2' THEN r.s2 ELSE r.sense END, 0) > 0
               OR COALESCE(CASE o.o WHEN 's1' THEN r.s2 WHEN 's2' THEN r.s1 ELSE r.antisense END, 0) > 0""")
        samples = [r[0] for r in conn.execute("SELECT DISTINCT sample_id FROM _new ORDER BY 1").fetchall()]
        rows = conn.execute("SELECT COUNT(*) FROM _new").fetchone()[0]
        conn.execute("BEGIN TRANSACTION")
        conn.execute("""DELETE FROM gene_counts USING _new
                        WHERE gene_counts.sample_id = _new.sample_id AND gene_counts.protein_id = _new.protein_id""")
        conn.execute("INSERT INTO gene_counts SELECT sample_id, protein_id, bin_id, gene_length, sense, antisense FROM _new")
        for sample in samples:
            s, a, n_vote, n_sense = conn.execute(
                f"""SELECT SUM(sense), SUM(antisense),
                           COUNT(*) FILTER (WHERE sense + antisense >= {VOTE_MIN_FRAGMENTS}),
                           COUNT(*) FILTER (WHERE sense + antisense >= {VOTE_MIN_FRAGMENTS} AND sense > antisense)
                    FROM gene_counts WHERE sample_id = ?""", [sample]).fetchone()
            if a is None or s is None:
                spec = None
            elif n_vote >= VOTE_MIN_GENES:
                spec = n_sense / n_vote
            else:
                spec = (s / (s + a)) if (s + a) else None
            conn.execute(
                """INSERT INTO rna_samples VALUES (?, ?, ?, ?, ?, ?, ?, ?) ON CONFLICT (sample_id) DO UPDATE SET
                   metadata = COALESCE(excluded.metadata, rna_samples.metadata),
                   sense_orientation = excluded.sense_orientation, strand_specificity = excluded.strand_specificity,
                   sense_fragments = excluded.sense_fragments, antisense_fragments = excluded.antisense_fragments,
                   source = excluded.source, imported_at = excluded.imported_at""",
                [sample, json.dumps(metadata[sample]) if sample in metadata else None,
                 sense_of.get(sample) or orientation, spec, s, a, source or fmt, now])
        for path in paths:
            conn.execute("INSERT INTO imports VALUES (?, ?, ?, ?, ?, ?)",
                         [str(path), _sha256(path), fmt, samples, rows, now])
        conn.execute("COMMIT")
    except Exception:
        with contextlib.suppress(duckdb.Error):
            conn.execute("ROLLBACK")
        raise
    finally:
        conn.close()
    return {"sidecar": str(sidecar), "samples": samples, "rows": rows, "unknown_proteins": unknown,
            "sense_orientation": sense_of}


def import_gene_flags(sidecar: Path, path: Path, *, flag: str = "structural_rna_overlap",
                      min_overlap_bp: int = 60) -> dict[str, Any]:
    """Flag CDS from an overlap table (protein_id, overlap_bp and optional rna_name/rfam_acc columns).

    Rows with overlap_bp >= min_overlap_bp flag their protein; earlier flags of the same name are replaced.
    """
    conn = connect(sidecar)
    try:
        cols = {c.lower() for c in conn.execute(f"DESCRIBE SELECT * FROM {_relation(conn, path)}").df()["column_name"]}
        if not {"protein_id", "overlap_bp"} <= cols:
            raise ValueError(f"{path}: needs protein_id and overlap_bp columns")
        label = "rna_name" if "rna_name" in cols else ("rfam_acc" if "rfam_acc" in cols else "NULL")
        conn.execute("BEGIN TRANSACTION")
        conn.execute("DELETE FROM gene_flags WHERE flag = ?", [flag])
        conn.execute(f"""INSERT INTO gene_flags
            SELECT protein_id, ?, string_agg(DISTINCT {label} || ':' || overlap_bp, ',')
            FROM {_relation(conn, path)} WHERE overlap_bp >= ? GROUP BY protein_id""", [flag, min_overlap_bp])
        n = conn.execute("SELECT COUNT(*) FROM gene_flags WHERE flag = ?", [flag]).fetchone()[0]
        conn.execute("INSERT INTO imports VALUES (?, ?, ?, ?, ?, ?)",
                     [str(path), _sha256(path), f"gene_flags:{flag}", [], n, datetime.now(timezone.utc)])
        conn.execute("COMMIT")
    except Exception:
        with contextlib.suppress(duckdb.Error):
            conn.execute("ROLLBACK")
        raise
    finally:
        conn.close()
    return {"sidecar": str(sidecar), "flag": flag, "flagged_proteins": n, "min_overlap_bp": min_overlap_bp}


# --------------------------------------------------------------------------- #
# Queries
# --------------------------------------------------------------------------- #


def _unflagged(include_flagged: bool, alias: str = "") -> str:
    col = f"{alias}.protein_id" if alias else "protein_id"
    return "" if include_flagged else f" AND {col} NOT IN (SELECT protein_id FROM gene_flags)"


def _flagged_ids(conn) -> set[str]:
    return {r[0] for r in conn.execute("SELECT DISTINCT protein_id FROM gene_flags").fetchall()}


def samples(sidecar: Path) -> list[dict[str, Any]]:
    with connect(sidecar, read_only=True) as conn:
        return [{"sample_id": s, "metadata": json.loads(m) if m else {}, "sense_orientation": o,
                 "strand_specificity": p, "sense_fragments": n, "flagged_sense_fragments": f, "genes": g}
                for s, m, o, p, n, f, g in conn.execute(
                    """SELECT s.sample_id, s.metadata, s.sense_orientation, s.strand_specificity, s.sense_fragments,
                              COALESCE(SUM(c.sense) FILTER (WHERE c.protein_id IN (SELECT protein_id FROM gene_flags)), 0),
                              COUNT(c.protein_id) FILTER (WHERE c.sense > 0)
                       FROM rna_samples s LEFT JOIN gene_counts c USING (sample_id) GROUP BY ALL ORDER BY 1""").fetchall()]


def _sample_clause(sample_ids: list[str] | None, column: str = "sample_id") -> tuple[str, list[str]]:
    if not sample_ids:
        return "", []
    return f" AND {column} IN ({','.join('?' * len(sample_ids))})", list(sample_ids)


def genome_expression(sidecar: Path, *, sample_ids: list[str] | None = None,
                      bins: list[str] | None = None, include_flagged: bool = False) -> list[dict[str, Any]]:
    """Per genome and sample: sense fragments, genes with reads, share of the sample's sense fragments."""
    clause, params = _sample_clause(sample_ids)
    with connect(sidecar, read_only=True) as conn:
        rows = conn.execute(f"""
            WITH g AS (SELECT sample_id, bin_id, SUM(sense) AS fragments,
                              COUNT(*) FILTER (WHERE sense > 0) AS genes_with_reads
                       FROM gene_counts WHERE TRUE {clause}{_unflagged(include_flagged)} GROUP BY 1, 2)
            SELECT sample_id, bin_id, fragments, genes_with_reads,
                   fragments / NULLIF(SUM(fragments) OVER (PARTITION BY sample_id), 0) AS transcript_share
            FROM g ORDER BY sample_id, fragments DESC, bin_id""", params).fetchall()
    keep = set(bins) if bins else None
    return [{"sample_id": s, "bin_id": b, "fragments": f, "genes_with_reads": n, "transcript_share": x}
            for s, b, f, n, x in rows if keep is None or b in keep]


def _labels(store, protein_ids: list[str]) -> dict[str, str]:
    if not protein_ids:
        return {}
    marks = ",".join("?" * len(protein_ids))
    rows = store.execute(
        f"""SELECT protein_id, source, string_agg(DISTINCT COALESCE(name, accession), '+')
            FROM annotations WHERE protein_id IN ({marks}) AND source IN ('kegg', 'kofam', 'pfam')
            GROUP BY 1, 2""", protein_ids)
    out: dict[str, dict[str, str]] = {}
    for pid, src, label in rows:
        out.setdefault(pid, {})["kegg" if src in ("kegg", "kofam") else "pfam"] = label
    return {pid: "; ".join(f"{k}: {v}" for k, v in sorted(d.items())) for pid, d in out.items()}


def _within_genome_percentile(rate):
    """Share of the genome's genes at or below each gene's rate; genes without reads get 0.

    Ranking ties with ``method="max"`` alone would give every read-less gene the read-less fraction of the genome
    (often 0.8-0.9 in a sparsely covered genome), which reads as a highly expressed gene.
    """
    return rate.rank(pct=True, method="max").where(rate > 0, 0.0)


def genome_genes(store, sidecar: Path, bin_id: str, *, sample_ids: list[str] | None = None,
                 top: int = 20, min_fragments: int = 0, include_flagged: bool = False) -> dict[str, Any]:
    """Most expressed genes of one genome per sample, with within-genome TPM and percentile."""
    import pandas as pd  # noqa: PLC0415

    genes = pd.DataFrame(store.execute(
        "SELECT protein_id, (end_coord - start + 1) AS gene_length FROM proteins WHERE bin_id = ?",
        [bin_id]), columns=["protein_id", "gene_length"])
    if genes.empty:
        raise ValueError(f"No proteins for genome {bin_id}")
    clause, params = _sample_clause(sample_ids)
    with connect(sidecar, read_only=True) as conn:
        counts = conn.execute(f"SELECT sample_id, protein_id, sense FROM gene_counts WHERE bin_id = ? {clause}"
                              f"{_unflagged(include_flagged)}", [bin_id, *params]).df()
        flagged = set() if include_flagged else _flagged_ids(conn)
    genes = genes[~genes["protein_id"].isin(flagged)]
    result = []
    for sample, part in counts.groupby("sample_id"):
        if part["sense"].sum() < max(min_fragments, 1):
            continue
        frame = genes.merge(part[["protein_id", "sense"]], on="protein_id", how="left").fillna({"sense": 0})
        frame["rate"] = frame["sense"] / (frame["gene_length"] / 1000)
        frame["genome_tpm"] = frame["rate"] / frame["rate"].sum() * 1e6
        frame["percentile"] = _within_genome_percentile(frame["rate"])
        best = frame.sort_values("rate", ascending=False).head(top)
        result.append({"sample_id": sample, "fragments": int(part["sense"].sum()),
                       "genes_with_reads": int((frame["sense"] > 0).sum()), "genes": len(frame),
                       "top": best[["protein_id", "sense", "genome_tpm", "percentile"]].to_dict("records")})
    labels = _labels(store, sorted({g["protein_id"] for r in result for g in r["top"]}))
    for r in result:
        for g in r["top"]:
            g["label"] = labels.get(g["protein_id"], "")
    return {"bin_id": bin_id, "samples": sorted(result, key=lambda r: -r["fragments"])}


def feature_expression(store, sidecar: Path, *, annotation: str | None = None, predicate: str | None = None,
                       sample_ids: list[str] | None = None, min_genome_fragments: int = 100,
                       include_flagged: bool = False) -> dict[str, Any]:
    """Share of each sample's sense fragments on proteins carrying a feature, and where the feature's genes
    rank inside their own genomes (samples where the genome has >= min_genome_fragments)."""
    import pandas as pd  # noqa: PLC0415

    if (predicate is None) == (annotation is None):
        raise ValueError("Give exactly one of predicate or annotation")
    if predicate is not None:
        carriers = [r[0] for r in store.execute(
            "SELECT DISTINCT protein_id FROM protein_predicates WHERE list_contains(predicates, ?)", [predicate])]
        feature = {"predicate": predicate}
    else:
        carriers = [r[0] for r in store.execute(
            "SELECT DISTINCT protein_id FROM annotations WHERE accession = ? OR name = ?", [annotation, annotation])]
        feature = {"annotation": annotation}
    carrier_set = set(carriers)
    if not carriers:
        return {**feature, "carrier_proteins": 0, "samples": [], "genomes": []}
    bins = sorted({r[0] for r in store.execute(
        f"SELECT DISTINCT bin_id FROM proteins WHERE protein_id IN ({','.join('?' * len(carriers))})", carriers)})
    genes = pd.DataFrame(store.execute(
        f"SELECT protein_id, bin_id, (end_coord - start + 1) AS gene_length FROM proteins "
        f"WHERE bin_id IN ({','.join('?' * len(bins))})", bins),
        columns=["protein_id", "bin_id", "gene_length"])
    clause, params = _sample_clause(sample_ids)
    with connect(sidecar, read_only=True) as conn:
        totals = dict(conn.execute(f"SELECT sample_id, SUM(sense) FROM gene_counts WHERE TRUE {clause}"
                                   f"{_unflagged(include_flagged)} GROUP BY 1", params).fetchall())
        flagged = set() if include_flagged else _flagged_ids(conn)
        conn.register("_bins", pd.DataFrame({"bin_id": bins}))
        counts = conn.execute(f"SELECT sample_id, protein_id, sense FROM gene_counts JOIN _bins USING (bin_id) "
                              f"WHERE TRUE {clause}{_unflagged(include_flagged)}", params).df()
    genes = genes[~genes["protein_id"].isin(flagged)]
    carrier_set -= flagged
    per_sample, per_genome = [], {}
    for sample in sorted(totals):
        part = counts[counts["sample_id"] == sample]
        on_feature = int(part.loc[part["protein_id"].isin(carrier_set), "sense"].sum())
        per_sample.append({"sample_id": sample, "fragments": on_feature,
                           "share_of_sample": on_feature / totals[sample] if totals[sample] else None})
        frame = genes.merge(part[["protein_id", "sense"]], on="protein_id", how="left").fillna({"sense": 0})
        frame["rate"] = frame["sense"] / (frame["gene_length"] / 1000)
        for bin_id, g in frame.groupby("bin_id"):
            if g["sense"].sum() < min_genome_fragments:
                continue
            pct = _within_genome_percentile(g["rate"])
            mine = pct[g["protein_id"].isin(carrier_set)]
            if mine.empty:
                continue
            per_genome.setdefault(bin_id, []).append({"sample_id": sample, "genome_fragments": int(g["sense"].sum()),
                                                      "best_percentile": float(mine.max())})
    genomes = [{"bin_id": b, "samples": len(v),
                "median_best_percentile": float(pd.Series([x["best_percentile"] for x in v]).median()),
                "max_genome_fragments": max(x["genome_fragments"] for x in v)} for b, v in per_genome.items()]
    genomes.sort(key=lambda g: (-g["median_best_percentile"], g["bin_id"]))
    return {**feature, "carrier_proteins": len(carriers), "carrier_genomes": len(bins),
            "min_genome_fragments": min_genome_fragments, "samples": per_sample, "genomes": genomes}


# --------------------------------------------------------------------------- #
# Markdown
# --------------------------------------------------------------------------- #


def genome_expression_markdown(rows: list[dict[str, Any]], top: int = 10) -> str:
    lines = ["# Genome transcript share"]
    by_sample: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        by_sample.setdefault(row["sample_id"], []).append(row)
    for sample, items in by_sample.items():
        lines.append(f"\n## {sample} ({len(items)} genomes with reads)")
        for r in items[:top]:
            share = f"{r['transcript_share']:.2%}" if r["transcript_share"] is not None else "n/a"
            lines.append(f"- {r['bin_id']}: {share} of sense fragments ({r['fragments']:,}; "
                         f"{r['genes_with_reads']:,} genes with reads)")
    return "\n".join(lines)


def genome_genes_markdown(result: dict[str, Any]) -> str:
    lines = [f"# Most expressed genes: {result['bin_id']}"]
    for s in result["samples"]:
        lines.append(f"\n## {s['sample_id']} ({s['fragments']:,} fragments; "
                     f"{s['genes_with_reads']:,} of {s['genes']:,} genes with reads)")
        for g in s["top"]:
            lines.append(f"- {g['protein_id']}: {int(g['sense']):,} fragments, genome TPM {g['genome_tpm']:,.0f}, "
                         f"percentile {g['percentile']:.3f}" + (f" — {g['label']}" if g["label"] else ""))
    return "\n".join(lines)


def feature_expression_markdown(result: dict[str, Any], top: int = 15) -> str:
    feature = result.get("predicate") or result.get("annotation")
    lines = [f"# Expression of {feature}",
             f"{result['carrier_proteins']:,} proteins in {result.get('carrier_genomes', 0):,} genomes carry it."]
    for s in result["samples"]:
        share = f"{s['share_of_sample']:.3%}" if s["share_of_sample"] is not None else "n/a"
        lines.append(f"- {s['sample_id']}: {s['fragments']:,} fragments ({share} of sense fragments)")
    if result["genomes"]:
        lines.append(f"\nWithin-genome rank (samples with >= {result['min_genome_fragments']} genome fragments):")
        for g in result["genomes"][:top]:
            lines.append(f"- {g['bin_id']}: median best percentile {g['median_best_percentile']:.3f} "
                         f"over {g['samples']} samples")
    return "\n".join(lines)
