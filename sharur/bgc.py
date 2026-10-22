"""Biosynthetic gene clusters: load GECCO calls into an existing database.

GECCO's cluster rows become ``loci`` with ``locus_type='bgc'`` (the
``bgc_loci`` view) and their member proteins become ``locus_proteins``.
The runner (``scripts/run_gecco_bgcs.py``) hands GECCO the database's own
gene calls, so every cluster member is a stored protein. It names GenBank
records with short IDs, and each shard's ``record_map.tsv`` maps an ID back to
``(bin_id, contig_id)``.

Loading replaces every existing BGC locus, so a rerun is idempotent.
"""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any

CALLER = "gecco"
PROBABILITY_COLUMNS = (
    "alkaloid_probability", "nrp_probability", "polyketide_probability",
    "ripp_probability", "saccharide_probability", "terpene_probability",
)


def _float(value: str | None) -> float | None:
    try:
        return float(value) if value not in (None, "") else None
    except ValueError:
        return None


def read_record_map(path: Path) -> dict[str, tuple[str, str]]:
    """record name -> (bin_id, contig_id) from a runner shard's record_map.tsv."""
    with path.open() as handle:
        return {row["record"]: (row["bin_id"], row["contig_id"])
                for row in csv.DictReader(handle, delimiter="\t")}


def collect_gecco_clusters(gecco_dir: Path) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Parse every shard under ``gecco_dir`` into cluster dicts.

    A shard is a directory holding ``record_map.tsv`` and one GECCO
    ``*.clusters.tsv``. Shards without a finished clusters table are reported
    in ``meta["missing_shards"]``.
    """
    clusters: list[dict[str, Any]] = []
    meta: dict[str, Any] = {"shards": 0, "missing_shards": [], "gecco_version": None}
    for record_map in sorted(gecco_dir.glob("*/record_map.tsv")):
        shard = record_map.parent
        meta["shards"] += 1
        tables = sorted(shard.glob("*.clusters.tsv"))
        if not tables:
            meta["missing_shards"].append(shard.name)
            continue
        version_file = shard / "gecco_version.txt"
        if version_file.is_file() and meta["gecco_version"] is None:
            meta["gecco_version"] = version_file.read_text().strip()
        records = read_record_map(record_map)
        with tables[0].open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                record = row.get("sequence_id", "")
                if record not in records:
                    raise ValueError(f"{tables[0]}: record {record!r} is absent from {record_map}")
                bin_id, contig_id = records[record]
                proteins = [p for p in (row.get("proteins") or "").split(";") if p]
                clusters.append({
                    "bin_id": bin_id,
                    "contig_id": contig_id,
                    "gecco_id": row.get("bgc_id") or row.get("cluster_id") or "",
                    "start": int(row["start"]),
                    "end": int(row["end"]),
                    "bgc_type": row.get("type") or "Unknown",
                    "average_p": _float(row.get("average_p")),
                    "max_p": _float(row.get("max_p")),
                    "probabilities": {c.replace("_probability", ""): _float(row.get(c))
                                      for c in PROBABILITY_COLUMNS if c in row},
                    "domains": [d for d in (row.get("domains") or "").split(";") if d],
                    "proteins": proteins,
                })
    return clusters, meta


def load_gecco_bgcs(conn, gecco_dir: Path) -> dict[str, Any]:
    """Replace the database's BGC loci with the GECCO calls under ``gecco_dir``.

    Locus IDs are ``{bin_id}::bgc{n}``, numbered per genome in contig and
    coordinate order. Member proteins absent from ``proteins`` (or stored under
    another genome) are dropped and counted; a cluster left with no member is
    skipped and counted.
    """
    import pandas as pd  # noqa: PLC0415

    clusters, meta = collect_gecco_clusters(Path(gecco_dir))
    if meta["missing_shards"]:
        raise ValueError(f"GECCO shards without a clusters table: {meta['missing_shards']}")
    member_ids = sorted({p for c in clusters for p in c["proteins"]})
    owner: dict[str, tuple[str, str]] = {}
    if member_ids:
        conn.register("_bgc_members", pd.DataFrame({"protein_id": member_ids}))
        owner = {pid: (bin_id, contig_id) for pid, bin_id, contig_id in conn.execute(
            """SELECT p.protein_id, p.bin_id, p.contig_id FROM proteins p
               JOIN _bgc_members m USING (protein_id)""").fetchall()}
        conn.unregister("_bgc_members")

    loci_rows, member_rows = [], []
    dropped_members = skipped_clusters = 0
    per_bin: Counter[str] = Counter()
    types: Counter[str] = Counter()
    for cluster in sorted(clusters, key=lambda c: (c["bin_id"], c["contig_id"], c["start"])):
        members = [p for p in cluster["proteins"]
                   if owner.get(p) == (cluster["bin_id"], cluster["contig_id"])]
        dropped_members += len(cluster["proteins"]) - len(members)
        if not members:
            skipped_clusters += 1
            continue
        per_bin[cluster["bin_id"]] += 1
        locus_id = f"{cluster['bin_id']}::bgc{per_bin[cluster['bin_id']]}"
        types[cluster["bgc_type"]] += 1
        metadata = {
            "caller": CALLER, "gecco_version": meta["gecco_version"],
            "gecco_id": cluster["gecco_id"], "bin_id": cluster["bin_id"],
            "bgc_type": cluster["bgc_type"], "average_p": cluster["average_p"],
            "max_p": cluster["max_p"], "type_probabilities": cluster["probabilities"],
            "domains": cluster["domains"], "n_proteins": len(members),
        }
        loci_rows.append((locus_id, "bgc", cluster["contig_id"], cluster["start"], cluster["end"],
                          cluster["average_p"], json.dumps(metadata)))
        member_rows.extend((locus_id, pid, idx) for idx, pid in enumerate(members))

    previous = conn.execute("SELECT COUNT(*) FROM loci WHERE locus_type = 'bgc'").fetchone()[0]
    # DuckDB checks a foreign key against the committed child table, so member
    # rows are deleted and committed before their parent loci. A failure after
    # this point leaves old loci without members; rerunning the load repairs it.
    conn.execute("""DELETE FROM locus_proteins WHERE locus_id IN
                    (SELECT locus_id FROM loci WHERE locus_type = 'bgc')""")
    conn.execute("BEGIN TRANSACTION")
    try:
        conn.execute("DELETE FROM loci WHERE locus_type = 'bgc'")
        if loci_rows:
            conn.register("_bgc_loci", pd.DataFrame(loci_rows, columns=[
                "locus_id", "locus_type", "contig_id", "start", "end_coord", "confidence", "metadata"]))
            conn.execute("""INSERT INTO loci (locus_id, locus_type, contig_id, start, end_coord,
                                              confidence, metadata)
                            SELECT * FROM _bgc_loci""")
            conn.unregister("_bgc_loci")
            conn.register("_bgc_lp", pd.DataFrame(member_rows, columns=["locus_id", "protein_id", "position"]))
            conn.execute("INSERT INTO locus_proteins (locus_id, protein_id, position) SELECT * FROM _bgc_lp")
            conn.unregister("_bgc_lp")
        # the same view Stage 07 defines, so databases built without it gain the query surface
        conn.execute("CREATE OR REPLACE VIEW bgc_loci AS SELECT * FROM loci WHERE locus_type = 'bgc'")
        conn.execute("COMMIT")
    except Exception:
        conn.execute("ROLLBACK")
        raise
    return {
        "clusters": len(loci_rows), "replaced": previous, "genomes": len(per_bin),
        "member_proteins": len(member_rows), "dropped_members": dropped_members,
        "skipped_clusters": skipped_clusters, "by_type": dict(types.most_common()),
        "shards": meta["shards"], "gecco_version": meta["gecco_version"],
    }
