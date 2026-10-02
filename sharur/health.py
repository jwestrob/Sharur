"""Dataset health: fast checks for the problems that quietly distort analyses.

Each check reads the dataset (never writes), runs in well under a second on a
few million proteins, and returns a :class:`Finding`: a status, one line with
exact counts, a few example IDs, and the command that fixes it.

Statuses: ``ok`` (nothing to do), ``info`` (worth knowing), ``warn`` (results
that depend on it need care), ``fail`` (fix before analysing).

Use :func:`run_checks` with a database path, or with a query function when a
connection is already open (the browser passes one that takes its lock per
query).
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

STATUSES = ("ok", "info", "warn", "fail")
AREAS = ("Database", "Genomes", "Genes", "Annotations", "Callers", "References")
EXAMPLES = 10
BROAD_SOURCE = 0.10   # sources hitting at least this share of a typical genome's proteins
THIN_SOURCE = 0.25    # under this fraction of a broad source's typical share, a genome looks partly annotated
GAP_EXPECTED = 7.0    # hits a genome's size predicts from a source; with none, P(0) = e^-7 < 0.1%: never searched
# Rows written by system callers and classifiers: a genome without them holds no call, a result in itself.
CALLER_SOURCES = ("defensefinder_system", "txsscan_system", "hyddb_subgroup")


@dataclass
class Finding:
    key: str
    area: str
    title: str
    status: str
    summary: str
    fix: str = ""
    examples: list[dict[str, str]] = field(default_factory=list)   # {"kind": genome|protein|contig, "id", "note"}
    table: dict[str, Any] | None = None                            # {"columns": [...], "rows": [[...]]}

    def to_dict(self) -> dict[str, Any]:
        return {"key": self.key, "area": self.area, "title": self.title, "status": self.status,
                "summary": self.summary, "fix": self.fix, "examples": self.examples, "table": self.table}


@dataclass
class HealthReport:
    db_path: str
    findings: list[Finding]
    seconds: float

    @property
    def counts(self) -> dict[str, int]:
        return {s: sum(1 for f in self.findings if f.status == s) for s in STATUSES}

    @property
    def worst(self) -> str:
        present = {f.status for f in self.findings}
        return next((s for s in reversed(STATUSES) if s in present), "ok")

    def by_area(self) -> list[tuple[str, list[Finding]]]:
        return [(a, [f for f in self.findings if f.area == a]) for a in AREAS
                if any(f.area == a for f in self.findings)]

    def to_dict(self) -> dict[str, Any]:
        return {"db_path": self.db_path, "seconds": round(self.seconds, 2), "worst": self.worst,
                "counts": self.counts, "findings": [f.to_dict() for f in self.findings]}

    def to_markdown(self) -> str:
        mark = {"ok": "✓", "info": "·", "warn": "!", "fail": "✗"}
        c = self.counts
        lines = [f"# Dataset health: {Path(self.db_path).parent.name or self.db_path}",
                 f"{c['fail']} fail · {c['warn']} warn · {c['info']} info · {c['ok']} ok "
                 f"({self.seconds:.1f} s)"]
        for area, findings in self.by_area():
            lines += ["", f"## {area}"]
            for f in findings:
                lines.append(f"- [{mark[f.status]} {f.status}] **{f.title}**: {f.summary}")
                if f.examples:
                    lines.append("  e.g. " + ", ".join(
                        e["id"] + (f" ({e['note']})" if e.get("note") else "") for e in f.examples[:5]))
                if f.fix and f.status != "ok":
                    lines.append(f"  fix: {f.fix}")
        return "\n".join(lines)


class Context:
    """What the checks share: a query function, the table list, and the dataset folder."""

    def __init__(self, query: Callable[..., list], db_path: Path, assemblies: dict[str, Path] | None = None):
        self.q = query
        self.db_path = Path(db_path)
        self.dataset_dir = self.db_path.resolve().parent
        self.tables = {r[0] for r in query(
            "SELECT table_name FROM information_schema.tables WHERE table_schema = 'main'")}
        self._assemblies = assemblies
        self._genomes: dict[str, dict[str, Any]] | None = None

    def scalar(self, sql: str, params: list | None = None) -> Any:
        rows = self.q(sql, params or [])
        return rows[0][0] if rows else None

    @property
    def assemblies(self) -> dict[str, Path]:
        if self._assemblies is None:
            from sharur.assemblies import find_assemblies  # noqa: PLC0415

            self._assemblies = find_assemblies(self.dataset_dir)
        return self._assemblies

    @property
    def genomes(self) -> dict[str, dict[str, Any]]:
        """bin_id -> {proteins, completeness, contamination}, every genome in ``bins``."""
        if self._genomes is None:
            rows = self.q("""
                SELECT b.bin_id, b.completeness, b.contamination, COALESCE(p.n, 0)
                FROM bins b LEFT JOIN (SELECT bin_id, COUNT(*) AS n FROM proteins GROUP BY 1) p USING (bin_id)""")
            self._genomes = {b: {"completeness": comp, "contamination": cont, "proteins": n}
                             for b, comp, cont, n in rows}
        return self._genomes


def _ex(kind: str, ids, note: Callable[[str], str] | None = None) -> list[dict[str, str]]:
    return [{"kind": kind, "id": str(i), "note": note(i) if note else ""} for i in list(ids)[:EXAMPLES]]


def _count(n: int, noun: str = "genome") -> str:
    return f"{n:,} {noun}{'' if n == 1 else 's'}"


def _verb(n: int, singular: str, plural: str) -> str:
    return singular if n == 1 else plural


def _pct(n: int, d: int) -> str:
    return f"{100 * n / d:.1f}%" if d else "0%"


# --------------------------------------------------------------------------- #
# Database
# --------------------------------------------------------------------------- #


def check_schema(ctx: Context) -> Finding:
    from sharur.storage.schema import SCHEMA_VERSION  # noqa: PLC0415

    version = ctx.scalar("SELECT MAX(version) FROM schema_version") if "schema_version" in ctx.tables else None
    if version is None:
        return Finding("schema", "Database", "Schema version", "fail",
                       "The database records no schema version.", "Rebuild it with `sharur-ingest`.")
    if version < SCHEMA_VERSION:
        return Finding("schema", "Database", "Schema version", "warn",
                       f"Schema {version}; this Sharur expects {SCHEMA_VERSION}. Reading works; writes migrate first.",
                       "Review the pending migrations, then run `sharur migrate --db DATASET/sharur.duckdb`.")
    if version > SCHEMA_VERSION:
        return Finding("schema", "Database", "Schema version", "warn",
                       f"Schema {version} is newer than this Sharur ({SCHEMA_VERSION}).", "Update Sharur.")
    return Finding("schema", "Database", "Schema version", "ok", f"Schema {version}, as this Sharur expects.")


def check_seal(ctx: Context) -> Finding:
    """Structural comparison only: table row counts and schema version against the seal."""
    from sharur.dataset_seal import DEFAULT_SEAL_NAME  # noqa: PLC0415

    seal_path = ctx.dataset_dir / DEFAULT_SEAL_NAME
    if not seal_path.is_file():
        return Finding("seal", "Database", "Dataset seal", "info", "The dataset has no seal.",
                       "Record its state with `sharur seal --db DATASET/sharur.duckdb` once writes are done.")
    try:
        seal = json.loads(seal_path.read_text())
        sealed = seal["identity"]["database"]
    except (ValueError, KeyError, OSError) as exc:
        return Finding("seal", "Database", "Dataset seal", "warn", f"The seal could not be read ({exc}).",
                       "Reseal with `sharur seal --db DATASET/sharur.duckdb --force`.")
    when = str(seal.get("generated_at", ""))[:10]
    changed = []
    for row in sealed.get("tables", []):
        table = row.get("table")
        if table not in ctx.tables:
            changed.append((table, row.get("rows"), None))
            continue
        now = ctx.scalar(f'SELECT COUNT(*) FROM "{table}"')
        if now != row.get("rows"):
            changed.append((table, row.get("rows"), now))
    base = {r[0] for r in ctx.q("SELECT table_name FROM information_schema.tables "
                                "WHERE table_schema = 'main' AND table_type = 'BASE TABLE'")}
    new_tables = base - {r.get("table") for r in sealed.get("tables", [])}
    current_version = ctx.scalar("SELECT MAX(version) FROM schema_version") if "schema_version" in ctx.tables else None
    if sealed.get("schema_version") != current_version:
        changed.append(("schema_version", sealed.get("schema_version"), current_version))
    if changed or new_tables:
        parts = [f"{t}: {a:,} → {b:,}" if isinstance(a, int) and isinstance(b, int) else f"{t}: {a} → {b}"
                 for t, a, b in changed]
        return Finding("seal", "Database", "Dataset seal", "warn",
                       f"{len(changed) + len(new_tables)} tables differ from the seal of {when} "
                       f"(row counts or presence).",
                       "After deliberate writes, reseal with `sharur seal --db DATASET/sharur.duckdb --force`; "
                       "`sharur verify-seal` gives the full comparison.",
                       examples=[{"kind": "text", "id": p, "note": ""} for p in parts[:EXAMPLES]]
                       + [{"kind": "text", "id": f"{t}: new", "note": ""} for t in sorted(new_tables)[:EXAMPLES]])
    return Finding("seal", "Database", "Dataset seal", "ok",
                   f"Table row counts and schema match the {seal.get('seal_strength', '')} seal of {when}.".replace("  ", " "))



# --------------------------------------------------------------------------- #
# Genomes
# --------------------------------------------------------------------------- #


def check_completeness(ctx: Context) -> Finding:
    genomes = ctx.genomes
    n = len(genomes)
    known = [g["completeness"] for g in genomes.values() if g["completeness"] is not None]
    fix = ("Import CheckM, CheckM2 or GTDB estimates with `sharur import-quality TABLE --db DATASET/sharur.duckdb`.")
    if not n:
        return Finding("completeness", "Genomes", "Completeness estimates", "fail", "The dataset has no genomes.")
    if not known:
        return Finding("completeness", "Genomes", "Completeness estimates", "warn",
                       f"None of {n:,} genomes has a completeness estimate, so absences can't be weighed "
                       "against incompleteness.", fix)
    known.sort()
    med = known[len(known) // 2]
    missing = [b for b, g in genomes.items() if g["completeness"] is None]
    if missing:
        return Finding("completeness", "Genomes", "Completeness estimates", "warn",
                       f"{len(known):,} of {n:,} genomes have completeness (median {med:.1f}%); "
                       f"{len(missing):,} {_verb(len(missing), 'has', 'have')} none.", fix, _ex("genome", sorted(missing)))
    low = sum(1 for c in known if c < 50)
    return Finding("completeness", "Genomes", "Completeness estimates", "ok",
                   f"All {n:,} genomes have completeness: median {med:.1f}%, {low:,} below 50%.")


def check_contamination(ctx: Context) -> Finding:
    genomes = ctx.genomes
    known = {b: g["contamination"] for b, g in genomes.items() if g["contamination"] is not None}
    if not known:
        return Finding("contamination", "Genomes", "Contamination", "info", "No contamination estimates recorded.",
                       "`sharur import-quality` fills them with completeness.")
    high = sorted((b for b, c in known.items() if c > 10), key=lambda b: -known[b])
    if high:
        return Finding("contamination", "Genomes", "Contamination", "warn",
                       f"{len(high):,} of {len(known):,} genomes {_verb(len(high), 'exceeds', 'exceed')} 10% contamination; their gene "
                       "content mixes organisms.",
                       "Check them before attributing their genes to their lineage.",
                       _ex("genome", high, lambda b: f"{known[b]:.1f}%"))
    return Finding("contamination", "Genomes", "Contamination", "ok",
                   f"All {len(known):,} estimated genomes are at or below 10% contamination.")


def check_empty_genomes(ctx: Context) -> Finding:
    empty = sorted(b for b, g in ctx.genomes.items() if not g["proteins"])
    if empty:
        return Finding("empty_genomes", "Genomes", "Genomes without proteins", "fail",
                       f"{_count(len(empty))} {_verb(len(empty), 'has', 'have')} no proteins.",
                       "Re-run gene calling (stage 03) for them and rebuild.", _ex("genome", empty))
    return Finding("empty_genomes", "Genomes", "Genomes without proteins", "ok",
                   f"Every one of {len(ctx.genomes):,} genomes has proteins.")


def check_gene_calls(ctx: Context) -> Finding:
    from sharur.assemblies import DEFICIT, assembly_kb, gene_call_deficits  # noqa: PLC0415

    genomes = ctx.genomes
    kb = assembly_kb({b: p for b, p in ctx.assemblies.items() if b in genomes})
    proteins = {b: g["proteins"] for b, g in genomes.items() if g["proteins"]}
    flagged = gene_call_deficits(proteins, {b: g["completeness"] for b, g in genomes.items()}, kb)
    judged_by_completeness = sum(1 for b, g in genomes.items() if b not in kb and g["completeness"])
    if not kb and not judged_by_completeness:
        return Finding("gene_calls", "Genes", "Gene calls per genome", "info",
                       "Gene calls can't be weighed: no assembly files sit beside the dataset and no genome has "
                       "a completeness estimate.",
                       "Place assemblies in the dataset folder (genomes_fna/ or source/) or run "
                       "`sharur import-quality`.")
    basis = (f"{len(kb):,} of {len(genomes):,} genomes judged against their assembly"
             + (f", {judged_by_completeness:,} against completeness" if judged_by_completeness else ""))
    if flagged:
        order = sorted(flagged, key=lambda b: (proteins[b] / kb[b]) if kb.get(b) else 1.0)
        return Finding("gene_calls", "Genes", "Gene calls per genome", "warn",
                       f"{_count(len(flagged))} {_verb(len(flagged), 'carries', 'carry')} under {DEFICIT} proteins per kb of assembly (or half what "
                       f"completeness implies), so most of their gene calls are missing. {basis}.",
                       "Re-run gene calling (stage 03) on their assemblies and re-annotate; meanwhile treat "
                       "their absences as missing data.",
                       _ex("genome", order, lambda b: (f"{proteins[b]:,} proteins on {kb[b]:,.0f} kb" if kb.get(b)
                                                       else f"{proteins[b]:,} proteins")))
    return Finding("gene_calls", "Genes", "Gene calls per genome", "ok",
                   f"No genome falls below {DEFICIT} proteins per kb. {basis}.")


# --------------------------------------------------------------------------- #
# Genes
# --------------------------------------------------------------------------- #


def check_strand(ctx: Context) -> Finding:
    rows = ctx.q("SELECT strand, COUNT(*) FROM proteins GROUP BY 1 ORDER BY 2 DESC")
    counts = {str(s): n for s, n in rows}
    symbolic = {s for s in counts if s in ("+", "-")}
    numeric = {s for s in counts if s in ("1", "-1")}
    other = {s for s in counts if s not in ("+", "-", "1", "-1")}
    detail = ", ".join(f"'{s}' {n:,}" for s, n in counts.items())
    if other:
        return Finding("strand", "Genes", "Strand encoding", "warn",
                       f"Unexpected strand values: {detail}.", "Gene order and operon logic assume '+' and '-'.")
    if symbolic and numeric:
        return Finding("strand", "Genes", "Strand encoding", "warn",
                       f"Two encodings mixed: {detail}. Code comparing strands must accept both.",
                       "Treat '1' as '+' and '-1' as '-' when reading; stage 07 writes '+'/'-'.")
    return Finding("strand", "Genes", "Strand encoding", "ok", f"One encoding: {detail}.")


def check_positionless(ctx: Context) -> Finding:
    self_rows = ctx.q("""SELECT COUNT(*), COUNT(DISTINCT bin_id) FROM proteins WHERE contig_id = protein_id""")[0]
    shared = ctx.q("""
        SELECT contig_id, COUNT(DISTINCT bin_id) AS g, COUNT(*) AS n FROM proteins
        GROUP BY 1 HAVING COUNT(DISTINCT bin_id) > 1 ORDER BY n DESC""")
    # stacked at 0 on contigs that are otherwise ordinary (shared and self contigs counted above)
    stacked = ctx.scalar("""
        WITH shared AS (SELECT contig_id FROM proteins GROUP BY 1 HAVING COUNT(DISTINCT bin_id) > 1)
        SELECT COALESCE(SUM(n), 0) FROM (
            SELECT COUNT(*) AS n FROM proteins
            WHERE start = 0 AND contig_id <> protein_id AND contig_id NOT IN (SELECT contig_id FROM shared)
            GROUP BY contig_id, start, end_coord HAVING COUNT(*) > 1)""")
    total = ctx.scalar("SELECT COUNT(*) FROM proteins")
    n_self, g_self = self_rows
    n_shared = sum(r[2] for r in shared)
    lost = n_self + n_shared
    if not lost and not stacked:
        return Finding("positionless", "Genes", "Genes without a genomic position", "ok",
                       f"Every one of {total:,} proteins sits on its own genome's contig.")
    parts = []
    if n_self:
        parts.append(f"{_count(n_self, 'protein')} in {_count(g_self)} {_verb(n_self, 'is its', 'are their')} "
                     "own contig")
    if shared:
        parts.append(f"{_count(n_shared, 'protein')} share {_count(len(shared), 'contig ID')} across genomes")
    if stacked:
        parts.append(f"{_count(stacked, 'protein')} {_verb(stacked, 'is', 'are')} stacked at coordinate 0")
    examples = [{"kind": "contig", "id": c, "note": f"{n:,} proteins in {g} genomes"} for c, g, n in shared[:5]]
    if n_self:
        examples += _ex("protein", [r[0] for r in ctx.q(
            "SELECT protein_id FROM proteins WHERE contig_id = protein_id ORDER BY protein_id LIMIT 5")],
            lambda _: "own contig")
    return Finding("positionless", "Genes", "Genes without a genomic position", "warn",
                   "; ".join(parts) + f" ({_pct(lost + stacked, total)} of proteins). Neighborhoods, operons "
                   "and contig-edge logic skip or mistreat them.",
                   "Map these proteins to their assembly coordinates (stage 03 gene calls) and rebuild; "
                   "gene-order features treat them as unplaced meanwhile.", examples)


def check_duplicate_coordinates(ctx: Context) -> Finding:
    rows = ctx.q("""
        SELECT contig_id, start, end_coord, COUNT(*) AS n FROM proteins WHERE start > 0
        GROUP BY 1, 2, 3 HAVING COUNT(*) > 1 ORDER BY n DESC""")
    if not rows:
        return Finding("duplicate_coordinates", "Genes", "Proteins sharing coordinates", "ok",
                       "No two placed proteins share contig coordinates.")
    proteins = sum(r[3] for r in rows)
    return Finding("duplicate_coordinates", "Genes", "Proteins sharing coordinates", "info",
                   f"{proteins:,} proteins share {len(rows):,} coordinate spans (several accessions mapped to one "
                   "locus); counts per genome include each copy.",
                   "Deduplicate by coordinates when counting genes.",
                   [{"kind": "contig", "id": c, "note": f"{n} proteins at {s:,}–{e:,}"} for c, s, e, n in rows[:EXAMPLES]])


def check_contigs(ctx: Context) -> Finding:
    total = ctx.scalar("SELECT COUNT(*) FROM contigs")
    empty, empty_median = ctx.q("SELECT COUNT(*), MEDIAN(length) FROM contigs c ANTI JOIN "
                                "(SELECT DISTINCT contig_id FROM proteins) p USING (contig_id)")[0]
    sources = {}
    if "length_source" in {r[0] for r in ctx.q(
            "SELECT column_name FROM information_schema.columns WHERE table_name = 'contigs'")}:
        sources = dict(ctx.q("SELECT COALESCE(length_source, 'unknown'), COUNT(*) FROM contigs GROUP BY 1"))
    gene_span = sources.get("gene_span", 0) + sources.get("unknown", 0)
    parts = [f"{total:,} contigs"]
    if empty:
        parts.append(f"{empty:,} without genes" + (f" (median {empty_median:,.0f} bp)" if empty_median else ""))
    if gene_span:
        parts.append(f"{gene_span:,} marked with gene-span lengths (the last gene's end) instead of assembly "
                     "lengths")
    status = "info" if (empty or gene_span) else "ok"
    fix = ("Record assembly contig lengths with `sharur backfill-contig-context --db DATASET/sharur.duckdb` "
           "(writes; reseal after)." if gene_span else "")
    return Finding("contigs", "Genes", "Contig lengths", status, "; ".join(parts) + ".", fix)


# --------------------------------------------------------------------------- #
# Annotations
# --------------------------------------------------------------------------- #


def check_annotation_coverage(ctx: Context) -> Finding:
    n_genomes = len(ctx.genomes)
    rows = ctx.q("""
        SELECT LOWER(a.source), COUNT(DISTINCT p.bin_id), COUNT(DISTINCT a.protein_id)
        FROM annotations a JOIN proteins p USING (protein_id) GROUP BY 1 ORDER BY 2 DESC, 3 DESC""")
    total_proteins = ctx.scalar("SELECT COUNT(*) FROM proteins") or 0
    annotated = ctx.scalar("SELECT COUNT(DISTINCT protein_id) FROM annotations") or 0
    if not rows:
        return Finding("annotation_coverage", "Annotations", "Annotation coverage", "fail",
                       "No annotations.", "Run stage 04 (`sharur-ingest`) to annotate proteins.")
    # Broad sources annotate a sizeable share of a typical genome's proteins (Pfam, KOfam, VOGdb); narrow ones
    # (hydrogenases, defense, secretion) hit few genes, so their absence from a genome is expected.
    share = dict(ctx.q("""
        WITH per AS (SELECT LOWER(a.source) AS source, p.bin_id, COUNT(DISTINCT a.protein_id) AS hits
                     FROM annotations a JOIN proteins p USING (protein_id) GROUP BY 1, 2),
        sizes AS (SELECT bin_id, COUNT(*) AS n FROM proteins GROUP BY 1)
        SELECT source, MEDIAN(hits * 1.0 / n) FROM per JOIN sizes USING (bin_id) GROUP BY 1"""))
    broad = {src for src, frac in share.items() if frac is not None and frac >= BROAD_SOURCE}
    table = {"columns": ["Source", "Kind", "Genomes", "Share of genomes", "Proteins", "Typical share of a genome"],
             "rows": [[s, "broad" if s in broad else "narrow", g, _pct(g, n_genomes), p,
                       f"{(share.get(s) or 0) * 100:.0f}%"] for s, g, p in rows]}
    genomes_with = {s: g for s, g, _ in rows}
    subset = sorted(s for s in broad if genomes_with[s] < 0.5 * n_genomes)
    gaps, thin = [], []
    # A genome large enough to expect several hits from a source, yet with none, was most likely never
    # searched with it (an interrupted run leaves exactly this pattern, narrow sources included).
    for source in sorted(set(share) - set(subset) - set(CALLER_SOURCES)):
        missing = [r[0] for r in ctx.q("""
            WITH sizes AS (SELECT bin_id, COUNT(*) AS n FROM proteins GROUP BY 1)
            SELECT z.bin_id FROM sizes z ANTI JOIN (
                SELECT DISTINCT p.bin_id FROM annotations a JOIN proteins p USING (protein_id)
                WHERE LOWER(a.source) = ?) s USING (bin_id)
            WHERE z.n * ? > ? ORDER BY z.bin_id""", [source, share[source] or 0, GAP_EXPECTED])]
        if missing:
            gaps.append((source, missing))
        if source not in broad:
            continue
        thin += [(source, b, f) for b, f in ctx.q("""
            WITH hits AS (SELECT p.bin_id, COUNT(DISTINCT a.protein_id) AS hits FROM annotations a
                          JOIN proteins p USING (protein_id) WHERE LOWER(a.source) = ? GROUP BY 1),
            sizes AS (SELECT bin_id, COUNT(*) AS n FROM proteins GROUP BY 1 HAVING COUNT(*) >= 50)
            SELECT bin_id, hits * 1.0 / n FROM hits JOIN sizes USING (bin_id)
            WHERE hits * 1.0 / n < ? ORDER BY 2, 1""", [source, THIN_SOURCE * share[source]])]
    unannotated_genomes = [r[0] for r in ctx.q("""
        SELECT b.bin_id FROM bins b ANTI JOIN (
            SELECT DISTINCT p.bin_id FROM annotations a JOIN proteins p USING (protein_id)) s USING (bin_id)
        ORDER BY 1""")]
    summary = (f"{len(rows)} sources ({len(broad)} broad); {_pct(annotated, total_proteins)} of "
               f"{total_proteins:,} proteins carry at least one hit.")
    if subset:
        summary += " " + "; ".join(f"{s} was run on {genomes_with[s]:,} of {n_genomes:,} genomes"
                                   for s in subset) + "."
    if unannotated_genomes or gaps or thin:
        parts = [f"{len(m):,} {'genome lacks' if len(m) == 1 else 'genomes lack'} {s} hits that "
                 f"{n_genomes - len(m):,} others have" for s, m in gaps]
        for source in sorted({s for s, _, _ in thin}):
            n_thin = sum(1 for s, _, _ in thin if s == source)
            parts.append(f"{_count(n_thin)} {_verb(n_thin, 'carries', 'carry')} {source} hits on under a quarter of the usual share of their "
                         "proteins")
        if unannotated_genomes:
            parts.append(f"{_count(len(unannotated_genomes))} {_verb(len(unannotated_genomes), 'has', 'have')} no annotation at all")
        examples = _ex("genome", unannotated_genomes, lambda _: "no annotation")
        for s, m in gaps:
            examples += _ex("genome", m[:3], lambda _, s=s: f"no {s}")
        examples += [{"kind": "genome", "id": b, "note": f"{s} on {f:.0%} of proteins"} for s, b, f in thin[:4]]
        return Finding("annotation_coverage", "Annotations", "Annotation coverage", "warn",
                       summary + " " + "; ".join(parts) + ".",
                       "Annotate the missing genomes with that source (stage 04) so absences compare like "
                       "with like.", examples[:EXAMPLES], table)
    return Finding("annotation_coverage", "Annotations", "Annotation coverage", "ok", summary, table=table)


def check_predicate_maps(ctx: Context) -> Finding:
    fix = "Regenerate predicates with `sharur compute-predicates --db DATASET/sharur.duckdb` (writes; reseal after)."
    if "predicate_provenance" not in ctx.tables:
        return Finding("predicate_maps", "Annotations", "Function labels vs installed maps", "warn",
                       "Function labels carry no map stamp (generated before stamping).", fix)
    try:
        from sharur.predicates.provenance import map_status  # noqa: PLC0415

        status = map_status(_StoreShim(ctx.q))
    except Exception as exc:  # noqa: BLE001 - advisory
        return Finding("predicate_maps", "Annotations", "Function labels vs installed maps", "info",
                       f"Map provenance could not be read ({type(exc).__name__}).")
    if status.state == "current":
        return Finding("predicate_maps", "Annotations", "Function labels vs installed maps", "ok",
                       "Function labels match the installed maps.")
    if status.state == "stale":
        return Finding("predicate_maps", "Annotations", "Function labels vs installed maps", "warn",
                       f"Function labels predate the installed maps ({', '.join(status.changed)} changed).", fix)
    return Finding("predicate_maps", "Annotations", "Function labels vs installed maps", "warn",
                   "Function labels carry no map stamp.", fix)


class _StoreShim:
    """``store.execute`` over a query function, for helpers written against DuckDBStore."""

    def __init__(self, query: Callable[..., list]):
        self._q = query

    def execute(self, sql: str, params: list | None = None) -> list:
        return self._q(sql, params or [])


# --------------------------------------------------------------------------- #
# Callers
# --------------------------------------------------------------------------- #


def check_crispr_scan(ctx: Context) -> Finding:
    reports = {f.name[: -len("_crispr.txt")] for d in ctx.dataset_dir.glob("stage05c*") if d.is_dir()
               for f in d.glob("*_crispr.txt")}
    genomes = set(ctx.genomes)
    scannable = genomes & set(ctx.assemblies)
    loci = ctx.scalar("SELECT COUNT(*) FROM loci WHERE LOWER(locus_type) LIKE '%crispr%'") if "loci" in ctx.tables else 0
    if not reports:
        if loci:
            return Finding("crispr_scan", "Callers", "CRISPR array scan", "info",
                           f"{loci:,} CRISPR arrays are loaded; no MinCED reports sit beside the dataset.")
        return Finding("crispr_scan", "Callers", "CRISPR array scan", "info",
                       "No CRISPR array scan recorded.",
                       "Run stage 05c (MinCED) through `sharur-ingest`, then load the arrays.")
    scanned = reports & genomes
    unscanned = sorted(scannable - reports)
    if unscanned:
        return Finding("crispr_scan", "Callers", "CRISPR array scan", "warn",
                       f"MinCED scanned {len(scanned):,} of {len(scannable):,} genomes with assemblies; "
                       f"{_count(len(unscanned), 'unscanned genome')} can't contribute arrays. {loci:,} arrays loaded.",
                       "Run MinCED on the rest (stage 05c) and load their arrays.", _ex("genome", unscanned))
    return Finding("crispr_scan", "Callers", "CRISPR array scan", "ok",
                   f"MinCED scanned all {len(scanned):,} genomes with assemblies; {loci:,} arrays loaded.")


def check_callers(ctx: Context) -> Finding:
    """Curated system callers: tables present and populated, given the raw hits they need."""
    sources = {s for (s,) in ctx.q("SELECT DISTINCT LOWER(source) FROM annotations")}
    callers = [
        ("defense_systems", "defense systems", "defensefinder",
         "Stage 07 calls defense systems from DefenseFinder hits; rebuild it (`sharur-ingest` resumes)."),
        ("secretion_systems", "secretion systems", "txsscan",
         "Stage 07 calls secretion systems from TXSScan hits; rebuild it (`sharur-ingest` resumes)."),
        ("crispr_cas_systems", "CRISPR-Cas subtypes", None, "Run `sharur cas-type --db DATASET/sharur.duckdb`."),
    ]
    rows, missing = [], []
    for table, label, needs, fix in callers:
        n = ctx.scalar(f'SELECT COUNT(*) FROM "{table}"') if table in ctx.tables else None
        g = ctx.scalar(f'SELECT COUNT(DISTINCT genome_id) FROM "{table}"') if n else 0
        rows.append([label, "absent" if n is None else f"{n:,}", f"{g:,}"])
        if not n and (needs is None or needs in sources):
            missing.append((label, fix))
    table = {"columns": ["Caller", "Calls", "Genomes"], "rows": rows}
    if missing:
        return Finding("callers", "Callers", "Curated system callers", "info",
                       "No calls yet for " + ", ".join(m[0] for m in missing) + "; named systems of those kinds "
                       "need them.", " ".join(m[1] for m in missing), table=table)
    return Finding("callers", "Callers", "Curated system callers", "ok",
                   "Every curated caller whose inputs exist has calls.", table=table)


# --------------------------------------------------------------------------- #
# References
# --------------------------------------------------------------------------- #


def check_kegg(ctx: Context) -> Finding:
    try:
        from sharur.predicates.mappings.kegg_map import default_kegg_dir  # noqa: PLC0415

        kegg = default_kegg_dir()
    except Exception:  # noqa: BLE001
        kegg = None
    if kegg and (Path(kegg) / "kegg_modules.tsv").is_file():
        return Finding("kegg", "References", "Local KEGG build", "ok",
                       f"KEGG modules and KO names available ({kegg}).")
    return Finding("kegg", "References", "Local KEGG build", "info",
                   "No local KEGG build: module completeness and KO names are unavailable.",
                   "Build it with `sharur setup-kegg` (KEGG data stays on this machine).")


CHECKS: list[Callable[[Context], Finding]] = [
    check_schema, check_seal,
    check_completeness, check_contamination, check_empty_genomes,
    check_gene_calls, check_strand, check_positionless, check_duplicate_coordinates, check_contigs,
    check_annotation_coverage, check_predicate_maps,
    check_crispr_scan, check_callers,
    check_kegg,
]


def run_checks(db_path: str | Path, *, query: Callable[..., list] | None = None,
               assemblies: dict[str, Path] | None = None) -> HealthReport:
    """Every check against one dataset. Opens the database read-only unless ``query`` is given."""
    started = time.perf_counter()
    store = None
    if query is None:
        from sharur.storage.duckdb_store import DuckDBStore  # noqa: PLC0415

        store = DuckDBStore(str(db_path), read_only=True)

        def query(sql: str, params: list | None = None) -> list:  # noqa: E306
            return store.execute(sql, params or [])
    try:
        ctx = Context(query, Path(db_path), assemblies)
        findings = []
        for check in CHECKS:
            try:
                findings.append(check(ctx))
            except Exception as exc:  # noqa: BLE001 - one broken check must not hide the rest
                findings.append(Finding(check.__name__.removeprefix("check_"), "Database",
                                        check.__name__.removeprefix("check_").replace("_", " ").capitalize(),
                                        "info", f"Check could not run ({type(exc).__name__}: {exc})."))
    finally:
        if store is not None:
            store.close()
    return HealthReport(str(db_path), findings, time.perf_counter() - started)


__all__ = ["CHECKS", "Finding", "HealthReport", "run_checks"]
