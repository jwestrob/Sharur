"""CRISPR arrays: catalog, per-array page with repeats and spacers, genomic context.

Array calls (``loci`` rows with ``locus_type = 'crispr'``) carry their span and
repeat consensus. Repeats and spacers come from MinCED's own text report when
Stage 05c kept it beside the dataset. Without it, they are recovered from the
assembly:
the consensus is matched along the array span (Hamming distance up to about a
tenth of its length), diverged copies (up to 30%) are then sought at the
expected spacing inside the span, and spacers are the stretches between
consecutive repeats. On DPANN (474 arrays), recovered repeats span at least 95%
of the recorded array in 440 arrays.
Nearby Cas-domain genes come from the evidence-backed labels and from domain
names; these are observed domains, named here as "Cas-domain genes".
"""

from __future__ import annotations

import gzip
import hashlib
import json
import re
import threading
from collections import Counter, OrderedDict
from pathlib import Path
from typing import Any

from fastapi import HTTPException, Query, Request
from fastapi.responses import HTMLResponse, PlainTextResponse
from markupsafe import Markup

from sharur.browser import charts
from sharur.crispr import annotate_repeats, parse_minced_text

CONTEXT_BP = 15000
CAS_LABELS = {"cas_domain", "crispr_associated", "cas_nuclease", "crispr_adaptation", "crispr_accessory",
              "crispr_class1", "crispr_class2"}
_CAS_NAME = re.compile(r"^(Cas|cas|Csm|Cmr|Csx|Csa|Csb|Csc|Csd|Cse|Csf|Csn|Cpf|Csy|Cst|Csh|Cmx|DinG_cas|CRISPR)")
_SUFFIXES = (".fna", ".fa", ".fasta", ".fas", ".fna.gz", ".fa.gz", ".fasta.gz")


# --------------------------------------------------------------------------- #
# Assemblies
# --------------------------------------------------------------------------- #


class Assemblies:
    """Find a genome's assembly FASTA and read single contigs (small LRU cache)."""

    def __init__(self, dataset_dir: Path, extra: list[Path] | None = None):
        self.paths: dict[str, Path] = {}
        dirs = list(extra or []) + [dataset_dir / "stage00_prepared" / "genomes", dataset_dir / "genomes_fna",
                                    dataset_dir / "genomes_fna_new", dataset_dir / "source", dataset_dir / "assemblies"]
        manifest = dataset_dir / "stage00_prepared" / "processing_manifest.json"
        if manifest.is_file():
            try:
                for entry in json.loads(manifest.read_text()).get("genomes", []):
                    if entry.get("genome_id") and entry.get("output_path") and Path(entry["output_path"]).is_file():
                        self.paths.setdefault(entry["genome_id"], Path(entry["output_path"]))
            except (ValueError, OSError):
                pass
        for d in dirs:
            if not d.is_dir():
                continue
            for f in d.iterdir():
                for suffix in _SUFFIXES:
                    if f.name.endswith(suffix):
                        self.paths.setdefault(f.name[: -len(suffix)], f)
                        break
        self._cache: OrderedDict[tuple[str, str], str] = OrderedDict()
        self._lock = threading.Lock()

    def contig(self, bin_id: str, contig_id: str) -> str | None:
        key = (bin_id, contig_id)
        with self._lock:
            if key in self._cache:
                self._cache.move_to_end(key)
                return self._cache[key]
        path = self.paths.get(bin_id)
        if path is None:
            return None
        opener = gzip.open if path.suffix == ".gz" else open
        wanted = {contig_id, contig_id.rsplit("|", 1)[-1]}
        chunks: list[str] = []
        found = False
        with opener(path, "rt") as handle:
            for line in handle:
                if line.startswith(">"):
                    if found:
                        break
                    found = line[1:].split(None, 1)[0] in wanted if line[1:].strip() else False
                elif found:
                    chunks.append(line.strip())
        if not found:
            return None
        seq = "".join(chunks).upper()
        with self._lock:
            self._cache[key] = seq
            while len(self._cache) > 48:
                self._cache.popitem(last=False)
        return seq


# --------------------------------------------------------------------------- #
# Repeats and spacers
# --------------------------------------------------------------------------- #


def _best_window(seq: str, consensus: str, lo: int, hi: int) -> tuple[int, int] | None:
    """(mismatches, 0-based start) of the closest window to ``consensus`` starting in [lo, hi]."""
    n = len(consensus)
    best = None
    for i in range(max(0, lo), min(len(seq) - n, hi) + 1):
        mism = sum(1 for a, b in zip(seq[i:i + n], consensus) if a != b)
        if best is None or mism < best[0]:
            best = (mism, i)
    return best


def find_repeats(seq: str, consensus: str, start: int, end: int) -> list[dict[str, Any]]:
    """Repeat copies of ``consensus`` within [start, end] (1-based).

    First pass: copies within about a tenth of the length of the consensus,
    greedy and non-overlapping. Second pass, inside the recorded span only:
    diverged copies (up to 30% mismatches) at the expected spacing beyond the
    outermost copies and inside gaps much longer than the typical spacer;
    arrays commonly end in such a degenerate repeat.
    """
    consensus = consensus.upper()
    n = len(consensus)
    if not n or not seq:
        return []
    budget, loose = max(2, n // 10), max(4, int(0.3 * n))
    lo, hi = max(0, start - 1 - n), min(len(seq), end + n)
    hits = []
    for i in range(lo, hi - n + 1):
        mism = sum(1 for a, b in zip(seq[i:i + n], consensus) if a != b)
        if mism <= budget:
            hits.append((mism, i))
    chosen: list[tuple[int, int]] = []
    for mism, i in sorted(hits):
        if all(abs(i - j) >= n for _, j in chosen):
            chosen.append((mism, i))
    chosen.sort(key=lambda h: h[1])
    if chosen:
        gaps = [b - a - n for (_, a), (_, b) in zip(chosen, chosen[1:])]
        spacer = sorted(gaps)[len(gaps) // 2] if gaps else 35
        span_lo, span_hi = start - 1, end - n  # 0-based limits for a repeat start inside the span
        added = True
        while added:  # extend outward while diverged copies keep turning up at the expected spacing
            added = False
            first, last = chosen[0][1], chosen[-1][1]
            for target, bound in ((first - n - spacer, first - n - 15), (last + n + spacer, None)):
                a_lo, a_hi = target - 12, target + 12
                if bound is not None:
                    a_hi = min(a_hi, bound)
                a_lo, a_hi = max(a_lo, span_lo), min(a_hi, span_hi)
                if a_lo > a_hi:
                    continue
                best = _best_window(seq, consensus, a_lo, a_hi)
                if best and best[0] <= loose and all(abs(best[1] - j) >= n for _, j in chosen):
                    chosen.append(best)
                    chosen.sort(key=lambda h: h[1])
                    added = True
        filled = []
        for (_, a), (_, b) in zip(chosen, chosen[1:]):  # long internal gaps: look for a diverged copy midway
            if b - a - n > 1.8 * spacer + n:
                best = _best_window(seq, consensus, a + n + spacer - 12, a + n + spacer + 12)
                if best and best[0] <= loose and best[1] + n <= b:
                    filled.append(best)
        chosen = sorted(chosen + filled, key=lambda h: h[1])
    out = []
    for mism, i in chosen:
        window = seq[i:i + n]
        out.append({"start": i + 1, "end": i + n, "seq": window, "mismatches": mism, "diverged": mism > budget,
                    "diff": [k for k, (a, b) in enumerate(zip(window, consensus)) if a != b]})
    return out


def spacers_between(seq: str, repeats: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for a, b in zip(repeats, repeats[1:]):
        s, e = a["end"] + 1, b["start"] - 1
        if e >= s:
            out.append({"start": s, "end": e, "seq": seq[s - 1:e], "length": e - s + 1})
    return out


def _spacer_color(text: str) -> str:
    return charts.PALETTE[int(hashlib.md5(text.encode()).hexdigest()[:6], 16) % len(charts.PALETTE)]


def array_svg(repeats: list[dict[str, Any]], spacers: list[dict[str, Any]], start: int, end: int,
              width: int = 1100) -> Markup:
    """Repeats as diamonds, spacers as colored blocks, numbered every fifth spacer."""
    if not repeats:
        return Markup("")
    lo, hi = repeats[0]["start"], repeats[-1]["end"]
    span = max(1, hi - lo)
    pad = 16
    scale = (width - 2 * pad) / span
    parts = [f'<line x1="{pad}" y1="30" x2="{width - pad}" y2="30" class="backbone-line"/>']
    for k, s in enumerate(spacers, 1):
        x1, x2 = pad + (s["start"] - lo) * scale, pad + (s["end"] - lo) * scale
        parts.append(f'<g><title>Spacer {k} · {s["length"]} bp · {s["start"]:,}–{s["end"]:,}</title>'
                     f'<rect x="{x1:.1f}" y="21" width="{max(x2 - x1, 1.5):.1f}" height="18" rx="3" '
                     f'fill="{_spacer_color(s["seq"])}"/></g>')
        if k == 1 or k % 5 == 0:
            parts.append(f'<text x="{(x1 + x2) / 2:.1f}" y="54" class="axis" text-anchor="middle">{k}</text>')
    for r in repeats:
        cx = pad + ((r["start"] + r["end"]) / 2 - lo) * scale
        w = max(3.0, (r["end"] - r["start"]) * scale / 2)
        cls = "repeat" + (" variant" if r["mismatches"] else "")
        parts.append(f'<g><title>Repeat {r["start"]:,}–{r["end"]:,}'
                     f'{" · " + str(r["mismatches"]) + " mismatches" if r["mismatches"] else ""}</title>'
                     f'<polygon class="{cls}" points="{cx - w:.1f},30 {cx:.1f},18 {cx + w:.1f},30 {cx:.1f},42"/></g>')
    parts.append(f'<text x="{pad}" y="12" class="axis">{lo:,}</text>'
                 f'<text x="{width - pad}" y="12" class="axis" text-anchor="end">{hi:,}</text>')
    return Markup(f'<svg class="crispr-array" viewBox="0 0 {width} 60" role="img" aria-label="CRISPR array">'
                  f'{"".join(parts)}</svg>')


def context_svg(genes: list[dict[str, Any]], arrays: list[dict[str, Any]], start: int, end: int,
                width: int = 1100) -> Markup:
    """Genes around the array on both strands, Cas-domain genes highlighted, arrays as striped blocks."""
    span = max(1, end - start)
    pad = 12
    scale = (width - 2 * pad) / span
    parts = [f'<line x1="0" y1="44" x2="{width}" y2="44" class="backbone-line"/>']
    step = next(s for s in (1000, 2000, 5000, 10000, 20000) if span / s <= 12)
    for pos in range(((start + step - 1) // step) * step, end + 1, step):
        x = pad + (pos - start) * scale
        parts.append(f'<line x1="{x:.1f}" y1="10" x2="{x:.1f}" y2="78" class="grid-line"/>'
                     f'<text x="{x + 3:.1f}" y="10" class="axis">{pos / 1000:,.0f} kb</text>')
    for a in arrays:
        x1, x2 = pad + (max(a["start"], start) - start) * scale, pad + (min(a["end"], end) - start) * scale
        parts.append(f'<a href="/crispr/{charts.quote(a["locus_id"], safe="")}"><g><title>CRISPR array '
                     f'{a["start"]:,}–{a["end"]:,}</title><rect x="{x1:.1f}" y="32" width="{max(x2 - x1, 4):.1f}" '
                     f'height="24" rx="4" class="array-block{" current" if a.get("current") else ""}"/></g></a>')
    for g in genes:
        x1, x2 = pad + (max(g["start"], start) - start) * scale, pad + (min(g["end"], end) - start) * scale
        y, h = (26 if g["strand"] != "-" else 62), 9
        head = min(8.0, (x2 - x1) * 0.45)
        if g["strand"] == "-":
            pts = f"{x2:.1f},{y - h} {x1 + head:.1f},{y - h} {x1:.1f},{y} {x1 + head:.1f},{y + h} {x2:.1f},{y + h}"
        else:
            pts = f"{x1:.1f},{y - h} {x2 - head:.1f},{y - h} {x2:.1f},{y} {x2 - head:.1f},{y + h} {x1:.1f},{y + h}"
        cls = "gene cas" if g["cas"] else ("gene" if g["label"] else "gene dark")
        label = g["label"] or "no annotation"
        text = ""
        short = label.split(" (")[0]
        if x2 - x1 > 6.0 * len(short) + 10:
            text = f'<text x="{(x1 + x2) / 2:.1f}" y="{y + 4}" class="gene-label">{charts._e(short)}</text>'
        parts.append(f'<a href="/protein/{charts.quote(g["protein_id"], safe="")}"><g><title>{charts._e(label)}'
                     f'{" · Cas-domain gene" if g["cas"] else ""}</title><polygon class="{cls}" points="{pts}"/>'
                     f'{text}</g></a>')
    return Markup(f'<svg class="crispr-context" viewBox="0 0 {width} 84" role="img" aria-label="Genomic context">'
                  f'{"".join(parts)}</svg>')


# --------------------------------------------------------------------------- #
# Routes
# --------------------------------------------------------------------------- #


def register(app, ctx, assembly_dirs: list[Path] | None = None) -> None:
    catalog = ctx.catalog
    dataset_dir = ctx.db_path.resolve().parent
    assemblies = Assemblies(dataset_dir, assembly_dirs)
    ctx.assemblies = assemblies
    # MinCED text reports (exact repeat and spacer calls) from Stage 05c, when kept beside the dataset
    reports = {}
    for d in sorted(dataset_dir.glob("stage05c*")):
        for f in d.glob("*_crispr.txt"):
            reports.setdefault(f.name[: -len("_crispr.txt")], f)
    parsed: dict[str, list[dict[str, Any]]] = {}

    def minced_array(bin_id: str, contig: str, start: int) -> dict[str, Any] | None:
        if bin_id not in reports:
            return None
        if bin_id not in parsed:
            parsed[bin_id] = parse_minced_text(reports[bin_id])
        names = {contig, contig.rsplit("|", 1)[-1]}
        return next((a for a in parsed[bin_id] if a["contig"] in names and a["start"] == start), None)

    def arrays() -> list[dict[str, Any]]:
        if not hasattr(ctx, "_crispr_arrays"):
            with ctx.lock:
                rows = ctx.store.execute(
                    """SELECT l.locus_id, l.contig_id, l.start, l.end_coord, l.metadata, c.bin_id, c.length
                       FROM loci l LEFT JOIN contigs c USING (contig_id)
                       WHERE LOWER(l.locus_type) LIKE '%crispr%' ORDER BY c.bin_id, l.contig_id, l.start""")
            out = []
            for locus_id, contig, start, end, metadata, bin_id, contig_length in rows:
                consensus, stored = "", None
                try:
                    meta = json.loads(metadata) if isinstance(metadata, str) else (metadata or {})
                    inner = meta.get("metadata", meta)
                    consensus = inner.get("rpt_unit_seq") or inner.get("repeat_consensus") or ""
                    if meta.get("repeats"):
                        stored = {"repeats": meta["repeats"], "spacers": meta.get("spacers", [])}
                except (ValueError, AttributeError):
                    pass
                g = catalog.by_bin.get(bin_id)
                out.append({"stored": stored, "locus_id": locus_id, "contig_id": contig, "start": start, "end": end,
                            "length": (end or 0) - (start or 0) + 1, "consensus": consensus.upper(),
                            "bin_id": bin_id, "lineage": g.label if g else "", "contig_length": contig_length})
            ctx._crispr_arrays = out
            ctx._crispr_by_id = {a["locus_id"]: a for a in out}
        return ctx._crispr_arrays

    def genes_near(contig: str, lo: int, hi: int) -> list[dict[str, Any]]:
        with ctx.lock:
            rows = ctx.store.execute(
                """WITH p AS (SELECT protein_id, start, end_coord, strand, sequence_length FROM proteins
                              WHERE contig_id = ? AND end_coord >= ? AND start <= ?),
                        best AS (SELECT a.protein_id, ARG_MIN(COALESCE(NULLIF(a.name, ''), a.accession),
                                                              COALESCE(a.evalue, 1)) AS top,
                                        LIST(DISTINCT COALESCE(NULLIF(a.name, ''), a.accession)) AS names
                                 FROM annotations a JOIN p USING (protein_id)
                                 WHERE LOWER(a.source) NOT IN ('defensefinder_system', 'txsscan_system')
                                 GROUP BY 1)
                   SELECT p.protein_id, p.start, p.end_coord, p.strand, p.sequence_length, best.top, best.names,
                          pp.predicates
                   FROM p LEFT JOIN best USING (protein_id) LEFT JOIN protein_predicates pp USING (protein_id)
                   ORDER BY p.start""", [contig, lo, hi])
        out = []
        for pid, s, e, strand, n, top, names, predicates in rows:
            labels = set(predicates or [])
            cas_labels = sorted(labels & CAS_LABELS | {p for p in labels if p.startswith("crispr_type")})
            cas_names = sorted({x for x in (names or []) if _CAS_NAME.match(str(x))})
            out.append({"protein_id": pid, "start": s, "end": e, "strand": strand, "length_aa": n,
                        "label": ctx.describe_hit(top), "cas": bool(cas_labels or cas_names),
                        "cas_labels": cas_labels, "cas_names": cas_names})
        return out

    @app.get("/crispr", response_class=HTMLResponse)
    def crispr_catalog(request: Request, rank: str = Query("class")):
        rows = arrays()
        carriers = {catalog.by_bin[a["bin_id"]].index: 1 for a in rows if a["bin_id"] in catalog.by_bin}
        per_genome = Counter(a["bin_id"] for a in rows)
        consensus = Counter(a["consensus"] for a in rows if a["consensus"])
        return ctx.render(request, "crispr_catalog.html", "systems", arrays=rows, genomes=len(per_genome),
                          per_genome=per_genome.most_common(10), shared=consensus.most_common(8),
                          prevalence=[r for r in catalog.prevalence_by(carriers, rank) if r["genomes"] >= 3][:25],
                          rank=rank, assemblies=len(assemblies.paths))

    @app.get("/crispr/{locus_id:path}/spacers.fasta")
    def spacers_fasta(locus_id: str):
        detail = _detail(locus_id)
        body = "".join(f">{locus_id}_spacer{k} {detail['array']['contig_id']}:{s['start']}-{s['end']}\n{s['seq']}\n"
                       for k, s in enumerate(detail["spacers"], 1))
        name = re.sub(r"[^\w.-]+", "_", locus_id)[:100]
        return PlainTextResponse(body, headers={"Content-Disposition": f'attachment; filename="{name}_spacers.fna"'})

    def _detail(locus_id: str) -> dict[str, Any]:
        arrays()
        a = ctx._crispr_by_id.get(locus_id)
        if a is None:
            raise HTTPException(404, "CRISPR array not found")
        called = None
        if a.get("stored"):
            called = {"repeats": [dict(r) for r in a["stored"]["repeats"]], "spacers": a["stored"]["spacers"]}
        elif a["bin_id"]:
            called = minced_array(a["bin_id"], a["contig_id"], a["start"])
        if called and called["repeats"]:
            annotate_repeats(called, a["consensus"] or None)
            return {"array": a, "seq_available": True, "repeats": called["repeats"], "spacers": called["spacers"],
                    "contig_length": a["contig_length"], "source": "MinCED"}
        seq = assemblies.contig(a["bin_id"], a["contig_id"]) if a["bin_id"] else None
        repeats = find_repeats(seq, a["consensus"], a["start"], a["end"]) if seq and a["consensus"] else []
        return {"array": a, "seq_available": seq is not None, "repeats": repeats,
                "spacers": spacers_between(seq, repeats) if seq else [],
                "contig_length": len(seq) if seq else a["contig_length"], "source": "sequence scan"}

    @app.get("/crispr/{locus_id:path}", response_class=HTMLResponse)
    def crispr_page(request: Request, locus_id: str, flank: int = Query(CONTEXT_BP, ge=2000, le=60000)):
        detail = _detail(locus_id)
        a = detail["array"]
        lo, hi = max(1, a["start"] - flank), a["end"] + flank
        if detail["contig_length"]:
            hi = min(hi, detail["contig_length"])
        genes = genes_near(a["contig_id"], lo, hi)
        neighbors = [dict(x, current=x["locus_id"] == locus_id) for x in arrays()
                     if x["contig_id"] == a["contig_id"] and x["end"] >= lo and x["start"] <= hi]
        cas = [g for g in genes if g["cas"]]
        nearest = min((min(abs(g["start"] - a["end"]), abs(a["start"] - g["end"])) for g in cas), default=None)
        spacer_lengths = [s["length"] for s in detail["spacers"]]
        g = catalog.by_bin.get(a["bin_id"])
        duplicates = Counter(s["seq"] for s in detail["spacers"])
        return ctx.render(
            request, "crispr.html", "systems", a=a, g=g, detail=detail, genes=genes, cas=cas, nearest=nearest,
            context=context_svg(genes, neighbors, lo, hi), array=array_svg(detail["repeats"], detail["spacers"],
                                                                           a["start"], a["end"]),
            lo=lo, hi=hi, flank=flank, spacer_lengths=spacer_lengths,
            repeated_spacers=sum(1 for n in duplicates.values() if n > 1),
            at_contig_start=lo <= 1, at_contig_end=bool(detail["contig_length"]) and hi >= detail["contig_length"],
            spacer_color=_spacer_color)

    # ------------------------------------------------------------------ #
    # CRISPR-Cas loci: arrays and Cas-domain genes merged along each contig
    # ------------------------------------------------------------------ #

    def cas_loci() -> list[dict[str, Any]]:
        if hasattr(ctx, "_cas_loci"):
            return ctx._cas_loci
        with ctx.lock:
            rows = ctx.store.execute(
                """WITH by_label AS (
                       SELECT protein_id FROM protein_predicates
                       WHERE list_has_any(predicates, ['cas_domain', 'crispr_associated', 'cas_nuclease',
                                                       'crispr_adaptation', 'crispr_accessory'])),
                        by_name AS (
                       SELECT DISTINCT protein_id FROM annotations
                       WHERE regexp_matches(name, ?) OR description ILIKE '%CRISPR-associated%'),
                        cas AS (SELECT protein_id FROM by_label UNION SELECT protein_id FROM by_name),
                        named AS (
                       SELECT a.protein_id, ARG_MIN(COALESCE(NULLIF(a.name, ''), a.accession), COALESCE(a.evalue, 1)) AS top
                       FROM annotations a JOIN cas USING (protein_id)
                       WHERE regexp_matches(a.name, ?) OR a.description ILIKE '%CRISPR-associated%' GROUP BY 1)
                   SELECT p.protein_id, p.contig_id, p.bin_id, p.start, p.end_coord, p.strand, named.top
                   FROM cas JOIN proteins p USING (protein_id) LEFT JOIN named USING (protein_id)
                   ORDER BY p.contig_id, p.start""", [_CAS_NAME.pattern, _CAS_NAME.pattern])
        features: dict[str, list[dict[str, Any]]] = {}
        for pid, contig, bin_id, start, end, strand, top in rows:
            features.setdefault(contig, []).append({"kind": "cas", "protein_id": pid, "bin_id": bin_id,
                                                    "start": start, "end": end, "strand": strand, "name": top})
        for a in arrays():
            features.setdefault(a["contig_id"], []).append({"kind": "array", "locus_id": a["locus_id"],
                                                            "bin_id": a["bin_id"], "start": a["start"],
                                                            "end": a["end"]})
        loci = []
        for contig, items in features.items():
            items.sort(key=lambda f: f["start"])
            group: list[dict[str, Any]] = []
            for f in items + [None]:
                if f is not None and (not group or f["start"] - max(g["end"] for g in group) <= 10000):
                    group.append(f)
                    continue
                if group and (any(g["kind"] == "array" for g in group) or
                              sum(g["kind"] == "cas" for g in group) >= 2):
                    # a lone Cas-domain gene without an array is a domain hit, not a locus
                    cas = [g for g in group if g["kind"] == "cas"]
                    arr = [g for g in group if g["kind"] == "array"]
                    kind = "array + Cas" if cas and arr else ("Cas genes only" if cas else "array only")
                    bin_id = group[0]["bin_id"]
                    g0 = catalog.by_bin.get(bin_id)
                    loci.append({"id": f"{contig}:{group[0]['start']}", "contig_id": contig, "bin_id": bin_id,
                                 "lineage": g0.label if g0 else "", "start": min(g["start"] for g in group),
                                 "end": max(g["end"] for g in group), "cas": cas, "arrays": arr, "kind": kind})
                group = [f] if f is not None else []
        loci.sort(key=lambda l: (-len(l["cas"]) - 3 * len(l["arrays"]), l["bin_id"] or "", l["start"]))
        ctx._cas_loci = loci
        return loci

    def locus_row_svg(locus: dict[str, Any], genes: list[dict[str, Any]], lo: int, hi: int, flip: bool,
                      width: int = 1040) -> Markup:
        span = max(1, hi - lo)
        pad = 12
        scale = (width - 2 * pad) / span

        def x(pos: int) -> float:
            return pad + ((hi - pos) if flip else (pos - lo)) * scale

        parts = [f'<line x1="0" y1="22" x2="{width}" y2="22" class="backbone-line"/>']
        for a in locus["arrays"]:
            x1, x2 = sorted((x(max(a["start"], lo)), x(min(a["end"], hi))))
            parts.append(f'<a href="/crispr/{charts.quote(a["locus_id"], safe="")}"><g><title>CRISPR array '
                         f'{a["start"]:,}–{a["end"]:,}</title><rect x="{x1:.1f}" y="10" width="{max(x2 - x1, 4):.1f}" '
                         f'height="24" rx="4" class="array-block current"/></g></a>')
        for g in genes:
            x1, x2 = sorted((x(max(g["start"], lo)), x(min(g["end"], hi))))
            forward = (g["strand"] != "-") != flip
            y, h = 22, 10 if g["cas"] else 7
            head = min(8.0, (x2 - x1) * 0.45)
            if forward:
                pts = f"{x1:.1f},{y - h} {x2 - head:.1f},{y - h} {x2:.1f},{y} {x2 - head:.1f},{y + h} {x1:.1f},{y + h}"
            else:
                pts = f"{x2:.1f},{y - h} {x1 + head:.1f},{y - h} {x1:.1f},{y} {x1 + head:.1f},{y + h} {x2:.1f},{y + h}"
            cls = "gene cas" if g["cas"] else ("gene" if g["label"] else "gene dark")
            short = (g["cas_name"] or g["label"] or "").split(" (")[0]
            text = (f'<text x="{(x1 + x2) / 2:.1f}" y="{y + 4}" class="gene-label">{charts._e(short)}</text>'
                    if short and x2 - x1 > 6.0 * len(short) + 10 else "")
            parts.append(f'<a href="/protein/{charts.quote(g["protein_id"], safe="")}"><g><title>'
                         f'{charts._e(short or "no annotation")}</title><polygon class="{cls}" points="{pts}"/>{text}</g></a>')
        return Markup(f'<svg class="stack-row" viewBox="0 0 {width} 44" role="img" aria-label="CRISPR-Cas locus">'
                      f'{"".join(parts)}</svg>')

    @app.get("/crispr-cas", response_class=HTMLResponse)
    def crispr_cas(request: Request, kind: str = Query(""), page: int = Query(1, ge=1),
                   flank: int = Query(3000, ge=0, le=20000)):
        loci = cas_loci()
        counts = Counter(l["kind"] for l in loci)
        selected = [l for l in loci if not kind or l["kind"] == kind]
        per_page = 25
        pages = max(1, (len(selected) + per_page - 1) // per_page)
        page = min(page, pages)
        shown = selected[(page - 1) * per_page: page * per_page]
        windows = [(l, max(1, l["start"] - flank), l["end"] + flank) for l in shown]
        span = max((hi - lo for _, lo, hi in windows), default=1)
        rows = []
        for l, lo, hi in windows:
            genes = genes_near(l["contig_id"], lo, hi)
            names = {c["protein_id"]: c["name"] for c in l["cas"]}
            for g in genes:
                g["cas_name"] = names.get(g["protein_id"])
                g["cas"] = g["cas"] or g["protein_id"] in names
            # orient on the Cas1 gene when present (or the first Cas gene); shared scale = widest window
            key = next((c for c in l["cas"] if re.search(r"cas1(?![0-9])", str(c["name"] or ""), re.I)),
                       l["cas"][0] if l["cas"] else None)
            flip = bool(key and key["strand"] == "-")
            center = (lo + hi) // 2
            rows.append((l, locus_row_svg(l, genes, center - span // 2, center + span // 2, flip), flip))
        return ctx.render(request, "crispr_cas.html", "systems", rows=rows, kind=kind, counts=counts,
                          total=len(selected), page=page, pages=pages, flank=flank,
                          genomes=len({l["bin_id"] for l in loci}), arrays_total=len(arrays()))

    ctx.crispr_arrays = arrays
    ctx.cas_loci = cas_loci
