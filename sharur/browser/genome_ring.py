"""A genome as a ring: contigs around the circle, genes by strand and function, and marked features.

Contigs sit end to end around the circle in order of length, at their recorded
lengths. Genes on the forward and reverse strands take their function
category's colour (unannotated genes stay grey). The inner ring marks curated
system members, CRISPR arrays and CRISPR-Cas calls, prophage and island
regions, and proteins of at least 3,000 aa. The SVG links each mark to its
page and carries a title for hovering.
"""

from __future__ import annotations

import math
from html import escape
from typing import Any
from urllib.parse import quote

from markupsafe import Markup

from sharur.browser.catalog import CATEGORY_LABELS
from sharur.browser.charts import CATEGORY_COLORS, OVERLAY_COLORS

SIZE = 720
CENTER = SIZE / 2
R_CONTIG = (306, 318)
R_FORWARD = (290, 302)
R_REVERSE = (274, 286)
R_FEATURE = (250, 268)
GIANT_AA = 3000
MARK_COLORS = {**OVERLAY_COLORS, "giant": "#e0895c", "cas": "#45aec4"}
MARK_LABELS = {"defense": "Defense system", "secretion": "Secretion system", "crispr": "CRISPR array",
               "cas": "CRISPR-Cas call", "prophage": "Prophage region", "island": "Island",
               "giant": f"Protein ≥ {GIANT_AA:,} aa"}


def _url(*parts: str) -> str:
    return "/" + "/".join(quote(str(p), safe="") for p in parts)


def _polar(angle: float, r: float) -> tuple[float, float]:
    return CENTER + r * math.sin(angle), CENTER - r * math.cos(angle)


def _arc(a0: float, a1: float, r0: float, r1: float) -> str:
    """Annular sector path from angle a0 to a1 (radians, clockwise from 12 o'clock)."""
    a1 = max(a1, a0 + 1e-4)
    large = 1 if a1 - a0 > math.pi else 0
    x0, y0 = _polar(a0, r1)
    x1, y1 = _polar(a1, r1)
    x2, y2 = _polar(a1, r0)
    x3, y3 = _polar(a0, r0)
    return (f"M{x0:.2f},{y0:.2f} A{r1},{r1} 0 {large} 1 {x1:.2f},{y1:.2f} "
            f"L{x2:.2f},{y2:.2f} A{r0},{r0} 0 {large} 0 {x3:.2f},{y3:.2f} Z")


def ring_data(store, catalog, bin_id: str, *, top_categories, cctyper_calls=()) -> dict[str, Any] | None:
    contigs = store.execute("SELECT contig_id, length FROM contigs WHERE bin_id = ? ORDER BY length DESC, contig_id",
                            [bin_id])
    contigs = [(c, n) for c, n in contigs if n and n > 0]
    if not contigs:
        return None
    genes = store.execute(
        """SELECT protein_id, contig_id, start, end_coord, strand, sequence_length FROM proteins
           WHERE bin_id = ? AND contig_id <> protein_id""", [bin_id])
    categories = top_categories([g[0] for g in genes])
    where = {pid: (contig, start, end) for pid, contig, start, end, *_ in genes}
    marks: list[dict[str, Any]] = []
    for s in catalog.systems:
        if s["bin_id"] != bin_id:
            continue
        spans = [where[p] for p in s["proteins"] if p in where]
        if spans:
            contig = spans[0][0]
            marks.append({"kind": s["kind"], "contig": contig, "start": min(x[1] for x in spans if x[0] == contig),
                          "end": max(x[2] for x in spans if x[0] == contig), "label": f"{s['type']} ({s['kind']})",
                          "url": f"/call/{quote(s['system_id'], safe='')}"})
    for locus in catalog.loci:
        if locus["bin_id"] == bin_id and locus["type"] in ("crispr", "prophage", "island"):
            kind = locus["type"]
            url = _url("crispr", locus["locus_id"]) if kind == "crispr" else \
                f"{_url('contig', locus['contig_id'])}?start={max(1, locus['start'] - 2000)}&span={locus['end'] - locus['start'] + 4000}"
            marks.append({"kind": kind, "contig": locus["contig_id"], "start": locus["start"], "end": locus["end"],
                          "label": MARK_LABELS[kind], "url": url})
    for call in cctyper_calls:
        if call.get("confident"):
            marks.append({"kind": "cas", "contig": call["contig_id"], "start": call["start"], "end": call["end"],
                          "label": f"CRISPR-Cas {call['prediction']}", "url": f"/cas-system/{quote(call['system_id'], safe='')}"})
    for pid, contig, start, end, _strand, length in genes:
        if length and length >= GIANT_AA:
            marks.append({"kind": "giant", "contig": contig, "start": start, "end": end,
                          "label": f"{pid} · {length:,} aa", "url": _url("protein", pid)})
    return {"contigs": contigs, "genes": genes, "categories": categories, "marks": marks}


def render(data: dict[str, Any], *, title: str = "", subtitle: str = "") -> Markup:
    contigs = data["contigs"]
    total = sum(n for _, n in contigs)
    n = len(contigs)
    gap = 0.0 if n > 240 else min(0.012, 1.2 / max(n, 1))
    usable = 2 * math.pi - gap * n
    offset: dict[str, tuple[float, float]] = {}   # contig -> (start angle, radians per bp)
    angle = 0.0
    parts = []
    for i, (contig, length) in enumerate(contigs):
        span = usable * length / total
        offset[contig] = (angle, span / length)
        shade = "ring-contig" if i % 2 == 0 else "ring-contig alt"
        parts.append(f'<a href="{_url("contig", contig)}"><path d="{_arc(angle, angle + span, *R_CONTIG)}" '
                     f'class="{shade}"><title>{escape(contig)} · {length:,} bp</title></path></a>')
        angle += span + gap

    def at(contig: str, pos: int) -> float | None:
        if contig not in offset:
            return None
        a0, per_bp = offset[contig]
        return a0 + max(0, pos) * per_bp

    cats = data["categories"]
    used_categories = set()
    for pid, contig, start, end, strand, _length in data["genes"]:
        a0, a1 = at(contig, start), at(contig, end)
        if a0 is None:
            continue
        a1 = max(a1, a0 + 0.0025)
        category = cats.get(pid)
        if category:
            used_categories.add(category)
        fill = CATEGORY_COLORS.get(category, "")
        ring = R_FORWARD if strand in ("+", "1") else R_REVERSE
        style = f' style="fill:{fill}"' if fill else ""
        parts.append(f'<path d="{_arc(a0, a1, *ring)}" class="ring-gene"{style}/>')
    used_marks = set()
    for m in data["marks"]:
        a0, a1 = at(m["contig"], m["start"]), at(m["contig"], m["end"])
        if a0 is None:
            continue
        used_marks.add(m["kind"])
        a1 = max(a1, a0 + 0.012)
        parts.append(f'<a href="{escape(m["url"])}"><path d="{_arc(a0, a1, *R_FEATURE)}" class="ring-mark" '
                     f'style="fill:{MARK_COLORS[m["kind"]]}"><title>{escape(m["label"])}</title></path></a>')
    parts.append(f'<circle cx="{CENTER}" cy="{CENTER}" r="{R_REVERSE[0] - 1}" class="ring-guide"/>')
    parts.append(f'<text x="{CENTER}" y="{CENTER - 6}" class="ring-title">{escape(title)}</text>')
    if subtitle:
        parts.append(f'<text x="{CENTER}" y="{CENTER + 16}" class="ring-sub">{escape(subtitle)}</text>')
    svg = (f'<svg class="genome-ring" viewBox="0 0 {SIZE} {SIZE}" role="img" aria-label="Genome ring">'
           + "".join(parts) + "</svg>")
    legend = "".join(f'<span><i style="background:{CATEGORY_COLORS[c]}"></i>{escape(CATEGORY_LABELS.get(c, c))}</span>'
                     for c in CATEGORY_COLORS if c in used_categories)
    marks = "".join(f'<span><i style="background:{MARK_COLORS[k]}"></i>{escape(MARK_LABELS[k])}</span>'
                    for k in MARK_LABELS if k in used_marks)
    return Markup(svg), Markup(legend), Markup(marks)


def register(templates, ctx) -> None:
    """Expose ``genome_ring(genome)`` to templates: (svg, category legend, mark legend) or None."""
    cache: dict[str, Any] = {}

    def for_template(genome):
        if genome.bin_id in cache:
            return cache[genome.bin_id]
        calls = ctx.cctyper_systems() if hasattr(ctx, "cctyper_systems") else []
        with ctx.lock:
            data = ring_data(ctx.store, ctx.catalog, genome.bin_id, top_categories=ctx.top_categories,
                             cctyper_calls=[c for c in calls if c["bin_id"] == genome.bin_id])
        if data is None:
            out = None
        else:
            sub = f"{len(data['contigs']):,} contigs · {len(data['genes']):,} genes"
            out = render(data, title=genome.label, subtitle=sub)
        if len(cache) > 64:
            cache.clear()
        cache[genome.bin_id] = out
        return out

    templates.env.globals["genome_ring"] = for_template
