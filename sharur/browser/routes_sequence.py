"""Protein-page sequence properties, loaded after the page: hydropathy, segments, composition, self-similarity.

Registered by :func:`sharur.browser.app.create_app` before the catch-all protein
route. The computations live in :mod:`sharur.browser.sequence_properties`;
this module draws them on the domain-architecture scale (the same 1000-unit
frame and padding as :func:`sharur.browser.charts.domain_track`).
"""

from __future__ import annotations

from collections import OrderedDict
from html import escape
from threading import Lock
from types import SimpleNamespace
from typing import Any

import numpy as np
from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import HTMLResponse
from markupsafe import Markup

from sharur.browser import charts
from sharur.browser import sequence_properties as sp

WIDTH, PAD = 1000, 8            # matches charts.domain_track
KD_LO, KD_HI = -4.5, 4.5
CACHE = 24


def _x(pos: float, length: int) -> float:
    return PAD + (pos - 1) * (WIDTH - 2 * PAD) / length


def profile_svg(a: dict[str, Any]) -> Markup:
    """Hydropathy (band = window range per pixel, line = mean), threshold, and segment rows."""
    n, profile = a["length"], a["profile"]
    top, base = 4, 60
    span = base - top

    def y(v: float) -> float:
        return top + span * (KD_HI - min(max(v, KD_LO), KD_HI)) / (KD_HI - KD_LO)

    parts = [f'<line x1="{PAD}" x2="{WIDTH - PAD}" y1="{y(0):.1f}" y2="{y(0):.1f}" class="seqp-zero"/>',
             f'<line x1="{PAD}" x2="{WIDTH - PAD}" y1="{y(sp.TM_THRESHOLD):.1f}" y2="{y(sp.TM_THRESHOLD):.1f}" '
             f'class="seqp-threshold"><title>Kyte-Doolittle 1.6</title></line>']
    if len(profile):
        cols = min(WIDTH - 2 * PAD, len(profile))
        edges = np.linspace(0, len(profile), cols + 1).astype(int)
        lo = np.minimum.reduceat(profile, edges[:-1])
        hi = np.maximum.reduceat(profile, edges[:-1])
        mean = np.add.reduceat(profile, edges[:-1]) / np.maximum(np.diff(edges), 1)
        centre = (edges[:-1] + np.diff(edges) / 2) + 1 + sp.TM_WINDOW // 2      # window centre, 1-based residue
        xs = [_x(c, n) for c in centre]
        band = " ".join(f"{x:.1f},{y(v):.1f}" for x, v in zip(xs, hi)) + " " + \
            " ".join(f"{x:.1f},{y(v):.1f}" for x, v in zip(reversed(xs), reversed(lo)))
        parts.append(f'<polygon points="{band}" class="seqp-band"/>')
        parts.append('<polyline points="' + " ".join(f"{x:.1f},{y(v):.1f}" for x, v in zip(xs, mean))
                     + '" class="seqp-line"/>')
    for row, (items, cls, label) in enumerate(((a["segments"], "seqp-tm", "hydrophobic segment"),
                                               (a["low_complexity"], "seqp-lc", "low complexity"))):
        y0 = 66 + row * 12
        parts.append(f'<rect x="{PAD}" y="{y0 + 4}" width="{WIDTH - 2 * PAD}" height="1" class="seqp-rail"/>')
        for it in items:
            x0, x1 = _x(it["start"], n), _x(it["end"] + 1, n)
            extra = f' · {it["top"]} {it["top_share"]:.0%}' if "top" in it else f' · peak {it["peak"]}'
            parts.append(f'<rect x="{x0:.1f}" y="{y0}" width="{max(1.5, x1 - x0):.1f}" height="9" rx="2" class="{cls}">'
                         f'<title>{label}: aa {it["start"]:,}–{it["end"]:,}{extra}</title></rect>')
    nt = a["n_terminal"]
    if nt:
        x0, x1 = _x(nt["start"], n), _x(nt["end"] + 1, n)
        parts.append(f'<rect x="{x0:.1f}" y="{top}" width="{max(3.0, x1 - x0):.1f}" height="{span}" class="seqp-nterm">'
                     f'<title>N-terminal hydrophobic stretch aa {nt["start"]}–{nt["end"]} (mean {nt["mean"]})</title></rect>')
    return Markup(f'<svg class="track seqp-fig" viewBox="0 0 {WIDTH} 90" role="img" '
                  f'aria-label="Hydropathy and sequence segments">{"".join(parts)}</svg>')


def dotplot_svg(a: dict[str, Any], domains: list[dict[str, Any]]) -> Markup | None:
    s = a["self"]
    if not s:
        return None
    n, bins = a["length"], s["bins"]
    off, size = 46, 560
    k = size / n
    png = sp.dotplot_png(s["cells"])
    parts = [f'<rect x="{off}" y="{off}" width="{size}" height="{size}" class="seqp-dot-bg"/>',
             f'<image x="{off}" y="{off}" width="{size}" height="{size}" preserveAspectRatio="none" '
             f'style="image-rendering:pixelated" href="data:image/png;base64,{png}"/>',
             f'<line x1="{off}" y1="{off}" x2="{off + size}" y2="{off + size}" class="seqp-diag"/>']
    for d in domains:
        if d.get("start_aa") is None:
            continue
        a0, a1 = (d["start_aa"] - 1) * k, d["end_aa"] * k
        color = charts.color_for(d["name"])
        title = f'<title>{escape(d["name"])} aa {d["start_aa"]}–{d["end_aa"]}</title>'
        parts.append(f'<rect x="{off + a0:.1f}" y="{off - 14}" width="{max(1.5, a1 - a0):.1f}" height="9" rx="2" '
                     f'fill="{color}">{title}</rect>')
        parts.append(f'<rect x="{off - 14}" y="{off + a0:.1f}" width="9" height="{max(1.5, a1 - a0):.1f}" rx="2" '
                     f'fill="{color}">{title}</rect>')
    p = s["period"]
    if p:
        a0, a1 = (p["start"] - 1) * k, p["end"] * k
        parts.append(f'<line x1="{off + a0:.1f}" x2="{off + a1:.1f}" y1="{off + size + 8}" y2="{off + size + 8}" '
                     f'class="seqp-repeat"><title>Repeat stretch aa {p["start"]:,}–{p["end"]:,}</title></line>')
    for frac in (0, 0.5, 1):
        pos = max(1, round(frac * n))
        parts.append(f'<text x="{off + frac * size:.1f}" y="{off + size + 24}" class="seqp-axis" '
                     f'text-anchor="{"start" if frac == 0 else "end" if frac == 1 else "middle"}">{pos:,}</text>')
    return Markup(f'<svg class="seqp-fig seqp-dotplot" viewBox="0 0 {off + size + 8} {off + size + 30}" role="img" '
                  f'aria-label="Self-similarity dot plot">{"".join(parts)}</svg>')


def _covered(regions: list[dict[str, Any]], n: int) -> float:
    mask = np.zeros(n + 1, dtype=bool)
    for r in regions:
        mask[r["start"]:r["end"] + 1] = True
    return float(mask.sum()) / n if n else 0.0


def register(app: FastAPI, ctx: SimpleNamespace) -> None:
    """Add ``/protein/{id}/properties`` (an HTML fragment). ``ctx``: store, lock, render."""
    store = ctx.store
    cache: OrderedDict[str, Any] = OrderedDict()
    cache_lock = Lock()

    @app.get("/protein/{protein_id:path}/properties", response_class=HTMLResponse)
    def protein_properties(request: Request, protein_id: str):
        from sharur.architecture import architecture  # noqa: PLC0415

        with cache_lock:
            hit = cache.get(protein_id)
        if hit is None:
            with ctx.lock:
                rows = store.execute("SELECT sequence FROM proteins WHERE protein_id = ?", [protein_id])
                if not rows:
                    raise HTTPException(404, "Protein not found")
                domains = [d.to_dict() for d in architecture(store, protein_id)]
            analysis = sp.analyse(rows[0][0] or "")
            hit = (analysis, domains)
            with cache_lock:
                cache[protein_id] = hit
                while len(cache) > CACHE:
                    cache.popitem(last=False)
        analysis, domains = hit
        if analysis is None:
            return HTMLResponse('<div class="panel"><p class="empty">No sequence stored for this protein.</p></div>')
        n = analysis["length"]
        return ctx.render(request, "_sequence_properties.html", "taxa", a=analysis, comp=analysis["composition"],
                          profile=profile_svg(analysis), dot=dotplot_svg(analysis, domains),
                          lc_share=_covered(analysis["low_complexity"], n),
                          tm_share=_covered(analysis["segments"], n))
