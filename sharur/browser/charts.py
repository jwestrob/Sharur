"""Server-rendered charts (HTML/CSS and inline SVG) for the browser."""

from __future__ import annotations

import html
from typing import Any
from urllib.parse import quote

from markupsafe import Markup

# Categorical palette: distinct in light and dark themes, assigned in fixed order.
PALETTE = ["#4e79a7", "#f28e2b", "#59a14f", "#e15759", "#76b7b2", "#edc948",
           "#b07aa1", "#ff9da7", "#9c755f", "#86bcb6", "#d37295", "#a0cbe8"]
CATEGORY_COLORS = {
    "metabolism": "#59a14f", "enzyme": "#4e79a7", "transport": "#76b7b2", "binding": "#edc948",
    "info_processing": "#b07aa1", "regulation": "#9c755f", "envelope": "#f28e2b", "stress": "#e15759",
    "mobile": "#d37295", "viral": "#ff9da7", "cazy": "#86bcb6", "division": "#a0cbe8", "structure": "#bab0ac",
}


def _e(value: Any) -> str:
    return html.escape("" if value is None else str(value))


def color_for(name: str) -> str:
    return PALETTE[sum(map(ord, name)) % len(PALETTE)]


# --------------------------------------------------------------------------- #
# Treemap (squarified), rendered as positioned divs so it scales with the page
# --------------------------------------------------------------------------- #


def _squarify(values: list[float], x: float, y: float, w: float, h: float) -> list[tuple[float, float, float, float]]:
    rects: list[tuple[float, float, float, float]] = []
    values = list(values)
    while values:
        short = min(w, h)
        total = sum(values)
        row: list[float] = []
        best = float("inf")
        for v in values:
            candidate = row + [v]
            s = sum(candidate) * w * h / total
            side = s / short
            worst = max(max((c * w * h / total) / side ** 2, side ** 2 / (c * w * h / total)) for c in candidate)
            if worst > best:
                break
            row, best = candidate, worst
        s = sum(row) * w * h / total
        side = s / short
        offset = 0.0
        for v in row:
            length = (v * w * h / total) / side
            if w >= h:
                rects.append((x, y + offset, side, length))
            else:
                rects.append((x + offset, y, length, side))
            offset += length
        if w >= h:
            x, w = x + side, w - side
        else:
            y, h = y + side, h - side
        values = values[len(row):]
    return rects


def treemap(items: list[dict[str, Any]], height: int = 320) -> Markup:
    """``items``: name, value, url, sub (optional caption)."""
    items = [i for i in items if i["value"] > 0]
    if not items:
        return Markup("")
    items.sort(key=lambda i: -i["value"])
    rects = _squarify([i["value"] for i in items], 0, 0, 100.0, 100.0)
    cells = []
    for item, (x, y, w, h) in zip(items, rects, strict=True):
        label = (f'<span class="tm-name">{_e(item["name"])}</span>'
                 f'<span class="tm-sub">{_e(item.get("sub", ""))}</span>') if w * h > 18 else ""
        cells.append(
            f'<a class="tm-cell" href="{_e(item["url"])}" title="{_e(item["name"])} · {_e(item.get("sub", ""))}" '
            f'style="left:{x:.3f}%;top:{y:.3f}%;width:{w:.3f}%;height:{h:.3f}%;'
            f'--tm:{item.get("color") or color_for(item["name"])}">{label}</a>')
    return Markup(f'<div class="treemap" style="height:{height}px">{"".join(cells)}</div>')


# --------------------------------------------------------------------------- #
# Genome-scale drawings
# --------------------------------------------------------------------------- #


def contig_strip(lengths: list[int], width: int = 1000, height: int = 28) -> Markup:
    """All contigs of a genome, longest first, as one proportional strip."""
    total = sum(lengths)
    if not total:
        return Markup("")
    x, parts = 0.0, []
    for i, length in enumerate(sorted(lengths, reverse=True)):
        w = length / total * width
        shade = "var(--strip-a)" if i % 2 == 0 else "var(--strip-b)"
        parts.append(f'<rect x="{x:.2f}" y="0" width="{max(w, 0.3):.2f}" height="{height}" fill="{shade}">'
                     f'<title>{length:,} bp</title></rect>')
        x += w
    return Markup(f'<svg class="strip" viewBox="0 0 {width} {height}" preserveAspectRatio="none" '
                  f'role="img" aria-label="Contigs by length">{"".join(parts)}</svg>')


# --------------------------------------------------------------------------- #
# Protein-scale drawings
# --------------------------------------------------------------------------- #


def domain_track(length: int | None, domains: list[dict[str, Any]], width: int = 1000) -> Markup:
    """Domains to scale; labels in lanes below so none overlap."""
    if not length:
        return Markup("")
    pad, scale = 8, (width - 16) / length
    lanes: list[float] = []
    labels, boxes = [], []
    for d in domains:
        x = pad + (d["start_aa"] - 1) * scale
        w = max(3.0, (d["end_aa"] - d["start_aa"] + 1) * scale)
        color = color_for(d["name"])
        boxes.append(f'<g><title>{_e(d["name"])} ({_e(d["accession"])}) aa {d["start_aa"]}–{d["end_aa"]}</title>'
                     f'<rect x="{x:.1f}" y="14" width="{w:.1f}" height="22" rx="5" fill="{color}"/></g>')
        text_w = 6.6 * len(d["name"]) + 6
        for lane, free in enumerate(lanes):
            if x >= free:
                lanes[lane] = x + text_w
                break
        else:
            lane = len(lanes)
            lanes.append(x + text_w)
        if lane < 4:
            y = 52 + lane * 15
            labels.append(f'<line x1="{x + 1:.1f}" y1="36" x2="{x + 1:.1f}" y2="{y - 10}" class="tick"/>'
                          f'<text x="{x + 1:.1f}" y="{y}">{_e(d["name"])}</text>')
    height = 52 + min(len(lanes), 4) * 15
    ticks = []
    step = next(s for s in (50, 100, 250, 500, 1000, 2500, 5000, 10000, 25000) if length / s <= 10)
    for aa in range(0, length + 1, step):
        x = pad + aa * scale
        ticks.append(f'<text x="{x:.1f}" y="9" class="axis">{aa:,}</text>')
    return Markup(f'<svg class="track" viewBox="0 0 {width} {height}" role="img" aria-label="Domain architecture">'
                  f'{"".join(ticks)}<rect x="{pad}" y="23" width="{length * scale:.1f}" height="4" rx="2" '
                  f'class="backbone"/>{"".join(boxes)}{"".join(labels)}</svg>')


def neighborhood(genes: list[dict[str, Any]], *, start_edge: bool, end_edge: bool, width: int = 1000) -> Markup:
    """Genes to scale with strand arrows, numbered to match the legend table."""
    if not genes:
        return Markup("")
    lo, hi = min(g["start"] for g in genes), max(g["end"] for g in genes)
    pad = 24
    scale = (width - 2 * pad) / max(1, hi - lo)
    parts = ['<line x1="0" y1="40" x2="{0}" y2="40" class="backbone-line"/>'.format(width)]
    if start_edge:
        parts.append(f'<line x1="{pad - 10}" y1="14" x2="{pad - 10}" y2="66" class="contig-end"/>'
                     f'<text x="{pad - 6}" y="12" class="axis">contig start</text>')
    if end_edge:
        parts.append(f'<line x1="{width - pad + 10}" y1="14" x2="{width - pad + 10}" y2="66" class="contig-end"/>'
                     f'<text x="{width - pad - 62}" y="12" class="axis">contig end</text>')
    for g in genes:
        x1, x2 = pad + (g["start"] - lo) * scale, pad + (g["end"] - lo) * scale
        head = min(10.0, (x2 - x1) * 0.45)
        y, h = 40, 12
        if g["strand"] == "-":
            points = f"{x2:.1f},{y - h} {x1 + head:.1f},{y - h} {x1:.1f},{y} {x1 + head:.1f},{y + h} {x2:.1f},{y + h}"
        else:
            points = f"{x1:.1f},{y - h} {x2 - head:.1f},{y - h} {x2:.1f},{y} {x2 - head:.1f},{y + h} {x1:.1f},{y + h}"
        cls = "gene anchor" if g.get("is_anchor") else ("gene truncated" if g.get("edge_status") == "truncated"
                                                       else ("gene" if g.get("annotation") else "gene dark"))
        fill = "" if g.get("is_anchor") else (f' style="fill:{g["color"]}"' if g.get("color") else "")
        parts.append(
            f'<a href="/protein/{quote(g["protein_id"], safe="")}"><g><title>{_e(g["number"])}. '
            f'{_e(g.get("annotation") or "no annotation")}</title><polygon class="{cls}" points="{points}"{fill}/>'
            f'<text x="{(x1 + x2) / 2:.1f}" y="{y + 4}" class="gene-num">{g["number"]}</text></g></a>')
    return Markup(f'<svg class="hood" viewBox="0 0 {width} 72" role="img" aria-label="Gene neighborhood">'
                  f'{"".join(parts)}</svg>')
