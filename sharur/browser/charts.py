"""Server-rendered charts (HTML/CSS and inline SVG) for the browser."""

from __future__ import annotations

import html
from typing import Any
from urllib.parse import quote

from markupsafe import Markup

# Categorical palette: distinct in light and dark themes, assigned in fixed order.
# Hues at matched lightness and chroma, so no category shouts over another.
PALETTE = ["#5b8def", "#e0895c", "#4fb38d", "#d9667e", "#9a7fe0", "#d9ae45",
           "#45aec4", "#c97aa4", "#88a85a", "#d1834c", "#6f93c4", "#b48cd6"]
CATEGORY_COLORS = {
    "metabolism": "#4fb38d", "enzyme": "#5b8def", "transport": "#45aec4", "binding": "#d9ae45",
    "info_processing": "#9a7fe0", "regulation": "#b0896a", "envelope": "#e0895c", "stress": "#d9667e",
    "mobile": "#c97aa4", "viral": "#e88fa8", "cazy": "#88a85a", "division": "#6f93c4", "structure": "#9aa3ad",
}
VOG_COLORS = {"Xr": "#9a7fe0", "Xs": "#45aec4", "Xh": "#4fb38d", "Xp": "#c97aa4", "Xu": "#9aa3ad"}
OVERLAY_COLORS = {"prophage": "#c97aa4", "island": "#d9ae45", "defense": "#d9667e", "secretion": "#5b8def",
                  "crispr": "#45aec4"}


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


def contig_strip(contigs: list[tuple[str, int]], url_for, width: int = 1000, height: int = 34) -> Markup:
    """All contigs of a genome, longest first, as one proportional strip; each links to its viewer."""
    total = sum(n for _, n in contigs)
    if not total:
        return Markup("")
    x, parts = 0.0, []
    for i, (contig, length) in enumerate(sorted(contigs, key=lambda c: -c[1])):
        w = length / total * width
        shade = "var(--strip-a)" if i % 2 == 0 else "var(--strip-b)"
        parts.append(f'<a href="{_e(url_for(contig))}"><rect x="{x:.2f}" y="0" width="{max(w, 0.35):.2f}" '
                     f'height="{height}" fill="{shade}"><title>{_e(contig)} · {length:,} bp</title></rect></a>')
        x += w
    return Markup(f'<svg class="strip" viewBox="0 0 {width} {height}" preserveAspectRatio="none" '
                  f'role="img" aria-label="Contigs by length">{"".join(parts)}</svg>')


def contig_overview(length: int, start: int, end: int, genes: list[tuple[int, int]], base_url: str,
                    span: int, width: int = 1000) -> Markup:
    """Whole contig with gene density ticks and the visible window; click a segment to jump there."""
    if not length:
        return Markup("")
    scale = width / length
    ticks = "".join(f'<rect x="{s * scale:.1f}" y="10" width="{max((e - s) * scale, .6):.1f}" height="10" '
                    f'class="ov-gene"/>' for s, e in genes)
    segments = []
    n = 40
    for k in range(n):
        seg_start = int(k * length / n) + 1
        target = max(1, min(seg_start - span // 2 + int(length / n / 2), max(1, length - span + 1)))
        segments.append(f'<a href="{_e(base_url)}?start={target}&amp;span={span}"><rect x="{k * width / n:.1f}" '
                        f'y="0" width="{width / n:.1f}" height="30" class="ov-hit"><title>Jump to '
                        f'{seg_start:,} bp</title></rect></a>')
    window = (f'<rect x="{(start - 1) * scale:.1f}" y="1" width="{max((end - start + 1) * scale, 3):.1f}" '
              f'height="28" rx="4" class="ov-window"/>')
    return Markup(f'<svg class="overview" viewBox="0 0 {width} 30" preserveAspectRatio="none" role="img" '
                  f'aria-label="Contig overview">{ticks}{window}{"".join(segments)}</svg>')


def contig_track(genes: list[dict[str, Any]], overlays: list[dict[str, Any]], start: int, end: int,
                 width: int = 1200, *, contig_start: bool = False, contig_end: bool = False) -> Markup:
    """Genes in the window on two strand lanes, with feature bands above.

    Genes marked ``member`` are outlined and labeled with their ``profile``;
    ``contig_start``/``contig_end`` draw dashed markers where the window meets a contig end.
    """
    span = max(1, end - start + 1)
    scale = (width - 20) / span
    lanes: list[int] = []
    bands = []
    for o in sorted(overlays, key=lambda o: o["start"] or 0):
        if o["end"] is None or o["start"] is None or o["end"] < start or o["start"] > end:
            continue
        x1, x2 = 10 + (max(o["start"], start) - start) * scale, 10 + (min(o["end"], end) - start) * scale
        for lane, free in enumerate(lanes):
            if x1 >= free:
                lanes[lane] = int(x2 + 70)
                break
        else:
            lane = len(lanes)
            lanes.append(int(x2 + 70))
        y = 4 + lane * 18
        color = OVERLAY_COLORS.get(o["label"], OVERLAY_COLORS.get(o["kind"], "#9aa3ad"))
        bands.append(f'<g><title>{_e(o["label"])} {o["start"]:,}–{o["end"]:,}</title><rect x="{x1:.1f}" y="{y}" '
                     f'width="{max(x2 - x1, 3):.1f}" height="14" rx="4" fill="{color}" opacity=".85"/>'
                     f'<text x="{x1 + 4:.1f}" y="{y + 11}" class="band-label">{_e(o["label"])}</text></g>')
    top = 8 + len(lanes) * 18
    axis_y = top + 6
    ticks = []
    step = next(s for s in (500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 250000, 500000, 10**6, 10**7) if span / s <= 12)
    first = ((start + step - 1) // step) * step
    for pos in range(first, end + 1, step):
        x = 10 + (pos - start) * scale
        ticks.append(f'<line x1="{x:.1f}" y1="{axis_y}" x2="{x:.1f}" y2="{axis_y + 112}" class="grid-line"/>'
                     f'<text x="{x + 3:.1f}" y="{axis_y + 10}" class="axis">{pos / 1000:,.0f} kb</text>')
    parts = []
    for g in genes:
        x1, x2 = 10 + (max(g["start"], start) - start) * scale, 10 + (min(g["end"], end) - start) * scale
        forward = g["strand"] != "-"
        y = axis_y + (36 if forward else 76)
        h = 13
        head = min(9.0, (x2 - x1) * 0.45)
        if forward:
            points = f"{x1:.1f},{y - h} {x2 - head:.1f},{y - h} {x2:.1f},{y} {x2 - head:.1f},{y + h} {x1:.1f},{y + h}"
        else:
            points = f"{x2:.1f},{y - h} {x1 + head:.1f},{y - h} {x1:.1f},{y} {x1 + head:.1f},{y + h} {x2:.1f},{y + h}"
        cls = ("gene" if g.get("annotation") or g.get("member") else "gene dark") + (" member" if g.get("member") else "")
        fill = f' style="fill:{g["color"]}"' if g.get("color") else ""
        label = g.get("profile") if g.get("member") and g.get("profile") else (g.get("annotation") or "").split(" (")[0]
        text = ""
        if label and x2 - x1 > 6.2 * len(label) + 12:
            text = f'<text x="{(x1 + x2) / 2:.1f}" y="{y + 4}" class="gene-label">{_e(label)}</text>'
        parts.append(f'<a href="/protein/{quote(g["protein_id"], safe="")}"><g><title>{_e(label or "no annotation")}'
                     f' · {g["length_aa"]} aa · {g["start"]:,}–{g["end"]:,} ({g["strand"]})</title>'
                     f'<polygon class="{cls}" points="{points}"{fill}/>{text}</g></a>')
    if contig_start:
        parts.append(f'<line x1="8" y1="{axis_y + 14}" x2="8" y2="{axis_y + 98}" class="contig-end"/>')
    if contig_end:
        parts.append(f'<line x1="{width - 12}" y1="{axis_y + 14}" x2="{width - 12}" y2="{axis_y + 98}" class="contig-end"/>')
    height = axis_y + 100
    strands = (f'<text x="{width - 6}" y="{axis_y + 26}" class="axis" text-anchor="end">+ strand</text>'
               f'<text x="{width - 6}" y="{axis_y + 98}" class="axis" text-anchor="end">− strand</text>')
    return Markup(f'<svg class="contig-track" viewBox="0 0 {width} {height}" role="img" aria-label="Contig genes">'
                  f'{"".join(ticks)}{"".join(bands)}<line x1="10" y1="{axis_y + 56}" x2="{width - 10}" '
                  f'y2="{axis_y + 56}" class="backbone-line"/>{"".join(parts)}{strands}</svg>')


def histogram(values: list[float], label: str, bins: int = 24, width: int = 420, height: int = 210) -> Markup:
    """Log-free histogram with median marker."""
    if not values:
        return Markup("")
    values = sorted(values)
    hi = values[int(len(values) * 0.99) - 1] if len(values) > 100 else values[-1]
    hi = max(hi, 1)
    counts = [0] * bins
    for v in values:
        counts[min(bins - 1, int(v / hi * bins))] += 1
    top = max(counts) or 1
    bw = (width - 20) / bins
    bars = "".join(f'<rect x="{10 + i * bw:.1f}" y="{height - 24 - c / top * (height - 40):.1f}" '
                   f'width="{bw - 2:.1f}" height="{c / top * (height - 40):.1f}" rx="2" class="hist-bar">'
                   f'<title>{int(i * hi / bins):,}–{int((i + 1) * hi / bins):,}: {c:,}</title></rect>'
                   for i, c in enumerate(counts))
    median = values[len(values) // 2]
    mx = 10 + min(median / hi, 1) * (width - 20)
    return Markup(f'<svg class="hist" viewBox="0 0 {width} {height}" role="img" aria-label="{_e(label)}">{bars}'
                  f'<line x1="{mx:.1f}" y1="8" x2="{mx:.1f}" y2="{height - 22}" class="median"/>'
                  f'<text x="{mx + 4:.1f}" y="16" class="axis">median {median:,.0f}</text>'
                  f'<text x="10" y="{height - 6}" class="axis">0</text>'
                  f'<text x="{width - 10}" y="{height - 6}" class="axis" text-anchor="end">{hi:,.0f}</text>'
                  f'<text x="{width / 2}" y="{height - 6}" class="axis" text-anchor="middle">{_e(label)}</text></svg>')


# --------------------------------------------------------------------------- #
# Protein-scale drawings
# --------------------------------------------------------------------------- #


def _domain_href(d: dict[str, Any]) -> str | None:
    if d.get("href", "") != "":
        return d.get("href")
    return f'/domain/{quote(d["accession"].split(".")[0], safe="")}'


def domain_track(length: int | None, domains: list[dict[str, Any]], width: int = 1000) -> Markup:
    """Domains to scale; labels in lanes below so none overlap.

    A domain's ``label`` replaces its name on the track, and ``href`` its link
    (default: the Pfam domain page; None for no link).
    """
    if not length:
        return Markup("")
    pad, scale = 8, (width - 16) / length
    lanes: list[float] = []
    labels, boxes = [], []
    for d in domains:
        x = pad + (d["start_aa"] - 1) * scale
        w = max(3.0, (d["end_aa"] - d["start_aa"] + 1) * scale)
        color = color_for(d["name"])
        label = d.get("label") or d["name"]
        href = _domain_href(d)
        box = (f'<g class="dom"><title>{_e(d.get("title") or d["name"])} ({_e(d["accession"])}) '
               f'aa {d["start_aa"]}–{d["end_aa"]}</title>'
               f'<rect x="{x:.1f}" y="14" width="{w:.1f}" height="22" rx="6" fill="{color}"/></g>')
        external = ' target="_blank" rel="noopener"' if href and href.startswith("http") else ""
        boxes.append(f'<a href="{_e(href)}"{external}>{box}</a>' if href else box)
        text_w = 6.6 * len(label) + 6
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
                          f'<text x="{x + 1:.1f}" y="{y}">{_e(label)}</text>')
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


# --------------------------------------------------------------------------- #
# Small drawings for lists (Discover)
# --------------------------------------------------------------------------- #


def mini_architecture(length: int, domains: list[dict[str, Any]], scale: int, *, width: int = 600,
                      height: int = 16, highlight: tuple[int, int] | None = None, partial: bool = False) -> Markup:
    """A protein to a shared scale: backbone, Pfam domains, an optional highlighted span.

    ``scale`` is the length that fills ``width``; ``partial`` adds a dashed tail
    for a gene that runs off its contig.
    """
    if not length or not scale:
        return Markup("")
    usable = width - 10
    w = max(6.0, usable * min(length, scale) / scale)
    k = w / length
    mid = height / 2
    parts = [f'<rect x="1" y="{mid - 1.5:.1f}" width="{w:.1f}" height="3" rx="1.5" class="mini-backbone"/>']
    if highlight:
        lo, hi = highlight
        parts.append(f'<rect x="{1 + (lo - 1) * k:.1f}" y="0.5" width="{max(2.0, (hi - lo + 1) * k):.1f}" '
                     f'height="{height - 1}" rx="3" class="mini-highlight"/>')
    for d in domains:
        if d.get("start_aa") is None:
            continue
        x = 1 + (d["start_aa"] - 1) * k
        dw = max(1.5, (d["end_aa"] - d["start_aa"] + 1) * k)
        parts.append(f'<rect x="{x:.1f}" y="3" width="{dw:.1f}" height="{height - 6}" rx="2" '
                     f'fill="{color_for(d["name"])}"><title>{_e(d["name"])} aa {d["start_aa"]}–{d["end_aa"]}</title></rect>')
    if partial:
        parts.append(f'<line x1="{1 + w:.1f}" y1="{mid:.1f}" x2="{min(width - 1, w + 9):.1f}" y2="{mid:.1f}" '
                     f'class="mini-partial"/>')
    # stretches across its cell; only horizontal extents carry meaning
    return Markup(f'<svg class="mini mini-fill" viewBox="0 0 {width} {height}" width="100%" height="{height}" '
                  f'preserveAspectRatio="none" role="img" aria-label="{_e(length)} aa">{"".join(parts)}</svg>')


def gene_arrows(strands: str, *, gene: int = 9, gap: int = 2, height: int = 14, max_genes: int = 40) -> Markup:
    """Consecutive genes as small strand arrows (``strands``: one '+' or '-' per gene)."""
    shown = strands[:max_genes]
    parts = []
    for i, s in enumerate(shown):
        x = i * (gene + gap)
        tip = 3
        if s == "+":
            pts = f"{x},2 {x + gene - tip},2 {x + gene},{height / 2} {x + gene - tip},{height - 2} {x},{height - 2}"
        else:
            pts = f"{x + gene},2 {x + tip},2 {x},{height / 2} {x + tip},{height - 2} {x + gene},{height - 2}"
        parts.append(f'<polygon points="{pts}" class="mini-gene"/>')
    width = len(shown) * (gene + gap) + (14 if len(strands) > max_genes else 0)
    if len(strands) > max_genes:
        parts.append(f'<text x="{width - 12}" y="{height - 3}" class="mini-more">…</text>')
    return Markup(f'<svg class="mini" viewBox="0 0 {width} {height}" width="{width}" height="{height}" '
                  f'role="img" aria-label="{len(strands)} genes">{"".join(parts)}</svg>')
