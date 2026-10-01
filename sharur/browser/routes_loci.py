"""Stacked locus views: every call of a system type, or every neighborhood of a family.

Rows are aligned on an anchor gene (a system's core component, or the carrier
of a domain, VOG or function) and flipped so the anchor always points right.
Genes take a color from their family, consistent across rows, so conserved
architecture reads as vertical bands.
"""

from __future__ import annotations

import hashlib
from collections import Counter, defaultdict
from typing import Any

from fastapi import HTTPException, Query, Request
from fastapi.responses import HTMLResponse
from markupsafe import Markup

from sharur.browser import charts
from sharur.browser.catalog import RANKS

WIDTH = 1040


def _family_color(key: str) -> str:
    digest = int(hashlib.md5(key.encode()).hexdigest()[:6], 16)
    return charts.PALETTE[digest % len(charts.PALETTE)]


def _family_key(label: str | None) -> str | None:
    """Collapse a hit label to its family: 'GH161 (GH161)' -> 'GH161'."""
    if not label:
        return None
    return label.split(" (")[0].strip() or None


def build_rows(ctx, anchors: list[dict[str, Any]], flank: int) -> list[dict[str, Any]]:
    """Gene windows around anchors.

    ``anchors``: dicts with ``anchor`` (protein id), ``members`` (ids to highlight),
    ``profiles`` (id -> label for highlighted genes), ``bin_id`` and ``title``.
    """
    store, lock = ctx.store, ctx.lock
    ids = list({a["anchor"] for a in anchors} | {m for a in anchors for m in a["members"]})
    if not ids:
        return []
    with lock:
        located = dict((pid, contig) for pid, contig in store.execute(
            "SELECT protein_id, contig_id FROM proteins WHERE protein_id IN (SELECT UNNEST(?::VARCHAR[]))", [ids]))
        contigs = sorted({c for c in located.values() if c})
        genes_by_contig: dict[str, list[tuple]] = defaultdict(list)
        for row in store.execute(
                """SELECT contig_id, protein_id, start, end_coord, strand, sequence_length FROM proteins
                   WHERE contig_id IN (SELECT UNNEST(?::VARCHAR[])) ORDER BY contig_id, start, end_coord""",
                [contigs]):
            genes_by_contig[row[0]].append(row[1:])
    windows = []
    needed: set[str] = set()
    for a in anchors:
        contig = located.get(a["anchor"])
        genes = genes_by_contig.get(contig or "", [])
        index = {g[0]: i for i, g in enumerate(genes)}
        if a["anchor"] not in index:
            continue
        member_idx = [index[m] for m in a["members"] if m in index] or [index[a["anchor"]]]
        lo, hi = max(0, min(member_idx) - flank), min(len(genes) - 1, max(member_idx) + flank)
        window = genes[lo:hi + 1]
        needed.update(g[0] for g in window)
        windows.append((a, contig, window, lo == 0, hi == len(genes) - 1))
    with lock:
        best = dict(store.execute(
            """SELECT protein_id, ARG_MIN(COALESCE(NULLIF(name, ''), accession), COALESCE(evalue, 1))
               FROM annotations WHERE protein_id IN (SELECT UNNEST(?::VARCHAR[]))
               AND LOWER(source) NOT IN ('defensefinder_system', 'txsscan_system', 'hyddb_subgroup')
               GROUP BY 1""", [list(needed)])) if needed else {}
    rows = []
    for a, contig, window, at_start, at_end in windows:
        anchor = next(g for g in window if g[0] == a["anchor"])
        flip = anchor[3] == "-"
        origin = anchor[1] if not flip else anchor[2]
        genes = []
        for pid, start, end, strand, length in window:
            if flip:
                x1, x2, s = origin - end, origin - start, ("+" if strand == "-" else "-")
            else:
                x1, x2, s = start - origin, end - origin, strand
            label = ctx.describe_hit(best.get(pid))
            member = pid in a["members"]
            key = a["profiles"].get(pid) if member and a["profiles"].get(pid) else _family_key(label)
            genes.append({"protein_id": pid, "x1": x1, "x2": x2, "strand": s, "length_aa": length,
                          "label": label, "key": key, "member": member, "anchor": pid == a["anchor"]})
        left_edge, right_edge = (at_end, at_start) if flip else (at_start, at_end)
        rows.append({**a, "contig": contig, "genes": genes, "flipped": flip,
                     "left_edge": left_edge, "right_edge": right_edge})
    return rows


def render_stack(rows: list[dict[str, Any]], *, member_label: str = "system gene") -> tuple[Markup, list]:
    """One SVG per row on a shared bp scale; returns markup and the family legend."""
    if not rows:
        return Markup(""), []
    lo = min(g["x1"] for r in rows for g in r["genes"])
    hi = max(g["x2"] for r in rows for g in r["genes"])
    span = max(1, hi - lo)
    pad = 14
    scale = (WIDTH - 2 * pad) / span
    counts = Counter(g["key"] for r in rows for g in r["genes"] if g["key"])
    legend = [(k, _family_color(k), n) for k, n in counts.most_common(14)]
    out = []
    for r in rows:
        parts = [f'<line x1="0" y1="22" x2="{WIDTH}" y2="22" class="backbone-line"/>']
        if r["left_edge"]:
            x = pad + (min(g["x1"] for g in r["genes"]) - lo) * scale - 6
            parts.append(f'<line x1="{x:.1f}" y1="6" x2="{x:.1f}" y2="38" class="contig-end"/>')
        if r["right_edge"]:
            x = pad + (max(g["x2"] for g in r["genes"]) - lo) * scale + 6
            parts.append(f'<line x1="{x:.1f}" y1="6" x2="{x:.1f}" y2="38" class="contig-end"/>')
        for g in r["genes"]:
            x1, x2 = pad + (g["x1"] - lo) * scale, pad + (g["x2"] - lo) * scale
            head = min(8.0, (x2 - x1) * 0.45)
            y, h = 22, 11 if (g["member"] or g["anchor"]) else 8
            if g["strand"] == "-":
                pts = f"{x2:.1f},{y - h} {x1 + head:.1f},{y - h} {x1:.1f},{y} {x1 + head:.1f},{y + h} {x2:.1f},{y + h}"
            else:
                pts = f"{x1:.1f},{y - h} {x2 - head:.1f},{y - h} {x2:.1f},{y} {x2 - head:.1f},{y + h} {x1:.1f},{y + h}"
            color = _family_color(g["key"]) if g["key"] else None
            cls = "gene" + (" member" if g["member"] else "") + (" anchor-ring" if g["anchor"] else "") + (
                "" if color else " dark")
            style = f' style="fill:{color}"' if color else ""
            title = charts._e(f'{g["key"] or "no annotation"}{" · " + member_label if g["member"] else ""}'
                              f' · {g["length_aa"]} aa')
            text = ""
            if g["key"] and x2 - x1 > 6.0 * len(g["key"]) + 10:
                text = f'<text x="{(x1 + x2) / 2:.1f}" y="{y + 4}" class="gene-label">{charts._e(g["key"])}</text>'
            parts.append(f'<a href="/protein/{charts.quote(g["protein_id"], safe="")}"><g><title>{title}</title>'
                         f'<polygon class="{cls}" points="{pts}"{style}/>{text}</g></a>')
        out.append(Markup(f'<svg class="stack-row" viewBox="0 0 {WIDTH} 44" role="img" '
                          f'aria-label="Locus">{"".join(parts)}</svg>'))
    return out, legend


def register(app, ctx) -> None:
    catalog = ctx.catalog

    def clade_filter(genomes_rank: str | None, name: str | None):
        if not genomes_rank or genomes_rank not in RANKS or not name:
            return None
        return {g.bin_id for g in catalog.clade(genomes_rank, name)}

    @app.get("/system/{kind}/{system_type:path}/loci", response_class=HTMLResponse)
    def system_loci(request: Request, kind: str, system_type: str, subtype: str = Query(""),
                    rank: str = Query(""), clade: str = Query(""), page: int = Query(1, ge=1),
                    flank: int = Query(3, ge=0, le=12)):
        calls = [s for s in catalog.systems if s["kind"] == kind and s["type"] == system_type]
        if not calls:
            raise HTTPException(404, "System type not found")
        keep = clade_filter(rank, clade)
        selected = [s for s in calls if (not subtype or (s["subtype"] or "") == subtype)
                    and (keep is None or s["bin_id"] in keep)]
        per_page = 25
        pages = max(1, (len(selected) + per_page - 1) // per_page)
        page = min(page, pages)
        shown = selected[(page - 1) * per_page: page * per_page]
        with ctx.lock:
            profiles = defaultdict(dict)
            positions: dict[str, list[int]] = defaultdict(list)
            for system_id, protein_id, profile, position in ctx.store.execute(
                    "SELECT system_id, protein_id, profile_name, position FROM system_proteins "
                    "WHERE system_id IN (SELECT UNNEST(?::VARCHAR[]))", [[s["system_id"] for s in shown]]):
                name = (profile or "").split("__")[-1] or None
                profiles[system_id][protein_id] = name
                if name and position is not None:
                    positions[name].append(position)
        # the core component: present in the most calls on the page, then earliest in the system
        core_counts = Counter(p for s in shown for p in set(profiles[s["system_id"]].values()) if p)
        core = min(core_counts, key=lambda p: (-core_counts[p], sum(positions[p]) / max(1, len(positions[p])), p)) \
            if core_counts else None
        anchors = []
        for s in shown:
            prof = profiles[s["system_id"]]
            members = [p for p in s["proteins"]] or list(prof)
            anchor = next((p for p, name in prof.items() if name == core), members[0] if members else None)
            if anchor:
                g = catalog.by_bin.get(s["bin_id"])
                anchors.append({"anchor": anchor, "members": set(members), "profiles": prof, "bin_id": s["bin_id"],
                                "title": s["subtype"] or s["type"], "system_id": s["system_id"],
                                "genome_label": g.label if g else ""})
        rows = build_rows(ctx, anchors, flank)
        svgs, legend = render_stack(rows)
        subtypes = Counter(s["subtype"] or "–" for s in calls).most_common()
        return ctx.render(request, "loci.html", "systems", mode="system", kind=kind, title=system_type,
                          rows=list(zip(rows, svgs)), legend=legend, total=len(selected), page=page, pages=pages,
                          subtype=subtype, subtypes=subtypes, rank=rank, clade=clade, flank=flank, core=core,
                          base=ctx.url("system", kind, system_type) + "/loci", back=ctx.url("system", kind, system_type))

    @app.get("/stack/{kind}/{ident:path}", response_class=HTMLResponse)
    def family_stack(request: Request, kind: str, ident: str, rank: str = Query(""), clade: str = Query(""),
                     limit: int = Query(30, ge=5, le=100), flank: int = Query(6, ge=1, le=15)):
        with ctx.lock:
            if kind == "domain":
                rows = ctx.store.execute(
                    """SELECT DISTINCT a.protein_id, p.bin_id FROM annotations a JOIN proteins p USING (protein_id)
                       WHERE LOWER(a.source) = 'pfam' AND split_part(a.accession, '.', 1) = ?""", [ident])
                name = catalog.domains.get(ident, {}).get("name", ident)
                back = ctx.url("domain", ident)
            elif kind == "vog":
                rows = ctx.store.execute(
                    """SELECT DISTINCT a.protein_id, p.bin_id FROM annotations a JOIN proteins p USING (protein_id)
                       WHERE LOWER(a.source) IN ('vogdb', 'vog') AND a.accession = ?""", [ident])
                name = catalog.vogs.get(ident, {}).get("name", ident)
                back = ctx.url("vog", ident)
            elif kind == "function":
                rows = ctx.store.execute(
                    """SELECT pp.protein_id, p.bin_id FROM protein_predicates pp JOIN proteins p USING (protein_id)
                       WHERE list_contains(pp.predicates, ?)""", [ident])
                definition = ctx.predicates.get(ident)
                name = definition.name if definition else ident
                back = ctx.url("function", ident)
            else:
                raise HTTPException(404, "Unknown stack kind")
        if not rows:
            raise HTTPException(404, "No carriers")
        keep = clade_filter(rank, clade)
        rows = [r for r in rows if keep is None or r[1] in keep]
        # spread the sample across genomes: one carrier per genome first, in a stable order
        by_genome: dict[str, list[str]] = defaultdict(list)
        for pid, bin_id in sorted(rows, key=lambda r: hashlib.md5(r[0].encode()).hexdigest()):
            by_genome[bin_id].append(pid)
        sample: list[tuple[str, str]] = []
        depth = 0
        while len(sample) < limit and any(len(v) > depth for v in by_genome.values()):
            for bin_id, pids in sorted(by_genome.items(), key=lambda kv: hashlib.md5(kv[0].encode()).hexdigest()):
                if len(pids) > depth and len(sample) < limit:
                    sample.append((pids[depth], bin_id))
            depth += 1
        anchors = []
        for pid, bin_id in sample:
            g = catalog.by_bin.get(bin_id)
            anchors.append({"anchor": pid, "members": {pid}, "profiles": {pid: name}, "bin_id": bin_id,
                            "title": "", "genome_label": g.label if g else ""})
        built = build_rows(ctx, anchors, flank)
        svgs, legend = render_stack(built, member_label="carrier")
        return ctx.render(request, "loci.html", "domains" if kind in ("domain", "vog") else "functions",
                          mode="family", kind=kind, ident=ident, title=name, rows=list(zip(built, svgs)),
                          legend=legend, total=len(rows), shown=len(built), genomes=len(by_genome), rank=rank,
                          clade=clade, flank=flank, limit=limit, base=ctx.url("stack", kind, ident), back=back,
                          page=1, pages=1)
