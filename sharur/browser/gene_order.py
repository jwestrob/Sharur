"""Gene order between two genomes: a dot plot and the largest collinear blocks.

Genes are matched by family: the protein's best KEGG ortholog, else its exact
resolved Pfam architecture. A family carried once by each genome gives one
unique pair; families carried several times give every pairing, up to
``MAX_PAIRS_PER_FAMILY``. Homology is inferred from these shared annotations,
so unannotated genes take no part.

Collinear blocks are runs of matched genes in the same order on one contig of
each genome (either orientation), allowing up to ``GAP`` unmatched genes
between consecutive matches on either side.
"""

from __future__ import annotations

import math
from collections import defaultdict
from html import escape
from typing import Any
from urllib.parse import quote

from markupsafe import Markup

MAX_PAIRS_PER_FAMILY = 16
GAP = 3
MIN_BLOCK = 3
RIBBONS = 6


def _url(*parts: str) -> str:
    return "/" + "/".join(quote(str(p), safe="") for p in parts)


def load_genome(store, bin_id: str) -> dict[str, Any]:
    """Positioned genes of one genome with their family, and its contigs in plot order."""
    from sharur.architecture import _hits, compact, resolve  # noqa: PLC0415

    rows = store.execute(
        """WITH shared AS (
               SELECT contig_id FROM proteins
               WHERE contig_id IN (SELECT DISTINCT contig_id FROM proteins WHERE bin_id = ?)
               GROUP BY 1 HAVING COUNT(DISTINCT bin_id) > 1)
           SELECT protein_id, contig_id, start, end_coord, strand FROM proteins
           WHERE bin_id = ? AND contig_id <> protein_id AND contig_id NOT IN (SELECT contig_id FROM shared)
           ORDER BY contig_id, start, end_coord""", [bin_id, bin_id])
    ko: dict[str, tuple[str, float]] = {}
    for pid, acc, evalue in store.execute(
            """SELECT a.protein_id, a.accession, a.evalue FROM annotations a JOIN proteins p USING (protein_id)
               WHERE p.bin_id = ? AND LOWER(a.source) IN ('kofam', 'kegg')
                 AND regexp_matches(a.accession, '^K[0-9]{5}$')""", [bin_id]):
        e = evalue if evalue is not None else math.inf
        if pid not in ko or e < ko[pid][1]:
            ko[pid] = (acc, e)
    pfam = _hits(store, ["pfam"], "AND a.protein_id IN (SELECT protein_id FROM proteins WHERE bin_id = ?)", [bin_id])
    lengths = dict(store.execute("SELECT contig_id, length FROM contigs WHERE bin_id = ?", [bin_id]))

    genes: list[dict[str, Any]] = []
    by_contig: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for pid, contig, start, end, strand in rows:
        if pid in ko:
            family = ko[pid][0]
        else:
            domains = resolve(pfam.get(pid, []))
            family = "pfam:" + compact([d.name for d in domains]) if domains else None
        g = {"protein_id": pid, "contig": contig, "start": int(start or 0), "end": int(end or 0),
             "strand": 1 if strand in ("+", "1") else -1, "family": family}
        by_contig[contig].append(g)
        genes.append(g)
    contigs = []
    for contig, members in by_contig.items():
        members.sort(key=lambda g: (g["start"], g["end"]))
        for i, g in enumerate(members):
            g["pos"] = i
        length = max(int(lengths.get(contig) or 0), max(g["end"] for g in members))
        contigs.append((contig, length))
    contigs.sort(key=lambda c: (-c[1], c[0]))
    offset, at = {}, 0
    for contig, length in contigs:
        offset[contig] = at
        at += length
    return {"bin_id": bin_id, "genes": genes, "by_contig": dict(by_contig), "contigs": contigs, "offset": offset,
            "total": at}


def match(a: dict[str, Any], b: dict[str, Any], *, cap: int = MAX_PAIRS_PER_FAMILY) -> dict[str, Any]:
    """Gene pairs sharing a family; ``unique`` marks single-copy families on both sides."""
    fa: dict[str, list] = defaultdict(list)
    fb: dict[str, list] = defaultdict(list)
    for g in a["genes"]:
        if g["family"]:
            fa[g["family"]].append(g)
    for g in b["genes"]:
        if g["family"]:
            fb[g["family"]].append(g)
    shared = sorted(set(fa) & set(fb))
    pairs, skipped = [], 0
    for fam in shared:
        ga, gb = fa[fam], fb[fam]
        if len(ga) * len(gb) > cap:
            skipped += 1
            continue
        unique = len(ga) == 1 and len(gb) == 1
        for x in ga:
            for y in gb:
                pairs.append({"a": x, "b": y, "family": fam, "unique": unique})
    return {"pairs": pairs, "shared": len(shared), "skipped_families": skipped,
            "families_a": len(fa), "families_b": len(fb)}


def blocks(pairs: list[dict[str, Any]], *, gap: int = GAP, min_size: int = MIN_BLOCK) -> list[dict[str, Any]]:
    """Collinear runs: same contig pair, positions advancing together (or B in reverse)."""
    by_contigs: dict[tuple[str, str], list] = defaultdict(list)
    for p in pairs:
        by_contigs[(p["a"]["contig"], p["b"]["contig"])].append(p)
    found = []
    for (ca, cb), group in by_contigs.items():
        group.sort(key=lambda p: (p["a"]["pos"], p["b"]["pos"]))
        chains: list[dict[str, Any]] = []
        used_b: set[int] = set()
        for p in group:
            pa, pb = p["a"]["pos"], p["b"]["pos"]
            if pb in used_b:
                continue
            best, best_cost = None, None
            for ch in chains:
                last = ch["pairs"][-1]
                da = pa - last["a"]["pos"]
                db = pb - last["b"]["pos"]
                if not 0 < da <= gap + 1 or db == 0 or abs(db) > gap + 1:
                    continue
                direction = 1 if db > 0 else -1
                if ch["dir"] not in (0, direction):
                    continue
                cost = da + abs(db)
                if best_cost is None or cost < best_cost:
                    best, best_cost = ch, cost
            if best is None:
                chains.append({"pairs": [p], "dir": 0})
            else:
                last = best["pairs"][-1]
                best["dir"] = 1 if pb > last["b"]["pos"] else -1
                best["pairs"].append(p)
            used_b.add(pb)
        for ch in chains:
            if len(ch["pairs"]) >= min_size:
                ps = ch["pairs"]
                found.append({"contig_a": ca, "contig_b": cb, "pairs": ps, "size": len(ps),
                              "reverse": ch["dir"] < 0,
                              "a_lo": min(p["a"]["pos"] for p in ps), "a_hi": max(p["a"]["pos"] for p in ps),
                              "b_lo": min(p["b"]["pos"] for p in ps), "b_hi": max(p["b"]["pos"] for p in ps)})
    found.sort(key=lambda b: (-b["size"], b["contig_a"], b["a_lo"]))
    return found


def compare_order(a: dict[str, Any], b: dict[str, Any]) -> dict[str, Any]:
    m = match(a, b)
    bl = blocks(m["pairs"])
    in_blocks_a = {p["a"]["protein_id"] for blk in bl for p in blk["pairs"]}
    positioned = len(a["genes"]) or 1
    return {**m, "blocks": bl, "in_blocks": len(in_blocks_a), "share_a": len(in_blocks_a) / positioned,
            "longest": bl[0]["size"] if bl else 0, "unique_pairs": sum(1 for p in m["pairs"] if p["unique"]),
            "genes_a": len(a["genes"]), "genes_b": len(b["genes"])}


# --------------------------------------------------------------------------- #
# Drawing
# --------------------------------------------------------------------------- #

PLOT = 640
PAD = 44


def dotplot(a: dict[str, Any], b: dict[str, Any], result: dict[str, Any], label_a: str, label_b: str) -> Markup:
    size = PLOT
    span = size - PAD - 8
    sx = span / (a["total"] or 1)
    sy = span / (b["total"] or 1)
    parts = [f'<rect x="{PAD}" y="8" width="{span}" height="{span}" class="dp-frame"/>']
    # contig boundaries: lines when few, alternating bands when many
    for genome, horizontal in ((a, False), (b, True)):
        scale = sx if not horizontal else sy
        many = len(genome["contigs"]) > 120
        for i, (contig, length) in enumerate(genome["contigs"]):
            pos = genome["offset"][contig] * scale
            if many:
                w = length * scale
                if i % 2 and w >= 1:            # bands under a pixel wide read as noise
                    parts.append(f'<rect x="{PAD}" y="{8 + pos:.2f}" width="{span}" height="{w:.2f}" class="dp-band"/>'
                                 if horizontal else
                                 f'<rect x="{PAD + pos:.2f}" y="8" width="{w:.2f}" height="{span}" class="dp-band"/>')
            elif i:
                parts.append(f'<line x1="{PAD}" y1="{8 + pos:.2f}" x2="{PAD + span}" y2="{8 + pos:.2f}" class="dp-grid"/>'
                             if horizontal else
                             f'<line x1="{PAD + pos:.2f}" y1="8" x2="{PAD + pos:.2f}" y2="{8 + span}" class="dp-grid"/>')
    for blk in result["blocks"]:
        ps = sorted(blk["pairs"], key=lambda p: p["a"]["pos"])
        pts = " ".join(f'{PAD + (a["offset"][p["a"]["contig"]] + p["a"]["start"]) * sx:.1f},'
                       f'{8 + (b["offset"][p["b"]["contig"]] + p["b"]["start"]) * sy:.1f}' for p in ps)
        parts.append(f'<polyline points="{pts}" class="dp-block"/>')
    for p in result["pairs"]:
        x = PAD + (a["offset"][p["a"]["contig"]] + p["a"]["start"]) * sx
        y = 8 + (b["offset"][p["b"]["contig"]] + p["b"]["start"]) * sy
        same = p["a"]["strand"] == p["b"]["strand"]
        cls = ("dp-dot same" if same else "dp-dot flip") + ("" if p["unique"] else " multi")
        title = f'{p["a"]["protein_id"]} ↔ {p["b"]["protein_id"]} · {p["family"]}'
        parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{2.2 if p["unique"] else 1.5}" class="{cls}" '
                     f'data-a="{escape(p["a"]["protein_id"])}" data-b="{escape(p["b"]["protein_id"])}" '
                     f'data-f="{escape(p["family"])}"><title>{escape(title)}</title></circle>')
    parts.append(f'<text x="{PAD + span / 2}" y="{size - 6}" class="dp-axis" text-anchor="middle">'
                 f'{escape(label_a)} · {len(a["contigs"]):,} contigs, longest first</text>')
    parts.append(f'<text x="14" y="{8 + span / 2}" class="dp-axis" text-anchor="middle" '
                 f'transform="rotate(-90 14 {8 + span / 2})">{escape(label_b)} · {len(b["contigs"]):,} contigs</text>')
    return Markup(f'<svg class="dotplot" viewBox="0 0 {size} {size}" role="img" aria-label="Gene order dot plot">'
                  + "".join(parts) + "</svg>")


def _region(genome: dict[str, Any], contig: str, lo: int, hi: int, flank: int = 1) -> list[dict[str, Any]]:
    genes = genome["by_contig"][contig]
    return genes[max(0, lo - flank): hi + flank + 1]


def ribbon(a: dict[str, Any], b: dict[str, Any], blk: dict[str, Any], width: int = 1000) -> Markup:
    """One block: A's region above, B's below (flipped when the block runs in reverse), matches joined."""
    from sharur.browser.charts import color_for  # noqa: PLC0415

    ra = _region(a, blk["contig_a"], blk["a_lo"], blk["a_hi"])
    rb = _region(b, blk["contig_b"], blk["b_lo"], blk["b_hi"])
    pad = 10
    span = width - 2 * pad

    def placer(region, flip):
        lo = min(g["start"] for g in region)
        hi = max(g["end"] for g in region)
        k = span / max(hi - lo, 1)

        def x(pos):
            v = (pos - lo) * k
            return pad + (span - v if flip else v)
        return x

    xa, xb = placer(ra, False), placer(rb, blk["reverse"])
    ya, yb, h = 14, 92, 16
    parts = []
    matched = {(p["a"]["protein_id"]): p for p in blk["pairs"]}
    for p in blk["pairs"]:
        a0, a1 = sorted((xa(p["a"]["start"]), xa(p["a"]["end"])))
        b0, b1 = sorted((xb(p["b"]["start"]), xb(p["b"]["end"])))
        parts.append(f'<polygon points="{a0:.1f},{ya + h} {a1:.1f},{ya + h} {b1:.1f},{yb} {b0:.1f},{yb}" '
                     f'class="rb-link" style="fill:{color_for(p["family"])}"/>')

    def gene(g, x, y, flip, partner):
        x0, x1 = sorted((x(g["start"]), x(g["end"])))
        forward = (g["strand"] > 0) != flip
        tip = min(6.0, (x1 - x0) / 2)
        pts = (f"{x0:.1f},{y} {x1 - tip:.1f},{y} {x1:.1f},{y + h / 2} {x1 - tip:.1f},{y + h} {x0:.1f},{y + h}"
               if forward else
               f"{x1:.1f},{y} {x0 + tip:.1f},{y} {x0:.1f},{y + h / 2} {x0 + tip:.1f},{y + h} {x1:.1f},{y + h}")
        fill = f' style="fill:{color_for(partner["family"])}"' if partner else ""
        title = f'{g["protein_id"]}' + (f' · {partner["family"]}' if partner else (f' · {g["family"]}' if g["family"] else ""))
        return (f'<a href="{_url("protein", g["protein_id"])}"><polygon points="{pts}" '
                f'class="rb-gene{" on" if partner else ""}"{fill}><title>{escape(title)}</title></polygon></a>')

    by_b = {p["b"]["protein_id"]: p for p in blk["pairs"]}
    for g in ra:
        parts.append(gene(g, xa, ya, False, matched.get(g["protein_id"])))
    for g in rb:
        parts.append(gene(g, xb, yb, blk["reverse"], by_b.get(g["protein_id"])))
    return Markup(f'<svg class="ribbon" viewBox="0 0 {width} 122" role="img" aria-label="Collinear block">'
                  + "".join(parts) + "</svg>")


def region_link(genome: dict[str, Any], contig: str, lo: int, hi: int) -> str:
    genes = genome["by_contig"][contig][lo: hi + 1]
    start, end = min(g["start"] for g in genes), max(g["end"] for g in genes)
    return f'{_url("contig", contig)}?start={max(1, start - 3000)}&span={end - start + 6000}'


def register(templates, ctx) -> None:
    """Expose ``gene_order(a, b)`` to templates for two genome sides; cached per pair."""
    genomes: dict[str, Any] = {}
    results: dict[tuple[str, str], Any] = {}

    def genome(bin_id):
        if bin_id not in genomes:
            if len(genomes) > 24:
                genomes.clear()
            with ctx.lock:
                genomes[bin_id] = load_genome(ctx.store, bin_id)
        return genomes[bin_id]

    def for_template(side_a, side_b):
        key = (side_a.key, side_b.key)
        if key not in results:
            ga, gb = genome(side_a.key), genome(side_b.key)
            if not ga["genes"] or not gb["genes"]:
                results[key] = None
            else:
                r = compare_order(ga, gb)
                r["plot"] = dotplot(ga, gb, r, side_a.label, side_b.label)
                r["ribbons"] = [{"block": blk, "svg": ribbon(ga, gb, blk),
                                 "link_a": region_link(ga, blk["contig_a"], blk["a_lo"], blk["a_hi"]),
                                 "link_b": region_link(gb, blk["contig_b"], blk["b_lo"], blk["b_hi"])}
                                for blk in r["blocks"][:RIBBONS]]
                if len(results) > 32:
                    results.clear()
                results[key] = r
        return results[key]

    templates.env.globals["gene_order"] = for_template
