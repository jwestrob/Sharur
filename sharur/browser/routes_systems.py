"""One page per system call: its genes in genomic context.

The call's member genes are outlined and labeled with their profile; the rest
of the window takes function-category colors. Other systems, prophages,
islands and CRISPR arrays overlapping the window appear as bands above.
"""

from __future__ import annotations

from typing import Any

from fastapi import HTTPException, Query, Request
from fastapi.responses import HTMLResponse

from sharur.browser import charts

FLANKS = (5000, 10000, 20000, 40000)


def register(app, ctx) -> None:
    catalog = ctx.catalog

    @app.get("/call/{system_id:path}", response_class=HTMLResponse)
    def system_call(request: Request, system_id: str, flank: int = Query(10000, ge=1000, le=100000)):
        call = next((s for s in catalog.systems if s["system_id"] == system_id), None)
        if call is None:
            raise HTTPException(404, "System call not found")
        with ctx.lock:
            profiles = {pid: ((name or "").split("__")[-1] or None, position) for pid, name, position in ctx.store.execute(
                "SELECT protein_id, profile_name, position FROM system_proteins WHERE system_id = ?", [system_id])}
            members = list(dict.fromkeys(list(call["proteins"]) + list(profiles)))
            located = ctx.store.execute(
                "SELECT protein_id, contig_id, start, end_coord FROM proteins "
                "WHERE protein_id IN (SELECT UNNEST(?::VARCHAR[]))", [members]) if members else []
            if not located:
                raise HTTPException(404, "System members not found in the protein table")
            contig = located[0][1]
            lo_member = min(r[2] for r in located if r[1] == contig)
            hi_member = max(r[3] for r in located if r[1] == contig)
            info = ctx.store.execute("SELECT length FROM contigs WHERE contig_id = ?", [contig])
            genes_on = ctx.store.execute("SELECT MIN(start), MAX(end_coord) FROM proteins WHERE contig_id = ?", [contig])
            contig_length = max(info[0][0] or 0, genes_on[0][1] or 0) if info else (genes_on[0][1] or 0)
            lo, hi = max(1, lo_member - flank), min(contig_length or hi_member + flank, hi_member + flank)
            rows = ctx.store.execute(
                """WITH p AS (SELECT protein_id, start, end_coord, strand, sequence_length FROM proteins
                              WHERE contig_id = ? AND end_coord >= ? AND start <= ?),
                        best AS (SELECT a.protein_id, ARG_MIN(COALESCE(NULLIF(a.name, ''), a.accession),
                                                              COALESCE(a.evalue, 1)) AS top
                                 FROM annotations a JOIN p USING (protein_id)
                                 WHERE LOWER(a.source) NOT IN ('defensefinder_system', 'txsscan_system',
                                                               'hyddb_subgroup')
                                 GROUP BY 1)
                   SELECT p.protein_id, p.start, p.end_coord, p.strand, p.sequence_length, best.top
                   FROM p LEFT JOIN best USING (protein_id) ORDER BY p.start, p.end_coord""", [contig, lo, hi])
            categories = ctx.top_categories([r[0] for r in rows])
        member_set = set(members)
        genes = []
        for pid, start, end, strand, length, top in rows:
            category = categories.get(pid)
            profile = profiles.get(pid, (None, None))[0]
            genes.append({"protein_id": pid, "start": start, "end": end, "strand": strand, "length_aa": length,
                          "annotation": ctx.describe_hit(top), "member": pid in member_set, "profile": profile,
                          "position": profiles.get(pid, (None, None))[1], "category": category,
                          "color": charts.CATEGORY_COLORS.get(category) if category else None})
        overlays: list[dict[str, Any]] = []
        positions = {g["protein_id"]: (g["start"], g["end"]) for g in genes}
        for other in catalog.systems:
            if other["system_id"] == system_id or other["bin_id"] != call["bin_id"]:
                continue
            inside = [positions[p] for p in other["proteins"] if p in positions]
            if inside:
                overlays.append({"label": other["type"], "kind": other["kind"],
                                 "start": min(a for a, _ in inside), "end": max(b for _, b in inside),
                                 "url": ctx.url("call", other["system_id"])})
        for locus in catalog.loci:
            if locus["contig_id"] == contig and locus["start"] is not None and locus["end"] is not None \
                    and locus["end"] >= lo and locus["start"] <= hi:
                overlays.append({"label": locus["type"], "kind": locus["type"], "start": locus["start"],
                                 "end": locus["end"]})
        overlays.insert(0, {"label": call["type"], "kind": call["kind"], "start": lo_member, "end": hi_member})
        member_rows = sorted((g for g in genes if g["member"]), key=lambda g: (g["position"] or 0, g["start"]))
        flanking = [g for g in genes if not g["member"]]
        g = catalog.by_bin.get(call["bin_id"])
        span = hi - lo + 1
        return ctx.render(
            request, "call.html", "systems", call=call, g=g, contig=contig, lo=lo, hi=hi, flank=flank,
            flanks=FLANKS, lo_member=lo_member, hi_member=hi_member, members=member_rows, flanking=flanking,
            others=[o for o in overlays[1:] if o.get("url")], other_loci=[o for o in overlays[1:] if not o.get("url")],
            track=charts.contig_track(genes, overlays, lo, hi, contig_start=lo <= 1,
                                      contig_end=bool(contig_length) and hi >= contig_length),
            contig_length=contig_length, viewer_start=max(1, lo_member - max(span, 20000) // 2 + (hi_member - lo_member) // 2),
            viewer_span=max(span, 20000), categories=sorted({x["category"] for x in genes if x["category"]}))
