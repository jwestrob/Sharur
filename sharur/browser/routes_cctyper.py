"""CRISPR-Cas subtype calls from ``sharur cas-type`` (CRISPRCasTyper scoring).

Calls live in ``crispr_cas_systems`` (one row per Cas operon) and array repeat
predictions in ``crispr_array_types``; both arrive with schema v9. Statuses
``crispr_cas`` and ``cas`` are the caller's confident subtype calls; the
``*_putative`` rows are candidates and are shown as such. Datasets without the
tables (or without rows) show nothing from this module.
"""

from __future__ import annotations

import re
from collections import Counter
from typing import Any

from fastapi import HTTPException, Query, Request
from fastapi.responses import HTMLResponse

from sharur.browser import routes_crispr

CONFIDENT = ("crispr_cas", "cas")
STATUS_LABELS = {
    "crispr_cas": "CRISPR-Cas (operon + array)",
    "cas": "Cas operon, no array within 10 kb",
    "crispr_cas_putative": "candidate: array nearby, subtype unresolved",
    "cas_putative": "candidate: too little Cas evidence or ambiguous",
}


def _split(value: str | None) -> list[str]:
    return [x for x in (value or "").split(",") if x]


def gene_name(profile: str | None) -> str:
    """Gene a CRISPRCasTyper profile belongs to: ``Cas7_0_IB`` -> ``Cas7``."""
    return (profile or "").split("_", 1)[0]


def register(app, ctx) -> None:
    catalog = ctx.catalog

    def _table(name: str) -> bool:
        with ctx.lock:
            return bool(ctx.store.execute(
                "SELECT COUNT(*) FROM information_schema.tables WHERE table_name = ?", [name])[0][0])

    def systems() -> list[dict[str, Any]]:
        if hasattr(ctx, "_cctyper_systems"):
            return ctx._cctyper_systems
        out: list[dict[str, Any]] = []
        if _table("crispr_cas_systems"):
            with ctx.lock:
                rows = ctx.store.execute(
                    """SELECT system_id, genome_id, contig_id, start, end_coord, status, prediction, prediction_cas,
                              best_type, best_score, complete_interference, complete_adaptation, genes_count,
                              protein_ids, profile_names, crispr_locus_ids, crispr_distances, caller
                       FROM crispr_cas_systems ORDER BY genome_id, contig_id, start""")
            for (sid, bin_id, contig, start, end, status, prediction, prediction_cas, best_type, best_score,
                 interf, adapt, n_genes, proteins, profiles, loci, distances, caller) in rows:
                g = catalog.by_bin.get(bin_id)
                out.append({
                    "system_id": sid, "bin_id": bin_id, "lineage": g.label if g else "", "contig_id": contig,
                    "start": start, "end": end, "status": status, "confident": status in CONFIDENT,
                    "prediction": prediction, "prediction_cas": prediction_cas, "best_type": best_type,
                    "best_score": best_score, "complete_interference": interf, "complete_adaptation": adapt,
                    "genes_count": n_genes, "proteins": _split(proteins), "profiles": _split(profiles),
                    "arrays": _split(loci), "distances": _split(distances), "caller": caller})
        ctx._cctyper_systems = out
        return out

    def array_types() -> dict[str, dict[str, Any]]:
        if hasattr(ctx, "_cctyper_arrays"):
            return ctx._cctyper_arrays
        out: dict[str, dict[str, Any]] = {}
        if _table("crispr_array_types"):
            with ctx.lock:
                rows = ctx.store.execute(
                    """SELECT locus_id, subtype, probability, prediction, repeat_identity, spacer_identity,
                              spacer_sem, trusted, near_cas FROM crispr_array_types""")
            for locus_id, subtype, prob, prediction, r_id, s_id, sem, trusted, near in rows:
                out[locus_id] = {"subtype": subtype, "probability": prob, "prediction": prediction,
                                 "repeat_identity": r_id, "spacer_identity": s_id, "spacer_sem": sem,
                                 "trusted": trusted, "near_cas": near}
        ctx._cctyper_arrays = out
        return out

    def systems_for_array(locus_id: str) -> list[dict[str, Any]]:
        return [s for s in systems() if locus_id in s["arrays"]]

    def summary() -> dict[str, Any]:
        """Systems-page rows: confident subtypes, then the candidate statuses."""
        calls = systems()
        n = len(catalog.genomes) or 1
        by_subtype: dict[str, list[dict[str, Any]]] = {}
        for s in calls:
            if s["confident"]:
                by_subtype.setdefault(s["prediction"], []).append(s)
        subtypes = sorted(({"subtype": k, "count": len(v), "genomes": len({s["bin_id"] for s in v}),
                            "share": len({s["bin_id"] for s in v}) / n} for k, v in by_subtype.items()),
                          key=lambda r: -r["count"])
        putative = []
        for status in ("crispr_cas_putative", "cas_putative"):
            members = [s for s in calls if s["status"] == status]
            if members:
                putative.append({"status": status, "label": STATUS_LABELS[status], "count": len(members),
                                 "share": len({s["bin_id"] for s in members}) / n})
        return {"subtypes": subtypes, "putative": putative, "total": len(calls)}

    ctx.cctyper_systems = systems
    ctx.cctyper_array_types = array_types
    ctx.cctyper_systems_for_array = systems_for_array
    ctx.cctyper_summary = summary

    def _member_genes(call: dict[str, Any], lo: int, hi: int) -> list[dict[str, Any]]:
        profiles = dict(zip(call["proteins"], call["profiles"]))
        genes = ctx.genes_near(call["contig_id"], lo, hi)
        for g in genes:
            profile = profiles.get(g["protein_id"])
            g["member"] = profile is not None
            g["profile"] = profile
            g["cas_name"] = gene_name(profile) if profile else None
            # the context drawing highlights members only: Cas-domain hits outside the call stay plain
            g["cas"] = g["member"]
            if profile:
                g["label"] = gene_name(profile)
        return genes

    def _linked_arrays(call: dict[str, Any]) -> list[dict[str, Any]]:
        ctx.crispr_arrays()
        by_id = getattr(ctx, "_crispr_by_id", {})
        return [dict(by_id[a], current=True) for a in call["arrays"] if a in by_id]

    @app.get("/crispr-cas/calls", response_class=HTMLResponse)
    def cas_calls(request: Request, subtype: str = Query(""), status: str = Query(""),
                  page: int = Query(1, ge=1), flank: int = Query(3000, ge=0, le=20000)):
        calls = systems()
        if status:
            selected = [s for s in calls if s["status"] == status]
            title = STATUS_LABELS.get(status, status)
        else:
            selected = [s for s in calls if s["confident"] and (not subtype or s["prediction"] == subtype)]
            title = subtype or "all subtypes"
        per_page = 25
        pages = max(1, (len(selected) + per_page - 1) // per_page)
        page = min(page, pages)
        shown = selected[(page - 1) * per_page: page * per_page]
        windows = []
        for s in shown:
            arrays = _linked_arrays(s)
            lo = min([s["start"]] + [a["start"] for a in arrays]) - flank
            hi = max([s["end"]] + [a["end"] for a in arrays]) + flank
            windows.append((s, arrays, max(1, lo), hi))
        span = max((hi - lo for *_, lo, hi in windows), default=1)
        rows = []
        for s, arrays, lo, hi in windows:
            center = (lo + hi) // 2
            genes = _member_genes(s, center - span // 2, center + span // 2)
            members = [g for g in genes if g["member"]]
            key = next((g for g in members if re.fullmatch(r"cas1", g["cas_name"] or "", re.I)),
                       members[0] if members else None)
            flip = bool(key and key["strand"] == "-")
            svg = ctx.locus_row_svg({"arrays": arrays}, genes, center - span // 2, center + span // 2, flip)
            rows.append((s, svg, flip))
        subtypes = Counter(s["prediction"] for s in calls if s["confident"])
        return ctx.render(request, "cas_calls.html", "systems", rows=rows, subtype=subtype, status=status,
                          title=title, total=len(selected), page=page, pages=pages, flank=flank,
                          subtypes=subtypes.most_common(), confident_total=sum(subtypes.values()), status_labels=STATUS_LABELS,
                          statuses=Counter(s["status"] for s in calls))

    @app.get("/cas-system/{system_id:path}", response_class=HTMLResponse)
    def cas_system(request: Request, system_id: str, flank: int = Query(8000, ge=1000, le=60000)):
        call = next((s for s in systems() if s["system_id"] == system_id), None)
        if call is None:
            raise HTTPException(404, "CRISPR-Cas call not found")
        arrays = _linked_arrays(call)
        lo = max(1, min([call["start"]] + [a["start"] for a in arrays]) - flank)
        hi = max([call["end"]] + [a["end"] for a in arrays]) + flank
        with ctx.lock:
            length = ctx.store.execute("SELECT length FROM contigs WHERE contig_id = ?", [call["contig_id"]])
            scores = {pid: score for pid, score in ctx.store.execute(
                "SELECT protein_id, score FROM system_proteins WHERE system_id = ? AND system_source = 'cctyper'",
                [system_id])}
        contig_length = length[0][0] if length and length[0][0] else None
        if contig_length:
            hi = min(hi, contig_length)
        genes = _member_genes(call, lo, hi)
        members = [dict(g, score=scores.get(g["protein_id"])) for g in genes if g["member"]]
        array_info = []
        types = array_types()
        for a, d in zip(arrays, call["distances"] + [""] * len(arrays)):
            array_info.append(dict(a, distance=d, typing=types.get(a["locus_id"])))
        g = catalog.by_bin.get(call["bin_id"])
        return ctx.render(request, "cas_system.html", "systems", call=call, g=g, members=members,
                          arrays=array_info, lo=lo, hi=hi, flank=flank, status_label=STATUS_LABELS.get(call["status"]),
                          context=routes_crispr.context_svg(genes, arrays, lo, hi),
                          at_contig_start=lo <= 1, at_contig_end=bool(contig_length) and hi >= contig_length)
