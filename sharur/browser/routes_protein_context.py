"""Protein-page context panels, loaded after the page: architecture elsewhere, paralogs, usual neighbours.

Registered by :func:`sharur.browser.app.create_app` before the catch-all protein route.
See :mod:`sharur.browser.protein_context` for what each panel counts.
"""

from __future__ import annotations

import re
from types import SimpleNamespace

from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import HTMLResponse

from sharur.browser import charts
from sharur.browser import protein_context as pc

_PLAIN_NAME = re.compile(r"^[A-Za-z0-9_.\-]+$")


def exact_pattern(domains: list[dict]) -> str:
    """Architecture-search pattern matching exactly this domain order, e.g. ``^ Big_2{3} VWA $``."""
    runs: list[list] = []
    for d in domains:
        token = d["name"] if _PLAIN_NAME.match(d["name"] or "") else d["accession"].split(".")[0]
        if runs and runs[-1][0] == token:
            runs[-1][1] += 1
        else:
            runs.append([token, 1])
    return "^ " + " ".join(t if k == 1 else f"{t}{{{k}}}" for t, k in runs) + " $"


def register(app: FastAPI, ctx: SimpleNamespace) -> None:
    """Add ``/protein/{id}/context`` (an HTML fragment). ``ctx``: store, lock, catalog, render, ko_names."""
    store, catalog = ctx.store, ctx.catalog

    def ko_label(ko: str) -> tuple[str, str]:
        known = ctx.ko_names.get(ko) if getattr(ctx, "ko_names", None) else None
        if not known:
            return "", ""
        symbol = known[0].split(",")[0].strip()
        return symbol, re.sub(r"\s*\[EC:[^\]]*\]", "", known[1])

    @app.get("/protein/{protein_id:path}/context", response_class=HTMLResponse)
    def protein_context(request: Request, protein_id: str):
        from sharur.architecture import architecture  # noqa: PLC0415

        with ctx.lock:
            rows = store.execute("""SELECT p.bin_id, p.contig_id, p.sequence_length, COALESCE(p.partial, '00'), c.length
                                    FROM proteins p LEFT JOIN contigs c USING (contig_id) WHERE p.protein_id = ?""",
                                 [protein_id])
            if not rows:
                raise HTTPException(404, "Protein not found")
            bin_id, contig_id, length, partial, contig_length = rows[0]
            domains = [d.to_dict() for d in architecture(store, protein_id)]
            ko = pc.best_ko(store, protein_id)
            family = pc.family_of(store, protein_id, domains, ko)
            same = pc.same_architecture(store, catalog, protein_id, domains, length or 0, ko=ko)
            paralogs = pc.paralogs(store, protein_id, bin_id, ko=ko, domains=domains)
            hood = pc.neighbourhood(store, family) if family else None
            position = pc.contig_position(store, protein_id, contig_id)
        if family and family["kind"] == "ko":
            family["symbol"], family["definition"] = ko_label(family["id"])
        if hood:
            for row in hood["rows"]:
                if row["kind"] == "ko":
                    row["symbol"], row["definition"] = ko_label(row["id"])
        minimap = None
        if position:
            span = max(6000, (position["end"] - position["start"]) * 4)
            minimap = charts.contig_overview(contig_length or position["genes"][-1][1], position["start"],
                                             position["end"], position["genes"], ctx.url("contig", contig_id), span)
            position["to_start"] = position["start"] - 1
            position["to_end"] = max(0, (contig_length or position["end"]) - position["end"])
            position["span"] = span
        return ctx.render(request, "_protein_context.html", "taxa", protein_id=protein_id, bin_id=bin_id,
                          contig_id=contig_id, length=length, partial=partial != "00", domains=domains,
                          pattern=exact_pattern(domains) if domains else "", ko=ko,
                          ko_symbol=ko_label(ko)[0] if ko else "", family=family, same=same, paralogs=paralogs,
                          hood=hood, position=position, minimap=minimap)
