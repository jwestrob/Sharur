"""A read-only JSON API over the browsed dataset, and a page that documents it.

Every endpoint answers GET with JSON; lists take ``limit`` and ``offset``.
Proteins come without sequences unless ``sequence=1``. The same data backs
the HTML pages, so a script and a page agree.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.responses import HTMLResponse, JSONResponse

from sharur.browser.catalog import RANKS, UNCLASSIFIED

MAX_LIMIT = 5000


def _clean(value: Any) -> Any:
    """JSON-safe copy: numpy scalars and arrays to Python, NaN to None, sets to sorted lists."""
    if isinstance(value, dict):
        return {str(k): _clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_clean(v) for v in value]
    if isinstance(value, set):
        return sorted(_clean(v) for v in value)
    if isinstance(value, np.ndarray):
        return _clean(value.tolist())
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return None
    return value


def _genome(g) -> dict[str, Any]:
    return {"bin_id": g.bin_id, "label": g.label, "taxonomy": {r: g.taxonomy.get(r) for r in RANKS},
            "completeness": g.completeness, "contamination": g.contamination, "contigs": g.contigs,
            "length": g.length, "n50": g.n50, "proteins": g.proteins, "annotated": g.annotated,
            "longest_protein": g.longest_protein}


def _page(rows: list, limit: int, offset: int) -> dict[str, Any]:
    limit = max(1, min(limit, MAX_LIMIT))
    return {"total": len(rows), "offset": offset, "limit": limit, "rows": rows[offset:offset + limit]}


def endpoints(catalog) -> list[dict[str, str]]:
    g = catalog.genomes[0] if catalog.genomes else None
    clade = next(((r, g.taxonomy[r]) for r in ("order", "class", "phylum") if g and g.taxonomy.get(r)
                  and g.taxonomy[r] != UNCLASSIFIED), None)
    gid = g.bin_id if g else "GENOME"
    clade_path = f"{clade[0]}/{clade[1]}" if clade else "order/NAME"
    clade_token = f"{clade[0]}:{clade[1]}" if clade else "order:NAME"
    return [
        {"path": "/api/v1", "what": "This list, as JSON."},
        {"path": "/api/v1/genomes", "what": "Genomes with taxonomy and assembly statistics. Filters: clade=rank:name, "
                                            "min_completeness, max_contamination.",
         "example": f"/api/v1/genomes?clade={clade_token}&limit=5"},
        {"path": "/api/v1/genome/{bin_id}", "what": "One genome: statistics, systems, CRISPR-Cas calls, pathways at "
                                                    "least half complete.", "example": f"/api/v1/genome/{gid}"},
        {"path": "/api/v1/protein/{protein_id}", "what": "A protein card: location, genome, annotations by source, "
                                                         "labels, systems, neighbourhood. sequence=1 adds the sequence."},
        {"path": "/api/v1/clade/{rank}/{name}", "what": "A clade: genomes, child clades, gene content and signature "
                                                        "families.", "example": f"/api/v1/clade/{clade_path}"},
        {"path": "/api/v1/matrix", "what": "The presence/absence matrix, same parameters as the /matrix page.",
         "example": f"/api/v1/matrix?a={clade_token}&kind=ko&n=10"},
        {"path": "/api/v1/discover/{list}", "what": "A Discover list: giants, dark, repeats, fusions, islands, "
                                                     "clades, systems.", "example": "/api/v1/discover/fusions?limit=5"},
        {"path": "/api/v1/search", "what": "Search suggestions; 'term in clade-or-genome' runs a scoped protein "
                                           "search.", "example": "/api/v1/search?q=rubisco"},
    ]


def register(app: FastAPI, ctx) -> None:
    catalog = ctx.catalog

    @app.get("/api", response_class=HTMLResponse)
    def api_page(request: Request):
        return ctx.render(request, "api.html", "discover", endpoints=endpoints(catalog))

    @app.get("/api/v1")
    def api_index():
        return JSONResponse({"dataset": getattr(ctx, "dataset", None), "endpoints": endpoints(catalog)})

    @app.get("/api/v1/genomes")
    def api_genomes(clade: str = Query("", max_length=500), min_completeness: float = Query(0.0, ge=0, le=100),
                    max_contamination: float = Query(100.0, ge=0), limit: int = Query(500), offset: int = Query(0, ge=0)):
        genomes = catalog.genomes
        if clade:
            if ":" not in clade or clade.split(":", 1)[0] not in RANKS:
                raise HTTPException(400, "clade takes rank:name, e.g. order:Micrarchaeales")
            rank, name = clade.split(":", 1)
            genomes = catalog.clade(rank, name)
        rows = [_genome(g) for g in genomes
                if (g.completeness or 0) >= min_completeness
                and (g.contamination is None or g.contamination <= max_contamination)]
        return JSONResponse(_clean(_page(rows, limit, offset)))

    @app.get("/api/v1/genome/{bin_id:path}")
    def api_genome(bin_id: str):
        g = catalog.by_bin.get(bin_id)
        if g is None:
            raise HTTPException(404, "Genome not found")
        out = _genome(g)
        out["systems"] = [{k: s[k] for k in ("kind", "system_id", "type", "subtype", "genes", "proteins")}
                          for s in catalog.systems if s["bin_id"] == bin_id]
        calls = ctx.cctyper_systems() if hasattr(ctx, "cctyper_systems") else []
        out["crispr_cas"] = [{k: c[k] for k in ("system_id", "contig_id", "start", "end", "status", "prediction",
                                                "confident", "proteins", "arrays")}
                             for c in calls if c["bin_id"] == bin_id]
        out["loci"] = [{k: v for k, v in locus.items() if k != "bin_id"} for locus in catalog.loci
                       if locus["bin_id"] == bin_id]
        modules = []
        if catalog.module_completeness is not None:
            row = catalog.module_completeness[g.index]
            for j in np.argsort(-row):
                if row[j] < 0.5:
                    break
                d = catalog.modules[catalog.module_ids[j]]
                modules.append({"module": d.module, "name": d.name, "completeness": float(row[j])})
        out["modules"] = modules
        return JSONResponse(_clean(out))

    @app.get("/api/v1/protein/{protein_id:path}")
    def api_protein(protein_id: str, sequence: int = Query(0), window: int = Query(5, ge=0, le=20)):
        from sharur.operators.cards import card  # noqa: PLC0415

        with ctx.lock:
            c = card(ctx.store, protein_id, window=window)
            if c.get("found") and sequence:
                rows = ctx.store.execute("SELECT sequence FROM proteins WHERE protein_id = ?", [protein_id])
                c["sequence"] = rows[0][0] if rows else None
        if not c.get("found"):
            raise HTTPException(404, "Protein not found")
        return JSONResponse(_clean(c))

    @app.get("/api/v1/clade/{rank}/{name:path}")
    def api_clade(rank: str, name: str):
        if rank not in RANKS:
            raise HTTPException(404, "Unknown rank")
        genomes = catalog.clade(rank, name)
        if not genomes:
            raise HTTPException(404, "Clade not found")
        child, counts = catalog.children(genomes, rank)
        out: dict[str, Any] = {"rank": rank, "name": name, "genomes": len(genomes),
                               "lineage": catalog.lineage_of(rank, name),
                               "children": {"rank": child, "counts": dict(counts)}}
        sets = getattr(ctx, "feature_sets", None)
        if sets is not None:
            from sharur.browser.clade_content import clade_content  # noqa: PLC0415

            try:
                out["gene_content"] = clade_content(ctx, sets, genomes)
            except LookupError:
                out["gene_content"] = None
        if catalog.module_completeness is not None:
            out["modules"] = catalog.clade_modules(genomes)[:50]
        return JSONResponse(_clean(out))

    @app.get("/api/v1/matrix")
    def api_matrix(a: str = Query("", max_length=500), b: str = Query("", max_length=500), kind: str = Query("ko"),
                   features: str = Query("", max_length=4000), module: str = Query("", max_length=20),
                   pick: str = Query(""), n: int = Query(40), order: str = Query("taxonomy"),
                   min_completeness: float = Query(0.0, ge=0, le=100)):
        compute = getattr(ctx, "matrix_compute", None)
        if compute is None:
            raise HTTPException(404, "Matrix unavailable")
        _, result, message = compute(a, b, kind, features, module, pick, n, order, min_completeness)
        if not result or message:
            raise HTTPException(400, message or "Choose a genome or clade (a=...)")
        out = {k: v for k, v in result.items() if k != "values"}
        out["values"] = result["values"]
        return JSONResponse(_clean(out))

    @app.get("/api/v1/discover/{feed}")
    def api_discover(feed: str, limit: int = Query(100), offset: int = Query(0, ge=0)):
        notable = catalog.notable or {}
        if feed == "systems":
            from sharur.browser import discover  # noqa: PLC0415

            rows = discover.rare_systems(catalog, getattr(ctx, "cctyper_systems", None))
        elif feed == "clades":
            rows = (notable.get("giant_clades") or {}).get("rows", [])
        elif feed in ("giants", "dark", "repeats", "fusions", "islands"):
            if feed not in notable:
                raise HTTPException(503, "Discovery lists are still being computed")
            rows = notable[feed]
        else:
            raise HTTPException(404, "Unknown list")
        return JSONResponse(_clean(_page(list(rows), limit, offset)))

    @app.get("/api/v1/search")
    def api_search(q: str = Query("", max_length=500)):
        from sharur.browser.routes_search import SCOPED, resolve_scope, scoped_search  # noqa: PLC0415

        q = q.strip()
        m = SCOPED.match(q)
        if m:
            scope = resolve_scope(ctx, m.group(2))
            if scope:
                kind, label, bins = scope
                found = scoped_search(ctx, m.group(1), bins)   # takes the store lock itself
                return JSONResponse(_clean({"query": q, "scope": {"kind": kind, "label": label, "genomes": len(bins)},
                                            **found}))
        suggest = getattr(ctx, "suggest", None)
        return JSONResponse(_clean({"query": q, "suggestions": suggest(q, limit=50) if suggest and len(q) >= 2 else []}))
