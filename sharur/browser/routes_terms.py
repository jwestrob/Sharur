"""V2 semantic-term search and per-protein stored term rows.

``/terms`` finds proteins by V2 term membership: all of (AND), any of (OR), none of (NOT),
optionally within a genome or clade. Membership is the materialized term table's active
rows, ``(term_kind != 'atom' OR relation != 'excludes')``; results page in canonical
protein-ID order with an exact total. ``/protein/{id}/terms`` lists every stored term
row for one protein exactly as stored. Both read through the dataset's semantic provider:
SQL over ``semantic_terms`` by default, or a compact index generation chosen at launch
(``sharur browse --semantic-index``). Every response names the backend that served it.

V2 terms and the browser's function labels are separate vocabularies; these pages show
V2 term IDs and their stored rows as recorded.
"""

from __future__ import annotations

import re
import time
from contextlib import nullcontext
from typing import Any
from urllib.parse import urlencode

from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.responses import HTMLResponse, JSONResponse

from sharur.browser.routes_search import resolve_scope
from sharur.semantic_index.provider import MAX_TERMS, SearchRequest


PAGE = 50
MAX_LIMIT = 500
SPLIT = re.compile(r"[\s,]+")
FIELDS = ("term_id", "term_kind", "facet", "relation", "source_db", "source_accession")


def split_terms(values: list[str]) -> list[str]:
    """Term IDs from repeated parameters, each holding IDs separated by commas or whitespace."""
    return [t for v in values for t in SPLIT.split(v.strip()) if t]


def is_active(term_kind: str | None, relation: str | None) -> bool:
    """SQL ``(term_kind != 'atom' OR relation != 'excludes')`` with NULL semantics."""
    return (term_kind is not None and term_kind != "atom") or (relation is not None and relation != "excludes")


def _guard(ctx):
    return ctx.lock if ctx.semantic.needs_store_lock else nullcontext()


def backend_label(provider) -> str:
    identity = provider.identity()
    return f"{identity['backend']} {identity['generation_id']}" if identity.get("generation_id") else identity["backend"]


def _stamp(response, provider, timing: dict[str, float] | None = None):
    response.headers["X-Sharur-Semantic-Backend"] = backend_label(provider)
    if timing:
        response.headers["Server-Timing"] = ", ".join(f"{k};dur={v:.3f}" for k, v in timing.items())
    return response


def _resolve_genomes(ctx, request: Request, scope: str) -> tuple[tuple[str, ...] | None, dict | None, str | None]:
    """Genome scope from ``scope`` (a genome or clade name) or explicit ``genomes`` IDs; None means unscoped."""
    explicit = request.query_params.getlist("genomes")
    if scope.strip() and explicit:
        return None, None, "Give either a scope name or explicit genomes, not both."
    if scope.strip():
        found = resolve_scope(ctx, scope)
        if found is None:
            return None, None, f"No genome or clade matches “{scope.strip()}”."
        kind, label, bins = found
        return tuple(bins), {"kind": kind, "label": label, "genomes": len(bins)}, None
    if explicit:
        genomes = tuple(split_terms(explicit))
        return genomes, {"kind": "genomes", "label": ", ".join(genomes[:3]) + ("…" if len(genomes) > 3 else ""),
                         "genomes": len(genomes)}, None
    return None, None, None


def run_search(ctx, request: SearchRequest, offset: int, limit: int) -> dict[str, Any]:
    """One page of term search: IDs from the provider, matched stored rows, scalar metadata from the dataset."""
    provider = ctx.semantic
    started = time.perf_counter()
    with _guard(ctx):
        page = provider.search(request, offset=offset, limit=limit)
        t1 = time.perf_counter()
        positives = set(request.has + request.any_of)
        evidence = provider.rich_rows_many(page.protein_ids) if positives and page.protein_ids else {}
    t2 = time.perf_counter()
    info = {}
    if page.protein_ids:
        with ctx.lock:
            info = {r[0]: r for r in ctx.store.execute(
                "SELECT protein_id, bin_id, contig_id, start, end_coord, strand, sequence_length FROM proteins "
                "WHERE protein_id IN (SELECT UNNEST(?::VARCHAR[]))", [page.protein_ids])}
    t3 = time.perf_counter()
    rows = []
    for pid in page.protein_ids:
        _, bin_id, contig, start, end, strand, length = info.get(pid, (pid, None, None, None, None, None, None))
        matched, seen = [], set()
        for term_id, kind, facet, relation, source_db, accession in evidence.get(pid, ()):
            key = (term_id, relation, source_db, accession)
            if term_id in positives and is_active(kind, relation) and key not in seen:
                seen.add(key)
                matched.append({"term_id": term_id, "relation": relation, "facet": facet,
                                "source_db": source_db, "source_accession": accession})
        rows.append({"protein_id": pid, "bin_id": bin_id, "contig_id": contig, "start": start, "end": end,
                     "strand": strand, "length": length, "matched": matched})
    with _guard(ctx):
        counts = provider.term_counts(request.terms())
    timing = {"search": page.elapsed_ms, "evidence": (t2 - t1) * 1000, "hydrate": (t3 - t2) * 1000,
              "total": (time.perf_counter() - started) * 1000}
    return {"total": page.total, "offset": offset, "limit": limit, "rows": rows, "term_counts": counts,
            "timing_ms": timing}


def protein_terms(ctx, protein_id: str) -> dict[str, Any] | None:
    with ctx.lock:
        found = ctx.store.execute("SELECT bin_id FROM proteins WHERE protein_id = ?", [protein_id])
    if not found:
        return None
    started = time.perf_counter()
    with _guard(ctx):
        rows = ctx.semantic.rich_rows(protein_id)
        elapsed = (time.perf_counter() - started) * 1000
        counts = ctx.semantic.term_counts(dict.fromkeys(r[0] for r in rows))
    groups: list[dict[str, Any]] = []
    for row in rows:
        if not groups or groups[-1]["term_id"] != row[0] or groups[-1]["term_kind"] != row[1]:
            groups.append({"term_id": row[0], "term_kind": row[1], "rows": [], "active": False,
                           "proteins": counts.get(row[0], 0)})
        groups[-1]["rows"].append(row)
        groups[-1]["active"] |= is_active(row[1], row[3])
    return {"protein_id": protein_id, "bin_id": found[0][0], "rows": rows, "groups": groups,
            "timing_ms": {"rich_rows": elapsed}}


def register(app: FastAPI, ctx) -> None:
    def parse(has, any_, lacks, request: Request, scope: str):
        has_t, any_t, lacks_t = split_terms(has), split_terms(any_), split_terms(lacks)
        genomes, scope_info, error = _resolve_genomes(ctx, request, scope)
        if error:
            return None, scope_info, error
        try:
            return SearchRequest.build(has_t, any_t, lacks_t, genomes), scope_info, None
        except ValueError as exc:
            return None, scope_info, str(exc)

    @app.get("/terms", response_class=HTMLResponse)
    def terms_page(request: Request, has: list[str] = Query([]), any: list[str] = Query([]),
                   lacks: list[str] = Query([]), scope: str = Query("", max_length=500),
                   offset: int = Query(0, ge=0), limit: int = Query(PAGE, ge=1, le=MAX_LIMIT)):
        req, scope_info, error = parse(has, any, lacks, request, scope)
        result = None
        if req is not None and not req.empty:
            result = run_search(ctx, req, offset, limit)
        names = {}
        if result:
            for term in req.terms():
                names[term] = term_family_name(ctx, term)
            for r in result["rows"]:
                g = ctx.catalog.by_bin.get(r["bin_id"])
                r["lineage"] = g.label if g else ""
        pages = {}
        if result:
            def page_url(off):
                params = [(k, v) for k, v in request.query_params.multi_items() if k != "offset"]
                return "/terms?" + urlencode([*params, ("offset", str(off))])
            pages = {"prev": page_url(max(0, offset - limit)) if offset else None,
                     "next": page_url(offset + limit) if offset + limit < result["total"] else None}
        response = ctx.render(request, "terms.html", "terms", pages=pages, form={"has": " ".join(split_terms(has)),
                              "any": " ".join(split_terms(any)), "lacks": " ".join(split_terms(lacks)),
                              "scope": scope}, req=req, result=result, scope_info=scope_info, error=error,
                              names=names, backend=ctx.semantic.identity(), max_terms=MAX_TERMS,
                              query_string=str(request.query_params))
        return _stamp(response, ctx.semantic, result["timing_ms"] if result else None)

    @app.get("/protein/{protein_id:path}/terms", response_class=HTMLResponse)
    def protein_terms_page(request: Request, protein_id: str):
        out = protein_terms(ctx, protein_id)
        if out is None:
            raise HTTPException(404, "Protein not found")
        response = ctx.render(request, "protein_terms.html", "taxa", t=out, g=ctx.catalog.by_bin.get(out["bin_id"]),
                              backend=ctx.semantic.identity(),
                              names={grp["term_id"]: term_family_name(ctx, grp["term_id"]) for grp in out["groups"]})
        return _stamp(response, ctx.semantic, out["timing_ms"])

    @app.get("/api/v1/terms")
    def api_terms(request: Request, has: list[str] = Query([]), any: list[str] = Query([]),
                  lacks: list[str] = Query([]), scope: str = Query("", max_length=500),
                  offset: int = Query(0, ge=0), limit: int = Query(PAGE, ge=1, le=MAX_LIMIT)):
        req, scope_info, error = parse(has, any, lacks, request, scope)
        if error:
            raise HTTPException(400, error)
        result = run_search(ctx, req, offset, limit)
        body = {"backend": ctx.semantic.identity(),
                "request": {"has": list(req.has), "any": list(req.any_of), "lacks": list(req.lacks),
                            # a named scope is described by kind/label/count; explicit IDs are echoed
                            "genomes": list(req.genomes) if scope_info and scope_info["kind"] == "genomes" else None,
                            "scope": scope_info},
                **result}
        return _stamp(JSONResponse(body), ctx.semantic, result["timing_ms"])

    @app.get("/api/v1/terms/suggest")
    def api_terms_suggest(q: str = Query("", max_length=200), limit: int = Query(20, ge=1, le=100)):
        with _guard(ctx):
            found = ctx.semantic.suggest_terms(q, limit) if len(q.strip()) >= 2 else []
        return _stamp(JSONResponse([{"term_id": t, "proteins": n} for t, n in found]), ctx.semantic)

    @app.get("/api/v1/protein/{protein_id:path}/terms")
    def api_protein_terms(protein_id: str):
        out = protein_terms(ctx, protein_id)
        if out is None:
            raise HTTPException(404, "Protein not found")
        body = {"backend": ctx.semantic.identity(), "protein_id": protein_id, "fields": list(FIELDS),
                "count": len(out["rows"]), "rows": [list(r) for r in out["rows"]], "timing_ms": out["timing_ms"]}
        return _stamp(JSONResponse(body), ctx.semantic, out["timing_ms"])

    @app.get("/api/v1/semantic")
    def api_semantic():
        provider = ctx.semantic
        body = {"backend": provider.identity(), "calls": provider.counters.snapshot()}
        generation = getattr(provider, "generation", None)
        if generation is not None:
            body["counts"] = generation.record["counts"]
            body["adopted_at"] = generation.record.get("adopted_at")
        return _stamp(JSONResponse(body), provider)


def term_family_name(ctx, term: str) -> str:
    """The annotation family name behind a ``pfam:`` or ``kegg:`` term ID (the family's own name only)."""
    source, _, accession = term.partition(":")
    if source == "pfam":
        return (ctx.catalog.domains.get(accession) or {}).get("name", "")
    if source == "kegg" and ctx.ko_names is not None:
        known = ctx.ko_names.get(accession)
        return known[0].split(",")[0].strip() if known else ""
    return ""
