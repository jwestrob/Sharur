"""Routes for the dataset browser."""

from __future__ import annotations

import secrets
import threading
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any
from urllib.parse import quote

import numpy as np
from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.responses import HTMLResponse, JSONResponse, RedirectResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

from sharur.browser import charts
from sharur.browser.catalog import (
    BIOLOGICAL,
    CATEGORY_LABELS,
    RANKS,
    UNCLASSIFIED,
    Catalog,
    load_background,
    load_catalog,
)
from sharur.predicates.vocabulary import PREDICATE_BY_ID
from sharur.storage.duckdb_store import DuckDBStore

HERE = Path(__file__).parent
COOKIE = "sharur_browse"


def _url(*parts: str) -> str:
    return "/" + "/".join(quote(str(p), safe="") for p in parts)


def _bp(n: int | float | None) -> str:
    if not n:
        return "0 bp"
    for unit, size in (("Gb", 1e9), ("Mb", 1e6), ("kb", 1e3)):
        if n >= size:
            return f"{n / size:.1f} {unit}"
    return f"{int(n)} bp"


def _num(n: Any) -> str:
    if n is None:
        return "–"
    if isinstance(n, float) and not n.is_integer():
        return f"{n:,.1f}"
    return f"{int(n):,}"


def _pct(x: float | None, digits: int = 0) -> str:
    return "–" if x is None else f"{x * 100:.{digits}f}%"


def _short(pid: str) -> str:
    """Readable protein label: the part after a genome prefix, middle-elided when long."""
    tail = pid.rsplit("|", 1)[-1]
    return tail if len(tail) <= 34 else tail[:16] + "…" + tail[-14:]


def _evalue(x: Any) -> str:
    return f"{x:.1e}" if isinstance(x, float) else ("–" if x is None else str(x))


def create_app(db_path: str | Path, *, token: str | None = None, background: bool = True) -> FastAPI:
    """Browser over one dataset, opened read-only."""
    from sharur.abundance import default_path  # noqa: PLC0415

    store = DuckDBStore(str(db_path), read_only=True)
    lock = threading.Lock()
    catalog = load_catalog(store)
    if background:
        threading.Thread(target=load_background, args=(store, catalog, lock), daemon=True).start()
    else:
        load_background(store, catalog, lock)
    sidecar = default_path(db_path)
    dataset_name = Path(db_path).resolve().parent.name

    templates = Jinja2Templates(directory=str(HERE / "templates"))
    templates.env.globals.update(url=_url, bp=_bp, num=_num, pct=_pct, evalue=_evalue, quote=quote, short=_short,
                                 PREDICATE_BY_ID=PREDICATE_BY_ID,
                                 charts=charts, CATEGORY_LABELS=CATEGORY_LABELS, dataset=dataset_name,
                                 catalog=catalog, CATEGORY_COLORS=charts.CATEGORY_COLORS)

    app = FastAPI(title="Sharur browser", docs_url=None, redoc_url=None, openapi_url=None)
    app.mount("/static", StaticFiles(directory=str(HERE / "static")), name="static")
    app.state.store, app.state.catalog = store, catalog

    def render(request: Request, template: str, section: str, /, **context: Any) -> HTMLResponse:
        return templates.TemplateResponse(request, template, {"section": section, **context})

    @app.middleware("http")
    async def require_token(request: Request, call_next):
        if token is None:
            return await call_next(request)
        supplied = request.query_params.get("token")
        if supplied is not None and secrets.compare_digest(supplied, token):
            response = RedirectResponse(str(request.url.remove_query_params("token")), status_code=303)
            response.set_cookie(COOKIE, token, httponly=True, samesite="strict")
            return response
        cookie = request.cookies.get(COOKIE)
        if cookie is None or not secrets.compare_digest(cookie, token):
            return HTMLResponse("Open the shared link (with ?token=…) to use this browser.", status_code=401)
        return await call_next(request)

    # ------------------------------------------------------------------ #
    # Overview and taxonomy
    # ------------------------------------------------------------------ #

    def clade_treemap(genomes, rank):
        child_rank, counts = catalog.children(genomes, rank)
        if child_rank is None:
            return None, []
        items = [{"name": name, "value": n, "sub": f"{n:,} genomes", "url": _url("taxa", child_rank, name)}
                 for name, n in counts]
        return child_rank, items

    @app.get("/", response_class=HTMLResponse)
    def home(request: Request):
        child_rank, items = clade_treemap(catalog.genomes, None)
        systems = Counter(s["kind"] for s in catalog.systems)
        loci = Counter(l["type"] for l in catalog.loci)
        return render(request, "home.html", "home", child_rank=child_rank, tree=charts.treemap(items, 340),
                      children=items,
                      systems=systems, loci=loci,
                      notable=catalog.notable.get("giants", [])[:6])

    @app.get("/taxa", response_class=HTMLResponse)
    @app.get("/taxa/{rank}/{name:path}", response_class=HTMLResponse)
    def taxa(request: Request, rank: str | None = None, name: str | None = None):
        if rank is not None and rank not in RANKS:
            raise HTTPException(404, "Unknown rank")
        genomes = catalog.clade(rank, name)
        if not genomes:
            raise HTTPException(404, "No genomes in this clade")
        child_rank, items = clade_treemap(genomes, rank)
        enriched = catalog.enrichment(genomes)
        sizes = [g.length for g in genomes if g.length]
        return render(request, "taxa.html", "taxa", rank=rank, name=name, genomes=genomes,
                      lineage=catalog.lineage_of(rank, name) if rank else [], child_rank=child_rank,
                      tree=charts.treemap(items, 300) if len(items) > 1 else None, children=items,
                      enriched=enriched, median_size=float(np.median(sizes)) if sizes else 0,
                      median_proteins=float(np.median([g.proteins for g in genomes])),
                      systems=_clade_systems(genomes))

    def _clade_systems(genomes):
        members = {g.bin_id for g in genomes}
        counts = Counter((s["kind"], s["type"]) for s in catalog.systems if s["bin_id"] in members)
        carriers = defaultdict(set)
        for s in catalog.systems:
            if s["bin_id"] in members:
                carriers[(s["kind"], s["type"])].add(s["bin_id"])
        return [{"kind": k, "type": t, "count": n, "share": len(carriers[(k, t)]) / len(members)}
                for (k, t), n in counts.most_common(12)]

    @app.get("/genomes", response_class=HTMLResponse)
    def genomes_page(request: Request):
        return render(request, "genomes.html", "taxa", genomes=catalog.genomes)

    # ------------------------------------------------------------------ #
    # Genome
    # ------------------------------------------------------------------ #

    @app.get("/genome/{bin_id:path}/proteins", response_class=HTMLResponse)
    def genome_proteins(request: Request, bin_id: str):
        genome = catalog.by_bin.get(bin_id)
        if genome is None:
            raise HTTPException(404, "Genome not found")
        with lock:
            rows = store.execute(
                """WITH p AS (SELECT protein_id, contig_id, start, sequence_length FROM proteins WHERE bin_id = ?),
                        best AS (SELECT a.protein_id, ARG_MIN(COALESCE(NULLIF(a.name, ''), a.accession),
                                                              COALESCE(a.evalue, 1)) AS top, COUNT(*) AS hits
                                 FROM annotations a JOIN p USING (protein_id) GROUP BY 1)
                   SELECT p.protein_id, p.contig_id, p.start, p.sequence_length, best.top, COALESCE(best.hits, 0)
                   FROM p LEFT JOIN best USING (protein_id) ORDER BY p.contig_id, p.start""", [bin_id])
        return render(request, "genome_proteins.html", "taxa", g=genome, rows=rows)

    @app.get("/genome/{bin_id:path}", response_class=HTMLResponse)
    def genome_page(request: Request, bin_id: str):
        genome = catalog.by_bin.get(bin_id)
        if genome is None:
            raise HTTPException(404, "Genome not found")
        with lock:
            lengths = [r[0] for r in store.execute("SELECT length FROM contigs WHERE bin_id = ?", [bin_id])]
            largest = store.execute(
                """WITH top AS (SELECT protein_id, sequence_length FROM proteins WHERE bin_id = ?
                                ORDER BY sequence_length DESC LIMIT 12)
                   SELECT top.protein_id, top.sequence_length, COUNT(a.protein_id) FROM top
                   LEFT JOIN annotations a USING (protein_id) GROUP BY ALL ORDER BY 2 DESC""", [bin_id])
        shares = catalog.category_share.get(bin_id, {})
        profile = [{"category": c, "label": CATEGORY_LABELS.get(c, c), "genome": shares.get(c, 0.0),
                    "dataset": catalog.dataset_category_share.get(c, 0.0)} for c in BIOLOGICAL
                   if catalog.dataset_category_share.get(c)]
        modules = []
        column_ids = catalog.module_ids
        if catalog.module_completeness is not None:
            row = catalog.module_completeness[genome.index]
            for j in np.argsort(-row):
                if row[j] < 0.5:
                    break
                d = catalog.modules[column_ids[j]]
                modules.append({"module": d.module, "name": d.name, "completeness": float(row[j]),
                                "class": d.module_class.split(";")[-1].strip()})
        systems = [s for s in catalog.systems if s["bin_id"] == bin_id]
        loci = [l for l in catalog.loci if l["bin_id"] == bin_id]
        abundance = []
        if sidecar.is_file():
            from sharur.abundance import genome_abundance  # noqa: PLC0415

            abundance = genome_abundance(sidecar, bins=[bin_id])
        return render(request, "genome.html", "taxa", g=genome, strip=charts.contig_strip(lengths),
                      longest_contig=max(lengths) if lengths else 0, largest=largest, profile=profile,
                      modules=modules, systems=systems, loci=loci, abundance=abundance)

    # ------------------------------------------------------------------ #
    # Protein
    # ------------------------------------------------------------------ #

    @app.get("/protein/{protein_id:path}/why/{predicate}", response_class=HTMLResponse)
    def why_page(request: Request, protein_id: str, predicate: str):
        from sharur.operators.cards import why  # noqa: PLC0415

        with lock:
            result = why(store, protein_id, predicate)
        bin_id = None
        with lock:
            row = store.execute("SELECT bin_id FROM proteins WHERE protein_id = ?", [protein_id])
            bin_id = row[0][0] if row else None
        return render(request, "why.html", "functions", protein_id=protein_id, predicate=predicate,
                      result=result, g=catalog.by_bin.get(bin_id))

    @app.get("/protein/{protein_id:path}", response_class=HTMLResponse)
    def protein_page(request: Request, protein_id: str):
        from sharur.architecture import architecture  # noqa: PLC0415
        from sharur.contig_context import EdgeContext, describe_edge  # noqa: PLC0415
        from sharur.operators.cards import card  # noqa: PLC0415
        from sharur.operators.navigation import get_neighborhood  # noqa: PLC0415

        with lock:
            c = card(store, protein_id, window=0)
            if not c.get("found"):
                raise HTTPException(404, "Protein not found")
            domains = [d.to_dict() for d in architecture(store, protein_id)]
            hood = get_neighborhood(store, protein_id, window=8).raw or {}
            ids = [g["protein_id"] for g in hood.get("proteins", [])]
            categories = _top_categories(ids)
            nearby = []
            try:
                from sharur.modules import locus_modules  # noqa: PLC0415

                nearby = [m for m in locus_modules(store, protein_id, window=10) if m["steps_complete"] > 0][:8]
            except FileNotFoundError:
                pass
        genes = []
        for i, g in enumerate(hood.get("proteins", []), 1):
            category = categories.get(g["protein_id"])
            genes.append({**g, "number": i, "category": category,
                          "color": charts.CATEGORY_COLORS.get(category) if category else None})
        edge = c.get("contig_edge")
        genome = catalog.by_bin.get(c["genome"]["bin_id"])
        systems_here = [s for s in catalog.systems if protein_id in s["proteins"]]
        return render(request, "protein.html", "taxa", c=c, g=genome, domains=domains,
                      track=charts.domain_track(c["location"]["length_aa"], domains),
                      hood=charts.neighborhood(genes, start_edge=hood.get("contig_start_in_window", False),
                                               end_edge=hood.get("contig_end_in_window", False)),
                      genes=genes, edge=edge, edge_text=describe_edge(EdgeContext(**edge)) if edge else "",
                      nearby=nearby, systems_here=systems_here)

    def _top_categories(ids: list[str]) -> dict[str, str]:
        """Most specific biological category per protein, for neighborhood colors."""
        if not ids:
            return {}
        rows = store.execute(
            "SELECT protein_id, predicates FROM protein_predicates WHERE protein_id IN "
            f"({','.join('?' * len(ids))})", ids)
        out = {}
        for pid, predicates in rows:
            counts = Counter(PREDICATE_BY_ID[p].category for p in predicates or []
                             if p in PREDICATE_BY_ID and PREDICATE_BY_ID[p].category in BIOLOGICAL)
            if counts:
                out[pid] = counts.most_common(1)[0][0]
        return out

    # ------------------------------------------------------------------ #
    # Functions
    # ------------------------------------------------------------------ #

    @app.get("/functions", response_class=HTMLResponse)
    def functions(request: Request):
        groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
        n = len(catalog.genomes) or 1
        for k, pred in enumerate(catalog.predicates):
            definition = PREDICATE_BY_ID.get(pred)
            if definition is None:
                continue
            groups[definition.category].append({
                "id": pred, "name": definition.name, "genomes": int(catalog.predicate_genomes[k]),
                "share": catalog.predicate_genomes[k] / n, "proteins": int(catalog.predicate_proteins[k]),
                "level": definition.level})
        ordered = [(c, sorted(groups[c], key=lambda r: -r["genomes"])) for c in (*BIOLOGICAL, *sorted(
            set(groups) - set(BIOLOGICAL))) if c in groups]
        # labels that separate genomes: carried by some and missing from others
        variable = sorted((r | {"category": c} for c, rows in ordered if c in BIOLOGICAL for r in rows
                           if 0.15 <= r["share"] <= 0.85 and PREDICATE_BY_ID[r["id"]].parent is not None),
                          key=lambda r: abs(r["share"] - 0.5))[:24]
        return render(request, "functions.html", "functions", groups=ordered,
                      variable=sorted(variable, key=lambda r: -r["share"]))

    @app.get("/function/{predicate}", response_class=HTMLResponse)
    def function_page(request: Request, predicate: str, rank: str = Query("class")):
        definition = PREDICATE_BY_ID.get(predicate)
        if definition is None and predicate not in catalog.predicate_index:
            raise HTTPException(404, "Unknown predicate")
        rank = rank if rank in RANKS else "class"
        carriers = catalog.carriers(predicate)
        prevalence = [r for r in catalog.prevalence_by(carriers, rank) if r["genomes"] >= 3][:30]
        top_genomes = sorted(carriers.items(), key=lambda kv: -kv[1])[:15]
        with lock:
            examples = store.execute(
                """SELECT pp.protein_id, p.bin_id, p.sequence_length FROM protein_predicates pp
                   JOIN proteins p USING (protein_id) WHERE list_contains(pp.predicates, ?)
                   ORDER BY p.sequence_length DESC LIMIT 25""", [predicate])
        children = [p for p in PREDICATE_BY_ID.values() if p.parent == predicate]
        parent = PREDICATE_BY_ID.get(definition.parent) if definition and definition.parent else None
        return render(request, "function.html", "functions", predicate=predicate, d=definition,
                      carriers=carriers, prevalence=prevalence, rank=rank,
                      top_genomes=[(catalog.genomes[i], n) for i, n in top_genomes], examples=examples,
                      children=children, parent=parent)

    # ------------------------------------------------------------------ #
    # Pathways
    # ------------------------------------------------------------------ #

    @app.get("/pathways", response_class=HTMLResponse)
    def pathways(request: Request):
        if catalog.module_completeness is None:
            return render(request, "pending.html", "pathways", what="Pathway summaries",
                          missing=catalog.ready.is_set())
        matrix = catalog.module_completeness
        groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for j, module in enumerate(catalog.module_ids):
            column = matrix[:, j]
            if not column.any():
                continue
            d = catalog.modules[module]
            parts = [p.strip() for p in d.module_class.split(";")]
            groups[parts[1] if len(parts) > 1 else parts[0]].append({
                "module": module, "name": d.name, "subclass": parts[-1],
                "complete": float((column >= 0.75).mean()), "partial": float(((column > 0) & (column < 0.75)).mean())})
        ordered = sorted(((k, sorted(v, key=lambda r: -r["complete"])) for k, v in groups.items()),
                         key=lambda kv: -max(r["complete"] for r in kv[1]))
        return render(request, "pathways.html", "pathways", groups=ordered)

    @app.get("/pathway/{module}", response_class=HTMLResponse)
    def pathway_page(request: Request, module: str, rank: str = Query("order")):
        from sharur.modules import evaluate  # noqa: PLC0415

        column = catalog.module_column(module)
        if column is None:
            raise HTTPException(404, "Module not found or pathways still loading")
        rank = rank if rank in RANKS else "order"
        d = catalog.modules[module]
        by_taxon: dict[str, list[float]] = defaultdict(list)
        for g in catalog.genomes:
            by_taxon[g.taxonomy.get(rank) or UNCLASSIFIED].append(float(column[g.index]))
        clades = sorted(({"taxon": t, "genomes": len(v), "mean": float(np.mean(v)),
                          "complete": float(np.mean([x >= 0.75 for x in v]))} for t, v in by_taxon.items()
                         if len(v) >= 3), key=lambda r: -r["genomes"])[:30]
        # step-level: how often each step is satisfied among genomes with the module partly present
        partial = [g for g in catalog.genomes if column[g.index] > 0]
        step_hits: Counter = Counter()
        step_missing: dict[int, Counter] = defaultdict(Counter)
        n_steps = 0
        for g in partial:
            result = evaluate(d, catalog.ko_sets.get(g.bin_id, set()), catalog.modules)
            n_steps = max(n_steps, result.steps_total)
            for s in result.steps:
                if s.score >= 1.0:
                    step_hits[s.index] += 1
                else:
                    for ko in s.missing[:3]:
                        step_missing[s.index][ko] += 1
        steps = [{"index": i, "share": step_hits[i] / len(partial) if partial else 0,
                  "missing": step_missing[i].most_common(3)}
                 for i in sorted({s.index for s in evaluate(d, set(), catalog.modules).steps})]
        genomes = sorted(((catalog.genomes[i], float(column[i])) for i in np.nonzero(column)[0]),
                         key=lambda gc: -gc[1])[:40]
        histogram = np.histogram(column[column > 0], bins=10, range=(0, 1))[0].tolist() if partial else []
        return render(request, "pathway.html", "pathways", d=d, clades=clades, steps=steps, rank=rank,
                      genomes=genomes, partial=len(partial), complete=int((column >= 0.75).sum()),
                      histogram=histogram)

    # ------------------------------------------------------------------ #
    # Systems
    # ------------------------------------------------------------------ #

    @app.get("/systems", response_class=HTMLResponse)
    def systems(request: Request):
        kinds: dict[str, list[dict[str, Any]]] = defaultdict(list)
        grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
        for s in catalog.systems:
            grouped[(s["kind"], s["type"])].append(s)
        n = len(catalog.genomes) or 1
        for (kind, system_type), members in grouped.items():
            kinds[kind].append({"type": system_type, "count": len(members),
                                "share": len({m["bin_id"] for m in members}) / n})
        loci = defaultdict(set)
        locus_counts = Counter()
        for l in catalog.loci:
            loci[l["type"]].add(l["bin_id"])
            locus_counts[l["type"]] += 1
        return render(request, "systems.html", "systems",
                      kinds={k: sorted(v, key=lambda r: -r["count"]) for k, v in kinds.items()},
                      loci=[{"type": t, "count": locus_counts[t], "share": len(b) / n} for t, b in loci.items()])

    @app.get("/system/{kind}/{system_type:path}", response_class=HTMLResponse)
    def system_page(request: Request, kind: str, system_type: str, rank: str = Query("class")):
        members = [s for s in catalog.systems if s["kind"] == kind and s["type"] == system_type]
        if not members:
            raise HTTPException(404, "System type not found")
        rank = rank if rank in RANKS else "class"
        carriers = {catalog.by_bin[s["bin_id"]].index: 1 for s in members if s["bin_id"] in catalog.by_bin}
        prevalence = [r for r in catalog.prevalence_by(carriers, rank) if r["genomes"] >= 3][:30]
        subtypes = Counter(s["subtype"] or "–" for s in members)
        return render(request, "system.html", "systems", kind=kind, system_type=system_type, members=members[:300],
                      total=len(members), prevalence=prevalence, rank=rank, subtypes=subtypes.most_common(),
                      genomes=len(carriers))

    @app.get("/loci/{locus_type}", response_class=HTMLResponse)
    def loci_page(request: Request, locus_type: str):
        members = [l for l in catalog.loci if l["type"] == locus_type]
        if not members:
            raise HTTPException(404, "Locus type not found")
        per_genome = Counter(l["bin_id"] for l in members)
        return render(request, "loci.html", "systems", locus_type=locus_type, members=members[:500],
                      total=len(members), per_genome=per_genome.most_common(20))

    # ------------------------------------------------------------------ #
    # Discover, search
    # ------------------------------------------------------------------ #

    @app.get("/discover", response_class=HTMLResponse)
    def discover(request: Request):
        if not catalog.notable:
            return render(request, "pending.html", "discover", what="Notable proteins",
                          missing=catalog.ready.is_set())
        return render(request, "discover.html", "discover", notable=catalog.notable)

    @app.get("/architecture", response_class=HTMLResponse)
    def architecture_page(request: Request, pattern: str = Query("", max_length=500),
                          limit: int = Query(200, ge=1, le=2000)):
        from sharur.architecture import PatternError, search_architecture  # noqa: PLC0415

        result, error = None, None
        if pattern.strip():
            try:
                with lock:
                    result = search_architecture(store, pattern, limit=limit)
            except PatternError as exc:
                error = str(exc)
        return render(request, "architecture.html", "discover", pattern=pattern, result=result, error=error)

    @app.get("/search")
    def search(request: Request, q: str = Query("", max_length=500)):
        q = q.strip()
        if not q:
            return RedirectResponse("/", status_code=303)
        exact = _exact(q)
        if exact:
            return RedirectResponse(exact, status_code=303)
        results = _suggest(q, limit=60)
        return render(request, "search.html", "home", q=q, results=results)

    def _exact(q: str) -> str | None:
        if q in catalog.by_bin:
            return _url("genome", q)
        if q in PREDICATE_BY_ID:
            return _url("function", q)
        if q in catalog.modules:
            return _url("pathway", q)
        with lock:
            if store.execute("SELECT 1 FROM proteins WHERE protein_id = ?", [q]):
                return _url("protein", q)
        return None

    def _suggest(q: str, limit: int = 12) -> list[dict[str, str]]:
        needle = q.lower()
        out: list[dict[str, str]] = []

        def add(kind, label, sub, url):
            out.append({"kind": kind, "label": label, "sub": sub, "url": url})

        seen = set()
        for g in catalog.genomes:
            for rank in RANKS:
                name = g.taxonomy.get(rank)
                if name and name != UNCLASSIFIED and needle in name.lower() and (rank, name) not in seen:
                    seen.add((rank, name))
                    add("taxon", name, rank, _url("taxa", rank, name))
        for pred, d in PREDICATE_BY_ID.items():
            if pred in catalog.predicate_index and (needle in pred or needle in d.name.lower()):
                add("function", d.name, pred, _url("function", pred))
        for module, d in catalog.modules.items():
            if needle in module.lower() or needle in d.name.lower():
                add("pathway", d.name, module, _url("pathway", module))
        for kind, system_type in sorted({(s["kind"], s["type"]) for s in catalog.systems}):
            if needle in str(system_type).lower():
                add("system", system_type, kind, _url("system", kind, system_type))
        for g in catalog.genomes:
            if needle in g.bin_id.lower():
                add("genome", g.bin_id, g.label, _url("genome", g.bin_id))
                if sum(1 for o in out if o["kind"] == "genome") >= 8:
                    break
        with lock:
            for (pid,) in store.execute(
                    "SELECT protein_id FROM proteins WHERE protein_id ILIKE ? LIMIT 6", [f"%{q}%"]):
                add("protein", pid, "", _url("protein", pid))
        exact_first = sorted(out, key=lambda o: (o["label"].lower() != needle,
                                                 not o["label"].lower().startswith(needle)))
        return exact_first[:limit]

    @app.get("/api/suggest")
    def suggest(q: str = Query("", max_length=200)):
        return JSONResponse(_suggest(q.strip()) if len(q.strip()) >= 2 else [])

    @app.get("/api/status")
    def status():
        return {"status": catalog.status, "ready": catalog.ready.is_set()}

    return app


__all__ = ["Catalog", "create_app"]
