"""Routes for the dataset browser."""

from __future__ import annotations

import secrets
import threading
from collections import Counter, defaultdict
from pathlib import Path
from types import SimpleNamespace
from typing import Any
from urllib.parse import quote

import numpy as np
from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.responses import HTMLResponse, JSONResponse, PlainTextResponse, RedirectResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

from sharur.browser import (
    charts,
    clade_content,
    routes_compare,
    routes_cctyper,
    routes_crispr,
    routes_curation,
    routes_loci,
    routes_matrix,
    routes_search,
    routes_systems,
)
from sharur.browser.routes_insight import register as register_insight
from sharur.browser.catalog import (
    BIOLOGICAL,
    CATEGORY_LABELS,
    RANKS,
    UNCLASSIFIED,
    Catalog,
    load_background,
    load_catalog,
)
from sharur.predicates.mappings.pfam_map import PFAM_EVIDENCE, PFAM_TO_PREDICATES
from sharur.predicates.mappings.vog_map import VOG_CATEGORY_NAMES
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


def _gshort(bin_id: str, limit: int = 32) -> str:
    """Genome names share long prefixes; keep both ends so names stay distinguishable."""
    return bin_id if len(bin_id) <= limit else bin_id[:11] + "…" + bin_id[-(limit - 12):]


_RESIDUE_MASS = {  # average residue masses (Da), water added once per chain
    "A": 71.0788, "R": 156.1875, "N": 114.1038, "D": 115.0886, "C": 103.1388, "E": 129.1155, "Q": 128.1307,
    "G": 57.0519, "H": 137.1411, "I": 113.1594, "L": 113.1594, "K": 128.1741, "M": 131.1926, "F": 147.1766,
    "P": 97.1167, "S": 87.0782, "T": 101.1051, "W": 186.2132, "Y": 163.1760, "V": 99.1326,
}


def _sequence_view(sequence: str | None) -> dict[str, Any] | None:
    """Sequence split into numbered 60-residue lines of 10-residue blocks, with simple stats."""
    if not sequence:
        return None
    seq = sequence.strip().rstrip("*").upper()
    lines = [(i + 1, [seq[j:j + 10] for j in range(i, min(i + 60, len(seq)), 10)]) for i in range(0, len(seq), 60)]
    mass = sum(_RESIDUE_MASS.get(a, 110.0) for a in seq) + 18.015
    return {"raw": seq, "lines": lines, "length": len(seq), "kda": mass / 1000}


def _evalue(x: Any) -> str:
    return f"{x:.1e}" if isinstance(x, float) else ("–" if x is None else str(x))


def create_app(db_path: str | Path, *, token: str | None = None, background: bool = True,
               notes_path: str | Path | None = None, assemblies: list[str | Path] | None = None) -> FastAPI:
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
    templates.env.globals.update(url=_url, bp=_bp, num=_num, pct=_pct, evalue=_evalue, quote=quote, short=_short, gshort=_gshort,
                                 PREDICATE_BY_ID=PREDICATE_BY_ID, VOG_CATEGORY_NAMES=VOG_CATEGORY_NAMES,
                                 charts=charts, CATEGORY_LABELS=CATEGORY_LABELS, dataset=dataset_name,
                                 catalog=catalog, CATEGORY_COLORS=charts.CATEGORY_COLORS,
                                 # static assets change with the package; bust browser caches on upgrade
                                 asset_version=int(max(f.stat().st_mtime for f in (HERE / "static").iterdir())))

    app = FastAPI(title="Sharur browser", docs_url=None, redoc_url=None, openapi_url=None)
    app.mount("/static", StaticFiles(directory=str(HERE / "static")), name="static")
    app.state.store, app.state.catalog = store, catalog

    from sharur.browser.annotations import SOURCE_LABELS, KoNames, describe_ko, domain_lanes, external_url, hit_view  # noqa: PLC0415

    ko_names = KoNames(store)
    with lock:
        ko_names.get("K00001")  # load once, before handlers share the connection

    def describe_hit(label: str | None) -> str | None:
        """Name bare VOG and KO ids in a hit label (VOGdb consensus, KEGG symbol and definition)."""
        if not label:
            return label
        head = label.split(" ", 1)[0].strip("()")
        vog = catalog.vogs.get(head)
        if vog and vog["description"]:
            return f'{vog["description"]} ({head})'
        return describe_ko(label, ko_names)

    def hit_rows(annotations: dict[str, list[dict[str, Any]]]) -> list[dict[str, Any]]:
        rows = []
        for src, hits in annotations.items():
            for h in hits:
                rows.append({**h, **hit_view(src, h, ko_names), "source": src})
        return rows

    templates.env.globals.update(describe_hit=describe_hit, hit_rows=hit_rows, SOURCE_LABELS=SOURCE_LABELS)

    def render(request: Request, template: str, section: str, /, **context: Any) -> HTMLResponse:
        return templates.TemplateResponse(request, template, {"section": section, **context})

    # Feature modules register before the catch-all path routes below.
    ctx = SimpleNamespace(store=store, lock=lock, catalog=catalog, render=render, url=_url,
                          describe_hit=describe_hit, predicates=PREDICATE_BY_ID, db_path=Path(db_path),
                          notes_path=Path(notes_path) if notes_path else None, ko_names=ko_names)
    app.state.ctx = ctx
    routes_loci.register(app, ctx)
    routes_curation.register(app, ctx)
    routes_systems.register(app, ctx)
    routes_crispr.register(app, ctx, [Path(p) for p in assemblies or []])
    routes_cctyper.register(app, ctx)

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
    def taxa(request: Request, rank: str | None = None, name: str | None = None,
             pathways: str = Query("median")):
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
                      systems=_clade_systems(genomes),
                      pathway_mode="present" if pathways == "present" else "median",
                      clade_modules=catalog.clade_modules(genomes, mode="present" if pathways == "present"
                                                          else "median"),
                      distinctive_modules=catalog.distinctive_modules(genomes))

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

    @app.get("/genome/{bin_id:path}/contigs", response_class=HTMLResponse)
    def genome_contigs(request: Request, bin_id: str):
        genome = catalog.by_bin.get(bin_id)
        if genome is None:
            raise HTTPException(404, "Genome not found")
        with lock:
            rows = store.execute(
                """WITH p AS (SELECT protein_id, contig_id FROM proteins WHERE bin_id = ?),
                        ann AS (SELECT DISTINCT a.protein_id FROM annotations a JOIN p USING (protein_id))
                   SELECT c.contig_id, c.length, COUNT(p.protein_id), COUNT(ann.protein_id)
                   FROM contigs c LEFT JOIN p USING (contig_id) LEFT JOIN ann USING (protein_id)
                   WHERE c.bin_id = ? GROUP BY 1, 2 ORDER BY 2 DESC""", [bin_id, bin_id])
        features: dict[str, list[str]] = defaultdict(list)
        for l in catalog.loci:
            if l["bin_id"] == bin_id:
                features[l["contig_id"]].append(l["type"])
        contig_of = _protein_contigs([p for s in catalog.systems if s["bin_id"] == bin_id for p in s["proteins"]])
        for s in catalog.systems:
            if s["bin_id"] == bin_id:
                for contig in {contig_of.get(p) for p in s["proteins"]} - {None}:
                    features[contig].append(s["type"])
        return render(request, "genome_contigs.html", "taxa", g=genome, rows=rows, features=features)

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
        rows = [(pid, contig, start, length, describe_hit(top), hits) for pid, contig, start, length, top, hits in rows]
        return render(request, "genome_proteins.html", "taxa", g=genome, rows=rows)

    @app.get("/genome/{bin_id:path}", response_class=HTMLResponse)
    def genome_page(request: Request, bin_id: str):
        genome = catalog.by_bin.get(bin_id)
        if genome is None:
            raise HTTPException(404, "Genome not found")
        with lock:
            contigs = store.execute("SELECT contig_id, length FROM contigs WHERE bin_id = ?", [bin_id])
            lengths = [n for _, n in contigs]
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
        cas_calls = [c for c in ctx.cctyper_systems() if c["bin_id"] == bin_id]
        abundance = []
        if sidecar.is_file():
            from sharur.abundance import genome_abundance  # noqa: PLC0415

            abundance = genome_abundance(sidecar, bins=[bin_id])
        return render(request, "genome.html", "taxa", g=genome,
                      strip=charts.contig_strip(contigs, lambda c: _url("contig", c)),
                      longest_contig=max(lengths) if lengths else 0, largest=largest, profile=profile,
                      modules=modules, systems=systems, loci=loci, abundance=abundance,
                      cas_calls=cas_calls)

    # ------------------------------------------------------------------ #
    # Protein
    # ------------------------------------------------------------------ #

    @app.get("/fasta/{protein_id:path}")
    def protein_fasta(protein_id: str):
        with lock:
            rows = store.execute("SELECT bin_id, sequence FROM proteins WHERE protein_id = ?", [protein_id])
        if not rows or not rows[0][1]:
            raise HTTPException(404, "Protein or sequence not found")
        bin_id, seq = rows[0]
        seq = seq.strip().rstrip("*")
        body = f">{protein_id} genome={bin_id}\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n"
        filename = "".join(c if c.isalnum() or c in "._-" else "_" for c in protein_id)[:120] + ".faa"
        return PlainTextResponse(body, headers={"Content-Disposition": f'inline; filename="{filename}"'})

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
            lanes = domain_lanes(store, protein_id)
            sequence_rows = store.execute("SELECT sequence FROM proteins WHERE protein_id = ?", [protein_id])
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
            genes.append({**g, "annotation": describe_hit(g.get("annotation")), "number": i, "category": category,
                          "color": charts.CATEGORY_COLORS.get(category) if category else None})
        edge = c.get("contig_edge")
        genome = catalog.by_bin.get(c["genome"]["bin_id"])
        order = [g["protein_id"] for g in genes]
        here = order.index(protein_id) if protein_id in order else -1
        prev_gene = order[here - 1] if here > 0 else None
        next_gene = order[here + 1] if 0 <= here < len(order) - 1 else None
        systems_here = [s for s in catalog.systems if protein_id in s["proteins"]]
        length = c["location"]["length_aa"]
        tracks = []
        for source, placed in lanes:
            for d in placed:
                ko = d["accession"].split(".")[0]
                known = ko_names.get(ko) if ko[:1] == "K" else None
                if known:
                    d["label"] = known[0].split(",")[0].strip() or ko
                    d["title"] = f"{ko} {known[1]}"
                if source == "pfam":
                    continue
                d["href"] = external_url(source, d["accession"])
            tracks.append((source, len(placed), charts.domain_track(length, placed)))
        return render(request, "protein.html", "taxa", c=c, g=genome, domains=domains,
                      track=charts.domain_track(length, domains), tracks=tracks,
                      hood=charts.neighborhood(genes, start_edge=hood.get("contig_start_in_window", False),
                                               end_edge=hood.get("contig_end_in_window", False)),
                      genes=genes, edge=edge, edge_text=describe_edge(EdgeContext(**edge)) if edge else "",
                      nearby=nearby, systems_here=systems_here, prev_gene=prev_gene, next_gene=next_gene,
                      sequence=_sequence_view(sequence_rows[0][0] if sequence_rows else None))

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

    ctx.top_categories = _top_categories

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
    def function_page(request: Request, predicate: str, rank: str = Query("class"),
                      clade: str = Query("", max_length=500), genome: str = Query("", max_length=500)):
        from sharur.browser.routes_compare import clade_filter  # noqa: PLC0415

        definition = PREDICATE_BY_ID.get(predicate)
        if definition is None and predicate not in catalog.predicate_index:
            raise HTTPException(404, "Unknown predicate")
        rank = rank if rank in RANKS else "class"
        restrict = clade_filter(catalog, genome or clade)
        carriers = catalog.carriers(predicate)
        if restrict:
            carriers = {i: n for i, n in carriers.items() if catalog.genomes[i].bin_id in restrict[1]}
        prevalence = [r for r in catalog.prevalence_by(carriers, rank) if r["genomes"] >= 3][:30]
        top_genomes = sorted(carriers.items(), key=lambda kv: -kv[1])[:15]
        with lock:
            if restrict:
                examples = store.execute(
                    """SELECT pp.protein_id, p.bin_id, p.sequence_length FROM protein_predicates pp
                       JOIN proteins p USING (protein_id) WHERE list_contains(pp.predicates, ?)
                         AND p.bin_id IN (SELECT UNNEST(?::VARCHAR[]))
                       ORDER BY p.sequence_length DESC LIMIT 50""", [predicate, sorted(restrict[1])])
            else:
                examples = store.execute(
                    """SELECT pp.protein_id, p.bin_id, p.sequence_length FROM protein_predicates pp
                       JOIN proteins p USING (protein_id) WHERE list_contains(pp.predicates, ?)
                       ORDER BY p.sequence_length DESC LIMIT 25""", [predicate])
        children = [p for p in PREDICATE_BY_ID.values() if p.parent == predicate]
        parent = PREDICATE_BY_ID.get(definition.parent) if definition and definition.parent else None
        return render(request, "function.html", "functions", predicate=predicate, d=definition,
                      carriers=carriers, prevalence=prevalence, rank=rank,
                      top_genomes=[(catalog.genomes[i], n) for i, n in top_genomes], examples=examples,
                      children=children, parent=parent, restrict=restrict, restrict_key=genome or clade)

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
            if "crispr" in (l["type"] or "").lower():
                continue  # CRISPR has its own panel
            loci[l["type"]].add(l["bin_id"])
            locus_counts[l["type"]] += 1
        crispr = []
        arrays = ctx.crispr_arrays()
        if arrays:
            crispr.append({"label": "CRISPR arrays", "url": "/crispr", "count": len(arrays),
                           "share": len({a["bin_id"] for a in arrays}) / n})
        cas = ctx.cas_loci()
        for kind_name in ("array + Cas", "Cas genes only", "array only"):
            members = [l for l in cas if l["kind"] == kind_name]
            if members:
                crispr.append({"label": f"Cas-domain loci: {kind_name}",
                               "url": "/crispr-cas?kind=" + quote(kind_name), "count": len(members),
                               "share": len({l["bin_id"] for l in members}) / n})
        return render(request, "systems.html", "systems",
                      kinds={k: sorted(v, key=lambda r: -r["count"]) for k, v in kinds.items()},
                      loci=[{"type": t, "count": locus_counts[t], "share": len(b) / n} for t, b in loci.items()],
                      crispr=crispr, cctyper=ctx.cctyper_summary())

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

    from sharur.browser import discover as discover_feeds  # noqa: PLC0415

    def rare_systems() -> list[dict[str, Any]]:
        if not hasattr(ctx, "_rare_systems"):
            ctx._rare_systems = discover_feeds.rare_systems(catalog, getattr(ctx, "cctyper_systems", None))
        return ctx._rare_systems

    def feeds_or_pending(request: Request):
        if not catalog.notable or "fusions" not in catalog.notable:
            return render(request, "pending.html", "discover", what="Discovery lists",
                          missing=catalog.ready.is_set())
        return None

    @app.get("/discover", response_class=HTMLResponse)
    def discover(request: Request):
        return feeds_or_pending(request) or render(request, "discover.html", "discover", f=catalog.notable,
                                                   rare=rare_systems(), feeds=discover_feeds.FEEDS)

    @app.get("/discover/random")
    def discover_random():
        import random  # noqa: PLC0415

        f = catalog.notable or {}
        picks = [_url("protein", p["protein_id"]) for p in f.get("giants", [])[:100] + f.get("dark", [])[:100]]
        picks += [_url("protein", x["example"]) for x in f.get("fusions", [])]
        picks += [f'{_url("contig", i["contig_id"])}?start={max(1, i["start"] - 3000)}&span={i["end"] - i["start"] + 6000}'
                  for i in f.get("islands", [])[:100]]
        return RedirectResponse(random.choice(picks) if picks else "/discover", status_code=303)

    @app.get("/discover/{feed}", response_class=HTMLResponse)
    def discover_feed(request: Request, feed: str):
        if feed not in discover_feeds.FEEDS:
            raise HTTPException(404, "Unknown list")
        pending = feeds_or_pending(request)
        if pending:
            return pending
        rows = rare_systems() if feed == "systems" else catalog.notable.get(feed, [])
        title, why = discover_feeds.FEEDS[feed]
        return render(request, "discover_feed.html", "discover", feed=feed, rows=rows, f=catalog.notable,
                      title=title, why=why)

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
        scoped = routes_search.SCOPED.match(q)
        if scoped:
            scope = routes_search.resolve_scope(ctx, scoped.group(2))
            if scope is not None:
                kind, label, bins = scope
                found = routes_search.scoped_search(ctx, scoped.group(1), bins)
                scope_url = _url("genome", label) if kind == "genome" else (
                    _url("taxa", kind, label) if kind in RANKS else None)
                return render(request, "search_scoped.html", "home", q=q, term=scoped.group(1), scope_kind=kind,
                              scope_label=label, scope_url=scope_url, scope_genomes=len(bins),
                              other_ranks=[r for r in routes_search.same_name_ranks(ctx, label) if r[0] != kind]
                              if kind in RANKS else [], **found)
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
        if q.upper() in catalog.vogs:
            return _url("vog", q.upper())
        if q.split(".")[0].upper() in catalog.domains:
            return _url("domain", q.split(".")[0].upper())
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
        for d in catalog.domains.values():
            if needle == d["accession"].lower() or needle in d["name"].lower():
                add("domain", d["name"], f'{d["accession"]} · {d["genomes"]:,} genomes', _url("domain", d["accession"]))
                if sum(1 for o in out if o["kind"] == "domain") >= 8:
                    break
        for v in catalog.vogs.values():
            if needle == v["accession"].lower() or (v["description"] and needle in v["description"].lower()):
                add("VOG", v["name"], f'{v["accession"]} · {v["genomes"]:,} genomes', _url("vog", v["accession"]))
                if sum(1 for o in out if o["kind"] == "VOG") >= 6:
                    break
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

    # ------------------------------------------------------------------ #
    # Pfam domains
    # ------------------------------------------------------------------ #

    @app.get("/domains", response_class=HTMLResponse)
    def domains(request: Request, q: str = Query("", max_length=200), sort: str = Query("genomes"),
                page: int = Query(1, ge=1)):
        if not catalog.domains:
            return render(request, "pending.html", "domains", what="Pfam domains", missing=catalog.ready.is_set())
        n = len(catalog.genomes) or 1
        rows = list(catalog.domains.values())
        needle = q.strip().lower()
        if needle:
            rows = [r for r in rows if needle in r["name"].lower() or needle in r["accession"].lower()
                    or needle in r["description"].lower()]
        key = {"genomes": lambda r: -r["genomes"], "proteins": lambda r: -r["proteins"],
               "name": lambda r: r["name"].lower(), "rare": lambda r: (r["genomes"], r["name"].lower()),
               "patchy": lambda r: abs(r["genomes"] / n - 0.5)}.get(sort, lambda r: -r["genomes"])
        rows.sort(key=key)
        per_page = 100
        pages = max(1, (len(rows) + per_page - 1) // per_page)
        page = min(page, pages)
        shown = rows[(page - 1) * per_page: page * per_page]
        mapped = {r["accession"]: PFAM_TO_PREDICATES.get(r["accession"], []) for r in shown}
        return render(request, "domains.html", "domains", rows=shown, total=len(rows), q=q, sort=sort, page=page,
                      pages=pages, n=n, mapped=mapped, families=len(catalog.domains))

    @app.get("/domain/{accession}", response_class=HTMLResponse)
    def domain_page(request: Request, accession: str, rank: str = Query("class")):
        from sharur.architecture import compact, resolve  # noqa: PLC0415
        from sharur.architecture import _hits as architecture_hits  # noqa: PLC0415

        acc = accession.split(".")[0]
        summary = catalog.domains.get(acc)
        if summary is None:
            raise HTTPException(404, "Pfam family not found in this dataset")
        rank = rank if rank in RANKS else "class"
        with lock:
            per_genome = store.execute(
                """SELECT p.bin_id, COUNT(DISTINCT a.protein_id) FROM annotations a JOIN proteins p USING (protein_id)
                   WHERE LOWER(a.source) = 'pfam' AND split_part(a.accession, '.', 1) = ? GROUP BY 1""", [acc])
            carriers_ids = [r[0] for r in store.execute(
                """SELECT DISTINCT protein_id FROM annotations
                   WHERE LOWER(source) = 'pfam' AND split_part(accession, '.', 1) = ?
                   ORDER BY hash(protein_id) LIMIT 2500""", [acc])]
            lengths = [r[0] for r in store.execute(
                """SELECT p.sequence_length FROM proteins p JOIN (SELECT DISTINCT protein_id FROM annotations
                   WHERE LOWER(source) = 'pfam' AND split_part(accession, '.', 1) = ?) c USING (protein_id)""", [acc])]
            hits = architecture_hits(store, ["pfam"], "AND a.protein_id IN (SELECT UNNEST(?::VARCHAR[]))",
                                     [carriers_ids]) if carriers_ids else {}
            examples = store.execute(
                """SELECT p.protein_id, p.bin_id, p.sequence_length FROM proteins p
                   JOIN (SELECT DISTINCT protein_id FROM annotations
                         WHERE LOWER(source) = 'pfam' AND split_part(accession, '.', 1) = ?) c USING (protein_id)
                   ORDER BY p.sequence_length DESC LIMIT 12""", [acc])
        carriers = {catalog.by_bin[b].index: n for b, n in per_genome if b in catalog.by_bin}
        architectures: Counter = Counter()
        patterns: dict[str, str] = {}
        partners: Counter = Counter()
        for pid, protein_hits in hits.items():
            resolved = resolve(protein_hits)
            names = [d.name for d in resolved]
            key = compact(names)
            architectures[key] += 1
            if key not in patterns:
                runs: list[list[Any]] = []
                for n in names:
                    if runs and runs[-1][0] == n:
                        runs[-1][1] += 1
                    else:
                        runs.append([n, 1])
                patterns[key] = "^ " + " ".join(n if k == 1 else f"{n} {{{k}}}" for n, k in runs) + " $"
            for other in {d.accession.split(".")[0]: d.name for d in resolved if d.accession.split(".")[0] != acc}.items():
                partners[other] += 1
        sampled = len(hits)
        labels = sorted(PFAM_EVIDENCE.get(acc, {}).items())
        return render(request, "domain.html", "domains", d=summary, carriers=carriers,
                      prevalence=[r for r in catalog.prevalence_by(carriers, rank) if r["genomes"] >= 3][:30],
                      rank=rank, architectures=[(a, n, patterns[a]) for a, n in architectures.most_common(12)],
                      partners=partners.most_common(16),
                      sampled=sampled, labels=labels, examples=examples,
                      histogram=charts.histogram(lengths, "Protein length (aa)"),
                      top_genomes=sorted(((catalog.genomes[i], n) for i, n in carriers.items()), key=lambda x: -x[1])[:10])

    @app.get("/vogs", response_class=HTMLResponse)
    def vogs(request: Request, q: str = Query("", max_length=200), sort: str = Query("genomes"),
             category: str = Query(""), page: int = Query(1, ge=1)):
        if not catalog.vogs:
            return render(request, "pending.html", "domains", what="VOG families", missing=catalog.ready.is_set())
        n = len(catalog.genomes) or 1
        rows = list(catalog.vogs.values())
        needle = q.strip().lower()
        if needle:
            rows = [r for r in rows if needle in r["accession"].lower() or needle in r["description"].lower()]
        if category:
            rows = [r for r in rows if category in (r["category"] or "")]
        key = {"genomes": lambda r: -r["genomes"], "proteins": lambda r: -r["proteins"],
               "name": lambda r: r["name"].lower(), "rare": lambda r: (r["genomes"], r["name"].lower()),
               "patchy": lambda r: abs(r["genomes"] / n - 0.5)}.get(sort, lambda r: -r["genomes"])
        rows.sort(key=key)
        pages = max(1, (len(rows) + 99) // 100)
        page = min(page, pages)
        counts = Counter(code for r in catalog.vogs.values() for code in VOG_CATEGORY_NAMES if code in (r["category"] or ""))
        return render(request, "vogs.html", "domains", rows=rows[(page - 1) * 100: page * 100], total=len(rows),
                      q=q, sort=sort, category=category, page=page, pages=pages, n=n, counts=counts,
                      families=len(catalog.vogs), has_reference=any(v["annotated"] for v in catalog.vogs.values()))

    @app.get("/vog/{vog_id}", response_class=HTMLResponse)
    def vog_page(request: Request, vog_id: str, rank: str = Query("class")):
        from sharur.predicates.mappings.vog_map import vog_evidence  # noqa: PLC0415

        vog = catalog.vogs.get(vog_id)
        if vog is None:
            raise HTTPException(404, "VOG not found in this dataset")
        rank = rank if rank in RANKS else "class"
        with lock:
            carriers_rows = store.execute(
                """SELECT p.protein_id, p.bin_id, p.contig_id, p.start, p.end_coord, p.sequence_length
                   FROM annotations a JOIN proteins p USING (protein_id)
                   WHERE LOWER(a.source) IN ('vogdb', 'vog') AND a.accession = ?""", [vog_id])
            ids = list({r[0] for r in carriers_rows})
            partners = store.execute(
                """SELECT COALESCE(NULLIF(name, ''), accession), split_part(accession, '.', 1), COUNT(DISTINCT protein_id)
                   FROM annotations WHERE LOWER(source) = 'pfam' AND protein_id IN (SELECT UNNEST(?::VARCHAR[]))
                   GROUP BY 1, 2 ORDER BY 3 DESC LIMIT 15""", [ids]) if ids else []
        per_genome = Counter(r[1] for r in {(r[0], r[1]) for r in carriers_rows})
        carriers = {catalog.by_bin[b].index: n for b, n in per_genome.items() if b in catalog.by_bin}
        # where carriers sit: prophage regions, islands, or elsewhere
        regions: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
        for l in catalog.loci:
            if l["start"] is not None and l["end"] is not None:
                regions[l["contig_id"]].append((l["type"], l["start"], l["end"]))
        context: Counter = Counter()
        seen = set()
        for pid, _bin, contig, start, end, _n in carriers_rows:
            if pid in seen:
                continue
            seen.add(pid)
            kinds = {t for t, a, b in regions.get(contig, []) if start <= b and end >= a}
            context["prophage" if "prophage" in kinds else ("island" if kinds else "elsewhere")] += 1
        lengths = [n for pid, *_, n in {r[0]: r for r in carriers_rows}.values() if n]
        examples = sorted({r[0]: r for r in carriers_rows}.values(), key=lambda r: -(r[5] or 0))[:12]
        labels = sorted(vog_evidence(vog_id, vog["category"], vog["description"]).items())
        return render(request, "vog.html", "domains", v=vog, carriers=carriers, context=context, total=len(seen),
                      prevalence=[r for r in catalog.prevalence_by(carriers, rank) if r["genomes"] >= 3][:30],
                      rank=rank, partners=partners, labels=labels, examples=examples,
                      histogram=charts.histogram(lengths, "Protein length (aa)"),
                      top_genomes=sorted(((catalog.genomes[i], n) for i, n in carriers.items()), key=lambda x: -x[1])[:10])

    # ------------------------------------------------------------------ #
    # Genome browser
    # ------------------------------------------------------------------ #

    def _protein_contigs(protein_ids: list[str]) -> dict[str, str]:
        if not protein_ids:
            return {}
        with lock:
            return dict(store.execute("SELECT protein_id, contig_id FROM proteins WHERE protein_id IN "
                                      "(SELECT UNNEST(?::VARCHAR[]))", [list(set(protein_ids))]))

    @app.get("/contig/{contig_id:path}", response_class=HTMLResponse)
    def contig_page(request: Request, contig_id: str, start: int = Query(1, ge=1), span: int = Query(0, ge=0)):
        with lock:
            info = store.execute("SELECT bin_id, length FROM contigs WHERE contig_id = ?", [contig_id])
            if not info:
                raise HTTPException(404, "Contig not found")
            bin_id, length = info[0]
            genes = store.execute(
                """WITH p AS (SELECT protein_id, start, end_coord, strand, sequence_length FROM proteins
                              WHERE contig_id = ?),
                        best AS (SELECT a.protein_id, ARG_MIN(COALESCE(NULLIF(a.name, ''), a.accession),
                                                              COALESCE(a.evalue, 1)) AS top
                                 FROM annotations a JOIN p USING (protein_id) GROUP BY 1)
                   SELECT p.protein_id, p.start, p.end_coord, p.strand, p.sequence_length, best.top
                   FROM p LEFT JOIN best USING (protein_id) ORDER BY p.start""", [contig_id])
            categories = _top_categories([g[0] for g in genes])
        length = max(length or 0, max((g[2] for g in genes), default=0))
        span = span or min(length, 60000)
        span = max(2000, min(span, length))
        start = max(1, min(start, max(1, length - span + 1)))
        end = start + span - 1
        window = [{"protein_id": pid, "start": s, "end": e, "strand": st, "length_aa": n,
                   "annotation": describe_hit(top),
                   "category": categories.get(pid), "color": charts.CATEGORY_COLORS.get(categories.get(pid))}
                  for pid, s, e, st, n, top in genes if e >= start and s <= end]
        members = {g[0] for g in genes}
        overlays = [{"label": l["type"], "start": l["start"], "end": l["end"], "kind": "locus"}
                    for l in catalog.loci if l["contig_id"] == contig_id]
        positions = {g[0]: (g[1], g[2]) for g in genes}
        for s in catalog.systems:
            inside = [positions[p] for p in s["proteins"] if p in members]
            if inside:
                overlays.append({"label": s["type"], "start": min(a for a, _ in inside),
                                 "end": max(b for _, b in inside), "kind": s["kind"]})
        genome = catalog.by_bin.get(bin_id)
        return render(request, "contig.html", "taxa", contig_id=contig_id, g=genome, length=length, start=start,
                      end=end, span=span, genes=window, total_genes=len(genes),
                      overview=charts.contig_overview(length, start, end, [(g[1], g[2]) for g in genes],
                                                      _url("contig", contig_id), span),
                      track=charts.contig_track(window, overlays, start, end))

    @app.get("/api/suggest")
    def suggest(q: str = Query("", max_length=200), kind: str = Query("", max_length=100)):
        if len(q.strip()) < 2:
            return JSONResponse([])
        if not kind:
            return JSONResponse(_suggest(q.strip()))
        kinds = {k.strip() for k in kind.split(",")}
        return JSONResponse([r for r in _suggest(q.strip(), limit=200) if r["kind"] in kinds][:12])

    register_insight(app, SimpleNamespace(store=store, lock=lock, catalog=catalog, templates=templates, render=render,
                                          db_path=db_path, url=_url, describe_hit=describe_hit))

    @app.get("/api/status")
    def status():
        return {"status": catalog.status, "ready": catalog.ready.is_set()}

    ctx.background = background
    routes_compare.register(app, ctx)
    routes_matrix.register(app, ctx)
    clade_content.register(templates, ctx)
    return app


__all__ = ["Catalog", "create_app"]
