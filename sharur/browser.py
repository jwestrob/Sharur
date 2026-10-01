"""Read-only web browser for a Sharur dataset: cards, neighborhoods, evidence, modules.

``create_app(db_path)`` returns a FastAPI app; ``sharur browse`` serves it.
Pages are server-rendered HTML with inline SVG and no scripts. Every value is
HTML-escaped and no page shows a sequence. Links are stable
(``/protein/<id>``, ``/genome/<id>``), so a URL can be shared with anyone who
can reach the server.

The database is opened read-only and requests are serialized on one
connection, which suits interactive browsing by a few people.
"""

from __future__ import annotations

import html
import secrets
import threading
from pathlib import Path
from typing import Any
from urllib.parse import quote

from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.responses import HTMLResponse, RedirectResponse

from sharur.storage.duckdb_store import DuckDBStore

COOKIE = "sharur_browse"

CSS = """
:root { --bg:#fbfaf7; --fg:#1d232b; --muted:#5d6672; --line:#d9d6cf; --card:#fff; --accent:#2d6a8f;
        --pos:#3f7f5f; --neg:#8a5a2b; --edge:#b23a3a; --chip:#eef2f5; }
@media (prefers-color-scheme: dark) {
  :root { --bg:#14171b; --fg:#e4e6e9; --muted:#9aa3ad; --line:#2c3238; --card:#1b1f24; --accent:#7fb3d5;
          --pos:#7dbf9a; --neg:#d6a46b; --edge:#e27b7b; --chip:#232a31; } }
* { box-sizing:border-box; }
body { margin:0; background:var(--bg); color:var(--fg); font:15px/1.5 system-ui, -apple-system, sans-serif; }
header { border-bottom:1px solid var(--line); padding:10px 16px; display:flex; gap:16px; align-items:center;
         flex-wrap:wrap; }
header a.home { font-weight:700; color:var(--fg); text-decoration:none; }
header form { display:flex; gap:6px; flex:1; min-width:240px; }
input[type=text] { flex:1; padding:6px 8px; border:1px solid var(--line); border-radius:6px; background:var(--card);
                   color:var(--fg); font:inherit; min-width:0; }
button { padding:6px 10px; border:1px solid var(--line); border-radius:6px; background:var(--chip); color:var(--fg);
         font:inherit; cursor:pointer; }
main { max-width:1100px; margin:0 auto; padding:16px; }
h1 { font-size:1.35rem; margin:0.2rem 0 0.6rem; word-break:break-all; }
h2 { font-size:1.05rem; margin:1.4rem 0 0.5rem; }
.card { background:var(--card); border:1px solid var(--line); border-radius:8px; padding:12px 14px; }
.muted { color:var(--muted); }
a { color:var(--accent); }
table { border-collapse:collapse; width:100%; font-size:0.92rem; }
th, td { text-align:left; padding:4px 8px; border-bottom:1px solid var(--line); vertical-align:top; }
th { color:var(--muted); font-weight:600; }
.scroll { overflow-x:auto; }
.chip { display:inline-block; background:var(--chip); border-radius:999px; padding:1px 8px; margin:2px 3px 2px 0;
        font-size:0.85rem; }
.edge { color:var(--edge); font-weight:600; }
code { font-size:0.9em; }
svg text { fill:var(--fg); font:11px system-ui, sans-serif; }
"""


def _e(value: Any) -> str:
    return html.escape("" if value is None else str(value))


def _evalue(value: Any) -> str:
    return f"{value:.1e}" if isinstance(value, float) else _e(value)


def _link(kind: str, ident: str, label: str | None = None) -> str:
    return f'<a href="/{kind}/{quote(ident, safe="")}">{_e(label or ident)}</a>'


def _page(title: str, body: str) -> HTMLResponse:
    return HTMLResponse(f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1"><title>{_e(title)} · Sharur</title>
<style>{CSS}</style></head><body><header><a class="home" href="/">Sharur</a>
<form action="/go" method="get"><input type="text" name="q" placeholder="protein ID, genome ID, predicate, or domain pattern"
 aria-label="Search"><button type="submit">Go</button></form></header><main>{body}</main></body></html>""")


# --------------------------------------------------------------------------- #
# SVG drawings
# --------------------------------------------------------------------------- #


def _domains_svg(length: int | None, domains: list[dict[str, Any]], width: int = 900) -> str:
    if not length or not domains:
        return ""
    scale = (width - 20) / length
    parts = [f'<line x1="10" y1="20" x2="{10 + length * scale:.1f}" y2="20" stroke="currentColor" '
             'stroke-width="2" opacity="0.4"/>']
    for i, d in enumerate(domains):
        x, w = 10 + (d["start_aa"] - 1) * scale, max(2.0, (d["end_aa"] - d["start_aa"] + 1) * scale)
        hue = (sum(map(ord, d["name"])) * 37) % 360
        parts.append(f'<g><title>{_e(d["name"])} {d["start_aa"]}-{d["end_aa"]}</title>'
                     f'<rect x="{x:.1f}" y="10" width="{w:.1f}" height="20" rx="4" fill="hsl({hue},45%,60%)"/></g>')
        if w > 7 * len(d["name"]):
            parts.append(f'<text x="{x + 3:.1f}" y="24">{_e(d["name"])}</text>')
        elif i % 2 == 0 and len(domains) <= 30:
            parts.append(f'<text x="{x:.1f}" y="44">{_e(d["name"])}</text>')
    return (f'<div class="scroll"><svg viewBox="0 0 {width} 50" width="100%" style="min-width:600px" '
            f'role="img" aria-label="Domain architecture">{"".join(parts)}</svg></div>')


def _neighborhood_svg(raw: dict[str, Any], width: int = 1000) -> str:
    genes = raw.get("proteins") or []
    if not genes:
        return ""
    lo, hi = min(g["start"] for g in genes), max(g["end"] for g in genes)
    span = max(1, hi - lo)
    scale = (width - 40) / span
    parts = []
    if raw.get("contig_start_in_window"):
        parts.append('<line x1="14" y1="8" x2="14" y2="58" class="edge" stroke="currentColor" stroke-dasharray="3 3"/>'
                     '<text x="16" y="12" class="edge">contig start</text>')
    if raw.get("contig_end_in_window"):
        parts.append(f'<line x1="{width - 14}" y1="8" x2="{width - 14}" y2="58" stroke="currentColor" '
                     f'stroke-dasharray="3 3"/><text x="{width - 80}" y="12">contig end</text>')
    for g in genes:
        x1, x2 = 20 + (g["start"] - lo) * scale, 20 + (g["end"] - lo) * scale
        head = min(8.0, (x2 - x1) / 2)
        y, h = 34, 9
        if g["strand"] == "-":
            points = f"{x2},{y - h} {x1 + head},{y - h} {x1},{y} {x1 + head},{y + h} {x2},{y + h}"
        else:
            points = f"{x1},{y - h} {x2 - head},{y - h} {x2},{y} {x2 - head},{y + h} {x1},{y + h}"
        fill = "var(--accent)" if g.get("is_anchor") else (
            "var(--edge)" if g.get("edge_status") == "truncated" else "var(--chip)")
        label = g.get("annotation") or ""
        parts.append(f'<a href="/protein/{quote(g["protein_id"], safe="")}"><g><title>{_e(g["protein_id"])} '
                     f'{_e(label)}</title><polygon points="{points}" fill="{fill}" stroke="currentColor" '
                     f'stroke-width="0.6"/></g></a>')
        if x2 - x1 > 40 and label:
            short = label.split(":")[-1].split(" (")[0][:max(3, int((x2 - x1) / 7))]
            parts.append(f'<text x="{x1:.1f}" y="{y + 24}">{_e(short)}</text>')
    return (f'<div class="scroll"><svg viewBox="0 0 {width} 66" width="100%" style="min-width:700px" '
            f'role="img" aria-label="Gene neighborhood">{"".join(parts)}</svg></div>')


# --------------------------------------------------------------------------- #
# App
# --------------------------------------------------------------------------- #


def create_app(db_path: str | Path, *, token: str | None = None) -> FastAPI:
    """Browser app over one dataset, opened read-only."""
    from sharur.abundance import default_path  # noqa: PLC0415

    store = DuckDBStore(str(db_path), read_only=True)
    lock = threading.Lock()
    sidecar = default_path(db_path)
    app = FastAPI(title="Sharur browser", docs_url=None, redoc_url=None, openapi_url=None)

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
            return HTMLResponse("Open the link with ?token=… to use this browser.", status_code=401)
        return await call_next(request)

    @app.get("/", response_class=HTMLResponse)
    def home():
        from sharur.operators.introspection import describe_dataset  # noqa: PLC0415

        with lock:
            info = describe_dataset(store)
            top_bins = store.execute(
                "SELECT bin_id, COUNT(*) FROM proteins GROUP BY 1 ORDER BY 2 DESC, 1 LIMIT 25")
        sources = "".join(f"<tr><td>{_e(s['source'])}</td><td>{s['proteins']:,}</td>"
                          f"<td>{s['protein_fraction']:.1%}</td></tr>" for s in info["annotation_sources"])
        callers = "".join(f"<tr><td>{_e(c['table'])}</td><td>{_e(c['meaning'])}</td><td>{c['rows']:,}</td></tr>"
                          for c in info["curated_callers"])
        genomes = "".join(f"<span class='chip'>{_link('genome', b)} · {n:,}</span>" for b, n in top_bins)
        return _page("Dataset", f"""
<h1>Dataset overview</h1>
<p class="card">{info['proteins']:,} proteins in {info['bins']['count']:,} genomes ·
predicate maps: {_e(info['predicates']['map_status'])} ·
KEGG modules: {'available' if info.get('local_kegg_build') else 'not built (sharur setup-kegg)'}</p>
<p class="muted">Search above by protein or genome ID, a predicate such as <code>nife_group3</code>, or a domain
pattern such as <code>TPR_* {{10,}}</code>.</p>
<h2>Annotation sources</h2><div class="scroll"><table><tr><th>Source</th><th>Proteins</th><th>Coverage</th></tr>
{sources}</table></div>
<h2>Curated callers</h2><div class="scroll"><table><tr><th>Table</th><th>Meaning</th><th>Rows</th></tr>{callers}</table></div>
<h2>Largest genomes</h2><p>{genomes}</p>""")

    @app.get("/go")
    def go(q: str = Query("", max_length=500)):
        q = q.strip()
        if not q:
            return RedirectResponse("/", status_code=303)
        with lock:
            if store.execute("SELECT 1 FROM proteins WHERE protein_id = ?", [q]):
                return RedirectResponse(f"/protein/{quote(q, safe='')}", status_code=303)
            if store.execute("SELECT 1 FROM bins WHERE bin_id = ?", [q]):
                return RedirectResponse(f"/genome/{quote(q, safe='')}", status_code=303)
        from sharur.predicates.vocabulary import PREDICATE_BY_ID  # noqa: PLC0415

        if q in PREDICATE_BY_ID:
            return RedirectResponse(f"/predicate/{quote(q, safe='')}", status_code=303)
        return RedirectResponse(f"/architecture?pattern={quote(q, safe='')}", status_code=303)

    @app.get("/protein/{protein_id:path}/why/{predicate}", response_class=HTMLResponse)
    def why_page(protein_id: str, predicate: str):
        from sharur.operators.cards import why  # noqa: PLC0415

        with lock:
            result = why(store, protein_id, predicate)
        rows = []
        for p in result.get("paths", []):
            mapping = p.get("mapping") or {}
            chain = " › ".join(mapping.get("chain", []))
            rows.append(f"<tr><td>{_e(p.get('source_db'))}</td><td>{_e(p.get('accession'))} "
                        f"{_e(p.get('annotation'))}</td><td>{_e(p.get('relation'))}</td>"
                        f"<td>{_e(mapping.get('evidence'))}</td><td>{_e(chain)}</td></tr>")
        definition = result.get("definition") or {}
        status = "carries" if result.get("present") else "does not carry"
        table = ("<div class='scroll'><table><tr><th>Source</th><th>Hit</th><th>Relation</th><th>Map evidence</th>"
                 f"<th>Is-a chain</th></tr>{''.join(rows)}</table></div>") if rows else ""
        return _page(f"{predicate} · {protein_id}", f"""
<h1>Why {_link('protein', protein_id)} {status} <code>{_e(predicate)}</code></h1>
<p class="muted">{_e(definition.get('description'))}</p>{table}""")

    @app.get("/protein/{protein_id:path}", response_class=HTMLResponse)
    def protein_page(protein_id: str):
        from sharur.architecture import architecture  # noqa: PLC0415
        from sharur.contig_context import EdgeContext, describe_edge  # noqa: PLC0415
        from sharur.operators.cards import card  # noqa: PLC0415
        from sharur.operators.navigation import get_neighborhood  # noqa: PLC0415

        with lock:
            c = card(store, protein_id)
            if not c.get("found"):
                raise HTTPException(404, "Protein not found")
            domains = [d.to_dict() for d in architecture(store, protein_id)]
            hood = get_neighborhood(store, protein_id, window=8).raw or {}
            nearby = []
            try:
                from sharur.modules import locus_modules  # noqa: PLC0415

                nearby = [m for m in locus_modules(store, protein_id, window=10) if m["steps_complete"] > 0][:8]
            except FileNotFoundError:
                pass
        loc, genome = c["location"], c["genome"]
        edge = describe_edge(EdgeContext(**c["contig_edge"])) if c.get("contig_edge") else ""
        quality = (f" · {genome['completeness']}% complete, {genome['contamination']}% contamination"
                   if genome.get("completeness") is not None else "")
        annotations = "".join(
            f"<tr><td>{_e(src)}</td><td>{_e(h['accession'])}</td><td>{_e(h['name'])}</td>"
            f"<td>{_e(h.get('description'))}</td><td>{_evalue(h.get('evalue'))}</td></tr>"
            for src, hits in c["annotations"].items() for h in hits)
        predicates = "".join(
            f"<p><b>{_e(facet)}</b>: " + "".join(
                f"<span class='chip'><a href='/protein/{quote(protein_id, safe='')}/why/{quote(p, safe='')}'>"
                f"{_e(p)}</a></span>" for p in values) + "</p>"
            for facet, values in c["predicates"].items() if values)
        systems = "".join(f"<span class='chip'>{_e(s['source'])}: {_e(s['system_id'])}</span>"
                          for s in c["validated_systems"])
        modules = "".join(f"<tr><td>{_e(m['module'])}</td><td>{_e(m['name'])}</td>"
                          f"<td>{m['steps_complete']}/{m['steps_total']}</td></tr>" for m in nearby)
        hydrogenase = c.get("hydrogenase_classification")
        return _page(protein_id, f"""
<h1>{_e(protein_id)}</h1>
<div class="card">{_e(loc['length_aa'])} aa · {_e(loc['contig_id'])}:{_e(loc['start'])}-{_e(loc['end'])}
({_e(loc['strand'])}) · genome {_link('genome', genome['bin_id'])}{_e(': ' + genome['taxonomy'] if genome.get('taxonomy') else '')}{_e(quality)}
<br><span class="{'edge' if c.get('contig_edge') and c['contig_edge']['edge_status'] in ('truncated', 'contig_edge') else 'muted'}">{_e(edge)}</span></div>
<h2>Domains</h2>{_domains_svg(loc['length_aa'], domains) or '<p class="muted">No placed domain hits.</p>'}
<p class="muted">{_e(c.get('architecture'))}</p>
<h2>Neighborhood</h2>{_neighborhood_svg(hood)}
<p class="muted">Anchor in blue; genes that run off a contig end in red. Click a gene to open it.</p>
{('<h2>Validated systems</h2><p>' + systems + '</p>') if systems else ''}
{('<h2>Hydrogenase classification</h2><p>' + _e(hydrogenase.get('outcome')) + ': ' + _e(hydrogenase.get('reference_label')) + '</p>') if hydrogenase else ''}
<h2>Predicates</h2><div class="card">{predicates or '<span class="muted">None</span>'}</div>
<p class="muted">Click a predicate for the evidence behind it.</p>
<h2>Annotations</h2><div class="scroll"><table><tr><th>Source</th><th>Accession</th><th>Name</th><th>Description</th>
<th>E-value</th></tr>{annotations}</table></div>
{('<h2>KEGG module steps nearby</h2><div class="scroll"><table><tr><th>Module</th><th>Name</th><th>Steps found within 10 genes</th></tr>' + modules + '</table></div>') if modules else ''}""")

    @app.get("/genome/{bin_id:path}", response_class=HTMLResponse)
    def genome_page(bin_id: str):
        with lock:
            bins = store.execute("SELECT * FROM bins WHERE bin_id = ?", [bin_id])
            if not bins:
                raise HTTPException(404, "Genome not found")
            columns = [r[0] for r in store.execute(
                "SELECT column_name FROM information_schema.columns WHERE table_name = 'bins' ORDER BY ordinal_position")]
            row = dict(zip(columns, bins[0], strict=True))
            stats = store.execute(
                """SELECT COUNT(DISTINCT p.contig_id), COUNT(*), COUNT(annotated.protein_id)
                   FROM proteins p
                   LEFT JOIN (SELECT DISTINCT a.protein_id FROM annotations a JOIN proteins q USING (protein_id)
                              WHERE q.bin_id = ?) annotated USING (protein_id)
                   WHERE p.bin_id = ?""", [bin_id, bin_id])[0]
            giants = store.execute(
                "SELECT protein_id, sequence_length FROM proteins WHERE bin_id = ? ORDER BY sequence_length DESC LIMIT 10",
                [bin_id])
            modules: list[dict[str, Any]] = []
            try:
                from sharur.modules import genome_modules  # noqa: PLC0415

                modules = genome_modules(store, bins=[bin_id], min_completeness=0.5)
            except FileNotFoundError:
                pass
        abundance_rows = []
        if sidecar.is_file():
            from sharur.abundance import genome_abundance  # noqa: PLC0415

            abundance_rows = genome_abundance(sidecar, bins=[bin_id])
        facts = " · ".join(f"{k}: {_e(row[k])}" for k in ("taxonomy", "completeness", "contamination")
                           if row.get(k) is not None)
        module_rows = "".join(
            f"<tr><td>{_e(m['module'])}</td><td>{_e(m['name'])}</td><td>{m['steps_complete']}/{m['steps_total']}</td>"
            f"<td>{'<span class=edge>found genes at a contig end</span>' if m.get('found_at_contig_edge') else ''}"
            f"</td></tr>" for m in sorted(modules, key=lambda m: -m["completeness"]))
        abundance = "".join(
            f"<tr><td>{_e(a['sample_id'])}</td><td>{(a['relative_abundance'] or 0):.2%}</td>"
            f"<td>{(a['mean_depth'] or 0):.1f}x</td></tr>" for a in abundance_rows)
        largest = "".join(f"<span class='chip'>{_link('protein', p)} · {n:,} aa</span>" for p, n in giants)
        return _page(bin_id, f"""
<h1>{_e(bin_id)}</h1>
<div class="card">{facts or '<span class="muted">No genome metadata</span>'}<br>
{stats[0]:,} contigs · {stats[1]:,} proteins · {stats[2] or 0:,} with annotation hits</div>
<h2>Largest proteins</h2><p>{largest}</p>
{('<h2>Abundance</h2><div class="scroll"><table><tr><th>Sample</th><th>Share of mapped reads</th><th>Mean depth</th></tr>' + abundance + '</table></div>') if abundance else ''}
<h2>KEGG modules ≥ 50% complete</h2>
{('<div class="scroll"><table><tr><th>Module</th><th>Name</th><th>Steps</th><th></th></tr>' + module_rows + '</table></div>') if module_rows else '<p class="muted">None, or KEGG modules not built (sharur setup-kegg).</p>'}""")

    @app.get("/predicate/{predicate}", response_class=HTMLResponse)
    def predicate_page(predicate: str, limit: int = Query(100, ge=1, le=1000)):
        with lock:
            rows = store.execute(
                """SELECT pp.protein_id, p.bin_id, p.sequence_length FROM protein_predicates pp
                   JOIN proteins p USING (protein_id) WHERE list_contains(pp.predicates, ?)
                   ORDER BY p.bin_id, pp.protein_id LIMIT ?""", [predicate, limit])
            total = store.execute(
                "SELECT COUNT(*) FROM protein_predicates WHERE list_contains(predicates, ?)", [predicate])[0][0]
        items = "".join(f"<tr><td>{_link('protein', pid)}</td><td>{_link('genome', b)}</td><td>{n}</td></tr>"
                        for pid, b, n in rows)
        return _page(predicate, f"""<h1>Proteins with <code>{_e(predicate)}</code></h1>
<p class="muted">{total:,} proteins; showing {len(rows):,}.</p>
<div class="scroll"><table><tr><th>Protein</th><th>Genome</th><th>Length (aa)</th></tr>{items}</table></div>""")

    @app.get("/architecture", response_class=HTMLResponse)
    def architecture_page(pattern: str = Query(..., max_length=500), limit: int = Query(100, ge=1, le=1000)):
        from sharur.architecture import PatternError, search_architecture  # noqa: PLC0415

        try:
            with lock:
                result = search_architecture(store, pattern, limit=limit)
        except PatternError as exc:
            return _page("Pattern", f"<h1>Domain pattern</h1><p class='edge'>{_e(exc)}</p>"
                         "<p class='muted'>See the architecture search guide for the pattern language.</p>")
        items = "".join(
            f"<tr><td>{_link('protein', r['protein_id'])}</td><td>{_link('genome', r['bin_id'] or '')}</td>"
            f"<td>{_e(r['length_aa'])}</td><td>{_e(r['architecture'])}</td></tr>" for r in result["records"])
        return _page("Architecture search", f"""<h1>Domain pattern <code>{_e(pattern)}</code></h1>
<p class="muted">{result['total']:,} proteins match; showing {len(result['records']):,}.</p>
<div class="scroll"><table><tr><th>Protein</th><th>Genome</th><th>aa</th><th>Architecture</th></tr>{items}</table></div>""")

    app.state.store = store
    return app
