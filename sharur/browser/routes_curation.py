"""Curation: notes and flags, triage, the collection cart, and hover previews."""

from __future__ import annotations

import json
from typing import Any
from urllib.parse import quote, unquote

from fastapi import HTTPException, Query, Request
from fastapi.responses import HTMLResponse, JSONResponse, PlainTextResponse, RedirectResponse

from sharur.browser.notes import FLAGS, KINDS, NotesStore, default_notes_path
from sharur.browser.routes_loci import build_rows, render_stack

USER_COOKIE = "sharur_user"


def _author(request: Request) -> str:
    name = (request.cookies.get(USER_COOKIE) or "").strip()
    return name[:60] or "anonymous"


def entity_url(ctx, kind: str, entity: str) -> str:
    if kind == "system":
        call = ctx.system_by_id.get(entity)
        return ctx.url("system", call["kind"], call["type"]) + "/loci" if call else "#"
    if kind == "clade":
        rank, _, name = entity.partition(":")
        return ctx.url("taxa", rank, name)
    if kind == "module":
        return ctx.url("pathway", entity)
    if kind == "crispr":
        return ctx.url("crispr", entity)
    return ctx.url(kind, entity)


def entity_label(ctx, kind: str, entity: str) -> str:
    catalog = ctx.catalog
    if kind == "domain" and entity in catalog.domains:
        return f'{catalog.domains[entity]["name"]} ({entity})'
    if kind == "vog" and entity in catalog.vogs:
        return f'{catalog.vogs[entity]["name"]} ({entity})'
    if kind == "function" and entity in ctx.predicates:
        return ctx.predicates[entity].name
    if kind == "system" and entity in ctx.system_by_id:
        call = ctx.system_by_id[entity]
        return f'{call["type"]} in {call["bin_id"]}'
    return entity


def register(app, ctx) -> None:
    notes = NotesStore(ctx.notes_path or default_notes_path(ctx.db_path))
    ctx.notes = notes
    ctx.system_by_id = {s["system_id"]: s for s in ctx.catalog.systems}
    catalog = ctx.catalog

    # ------------------------------------------------------------------ #
    # Notes API
    # ------------------------------------------------------------------ #

    @app.get("/api/notes")
    def get_notes(kind: str, id: str):  # noqa: A002 - query parameter name
        if kind not in KINDS:
            raise HTTPException(400, "unknown kind")
        return notes.state(kind, id)

    @app.post("/api/notes")
    async def post_note(request: Request):
        body = await request.json()
        kind, entity = body.get("kind", ""), str(body.get("id", ""))
        author = _author(request)
        try:
            if body.get("flag"):
                return notes.toggle_flag(kind, entity, body["flag"], author)
            return notes.add_note(kind, entity, str(body.get("text", "")), author)
        except ValueError as exc:
            raise HTTPException(400, str(exc)) from exc

    @app.post("/api/notes/{note_id}/delete")
    def delete_note(request: Request, note_id: int):
        target = notes.delete(note_id, _author(request))
        if target is None:
            raise HTTPException(403, "Only the author can delete a note")
        return notes.state(*target)

    @app.post("/api/me")
    async def set_me(request: Request):
        body = await request.json()
        name = str(body.get("name", "")).strip()[:60]
        response = JSONResponse({"name": name or "anonymous"})
        response.set_cookie(USER_COOKIE, name, max_age=60 * 60 * 24 * 365, samesite="strict")
        return response

    # ------------------------------------------------------------------ #
    # Flags overview and exports
    # ------------------------------------------------------------------ #

    def _rows(flag: str, kind: str, author: str) -> list[dict[str, Any]]:
        rows = notes.listing(flag=flag, kind=kind, author=author)
        for r in rows:
            r["url"] = entity_url(ctx, r["kind"], r["entity"])
            r["label"] = entity_label(ctx, r["kind"], r["entity"])
        return rows

    @app.get("/flags", response_class=HTMLResponse)
    def flags_page(request: Request, flag: str = Query(""), kind: str = Query(""), author: str = Query("")):
        rows = _rows(flag, kind, author)
        authors = sorted({r["author"] for r in notes.listing()})
        counts = {f: len(notes.flagged_entities(f)) for f in FLAGS}
        return ctx.render(request, "flags.html", "flags", rows=rows[:1000], total=len(rows), flag=flag, kind=kind,
                          author=author, authors=authors, counts=counts, FLAGS=FLAGS, KINDS=KINDS,
                          notes_path=str(notes.path))

    @app.get("/flags.{fmt}")
    def flags_export(fmt: str, flag: str = Query(""), kind: str = Query(""), author: str = Query("")):
        rows = _rows(flag, kind, author)
        fields = ["id", "kind", "entity", "label", "flag", "text", "author", "created_at"]
        if fmt == "jsonl":
            return PlainTextResponse("".join(json.dumps({k: r.get(k) for k in fields}) + "\n" for r in rows),
                                     headers={"Content-Disposition": 'attachment; filename="sharur_notes.jsonl"'})
        if fmt == "tsv":
            lines = ["\t".join(fields)] + ["\t".join(str(r.get(k) or "").replace("\t", " ").replace("\n", " ")
                                                     for k in fields) for r in rows]
            return PlainTextResponse("\n".join(lines) + "\n",
                                     headers={"Content-Disposition": 'attachment; filename="sharur_notes.tsv"'})
        raise HTTPException(404, "Use .tsv or .jsonl")

    # ------------------------------------------------------------------ #
    # Triage
    # ------------------------------------------------------------------ #

    def triage_items(source: str, kind: str, ident: str, subtype: str, flag: str) -> tuple[list[tuple[str, str]], str, str]:
        """Ordered (entity_kind, entity) items, a title and a back link."""
        if source == "system":
            calls = [s for s in catalog.systems if s["kind"] == kind and s["type"] == ident
                     and (not subtype or (s["subtype"] or "") == subtype)]
            return [("system", s["system_id"]) for s in calls], f"{ident} calls", ctx.url("system", kind, ident)
        if source == "family":
            with ctx.lock:
                if kind == "domain":
                    rows = ctx.store.execute(
                        "SELECT DISTINCT protein_id FROM annotations WHERE LOWER(source) = 'pfam' "
                        "AND split_part(accession, '.', 1) = ? ORDER BY 1", [ident])
                elif kind == "vog":
                    rows = ctx.store.execute(
                        "SELECT DISTINCT protein_id FROM annotations WHERE LOWER(source) IN ('vogdb', 'vog') "
                        "AND accession = ? ORDER BY 1", [ident])
                elif kind == "function":
                    rows = ctx.store.execute(
                        "SELECT protein_id FROM protein_predicates WHERE list_contains(predicates, ?) ORDER BY 1",
                        [ident])
                else:
                    raise HTTPException(404, "Unknown family kind")
            back = ctx.url({"domain": "domain", "vog": "vog", "function": "function"}[kind], ident)
            return [("protein", r[0]) for r in rows], f"{entity_label(ctx, kind, ident)} carriers", back
        if source == "flags":
            if flag not in FLAGS:
                raise HTTPException(404, "Unknown flag")
            return notes.flagged_entities(flag, kind), f"Flagged {FLAGS[flag][0].lower()}", "/flags"
        raise HTTPException(404, "Unknown triage source")

    @app.get("/triage", response_class=HTMLResponse)
    def triage(request: Request, source: str, kind: str = Query(""), id: str = Query(""),  # noqa: A002
               subtype: str = Query(""), flag: str = Query(""), i: int = Query(0, ge=0)):
        items, title, back = triage_items(source, kind, id, subtype, flag)
        if not items:
            return ctx.render(request, "triage.html", "flags", title=title, back=back, item=None, total=0, i=0,
                              FLAGS=FLAGS, base="")
        i = min(i, len(items) - 1)
        entity_kind, entity = items[i]
        base = f"/triage?source={quote(source)}&kind={quote(kind)}&id={quote(id)}&subtype={quote(subtype)}&flag={quote(flag)}"
        detail = _triage_detail(entity_kind, entity)
        return ctx.render(request, "triage.html", "flags", title=title, back=back, total=len(items), i=i,
                          item={"kind": entity_kind, "entity": entity, "label": entity_label(ctx, entity_kind, entity),
                                "url": entity_url(ctx, entity_kind, entity), **detail},
                          base=base, FLAGS=FLAGS)

    def _triage_detail(entity_kind: str, entity: str) -> dict[str, Any]:
        from sharur.architecture import architecture, compact  # noqa: PLC0415

        if entity_kind == "system" and entity in ctx.system_by_id:
            call = ctx.system_by_id[entity]
            with ctx.lock:
                profiles = {p: (name or "").split("__")[-1] or None for _, p, name in ctx.store.execute(
                    "SELECT system_id, protein_id, profile_name FROM system_proteins WHERE system_id = ?", [entity])}
            members = set(call["proteins"]) or set(profiles)
            anchor = next(iter(sorted(members)), None)
            g = catalog.by_bin.get(call["bin_id"])
            rows = build_rows(ctx, [{"anchor": anchor, "members": members, "profiles": profiles,
                                     "bin_id": call["bin_id"], "title": call["subtype"] or call["type"],
                                     "genome_label": g.label if g else ""}], 6) if anchor else []
            svgs, legend = render_stack(rows)
            return {"svg": svgs[0] if svgs else "", "legend": legend,
                    "facts": [("Genome", call["bin_id"]), ("Lineage", g.label if g else ""),
                              ("Subtype", call["subtype"] or "–"), ("Genes", str(call["genes"] or len(members)))],
                    "genome": call["bin_id"]}
        if entity_kind == "protein":
            with ctx.lock:
                row = ctx.store.execute("SELECT bin_id, sequence_length FROM proteins WHERE protein_id = ?", [entity])
                arch = compact([d.name for d in architecture(ctx.store, entity)])
            bin_id = row[0][0] if row else None
            g = catalog.by_bin.get(bin_id)
            rows = build_rows(ctx, [{"anchor": entity, "members": {entity}, "profiles": {}, "bin_id": bin_id,
                                     "title": "", "genome_label": g.label if g else ""}], 6)
            svgs, legend = render_stack(rows, member_label="this protein")
            return {"svg": svgs[0] if svgs else "", "legend": legend,
                    "facts": [("Genome", bin_id or ""), ("Lineage", g.label if g else ""),
                              ("Length", f"{row[0][1]:,} aa" if row and row[0][1] else "–"),
                              ("Pfam architecture", arch or "no placed domains")], "genome": bin_id}
        return {"svg": "", "legend": [], "facts": [], "genome": None}

    # ------------------------------------------------------------------ #
    # Collection, previews
    # ------------------------------------------------------------------ #

    @app.get("/collection", response_class=HTMLResponse)
    def collection(request: Request):
        return ctx.render(request, "collection.html", "collection")

    @app.post("/api/collection")
    async def collection_rows(request: Request):
        body = await request.json()
        proteins = [str(p) for p in body.get("proteins", [])][:5000]
        genomes = [str(g) for g in body.get("genomes", [])][:5000]
        rows = []
        if proteins:
            with ctx.lock:
                found = ctx.store.execute(
                    """WITH p AS (SELECT protein_id, bin_id, sequence_length FROM proteins
                                  WHERE protein_id IN (SELECT UNNEST(?::VARCHAR[]))),
                            best AS (SELECT a.protein_id, ARG_MIN(COALESCE(NULLIF(a.name, ''), a.accession),
                                                                  COALESCE(a.evalue, 1)) AS top
                                     FROM annotations a JOIN p USING (protein_id) GROUP BY 1)
                       SELECT p.protein_id, p.bin_id, p.sequence_length, best.top FROM p LEFT JOIN best USING (protein_id)""",
                    [proteins])
            for pid, bin_id, length, top in found:
                g = catalog.by_bin.get(bin_id)
                rows.append({"kind": "protein", "id": pid, "url": ctx.url("protein", pid), "genome": bin_id,
                             "lineage": g.label if g else "", "length": length, "best": ctx.describe_hit(top) or ""})
        for bin_id in genomes:
            g = catalog.by_bin.get(bin_id)
            if g:
                rows.append({"kind": "genome", "id": bin_id, "url": ctx.url("genome", bin_id), "genome": bin_id,
                             "lineage": g.label, "length": g.length, "best": f"{g.proteins:,} proteins"})
        return rows

    @app.post("/api/fasta")
    async def fasta_batch(request: Request):
        form = await request.form()
        ids = [i for i in str(form.get("ids", "")).split("\n") if i.strip()][:5000]
        with ctx.lock:
            rows = ctx.store.execute("SELECT protein_id, bin_id, sequence FROM proteins "
                                     "WHERE protein_id IN (SELECT UNNEST(?::VARCHAR[]))", [[i.strip() for i in ids]])
        body = "".join(f">{pid} genome={bin_id}\n" + "\n".join((seq or "").rstrip("*")[i:i + 60]
                                                               for i in range(0, len((seq or "").rstrip("*")), 60)) + "\n"
                       for pid, bin_id, seq in rows if seq)
        return PlainTextResponse(body, headers={"Content-Disposition": 'attachment; filename="sharur_collection.faa"'})

    @app.get("/api/preview", response_class=HTMLResponse)
    def preview(request: Request, href: str):
        from sharur.architecture import architecture, compact  # noqa: PLC0415

        path = unquote(href.split("?")[0])
        if path.startswith("/protein/") and "/why/" not in path:
            pid = path[len("/protein/"):]
            with ctx.lock:
                row = ctx.store.execute("SELECT bin_id, sequence_length FROM proteins WHERE protein_id = ?", [pid])
                if not row:
                    raise HTTPException(404)
                arch = compact([d.name for d in architecture(ctx.store, pid)])
                labels = ctx.store.execute("SELECT predicates FROM protein_predicates WHERE protein_id = ?", [pid])
            g = catalog.by_bin.get(row[0][0])
            names = [ctx.predicates[p].name for p in (labels[0][0] if labels else []) or []
                     if p in ctx.predicates and ctx.predicates[p].category not in ("annotation", "size", "composition")][:6]
            return ctx.render(request, "_preview.html", "", kind="protein", title=pid.rsplit("|", 1)[-1],
                              sub=g.label if g else "", facts=[("Length", f"{row[0][1]:,} aa" if row[0][1] else "–"),
                                                               ("Architecture", arch or "no placed domains")],
                              chips=names, flags=notes.state("protein", pid)["flags"])
        if path.startswith("/genome/"):
            bin_id = path[len("/genome/"):].split("/")[0]
            g = catalog.by_bin.get(bin_id)
            if g is None:
                raise HTTPException(404)
            return ctx.render(request, "_preview.html", "", kind="genome", title=bin_id,
                              sub=" › ".join(n for _, n in g.lineage[1:]),
                              facts=[("Size", f"{g.length / 1e6:.2f} Mb"), ("Contigs", f"{g.contigs:,}"),
                                     ("Proteins", f"{g.proteins:,}")], chips=[],
                              flags=notes.state("genome", bin_id)["flags"])
        raise HTTPException(404)

    @app.get("/me")
    def me(request: Request):
        return RedirectResponse("/flags")
