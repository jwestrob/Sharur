"""Dataset health in the browser: ``/health``, ``/api/health`` and a rail pill when checks need care.

The report comes from :func:`sharur.health.run_checks`, run once in the
background after the catalog loads (about a second on a few million proteins)
and kept; ``/health?refresh=1`` runs it again.
"""

from __future__ import annotations

import re
import threading
from types import SimpleNamespace

from fastapi import FastAPI, Query, Request
from fastapi.responses import HTMLResponse, JSONResponse
from markupsafe import Markup, escape

_CODE = re.compile(r"`([^`]+)`")


def code_html(text: str) -> Markup:
    """Escape ``text`` and set its backticked spans as code."""
    return Markup(_CODE.sub(lambda m: f"<code>{m.group(1)}</code>", str(escape(text or ""))))


def example_href(url, example: dict[str, str]) -> str | None:
    kind = example.get("kind")
    if kind in ("genome", "protein", "contig"):
        return url(kind, example["id"])
    return None


def register(app: FastAPI, ctx: SimpleNamespace) -> None:
    """Add /health and /api/health. ``ctx``: store, lock, catalog, render, url, db_path, background."""
    from sharur.health import run_checks  # noqa: PLC0415

    state: dict[str, object] = {"report": None}
    running = threading.Lock()

    def query(sql: str, params: list | None = None) -> list:
        with ctx.lock:
            return ctx.store.execute(sql, params or [])

    def compute(force: bool = False):
        with running:
            if state["report"] is None or force:
                assemblies = getattr(getattr(ctx, "assemblies", None), "paths", None)
                state["report"] = run_checks(ctx.db_path, query=query, assemblies=assemblies)
        return state["report"]

    def after_catalog() -> None:
        ctx.catalog.ready.wait()
        try:
            compute()
        except Exception:  # pragma: no cover - the page reruns it and shows the error
            pass

    if getattr(ctx, "background", True):
        threading.Thread(target=after_catalog, daemon=True).start()

    ctx.health_report = lambda: state["report"]
    templates = getattr(ctx, "templates", None)
    if templates is not None:
        templates.env.globals["health_report"] = ctx.health_report

    @app.get("/health", response_class=HTMLResponse)
    def health_page(request: Request, refresh: int = Query(0)):
        report = compute(force=bool(refresh))
        return ctx.render(request, "health.html", "health", report=report, code_html=code_html,
                          example_href=lambda e: example_href(ctx.url, e))

    @app.get("/api/health")
    def health_api(refresh: int = Query(0)):
        return JSONResponse(compute(force=bool(refresh)).to_dict())
