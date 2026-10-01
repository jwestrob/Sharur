"""Functional landscape: every genome placed by what it encodes.

The map is a principal-component analysis of genome × feature presence (KEGG
orthologs by default; Pfam families or function labels on request). The
presence matrix is centred per feature, and its two leading components come
from a truncated SVD (``scipy.sparse.linalg.svds``) on the sparse matrix, so
genomes sharing many features sit close together. Features carried by a
single genome, or by every genome, carry no contrast and are left out.

Presence counts grow with genome size and completeness, so the page reports
how strongly each axis follows the number of features per genome and genome
completeness, and lists the features that load each end of each axis.

The genome page's closest genomes rank every other genome by the Jaccard
similarity of their feature sets (shared / either), computed from the same
matrix.
"""

from __future__ import annotations

import threading
import time
from dataclasses import dataclass, field
from types import SimpleNamespace
from typing import Any

import numpy as np
from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.responses import HTMLResponse, JSONResponse

from sharur.browser.catalog import RANKS, UNCLASSIFIED

KINDS = {"ko": "KEGG orthologs", "pfam": "Pfam families", "function": "Function labels"}
FEATURE_URLS = {"ko": "/search?q={id}", "pfam": "/domain/{id}", "function": "/function/{id}"}
IN_TEXT = {"ko": "KEGG orthologs", "pfam": "Pfam families", "function": "function labels"}
INLINE_SELECTION = 450   # characters: longer selections are stored and passed as a short key
TOP_LOADINGS = 6


@dataclass
class Embedding:
    kind: str
    coords: np.ndarray                    # genomes × 2; NaN for genomes without features
    explained: list[float]                # share of total variance per axis
    n_features: np.ndarray                # features per genome
    used_features: int                    # features with contrast (in 2+ genomes, not all)
    axis_size_r: list[float]              # Pearson r of each axis with features per genome
    axis_completeness_r: list[float | None]
    loadings: list[dict[str, list[dict[str, Any]]]]
    seconds: float
    matrix: Any = field(repr=False)       # genomes × features, binary CSR (every feature of the kind)
    note: str = ""


def presence(fs, n_genomes: int):
    """Binary genome × feature CSR matrix from a matrix feature set."""
    from scipy.sparse import csr_matrix  # noqa: PLC0415

    keep = fs.value >= fs.threshold
    rows, cols = fs.bin_idx[keep], fs.feat_idx[keep]
    m = csr_matrix((np.ones(len(rows), dtype=np.float64), (rows, cols)), shape=(n_genomes, len(fs.ids)))
    m.sum_duplicates()
    m.data[:] = 1.0
    return m


def _pearson(a: np.ndarray, b: np.ndarray) -> float | None:
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3 or np.std(a[ok]) == 0 or np.std(b[ok]) == 0:
        return None
    return float(np.corrcoef(a[ok], b[ok])[0, 1])


def embed(fs, catalog) -> Embedding:
    """PCA (two components) of centred genome × feature presence."""
    from scipy.sparse.linalg import LinearOperator, svds  # noqa: PLC0415

    t0 = time.perf_counter()
    n = len(catalog.genomes)
    full = presence(fs, n)
    sizes = np.asarray(full.sum(axis=1)).ravel()
    rows = np.flatnonzero(sizes > 0)
    a = full[rows]
    counts = np.asarray(a.sum(axis=0)).ravel()
    cols = np.flatnonzero((counts >= 2) & (counts < len(rows)))
    coords = np.full((n, 2), np.nan)
    if len(rows) < 4 or len(cols) < 3:
        return Embedding(fs.kind, coords, [0.0, 0.0], sizes, int(len(cols)), [None, None], [None, None],
                         [{"high": [], "low": []}, {"high": [], "low": []}], time.perf_counter() - t0, full,
                         note="Too few genomes or shared features to place on a map.")
    a = a[:, cols].tocsr()
    at = a.T.tocsr()
    m, q = a.shape
    mu = np.asarray(a.mean(axis=0)).ravel()

    def matvec(v):
        v = np.ravel(v)
        return a @ v - float(mu @ v)

    def rmatvec(u):
        u = np.ravel(u)
        return at @ u - mu * float(u.sum())

    op = LinearOperator((m, q), matvec=matvec, rmatvec=rmatvec, dtype=np.float64)
    k = min(3, min(m, q) - 1)
    v0 = np.random.default_rng(0).standard_normal(min(m, q))
    # macOS Accelerate BLAS raises spurious floating-point warnings inside sparse matmul;
    # the result matches a dense SVD, and is checked for finite values below
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        u, s, vt = svds(op, k=k, v0=v0)
    if not (np.all(np.isfinite(u)) and np.all(np.isfinite(s))):
        raise ArithmeticError("The SVD of the presence matrix did not converge to finite values.")
    order = np.argsort(-s)[:2]
    u, s, vt = u[:, order], s[order], vt[order]
    # sign: the feature with the largest loading points to positive
    for i in range(len(s)):
        if vt[i, np.argmax(np.abs(vt[i]))] < 0:
            u[:, i] *= -1
            vt[i] *= -1
    col_counts = counts[cols]
    total = float(np.sum(col_counts - col_counts ** 2 / m))
    explained = [float(x ** 2 / total) if total else 0.0 for x in s]
    coords[rows] = u * s
    completeness = np.array([g.completeness if g.completeness is not None else np.nan for g in catalog.genomes])
    size_r = [_pearson(coords[rows, i], sizes[rows]) for i in range(2)]
    comp_r = [_pearson(coords[rows, i], completeness[rows]) for i in range(2)]
    loadings = []
    for i in range(2):
        idx = np.argsort(vt[i])
        def item(j):
            f = int(cols[j])
            return {"id": fs.ids[f], "label": fs.labels[f], "name": fs.names[f], "weight": float(vt[i, j])}
        loadings.append({"high": [item(j) for j in idx[::-1][:TOP_LOADINGS]],
                         "low": [item(j) for j in idx[:TOP_LOADINGS]]})
    return Embedding(fs.kind, coords, explained, sizes, int(q), size_r, comp_r, loadings,
                     time.perf_counter() - t0, full)


def neighbours(emb: Embedding, catalog, genome_index: int, top: int = 10) -> list[dict[str, Any]]:
    """Genomes ranked by Jaccard similarity of feature sets to one genome."""
    x = emb.matrix[genome_index]
    own = float(emb.n_features[genome_index])
    if own == 0:
        return []
    shared = np.asarray((emb.matrix @ x.T).todense()).ravel()
    union = own + emb.n_features - shared
    with np.errstate(divide="ignore", invalid="ignore"):
        jac = np.where(union > 0, shared / union, 0.0)
    jac[genome_index] = -1.0
    me = catalog.genomes[genome_index]
    out = []
    for j in np.argsort(-jac)[:top]:
        if jac[j] <= 0:
            break
        g = catalog.genomes[int(j)]
        differs = next((r for r in ("class", "order") if (me.taxonomy.get(r) or UNCLASSIFIED) != UNCLASSIFIED
                        and (g.taxonomy.get(r) or UNCLASSIFIED) != UNCLASSIFIED
                        and me.taxonomy.get(r) != g.taxonomy.get(r)), None)
        names = [g.taxonomy.get(r) for r in ("class", "order", "family")]
        names = [n for n in names if n and n != UNCLASSIFIED]
        out.append({"genome": g, "similarity": float(jac[j]), "shared": int(shared[j]),
                    "features": int(emb.n_features[j]), "differs_at": differs,
                    "lineage": "; ".join(n for k, n in enumerate(names) if k == 0 or n != names[k - 1])})
    return out


class Landscapes:
    """Embeddings per feature kind, computed once in a background thread and kept."""

    def __init__(self, ctx, background: bool = True) -> None:
        self.ctx, self.background = ctx, background
        self._done: dict[str, Embedding] = {}
        self._errors: dict[str, str] = {}
        self._running: set[str] = set()
        self._lock = threading.Lock()

    def _compute(self, kind: str) -> None:
        try:
            fs = self.ctx.feature_sets.get(kind)
            note = ""
            if kind == "ko" and not len(fs.ids):
                fs, note = self.ctx.feature_sets.get("pfam"), "This dataset has no KEGG orthologs; showing Pfam families."
            emb = embed(fs, self.ctx.catalog)
            emb.note = emb.note or note
            with self._lock:
                self._done[kind] = emb
        except Exception as exc:  # surfaced on the page
            with self._lock:
                self._errors[kind] = str(exc) or exc.__class__.__name__
        finally:
            with self._lock:
                self._running.discard(kind)

    def get(self, kind: str) -> Embedding | None:
        """The embedding, or None while it is computed. Raises LookupError on failure."""
        with self._lock:
            if kind in self._done:
                return self._done[kind]
            if kind in self._errors:
                raise LookupError(self._errors[kind])
            start = kind not in self._running
            self._running.add(kind)
        if start:
            if self.background:
                threading.Thread(target=self._compute, args=(kind,), daemon=True).start()
            else:
                self._compute(kind)
                return self.get(kind)
        return None


def payload(emb: Embedding, catalog) -> dict[str, Any]:
    """Points, taxonomy as indexed name tables, and the axis summary for the page."""
    from sharur.browser.charts import PALETTE  # noqa: PLC0415

    keep = np.flatnonzero(np.isfinite(emb.coords[:, 0]))
    genomes = [catalog.genomes[i] for i in keep]
    ranks = {}
    for rank in RANKS[1:]:
        names: dict[str, int] = {}
        idx = [names.setdefault(g.taxonomy.get(rank) or UNCLASSIFIED, len(names)) for g in genomes]
        ranks[rank] = {"names": list(names), "idx": idx}
    default_rank = next((r for r in ("order", "class", "family", "phylum")
                         if 4 <= len(set(ranks[r]["names"]) - {UNCLASSIFIED}) <= 40), "order")
    return {
        "kind": emb.kind, "label": IN_TEXT.get(emb.kind, emb.kind), "note": emb.note,
        "feature_url": FEATURE_URLS.get(emb.kind, "/search?q={id}"),
        "ids": [g.bin_id for g in genomes],
        "x": [round(float(v), 4) for v in emb.coords[keep, 0]],
        "y": [round(float(v), 4) for v in emb.coords[keep, 1]],
        "completeness": [None if g.completeness is None else round(float(g.completeness), 1) for g in genomes],
        "features": [int(emb.n_features[i]) for i in keep],
        "ranks": ranks, "default_rank": default_rank, "palette": PALETTE,
        "explained": emb.explained, "size_r": emb.axis_size_r, "completeness_r": emb.axis_completeness_r,
        "loadings": emb.loadings, "used_features": emb.used_features, "seconds": round(emb.seconds, 2),
        "left_out": int(len(catalog.genomes) - len(keep)), "inline_selection": INLINE_SELECTION,
    }


def register(app: FastAPI, ctx: SimpleNamespace) -> None:
    """Add /landscape, /api/landscape and /api/neighbors. Needs ``ctx.feature_sets`` (routes_matrix)."""
    catalog = ctx.catalog
    maps = Landscapes(ctx, background=getattr(ctx, "background", True))
    app.state.landscapes = maps

    @app.get("/landscape", response_class=HTMLResponse)
    def landscape_page(request: Request, kind: str = Query("ko"), color: str = Query(""),
                       q: str = Query("", max_length=200)):
        kind = kind if kind in KINDS else "ko"
        maps.get(kind) if not maps.background else None  # tests: compute up front
        return ctx.render(request, "landscape.html", "taxa", kind=kind, color=color if color in RANKS else "",
                          q=q, KINDS=KINDS, in_text=IN_TEXT[kind])

    @app.get("/api/landscape")
    def landscape_data(kind: str = Query("ko")):
        if kind not in KINDS:
            raise HTTPException(404, "Unknown feature kind")
        try:
            emb = maps.get(kind)
        except LookupError as exc:
            return JSONResponse({"status": "error", "message": str(exc)}, status_code=200)
        if emb is None:
            return JSONResponse({"status": "computing"}, status_code=202)
        return {"status": "ready", **payload(emb, catalog)}

    @app.post("/api/selection")
    async def store_selection(request: Request):
        """Keep a set of genomes under a short key usable wherever a genome or clade is accepted."""
        import hashlib  # noqa: PLC0415

        body = await request.json()
        ids = sorted({str(i) for i in (body.get("genomes") or []) if str(i) in catalog.by_bin})
        if not ids:
            raise HTTPException(400, "No known genomes in the selection")
        key = "selection:" + hashlib.sha1("\n".join(ids).encode()).hexdigest()[:12]
        catalog.selections[key] = ids
        return {"token": key, "genomes": len(ids)}

    @app.get("/api/neighbors/{bin_id:path}", response_class=HTMLResponse)
    def genome_neighbours(request: Request, bin_id: str, kind: str = Query("ko")):
        genome = catalog.by_bin.get(bin_id)
        if genome is None:
            raise HTTPException(404, "Genome not found")
        kind = kind if kind in KINDS else "ko"
        try:
            emb = maps.get(kind)
        except LookupError as exc:
            return HTMLResponse(f'<p class="empty">Similarity unavailable: {exc}</p>')
        if emb is None:
            return HTMLResponse('<p class="empty" data-retry>Comparing feature sets across genomes…</p>',
                                status_code=202)
        rows = neighbours(emb, catalog, genome.index)
        return ctx.render(request, "_neighbors.html", "taxa", g=genome, rows=rows, kind=emb.kind,
                          kind_label=IN_TEXT.get(emb.kind, emb.kind), own=int(emb.n_features[genome.index]))
