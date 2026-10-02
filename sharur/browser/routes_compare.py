"""Comparison views for the browser: genomes or clades side by side, and a function × clade heatmap.

Registered by :func:`sharur.browser.app.create_app` through :func:`register`.
Both views read the catalog's genome × predicate arrays; Pfam family presence
per genome is loaded once (in the background) into the same compact form.
"""

from __future__ import annotations

import threading
from collections import Counter
from dataclasses import dataclass
from types import SimpleNamespace
from typing import Any

import numpy as np
from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.responses import HTMLResponse

from sharur.browser.catalog import BIOLOGICAL, RANKS, UNCLASSIFIED, Genome
from sharur.predicates.vocabulary import PREDICATE_BY_ID

MIN_CLADE = 3


@dataclass
class Side:
    """One side of a comparison: a single genome or a clade."""

    key: str            # the URL token: a genome id or ``rank:name``
    label: str
    kind: str           # "genome" or the rank
    genomes: list[Genome]
    url: str

    @property
    def is_genome(self) -> bool:
        return self.kind == "genome"

    @property
    def filter_param(self) -> str:
        """Query string restricting a function page to this side."""
        from urllib.parse import quote  # noqa: PLC0415

        return f"genome={quote(self.key, safe='')}" if self.is_genome else f"clade={quote(self.key, safe='')}"


GENOME_SET = "genomes:"
SELECTION = "selection:"


def genome_set(catalog, token: str | None) -> list[Genome] | None:
    """Genomes named by ``genomes:ID1,ID2,...`` or by a stored ``selection:<key>`` (see
    ``routes_landscape``); unknown IDs drop out. None for any other token."""
    if not token:
        return None
    if token.startswith(SELECTION):
        ids = getattr(catalog, "selections", {}).get(token, [])
    elif token.startswith(GENOME_SET):
        ids = [t for t in dict.fromkeys(x.strip() for x in token[len(GENOME_SET):].split(",")) if t]
    else:
        return None
    return [catalog.by_bin[i] for i in ids if i in catalog.by_bin]


def resolve_side(catalog, token: str, url) -> Side | None:
    """``token`` is a genome id, ``rank:name`` (e.g. ``class:Omnitrophia``) or ``genomes:ID1,ID2,...``."""
    token = (token or "").strip()
    if not token:
        return None
    genome = catalog.by_bin.get(token)
    if genome is not None:
        return Side(token, token, "genome", [genome], url("genome", token))
    selected = genome_set(catalog, token)
    if selected is not None:
        if not selected:
            return None
        if len(selected) == 1:
            g = selected[0]
            return Side(g.bin_id, g.bin_id, "genome", selected, url("genome", g.bin_id))
        return Side(token, f"{len(selected)} selected genomes", "selection", selected, "/landscape")
    if ":" in token:
        rank, name = token.split(":", 1)
        if rank in RANKS:
            members = catalog.clade(rank, name)
            if members:
                return Side(token, name, rank, members, url("taxa", rank, name))
    return None


def clade_filter(catalog, value: str | None) -> tuple[str, set[str]] | None:
    """Parse a ``clade=rank:name`` or genome id filter into (label, bin ids)."""
    if not value:
        return None
    if value in catalog.by_bin:
        return value, {value}
    selected = genome_set(catalog, value)
    if selected:
        return f"{len(selected)} selected genomes", {g.bin_id for g in selected}
    if ":" in value:
        rank, name = value.split(":", 1)
        if rank in RANKS:
            members = catalog.clade(rank, name)
            if members:
                return f"{name} ({rank})", {g.bin_id for g in members}
    return None


def descendants(root: str) -> list[str]:
    """``root`` and every vocabulary predicate below it (is-a), breadth first."""
    out, frontier = [], [root]
    while frontier:
        current = frontier.pop(0)
        if current in PREDICATE_BY_ID and current not in out:
            out.append(current)
            frontier.extend(p.predicate_id for p in PREDICATE_BY_ID.values() if p.parent == current)
    return out


def presets(catalog) -> list[tuple[str, list[str]]]:
    """Ready-made function sets for the heatmap, limited to labels present in the dataset."""
    groups = [
        ("Hydrogenases", descendants("hydrogenase")),
        ("Terminal oxidases", descendants("terminal_oxidase")),
        ("Nitrogen metabolism", descendants("nitrogen_metabolism")),
        ("Carbon fixation", descendants("carbon_fixation")),
        ("CRISPR-Cas types", sorted(p for p in PREDICATE_BY_ID if p.startswith("crispr_type"))),
        ("Secretion systems", descendants("secretion_system")),
    ]
    present = set(catalog.predicate_index)
    out = []
    for name, ids in groups:
        ids = [p for p in ids if p in present]
        if len(ids) >= 2:
            out.append((name, ids))
    return out


# --------------------------------------------------------------------------- #
# Pfam presence per genome (compact pairs, like the predicate pairs)
# --------------------------------------------------------------------------- #


class PfamPairs:
    def __init__(self) -> None:
        self.ready = threading.Event()
        self.accessions: list[str] = []
        self.names: dict[str, str] = {}
        self.bin_idx: np.ndarray | None = None
        self.acc_idx: np.ndarray | None = None

    def load(self, store, catalog, lock) -> None:
        # Families come back as integer codes per genome (one row per genome), not one row per pair:
        # millions of (genome, accession) tuples cost far more memory than the arrays they become.
        try:
            with lock:
                families = store.execute("""
                    SELECT split_part(accession, '.', 1) AS acc, MIN(COALESCE(NULLIF(name, ''), accession))
                    FROM annotations WHERE LOWER(source) = 'pfam' GROUP BY 1 ORDER BY 1""")
                per_genome = store.execute("""
                    WITH codes AS (SELECT acc, (ROW_NUMBER() OVER (ORDER BY acc) - 1)::INTEGER AS k FROM (
                        SELECT DISTINCT split_part(accession, '.', 1) AS acc FROM annotations WHERE LOWER(source) = 'pfam'))
                    SELECT p.bin_id, LIST(DISTINCT c.k) FROM annotations a JOIN proteins p USING (protein_id)
                    JOIN codes c ON c.acc = split_part(a.accession, '.', 1)
                    WHERE LOWER(a.source) = 'pfam' GROUP BY 1""")
            self.names = dict(families)
            bins, accs = [], []
            for bin_id, codes in per_genome:
                genome = catalog.by_bin.get(bin_id)
                if genome is not None:
                    bins.append(np.full(len(codes), genome.index, dtype=np.int32))
                    accs.append(np.asarray(codes, dtype=np.int32))
            self.accessions = [acc for acc, _ in families]
            self.bin_idx = np.concatenate(bins) if bins else np.zeros(0, dtype=np.int32)
            self.acc_idx = np.concatenate(accs) if accs else np.zeros(0, dtype=np.int32)
        finally:
            self.ready.set()

    def shares(self, genome_indices: list[int], n_genomes: int) -> np.ndarray:
        inside = np.zeros(n_genomes, dtype=bool)
        inside[genome_indices] = True
        hits = np.bincount(self.acc_idx[inside[self.bin_idx]], minlength=len(self.accessions))
        return hits / max(len(genome_indices), 1)


# --------------------------------------------------------------------------- #
# Comparison
# --------------------------------------------------------------------------- #


def _predicate_shares(catalog, side: Side) -> tuple[np.ndarray, np.ndarray]:
    """Per predicate: share of the side's genomes carrying it, and total proteins."""
    n = len(catalog.genomes)
    inside = np.zeros(n, dtype=bool)
    inside[[g.index for g in side.genomes]] = True
    mask = inside[catalog.pair_bin]
    genomes = np.bincount(catalog.pair_pred[mask], minlength=len(catalog.predicates))
    proteins = np.bincount(catalog.pair_pred[mask], weights=catalog.pair_count[mask], minlength=len(catalog.predicates))
    return genomes / len(side.genomes), proteins


def _differences(names: list[str], a: np.ndarray, b: np.ndarray, *, top: int, threshold: float,
                 keep=lambda name: True) -> tuple[list[tuple[str, float, float]], list[tuple[str, float, float]]]:
    """Items most enriched on each side (difference in share ≥ ``threshold``)."""
    diff = a - b
    more_a = [(names[k], float(a[k]), float(b[k])) for k in np.argsort(-diff)
              if diff[k] >= threshold and keep(names[k])][:top]
    more_b = [(names[k], float(a[k]), float(b[k])) for k in np.argsort(diff)
              if -diff[k] >= threshold and keep(names[k])][:top]
    return more_a, more_b


def _side_stats(side: Side) -> dict[str, Any]:
    genomes = side.genomes
    sizes = [g.length for g in genomes if g.length]
    return {
        "genomes": len(genomes),
        "size": float(np.median(sizes)) if sizes else 0.0,
        "proteins": float(np.median([g.proteins for g in genomes])),
        "contigs": float(np.median([g.contigs for g in genomes])),
        "n50": float(np.median([g.n50 for g in genomes])),
        "annotated": float(np.median([g.annotated / (g.proteins or 1) for g in genomes])),
    }


def _system_shares(catalog, side: Side) -> dict[tuple[str, str], float]:
    members = {g.bin_id for g in side.genomes}
    carriers: dict[tuple[str, str], set[str]] = {}
    for s in catalog.systems:
        if s["bin_id"] in members:
            carriers.setdefault((s["kind"], s["type"]), set()).add(s["bin_id"])
    return {k: len(v) / len(members) for k, v in carriers.items()}


def compare(catalog, pfam: PfamPairs, a: Side, b: Side) -> dict[str, Any]:
    n = len(catalog.genomes)
    single = a.is_genome and b.is_genome
    threshold = 0.5 if single else 0.25
    out: dict[str, Any] = {"stats": (_side_stats(a), _side_stats(b))}

    if catalog.pair_bin is not None:
        share_a, prot_a = _predicate_shares(catalog, a)
        share_b, prot_b = _predicate_shares(catalog, b)
        biological = {p for p in catalog.predicates if p in PREDICATE_BY_ID
                      and PREDICATE_BY_ID[p].category in BIOLOGICAL}
        more_a, more_b = _differences(catalog.predicates, share_a, share_b, top=30, threshold=threshold,
                                      keep=lambda p: p in biological)
        k = catalog.predicate_index
        out["functions"] = tuple([{"id": p, "name": PREDICATE_BY_ID[p].name, "category": PREDICATE_BY_ID[p].category,
                                   "a": sa, "b": sb, "copies_a": int(prot_a[k[p]]), "copies_b": int(prot_b[k[p]])}
                                  for p, sa, sb in rows] for rows in (more_a, more_b))
        shared = int(((share_a > 0) & (share_b > 0)).sum())
        out["function_counts"] = {"a": int((share_a > 0).sum()), "b": int((share_b > 0).sum()), "shared": shared}

    if catalog.module_completeness is not None:
        mean_a = catalog.module_completeness[[g.index for g in a.genomes]].mean(axis=0)
        mean_b = catalog.module_completeness[[g.index for g in b.genomes]].mean(axis=0)
        more_a, more_b = _differences(catalog.module_ids, mean_a, mean_b, top=20, threshold=0.25)
        out["modules"] = tuple([{"module": m, "name": catalog.modules[m].name, "a": sa, "b": sb} for m, sa, sb in rows]
                               for rows in (more_a, more_b))

    sys_a, sys_b = _system_shares(catalog, a), _system_shares(catalog, b)
    keys = set(sys_a) | set(sys_b)
    systems = sorted(({"kind": k, "type": t, "a": sys_a.get((k, t), 0.0), "b": sys_b.get((k, t), 0.0)}
                      for k, t in keys), key=lambda r: -abs(r["a"] - r["b"]))
    out["systems"] = [r for r in systems if abs(r["a"] - r["b"]) >= (0.5 if single else 0.15)][:24]

    if pfam.ready.is_set() and pfam.bin_idx is not None:
        pa = pfam.shares([g.index for g in a.genomes], n)
        pb = pfam.shares([g.index for g in b.genomes], n)
        more_a, more_b = _differences(pfam.accessions, pa, pb, top=30, threshold=threshold)
        out["domains"] = tuple([{"accession": acc, "name": pfam.names.get(acc, acc), "a": sa, "b": sb}
                                for acc, sa, sb in rows] for rows in (more_a, more_b))
        out["domain_counts"] = {"a": int((pa > 0).sum()), "b": int((pb > 0).sum()),
                                "shared": int(((pa > 0) & (pb > 0)).sum())}
    return out


# --------------------------------------------------------------------------- #
# Heatmap
# --------------------------------------------------------------------------- #


def heatmap(catalog, functions: list[str], rank: str, include_small: bool) -> dict[str, Any]:
    clades = Counter(g.taxonomy.get(rank) or UNCLASSIFIED for g in catalog.genomes)
    names = [name for name, n in clades.most_common() if include_small or n >= MIN_CLADE]
    column = {name: j for j, name in enumerate(names)}
    clade_of = np.array([column.get(g.taxonomy.get(rank) or UNCLASSIFIED, -1) for g in catalog.genomes],
                        dtype=np.int32)
    sizes = np.array([clades[name] for name in names], dtype=float)
    rows = []
    for p in functions:
        k = catalog.predicate_index.get(p)
        if k is None or catalog.pair_pred is None:
            counts = np.zeros(len(names))
        else:
            bins = catalog.pair_bin[catalog.pair_pred == k]
            cols = clade_of[bins]
            counts = np.bincount(cols[cols >= 0], minlength=len(names)).astype(float)
        d = PREDICATE_BY_ID.get(p)
        rows.append({"id": p, "name": d.name if d else p, "category": d.category if d else "",
                     "shares": (counts / sizes).tolist() if len(names) else [], "carriers": int(counts.sum())})
    hidden = sum(1 for n in clades.values() if n < MIN_CLADE) if not include_small else 0
    return {"clades": [(name, int(clades[name])) for name in names], "rows": rows, "hidden": hidden}


# --------------------------------------------------------------------------- #
# Routes
# --------------------------------------------------------------------------- #


def register(app: FastAPI, ctx: SimpleNamespace) -> None:
    """Add /compare and /heatmap. ``ctx``: store, lock, catalog, render, url."""
    catalog, render, url = ctx.catalog, ctx.render, ctx.url
    pfam = PfamPairs()
    app.state.pfam_pairs = pfam
    summaries = getattr(ctx, "summaries", ctx.store)
    if ctx.background:
        threading.Thread(target=pfam.load, args=(summaries, catalog, ctx.lock), daemon=True).start()
    else:
        pfam.load(summaries, catalog, ctx.lock)

    @app.get("/compare", response_class=HTMLResponse)
    def compare_page(request: Request, a: str = Query("", max_length=500), b: str = Query("", max_length=500)):
        side_a, side_b = resolve_side(catalog, a, url), resolve_side(catalog, b, url)
        if (a and side_a is None) or (b and side_b is None):
            raise HTTPException(404, "Unknown genome or clade; use a genome id or rank:name")
        result = compare(catalog, pfam, side_a, side_b) if side_a and side_b else None
        return render(request, "compare.html", "taxa", a=side_a, b=side_b, a_key=a, b_key=b, result=result,
                      pfam_ready=pfam.ready.is_set())

    @app.get("/heatmap", response_class=HTMLResponse)
    def heatmap_page(request: Request, functions: str = Query("", max_length=4000), rank: str = Query("order"),
                     small: int = Query(0), preset: str = Query("", max_length=100)):
        rank = rank if rank in RANKS[1:] else "order"
        options = presets(catalog)
        chosen = [f for f in dict.fromkeys(x.strip() for x in functions.split(",")) if f]
        if preset and not chosen:
            chosen = next((ids for name, ids in options if name == preset), [])
        if not chosen and options:
            chosen = options[0][1]
        chosen = [f for f in chosen if f in PREDICATE_BY_ID or f in catalog.predicate_index][:60]
        data = heatmap(catalog, chosen, rank, bool(small))
        return render(request, "heatmap.html", "functions", data=data, chosen=chosen, rank=rank, small=bool(small),
                      presets=options, ranks=RANKS[1:])
