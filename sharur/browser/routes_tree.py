"""Taxonomy tree with features mapped onto it.

``/tree`` draws the dataset's GTDB taxonomy as a cladogram (radial or
rectangular) down to genus. Clades open and close in the page; a clade's
name zooms the tree to it (``/tree?root=rank:name``). Up to eight features,
each a KO, Pfam family, function label, KEGG module (completeness ≥ 0.75) or
curated system type, are drawn as rings (radial) or columns (rectangular)
holding the share of each clade's genomes that carry the feature.

Features are written ``kind:id``: ``ko:K00001``, ``pfam:PF00005``,
``function:hydrogenase``, ``module:M00001``, ``system:defense:CBASS``,
``system:crispr:I-B``. Prevalence comes from the presence/absence matrix's
feature tables (``routes_matrix.FeatureSets``).
"""

from __future__ import annotations

import re
from typing import Any

import numpy as np
from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.responses import HTMLResponse, JSONResponse

from sharur.browser.catalog import RANKS, UNCLASSIFIED, Genome
from sharur.browser.charts import PALETTE
from sharur.browser.routes_matrix import KINDS, FeatureSets
from sharur.predicates.vocabulary import PREDICATE_BY_ID

MAX_FEATURES = 8
_EC = re.compile(r"\s*\[EC:[^]]*\]")
FEATURE_KINDS = ("ko", "pfam", "function", "module", "system")


# --------------------------------------------------------------------------- #
# Features
# --------------------------------------------------------------------------- #


def parse_features(text: str) -> list[tuple[str, str]]:
    """``kind:id`` tokens, comma separated, in order and without repeats."""
    out: list[tuple[str, str]] = []
    for token in (t.strip() for t in (text or "").split(",")):
        if ":" not in token:
            continue
        kind, ident = token.split(":", 1)
        kind = kind.strip().lower()
        if kind in FEATURE_KINDS and ident.strip() and (kind, ident.strip()) not in out:
            out.append((kind, ident.strip()))
    return out[:MAX_FEATURES]


def feature_mask(fs, ident: str, n_genomes: int) -> np.ndarray | None:
    """Genomes carrying feature ``ident`` (present at the set's threshold)."""
    k = fs.index.get(ident)
    if k is None:
        return None
    keep = (fs.feat_idx == k) & (fs.value >= fs.threshold)
    mask = np.zeros(n_genomes, dtype=bool)
    mask[fs.bin_idx[keep]] = True
    return mask


# --------------------------------------------------------------------------- #
# Tree
# --------------------------------------------------------------------------- #


def build_tree(genomes: list[Genome], root_rank: str | None, root_name: str | None,
               masks: list[np.ndarray]) -> dict[str, Any]:
    """Nested clades below the root, down to genus, with genome and carrier counts.

    A genome unclassified at a rank stops there, under an ``Unclassified``
    clade of that rank.
    """
    nfeat = len(masks)
    start = RANKS.index(root_rank) + 1 if root_rank else 0

    def node(name: str, rank: str) -> dict[str, Any]:
        return {"name": name, "rank": rank, "n": 0, "k": [0] * nfeat, "children": {}}

    root = node(root_name or "All genomes", root_rank or "root")
    for g in genomes:
        hits = [int(m[g.index]) for m in masks]
        path = [root]
        current = root
        for rank in RANKS[start:]:
            name = g.taxonomy.get(rank) or UNCLASSIFIED
            child = current["children"].get(name)
            if child is None:
                child = current["children"][name] = node(name, rank)
            path.append(child)
            if name == UNCLASSIFIED:
                break
            current = child
        for n in path:
            n["n"] += 1
            for i, h in enumerate(hits):
                n["k"][i] += h

    def finish(n: dict[str, Any]) -> dict[str, Any]:
        children = sorted(n["children"].values(), key=lambda c: (c["name"] == UNCLASSIFIED, -c["n"], c["name"]))
        n["children"] = [finish(c) for c in children]
        return n

    return finish(root)


def effective_root(tree: dict[str, Any]) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Descend through single-child clades; returns the root shown and the clades passed."""
    passed = []
    while len(tree["children"]) == 1 and tree["children"][0]["children"]:
        passed.append(tree)
        tree = tree["children"][0]
    return tree, passed


def depth(tree: dict[str, Any]) -> int:
    return 1 + max((depth(c) for c in tree["children"]), default=0)


def count_nodes(tree: dict[str, Any]) -> int:
    return 1 + sum(count_nodes(c) for c in tree["children"])


# --------------------------------------------------------------------------- #
# Feature search
# --------------------------------------------------------------------------- #


def feature_options(ctx) -> list[dict[str, Any]]:
    """Every mappable feature in the dataset: token, kind, id, label, genomes (when known)."""
    cached = getattr(ctx, "_tree_feature_options", None)
    if cached is not None and ctx.catalog.ready.is_set():
        return cached
    catalog = ctx.catalog
    out: list[dict[str, Any]] = []
    for p in catalog.predicates:
        d = PREDICATE_BY_ID.get(p)
        k = catalog.predicate_index.get(p)
        out.append({"token": f"function:{p}", "kind": "function", "id": p, "label": d.name if d else p,
                    "genomes": int(catalog.predicate_genomes[k]) if catalog.predicate_genomes is not None else None})
    for acc, info in catalog.domains.items():
        out.append({"token": f"pfam:{acc}", "kind": "pfam", "id": acc, "label": info["name"],
                    "detail": info.get("description", ""), "genomes": info.get("genomes")})
    for m, definition in catalog.modules.items():
        out.append({"token": f"module:{m}", "kind": "module", "id": m, "label": definition.name, "genomes": None})
    kos = set()
    for present in catalog.ko_sets.values():
        kos |= present
    names = getattr(getattr(ctx, "ko_names", None), "_names", None) or {}
    for ko in sorted(kos):
        symbols, definition = names.get(ko, ("", ""))
        definition = _EC.sub("", definition)
        label = f"{symbols.split(',')[0].strip()} · {definition}".strip(" ·") or ko
        out.append({"token": f"ko:{ko}", "kind": "ko", "id": ko, "label": label, "genomes": None})
    systems: dict[str, set[str]] = {}
    for s in catalog.systems:
        systems.setdefault(f"{s['kind']}:{s['type']}", set()).add(s["bin_id"])
    cctyper = getattr(ctx, "cctyper_systems", None)
    for s in (cctyper() if cctyper else []):
        if s.get("confident"):
            systems.setdefault(f"crispr:{s['prediction']}", set()).add(s["bin_id"])
    for key, bins in systems.items():
        family, name = key.split(":", 1)
        label = f"CRISPR-Cas {name}" if family == "crispr" else name
        out.append({"token": f"system:{key}", "kind": "system", "id": key, "label": label, "detail": f"{family} system",
                    "genomes": len(bins)})
    if catalog.ready.is_set():
        ctx._tree_feature_options = out
    return out


def search_features(options: list[dict[str, Any]], q: str, limit: int = 15) -> list[dict[str, Any]]:
    q = q.strip().lower()
    if len(q) < 2:
        return []
    scored = []
    for o in options:
        ident, label = o["id"].lower(), o["label"].lower()
        if q == ident or q == label:
            score = 0
        elif ident.startswith(q) or label.startswith(q):
            score = 1
        elif q in label or q in ident or q in (o.get("detail") or "").lower():
            score = 2
        else:
            continue
        scored.append((score, -(o.get("genomes") or 0), o["label"], o))
    scored.sort(key=lambda t: t[:3])
    return [o for *_, o in scored[:limit]]


def presets(ctx) -> list[tuple[str, list[str]]]:
    """Ready-made feature sets present in the dataset."""
    from sharur.browser.routes_compare import presets as function_presets  # noqa: PLC0415

    out = [(name, [f"function:{p}" for p in ids[:MAX_FEATURES]]) for name, ids in function_presets(ctx.catalog)]
    counts: dict[str, set[str]] = {}
    for s in ctx.catalog.systems:
        if s["kind"] == "defense":
            counts.setdefault(s["type"], set()).add(s["bin_id"])
    top = sorted(counts, key=lambda t: -len(counts[t]))[:6]
    if len(top) >= 2:
        out.append(("Most common defense systems", [f"system:defense:{t}" for t in top]))
    cctyper = getattr(ctx, "cctyper_systems", None)
    subtypes: dict[str, set[str]] = {}
    for s in (cctyper() if cctyper else []):
        if s.get("confident") and not s["prediction"].startswith("Hybrid"):
            subtypes.setdefault(s["prediction"], set()).add(s["bin_id"])
    top = sorted(subtypes, key=lambda t: -len(subtypes[t]))[:6]
    if len(top) >= 2:
        out.append(("CRISPR-Cas subtypes", [f"system:crispr:{t}" for t in top]))
    return out


# --------------------------------------------------------------------------- #
# Routes
# --------------------------------------------------------------------------- #


def register(app: FastAPI, ctx) -> None:
    """Add /tree and /api/tree/features. ``ctx``: catalog, render, url, and the matrix's ``feature_sets``."""
    catalog = ctx.catalog

    def sets() -> FeatureSets:
        fs = getattr(ctx, "feature_sets", None)
        if fs is None:
            fs = ctx.feature_sets = FeatureSets(ctx)
        return fs

    masks_cache: dict[tuple[str, str], np.ndarray | None] = {}

    def resolve(kind: str, ident: str) -> tuple[dict[str, Any], np.ndarray | None, str]:
        """Feature description, genome mask, and a message when it can't be mapped."""
        try:
            fs = sets().get(kind)
        except LookupError as exc:
            return {"kind": kind, "id": ident, "label": ident}, None, str(exc)
        k = fs.index.get(ident)
        info = {"kind": kind, "id": ident, "token": f"{kind}:{ident}",
                "label": fs.labels[k] if k is not None else ident, "name": fs.names[k] if k is not None else ident,
                "kind_label": KINDS[kind]}
        if k is None:
            return info, None, f"{ident} is not in this dataset"
        key = (kind, ident)
        if key not in masks_cache or (kind in ("module", "pfam") and not catalog.ready.is_set()):
            masks_cache[key] = feature_mask(fs, ident, len(catalog.genomes))
        return info, masks_cache[key], ""

    @app.get("/api/tree/features")
    def tree_features(q: str = Query("", max_length=100)):
        return JSONResponse(search_features(feature_options(ctx), q))

    @app.get("/tree", response_class=HTMLResponse)
    def tree_page(request: Request, root: str = Query("", max_length=300), features: str = Query("", max_length=2000),
                  layout: str = Query("radial"), open: str = Query("", max_length=20)):  # noqa: A002
        root_rank = root_name = None
        if root:
            if ":" not in root:
                raise HTTPException(404, "Use rank:name for the root, e.g. order:Woesearchaeales")
            root_rank, root_name = root.split(":", 1)
            if root_rank not in RANKS:
                raise HTTPException(404, "Unknown rank")
        genomes = catalog.clade(root_rank, root_name)
        if not genomes:
            raise HTTPException(404, "No genomes in that clade")
        chosen, masks, notes = [], [], []
        for kind, ident in parse_features(features):
            info, mask, message = resolve(kind, ident)
            if mask is None:
                notes.append(message)
                continue
            chosen.append(info)
            masks.append(mask)
        tree = build_tree(genomes, root_rank, root_name, masks)
        shown, passed = effective_root(tree)
        lineage = catalog.lineage_of(root_rank, root_name) if root_rank else []
        for i, f in enumerate(chosen):
            f["color"] = PALETTE[i % len(PALETTE)]
            f["carriers"] = int(tree["k"][i])
        payload = {"tree": shown, "features": chosen, "layout": layout if layout in ("radial", "rect") else "radial",
                   "open": open if open in RANKS else "", "depth": depth(shown) - 1, "root": root,
                   "ranks": list(RANKS)}
        return ctx.render(request, "tree.html", "taxa", payload=payload, tree=shown, passed=passed, lineage=lineage,
                          root=root, root_rank=root_rank, root_name=root_name, chosen=chosen, notes=notes,
                          features=",".join(f["token"] for f in chosen), presets=presets(ctx),
                          layout=payload["layout"], n_genomes=len(genomes), nodes=count_nodes(shown),
                          max_features=MAX_FEATURES,
                          open_ranks=[r for r in RANKS if r != "domain" and
                                      (shown["rank"] == "root" or RANKS.index(r) > RANKS.index(shown["rank"]))])
