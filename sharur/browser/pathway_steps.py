"""Where a KEGG module's steps are found, and where they go missing.

For one module, every genome's step-by-step result (``sharur.modules.evaluate``)
is computed once and cached. From it the pathway page draws:

- a step diagram: the module definition left to right, alternatives stacked
  within a step, complex components joined with ``+`` (optional ones dashed),
  nested modules as links; each KO box shaded by the share of genomes in scope
  carrying that KO;
- a clade × step heatmap: the share of each clade's genomes completing each
  step, with each clade's usual gap outlined;
- genomes one step short of complete, grouped by the step they lack.

A step can look missing because its gene was not annotated (a divergent
homolog below the KOfam threshold), sits in an unassembled part of the genome,
or is absent; the page says so where it reports a gap.
"""

from __future__ import annotations

import threading
from collections import Counter, OrderedDict
from dataclasses import dataclass
from html import escape
from typing import Any, Callable
from urllib.parse import quote

import numpy as np
from markupsafe import Markup

from sharur.browser.catalog import RANKS, UNCLASSIFIED
from sharur.modules import ModuleDefinition, Node, evaluate, kos_in

KO_W, KO_H = 66, 32
GAP, PLUS_GAP, ARROW_GAP = 4, 14, 22
HEADER_H = 34
SHADE_MIN, SHADE_SPAN = 6, 64


# --------------------------------------------------------------------------- #
# Per-genome step results, cached per module
# --------------------------------------------------------------------------- #


@dataclass
class StepTable:
    module: str
    steps: list[int]                 # step numbers as evaluate() reports them (1-based, undefined steps skipped)
    complete: np.ndarray             # genomes × steps, bool
    missing: dict[tuple[int, int], list[str]]   # (genome index, step column) -> KOs that would complete it
    relevant: np.ndarray             # genomes with any KO of the module (or of modules it refers to)

    @property
    def n_steps(self) -> int:
        return len(self.steps)


def _all_kos(definition: ModuleDefinition, modules: dict[str, ModuleDefinition], depth: int = 0) -> set[str]:
    kos = set(kos_in(definition.tree))
    if depth < 5:
        for ref in _module_refs(definition.tree):
            if ref in modules:
                kos |= _all_kos(modules[ref], modules, depth + 1)
    return kos


def _module_refs(node: Node) -> set[str]:
    if node.kind == "module":
        return {node.value}
    return {m for child in node.children for m in _module_refs(child)}


def step_table(catalog, definition: ModuleDefinition) -> StepTable:
    modules = catalog.modules
    empty = evaluate(definition, set(), modules)
    steps = [s.index for s in empty.steps]
    col = {s: j for j, s in enumerate(steps)}
    n = len(catalog.genomes)
    complete = np.zeros((n, len(steps)), dtype=bool)
    relevant = np.zeros(n, dtype=bool)
    missing: dict[tuple[int, int], list[str]] = {}
    wanted = _all_kos(definition, modules)
    for g in catalog.genomes:
        present = catalog.ko_sets.get(g.bin_id)
        if not present or not (present & wanted):
            continue
        relevant[g.index] = True
        result = evaluate(definition, present, modules)
        for s in result.steps:
            j = col.get(s.index)
            if j is None:
                continue
            if s.score >= 1.0:
                complete[g.index, j] = True
            else:
                missing[(g.index, j)] = s.missing
    return StepTable(definition.module, steps, complete, missing, relevant)


class StepCache:
    """Step tables per module, recomputed when the catalog's module data change."""

    def __init__(self, size: int = 48) -> None:
        self._tables: OrderedDict[str, StepTable] = OrderedDict()
        self._lock = threading.Lock()
        self._size = size

    def get(self, catalog, definition: ModuleDefinition) -> StepTable:
        with self._lock:
            table = self._tables.get(definition.module)
            if table is not None:
                self._tables.move_to_end(definition.module)
                return table
        table = step_table(catalog, definition)
        with self._lock:
            self._tables[definition.module] = table
            while len(self._tables) > self._size:
                self._tables.popitem(last=False)
        return table


# --------------------------------------------------------------------------- #
# Aggregates
# --------------------------------------------------------------------------- #


def ko_shares(catalog, kos: set[str], genomes) -> dict[str, float]:
    if not genomes:
        return dict.fromkeys(kos, 0.0)
    counts: Counter = Counter()
    for g in genomes:
        present = catalog.ko_sets.get(g.bin_id)
        if present:
            counts.update(present & kos)
    return {k: counts[k] / len(genomes) for k in kos}


def step_shares(table: StepTable, genomes) -> list[float]:
    if not genomes or not table.n_steps:
        return [0.0] * table.n_steps
    idx = [g.index for g in genomes]
    return table.complete[idx].mean(axis=0).tolist()


def clade_rows(catalog, table: StepTable, genomes, rank: str, *, min_genomes: int = 3,
               limit: int = 30) -> list[dict[str, Any]]:
    """Per clade at ``rank`` (among ``genomes``): share completing each step and the whole module."""
    groups: dict[str, list[int]] = {}
    for g in genomes:
        groups.setdefault(g.taxonomy.get(rank) or UNCLASSIFIED, []).append(g.index)
    rows = []
    for clade, idx in groups.items():
        if len(idx) < min_genomes:
            continue
        sub = table.complete[idx]
        carriers = int(table.relevant[idx].sum())
        if not carriers:
            continue
        shares = sub.mean(axis=0).tolist() if table.n_steps else []
        whole = float((sub.mean(axis=1) >= 0.75).mean()) if table.n_steps else 0.0
        gap = None
        if shares:
            j = int(np.argmin(shares))
            # the usual gap: one step well below the clade's typical step
            if shares[j] < 0.5 and float(np.median(shares)) >= 0.5:
                gap = j
        rows.append({"clade": clade, "genomes": len(idx), "carriers": carriers, "shares": shares,
                     "complete": whole, "gap": gap})
    rows.sort(key=lambda r: (-r["carriers"], r["clade"]))
    return rows[:limit]


def nearly_complete(catalog, table: StepTable, genomes, *, per_group: int = 12) -> list[dict[str, Any]]:
    """Genomes completing every step but one, grouped by the step they lack."""
    if table.n_steps < 2:
        return []
    groups: dict[int, list] = {}
    for g in genomes:
        row = table.complete[g.index]
        if int(row.sum()) == table.n_steps - 1:
            j = int(np.flatnonzero(~row)[0])
            groups.setdefault(j, []).append(g)
    out = []
    for j, members in groups.items():
        kos: Counter = Counter()
        for g in members:
            kos.update(table.missing.get((g.index, j), []))
        comp = [g.completeness for g in members if g.completeness is not None]
        members.sort(key=lambda g: (-(g.completeness or 0), g.bin_id))
        out.append({"column": j, "step": table.steps[j], "genomes": members[:per_group], "count": len(members),
                    "ids": [g.bin_id for g in members], "missing": kos.most_common(4),
                    "median_completeness": float(np.median(comp)) if comp else None})
    out.sort(key=lambda r: -r["count"])
    return out


# --------------------------------------------------------------------------- #
# Drawing
# --------------------------------------------------------------------------- #


def _shade(share: float) -> str:
    return f"color-mix(in srgb, var(--accent) {SHADE_MIN + SHADE_SPAN * max(0.0, min(1.0, share)):.0f}%, var(--surface-2))"


def _layout(node: Node, ctx: dict[str, Any]):
    """(width, height, draw(x, y) -> list of SVG strings) for one definition node."""
    if node.kind == "ko":
        ko = node.value
        share = ctx["shares"].get(ko, 0.0)
        known = ctx["names"](ko)
        symbol = (known[0].split(",")[0].strip() if known and known[0] else "") or ko
        title = f"{ko}{(' ' + known[0]) if known and known[0] else ''}: {known[1] if known else 'no KEGG name here'} · {share * 100:.0f}% of genomes"

        def draw(x, y, ko=ko, symbol=symbol, title=title, share=share):
            label = symbol if len(symbol) <= 9 else symbol[:8] + "…"
            return [f'<a href="{escape(ctx["ko_url"](ko))}"><g class="pw-ko{" pw-absent" if share == 0 else ""}">'
                    f'<title>{escape(title)}</title>'
                    f'<rect x="{x}" y="{y}" width="{KO_W}" height="{KO_H}" rx="6" style="fill:{_shade(share)}"/>'
                    f'<text x="{x + KO_W / 2}" y="{y + 14}" class="pw-sym">{escape(label)}</text>'
                    f'<text x="{x + KO_W / 2}" y="{y + 26}" class="pw-id">{ko}</text></g></a>']
        return KO_W, KO_H, draw
    if node.kind == "module":
        ref = node.value

        def draw(x, y, ref=ref):
            return [f'<a href="/pathway/{quote(ref)}"><g class="pw-ref"><title>Module {ref}</title>'
                    f'<rect x="{x}" y="{y}" width="{KO_W}" height="{KO_H}" rx="6"/>'
                    f'<text x="{x + KO_W / 2}" y="{y + 20}" class="pw-sym">{ref}</text></g></a>']
        return KO_W, KO_H, draw
    if node.kind == "undefined":
        def draw(x, y):
            return [f'<g class="pw-undef"><title>A step without a KEGG ortholog</title>'
                    f'<rect x="{x}" y="{y}" width="34" height="{KO_H}" rx="6"/>'
                    f'<text x="{x + 17}" y="{y + 20}" class="pw-id">—</text></g>']
        return 34, KO_H, draw
    if node.kind == "alt":
        parts = [_layout(c, ctx) for c in node.children]
        w = max(p[0] for p in parts) + 8
        h = sum(p[1] for p in parts) + GAP * (len(parts) - 1) + 8

        def draw(x, y, parts=parts, w=w, h=h):
            out = [f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="9" class="pw-alt"><title>Alternatives: any one completes this part</title></rect>']
            cy = y + 4
            for pw, ph, pd in parts:
                out += pd(x + 4 + (w - 8 - pw) / 2, cy)
                cy += ph + GAP
            return out
        return w, h, draw
    if node.kind in ("complex", "seq"):
        parts = [_layout(c, ctx) for c in node.children]
        optional = list(node.optional) if node.kind == "complex" else [False] * len(parts)
        gap = PLUS_GAP if node.kind == "complex" else ARROW_GAP
        w = sum(p[0] for p in parts) + gap * (len(parts) - 1)
        h = max(p[1] for p in parts)

        def draw(x, y, parts=parts, optional=optional, gap=gap, h=h, kind=node.kind):
            out, cx = [], x
            for k, (pw, ph, pd) in enumerate(parts):
                if k:
                    mid = y + h / 2
                    if kind == "complex":
                        out.append(f'<text x="{cx - gap / 2}" y="{mid + 4}" class="pw-plus">{"−" if optional[k] else "+"}</text>')
                    else:
                        out.append(f'<path d="M{cx - gap + 3},{mid} L{cx - 4},{mid}" class="pw-arrow"/>')
                drawn = pd(cx, y + (h - ph) / 2)
                if optional[k]:
                    drawn = [f'<g class="pw-optional"><title>Optional component</title>'] + drawn + ["</g>"]
                out += drawn
                cx += pw + gap
            return out
        return w, h, draw
    return 0, 0, lambda x, y: []


def diagram(definition: ModuleDefinition, *, shares: dict[str, float], step_complete: list[float],
            steps: list[int], names: Callable[[str], Any], ko_url: Callable[[str], str]) -> Markup:
    """The module definition as an SVG step diagram."""
    ctx = {"shares": shares, "names": names, "ko_url": ko_url}
    share_of = dict(zip(steps, step_complete))
    top = definition.tree.children if definition.tree.kind == "seq" else (definition.tree,)
    laid = [(i, _layout(node, ctx)) for i, node in enumerate(top, 1)]
    height = max((p[1] for _, p in laid), default=KO_H) + HEADER_H + 10
    x, parts = 10, []
    for k, (i, (w, h, draw)) in enumerate(laid):
        col_w = max(w, 58)
        if k:
            mid = HEADER_H + (height - HEADER_H) / 2
            parts.append(f'<path d="M{x - ARROW_GAP + 4},{mid} L{x - 6},{mid}" class="pw-arrow strong"/>')
        share = share_of.get(i)
        head = f"Step {i}"
        parts.append(f'<text x="{x + col_w / 2}" y="13" class="pw-step">{head}</text>')
        if share is not None:
            parts.append(f'<rect x="{x + col_w / 2 - 22}" y="19" width="44" height="5" rx="2.5" class="pw-track"/>'
                         f'<rect x="{x + col_w / 2 - 22}" y="19" width="{44 * share:.1f}" height="5" rx="2.5" class="pw-fill">'
                         f'<title>{share * 100:.0f}% of genomes complete step {i}</title></rect>')
        parts += draw(x + (col_w - w) / 2, HEADER_H + (height - HEADER_H - 10 - h) / 2)
        x += col_w + ARROW_GAP
    width = x - ARROW_GAP + 10
    return Markup(f'<svg class="pw-diagram" viewBox="0 0 {width:.0f} {height:.0f}" width="{width:.0f}" height="{height:.0f}" '
                  f'role="img" aria-label="Module step diagram">{"".join(parts)}</svg>')


def heatmap(rows: list[dict[str, Any]], steps: list[int], rank: str, url) -> Markup:
    """Clade × step SVG: the share of each clade's genomes completing each step."""
    if not rows or not steps:
        return Markup("")
    label_w, cell, gap = 190, 34, 2
    head = 26
    cols = len(steps) + 1
    width = label_w + cols * (cell + gap) + 60
    height = head + len(rows) * (cell - 8 + gap) + 6
    rh = cell - 8
    parts = []
    for j, s in enumerate(steps):
        parts.append(f'<text x="{label_w + j * (cell + gap) + cell / 2}" y="16" class="pw-col">{s}</text>')
    xa = label_w + len(steps) * (cell + gap) + 6
    parts.append(f'<text x="{xa + cell / 2}" y="16" class="pw-col">≥75%</text>')
    for r, row in enumerate(rows):
        y = head + r * (rh + gap)
        name = row["clade"] if len(row["clade"]) <= 26 else row["clade"][:25] + "…"
        parts.append(f'<a href="{escape(url("taxa", rank, row["clade"]))}"><text x="{label_w - 8}" y="{y + rh / 2 + 4}" '
                     f'class="pw-row">{escape(name)}<title>{escape(row["clade"])}: {row["genomes"]} genomes</title></text></a>')
        for j, share in enumerate(row["shares"]):
            x = label_w + j * (cell + gap)
            cls = "pw-cell gap" if row["gap"] == j else "pw-cell"
            parts.append(f'<rect x="{x}" y="{y}" width="{cell}" height="{rh}" rx="3" class="{cls}" style="fill:{_shade(share)}">'
                         f'<title>{escape(row["clade"])} · step {steps[j]}: {share * 100:.0f}% of {row["genomes"]} genomes</title></rect>')
            if share >= 0.995 or share == 0:
                continue
            parts.append(f'<text x="{x + cell / 2}" y="{y + rh / 2 + 4}" class="pw-val">{share * 100:.0f}</text>')
        parts.append(f'<rect x="{xa}" y="{y}" width="{cell}" height="{rh}" rx="3" class="pw-cell whole" style="fill:{_shade(row["complete"])}">'
                     f'<title>{escape(row["clade"])}: module ≥ 75% complete in {row["complete"] * 100:.0f}% of genomes</title></rect>'
                     f'<text x="{xa + cell / 2}" y="{y + rh / 2 + 4}" class="pw-val">{row["complete"] * 100:.0f}</text>'
                     f'<text x="{xa + cell + 8}" y="{y + rh / 2 + 4}" class="pw-n">{row["genomes"]}</text>')
    return Markup(f'<svg class="pw-heat" viewBox="0 0 {width} {height}" width="{width}" height="{height}" role="img" '
                  f'aria-label="Step completion by {escape(rank)}">{"".join(parts)}</svg>')


def scope_genomes(catalog, clade: str):
    """(label, genomes) for ``rank:name``, or the whole dataset."""
    if clade and ":" in clade:
        rank, name = clade.split(":", 1)
        if rank in RANKS:
            members = catalog.clade(rank, name)
            if members:
                return rank, name, members
    return None, None, list(catalog.genomes)
