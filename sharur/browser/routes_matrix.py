"""Presence/absence matrix: genomes × features, for one clade or two.

Rows are the genomes of one or two groups (a clade, a genome, or the whole
dataset); columns are features of one kind: KEGG orthologs, Pfam families,
KEGG modules (cells hold completeness), curated systems, or function labels.
Columns are picked explicitly, as a KEGG module's KOs in step order, or
automatically:

- ``variable``: features present in about half of one group (largest p(1-p)).
- ``absences``: the most common features (in at least half of a group) whose
  absences exceed what genome incompleteness explains, at BH q < 0.05 across
  the features tested (needs completeness; see below).
- ``differential``: largest prevalence differences between two groups, with
  two-sided Fisher's exact p and Benjamini-Hochberg q over every feature
  carried by either group.

Completeness control: if a feature were in every genome of a group and each
genome reported it with probability equal to its completeness, the number of
carriers would follow a Poisson-binomial distribution over those
completeness values. Its lower tail at the observed count says whether
incompleteness explains the absences. Absences beyond it come from gene loss,
gene calls missing from the dataset, or annotation thresholds that miss
divergent homologs. The test fits single-copy genes; it is reported for KOs,
Pfam families and function labels. Genomes whose gene calls are mostly
missing (``protein_deficits``) are left out of it.
"""

from __future__ import annotations

import csv
import io
import re
import threading
from dataclasses import dataclass, field
from types import SimpleNamespace
from typing import Any

import numpy as np
from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.responses import HTMLResponse, Response

from sharur.browser.catalog import BIOLOGICAL, RANKS, UNCLASSIFIED, Genome
from sharur.predicates.vocabulary import PREDICATE_BY_ID

KINDS = {"ko": "KEGG orthologs", "pfam": "Pfam families", "module": "KEGG modules",
         "system": "Curated systems", "function": "Function labels"}
PICKS = {"variable": "Most variable", "absences": "Absences beyond incompleteness",
         "differential": "Largest differences"}
TESTABLE = ("ko", "pfam", "function")   # single-gene features: completeness test applies
MODULE_PRESENT = 0.75
MAX_FEATURES = 150
ALPHA = 0.05


# --------------------------------------------------------------------------- #
# Feature sets: one sparse genome × feature table per kind
# --------------------------------------------------------------------------- #


@dataclass
class FeatureSet:
    kind: str
    ids: list[str]
    labels: list[str]               # short column labels
    names: list[str]                # full names
    bin_idx: np.ndarray             # genome index per (genome, feature) pair
    feat_idx: np.ndarray
    value: np.ndarray               # copies, calls, or module completeness
    threshold: float = 1e-9         # value at which a feature counts as present
    selectable: np.ndarray | None = None  # features offered by automatic picks
    index: dict[str, int] = field(init=False)

    def __post_init__(self) -> None:
        self.index = {f: i for i, f in enumerate(self.ids)}
        if self.selectable is None:
            self.selectable = np.ones(len(self.ids), dtype=bool)

    def carriers(self, genome_mask: np.ndarray) -> np.ndarray:
        """Genomes carrying each feature, among ``genome_mask``."""
        keep = genome_mask[self.bin_idx] & (self.value >= self.threshold)
        return np.bincount(self.feat_idx[keep], minlength=len(self.ids))


def _sparse(catalog, rows) -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray]:
    ids = sorted({f for _, f, _ in rows})
    index = {f: i for i, f in enumerate(ids)}
    by_bin = {g.bin_id: g.index for g in catalog.genomes}
    bins = np.fromiter((by_bin.get(b, -1) for b, _, _ in rows), dtype=np.int32, count=len(rows))
    feats = np.fromiter((index[f] for _, f, _ in rows), dtype=np.int32, count=len(rows))
    values = np.fromiter((v for _, _, v in rows), dtype=np.float32, count=len(rows))
    keep = bins >= 0
    return ids, bins[keep], feats[keep], values[keep]


def _load_ko(ctx) -> FeatureSet:
    with ctx.lock:
        rows = ctx.store.execute(
            """SELECT p.bin_id, a.accession, COUNT(DISTINCT a.protein_id)
               FROM annotations a JOIN proteins p USING (protein_id)
               WHERE LOWER(a.source) IN ('kofam', 'kegg') AND regexp_matches(a.accession, '^K[0-9]{5}$')
               GROUP BY 1, 2""")
    ids, bins, feats, values = _sparse(ctx.catalog, rows)
    labels, names = [], []
    for ko in ids:
        known = ctx.ko_names.get(ko) if getattr(ctx, "ko_names", None) else None
        symbol = known[0].split(",")[0].strip() if known else ""
        labels.append(f"{ko} {symbol}".strip())
        names.append(re.sub(r"\s*\[EC:[^\]]*\]", "", known[1]) if known else ko)
    return FeatureSet("ko", ids, labels, names, bins, feats, values)


def _load_pfam(ctx) -> FeatureSet:
    with ctx.lock:
        rows = ctx.store.execute(
            """SELECT p.bin_id, split_part(a.accession, '.', 1), COUNT(DISTINCT a.protein_id)
               FROM annotations a JOIN proteins p USING (protein_id)
               WHERE LOWER(a.source) = 'pfam' GROUP BY 1, 2""")
        named = dict(ctx.store.execute(
            """SELECT split_part(accession, '.', 1), ANY_VALUE(COALESCE(NULLIF(name, ''), accession))
               FROM annotations WHERE LOWER(source) = 'pfam' GROUP BY 1"""))
    ids, bins, feats, values = _sparse(ctx.catalog, rows)
    domains = ctx.catalog.domains
    labels = [named.get(acc, acc) for acc in ids]
    names = [f"{named.get(acc, acc)}: {domains[acc]['description']}" if acc in domains and domains[acc]["description"]
             else named.get(acc, acc) for acc in ids]
    return FeatureSet("pfam", ids, labels, names, bins, feats, values)


def _load_function(ctx) -> FeatureSet:
    catalog = ctx.catalog
    if catalog.pair_bin is None:
        raise LookupError("This dataset has no function labels.")
    ids = list(catalog.predicates)
    defs = [PREDICATE_BY_ID.get(p) for p in ids]
    selectable = np.array([d is not None and d.category in BIOLOGICAL for d in defs], dtype=bool)
    return FeatureSet("function", ids, [d.name if d else p for d, p in zip(defs, ids)],
                      [d.description if d and getattr(d, "description", "") else (d.name if d else p)
                       for d, p in zip(defs, ids)],
                      catalog.pair_bin, catalog.pair_pred, catalog.pair_count.astype(np.float32),
                      selectable=selectable)


def _load_module(ctx) -> FeatureSet:
    catalog = ctx.catalog
    if catalog.module_completeness is None:
        if not catalog.ready.is_set():
            raise LookupError("Pathway completeness is still being computed; reload in a moment.")
        raise LookupError("No KEGG module definitions (run `sharur setup-kegg`).")
    matrix = catalog.module_completeness
    rows, cols = np.nonzero(matrix)
    ids = list(catalog.module_ids)
    return FeatureSet("module", ids, ids, [catalog.modules[m].name for m in ids], rows.astype(np.int32),
                      cols.astype(np.int32), matrix[rows, cols].astype(np.float32), threshold=MODULE_PRESENT)


def _load_system(ctx) -> FeatureSet:
    calls: dict[tuple[str, str], int] = {}
    labels: dict[str, tuple[str, str]] = {}
    for s in ctx.catalog.systems:
        key = f"{s['kind']}:{s['type']}"
        calls[(s["bin_id"], key)] = calls.get((s["bin_id"], key), 0) + 1
        labels[key] = (s["type"], f"{s['type']} ({s['kind']} system)")
    cctyper = getattr(ctx, "cctyper_systems", None)
    for s in (cctyper() if cctyper else []):
        if s["confident"]:
            key = f"crispr:{s['prediction']}"
            calls[(s["bin_id"], key)] = calls.get((s["bin_id"], key), 0) + 1
            labels[key] = (f"CRISPR-Cas {s['prediction']}", f"CRISPR-Cas {s['prediction']} (CRISPRCasTyper)")
    if not calls:
        raise LookupError("This dataset has no curated system calls.")
    rows = [(b, k, float(n)) for (b, k), n in calls.items()]
    ids, bins, feats, values = _sparse(ctx.catalog, rows)
    return FeatureSet("system", ids, [labels[k][0] for k in ids], [labels[k][1] for k in ids], bins, feats, values)


LOADERS = {"ko": _load_ko, "pfam": _load_pfam, "function": _load_function, "module": _load_module,
           "system": _load_system}


class FeatureSets:
    """Feature sets loaded on first use and kept."""

    def __init__(self, ctx) -> None:
        self.ctx, self._sets, self._lock = ctx, {}, threading.Lock()

    def get(self, kind: str) -> FeatureSet:
        with self._lock:
            if kind not in self._sets:
                fs = LOADERS[kind](self.ctx)
                # modules and Pfam names fill in after the background load; keep those only once complete
                if kind not in ("module", "pfam") or self.ctx.catalog.ready.is_set():
                    self._sets[kind] = fs
                return fs
            return self._sets[kind]


# --------------------------------------------------------------------------- #
# Groups
# --------------------------------------------------------------------------- #


@dataclass
class Group:
    key: str
    label: str
    rank: str | None            # None for a genome or the whole dataset
    genomes: list[Genome]


def resolve_group(ctx, token: str) -> Group | None:
    """A genome id, ``rank:name``, a taxon name at any rank (broadest wins), or ``all``."""
    from sharur.browser.routes_search import resolve_scope  # noqa: PLC0415

    catalog = ctx.catalog
    token = (token or "").strip()
    if not token:
        return None
    if token.lower() == "all":
        return Group("all", "All genomes", None, list(catalog.genomes))
    if token in catalog.by_bin:
        return Group(token, token, None, [catalog.by_bin[token]])
    if ":" in token:
        rank, name = token.split(":", 1)
        if rank in RANKS:
            members = catalog.clade(rank, name)
            return Group(token, name, rank, members) if members else None
    scope = resolve_scope(ctx, token)
    if scope is None:
        return None
    kind, label, bins = scope
    genomes = [catalog.by_bin[b] for b in bins if b in catalog.by_bin]
    if kind in RANKS:
        return Group(f"{kind}:{label}", label, kind, genomes)
    return Group(token, label, None, genomes)


# --------------------------------------------------------------------------- #
# Statistics
# --------------------------------------------------------------------------- #


def poisson_binomial_cdf(probs: np.ndarray) -> np.ndarray:
    """P(K <= k) for k = 0..n, K the number of successes of independent trials."""
    pmf = np.ones(1)
    for p in probs:
        pmf = np.convolve(pmf, (1.0 - p, p))
    return np.minimum(np.cumsum(pmf), 1.0)


DEFICIT = 0.5   # proteins per assembly kb (prokaryotes carry about 1), or half the completeness-implied count


def assembly_kb(ctx) -> dict[str, float]:
    """Assembly length per genome (kb) from uncompressed FASTA file sizes, cached.

    A FASTA file holds about 1.3% newlines and headers beyond its bases.
    """
    cached = getattr(ctx, "_assembly_kb", None)
    if cached is None:
        import os  # noqa: PLC0415

        cached = {}
        for bin_id, path in getattr(getattr(ctx, "assemblies", None), "paths", {}).items():
            if not str(path).endswith(".gz"):
                try:
                    cached[bin_id] = os.path.getsize(path) / 1013.0
                except OSError:
                    pass
        ctx._assembly_kb = cached
    return cached


def protein_deficits(genomes: list[Genome], reference: list[Genome], kb: dict[str, float] | None = None) -> set[str]:
    """Genomes whose gene calls are mostly missing.

    With the assembly at hand: under ``DEFICIT`` proteins per kb. Otherwise:
    under ``DEFICIT`` × the proteins its completeness implies, from the median
    proteins-per-completeness of ``reference``. Completeness estimated on the
    assembly says little about absences in such a genome.
    """
    kb = kb or {}

    def ratio(g: Genome) -> float | None:
        return g.proteins / (g.completeness / 100.0) if g.completeness and g.proteins else None

    ratios = [r for r in map(ratio, reference) if r]
    per_complete = float(np.median(ratios)) if len(ratios) >= 5 else None
    flagged = set()
    for g in genomes:
        if kb.get(g.bin_id):
            if g.proteins < DEFICIT * kb[g.bin_id]:
                flagged.add(g.bin_id)
        elif per_complete and g.completeness and g.proteins < DEFICIT * per_complete * g.completeness / 100.0:
            flagged.add(g.bin_id)
    return flagged


def completeness_test(genomes: list[Genome], carriers_known: np.ndarray) -> dict[str, Any] | None:
    """Per feature: expected carriers if universal, and P(carriers <= observed).

    Uses the genomes with a completeness estimate; ``carriers_known`` counts
    carriers among those genomes.
    """
    probs = np.array([min(max(g.completeness / 100.0, 0.0), 1.0) for g in genomes if g.completeness is not None])
    if not len(probs):
        return None
    cdf = poisson_binomial_cdf(probs)
    return {"n": len(probs), "expected": float(probs.sum()), "p": cdf[np.minimum(carriers_known, len(probs))]}


def pval(p: float | None) -> str:
    if p is None or p != p:
        return "–"
    return f"{p:.3f}" if p >= 0.001 else f"{p:.1e}"


def fisher_two_sided(ka: int, na: int, kb: int, nb: int) -> float:
    from scipy.stats import fisher_exact  # noqa: PLC0415

    return float(fisher_exact([[ka, na - ka], [kb, nb - kb]])[1])


def bh_qvalues(pvalues: np.ndarray, m: int) -> np.ndarray:
    """Benjamini-Hochberg q for ``pvalues`` within ``m`` tests (the untested count as p = 1)."""
    if not len(pvalues):
        return pvalues
    order = np.argsort(pvalues)
    ranked = pvalues[order] * m / np.arange(1, len(pvalues) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty_like(ranked)
    q[order] = np.minimum(ranked, 1.0)
    return q


# --------------------------------------------------------------------------- #
# Feature choice and ordering
# --------------------------------------------------------------------------- #


def module_kos(catalog, module: str) -> tuple[list[str], list[int]]:
    """A module's KOs in definition order, with the top-level step of each."""
    from sharur.modules import kos_in  # noqa: PLC0415

    definition = catalog.modules.get(module)
    if definition is None:
        return [], []
    tree = definition.tree
    steps = tree.children if tree.kind == "seq" else (tree,)
    kos, step_of = [], []
    for i, step in enumerate(steps):
        for ko in sorted(kos_in(step), key=lambda k: definition.definition.find(k)):
            if ko not in kos:
                kos.append(ko)
                step_of.append(i)
    return kos, step_of


def _mask(n: int, genomes: list[Genome]) -> np.ndarray:
    mask = np.zeros(n, dtype=bool)
    mask[[g.index for g in genomes]] = True
    return mask


def _jaccard_order(presence: np.ndarray) -> list[int]:
    """Row order from average-linkage clustering on Jaccard distance."""
    n = presence.shape[0]
    if n < 3:
        return list(range(n))
    from scipy.cluster.hierarchy import leaves_list, linkage  # noqa: PLC0415
    from scipy.spatial.distance import pdist  # noqa: PLC0415

    dist = np.nan_to_num(pdist(presence.astype(bool), "jaccard"), nan=0.0)
    return [int(i) for i in leaves_list(linkage(dist, "average"))]


def _taxonomy_key(g: Genome) -> tuple:
    return (*[(g.taxonomy.get(r) or "~") for r in RANKS], -(g.completeness or 0.0), g.bin_id)


def build_matrix(ctx, fs: FeatureSet, groups: list[Group], *, features: list[str] | None = None,
                 module: str = "", pick: str = "", n: int = 40, order: str = "taxonomy",
                 min_completeness: float = 0.0) -> dict[str, Any]:
    catalog = ctx.catalog
    n_all = len(catalog.genomes)
    notes: list[str] = []
    dropped = 0
    if min_completeness > 0:
        kept_groups = []
        for grp in groups:
            kept = [g for g in grp.genomes if g.completeness is not None and g.completeness >= min_completeness]
            dropped += len(grp.genomes) - len(kept)
            kept_groups.append(Group(grp.key, grp.label, grp.rank, kept))
        groups = kept_groups
    groups = [grp for grp in groups if grp.genomes]
    if not groups:
        return {"empty": f"No genomes left at completeness ≥ {min_completeness:g}%."}

    masks = [_mask(n_all, grp.genomes) for grp in groups]
    reference = [g for grp in groups for g in grp.genomes]
    if len(reference) < 20:
        reference = list(catalog.genomes)
    deficit = protein_deficits([g for grp in groups for g in grp.genomes], reference, assembly_kb(ctx))
    testable = [[g for g in grp.genomes if g.completeness is not None and g.bin_id not in deficit] for grp in groups]
    known = [_mask(n_all, gs) for gs in testable]
    carriers = [fs.carriers(m) for m in masks]
    carriers_known = [fs.carriers(m) for m in known]
    sizes = [len(grp.genomes) for grp in groups]
    prevalence = [c / s for c, s in zip(carriers, sizes)]
    tests = [completeness_test(gs, ck) if fs.kind in TESTABLE else None for gs, ck in zip(testable, carriers_known)]
    for test, prev in zip(tests, prevalence):
        if test is not None:
            # BH across every feature tested in the group (those in at least half its genomes)
            tested = np.flatnonzero(prev >= 0.5)
            test["q"] = np.full(len(prev), np.nan)
            test["q"][tested] = bh_qvalues(test["p"][tested], len(tested))

    steps: list[int] | None = None
    fisher: dict[int, tuple[float, float]] = {}
    if module and fs.kind == "ko":
        kos, step_of = module_kos(catalog, module)
        if not kos:
            return {"empty": f"Unknown KEGG module {module}."}
        chosen = [fs.index[k] if k in fs.index else -1 for k in kos]
        absent_ids = [k for k, i in zip(kos, chosen) if i < 0]
        if absent_ids:
            notes.append(f"{len(absent_ids)} of the module's {len(kos)} KOs occur in no genome of this dataset: "
                         + ", ".join(absent_ids[:12]) + ("…" if len(absent_ids) > 12 else ""))
        keep = [k for k, i in enumerate(chosen) if i >= 0]
        chosen, steps = [chosen[k] for k in keep], [step_of[k] for k in keep]
        pick = "module"
    elif features:
        chosen = [fs.index[f] for f in features if f in fs.index]
        missing = [f for f in features if f not in fs.index]
        if missing:
            notes.append(f"Not in this dataset: {', '.join(missing[:12])}")
        pick = "explicit"
    else:
        selectable = fs.selectable & (sum(carriers) > 0)
        if len(groups) == 2:
            pick = "differential"
            delta = prevalence[0] - prevalence[1]
            candidates = np.flatnonzero(selectable)
            ranked = candidates[np.argsort(-np.abs(delta[candidates]), kind="stable")]
            tested = ranked[: max(4 * n, 200)]
            p = np.array([fisher_two_sided(int(carriers[0][k]), sizes[0], int(carriers[1][k]), sizes[1])
                          for k in tested])
            q = bh_qvalues(p, len(candidates))
            fisher = {int(k): (float(pk), float(qk)) for k, pk, qk in zip(tested, p, q)}
            top = list(ranked[:n])
            chosen = sorted((k for k in top if delta[k] > 0), key=lambda k: -delta[k]) + \
                sorted((k for k in top if delta[k] <= 0), key=lambda k: delta[k])
        elif pick == "absences":
            test = tests[0]
            if test is None:
                return {"empty": "No completeness estimates for these genomes, so absences can't be weighed "
                                 "against incompleteness. Import them with `sharur import-quality`."}
            # the most common features whose absences still exceed incompleteness
            p0 = prevalence[0]
            candidates = np.flatnonzero(selectable & (p0 >= 0.5) & (test["q"] < ALPHA))
            chosen = [int(k) for k in candidates[np.argsort(-p0[candidates], kind="stable")][:n]]
        else:
            pick = "variable"
            p0 = prevalence[0]
            candidates = np.flatnonzero(selectable & (p0 > 0) & (p0 < 1))
            chosen = [int(k) for k in candidates[np.argsort(-(p0[candidates] * (1 - p0[candidates])), kind="stable")][:n]]
    chosen = chosen[:MAX_FEATURES]
    if not chosen:
        return {"empty": "No features to show for this selection."}
    if len(groups) == 2 and not fisher:
        # explicit or module columns: correct over the features shown
        p = np.array([fisher_two_sided(int(carriers[0][k]), sizes[0], int(carriers[1][k]), sizes[1]) for k in chosen])
        fisher = {int(k): (float(pk), float(qk)) for k, pk, qk in zip(chosen, p, bh_qvalues(p, len(chosen)))}

    # dense values for the displayed genomes and features
    genome_rows = [g for grp in groups for g in grp.genomes]
    group_of = [gi for gi, grp in enumerate(groups) for _ in grp.genomes]
    row_of = np.full(n_all, -1, dtype=np.int64)
    # a genome in both groups appears in each; index rows per group
    values = np.zeros((len(genome_rows), len(chosen)), dtype=np.float32)
    col_of = np.full(len(fs.ids), -1, dtype=np.int64)
    col_of[chosen] = np.arange(len(chosen))
    start = 0
    for grp in groups:
        row_of[:] = -1
        row_of[[g.index for g in grp.genomes]] = np.arange(start, start + len(grp.genomes))
        keep = (row_of[fs.bin_idx] >= 0) & (col_of[fs.feat_idx] >= 0)
        values[row_of[fs.bin_idx[keep]], col_of[fs.feat_idx[keep]]] = fs.value[keep]
        start += len(grp.genomes)
    presence = values >= fs.threshold

    # feature order: clustered for automatic single-group picks
    if pick in ("variable", "absences") and len(chosen) > 2:
        col_order = _jaccard_order(presence.T)
        chosen = [chosen[k] for k in col_order]
        values, presence = values[:, col_order], presence[:, col_order]

    # genome order within each group
    order_idx: list[int] = []
    start = 0
    for grp in groups:
        rows = list(range(start, start + len(grp.genomes)))
        if order == "cluster":
            local = _jaccard_order(presence[rows])
            rows = [rows[k] for k in local]
        elif order == "completeness":
            rows.sort(key=lambda r: (-(genome_rows[r].completeness or 0.0), genome_rows[r].bin_id))
        else:
            rows.sort(key=lambda r: _taxonomy_key(genome_rows[r]))
        order_idx.extend(rows)
        start += len(grp.genomes)
    genome_rows = [genome_rows[r] for r in order_idx]
    group_of = [group_of[r] for r in order_idx]
    values = values[order_idx]

    color_rank = _color_rank(catalog, groups)
    table = []
    for j, k in enumerate(chosen):
        row: dict[str, Any] = {"id": fs.ids[k], "label": fs.labels[k], "name": fs.names[k], "groups": []}
        for gi in range(len(groups)):
            col = values[np.array(group_of) == gi, j]
            present = col >= fs.threshold
            entry = {"carriers": int(carriers[gi][k]), "n": sizes[gi], "share": float(prevalence[gi][k]),
                     "median": float(np.median(col[present])) if present.any() else None}
            test = tests[gi]
            if test is not None and prevalence[gi][k] >= 0.5:
                entry["expected"] = test["expected"]
                entry["observed"] = int(carriers_known[gi][k])
                entry["tested_n"] = test["n"]
                entry["p"] = float(test["p"][k])
                entry["q"] = float(test["q"][k])
            row["groups"].append(entry)
        if len(groups) == 2:
            row["delta"] = float(prevalence[0][k] - prevalence[1][k])
            row["fisher_p"], row["fisher_q"] = fisher[int(k)]
        if steps is not None:
            row["step"] = steps[j]
        table.append(row)

    if dropped:
        notes.append(f"{dropped} genomes below {min_completeness:g}% completeness left out.")
    missing_comp = sum(1 for g in genome_rows if g.completeness is None)
    return {
        "kind": fs.kind, "pick": pick, "order": order, "threshold": fs.threshold, "color_rank": color_rank,
        "groups": [{"key": grp.key, "label": grp.label, "rank": grp.rank, "n": len(grp.genomes)} for grp in groups],
        "genomes": [{"id": g.bin_id, "label": g.label, "group": gi, "completeness": g.completeness,
                     "proteins": g.proteins, "deficit": g.bin_id in deficit,
                     "contamination": g.contamination, "clade": g.taxonomy.get(color_rank) or UNCLASSIFIED,
                     "lineage": "; ".join(n for _, n in g.lineage if n != UNCLASSIFIED)}
                    for g, gi in zip(genome_rows, group_of)],
        "features": table, "values": values, "notes": notes, "tested": fs.kind in TESTABLE,
        "missing_completeness": missing_comp,
        "deficit": sorted(deficit & {g.bin_id for g in genome_rows}),
    }


def _color_rank(catalog, groups: list[Group]) -> str:
    """The broadest rank that splits the genomes within a group."""
    found = []
    for grp in groups:
        child, _ = catalog.children(grp.genomes, grp.rank)
        if child:
            found.append(child)
    return min(found, key=RANKS.index) if found else "order"


# --------------------------------------------------------------------------- #
# Links
# --------------------------------------------------------------------------- #


def feature_url(url, kind: str, feature: str) -> str:
    if kind == "pfam":
        return url("domain", feature)
    if kind == "function":
        return url("function", feature)
    if kind == "module":
        return url("pathway", feature)
    if kind == "system":
        family, name = feature.split(":", 1)
        return f"/crispr-cas/calls?subtype={name}" if family == "crispr" else url("system", family, name)
    return f"/search?q={feature}"


def cell_url_template(kind: str) -> str:
    """URL for one genome × feature cell; {g} and {f} are filled in the page."""
    if kind == "function":
        return "/function/{f}?genome={g}"
    if kind in ("module", "system"):
        return "/genome/{g}"
    return "/search?q={f}%20in%20{g}"


def external_url(kind: str, feature: str) -> str | None:
    if kind == "ko":
        return f"https://www.kegg.jp/entry/{feature}"
    if kind == "module":
        return f"https://www.kegg.jp/entry/{feature}"
    if kind == "pfam":
        return f"https://www.ebi.ac.uk/interpro/entry/pfam/{feature}/"
    return None


# --------------------------------------------------------------------------- #
# Routes
# --------------------------------------------------------------------------- #


def _tsv(result: dict[str, Any]) -> str:
    out = io.StringIO()
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(["genome", "group", "completeness", "contamination", "lineage",
                     *[f["id"] for f in result["features"]]])
    labels = [grp["label"] for grp in result["groups"]]
    for g, row in zip(result["genomes"], result["values"]):
        writer.writerow([g["id"], labels[g["group"]], "" if g["completeness"] is None else g["completeness"],
                         "" if g["contamination"] is None else g["contamination"], g["lineage"],
                         *[f"{v:g}" for v in row]])
    return out.getvalue()


def register(app: FastAPI, ctx: SimpleNamespace) -> None:
    """Add /matrix and /matrix.tsv. ``ctx``: store, lock, catalog, render, url, ko_names, cctyper_systems."""
    sets = FeatureSets(ctx)
    ctx.feature_sets = sets

    def compute(a: str, b: str, kind: str, features: str, module: str, pick: str, n: int, order: str,
                min_completeness: float) -> tuple[list[Group], dict[str, Any] | None, str]:
        groups = []
        for token in (a, b):
            if token:
                grp = resolve_group(ctx, token)
                if grp is None:
                    raise HTTPException(404, f"Unknown genome or clade: {token}")
                groups.append(grp)
        if not groups:
            return [], None, ""
        try:
            fs = sets.get(kind)
        except LookupError as exc:
            return groups, None, str(exc)
        wanted = [f for f in dict.fromkeys(x.strip() for x in re.split(r"[,\s]+", features)) if f]
        result = build_matrix(ctx, fs, groups, features=wanted, module=module.strip(), pick=pick, n=n, order=order,
                              min_completeness=min_completeness)
        return groups, result, result.get("empty", "")

    def params(kind, pick, order, n):
        kind = kind if kind in KINDS else "ko"
        order = order if order in ("taxonomy", "cluster", "completeness") else "taxonomy"
        return kind, pick if pick in PICKS else "", order, max(5, min(n, MAX_FEATURES))

    @app.get("/matrix", response_class=HTMLResponse)
    def matrix_page(request: Request, a: str = Query("", max_length=500), b: str = Query("", max_length=500),
                    kind: str = Query("ko"), features: str = Query("", max_length=4000),
                    module: str = Query("", max_length=20), pick: str = Query(""), n: int = Query(40),
                    order: str = Query("taxonomy"), min_completeness: float = Query(0.0, ge=0, le=100)):
        kind, pick, order, n = params(kind, pick, order, n)
        groups, result, message = compute(a, b, kind, features, module, pick, n, order, min_completeness)
        payload = None
        if result and not message:
            from sharur.browser.charts import color_for  # noqa: PLC0415

            genomes = [{**g, "color": color_for(g["clade"]) if g["clade"] != UNCLASSIFIED else None}
                       for g in result["genomes"]]
            payload = {
                "kind": kind, "threshold": result["threshold"],
                "groups": result["groups"], "genomes": genomes,
                "features": [{"id": f["id"], "label": f["label"], "name": f["name"], "step": f.get("step"),
                              "url": feature_url(ctx.url, kind, f["id"])} for f in result["features"]],
                "values": [[round(float(v), 3) for v in row] for row in result["values"]],
                "cell_url": cell_url_template(kind),
            }
        return ctx.render(request, "matrix.html", "taxa", a=a, b=b, kind=kind, features=features, module=module,
                          pick=pick, n=n, order=order, min_completeness=min_completeness, groups=groups,
                          result=result, message=message, payload=payload, KINDS=KINDS, PICKS=PICKS,
                          external_url=external_url, feature_url=lambda f: feature_url(ctx.url, kind, f),
                          ALPHA=ALPHA, pval=pval)

    @app.get("/matrix.tsv")
    def matrix_tsv(a: str = Query("", max_length=500), b: str = Query("", max_length=500), kind: str = Query("ko"),
                   features: str = Query("", max_length=4000), module: str = Query("", max_length=20),
                   pick: str = Query(""), n: int = Query(40), order: str = Query("taxonomy"),
                   min_completeness: float = Query(0.0, ge=0, le=100)):
        kind, pick, order, n = params(kind, pick, order, n)
        _, result, message = compute(a, b, kind, features, module, pick, n, order, min_completeness)
        if not result or message:
            raise HTTPException(400, message or "Choose a genome or clade (a=...)")
        return Response(_tsv(result), media_type="text/tab-separated-values",
                        headers={"Content-Disposition": f'attachment; filename="sharur_matrix_{kind}.tsv"'})
