"""What a clade's gene content looks like, and what sets it apart.

For a clade's genomes and one feature kind (KEGG orthologs, else Pfam
families):

- ``spectrum``: features binned by the share of the clade's genomes carrying
  them (the gene-content frequency spectrum).
- ``core``: features in every genome with gene calls (genomes whose gene
  calls are mostly missing left out), and features consistent with every
  genome given completeness (the Poisson-binomial test of routes_matrix, at
  BH q >= 0.05, on genomes with completeness and gene calls).
- ``signature``: features in most of the clade (>= 75%) and rare elsewhere
  (<= 10%), ranked by the difference.
- ``absent``: features common elsewhere (>= 75%) and rare in the clade
  (<= 10%). Absence of an annotation can also reflect divergent homologs that
  annotation thresholds miss.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from sharur.browser.routes_matrix import (ALPHA, bh_qvalues, completeness_test, protein_deficits,
                                          assembly_kb)

BINS = 10
SIGNATURE_IN, SIGNATURE_OUT = 0.75, 0.10
TOP = 15


def clade_content(ctx, sets, genomes, *, top: int = TOP) -> dict[str, Any] | None:
    catalog = ctx.catalog
    if len(genomes) < 3:
        return None
    fs = None
    for kind in ("ko", "pfam"):
        try:
            fs = sets.get(kind)
            break
        except LookupError:
            continue
    if fs is None or not len(fs.ids):
        return None
    n_all = len(catalog.genomes)
    inside = np.zeros(n_all, dtype=bool)
    inside[[g.index for g in genomes]] = True
    outside = ~inside
    n_in, n_out = int(inside.sum()), int(outside.sum())
    c_in, c_out = fs.carriers(inside), fs.carriers(outside)
    p_in = c_in / n_in
    p_out = c_out / n_out if n_out else np.zeros_like(p_in)
    present = c_in > 0

    edges = np.linspace(0, 1, BINS + 1)
    spectrum = []
    for k in range(BINS):
        lo, hi = edges[k], edges[k + 1]
        sel = present & (p_in > lo) & (p_in <= hi)
        spectrum.append({"lo": float(lo), "hi": float(hi), "count": int(sel.sum())})

    consistent = None
    deficit = protein_deficits(genomes, genomes if len(genomes) >= 20 else catalog.genomes, assembly_kb(ctx))
    called = np.zeros(n_all, dtype=bool)
    called[[g.index for g in genomes if g.bin_id not in deficit]] = True
    every = int((fs.carriers(called) == int(called.sum())).sum()) if called.any() else 0
    testable = [g for g in genomes if g.completeness is not None and g.bin_id not in deficit]
    if len(testable) >= 3:
        known = np.zeros(n_all, dtype=bool)
        known[[g.index for g in testable]] = True
        ck = fs.carriers(known)
        test = completeness_test(testable, ck)
        tested = np.flatnonzero(p_in >= 0.5)
        q = bh_qvalues(test["p"][tested], len(tested))
        consistent = int((q >= ALPHA).sum())

    def rows(mask, order):
        idx = np.flatnonzero(mask)
        idx = idx[np.argsort(order[idx], kind="stable")][:top]
        return [{"id": fs.ids[k], "label": fs.labels[k], "name": fs.names[k], "in": float(p_in[k]),
                 "out": float(p_out[k]), "carriers": int(c_in[k])} for k in idx]

    selectable = fs.selectable
    signature = rows(selectable & (p_in >= SIGNATURE_IN) & (p_out <= SIGNATURE_OUT), -(p_in - p_out)) if n_out else []
    absent = rows(selectable & (p_out >= SIGNATURE_IN) & (p_in <= SIGNATURE_OUT), -(p_out - p_in)) if n_out else []
    comp = [g.completeness for g in genomes if g.completeness is not None]
    return {
        "kind": fs.kind, "genomes": n_in, "rest": n_out, "observed": int(present.sum()), "spectrum": spectrum,
        "every": every, "consistent": consistent, "tested_genomes": len(testable), "deficit": len(deficit),
        "rare": int((present & (p_in <= 0.1)).sum()), "signature": signature, "absent": absent,
        "median_completeness": float(np.median(comp)) if comp else None,
    }


def register(templates, ctx) -> None:
    """Expose ``clade_content(genomes, rank, name)`` to templates, cached per clade."""
    cache: dict[tuple[str, str], Any] = {}

    def for_template(genomes, rank, name):
        sets = getattr(ctx, "feature_sets", None)
        if sets is None or not rank:
            return None
        key = (rank, name)
        if key not in cache:
            try:
                cache[key] = clade_content(ctx, sets, genomes)
            except LookupError:
                cache[key] = None
        return cache[key]

    templates.env.globals["clade_content"] = for_template
