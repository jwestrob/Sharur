"""The overview page: the dataset at a glance, from the in-memory catalog.

Everything here reads the catalog (and, when present, the Discover feeds, the
CRISPR-Cas summary and the curation notes), so the page renders in
milliseconds once the catalog has loaded and fills in as background summaries
finish.
"""

from __future__ import annotations

from collections import Counter
from statistics import median
from typing import Any

from sharur.browser.catalog import UNCLASSIFIED

# Genome quality tiers (MIMAG): high quality > 90% complete and < 5% contamination;
# medium quality >= 50% complete and < 10% contamination.
HQ = (90.0, 5.0)
MQ = (50.0, 10.0)
SCATTER_GROUPS = 3   # the scatter colours this many phyla; the rest read as "Other"


def quality(catalog) -> dict[str, Any]:
    """Completeness/contamination points, MIMAG tiers and the phyla that colour the scatter."""
    scored = [g for g in catalog.genomes if g.completeness is not None]
    tiers = Counter()
    for g in scored:
        contamination = g.contamination if g.contamination is not None else 0.0
        if g.completeness > HQ[0] and contamination < HQ[1]:
            tiers["high"] += 1
        elif g.completeness >= MQ[0] and contamination < MQ[1]:
            tiers["medium"] += 1
        else:
            tiers["low"] += 1
    phyla = Counter(g.taxonomy.get("phylum") or UNCLASSIFIED for g in scored)
    groups = [p for p, _ in phyla.most_common() if p != UNCLASSIFIED][:SCATTER_GROUPS]
    slot = {p: i for i, p in enumerate(groups)}
    points = [{"id": g.bin_id, "x": float(g.completeness),
               "y": float(g.contamination) if g.contamination is not None else 0.0,
               "group": slot.get(g.taxonomy.get("phylum"), SCATTER_GROUPS), "label": g.label}
              for g in scored]
    legend = [(p, phyla[p]) for p in groups]
    other = len(scored) - sum(n for _, n in legend)
    return {"points": points, "tiers": tiers, "scored": len(scored), "total": len(catalog.genomes),
            "legend": legend, "other": other,
            "median": median(g.completeness for g in scored) if scored else None}


def sizes(catalog) -> dict[str, list[float]]:
    return {"mb": [g.length / 1e6 for g in catalog.genomes if g.length],
            "proteins": [float(g.proteins) for g in catalog.genomes if g.proteins]}


def systems(catalog, cctyper_summary=None) -> dict[str, Any]:
    """Curated system calls by kind, the commonest types per kind, and CRISPR-Cas subtypes."""
    by_kind: dict[str, Counter] = {}
    genomes: dict[str, set[str]] = {}
    for s in catalog.systems:
        by_kind.setdefault(s["kind"], Counter())[s["type"]] += 1
        genomes.setdefault(s["kind"], set()).add(s["bin_id"])
    kinds = [{"kind": k, "calls": sum(c.values()), "types": len(c), "genomes": len(genomes[k]),
              "top": c.most_common(6)} for k, c in sorted(by_kind.items(), key=lambda kv: -sum(kv[1].values()))]
    crispr = None
    if cctyper_summary is not None:
        try:
            summary = cctyper_summary()
        except Exception:  # noqa: BLE001 - the panel is optional
            summary = None
        if summary and summary.get("subtypes"):
            crispr = {"calls": sum(r["count"] for r in summary["subtypes"]), "subtypes": summary["subtypes"][:8],
                      "n_subtypes": len(summary["subtypes"])}
    return {"kinds": kinds, "crispr": crispr}


def highlights(catalog) -> list[dict[str, Any]]:
    """A few standouts from the Discover feeds, once they are ready."""
    feeds = catalog.notable or {}
    out = []
    giants = feeds.get("giants") or []
    if giants:
        g = giants[0]
        out.append({"kind": "giant", "title": "Longest protein", "value": f"{g['length']:,} aa",
                    "protein_id": g["protein_id"], "bin_id": g["bin_id"], "domains": g.get("domains", []),
                    "length": g["length"], "partial": g.get("partial", False),
                    "detail": g.get("architecture") or "no placed Pfam domains", "href": "/discover/giants"})
    fusions = feeds.get("fusions") or []
    if fusions:
        f = fusions[0]
        out.append({"kind": "fusion", "title": "Clade-specific domain fusion",
                    "value": " + ".join(name for _, name in f["domains"]),
                    "detail": f"{f['genomes']:,} of {f['clade_size']:,} {f['clade']} genomes",
                    "protein_id": f["example"], "href": "/discover/fusions"})
    islands = feeds.get("islands") or []
    if islands:
        i = islands[0]
        out.append({"kind": "island", "title": "Longest unannotated stretch", "value": f"{i['genes']} genes",
                    "detail": f"{(i['end'] - i['start']) / 1000:,.1f} kb on {i['contig_id']}",
                    "strands": i["strands"], "contig_id": i["contig_id"], "start": i["start"], "end": i["end"],
                    "bin_id": i["bin_id"], "href": "/discover/islands"})
    clades = (feeds.get("giant_clades") or {}).get("rows") or []
    if clades:
        c = clades[0]
        rank = feeds["giant_clades"]["rank"]
        out.append({"kind": "clade", "title": "Richest in giant proteins", "value": c["clade"],
                    "detail": f"{c['per_genome']:.1f} proteins ≥ {feeds['giant_clades']['threshold']:,} aa per genome",
                    "rank": rank, "href": "/discover/clades"})
    return out


def recent_curation(ctx, limit: int = 6) -> list[dict[str, Any]]:
    notes = getattr(ctx, "notes", None)
    if notes is None:
        return []
    from sharur.browser.routes_curation import entity_label, entity_url  # noqa: PLC0415

    rows = []
    for r in notes.listing()[:limit]:
        try:
            label = entity_label(ctx, r["kind"], r["entity"])
        except Exception:  # noqa: BLE001 - a stale entity still lists by id
            label = r["entity"]
        try:
            url = entity_url(ctx, r["kind"], r["entity"])
        except Exception:  # noqa: BLE001
            url = "#"
        rows.append({**r, "url": url, "label": label})
    return rows


def length_basis(ctx) -> str:
    """How contig lengths were stored: 'assembly', 'gene_span' or 'mixed' (cached)."""
    cached = getattr(ctx, "_length_basis", None)
    if cached is None:
        cached = "assembly"
        try:
            with ctx.lock:
                cols = {r[0] for r in ctx.store.execute(
                    "SELECT column_name FROM information_schema.columns WHERE table_name = 'contigs'")}
                if "length_source" in cols:
                    kinds = {r[0] for r in ctx.store.execute("SELECT DISTINCT length_source FROM contigs")}
                    kinds.discard(None)
                    cached = "gene_span" if kinds == {"gene_span"} else "mixed" if "gene_span" in kinds else "assembly"
        except Exception:  # noqa: BLE001 - the label falls back to "assembled"
            cached = "assembly"
        ctx._length_basis = cached
    return cached


SEQUENCE_LABELS = {"assembly": "assembled sequence", "gene_span": "sequence spanned by gene calls",
                   "mixed": "contig sequence"}
LOCUS_LABELS = {"crispr": "CRISPR arrays", "prophage": "prophage regions", "island": "genomic islands"}


def build(catalog, ctx) -> dict[str, Any]:
    basis = length_basis(ctx)
    return {"basis": basis, "sequence_label": SEQUENCE_LABELS[basis], "locus_labels": LOCUS_LABELS,
            "quality": quality(catalog), "sizes": sizes(catalog),
            "systems": systems(catalog, getattr(ctx, "cctyper_summary", None)),
            "highlights": highlights(catalog), "recent": recent_curation(ctx)}
