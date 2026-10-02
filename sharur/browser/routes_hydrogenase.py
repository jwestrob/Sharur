"""Hydrogenases: the nearest-HydDB-reference classifier's calls, each with its evidence layers.

Calls come from ``hydrogenase_classifications`` (``scripts/classify_hydrogenases.py``). The page keeps
apart the layers ``.claude/skills/hydrogenase.md`` names: a HydDB HMM hit, the subgroup of the nearest
HydDB reference, the catalytic-domain check, KOfam corroboration, and the gene neighborhood
(``sharur.hydrogenase.neighborhood``). Datasets that predate the classifier show the subgroup labels
they store and the command that classifies them.
"""

from __future__ import annotations

import re
from collections import Counter, defaultdict
from typing import Any

from fastapi import HTTPException, Query, Request
from fastapi.responses import HTMLResponse

from sharur.browser.catalog import RANKS
from sharur.hydrogenase.neighborhood import VERDICT_LABELS, VERDICTS, neighborhood_contexts

CLASS_LABELS = {"NiFe": "[NiFe]", "FeFe": "[FeFe]", "Fe": "[Fe]"}
GROUPS = (("NiFe", "1"), ("NiFe", "2"), ("NiFe", "3"), ("NiFe", "4"), ("FeFe", "A"), ("FeFe", "B"), ("FeFe", "C"),
          ("Fe", ""))
EVIDENCE = {
    "curated": "catalytic domain observed, or hydrogenase genes nearby",
    "assigned": "every nearest-reference assignment",
}
INTERPRETATION = {
    "characterized": "Søndergaard et al. 2016 Table 1 characterizes this subgroup",
    "putative": "Table 1 role is putative",
    "unresolved": "Table 1 leaves the role unresolved",
    "unverified": "label present in the installed HydDB reference; role unverified",
}
KO_SUPPORT = {
    "subgroup": "a KOfam hit captures this subgroup's HydDB references (≥80%)",
    "compatible": "a KOfam hit captures some references of this subgroup",
    "group": "KOfam hits capture this group, other subgroups",
    "conflict": "KOfam hits capture other groups",
    "none": "no KOfam hit associated with HydDB references",
}


def group_of(hyd_class: str | None, subgroup: str | None) -> str:
    """``Group_4g`` -> ``4``; ``Group_A3`` -> ``A``; [Fe] has no groups."""
    if hyd_class == "Fe":
        return ""
    m = re.match(r"Group_([0-9A-Z])", subgroup or "")
    return m.group(1) if m else ""


def group_label(hyd_class: str, group: str) -> str:
    return f"{CLASS_LABELS.get(hyd_class, hyd_class)} Group {group}" if group else CLASS_LABELS.get(hyd_class, hyd_class)


def _subgroup_key(row: dict[str, Any]) -> tuple:
    order = [c for c, _ in GROUPS]
    sub = row["subgroup"] or ""
    m = re.match(r"Group_([0-9]+)([a-z]*)", sub)
    natural = (int(m.group(1)), m.group(2)) if m else (99, sub)
    return (order.index(row["class"]) if row["class"] in order else 9, natural)


def register(app, ctx) -> None:
    catalog = ctx.catalog

    def _has_table() -> bool:
        with ctx.lock:
            return bool(ctx.store.execute(
                "SELECT COUNT(*) FROM information_schema.tables WHERE table_name = 'hydrogenase_classifications'")[0][0])

    def calls() -> list[dict[str, Any]] | None:
        """One record per HydDB-hit protein, or None when the dataset predates the classifier."""
        if hasattr(ctx, "_hydrogenase_calls"):
            return ctx._hydrogenase_calls
        out = None
        if _has_table():
            with ctx.lock:
                rows = ctx.store.execute(
                    """SELECT h.protein_id, p.bin_id, p.contig_id, p.start, p.end_coord, p.strand, h.outcome,
                              h.discovery_classes, h.reference_label, h.reference_class, h.reference_subgroup,
                              h.interpretation_status, h.reference_role, h.pident, h.has_nifese_hases, h.has_fe_hyd,
                              h.has_complex1, h.has_hmd, h.curation_status, h.curation_reason, h.ko_support,
                              h.reference_release, h.classifier_version
                       FROM hydrogenase_classifications h JOIN proteins p USING (protein_id)
                       ORDER BY p.bin_id, h.protein_id""")
                contexts = neighborhood_contexts(ctx.store, [r[0] for r in rows])
            out = [] if rows else None   # an empty table reads as unclassified
            for (pid, bin_id, contig, start, end, strand, outcome, discovery, label, hyd_class, subgroup, interp, role,
                 pident, nifese, fe_hyd, complex1, hmd, curation, reason, ko, release, version) in rows or []:
                g = catalog.by_bin.get(bin_id)
                nb = contexts.get(pid)
                cleared = curation == "domain_check_cleared"
                verdict = nb.verdict if nb else "no_position"
                domains = [d for d, present in (("NiFeSe_Hases", nifese), ("Fe_hyd", fe_hyd),
                                                ("Complex I superfamily", complex1), ("HMD", hmd)) if present]
                out.append({
                    "protein_id": pid, "bin_id": bin_id, "lineage": g.label if g else "", "contig_id": contig,
                    "start": start, "end": end, "strand": strand, "outcome": outcome, "discovery": discovery,
                    "label": label, "class": hyd_class, "subgroup": subgroup, "group": group_of(hyd_class, subgroup),
                    "interpretation": interp, "role": role, "pident": pident, "domains": domains,
                    "cleared": cleared, "curation_reason": reason, "ko_support": ko or "none",
                    "verdict": verdict, "hydrogenase_markers": sorted(nb.hydrogenase) if nb else [],
                    "complex_i_markers": sorted(nb.complex_i) if nb else [],
                    "hydrogenase_genes": nb.hydrogenase_genes if nb else set(),
                    "complex_i_genes": nb.complex_i_genes if nb else set(),
                    "curated": outcome == "assigned" and (cleared or verdict == "supported"),
                    "release": release, "version": version})
        ctx._hydrogenase_calls = out
        return out

    def by_protein() -> dict[str, dict[str, Any]]:
        if not hasattr(ctx, "_hydrogenase_by_protein"):
            ctx._hydrogenase_by_protein = {c["protein_id"]: c for c in calls() or []}
        return ctx._hydrogenase_by_protein

    def stored_labels() -> dict[str, Any]:
        """Datasets without the classifier table: their ``hyddb_subgroup`` labels and raw HydDB hits."""
        with ctx.lock:
            labels = ctx.store.execute(
                """SELECT a.accession, COUNT(DISTINCT a.protein_id), COUNT(DISTINCT p.bin_id)
                   FROM annotations a JOIN proteins p USING (protein_id) WHERE LOWER(a.source) = 'hyddb_subgroup'
                   GROUP BY 1 ORDER BY 2 DESC, 1""")
            hits = ctx.store.execute(
                """SELECT COALESCE(NULLIF(a.name, ''), a.accession), COUNT(DISTINCT a.protein_id), COUNT(DISTINCT p.bin_id)
                   FROM annotations a JOIN proteins p USING (protein_id) WHERE LOWER(a.source) = 'hyddb'
                   GROUP BY 1 ORDER BY 2 DESC, 1""")
        return {"labels": labels, "hits": hits}

    def summary() -> dict[str, Any] | None:
        """Evidence ladder, subgroup table and per-group totals (None without the classifier table)."""
        rows = calls()
        if rows is None:
            return None
        n = len(catalog.genomes) or 1
        assigned = [c for c in rows if c["outcome"] == "assigned"]
        flagged = [c for c in assigned if not c["cleared"]]
        subgroups: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
        for c in assigned:
            subgroups[(c["class"], c["subgroup"])].append(c)
        table = []
        for (hyd_class, subgroup), members in subgroups.items():
            first = members[0]
            table.append({
                "class": hyd_class, "subgroup": subgroup, "label": first["label"], "group": first["group"],
                "role": first["role"], "interpretation": first["interpretation"], "calls": len(members),
                "cleared": sum(c["cleared"] for c in members),
                "flagged": Counter(c["verdict"] for c in members if not c["cleared"]),
                "ko_subgroup": sum(c["ko_support"] == "subgroup" for c in members),
                "genomes": len({c["bin_id"] for c in members}),
                "curated_genomes": len({c["bin_id"] for c in members if c["curated"]}),
                "share": len({c["bin_id"] for c in members if c["curated"]}) / n})
        table.sort(key=_subgroup_key)
        groups = []
        for hyd_class, group in GROUPS:
            members = [c for c in assigned if c["class"] == hyd_class and c["group"] == group]
            if members:
                curated = {c["bin_id"] for c in members if c["curated"]}
                groups.append({"class": hyd_class, "group": group, "label": group_label(hyd_class, group),
                               "calls": len(members), "curated": sum(c["curated"] for c in members),
                               "genomes": len(curated), "share": len(curated) / n})
        return {
            "proteins": len(rows), "assigned": len(assigned), "conflicts": sum(c["outcome"] == "class_conflict" for c in rows),
            "other_outcomes": Counter(c["outcome"] for c in rows if c["outcome"] not in ("assigned", "class_conflict")),
            "cleared": len(assigned) - len(flagged), "flagged": len(flagged),
            "flagged_verdicts": Counter(c["verdict"] for c in flagged),
            "cleared_verdicts": Counter(c["verdict"] for c in assigned if c["cleared"]),
            "ko": Counter(c["ko_support"] for c in assigned),
            "curated": sum(c["curated"] for c in assigned),
            "genomes": len({c["bin_id"] for c in assigned if c["curated"]}), "table": table, "groups": groups,
            "release": rows[0]["release"] if rows else None, "version": rows[0]["version"] if rows else None}

    def distribution(rank: str, evidence: str) -> dict[str, Any]:
        """Share of each clade's genomes with a call per group."""
        rows = [c for c in calls() or [] if c["outcome"] == "assigned" and (evidence == "assigned" or c["curated"])]
        carriers: dict[tuple[str, str], set[int]] = defaultdict(set)
        for c in rows:
            g = catalog.by_bin.get(c["bin_id"])
            if g is not None:
                carriers[(c["class"], c["group"])].add(g.index)
        columns = [(cls, grp) for cls, grp in GROUPS if carriers.get((cls, grp))]
        totals = Counter(g.taxonomy.get(rank) or "Unclassified" for g in catalog.genomes)
        clade_of = {g.index: g.taxonomy.get(rank) or "Unclassified" for g in catalog.genomes}
        counts = {col: Counter(clade_of[i] for i in carriers[col]) for col in columns}
        clades = [t for t, k in sorted(totals.items(), key=lambda kv: (-kv[1], kv[0])) if k >= 3][:40]
        return {"columns": [{"class": cls, "group": grp, "label": group_label(cls, grp),
                             "short": f"{CLASS_LABELS.get(cls, cls)} {grp}".strip()} for cls, grp in columns],
                "rows": [{"taxon": t, "genomes": totals[t],
                          "cells": [{"n": counts[col][t], "share": counts[col][t] / totals[t]} for col in columns]}
                         for t in clades]}

    ctx.hydrogenase_calls = calls
    ctx.hydrogenase_summary = summary
    ctx.hydrogenase_call = lambda pid: by_protein().get(pid)
    ctx.hydrogenases_in = lambda bin_id: [c for c in calls() or [] if c["bin_id"] == bin_id and c["outcome"] == "assigned"]
    ctx.templates.env.globals["hydrogenase_call"] = ctx.hydrogenase_call
    ctx.templates.env.globals["VERDICT_LABELS"] = VERDICT_LABELS

    @app.get("/hydrogenases", response_class=HTMLResponse)
    def hydrogenases(request: Request, rank: str = Query("class"), evidence: str = Query("curated")):
        rank = rank if rank in RANKS else "class"
        evidence = evidence if evidence in EVIDENCE else "curated"
        s = summary()
        return ctx.render(request, "hydrogenases.html", "systems", s=s, rank=rank, ranks=RANKS[1:],
                          evidence=evidence, evidence_labels=EVIDENCE,
                          dist=distribution(rank, evidence) if s else None,
                          stored=None if s else stored_labels(), verdicts=VERDICTS, verdict_labels=VERDICT_LABELS,
                          interpretation=INTERPRETATION, ko_labels=KO_SUPPORT, class_labels=CLASS_LABELS)

    @app.get("/hydrogenases/calls", response_class=HTMLResponse)
    def hydrogenase_calls(request: Request, hyd_class: str = Query("", alias="class"), subgroup: str = Query(""),
                          group: str = Query(""), evidence: str = Query(""), verdict: str = Query(""),
                          genome: str = Query(""), clade: str = Query(""), rank: str = Query("class"),
                          page: int = Query(1, ge=1), flank: int = Query(6000, ge=0, le=30000)):
        rows = calls()
        if rows is None:
            raise HTTPException(404, "This dataset has no hydrogenase classifier calls")
        selected = [c for c in rows if c["outcome"] == "assigned"
                    and (not hyd_class or c["class"] == hyd_class) and (not subgroup or c["subgroup"] == subgroup)
                    and (not group or c["group"] == group) and (not genome or c["bin_id"] == genome)
                    and (not verdict or c["verdict"] == verdict)
                    and (evidence != "curated" or c["curated"])
                    and (evidence != "cleared" or c["cleared"]) and (evidence != "flagged" or not c["cleared"])]
        if clade:
            r = rank if rank in RANKS else "class"
            selected = [c for c in selected if c["bin_id"] in catalog.by_bin
                        and (catalog.by_bin[c["bin_id"]].taxonomy.get(r) or "Unclassified") == clade]
        per_page = 25
        pages = max(1, (len(selected) + per_page - 1) // per_page)
        page = min(page, pages)
        shown = selected[(page - 1) * per_page: page * per_page]
        placed = {c["protein_id"] for c in shown if c["verdict"] != "no_position"}
        span = max((c["end"] - c["start"] + 2 * flank for c in shown if c["protein_id"] in placed), default=1)
        drawn = []
        for c in shown:
            svg = None
            if c["protein_id"] in placed:
                center = (c["start"] + c["end"]) // 2
                lo, hi = center - span // 2, center + span // 2
                genes = ctx.genes_near(c["contig_id"], lo, hi)
                for gene in genes:
                    pid = gene["protein_id"]
                    if pid == c["protein_id"]:
                        gene["cls"], gene["cas_name"] = "gene hyd-focal", (c["subgroup"] or "").replace("Group_", "")
                    elif pid in c["hydrogenase_genes"]:
                        gene["cls"] = "gene hyd-marker"
                    elif pid in c["complex_i_genes"]:
                        gene["cls"] = "gene c1-marker"
                    else:
                        gene["cls"] = "gene" if gene["label"] else "gene dark"
                svg = ctx.locus_row_svg({"arrays": []}, genes, lo, hi, c["strand"] == "-", label="Hydrogenase locus")
            drawn.append((c, svg))
        title = (next((c["label"] for c in selected if c["subgroup"] == subgroup), subgroup) if subgroup
                 else group_label(hyd_class, group) if hyd_class else "All assignments")
        params = {"class": hyd_class, "subgroup": subgroup, "group": group, "evidence": evidence, "verdict": verdict,
                  "genome": genome, "clade": clade, "rank": rank, "flank": flank}
        return ctx.render(request, "hydrogenase_calls.html", "systems", rows=drawn, total=len(selected), page=page,
                          pages=pages, title=title, params=params, verdict_labels=VERDICT_LABELS, verdicts=VERDICTS,
                          ko_labels=KO_SUPPORT, interpretation=INTERPRETATION,
                          subgroups=sorted({(c["class"], c["subgroup"], c["label"]) for c in rows if c["outcome"] == "assigned"},
                                           key=lambda t: _subgroup_key({"class": t[0], "subgroup": t[1]})))
