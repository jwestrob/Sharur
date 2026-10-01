"""Synteny views over a dataset's synteny sidecar: conserved gene-order clusters.

Registered only when the dataset has a readable synteny sidecar whose genes
resolve to the dataset's proteins; without one, nothing here exists (no
routes, no panels, no links). A sidecar built on an earlier dataset version is
shown with a note, provided its sampled gene IDs still resolve.

- Protein pages load a panel listing the clusters the protein belongs to.
- ``/synteny/cluster/{key}`` draws a cluster's loci as stacked, aligned rows.
- ``/synteny`` summarizes the run and lists the largest clusters, those spread
  across the most lineages, and those carried by few genomes in distant
  lineages, plus search by protein or genome.

Per-request queries use the sidecar's (run_id, protein_id) and
(run_id, cluster_key) indexes; summaries are computed once in the background.
"""

from __future__ import annotations

import json
import logging
import re
import statistics
import threading
from collections import Counter
from pathlib import Path
from typing import Any

from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.responses import HTMLResponse, RedirectResponse

from sharur.browser.catalog import RANKS, UNCLASSIFIED
from sharur.browser.routes_loci import _family_key, build_rows, render_stack

logger = logging.getLogger(__name__)

PER_PAGE = 25
MIN_RESOLVED = 0.98      # share of sampled sidecar genes that must exist in the dataset
RESOLVE_SAMPLE = 2000
LIST_LIMIT = 300
PATCHY_MAX_GENOMES = 10  # "few genomes" for the distant-lineage list
PATCHY_MIN_GENES = 4     # median genes per locus, so pairs of genes don't dominate


def lca(genomes) -> tuple[str, str] | None:
    """Deepest rank every genome shares."""
    found = None
    for rank in RANKS:
        values = {g.taxonomy.get(rank) for g in genomes}
        if len(values) == 1:
            value = next(iter(values))
            if value and value != UNCLASSIFIED:
                found = (rank, value)
                continue
        break
    return found


def spread(catalog, genome_ids) -> dict[str, Any]:
    """Taxonomic spread of a set of genomes: LCA, distinct clades per rank, clades below the LCA."""
    genomes = [catalog.by_bin[b] for b in dict.fromkeys(genome_ids) if b in catalog.by_bin]
    common = lca(genomes)
    distinct = {r: len({g.taxonomy.get(r) for g in genomes} - {None, UNCLASSIFIED}) for r in RANKS}
    below = None
    start = RANKS.index(common[0]) + 1 if common else 0
    for rank in RANKS[start:]:
        if distinct[rank] > 1:
            below = rank
            break
    clades = Counter(g.taxonomy.get(below) or UNCLASSIFIED for g in genomes).most_common() if below else []
    return {"genomes": len(genomes), "lca": common, "distinct": distinct, "rank": below, "clades": clades}


def _flatten(prefix: str, value: Any, out: list[tuple[str, str]]) -> None:
    """Run parameters as dotted key/value pairs, keeping the method's name out of what is shown."""
    if isinstance(value, dict):
        for k, v in value.items():
            key = re.sub(r"(?i)elsa[_ ]?", "", str(k)) or str(k)
            _flatten(f"{prefix}.{key}" if prefix else key, v, out)
        return
    text = (", ".join(map(str, value[:12])) + ("…" if len(value) > 12 else "")) if isinstance(value, list) else str(value)
    if "elsa" not in text.lower():
        out.append((prefix, text))


class SyntenyView:
    """One open synteny sidecar, its active run, and cached summaries."""

    def __init__(self, ctx, path: Path):
        from sharur.synteny import (  # noqa: PLC0415
            SyntenyStore,
            inspect_synteny_dataset_identity,
            inspect_synteny_sidecar,
        )

        self.ctx, self.path = ctx, path
        self.inspection = inspect_synteny_sidecar(path)
        if self.inspection.state != "available":
            raise RuntimeError(self.inspection.error or self.inspection.state)
        self.identity = inspect_synteny_dataset_identity(ctx.db_path, self.inspection)
        self.store = SyntenyStore(path, core_db_path=ctx.db_path, allow_stale=True)
        self.conn = self.store.connection
        self.lock = threading.Lock()
        self.run_id = self.inspection.active_run_id
        self.run = next(r for r in self.store.runs() if r["run_id"] == self.run_id)
        self.resolved = self._resolve_share()
        self.summary: dict[str, Any] | None = None
        self.summary_error: str | None = None

    # ---- state -------------------------------------------------------------------------------------------

    def _resolve_share(self) -> float:
        with self.lock:
            ids = [r[0] for r in self.conn.execute(
                f"SELECT protein_id FROM (SELECT protein_id FROM elsa_genes WHERE run_id = ?) "
                f"USING SAMPLE {RESOLVE_SAMPLE} ROWS", [self.run_id]).fetchall()]
        if not ids:
            return 0.0
        with self.ctx.lock:
            found = self.ctx.store.execute(
                "SELECT COUNT(*) FROM proteins WHERE protein_id IN (SELECT UNNEST(?::VARCHAR[]))", [ids])[0][0]
        return found / len(ids)

    @property
    def usable(self) -> bool:
        return self.identity.compatible or self.resolved >= MIN_RESOLVED

    @property
    def stale(self) -> bool:
        return not self.identity.compatible

    def run_info(self) -> dict[str, Any]:
        params: list[tuple[str, str]] = []
        try:
            _flatten("", json.loads(self.run.get("parameters_json") or "{}"), params)
        except ValueError:
            pass
        return {**{k: self.run.get(k) for k in ("run_id", "run_label", "created_at", "elsa_version", "status",
                                                  "gene_count", "block_count", "cluster_count", "singleton_count",
                                                  "locus_count", "member_count")},
                "version": self.run.get("elsa_version"), "parameters": params, "stale": self.stale,
                "resolved": self.resolved}

    # ---- per-request lookups --------------------------------------------------------------------------------

    def clusters_for_protein(self, protein_id: str, limit: int = 40) -> list[dict[str, Any]]:
        with self.lock:
            members = self.conn.execute(
                "SELECT cluster_key, member_role FROM elsa_cluster_members WHERE run_id = ? AND protein_id = ?",
                [self.run_id, protein_id]).fetchall()
            keys = list(dict.fromkeys(k for k, _ in members))
            if not keys:
                return []
            rows = self.conn.execute(
                """SELECT cluster_key, source_cluster_id, cluster_kind, size, genome_support, locus_count,
                          member_count, mean_chain_score
                   FROM elsa_clusters WHERE run_id = ? AND cluster_key IN (SELECT UNNEST(?::VARCHAR[]))""",
                [self.run_id, keys]).fetchall()
        roles: dict[str, set[str]] = {}
        for key, role in members:
            roles.setdefault(key, set()).add(role)
        out = [{"cluster_key": k, "source_cluster_id": sid, "kind": kind, "blocks": size, "genomes": support,
                "loci": loci, "members": n_members, "chain_score": score,
                "role": "anchor" if "anchor" in roles.get(k, ()) else "context"}
               for k, sid, kind, size, support, loci, n_members, score in rows]
        out.sort(key=lambda c: (c["kind"] != "cluster", -(c["genomes"] or 0), -(c["blocks"] or 0), c["cluster_key"]))
        out = out[:limit]
        for c in out[:12]:
            if (c["genomes"] or 0) >= 2:
                c["spread"] = spread(self.ctx.catalog, self.genomes_of(c["cluster_key"]))
        return out

    def genomes_of(self, cluster_key: str) -> list[str]:
        with self.lock:
            return [r[0] for r in self.conn.execute(
                "SELECT DISTINCT genome_id FROM elsa_cluster_loci WHERE run_id = ? AND cluster_key = ?",
                [self.run_id, cluster_key]).fetchall()]

    def cluster(self, key: str) -> dict[str, Any] | None:
        with self.lock:
            cursor = self.conn.execute(
                "SELECT * FROM elsa_clusters WHERE run_id = ? AND cluster_key = ?", [self.run_id, key])
            row = cursor.fetchone()
            if row is None and key.isdigit():
                cursor = self.conn.execute(
                    "SELECT * FROM elsa_clusters WHERE run_id = ? AND source_cluster_id = ?", [self.run_id, int(key)])
                row = cursor.fetchone()
            if row is None:
                return None
            summary = dict(zip([d[0] for d in cursor.description], row))
            cursor = self.conn.execute(
                """SELECT locus_key, genome_id, contig_id, start_position_index, end_position_index, start_bp, end_bp,
                          n_genes, block_support
                   FROM elsa_cluster_loci WHERE run_id = ? AND cluster_key = ?""",
                [self.run_id, summary["cluster_key"]])
            columns = [d[0] for d in cursor.description]
            summary["loci"] = [dict(zip(columns, r)) for r in cursor.fetchall()]
        return summary

    def members(self, cluster_key: str, locus_keys: list[str]) -> list[tuple[str, str, str, int]]:
        with self.lock:
            return self.conn.execute(
                """SELECT locus_key, protein_id, member_role, position_index FROM elsa_cluster_members
                   WHERE run_id = ? AND cluster_key = ? AND locus_key IN (SELECT UNNEST(?::VARCHAR[]))
                   ORDER BY locus_key, position_index""", [self.run_id, cluster_key, locus_keys]).fetchall()

    def clusters_in_genome(self, genome_id: str, limit: int = 200) -> list[dict[str, Any]]:
        with self.lock:
            rows = self.conn.execute(
                """WITH g AS (SELECT DISTINCT cluster_key FROM elsa_cluster_loci WHERE run_id = ? AND genome_id = ?)
                   SELECT c.cluster_key, c.source_cluster_id, c.cluster_kind, c.size, c.genome_support, c.locus_count
                   FROM elsa_clusters c JOIN g USING (cluster_key) WHERE c.run_id = ?
                   ORDER BY c.cluster_kind <> 'cluster', c.genome_support DESC, c.size DESC, c.cluster_key LIMIT ?""",
                [self.run_id, genome_id, self.run_id, limit]).fetchall()
        return [{"cluster_key": k, "source_cluster_id": sid, "kind": kind, "blocks": size, "genomes": support,
                 "loci": loci} for k, sid, kind, size, support, loci in rows]

    # ---- background summaries ---------------------------------------------------------------------------

    def summarize(self) -> None:
        try:
            self.summary = self._summarize()
        except Exception as exc:  # noqa: BLE001 - the overview shows the reason
            logger.exception("synteny summary failed")
            self.summary_error = f"{type(exc).__name__}: {exc}"

    def _summarize(self) -> dict[str, Any]:
        import pandas as pd  # noqa: PLC0415

        catalog = self.ctx.catalog
        tax = pd.DataFrame([{"bin_id": g.bin_id, **{r: (None if g.taxonomy.get(r) in (None, UNCLASSIFIED)
                                                        else g.taxonomy.get(r)) for r in RANKS}}
                            for g in catalog.genomes])
        with self.lock:
            counts = dict(self.conn.execute(
                "SELECT cluster_kind, COUNT(*) FROM elsa_clusters WHERE run_id = ? GROUP BY 1", [self.run_id]).fetchall())
            largest = self.conn.execute(
                """SELECT cluster_key, source_cluster_id, size, genome_support, locus_count, member_count
                   FROM elsa_clusters WHERE run_id = ? AND cluster_kind = 'cluster'
                   ORDER BY genome_support DESC, size DESC, cluster_key LIMIT ?""", [self.run_id, LIST_LIMIT]).fetchall()
            self.conn.register("sharur_synteny_tax", tax)
            try:
                lineages = self.conn.execute(
                    """WITH l AS (SELECT cluster_key, genome_id, MAX(n_genes) AS n_genes FROM elsa_cluster_loci
                                  WHERE run_id = ? GROUP BY 1, 2),
                       c AS (SELECT cluster_key, source_cluster_id, size, genome_support, locus_count FROM elsa_clusters
                             WHERE run_id = ? AND cluster_kind = 'cluster' AND genome_support >= 3)
                       SELECT c.cluster_key, c.source_cluster_id, c.size, c.genome_support, c.locus_count,
                              COUNT(DISTINCT t.phylum), COUNT(DISTINCT t.class), COUNT(DISTINCT t."order"),
                              MEDIAN(l.n_genes), LIST(DISTINCT l.genome_id)
                       FROM c JOIN l USING (cluster_key) JOIN sharur_synteny_tax t ON t.bin_id = l.genome_id
                       GROUP BY ALL HAVING COUNT(DISTINCT t.class) >= 2""", [self.run_id, self.run_id]).fetchall()
            finally:
                self.conn.unregister("sharur_synteny_tax")

        def entry(key, sid, size, support, loci, phyla=None, classes=None, orders=None, genes=None, genomes=None):
            row = {"cluster_key": key, "source_cluster_id": sid, "blocks": size, "genomes": support, "loci": loci,
                   "phyla": phyla, "classes": classes, "orders": orders, "median_genes": genes}
            if genomes is not None:
                row["spread"] = spread(catalog, genomes)
            return row

        lineage_rows = [entry(*r) for r in lineages]
        widespread = sorted(lineage_rows, key=lambda r: (-r["phyla"], -r["classes"], -r["genomes"], r["cluster_key"]))
        patchy = sorted((r for r in lineage_rows if r["genomes"] <= PATCHY_MAX_GENOMES
                         and (r["median_genes"] or 0) >= PATCHY_MIN_GENES),
                        key=lambda r: (-r["phyla"], -r["classes"], -(r["median_genes"] or 0), r["genomes"],
                                       r["cluster_key"]))
        largest_rows = [entry(*r[:5]) for r in largest]
        for r in largest_rows[:60]:
            r["spread"] = spread(catalog, self.genomes_of(r["cluster_key"]))
        return {"clusters": counts.get("cluster", 0), "pairs": counts.get("singleton", 0),
                "largest": largest_rows, "widespread": widespread[:LIST_LIMIT], "patchy": patchy[:LIST_LIMIT],
                "lineage_total": len(lineage_rows), "patchy_total": sum(1 for r in lineage_rows
                                                                         if r["genomes"] <= PATCHY_MAX_GENOMES
                                                                         and (r["median_genes"] or 0) >= PATCHY_MIN_GENES)}


def open_view(ctx) -> SyntenyView | None:
    """The dataset's synteny view, or None when it has no usable sidecar."""
    from sharur.synteny import discover_synteny_sidecar  # noqa: PLC0415

    path = discover_synteny_sidecar(ctx.db_path)
    if path is None:
        return None
    try:
        view = SyntenyView(ctx, path)
    except Exception:  # noqa: BLE001 - an unreadable sidecar leaves the browser as if it were absent
        logger.exception("synteny sidecar %s could not be opened", path)
        return None
    return view if view.usable else None


def register(app: FastAPI, ctx, templates) -> None:
    """Add synteny routes and template hooks when the dataset has a usable sidecar; otherwise add nothing."""
    templates.env.globals["synteny"] = None
    view = open_view(ctx)
    if view is None:
        return
    templates.env.globals["synteny"] = view
    ctx.synteny = view
    catalog, url = ctx.catalog, ctx.url
    threading.Thread(target=view.summarize, daemon=True).start()

    def patchy_rows() -> list[dict[str, Any]] | None:
        return view.summary["patchy"] if view.summary else None

    ctx.synteny_patchy = patchy_rows

    @app.get("/synteny", response_class=HTMLResponse)
    def synteny_overview(request: Request, q: str = Query("", max_length=300)):
        q = q.strip()
        genome_hits = None
        if q:
            with ctx.lock:
                is_protein = bool(ctx.store.execute("SELECT 1 FROM proteins WHERE protein_id = ?", [q]))
            if is_protein:
                return RedirectResponse(url("protein", q) + "#synteny", status_code=303)
            if q in catalog.by_bin:
                genome_hits = view.clusters_in_genome(q)
        return ctx.render(request, "synteny.html", "discover", run=view.run_info(), summary=view.summary,
                          summary_error=view.summary_error, q=q, genome_hits=genome_hits,
                          patchy_max=PATCHY_MAX_GENOMES, patchy_genes=PATCHY_MIN_GENES)

    @app.get("/synteny/protein/{protein_id:path}", response_class=HTMLResponse)
    def synteny_protein(request: Request, protein_id: str):
        return ctx.render(request, "_synteny_clusters.html", "taxa", clusters=view.clusters_for_protein(protein_id),
                          protein_id=protein_id)

    @app.get("/synteny/cluster/{cluster_key:path}", response_class=HTMLResponse)
    def synteny_cluster(request: Request, cluster_key: str, page: int = Query(1, ge=1),
                        flank: int = Query(2, ge=0, le=10), rank: str = Query(""), clade: str = Query("")):
        summary = view.cluster(cluster_key)
        if summary is None:
            raise HTTPException(404, "Syntenic cluster not found")
        key = summary["cluster_key"]
        loci = summary.pop("loci")
        taxa = spread(catalog, [l["genome_id"] for l in loci])

        def order(locus):
            g = catalog.by_bin.get(locus["genome_id"])
            lineage = tuple((g.taxonomy.get(r) or "~") for r in RANKS) if g else ("~",) * len(RANKS)
            return (*lineage, locus["genome_id"], locus["contig_id"] or "", locus["start_position_index"] or 0)

        loci.sort(key=order)
        if rank in RANKS and clade:
            loci = [l for l in loci if (catalog.by_bin.get(l["genome_id"]) and
                                        catalog.by_bin[l["genome_id"]].taxonomy.get(rank) == clade)]
        pages = max(1, (len(loci) + PER_PAGE - 1) // PER_PAGE)
        page = min(page, pages)
        shown = loci[(page - 1) * PER_PAGE: page * PER_PAGE]
        by_locus: dict[str, list[tuple[str, str, int]]] = {}
        for locus_key, pid, role, pos in view.members(key, [l["locus_key"] for l in shown]):
            by_locus.setdefault(locus_key, []).append((pid, role, pos))
        # orient every row on the same gene family: the commonest family among anchor genes on this page
        anchor_ids = [pid for rows in by_locus.values() for pid, role, _ in rows if role == "anchor"]
        with ctx.lock:
            best = dict(ctx.store.execute(
                """SELECT protein_id, ARG_MIN(COALESCE(NULLIF(name, ''), accession), COALESCE(evalue, 1))
                   FROM annotations WHERE protein_id IN (SELECT UNNEST(?::VARCHAR[]))
                   AND LOWER(source) NOT IN ('defensefinder_system', 'txsscan_system', 'hyddb_subgroup')
                   GROUP BY 1""", [anchor_ids])) if anchor_ids else {}
        family = {pid: _family_key(ctx.describe_hit(best.get(pid))) for pid in anchor_ids}
        reference = Counter(f for f in family.values() if f).most_common(1)
        reference = reference[0][0] if reference else None
        anchors = []
        for locus in shown:
            rows = by_locus.get(locus["locus_key"], [])
            if not rows:
                continue
            anchored = [r for r in rows if r[1] == "anchor"] or rows
            pick = next((r for r in anchored if reference and family.get(r[0]) == reference),
                        anchored[len(anchored) // 2])
            g = catalog.by_bin.get(locus["genome_id"])
            anchors.append({"anchor": pick[0], "members": {r[0] for r in rows}, "profiles": {},
                            "bin_id": locus["genome_id"], "title": f"{locus['n_genes']} genes",
                            "genome_label": g.label if g else "", "locus": locus})
        built = build_rows(ctx, anchors, flank)
        svgs, legend = render_stack(built, member_label="cluster gene")
        genes = [l["n_genes"] for l in loci if l["n_genes"]]
        return ctx.render(request, "synteny_cluster.html", "discover", c=summary, taxa=taxa, rows=list(zip(built, svgs)),
                          legend=legend, total=len(loci), page=page, pages=pages, flank=flank, rank=rank, clade=clade,
                          reference=reference, median_genes=statistics.median(genes) if genes else None,
                          run=view.run_info(), base=url("synteny", "cluster", str(summary["source_cluster_id"])))
