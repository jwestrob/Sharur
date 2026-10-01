"""Scoped search: "<term> in <genome or clade>", e.g. "rubisco in EX4484-52".

The scope resolves to a genome (exact ID, or a unique part of one) or a taxon
at any rank. The term matches functional labels (by ID or name), annotation
names, accessions and descriptions (Pfam, KEGG, CAZy, ...) and VOG consensus
descriptions, wherever a word starts with it. Each hit lists what it matched.
"""

from __future__ import annotations

import re
from collections import Counter, defaultdict
from typing import Any

from sharur.browser.catalog import RANKS, UNCLASSIFIED

SCOPED = re.compile(r"^\s*(.+?)\s+in\s+(.+?)\s*$", re.I)


def resolve_scope(ctx, text: str) -> tuple[str, str, list[str]] | None:
    """(kind, label, bin_ids) for a genome or taxon named by ``text``."""
    catalog = ctx.catalog
    needle = text.strip()
    low = needle.lower()
    if needle in catalog.by_bin:
        return "genome", needle, [needle]
    # GTDB placeholder names recur across ranks (a class and a genus can share a name):
    # take the broadest clade by that name; the page links the narrower ones
    for rank in RANKS:
        members = [g.bin_id for g in catalog.genomes if (g.taxonomy.get(rank) or "").lower() == low]
        if members and low != UNCLASSIFIED.lower():
            name = catalog.by_bin[members[0]].taxonomy[rank]
            return rank, name, members
    partial = [g.bin_id for g in catalog.genomes if low in g.bin_id.lower()]
    if len(partial) == 1:
        return "genome", partial[0], partial
    if partial and len(partial) <= 50:
        return "genomes matching", needle, partial
    return None


def scoped_search(ctx, term: str, bins: list[str], limit: int = 2000) -> dict[str, Any]:
    """Proteins in ``bins`` matching ``term``, with what each matched."""
    term_l = term.strip().lower()
    # the term starts a word: "hydrogenase" finds hydrogenases, not dehydrogenases
    word = re.compile(r"(?<![a-z0-9])" + re.escape(term_l))
    labels = [p for p, d in ctx.predicates.items()
              if term_l == p or word.search(d.name.lower()) or word.search(p.replace("_", " "))]
    vogs = [v for v, info in ctx.catalog.vogs.items() if word.search((info.get("description") or "").lower())]
    boundary = "(^|[^A-Za-z0-9])" + re.escape(term.strip())
    matched: dict[str, dict[str, Any]] = defaultdict(lambda: {"labels": set(), "hits": set()})
    with ctx.lock:
        if labels:
            for pid, predicates in ctx.store.execute(
                    """SELECT pp.protein_id, pp.predicates FROM protein_predicates pp JOIN proteins p USING (protein_id)
                       WHERE p.bin_id IN (SELECT UNNEST(?::VARCHAR[])) AND list_has_any(pp.predicates, ?::VARCHAR[])""",
                    [bins, labels]):
                matched[pid]["labels"].update(set(predicates) & set(labels))
        pattern = f"%{term.strip()}%"
        for pid, source, accession, name, description in ctx.store.execute(
                """SELECT a.protein_id, LOWER(a.source), a.accession, a.name, a.description
                   FROM annotations a JOIN proteins p USING (protein_id)
                   WHERE p.bin_id IN (SELECT UNNEST(?::VARCHAR[]))
                     AND (((a.name ILIKE ? OR a.description ILIKE ? OR a.accession ILIKE ?)
                           AND regexp_matches(concat_ws(' ', a.name, a.description, a.accession), ?, 'i'))
                          OR a.accession IN (SELECT UNNEST(?::VARCHAR[])))""",
                [bins, pattern, pattern, pattern, boundary, vogs]):
            text = ctx.describe_hit(name or accession) if source in ("vogdb", "vog") else (name or accession)
            matched[pid]["hits"].add(f"{source}: {text}" + (f" — {description}" if description and source not in ("vogdb", "vog") else ""))
        ids = list(matched)[:limit]
        info = {pid: (b, n, contig, start) for pid, b, n, contig, start in ctx.store.execute(
            "SELECT protein_id, bin_id, sequence_length, contig_id, start FROM proteins "
            "WHERE protein_id IN (SELECT UNNEST(?::VARCHAR[]))", [ids])} if ids else {}
        rest = [pid for pid in matched if pid not in info]
        bin_of = {pid: v[0] for pid, v in info.items()}
        if rest:
            bin_of.update(dict(ctx.store.execute(
                "SELECT protein_id, bin_id FROM proteins WHERE protein_id IN (SELECT UNNEST(?::VARCHAR[]))", [rest])))
        # the KOs and Pfam families the matched proteins carry, however they matched
        families: dict[tuple[str, str], set[str]] = defaultdict(set)
        sample = list(matched)[:5000]
        if sample:
            for pid, source, accession in ctx.store.execute(
                    """SELECT protein_id, LOWER(source), accession FROM annotations
                       WHERE protein_id IN (SELECT UNNEST(?::VARCHAR[]))
                         AND LOWER(source) IN ('kofam', 'kegg', 'pfam')""", [sample]):
                acc = (accession or "").split(".")[0]
                if source in ("kofam", "kegg") and re.match(r"^K\d{5}$", acc):
                    families[("ko", acc)].add(pid)
                elif source == "pfam" and re.match(r"^PF\d{5}$", acc):
                    families[("pfam", acc)].add(pid)
    per_genome = Counter(bin_of[pid] for pid in matched if bin_of.get(pid))
    rows = []
    for pid in ids:
        bin_id, length, contig, start = info.get(pid, (None, None, None, None))
        g = ctx.catalog.by_bin.get(bin_id)
        rows.append({"protein_id": pid, "bin_id": bin_id, "lineage": g.label if g else "", "length": length,
                     "contig": contig, "start": start,
                     "labels": sorted(ctx.predicates[p].name for p in matched[pid]["labels"] if p in ctx.predicates),
                     "hits": sorted(matched[pid]["hits"])[:4]})
    rows.sort(key=lambda r: (r["bin_id"] or "", r["contig"] or "", r["start"] or 0))
    top_families = sorted(families.items(), key=lambda kv: (-len(kv[1]), kv[0][1]))
    return {"rows": rows, "total": len(matched), "labels": [ctx.predicates[p].name for p in labels[:8]],
            "genomes": len(per_genome), "per_genome": per_genome,
            "families": {kind: [(acc, len(pids)) for (k, acc), pids in top_families if k == kind][:20]
                         for kind in ("ko", "pfam")}}


def same_name_ranks(ctx, name: str) -> list[tuple[str, int]]:
    """(rank, genome count) for every rank where a taxon carries ``name``."""
    low = name.lower()
    out = []
    for rank in RANKS:
        n = sum(1 for g in ctx.catalog.genomes if (g.taxonomy.get(rank) or "").lower() == low)
        if n:
            out.append((rank, n))
    return out
