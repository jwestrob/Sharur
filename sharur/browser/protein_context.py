"""Context for one protein: its architecture elsewhere, its paralogs, its usual neighbours.

Everything here is observed: domain arrangements, shared KOs, genes that sit
nearby. Queries start from indexed accession lookups and touch at most a few
thousand proteins, so the panels fill in a fraction of a second after the page.

- :func:`family_of` picks the protein's family: its best KO, or else its
  best-scoring Pfam domain.
- :func:`same_architecture` counts proteins with exactly this resolved Pfam
  architecture: every protein carrying all of its domains is a candidate, and
  up to ``MAX_CANDIDATES`` of them (a stable sample) are resolved and compared,
  so very common domain sets give an estimate. Without Pfam domains it counts
  proteins sharing the best KO.
- :func:`paralogs` lists proteins in the same genome with the same KO or the
  same Pfam domains, flagging tandem copies.
- :func:`neighbourhood` samples up to ``max_genomes`` carriers of the family
  (one per genome) and ranks the families found within ``flank`` genes of
  them, with their usual side relative to the carrier's direction.
"""

from __future__ import annotations

import hashlib
import statistics
from collections import Counter, defaultdict
from typing import Any

from sharur.browser.catalog import RANKS, UNCLASSIFIED

MAX_CANDIDATES = 1000      # same-domain proteins resolved exactly; more are sampled
CARRIER_SAMPLE = 800       # carriers looked up when choosing one per genome
MAX_GENOMES = 150
FLANK = 5


def _rows(store, sql: str, params: list | None = None) -> list:
    result = store.execute(sql, params or [])
    return result.fetchall() if hasattr(result, "fetchall") else list(result)


def _in(column: str, ids: list[str]) -> tuple[str, list[Any]]:
    """``column IN (...)`` with literal parameters, which DuckDB answers from the index."""
    if not ids:
        return "FALSE", []
    if len(ids) <= 5000:
        return f"{column} IN ({','.join('?' * len(ids))})", list(ids)
    return f"{column} IN (SELECT UNNEST(?::VARCHAR[]))", [list(ids)]


def _stable(key: str) -> str:
    return hashlib.md5(key.encode()).hexdigest()


def best_ko(store, protein_id: str) -> str | None:
    rows = _rows(store, """SELECT accession FROM annotations
                           WHERE protein_id = ? AND LOWER(source) IN ('kofam', 'kegg')
                             AND regexp_matches(accession, '^K[0-9]{5}$')
                           ORDER BY evalue NULLS LAST, score DESC NULLS LAST LIMIT 1""", [protein_id])
    return rows[0][0] if rows else None


def family_of(store, protein_id: str, domains: list[dict[str, Any]], ko: str | None = None) -> dict[str, str] | None:
    """{'kind': 'ko'|'domain', 'id', 'name', 'accession'} for the best KO, else the best Pfam domain."""
    ko = ko or best_ko(store, protein_id)
    if ko:
        return {"kind": "ko", "id": ko, "name": ko, "accession": ko}
    placed = [d for d in domains if d.get("accession")]
    if not placed:
        return None
    best = min(placed, key=lambda d: d["evalue"] if d.get("evalue") is not None else 1.0)
    return {"kind": "domain", "id": best["accession"].split(".")[0], "name": best["name"],
            "accession": best["accession"]}


def carrier_ids(store, family: dict[str, str]) -> list[str]:
    if family["kind"] == "ko":
        sql = "SELECT DISTINCT protein_id FROM annotations WHERE LOWER(source) IN ('kofam', 'kegg') AND accession = ?"
    else:
        sql = "SELECT DISTINCT protein_id FROM annotations WHERE source = 'pfam' AND accession = ?"
    return [r[0] for r in _rows(store, sql, [family["accession"]])]


def _clades(catalog, bins: list[str], top: int = 6) -> dict[str, Any]:
    """Genomes over the first rank that splits them."""
    genomes = [catalog.by_bin[b] for b in set(bins) if b in catalog.by_bin]
    rank = next((r for r in RANKS[1:] if len({g.taxonomy.get(r) or UNCLASSIFIED for g in genomes}) > 1), "genus")
    counts = Counter(g.taxonomy.get(rank) or UNCLASSIFIED for g in genomes)
    rows = counts.most_common(top)
    return {"rank": rank, "rows": rows, "other": sum(counts.values()) - sum(n for _, n in rows),
            "total": sum(counts.values()), "max": rows[0][1] if rows else 1}


def _length_summary(lengths: list[int], own: int) -> dict[str, Any]:
    if not lengths:
        return {}
    median = statistics.median(lengths)
    out: dict[str, Any] = {"min": min(lengths), "median": median, "max": max(lengths),
                           "percentile": sum(1 for n in lengths if n < own) / len(lengths)}
    if median and own < 0.7 * median:
        out["note"] = "shorter"
    elif median and own > 1.3 * median:
        out["note"] = "longer"
    return out


def same_architecture(store, catalog, protein_id: str, domains: list[dict[str, Any]], own_length: int, *,
                      ko: str | None = None) -> dict[str, Any] | None:
    """Proteins with this one's exact resolved Pfam architecture (or, with no domains, its best KO)."""
    from sharur.architecture import _hits, compact, resolve  # noqa: PLC0415

    accessions = sorted({d["accession"] for d in domains if d.get("accession")})
    if accessions:
        sizes = {}
        for acc in accessions:
            known = catalog.domains.get(acc.split(".")[0], {}).get("proteins")
            sizes[acc] = known if known is not None else _rows(
                store, "SELECT COUNT(*) FROM annotations WHERE source = 'pfam' AND accession = ?", [acc])[0][0]
        rarest = min(accessions, key=lambda a: sizes[a])
        marks = ",".join("?" * len(accessions))
        # every protein carrying all of these domains; its raw hits may hold extras that resolve away
        same_set = sorted((r[0] for r in _rows(store, f"""
            WITH c AS (SELECT DISTINCT protein_id FROM annotations WHERE source = 'pfam' AND accession = ?)
            SELECT a.protein_id FROM annotations a JOIN c USING (protein_id)
            WHERE a.source = 'pfam' AND a.accession IN ({marks}) GROUP BY 1
            HAVING COUNT(DISTINCT a.accession) = {len(accessions)}""", [rarest, *accessions])), key=_stable)
        checked = same_set[:MAX_CANDIDATES]
        if protein_id not in checked:
            checked.append(protein_id)
        target = compact([d["name"] for d in domains])
        where, params = _in("a.protein_id", checked)
        hits = _hits(store, ["pfam"], f"AND {where}", params)
        exact = [pid for pid in checked if compact([d.name for d in resolve(hits.get(pid, []))]) == target]
        if protein_id not in exact:
            exact.append(protein_id)
        sampled = len(same_set) > MAX_CANDIDATES
        basis = {"kind": "architecture", "label": target, "same_domains": len(same_set), "sampled": sampled,
                 "checked": len(checked)}
    elif ko:
        exact = carrier_ids(store, {"kind": "ko", "accession": ko})
        if protein_id not in exact:
            exact.append(protein_id)
        basis = {"kind": "ko", "label": ko, "same_domains": None, "sampled": False, "checked": len(exact)}
    else:
        return None
    where, params = _in("protein_id", exact)
    info = _rows(store, f"""SELECT protein_id, bin_id, sequence_length, COALESCE(partial, '00')
                            FROM proteins WHERE {where}""", params)
    others = [(b, n, part) for pid, b, n, part in info if pid != protein_id]
    estimated = None
    if basis["sampled"] and basis["checked"]:
        estimated = round(len(exact) * basis["same_domains"] / basis["checked"])
    other_bins = [b for b, _, _ in others]
    return {**basis, "proteins": len(exact), "estimated": estimated, "others": len(others),
            "genomes": len(set(other_bins)), "clades": _clades(catalog, other_bins),
            "lengths": _length_summary([n for _, n, _ in others if n], own_length),
            "partial_others": sum(1 for _, _, part in others if part != "00")}


def paralogs(store, protein_id: str, bin_id: str, *, ko: str | None,
             domains: list[dict[str, Any]] | None = None) -> list[dict[str, Any]]:
    """Other proteins in this genome sharing the best KO or all of the Pfam domains."""
    shared: dict[str, set[str]] = defaultdict(set)
    if ko:
        for (pid,) in _rows(store, """SELECT DISTINCT a.protein_id FROM annotations a JOIN proteins p USING (protein_id)
                                      WHERE LOWER(a.source) IN ('kofam', 'kegg') AND a.accession = ? AND p.bin_id = ?""",
                            [ko, bin_id]):
            shared[pid].add(f"same KO ({ko})")
    accessions = sorted({d["accession"] for d in (domains or []) if d.get("accession")})
    if accessions:
        marks = ",".join("?" * len(accessions))
        for (pid,) in _rows(store, f"""
                SELECT a.protein_id FROM annotations a JOIN proteins p USING (protein_id)
                WHERE a.source = 'pfam' AND p.bin_id = ? AND a.accession IN ({marks}) GROUP BY 1
                HAVING COUNT(DISTINCT a.accession) = {len(accessions)}""", [bin_id, *accessions]):
            shared[pid].add("same Pfam domains")
    shared.pop(protein_id, None)
    if not shared:
        return []
    where, params = _in("protein_id", list(shared) + [protein_id])
    rows = {pid: (contig, start, gi, n, part) for pid, contig, start, gi, n, part in _rows(
        store, f"""SELECT protein_id, contig_id, start, gene_index, sequence_length, COALESCE(partial, '00')
                   FROM proteins WHERE {where}""", params)}
    me = rows.get(protein_id)
    out = []
    for pid, why in shared.items():
        contig, start, gi, n, part = rows.get(pid, (None, None, None, None, "00"))
        gap = None
        if me and contig == me[0] and contig != pid and start and me[1] and gi is not None and me[2] is not None:
            gap = abs(gi - me[2])
        out.append({"protein_id": pid, "contig_id": contig, "start": start, "length": n, "partial": part != "00",
                    "why": sorted(why), "gap": gap,
                    "tandem": "adjacent" if gap == 1 else ("nearby" if gap is not None and gap <= 3 else None)})
    out.sort(key=lambda r: (r["gap"] is None, r["gap"] or 0, r["contig_id"] or "", r["start"] or 0))
    return out


def neighbourhood(store, family: dict[str, str], *, max_genomes: int = MAX_GENOMES, flank: int = FLANK,
                  top: int = 12) -> dict[str, Any] | None:
    """Families within ``flank`` genes of the family's carriers, ranked by how many carriers have them nearby."""
    ids = carrier_ids(store, family)
    if not ids:
        return None
    sample = sorted(ids, key=_stable)[:CARRIER_SAMPLE]
    where, params = _in("protein_id", sample)
    rows = _rows(store, f"""SELECT protein_id, bin_id, contig_id, gene_index, strand FROM proteins
                            WHERE {where} AND start > 0 AND contig_id <> protein_id AND gene_index IS NOT NULL""",
                 params)
    by_genome: dict[str, tuple] = {}
    for row in sorted(rows, key=lambda r: _stable(r[0])):
        by_genome.setdefault(row[1], row)
    picked = [by_genome[b] for b in sorted(by_genome, key=_stable)][:max_genomes]
    base = {"family": family, "carriers": len(ids), "sampled": len(picked), "flank": flank, "rows": [],
            "vocabulary": None}
    if not picked:
        return base
    strand = {pid: 1 if s in ("+", "1") else -1 for pid, _, _, _, s in picked}
    window = _rows(store, """
        WITH c AS (SELECT UNNEST(?::VARCHAR[]) AS anchor, UNNEST(?::VARCHAR[]) AS contig_id,
                          UNNEST(?::INTEGER[]) AS gi)
        SELECT c.anchor, p.protein_id, p.gene_index - c.gi
        FROM c JOIN proteins p ON p.contig_id = c.contig_id
        WHERE p.gene_index BETWEEN c.gi - ? AND c.gi + ? AND p.protein_id <> c.anchor""",
                   [[p[0] for p in picked], [p[2] for p in picked], [p[3] for p in picked], flank, flank])
    neighbours = sorted({pid for _, pid, _ in window})
    ko_of: dict[str, str] = {}
    pfam_of: dict[str, tuple[str, str]] = {}
    if neighbours:
        where, params = _in("protein_id", neighbours)
        for pid, source, acc, name in _rows(store, f"""
                SELECT protein_id, LOWER(source), accession, COALESCE(NULLIF(name, ''), accession)
                FROM annotations WHERE {where} AND LOWER(source) IN ('kofam', 'kegg', 'pfam')
                ORDER BY evalue NULLS LAST""", params):
            if source == "pfam":
                pfam_of.setdefault(pid, (acc.split(".")[0], name))
            elif len(acc) == 6 and acc.startswith("K"):
                ko_of.setdefault(pid, acc)
    # name neighbours in whichever vocabulary covers more of them
    use_ko = len(ko_of) > len(pfam_of)
    seen: dict[str, set[str]] = defaultdict(set)
    offsets: dict[str, list[int]] = defaultdict(list)
    names: dict[str, str] = {}
    for anchor, pid, offset in window:
        if use_ko:
            key = label = ko_of.get(pid)
        else:
            key, label = pfam_of.get(pid, (None, None))
        if not key:
            continue
        names[key] = label
        if anchor not in seen[key]:
            seen[key].add(anchor)
            offsets[key].append(offset * strand[anchor])   # + downstream of the carrier, - upstream
    out = []
    for key, anchors in sorted(seen.items(), key=lambda kv: (-len(kv[1]), kv[0]))[:top]:
        median = statistics.median(offsets[key])
        out.append({"id": key, "name": names[key], "kind": "ko" if use_ko else "domain", "carriers": len(anchors),
                    "share": len(anchors) / len(picked), "offset": median,
                    "side": "downstream" if median > 0 else ("upstream" if median < 0 else "either side")})
    return {**base, "rows": out, "vocabulary": "KO" if use_ko else "Pfam domain",
            "named": len(ko_of if use_ko else pfam_of), "neighbours": len(neighbours)}


def contig_position(store, protein_id: str, contig_id: str) -> dict[str, Any] | None:
    """Where the gene sits among the genes of its contig (None without a genomic position)."""
    if not contig_id or contig_id == protein_id:
        return None
    genes = _rows(store, "SELECT protein_id, start, end_coord FROM proteins WHERE contig_id = ? ORDER BY start, end_coord",
                  [contig_id])
    index = next((i for i, (pid, _, _) in enumerate(genes) if pid == protein_id), None)
    if index is None or genes[index][1] <= 0:
        return None
    return {"genes": [(s, e) for _, s, e in genes], "index": index + 1, "count": len(genes),
            "start": genes[index][1], "end": genes[index][2]}
