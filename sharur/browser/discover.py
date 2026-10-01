"""Discovery feeds: extremes and oddities worth a second look.

Computed once in the background after the catalog loads (about a second per
feed on a few million proteins):

- ``giants``: the longest proteins, with their Pfam architecture and whether
  the gene runs off a contig end (Prodigal ``partial``), plus giants per
  genome by clade.
- ``dark``: the longest proteins with no annotation from any source.
- ``repeats``: the longest uninterrupted runs of one Pfam domain.
- ``fusions``: Pfam domain pairs on one protein that recur in at least three
  genomes, all within one order, family or genus, while each domain is common
  outside that clade and mostly found apart. Pairs from the same clade with
  shared domains and carriers merge into one entry.
- ``islands``: the longest runs of consecutive genes with no annotation, with
  whether they span a whole contig or reach its end.
- ``rare_systems``: curated system types called in at most three genomes.

Everything here is observed: domain co-occurrence, length, annotation gaps.
Names for systems come only from the curated callers.
"""

from __future__ import annotations

from collections import Counter
from typing import Any

from sharur.browser.catalog import RANKS, UNCLASSIFIED

GIANT_AA = 3000
LIMIT = 300


def _positioned(where_alias: str = "p") -> str:
    """Proteins with a genomic position (excludes per-protein and shared pseudo-contigs)."""
    a = where_alias
    return f"{a}.contig_id <> {a}.protein_id AND {a}.start > 0"


def giants(store) -> list[dict[str, Any]]:
    rows = store.execute(f"""
        SELECT protein_id, bin_id, contig_id, sequence_length, COALESCE(partial, '00')
        FROM proteins ORDER BY sequence_length DESC, protein_id LIMIT {LIMIT}""")
    return [{"protein_id": p, "bin_id": b, "contig_id": c, "length": n, "partial": part != "00",
             "partial_code": part} for p, b, c, n, part in rows]


def dark(store) -> list[dict[str, Any]]:
    rows = store.execute(f"""
        SELECT p.protein_id, p.bin_id, p.sequence_length, COALESCE(p.partial, '00')
        FROM proteins p ANTI JOIN annotations a USING (protein_id)
        ORDER BY p.sequence_length DESC, p.protein_id LIMIT {LIMIT}""")
    return [{"protein_id": p, "bin_id": b, "length": n, "partial": part != "00"} for p, b, n, part in rows]


def repeats(store) -> list[dict[str, Any]]:
    rows = store.execute(f"""
        WITH d AS (
            SELECT protein_id, COALESCE(NULLIF(name, ''), accession) AS domain, start_aa, end_aa,
                   -- annotation_id breaks ties between hits sharing a start, so runs are deterministic
                   ROW_NUMBER() OVER (PARTITION BY protein_id ORDER BY start_aa, accession, annotation_id)
                 - ROW_NUMBER() OVER (PARTITION BY protein_id, COALESCE(NULLIF(name, ''), accession)
                                      ORDER BY start_aa, accession, annotation_id) AS island
            FROM annotations WHERE LOWER(source) = 'pfam' AND start_aa IS NOT NULL
        ), runs AS (
            SELECT protein_id, domain, COUNT(*) AS run, MIN(start_aa) AS lo, MAX(end_aa) AS hi
            FROM d GROUP BY protein_id, domain, island
        ), best AS (
            SELECT *, ROW_NUMBER() OVER (PARTITION BY protein_id ORDER BY run DESC, domain) AS pick FROM runs
        )
        SELECT r.protein_id, p.bin_id, p.sequence_length, r.domain, r.run, r.lo, r.hi
        FROM best r JOIN proteins p USING (protein_id)
        WHERE r.pick = 1 AND r.run >= 6
        ORDER BY r.run DESC, p.sequence_length DESC, r.protein_id LIMIT {LIMIT}""")
    return [{"protein_id": p, "bin_id": b, "length": n, "domain": d, "run": r, "lo": lo, "hi": hi}
            for p, b, n, d, r, lo, hi in rows]


def _lca(genomes) -> tuple[str, str] | None:
    """Deepest rank shared by every genome."""
    lca = None
    for rank in RANKS:
        values = {g.taxonomy.get(rank) for g in genomes}
        if len(values) == 1:
            value = next(iter(values))
            if value and value != UNCLASSIFIED:
                lca = (rank, value)
                continue
        break
    return lca


def fusions(store, catalog, *, min_genomes: int = 3, outside: int = 20,
            max_share: float = 0.2) -> tuple[list[dict[str, Any]], int]:
    """Clade-restricted domain fusions (merged), and how many there are before the display limit."""
    rows = store.execute(f"""
        WITH d AS (
            SELECT DISTINCT a.protein_id, p.bin_id, split_part(a.accession, '.', 1) AS acc
            FROM annotations a JOIN proteins p USING (protein_id) WHERE LOWER(a.source) = 'pfam'),
        names AS (SELECT split_part(accession, '.', 1) AS acc, ANY_VALUE(COALESCE(NULLIF(name, ''), accession)) AS nm
                  FROM annotations WHERE LOWER(source) = 'pfam' GROUP BY 1),
        n AS (SELECT protein_id, COUNT(*) AS k FROM d GROUP BY 1),
        dd AS (SELECT d.* FROM d JOIN n USING (protein_id) WHERE n.k BETWEEN 2 AND 15),
        freq AS (SELECT acc, COUNT(DISTINCT protein_id) AS proteins, COUNT(DISTINCT bin_id) AS genomes FROM d GROUP BY 1),
        pairs AS (
            SELECT x.acc AS a, y.acc AS b, COUNT(DISTINCT x.protein_id) AS proteins, LIST(DISTINCT x.bin_id) AS bins,
                   MIN(x.protein_id) AS example
            FROM dd x JOIN dd y ON x.protein_id = y.protein_id AND x.acc < y.acc
            GROUP BY 1, 2 HAVING COUNT(DISTINCT x.bin_id) >= {int(min_genomes)})
        SELECT p.a, na.nm, p.b, nb.nm, p.proteins, p.bins, p.example, fa.proteins, fa.genomes, fb.proteins, fb.genomes
        FROM pairs p JOIN freq fa ON fa.acc = p.a JOIN freq fb ON fb.acc = p.b
                     JOIN names na ON na.acc = p.a JOIN names nb ON nb.acc = p.b""")
    clade_size: dict[tuple[str, str], int] = {}
    kept = []
    for a, a_name, b, b_name, proteins, bins, example, a_prot, a_gen, b_prot, b_gen in rows:
        genomes = [catalog.by_bin[x] for x in bins if x in catalog.by_bin]
        lca = _lca(genomes)
        if lca is None or lca[0] not in ("order", "family", "genus"):
            continue
        if a_gen - len(genomes) < outside or b_gen - len(genomes) < outside:
            continue
        if proteins > max_share * a_prot or proteins > max_share * b_prot:
            continue
        if lca not in clade_size:
            clade_size[lca] = sum(1 for g in catalog.genomes if g.taxonomy.get(lca[0]) == lca[1])
        kept.append({"domains": [(a, a_name), (b, b_name)], "genomes": len(genomes), "bins": set(bins),
                     "proteins": proteins, "example": example, "rank": lca[0], "clade": lca[1],
                     "clade_size": clade_size[lca], "domain_proteins": {a: a_prot, b: b_prot}})
    kept.sort(key=lambda f: (-f["genomes"], -f["genomes"] / f["clade_size"], f["example"]))
    merged: list[dict[str, Any]] = []
    for f in kept:
        accs = {acc for acc, _ in f["domains"]}
        for m in merged:
            if (m["rank"], m["clade"]) != (f["rank"], f["clade"]):
                continue
            if not accs & {acc for acc, _ in m["domains"]}:
                continue
            if len(m["bins"] & f["bins"]) / len(m["bins"] | f["bins"]) >= 0.6:
                m["domains"] += [d for d in f["domains"] if d[0] not in {acc for acc, _ in m["domains"]}]
                m["domain_proteins"].update(f["domain_proteins"])
                break
        else:
            merged.append(f)
    for m in merged:
        m["share"] = m["genomes"] / m["clade_size"] if m["clade_size"] else 0.0
        del m["bins"]
    return merged[:LIMIT], len(merged)


def islands(store, *, min_genes: int = 6) -> list[dict[str, Any]]:
    rows = store.execute(f"""
        WITH ann AS (SELECT DISTINCT protein_id FROM annotations),
        shared AS (SELECT contig_id FROM proteins GROUP BY 1 HAVING COUNT(DISTINCT bin_id) > 1),
        p AS (
            SELECT p.protein_id, p.contig_id, p.bin_id, p.start, p.end_coord, p.strand,
                   (a.protein_id IS NULL) AS dark,
                   ROW_NUMBER() OVER (PARTITION BY p.contig_id ORDER BY p.start, p.end_coord) AS pos,
                   COUNT(*) OVER (PARTITION BY p.contig_id) AS genes_on_contig
            FROM proteins p LEFT JOIN ann a USING (protein_id)
            WHERE {_positioned()} AND p.contig_id NOT IN (SELECT contig_id FROM shared)),
        r AS (SELECT *, pos - ROW_NUMBER() OVER (PARTITION BY contig_id, dark ORDER BY pos) AS grp FROM p),
        runs AS (
            SELECT contig_id, ANY_VALUE(bin_id) AS bin_id, COUNT(*) AS genes, MIN(start) AS lo, MAX(end_coord) AS hi,
                   MIN(pos) AS first_pos, MAX(pos) AS last_pos, ANY_VALUE(genes_on_contig) AS contig_genes,
                   ARG_MIN(protein_id, pos) AS first_protein,
                   STRING_AGG(CASE WHEN strand IN ('+', '1') THEN '+' ELSE '-' END, '' ORDER BY pos) AS strands
            FROM r WHERE dark GROUP BY contig_id, grp HAVING COUNT(*) >= {int(min_genes)})
        SELECT runs.*, c.length FROM runs LEFT JOIN contigs c USING (contig_id)
        ORDER BY genes DESC, hi - lo DESC, contig_id LIMIT {LIMIT}""")
    out = []
    for contig, bin_id, genes, lo, hi, first, last, contig_genes, first_protein, strands, length in rows:
        out.append({"contig_id": contig, "bin_id": bin_id, "genes": genes, "start": lo, "end": hi,
                    "whole_contig": first == 1 and last == contig_genes,
                    "at_edge": first == 1 or last == contig_genes, "contig_genes": contig_genes,
                    "first_protein": first_protein, "strands": strands, "contig_length": length})
    return out


def island_count(store, *, min_genes: int = 10) -> int:
    return int(store.execute(f"""
        WITH ann AS (SELECT DISTINCT protein_id FROM annotations),
        shared AS (SELECT contig_id FROM proteins GROUP BY 1 HAVING COUNT(DISTINCT bin_id) > 1),
        p AS (SELECT p.contig_id, (a.protein_id IS NULL) AS dark,
                     ROW_NUMBER() OVER (PARTITION BY p.contig_id ORDER BY p.start, p.end_coord) AS pos
              FROM proteins p LEFT JOIN ann a USING (protein_id)
              WHERE {_positioned()} AND p.contig_id NOT IN (SELECT contig_id FROM shared)),
        r AS (SELECT *, pos - ROW_NUMBER() OVER (PARTITION BY contig_id, dark ORDER BY pos) AS grp FROM p)
        SELECT COUNT(*) FROM (SELECT 1 FROM r WHERE dark GROUP BY contig_id, grp HAVING COUNT(*) >= {int(min_genes)})
        """)[0][0])


def giant_clades(store, catalog, *, threshold: int = GIANT_AA, min_genomes: int = 3) -> dict[str, Any]:
    """Proteins of at least ``threshold`` aa per genome, by the rank that best splits the dataset."""
    counts = dict(store.execute(
        f"SELECT bin_id, COUNT(*) FROM proteins WHERE sequence_length >= {int(threshold)} GROUP BY 1"))
    rank, _ = catalog.children(catalog.genomes, None)
    for candidate in ("order", "class", "family"):
        named = {g.taxonomy.get(candidate) for g in catalog.genomes} - {None, UNCLASSIFIED}
        if 4 <= len(named) <= 200:
            rank = candidate
            break
    genomes: Counter = Counter()
    giant_total: Counter = Counter()
    carriers: Counter = Counter()
    for g in catalog.genomes:
        clade = g.taxonomy.get(rank) or UNCLASSIFIED
        genomes[clade] += 1
        n = counts.get(g.bin_id, 0)
        giant_total[clade] += n
        carriers[clade] += n > 0
    rows = [{"clade": c, "genomes": genomes[c], "giants": giant_total[c], "per_genome": giant_total[c] / genomes[c],
             "carriers": carriers[c] / genomes[c]}
            for c in genomes if genomes[c] >= min_genomes and c != UNCLASSIFIED]
    rows.sort(key=lambda r: -r["per_genome"])
    return {"rank": rank, "threshold": threshold, "rows": rows,
            "total": sum(counts.values()), "genomes_with": sum(1 for v in counts.values() if v)}


def rare_systems(catalog, cctyper=None, *, max_genomes: int = 3) -> list[dict[str, Any]]:
    carriers: dict[tuple[str, str], set[str]] = {}
    for s in catalog.systems:
        carriers.setdefault((s["kind"], s["type"]), set()).add(s["bin_id"])
    for s in (cctyper() if cctyper else []):
        if s.get("confident"):
            carriers.setdefault(("crispr", s["prediction"]), set()).add(s["bin_id"])
    out = []
    for (kind, kind_type), bins in carriers.items():
        if len(bins) <= max_genomes:
            genomes = [catalog.by_bin[b] for b in sorted(bins) if b in catalog.by_bin]
            out.append({"kind": kind, "type": kind_type, "genomes": [g.bin_id for g in genomes],
                        "labels": [g.label for g in genomes]})
    out.sort(key=lambda r: (len(r["genomes"]), r["kind"], r["type"]))
    return out


def architectures(store, protein_ids: list[str]) -> dict[str, list[dict[str, Any]]]:
    """Resolved Pfam architecture for each protein, in one query."""
    from sharur.architecture import _hits, resolve  # noqa: PLC0415

    if not protein_ids:
        return {}
    hits = _hits(store, ["pfam"], "AND a.protein_id IN (SELECT UNNEST(?::VARCHAR[]))", [protein_ids])
    return {pid: [d.to_dict() for d in resolve(found)] for pid, found in hits.items()}


def compute(store, catalog) -> dict[str, Any]:
    """Every feed but ``rare_systems`` (which needs the CRISPR-Cas calls; see the route)."""
    feeds: dict[str, Any] = {"giants": giants(store), "dark": dark(store), "repeats": repeats(store)}
    arch = architectures(store, [p["protein_id"] for p in feeds["giants"]] +
                         [p["protein_id"] for p in feeds["repeats"]])
    from sharur.architecture import compact  # noqa: PLC0415

    for row in feeds["giants"] + feeds["repeats"]:
        row["domains"] = arch.get(row["protein_id"], [])
        row["architecture"] = compact([d["name"] for d in row["domains"]])
    feeds["fusions"], fusion_total = fusions(store, catalog)
    feeds["islands"] = islands(store)
    feeds["giant_clades"] = giant_clades(store, catalog)
    feeds["stats"] = {
        "giants": feeds["giant_clades"]["total"],
        "dark_1000": int(store.execute(
            "SELECT COUNT(*) FROM proteins p ANTI JOIN annotations a USING (protein_id) "
            "WHERE p.sequence_length >= 1000")[0][0]),
        "fusions": fusion_total,
        "islands_10": island_count(store),
    }
    return feeds


FEEDS = {
    "giants": ("Giant proteins",
               "The longest proteins, drawn to one scale with their Pfam domains. A dashed tail marks a gene that "
               "runs off its contig: that protein is longer than its sequence here."),
    "clades": ("Where the giants are",
               f"Proteins of {GIANT_AA:,} aa or more per genome, by clade."),
    "fusions": ("Clade-specific domain fusions",
                "Pfam domains found together on one protein in at least three genomes, all within one order, family "
                "or genus, while each domain is common elsewhere and usually found apart. Domain co-occurrence is "
                "observed here; what the fused protein does is a hypothesis."),
    "islands": ("Unannotated stretches",
                "The longest runs of consecutive genes with no annotation from any source. Whole unannotated contigs "
                "are worth checking for phage, plasmid or contamination."),
    "dark": ("Largest unannotated proteins",
             "The longest proteins with no hit in any annotation source."),
    "repeats": ("Longest domain repeats",
                "Proteins with the longest uninterrupted runs of one Pfam domain; the run is outlined."),
    "systems": ("Rare systems",
                "Curated system types called in three genomes or fewer."),
}
