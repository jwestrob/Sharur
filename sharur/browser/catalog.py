"""Dataset summaries for browsing, computed once when the browser starts.

Fast summaries (genomes, taxonomy, predicate presence, systems, loci) load
synchronously; slower ones (KEGG module completeness, notable proteins) load in
a background thread and pages show them once ready.

Genome x predicate presence is held as three numpy arrays (genome index,
predicate index, protein count) so clade enrichment and per-clade prevalence
are vectorized lookups.
"""

from __future__ import annotations

import logging
import re
import threading
import time
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from sharur.predicates.vocabulary import PREDICATE_BY_ID

logger = logging.getLogger(__name__)

RANKS = ("domain", "phylum", "class", "order", "family", "genus")
_GTDB_PREFIX = {"d": "domain", "p": "phylum", "c": "class", "o": "order", "f": "family", "g": "genus"}
UNCLASSIFIED = "Unclassified"

CATEGORY_LABELS = {
    "enzyme": "Enzymes", "metabolism": "Metabolism", "transport": "Transport", "binding": "Cofactors & binding",
    "info_processing": "Information processing", "regulation": "Regulation", "envelope": "Cell envelope",
    "mobile": "Mobile elements", "viral": "Viral", "stress": "Stress & defense", "cazy": "Carbohydrate-active",
    "structure": "Structure", "division": "Cell division", "topology": "Topology & localization",
    "annotation": "Annotation status", "size": "Size", "composition": "Composition",
}
# categories that describe biology (the rest describe the annotation itself)
BIOLOGICAL = ("metabolism", "enzyme", "transport", "binding", "info_processing", "regulation", "envelope",
              "stress", "mobile", "viral", "cazy", "division", "structure")


@dataclass
class Genome:
    index: int
    bin_id: str
    taxonomy: dict[str, str]
    completeness: float | None
    contamination: float | None
    contigs: int
    length: int
    n50: int
    proteins: int
    annotated: int
    longest_protein: int

    @property
    def lineage(self) -> list[tuple[str, str]]:
        return [(r, self.taxonomy[r]) for r in RANKS if self.taxonomy.get(r)]

    @property
    def label(self) -> str:
        for rank in reversed(RANKS):
            name = self.taxonomy.get(rank)
            if name and name != UNCLASSIFIED:
                return name
        return UNCLASSIFIED


@dataclass
class Catalog:
    genomes: list[Genome] = field(default_factory=list)
    by_bin: dict[str, Genome] = field(default_factory=dict)
    predicates: list[str] = field(default_factory=list)
    predicate_index: dict[str, int] = field(default_factory=dict)
    predicate_proteins: np.ndarray | None = None
    predicate_genomes: np.ndarray | None = None
    pair_bin: np.ndarray | None = None
    pair_pred: np.ndarray | None = None
    pair_count: np.ndarray | None = None
    category_share: dict[str, dict[str, float]] = field(default_factory=dict)  # bin -> category -> share
    dataset_category_share: dict[str, float] = field(default_factory=dict)
    systems: list[dict[str, Any]] = field(default_factory=list)
    loci: list[dict[str, Any]] = field(default_factory=list)
    sources: list[dict[str, Any]] = field(default_factory=list)
    totals: dict[str, int] = field(default_factory=dict)
    # background
    modules: dict[str, Any] = field(default_factory=dict)          # module -> definition
    module_completeness: np.ndarray | None = None                  # genomes x modules
    module_ids: list[str] = field(default_factory=list)
    ko_sets: dict[str, set[str]] = field(default_factory=dict)
    notable: dict[str, list[dict[str, Any]]] = field(default_factory=dict)
    domains: dict[str, dict[str, Any]] = field(default_factory=dict)  # Pfam accession -> summary
    vogs: dict[str, dict[str, Any]] = field(default_factory=dict)     # VOG id -> summary (+ VOGdb annotation)
    ready: threading.Event = field(default_factory=threading.Event)
    status: str = "loading"
    map_state: str = "unknown"

    # ------------------------------------------------------------------ #
    # Taxonomy
    # ------------------------------------------------------------------ #

    def clade(self, rank: str | None = None, name: str | None = None) -> list[Genome]:
        if rank is None:
            return self.genomes
        return [g for g in self.genomes if g.taxonomy.get(rank) == name]

    def children(self, genomes: list[Genome], rank: str | None) -> tuple[str | None, list[tuple[str, int]]]:
        """Next rank below ``rank`` that splits these genomes, with counts."""
        start = 0 if rank is None else RANKS.index(rank) + 1
        for child in RANKS[start:]:
            counts = Counter(g.taxonomy.get(child) or UNCLASSIFIED for g in genomes)
            named = [n for n in counts if n != UNCLASSIFIED]
            if len(named) > 1 or child == RANKS[-1]:
                return child, counts.most_common()
        return None, []

    def lineage_of(self, rank: str, name: str) -> list[tuple[str, str]]:
        members = self.clade(rank, name)
        if not members:
            return []
        lineage = []
        for r in RANKS[: RANKS.index(rank) + 1]:
            values = {g.taxonomy.get(r) for g in members}
            if len(values) == 1 and next(iter(values)):
                lineage.append((r, next(iter(values))))
        return lineage

    # ------------------------------------------------------------------ #
    # Predicates
    # ------------------------------------------------------------------ #

    def carriers(self, predicate: str) -> dict[int, int]:
        """genome index -> proteins carrying the predicate."""
        k = self.predicate_index.get(predicate)
        if k is None or self.pair_pred is None:
            return {}
        mask = self.pair_pred == k
        return dict(zip(self.pair_bin[mask].tolist(), self.pair_count[mask].tolist(), strict=True))

    def enrichment(self, genomes: list[Genome], *, min_genomes: int = 3, top: int = 20) -> list[dict[str, Any]]:
        """Predicates over-represented in these genomes relative to the rest of the dataset."""
        if self.pair_bin is None or not genomes or len(genomes) == len(self.genomes):
            return []
        inside = np.zeros(len(self.genomes), dtype=bool)
        inside[[g.index for g in genomes]] = True
        n_in, n_out = int(inside.sum()), len(self.genomes) - int(inside.sum())
        hit_in = np.bincount(self.pair_pred[inside[self.pair_bin]], minlength=len(self.predicates))
        hit_out = self.predicate_genomes - hit_in
        frac_in, frac_out = hit_in / n_in, hit_out / max(n_out, 1)
        rows = []
        for k in np.argsort(-(frac_in - frac_out)):
            pred = self.predicates[k]
            definition = PREDICATE_BY_ID.get(pred)
            if hit_in[k] < min_genomes or definition is None or definition.category not in BIOLOGICAL:
                continue
            if frac_in[k] - frac_out[k] <= 0.05:
                break
            rows.append({"predicate": pred, "name": definition.name, "category": definition.category,
                         "inside": float(frac_in[k]), "outside": float(frac_out[k]), "genomes": int(hit_in[k])})
            if len(rows) >= top:
                break
        return rows

    def prevalence_by(self, genome_indices: dict[int, Any], rank: str) -> list[dict[str, Any]]:
        """Share of genomes per taxon at ``rank`` that appear in ``genome_indices``."""
        totals, hits = Counter(), Counter()
        for g in self.genomes:
            taxon = g.taxonomy.get(rank) or UNCLASSIFIED
            totals[taxon] += 1
            if g.index in genome_indices:
                hits[taxon] += 1
        rows = [{"taxon": t, "genomes": n, "with": hits[t], "share": hits[t] / n} for t, n in totals.items()]
        if sum(hits.values()) < 0.05 * len(self.genomes):
            # rare features: list the clades that carry them, most affected first
            return sorted((r for r in rows if r["with"]), key=lambda r: (-r["share"], -r["genomes"], r["taxon"]))
        return sorted(rows, key=lambda r: (-r["genomes"], r["taxon"]))

    # ------------------------------------------------------------------ #
    # Modules
    # ------------------------------------------------------------------ #

    def clade_modules(self, genomes: list[Genome], *, mode: str = "median",
                      min_median: float = 0.5, min_present: float = 0.25,
                      complete: float = 0.75) -> list[dict[str, Any]]:
        """KEGG module completeness aggregated over ``genomes``.

        Per module: median and mean completeness across the clade, the share of
        genomes with the module at least ``complete`` (prevalence), the share
        with any step, and the same median and prevalence over the rest of the
        dataset for comparison. ``mode='median'`` keeps modules whose clade
        median is at least ``min_median``; ``mode='present'`` keeps modules with
        any step in at least ``min_present`` of the clade's genomes.
        """
        matrix = self.module_completeness
        if matrix is None or not genomes:
            return []
        inside = np.zeros(matrix.shape[0], dtype=bool)
        inside[[g.index for g in genomes]] = True
        clade, rest = matrix[inside], matrix[~inside]
        median = np.median(clade, axis=0)
        mean = clade.mean(axis=0)
        prevalence = (clade >= complete).mean(axis=0)
        present = (clade > 0).mean(axis=0)
        overall_median = np.median(matrix, axis=0)
        rest_prevalence = (rest >= complete).mean(axis=0) if len(rest) else np.zeros(matrix.shape[1])
        keep = median >= min_median if mode == "median" else present >= min_present
        rows = []
        for j in np.nonzero(keep)[0]:
            d = self.modules[self.module_ids[j]]
            rows.append({"module": d.module, "name": d.name, "class": d.module_class.split(";")[-1].strip(),
                         "median": float(median[j]), "mean": float(mean[j]), "prevalence": float(prevalence[j]),
                         "present": float(present[j]), "dataset_median": float(overall_median[j]),
                         "rest_prevalence": float(rest_prevalence[j])})
        rows.sort(key=lambda r: (-r["median"], -r["prevalence"], r["module"]))
        return rows

    def distinctive_modules(self, genomes: list[Genome], *, margin: float = 0.2,
                            complete: float = 0.75, top: int = 15) -> list[dict[str, Any]]:
        """Modules complete (>= ``complete``) in a larger share of the clade than of the rest, by >= ``margin``."""
        if self.module_completeness is None or not genomes or len(genomes) == len(self.genomes):
            return []
        rows = [r for r in self.clade_modules(genomes, mode="present", min_present=0.0, complete=complete)
                if r["prevalence"] - r["rest_prevalence"] >= margin]
        rows.sort(key=lambda r: (-(r["prevalence"] - r["rest_prevalence"]), r["module"]))
        return rows[:top]

    def module_column(self, module: str) -> np.ndarray | None:
        if self.module_completeness is None or module not in self.module_ids:
            return None
        return self.module_completeness[:, self.module_ids.index(module)]


# ---------------------------------------------------------------------- #
# Loading
# ---------------------------------------------------------------------- #


def _tables(store) -> set[str]:
    return {r[0] for r in store.execute("SELECT table_name FROM information_schema.tables")}


def _columns(store, table: str) -> set[str]:
    return {r[0] for r in store.execute(
        "SELECT column_name FROM information_schema.columns WHERE table_name = ?", [table])}


def _taxonomy(store, tables: set[str]) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = defaultdict(dict)
    for bin_id, text in store.execute("SELECT bin_id, taxonomy FROM bins"):
        for part in (text or "").split(";"):
            if "__" in part:
                prefix, name = part.split("__", 1)
                if prefix in _GTDB_PREFIX and name:
                    result[bin_id][_GTDB_PREFIX[prefix]] = name
    if "taxonomy" in tables:
        columns = _columns(store, "taxonomy")
        wanted = [(r, "order_" if r == "order" and "order_" in columns else r) for r in RANKS]
        if "bin_id" in columns and all(c in columns for _, c in wanted):
            select = ", ".join(f'"{c}"' for _, c in wanted)
            for row in store.execute(f"SELECT bin_id, {select} FROM taxonomy"):
                for (rank, _), value in zip(wanted, row[1:], strict=True):
                    if value:
                        result[row[0]][rank] = value
    return result


def load_catalog(store) -> Catalog:
    t0 = time.time()
    tables = _tables(store)
    catalog = Catalog()
    taxonomy = _taxonomy(store, tables)
    bin_columns = _columns(store, "bins")
    quality = {}
    if {"completeness", "contamination"} <= bin_columns:
        quality = {b: (c, x) for b, c, x in store.execute("SELECT bin_id, completeness, contamination FROM bins")}
    lengths: dict[str, list[int]] = defaultdict(list)
    for bin_id, length in store.execute("SELECT bin_id, length FROM contigs"):
        lengths[bin_id].append(length or 0)
    contig_stats = {b: (len(v), sum(v), _n50(v)) for b, v in lengths.items()}
    protein_stats = {b: (n, longest) for b, n, longest in store.execute(
        "SELECT bin_id, COUNT(*), MAX(sequence_length) FROM proteins GROUP BY 1")}
    annotated = dict(store.execute(
        """SELECT p.bin_id, COUNT(*) FROM proteins p
           JOIN (SELECT DISTINCT protein_id FROM annotations) a USING (protein_id) GROUP BY 1"""))
    for i, bin_id in enumerate(sorted(r[0] for r in store.execute("SELECT bin_id FROM bins"))):
        n_contigs, length, n50 = contig_stats.get(bin_id, (0, 0, 0))
        proteins, longest = protein_stats.get(bin_id, (0, 0))
        completeness, contamination = quality.get(bin_id, (None, None))
        tax = {r: taxonomy.get(bin_id, {}).get(r) or UNCLASSIFIED for r in RANKS}
        genome = Genome(i, bin_id, tax, completeness, contamination, n_contigs or 0, length or 0, n50 or 0,
                        proteins or 0, annotated.get(bin_id, 0), longest or 0)
        catalog.genomes.append(genome)
        catalog.by_bin[bin_id] = genome
    catalog.totals = {"genomes": len(catalog.genomes), "proteins": sum(g.proteins for g in catalog.genomes),
                      "contigs": sum(g.contigs for g in catalog.genomes),
                      "bases": sum(g.length for g in catalog.genomes),
                      "annotated": sum(g.annotated for g in catalog.genomes)}

    if "protein_predicates" in tables:
        pairs = store.execute("""
            SELECT pr.bin_id, x.p, COUNT(*) FROM
                (SELECT protein_id, UNNEST(predicates) AS p FROM protein_predicates) x
            JOIN proteins pr USING (protein_id) WHERE x.p NOT LIKE '%:%' GROUP BY 1, 2""")
        names = sorted({p for _, p, _ in pairs})
        catalog.predicates = names
        catalog.predicate_index = {p: k for k, p in enumerate(names)}
        bins = np.array([catalog.by_bin[b].index if b in catalog.by_bin else -1 for b, _, _ in pairs], dtype=np.int32)
        preds = np.array([catalog.predicate_index[p] for _, p, _ in pairs], dtype=np.int32)
        counts = np.array([n for _, _, n in pairs], dtype=np.int32)
        keep = bins >= 0
        catalog.pair_bin, catalog.pair_pred, catalog.pair_count = bins[keep], preds[keep], counts[keep]
        catalog.predicate_genomes = np.bincount(catalog.pair_pred, minlength=len(names))
        catalog.predicate_proteins = np.bincount(catalog.pair_pred, weights=catalog.pair_count, minlength=len(names))
        _category_shares(store, catalog)

    if "defense_systems" in tables:
        for system_id, genome_id, system_type, subtype, genes, protein_ids in store.execute(
                "SELECT system_id, genome_id, system_type, system_subtype, genes_count, protein_ids "
                "FROM defense_systems"):
            catalog.systems.append({"kind": "defense", "system_id": system_id, "bin_id": genome_id,
                                    "type": system_type, "subtype": subtype, "genes": genes,
                                    "proteins": _as_list(protein_ids)})
    if "secretion_systems" in tables:
        for system_id, genome_id, system_type, subtype, genes, protein_ids in store.execute(
                "SELECT system_id, genome_id, system_type, system_subtype, genes_count, protein_ids "
                "FROM secretion_systems"):
            catalog.systems.append({"kind": "secretion", "system_id": system_id, "bin_id": genome_id,
                                    "type": system_type, "subtype": subtype, "genes": genes,
                                    "proteins": _as_list(protein_ids)})
    if "loci" in tables:
        for locus_id, locus_type, contig_id, start, end, bin_id in store.execute(
                """SELECT l.locus_id, l.locus_type, l.contig_id, l.start, l.end_coord, c.bin_id
                   FROM loci l LEFT JOIN contigs c USING (contig_id)"""):
            catalog.loci.append({"locus_id": locus_id, "type": locus_type, "contig_id": contig_id,
                                 "start": start, "end": end, "bin_id": bin_id})
    try:
        from sharur.predicates.provenance import map_status  # noqa: PLC0415

        catalog.map_state = map_status(store).state if "predicate_provenance" in tables else "unstamped"
    except Exception:  # pragma: no cover - provenance is advisory here
        catalog.map_state = "unknown"
    catalog.sources = [{"source": s, "proteins": n} for s, n in store.execute(
        "SELECT LOWER(source), COUNT(DISTINCT protein_id) FROM annotations GROUP BY 1 ORDER BY 2 DESC")]
    logger.info("catalog loaded in %.1fs", time.time() - t0)
    return catalog


def _n50(lengths: list[int]) -> int:
    half, running = sum(lengths) / 2, 0
    for length in sorted(lengths, reverse=True):
        running += length
        if running >= half:
            return length
    return 0


def _as_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, list):
        return [str(v) for v in value]
    return [v for v in str(value).replace(";", ",").split(",") if v]


def _category_shares(store, catalog: Catalog) -> None:
    import pandas as pd  # noqa: PLC0415

    mapping = pd.DataFrame([(p, PREDICATE_BY_ID[p].category) for p in catalog.predicates if p in PREDICATE_BY_ID],
                           columns=["predicate", "category"])
    store.conn.register("_browser_categories", mapping)
    try:
        rows = store.execute("""
            SELECT pr.bin_id, m.category, COUNT(DISTINCT x.protein_id) FROM
                (SELECT protein_id, UNNEST(predicates) AS predicate FROM protein_predicates) x
            JOIN _browser_categories m USING (predicate)
            JOIN proteins pr USING (protein_id) GROUP BY 1, 2""")
    finally:
        store.conn.unregister("_browser_categories")
    totals: Counter = Counter()
    for bin_id, category, n in rows:
        genome = catalog.by_bin.get(bin_id)
        if genome and genome.proteins:
            catalog.category_share.setdefault(bin_id, {})[category] = n / genome.proteins
            totals[category] += n
    proteins = catalog.totals["proteins"] or 1
    catalog.dataset_category_share = {c: n / proteins for c, n in totals.items()}


def load_background(store, catalog: Catalog, lock: threading.Lock) -> None:
    """Module completeness and notable proteins (slow); sets ``catalog.ready``."""
    try:
        from sharur.modules import evaluate, load_modules  # noqa: PLC0415

        catalog.status = "computing pathways"
        definitions = load_modules()
        if definitions:
            with lock:
                hits = store.execute(
                    """SELECT p.bin_id, a.accession FROM annotations a JOIN proteins p USING (protein_id)
                       WHERE LOWER(a.source) IN ('kofam', 'kegg') AND regexp_matches(a.accession, '^K[0-9]{5}$')
                       GROUP BY 1, 2""")
            for bin_id, ko in hits:
                catalog.ko_sets.setdefault(bin_id, set()).add(ko)
            ids = sorted(definitions)
            matrix = np.zeros((len(catalog.genomes), len(ids)), dtype=np.float32)
            module_kos = {m: _kos(definitions[m]) for m in ids}
            for genome in catalog.genomes:
                present = catalog.ko_sets.get(genome.bin_id, set())
                if not present:
                    continue
                for j, m in enumerate(ids):
                    if module_kos[m] & present:
                        matrix[genome.index, j] = evaluate(definitions[m], present, definitions).completeness
            catalog.modules, catalog.module_ids, catalog.module_completeness = definitions, ids, matrix
        catalog.status = "summarizing Pfam domains"
        with lock:
            catalog.domains = _domains(store)
            catalog.vogs = _vogs(store)
        catalog.status = "finding notable proteins"
        with lock:
            notable = _notable(store)
            from sharur.architecture import architecture, compact  # noqa: PLC0415

            for row in notable["giants"][:60]:
                row["architecture"] = compact([d.name for d in architecture(store, row["protein_id"])])
        catalog.notable = notable
        catalog.status = "ready"
    except Exception:  # pragma: no cover - logged and surfaced as status
        logger.exception("background catalog load failed")
        catalog.status = "background summaries failed; see server log"
    finally:
        catalog.ready.set()


def _kos(definition) -> set[str]:
    from sharur.modules import kos_in  # noqa: PLC0415

    return kos_in(definition.tree)


def _domains(store) -> dict[str, dict[str, Any]]:
    rows = store.execute("""
        SELECT split_part(a.accession, '.', 1) AS acc, ANY_VALUE(NULLIF(a.name, '')), ANY_VALUE(NULLIF(a.description, '')),
               COUNT(DISTINCT a.protein_id), COUNT(DISTINCT p.bin_id), COUNT(*)
        FROM annotations a JOIN proteins p USING (protein_id)
        WHERE LOWER(a.source) = 'pfam' GROUP BY 1""")
    return {acc: {"accession": acc, "name": name or acc, "description": desc or "", "proteins": proteins,
                  "genomes": genomes, "hits": hits}
            for acc, name, desc, proteins, genomes, hits in rows}


def _vogs(store) -> dict[str, dict[str, Any]]:
    from sharur.predicates.mappings.vog_map import load_vog_annotations, vog_category_names  # noqa: PLC0415

    reference = load_vog_annotations()
    rows = store.execute("""
        SELECT a.accession, COUNT(DISTINCT a.protein_id), COUNT(DISTINCT p.bin_id), COUNT(*)
        FROM annotations a JOIN proteins p USING (protein_id)
        WHERE LOWER(a.source) IN ('vogdb', 'vog') GROUP BY 1""")
    out = {}
    for vog_id, proteins, genomes, hits in rows:
        ref = reference.get(vog_id, {})
        description = vog_description(str(ref.get("description") or ""))
        out[vog_id] = {"accession": vog_id, "name": description or vog_id, "description": description,
                       "category": ref.get("category") or "", "categories": vog_category_names(ref.get("category")),
                       "vogdb_proteins": ref.get("proteins"), "vogdb_species": ref.get("species"),
                       "proteins": proteins, "genomes": genomes, "hits": hits, "annotated": bool(ref)}
    return out


_UNIPROT_HEADER = re.compile(r"^(?:sp|tr)\|[^|\s]+\|\S+\s+")


def vog_description(text: str) -> str:
    """VOGdb consensus descriptions carry a source prefix (``REFSEQ ...``, ``sp|ACC|ENTRY ...``); drop it."""
    text = text[7:] if text.startswith("REFSEQ ") else text
    return _UNIPROT_HEADER.sub("", text).strip()


def _notable(store) -> dict[str, list[dict[str, Any]]]:
    giants = store.execute("""
        WITH top AS (SELECT protein_id, bin_id, sequence_length FROM proteins
                     ORDER BY sequence_length DESC, protein_id LIMIT 150)
        SELECT top.protein_id, top.bin_id, top.sequence_length, COUNT(a.protein_id)
        FROM top LEFT JOIN annotations a USING (protein_id)
        GROUP BY ALL ORDER BY top.sequence_length DESC, top.protein_id""")
    dark = store.execute("""
        SELECT p.protein_id, p.bin_id, p.sequence_length FROM proteins p
        ANTI JOIN annotations a USING (protein_id)
        ORDER BY p.sequence_length DESC, p.protein_id LIMIT 150""")
    repeats = store.execute("""
        WITH d AS (
            SELECT protein_id, COALESCE(NULLIF(name, ''), accession) AS domain, start_aa,
                   -- annotation_id breaks ties between hits sharing a start, so runs are deterministic
                   ROW_NUMBER() OVER (PARTITION BY protein_id ORDER BY start_aa, accession, annotation_id)
                 - ROW_NUMBER() OVER (PARTITION BY protein_id, COALESCE(NULLIF(name, ''), accession)
                                      ORDER BY start_aa, accession, annotation_id) AS island
            FROM annotations WHERE LOWER(source) = 'pfam' AND start_aa IS NOT NULL
        ), runs AS (
            SELECT protein_id, domain, COUNT(*) AS run FROM d GROUP BY protein_id, domain, island
        )
        , best AS (
            SELECT protein_id, domain, run,
                   ROW_NUMBER() OVER (PARTITION BY protein_id ORDER BY run DESC, domain) AS pick
            FROM runs
        )
        SELECT r.protein_id, p.bin_id, p.sequence_length, r.domain, r.run
        FROM best r JOIN proteins p USING (protein_id)
        WHERE r.pick = 1 AND r.run >= 6
        ORDER BY r.run DESC, p.sequence_length DESC, r.protein_id LIMIT 150""")
    return {
        "giants": [{"protein_id": p, "bin_id": b, "length": n, "hits": h} for p, b, n, h in giants],
        "dark": [{"protein_id": p, "bin_id": b, "length": n} for p, b, n in dark],
        "repeats": [{"protein_id": p, "bin_id": b, "length": n, "domain": d, "run": r} for p, b, n, d, r in repeats],
    }
