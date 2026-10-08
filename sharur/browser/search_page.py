"""The search results page: matches grouped by kind, the best match previewed, near misses offered.

An index of everything nameable in the dataset (taxa, genomes, function labels,
KEGG modules, KEGG orthologs with their KEGG names and EC numbers, Pfam and VOG
families, curated system types) is built once the catalog is ready and kept.
Each query scores entries by how closely the name matches (exact, prefix,
word, substring, description) and breaks ties by how many genomes carry them.
Proteins match by identifier straight from the store.
"""

from __future__ import annotations

import difflib
import re
import threading
from collections import Counter
from dataclasses import dataclass, field
from typing import Any
from urllib.parse import quote

import numpy as np
from markupsafe import Markup, escape

from sharur.browser.catalog import RANKS, UNCLASSIFIED
from sharur.predicates.vocabulary import PREDICATE_BY_ID

KIND_ORDER = ("taxon", "genome", "function", "pathway", "ko", "pfam", "vog", "system", "term", "protein")
KIND_LABELS = {"taxon": "Taxa", "genome": "Genomes", "function": "Function labels", "pathway": "Pathways",
               "ko": "KEGG orthologs", "pfam": "Pfam families", "vog": "VOG families", "system": "Systems",
               "term": "V2 semantic terms", "protein": "Proteins"}
KIND_TAGS = {"taxon": "taxon", "genome": "genome", "function": "function", "pathway": "pathway", "ko": "KO",
             "pfam": "Pfam", "vog": "VOG", "system": "system", "term": "V2 term", "protein": "protein"}
KO_RE = re.compile(r"^K\d{5}$", re.I)
EC_RE = re.compile(r"^(?:EC[:\s]*)?(\d+\.(?:\d+|-)\.(?:\d+|-)\.(?:n?\d+|-))$", re.I)
EC_IN_TEXT = re.compile(r"\[EC:([^\]]+)\]")
SHOWN = 60          # rows rendered per group
PROTEINS = 50


@dataclass
class Entry:
    kind: str
    id: str
    label: str
    sub: str
    url: str
    name_l: str                 # lowercased label + id, for name matches
    extra_l: str = ""           # lowercased description, for weaker matches
    weight: int = 0             # genomes carrying it (tie-break)
    share: float | None = None  # share of genomes carrying it
    data: dict[str, Any] = field(default_factory=dict)


def _url(*parts: str) -> str:
    return "/" + "/".join(quote(str(p), safe="") for p in parts)


class SearchIndex:
    """Everything nameable in the dataset, built once the catalog is ready."""

    def __init__(self, ctx) -> None:
        self.ctx = ctx
        self.entries: list[Entry] = []
        self.ec_to_kos: dict[str, list[str]] = {}
        self.by_key: dict[tuple[str, str], Entry] = {}
        self.vocab: list[str] = []
        self.complete = False
        self._lock = threading.Lock()

    def ensure(self) -> "SearchIndex":
        with self._lock:
            if not self.complete:
                self._build()
        return self

    def _add(self, e: Entry) -> None:
        self.entries.append(e)
        self.by_key[(e.kind, e.id)] = e

    def _build(self) -> None:
        ctx, catalog = self.ctx, self.ctx.catalog
        self.entries, self.by_key, self.ec_to_kos = [], {}, {}
        n = len(catalog.genomes) or 1
        ready = catalog.ready.is_set()

        taxa: Counter = Counter()
        for g in catalog.genomes:
            for rank in RANKS:
                name = g.taxonomy.get(rank)
                if name and name != UNCLASSIFIED:
                    taxa[(rank, name)] += 1
        for (rank, name), count in taxa.items():
            self._add(Entry("taxon", f"{rank}:{name}", name, f"{rank} · {count:,} genome{'s' if count != 1 else ''}",
                            _url("taxa", rank, name), name.lower(), "", count, count / n,
                            {"rank": rank, "name": name}))
        for g in catalog.genomes:
            lineage = "; ".join(v for _, v in g.lineage if v != UNCLASSIFIED)
            # genomes match by identifier; their lineage names match through the taxa
            # weights elsewhere count genomes; a genome carries none, so equal names favour clades and families
            self._add(Entry("genome", g.bin_id, g.bin_id, g.label, _url("genome", g.bin_id), g.bin_id.lower(),
                            "", 0, None, {"genome": g, "lineage": lineage}))

        if catalog.predicate_genomes is not None:
            for pred, k in catalog.predicate_index.items():
                d = PREDICATE_BY_ID.get(pred)
                carriers = int(catalog.predicate_genomes[k])
                label = d.name if d else pred
                self._add(Entry("function", pred, label, pred, _url("function", pred), f"{label} {pred}".lower(),
                                (d.description if d else "").lower(), carriers, carriers / n,
                                {"proteins": int(catalog.predicate_proteins[k]) if catalog.predicate_proteins is not None else None}))

        if catalog.modules:
            complete = None
            if catalog.module_completeness is not None:
                complete = (catalog.module_completeness >= 0.75).sum(axis=0)
            index = {m: j for j, m in enumerate(catalog.module_ids)}
            for module, d in catalog.modules.items():
                carriers = int(complete[index[module]]) if complete is not None and module in index else 0
                self._add(Entry("pathway", module, d.name, module, _url("pathway", module),
                                f"{d.name} {module}".lower(), d.module_class.lower(), carriers,
                                carriers / n if complete is not None else None,
                                {"class": d.module_class.split(";")[-1].strip()}))

        sets = getattr(ctx, "feature_sets", None)
        ko_names = getattr(ctx, "ko_names", None)
        if sets is not None:
            try:
                fs = sets.get("ko")
            except LookupError:
                fs = None
            if fs is not None and len(fs.ids):
                carriers = fs.carriers(np.ones(len(catalog.genomes), dtype=bool))
                for k, ko in enumerate(fs.ids):
                    known = ko_names.get(ko) if ko_names else None
                    symbols, definition = known if known else ("", "")
                    symbol = symbols.split(",")[0].strip()
                    plain = EC_IN_TEXT.sub("", definition).strip()
                    for block in EC_IN_TEXT.findall(definition):
                        for ec in block.split():
                            self.ec_to_kos.setdefault(ec, []).append(ko)
                    self._add(Entry("ko", ko, f"{ko} {symbol}".strip(), plain, f"/search?q={ko}",
                                    f"{ko} {symbols}".lower(), plain.lower(), int(carriers[k]), int(carriers[k]) / n,
                                    {"symbols": symbols, "definition": definition}))

        for d in catalog.domains.values():
            self._add(Entry("pfam", d["accession"], d["name"], f'{d["accession"]} · {d["description"]}',
                            _url("domain", d["accession"]), f'{d["name"]} {d["accession"]}'.lower(),
                            d["description"].lower(), d["genomes"], d["genomes"] / n,
                            {"proteins": d["proteins"], "description": d["description"]}))
        for v in catalog.vogs.values():
            desc = v.get("description") or ""
            self._add(Entry("vog", v["accession"], v["name"], desc, _url("vog", v["accession"]),
                            f'{v["name"]} {v["accession"]}'.lower(), desc.lower(), v["genomes"], v["genomes"] / n,
                            {"description": desc}))

        systems: dict[tuple[str, str], set[str]] = {}
        for s in catalog.systems:
            systems.setdefault((s["kind"], str(s["type"])), set()).add(s["bin_id"])
        cctyper = getattr(ctx, "cctyper_systems", None)
        for s in (cctyper() if cctyper else []):
            if s.get("confident"):
                systems.setdefault(("crispr", s["prediction"]), set()).add(s["bin_id"])
        for (kind, system_type), bins in systems.items():
            if kind == "crispr":
                label, url = f"CRISPR-Cas {system_type}", f"/crispr-cas/calls?subtype={quote(system_type, safe='')}"
            else:
                label, url = system_type, _url("system", kind, system_type)
            self._add(Entry("system", f"{kind}:{system_type}", label, f"{kind} system · {len(bins):,} genomes", url,
                            f"{label} {system_type}".lower(), kind, len(bins), len(bins) / n, {"kind": kind}))

        words = set()
        for e in self.entries:
            if e.kind in ("taxon", "function", "pathway", "pfam", "system"):
                words.add(e.label.lower())
            elif e.kind == "ko" and e.data.get("symbols"):
                words.update(s.strip().lower() for s in e.data["symbols"].split(",") if s.strip())
        self.vocab = sorted(words)
        self.complete = ready and bool(catalog.notable)

    # ------------------------------------------------------------------ #

    @staticmethod
    def _score(e: Entry, q: str, word: re.Pattern, tokens: list[str]) -> int:
        label = e.label.lower()
        idl = label if e.kind == "taxon" else e.id.lower()   # a taxon's id carries its rank
        if q == idl or q == label:
            return 100
        if label.startswith(q) or idl.startswith(q):
            return 80
        if word.search(e.name_l):
            return 65
        if q in e.name_l:
            return 50
        if len(tokens) > 1 and q in e.extra_l:              # the whole phrase in the description
            return 40
        if len(tokens) > 1 and all(t in e.name_l for t in tokens):
            return 38
        if len(tokens) > 1 and all(t in e.name_l or t in e.extra_l for t in tokens):
            return 30
        if len(q) >= 3 and q in e.extra_l:
            return 25
        return 0

    def search(self, q: str) -> list[tuple[int, Entry]]:
        ql = q.strip().lower()
        if not ql:
            return []
        word = re.compile(r"(?<![a-z0-9])" + re.escape(ql))
        tokens = [t for t in ql.split() if len(t) >= 2]
        ec = EC_RE.match(q.strip())
        found = []
        if ec:
            # an EC number lists the KOs KEGG assigns it; "1.2.7.-" covers every 1.2.7.x
            target = ec.group(1)
            prefix = target.rstrip("-").rstrip(".") + "." if target.endswith("-") else None
            for ec_id, kos in self.ec_to_kos.items():
                if ec_id == target or (prefix and ec_id.startswith(prefix)):
                    found.extend((90, self.by_key[("ko", ko)]) for ko in kos if ("ko", ko) in self.by_key)
        for e in self.entries:
            s = self._score(e, ql, word, tokens)
            if s:
                found.append((s, e))
        found.sort(key=lambda se: (-se[0], -se[1].weight, se[1].label.lower()))
        seen, out = set(), []
        for s, e in found:
            if (e.kind, e.id) not in seen:
                seen.add((e.kind, e.id))
                out.append((s, e))
        return out

    def near(self, q: str, n: int = 5) -> list[str]:
        ql = q.strip().lower()
        if len(ql) < 3:
            return []
        pool = [w for w in self.vocab if abs(len(w) - len(ql)) <= max(2, len(ql) // 3)]
        return difflib.get_close_matches(ql, pool, n=n, cutoff=0.72)


def index_for(ctx) -> SearchIndex:
    idx = getattr(ctx, "_search_index", None)
    if idx is None:
        idx = SearchIndex(ctx)
        ctx._search_index = idx
    return idx.ensure()


def highlight(text: Any, q: str) -> Markup:
    """Escape ``text`` and mark every occurrence of the query's words."""
    text = "" if text is None else str(text)
    tokens = sorted({t for t in re.split(r"\s+", q.strip()) if len(t) >= 2}, key=len, reverse=True)
    if not tokens:
        return Markup(escape(text))
    pattern = re.compile("|".join(re.escape(t) for t in tokens), re.I)
    out, last = [], 0
    for m in pattern.finditer(text):
        out.append(escape(text[last:m.start()]))
        out.append(Markup("<mark>") + escape(m.group(0)) + Markup("</mark>"))
        last = m.end()
    out.append(escape(text[last:]))
    return Markup("").join(out)


def preview(ctx, e: Entry) -> dict[str, Any]:
    """What the best match shows: a few numbers and the next places to go."""
    catalog = ctx.catalog
    n = len(catalog.genomes) or 1
    out: dict[str, Any] = {"stats": [], "links": [], "text": ""}
    tree = lambda token: f"/tree?features={quote(token, safe=':')}"  # noqa: E731
    if e.kind == "function":
        d = PREDICATE_BY_ID.get(e.id)
        out["text"] = d.description if d else ""
        out["stats"] = [(f"{e.weight:,}", "genomes"), (f"{(e.data.get('proteins') or 0):,}", "proteins")]
        if d:
            out["stats"].append((d.category.replace("_", " "), "system-level" if d.level == "system" else "component-level"))
        out["links"] = [("Function page", e.url), ("Map on tree", tree(f"function:{e.id}")),
                        ("Matrix", f"/matrix?a=all&kind=function&features={quote(e.id)}")]
    elif e.kind == "taxon":
        rank, name = e.data["rank"], e.data["name"]
        genomes = catalog.clade(rank, name)
        child, counts = catalog.children(genomes, rank)
        comp = [g.completeness for g in genomes if g.completeness is not None]
        out["text"] = " › ".join(v for _, v in catalog.lineage_of(rank, name))
        out["stats"] = [(f"{len(genomes):,}", "genomes"), (rank, "rank")]
        if child and counts:
            out["stats"].append((f"{len([c for c in counts if c[0] != UNCLASSIFIED]):,}", f"{child} groups"))
        if comp:
            out["stats"].append((f"{float(np.median(comp)):.0f}%", "median completeness"))
        token = f"{rank}:{name}"
        out["links"] = [("Clade page", e.url), ("Tree", f"/tree?root={quote(token, safe=':')}"),
                        ("Matrix", f"/matrix?a={quote(token, safe=':')}"),
                        ("Functional landscape", f"/landscape?q={quote(name)}")]
    elif e.kind == "genome":
        g = e.data["genome"]
        out["text"] = "; ".join(v for _, v in g.lineage if v != UNCLASSIFIED)
        out["stats"] = [(f"{g.proteins:,}", "proteins"), (f"{g.contigs:,}", "contigs")]
        if g.completeness is not None:
            out["stats"].append((f"{g.completeness:.1f}%", "complete"))
        out["links"] = [("Genome page", e.url), ("Matrix", f"/matrix?a={quote(g.bin_id)}")]
    elif e.kind == "ko":
        from sharur.browser.annotations import link_ec  # noqa: PLC0415

        out["text"] = link_ec(e.data.get("definition") or "")
        out["stats"] = [(f"{e.weight:,}", "genomes"), (f"{e.share:.0%}" if e.share is not None else "–", "of the dataset")]
        try:
            from sharur.predicates.mappings.kegg_map import KEGG_TO_PREDICATES  # noqa: PLC0415

            labels = [p for p in KEGG_TO_PREDICATES.get(e.id, []) if p in PREDICATE_BY_ID][:8]
        except Exception:  # noqa: BLE001 - labels are a nicety
            labels = []
        out["labels"] = [(p, PREDICATE_BY_ID[p].name) for p in labels]
        out["links"] = [("KEGG ↗", f"https://www.kegg.jp/entry/{e.id}"), ("Map on tree", tree(f"ko:{e.id}")),
                        ("Matrix", f"/matrix?a=all&kind=ko&features={e.id}")]
        out["hint"] = f"{e.id} in <clade or genome> lists its proteins there"
    elif e.kind == "pfam":
        out["text"] = e.data.get("description", "")
        out["stats"] = [(f"{e.weight:,}", "genomes"), (f"{(e.data.get('proteins') or 0):,}", "proteins")]
        out["links"] = [("Family page", e.url), ("Map on tree", tree(f"pfam:{e.id}")),
                        ("Matrix", f"/matrix?a=all&kind=pfam&features={e.id}"),
                        ("InterPro ↗", f"https://www.ebi.ac.uk/interpro/entry/pfam/{e.id}/")]
    elif e.kind == "vog":
        out["text"] = e.data.get("description", "")
        out["stats"] = [(f"{e.weight:,}", "genomes")]
        out["links"] = [("Family page", e.url)]
    elif e.kind == "pathway":
        out["text"] = e.data.get("class", "")
        if e.share is not None:
            out["stats"] = [(f"{e.weight:,}", "genomes ≥ 75% complete")]
        out["links"] = [("Pathway page", e.url), ("Steps across genomes", f"/matrix?a=all&module={e.id}&order=cluster"),
                        ("Map on tree", tree(f"module:{e.id}")), ("KEGG ↗", f"https://www.kegg.jp/entry/{e.id}")]
    elif e.kind == "system":
        kind = e.data["kind"]
        out["stats"] = [(f"{e.weight:,}", "genomes"), (kind, "curated caller")]
        out["links"] = [("Calls", e.url), ("Map on tree", tree(f"system:{e.id}"))]
    elif e.kind == "term":
        out["text"] = "V2 term stored in the semantic term table; its proteins and their stored rows are one click away."
        out["stats"] = [(f"{e.weight:,}", "proteins")]
        out["links"] = [("Proteins with this term", e.url)]
    out["share"] = e.share
    return out


def run(ctx, q: str) -> dict[str, Any]:
    """Grouped results, the best match and near misses for one query."""
    idx = index_for(ctx)
    scored = idx.search(q)
    groups: dict[str, list[tuple[int, Entry]]] = {k: [] for k in KIND_ORDER}
    for s, e in scored:
        groups[e.kind].append((s, e))
    stripped = q.strip()
    proteins, protein_count = [], 0
    if len(stripped) >= 3 and " " not in stripped:
        with ctx.lock:
            rows = ctx.store.execute(
                f"SELECT protein_id, bin_id, sequence_length FROM proteins WHERE protein_id ILIKE ? "
                f"ORDER BY protein_id LIMIT {PROTEINS + 1}", [f"%{stripped}%"])
        protein_count = len(rows)
        for pid, bin_id, length in rows[:PROTEINS]:
            g = ctx.catalog.by_bin.get(bin_id)
            score = 100 if pid.lower() == stripped.lower() else 45
            proteins.append((score, Entry("protein", pid, pid, f"{g.label if g else bin_id} · {length or 0:,} aa",
                                          _url("protein", pid), pid.lower())))
    groups["protein"] = proteins
    semantic = getattr(ctx, "semantic", None)
    if semantic is not None and semantic.cheap_catalog and len(stripped) >= 2:
        # V2 terms are their own vocabulary; they are listed beside function labels, never merged with them
        for term, n in semantic.suggest_terms(stripped, SHOWN):
            low = term.lower()
            score = 99 if low == stripped.lower() else 79 if low.startswith(stripped.lower()) else 49
            groups["term"].append((score, Entry("term", term, term, f"{n:,} proteins",
                                                "/terms?has=" + quote(term, safe=":"), low, "", n)))
    out_groups = []
    for kind in KIND_ORDER:
        rows = groups[kind]
        if not rows:
            continue
        count = len(rows) if kind != "protein" else protein_count
        more = None
        if kind == "pfam":
            more = f"/domains?q={quote(stripped)}"
        out_groups.append({"kind": kind, "label": KIND_LABELS[kind], "count": count,
                           "count_label": f"{PROTEINS}+" if kind == "protein" and count > PROTEINS else f"{count:,}",
                           "rows": [e for _, e in rows[:SHOWN]], "more": more if count > SHOWN else None})
    best = None
    candidates = [(s, e) for g in KIND_ORDER for s, e in groups[g][:1]]
    if candidates:
        candidates.sort(key=lambda se: (-se[0], -se[1].weight))
        s, e = candidates[0]
        if s >= 50:
            best = {"entry": e, "score": s, "tag": KIND_TAGS[e.kind], **preview(ctx, e)}
    total = sum(len(groups[k]) for k in KIND_ORDER if k != "protein") + protein_count
    missing_ko = None
    if KO_RE.match(stripped) and ("ko", stripped.upper()) not in idx.by_key:
        known = ctx.ko_names.get(stripped.upper()) if getattr(ctx, "ko_names", None) else None
        missing_ko = {"id": stripped.upper(), "symbols": known[0] if known else "", "definition": known[1] if known else ""}
    return {"groups": out_groups, "best": best, "total": total,
            "near": idx.near(stripped) if not any(s >= 65 for s, _ in scored) and total < 3 else [],
            "missing_ko": missing_ko,
            "indexed": idx.complete}


# --------------------------------------------------------------------------- #
# Scoped search additions
# --------------------------------------------------------------------------- #


def scope_token(kind: str, label: str, bins: list[str]) -> str:
    if kind in RANKS:
        return f"{kind}:{label}"
    if kind == "genome":
        return label
    return "genomes:" + ",".join(bins)


def breakdown(ctx, kind: str, label: str, bins: list[str], per_genome: Counter) -> dict[str, Any]:
    """Matches by genome, and by the next rank down when the scope is a clade."""
    catalog = ctx.catalog
    genomes = [catalog.by_bin[b] for b in bins if b in catalog.by_bin]
    top = per_genome.most_common(10)
    out: dict[str, Any] = {"genomes": [(b, catalog.by_bin[b].label if b in catalog.by_bin else "", n) for b, n in top],
                           "max": top[0][1] if top else 1, "children": [], "child_rank": None}
    if kind in RANKS and len(genomes) > 1:
        child, counts = catalog.children(genomes, kind)
        if child:
            hit = Counter((catalog.by_bin[b].taxonomy.get(child) or UNCLASSIFIED) for b in per_genome if b in catalog.by_bin)
            size = Counter(g.taxonomy.get(child) or UNCLASSIFIED for g in genomes)
            rows = [(name, hit.get(name, 0), n) for name, n in size.most_common()]
            out["children"] = rows[:12]
            out["child_rank"] = child
    return out
