"""CRISPR-Cas subtyping: a port of CRISPRCasTyper's typing logic onto Sharur data.

CRISPRCasTyper (Russel et al. 2020, The CRISPR Journal 3:462; MIT) types Cas
operons by searching ~700 Cas profiles, clustering hits into operons, and
scoring each operon against a subtype table. It also classifies CRISPR
repeats with an XGBoost model and joins operons to arrays within 10 kb. This
module reproduces those steps on Sharur's own gene calls and MinCED arrays,
so calls attach to stored protein and locus IDs:

1. :func:`search_genome` runs the profiles against one genome's proteins with
   pyhmmer. E-values are per genome, as in a per-assembly cctyper run, and
   per (profile, protein) hits are merged into union coverages.
2. :func:`filter_hits` keeps each protein's best profile, then applies the
   effector-specific cutoffs (``cutoffs.tab``) and the overall E-value and
   coverage thresholds.
3. :func:`cluster_operons` joins hits separated by at most ``dist`` other
   genes; :func:`type_operon` scores each operon against ``CasScoring.csv``
   and reports interference and adaptation completeness.
4. :func:`type_repeats` predicts an array subtype from its consensus repeat;
   :func:`array_trust` applies cctyper's repeat/spacer identity rules.
5. :func:`link_arrays` joins operons and arrays and reconciles their subtypes.

Data comes from Aksha's ``CCTyper`` database (``aksha`` installs the profiles
with the scoring tables and repeat model beside them). Steps that need the
contig sequence (cctyper's search for known repeats near orphan operons) and
circular replicons are outside this port.
"""

from __future__ import annotations

import json
import math
import os
import re
import statistics
from dataclasses import dataclass, field
from itertools import combinations, product
from pathlib import Path
from typing import Any, Iterable

CCTYPER_VERSION = "1.9.0"
CALLER = f"cctyper-port {CCTYPER_VERSION}"

# cctyper's command-line defaults
DEFAULTS = {
    "dist": 3,             # unknown genes allowed between Cas genes of one operon
    "overall_eval": 0.01,
    "overall_cov_seq": 0.3,
    "overall_cov_hmm": 0.3,
    "ccd": 10000,          # bp between an operon and an array it is linked to
    "pred_prob": 0.75,     # repeat-subtype probability for a named array prediction
    "repeat_id": 70,       # trusted arrays: mean repeat identity above this (%)
    "spacer_id": 55,       # ... mean spacer identity below this (%)
    "spacer_sem": 3.5,     # ... spacer length SEM below this
}
PROFILES = "cctyper_profiles.hmm"
TABLES = ("CasScoring.csv", "cutoffs.tab", "interference.json", "adaptation.json")


def gene_name(hmm: str) -> str:
    """Gene a profile belongs to: ``Cas1_0_IA`` -> ``Cas1``."""
    return hmm.split("_", 1)[0]


def find_database(explicit: str | Path | None = None) -> Path:
    """Directory holding the CCTyper profiles and tables.

    Order: ``explicit``, ``$SHARUR_CCTYPER_DB``, Aksha's install record, then
    Aksha's default store ``~/.config/Astra/CCTyper``.
    """
    candidates: list[Path] = []
    if explicit:
        candidates.append(Path(explicit).expanduser())
    if os.environ.get("SHARUR_CCTYPER_DB"):
        candidates.append(Path(os.environ["SHARUR_CCTYPER_DB"]).expanduser())
    for config in (Path.home() / ".config" / "Aksha" / "hmm_databases.json",
                   Path.home() / "Library" / "Application Support" / "Aksha" / "hmm_databases.json"):
        if config.exists():
            try:
                parsed = json.loads(config.read_text())
            except ValueError:
                continue
            for db in parsed.get("db_urls", []):
                if db.get("name") == "CCTyper" and db.get("installation_dir"):
                    candidates.append(Path(os.path.expandvars(db["installation_dir"])).expanduser())
            root = os.path.expandvars(os.path.expanduser(parsed.get("db_path", "")))
            if root:
                candidates.append(Path(root) / "CCTyper")
    # Aksha's default database store (shared with Astra installs)
    candidates.append(Path.home() / ".config" / "Astra" / "CCTyper")
    for path in candidates:
        if (path / PROFILES).exists() and all((path / t).exists() for t in TABLES):
            return path
    raise FileNotFoundError(
        "CCTyper database not found. Install it with `aksha` (database 'CCTyper') "
        "or point SHARUR_CCTYPER_DB at a directory with "
        f"{PROFILES} and {', '.join(TABLES)}.")


@dataclass
class TypingData:
    """CasScoring table, cutoffs and gene sets from a CCTyper database."""

    root: Path
    types: list[str]
    scores: dict[str, dict[str, float]]          # profile -> {type: score}
    cutoffs: list[tuple[str, float, float, float]]  # (key, evalue, cov_seq, cov_hmm), file order
    interference: dict[str, list[list[str]]]
    adaptation: dict[str, list[list[str]]]
    signature: set[str] = field(default_factory=set)
    single_effector: set[str] = field(default_factory=set)

    @classmethod
    def load(cls, root: str | Path) -> "TypingData":
        root = Path(root)
        lines = (root / "CasScoring.csv").read_text().splitlines()
        types = lines[0].split(",")[1:]
        scores: dict[str, dict[str, float]] = {}
        for line in lines[1:]:
            cells = line.split(",")
            scores[cells[0]] = {t: float(v) for t, v in zip(types, cells[1:]) if v != ""}
        cutoffs = []
        for line in (root / "cutoffs.tab").read_text().splitlines():
            if line.strip():
                key, values = line.split(":", 1)
                ev, cs, ch = (float(x) for x in values.split(","))
                cutoffs.append((key, ev, cs, ch))
        data = cls(root, types, scores, cutoffs,
                   json.loads((root / "interference.json").read_text()),
                   json.loads((root / "adaptation.json").read_text()))
        # Single-gene signatures: profiles with a specific cutoff, and the types
        # those profiles score for. cctyper derives both from the profiles that
        # hit in the current input; the full table is the run-independent form.
        specific = [h for h in scores if data.specific_key(h) is not None]
        data.signature = {gene_name(h) for h in specific}
        data.single_effector = {t for t in types
                                if sum(max(scores[h].get(t, 0.0), 0.0) for h in specific) > 0}
        return data

    def specific_key(self, hmm: str) -> str | None:
        low = hmm.lower()
        for key, *_ in self.cutoffs:
            if key.lower() in low:
                return key
        return None


# ---------------------------------------------------------------------------
# 1. profile search
# ---------------------------------------------------------------------------

def merge_domains(domains: list[dict[str, Any]], tlen: int, qlen: int) -> dict[str, Any]:
    """One (profile, protein) hit from its domains: best domain's E-value and
    score, coverage as the union of aligned target and profile positions."""
    seq_span: set[int] = set()
    hmm_span: set[int] = set()
    for d in domains:
        seq_span.update(range(d["ali_from"], d["ali_to"] + 1))
        if d["hmm_from"] and d["hmm_to"]:
            hmm_span.update(range(d["hmm_from"], d["hmm_to"] + 1))
    best = max(domains, key=lambda d: d["score"])
    return {"evalue": best["evalue"], "score": best["score"],
            "cov_seq": len(seq_span) / tlen if tlen else 0.0,
            "cov_hmm": len(hmm_span) / qlen if qlen else 0.0}


_HMMS: list | None = None


def _load_hmms(path: str) -> list:
    global _HMMS
    if _HMMS is None:
        import pyhmmer
        with pyhmmer.plan7.HMMFile(path) as handle:
            _HMMS = list(handle)
    return _HMMS


def search_genome(proteins: list[tuple[str, str]], hmm_path: str | Path, cpus: int = 1) -> list[dict[str, Any]]:
    """Every included (profile, protein) hit for one genome's proteins."""
    import pyhmmer

    alphabet = pyhmmer.easel.Alphabet.amino()
    seqs = []
    for pid, seq in proteins:
        if seq:
            text = pyhmmer.easel.TextSequence(name=pid.encode(), sequence=seq.rstrip("*").replace("*", "X"))
            seqs.append(text.digitize(alphabet))
    if not seqs:
        return []
    lengths = {s.name: len(s) for s in seqs}
    block = pyhmmer.easel.DigitalSequenceBlock(alphabet, seqs)
    hits_out = []
    for hits in pyhmmer.hmmsearch(_load_hmms(str(hmm_path)), block, cpus=cpus):
        hmm = hits.query.name
        hmm = hmm.decode() if isinstance(hmm, bytes) else hmm
        qlen = hits.query.M
        for hit in hits.included:
            name = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
            domains = []
            for dom in hit.domains.included:
                aln = dom.alignment
                domains.append({"evalue": dom.i_evalue, "score": dom.score,
                                "ali_from": aln.target_from, "ali_to": aln.target_to,
                                "hmm_from": aln.hmm_from, "hmm_to": aln.hmm_to})
            if domains:
                tlen = lengths[hit.name]
                hits_out.append({"hmm": hmm, "protein_id": name, **merge_domains(domains, tlen, qlen)})
    return hits_out


# ---------------------------------------------------------------------------
# 2. thresholds
# ---------------------------------------------------------------------------

def filter_hits(hits: list[dict[str, Any]], data: TypingData, *, overall_eval: float = DEFAULTS["overall_eval"],
                overall_cov_seq: float = DEFAULTS["overall_cov_seq"],
                overall_cov_hmm: float = DEFAULTS["overall_cov_hmm"]) -> list[dict[str, Any]]:
    """Each protein's best-scoring profile, kept if it passes its thresholds.

    The best profile is chosen before thresholds, as in cctyper: a protein
    whose best profile fails is dropped even when a weaker profile would pass.
    """
    best: dict[str, dict[str, Any]] = {}
    for h in hits:
        cur = best.get(h["protein_id"])
        if cur is None or h["score"] > cur["score"]:
            best[h["protein_id"]] = h
    kept = []
    for h in best.values():
        key = data.specific_key(h["hmm"])
        if key is not None:
            _, ev, cs, ch = next(c for c in data.cutoffs if c[0] == key)
            ok = h["evalue"] < ev and h["cov_seq"] >= cs and h["cov_hmm"] >= ch
        else:
            ok = (h["evalue"] < overall_eval and h["cov_seq"] >= overall_cov_seq
                  and h["cov_hmm"] >= overall_cov_hmm)
        if ok:
            kept.append(h)
    return kept


# ---------------------------------------------------------------------------
# 3. operons and subtypes
# ---------------------------------------------------------------------------

def cluster_operons(positions: Iterable[int], dist: int = DEFAULTS["dist"]) -> list[list[int]]:
    """Group gene ranks on one contig: neighbours with at most ``dist`` genes between them."""
    groups: list[list[int]] = []
    for pos in sorted(set(positions)):
        if groups and pos - groups[-1][-1] - 1 <= dist:
            groups[-1].append(pos)
        else:
            groups.append([pos])
    return groups


def _completeness(best_types: list[str], sets: dict[str, list[list[str]]], hmms: list[str]) -> list[str]:
    out = []
    for t in best_types:
        if t in sets:
            groups = sets[t]
            found = sum(1 for g in groups if any(x in h for x in g for h in hmms))
            out.append(f"{round(found / len(groups) * 100)}%")
        else:
            out.append("NA")
    return out


def _strand(rows: list[dict[str, Any]], genes: list[str]) -> Any:
    strands = {r["strand"] for r in rows if any(g in r["hmm"] for g in genes)}
    if not strands:
        return "NA"
    return strands.pop() if len(strands) == 1 else 0


def type_operon(rows: list[dict[str, Any]], data: TypingData) -> dict[str, Any]:
    """Subtype one operon: ``rows`` are its filtered hits with ``hmm``,
    ``score``, ``start``, ``end``, ``strand`` (+1/-1) and gene ``pos``.

    Only profiles in the scoring table count, as in cctyper's inner merge.
    """
    rows = sorted((r for r in rows if r["hmm"] in data.scores), key=lambda r: r["pos"])
    if not rows:
        return {}
    # one row per gene name: the best-scoring profile of that gene
    per_gene: dict[str, dict[str, Any]] = {}
    for r in sorted(rows, key=lambda r: -r["score"]):
        per_gene.setdefault(gene_name(r["hmm"]), r)
    type_scores = {t: 0.0 for t in data.types}
    for r in per_gene.values():
        for t, v in data.scores[r["hmm"]].items():
            type_scores[t] += v
    best_score = max(type_scores.values())
    best_type: Any = next(t for t in data.types if type_scores[t] == best_score)
    ties = [t for t in data.types if type_scores[t] == best_score]
    genes = list(per_gene)
    n = len(per_gene)

    if n >= 3:
        if best_score <= 5:
            if any(g in data.signature for g in genes):
                if len(ties) > 1:
                    prediction, best_type = "Ambiguous", ties
                else:
                    prediction = best_type
            else:
                prediction = "False"
        elif n >= 6:
            # adjacent systems: types with a specific (>=3) profile whose
            # type-unique genes sum to at least 6 make a hybrid
            matrix = {g: data.scores[r["hmm"]] for g, r in per_gene.items()}
            strong = [t for t in data.types if any(matrix[g].get(t, 0) >= 3 for g in matrix)]
            unique_genes = [g for g in matrix if sum(1 for t in strong if matrix[g].get(t, 0) > 0) == 1]
            sums = {t: sum(max(matrix[g].get(t, 0), 0) for g in unique_genes) for t in strong}
            hybrid = [t for t in strong if sums[t] >= 6]
            if len(hybrid) > 1:
                prediction, best_type = f"Hybrid({','.join(hybrid)})", hybrid
            else:
                prediction = best_type
        elif len(ties) > 1:
            prediction, best_type = "Ambiguous", ties
        else:
            prediction = best_type
    else:
        if any(g in data.signature for g in genes):
            if len(ties) > 1:
                prediction, best_type = "Ambiguous", ties
            elif best_type not in data.single_effector:
                prediction = "Ambiguous"
                best_type = [t for t in data.types if type_scores[t] > 0]
            else:
                prediction = best_type
        else:
            prediction = "False"

    best_list = [best_type] if isinstance(best_type, str) else list(best_type)
    hmms = [r["hmm"] for r in rows]
    interf = _completeness(best_list, data.interference, hmms)
    adapt = _completeness(best_list, data.adaptation, hmms)
    interf_genes = [g for t in best_list for grp in data.interference.get(t, []) for g in grp]
    adapt_genes = [g for t in best_list for grp in data.adaptation.get(t, []) for g in grp]
    return {
        "start": min(min(r["start"], r["end"]) for r in rows),
        "end": max(max(r["start"], r["end"]) for r in rows),
        "prediction": prediction,
        "best_type": best_type,
        "best_score": best_score,
        "complete_interference": interf[0] if isinstance(best_type, str) else interf,
        "complete_adaptation": adapt[0] if isinstance(best_type, str) else adapt,
        "strand_interference": _strand(rows, interf_genes),
        "strand_adaptation": _strand(rows, adapt_genes),
        "genes": rows,
    }


def type_contig(hits: list[dict[str, Any]], data: TypingData, dist: int = DEFAULTS["dist"]) -> list[dict[str, Any]]:
    """Operons on one contig from its filtered hits (each with gene ``pos``)."""
    by_pos: dict[int, list[dict[str, Any]]] = {}
    for h in hits:
        by_pos.setdefault(h["pos"], []).append(h)
    operons = []
    for group in cluster_operons(by_pos, dist):
        typed = type_operon([h for p in group for h in by_pos[p]], data)
        if typed:
            operons.append(typed)
    return operons


# ---------------------------------------------------------------------------
# 4. arrays
# ---------------------------------------------------------------------------

_COMP = str.maketrans("ACGT", "TGCA")


def canonical_kmers(k: int = 4) -> list[str]:
    kmers = ["".join(p) for p in product("ACGT", repeat=k)]
    return sorted({min(x, x.translate(_COMP)[::-1]) for x in kmers})


def repeat_features(repeat: str, k: int = 4) -> dict[str, float]:
    """cctyper's repeat features: canonical k-mer counts, length and GC."""
    seq = repeat.upper()
    counts = dict.fromkeys(canonical_kmers(k), 0.0)
    for i in range(len(seq) - k + 1):
        fwd = seq[i:i + k]
        rev = fwd.translate(_COMP)[::-1]
        kmer = fwd if fwd < rev else rev
        if kmer in counts:
            counts[kmer] += 1
    counts["Length"] = float(len(seq))
    counts["GC"] = (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0
    return counts


def type_repeats(repeats: list[str], root: str | Path) -> list[tuple[str, float]]:
    """(subtype, probability) per consensus repeat, from cctyper's XGBoost model."""
    if not repeats:
        return []
    import numpy as np
    import xgboost as xgb

    root = Path(root)
    model = next((root / n for n in ("xgb_repeats.json", "xgb_repeats.ubj") if (root / n).exists()), None)
    if model is None:
        raise FileNotFoundError(f"no xgb_repeats.json in {root}")
    bst = xgb.Booster()
    bst.load_model(str(model))
    labels = {}
    for line in (root / "type_dict.tab").read_text().splitlines():
        if line.strip():
            name, idx = line.split(":")
            labels[int(idx)] = name
    rows = [repeat_features(r) for r in repeats]
    columns = sorted(rows[0])
    matrix = np.array([[row[c] for c in columns] for row in rows], dtype=float)
    best = bst.attr("best_iteration")
    pred = bst.predict(xgb.DMatrix(matrix, feature_names=columns),
                       iteration_range=(0, int(best)) if best is not None else (0, 0))
    return [(labels[int(p.argmax())], float(p.max())) for p in pred]


def _identity(a: str, b: str) -> float:
    """cctyper's pairwise identity: LCS length over the longer sequence (%)."""
    if not a or not b:
        return 0.0
    prev = [0] * (len(b) + 1)
    for ca in a:
        curr = [0]
        for j, cb in enumerate(b, start=1):
            curr.append(prev[j - 1] + 1 if ca == cb else max(prev[j], curr[-1]))
        prev = curr
    return 100.0 * prev[-1] / max(len(a), len(b))


def consensus_repeat(repeats: list[str]) -> str:
    """cctyper's consensus: the most frequent repeat (ties: lexicographic first)."""
    return max(sorted(set(repeats)), key=repeats.count) if repeats else ""


def _mean_identity(seqs: list[str], sample: int = 10) -> float:
    picked = seqs if len(seqs) <= sample else seqs[:: max(1, len(seqs) // sample)][:sample]
    pairs = list(combinations(picked, 2))
    return statistics.mean(_identity(a, b) for a, b in pairs) if pairs else 100.0


def array_trust(repeats: list[str], spacers: list[str], *, repeat_id: float = DEFAULTS["repeat_id"],
                spacer_id: float = DEFAULTS["spacer_id"], spacer_sem: float = DEFAULTS["spacer_sem"]) -> dict[str, Any]:
    """cctyper's array statistics and its 'trusted' rule (conserved repeats,
    diverse spacers of even length)."""
    if len(spacers) > 1:
        s_ident = _mean_identity(spacers)
        s_sem = statistics.stdev(map(len, spacers)) / math.sqrt(len(spacers))
    else:
        s_ident, s_sem = 0.0, 0.0
    r_ident = _mean_identity(repeats)
    return {"repeat_identity": round(r_ident, 1), "spacer_identity": round(s_ident, 1),
            "spacer_sem": round(s_sem, 1),
            "trusted": r_ident > repeat_id and s_ident < spacer_id and s_sem < spacer_sem}


# ---------------------------------------------------------------------------
# 5. operon-array linking
# ---------------------------------------------------------------------------

def interval_distance(a: tuple[int, int], b: tuple[int, int]) -> int:
    """Gap in bp between two intervals; negative when they overlap."""
    return b[0] - a[1] if b[0] > a[1] else a[0] - b[1]


def _major(subtype: str) -> str:
    return re.sub("-.*$", "", subtype)


def reconcile(prediction_cas: str, best_cas: Any, prediction_crispr: str) -> str:
    """cctyper's joint call from an operon and its nearest array."""
    if prediction_cas == prediction_crispr:
        return prediction_cas
    if prediction_cas == "Ambiguous":
        best = list(best_cas) if not isinstance(best_cas, str) else [best_cas]
        if prediction_crispr in best:
            return prediction_crispr
        if _major(prediction_crispr) in [_major(x) for x in best]:
            return _major(prediction_crispr)
        return "Unknown"
    if prediction_cas not in ("False", "Partial"):
        return prediction_cas
    if best_cas == prediction_crispr:
        return f"{best_cas}(Putative)"
    if isinstance(best_cas, str) and _major(best_cas) == _major(prediction_crispr):
        return f"{_major(prediction_crispr)}(Putative)"
    return "Unknown"


def link_arrays(operons: list[dict[str, Any]], arrays: list[dict[str, Any]],
                ccd: int = DEFAULTS["ccd"]) -> None:
    """Attach arrays within ``ccd`` bp to each non-False operon on the same
    contig, set the joint ``prediction`` and ``status``, and mark linked arrays
    trusted. Operons and arrays carry ``contig_id``, ``start``, ``end``; arrays
    carry ``locus_id``, ``prediction`` and ``trusted``."""
    for op in operons:
        op["arrays"], op["distances"] = [], []
        if op["prediction"] != "False" and op.get("positioned", True):
            for arr in arrays:
                if arr["contig_id"] != op["contig_id"]:
                    continue
                d = interval_distance((op["start"], op["end"]), (arr["start"], arr["end"]))
                if d <= ccd:
                    op["arrays"].append(arr)
                    op["distances"].append(max(d, 0))
        if op["arrays"]:
            nearest = op["arrays"][op["distances"].index(min(op["distances"]))]
            op["joint_prediction"] = reconcile(op["prediction"], op["best_type"], nearest["prediction"])
            putative = "Unknown" in op["joint_prediction"] or "Putative" in op["joint_prediction"]
            op["status"] = "crispr_cas_putative" if putative else "crispr_cas"
            for arr in op["arrays"]:
                arr["trusted"] = True
                arr["near_cas"] = True
        else:
            op["joint_prediction"] = op["prediction"]
            op["status"] = "cas_putative" if op["prediction"] in ("False", "Ambiguous") else "cas"


# ---------------------------------------------------------------------------
# dataset runner
# ---------------------------------------------------------------------------

def _search_worker(args: tuple[str, list[tuple[str, str]], str]) -> tuple[str, list[dict[str, Any]]]:
    bin_id, proteins, hmm_path = args
    return bin_id, search_genome(proteins, hmm_path, cpus=1)


def _minced_reports(dataset_dir: Path) -> dict[str, Path]:
    reports: dict[str, Path] = {}
    for d in sorted(dataset_dir.glob("stage05c*")):
        for f in d.glob("*_crispr.txt"):
            reports.setdefault(f.name[: -len("_crispr.txt")], f)
    return reports


def load_arrays(conn, dataset_dir: Path) -> list[dict[str, Any]]:
    """CRISPR array loci with consensus repeat and, when MinCED's report is
    present, every repeat and spacer."""
    from sharur.crispr import parse_minced_text

    rows = conn.execute(
        """SELECT l.locus_id, l.contig_id, l.start, l.end_coord, l.metadata, c.bin_id
           FROM loci l JOIN contigs c USING (contig_id)
           WHERE LOWER(l.locus_type) LIKE '%crispr%'""").fetchall()
    reports = _minced_reports(dataset_dir)
    parsed: dict[str, dict[tuple[str, int], dict[str, Any]]] = {}
    arrays = []
    for locus_id, contig, start, end, metadata, bin_id in rows:
        meta = json.loads(metadata) if isinstance(metadata, str) else (metadata or {})
        inner = meta.get("metadata") or {}
        repeats, spacers = [r.get("seq", "") for r in meta.get("repeats", [])], \
            [s.get("seq", "") for s in meta.get("spacers", [])]
        if not repeats and bin_id in reports:
            if bin_id not in parsed:
                parsed[bin_id] = {(a["contig"], a["start"]): a for a in parse_minced_text(reports[bin_id])}
            hit = parsed[bin_id].get((contig, start))
            if hit:
                repeats = [r["seq"] for r in hit["repeats"]]
                spacers = [s["seq"] for s in hit["spacers"]]
        consensus = consensus_repeat(repeats) or (inner.get("rpt_unit_seq") or "").upper()
        arrays.append({"locus_id": locus_id, "contig_id": contig, "genome_id": bin_id, "start": start,
                       "end": end, "repeats": repeats, "spacers": spacers, "consensus": consensus})
    return arrays


def type_arrays(arrays: list[dict[str, Any]], root: Path, pred_prob: float = DEFAULTS["pred_prob"]) -> None:
    """Repeat subtype and trust for each array, in place."""
    typed = type_repeats([a["consensus"] for a in arrays if a["consensus"]], root)
    it = iter(typed)
    for a in arrays:
        subtype, prob = next(it) if a["consensus"] else ("Unknown", 0.0)
        a["subtype"], a["probability"] = subtype, prob
        a["prediction"] = subtype if prob >= pred_prob else "Unknown"
        stats = array_trust(a["repeats"], a["spacers"]) if len(a["repeats"]) > 1 and a["spacers"] else {
            "repeat_identity": None, "spacer_identity": None, "spacer_sem": None, "trusted": False}
        a.update(stats)
        # cctyper trusts arrays whose repeat matches a known subtype confidently
        if prob >= 0.9:
            a["trusted"] = True
        a["near_cas"] = False


def run_dataset(db_path: str | Path, *, cctyper_db: str | Path | None = None, workers: int | None = None,
                genomes: list[str] | None = None, progress=None) -> dict[str, Any]:
    """Type every genome's Cas operons and arrays (read-only on the dataset).

    Returns ``{"systems": [...], "arrays": [...], "genomes": n}``; write them
    with :func:`write_results`.
    """
    from concurrent.futures import ProcessPoolExecutor

    import duckdb

    root = find_database(cctyper_db)
    data = TypingData.load(root)
    db_path = Path(db_path)
    conn = duckdb.connect(str(db_path), read_only=True)
    try:
        where, params = ("WHERE bin_id IN (SELECT UNNEST(?::VARCHAR[]))", [genomes]) if genomes else ("", [])
        # Proteins without a genomic position (a contig shared across genomes, or
        # stacked at coordinate 0) are typed alone, as cctyper treats protein-only input.
        rows = conn.execute(
            f"""WITH p AS (
                    SELECT *, COUNT(DISTINCT bin_id) OVER (PARTITION BY contig_id) AS n_bins,
                           COUNT(*) OVER (PARTITION BY contig_id, start, end_coord) AS n_same
                    FROM proteins {where})
                SELECT bin_id, protein_id, contig_id, start, end_coord, strand, sequence,
                       CASE WHEN n_bins > 1 OR (start = 0 AND n_same > 1) THEN protein_id ELSE contig_id END AS unit,
                       ROW_NUMBER() OVER (PARTITION BY bin_id, contig_id ORDER BY start, end_coord, protein_id) AS pos
                FROM p ORDER BY bin_id""", params).fetchall()
        arrays = load_arrays(conn, db_path.parent)
    finally:
        conn.close()
    if genomes:
        arrays = [a for a in arrays if a["genome_id"] in set(genomes)]

    by_bin: dict[str, list[tuple[str, str]]] = {}
    gene: dict[str, tuple[str, str, str, int, int, int, int]] = {}
    for bin_id, pid, contig, start, end, strand, seq, unit, pos in rows:
        by_bin.setdefault(bin_id, []).append((pid, seq))
        gene[pid] = (bin_id, contig, unit, start, end, 1 if strand in ("+", "1") else -1, pos)
    del rows

    jobs = [(b, prots, str(root / PROFILES)) for b, prots in by_bin.items()]
    workers = workers or os.cpu_count() or 1
    systems: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for done, (bin_id, hits) in enumerate(pool.map(_search_worker, jobs, chunksize=4), 1):
            kept = filter_hits(hits, data)
            per_unit: dict[tuple[str, str], list[dict[str, Any]]] = {}
            for h in kept:
                _, contig, unit, start, end, strand, pos = gene[h["protein_id"]]
                per_unit.setdefault((contig, unit), []).append(
                    {**h, "start": start, "end": end, "strand": strand, "pos": pos})
            for (contig, unit), unit_hits in sorted(per_unit.items()):
                for op in type_contig(unit_hits, data):
                    op.update(genome_id=bin_id, contig_id=contig, positioned=unit == contig)
                    systems.append(op)
            if progress:
                progress(done, len(jobs))
    type_arrays(arrays, root)
    link_arrays(systems, arrays)
    return {"systems": systems, "arrays": arrays, "genomes": len(jobs)}


def _fmt(value: Any) -> str:
    return ",".join(map(str, value)) if isinstance(value, (list, tuple)) else str(value)


def system_rows(systems: list[dict[str, Any]]) -> list[tuple]:
    rows = []
    for op in sorted(systems, key=lambda o: (o["genome_id"], o["contig_id"], o["start"])):
        where = (f"{op['contig_id']}:{op['start']}-{op['end']}" if op.get("positioned", True)
                 else op["genes"][0]["protein_id"])
        sid = f"cctyper:{op['genome_id']}:{where}"
        op["system_id"] = sid
        rows.append((sid, op["genome_id"], op["contig_id"], int(op["start"]), int(op["end"]), op["status"],
                     op["joint_prediction"], op["prediction"], _fmt(op["best_type"]), float(op["best_score"]),
                     _fmt(op["complete_interference"]), _fmt(op["complete_adaptation"]),
                     str(op["strand_interference"]), str(op["strand_adaptation"]), len(op["genes"]),
                     ",".join(g["protein_id"] for g in op["genes"]), ",".join(g["hmm"] for g in op["genes"]),
                     ",".join(a["locus_id"] for a in op["arrays"]), _fmt(op["distances"]), CALLER))
    return rows


def write_results(db_path: str | Path, result: dict[str, Any]) -> dict[str, int]:
    """Replace this caller's rows in a dataset already at the current schema."""
    import duckdb

    from sharur.storage.migrations import get_current_version, run_migrations
    from sharur.storage.schema import SCHEMA_VERSION

    conn = duckdb.connect(str(db_path))
    try:
        version = get_current_version(conn)
        if version < SCHEMA_VERSION - 1:
            raise RuntimeError(f"{db_path} is at schema {version}; run `sharur migrate` (and review what "
                               f"it changes) before writing CRISPR-Cas calls.")
        if version < SCHEMA_VERSION:
            run_migrations(conn)   # the one pending step only adds the CRISPR-Cas tables
        sys_rows = system_rows(result["systems"])
        arr_rows = [(a["locus_id"], a["genome_id"], a["contig_id"], len(a["consensus"]), len(a["repeats"]) or None,
                     a["subtype"], a["probability"], a["prediction"], a["repeat_identity"], a["spacer_identity"],
                     a["spacer_sem"], bool(a["trusted"]), bool(a["near_cas"]), CALLER) for a in result["arrays"]]
        members = [(op["system_id"], g["protein_id"], "cctyper", i, g["hmm"], float(g["score"]))
                   for op in result["systems"] for i, g in enumerate(op["genes"])]
        conn.execute("BEGIN TRANSACTION")
        conn.execute("DELETE FROM crispr_cas_systems WHERE caller LIKE 'cctyper-port%'")
        conn.execute("DELETE FROM crispr_array_types WHERE caller LIKE 'cctyper-port%'")
        conn.execute("DELETE FROM system_proteins WHERE system_source = 'cctyper'")
        if sys_rows:
            conn.executemany("INSERT INTO crispr_cas_systems (system_id, genome_id, contig_id, start, end_coord, "
                             "status, prediction, prediction_cas, best_type, best_score, complete_interference, "
                             "complete_adaptation, strand_interference, strand_adaptation, genes_count, protein_ids, "
                             "profile_names, crispr_locus_ids, crispr_distances, caller) "
                             "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", sys_rows)
        if arr_rows:
            conn.executemany("INSERT INTO crispr_array_types VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", arr_rows)
        if members:
            conn.executemany("INSERT INTO system_proteins (system_id, protein_id, system_source, position, "
                             "profile_name, score) VALUES (?, ?, ?, ?, ?, ?)", members)
        conn.commit()
    finally:
        conn.close()
    return {"systems": len(sys_rows), "arrays": len(arr_rows), "members": len(members)}
