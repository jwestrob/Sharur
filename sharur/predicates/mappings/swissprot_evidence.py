"""Reviewed-protein (UniProtKB/Swiss-Prot) consensus evidence for Pfam and KO predicates.

Each reviewed protein carries the predicates of its curator-assigned EC numbers
and experimentally evidenced GO terms (closed over is_a/part_of ancestors and
matched to the GO anchors in :mod:`sharur.predicates.mappings.pfam_evidence_spec`).
For Pfam consensus a protein also carries the KEGG-evidenced predicates of its
KOs (KEGG's gene -> KO links for the protein's ``DR KEGG`` genes).

Thresholds live beside ``SWISSPROT_COVERAGE`` in
:mod:`sharur.predicates.mappings.pfam_evidence`. Used at build time by
``scripts/build_pfam_predicate_map.py`` and ``scripts/build_kegg_swissprot_consensus.py``,
and by the integrity tests' optional recounts.
"""

from __future__ import annotations

import gzip
import re
from collections import defaultdict
from dataclasses import dataclass
from functools import cache
from pathlib import Path

from sharur.predicates.mappings.kegg_map import get_predicates_for_ec
from sharur.predicates.mappings.pfam_evidence import (
    EVIDENCE,
    SWISSPROT_CODOMAIN_FRACTION,
    SWISSPROT_COVERAGE,
    SWISSPROT_EXCLUDED,
    SWISSPROT_MIN_LOWER_BOUND,
    SWISSPROT_SINGLE,
    wilson_lower_bound,
)
from sharur.predicates.vocabulary import PREDICATE_BY_ID


EXPERIMENTAL_GO = {"EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP"}
KO_CONSENSUS = re.compile(r"swissprot:(\d+)/(\d+)")


@dataclass(frozen=True)
class Reviewed:
    accession: str
    pfams: frozenset[str]
    predicates: frozenset[str]  # from curated EC and experimental GO
    kegg_genes: frozenset[str]


def go_ancestors(obo: Path):
    """Callable GO id -> frozenset of itself and its is_a/part_of ancestors."""
    parents: dict[str, list[str]] = defaultdict(list)
    term = None
    with open(obo) as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line == "[Term]":
                term = None
            elif line.startswith("id: GO:"):
                term = line[4:]
            elif term and line.startswith("is_a: "):
                parents[term].append(line[6:16])
            elif term and line.startswith("relationship: part_of "):
                parents[term].append(line[22:32])

    @cache
    def ancestors(go: str) -> frozenset[str]:
        found = {go}
        for parent in parents.get(go, ()):
            found |= ancestors(parent)
        return frozenset(found)

    return ancestors


def read_swissprot(path: Path, ancestors) -> list[Reviewed]:
    """Every reviewed protein with a Pfam or KEGG cross-reference."""
    anchor_preds: dict[str, set[str]] = defaultdict(set)
    for pred, spec in EVIDENCE.items():
        for go in spec["go"]:
            anchor_preds[go].add(pred)
    vocabulary = set(PREDICATE_BY_ID)

    @cache
    def go_preds(go: str) -> frozenset[str]:
        return frozenset(p for a in ancestors(go) for p in anchor_preds.get(a, ()))

    @cache
    def ec_preds(ec: str) -> frozenset[str]:
        return frozenset(get_predicates_for_ec(ec))

    proteins: list[Reviewed] = []
    accession, pfams, preds, genes = None, set(), set(), set()
    with gzip.open(path, "rt", encoding="latin-1") as handle:
        for line in handle:
            tag = line[:2]
            if tag == "AC" and accession is None:
                accession = line[5:].split(";")[0].strip()
            elif tag == "DR":
                if line.startswith("DR   Pfam; "):
                    pfams.add(line[11:18])
                elif line.startswith("DR   KEGG; "):
                    genes.add(line[11:].split(";")[0])
                elif line.startswith("DR   GO; "):
                    fields = line[5:].rstrip(".\n").split("; ")
                    if len(fields) >= 4 and fields[3].split(":")[0] in EXPERIMENTAL_GO:
                        preds |= go_preds(fields[1])
            elif tag == "DE" and "EC=" in line:
                preds |= ec_preds(line.split("EC=", 1)[1].split()[0].rstrip(";"))
            elif tag == "//":
                if pfams or genes:
                    proteins.append(Reviewed(accession, frozenset(pfams), frozenset(preds & vocabulary),
                                             frozenset(genes)))
                accession, pfams, preds, genes = None, set(), set(), set()
    return proteins


def read_gene_ko_links(path: Path) -> dict[str, frozenset[str]]:
    """KEGG gene -> KOs, from ``gene<TAB>KO[,KO]`` rows (see scripts/link_swissprot_kegg.py)."""
    links = {}
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            gene, kos = line.rstrip("\n").split("\t")
            links[gene] = frozenset(kos.split(","))
    return links


def protein_kos(protein: Reviewed, gene_ko: dict[str, frozenset[str]]) -> frozenset[str]:
    return frozenset(ko for gene in protein.kegg_genes for ko in gene_ko.get(gene, ()))


def ko_consensus(proteins: list[Reviewed], gene_ko, candidates: dict[str, set[str]]) -> dict[tuple[str, str], str]:
    """(KO, predicate) pairs on which the KO's reviewed proteins agree (coverage rule only).

    A KO is whole-protein orthology, so attribution to a co-occurring domain does
    not arise; per-protein predicates are EC and experimental GO only.
    """
    members: dict[str, list[Reviewed]] = defaultdict(list)
    for protein in proteins:
        for ko in protein_kos(protein, gene_ko):
            members[ko].append(protein)
    min_n, min_frac = SWISSPROT_COVERAGE
    result = {}
    for ko, preds in candidates.items():
        group = members.get(ko, ())
        if len(group) < min_n:
            continue
        for pred in preds - SWISSPROT_EXCLUDED:
            k = sum(pred in p.predicates for p in group)
            if k >= min_frac * len(group) and wilson_lower_bound(k, len(group)) >= SWISSPROT_MIN_LOWER_BOUND:
                result[(ko, pred)] = f"swissprot:{k}/{len(group)}"
    return result


def ortholog_support(members, pred) -> tuple[int, int]:
    """(KOs whose carriers mostly have ``pred``, distinct KOs) among ``members``.

    Swiss-Prot over-represents some orthologs (one reviewed entry per well-studied
    proteome), so a family's reviewed carriers can agree because one KO
    dominates; counting distinct KOs corrects for that redundancy.
    """
    by_ko: dict[str, list[bool]] = defaultdict(list)
    for preds, kos in members:
        for ko in kos:
            by_ko[ko].append(pred in preds)
    agreeing = sum(1 for hits in by_ko.values() if 2 * sum(hits) > len(hits))
    return agreeing, len(by_ko)


def family_consensus(proteins, candidates, supported) -> dict[tuple[str, str], str]:
    """(Pfam family, predicate) pairs supported by reviewed-protein consensus.

    ``proteins`` are (Pfam families, predicates, KOs) per reviewed protein.
    Coverage: the predicate holds for enough of the family's carriers and, when
    the carriers span two or more KOs, for enough of those KOs. Attribution:
    single-domain carriers agree, or (with too few of them) no co-occurring
    family that carries the predicate on its own evidence explains it.
    """
    by_family: dict[str, list[int]] = defaultdict(list)
    for i, (fams, _, _) in enumerate(proteins):
        for f in fams:
            by_family[f].append(i)
    min_n, min_frac = SWISSPROT_COVERAGE
    min_single, single_frac = SWISSPROT_SINGLE

    covered: dict[tuple[str, str], list[int]] = {}
    for fam, pred_set in candidates.items():
        members = by_family.get(fam, ())
        if len(members) < min_n:
            continue
        for pred in pred_set - SWISSPROT_EXCLUDED:
            hits = [i for i in members if pred in proteins[i][1]]
            if len(hits) < min_frac * len(members) \
                    or wilson_lower_bound(len(hits), len(members)) < SWISSPROT_MIN_LOWER_BOUND:
                continue
            k_kos, n_kos = ortholog_support([proteins[i][1:] for i in members], pred)
            if n_kos >= 2 and k_kos < min_frac * n_kos:
                continue
            covered[(fam, pred)] = (hits, k_kos, n_kos)

    result: dict[tuple[str, str], str] = {}
    unattributed = []
    for (fam, pred), (hits, k_kos, n_kos) in covered.items():
        coverage = f"swissprot:{len(hits)}/{len(by_family[fam])}"
        if n_kos >= 2:
            coverage += f" KOs {k_kos}/{n_kos}"
        single = [i for i in by_family[fam] if len(proteins[i][0]) == 1]
        if len(single) >= min_single:
            k = sum(pred in proteins[i][1] for i in single)
            if k >= single_frac * len(single):
                result[(fam, pred)] = f"{coverage} single-domain {k}/{len(single)}"
        else:
            unattributed.append((fam, pred, hits, coverage))

    def owns(fam, pred):
        return (fam, pred) in result or pred in supported.get(fam, ())

    for fam, pred, hits, coverage in unattributed:
        co: dict[str, int] = defaultdict(int)
        for i in hits:
            for g in proteins[i][0] - {fam}:
                co[g] += 1
        if not any(c >= SWISSPROT_CODOMAIN_FRACTION * len(hits) and owns(g, pred) for g, c in co.items()):
            result[(fam, pred)] = f"{coverage} co-domains excluded"
    return result


def pfam_view(proteins: list[Reviewed], gene_ko=None, ko_predicates=None) -> list[tuple[frozenset, frozenset, frozenset]]:
    """(Pfam families, predicates, KOs) per reviewed protein with Pfam, adding KEGG-evidenced KO predicates."""
    view = []
    for p in proteins:
        if not p.pfams:
            continue
        preds = set(p.predicates)
        kos = protein_kos(p, gene_ko) if gene_ko else frozenset()
        for ko in kos:
            preds |= (ko_predicates or {}).get(ko, set())
        view.append((p.pfams, frozenset(preds), kos))
    return view
