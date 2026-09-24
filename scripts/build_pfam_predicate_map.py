#!/usr/bin/env python3
"""Build the shipped Pfam -> predicate map from proposals and Pfam's own information.

Proposals come from ``pfam_predicate_proposals.tsv`` (curated family rows) and
``pfam_predicate_proposal_patterns.tsv`` (regexes applied to every family's
name/description). Each proposed (family, predicate) pair ships only when
:func:`sharur.predicates.mappings.pfam_evidence.judge` finds support in the
family's InterPro GO annotation, its Pfam name/description, or an ENZYME name
in its description. With ``--swissprot``, a pair also ships when reviewed
Swiss-Prot proteins carrying the family agree on it (curated EC numbers and
experimentally evidenced GO terms; see :func:`swissprot_consensus`). Outputs, all under ``sharur/predicates/mappings/data``:

- ``pfam_predicates.tsv``: accession, name, predicates, per-predicate evidence.
- ``pfam_evidence_snapshot.tsv``: the Pfam/GO/ENZYME facts each shipped pair
  relies on, so the test suite re-verifies the map without network access.

Usage:
    python scripts/build_pfam_predicate_map.py \\
        --pfam-hmm ~/.config/Astra/PFAM/Pfam-A.hmm --pfam-clans Pfam-A.clans.tsv.gz \\
        --pfam2go pfam2go --go-obo go-basic.obo --enzyme-dat enzyme.dat \\
        --swissprot uniprot_sprot.dat.gz --swissprot-release 2026_03 \\
        --report dropped_pairs.tsv
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import os
import re
import subprocess
import sys
from collections import defaultdict
from functools import cache
from pathlib import Path


sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sharur.predicates.mappings.kegg_map import get_predicates_for_ec
from sharur.predicates.mappings.pfam_evidence import (
    EVIDENCE,
    SWISSPROT_CODOMAIN_FRACTION,
    SWISSPROT_COVERAGE,
    SWISSPROT_EXCLUDED,
    SWISSPROT_MIN_LOWER_BOUND,
    SWISSPROT_SINGLE,
    enzyme_phrases,
    judge,
    resolve,
    wilson_lower_bound,
)
from sharur.predicates.vocabulary import PREDICATE_BY_ID


DATA = Path(__file__).resolve().parents[1] / "sharur/predicates/mappings/data"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_hmm_headers(path: Path) -> dict[str, tuple[str, str]]:
    out = subprocess.run(["grep", "-E", "^(NAME|ACC|DESC) ", str(path)],
                         capture_output=True, text=True, check=True).stdout
    families, cur = {}, {}
    for line in out.splitlines():
        tag, _, val = line.partition(" ")
        val = val.strip()
        if tag == "NAME":
            cur = {"name": val}
        elif tag == "ACC":
            cur["acc"] = val.split(".")[0]
        elif tag == "DESC":
            families[cur["acc"]] = (cur["name"], val)
    return families


def read_clans(path: Path) -> dict[str, tuple[str, str]]:
    with gzip.open(path, "rt") as handle:
        return {row[0]: (row[3], row[4]) for row in csv.reader(handle, delimiter="\t")}


def read_go(pfam2go: Path, obo: Path) -> dict[str, set[str]]:
    parents: dict[str, list[str]] = defaultdict(list)
    term = None
    for line in open(obo):
        line = line.rstrip("\n")
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

    closure: dict[str, set[str]] = defaultdict(set)
    for line in open(pfam2go):
        m = re.match(r"Pfam:(PF\d+) \S+ > GO:.* ; (GO:\d+)", line)
        if m:
            closure[m.group(1)] |= ancestors(m.group(2))
    read_go.ancestors = ancestors
    return closure


def read_enzyme_names(path: Path) -> dict[str, list[str]]:
    names: dict[str, set[str]] = defaultdict(set)
    cur = None
    for line in open(path, encoding="latin-1"):
        tag, val = line[:2], line[5:].rstrip("\n").strip()
        if tag == "ID":
            cur = {"ec": val, "de": "", "an": []}
        elif tag == "DE" and cur is not None:
            cur["de"] = (cur["de"] + " " + val).strip()
        elif tag == "AN" and cur is not None:
            if cur["an"] and not cur["an"][-1].endswith("."):
                cur["an"][-1] += " " + val
            else:
                cur["an"].append(val)
        elif tag == "//" and cur is not None:
            if not cur["de"].startswith(("Transferred entry", "Deleted entry")):
                for n in [cur["de"], *cur["an"]]:
                    n = n.rstrip(".").strip().lower()
                    if n:
                        names[n].add(cur["ec"])
                        stripped = re.sub(r"\s+", " ", re.sub(r"\s*\([^)]*\)\s*-?", " ", n)).strip(" -")
                        if stripped != n and len(stripped) >= 7:
                            names[stripped].add(cur["ec"])
            cur = None
    return {k: sorted(v) for k, v in names.items()}


EXPERIMENTAL_GO = {"EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP"}


def read_swissprot(path: Path, ancestors) -> list[tuple[frozenset, frozenset]]:
    """Reviewed proteins as (Pfam accessions, predicates from curated EC and experimental GO)."""
    anchor_preds: dict[str, set[str]] = defaultdict(set)
    for pred, spec in EVIDENCE.items():
        for go in spec["go"]:
            anchor_preds[go].add(pred)

    @cache
    def go_preds(go: str) -> frozenset[str]:
        return frozenset(p for a in ancestors(go) for p in anchor_preds.get(a, ()))

    @cache
    def ec_preds(ec: str) -> frozenset[str]:
        return frozenset(get_predicates_for_ec(ec))

    proteins, pfams, preds = [], set(), set()
    with gzip.open(path, "rt", encoding="latin-1") as handle:
        for line in handle:
            tag = line[:2]
            if tag == "DR":
                if line.startswith("DR   Pfam; "):
                    pfams.add(line[11:18])
                elif line.startswith("DR   GO; "):
                    fields = line[5:].rstrip(".\n").split("; ")
                    if len(fields) >= 4 and fields[3].split(":")[0] in EXPERIMENTAL_GO:
                        preds |= go_preds(fields[1])
            elif tag == "DE" and "EC=" in line:
                preds |= ec_preds(line.split("EC=", 1)[1].split()[0].rstrip(";"))
            elif tag == "//":
                if pfams:
                    proteins.append((frozenset(pfams), frozenset(preds & set(PREDICATE_BY_ID))))
                pfams, preds = set(), set()
    return proteins


def swissprot_consensus(proteins, candidates, supported) -> dict[tuple[str, str], str]:
    """Evidence strings for (family, predicate) pairs supported by reviewed-protein consensus.

    Coverage and attribution thresholds are documented beside SWISSPROT_COVERAGE.
    """
    by_family: dict[str, list[int]] = defaultdict(list)
    for i, (fams, _) in enumerate(proteins):
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
            if len(hits) >= min_frac * len(members) \
                    and wilson_lower_bound(len(hits), len(members)) >= SWISSPROT_MIN_LOWER_BOUND:
                covered[(fam, pred)] = hits

    result: dict[tuple[str, str], str] = {}
    unattributed = []
    for (fam, pred), hits in covered.items():
        coverage = f"swissprot:{len(hits)}/{len(by_family[fam])}"
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
        co = defaultdict(int)
        for i in hits:
            for g in proteins[i][0] - {fam}:
                co[g] += 1
        if not any(c >= SWISSPROT_CODOMAIN_FRACTION * len(hits) and owns(g, pred) for g, c in co.items()):
            result[(fam, pred)] = f"{coverage} co-domains excluded"
    return result


def read_proposals(path: Path, by_name: dict[str, str]) -> dict[str, set[str]]:
    proposals: dict[str, set[str]] = defaultdict(set)
    for line in open(path):
        if not line.strip() or line.startswith("#"):
            continue
        key, preds = line.rstrip("\n").split("\t")[:2]
        acc = key if key.startswith("PF") else by_name.get(key)
        if acc is None:
            raise SystemExit(f"{path}: unknown Pfam family {key!r}")
        proposals[acc].update(p.strip() for p in preds.split(",") if p.strip())
    return proposals


def read_patterns(path: Path) -> list[tuple[re.Pattern, list[str]]]:
    patterns = []
    for line in open(path):
        if not line.strip() or line.startswith("#"):
            continue
        pattern, preds = line.rstrip("\n").split("\t")[:2]
        predicates = [p.strip() for p in preds.split(",") if p.strip()]
        if predicates:
            patterns.append((re.compile(pattern, re.IGNORECASE), predicates))
    return patterns


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--pfam-hmm", type=Path, default=Path.home() / ".config/Astra/PFAM/Pfam-A.hmm")
    parser.add_argument("--pfam-clans", type=Path, required=True, help="Pfam-A.clans.tsv.gz (current release)")
    parser.add_argument("--pfam2go", type=Path, required=True)
    parser.add_argument("--go-obo", type=Path, required=True)
    parser.add_argument("--enzyme-dat", type=Path, required=True)
    parser.add_argument("--swissprot", type=Path, help="uniprot_sprot.dat.gz (reviewed-protein consensus tier)")
    parser.add_argument("--swissprot-release", default="", help="Release label recorded in the output header")
    parser.add_argument("--report", type=Path, help="TSV of proposed pairs that lack evidence")
    args = parser.parse_args()

    installed = read_hmm_headers(args.pfam_hmm)
    current = read_clans(args.pfam_clans)
    families = {**current, **installed}  # the installed release's wording wins
    by_name = {name: acc for acc, (name, _) in families.items()}
    closure = read_go(args.pfam2go, args.go_obo)
    enzymes = read_enzyme_names(args.enzyme_dat)

    proposals = read_proposals(DATA / "pfam_predicate_proposals.tsv", by_name)
    for pattern, preds in read_patterns(DATA / "pfam_predicate_proposal_patterns.tsv"):
        for acc, (name, desc) in families.items():
            if pattern.search(f"{name} {desc}".lower()):
                proposals[acc].update(preds)

    shipped: dict[str, dict[str, str]] = {}
    pending: dict[str, dict[str, str]] = defaultdict(dict)  # acc -> predicate -> original proposal
    for acc in sorted(proposals):
        name, desc = families[acc]
        for proposed in sorted(proposals[acc]):
            predicate = resolve(proposed, PREDICATE_BY_ID)
            evidence = judge(predicate, name, desc, closure.get(acc, ()), enzymes) if predicate else None
            if evidence:
                shipped.setdefault(acc, {})[predicate] = evidence
            else:
                pending[acc][predicate or f"(no vocabulary equivalent: {proposed})"] = proposed

    swissprot_note = ""
    if args.swissprot:
        proteins = read_swissprot(args.swissprot, read_go.ancestors)
        candidates: dict[str, set[str]] = defaultdict(set)
        for acc, preds in pending.items():
            candidates[acc] |= {p for p in preds if p in PREDICATE_BY_ID}
        for fams, preds in proteins:  # consensus may also add predicates nobody proposed
            if len(fams) == 1:
                (fam,) = fams
                if fam in families:
                    candidates[fam] |= preds
        supported = {acc: set(v) for acc, v in shipped.items()}
        added = 0
        for (acc, pred), evidence in sorted(swissprot_consensus(proteins, candidates, supported).items()):
            if acc in families and pred not in shipped.get(acc, {}):
                shipped.setdefault(acc, {})[pred] = evidence
                pending.get(acc, {}).pop(pred, None)
                added += 1
        swissprot_note = f"# Swiss-Prot: {args.swissprot.name} {args.swissprot_release} sha256 {sha256(args.swissprot)[:16]} ({len(proteins):,} reviewed proteins with Pfam)\n"
        print(f"Swiss-Prot consensus added {added:,} pairs")

    dropped = [(acc, *families[acc], proposed, pred) for acc, preds in sorted(pending.items())
               for pred, proposed in sorted(preds.items())]

    sources = (
        f"# Pfam: {args.pfam_hmm.name} sha256 {sha256(args.pfam_hmm)[:16]} ({len(installed)} families); "
        f"{args.pfam_clans.name} ({len(current)} families)\n"
        f"# GO: {args.pfam2go.name} sha256 {sha256(args.pfam2go)[:16]}; {args.go_obo.name} sha256 {sha256(args.go_obo)[:16]}\n"
        f"# ENZYME: {args.enzyme_dat.name} sha256 {sha256(args.enzyme_dat)[:16]}\n"
        + swissprot_note
    )
    anchors = {g for spec in EVIDENCE.values() for g in spec["go"]}
    with open(DATA / "pfam_predicates.tsv", "w") as out, open(DATA / "pfam_evidence_snapshot.tsv", "w") as snap:
        out.write("# Generated by scripts/build_pfam_predicate_map.py; edit the proposals, not this file.\n" + sources)
        out.write("# accession\tname\tpredicates\tevidence\n")
        snap.write("# Pfam/GO/ENZYME facts behind pfam_predicates.tsv (generated).\n" + sources)
        snap.write("# accession\tname\tdescription\tgo_anchors\tenzymes\n")
        # Installed-release families first: they own a profile name that a later
        # release reassigned, because that is the name Astra reports.
        for acc in sorted(shipped, key=lambda a: (a not in installed, a)):
            name, desc = families[acc]
            preds = shipped[acc]
            out.write(f"{acc}\t{name}\t{','.join(sorted(preds))}\t"
                      + ";".join(f"{p}={preds[p]}" for p in sorted(preds)) + "\n")
            go = sorted(closure.get(acc, set()) & anchors)
            enz = ";".join(f"{ph}={'|'.join(ecs)}" for ph, ecs in sorted(enzyme_phrases(desc, enzymes).items()))
            snap.write(f"{acc}\t{name}\t{desc}\t{','.join(go)}\t{enz}\n")

    if args.report:
        with open(args.report, "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["accession", "name", "description", "proposed", "resolved"])
            writer.writerows(dropped)

    pairs = sum(len(v) for v in shipped.values())
    print(f"proposed families: {len(proposals):,}; shipped: {len(shipped):,} families, {pairs:,} pairs; "
          f"dropped pairs without evidence: {len(dropped):,}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
