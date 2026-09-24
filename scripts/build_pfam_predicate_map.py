#!/usr/bin/env python3
"""Build the shipped Pfam -> predicate map from proposals and Pfam's own information.

Proposals come from ``pfam_predicate_proposals.tsv`` (curated family rows) and
``pfam_predicate_proposal_patterns.tsv`` (regexes applied to every family's
name/description). Each proposed (family, predicate) pair ships only when
:func:`sharur.predicates.mappings.pfam_evidence.judge` finds support in the
family's InterPro GO annotation, its Pfam name/description, or an ENZYME name
in its description. Outputs, all under ``sharur/predicates/mappings/data``:

- ``pfam_predicates.tsv``: accession, name, predicates, per-predicate evidence.
- ``pfam_evidence_snapshot.tsv``: the Pfam/GO/ENZYME facts each shipped pair
  relies on, so the test suite re-verifies the map without network access.

Usage:
    python scripts/build_pfam_predicate_map.py \\
        --pfam-hmm ~/.config/Astra/PFAM/Pfam-A.hmm --pfam-clans Pfam-A.clans.tsv.gz \\
        --pfam2go pfam2go --go-obo go-basic.obo --enzyme-dat enzyme.dat \\
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

from sharur.predicates.mappings.pfam_evidence import (  # noqa: E402
    EVIDENCE,
    enzyme_phrases,
    judge,
    resolve,
)
from sharur.predicates.vocabulary import PREDICATE_BY_ID  # noqa: E402


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
    dropped = []
    for acc in sorted(proposals):
        name, desc = families[acc]
        for proposed in sorted(proposals[acc]):
            predicate = resolve(proposed, PREDICATE_BY_ID)
            evidence = judge(predicate, name, desc, closure.get(acc, ()), enzymes) if predicate else None
            if evidence:
                shipped.setdefault(acc, {})[predicate] = evidence
            else:
                dropped.append((acc, name, desc, proposed, predicate or "(no vocabulary equivalent)"))

    sources = (
        f"# Pfam: {args.pfam_hmm.name} sha256 {sha256(args.pfam_hmm)[:16]} ({len(installed)} families); "
        f"{args.pfam_clans.name} ({len(current)} families)\n"
        f"# GO: {args.pfam2go.name} sha256 {sha256(args.pfam2go)[:16]}; {args.go_obo.name} sha256 {sha256(args.go_obo)[:16]}\n"
        f"# ENZYME: {args.enzyme_dat.name} sha256 {sha256(args.enzyme_dat)[:16]}\n"
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
