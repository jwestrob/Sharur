#!/usr/bin/env python3
"""Build the shipped KO -> predicate map from KEGG's own information.

Every KO receives the predicates its KEGG-stated EC numbers, BRITE placements,
module memberships and (for KOs KEGG names as hydrogenases) HydDB reference
labels support (rules in ``kegg_brite_predicates.tsv`` and
``kegg_module_predicates.tsv``). Curated proposals
(``kegg_predicate_proposals.tsv``) and definition patterns
(``kegg_predicate_proposal_patterns.tsv``) add a pair only when one of those
sources or the KO's own symbols/name support it. With ``--swissprot``, a pair
also ships when the KO's reviewed Swiss-Prot proteins agree on it through
curated EC numbers and experimental GO
(:mod:`sharur.predicates.mappings.kegg_evidence`). Outputs, under
``sharur/predicates/mappings/data``:

- ``kegg_predicates.tsv``: KO, predicates, per-predicate evidence.
- ``kegg_evidence_snapshot.tsv``: KEGG's symbols, name, BRITE placements and
  module memberships for each shipped KO, so the test suite re-verifies the map
  without network access.

Inputs are KEGG REST snapshots: ``list/ko``, ``get/br:<hierarchy>/json`` for
each hierarchy named in the BRITE rules, and ``get/<module>`` entries for every
module (concatenated flat file).

Usage:
    python scripts/build_kegg_predicate_map.py --ko-list ko_list.tsv \\
        --kofam-ko-list data/reference/ko_list --brite-dir brite/ --modules modules.txt --kegg-release "KEGG REST 2026-09-24" \\
        --report dropped_kegg_pairs.tsv
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
from collections import defaultdict
from pathlib import Path


sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sharur.predicates.mappings.kegg_evidence import (
    PATH_SEP,
    brite_evidence,
    hyddb_evidence,
    ko_definition,
    ko_text_evidence,
    load_brite_rules,
    load_hyddb_snapshot,
    load_module_rules,
    module_evidence,
    parse_module_definition,
)
from sharur.predicates.mappings.kegg_map import get_predicates_for_ec, parse_ec_numbers
from sharur.predicates.mappings.pfam_evidence import resolve
from sharur.predicates.mappings.swissprot_evidence import (
    go_ancestors,
    ko_consensus,
    protein_kos,
    read_gene_ko_links,
    read_swissprot,
)
from sharur.predicates.vocabulary import PREDICATE_BY_ID


DATA = Path(__file__).resolve().parents[1] / "sharur/predicates/mappings/data"
KO = re.compile(r"K\d{5}")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_ko_list(path: Path) -> dict[str, tuple[str, str]]:
    """KEGG REST list/ko -> {KO: (symbols, name)}."""
    kos = {}
    with open(path) as handle:
        for line in handle:
            ko, _, definition = line.rstrip("\n").partition("\t")
            ko = ko.removeprefix("ko:")
            symbols, sep, name = definition.partition("; ")
            kos[ko] = (symbols, name) if sep else ("", definition)
    return kos


def read_kofam_ko_list(path: Path) -> dict[str, tuple[str, str]]:
    """KOfam ko_list -> {KO: ("", definition)} (KOfam carries no symbols)."""
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        ko_idx, def_idx = header.index("knum"), header.index("definition")
        return {row[ko_idx]: ("", row[def_idx])
                for row in (line.rstrip("\n").split("\t") for line in handle)}


def read_brite(brite_dir: Path, hierarchies) -> dict[str, list[tuple[str, tuple[str, ...]]]]:
    """KO -> [(hierarchy, node labels from the top level to the KO's parent)]."""
    placements: dict[str, list] = defaultdict(list)

    def walk(node, hierarchy, labels):
        for child in node.get("children", ()):
            if "children" in child:
                label = child["name"]
                if PATH_SEP in label or "|" in label:
                    raise SystemExit(f"{hierarchy}: node label {label!r} contains a reserved separator")
                walk(child, hierarchy, (*labels, label))
            elif KO.match(child["name"]):
                entry = (hierarchy, labels)
                ko = child["name"][:6]
                if entry not in placements[ko]:
                    placements[ko].append(entry)

    for hierarchy in sorted(hierarchies):
        with open(brite_dir / f"{hierarchy}.json") as handle:
            walk(json.load(handle), hierarchy, ())
    return placements


def read_modules(path: Path) -> tuple[dict[str, str], dict[str, list[tuple[str, bool]]]]:
    """Module flat file -> ({module: name}, {KO: [(module, is_complex_component)]})."""
    names: dict[str, str] = {}
    memberships: dict[str, list] = defaultdict(list)
    entry, field, definition = None, None, []

    def flush():
        if entry:
            for ko, in_complex in parse_module_definition(" ".join(definition)).items():
                memberships[ko].append((entry, in_complex))

    with open(path) as handle:
        for line in handle:
            tag, value = line[:12].strip(), line[12:].rstrip("\n")
            if tag == "ENTRY":
                flush()
                entry, definition = value.split()[0], []
            elif tag == "NAME":
                names[entry] = value.strip()
            elif tag == "///":
                flush()
                entry = None
            if tag:
                field = tag
            if field == "DEFINITION" and tag in ("DEFINITION", ""):
                definition.append(value.strip())
    flush()
    return names, memberships


def read_proposals() -> dict[str, set[str]]:
    proposals: dict[str, set[str]] = defaultdict(set)
    for line in (DATA / "kegg_predicate_proposals.tsv").read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            ko, preds = line.split("\t")[:2]
            proposals[ko].update(p for p in preds.split(",") if p)
    return proposals


def read_patterns() -> list[tuple[re.Pattern, list[str]]]:
    patterns = []
    for line in (DATA / "kegg_predicate_proposal_patterns.tsv").read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            pattern, preds = line.split("\t")[:2]
            patterns.append((re.compile(pattern, re.IGNORECASE), [p for p in preds.split(",") if p]))
    return patterns


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--ko-list", type=Path, required=True, help="KEGG REST list/ko snapshot")
    parser.add_argument("--kofam-ko-list", type=Path, default=Path("data/reference/ko_list"),
                        help="KOfam ko_list; supplies definitions for KOfam KOs retired from KEGG REST")
    parser.add_argument("--brite-dir", type=Path, required=True, help="Directory of <hierarchy>.json files")
    parser.add_argument("--modules", type=Path, required=True, help="Concatenated KEGG module entries")
    parser.add_argument("--kegg-release", required=True, help="Release/retrieval label recorded in the header")
    parser.add_argument("--swissprot", type=Path, help="uniprot_sprot.dat.gz (reviewed-protein consensus tier)")
    parser.add_argument("--swissprot-kegg", type=Path, help="KEGG gene -> KO links (scripts/link_swissprot_kegg.py)")
    parser.add_argument("--go-obo", type=Path, help="go-basic.obo, for experimental GO ancestor closure")
    parser.add_argument("--swissprot-release", default="", help="Release label recorded in the header")
    parser.add_argument("--report", type=Path, help="TSV of proposed pairs that lack evidence")
    args = parser.parse_args()

    kos = read_ko_list(args.ko_list)
    retired = {ko: v for ko, v in read_kofam_ko_list(args.kofam_ko_list).items() if ko not in kos}
    kos.update(retired)
    brite_rules = load_brite_rules()
    module_rules = load_module_rules()
    placements = read_brite(args.brite_dir, {r.hierarchy for r in brite_rules})
    module_names, memberships = read_modules(args.modules)
    hyddb = load_hyddb_snapshot()
    for module, (name, _) in module_rules.items():
        if module_names.get(module) != name:
            raise SystemExit(f"{module}: KEGG names it {module_names.get(module)!r}, rules say {name!r}")
    unknown = {p for r in brite_rules for p in r.predicates} | {p for _, ps in module_rules.values() for p in ps}
    unknown -= set(PREDICATE_BY_ID)
    if unknown:
        raise SystemExit(f"rules name predicates outside the vocabulary: {sorted(unknown)}")

    proposals = read_proposals()
    missing = sorted(set(proposals) - set(kos))
    if missing:
        raise SystemExit(f"proposed KOs absent from the KO list: {missing}")
    for pattern, preds in read_patterns():
        for ko, (symbols, name) in kos.items():
            if pattern.search(ko_definition(symbols, name)):
                proposals[ko].update(preds)

    shipped: dict[str, dict[str, str]] = {}
    dropped = []
    for ko in sorted(kos):
        symbols, name = kos[ko]
        definition = ko_definition(symbols, name)
        evidence: dict[str, str] = {}
        for ec in parse_ec_numbers(name):
            for pred in get_predicates_for_ec(ec):
                evidence.setdefault(pred, f"ec:{ec}")
        for pred, ev in hyddb_evidence(hyddb.get(ko, {})).items():
            evidence.setdefault(pred, ev)
        for pred, ev in {**module_evidence(module_rules, memberships.get(ko, ())),
                         **brite_evidence(brite_rules, placements.get(ko, ()), definition)}.items():
            evidence.setdefault(pred, ev)
        for proposed in sorted(proposals.get(ko, ())):
            pred = resolve(proposed, PREDICATE_BY_ID)
            if pred in evidence:
                continue
            match = ko_text_evidence(pred, symbols, name) if pred else None
            if match:
                evidence[pred] = f"text:{match}"
            else:
                dropped.append((ko, symbols, name, proposed, pred or "(no vocabulary equivalent)"))
        evidence = {p: e for p, e in evidence.items() if p in PREDICATE_BY_ID}
        if evidence:
            shipped[ko] = evidence

    swissprot_note = ""
    if args.swissprot:
        if not (args.swissprot_kegg and args.go_obo):
            raise SystemExit("--swissprot needs --swissprot-kegg and --go-obo")
        proteins = read_swissprot(args.swissprot, go_ancestors(args.go_obo))
        gene_ko = read_gene_ko_links(args.swissprot_kegg)
        candidates: dict[str, set[str]] = defaultdict(set)
        for ko, _, _, _, pred in dropped:
            if pred in PREDICATE_BY_ID:
                candidates[ko].add(pred)
        for protein in proteins:  # consensus may also add predicates nobody proposed
            for ko in protein_kos(protein, gene_ko):
                if ko in kos:
                    candidates[ko] |= protein.predicates
        added = 0
        for (ko, pred), ev in sorted(ko_consensus(proteins, gene_ko, candidates).items()):
            if pred not in shipped.get(ko, {}):
                shipped.setdefault(ko, {})[pred] = ev
                added += 1
        dropped = [row for row in dropped if row[4] not in shipped.get(row[0], {})]
        linked = sum(1 for p in proteins if protein_kos(p, gene_ko))
        swissprot_note = (f"# Swiss-Prot: {args.swissprot.name} {args.swissprot_release} sha256 {sha256(args.swissprot)[:16]} "
                          f"({linked:,} reviewed proteins linked to KOs); KEGG gene links {args.swissprot_kegg.name} "
                          f"sha256 {sha256(args.swissprot_kegg)[:16]}; {args.go_obo.name} sha256 {sha256(args.go_obo)[:16]}\n")
        print(f"Swiss-Prot consensus added {added:,} pairs")
    shipped = dict(sorted(shipped.items()))

    header = (
        "# Generated by scripts/build_kegg_predicate_map.py; edit the proposals or rules, not this file.\n"
        f"# KEGG: {args.kegg_release}; {args.ko_list.name} sha256 {sha256(args.ko_list)[:16]} ({len(kos)} KOs); "
        f"{args.modules.name} sha256 {sha256(args.modules)[:16]} ({len(module_names)} modules); "
        f"{args.kofam_ko_list.name} sha256 {sha256(args.kofam_ko_list)[:16]} ({len(retired)} KOfam KOs absent from KEGG REST)\n"
        + swissprot_note
    )
    with open(DATA / "kegg_predicates.tsv", "w") as out, open(DATA / "kegg_evidence_snapshot.tsv", "w") as snap:
        out.write(header + "# ko\tpredicates\tevidence\n")
        snap.write(header.replace("Generated by", "KEGG facts behind kegg_predicates.tsv, generated by", 1)
                   + "# ko\tsymbols\tname\tbrite\tmodules\n")
        for ko, evidence in shipped.items():
            symbols, name = kos[ko]
            out.write(f"{ko}\t{','.join(sorted(evidence))}\t{';'.join(f'{p}={e}' for p, e in sorted(evidence.items()))}\n")
            brite = "|".join(f"{h}:{PATH_SEP.join(labels)}" for h, labels in placements.get(ko, ()))
            modules = ",".join(m + ("+" if c else "") for m, c in memberships.get(ko, ()))
            snap.write(f"{ko}\t{symbols}\t{name}\t{brite}\t{modules}\n")

    if args.report:
        with open(args.report, "w") as rep:
            rep.write("ko\tsymbols\tname\tproposed\tresolved\n")
            for row in dropped:
                rep.write("\t".join(row) + "\n")
    n_pairs = sum(map(len, shipped.values()))
    print(f"KOs: {len(kos):,}; shipped: {len(shipped):,} KOs, {n_pairs:,} pairs; "
          f"dropped proposed pairs without evidence: {len(dropped):,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
