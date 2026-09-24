"""Fetch KEGG data and build the KO -> predicate map locally.

Sharur ships its own rules (``data/kegg_*``: BRITE and module rules, KO
proposals keyed by KO ID, definition patterns, HydDB x KOfam counts, and
reviewed-protein consensus counts). KEGG text (KO names, BRITE placements,
module definitions) is fetched by the user at setup and stays on the user's
machine, subject to KEGG's terms: KEGG REST is for academic use; other users
need a KEGG license and can build from their licensed copy with ``--inputs``.

``sharur setup-kegg`` writes, into the KEGG directory
(:func:`sharur.predicates.mappings.kegg_map.kegg_dir`):

- ``inputs/``: the fetched KEGG files (list/ko, BRITE JSON, module entries);
- ``kegg_predicates.tsv``: KO, predicates, per-predicate evidence;
- ``kegg_evidence_snapshot.tsv``: KEGG symbols, name, BRITE placements and
  module memberships per mapped KO (for the integrity tests);
- ``provenance.json``: KEGG release, retrieval time, input and rule checksums.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import re
import time
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from sharur.predicates.mappings.kegg_evidence import (
    DATA,
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
from sharur.predicates.vocabulary import PREDICATE_BY_ID, component_level


REST = "https://rest.kegg.jp"
KOFAM_KO_LIST_URL = "https://www.genome.jp/ftp/db/kofam/ko_list.gz"
REQUEST_INTERVAL = 0.4  # KEGG REST asks for at most 3 calls per second
KO = re.compile(r"K\d{5}")
RULE_FILES = (
    "kegg_brite_predicates.tsv",
    "kegg_module_predicates.tsv",
    "kegg_predicate_proposals.tsv",
    "kegg_predicate_proposal_patterns.tsv",
    "kegg_hyddb_snapshot.tsv",
    "kegg_swissprot_consensus.tsv",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# --------------------------------------------------------------------------- #
# Fetch
# --------------------------------------------------------------------------- #


def _get(url: str, retries: int = 3) -> bytes:
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=300) as response:
                data = response.read()
            time.sleep(REQUEST_INTERVAL)
            return data
        except OSError:
            if attempt == retries - 1:
                raise
            time.sleep(5 * (attempt + 1))
    raise AssertionError("unreachable")


def fetch_inputs(inputs: Path, kofam_ko_list: Path | None = None, log=print) -> None:
    """Download the KEGG REST files the build needs (about 100 requests)."""
    inputs.mkdir(parents=True, exist_ok=True)
    (inputs / "brite").mkdir(exist_ok=True)
    log("Fetching KEGG release information and KO list")
    (inputs / "kegg_info.txt").write_bytes(_get(f"{REST}/info/kegg"))
    (inputs / "ko_list.tsv").write_bytes(_get(f"{REST}/list/ko"))
    hierarchies = sorted({r.hierarchy for r in load_brite_rules()})
    for i, hierarchy in enumerate(hierarchies, 1):
        log(f"Fetching BRITE {hierarchy} ({i}/{len(hierarchies)})")
        (inputs / "brite" / f"{hierarchy}.json").write_bytes(_get(f"{REST}/get/br:{hierarchy}/json"))
    modules = [line.split("\t")[0].removeprefix("md:")
               for line in _get(f"{REST}/list/module").decode().splitlines() if line]
    entries = []
    for start in range(0, len(modules), 10):
        log(f"Fetching modules {start + 1}-{min(start + 10, len(modules))} of {len(modules)}")
        entries.append(_get(f"{REST}/get/{'+'.join(modules[start:start + 10])}").decode())
    (inputs / "modules.txt").write_text("".join(entries))
    if kofam_ko_list is None or not kofam_ko_list.exists():
        log("Fetching KOfam ko_list")
        (inputs / "kofam_ko_list").write_bytes(gzip.decompress(_get(KOFAM_KO_LIST_URL)))


# --------------------------------------------------------------------------- #
# Readers
# --------------------------------------------------------------------------- #


def read_ko_list(path: Path) -> dict[str, tuple[str, str]]:
    """KEGG REST list/ko -> {KO: (symbols, name)}."""
    kos = {}
    with open(path) as handle:
        for line in handle:
            ko, _, definition = line.rstrip("\n").partition("\t")
            symbols, sep, name = definition.partition("; ")
            kos[ko.removeprefix("ko:")] = (symbols, name) if sep else ("", definition)
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
                    raise ValueError(f"{hierarchy}: node label {label!r} contains a reserved separator")
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


def read_modules(path: Path) -> tuple[set[str], dict[str, list[tuple[str, bool]]]]:
    """Module flat file -> (module IDs, {KO: [(module, is_complex_component)]})."""
    modules: set[str] = set()
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
                modules.add(entry)
            elif tag == "///":
                flush()
                entry = None
            if tag:
                field = tag
            if field == "DEFINITION" and tag in ("DEFINITION", ""):
                definition.append(value.strip())
    flush()
    return modules, memberships


def kegg_release(info: Path) -> str:
    """``ko 2026/09/24; brite 2026/09/24; module 2026/09/11`` from KEGG REST info/kegg."""
    if not info.exists():
        return "unknown (no info/kegg in inputs)"
    dates = {}
    for line in info.read_text().splitlines()[1:]:
        fields = line.split()
        if len(fields) == 3 and fields[0] in ("ko", "brite", "module"):
            dates[fields[0]] = fields[2]
    return "; ".join(f"{db} {dates[db]}" for db in ("ko", "brite", "module") if db in dates) or "unknown"


def _rows(name: str) -> list[list[str]]:
    return [line.split("\t") for line in (DATA / name).read_text().splitlines()
            if line.strip() and not line.startswith("#")]


def read_proposals() -> dict[str, set[str]]:
    proposals: dict[str, set[str]] = defaultdict(set)
    for ko, preds, *_ in _rows("kegg_predicate_proposals.tsv"):
        proposals[ko].update(p for p in preds.split(",") if p)
    return proposals


def read_patterns() -> list[tuple[re.Pattern, list[str]]]:
    return [(re.compile(pattern, re.IGNORECASE), [p for p in preds.split(",") if p])
            for pattern, preds, *_ in _rows("kegg_predicate_proposal_patterns.tsv")]


def load_swissprot_consensus() -> dict[str, dict[str, str]]:
    """KO -> {predicate: 'swissprot:k/n'} from the shipped reviewed-protein consensus table."""
    table: dict[str, dict[str, str]] = defaultdict(dict)
    for ko, pred, k, n in _rows("kegg_swissprot_consensus.tsv"):
        table[ko][pred] = f"swissprot:{k}/{n}"
    return table


# --------------------------------------------------------------------------- #
# Build
# --------------------------------------------------------------------------- #


def build(inputs: Path, out_dir: Path, kofam_ko_list: Path | None = None,
          report: Path | None = None, log=print) -> dict:
    """Build ``kegg_predicates.tsv`` and its evidence snapshot from fetched inputs."""
    kofam_ko_list = kofam_ko_list if kofam_ko_list and kofam_ko_list.exists() else inputs / "kofam_ko_list"
    kos = read_ko_list(inputs / "ko_list.tsv")
    retired = {ko: v for ko, v in read_kofam_ko_list(kofam_ko_list).items() if ko not in kos}
    kos.update(retired)
    brite_rules = load_brite_rules()
    module_rules = load_module_rules()
    placements = read_brite(inputs / "brite", {r.hierarchy for r in brite_rules})
    modules, memberships = read_modules(inputs / "modules.txt")
    missing_modules = sorted(set(module_rules) - modules)
    if missing_modules:
        raise ValueError(f"module rules name modules absent from this KEGG release: {missing_modules}")
    hyddb = load_hyddb_snapshot()
    consensus = load_swissprot_consensus()

    proposals = read_proposals()
    for pattern, preds in read_patterns():
        for ko, (symbols, name) in kos.items():
            if pattern.search(ko_definition(symbols, name)):
                proposals[ko].update(preds)

    mapped: dict[str, dict[str, str]] = {}
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
            elif pred not in consensus.get(ko, {}):
                dropped.append((ko, symbols, name, proposed, pred or "(no vocabulary equivalent)"))
        for pred, ev in consensus.get(ko, {}).items():
            evidence.setdefault(pred, ev)
        component: dict[str, str] = {}
        for pred, ev in evidence.items():  # maps carry component-level predicates only
            if pred in PREDICATE_BY_ID:
                component.setdefault(component_level(pred), ev)
        if component:
            mapped[ko] = component

    provenance = {
        "kegg_release": kegg_release(inputs / "kegg_info.txt"),
        "built": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "inputs": {
            "ko_list.tsv": sha256(inputs / "ko_list.tsv"),
            "modules.txt": sha256(inputs / "modules.txt"),
            "kofam_ko_list": sha256(kofam_ko_list),
            **{f"brite/{p.name}": sha256(p) for p in sorted((inputs / "brite").glob("*.json"))},
        },
        "rules": {name: sha256(DATA / name) for name in RULE_FILES},
        "kos": len(kos),
        "kofam_only_kos": len(retired),
        "mapped_kos": len(mapped),
        "pairs": sum(map(len, mapped.values())),
        "dropped_proposals": len(dropped),
    }

    out_dir.mkdir(parents=True, exist_ok=True)
    header = (f"# Built locally by sharur setup-kegg from KEGG data ({provenance['kegg_release']}); "
              "subject to KEGG's terms of use. Do not redistribute.\n")
    with open(out_dir / "kegg_predicates.tsv", "w") as out, \
            open(out_dir / "kegg_evidence_snapshot.tsv", "w") as snap:
        out.write(header + "# ko\tpredicates\tevidence\n")
        snap.write(header + "# ko\tsymbols\tname\tbrite\tmodules\n")
        for ko, evidence in mapped.items():
            symbols, name = kos[ko]
            out.write(f"{ko}\t{','.join(sorted(evidence))}\t"
                      f"{';'.join(f'{p}={e}' for p, e in sorted(evidence.items()))}\n")
            brite = "|".join(f"{h}:{PATH_SEP.join(labels)}" for h, labels in placements.get(ko, ()))
            modules_col = ",".join(m + ("+" if c else "") for m, c in memberships.get(ko, ()))
            snap.write(f"{ko}\t{symbols}\t{name}\t{brite}\t{modules_col}\n")
    provenance["map_sha256"] = sha256(out_dir / "kegg_predicates.tsv")
    (out_dir / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    if report:
        with open(report, "w") as rep:
            rep.write("ko\tsymbols\tname\tproposed\tresolved\n")
            for row in dropped:
                rep.write("\t".join(row) + "\n")
    log(f"KOs: {len(kos):,}; mapped: {len(mapped):,} KOs, {provenance['pairs']:,} pairs; "
        f"proposed pairs without evidence: {len(dropped):,}")
    return provenance


def setup(out_dir: Path, inputs: Path | None = None, kofam_ko_list: Path | None = None,
          report: Path | None = None, log=print) -> dict:
    """Fetch (unless ``inputs`` already holds KEGG files) and build into ``out_dir``."""
    fetched = inputs is None
    inputs = inputs or out_dir / "inputs"
    if fetched:
        fetch_inputs(inputs, kofam_ko_list, log=log)
    return build(inputs, out_dir, kofam_ko_list=kofam_ko_list, report=report, log=log)
