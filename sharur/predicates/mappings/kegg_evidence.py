"""Evidence rules tying KO-derived predicates to KEGG's own information.

A (KO, predicate) pair is supported by one of:

- ``ec:<EC>``: KEGG's KO definition lists the EC number, and its class maps to
  the predicate through :data:`sharur.predicates.mappings.kegg_map.EC_TO_PREDICATES`;
- ``brite:<hierarchy> <node path>``: KEGG places the KO under a BRITE node that
  ``data/kegg_brite_predicates.tsv`` maps to the predicate;
- ``module:<module>``: the KO is in a KEGG module that
  ``data/kegg_module_predicates.tsv`` maps to the predicate, or, for
  ``complex_subunit``, is joined to other components by '+' or '-' in the
  module definition;
- ``hyddb:k/n``: KEGG names the KO a hydrogenase (EC 1.12, or "hydrogenase"/
  "hydrogenlyase" in its name), and k of the n HydDB reference hydrogenases
  that score at or above the KO's KOfam threshold carry a HydDB label whose
  Table 1 interpretation (:mod:`sharur.hydrogenase.subgroups`) gives the
  predicate; consensus thresholds match the Swiss-Prot tier;
- ``text:<match>``: the KO's symbols or name state the predicate (the Pfam
  text conventions in :mod:`sharur.predicates.mappings.pfam_evidence_spec`).

``scripts/build_kegg_predicate_map.py`` applies these rules to build the shipped
map; ``tests/test_kegg_map_integrity.py`` re-verifies every shipped pair.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from sharur.hydrogenase.classifier import parse_reference_id
from sharur.hydrogenase.subgroups import SUBGROUPS
from sharur.predicates.mappings.pfam_evidence import (
    SWISSPROT_COVERAGE,
    SWISSPROT_MIN_LOWER_BOUND,
    text_evidence,
    wilson_lower_bound,
)


DATA = Path(__file__).with_name("data")
PATH_SEP = " > "
EC_SUFFIX = re.compile(r"\s*\[EC:[^\]]*\]")


@dataclass(frozen=True)
class BriteRule:
    hierarchy: str
    levels: tuple[str, ...]
    definition: re.Pattern | None
    predicates: tuple[str, ...]

    def matches(self, hierarchy: str, labels: tuple[str, ...], definition: str) -> bool:
        if hierarchy != self.hierarchy or len(labels) < len(self.levels):
            return False
        for level, label in zip(self.levels, labels, strict=False):
            if level.startswith("re:"):
                if not re.fullmatch(level[3:], label):
                    return False
            elif not label.startswith(level):
                return False
        return self.definition is None or bool(self.definition.search(definition))


def _rows(path: Path) -> list[list[str]]:
    return [line.rstrip("\n").split("\t") for line in path.read_text().splitlines()
            if line.strip() and not line.startswith("#")]


def load_brite_rules(path: Path = DATA / "kegg_brite_predicates.tsv") -> list[BriteRule]:
    rules = []
    for hierarchy, node_path, definition, preds in _rows(path):
        rules.append(BriteRule(
            hierarchy,
            tuple(node_path.split(PATH_SEP)) if node_path else (),
            re.compile(definition) if definition else None,
            tuple(preds.split()),
        ))
    return rules


def load_module_rules(path: Path = DATA / "kegg_module_predicates.tsv") -> dict[str, tuple[str, tuple[str, ...]]]:
    return {module: (name, tuple(preds.split())) for module, name, preds in _rows(path)}


def ko_definition(symbols: str, name: str) -> str:
    """KEGG's KO definition as written in KEGG REST ``list/ko``."""
    return f"{symbols}; {name}" if symbols else name


def brite_evidence(rules: list[BriteRule], placements, definition: str) -> dict[str, str]:
    """Predicates for a KO from its BRITE placements [(hierarchy, labels)], with evidence."""
    found: dict[str, str] = {}
    for hierarchy, labels in placements:
        for rule in rules:
            if rule.matches(hierarchy, labels, definition):
                for pred in rule.predicates:
                    found.setdefault(pred, f"brite:{hierarchy} {PATH_SEP.join(labels)}")
    return found


def module_evidence(rules, memberships) -> dict[str, str]:
    """Predicates for a KO from its module memberships [(module, is_complex_component)]."""
    found: dict[str, str] = {}
    for module, in_complex in memberships:
        for pred in rules.get(module, ("", ()))[1]:
            found.setdefault(pred, f"module:{module}")
        if in_complex:
            found.setdefault("complex_subunit", f"module:{module}")
    return found


def ko_text_evidence(predicate: str, symbols: str, name: str) -> str | None:
    return text_evidence(predicate, symbols, EC_SUFFIX.sub("", name))


HYDROGENASE_NAME = re.compile(r"(?<![a-z])(hydrogenase|hydrogenlyase)|EC:[^\]]*\b1\.12\.", re.IGNORECASE)
HYDDB_TYPE_PREDICATES = {"NiFe": "nife_hydrogenase", "FeFe": "fefe_hydrogenase", "Fe": "fe_only_hydrogenase"}


def names_hydrogenase(name: str) -> bool:
    """KEGG's KO name states hydrogenase activity (dehydrogenases and transhydrogenases excluded)."""
    return bool(HYDROGENASE_NAME.search(name))


def hyddb_label_predicates(label: str) -> tuple[str, ...]:
    """Predicates HydDB Table 1 supports for a reference label such as ``[NiFe]_Group_4a``."""
    _, _, _, hyd_type, subgroup = parse_reference_id(f"x|x|{label}")
    found = SUBGROUPS.get((hyd_type, subgroup))
    if found is None:
        return ()
    return ("hydrogenase", HYDDB_TYPE_PREDICATES[hyd_type], *found.predicates)


def hyddb_evidence(label_counts: dict[str, int]) -> dict[str, str]:
    """Predicates on which the HydDB references inside a KO agree, with k/n evidence."""
    n = sum(label_counts.values())
    min_n, min_frac = SWISSPROT_COVERAGE
    if n < min_n:
        return {}
    support: dict[str, int] = {}
    for label, count in label_counts.items():
        for pred in set(hyddb_label_predicates(label)):
            support[pred] = support.get(pred, 0) + count
    return {pred: f"hyddb:{k}/{n}" for pred, k in sorted(support.items())
            if k >= min_frac * n and wilson_lower_bound(k, n) >= SWISSPROT_MIN_LOWER_BOUND}


def load_hyddb_snapshot(path: Path = DATA / "kegg_hyddb_snapshot.tsv") -> dict[str, dict[str, int]]:
    snapshot = {}
    for ko, _, labels in _rows(path):
        snapshot[ko] = {lab: int(n) for lab, n in (x.rsplit("=", 1) for x in labels.split("|") if x)}
    return snapshot


_TOKEN = re.compile(r"K\d{5}|M\d{5}|--|[()+\-, ]")


def parse_module_definition(definition: str) -> dict[str, bool]:
    """KOs in a KEGG module definition -> whether each is a component of a complex.

    Grammar: steps separated by spaces, alternatives by ',', complex components
    by '+' (or '-' for optional components); parentheses group. A KO is a
    complex component when it sits, at any depth, inside a '+'/'-' joined term
    with two or more components.
    """
    tokens = _TOKEN.findall(definition)
    pos = 0
    result: dict[str, bool] = {}

    def expr(stop):
        nonlocal pos
        kos: list[list[str]] = []
        while pos < len(tokens) and tokens[pos] not in stop:
            if tokens[pos] in (" ", ","):
                pos += 1
                continue
            kos.append(term())
        return [k for group in kos for k in group]

    def term():
        nonlocal pos
        components = []
        if tokens[pos] == "-":  # leading optional component
            pos += 1
        components.append(component())
        while pos < len(tokens) and tokens[pos] in ("+", "-"):
            pos += 1
            components.append(component())
        found = [k for c in components for k in c]
        if len(components) > 1:
            for k in found:
                result[k] = True
        return found

    def component():
        nonlocal pos
        tok = tokens[pos]
        pos += 1
        if tok == "(":
            inner = expr({")"})
            pos += 1
            return inner
        if tok.startswith("K"):
            result.setdefault(tok, False)
            return [tok]
        return []  # module reference or '--' placeholder

    expr(set())
    return result
