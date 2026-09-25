"""KEGG module completeness for genomes and loci.

Module definitions come from the local KEGG build (``sharur setup-kegg`` writes
``kegg_modules.tsv``; KEGG text stays on the user's machine). A definition is
a space-separated sequence of steps; within a step ``,`` separates
alternatives, ``+`` joins required complex components and ``-`` marks optional
ones; parentheses nest; ``Mxxxxx`` refers to another module and ``--`` marks a
step without a KO (left out of the count).

Scores run from 0 to 1: a KO is 1 when present; a complex scores the fraction
of its required components present; alternatives take their best option; a
sequence averages its steps. A step is complete at 1. Module completeness is
the fraction of complete top-level steps; ``score`` keeps partial credit.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from functools import cache
from pathlib import Path
from typing import Any

from sharur.predicates.mappings.kegg_map import kegg_dir


_TOKEN = re.compile(r"K\d{5}|M\d{5}|--|[()+\-, ]")


# --------------------------------------------------------------------------- #
# Grammar
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Node:
    kind: str  # ko | module | undefined | seq | alt | complex
    value: str = ""
    children: tuple[Node, ...] = ()
    optional: tuple[bool, ...] = ()  # complex: per-child optional flag


def parse(definition: str) -> Node:
    """Parse a KEGG module definition into a sequence of steps."""
    tokens = _TOKEN.findall(definition)
    pos = 0

    def seq(stop: str | None) -> Node:
        nonlocal pos
        items = []
        while pos < len(tokens) and tokens[pos] != stop:
            if tokens[pos] == " ":
                pos += 1
                continue
            items.append(alt(stop))
        return Node("seq", children=tuple(items))

    def alt(stop: str | None) -> Node:
        nonlocal pos
        options = [complex_()]
        while pos < len(tokens) and tokens[pos] == ",":
            pos += 1
            options.append(complex_())
        return options[0] if len(options) == 1 else Node("alt", children=tuple(options))

    def complex_() -> Node:
        nonlocal pos
        parts, optional = [], []
        leading_optional = pos < len(tokens) and tokens[pos] == "-" and (pos + 1 < len(tokens) and tokens[pos + 1] != "-")
        if leading_optional:
            pos += 1
        parts.append(unit())
        optional.append(leading_optional)
        while pos < len(tokens) and tokens[pos] in ("+", "-"):
            optional.append(tokens[pos] == "-")
            pos += 1
            parts.append(unit())
        if len(parts) == 1 and not optional[0]:
            return parts[0]
        return Node("complex", children=tuple(parts), optional=tuple(optional))

    def unit() -> Node:
        nonlocal pos
        token = tokens[pos]
        pos += 1
        if token == "(":
            inner = seq(")")
            pos += 1
            return inner.children[0] if len(inner.children) == 1 else inner
        if token.startswith("K"):
            return Node("ko", token)
        if token.startswith("M"):
            return Node("module", token)
        return Node("undefined")

    return seq(None)


def kos_in(node: Node) -> set[str]:
    if node.kind == "ko":
        return {node.value}
    return {k for child in node.children for k in kos_in(child)}


# --------------------------------------------------------------------------- #
# Scoring
# --------------------------------------------------------------------------- #


@dataclass
class ModuleDefinition:
    module: str
    name: str
    module_class: str
    definition: str

    @property
    def tree(self) -> Node:
        return _parsed(self.definition)


@cache
def _parsed(definition: str) -> Node:
    return parse(definition)


def score(node: Node, present: set[str], modules: dict[str, ModuleDefinition], depth: int = 0) -> float | None:
    """Fractional presence of a node (None: undefined step, left out)."""
    if node.kind == "ko":
        return 1.0 if node.value in present else 0.0
    if node.kind == "undefined":
        return None
    if node.kind == "module":
        ref = modules.get(node.value)
        if ref is None or depth > 5:
            return 0.0
        return module_score(ref, present, modules, depth + 1)
    if node.kind == "alt":
        values = [v for v in (score(c, present, modules, depth) for c in node.children) if v is not None]
        return max(values) if values else None
    if node.kind == "complex":
        required = [score(c, present, modules, depth) for c, opt in zip(node.children, node.optional, strict=True) if not opt]
        required = [v for v in required if v is not None]
        return sum(required) / len(required) if required else None
    values = [v for v in (score(c, present, modules, depth) for c in node.children) if v is not None]
    return sum(values) / len(values) if values else None


def module_score(module: ModuleDefinition, present: set[str], modules: dict[str, ModuleDefinition],
                 depth: int = 0) -> float:
    value = score(module.tree, present, modules, depth)
    return value or 0.0


def _missing(node: Node, present: set[str]) -> list[str]:
    """KOs of the best-scoring option that are absent (what would complete this step)."""
    if node.kind == "ko":
        return [] if node.value in present else [node.value]
    if node.kind == "alt":
        best = max(node.children, key=lambda c: score(c, present, {}) or 0.0)
        return _missing(best, present)
    if node.kind == "complex":
        return [k for c, opt in zip(node.children, node.optional, strict=True) if not opt for k in _missing(c, present)]
    return [k for c in node.children for k in _missing(c, present)]


@dataclass
class StepResult:
    index: int
    score: float
    found: list[str]
    missing: list[str]


@dataclass
class ModuleResult:
    module: str
    name: str
    module_class: str
    completeness: float
    score: float
    steps_complete: int
    steps_total: int
    steps: list[StepResult] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "module": self.module, "name": self.name, "class": self.module_class,
            "completeness": round(self.completeness, 4), "score": round(self.score, 4),
            "steps_complete": self.steps_complete, "steps_total": self.steps_total,
            "steps": [{"step": s.index, "score": round(s.score, 4), "found": s.found, "missing": s.missing}
                      for s in self.steps],
        }


def evaluate(module: ModuleDefinition, present: set[str], modules: dict[str, ModuleDefinition]) -> ModuleResult:
    """Stepwise completeness of one module given the KOs present."""
    steps: list[StepResult] = []
    for i, step in enumerate(module.tree.children, 1):
        value = score(step, present, modules)
        if value is None:
            continue
        found = sorted(kos_in(step) & present)
        steps.append(StepResult(i, value, found, [] if value >= 1.0 else sorted(set(_missing(step, present)))))
    complete = sum(1 for s in steps if s.score >= 1.0)
    total = len(steps)
    return ModuleResult(module.module, module.name, module.module_class,
                        complete / total if total else 0.0,
                        sum(s.score for s in steps) / total if total else 0.0,
                        complete, total, steps)


# --------------------------------------------------------------------------- #
# Data
# --------------------------------------------------------------------------- #


def load_modules(directory: Path | None = None) -> dict[str, ModuleDefinition]:
    """Module definitions from the local KEGG build (empty when absent)."""
    directory = directory or kegg_dir()
    path = (directory / "kegg_modules.tsv") if directory else None
    if path is None or not path.exists():
        return {}
    modules = {}
    for line in path.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        module, name, module_class, definition = line.split("\t")
        modules[module] = ModuleDefinition(module, name, module_class, definition)
    return modules


def _ko_hits(store, where: str = "", params: list | None = None) -> list[tuple[str, str, str]]:
    """(bin_id, KO, protein_id) for KOfam/KEGG annotations."""
    return store.execute(
        f"""SELECT p.bin_id, a.accession, a.protein_id FROM annotations a
            JOIN proteins p USING (protein_id)
            WHERE LOWER(a.source) IN ('kofam', 'kegg') AND regexp_matches(a.accession, '^K[0-9]{{5}}$') {where}""",
        params or [])


def genome_modules(store, *, bins: list[str] | None = None, modules: list[str] | None = None,
                   min_completeness: float = 0.0, directory: Path | None = None) -> list[dict[str, Any]]:
    """Module completeness per genome, with the proteins behind each found KO."""
    definitions = load_modules(directory)
    if not definitions:
        raise FileNotFoundError("No local KEGG module definitions; run `sharur setup-kegg`.")
    where, params = "", []
    if bins:
        where = f"AND p.bin_id IN ({','.join('?' * len(bins))})"
        params = list(bins)
    by_bin: dict[str, dict[str, list[str]]] = {}
    for bin_id, ko, protein in _ko_hits(store, where, params):
        by_bin.setdefault(bin_id, {}).setdefault(ko, []).append(protein)
    bin_quality = {}
    columns = {r[0] for r in store.execute(
        "SELECT column_name FROM information_schema.columns WHERE table_name = 'bins'")}
    if {"completeness", "contamination"} <= columns:
        bin_quality = {b: (c, x) for b, c, x in store.execute("SELECT bin_id, completeness, contamination FROM bins")}
    wanted = modules or sorted(definitions)
    rows = []
    for bin_id in sorted(set(bins or []) | set(by_bin)):
        present = set(by_bin.get(bin_id, {}))
        for module in wanted:
            if module not in definitions:
                continue
            result = evaluate(definitions[module], present, definitions)
            if result.completeness < min_completeness or (result.score == 0 and not modules):
                continue
            row = result.to_dict()
            for step in row["steps"]:
                step["proteins"] = {ko: sorted(by_bin[bin_id][ko]) for ko in step["found"]}
            completeness, contamination = bin_quality.get(bin_id, (None, None))
            row.update({"bin_id": bin_id, "bin_completeness": completeness, "bin_contamination": contamination})
            rows.append(row)
    return rows


def locus_modules(store, protein_id: str, window: int = 10, directory: Path | None = None) -> list[dict[str, Any]]:
    """Modules with steps encoded within ``window`` genes of a protein on the same contig."""
    definitions = load_modules(directory)
    if not definitions:
        raise FileNotFoundError("No local KEGG module definitions; run `sharur setup-kegg`.")
    anchor = store.execute("SELECT contig_id, gene_index FROM proteins WHERE protein_id = ?", [protein_id])
    if not anchor or anchor[0][1] is None:
        return []
    contig, gene_index = anchor[0]
    hits = store.execute(
        """SELECT a.accession, a.protein_id, p.gene_index FROM annotations a JOIN proteins p USING (protein_id)
           WHERE p.contig_id = ? AND p.gene_index BETWEEN ? AND ?
             AND LOWER(a.source) IN ('kofam', 'kegg') AND regexp_matches(a.accession, '^K[0-9]{5}$')""",
        [contig, gene_index - window, gene_index + window])
    present: dict[str, list[tuple[str, int]]] = {}
    for ko, pid, gi in hits:
        present.setdefault(ko, []).append((pid, gi - gene_index))
    rows = []
    for _, definition in sorted(definitions.items()):
        if not kos_in(definition.tree) & set(present):
            continue
        row = evaluate(definition, set(present), definitions).to_dict()
        for step in row["steps"]:
            step["proteins"] = {ko: [f"{pid} ({offset:+d})" for pid, offset in present[ko]] for ko in step["found"]}
        rows.append(row)
    return sorted(rows, key=lambda r: (-r["completeness"], -r["score"], r["module"]))


def modules_markdown(rows: list[dict[str, Any]], title: str) -> str:
    lines = [f"# {title}"]
    if not rows:
        lines.append("No module steps found.")
    for r in rows:
        quality = ""
        if r.get("bin_completeness") is not None:
            quality = f" [bin {r['bin_completeness']}% complete]"
        lines.append(f"- {r.get('bin_id', '') + ' ' if r.get('bin_id') else ''}{r['module']} {r['name']}: "
                     f"{r['steps_complete']}/{r['steps_total']} steps ({r['completeness']:.0%}; score {r['score']:.2f})"
                     + quality)
        incomplete = [s for s in r["steps"] if s["score"] < 1.0]
        if incomplete and r["completeness"] >= 0.5:
            lines.append("  missing: " + "; ".join(
                f"step {s['step']} needs {'/'.join(s['missing'][:4])}" for s in incomplete[:6]))
    return "\n".join(lines)
