"""Domain architectures: ordered domains per protein and pattern search over them.

A protein's architecture is its domain hits with protein coordinates from one
or more sources (Pfam by default), overlaps resolved (best E-value first; a hit overlapping an accepted one by more than
``max_overlap`` of the shorter domain is dropped), ordered by start.

Patterns are whitespace-separated tokens matched against that order:

- ``Big_2``, ``PF00092``: one domain, by name or accession (case-sensitive);
- ``Big_*``: a glob over names (``*`` any characters, ``?`` one character);
- ``.``: any one domain;
- quantifiers after a token or group: ``?``, ``*``, ``+``, ``{m}``, ``{m,}``, ``{m,n}``
  (``?`` and ``*`` need a space before them, since attached they are globs);
- groups with alternatives: ``( Cadherin | Big_2 )+``;
- anchors: ``^`` (N-terminal end) and ``$`` (C-terminal end).

Without anchors a pattern may match anywhere in the architecture, e.g.
``Big_2{5,} . VWA`` finds five or more Big_2 domains followed by any one domain
and a VWA domain.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any

MAX_OVERLAP = 0.5
_QUANTIFIER = re.compile(r"^(\?|\*|\+|\{\d+(,\d*)?\})$")
_TOKEN = re.compile(r"\(|\)|\||\{\d+(?:,\d*)?\}|[?*+](?=\s|$|\))|[^\s()|{}]+")


class PatternError(ValueError):
    """The architecture pattern cannot be parsed."""


@dataclass(frozen=True)
class Domain:
    name: str
    accession: str
    start: int | None
    end: int | None
    evalue: float | None
    source: str

    def to_dict(self) -> dict[str, Any]:
        return {"name": self.name, "accession": self.accession, "start_aa": self.start,
                "end_aa": self.end, "evalue": self.evalue, "source": self.source}


@dataclass(frozen=True)
class CompiledPattern:
    text: str
    regex: re.Pattern
    # Each clause is a set of alternatives (names, accessions or globs); a
    # matching protein carries at least one hit per clause.
    clauses: tuple[frozenset[str], ...]
    min_domains: int

    @property
    def required(self) -> frozenset[str]:
        return frozenset(next(iter(c)) for c in self.clauses if len(c) == 1)


# --------------------------------------------------------------------------- #
# Patterns
# --------------------------------------------------------------------------- #


def _tokenize(pattern: str) -> list[str]:
    # ``Big_2+`` means one or more Big_2 (Pfam names contain no '+'); ``*`` and
    # ``?`` attached to a name are globs, so quantify with a space: ``Big_2 *``.
    pattern = re.sub(r"(?<=[^\s(|])\+", " +", pattern)
    tokens = _TOKEN.findall(pattern)
    if "".join(tokens) != re.sub(r"\s+", "", pattern):
        raise PatternError(f"Cannot parse pattern: {pattern!r}")
    return tokens


# Each domain is encoded as ``NAME\x1fACCESSION;`` so a token can name either.
_SEP = "\x1f"


def _is_glob(token: str) -> bool:
    return any(ch in token for ch in "*?[")


def _atom(token: str) -> str:
    if token == ".":
        return "(?:[^;]+;)"
    if _is_glob(token):
        body = "".join("[^;\x1f]*" if ch == "*" else "[^;\x1f]" if ch == "?" else re.escape(ch) for ch in token)
        return f"(?:{body}\x1f[^;]*;)"
    t = re.escape(token)
    return f"(?:{t}\x1f[^;]*;|[^;\x1f]*\x1f{t}(?:\\.\\d+)?;)"


@dataclass
class _Item:
    regex: str
    min_len: int
    atom: str | None = None          # single token (None for groups)
    alternatives: list[list["_Item"]] | None = None
    low: int = 1


def _quantify(item: _Item, token: str) -> None:
    low = {"?": 0, "*": 0, "+": 1}.get(token)
    if low is None:
        low = int(token[1:-1].split(",")[0])
    item.low = low
    item.regex += token


def _parse(tokens: list[str]) -> tuple[list[_Item], bool, bool]:
    anchored_start = bool(tokens) and tokens[0] == "^"
    anchored_end = bool(tokens) and tokens[-1] == "$"
    body = tokens[1 if anchored_start else 0: len(tokens) - (1 if anchored_end else 0)]
    if "^" in body or "$" in body:
        raise PatternError("'^' must start the pattern and '$' must end it")
    position = 0

    def sequence(depth: int) -> list[list[_Item]]:
        nonlocal position
        branches: list[list[_Item]] = [[]]
        while position < len(body):
            token = body[position]
            if token == ")":
                if depth == 0:
                    raise PatternError("Unbalanced ')'")
                return branches
            position += 1
            if token == "|":
                if depth == 0:
                    raise PatternError("'|' is only allowed inside ( ... )")
                branches.append([])
            elif token == "(":
                alternatives = sequence(depth + 1)
                if position >= len(body) or body[position] != ")":
                    raise PatternError("Unbalanced '('")
                position += 1
                if any(not alt for alt in alternatives):
                    raise PatternError("Empty alternative in ( ... )")
                regex = "(?:" + "|".join("".join(i.regex for i in alt) for alt in alternatives) + ")"
                branches[-1].append(_Item(regex, min(sum(i.min_len * i.low for i in alt) for alt in alternatives),
                                          alternatives=alternatives))
            elif _QUANTIFIER.match(token):
                if not branches[-1] or branches[-1][-1].regex.endswith(("?", "*", "+", "}")):
                    raise PatternError(f"Quantifier {token!r} has nothing to repeat")
                _quantify(branches[-1][-1], token)
            else:
                branches[-1].append(_Item(_atom(token), 1, atom=token))
        if depth:
            raise PatternError("Unbalanced '('")
        return branches

    items = sequence(0)
    if position != len(body):
        raise PatternError("Unbalanced ')'")
    return items[0], anchored_start, anchored_end


def compile_pattern(pattern: str) -> CompiledPattern:
    """Compile an architecture pattern to a regex over encoded architectures."""
    tokens = _tokenize(pattern.strip())
    if not tokens:
        raise PatternError("Empty pattern")
    items, anchored_start, anchored_end = _parse(tokens)
    if not items:
        raise PatternError("Pattern has no domains")
    body = "".join(i.regex for i in items)
    body = ("^" if anchored_start else "(?:^|(?<=;))") + body + ("$" if anchored_end else "")
    try:
        regex = re.compile(body)
    except re.error as exc:
        raise PatternError(f"Invalid pattern {pattern!r}: {exc}") from exc
    clauses = []
    for item in items:
        if item.low < 1:
            continue
        if item.atom is not None and item.atom != ".":
            clauses.append(frozenset([item.atom]))
        elif item.alternatives and all(len(alt) == 1 and alt[0].atom not in (None, ".") and alt[0].low >= 1
                                       for alt in item.alternatives):
            clauses.append(frozenset(alt[0].atom for alt in item.alternatives))
    return CompiledPattern(pattern, regex, tuple(clauses), sum(i.min_len * i.low for i in items))


# --------------------------------------------------------------------------- #
# Architectures
# --------------------------------------------------------------------------- #


def resolve(hits: list[Domain], max_overlap: float = MAX_OVERLAP) -> list[Domain]:
    """Overlap-resolved domains in N- to C-terminal order (hits without coordinates are left out)."""
    placed = [h for h in hits if h.start is not None and h.end is not None]
    accepted: list[Domain] = []
    for hit in sorted(placed, key=lambda h: (h.evalue if h.evalue is not None else 1.0, h.start)):
        length = hit.end - hit.start + 1
        if all(min(hit.end, a.end) - max(hit.start, a.start) + 1
               <= max_overlap * min(length, a.end - a.start + 1) for a in accepted):
            accepted.append(hit)
    return sorted(accepted, key=lambda h: (h.start, h.end))


def encode(domains: list[Domain]) -> tuple[str, list[int]]:
    """``NAME\x1fACCESSION;...`` string and each domain's character offset."""
    parts, offsets, position = [], [], 0
    for d in domains:
        token = f"{d.name}{_SEP}{d.accession};"
        offsets.append(position)
        parts.append(token)
        position += len(token)
    return "".join(parts), offsets


def compact(names: list[str]) -> str:
    """Run-length display: ``Big_2 x12 - VWA``."""
    runs: list[list[Any]] = []
    for name in names:
        if runs and runs[-1][0] == name:
            runs[-1][1] += 1
        else:
            runs.append([name, 1])
    return " - ".join(f"{n} x{k}" if k > 1 else n for n, k in runs)


def _hits(store, sources: list[str], where: str = "", params: list | None = None) -> dict[str, list[Domain]]:
    placeholders = ",".join("?" * len(sources))
    rows = store.execute(
        f"""SELECT a.protein_id, COALESCE(NULLIF(a.name, ''), a.accession), a.accession, a.start_aa, a.end_aa,
                   a.evalue, LOWER(a.source)
            FROM annotations a WHERE LOWER(a.source) IN ({placeholders}) {where}""",
        [*sources, *(params or [])])
    by_protein: dict[str, list[Domain]] = {}
    for pid, name, accession, start, end, evalue, source in rows:
        by_protein.setdefault(pid, []).append(Domain(name, accession, start, end, evalue, source))
    return by_protein


def architecture(store, protein_id: str, sources: tuple[str, ...] = ("pfam",),
                 max_overlap: float = MAX_OVERLAP) -> list[Domain]:
    """One protein's resolved architecture."""
    hits = _hits(store, [s.lower() for s in sources], "AND a.protein_id = ?", [protein_id])
    return resolve(hits.get(protein_id, []), max_overlap)


def search_architecture(store, pattern: str, *, sources: tuple[str, ...] = ("pfam",),
                        max_overlap: float = MAX_OVERLAP, bins: list[str] | None = None,
                        limit: int | None = 100) -> dict[str, Any]:
    """Proteins whose resolved architecture matches ``pattern``.

    Returns ``{"pattern", "sources", "total", "records"}``; each record has the
    protein, genome, length, compact architecture and the matched domains with
    their coordinates. ``total`` counts all matches; ``records`` stops at ``limit``.
    """
    compiled = compile_pattern(pattern)
    sources_l = [s.lower() for s in sources]
    in_sources = f"LOWER(source) IN ({','.join('?' * len(sources_l))})"
    where, params = "", []
    for clause in compiled.clauses:
        exact = sorted(t for t in clause if not _is_glob(t))
        globs = sorted(t for t in clause if _is_glob(t))
        tests = []
        if exact:
            marks = ",".join("?" * len(exact))
            tests.append(f"name IN ({marks}) OR split_part(accession, '.', 1) IN ({marks})")
            params_clause = [*exact, *exact]
        else:
            params_clause = []
        for glob in globs:
            tests.append("name GLOB ?")
            params_clause.append(glob)
        where += (f" AND a.protein_id IN (SELECT protein_id FROM annotations WHERE {in_sources} "
                  f"AND ({' OR '.join(tests)}))")
        params += [*sources_l, *params_clause]
    if compiled.min_domains > 1:
        where += (f" AND a.protein_id IN (SELECT protein_id FROM annotations WHERE {in_sources} "
                  "AND start_aa IS NOT NULL GROUP BY protein_id HAVING COUNT(*) >= ?)")
        params += [*sources_l, compiled.min_domains]
    if bins:
        where += f" AND a.protein_id IN (SELECT protein_id FROM proteins WHERE bin_id IN ({','.join('?' * len(bins))}))"
        params += list(bins)
    hits = _hits(store, sources_l, where, params)

    matches = []
    for pid in sorted(hits):
        domains = resolve(hits[pid], max_overlap)
        text, offsets = encode(domains)
        found = compiled.regex.search(text)
        if found is None:
            continue
        matched = [d for d, off in zip(domains, offsets, strict=True) if found.start() <= off < found.end()]
        matches.append((pid, domains, matched))
    total = len(matches)
    shown = matches if limit is None else matches[:limit]
    meta = {}
    if shown:
        ids = [m[0] for m in shown]
        for pid, bin_id, length in store.execute(
                "SELECT protein_id, bin_id, sequence_length FROM proteins WHERE protein_id IN "
                f"({','.join('?' * len(ids))})", ids):
            meta[pid] = (bin_id, length)
    records = []
    for pid, domains, matched in shown:
        bin_id, length = meta.get(pid, (None, None))
        records.append({
            "protein_id": pid, "bin_id": bin_id, "length_aa": length,
            "n_domains": len(domains), "architecture": compact([d.name for d in domains]),
            "match": {"domains": [d.to_dict() for d in matched],
                      "start_aa": matched[0].start if matched else None,
                      "end_aa": matched[-1].end if matched else None,
                      "n_domains": len(matched)},
        })
    return {"pattern": pattern, "sources": list(sources_l), "max_overlap": max_overlap,
            "total": total, "records": records}


def architecture_markdown(result: dict[str, Any]) -> str:
    lines = [f"# Architecture search: `{result['pattern']}`",
             f"{result['total']:,} proteins match ({', '.join(result['sources'])}; overlaps above "
             f"{result['max_overlap']:.0%} resolved by E-value)"]
    for r in result["records"]:
        m = r["match"]
        span = f"{m['start_aa']}-{m['end_aa']}" if m["start_aa"] is not None else "unplaced"
        lines.append(f"- {r['protein_id']} ({r['bin_id']}, {r['length_aa']} aa): {r['architecture']}"
                     f"  [match {m['n_domains']} domains, aa {span}]")
    if result["total"] > len(result["records"]):
        lines.append(f"... {result['total'] - len(result['records']):,} more (raise --limit)")
    return "\n".join(lines)
