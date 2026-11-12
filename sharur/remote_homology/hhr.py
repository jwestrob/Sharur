"""Parse HH-suite result files (``.hhr``) and classify hits; standard library only.

An ``.hhr`` file holds a header (``Query``, ``Match_columns``, ``No_of_seqs``, ``Neff``,
``Searched_HMMs``, ``Date``, ``Command``), a summary table with one line per hit, and an
alignment block per hit (``No N``, ``>full name``, ``Probab=... E-value=... Aligned_cols=...``,
then query and template alignment lines). The summary line gives the query and template
ranges and the template length; the alignment block gives the untruncated hit name and the
identity. :func:`parse_hhr` merges both by hit number. Alignment residues are never read
into the result.

Hit names differ by database:

``pfam``    ``PF00004.33 AAA; ATPase family`` (accession with version, ID, description), the
            ``-name`` that ``hhmake`` received when the Pfam seed alignments were converted
``pdb``     ``1abc_A`` or ``1abc_A mol:protein length:154  TITLE``; titles come from a
            ``pdb id -> title`` table when the name has none
``uniref``  ``UniRef100_X Name n=1 Tax=... TaxID=... RepID=...``
other       first word accession, the rest description
"""

from __future__ import annotations

import re
from collections.abc import Iterable
from dataclasses import asdict, dataclass, field
from pathlib import Path


# Family-tier sources and MSA kinds (family_annotations.source / .msa_kind).
SOURCE_PFAM = "hhsearch_pfam"
SOURCE_PDB = "hhsearch_pdb"
SOURCE_PDB_HHBLITS = "hhblits_pdb"     # PDB searched with HHblits (prefilter, one iteration)
SOURCE_UNIREF = "hhblits_uniref"
FAMILY_SOURCES = (SOURCE_PFAM, SOURCE_PDB, SOURCE_PDB_HHBLITS, SOURCE_UNIREF)
RULE_SOURCES = (SOURCE_PFAM, SOURCE_PDB, SOURCE_PDB_HHBLITS)
MSA_ENRICHED = "enriched"
MSA_FAMILY_ONLY = "family_only"
MSA_KINDS = (MSA_ENRICHED, MSA_FAMILY_ONLY)
STORE_MIN_PROBABILITY = 50.0

# Hit names that describe a family of unknown function: Pfam DUF/UPF families and
# "uncharacterized", "hypothetical" or "unknown" proteins and structures.
UNKNOWN_FUNCTION_PATTERN = (r"DUF\d|UPF\d|unknown function|uncharacteri|hypothetical|function unknown"
                            r"|unknown protein")
UNKNOWN_FUNCTION = re.compile(UNKNOWN_FUNCTION_PATTERN, re.IGNORECASE)

_SUMMARY = re.compile(
    r"^\s*(?P<rank>\d+)\s(?P<name>.*?)\s+(?P<prob>\d+(?:\.\d+)?)\s+(?P<evalue>\S+)\s+(?P<pvalue>\S+)"
    r"\s+(?P<score>-?\d+(?:\.\d+)?)\s+(?P<ss>-?\d+(?:\.\d+)?)\s+(?P<cols>\d+)\s+(?P<qs>\d+)-(?P<qe>\d+)"
    r"\s+(?P<ts>\d+)-(?P<te>\d+)\s*\(\s*(?P<tlen>\d+)\s*\)\s*$")
_KV = re.compile(r"(\w+)=(\S+)")
_PDB_ID = re.compile(r"^(?P<code>[0-9][A-Za-z0-9]{3})_(?P<chain>\S+)$")
_PDB_META = re.compile(r"\bmol:\S+\s+length:\d+\s*")
_PDB_RANGE = re.compile(r"^\d+-\d+(?:,\s*\d+-\d+)*(?=\s|$)")


def is_unknown_function(*texts: str | None) -> bool:
    """True when any of ``texts`` names a family of unknown function."""
    return any(bool(t) and UNKNOWN_FUNCTION.search(t) is not None for t in texts)


@dataclass
class HHRHit:
    rank: int
    label: str                    # the full hit name as HH-suite wrote it
    probability: float            # percent
    evalue: float
    pvalue: float | None
    score: float
    aligned_cols: int
    q_start: int
    q_end: int
    t_start: int
    t_end: int
    t_length: int | None
    identities: float | None = None   # percent, from the alignment block
    similarity: float | None = None
    sum_probs: float | None = None
    template_neff: float | None = None

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class HHRResult:
    query: str | None = None
    match_columns: int | None = None
    no_of_seqs: str | None = None
    neff: float | None = None
    searched_hmms: int | None = None
    command: str | None = None
    hits: list[HHRHit] = field(default_factory=list)


def _float(text: str | None) -> float | None:
    if text is None:
        return None
    try:
        return float(text.rstrip("%"))
    except ValueError:
        return None


def parse_hhr_text(text: str) -> HHRResult:
    """Parse the text of one ``.hhr`` file."""
    result = HHRResult()
    lines = text.splitlines()
    i = 0
    # header
    while i < len(lines):
        line = lines[i]
        if line.startswith(" No Hit"):
            i += 1
            break
        key, _, value = line.partition(" ")
        value = value.strip()
        if key == "Query":
            result.query = value.split()[0] if value else None
        elif key == "Match_columns":
            result.match_columns = int(value)
        elif key == "No_of_seqs":
            result.no_of_seqs = value
        elif key == "Neff":
            result.neff = _float(value)
        elif key == "Searched_HMMs":
            result.searched_hmms = int(value)
        elif key == "Command":
            result.command = value
        i += 1
    # summary table
    hits: dict[int, HHRHit] = {}
    while i < len(lines):
        line = lines[i]
        if not line.strip():
            i += 1
            break
        m = _SUMMARY.match(line)
        if m is None:
            raise ValueError(f"unreadable .hhr summary line: {line[:60]!r}")
        rank = int(m["rank"])
        hits[rank] = HHRHit(
            rank=rank, label=m["name"].strip(), probability=float(m["prob"]), evalue=float(m["evalue"]),
            pvalue=_float(m["pvalue"]), score=float(m["score"]), aligned_cols=int(m["cols"]),
            q_start=int(m["qs"]), q_end=int(m["qe"]), t_start=int(m["ts"]), t_end=int(m["te"]),
            t_length=int(m["tlen"]))
        i += 1
    # alignment blocks: full names and identities
    current: HHRHit | None = None
    while i < len(lines):
        line = lines[i]
        if line.startswith("No ") and line[3:].strip().isdigit():
            current = hits.get(int(line[3:].strip()))
        elif current is not None and line.startswith(">"):
            current.label = line[1:].strip()
        elif current is not None and line.startswith("Probab="):
            kv = dict(_KV.findall(line))
            current.probability = _float(kv.get("Probab")) or current.probability
            current.evalue = _float(kv.get("E-value")) if kv.get("E-value") else current.evalue
            current.aligned_cols = int(kv.get("Aligned_cols", current.aligned_cols))
            current.identities = _float(kv.get("Identities"))
            current.similarity = _float(kv.get("Similarity"))
            current.sum_probs = _float(kv.get("Sum_probs"))
            current.template_neff = _float(kv.get("Template_Neff"))
            current = None
        i += 1
    result.hits = [hits[k] for k in sorted(hits)]
    return result


def parse_hhr(path: str | Path) -> HHRResult:
    """Parse one ``.hhr`` file."""
    return parse_hhr_text(Path(path).read_text(errors="replace"))


# --------------------------------------------------------------------------- #
# Hit names by database kind
# --------------------------------------------------------------------------- #
def split_label(label: str, kind: str, titles: dict[str, str] | None = None) -> tuple[str, str, str]:
    """``(accession, name, description)`` of a hit name for database ``kind``."""
    label = label.strip()
    if kind == "pfam":
        head, _, description = label.partition(";")
        parts = head.split(None, 1)
        accession = parts[0].split(".")[0] if parts else label
        name = parts[1].strip() if len(parts) > 1 else accession
        return accession, name, description.strip()
    if kind == "pdb":
        parts = label.split(None, 1)
        chain = parts[0] if parts else label
        rest = _PDB_META.sub("", parts[1]).strip() if len(parts) > 1 else ""
        rest = _PDB_RANGE.sub("", rest).strip()  # a residue range ("1-233") is no description
        if not rest and titles:
            m = _PDB_ID.match(chain)
            if m:
                rest = titles.get(m["code"].lower(), "")
        return chain, chain, rest
    parts = label.split(None, 1)
    accession = parts[0] if parts else label
    description = parts[1].strip() if len(parts) > 1 else ""
    if kind == "uniref":
        description = re.split(r"\s+n=\d+\s+Tax=", description)[0].strip()
    return accession, accession, description


def read_titles(path: str | Path | None) -> dict[str, str]:
    """``pdb code (lower case) -> title`` from a two-column TSV (``code<TAB>title``)."""
    if not path:
        return {}
    titles: dict[str, str] = {}
    with open(path, errors="replace") as handle:
        for line in handle:
            code, _, title = line.rstrip("\n").partition("\t")
            if code and title:
                titles.setdefault(code.lower()[:4], title)
    return titles


HIT_COLUMNS = ("source", "msa_kind", "hit_rank", "hit_accession", "hit_name", "hit_description", "probability",
               "evalue", "score", "aligned_cols", "identities", "q_start", "q_end", "t_start", "t_end", "t_length",
               "unknown_function", "hit_label")


def hit_rows(result: HHRResult, *, source: str, msa_kind: str, kind: str, min_probability: float,
             max_hits: int | None = None, titles: dict[str, str] | None = None) -> Iterable[tuple]:
    """Rows in :data:`HIT_COLUMNS` order for hits at ``min_probability`` or above (first ``max_hits``)."""
    kept = 0
    for hit in result.hits:
        if hit.probability < min_probability:
            continue
        if max_hits is not None and kept >= max_hits:
            break
        kept += 1
        accession, name, description = split_label(hit.label, kind, titles)
        yield (source, msa_kind, hit.rank, accession, name, description, hit.probability, hit.evalue, hit.score,
               hit.aligned_cols, hit.identities, hit.q_start, hit.q_end, hit.t_start, hit.t_end, hit.t_length,
               is_unknown_function(name, description), hit.label)


__all__ = [
    "FAMILY_SOURCES",
    "HIT_COLUMNS",
    "MSA_ENRICHED",
    "MSA_FAMILY_ONLY",
    "MSA_KINDS",
    "RULE_SOURCES",
    "SOURCE_PDB",
    "SOURCE_PDB_HHBLITS",
    "SOURCE_PFAM",
    "SOURCE_UNIREF",
    "STORE_MIN_PROBABILITY",
    "UNKNOWN_FUNCTION",
    "UNKNOWN_FUNCTION_PATTERN",
    "HHRHit",
    "HHRResult",
    "hit_rows",
    "is_unknown_function",
    "parse_hhr",
    "parse_hhr_text",
    "read_titles",
    "split_label",
]
