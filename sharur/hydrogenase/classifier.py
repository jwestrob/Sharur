"""Nearest-reference hydrogenase subgroup assignment against HydDB references.

Procedure: proteins with a HydDB HMM hit (Astra ``hyddb`` annotations) are
searched with DIAMOND (``--ultra-sensitive``, best reference per query) against
the installed HydDB reference sequences, and receive the subgroup of their best
reference. This is a Sharur nearest-reference assignment. The published HydDB
classifier additionally votes over k=4 neighbors, excludes homologous
non-hydrogenase families, and uses downstream genes for FeFe Group A subtypes.

Each protein yields one reconciled :class:`Classification` that keeps the
discovery HMM evidence, the reference match, Pfam domain observations, and the
curation reason separately. Raw ``hyddb`` rows are read and never modified.
"""

from __future__ import annotations

import hashlib
import re
import shutil
import subprocess
import tempfile
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path

from sharur.hydrogenase.subgroups import UNVERIFIED, Subgroup, lookup
from sharur.predicates.pfam_identity import (
    COMPLEX1,
    FE_HYD,
    NIFESE_HASES,
    has_pfam_domain,
    pfam_keys,
)


CLASSIFIER_VERSION = "sharur-hyddb-nearest-reference/2"

HMD = frozenset({"PF03201", "HMD"})

# Outcomes: every protein with a HydDB HMM hit receives exactly one.
ASSIGNED = "assigned"
CLASS_CONFLICT = "class_conflict"
NO_REFERENCE_HIT = "no_reference_hit"
MISSING_SEQUENCE = "missing_sequence"
UNPARSED_REFERENCE_LABEL = "unparsed_reference_label"

CLEARED = "domain_check_cleared"
NEEDS_CURATION = "needs_curation"

CURATION_FLAG = "hyddb_needs_curation"
CLASS_CONFLICT_FLAG = "hyddb_class_conflict"
REVIEW_FLAGS = frozenset({CURATION_FLAG, CLASS_CONFLICT_FLAG})

_DISCOVERY_CLASSES = {"nife": "NiFe", "fefe": "FeFe", "fe_only": "Fe", "fe-only": "Fe", "feonly": "Fe"}
_REFERENCE_LABEL = re.compile(r"^\[(NiFe|FeFe)\]_(Group_\w+)$")


class HydrogenaseSearchError(RuntimeError):
    """The reference search failed; no classification may be written."""


@dataclass(frozen=True)
class ReferenceInfo:
    directory: Path
    release: str
    sha256: str


@dataclass(frozen=True)
class ReferenceHit:
    reference_accession: str
    organism: str
    label: str
    hyd_type: str | None
    subgroup: str | None
    pident: float
    evalue: float
    bitscore: float


@dataclass(frozen=True)
class DomainEvidence:
    nifese_hases: bool = False
    fe_hyd: bool = False
    complex1: bool = False
    hmd: bool = False

    @classmethod
    def from_identifiers(cls, identifiers) -> DomainEvidence:
        ids = tuple(identifiers)
        return cls(
            nifese_hases=has_pfam_domain(ids, NIFESE_HASES),
            fe_hyd=has_pfam_domain(ids, FE_HYD),
            complex1=has_pfam_domain(ids, COMPLEX1),
            hmd=has_pfam_domain(ids, HMD),
        )


@dataclass
class Classification:
    protein_id: str
    outcome: str
    discovery_classes: tuple[str, ...]
    discovery_hit_count: int
    discovery_best_score: float | None
    domains: DomainEvidence
    hit: ReferenceHit | None = None
    subgroup: Subgroup | None = None
    curation_status: str | None = None
    curation_reason: str = ""

    @property
    def derived_labels(self) -> tuple[str, ...]:
        """Predicate names written as ``hyddb_subgroup`` annotations."""
        labels: list[str] = []
        if self.outcome == ASSIGNED and self.subgroup is not None:
            labels.extend(self.subgroup.predicates)
            if self.curation_status == NEEDS_CURATION:
                labels.append(CURATION_FLAG)
        elif self.outcome == CLASS_CONFLICT:
            labels.append(CLASS_CONFLICT_FLAG)
        return tuple(dict.fromkeys(labels))


def _discovery_class(accession: str | None) -> str | None:
    return _DISCOVERY_CLASSES.get((accession or "").strip().lower())


def parse_reference_id(sid: str) -> tuple[str, str, str, str | None, str | None]:
    """Split ``WP_xxx|Organism|[NiFe]_Group_1a`` into accession, organism, label, type, subgroup."""
    parts = sid.split("|")
    accession = parts[0]
    organism = parts[1] if len(parts) > 1 else ""
    label = parts[2] if len(parts) > 2 else ""
    if label == "[Fe]":
        return accession, organism, label, "Fe", "Fe_only"
    match = _REFERENCE_LABEL.match(label)
    if match:
        return accession, organism, label, match.group(1), match.group(2)
    return accession, organism, label, None, None


def find_reference(directory: Path | None = None) -> ReferenceInfo:
    """Locate the installed HydDB DIAMOND reference and fingerprint it."""
    candidates = (
        [Path(directory)]
        if directory is not None
        else [
            Path(__file__).resolve().parents[2] / "data/reference/hyddb",
            Path("data/reference/hyddb"),
            Path.home() / ".sharur/hyddb",
        ]
    )
    for path in candidates:
        dmnd = path / "HydDB_all.dmnd"
        if dmnd.exists():
            releases = sorted(
                {m.group(1) for f in path.glob("*.hmm") if (m := re.search(r"HydDB_(\w+)\.hmm$", f.name))}
            )
            digest = hashlib.sha256(dmnd.read_bytes()).hexdigest()
            return ReferenceInfo(path, ",".join(releases) or "unknown", digest)
    raise FileNotFoundError(
        "HydDB reference not found. Run: diamond makedb --in HydDB_all_hydrogenases.faa --db HydDB_all"
    )


def run_diamond(sequences: dict[str, str], reference: ReferenceInfo, threads: int) -> dict[str, ReferenceHit]:
    """Best HydDB reference per query. Raises :class:`HydrogenaseSearchError` on failure."""
    if not sequences:
        return {}
    if shutil.which("diamond") is None:
        raise HydrogenaseSearchError("diamond not found on PATH")

    with tempfile.TemporaryDirectory(prefix="sharur-hyddb-") as tmp:
        query = Path(tmp) / "query.faa"
        with query.open("w") as handle:
            for pid, seq in sequences.items():
                handle.write(f">{pid}\n{seq}\n")
        result = subprocess.run(
            [
                "diamond", "blastp",
                "--db", str(reference.directory / "HydDB_all.dmnd"),
                "--query", str(query),
                "--outfmt", "6", "qseqid", "sseqid", "pident", "evalue", "bitscore",
                "--max-target-seqs", "1",
                "--threads", str(threads),
                "--ultra-sensitive",
            ],
            capture_output=True,
            text=True,
            check=False,
        )
    if result.returncode != 0:
        raise HydrogenaseSearchError(f"DIAMOND exited {result.returncode}: {result.stderr.strip()[-2000:]}")

    hits: dict[str, ReferenceHit] = {}
    for line in result.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        qid, sid, pident, evalue, bitscore = parts[:5]
        accession, organism, label, hyd_type, subgroup = parse_reference_id(sid)
        hit = ReferenceHit(accession, organism, label, hyd_type, subgroup,
                           float(pident), float(evalue), float(bitscore))
        # Several HSPs can be reported for the one best reference; keep the strongest.
        if qid not in hits or hit.bitscore > hits[qid].bitscore:
            hits[qid] = hit
    return hits


def curation(hyd_type: str, domains: DomainEvidence) -> tuple[str, str]:
    """Domain check for the assigned hydrogenase type (HydDB Methods, CDD step)."""
    if hyd_type == "NiFe":
        if domains.nifese_hases:
            return CLEARED, "NiFeSe_Hases catalytic domain observed"
        if domains.complex1:
            return NEEDS_CURATION, (
                "Complex1_49kDa/30kDa-superfamily domain without NiFeSe_Hases; "
                "the superfamily includes respiratory Complex I subunits"
            )
        return NEEDS_CURATION, "no NiFe catalytic-domain Pfam observed"
    if hyd_type == "FeFe":
        if domains.fe_hyd:
            return CLEARED, "Fe_hyd_lg_C or Fe_hyd_SSU catalytic domain observed"
        return NEEDS_CURATION, "no FeFe catalytic-domain Pfam observed"
    if domains.hmd:
        return CLEARED, "HMD domain observed"
    return NEEDS_CURATION, "no HMD domain observed"


SearchFn = Callable[[dict[str, str], ReferenceInfo, int], dict[str, ReferenceHit]]


def classify(conn, reference: ReferenceInfo, threads: int = 4, search: SearchFn = run_diamond) -> list[Classification]:
    """One reconciled classification per protein with a HydDB HMM hit."""
    conn.execute("""
        CREATE OR REPLACE TEMP TABLE _hyd_proteins AS
        SELECT DISTINCT protein_id FROM annotations WHERE LOWER(source) = 'hyddb'
    """)
    discovery: dict[str, list[tuple[str | None, float | None]]] = {}
    for pid, accession, score in conn.execute("""
        SELECT protein_id, accession, score FROM annotations
        WHERE LOWER(source) = 'hyddb' ORDER BY protein_id
    """).fetchall():
        discovery.setdefault(pid, []).append((_discovery_class(accession), score))

    pfam: dict[str, set[str]] = {}
    for pid, accession, name in conn.execute("""
        SELECT a.protein_id, a.accession, a.name FROM annotations a
        JOIN _hyd_proteins h USING (protein_id) WHERE LOWER(a.source) = 'pfam'
    """).fetchall():
        pfam.setdefault(pid, set()).update(pfam_keys(accession, name))

    sequences = {
        pid: seq for pid, seq in conn.execute("""
            SELECT p.protein_id, p.sequence FROM proteins p JOIN _hyd_proteins h USING (protein_id)
            WHERE p.sequence IS NOT NULL AND LENGTH(p.sequence) > 0
        """).fetchall()
    }
    conn.execute("DROP TABLE IF EXISTS _hyd_proteins")

    hits = search(sequences, reference, threads)

    results = []
    for pid in sorted(discovery):
        rows = discovery[pid]
        classes = tuple(sorted({c for c, _ in rows if c}))
        scores = [s for _, s in rows if s is not None]
        item = Classification(
            protein_id=pid,
            outcome=NO_REFERENCE_HIT,
            discovery_classes=classes,
            discovery_hit_count=len(rows),
            discovery_best_score=max(scores) if scores else None,
            domains=DomainEvidence.from_identifiers(pfam.get(pid, ())),
        )
        hit = hits.get(pid)
        if pid not in sequences:
            item.outcome = MISSING_SEQUENCE
        elif hit is None:
            item.outcome = NO_REFERENCE_HIT
        elif hit.hyd_type is None:
            item.hit, item.outcome = hit, UNPARSED_REFERENCE_LABEL
        elif hit.hyd_type not in classes:
            item.hit, item.outcome = hit, CLASS_CONFLICT
            item.curation_status = NEEDS_CURATION
            item.curation_reason = (
                f"discovery HMM class {'/'.join(classes) or 'unrecognized'} "
                f"disagrees with reference class {hit.hyd_type}"
            )
        else:
            item.hit, item.outcome = hit, ASSIGNED
            item.subgroup = lookup(hit.hyd_type, hit.subgroup) or Subgroup(
                hit.hyd_type, hit.subgroup, "(not in reference table)",
                "Label present in the installed reference; interpretation unverified.", UNVERIFIED,
            )
            item.curation_status, item.curation_reason = curation(hit.hyd_type, item.domains)
        results.append(item)
    return results


CLASSIFICATION_COLUMNS = (
    "protein_id VARCHAR PRIMARY KEY",
    "outcome VARCHAR NOT NULL",
    "discovery_classes VARCHAR",
    "discovery_hit_count INTEGER",
    "discovery_best_score DOUBLE",
    "reference_accession VARCHAR",
    "reference_organism VARCHAR",
    "reference_label VARCHAR",
    "reference_class VARCHAR",
    "reference_subgroup VARCHAR",
    "interpretation_status VARCHAR",
    "reference_role VARCHAR",
    "pident DOUBLE",
    "evalue DOUBLE",
    "bitscore DOUBLE",
    "has_nifese_hases BOOLEAN",
    "has_fe_hyd BOOLEAN",
    "has_complex1 BOOLEAN",
    "has_hmd BOOLEAN",
    "curation_status VARCHAR",
    "curation_reason VARCHAR",
    "derived_labels VARCHAR[]",
    "reference_release VARCHAR",
    "reference_sha256 VARCHAR",
    "classifier_version VARCHAR",
)


def _row(item: Classification, reference: ReferenceInfo) -> tuple:
    hit, sub, dom = item.hit, item.subgroup, item.domains
    return (
        item.protein_id, item.outcome, ";".join(item.discovery_classes), item.discovery_hit_count,
        item.discovery_best_score,
        hit.reference_accession if hit else None, hit.organism if hit else None,
        hit.label if hit else None, hit.hyd_type if hit else None, hit.subgroup if hit else None,
        sub.status if sub else None, sub.reference_role if sub else None,
        hit.pident if hit else None, hit.evalue if hit else None, hit.bitscore if hit else None,
        dom.nifese_hases, dom.fe_hyd, dom.complex1, dom.hmd,
        item.curation_status, item.curation_reason or None, list(item.derived_labels),
        reference.release, reference.sha256, CLASSIFIER_VERSION,
    )


def _description(item: Classification, label: str, reference: ReferenceInfo) -> str:
    if label == CURATION_FLAG:
        return f"HydDB assignment needs curation: {item.curation_reason}"
    if label == CLASS_CONFLICT_FLAG:
        return f"HydDB class conflict: {item.curation_reason}"
    sub = item.subgroup
    return (
        f"Sharur nearest-reference match to HydDB {reference.release}: {sub.label} "
        f"({sub.status}; provisional)"
    )


def write_classifications(conn, results: list[Classification], reference: ReferenceInfo) -> dict[str, int]:
    """Replace ``hydrogenase_classifications`` and ``hyddb_subgroup`` rows in one transaction."""
    conn.execute("BEGIN TRANSACTION")
    try:
        conn.execute(f"CREATE TABLE IF NOT EXISTS hydrogenase_classifications ({', '.join(CLASSIFICATION_COLUMNS)})")
        conn.execute("DELETE FROM hydrogenase_classifications")
        if results:
            placeholders = ",".join(["?"] * len(CLASSIFICATION_COLUMNS))
            conn.executemany(
                f"INSERT INTO hydrogenase_classifications VALUES ({placeholders})",
                [_row(item, reference) for item in results],
            )

        conn.execute("DELETE FROM annotations WHERE source = 'hyddb_subgroup'")
        next_id = conn.execute("SELECT COALESCE(MAX(annotation_id), 0) FROM annotations").fetchone()[0]
        ann_rows = []
        for item in results:
            evalue = item.hit.evalue if item.hit else None
            bitscore = item.hit.bitscore if item.hit else None
            for label in item.derived_labels:
                next_id += 1
                ann_rows.append((next_id, item.protein_id, "hyddb_subgroup", label, label,
                                 _description(item, label, reference), evalue, bitscore))
        if ann_rows:
            conn.executemany(
                """INSERT INTO annotations
                   (annotation_id, protein_id, source, accession, name, description, evalue, score)
                   VALUES (?, ?, ?, ?, ?, ?, ?, ?)""",
                ann_rows,
            )
        conn.execute("COMMIT")
    except Exception:
        conn.execute("ROLLBACK")
        raise
    return {"classifications": len(results), "subgroup_annotations": len(ann_rows)}


def classify_database(db_path: str | Path, threads: int = 4, reference_dir: Path | None = None,
                      search: SearchFn = run_diamond) -> list[Classification]:
    """Classify and write in place. Use on a build or staging database only."""
    import duckdb

    reference = find_reference(reference_dir)
    conn = duckdb.connect(str(db_path))
    try:
        results = classify(conn, reference, threads=threads, search=search)
        write_classifications(conn, results, reference)
    finally:
        conn.close()
    return results


@dataclass
class Summary:
    outcomes: dict[str, int] = field(default_factory=dict)
    curation: dict[str, int] = field(default_factory=dict)


def summarize(results: list[Classification]) -> Summary:
    summary = Summary()
    for item in results:
        summary.outcomes[item.outcome] = summary.outcomes.get(item.outcome, 0) + 1
        if item.curation_status:
            summary.curation[item.curation_status] = summary.curation.get(item.curation_status, 0) + 1
    return summary
