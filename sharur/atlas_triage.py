"""Reduce an Atlas candidate set to a short, ranked digest.

Groups occurrences into accession-set groups, drops census-tier systems,
flags known annotation confounds, and ranks by prevalence. Deterministic and
read-only; makes no model calls.
"""

from __future__ import annotations

import json
import re
import sqlite3
from collections import defaultdict
from pathlib import Path
from typing import Any

# A system present in at least this fraction of genomes is treated as census.
CENSUS_FRACTION = 0.60

# (label, trigger, discriminator): a group matching `trigger` but lacking
# `discriminator` asserts more than its accessions support. Keyed on what is
# observed, so it fires without anyone suspecting the trap in advance.
CONFOUNDS: list[tuple[str, str, str]] = [
    ("hydrogenase-vs-ComplexI", r"nife_group4|energy_conserving|ech_hydrogenase", r"NiFeSe_Hases"),
    ("RuBisCO-without-PRK", r"RuBisCO_large|rbcL", r"\bPRK\b|K00855"),
    ("CODH-vs-xanthine-DH", r"Ald_Xan_dh|XdhC", r"K03520|CoxL|coxM|K03519"),
    ("T2SS-vs-Tad", r"T2SS|gspD|secretin", r"GspD|Secretin|K02453"),
    ("WoodLjungdahl-methyl-branch", r"FTHFS|MTHFR", r"AcsB|CdhC|K14138"),
]

# (system, marker accessions, required accessions, corrected label): when a
# group carries the markers but none of the required set, the scanner's label
# is not supported and the group is counted under the corrected one instead.
RELABEL_RULES: list[tuple[str, set[str], set[str], str]] = [
    (
        "natural-competence",
        {"Tad", "TadE", "TadB", "TadC", "K12513", "RcpC", "Flp"},
        {"Competence", "K02237", "K02238", "HHH_3", "DUF4131", "SLBB"},
        "pilus-assembly(relabelled)",
    ),
]


def _loads(raw: Any) -> dict:
    if not raw:
        return {}
    try:
        value = json.loads(raw)
    except (TypeError, ValueError):
        return {}
    return value if isinstance(value, dict) else {}


def load_occurrences(db_path: Path) -> list[dict]:
    con = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    try:
        rows = con.execute(
            "SELECT signature, genome_id, candidate_type, reduction_features, evidence "
            "FROM candidate_occurrences"
        ).fetchall()
    finally:
        con.close()

    out = []
    for sig_raw, genome, ctype, rf_raw, ev in rows:
        sig, rf = _loads(sig_raw), _loads(rf_raw)
        out.append(
            {
                "system": sig.get("system", "?"),
                "genome": genome,
                "candidate_type": (ctype or "").replace("atlas-scan:", ""),
                "accessions": tuple(rf.get("accessions") or []),
                "subtype": rf.get("subtype") or "",
                "evidence": ev or "",
            }
        )
    return out


def flag_confounds(accessions: set[str]) -> list[str]:
    blob = " ".join(sorted(accessions))
    return [
        label
        for label, trigger, discriminator in CONFOUNDS
        if re.search(trigger, blob, re.I) and not re.search(discriminator, blob, re.I)
    ]


def relabel(system: str, accessions: set[str]) -> str:
    for target, markers, required, corrected in RELABEL_RULES:
        if system == target and (accessions & markers) and not (accessions & required):
            return corrected
    return system


def build_digest(
    occurrences: list[dict],
    *,
    min_genomes: int = 2,
    system: str | None = None,
) -> dict:
    if not occurrences:
        return {"total": 0, "n_genomes": 0, "census": [], "entries": []}

    n_genomes = len({o["genome"] for o in occurrences})
    by_system: dict[str, set[str]] = defaultdict(set)
    for o in occurrences:
        by_system[o["system"]].add(o["genome"])
    census = {s for s, g in by_system.items() if len(g) / n_genomes >= CENSUS_FRACTION}

    groups: dict[tuple[str, tuple], list[dict]] = defaultdict(list)
    for o in occurrences:
        label = relabel(o["system"], set(o["accessions"]))
        groups[(label, o["accessions"])].append(o)

    entries = []
    for (label, accessions), members in groups.items():
        genome_set = {m["genome"] for m in members}
        if len(genome_set) < min_genomes or (system and label != system):
            continue
        entries.append(
            {
                "system": label,
                "subtype": members[0]["subtype"],
                "accessions": accessions,
                "occurrences": len(members),
                "genome_set": genome_set,
                "confounds": flag_confounds(set(accessions)),
                "members": members,
                "is_census": label in census,
            }
        )

    for e in entries:
        e["genomes"] = len(e["genome_set"])
        e["prevalence"] = e["genomes"] / n_genomes
        # Unusual but recurrent: peaks mid-range, damped hard for census systems.
        e["interest"] = e["genomes"] * (1.0 - e["prevalence"]) ** 2
        if e["is_census"]:
            e["interest"] *= 0.05
    entries.sort(key=lambda e: -e["interest"])

    census_occ = sum(1 for o in occurrences if o["system"] in census)
    return {
        "total": len(occurrences),
        "n_genomes": n_genomes,
        "n_groups": len(groups),
        "census": sorted(census),
        "census_occurrences": census_occ,
        "entries": entries,
    }


def render(digest: dict, *, top: int = 25, evidence: int = 0) -> str:
    if not digest["total"]:
        return "no candidates"
    lines = [
        "",
        f"{digest['total']:,} occurrences / {digest['n_genomes']} genomes / "
        f"{digest['n_groups']:,} accession-set groups",
        f"census-tier systems (>={CENSUS_FRACTION:.0%} of genomes): {digest['census']}",
        f"  -> {digest['census_occurrences']:,} occurrences "
        f"({digest['census_occurrences'] / digest['total']:.0%}) are inventory, not findings",
        "",
        f"{'system':<26}{'sub':<14}{'occ':>5}{'gen':>5}{'prev':>7}  {'flags':<26}accessions",
    ]
    for e in digest["entries"][:top]:
        acc = ",".join(e["accessions"][:6]) + ("…" if len(e["accessions"]) > 6 else "")
        lines.append(
            f"{e['system'][:25]:<26}{e['subtype'][:13]:<14}{e['occurrences']:>5}"
            f"{e['genomes']:>5}{e['prevalence']:>7.1%}  "
            f"{','.join(e['confounds'])[:24]:<26}{acc[:70]}"
        )
        for m in e["members"][:evidence]:
            ev = _loads(m["evidence"])
            text = ev.get("hypothesis") or ev.get("interpretation") or ev.get("claim") or ""
            if text:
                lines.append(f"      · {text[:150]}")

    flagged = [e for e in digest["entries"] if e["confounds"]]
    if flagged:
        lines += [
            "",
            f"{len(flagged)} groups trip a confound check and must not be read as named calls:",
        ]
        lines += [
            f"   {e['system']:<26}{','.join(e['confounds'])}" for e in flagged[:10]
        ]
    return "\n".join(lines) + "\n"
