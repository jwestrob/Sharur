"""KO -> HydDB subgroup associations, used as supporting evidence.

``kegg_hyddb_snapshot.tsv`` (``scripts/build_kegg_hyddb_snapshot.py``) records,
for each KO that KEGG names as a hydrogenase, the HydDB labels of the reference
hydrogenases that score at or above the KO's KOfam threshold. That is an
association between a KOfam profile and HydDB reference labels: it describes
which HydDB references a KO captures, and it supports or questions a
nearest-reference subgroup assignment. It establishes no subgroup by itself,
and it leaves the assignment unchanged.

Support status for an assignment, per associated KO the protein hits (the best
status across KOs is reported; ``conflict`` requires every associated KO to
conflict):

- ``subgroup``: at least ``SUBGROUP_AGREEMENT`` of the KO's references carry the
  assigned label;
- ``compatible``: some of the KO's references carry the assigned label;
- ``group``: none carry the label, some share its HydDB group;
- ``conflict``: none share its group or type;
- ``none``: the protein hits no KO with at least ``MIN_REFERENCES`` references.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from pathlib import Path
from typing import TYPE_CHECKING

from sharur.hydrogenase.subgroups import group_of, parse_label


if TYPE_CHECKING:
    from collections.abc import Iterable


SNAPSHOT = Path(__file__).resolve().parents[1] / "predicates/mappings/data/kegg_hyddb_snapshot.tsv"
SUBGROUP_AGREEMENT = 0.8
MIN_REFERENCES = 5

SUBGROUP = "subgroup"
COMPATIBLE = "compatible"
GROUP = "group"
CONFLICT = "conflict"
NONE = "none"
_RANK = {SUBGROUP: 0, COMPATIBLE: 1, GROUP: 2, CONFLICT: 3}


@dataclass(frozen=True)
class Association:
    ko: str
    label: str
    hyd_type: str | None
    subgroup: str | None
    count: int
    total: int

    @property
    def fraction(self) -> float:
        return self.count / self.total


@dataclass(frozen=True)
class KOSupport:
    status: str
    detail: str = ""


@cache
def load_associations(path: Path = SNAPSHOT) -> dict[str, tuple[Association, ...]]:
    """KO -> its HydDB label associations, most frequent first."""
    associations: dict[str, tuple[Association, ...]] = {}
    with open(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            ko, total, labels = line.rstrip("\n").split("\t")
            rows = []
            for item in filter(None, labels.split("|")):
                label, count = item.rsplit("=", 1)
                rows.append(Association(ko, label, *parse_label(label), int(count), int(total)))
            if rows:
                associations[ko] = tuple(sorted(rows, key=lambda a: (-a.count, a.label)))
    return associations


def associations(ko: str) -> tuple[Association, ...]:
    """HydDB label associations of one KO (empty when it captures no references)."""
    return load_associations().get(ko, ())


def describe(ko: str) -> str:
    """``K15830 [NiFe]_Group_4a 47/47``, listing every associated label."""
    rows = associations(ko)
    if not rows:
        return ""
    return f"{ko} " + ", ".join(f"{a.label} {a.count}/{a.total}" for a in rows)


def _status(rows: tuple[Association, ...], hyd_type: str, subgroup: str) -> str:
    label_fraction = sum(a.fraction for a in rows if (a.hyd_type, a.subgroup) == (hyd_type, subgroup))
    if label_fraction >= SUBGROUP_AGREEMENT:
        return SUBGROUP
    if label_fraction > 0:
        return COMPATIBLE
    group = group_of(hyd_type, subgroup)
    if any(a.hyd_type == hyd_type and a.subgroup and group_of(a.hyd_type, a.subgroup) == group for a in rows):
        return GROUP
    return CONFLICT


def ko_support(kos: Iterable[str], hyd_type: str | None, subgroup: str | None) -> KOSupport:
    """How the KOs a protein hits relate to its assigned HydDB subgroup."""
    if not hyd_type or not subgroup:
        return KOSupport(NONE)
    graded = []
    for ko in sorted(set(kos)):
        rows = associations(ko)
        if rows and rows[0].total >= MIN_REFERENCES:
            graded.append((_status(rows, hyd_type, subgroup), describe(ko)))
    if not graded:
        return KOSupport(NONE)
    best = min(graded, key=lambda g: _RANK[g[0]])[0]
    return KOSupport(best, "; ".join(detail for _, detail in graded))
