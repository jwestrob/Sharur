"""MinCED CRISPR output: per-array repeat and spacer tables.

MinCED writes a text report beside its GFF. For every array it lists each
repeat's position with the repeat and the following spacer::

    Sequence 'contig_1' (490885 bp)

    CRISPR 1   Range: 13653 - 14585
    POSITION	REPEAT				SPACER
    --------	------------	------------
    13653		<repeat>	<spacer>	[ 28, 38 ]
    ...
    14558		<repeat>
    --------	------------	------------
    Repeats: 14	Average Length: 28		Average Length: 41

:func:`parse_minced_text` returns those arrays with 1-based repeat and spacer
coordinates, keyed by contig and array start so they join to the GFF calls.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

_SEQUENCE = re.compile(r"^Sequence '([^']+)'")
_ARRAY = re.compile(r"^CRISPR\s+(\d+)\s+Range:\s*(\d+)\s*-\s*(\d+)")
_ROW = re.compile(r"^(\d+)\s+([ACGTUNacgtun]+)(?:\s+([ACGTUNacgtun]+))?")


def parse_minced_text(path: Path) -> list[dict[str, Any]]:
    """Arrays from a MinCED text report, each with ``repeats`` and ``spacers``."""
    arrays: list[dict[str, Any]] = []
    contig = None
    current: dict[str, Any] | None = None
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            m = _SEQUENCE.match(line)
            if m:
                contig = m.group(1).split()[0]
                continue
            m = _ARRAY.match(line)
            if m:
                current = {"contig": contig, "number": int(m.group(1)), "start": int(m.group(2)),
                           "end": int(m.group(3)), "repeats": [], "spacers": []}
                arrays.append(current)
                continue
            if current is None:
                continue
            m = _ROW.match(line)
            if m:
                pos, repeat, spacer = int(m.group(1)), m.group(2).upper(), (m.group(3) or "").upper()
                current["repeats"].append({"start": pos, "end": pos + len(repeat) - 1, "seq": repeat})
                if spacer:
                    s = pos + len(repeat)
                    current["spacers"].append({"start": s, "end": s + len(spacer) - 1, "seq": spacer,
                                               "length": len(spacer)})
            elif line.startswith("Repeats:"):
                current = None
    return arrays


def consensus(repeats: list[str]) -> str:
    """Column-wise majority of equal-length repeats (the longest length class wins)."""
    if not repeats:
        return ""
    length = max(set(map(len, repeats)), key=lambda n: (sum(len(r) == n for r in repeats), n))
    same = [r for r in repeats if len(r) == length]
    return "".join(max(set(col), key=col.count) for col in zip(*same))


def annotate_repeats(array: dict[str, Any], reference: str | None = None) -> dict[str, Any]:
    """Mark each repeat's differences from ``reference`` (default: the array's own consensus)."""
    ref = (reference or consensus([r["seq"] for r in array["repeats"]])).upper()
    budget = max(2, len(ref) // 10)
    for r in array["repeats"]:
        diff = [k for k, (a, b) in enumerate(zip(r["seq"], ref)) if a != b]
        diff += list(range(min(len(r["seq"]), len(ref)), max(len(r["seq"]), len(ref))))
        r["diff"], r["mismatches"] = diff, len(diff)
        r["diverged"] = len(diff) > budget
    array["consensus"] = ref
    return array
