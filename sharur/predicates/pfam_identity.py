"""Pfam domain identity: match a hit by accession or by profile name.

Stage 07 stores the Pfam accession (``PF00374``) when its reference map knows
the profile, and keeps the profile name (``NiFeSe_Hases``) in the accession
column otherwise, e.g. when the reference TSV is absent or predates the HMM
database. Code that looks up a specific domain therefore checks both forms.
"""

from __future__ import annotations

import re


_VERSIONED_ACCESSION = re.compile(r"^(PF\d{5})\.\d+$")

# Supporting domains for hydrogenase corroboration, as {accession, profile name}.
NIFESE_HASES = frozenset({"PF00374", "NiFeSe_Hases"})
COMPLEX1_49KDA = frozenset({"PF00346", "Complex1_49kDa"})
COMPLEX1_30KDA = frozenset({"PF00329", "Complex1_30kDa"})
FE_HYD_LG_C = frozenset({"PF02906", "Fe_hyd_lg_C"})
FE_HYD_SSU = frozenset({"PF02256", "Fe_hyd_SSU"})
COMPLEX1 = COMPLEX1_49KDA | COMPLEX1_30KDA
FE_HYD = FE_HYD_LG_C | FE_HYD_SSU


def normalize_pfam_accession(value: str | None) -> str:
    """Strip whitespace and a Pfam version suffix (``PF00374.26`` -> ``PF00374``)."""
    text = (value or "").strip()
    match = _VERSIONED_ACCESSION.match(text)
    return match.group(1) if match else text


def pfam_keys(accession: str | None, name: str | None = None) -> tuple[str, ...]:
    """Lookup keys for one Pfam hit: normalized accession first, then profile name."""
    keys = []
    for key in (normalize_pfam_accession(accession), (name or "").strip()):
        if key and key not in keys:
            keys.append(key)
    return tuple(keys)


def has_pfam_domain(identifiers, domain: frozenset[str]) -> bool:
    """True when any identifier (accession or profile name) names ``domain``."""
    return any(normalize_pfam_accession(value) in domain for value in identifiers)
