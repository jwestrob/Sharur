"""Neighborhood evidence for hydrogenase calls (step 3 of ``.claude/skills/hydrogenase.md``).

[NiFe] Group 4 large subunits share the Complex1_49kDa superfamily with respiratory Complex I
subunit D, so a call without the NiFeSe_Hases catalytic domain needs gene context. Each call's
+/-6-gene neighborhood is read for markers on two sides:

* hydrogenase: KOs KEGG names as hydrogenases (the KEGG-HydDB snapshot, which includes the hyp
  maturation genes), the Hyc/Hyf/Ech/maturation KOs, the NiFeSe_Hases / Fe_hyd catalytic domains,
  and Ech/Eha/Ehb/Hyc/Hyf/Coo gene names;
* Complex I: KOfam nuoA-N. Pfam families such as Oxidored_q6 or Complex1_30kDa are left out:
  hydrogenase complexes carry them too (Oxidored_q6 covers [NiFe] small subunits), and on DPANN
  they sit beside 699 of 1,258 calls that carry the catalytic domain.

Verdicts: ``supported`` (hydrogenase markers only), ``complex_i_context`` (Complex I markers only),
``mixed`` (both), ``no_markers`` (neither), ``no_position`` (no genomic position). A verdict is
evidence about context: a named complex (Ech, Hyc, ...) still needs its subunit genes, and a
physiological role still needs the genome's metabolic context.
"""

from __future__ import annotations

import re
from collections import defaultdict
from dataclasses import dataclass, field

WINDOW = 6
COMPLEX_I_KOS = frozenset(f"K{n:05d}" for n in range(330, 344))            # nuoA-N
HYC_KOS = frozenset(f"K{n:05d}" for n in range(15827, 15834))              # hycB-G, hycA (formate hydrogenlyase)
HYF_KOS = frozenset(f"K{n:05d}" for n in range(12136, 12146))              # hyfA-J (hydrogenase-4)
ECH_KOS = frozenset(f"K{n:05d}" for n in range(14086, 14092))              # echA-F
MATURATION_KOS = frozenset({"K04651", "K04652", "K04653", "K04654", "K04655", "K04656", "K03605"})  # hypA-F, hyaD/hybD
HYDROGENASE_PFAM = frozenset({"NiFeSe_Hases", "Fe_hyd_lg_C", "Fe_hyd_SSU"})
_COMPLEX_NAME = re.compile(r"^(ech|eha|ehb|hyc|hyf|coo)", re.I)

VERDICTS = ("supported", "complex_i_context", "mixed", "no_markers", "no_position")
VERDICT_LABELS = {
    "supported": "hydrogenase genes nearby",
    "complex_i_context": "Complex I genes nearby",
    "mixed": "both nearby",
    "no_markers": "no marker genes nearby",
    "no_position": "no genomic position",
}


def hydrogenase_kos() -> frozenset[str]:
    """KOs read as hydrogenase-side evidence: KEGG's hydrogenase KOs plus the skill's named sets."""
    from sharur.hydrogenase.ko_association import named_kos  # noqa: PLC0415

    return named_kos() | HYC_KOS | HYF_KOS | ECH_KOS | MATURATION_KOS


@dataclass
class Context:
    verdict: str
    hydrogenase: set[str] = field(default_factory=set)        # "accession name" of each hydrogenase-side marker
    complex_i: set[str] = field(default_factory=set)
    hydrogenase_genes: set[str] = field(default_factory=set)  # neighbor protein IDs carrying them
    complex_i_genes: set[str] = field(default_factory=set)


def _marker(accession: str, name: str, hyd_kos: frozenset[str]) -> str | None:
    if accession in COMPLEX_I_KOS:
        return "complex_i"
    if accession in hyd_kos or {accession, name} & HYDROGENASE_PFAM or _COMPLEX_NAME.match(name):
        return "hydrogenase"
    return None


def neighborhood_contexts(store, protein_ids, window: int = WINDOW) -> dict[str, Context]:
    """Verdict per protein from its +/-``window``-gene neighborhood on the same contig."""
    ids = list(dict.fromkeys(protein_ids))
    if not ids:
        return {}
    hyd_kos = hydrogenase_kos()
    placed = {pid for (pid,) in store.execute(
        "SELECT protein_id FROM proteins WHERE protein_id IN (SELECT UNNEST(?::VARCHAR[])) AND contig_id <> protein_id",
        [ids])}
    rows = store.execute(
        """WITH f AS (SELECT protein_id, contig_id, gene_index FROM proteins
                      WHERE protein_id IN (SELECT UNNEST(?::VARCHAR[])) AND contig_id <> protein_id),
                nb AS (SELECT f.protein_id AS focal, q.protein_id FROM f JOIN proteins q
                       ON q.contig_id = f.contig_id AND q.gene_index BETWEEN f.gene_index - ? AND f.gene_index + ?
                       AND q.protein_id <> f.protein_id)
           SELECT nb.focal, nb.protein_id, a.accession, COALESCE(a.name, '') FROM nb
           JOIN annotations a ON a.protein_id = nb.protein_id
           WHERE a.accession IN (SELECT UNNEST(?::VARCHAR[])) OR a.name IN (SELECT UNNEST(?::VARCHAR[]))
              OR regexp_matches(LOWER(COALESCE(a.name, '')), '^(ech|eha|ehb|hyc|hyf|coo)')""",
        [ids, window, window, sorted(COMPLEX_I_KOS | hyd_kos), sorted(HYDROGENASE_PFAM)])
    found: dict[str, Context] = defaultdict(lambda: Context("no_markers"))
    for focal, neighbor, accession, name in rows:
        side = _marker(accession or "", name or "", hyd_kos)
        if side is None:
            continue
        ctx = found[focal]
        if side == "hydrogenase":
            ctx.hydrogenase.add(f"{accession} {name}".strip())
            ctx.hydrogenase_genes.add(neighbor)
        else:
            ctx.complex_i.add(f"{accession} {name}".strip())
            ctx.complex_i_genes.add(neighbor)
    out: dict[str, Context] = {}
    for pid in ids:
        if pid not in placed:
            out[pid] = Context("no_position")
            continue
        ctx = found.get(pid) or Context("no_markers")
        if ctx.hydrogenase and ctx.complex_i:
            ctx.verdict = "mixed"
        elif ctx.hydrogenase:
            ctx.verdict = "supported"
        elif ctx.complex_i:
            ctx.verdict = "complex_i_context"
        else:
            ctx.verdict = "no_markers"
        out[pid] = ctx
    return out
