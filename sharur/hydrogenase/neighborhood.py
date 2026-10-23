"""Neighborhood evidence for hydrogenase calls (step 3 of ``.claude/skills/hydrogenase.md``).

[NiFe] Group 4 large subunits share the Complex1_49kDa superfamily with respiratory Complex I
subunit D, so a call without the NiFeSe_Hases catalytic domain needs gene context. Each call's
+/-6-gene neighborhood is read for markers of three kinds:

* hydrogenase-specific (``Context.hydrogenase``): evidence a Complex I-like module lacks. The
  NiFeSe_Hases / Fe_hyd catalytic domains; maturation genes (hypA-F, hypX, the HycI-family
  maturation proteases hyaD/hybD, hycI and frhD, HycH); and KEGG-named subunits of the
  [NiFe] Group 1-3, [FeFe] and [Fe] hydrogenases (``SPECIFIC_KOS``).
* Complex I-like hydrogenase subunits (``Context.complex_i_like``): subunits of the Complex
  I-related hydrogenase complexes (Hyc, Hyf, Ech, Eha, Ehb, Mbh and the FHL formate
  dehydrogenase) and the NuoE/F/G-type diaphorase subunits (``COMPLEX_I_LIKE_KOS``), plus the
  Ech/Eha/Ehb/Hyc/Hyf subunit Pfam families. Each is a Complex I homolog or sits in a complex
  built from them: hyfB/D/F proteins carry Proton_antipo_M (the NuoL/M/N family), hyfG
  captures Group 4a large subunits (NuoD), and the hox diaphorase subunits carry the NuoE/F/G
  families. A Complex I-like module can draw these KOs, so they distinguish nothing alone.
* Complex I (``Context.complex_i``): KOfam profiles of respiratory Complex I and its archaeal
  F420H2 relatives (``COMPLEX_I_KOS``): nuoA-N and the fused nuo KOs, ndhA-N, Fpo and Fqo
  without their FrhB-family F subunits, and the mitochondrial core subunits. Mrp/Mnh
  antiporter KOs count on neither side, because Mbh-type Group 4 hydrogenases carry Mrp modules.

Pfam Complex I-superfamily families (Oxidored_q6, Complex1_30kDa, ...) count on neither side:
hydrogenase complexes carry them too.

Verdicts:

* ``supported``: a hydrogenase-specific marker and no Complex I marker;
* ``complex_i_context``: a Complex I marker and no hydrogenase-specific marker;
* ``mixed``: both;
* ``complex_i_homolog_only``: Complex I-like hydrogenase subunits only, read as unresolved;
* ``no_markers``, ``no_position``.

Complex I-like markers set the verdict only when neither informative side is present. A verdict
is evidence about context: a named complex (Ech, Hyc, ...) still needs its subunit genes, and a
physiological role still needs the genome's metabolic context.
"""

from __future__ import annotations

import re
from collections import defaultdict
from dataclasses import dataclass, field

from sharur.predicates.pfam_identity import normalize_pfam_accession

WINDOW = 6

HYDROGENASE = "hydrogenase"
COMPLEX_I_LIKE = "complex_i_like"
COMPLEX_I = "complex_i"

# Every KO below carries the first KEGG symbol of its entry in the local KEGG snapshot
# (``kegg_evidence_snapshot.tsv``); tests/test_hydrogenase_neighborhood.py re-checks them.
COMPLEX_I_KOS: dict[str, str] = {
    # prokaryotic Complex I (M00144) and its fused subunits
    "K00330": "nuoA", "K00331": "nuoB", "K00332": "nuoC", "K00333": "nuoD", "K00334": "nuoE", "K00335": "nuoF",
    "K00336": "nuoG", "K00337": "nuoH", "K00338": "nuoI", "K00339": "nuoJ", "K00340": "nuoK", "K00341": "nuoL",
    "K00342": "nuoM", "K00343": "nuoN", "K13378": "nuoCD", "K13380": "nuoBCD", "K15863": "nuoLM",
    # NAD(P)H:quinone oxidoreductase, chloroplasts and cyanobacteria (M00145)
    "K05572": "ndhA", "K05573": "ndhB", "K05574": "ndhC", "K05575": "ndhD", "K05576": "ndhE", "K05577": "ndhF",
    "K05578": "ndhG", "K05579": "ndhH", "K05580": "ndhI", "K05581": "ndhJ", "K05582": "ndhK", "K05583": "ndhL",
    "K05584": "ndhM", "K05585": "ndhN",
    # F420H2 dehydrogenase (Fpo) and F420H2:quinone oxidoreductase (Fqo). fpoF (K22162) and fqoF (K22174)
    # stay out: they belong to the FrhB_FdhB family of the F420-reducing [NiFe] hydrogenase subunit FrhB.
    "K22158": "fpoA", "K22159": "fpoB", "K22160": "fpoC", "K22161": "fpoD", "K22163": "fpoH", "K22164": "fpoI",
    "K22165": "fpoJ", "K22166": "fpoK", "K22167": "fpoL", "K22168": "fpoM", "K22169": "fpoN", "K22170": "fpoO",
    "K22171": "fqoA", "K22172": "fqoBC", "K22173": "fqoD", "K22175": "fqoH", "K22176": "fqoI", "K22177": "fqoJ",
    "K22178": "fqoK", "K22179": "fqoL", "K22180": "fqoM", "K22181": "fqoN",
    # mitochondrial core subunits (M00142, M00143). In SR-VP bacteria every NDUFS2/3/7/8 and NDUFV2 hit
    # sits beside a nuo gene; the accessory NDUFS4/5/6 and NDUFV3 stay out (NDUFS4/6 hits: 0 of 62).
    "K03878": "ND1", "K03879": "ND2", "K03880": "ND3", "K03881": "ND4", "K03882": "ND4L", "K03883": "ND5",
    "K03884": "ND6", "K03934": "NDUFS1", "K03935": "NDUFS2", "K03936": "NDUFS3", "K03940": "NDUFS7",
    "K03941": "NDUFS8", "K03942": "NDUFV1", "K03943": "NDUFV2",
}

SPECIFIC_KOS: dict[str, str] = {
    # maturation: hyp genes, hypX, and the HycI-family maturation proteases (frhD proteins carry Pfam HycI)
    "K04651": "hypA", "K04652": "hypB", "K04653": "hypC", "K04654": "hypD", "K04655": "hypE", "K04656": "hypF",
    "K03605": "hyaD", "K08315": "hycI", "K00442": "frhD", "K19640": "hypX", "K03618": "hyaF", "K03619": "hyaE",
    # [NiFe] Groups 1-3: large, small, cytochrome b and electron-transfer subunits
    "K06281": "hyaB", "K06282": "hyaA", "K03620": "hyaC", "K05922": "hydB", "K05927": "hydA", "K00437": "hydB",
    "K18008": "hydA", "K14068": "vhoA", "K14069": "vhoC", "K14070": "vhoG", "K23549": "hupV", "K23548": "hupU",
    "K00436": "hoxH", "K18007": "hoxY", "K00440": "frhA", "K00441": "frhB", "K00443": "frhG", "K14126": "mvhA",
    "K14127": "mvhD", "K14128": "mvhG", "K17993": "hydA", "K17994": "hydD", "K17995": "hydG", "K17996": "hydB",
    # [FeFe] catalytic subunits and the [Fe] hydrogenase
    "K00532": "E1.12.7.2", "K00533": "E1.12.7.2L", "K00534": "E1.12.7.2S", "K17997": "hydA", "K18332": "hndD",
    "K25123": "hydA2", "K13942": "hmd",
}

COMPLEX_I_LIKE_KOS: dict[str, str] = {
    # Complex I-related (Group 4) hydrogenase complexes: Hyc (formate hydrogenlyase), Hyf, Ech, Eha, Ehb, Mbh
    "K15827": "hycB", "K15828": "hycC", "K15829": "hycD", "K15830": "hycE", "K15831": "hycF", "K15832": "hycG",
    "K12136": "hyfA", "K12137": "hyfB", "K12138": "hyfC", "K12139": "hyfD", "K12140": "hyfE", "K12141": "hyfF",
    "K12142": "hyfG", "K12143": "hyfH", "K12144": "hyfI", "K12145": "hyfJ", "K14086": "echA", "K14087": "echB",
    "K14088": "echC", "K14089": "echD", "K14090": "echE", "K14091": "echF", "K14092": "ehaA", "K14093": "ehaB",
    "K14094": "ehaC", "K14095": "ehaD", "K14096": "ehaE", "K14097": "ehaF", "K14098": "ehaG", "K14099": "ehaH",
    "K14100": "ehaI", "K14101": "ehaJ", "K14102": "ehaK", "K14103": "ehaL", "K14104": "ehaM", "K14105": "ehaN",
    "K14106": "ehaO", "K14107": "ehaP", "K14108": "ehaQ", "K14109": "ehaR", "K14110": "ehbA", "K14111": "ehbB",
    "K14112": "ehbC", "K14113": "ehbD", "K14114": "ehbE", "K14115": "ehbF", "K14116": "ehbG", "K14117": "ehbH",
    "K14118": "ehbI", "K14119": "ehbJ", "K14120": "ehbK", "K14121": "ehbL", "K14122": "ehbM", "K14123": "ehbN",
    "K14124": "ehbO", "K14125": "ehbP", "K06862": "ehbQ", "K18016": "mbhL", "K18017": "mbhK", "K18023": "mbhJ",
    # FHL formate dehydrogenase: Molybdopterin family, which nuoG proteins carry
    "K22015": "fdhF",
    # NuoE/F/G-type diaphorase subunits of bidirectional [NiFe] and multimeric [FeFe] hydrogenases
    # (hoxE/F/U carry the NuoE/F/G Pfam families; KEGG gives K05586-K05588 the Complex I EC 7.1.1.2;
    # in SR-VP bacteria 144 of 147 hndB proteins sit beside nuo genes)
    "K05586": "hoxE", "K05587": "hoxF", "K05588": "hoxU", "K18005": "hoxF", "K18006": "hoxU", "K18330": "hndA",
    "K18331": "hndC", "K17992": "hndB", "K17998": "hydB", "K17999": "hydC", "K06441": "E1.12.7.2G",
}

# KOs of ``kegg_hyddb_snapshot.tsv`` (KOs KEGG names as hydrogenases) left on neither side, with the reason.
UNUSED_KOS: dict[str, str] = {
    "K12146": "hyfR: transcriptional activator (KEGG BRITE ko03000)",
    "K15836": "fhlA: transcriptional activator (KEGG BRITE ko03000); its SR-VP hits are sigma-54 activators",
    "K07388": "AIRS-family 'hydrogenase expression/formation protein'; 2 of 89 SR-VP archaeal copies sit beside a hyp gene",
    "K15833": "absent from the current KEGG snapshot",
    "K15834": "absent from the current KEGG snapshot",
}

CATALYTIC_PFAM = {"NiFeSe_Hases": "PF00374", "Fe_hyd_lg_C": "PF02906", "Fe_hyd_SSU": "PF02256"}
MATURATION_PFAM = {"HycI": "PF01750", "HycH": "PF07450"}  # hydrogenase maturation protease; FHL maturation protein
SPECIFIC_PFAM = CATALYTIC_PFAM | MATURATION_PFAM
# Pfam subunit families of the Ech/Eha/Ehb/Hyc/Hyf complexes (EhaB, EhaL, ...). Case-sensitive, one subunit
# letter: Pfam's ECH_1/ECH_2 (enoyl-CoA hydratase), CooT and CooC_C carry no hydrogenase subunit.
COMPLEX_I_LIKE_PFAM_NAME = r"^(Ech|Eha|Ehb|Hyc|Hyf)[A-Z]$"
COMPLEX_I_LIKE_PFAM = {"Ehbp": "PF10622", "NiFe_hyd_3_EhaA": "PF17367"}
_COMPLEX_I_LIKE_NAME = re.compile(COMPLEX_I_LIKE_PFAM_NAME)
_SPECIFIC_PFAM_KEYS = frozenset(SPECIFIC_PFAM) | frozenset(SPECIFIC_PFAM.values())
_COMPLEX_I_LIKE_PFAM_KEYS = frozenset(COMPLEX_I_LIKE_PFAM) | frozenset(COMPLEX_I_LIKE_PFAM.values())

VERDICTS = ("supported", "complex_i_context", "mixed", "complex_i_homolog_only", "no_markers", "no_position")
VERDICT_LABELS = {
    "supported": "hydrogenase-specific genes nearby",
    "complex_i_context": "Complex I genes nearby",
    "mixed": "both nearby",
    "complex_i_homolog_only": "only Complex I-like hydrogenase subunits nearby",
    "no_markers": "no marker genes nearby",
    "no_position": "no genomic position",
}


@dataclass
class Context:
    verdict: str
    hydrogenase: set[str] = field(default_factory=set)        # "accession name" of each hydrogenase-specific marker
    complex_i_like: set[str] = field(default_factory=set)     # Complex I-like hydrogenase subunits
    complex_i: set[str] = field(default_factory=set)
    hydrogenase_genes: set[str] = field(default_factory=set)  # neighbor protein IDs carrying them
    complex_i_like_genes: set[str] = field(default_factory=set)
    complex_i_genes: set[str] = field(default_factory=set)


def marker_kind(accession: str, name: str) -> str | None:
    """``hydrogenase``, ``complex_i_like``, ``complex_i`` or None for one annotation row."""
    acc = normalize_pfam_accession(accession)
    if acc in COMPLEX_I_KOS:
        return COMPLEX_I
    if acc in SPECIFIC_KOS or acc in _SPECIFIC_PFAM_KEYS or name in _SPECIFIC_PFAM_KEYS:
        return HYDROGENASE
    if (acc in COMPLEX_I_LIKE_KOS or acc in _COMPLEX_I_LIKE_PFAM_KEYS or name in _COMPLEX_I_LIKE_PFAM_KEYS
            or _COMPLEX_I_LIKE_NAME.match(name) or _COMPLEX_I_LIKE_NAME.match(acc)):
        return COMPLEX_I_LIKE
    return None


_KO_SYMBOLS = COMPLEX_I_KOS | SPECIFIC_KOS | COMPLEX_I_LIKE_KOS


def _label(accession: str, name: str) -> str:
    """``K12137 hyfB`` for a KO, ``NiFeSe_Hases`` (or ``PF00374 NiFeSe_Hases``) for a Pfam hit."""
    symbol = _KO_SYMBOLS.get(accession)
    if symbol:
        return f"{accession} {symbol}"
    return f"{accession} {name}" if name and name != accession else (accession or name)


def verdict_of(ctx: Context) -> str:
    """Verdict from the marker sets; Complex I-like subunits decide only when nothing informative is nearby."""
    if ctx.hydrogenase and ctx.complex_i:
        return "mixed"
    if ctx.hydrogenase:
        return "supported"
    if ctx.complex_i:
        return "complex_i_context"
    if ctx.complex_i_like:
        return "complex_i_homolog_only"
    return "no_markers"


def neighborhood_contexts(store, protein_ids, window: int = WINDOW) -> dict[str, Context]:
    """Verdict per protein from its +/-``window``-gene neighborhood on the same contig."""
    ids = list(dict.fromkeys(protein_ids))
    if not ids:
        return {}
    placed = {pid for (pid,) in store.execute(
        "SELECT protein_id FROM proteins WHERE protein_id IN (SELECT UNNEST(?::VARCHAR[])) AND contig_id <> protein_id",
        [ids])}
    keys = sorted(set(COMPLEX_I_KOS) | set(SPECIFIC_KOS) | set(COMPLEX_I_LIKE_KOS)
                  | _SPECIFIC_PFAM_KEYS | _COMPLEX_I_LIKE_PFAM_KEYS)
    rows = store.execute(
        """WITH f AS (SELECT protein_id, contig_id, gene_index FROM proteins
                      WHERE protein_id IN (SELECT UNNEST(?::VARCHAR[])) AND contig_id <> protein_id),
                nb AS (SELECT f.protein_id AS focal, q.protein_id FROM f JOIN proteins q
                       ON q.contig_id = f.contig_id AND q.gene_index BETWEEN f.gene_index - ? AND f.gene_index + ?
                       AND q.protein_id <> f.protein_id)
           SELECT nb.focal, nb.protein_id, a.accession, COALESCE(a.name, '') FROM nb
           JOIN annotations a ON a.protein_id = nb.protein_id
           WHERE regexp_replace(a.accession, '\\.[0-9]+$', '') IN (SELECT UNNEST(?::VARCHAR[]))
              OR a.name IN (SELECT UNNEST(?::VARCHAR[]))
              OR regexp_matches(COALESCE(a.name, ''), ?) OR regexp_matches(COALESCE(a.accession, ''), ?)""",
        [ids, window, window, keys, keys, COMPLEX_I_LIKE_PFAM_NAME, COMPLEX_I_LIKE_PFAM_NAME])
    found: dict[str, Context] = defaultdict(lambda: Context("no_markers"))
    for focal, neighbor, accession, name in rows:
        kind = marker_kind(accession or "", name or "")
        if kind is None:
            continue
        ctx = found[focal]
        label = _label(accession or "", name or "")
        if kind == HYDROGENASE:
            ctx.hydrogenase.add(label)
            ctx.hydrogenase_genes.add(neighbor)
        elif kind == COMPLEX_I_LIKE:
            ctx.complex_i_like.add(label)
            ctx.complex_i_like_genes.add(neighbor)
        else:
            ctx.complex_i.add(label)
            ctx.complex_i_genes.add(neighbor)
    out: dict[str, Context] = {}
    for pid in ids:
        if pid not in placed:
            out[pid] = Context("no_position")
            continue
        ctx = found.get(pid) or Context("no_markers")
        ctx.verdict = verdict_of(ctx)
        out[pid] = ctx
    return out
