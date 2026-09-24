"""Evidence rules tying Pfam-derived predicates to Pfam's own information.

A (Pfam family, predicate) pair is supported by one of:

- ``go:<GO id>``: the family's InterPro GO annotation (pfam2go), closed over
  is_a/part_of ancestors, contains a GO anchor of the predicate;
- ``text:<match>``: the family's Pfam name or description states the predicate
  (see :mod:`sharur.predicates.mappings.pfam_evidence_spec` for conventions);
- ``enzyme:<EC> (<name>)``: the description names an enzyme listed in the
  Expasy ENZYME database whose EC class maps to the predicate through
  :data:`sharur.predicates.mappings.kegg_map.EC_TO_PREDICATES`;
- ``swissprot:k/n single-domain k/n`` or ``swissprot:k/n co-domains excluded``:
  reviewed UniProtKB/Swiss-Prot proteins containing the family agree on the
  predicate through curator-assigned EC numbers, experimentally supported GO
  terms, or the KEGG-evidenced predicates of their KOs (thresholds below;
  :mod:`sharur.predicates.mappings.swissprot_evidence`).

``scripts/build_pfam_predicate_map.py`` applies these rules to build the shipped
map; ``tests/test_pfam_map_integrity.py`` re-verifies every shipped pair.
"""

from __future__ import annotations

import re
from functools import cache

from sharur.predicates.mappings.pfam_evidence_spec import HOMONYMS
from sharur.predicates.mappings.pfam_evidence_spec import E as EVIDENCE


# Swiss-Prot consensus. Coverage: the predicate holds for >= SWISSPROT_COVERAGE[1]
# of >= SWISSPROT_COVERAGE[0] reviewed carriers of the family, with a Wilson 95%
# lower bound >= SWISSPROT_MIN_LOWER_BOUND. Attribution: single-domain carriers
# (when >= SWISSPROT_SINGLE[0] exist) agree at >= SWISSPROT_SINGLE[1]; otherwise
# no family that carries the predicate on its own evidence co-occurs in
# >= SWISSPROT_CODOMAIN_FRACTION of the supporting proteins.
SWISSPROT_COVERAGE = (5, 0.8)
SWISSPROT_MIN_LOWER_BOUND = 0.5
SWISSPROT_SINGLE = (3, 0.8)
SWISSPROT_CODOMAIN_FRACTION = 0.9


def wilson_lower_bound(k: int, n: int, z: float = 1.96) -> float:
    """Lower bound of the Wilson score interval for k successes in n trials."""
    p = k / n
    centre = p + z * z / (2 * n)
    spread = z * (p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5
    return (centre - spread) / (1 + z * z / n)


# Fold/architecture and bookkeeping predicates: function evidence never establishes them.
SWISSPROT_EXCLUDED = frozenset({
    "repeat_domain", "tpr_repeat", "wd40_repeat", "lrr_repeat", "kelch_repeat", "heat_repeat",
    "ankyrin_repeat", "sel1_repeat", "helix_turn_helix", "winged_helix", "ribbon_helix_helix",
    "helix_loop_helix", "zinc_finger", "coiled_coil", "beta_barrel", "beta_helix", "alpha_helical",
    "p_loop", "aaa_domain", "cbs_domain", "pin_domain", "binding", "hypothetical", "unannotated",
})

# Predicate names used by earlier mapping layers that are not in the
# vocabulary: vocabulary equivalent, or None when no equivalent exists.
ALIASES: dict[str, str | None] = {
 "unknown_function": "hypothetical", "iron_sulfur_cluster": "iron_sulfur",
 "aminoacyl_trna_synthetase": "trna_synthetase", "methyl_coenzyme_m_reductase": "mcr_complex",
 "polysaccharide_biosynthesis": "exopolysaccharide", "aldolase": "lyase",
 "fatty_acid_metabolism": "lipid_metabolism", "carbohydrate_metabolism": "sugar_metabolism",
 "dna_replication": "replication", "aaa_atpase": "aaa_domain", "hth_regulator": "helix_turn_helix",
 "quinone_biosynthesis": "cofactor_biosynthesis", "virulence": None, "hydrogenase_group4": "nife_group4",
 "hydrogenase_group3": "nife_group3", "archaellum": "flagellum", "zinc_ribbon": "zinc_binding",
 "cation_transport": "ion_transporter", "detoxification": None, "dna_ligase": "ligase_dna",
 "release_factor": "translation_factor", "tautomerase": "tautomerase", "elongation_factor": "translation_factor",
 "ribosome_biogenesis": None, "binding_protein": "abc_substrate_binding", "ribonuclease": "rnase",
 "ribonucleotide_reductase": "nucleotide_metabolism", "ferredoxin_dependent": "ferredoxin",
 "pyruvate_metabolism": "central_metabolism", "escrt_system": "vesicle_trafficking",
 "peptide_transport": "peptide_transporter", "acyltransferase": "transferase", "coa_binding": "coenzyme_a_binding",
 "flavoprotein": "flavin_binding", "signal_transduction": "signaling", "glyoxalase": None,
 "rna_degradation": "rnase", "nad_metabolism": "nad_biosynthesis", "cell_wall_biosynthesis": "cell_wall",
 "cyclase": "lyase", "polyketide_biosynthesis": "polyketide_synthesis", "selenium": "selenium_metabolism",
 "phosphopeptide_binding": "protein_binding", "regulatory_domain": "regulatory", "sugar_biosynthesis": "sugar_metabolism",
 "steroid_metabolism": "lipid_metabolism", "receptor": None, "polysaccharide_deacetylase": "carbohydrate_esterase",
 "aminopeptidase": "peptidase", "ribosome": "translation", "menaquinone_biosynthesis": "cofactor_biosynthesis",
 "pyridine_nucleotide": "nad_binding", "chromate_resistance": "heavy_metal_resistance", "glutathione": None,
 "sugar_nucleotide_metabolism": "sugar_metabolism", "ras_family": "gtpase", "motility": None,
 "quinone_metabolism": "cofactor_biosynthesis", "ammonia_oxidation": "nitrification", "biotin_metabolism": "cofactor_biosynthesis",
 "lysine_biosynthesis": "amino_acid_biosynthesis", "pas_domain": "signaling", "sensor": "signaling",
 "pseudouridine_synthase": "rna_processing", "helix_hairpin_helix": "dna_binding", "proteasome": None,
 "protein_degradation": "protease", "efflux": "efflux_pump", "amp_binding": None, "arginine_biosynthesis": "amino_acid_biosynthesis",
 "phosphoesterase": "phosphodiesterase", "pentose_phosphate_pathway": "pentose_phosphate", "pyrimidine_biosynthesis": "pyrimidine_metabolism",
 "polyhydroxybutyrate": None, "protein_insertion": "translocase", "initiation_factor": "translation_factor",
 "lipopolysaccharide": "lps_biosynthesis", "serine_protease": "protease", "dna_gyrase_inhibitor": None,
 "trna_binding": "rna_binding", "metal_resistance": "heavy_metal_resistance", "amidohydrolase": "amidase",
 "arginine_metabolism": "amino_acid_metabolism", "purine_biosynthesis": "purine_metabolism",
 "metabolic_detoxification": None, "redox": "oxidoreductase", "nadh_dehydrogenase": "respiration",
 "tubulin": "ftsz", "molybdenum_cofactor": "molybdenum_binding", "trna_processing": "rna_processing",
}


def split_alternatives(pattern):
    """Split a regex on top-level '|' (outside groups and classes)."""
    parts, depth, cls, buf, i = [], 0, False, "", 0
    while i < len(pattern):
        c = pattern[i]
        if c == "\\" and i + 1 < len(pattern):
            buf += pattern[i:i + 2]; i += 2; continue
        if cls:
            cls = c != "]"
        elif c == "[":
            cls = True
        elif c == "(":
            depth += 1
        elif c == ")":
            depth -= 1
        elif c == "|" and depth == 0:
            parts.append(buf); buf = ""; i += 1; continue
        buf += c; i += 1
    parts.append(buf)
    return [p for p in parts if p]

def compile_evidence(pattern):
    """Lowercase alternatives are words (case-insensitive); alternatives with capitals are
    symbols (case-sensitive). Every alternative starts at a word boundary; alternatives of
    four or fewer letters must also end at one."""
    wrapped = []
    for alt in split_alternatives(pattern):
        letters = re.sub(r"\\[a-zA-Z]|\[[^\]]*\]|[^A-Za-z]", "", alt)
        has_upper = bool(re.search(r"[A-Z]", re.sub(r"\\[A-Za-z]", "", alt)))
        body = f"(?:{alt})" if has_upper else f"(?i:{alt})"
        tail = "(?![a-z])" if len(letters) <= 4 and alt[-1:].isalpha() else ""
        if re.fullmatch(r"[a-z]+ase", alt):
            # Enzyme-class suffixes name the class at the end of compound words
            # (aminotransferase, bisphosphatase); exclusions stop known traps.
            guard = "".join(f"(?<!{x})" for x in SUFFIX_EXCLUSIONS.get(alt, ()))
            wrapped.append(f"{guard}{body}")
        else:
            wrapped.append(f"(?<![A-Za-z0-9]){body}{tail}")
    return re.compile("|".join(wrapped))

SUFFIX_EXCLUSIONS = {"peptidase": ("trans", "Trans"), "sortase": ("exo", "Exo"), "hydrogenase": ("de", "De"), "reductase": ("oxido", "Oxido"), "oxidase": ("per", "Per"), "ligase": ("dia",)}

# A match that is the object of a regulator/transporter/sensor names that
# relation ("activator of glycolytic enzymes", "permease for ... thiamine").
PRECEDED = re.compile(
    r"(activator|regulator|regulatory protein|inhibitor|repressor|sensor|receptor|binding protein|"
    r"transporter|permease|exporter|importer|chaperone)s?\s+(of|for)\s+(\S+\s+){0,4}$",
    re.IGNORECASE,
)

LIKE = re.compile(r"^[\w/]*-(like(?![a-z])|,|\s+and\b)|^[\w/]*\s+like(?![a-z])|^[\w/-]*(\s+\w+){0,2}\s+like\s*$")
# A match followed by one of these names a relation to the thing, not the thing
# ("GTPase-activating protein", "toxin-antitoxin system", "protease inhibitor").
RELATION = re.compile(
    r"^[\w/+]*[- ](activating|activator|inactivating|inhibitor|inhibiting|associated|interacting|interaction|antitoxin|binding)(?![a-z])",
    re.IGNORECASE,
)



@cache
def _compiled(predicate: str) -> tuple[re.Pattern, ...]:
    return tuple(compile_evidence(x) for x in EVIDENCE.get(predicate, {}).get("text", ()))


def text_evidence(predicate: str, name: str, desc: str) -> str | None:
    """Return the matched wording when the Pfam name/description states the predicate."""
    text = f"{name} {desc}"
    found = None
    for pattern in _compiled(predicate):
        for m in pattern.finditer(text):
            if PRECEDED.search(text[:m.start()]):
                continue
            rest = text[m.end():]
            relation = RELATION.match(rest)
            if relation:
                kind = relation.group(1).lower()
                if kind in ("binding", "associated") and predicate.endswith(f"_{kind}"):
                    pass
                elif kind in ("binding", "associated", "antitoxin", "interaction"):
                    continue
                else:
                    return None  # the family is defined by a relation to the predicate
            if found is None and not LIKE.match(rest):
                found = m.group(0)
    return found


def enzyme_phrases(desc: str, enzyme_names: dict[str, list[str]]) -> dict[str, list[str]]:
    """ENZYME database names (>= 7 characters) occurring as phrases in a description."""
    words = desc.lower().replace(",", " ").replace(";", " ").split()
    hits = {}
    for i in range(len(words)):
        for j in range(i + 1, min(len(words), i + 8) + 1):
            phrase = " ".join(words[i:j]).strip(" ,;:()")
            if len(phrase) >= 7 and phrase in enzyme_names:
                hits[phrase] = enzyme_names[phrase]
    return hits


def resolve(predicate: str, vocabulary) -> str | None:
    """Vocabulary predicate for a proposed name (aliases applied), or None."""
    target = predicate if predicate in vocabulary else ALIASES.get(predicate)
    return target if target in vocabulary else None


def judge(predicate: str, name: str, desc: str, go_closure, enzyme_names) -> str | None:
    """Evidence string supporting ``predicate`` for one family, or None."""
    from sharur.predicates.mappings.kegg_map import get_predicates_for_ec

    if (name, predicate) in HOMONYMS:
        return None

    for go in EVIDENCE.get(predicate, {}).get("go", ()):
        if go in go_closure:
            return f"go:{go}"
    match = text_evidence(predicate, name, desc)
    if match:
        return f"text:{match}"
    for phrase, ecs in enzyme_phrases(desc, enzyme_names).items():
        for ec in ecs:
            if predicate in get_predicates_for_ec(ec):
                return f"enzyme:{ec} ({phrase})"
    return None
