"""
KEGG ortholog to predicate mappings.

Maps KEGG KO accessions to semantic predicates through the generated,
evidence-checked map in ``data/kegg_predicates.tsv``, and EC numbers to
predicates through :data:`EC_TO_PREDICATES`.
"""
import re
from functools import lru_cache
from pathlib import Path


# ============================================================================
# EC NUMBER TO PREDICATE MAPPINGS
# ============================================================================
# Maps EC class prefixes to predicates

EC_TO_PREDICATES: dict[str, list[str]] = {
    # EC 1: Oxidoreductases
    "1": ["oxidoreductase"],
    "1.1": ["oxidoreductase"],  # Acting on CH-OH donors
    "1.1.1": ["oxidoreductase", "dehydrogenase", "nad_binding"],  # With NAD+ or NADP+ as acceptor
    "1.2": ["oxidoreductase"],  # Acting on aldehyde/oxo
    "1.3": ["oxidoreductase"],  # Acting on CH-CH
    "1.4": ["oxidoreductase"],  # Acting on CH-NH2 donors
    "1.5": ["oxidoreductase"],  # Acting on CH-NH
    "1.6": ["oxidoreductase", "nad_binding"],  # Acting on NADH/NADPH
    "1.7": ["oxidoreductase", "nitrogen_metabolism"],  # Acting on N compounds
    "1.8": ["oxidoreductase", "sulfur_metabolism"],  # Acting on S compounds
    "1.9": ["oxidoreductase", "heme_binding"],  # Acting on heme
    "1.10": ["oxidoreductase"],  # Acting on diphenols
    "1.11": ["oxidoreductase", "peroxidase"],  # Acting on peroxide
    "1.12": ["oxidoreductase", "hydrogenase", "hydrogen_metabolism"],  # Acting on H2
    "1.13": ["oxidoreductase", "oxygenase"],  # Acting with single O2
    "1.14": ["oxidoreductase", "oxygenase"],  # Paired donors (ENZYME); mono-/dioxygenase set per sub-subclass
    "1.15": ["oxidoreductase", "superoxide_dismutase", "oxidative_stress"],  # Acting on superoxide
    "1.16": ["oxidoreductase", "metal_binding"],  # Oxidizing metal ions
    "1.17": ["oxidoreductase"],  # Acting on CH or CH2
    "1.18": ["oxidoreductase"],  # Acting on iron-sulfur proteins as donors
    "1.97": ["oxidoreductase"],  # Other oxidoreductases

    # Cofactor classes from Expasy ENZYME enzclass.txt descriptions
    "1.2.1": ["oxidoreductase", "nad_binding"],  # ENZYME enzclass: With NAD(+) or NADP(+) as acceptor.
    "1.2.8": ["oxidoreductase", "flavin_binding"],  # ENZYME enzclass: With a flavin or flavoprotein as acceptor.
    "1.3.1": ["oxidoreductase", "nad_binding"],  # ENZYME enzclass: With NAD(+) or NADP(+) as acceptor.
    "1.3.8": ["oxidoreductase", "flavin_binding"],  # ENZYME enzclass: With a flavin as acceptor.
    "1.4.1": ["oxidoreductase", "nad_binding"],  # ENZYME enzclass: With NAD(+) or NADP(+) as acceptor.
    "1.5.1": ["oxidoreductase", "nad_binding"],  # ENZYME enzclass: With NAD(+) or NADP(+) as acceptor.
    "1.5.8": ["oxidoreductase", "flavin_binding"],  # ENZYME enzclass: With a flavin as acceptor.
    "1.6.1": ["nad_binding", "oxidoreductase"],  # ENZYME enzclass: With NAD(+) or NADP(+) as acceptor.
    "1.6.8": ["nad_binding", "oxidoreductase", "flavin_binding"],  # ENZYME enzclass: With a flavin as acceptor.
    "1.7.1": ["nitrogen_metabolism", "oxidoreductase", "nad_binding"],  # ENZYME enzclass: With NAD(+) or NADP(+) as acceptor.
    "1.8.1": ["oxidoreductase", "sulfur_metabolism", "nad_binding"],  # ENZYME enzclass: With NAD(+) or NADP(+) as acceptor.
    "1.10.1": ["oxidoreductase", "nad_binding"],  # ENZYME enzclass: With NAD(+) or NADP(+) as acceptor.
    "1.12.1": ["hydrogen_metabolism", "hydrogenase", "oxidoreductase", "nad_binding"],  # ENZYME enzclass: With NAD(+) or NADP(+) as acceptor.
    "1.13.1": ["oxidoreductase", "oxygenase", "nad_binding"],  # ENZYME enzclass: With NADH or NADPH as one donor.
    "1.14.1": ["oxidoreductase", "oxygenase", "nad_binding"],  # ENZYME enzclass: With NADH or NADPH as one donor.
    "1.14.12": ["oxidoreductase", "oxygenase", "dioxygenase", "nad_binding"],  # ENZYME enzclass: With NADH or NADPH as one donor, and incorporation of two atoms of oxygen into one donor.
    "1.14.13": ["oxidoreductase", "oxygenase", "monooxygenase", "nad_binding"],  # ENZYME enzclass: With NADH or NADPH as one donor, and incorporation of one atom of oxygen.
    "1.14.14": ["oxidoreductase", "oxygenase", "monooxygenase", "flavin_binding"],  # ENZYME enzclass: With reduced flavin or flavoprotein as one donor, and incorporation of one atom of oxygen.
    "1.14.21": ["oxidoreductase", "oxygenase", "nad_binding"],  # ENZYME enzclass: With NADH or NADPH as one donor, and the other dehydrogenated.
    "1.13.11": ["oxidoreductase", "oxygenase", "dioxygenase"],  # ENZYME enzclass: With incorporation of two atoms of oxygen.
    "1.13.12": ["oxidoreductase", "oxygenase", "monooxygenase"],  # ENZYME enzclass: With incorporation of one atom of oxygen (internal monooxygenases or internal mixed function oxidases).
    "1.14.11": ["oxidoreductase", "oxygenase", "dioxygenase"],  # ENZYME enzclass: With 2-oxoglutarate as one donor, and incorporation of one atom each of oxygen into both donors.
    "1.14.15": ["oxidoreductase", "oxygenase", "monooxygenase"],  # ENZYME enzclass: With reduced iron-sulfur protein as one donor, and incorporation of one atom of oxygen.
    "1.14.16": ["oxidoreductase", "oxygenase", "monooxygenase"],  # ENZYME enzclass: With reduced pteridine as one donor, and incorporation of one atom of oxygen.
    "1.14.17": ["oxidoreductase", "oxygenase", "monooxygenase"],  # ENZYME enzclass: With reduced ascorbate as one donor, and incorporation of one atom of oxygen.
    "1.14.18": ["oxidoreductase", "oxygenase", "monooxygenase"],  # ENZYME enzclass: With another compound as one donor, and incorporation of one atom of oxygen.
    "1.16.1": ["metal_binding", "oxidoreductase", "nad_binding"],  # ENZYME enzclass: With NAD(+) or NADP(+) as acceptor.
    "1.16.8": ["metal_binding", "oxidoreductase", "flavin_binding"],  # ENZYME enzclass: With a flavin as acceptor.
    "1.17.1": ["oxidoreductase", "nad_binding"],  # ENZYME enzclass: With NAD(+) or NADP(+) as acceptor.
    "1.17.8": ["oxidoreductase", "flavin_binding"],  # ENZYME enzclass: With a flavin as acceptor.
    "1.18.1": ["oxidoreductase", "nad_binding"],  # ENZYME enzclass: With NAD(+) or NADP(+) as acceptor.
    "1.19.1": ["oxidoreductase", "nad_binding"],  # ENZYME enzclass: With NAD(+) or NADP(+) as acceptor.
    "1.20.1": ["oxidoreductase", "nad_binding"],  # ENZYME enzclass: With NAD(+) or NADP(+) as acceptor.
    "1.21.1": ["oxidoreductase", "nad_binding"],  # ENZYME enzclass: With NAD(+) or NADP(+) as acceptor.
    "1.22.1": ["oxidoreductase", "nad_binding"],  # ENZYME enzclass: With NAD(+) or NADP(+) as acceptor.
    "1.23.1": ["oxidoreductase", "nad_binding"],  # ENZYME enzclass: With NADH or NADPH as donor.

    # EC 2: Transferases
    "2": ["transferase"],
    "2.1": ["transferase"],  # Transferring one-carbon groups
    "2.1.1": ["transferase", "methyltransferase"],  # Methyltransferases
    "2.1.1.37": ["transferase", "methyltransferase", "dna_methylase"],  # DNA (cytosine-5-)-methyltransferase
    "2.1.1.72": ["transferase", "methyltransferase", "dna_methylase"],  # Site-specific DNA-methyltransferase (adenine-specific)
    "2.1.1.113": ["transferase", "methyltransferase", "dna_methylase"],  # Site-specific DNA-methyltransferase (cytosine-N4-specific)
    "2.3": ["transferase"],  # Acyltransferases
    "2.4": ["transferase", "glycosyltransferase"],  # Glycosyltransferases
    "2.5": ["transferase"],  # Alkyl/aryl transferases
    "2.6": ["transferase"],  # Transferring nitrogenous groups
    "2.6.1": ["transferase", "aminotransferase", "plp_binding"],  # Transaminases
    "2.7": ["transferase"],  # Transferring phosphorus-containing groups
    "2.7.1": ["transferase", "kinase", "atp_binding"],  # Kinases (alcohol acceptors)
    "2.7.2": ["transferase", "kinase", "atp_binding"],  # Kinases (carboxyl acceptors)
    "2.7.4": ["transferase", "kinase", "atp_binding"],  # Phosphotransferases
    "2.7.3": ["transferase", "kinase", "atp_binding"],  # Phosphotransferases with a nitrogenous acceptor
    "2.7.6": ["transferase", "kinase", "atp_binding"],  # Diphosphotransferases
    "2.7.9": ["transferase", "kinase", "atp_binding"],  # Phosphotransferases with paired acceptors
    "2.7.10": ["transferase", "kinase", "atp_binding"],  # Protein-tyrosine kinases
    "2.7.7": ["transferase", "nucleotidyltransferase"],  # Nucleotidyltransferases
    "2.7.11": ["transferase", "serine_threonine_kinase", "kinase"],  # Ser/Thr kinases
    "2.7.12": ["transferase", "kinase"],  # Dual-specificity kinases
    "2.7.13": ["transferase", "sensor_kinase", "two_component"],  # Histidine kinases
    "2.8": ["transferase", "sulfurtransferase"],  # S-containing groups

    # EC 3: Hydrolases
    "3": ["hydrolase"],
    "3.1": ["hydrolase"],  # Ester bonds (esterase is claimed at 3.1.1)
    "3.1.1": ["hydrolase", "esterase"],  # Carboxylic ester hydrolases
    "3.1.3": ["hydrolase", "phosphatase"],  # Phosphoric monoester
    "3.1.4": ["hydrolase", "phosphodiesterase"],  # Phosphoric diester
    "3.1.11": ["hydrolase", "nuclease", "exonuclease"],  # Exodeoxyribonucleases
    "3.1.13": ["hydrolase", "nuclease", "exonuclease", "rnase"],  # Exoribonucleases
    "3.1.21": ["hydrolase", "nuclease", "endonuclease", "dnase"],  # Endodeoxyribonucleases
    "3.1.22": ["hydrolase", "nuclease", "endonuclease", "dnase"],  # Type I restriction
    "3.1.26": ["hydrolase", "nuclease", "endonuclease", "rnase"],  # Endoribonucleases
    "3.2": ["hydrolase", "glycosidase"],  # Glycosidic bonds
    "3.2.1": ["hydrolase", "glycosidase", "carbohydrate_active"],  # Glycosidases
    "3.4": ["hydrolase", "protease", "peptidase"],  # Peptide bonds
    "3.4.11": ["hydrolase", "protease", "peptidase"],  # Aminopeptidases
    "3.4.21": ["hydrolase", "protease", "peptidase"],  # Serine endopeptidases
    "3.4.22": ["hydrolase", "protease", "cysteine_protease"],  # Cysteine proteases
    "3.4.23": ["hydrolase", "protease", "peptidase"],  # Aspartic endopeptidases
    "3.4.24": ["hydrolase", "protease", "peptidase", "metal_binding"],  # Metalloendopeptidases
    "3.5": ["hydrolase"],  # Acting on C-N bonds other than peptide bonds
    "3.5.1": ["hydrolase", "amidase"],  # In linear amides
    "3.6": ["hydrolase"],  # Acting on acid anhydrides
    "3.6.1": ["hydrolase"],  # In phosphorus-containing anhydrides
    "3.6.3": ["hydrolase", "atpase", "atp_binding", "transporter"],  # Catalyzing transmembrane movement
    "3.6.4": ["hydrolase", "atpase", "atp_binding"],  # Acting on ATP; cellular movement
    "3.6.5": ["hydrolase", "gtp_binding"],  # Acting on GTP

    # EC 4: Lyases
    "4": ["lyase"],
    "4.1": ["lyase"],  # Carbon-carbon lyases
    "4.1.1": ["lyase", "decarboxylase"],  # Carboxy-lyases
    "4.1.2": ["lyase"],  # Aldehyde lyases
    "4.2": ["lyase"],  # C-O lyases
    "4.2.1": ["lyase", "dehydratase"],  # Hydro-lyases
    "4.3": ["lyase"],  # C-N lyases
    "4.4": ["lyase"],  # C-S lyases
    "4.6": ["lyase"],  # P-O lyases

    # EC 5: Isomerases
    "5": ["isomerase"],
    "5.1": ["isomerase", "racemase"],  # Racemases/epimerases
    "5.2": ["isomerase"],  # Cis-trans isomerases
    "5.3": ["isomerase"],  # Intramolecular oxidoreductases
    "5.4": ["isomerase", "mutase"],  # Intramolecular transferases
    "5.5": ["isomerase"],  # Intramolecular lyases
    "5.6": ["isomerase"],  # Isomerases altering macromolecular conformation

    # EC 6: Ligases
    "6": ["ligase", "synthetase"],
    "6.1": ["ligase", "synthetase"],  # C-O bonds
    "6.1.1": ["ligase", "trna_synthetase", "translation"],  # Aminoacyl-tRNA synthetases
    "6.2": ["ligase", "synthetase"],  # C-S bonds
    "6.3": ["ligase", "synthetase"],  # C-N bonds
    "6.3.1": ["ligase", "synthetase"],  # Acid-ammonia (or amine) ligases
    "6.3.2": ["ligase", "synthetase"],  # Acid-amino-acid ligases
    "6.4": ["ligase", "carboxylase", "biotin_binding"],  # C-C bonds
    "6.5": ["ligase", "ligase_dna"],  # Forming phosphoric ester bonds
    "6.6": ["ligase"],  # N-metal bonds

    # EC 7: Translocases
    "7": ["translocase", "transporter", "membrane"],
    "7.1": ["translocase", "transporter"],  # Translocation of hydrons
    "7.1.1": ["translocase", "transporter", "electron_transport"],  # Hydron translocation linked to oxidoreduction
    "7.2": ["translocase", "transporter", "ion_transporter"],  # Linked to hydrolysis
    "7.3": ["translocase", "transporter"],  # Linked to decarboxylation
    "7.4": ["translocase", "transporter"],  # Linked to redox
    "7.5": ["translocase", "transporter"],  # Translocation of carbohydrates
    "7.5.2": ["translocase", "transporter", "atp_binding"],  # Carbohydrate translocation linked to ATP hydrolysis
    "7.6": ["translocase", "transporter"],  # Translocation of other compounds
    "7.6.2": ["translocase", "transporter", "atp_binding"],  # Translocation linked to ATP hydrolysis
}


# ============================================================================
# KO -> PREDICATE MAP (generated)
# ============================================================================
# ``data/kegg_predicates.tsv`` is generated by ``scripts/build_kegg_predicate_map.py``:
# every (KO, predicate) pair in it is supported by KEGG's own information (EC,
# BRITE placement, module membership), HydDB reference labels for KOs KEGG names
# as hydrogenases, or the KO's KEGG symbols/name, and records that evidence.
# Proposals live in ``data/kegg_predicate_proposals.tsv`` and
# ``data/kegg_predicate_proposal_patterns.tsv``; edit those and rebuild.

KEGG_MAPPING_FILE = Path(__file__).with_name("data") / "kegg_predicates.tsv"


def _load_kegg_mapping_file(path: Path) -> tuple[dict[str, list[str]], dict[str, dict[str, str]]]:
    predicates: dict[str, list[str]] = {}
    evidence: dict[str, dict[str, str]] = {}
    with path.open() as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            if not raw_line.strip() or raw_line.startswith("#"):
                continue
            parts = raw_line.rstrip("\n").split("\t")
            if len(parts) != 3:
                raise ValueError(f"{path}:{line_number}: expected ko, predicates, evidence")
            ko, preds, ev = parts
            pred_list = sorted(p for p in preds.split(",") if p)
            ev_map = dict(item.split("=", 1) for item in ev.split(";") if item)
            if set(ev_map) != set(pred_list):
                raise ValueError(f"{path}:{line_number}: evidence does not cover every predicate")
            predicates[ko] = pred_list
            evidence[ko] = ev_map
    return predicates, evidence


KEGG_TO_PREDICATES, KEGG_EVIDENCE = _load_kegg_mapping_file(KEGG_MAPPING_FILE)


def parse_ec_numbers(definition: str) -> list[str]:
    """Extract EC numbers from KEGG definition."""
    # Pattern: [EC:1.2.3.4] or EC:1.2.3.4 or [EC:1.2.3.4 1.2.3.5]
    ec_pattern = r"\[?EC:([^\]]+)\]?"
    match = re.search(ec_pattern, definition)
    if match:
        ec_string = match.group(1)
        # Split on space for multiple EC numbers
        return [ec.strip() for ec in ec_string.split() if ec.strip()]
    return []


def get_predicates_for_ec(ec_number: str) -> list[str]:
    """Get predicates for an EC number."""
    predicates = set()

    # Try increasingly specific prefixes
    parts = ec_number.split(".")
    for i in range(len(parts), 0, -1):
        prefix = ".".join(parts[:i])
        if prefix in EC_TO_PREDICATES:
            predicates.update(EC_TO_PREDICATES[prefix])
            break  # Use most specific match

    return list(predicates)


def get_predicates_for_kegg(
    ko_id: str,
    definition: str = "",
) -> list[str]:
    """
    Get predicates for a KEGG ortholog.

    Args:
        ko_id: KEGG ortholog ID (e.g., "K00001")
        definition: KEGG definition text; used only for KOs absent from the
            generated map

    Returns:
        List of predicate IDs
    """
    if ko_id in KEGG_TO_PREDICATES:
        return list(KEGG_TO_PREDICATES[ko_id])

    # KOs newer than the generated map: KEGG's own EC numbers only.
    predicates = set()
    for ec in parse_ec_numbers(_resolve_ko_definition(ko_id, definition)):
        predicates.update(get_predicates_for_ec(ec))
    return sorted(predicates)


def _resolve_ko_definition(ko_id: str, definition: str = "") -> str:
    """Use local KO metadata when an annotation only stores score labels."""
    if not _is_informative_definition(ko_id, definition):
        return _load_ko_definitions().get(ko_id, definition or "")
    return definition


def _is_informative_definition(ko_id: str, definition: str = "") -> bool:
    """Return whether a KEGG definition can drive EC mapping."""
    text = (definition or "").strip()
    if not text:
        return False
    if text == ko_id:
        return False
    return not (text == "GA" or text.startswith("evalue_"))


@lru_cache(maxsize=1)
def _load_ko_definitions() -> dict[str, str]:
    """Load optional local KOFAM KO definitions."""
    ko_list = _ko_list_path()
    if not ko_list.exists():
        return {}

    definitions: dict[str, str] = {}
    with open(ko_list) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        try:
            ko_idx = header.index("knum")
            definition_idx = header.index("definition")
        except ValueError:
            return {}

        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(ko_idx, definition_idx):
                continue
            definitions[fields[ko_idx]] = fields[definition_idx]

    return definitions


def _ko_list_path() -> Path:
    """Return the repo-local KO list path."""
    return Path(__file__).resolve().parents[3] / "data" / "reference" / "ko_list"


__all__ = [
    "EC_TO_PREDICATES",
    "KEGG_EVIDENCE",
    "KEGG_TO_PREDICATES",
    "get_predicates_for_ec",
    "get_predicates_for_kegg",
    "parse_ec_numbers",
]
