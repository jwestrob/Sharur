"""
CAZy family to predicate mappings.

Maps CAZy carbohydrate-active enzyme families to semantic predicates.
CAZy families are organized into major classes:
- GH: Glycoside Hydrolases
- GT: Glycosyltransferases
- PL: Polysaccharide Lyases
- CE: Carbohydrate Esterases
- CBM: Carbohydrate-Binding Modules
- AA: Auxiliary Activities
"""

from typing import Optional
import re


# ============================================================================
# CAZY CLASS TO PREDICATE MAPPINGS
# ============================================================================
# Maps CAZy class prefix to general predicates

CAZY_CLASS_PREDICATES: dict[str, list[str]] = {
    "GH": ["carbohydrate_active", "glycoside_hydrolase", "hydrolase"],
    "GT": ["carbohydrate_active", "glycosyltransferase_cazy", "transferase"],
    "PL": ["carbohydrate_active", "polysaccharide_lyase", "lyase"],
    "CE": ["carbohydrate_active", "carbohydrate_esterase", "esterase", "hydrolase"],
    "CBM": ["carbohydrate_active", "carbohydrate_binding"],
    "AA": ["carbohydrate_active", "auxiliary_activity", "oxidoreductase"],
}


# ============================================================================
# SUBSTRATE-SPECIFIC FAMILY MAPPINGS
# ============================================================================
# Maps specific CAZy families to substrate predicates

CAZY_TO_PREDICATES: dict[str, list[str]] = {
    # -------------------------------------------------------------------------
    # CELLULASES (Cellulose degradation)
    # -------------------------------------------------------------------------
    "GH5": ["cellulase", "glycoside_hydrolase", "carbohydrate_active"],  # Major cellulase family
    "GH6": ["cellulase", "glycoside_hydrolase", "carbohydrate_active"],  # Cellobiohydrolase
    "GH7": ["cellulase", "glycoside_hydrolase", "carbohydrate_active"],  # Cellobiohydrolase
    "GH9": ["cellulase", "glycoside_hydrolase", "carbohydrate_active"],  # Endoglucanase
    "GH44": ["cellulase", "glycoside_hydrolase", "carbohydrate_active"],  # Endoglucanase
    "GH45": ["cellulase", "glycoside_hydrolase", "carbohydrate_active"],  # Endoglucanase
    "GH48": ["cellulase", "glycoside_hydrolase", "carbohydrate_active"],  # Cellobiohydrolase
    "GH12": ["cellulase", "xylanase", "glycoside_hydrolase", "carbohydrate_active"],  # Xyloglucanase
    "GH74": ["cellulase", "glycoside_hydrolase", "carbohydrate_active"],  # Xyloglucanase

    # Beta-glucosidases (complete cellulose degradation)
    "GH1": ["glycoside_hydrolase", "carbohydrate_active"],  # Beta-glucosidase
    "GH3": ["glycoside_hydrolase", "carbohydrate_active"],  # Beta-glucosidase

    # -------------------------------------------------------------------------
    # CHITINASES (Chitin degradation)
    # -------------------------------------------------------------------------
    "GH18": ["chitinase", "glycoside_hydrolase", "carbohydrate_active"],  # Major chitinase
    "GH19": ["chitinase", "glycoside_hydrolase", "carbohydrate_active"],  # Plant chitinase
    "GH20": ["chitinase", "glycoside_hydrolase", "carbohydrate_active"],  # Beta-hexosaminidase
    "GH23": ["glycoside_hydrolase", "carbohydrate_active", "lytic_transglycosylase"],  # Lysozyme/LT
    "GH46": ["chitinase", "glycoside_hydrolase", "carbohydrate_active"],  # Chitosanase

    # -------------------------------------------------------------------------
    # AMYLASES (Starch degradation)
    # -------------------------------------------------------------------------
    "GH13": ["amylase", "glycoside_hydrolase", "carbohydrate_active"],  # Alpha-amylase superfamily
    "GH14": ["amylase", "glycoside_hydrolase", "carbohydrate_active"],  # Beta-amylase
    "GH15": ["amylase", "glycoside_hydrolase", "carbohydrate_active"],  # Glucoamylase
    "GH57": ["amylase", "glycoside_hydrolase", "carbohydrate_active"],  # Amylopullulanase
    "GH77": ["amylase", "glycoside_hydrolase", "carbohydrate_active"],  # Amylomaltase
    "GH97": ["amylase", "glycoside_hydrolase", "carbohydrate_active"],  # Alpha-glucosidase
    "GH119": ["amylase", "glycoside_hydrolase", "carbohydrate_active"],  # Alpha-amylase

    # -------------------------------------------------------------------------
    # XYLANASES (Xylan/Hemicellulose degradation)
    # -------------------------------------------------------------------------
    "GH10": ["xylanase", "glycoside_hydrolase", "carbohydrate_active"],  # Endo-xylanase
    "GH11": ["xylanase", "glycoside_hydrolase", "carbohydrate_active"],  # Endo-xylanase
    "GH30": ["xylanase", "glycoside_hydrolase", "carbohydrate_active"],  # Glucuronoxylanase
    "GH8": ["xylanase", "glycoside_hydrolase", "carbohydrate_active"],  # Reducing-end xylanase
    "GH43": ["xylanase", "glycoside_hydrolase", "carbohydrate_active"],  # Beta-xylosidase
    "GH67": ["xylanase", "glycoside_hydrolase", "carbohydrate_active"],  # Alpha-glucuronidase

    # -------------------------------------------------------------------------
    # PECTINASES (Pectin degradation)
    # -------------------------------------------------------------------------
    "GH28": ["pectinase", "glycoside_hydrolase", "carbohydrate_active"],  # Polygalacturonase
    "GH78": ["pectinase", "glycoside_hydrolase", "carbohydrate_active"],  # Alpha-rhamnosidase
    "GH88": ["pectinase", "glycoside_hydrolase", "carbohydrate_active"],  # Unsaturated glucuronyl hydrolase
    "GH105": ["pectinase", "glycoside_hydrolase", "carbohydrate_active"],  # Unsaturated rhamnogalacturonyl hydrolase
    "PL1": ["pectinase", "polysaccharide_lyase", "carbohydrate_active"],  # Pectate lyase
    "PL2": ["pectinase", "polysaccharide_lyase", "carbohydrate_active"],  # Pectate lyase
    "PL3": ["pectinase", "polysaccharide_lyase", "carbohydrate_active"],  # Pectate lyase
    "PL9": ["pectinase", "polysaccharide_lyase", "carbohydrate_active"],  # Pectate lyase
    "PL10": ["pectinase", "polysaccharide_lyase", "carbohydrate_active"],  # Pectate lyase
    "PL22": ["pectinase", "polysaccharide_lyase", "carbohydrate_active"],  # Oligogalacturonate lyase
    "CE8": ["pectinase", "carbohydrate_esterase", "carbohydrate_active"],  # Pectin methylesterase
    "CE12": ["pectinase", "carbohydrate_esterase", "carbohydrate_active"],  # Pectin acetylesterase

    # -------------------------------------------------------------------------
    # MANNANASES (Mannan degradation)
    # -------------------------------------------------------------------------
    "GH26": ["mannanase", "glycoside_hydrolase", "carbohydrate_active"],  # Beta-mannanase
    "GH113": ["mannanase", "glycoside_hydrolase", "carbohydrate_active"],  # Beta-mannanase
    "GH130": ["mannanase", "glycoside_hydrolase", "carbohydrate_active"],  # Beta-mannooligosaccharide phosphorylase
    "GH2": ["mannanase", "glycoside_hydrolase", "carbohydrate_active"],  # Beta-mannosidase

    # -------------------------------------------------------------------------
    # LYTIC POLYSACCHARIDE MONOOXYGENASES (LPMOs)
    # -------------------------------------------------------------------------
    "AA9": ["lytic_polysaccharide_monooxygenase", "auxiliary_activity", "carbohydrate_active", "copper_binding"],  # LPMO cellulose
    "AA10": ["lytic_polysaccharide_monooxygenase", "auxiliary_activity", "carbohydrate_active", "copper_binding"],  # LPMO chitin/cellulose
    "AA11": ["lytic_polysaccharide_monooxygenase", "auxiliary_activity", "carbohydrate_active", "copper_binding"],  # LPMO chitin
    "AA13": ["lytic_polysaccharide_monooxygenase", "auxiliary_activity", "carbohydrate_active", "copper_binding"],  # LPMO starch
    "AA14": ["lytic_polysaccharide_monooxygenase", "auxiliary_activity", "carbohydrate_active", "copper_binding"],  # LPMO xylan
    "AA15": ["lytic_polysaccharide_monooxygenase", "auxiliary_activity", "carbohydrate_active", "copper_binding"],  # LPMO various
    "AA16": ["lytic_polysaccharide_monooxygenase", "auxiliary_activity", "carbohydrate_active", "copper_binding"],  # LPMO cellulose
    "AA17": ["lytic_polysaccharide_monooxygenase", "auxiliary_activity", "carbohydrate_active", "copper_binding"],  # LPMO pectin

    # -------------------------------------------------------------------------
    # OTHER AUXILIARY ACTIVITIES
    # -------------------------------------------------------------------------
    "AA1": ["auxiliary_activity", "carbohydrate_active", "oxidase", "copper_binding"],  # Laccase
    "AA2": ["auxiliary_activity", "carbohydrate_active", "peroxidase", "heme_binding"],  # Class II peroxidase
    "AA3": ["auxiliary_activity", "carbohydrate_active", "oxidase", "fad_binding"],  # GMC oxidoreductases
    "AA4": ["auxiliary_activity", "carbohydrate_active", "oxidase", "fad_binding"],  # Vanillyl-alcohol oxidase
    "AA5": ["auxiliary_activity", "carbohydrate_active", "oxidase", "copper_binding"],  # Glyoxal/galactose oxidase
    "AA6": ["auxiliary_activity", "carbohydrate_active", "reductase"],  # 1,4-benzoquinone reductase
    "AA7": ["auxiliary_activity", "carbohydrate_active", "oxidase", "fad_binding"],  # Glucooligosaccharide oxidase
    "AA8": ["auxiliary_activity", "carbohydrate_active", "iron_binding"],  # Iron reductase domain

    # -------------------------------------------------------------------------
    # GLYCOSYLTRANSFERASES (Selected families)
    # -------------------------------------------------------------------------
    "GT1": ["glycosyltransferase_cazy", "carbohydrate_active", "transferase"],  # UDP-glucuronosyltransferase
    "GT2": ["glycosyltransferase_cazy", "carbohydrate_active", "transferase"],  # Cellulose synthase
    "GT4": ["glycosyltransferase_cazy", "carbohydrate_active", "transferase"],  # Sucrose synthase
    "GT5": ["glycosyltransferase_cazy", "carbohydrate_active", "transferase"],  # Starch synthase
    "GT20": ["glycosyltransferase_cazy", "carbohydrate_active", "transferase"],  # Alpha-amylase superfamily
    "GT28": ["glycosyltransferase_cazy", "carbohydrate_active", "transferase", "peptidoglycan"],  # MurG
    "GT35": ["glycosyltransferase_cazy", "carbohydrate_active", "transferase"],  # Glycogen phosphorylase
    "GT51": ["glycosyltransferase_cazy", "carbohydrate_active", "transferase", "peptidoglycan"],  # Murein polymerase

    # -------------------------------------------------------------------------
    # CARBOHYDRATE ESTERASES
    # -------------------------------------------------------------------------
    "CE1": ["carbohydrate_esterase", "carbohydrate_active", "esterase"],  # Acetyl xylan esterase
    "CE2": ["carbohydrate_esterase", "carbohydrate_active", "esterase"],  # Acetyl xylan esterase
    "CE3": ["carbohydrate_esterase", "carbohydrate_active", "esterase"],  # Acetyl xylan esterase
    "CE4": ["carbohydrate_esterase", "carbohydrate_active", "esterase"],  # Acetyl xylan/chitin deacetylase
    "CE5": ["carbohydrate_esterase", "carbohydrate_active", "esterase"],  # Cutinase
    "CE6": ["carbohydrate_esterase", "carbohydrate_active", "esterase"],  # Acetyl xylan esterase
    "CE7": ["carbohydrate_esterase", "carbohydrate_active", "esterase"],  # Acetyl xylan esterase
    "CE9": ["carbohydrate_esterase", "carbohydrate_active", "esterase"],  # N-acetylglucosamine 6-phosphate deacetylase
    "CE10": ["carbohydrate_esterase", "carbohydrate_active", "esterase"],  # Arylesterase
    "CE11": ["carbohydrate_esterase", "carbohydrate_active", "esterase", "peptidoglycan"],  # UDP-3-O-acyl-N-acetylglucosamine deacetylase
    "CE14": ["carbohydrate_esterase", "carbohydrate_active", "esterase"],  # N-acetyl-1-D-myo-inosityl-2-amino-2-deoxy-alpha-D-glucopyranoside deacetylase
    "CE15": ["carbohydrate_esterase", "carbohydrate_active", "esterase"],  # 4-O-methyl-glucuronoyl methylesterase
    "CE16": ["carbohydrate_esterase", "carbohydrate_active", "esterase"],  # Acetylesterase

    # -------------------------------------------------------------------------
    # CARBOHYDRATE-BINDING MODULES (Selected)
    # -------------------------------------------------------------------------
    "CBM1": ["carbohydrate_binding", "carbohydrate_active"],  # Cellulose-binding
    "CBM2": ["carbohydrate_binding", "carbohydrate_active"],  # Cellulose/xylan-binding
    "CBM3": ["carbohydrate_binding", "carbohydrate_active"],  # Cellulose-binding (scaffoldin)
    "CBM4": ["carbohydrate_binding", "carbohydrate_active"],  # Amorphous cellulose-binding
    "CBM5": ["carbohydrate_binding", "carbohydrate_active", "chitinase"],  # Chitin-binding
    "CBM6": ["carbohydrate_binding", "carbohydrate_active"],  # Various substrates
    "CBM9": ["carbohydrate_binding", "carbohydrate_active"],  # Cellulose-binding
    "CBM10": ["carbohydrate_binding", "carbohydrate_active"],  # Cellulose-binding
    "CBM12": ["carbohydrate_binding", "carbohydrate_active", "chitinase"],  # Chitin-binding
    "CBM13": ["carbohydrate_binding", "carbohydrate_active"],  # Lectin-like
    "CBM14": ["carbohydrate_binding", "carbohydrate_active", "chitinase"],  # Chitin-binding
    "CBM16": ["carbohydrate_binding", "carbohydrate_active"],  # Cellulose/xylan-binding
    "CBM18": ["carbohydrate_binding", "carbohydrate_active", "chitinase"],  # Chitin-binding
    "CBM19": ["carbohydrate_binding", "carbohydrate_active", "chitinase"],  # Chitin-binding
    "CBM20": ["carbohydrate_binding", "carbohydrate_active", "amylase"],  # Starch-binding
    "CBM21": ["carbohydrate_binding", "carbohydrate_active", "amylase"],  # Starch-binding
    "CBM22": ["carbohydrate_binding", "carbohydrate_active", "xylanase"],  # Xylan-binding
    "CBM25": ["carbohydrate_binding", "carbohydrate_active", "amylase"],  # Starch-binding
    "CBM26": ["carbohydrate_binding", "carbohydrate_active", "amylase"],  # Starch-binding
    "CBM32": ["carbohydrate_binding", "carbohydrate_active"],  # Galactose-binding
    "CBM33": ["lytic_polysaccharide_monooxygenase", "carbohydrate_active", "copper_binding"],  # Now AA10
    "CBM34": ["carbohydrate_binding", "carbohydrate_active", "amylase"],  # Granular starch-binding
    "CBM35": ["carbohydrate_binding", "carbohydrate_active", "xylanase"],  # Xylan-binding
    "CBM42": ["carbohydrate_binding", "carbohydrate_active"],  # Arabinofuranose-binding
    "CBM46": ["carbohydrate_binding", "carbohydrate_active", "cellulase"],  # Cellulose-binding
    "CBM48": ["carbohydrate_binding", "carbohydrate_active", "amylase"],  # Glycogen-binding
    "CBM50": ["carbohydrate_binding", "carbohydrate_active", "peptidoglycan"],  # LysM
    "CBM54": ["carbohydrate_binding", "carbohydrate_active", "xylanase"],  # Xylan-binding
    "CBM57": ["carbohydrate_binding", "carbohydrate_active"],  # Starch/cellulose-binding
    "CBM61": ["carbohydrate_binding", "carbohydrate_active"],  # Beta-galactan-binding
    "CBM63": ["carbohydrate_binding", "carbohydrate_active", "cellulase"],  # Cellulose-binding

    # -------------------------------------------------------------------------
    # POLYSACCHARIDE LYASES (Additional)
    # -------------------------------------------------------------------------
    "PL4": ["polysaccharide_lyase", "carbohydrate_active"],  # Rhamnogalacturonan lyase
    "PL5": ["polysaccharide_lyase", "carbohydrate_active"],  # Alginate lyase
    "PL6": ["polysaccharide_lyase", "carbohydrate_active"],  # Alginate lyase
    "PL7": ["polysaccharide_lyase", "carbohydrate_active"],  # Alginate lyase
    "PL8": ["polysaccharide_lyase", "carbohydrate_active"],  # Hyaluronate lyase
    "PL11": ["polysaccharide_lyase", "carbohydrate_active"],  # Rhamnogalacturonan lyase
    "PL12": ["polysaccharide_lyase", "carbohydrate_active"],  # Heparin sulfate lyase
    "PL14": ["polysaccharide_lyase", "carbohydrate_active"],  # Alginate lyase
    "PL15": ["polysaccharide_lyase", "carbohydrate_active"],  # Oligo-alginate lyase
    "PL17": ["polysaccharide_lyase", "carbohydrate_active"],  # Alginate lyase
    "PL18": ["polysaccharide_lyase", "carbohydrate_active"],  # Alginate lyase
    "PL21": ["polysaccharide_lyase", "carbohydrate_active"],  # Heparin lyase

    # -------------------------------------------------------------------------
    # ADDITIONAL GH FAMILIES
    # -------------------------------------------------------------------------
    "GH4": ["glycoside_hydrolase", "carbohydrate_active", "nad_binding"],  # 6-phospho-beta-glucosidase
    "GH16": ["glycoside_hydrolase", "carbohydrate_active"],  # Various (xyloglucanase, laminarinase)
    "GH17": ["glycoside_hydrolase", "carbohydrate_active"],  # 1,3-beta-glucanase
    "GH19": ["chitinase", "glycoside_hydrolase", "carbohydrate_active"],  # Chitinase
    "GH22": ["glycoside_hydrolase", "carbohydrate_active"],  # Lysozyme
    "GH24": ["glycoside_hydrolase", "carbohydrate_active"],  # Lysozyme
    "GH25": ["glycoside_hydrolase", "carbohydrate_active", "lytic_transglycosylase"],  # N,O-diacetylmuramidase
    "GH29": ["glycoside_hydrolase", "carbohydrate_active"],  # Alpha-L-fucosidase
    "GH31": ["glycoside_hydrolase", "carbohydrate_active"],  # Alpha-glucosidase
    "GH32": ["glycoside_hydrolase", "carbohydrate_active"],  # Invertase
    "GH33": ["glycoside_hydrolase", "carbohydrate_active"],  # Sialidase
    "GH35": ["glycoside_hydrolase", "carbohydrate_active"],  # Beta-galactosidase
    "GH36": ["glycoside_hydrolase", "carbohydrate_active"],  # Alpha-galactosidase
    "GH37": ["glycoside_hydrolase", "carbohydrate_active"],  # Alpha,alpha-trehalase
    "GH38": ["glycoside_hydrolase", "carbohydrate_active"],  # Alpha-mannosidase
    "GH39": ["glycoside_hydrolase", "carbohydrate_active"],  # Beta-xylosidase
    "GH42": ["glycoside_hydrolase", "carbohydrate_active"],  # Beta-galactosidase
    "GH47": ["glycoside_hydrolase", "carbohydrate_active"],  # Alpha-mannosidase
    "GH51": ["glycoside_hydrolase", "carbohydrate_active"],  # Alpha-L-arabinofuranosidase
    "GH52": ["glycoside_hydrolase", "carbohydrate_active"],  # Beta-xylosidase
    "GH53": ["glycoside_hydrolase", "carbohydrate_active"],  # Endo-1,4-beta-galactanase
    "GH54": ["glycoside_hydrolase", "carbohydrate_active"],  # Alpha-L-arabinofuranosidase
    "GH55": ["glycoside_hydrolase", "carbohydrate_active"],  # 1,3-beta-glucanase
    "GH62": ["glycoside_hydrolase", "carbohydrate_active"],  # Alpha-L-arabinofuranosidase
    "GH63": ["glycoside_hydrolase", "carbohydrate_active"],  # Mannosyl-oligosaccharide glucosidase
    "GH64": ["glycoside_hydrolase", "carbohydrate_active"],  # 1,3-beta-glucanase
    "GH65": ["glycoside_hydrolase", "carbohydrate_active"],  # Maltose phosphorylase
    "GH68": ["glycoside_hydrolase", "carbohydrate_active"],  # Levansucrase
    "GH70": ["glycoside_hydrolase", "carbohydrate_active"],  # Glucansucrase
    "GH71": ["glycoside_hydrolase", "carbohydrate_active"],  # Mutanase
    "GH72": ["glycoside_hydrolase", "carbohydrate_active"],  # 1,3-beta-glucanosyltransglycosylase
    "GH73": ["glycoside_hydrolase", "carbohydrate_active", "peptidoglycan"],  # Lysozyme
    "GH79": ["glycoside_hydrolase", "carbohydrate_active"],  # Beta-glucuronidase
    "GH81": ["glycoside_hydrolase", "carbohydrate_active"],  # 1,3-beta-glucanase
    "GH84": ["glycoside_hydrolase", "carbohydrate_active"],  # N-acetyl-beta-glucosaminidase
    "GH85": ["glycoside_hydrolase", "carbohydrate_active"],  # Endo-beta-N-acetylglucosaminidase
    "GH92": ["glycoside_hydrolase", "carbohydrate_active"],  # Alpha-mannosidase
    "GH94": ["glycoside_hydrolase", "carbohydrate_active"],  # Cellobiose phosphorylase
    "GH95": ["glycoside_hydrolase", "carbohydrate_active"],  # Alpha-L-fucosidase
    "GH99": ["glycoside_hydrolase", "carbohydrate_active"],  # Glycoprotein endo-alpha-mannosidase
    "GH102": ["glycoside_hydrolase", "carbohydrate_active", "peptidoglycan"],  # Peptidoglycan lytic transglycosylase
    "GH103": ["glycoside_hydrolase", "carbohydrate_active", "peptidoglycan"],  # Peptidoglycan lytic transglycosylase
    "GH104": ["glycoside_hydrolase", "carbohydrate_active", "peptidoglycan"],  # Peptidoglycan lytic transglycosylase
    "GH108": ["glycoside_hydrolase", "carbohydrate_active"],  # N-acetylmuramidase
    "GH109": ["glycoside_hydrolase", "carbohydrate_active"],  # Alpha-N-acetylgalactosaminidase
    "GH110": ["glycoside_hydrolase", "carbohydrate_active"],  # Alpha-galactosidase
    "GH115": ["glycoside_hydrolase", "carbohydrate_active", "xylanase"],  # Alpha-glucuronidase
    "GH116": ["glycoside_hydrolase", "carbohydrate_active"],  # Beta-glucosidase
    "GH117": ["glycoside_hydrolase", "carbohydrate_active"],  # Alpha-neoagarobiose hydrolase
    "GH120": ["glycoside_hydrolase", "carbohydrate_active"],  # Beta-xylosidase
    "GH121": ["glycoside_hydrolase", "carbohydrate_active"],  # Beta-L-arabinobiosidase
    "GH123": ["glycoside_hydrolase", "carbohydrate_active"],  # Beta-N-acetylgalactosaminidase
    "GH125": ["glycoside_hydrolase", "carbohydrate_active"],  # Exo-alpha-1,6-mannosidase
    "GH127": ["glycoside_hydrolase", "carbohydrate_active"],  # Beta-L-arabinofuranosidase
    "GH128": ["glycoside_hydrolase", "carbohydrate_active"],  # Endo-1,3-beta-glucanase
    "GH129": ["glycoside_hydrolase", "carbohydrate_active"],  # Alpha-N-acetylgalactosaminidase
    "GH131": ["glycoside_hydrolase", "carbohydrate_active", "cellulase"],  # Exo-1,3/1,4-beta-glucanase
    "GH133": ["glycoside_hydrolase", "carbohydrate_active", "amylase"],  # Amylo-alpha-1,6-glucosidase
    "GH144": ["glycoside_hydrolase", "carbohydrate_active"],  # Endo-beta-1,2-glucanase
}


# ============================================================================
# PATTERN-BASED FAMILY MATCHING
# ============================================================================
# For families not explicitly mapped, infer predicates from family number patterns

CAZY_FAMILY_PATTERNS: list[tuple[str, list[str]]] = [
    # Class-level patterns (catch-all for unmapped families)
    (r"^GH\d+", ["glycoside_hydrolase", "carbohydrate_active", "hydrolase"]),
    (r"^GT\d+", ["glycosyltransferase_cazy", "carbohydrate_active", "transferase"]),
    (r"^PL\d+", ["polysaccharide_lyase", "carbohydrate_active", "lyase"]),
    (r"^CE\d+", ["carbohydrate_esterase", "carbohydrate_active", "esterase"]),
    (r"^CBM\d+", ["carbohydrate_binding", "carbohydrate_active"]),
    (r"^AA\d+", ["auxiliary_activity", "carbohydrate_active", "oxidoreductase"]),
]


def get_predicates_for_cazy(family: str) -> list[str]:
    """
    Get predicates for a CAZy family.

    Args:
        family: CAZy family (e.g., "GH5", "GT2", "AA9", "CBM1")

    Returns:
        List of predicate IDs
    """
    predicates = set()

    # Normalize family ID (remove spaces, uppercase)
    family = family.strip().upper()

    # Direct mapping
    if family in CAZY_TO_PREDICATES:
        predicates.update(CAZY_TO_PREDICATES[family])
    else:
        # Pattern matching fallback
        for pattern, preds in CAZY_FAMILY_PATTERNS:
            if re.match(pattern, family):
                predicates.update(preds)
                break

    return sorted(predicates)


def get_cazy_class(family: str) -> Optional[str]:
    """
    Get the CAZy class for a family.

    Args:
        family: CAZy family (e.g., "GH5", "AA10")

    Returns:
        Class name (GH, GT, PL, CE, CBM, AA) or None
    """
    match = re.match(r"^([A-Z]+)", family.upper())
    if match:
        class_prefix = match.group(1)
        if class_prefix in CAZY_CLASS_PREDICATES:
            return class_prefix
    return None


__all__ = [
    "CAZY_CLASS_PREDICATES",
    "CAZY_TO_PREDICATES",
    "CAZY_FAMILY_PATTERNS",
    "get_predicates_for_cazy",
    "get_cazy_class",
]
