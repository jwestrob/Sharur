"""
VOGdb → Predicate mappings.

VOGdb (Virus Orthologous Groups) provides functional annotations for viral proteins.
This module maps VOG functional categories and description patterns to semantic predicates.

Reference: https://vogdb.org/
Data source: vog.annotations.tsv from VOGdb release 232 (September 2025)

E-value thresholds:
- For predicate assignment, we use the same thresholds as other annotation sources
- confident_hit: e-value < 1e-10
- weak_hit: e-value > 1e-5
- Typical VOGdb hits with e-value < 1e-5 are reliable

VOGdb provides ~48,000 VOGs, of which ~13,000 have informative annotations.
"""

from __future__ import annotations

import re
from typing import Optional

# =============================================================================
# VOGdb Functional Category Mappings
# =============================================================================
# VOGdb uses letter codes for broad functional categories.
# These can combine (e.g., XhXs = both host-beneficial and structural).

VOG_CATEGORY_TO_PREDICATES: dict[str, list[str]] = {
    "Xr": ["viral_replication"],
    "Xs": ["viral_structure"],
    "Xh": ["host_beneficial"],  # Function beneficial for host
    "Xp": ["virus_beneficial"],  # Function beneficial for virus
    "Xu": [],  # Unknown - no predicates assigned
}


def predicates_from_vog_category(category: str) -> list[str]:
    """
    Map VOGdb functional category code to predicates.

    Args:
        category: VOGdb category code (e.g., "Xr", "XhXs", "Xu")

    Returns:
        List of predicate names
    """
    predicates = []

    # Handle combined categories (e.g., "XhXs", "XhXpXrXs")
    # Extract individual category codes
    for code in ["Xr", "Xs", "Xh", "Xp"]:
        if code in category:
            predicates.extend(VOG_CATEGORY_TO_PREDICATES.get(code, []))

    return list(set(predicates))


# =============================================================================
# Description Pattern Mappings
# =============================================================================
# These patterns match against the ConsensusFunctionalDescription field.
# Patterns are checked in order; a description can match multiple patterns.

VOG_DESCRIPTION_PATTERNS: list[tuple[re.Pattern, list[str]]] = [
    # -------------------------------------------------------------------------
    # Virion Structure - Capsid/Head
    # -------------------------------------------------------------------------
    (re.compile(r"\bcapsid\b", re.I), ["viral_capsid", "viral_structure"]),
    (re.compile(r"\bmajor\s+head\b", re.I), ["viral_capsid", "viral_structure"]),
    (re.compile(r"\bhead\s+protein\b", re.I), ["viral_capsid", "viral_structure"]),
    (re.compile(r"\bcoat\s+protein\b", re.I), ["viral_capsid", "viral_structure"]),
    (re.compile(r"\bprohead\b", re.I), ["viral_capsid", "viral_structure"]),
    (re.compile(r"\bprocapsid\b", re.I), ["viral_capsid", "viral_structure"]),
    (re.compile(r"\bscaffold\b.*\b(protein|assembly)\b", re.I), ["viral_scaffold", "viral_structure"]),

    # -------------------------------------------------------------------------
    # Virion Structure - Portal
    # -------------------------------------------------------------------------
    (re.compile(r"\bportal\b", re.I), ["viral_portal", "viral_structure"]),
    (re.compile(r"\bconnector\b", re.I), ["viral_portal", "viral_structure"]),

    # -------------------------------------------------------------------------
    # Virion Structure - Tail
    # -------------------------------------------------------------------------
    (re.compile(r"\btail\s+fiber\b", re.I), ["viral_tail_fiber", "viral_tail", "viral_structure"]),
    (re.compile(r"\btail\s+spike\b", re.I), ["viral_tail_fiber", "viral_tail", "viral_structure"]),
    (re.compile(r"\btape\s+measure\b", re.I), ["viral_tape_measure", "viral_tail", "viral_structure"]),
    (re.compile(r"\btail\s+length\b", re.I), ["viral_tape_measure", "viral_tail", "viral_structure"]),
    (re.compile(r"\btail\s+sheath\b", re.I), ["viral_tail_sheath", "viral_tail", "viral_structure"]),
    (re.compile(r"\btail\s+tube\b", re.I), ["viral_tail", "viral_structure"]),
    (re.compile(r"\bbaseplate\b", re.I), ["viral_baseplate", "viral_tail", "viral_structure"]),
    (re.compile(r"\bwedge\b.*\btail\b", re.I), ["viral_baseplate", "viral_tail", "viral_structure"]),
    (re.compile(r"\bhead-tail\s+adaptor\b", re.I), ["viral_tail", "viral_structure"]),
    (re.compile(r"\bhead-tail\s+connector\b", re.I), ["viral_tail", "viral_portal", "viral_structure"]),
    (re.compile(r"\bminor\s+tail\b", re.I), ["viral_tail", "viral_structure"]),
    (re.compile(r"\btail\s+protein\b", re.I), ["viral_tail", "viral_structure"]),
    (re.compile(r"\btail\s+assembly\b", re.I), ["viral_tail", "viral_structure"]),

    # -------------------------------------------------------------------------
    # DNA Packaging
    # -------------------------------------------------------------------------
    (re.compile(r"\bterminase\b", re.I), ["viral_terminase", "viral_packaging"]),
    (re.compile(r"\bpackaging\b", re.I), ["viral_packaging"]),
    (re.compile(r"\bdna\s+packaging\b", re.I), ["viral_packaging"]),
    (re.compile(r"\bmaturation\b.*\bprotease\b", re.I), ["viral_maturation_protease", "protease"]),
    (re.compile(r"\bhead\s+maturation\b", re.I), ["viral_maturation_protease"]),

    # -------------------------------------------------------------------------
    # Lysis
    # -------------------------------------------------------------------------
    (re.compile(r"\bholin\b", re.I), ["holin", "lysis"]),
    (re.compile(r"\bendolysin\b", re.I), ["endolysin", "lysis"]),
    (re.compile(r"\blysin\b", re.I), ["lysin", "lysis"]),
    (re.compile(r"\blysozyme\b", re.I), ["lysozyme", "lysis"]),
    (re.compile(r"\blysis\b", re.I), ["lysis"]),
    (re.compile(r"\bspanin\b", re.I), ["spanin", "lysis"]),
    (re.compile(r"\bantiholin\b", re.I), ["lysis_inhibitor", "lysis"]),

    # -------------------------------------------------------------------------
    # Replication
    # -------------------------------------------------------------------------
    (re.compile(r"\bdna\s+polymerase\b", re.I), ["dna_polymerase", "viral_replication"]),
    (re.compile(r"\brna.+polymerase\b", re.I), ["rna_polymerase", "viral_replication"]),
    (re.compile(r"\breplicase\b", re.I), ["viral_replication"]),
    (re.compile(r"\bhelicase\b", re.I), ["helicase", "viral_replication"]),
    (re.compile(r"\bprimase\b", re.I), ["primase", "viral_replication"]),
    (re.compile(r"\bribonucleotide\s+reductase\b", re.I), ["oxidoreductase", "nucleotide_metabolism"]),
    (re.compile(r"\breplication\s+factor\s+C\b", re.I), ["replication"]),
    (re.compile(r"\bligase\b", re.I), ["ligase"]),
    (re.compile(r"\bsingle.strand.+binding\b", re.I), ["ssb_protein", "viral_replication"]),
    (re.compile(r"\brecombinase\b", re.I), ["recombinase"]),
    (re.compile(r"\bholliday\b", re.I), ["recombinase"]),

    # -------------------------------------------------------------------------
    # Nucleases
    # -------------------------------------------------------------------------
    (re.compile(r"\bendonuclease\b", re.I), ["endonuclease", "nuclease"]),
    (re.compile(r"\bexonuclease\b", re.I), ["exonuclease", "nuclease"]),
    (re.compile(r"\bhoming\s+endonuclease\b", re.I), ["homing_endonuclease", "mobile_element"]),
    (re.compile(r"\bnuclease\b", re.I), ["nuclease"]),
    (re.compile(r"\brnase\b", re.I), ["nuclease"]),

    # -------------------------------------------------------------------------
    # Integration / Lysogeny
    # -------------------------------------------------------------------------
    (re.compile(r"\bintegrase\b", re.I), ["integrase", "lysogeny"]),
    (re.compile(r"\bexcisionase\b", re.I), ["excisionase", "lysogeny"]),
    (re.compile(r"\brepressor\b", re.I), ["repressor"]),
    (re.compile(r"\banti.?repressor\b", re.I), ["anti_repressor"]),
    (re.compile(r"\bcI\s+protein\b", re.I), ["repressor", "lysogeny"]),
    (re.compile(r"\bcro\b", re.I), ["lysogeny"]),

    # -------------------------------------------------------------------------
    # Transcription
    # -------------------------------------------------------------------------
    (re.compile(r"\btranscription\s+factor\b", re.I), ["transcription_factor"]),
    (re.compile(r"\btranscription\s+regulator\b", re.I), ["transcription_factor"]),
    (re.compile(r"\btranscription\s+initiation\s+factor\b", re.I), ["transcription", "transcription_factor", "dna_binding"]),
    (re.compile(r"\bsigma\s+factor\b", re.I), ["sigma_factor", "transcription_factor"]),
    (re.compile(r"\banti.?sigma\b", re.I), ["anti_sigma", "transcription_factor"]),

    # -------------------------------------------------------------------------
    # DNA Modification
    # -------------------------------------------------------------------------
    (re.compile(r"\bmethyltransferase\b", re.I), ["methyltransferase", "dna_modification"]),
    (re.compile(r"\bmethylase\b", re.I), ["methyltransferase", "dna_modification"]),
    (re.compile(r"\bdam\b.*\bmethyl", re.I), ["dam_methylase", "dna_modification"]),
    (re.compile(r"\bdcm\b.*\bmethyl", re.I), ["dcm_methylase", "dna_modification"]),

    # -------------------------------------------------------------------------
    # Host Interaction / Anti-defense
    # -------------------------------------------------------------------------
    (re.compile(r"\banti.?crispr\b", re.I), ["anti_crispr", "anti_defense"]),
    (re.compile(r"\banti.?restriction\b", re.I), ["anti_restriction", "anti_defense"]),
    (re.compile(r"\breceptor\s+binding\b", re.I), ["receptor_binding", "host_attachment"]),
    (re.compile(r"\badsorption\b", re.I), ["host_attachment"]),
    (re.compile(r"\bhost\s+specificity\b", re.I), ["host_attachment"]),
    (re.compile(r"\binjection\b", re.I), ["dna_injection"]),
    (re.compile(r"\bejection\b", re.I), ["dna_injection"]),

    # -------------------------------------------------------------------------
    # Ubiquitin system
    # -------------------------------------------------------------------------
    (re.compile(r"\bE3\s+ubiquitin\s+ligase\b|\bubiquitin[-\s]*protein\s+ligase\b", re.I), ["ubiquitin_ligase", "protein_modification"]),
    (re.compile(r"\bdeubiquitin|ubiquitin\s+protease|ubiquitin\s+hydrolase|ubiquitin-specific\s+protease", re.I), ["deubiquitinase", "protease"]),
    (re.compile(r"\bubiquitin\b", re.I), ["ubiquitin_like"]),

    # -------------------------------------------------------------------------
    # Chromatin / histone
    # -------------------------------------------------------------------------
    (re.compile(r"\bhistone\b", re.I), ["histone"]),
    (re.compile(r"\bchromatin\b", re.I), ["chromatin"]),

    # -------------------------------------------------------------------------
    # Toxin-Antitoxin
    # -------------------------------------------------------------------------
    (re.compile(r"\btoxin\b", re.I), ["toxin"]),
    (re.compile(r"\bantitoxin\b", re.I), ["antitoxin"]),

    # -------------------------------------------------------------------------
    # Enzymes
    # -------------------------------------------------------------------------
    (re.compile(r"\bprotease\b|\bpeptidase\b", re.I), ["protease"]),
    (re.compile(r"\bkinase\b", re.I), ["kinase"]),
    (re.compile(r"\bphosphatase\b", re.I), ["phosphatase"]),
    (re.compile(r"\btriphosphatase\b", re.I), ["phosphatase", "hydrolase"]),
    (re.compile(r"\bgtpase\b", re.I), ["gtpase"]),
    (re.compile(r"\batpase\b", re.I), ["atpase"]),
    (re.compile(r"\bftsH\b", re.I), ["protease", "atpase", "membrane"]),
    (re.compile(r"\bcell\s+division\s+protein\s+48\b|\bcdc48\b|\bp97\b", re.I), ["aaa_domain", "atpase"]),
    (re.compile(r"\bdeaminase\b", re.I), ["deaminase"]),
    (re.compile(r"\bradical\s+sam\b", re.I), ["radical_sam", "iron_sulfur"]),

    # -------------------------------------------------------------------------
    # Structural features
    # -------------------------------------------------------------------------
    (re.compile(r"\bstructural\s+protein\b", re.I), ["viral_structure"]),
    (re.compile(r"\bvirion\b", re.I), ["viral_structure"]),
    (re.compile(r"\bleucin(e)?\s+rich\s+repeat\b", re.I), ["lrr_repeat", "repeat_domain"]),

    # -------------------------------------------------------------------------
    # Mobile elements
    # -------------------------------------------------------------------------
    (re.compile(r"\btransposase\b", re.I), ["transposase", "mobile_element"]),
]


def predicates_from_vog_description(description: str) -> list[str]:
    """
    Extract predicates from VOGdb description using pattern matching.

    Args:
        description: VOGdb ConsensusFunctionalDescription field

    Returns:
        List of predicate names
    """
    predicates = []

    # Skip hypothetical proteins
    if "hypothetical" in description.lower():
        return []

    for pattern, preds in VOG_DESCRIPTION_PATTERNS:
        if pattern.search(description):
            predicates.extend(preds)

    return list(set(predicates))


def predicates_from_vog(
    vog_id: str,
    category: Optional[str] = None,
    description: Optional[str] = None,
) -> list[str]:
    """
    Get all predicates for a VOG entry.

    Args:
        vog_id: VOG identifier (e.g., "VOG12955")
        category: VOGdb functional category (e.g., "Xr", "XhXs")
        description: VOGdb consensus functional description

    Returns:
        List of predicate names (deduplicated)
    """
    predicates = []

    # Add category-based predicates
    if category:
        predicates.extend(predicates_from_vog_category(category))

    # Add description-based predicates
    if description:
        predicates.extend(predicates_from_vog_description(description))

    return list(set(predicates))


# =============================================================================
# Direct VOG → Predicate Mappings (for high-value specific VOGs)
# =============================================================================
# These are curated mappings for specific VOGs that are particularly important
# or where the description pattern matching might miss them.

VOG_DIRECT_MAPPINGS: dict[str, list[str]] = {
    # Anti-CRISPR proteins (from analysis)
    "VOG00563": ["anti_crispr", "anti_defense"],  # AcrF1
    "VOG07654": ["anti_crispr", "anti_defense"],  # AcrIIA1
    "VOG11938": ["anti_crispr", "anti_defense"],  # Type I-E anti-CRISPR

    # Add more curated mappings as discovered during analyses
}


def get_vog_predicates(
    vog_id: str,
    category: Optional[str] = None,
    description: Optional[str] = None,
) -> list[str]:
    """
    Get predicates for a VOG, checking direct mappings first.

    This is the main entry point for VOG predicate lookup.

    Args:
        vog_id: VOG identifier (e.g., "VOG12955")
        category: VOGdb functional category
        description: VOGdb consensus functional description

    Returns:
        List of predicate names
    """
    # Check direct mappings first
    if vog_id in VOG_DIRECT_MAPPINGS:
        return VOG_DIRECT_MAPPINGS[vog_id].copy()

    # Fall back to pattern-based mapping
    return predicates_from_vog(vog_id, category, description)


__all__ = [
    "VOG_CATEGORY_TO_PREDICATES",
    "VOG_DESCRIPTION_PATTERNS",
    "VOG_DIRECT_MAPPINGS",
    "predicates_from_vog_category",
    "predicates_from_vog_description",
    "predicates_from_vog",
    "get_vog_predicates",
]
