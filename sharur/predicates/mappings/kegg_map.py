"""
KEGG ortholog to predicate mappings.

Maps KEGG KO accessions to semantic predicates using:
1. EC number classification (top-level enzyme classes)
2. Direct KO mappings for key pathways
3. Pattern-based matching on definitions
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
    "1.1.1": ["oxidoreductase", "dehydrogenase"],  # With NAD+ or NADP+ as acceptor
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
    "1.14": ["oxidoreductase", "oxygenase", "monooxygenase"],  # Paired oxygenases
    "1.15": ["oxidoreductase", "superoxide_dismutase", "oxidative_stress"],  # Acting on superoxide
    "1.16": ["oxidoreductase", "metal_binding"],  # Oxidizing metal ions
    "1.17": ["oxidoreductase"],  # Acting on CH or CH2
    "1.18": ["oxidoreductase"],  # Acting on iron-sulfur proteins as donors
    "1.97": ["oxidoreductase"],  # Other oxidoreductases

    # EC 2: Transferases
    "2": ["transferase"],
    "2.1": ["transferase"],  # Transferring one-carbon groups
    "2.1.1": ["transferase", "methyltransferase"],  # Methyltransferases
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
# DIRECT KEGG ORTHOLOG MAPPINGS
# ============================================================================
# Key KEGG orthologs mapped to predicates

KEGG_TO_PREDICATES: dict[str, list[str]] = {
    # -------------------------------------------------------------------------
    # HYDROGENASES
    # -------------------------------------------------------------------------
    # Generic / unclassified hydrogenases
    "K00532": ["hydrogenase", "hydrogen_metabolism"],  # E1.12.7.2; ferredoxin hydrogenase
    "K00533": ["hydrogenase", "hydrogen_metabolism"],  # E1.12.7.2L; ferredoxin hydrogenase large subunit
    "K00534": ["hydrogenase", "hydrogen_metabolism"],  # E1.12.7.2S; ferredoxin hydrogenase small subunit
    "K00437": ["hydrogenase", "nife_hydrogenase", "hydrogen_metabolism"],  # hydB; [NiFe] hydrogenase large subunit
    "K18008": ["hydrogenase", "nife_hydrogenase", "hydrogen_metabolism"],  # hydA; [NiFe] hydrogenase small subunit

    # Group 1 - respiratory uptake [NiFe]-hydrogenases
    "K06281": ["hydrogenase", "nife_hydrogenase", "hydrogen_metabolism", "nife_group1", "uptake_hydrogenase"],  # hyaB; hydrogenase large subunit
    "K06282": ["hydrogenase", "nife_hydrogenase", "hydrogen_metabolism", "nife_group1", "uptake_hydrogenase"],  # hyaA; hydrogenase small subunit
    "K05922": ["hydrogenase", "nife_hydrogenase", "hydrogen_metabolism", "nife_group1", "uptake_hydrogenase"],  # hydB; quinone-reactive Ni/Fe-hydrogenase large subunit
    "K05927": ["hydrogenase", "nife_hydrogenase", "hydrogen_metabolism", "nife_group1", "uptake_hydrogenase"],  # hydA; quinone-reactive Ni/Fe-hydrogenase small subunit

    # Group 3 - cofactor-coupled bidirectional [NiFe]-hydrogenases (Hox 3d, Frh 3a, Mvh 3c)
    "K00436": ["hydrogenase", "nife_hydrogenase", "hydrogen_metabolism", "nife_group3", "bidirectional_hydrogenase", "nad_coupled"],  # hoxH; NAD-reducing hydrogenase large subunit
    "K00440": ["hydrogenase", "nife_hydrogenase", "hydrogen_metabolism", "nife_group3", "f420_dependent", "f420_reducing"],  # frhA; coenzyme F420 hydrogenase subunit alpha
    "K00441": ["hydrogenase", "nife_hydrogenase", "hydrogen_metabolism", "nife_group3", "f420_dependent", "f420_reducing"],  # frhB; coenzyme F420 hydrogenase subunit beta
    "K00443": ["hydrogenase", "nife_hydrogenase", "hydrogen_metabolism", "nife_group3", "f420_dependent", "f420_reducing"],  # frhG; coenzyme F420 hydrogenase subunit gamma
    "K14126": ["hydrogenase", "nife_hydrogenase", "hydrogen_metabolism", "nife_group3", "heterodisulfide_reductase_linked"],  # mvhA; F420-non-reducing hydrogenase large subunit
    "K14127": ["hydrogenase", "nife_hydrogenase", "hydrogen_metabolism", "nife_group3", "heterodisulfide_reductase_linked"],  # mvhD; F420-non-reducing hydrogenase iron-sulfur subunit
    "K14128": ["hydrogenase", "nife_hydrogenase", "hydrogen_metabolism", "nife_group3", "heterodisulfide_reductase_linked"],  # mvhG; F420-non-reducing hydrogenase small subunit

    # Group 4 - H2-evolving [NiFe]-hydrogenase complexes (Hyc 4a, Ech 4e, Mbh)
    "K15827": ["hydrogenase", "nife_hydrogenase", "hydrogen_metabolism", "nife_group4", "h2_evolving", "formate_coupled"],  # hycB; formate hydrogenlyase subunit 2
    "K15828": ["hydrogenase", "nife_hydrogenase", "hydrogen_metabolism", "nife_group4", "h2_evolving", "formate_coupled", "membrane"],  # hycC; formate hydrogenlyase subunit 3
    "K15829": ["hydrogenase", "nife_hydrogenase", "hydrogen_metabolism", "nife_group4", "h2_evolving", "formate_coupled", "membrane"],  # hycD; formate hydrogenlyase subunit 4
    "K15830": ["hydrogenase", "nife_hydrogenase", "hydrogen_metabolism", "nife_group4", "h2_evolving", "formate_coupled"],  # hycE; formate hydrogenlyase subunit 5
    "K15831": ["hydrogenase", "nife_hydrogenase", "hydrogen_metabolism", "nife_group4", "h2_evolving", "formate_coupled"],  # hycF; formate hydrogenlyase subunit 6
    "K15832": ["hydrogenase", "nife_hydrogenase", "hydrogen_metabolism", "nife_group4", "h2_evolving", "formate_coupled"],  # hycG; formate hydrogenlyase subunit 7
    "K14086": ["hydrogenase", "nife_hydrogenase", "hydrogen_metabolism", "nife_group4", "ech_hydrogenase", "membrane"],  # echA; ech hydrogenase subunit A
    "K14087": ["hydrogenase", "nife_hydrogenase", "hydrogen_metabolism", "nife_group4", "ech_hydrogenase", "membrane"],  # echB; ech hydrogenase subunit B
    "K14088": ["hydrogenase", "nife_hydrogenase", "hydrogen_metabolism", "nife_group4", "ech_hydrogenase"],  # echC; ech hydrogenase subunit C
    "K14089": ["hydrogenase", "nife_hydrogenase", "hydrogen_metabolism", "nife_group4", "ech_hydrogenase"],  # echD; ech hydrogenase subunit D
    "K14090": ["hydrogenase", "nife_hydrogenase", "hydrogen_metabolism", "nife_group4", "ech_hydrogenase"],  # echE; ech hydrogenase subunit E
    "K14091": ["hydrogenase", "nife_hydrogenase", "hydrogen_metabolism", "nife_group4", "ech_hydrogenase"],  # echF; ech hydrogenase subunit F
    "K18016": ["hydrogenase", "nife_hydrogenase", "hydrogen_metabolism", "nife_group4", "mbh_hydrogenase"],  # mbhL; membrane-bound hydrogenase subunit alpha
    "K18017": ["hydrogenase", "nife_hydrogenase", "hydrogen_metabolism", "nife_group4", "mbh_hydrogenase"],  # mbhK; membrane-bound hydrogenase subunit beta
    "K18023": ["hydrogenase", "nife_hydrogenase", "hydrogen_metabolism", "nife_group4", "mbh_hydrogenase", "membrane"],  # mbhJ; membrane-bound hydrogenase subunit mbhJ

    # Hydrogenase maturation
    "K04651": ["hydrogenase_maturation", "nickel_binding"],  # hypA; hydrogenase nickel incorporation protein HypA/HybF
    "K04652": ["hydrogenase_maturation"],  # hypB; hydrogenase nickel incorporation protein HypB
    "K04653": ["hydrogenase_maturation"],  # hypC; hydrogenase expression/formation protein HypC
    "K04654": ["hydrogenase_maturation"],  # hypD; hydrogenase expression/formation protein HypD
    "K04655": ["hydrogenase_maturation"],  # hypE; hydrogenase expression/formation protein HypE
    "K04656": ["hydrogenase_maturation"],  # hypF; hydrogenase maturation protein HypF

    # -------------------------------------------------------------------------
    # NITROGEN FIXATION
    # -------------------------------------------------------------------------
    # NOTE: Only nifHDK (core nitrogenase subunits) should trigger nitrogen_fixation.
    # Accessory genes (nifBNE) are involved in cofactor biosynthesis but organisms
    # can have these genes without actually fixing nitrogen (genes repurposed).
    "K02586": ["nitrogenase", "nitrogen_fixation", "molybdenum_binding", "iron_sulfur"],  # nifD; nitrogenase molybdenum-iron protein alpha chain
    "K02588": ["nitrogenase", "nitrogen_fixation", "iron_sulfur", "atp_binding"],  # nifH; nitrogenase iron protein NifH
    "K02591": ["nitrogenase", "nitrogen_fixation", "molybdenum_binding", "iron_sulfur"],  # nifK; nitrogenase molybdenum-iron protein beta chain
    "K02592": ["nitrogenase_maturation"],  # nifN; nitrogenase molybdenum-iron protein NifN
    # Accessory genes - do NOT indicate nitrogen fixation without nifHDK
    "K02587": ["nitrogenase_maturation"],  # nifE; nitrogenase molybdenum-cofactor synthesis protein NifE
    "K02585": ["nitrogenase_maturation", "iron_sulfur"],  # nifB; nitrogen fixation protein NifB
    "K02594": ["nitrogenase_maturation"],  # nifV; homocitrate synthase NifV
    "K02593": ["nitrogen_fixation"],  # nifT; nitrogen fixation protein NifT

    # -------------------------------------------------------------------------
    # DENITRIFICATION
    # -------------------------------------------------------------------------
    "K00370": ["denitrification", "nitrate_reduction", "molybdenum_binding"],  # narG; nitrate reductase / nitrite oxidoreductase, alpha subunit
    "K00371": ["denitrification", "nitrate_reduction", "iron_sulfur"],  # narH; nitrate reductase / nitrite oxidoreductase, beta subunit
    "K00374": ["denitrification", "nitrate_reduction"],  # narI; nitrate reductase gamma subunit
    "K02567": ["denitrification", "nitrite_reductase"],  # napA; nitrate reductase (cytochrome)
    "K00368": ["denitrification", "nitrite_reductase", "heme_binding"],  # nirK; nitrite reductase (NO-forming)
    "K15864": ["denitrification", "nitrite_reductase", "heme_binding"],  # nirS; nitrite reductase (NO-forming) / hydroxylamine reductase
    "K00376": ["denitrification", "nitrous_oxide_reductase", "copper_binding"],  # nosZ; nitrous-oxide reductase
    "K04561": ["denitrification", "nitric_oxide_reductase"],  # norB; nitric oxide reductase subunit B

    # -------------------------------------------------------------------------
    # SULFUR METABOLISM
    # -------------------------------------------------------------------------
    "K11180": ["sulfate_reduction", "sulfur_metabolism"],  # dsrA; dissimilatory sulfite reductase alpha subunit
    "K11181": ["sulfate_reduction", "sulfur_metabolism"],  # dsrB; dissimilatory sulfite reductase beta subunit
    "K00394": ["sulfate_reduction", "sulfur_metabolism", "atp_binding"],  # aprA; adenylylsulfate reductase, subunit A
    "K00395": ["sulfate_reduction", "sulfur_metabolism"],  # aprB; adenylylsulfate reductase, subunit B
    "K00958": ["sulfur_assimilation", "sulfur_metabolism", "atp_binding"],  # sat; sulfate adenylyltransferase
    "K00860": ["sulfur_assimilation", "sulfur_metabolism"],  # cysC; adenylylsulfate kinase
    "K17218": ["sulfur_oxidation", "sulfur_metabolism"],  # sqr; sulfide:quinone oxidoreductase
    "K17222": ["sulfur_oxidation", "sulfur_metabolism"],  # soxA; L-cysteine S-thiosulfotransferase
    "K17223": ["sulfur_oxidation", "sulfur_metabolism"],  # soxX; L-cysteine S-thiosulfotransferase
    "K17224": ["sulfur_oxidation", "sulfur_metabolism"],  # soxB; S-sulfosulfanyl-L-cysteine sulfohydrolase
    "K17225": ["sulfur_oxidation", "sulfur_metabolism"],  # soxC; sulfane dehydrogenase subunit SoxC
    "K17226": ["sulfur_oxidation", "sulfur_metabolism"],  # soxY; sulfur-oxidizing protein SoxY
    "K17227": ["sulfur_oxidation", "sulfur_metabolism"],  # soxZ; sulfur-oxidizing protein SoxZ

    # -------------------------------------------------------------------------
    # METHANOGENESIS
    # -------------------------------------------------------------------------
    # NOTE: Only MCR (methyl-CoM reductase) subunits are DEFINITIVE for methanogenesis.
    # Other enzymes (HDR, Fwd, Mtr, Mer) are shared with non-methanogenic archaea.
    "K00399": ["methanogenesis", "mcr_complex", "one_carbon_metabolism", "nickel_binding"],  # mcrA; methyl-coenzyme M reductase alpha subunit
    "K00401": ["methanogenesis", "mcr_complex", "one_carbon_metabolism"],  # mcrB; methyl-coenzyme M reductase beta subunit
    "K00402": ["methanogenesis", "mcr_complex", "one_carbon_metabolism"],  # mcrG; methyl-coenzyme M reductase gamma subunit
    # HDR, Fwd, Mtr, Mer - shared with non-methanogens (e.g., Altiarchaeota)
    "K03388": ["archaeal_one_carbon", "one_carbon_metabolism", "iron_sulfur"],  # hdrA2; heterodisulfide reductase subunit A2
    "K03389": ["archaeal_one_carbon", "one_carbon_metabolism"],  # hdrB2; heterodisulfide reductase subunit B2
    "K03390": ["archaeal_one_carbon", "one_carbon_metabolism"],  # hdrC2; heterodisulfide reductase subunit C2
    "K00200": ["archaeal_one_carbon", "one_carbon_metabolism", "wood_ljungdahl"],  # fwdA; formylmethanofuran dehydrogenase subunit A
    "K00201": ["archaeal_one_carbon", "one_carbon_metabolism", "wood_ljungdahl"],  # fwdB; formylmethanofuran dehydrogenase subunit B
    "K00577": ["archaeal_one_carbon", "one_carbon_metabolism", "cobalamin_binding"],  # mtrA; tetrahydromethanopterin S-methyltransferase subunit A
    "K00320": ["archaeal_one_carbon", "one_carbon_metabolism", "fad_binding"],  # mer; 5,10-methylenetetrahydromethanopterin reductase

    # Methanotrophy
    "K10944": ["methanotrophy", "one_carbon_metabolism", "copper_binding"],  # pmoA; methane monooxygenase subunit A
    "K10945": ["methanotrophy", "one_carbon_metabolism"],  # pmoB-amoB; methane/ammonia monooxygenase subunit B
    "K10946": ["methanotrophy", "one_carbon_metabolism"],  # pmoC-amoC; methane/ammonia monooxygenase subunit C
    "K16157": ["methanotrophy", "one_carbon_metabolism"],  # mmoX; methane monooxygenase component A alpha chain
    "K16158": ["methanotrophy", "one_carbon_metabolism"],  # mmoY; methane monooxygenase component A beta chain
    "K16159": ["methanotrophy", "one_carbon_metabolism"],  # mmoZ; methane monooxygenase component A gamma chain

    # -------------------------------------------------------------------------
    # WOOD-LJUNGDAHL PATHWAY
    # -------------------------------------------------------------------------
    "K00192": ["wood_ljungdahl", "one_carbon_metabolism", "nickel_binding"],  # cdhA; anaerobic carbon-monoxide dehydrogenase, CODH/ACS complex subunit a...
    "K00193": ["wood_ljungdahl", "one_carbon_metabolism"],  # cdhC; acetyl-CoA decarbonylase/synthase, CODH/ACS complex subunit beta
    "K00194": ["wood_ljungdahl", "one_carbon_metabolism"],  # cdhD; acetyl-CoA decarbonylase/synthase, CODH/ACS complex subunit delta
    "K00197": ["wood_ljungdahl", "one_carbon_metabolism", "iron_sulfur"],  # cdhE; acetyl-CoA decarbonylase/synthase, CODH/ACS complex subunit gamma
    "K00198": ["wood_ljungdahl", "co_oxidation", "one_carbon_metabolism"],  # cooS; anaerobic carbon-monoxide dehydrogenase catalytic subunit

    # -------------------------------------------------------------------------
    # PHOTOSYNTHESIS
    # -------------------------------------------------------------------------
    "K02703": ["photosynthesis", "photosystem_ii"],  # psbA; photosystem II P680 reaction center D1 protein
    "K02706": ["photosynthesis", "photosystem_ii"],  # psbD; photosystem II P680 reaction center D2 protein
    "K02689": ["photosynthesis", "photosystem_i"],  # psaA; photosystem I P700 chlorophyll a apoprotein A1
    "K02690": ["photosynthesis", "photosystem_i"],  # psaB; photosystem I P700 chlorophyll a apoprotein A2
    "K02634": ["photosynthesis", "electron_transport"],  # petA; apocytochrome f
    "K02635": ["photosynthesis", "electron_transport"],  # petB; cytochrome b6
    "K02636": ["photosynthesis", "electron_transport", "iron_sulfur"],  # petC; cytochrome b6-f complex iron-sulfur subunit
    # NOTE: RuBisCO alone does NOT mean Calvin cycle - RuBisCO-like proteins (RLPs) exist.
    # Only assign rubisco and carbon_fixation, NOT calvin_cycle (requires PRK confirmation).
    "K01601": ["rubisco", "carbon_fixation"],  # rbcL; ribulose-bisphosphate carboxylase large chain
    "K01602": ["rubisco", "carbon_fixation"],  # rbcS; ribulose-bisphosphate carboxylase small chain

    # -------------------------------------------------------------------------
    # ATP SYNTHESIS
    # -------------------------------------------------------------------------
    "K02111": ["atp_synthesis", "energy_metabolism", "atp_binding"],  # ATPF1A; F-type H+/Na+-transporting ATPase subunit alpha
    "K02112": ["atp_synthesis", "energy_metabolism", "atp_binding"],  # ATPF1B; F-type H+/Na+-transporting ATPase subunit beta
    "K02113": ["atp_synthesis", "energy_metabolism"],  # ATPF1D; F-type H+-transporting ATPase subunit delta
    "K02109": ["atp_synthesis", "energy_metabolism", "membrane"],  # ATPF0B; F-type H+-transporting ATPase subunit b
    "K02108": ["atp_synthesis", "energy_metabolism", "membrane"],  # ATPF0A; F-type H+-transporting ATPase subunit a
    "K02110": ["atp_synthesis", "energy_metabolism", "membrane"],  # ATPF0C; F-type H+-transporting ATPase subunit c
    "K02115": ["atp_synthesis", "energy_metabolism"],  # ATPF1G; F-type H+-transporting ATPase subunit gamma

    # -------------------------------------------------------------------------
    # CENTRAL CARBON METABOLISM
    # -------------------------------------------------------------------------
    # Glycolysis
    "K00844": ["glycolysis", "central_metabolism", "kinase"],  # HK; hexokinase
    "K01810": ["glycolysis", "central_metabolism", "isomerase"],  # GPI; glucose-6-phosphate isomerase
    "K00850": ["glycolysis", "central_metabolism", "kinase"],  # pfkA; 6-phosphofructokinase 1
    "K01623": ["glycolysis", "central_metabolism", "lyase"],  # ALDO; fructose-bisphosphate aldolase, class I
    "K01624": ["glycolysis", "central_metabolism", "lyase"],  # FBA; fructose-bisphosphate aldolase, class II
    "K00134": ["glycolysis", "central_metabolism", "dehydrogenase", "nad_binding"],  # GAPDH; glyceraldehyde 3-phosphate dehydrogenase (phosphorylating)
    "K00927": ["glycolysis", "central_metabolism", "kinase"],  # PGK; phosphoglycerate kinase
    "K01689": ["glycolysis", "central_metabolism"],  # ENO1_2_3; enolase 1/2/3
    "K00873": ["glycolysis", "central_metabolism", "kinase"],  # PK; pyruvate kinase

    # TCA cycle
    "K01647": ["tca_cycle", "central_metabolism"],  # CS; citrate synthase
    "K01681": ["tca_cycle", "central_metabolism"],  # ACO; aconitate hydratase
    "K00031": ["tca_cycle", "central_metabolism", "dehydrogenase", "nad_binding"],  # IDH1; isocitrate dehydrogenase
    "K00164": ["tca_cycle", "central_metabolism", "dehydrogenase"],  # OGDH; 2-oxoglutarate dehydrogenase E1 component
    "K01902": ["tca_cycle", "central_metabolism", "ligase"],  # sucD; succinyl-CoA synthetase alpha subunit
    "K00239": ["tca_cycle", "central_metabolism", "dehydrogenase", "fad_binding"],  # sdhA; succinate dehydrogenase flavoprotein subunit
    "K01676": ["tca_cycle", "central_metabolism"],  # E4.2.1.2A; fumarate hydratase, class I
    "K00024": ["tca_cycle", "central_metabolism", "dehydrogenase", "nad_binding"],  # mdh; malate dehydrogenase

    # -------------------------------------------------------------------------
    # TRANSPORTERS
    # -------------------------------------------------------------------------
    "K02003": ["transporter", "abc_transporter", "abc_atpase", "atp_binding"],  # ABC.CD.A; putative ABC transport system ATP-binding protein
    "K02004": ["transporter", "abc_transporter", "abc_permease"],  # ABC.CD.P; putative ABC transport system permease protein
    "K02020": ["transporter", "abc_transporter", "molybdenum_binding"],  # modA; molybdate transport system substrate-binding protein
    "K02035": ["transporter", "abc_transporter", "peptide_transporter"],  # ABC.PE.S; peptide/nickel transport system substrate-binding protein
    "K02036": ["transporter", "abc_transporter", "phosphate_transporter", "atp_binding"],  # pstB; phosphate transport system ATP-binding protein
    "K02037": ["transporter", "abc_transporter", "phosphate_transporter"],  # pstC; phosphate transport system permease protein
    "K02038": ["transporter", "abc_transporter", "phosphate_transporter"],  # pstA; phosphate transport system permease protein
    "K02040": ["transporter", "abc_transporter", "phosphate_transporter"],  # pstS; phosphate transport system substrate-binding protein
    "K02041": ["transporter", "abc_transporter", "atp_binding"],  # phnC; phosphonate transport system ATP-binding protein
    "K02046": ["transporter", "abc_transporter", "sulfate_transporter"],  # cysU; sulfate/thiosulfate transport system permease protein
    "K02055": ["transporter", "abc_transporter"],  # ABC.SP.S; putative spermidine/putrescine transport system substrate-binding p...
    "K15580": ["transporter", "abc_transporter", "peptide_transporter"],  # oppA; oligopeptide transport system substrate-binding protein
    "K23163": ["transporter", "abc_transporter", "sulfate_transporter"],  # sbp; sulfate/thiosulfate transport system substrate-binding protein
    "K02048": ["transporter", "abc_transporter", "sulfate_transporter"],  # cysP; sulfate/thiosulfate transport system substrate-binding protein
    "K10542": ["transporter", "abc_transporter", "sugar_transporter", "atp_binding"],  # mglA; methyl-galactoside transport system ATP-binding protein
    "K03088": ["sigma_factor", "regulator"],  # rpoE; RNA polymerase sigma-70 factor, ECF subfamily

    # -------------------------------------------------------------------------
    # TWO-COMPONENT SYSTEMS
    # -------------------------------------------------------------------------
    "K07638": ["two_component", "sensor_kinase", "regulator"],  # envZ; two-component system, OmpR family, osmolarity sensor histidine kina...
    "K07659": ["two_component", "response_regulator", "regulator"],  # ompR; two-component system, OmpR family, phosphate regulon response regul...
    "K07648": ["two_component", "sensor_kinase", "regulator"],  # arcB; two-component system, OmpR family, aerobic respiration control sens...
    "K07657": ["two_component", "response_regulator", "regulator"],  # phoB; two-component system, OmpR family, phosphate regulon response regul...
    "K07678": ["two_component", "sensor_kinase", "regulator"],  # barA; two-component system, NarL family, sensor histidine kinase BarA
    "K07636": ["two_component", "sensor_kinase", "regulator"],  # phoR; two-component system, OmpR family, phosphate regulon sensor histidi...
    "K07673": ["two_component", "sensor_kinase", "regulator"],  # narX; two-component system, NarL family, nitrate/nitrite sensor histidine...
    "K07684": ["two_component", "response_regulator", "regulator", "narl_family"],  # narL; two-component system, NarL family, nitrate/nitrite response regulat...
    "K02483": ["two_component", "response_regulator", "regulator"],  # K02483; two-component system, OmpR family, response regulator
    "K03407": ["two_component", "sensor_kinase", "chemotaxis"],  # cheA; two-component system, chemotaxis family, sensor kinase CheA
    "K03413": ["two_component", "response_regulator", "chemotaxis"],  # cheY; two-component system, chemotaxis family, chemotaxis protein CheY

    # -------------------------------------------------------------------------
    # SECRETION SYSTEMS
    # -------------------------------------------------------------------------
    "K03070": ["secretion_system", "sec_pathway", "translocase"],  # secA; preprotein translocase subunit SecA
    "K03071": ["secretion_system", "sec_pathway", "chaperone"],  # secB; preprotein translocase subunit SecB
    "K03072": ["secretion_system", "sec_pathway", "membrane"],  # secD; preprotein translocase subunit SecD
    "K03073": ["secretion_system", "sec_pathway", "membrane"],  # secE; preprotein translocase subunit SecE
    "K03076": ["secretion_system", "sec_pathway"],  # secY; preprotein translocase subunit SecY
    "K03075": ["secretion_system", "sec_pathway", "membrane"],  # secG; preprotein translocase subunit SecG
    "K03116": ["secretion_system", "tat_pathway"],  # tatA; sec-independent protein translocase protein TatA
    "K03117": ["secretion_system", "tat_pathway"],  # tatB; sec-independent protein translocase protein TatB
    "K03118": ["secretion_system", "tat_pathway", "membrane"],  # tatC; sec-independent protein translocase protein TatC
    "K03205": ["secretion_system", "type_iv_secretion", "atp_binding"],  # virD4; type IV secretion system protein VirD4
    "K03219": ["secretion_system", "type_iii_secretion"],  # yscC; type III secretion protein C
    "K03224": ["secretion_system", "type_iii_secretion", "atp_binding"],  # yscN; ATP synthase in type III secretion protein N
    "K03230": ["secretion_system", "type_iii_secretion", "membrane"],  # yscV; type III secretion protein V
    "K03194": ["secretion_system", "type_iv_secretion"],  # virB1; type IV secretion system protein VirB1
    "K03195": ["secretion_system", "type_iv_secretion", "membrane"],  # virB10; type IV secretion system protein VirB10
    "K11901": ["secretion_system", "type_vi_secretion"],  # impB; type VI secretion system protein ImpB
    "K11900": ["secretion_system", "type_vi_secretion"],  # impC; type VI secretion system protein ImpC

    # -------------------------------------------------------------------------
    # MOBILE ELEMENTS / DEFENSE
    # -------------------------------------------------------------------------
    "K07481": ["transposase", "mobile_element"],  # K07481; transposase, IS5 family
    "K01356": ["regulator", "repressor", "dna_repair"],  # lexA; repressor LexA
    "K03529": ["chromosome_partitioning", "atp_binding"],  # smc; chromosome segregation protein
    "K09951": ["crispr_associated", "defense_system"],  # cas2; CRISPR-associated protein Cas2
    "K07012": ["crispr_associated", "cas_nuclease", "defense_system"],  # cas3; CRISPR-associated endonuclease/helicase Cas3
    "K19086": ["crispr_associated", "defense_system"],  # csa4; CRISPR-associated protein Csa4
    "K19087": ["crispr_associated", "defense_system"],  # csa5; CRISPR-associated protein Csa5
    "K15342": ["crispr_associated", "defense_system"],  # cas1; CRISP-associated protein Cas1
    "K09952": ["crispr_associated", "cas_nuclease", "defense_system"],  # csn1; CRISPR-associated endonuclease Csn1

    # -------------------------------------------------------------------------
    # CELL DIVISION
    # -------------------------------------------------------------------------
    "K03531": ["cell_division", "ftsz", "gtp_binding"],  # ftsZ; cell division protein FtsZ
    "K03587": ["cell_division", "divisome"],  # ftsI; cell division protein FtsI (penicillin-binding protein 3)
    "K03589": ["cell_division", "divisome"],  # ftsQ; cell division protein FtsQ
    "K03590": ["cell_division", "divisome"],  # ftsA; cell division protein FtsA
    "K03591": ["cell_division", "divisome"],  # ftsN; cell division protein FtsN
    "K03592": ["protease"],  # pmbA; PmbA protein
    "K03593": ["atp_binding"],  # mrp; ATP-binding protein involved in chromosome partitioning
    "K03466": ["cell_division", "divisome", "atp_binding"],  # ftsK; DNA segregation ATPase FtsK/SpoIIIE, S-DNA-T family
    "K03586": ["cell_division", "divisome"],  # ftsL; cell division protein FtsL
    "K03588": ["cell_division", "divisome", "membrane"],  # ftsW; peptidoglycan glycosyltransferase
    "K03595": ["gtp_binding"],  # era; GTPase
    "K03496": ["cell_division", "chromosome_partitioning"],  # parA; chromosome partitioning protein
    "K03497": ["cell_division", "chromosome_partitioning"],  # parB; ParB family transcriptional regulator, chromosome partitioning protein

    # -------------------------------------------------------------------------
    # DNA REPLICATION
    # -------------------------------------------------------------------------
    "K02337": ["replication", "dna_polymerase"],  # dnaE; DNA polymerase III subunit alpha
    "K02338": ["replication", "dna_polymerase"],  # dnaN; DNA polymerase III subunit beta
    "K02340": ["replication", "dna_polymerase"],  # holA; DNA polymerase III subunit delta
    "K02314": ["replication", "helicase", "atp_binding"],  # dnaB; replicative DNA helicase
    "K02316": ["replication", "primase"],  # dnaG; DNA primase
    "K02335": ["replication", "dna_polymerase"],  # polA; DNA polymerase I
    "K02342": ["replication", "dna_polymerase"],  # dnaQ; DNA polymerase III subunit epsilon
    "K02313": ["replication", "dna_binding", "atp_binding"],  # dnaA; chromosomal replication initiator protein
    "K02469": ["replication", "topoisomerase"],  # gyrA; DNA gyrase subunit A
    "K02470": ["replication", "topoisomerase"],  # gyrB; DNA gyrase subunit B

    # -------------------------------------------------------------------------
    # TRANSLATION
    # -------------------------------------------------------------------------
    "K02863": ["ribosomal_protein", "translation"],  # RP-L1; large subunit ribosomal protein L1
    "K02886": ["ribosomal_protein", "translation"],  # RP-L2; large subunit ribosomal protein L2
    "K02906": ["ribosomal_protein", "translation"],  # RP-L3; large subunit ribosomal protein L3
    "K02945": ["ribosomal_protein", "translation"],  # RP-S1; small subunit ribosomal protein S1
    "K02935": ["ribosomal_protein", "translation"],  # RP-L7; large subunit ribosomal protein L7/L12
    "K02355": ["translation_factor", "translation", "gtp_binding"],  # fusA; elongation factor G
    "K02357": ["translation_factor", "translation"],  # tsf; elongation factor Ts
    "K02358": ["translation_factor", "translation", "gtp_binding"],  # tuf; elongation factor Tu
    "K02518": ["translation_factor", "translation"],  # infA; translation initiation factor IF-1
    "K02519": ["translation_factor", "translation", "gtp_binding"],  # infB; translation initiation factor IF-2
    "K02520": ["translation_factor", "translation"],  # infC; translation initiation factor IF-3
    "K01866": ["trna_synthetase", "translation"],  # YARS; tyrosyl-tRNA synthetase
    "K01867": ["trna_synthetase", "translation"],  # WARS; tryptophanyl-tRNA synthetase
    "K01869": ["trna_synthetase", "translation"],  # LARS; leucyl-tRNA synthetase
    "K01868": ["trna_synthetase", "translation"],  # TARS; threonyl-tRNA synthetase

    # -------------------------------------------------------------------------
    # STRESS RESPONSE
    # -------------------------------------------------------------------------
    "K04043": ["heat_shock", "chaperone", "stress_response", "atp_binding"],  # dnaK; molecular chaperone DnaK
    "K03686": ["heat_shock", "chaperone", "stress_response"],  # dnaJ; molecular chaperone DnaJ
    "K04077": ["heat_shock", "chaperone", "stress_response"],  # groEL; chaperonin GroEL
    "K04078": ["heat_shock", "chaperone", "stress_response"],  # groES; chaperonin GroES
    "K04769": ["regulator"],  # spoVT; AbrB family transcriptional regulator, stage V sporulation protein T
    "K04079": ["heat_shock", "stress_response", "chaperone"],  # HSP90A; molecular chaperone HtpG
    "K03544": ["clp_protease", "chaperone", "atp_binding"],  # clpX; ATP-dependent Clp protease ATP-binding subunit ClpX
    "K03694": ["clp_protease", "chaperone", "atp_binding"],  # clpA; ATP-dependent Clp protease ATP-binding subunit ClpA
    "K03695": ["clp_protease", "chaperone", "atp_binding"],  # clpB; ATP-dependent Clp protease ATP-binding subunit ClpB
    "K01358": ["clp_protease", "protease"],  # clpP; ATP-dependent Clp protease, protease subunit
    "K03671": ["oxidoreductase"],  # TXN; thioredoxin
    "K03386": ["oxidative_stress", "peroxiredoxin"],  # PRDX2_4; peroxiredoxin 2/4
    "K04564": ["oxidative_stress", "superoxide_dismutase", "metal_binding"],  # SOD2; superoxide dismutase, Fe-Mn family
    "K04565": ["oxidative_stress", "superoxide_dismutase", "copper_binding", "zinc_binding"],  # SOD1; superoxide dismutase, Cu-Zn family
    "K03781": ["oxidative_stress", "catalase", "heme_binding"],  # katE; catalase
    "K03782": ["oxidative_stress", "catalase", "heme_binding"],  # katG; catalase-peroxidase

    # -------------------------------------------------------------------------
    # ANTIBIOTIC RESISTANCE
    # -------------------------------------------------------------------------
    "K01467": ["antibiotic_resistance", "beta_lactamase"],  # ampC; beta-lactamase class C
    "K18698": ["antibiotic_resistance", "beta_lactamase"],  # blaTEM; beta-lactamase class A TEM
    "K17836": ["antibiotic_resistance", "beta_lactamase"],  # penP; beta-lactamase class A
    "K00984": ["antibiotic_resistance", "aminoglycoside_resistance"],  # aadA; streptomycin 3"-adenylyltransferase
    "K00663": ["antibiotic_resistance", "aminoglycoside_resistance"],  # aacA; aminoglycoside 6'-N-acetyltransferase
    "K10673": ["antibiotic_resistance", "aminoglycoside_resistance"],  # strA; streptomycin 3"-kinase
    "K18139": ["antibiotic_resistance", "efflux_pump", "multidrug_resistance"],  # oprM; outer membrane protein, multidrug efflux system
    "K18140": ["regulator", "repressor"],  # envR; TetR/AcrR family transcriptional regulator, acrEF/envCD operon repr...
    "K03585": ["antibiotic_resistance", "efflux_pump", "multidrug_resistance"],  # acrA; membrane fusion protein, multidrug efflux system
    "K18138": ["antibiotic_resistance", "efflux_pump", "multidrug_resistance", "membrane"],  # acrB; multidrug efflux pump

    # -------------------------------------------------------------------------
    # ADDITIONAL HIGH-ABUNDANCE KEGG ORTHOLOGS (from dataset analysis)
    # -------------------------------------------------------------------------
    # Glycosyltransferases - LPS/O-antigen biosynthesis
    "K20578": ["glycosyltransferase", "transferase", "carbohydrate_active"],  # tobM2; glycosyltransferase
    "K23102": ["glycosyltransferase", "transferase", "lps_biosynthesis", "carbohydrate_active"],  # wekW; O18-antigen biosynthesis alpha-1,3-galactosyltransferase
    "K24654": ["glycosyltransferase", "transferase", "lps_biosynthesis", "carbohydrate_active"],  # wclR; O3-antigen biosynthesis alpha-1,3-galactosyltransferase
    "K24649": ["glycosyltransferase", "transferase", "lps_biosynthesis", "carbohydrate_active"],  # wbdM; O111-antigen biosynthesis glucosyltransferase
    "K20438": ["glycosyltransferase", "transferase", "carbohydrate_active"],  # valG; validoxylamine A glucosyltransferase
    "K23101": ["glycosyltransferase", "transferase", "lps_biosynthesis", "carbohydrate_active"],  # wekV; O18-antigen biosynthesis glucosyltransferase
    "K21586": ["glycosyltransferase", "transferase", "lps_biosynthesis", "carbohydrate_active"],  # wbgM; O55-antigen biosynthesis alpha-1,3-galactosyltransferase
    "K23229": ["glycosyltransferase", "transferase", "lps_biosynthesis", "carbohydrate_active"],  # wbgP; O55-antigen biosynthesis glycosyltransferase
    "K23105": ["glycosyltransferase", "transferase", "lps_biosynthesis", "carbohydrate_active"],  # wbdL; O111-antigen biosynthesis colitosyltransferase
    "K23104": ["glycosyltransferase", "transferase", "lps_biosynthesis", "carbohydrate_active"],  # wbdH; O111-antigen biosynthesis galactosyltransferase
    "K26073": ["glycosyltransferase", "transferase", "lps_biosynthesis", "carbohydrate_active"],  # wejO; O78-antigen biosynthesis glycosyltransferase
    "K23103": ["glycosyltransferase", "transferase", "lps_biosynthesis", "carbohydrate_active"],  # wekU; O18-antigen biosynthesis glycosyltransferase
    "K23130": ["glycosyltransferase", "transferase", "lps_biosynthesis", "carbohydrate_active"],  # wbdN; O157-antigen biosynthesis beta-1,3-glucosyltransferase
    "K23199": ["glycosyltransferase", "transferase", "lps_biosynthesis", "carbohydrate_active"],  # wbbC; O7-antigen biosynthesis mannosyltransferase
    "K23131": ["glycosyltransferase", "transferase", "lps_biosynthesis", "carbohydrate_active"],  # wbdO; O157-antigen biosynthesis glycosyltransferase
    "K21216": ["glycosyltransferase", "transferase", "carbohydrate_active"],  # ncsC6; glycosyltransferase
    "K21365": ["glycosyltransferase", "transferase", "lps_biosynthesis", "carbohydrate_active"],  # wbnH; O86/O127-antigen biosynthesis alpha-1,3-N-acetylgalactosaminyltrans...
    "K23073": ["glycosyltransferase", "transferase", "lps_biosynthesis", "carbohydrate_active"],  # wbbK; O16-antigen biosynthesis glucosyltransferase
    "K23250": ["glycosyltransferase", "transferase", "lps_biosynthesis", "carbohydrate_active"],  # wbwB; O104-antigen biosynthesis alpha-1,4-galactosyltransferase
    "K23230": ["glycosyltransferase", "transferase", "lps_biosynthesis", "carbohydrate_active"],  # wbgO; O55-antigen biosynthesis beta-1,3-galactosyltransferase
    "K20573": ["glycosyltransferase", "transferase", "carbohydrate_active", "secondary_metabolism"],  # kanE; 2'-deamino-2'-hydroxyneamine 1-alpha-D-kanosaminyltransferase
    "K20589": ["glycosyltransferase", "transferase", "carbohydrate_active", "secondary_metabolism"],  # genM2; paromamine D-xylosyltransferase
    "K23251": ["glycosyltransferase", "transferase", "lps_biosynthesis", "carbohydrate_active"],  # wbwC; O104-antigen biosynthesis beta-1,3-galactosyltransferase

    # ABC transporters
    "K15558": ["transporter", "abc_transporter", "abc_atpase", "atp_binding"],  # ophH; phthalate transport system ATP-binding protein
    "K20386": ["transporter", "abc_transporter", "abc_atpase", "atp_binding"],  # cylB; ATP-binding cassette, subfamily B, bacterial CylB

    # Selenate/chlorate reductases - anaerobic respiration
    "K27875": ["oxidoreductase", "anaerobic_respiration", "selenium_metabolism"],  # srdB; selenate reductase (quinol) subunit B
    "K17051": ["oxidoreductase", "anaerobic_respiration", "selenium_metabolism"],  # serB; selenate/chlorate reductase subunit beta

    # Terminal oxidases - respiratory chain
    "K02274": ["cytochrome", "electron_transport", "respiration", "terminal_oxidase", "cytochrome_c_oxidase", "heme_binding"],  # coxA; cytochrome c oxidase subunit I
    "K02275": ["cytochrome", "electron_transport", "respiration", "terminal_oxidase", "cytochrome_c_oxidase", "copper_binding"],  # coxB; cytochrome c oxidase subunit II
    "K02276": ["cytochrome", "electron_transport", "respiration", "terminal_oxidase", "cytochrome_c_oxidase"],  # coxC; cytochrome c oxidase subunit III
    "K02277": ["cytochrome", "electron_transport", "respiration", "terminal_oxidase", "cytochrome_c_oxidase"],  # coxD; cytochrome c oxidase subunit IV
    "K15862": ["cytochrome", "electron_transport", "respiration", "terminal_oxidase", "cytochrome_c_oxidase"],  # ccoNO; cytochrome c oxidase cbb3-type subunit I/II
    "K00425": ["cytochrome", "electron_transport", "respiration", "terminal_oxidase", "cytochrome_bd_oxidase"],  # cydA; cytochrome bd ubiquinol oxidase subunit I
    "K00426": ["cytochrome", "electron_transport", "respiration", "terminal_oxidase", "cytochrome_bd_oxidase"],  # cydB; cytochrome bd ubiquinol oxidase subunit II

    # Halogenases - secondary metabolism
    "K21256": ["oxidoreductase", "fad_binding", "secondary_metabolism"],  # calO3; flavin-dependent halogenase
    "K12714": ["oxidoreductase", "fad_binding", "secondary_metabolism"],  # clohal; clorobiocin biosynthesis protein Clo-hal

    # Methyltransferases
    "K16414": ["methyltransferase", "transferase", "sam_binding"],  # stiK; methyltransferase
    "K19889": ["methyltransferase", "transferase", "sam_binding", "secondary_metabolism"],  # rebM; demethylrebeccamycin-D-glucose O-methyltransferase
    "K27779": ["methyltransferase", "transferase", "sam_binding", "secondary_metabolism"],  # NNMT; norajmaline N-methyltransferase
    "K14374": ["methyltransferase", "transferase", "sam_binding", "secondary_metabolism"],  # aveD; avermectin B 5-O-methyltransferase
    "K12914": ["methyltransferase", "transferase", "sam_binding"],  # phpK; P-methyltransferase
    "K27844": ["methyltransferase", "transferase", "sam_binding", "secondary_metabolism"],  # ANMT; ajmaline N-methyltransferase
    "K20592": ["methyltransferase", "transferase", "sam_binding"],  # genN; SAM-dependent 3''-N-methyltransferase
    "K16403": ["methyltransferase", "transferase", "sam_binding"],  # sorM; O-methyltransferase
    "K12902": ["methyltransferase", "transferase", "sam_binding"],  # fom3; phosphonoacetaldehyde methylase
    "K21896": ["methyltransferase", "transferase", "sam_binding", "secondary_metabolism"],  # NMT; 3-hydroxy-16-methoxy-2,3-dihydrotabersonine N-methyltransferase
    "K20594": ["methyltransferase", "transferase", "sam_binding", "secondary_metabolism"],  # genK; gentamicin X2 methyltransferase
    "K12633": ["methyltransferase", "transferase", "sam_binding"],  # pur5; N-methyl-transferase
    "K21541": ["methyltransferase", "transferase", "sam_binding", "terpene_synthesis"],  # TMT-1_2; squalene methyltransferase
    "K21542": ["methyltransferase", "transferase", "sam_binding", "terpene_synthesis"],  # TMT-3; botryococcene C-methyltransferase
    "K23153": ["methyltransferase", "transferase", "sam_binding", "secondary_metabolism"],  # ncmP; nocamycin O-methyltransferase
    "K21458": ["methyltransferase", "transferase", "sam_binding", "secondary_metabolism"],  # E2.1.1.300; pavine N-methyltransferase
    "K12705": ["methyltransferase", "transferase", "sam_binding", "secondary_metabolism"],  # novO; 8-demethylnovobiocic acid C8-methyltransferase
    "K00571": ["methyltransferase", "dna_methylase", "sam_binding"],  # E2.1.1.72; site-specific DNA-methyltransferase (adenine-specific)
    "K00558": ["methyltransferase", "dna_methylase", "sam_binding"],  # DNMT1; DNA (cytosine-5)-methyltransferase 1
    "K06223": ["methyltransferase", "dna_methylase", "sam_binding"],  # dam; DNA adenine methylase

    # Acetyltransferases
    "K27979": ["acetyltransferase", "transferase", "secondary_metabolism"],  # sbzA; altemicidin L-isoleucyltransferase
    "K20681": ["acetyltransferase", "transferase", "carbohydrate_active"],  # fdtC; dTDP-3-amino-3,6-dideoxy-alpha-D-galactopyranose 3-N-acetyltransferase
    "K24251": ["acetyltransferase", "transferase", "carbohydrate_active"],  # qdtC; dTDP-D-Quip3N acetylase
    "K17939": ["acetyltransferase", "transferase", "carbohydrate_active"],  # perB; GDP-perosamine N-acetyltransferase
    "K19313": ["acetyltransferase", "transferase", "secondary_metabolism"],  # apmA; aminocyclitol acetyltransferase
    "K20680": ["acetyltransferase", "transferase", "carbohydrate_active"],  # fdtB; dTDP-3-amino-3,6-dideoxy-alpha-D-galactopyranose transaminase
    "K26075": ["acetyltransferase", "transferase", "lps_biosynthesis"],  # wckD; O104-antigen biosynthesis O-acyltransferase

    # Dehydrogenases and reductases - carbohydrate metabolism
    "K13306": ["oxidoreductase", "dehydrogenase", "nad_binding", "carbohydrate_active"],  # fcd; dTDP-4-dehydro-6-deoxyglucose reductase
    "K18836": ["oxidoreductase", "reductase"],  # vcaE; reductase VcaE
    "K24122": ["oxidoreductase", "dehydrogenase", "terpene_synthesis"],  # ctmAB; limonene dehydrogenase
    "K15858": ["oxidoreductase", "reductase", "carbohydrate_active"],  # ascF; CDP-3, 6-dideoxy-D-glycero-L-glycero-4-hexulose-4-reductase
    "K21271": ["oxidoreductase", "dehydrogenase", "secondary_metabolism"],  # auaH; aurachin B dehydrogenase
    "K19652": ["oxidoreductase", "dehydrogenase", "nad_binding", "carbohydrate_active"],  # tal; dTDP-6-deoxy-L-talose 4-dehydrogenase [NAD(P)+]
    "K19859": ["oxidoreductase", "reductase", "carbohydrate_active"],  # ang14; dTDP-4-keto-6-deoxyhexose reductase
    "K19180": ["oxidoreductase", "dehydrogenase", "nad_binding", "carbohydrate_active"],  # tll; dTDP-6-deoxy-L-talose 4-dehydrogenase (NAD+)
    "K13319": ["oxidoreductase", "reductase", "secondary_metabolism"],  # oleU; 4-ketoreductase
    "K19857": ["oxidoreductase", "reductase", "carbohydrate_active"],  # aveBIV; dTDP-4-keto-6-deoxy-L-hexose 4-reductase
    "K22098": ["oxidoreductase", "dehydrogenase", "secondary_metabolism"],  # SDR1; noscapine synthase
    "K20152": ["oxidoreductase", "reductase", "secondary_metabolism"],  # kanK; 2'-dehydrokanamycin reductase
    "K18535": ["oxidoreductase"],  # thnA; putative oxidoreductase
    "K14373": ["oxidoreductase", "reductase"],  # aveF; C-5 ketoreductase

    # Isomerases
    "K20930": ["isomerase", "lyase"],  # pmi; phosphinomethylmalate isomerase
    "K21212": ["lyase", "dehydratase", "carbohydrate_active"],  # ncsC2; NDP-hexose 2,3-dehydratase

    # Monooxygenases - cytochrome P450
    "K16415": ["oxidoreductase", "monooxygenase", "heme_binding"],  # stiL; cytochrome P450 dependent monooxygenase
    "K15966": ["oxidoreductase", "monooxygenase"],  # mtmOIV; monooxygenase
    "K13608": ["oxidoreductase", "monooxygenase", "secondary_metabolism"],  # SNO; senecionine N-oxygenase
    "K27593": ["oxidoreductase", "monooxygenase", "secondary_metabolism"],  # asuE1; protoasukamycin 4-monooxygenase
    "K00492": ["oxidoreductase", "monooxygenase"],  # tmuM; 1,3,7-trimethyluric acid 5-monooxygenase
    "K20944": ["oxidoreductase", "monooxygenase"],  # tsdB; resorcinol 4-hydroxylase (NADH)

    # Aminotransferases
    "K21175": ["aminotransferase", "transferase", "plp_binding", "aromatic_aa_metabolism"],  # sgcD1; 2-amino-4-deoxychorismate synthase, glutamine amidotransferase comp...
    "K23397": ["aminotransferase", "transferase", "plp_binding"],  # cpkG; 5-hydroxydodecatetraenal 1-aminotransferase
    "K20442": ["aminotransferase", "transferase", "plp_binding"],  # cetM; 2-keto-5-epi-valolone aminotransferase
    "K20591": ["aminotransferase", "transferase", "plp_binding"],  # genS2; pyridoxal phosphate-dependent aminotransferase
    "K26217": ["aminotransferase", "transferase", "plp_binding", "lps_biosynthesis"],  # wbgX; UDP-4-amino-4-deoxy-L-arabinose-oxoglutarate aminotransferase
    "K24305": ["hydrolase", "carbohydrate_active"],  # K24305; UDP-2-acetamido-4-(D-alanylamino)-trideoxy-mannopyranose hydrolase
    "K21326": ["aminotransferase", "transferase", "plp_binding"],  # sgcA4; aminotransferase

    # Circadian clock and signaling
    "K08482": ["atpase", "atp_binding", "signaling"],  # kaiC; circadian clock protein KaiC

    # Lyases
    "K20615": ["lyase", "amino_acid_biosynthesis"],  # vioD; capreomycidine synthase
    "K22004": ["lyase", "decarboxylase"],  # TAD1; trans-aconitate decarboxylase

    # Hydrolases and deaminases
    "K05394": ["hydrolase", "dehalogenase"],  # atzA; atrazine chlorohydrolase
    "K21593": ["hydrolase", "amidase"],  # triA; melamine deaminase
    "K21287": ["hydrolase", "amidase"],  # atdA2; gamma-glutamylanilide hydrolase
    "K01471": ["hydrolase", "amidase"],  # nylA; 6-aminohexanoate-cyclic-dimer hydrolase
    "K03391": ["oxidoreductase", "monooxygenase"],  # pcpB; pentachlorophenol monooxygenase

    # Synthases - secondary metabolism
    "K21934": ["transferase", "lyase"],  # frbC; 2-phosphonomethylmalate synthase
    "K12910": ["transferase"],  # pmmS; 2-phosphinomethylmalate synthase

    # NRPS - nonribosomal peptide synthesis
    "K23646": ["nrps", "secondary_metabolism", "ligase"],  # DTPA; ditryptophenaline biosynthesis nonribosomal peptide synthase
    "K15384": ["nrps", "polyketide_synthesis", "secondary_metabolism"],  # APDA; aspyridone synthetase, hybrid polyketide synthase / nonribosomal pe...
    "K16106": ["nrps", "secondary_metabolism"],  # blmVI; nonribosomal peptide synthetase protein BlmVI
    "K16097": ["nrps", "secondary_metabolism", "antibiotic_biosynthesis"],  # grsB; gramicidin S synthase 2
    "K23888": ["nrps", "secondary_metabolism"],  # ROQA; histidyltryptophyldiketopiperazine synthetase
    "K16108": ["nrps", "secondary_metabolism"],  # blmX; nonribosomal peptide synthetase protein BlmX
    "K16109": ["nrps", "secondary_metabolism"],  # blmIX; nonribosomal peptide synthetase protein BlmIX
    "K23631": ["nrps", "secondary_metabolism"],  # ASAC; aspergillic acid biosynthesis nonribosomal peptide synthetase
    "K16112": ["nrps", "secondary_metabolism"],  # blmIV; nonribosomal peptide synthetase protein BlmIV
    "K16096": ["nrps", "secondary_metabolism", "antibiotic_biosynthesis"],  # grsA; gramicidin S synthase 1
    "K16117": ["nrps", "secondary_metabolism", "antibiotic_biosynthesis"],  # snbC; pristinamycin I synthase 2
    "K15391": ["nrps", "polyketide_synthesis", "secondary_metabolism"],  # CPAA; cyclopiazonic acid synthetase, hybrid polyketide synthase / nonribo...
    "K16116": ["nrps", "secondary_metabolism", "antibiotic_biosynthesis"],  # snbA; pristinamycin I synthetase 1
    "K23627": ["nrps", "secondary_metabolism"],  # ASQK; cyclopeptine synthase
    "K12913": ["nrps", "secondary_metabolism", "antibiotic_biosynthesis"],  # phsC; phosphinothricin tripeptide synthetase PhsC
    "K12912": ["nrps", "secondary_metabolism", "antibiotic_biosynthesis"],  # phsB; phosphinothricin tripeptide synthetase phsB
    "K22581": ["nrps", "secondary_metabolism"],  # TAS1; tenuazonic acid synthetase

    # PKS - polyketide synthesis
    "K16395": ["polyketide_synthesis", "secondary_metabolism"],  # epoB; epothilone synthetase B
    "K15393": ["oxidoreductase", "secondary_metabolism"],  # CPAO; beta-cyclopiazonate dehydrogenase
    "K27305": ["polyketide_synthesis", "secondary_metabolism"],  # aprG; polyketide synthase / 4-carboxy-3-alkylbut-2-enoyl-[acp] decarboxylase
    "K14245": ["polyketide_synthesis", "secondary_metabolism"],  # auaC; PKS ketosynthase (KS/KS alpha)
    "K21107": ["fatty_acid_synthesis", "lipid_metabolism"],  # jamA; medium-chain-fatty-acid---[acyl-carrier-protein] ligase

    # Archaeal flagella
    "K07332": ["flagellum", "atpase", "atp_binding"],  # flaI; archaeal flagellar protein FlaI
    "K07333": ["flagellum", "membrane"],  # flaJ; archaeal flagellar protein FlaJ

    # Ligases and related
    "K12719": ["ligase", "secondary_metabolism"],  # cloN4; L-proline---[L-prolyl-carrier protein] ligase
    "K14417": ["ligase", "secondary_metabolism"],  # fcbA; 4-chlorobenzoate-CoA ligase
    "K23108": ["ligase", "secondary_metabolism"],  # BZO1; benzoate---CoA ligase
    "K20419": ["ligase", "secondary_metabolism"],  # atdA1; gamma-glutamylanilide synthase

    # Ferredoxin components
    "K24245": ["iron_sulfur", "electron_transport", "metal_binding"],  # ligXd; 5,5'-dehydrodivanillate O-demethylase ferredoxin reductase subunit
    "K22362": ["iron_sulfur", "electron_transport", "reductase"],  # amoF; alkene monooxygenase ferredoxin reductase component
    "K20809": ["iron_sulfur", "electron_transport"],  # cbdC; 2-halobenzoate 1,2-dioxygenase electron transfer component

    # Radical SAM enzymes
    "K20575": ["radical_sam", "iron_sulfur", "lyase"],  # aprD4; radical SAM diol-dehydratase
    "K20593": ["radical_sam", "iron_sulfur", "cobalamin_binding", "methyltransferase"],  # genD1; cobalamin-dependent radical SAM methyltransferase

    # Miscellaneous enzymes
    "K17048": ["oxidoreductase", "iron_sulfur"],  # ebdB; ethylbenzene hydroxylase subunit beta
    "K20439": ["phosphatase", "isomerase"],  # salP; putative phosphohexomutase / phosphatase
    "K19632": ["oxidoreductase", "carbohydrate_active"],  # rfbJ; CDP-abequose synthase
    "K27986": ["transferase", "secondary_metabolism"],  # nasB; branched-chain 2-oxoacid:malonyl-[acyl-carrier protein] acyltransfe...
    "K19884": ["oxidoreductase", "secondary_metabolism"],  # rebO; 7-chloro-L-tryptophan oxidase
    "K13549": ["oxidoreductase", "dehydrogenase"],  # btrN; 2-deoxy-scyllo-inosamine dehydrogenase (SAM-dependent)
    "K27193": ["esterase", "hydrolase"],  # mheI; carbendazim hydrolysing esterase
    "K20586": ["oxidoreductase"],  # livW; oxidoreductase
    "K17055": ["lyase", "terpene_synthesis"],  # EGS1; eugenol synthase
    # WARNING: K05914 (luciferase) is often mis-assigned to NRPS proteins by KofamScan
    # due to AMP-binding domain similarity. If a protein has both K05914 and NRPS annotations,
    # the bioluminescence assignment is likely a FALSE POSITIVE.
    "K05914": ["oxidoreductase", "bioluminescence"],  # E1.13.12.7; photinus-luciferin 4-monooxygenase (ATP-hydrolysing)
    "K15635": ["isomerase", "glycolysis"],  # apgM; 2,3-bisphosphoglycerate-independent phosphoglycerate mutase

    # Transcription factors
    "K03124": ["transcription", "dna_binding"],  # TFIIB; transcription initiation factor TFIIB
    "K03120": ["transcription", "dna_binding"],  # TBP; transcription initiation factor TFIID TATA-box-binding protein

    # DNA replication and repair
    "K02684": ["primase", "replication"],  # PRI1; DNA primase small subunit
    "K04799": ["nuclease", "dna_repair", "base_excision_repair"],  # FEN1; flap endonuclease-1
    "K01520": ["hydrolase", "nucleotide_metabolism"],  # dut; dUTP diphosphatase

    # CRISPR-Cas proteins (direct mapping for missing definitions)
    "K19091": ["crispr_associated", "cas_domain", "nuclease", "defense_system"],  # cas6; CRISPR-associated endoribonuclease Cas6
    "K19144": ["crispr_associated", "cas_domain", "nuclease", "defense_system"],  # csx3; CRISPR-associated protein Csx3

    # High-volume generic KOFAM assignments with specific definitions
    "K00610": ["transferase", "pyrimidine_metabolism"],  # pyrI; aspartate carbamoyltransferase regulatory subunit
    "K01999": ["transporter", "amino_acid_transporter", "abc_substrate_binding"],  # livK; branched-chain amino acid transport system substrate-binding protein
    "K02600": ["transcription", "transcription_termination"],  # nusA; transcription termination/antitermination protein NusA
    "K03154": ["sulfur_metabolism"],  # thiS; sulfur carrier protein
    "K03498": ["transporter", "ion_transporter"],  # trkH; trk/ktr system potassium uptake protein
    "K03568": ["protease", "hydrolase"],  # tldD; TldD protein
    "K03701": ["dna_repair", "nucleotide_excision_repair"],  # uvrA; excinuclease ABC subunit A
    "K03702": ["dna_repair", "nucleotide_excision_repair"],  # uvrB; excinuclease ABC subunit B
    "K03703": ["dna_repair", "nucleotide_excision_repair"],  # uvrC; excinuclease ABC subunit C
    "K03926": ["metal_homeostasis", "periplasmic"],  # cutA; periplasmic divalent cation tolerance protein
    "K02238": ["transporter", "membrane"],  # comEC; competence protein ComEC
    "K03499": ["transporter", "ion_transporter"],  # trkA; trk/ktr system potassium uptake protein
    "K03555": ["dna_repair", "mismatch_repair"],  # mutS; DNA mismatch repair protein MutS
    "K03572": ["dna_repair", "mismatch_repair"],  # mutL; DNA mismatch repair protein MutL
    "K03699": ["transporter", "ion_transporter", "metal_transporter", "metal_homeostasis"],  # tlyC; magnesium and cobalt exporter, CNNM family
    "K04488": ["iron_sulfur_biosynthesis", "iron_sulfur"],  # iscU; nitrogen fixation protein NifU and related proteins
    "K04795": ["rna_processing", "rrna_modification"],  # flpA; fibrillarin-like pre-rRNA processing protein
    "K05770": ["translocase", "transporter", "membrane"],  # TSPO; translocator protein
    "K06196": ["cytochrome", "heme_biosynthesis"],  # ccdA; cytochrome c-type biogenesis protein
    "K07059": ["protease", "hydrolase", "membrane"],  # rho2; rhomboid family protease
    "K07341": ["toxin_antitoxin"],  # doc; death on curing protein
    "K07465": ["exonuclease", "nuclease", "hydrolase", "dna_repair"],  # K07465; putative RecB family exonuclease
    "K07463": ["exonuclease", "nuclease", "hydrolase", "dna_repair"],  # K07463; archaea-specific RecJ-like exonuclease
    "K07477": ["dna_binding", "rna_binding"],  # K07477; translin
    "K07579": ["methyltransferase", "transferase"],  # K07579; putative methylase
    "K09741": ["trna_modification"],  # pcc1; KEOPS complex subunit Pcc1
    "K09119": ["trna_modification"],  # cgi121; KEOPS complex subunit Cgi121
    "K13993": ["small_hsp", "chaperone", "stress_response"],  # HSP20; HSP20 family protein
    "K14623": ["dna_repair", "sos_response"],  # dinD; DNA-damage-inducible protein D
    "K15977": ["oxidoreductase"],  # K15977; putative oxidoreductase
    "K18882": ["primase", "replication"],  # priL; DNA primase large subunit
    "K12063": ["conjugation", "mobile_element", "atpase", "atp_binding"],  # traC; conjugal transfer ATP-binding protein TraC
    "K25156": ["transporter", "abc_transporter", "atp_binding"],  # evrA; viologen exporter family transport system ATP-binding protein
    "K26996": ["immune_related"],  # sdpI; immunity protein, SdpI family

    # ATPases
    "K03924": ["atpase", "atp_binding", "chaperone"],  # moxR; MoxR-like ATPase

    # Transposases
    "K07496": ["transposase", "mobile_element"],  # K07496; putative transposase

    # Asparagine synthase
    "K01953": ["ligase", "amino_acid_biosynthesis", "atp_binding"],  # asnB; asparagine synthase (glutamine-hydrolysing)

    # Pyruvate formate lyase
    "K04069": ["radical_sam", "iron_sulfur", "glycolysis", "fermentation"],  # pflA; pyruvate formate lyase activating enzyme

    # Bromoperoxidase
    "K05918": ["oxidoreductase", "halogenase"],  # bmp5; 4-hydroxybenzoate brominase (decarboxylating)

    # Uncharacterized proteins
    "K07133": ["hypothetical"],  # K07133; uncharacterized protein

    # Viral proteins (may appear in metagenomes)
    "K26370": ["rna_polymerase", "transcription", "phage_related"],  # IIV6-343L; Iridovirus probable DNA-directed RNA polymerase subunit
    "K21664": ["phage_related", "dna_binding"],  # vLANA; KSHV latency-associated nuclear antigen

    # Ubiquitin system
    "K15343": ["ubiquitin_ligase", "protein_modification"],  # sspH1; E3 ubiquitin-protein ligase SspH1
}


# ============================================================================
# PATTERN-BASED MAPPINGS FOR KEGG DEFINITIONS
# ============================================================================
KEGG_PATTERNS: list[tuple[str, list[str]]] = [
    # Enzymes - core classes
    (r"\bdehydrogenase\b", ["dehydrogenase", "oxidoreductase"]),
    (r"\breductase\b", ["reductase", "oxidoreductase"]),
    (r"\boxidase\b", ["oxidase", "oxidoreductase"]),
    (r"\bkinase\b", ["kinase", "transferase"]),
    (r"\bphosphatase\b", ["phosphatase", "hydrolase"]),
    (r"\btransferase\b", ["transferase"]),
    (r"\bhydrolase\b", ["hydrolase"]),
    (r"\blyase\b", ["lyase"]),
    (r"\bsynthase\b", ["synthase"]),
    (r"\bsynthetase\b", ["synthetase", "ligase"]),
    (r"\bprotease\b|\bpeptidase\b", ["protease", "hydrolase"]),
    (r"\bnuclease\b", ["nuclease", "hydrolase"]),
    (r"\bligase\b", ["ligase"]),
    (r"\bisomerase\b", ["isomerase"]),

    # More specific enzyme patterns
    (r"monooxygenase", ["oxidoreductase", "monooxygenase"]),
    (r"dioxygenase", ["oxidoreductase", "dioxygenase"]),
    (r"peroxidase", ["oxidoreductase", "peroxidase"]),
    (r"oxygenase", ["oxidoreductase", "oxygenase"]),
    (r"aminotransferase|transaminase", ["aminotransferase", "transferase", "plp_binding"]),
    (r"acetyltransferase", ["acetyltransferase", "transferase"]),
    (r"methyltransferase|O-methyltransferase|N-methyltransferase", ["methyltransferase", "transferase", "sam_binding"]),
    (r"glycosyltransferase|glucosyltransferase|galactosyltransferase|mannosyltransferase", ["glycosyltransferase", "transferase", "carbohydrate_active"]),
    (r"phosphotransferase", ["kinase", "transferase"]),
    (r"deaminase", ["hydrolase", "amidase"]),
    (r"dehalogenase|chlorohydrolase", ["hydrolase", "dehalogenase"]),
    (r"decarboxylase", ["lyase", "decarboxylase"]),
    (r"dehydratase", ["lyase", "dehydratase"]),
    (r"\bepimerase\b", ["isomerase", "epimerase"]),
    (r"\bmutase\b", ["isomerase", "mutase"]),
    (r"\bracemase\b", ["isomerase", "racemase"]),

    # Transporters
    (r"\btransporter\b|\bpermease\b", ["transporter"]),
    (r"ABC.*transporter|ABC.*transport system|ATP-binding cassette", ["transporter", "abc_transporter", "atp_binding"]),
    (r"\bchannel\b", ["transporter", "ion_channel"]),
    (r"\bsymporter\b", ["transporter", "symporter"]),
    (r"\bantiporter\b", ["transporter", "antiporter"]),
    (r"efflux", ["transporter", "efflux_pump"]),

    # Regulators
    (r"transcription.*regulator|regulator.*transcription|transcriptional repressor", ["regulator", "transcription_factor"]),
    (r"response.*regulator", ["response_regulator", "two_component"]),
    (r"sensor.*kinase|histidine.*kinase", ["sensor_kinase", "two_component"]),
    (r"sigma.*factor", ["sigma_factor"]),
    (r"transcription.*factor|TFIIB|TBP", ["transcription", "dna_binding"]),
    (r"HTH-type|DNA[-/ ].*binding|DNA/RNA-binding", ["dna_binding"]),
    (r"nucleotide binding protein", ["nucleotide_binding", "binding"]),

    # Metabolism keywords
    # NOTE: \b word boundary prevents "dehydrogenase" from matching "hydrogenase"
    (r"\bhydrogenase\b", ["hydrogenase", "hydrogen_metabolism"]),
    (r"\bnitrogenase\b", ["nitrogenase", "nitrogen_fixation"]),
    (r"nitrate.*reductase|nar[GHI]", ["nitrate_reduction", "denitrification"]),
    (r"nitrite.*reductase|nir[SK]", ["denitrification"]),
    (r"sulfite.*reductase|dsr[AB]", ["sulfate_reduction", "sulfur_metabolism"]),
    (r"methyl.*coenzyme.*M.*reductase|mcr[ABG]", ["methanogenesis", "one_carbon_metabolism"]),
    # NOTE: RuBisCO alone does NOT mean Calvin cycle - RuBisCO-like proteins (RLPs) exist.
    (r"RuBisCO|ribulose.*bisphosphate.*carboxylase|rbc[LS]", ["rubisco", "carbon_fixation"]),
    (r"photosystem|psa[AB]|psb[AD]", ["photosynthesis"]),

    # Carbohydrate metabolism
    (r"O-antigen.*biosynthesis|O\d+-antigen", ["glycosyltransferase", "lps_biosynthesis", "carbohydrate_active"]),
    (r"dTDP.*sugar|dTDP-\d-deoxy|dTDP-.*hexose", ["carbohydrate_active"]),
    (r"UDP.*sugar|CDP.*sugar|GDP.*sugar|NDP.*sugar", ["carbohydrate_active"]),
    (r"polysaccharide", ["carbohydrate_active"]),

    # Secondary metabolism / natural products
    (r"nonribosomal.*peptide|NRPS", ["nrps", "secondary_metabolism"]),
    (r"polyketide|PKS", ["polyketide_synthesis", "secondary_metabolism"]),
    (r"halogenase|brominase", ["oxidoreductase", "secondary_metabolism"]),
    (r"biosynthesis.*protein", ["biosynthesis"]),
    (r"avermectin|erythromycin|streptomycin|kanamycin|gentamicin|novobiocin|pristinamycin|gramicidin|bleomycin", ["antibiotic_biosynthesis", "secondary_metabolism"]),

    # Structural
    (r"flagell", ["flagellum"]),
    (r"pil[iu]s|fimbr", ["pilus"]),
    (r"\bribosom", ["ribosomal_protein", "translation"]),
    (r"\bchaperone\b|chaperonin|prefoldin", ["chaperone", "stress_response"]),
    (r"zinc finger|Zn-ribbon", ["zinc_finger", "zinc_binding"]),
    (r"uncharacterized protein|hypothetical protein|unknown function", ["hypothetical"]),

    # Information processing
    (r"(translation )?(initiation|elongation|release) factor|peptide chain release factor", ["translation"]),
    (r"tRNA synthetase|aminoacyl-tRNA", ["translation", "trna_synthetase", "synthetase", "ligase"]),
    (r"proliferating cell nuclear antigen|\bPCNA\b", ["replication", "sliding_clamp"]),
    (r"replication factor C", ["replication", "clamp_loader", "atp_binding"]),
    (r"DNA replication factor|archaeal cell division control protein 6|\bcdc6", ["replication", "atpase", "dna_binding"]),
    (r"DNA repair|\bRad[AB]\b|\bRad50\b|\bMre11\b|\bSbc[CD]\b|\bRecB\b|DNA polymerase", ["dna_repair"]),
    (r"endonuclease", ["endonuclease", "nuclease", "hydrolase"]),
    (r"exosome complex|ribonuclease|[mr]RNA .*processing|RNA-binding|ribonucleoprotein", ["rna_binding", "rnase"]),
    (r"signal recognition particle receptor|\bftsY\b", ["signal_recognition", "gtp_binding"]),
    (r"protein pelota|\bpelA\b", ["translation"]),
    (r"protein archease", ["chaperone"]),
    (r"\btranslin\b", ["dna_binding", "rna_binding"]),
    (r"segregation and condensation protein|chromosome.*partition", ["chromosome_partitioning"]),
    (r"cell division protein|septum site", ["cell_division"]),

    # Mobile elements
    (r"\btransposase\b", ["transposase", "mobile_element"]),
    (r"\bintegrase\b", ["integrase", "mobile_element"]),
    (r"\bphage\b", ["phage_related"]),
    (r"recombinase", ["recombinase", "mobile_element"]),

    # Signaling
    (r"circadian|clock.*protein|KaiC", ["signaling"]),
    (r"ATPase", ["atpase", "atp_binding"]),
    (r"GTPase", ["gtpase", "gtp_binding"]),
    (r"GTP-binding", ["gtp_binding"]),

    # Cofactor/coenzyme related
    (r"radical.*SAM", ["radical_sam", "iron_sulfur"]),
    (r"Fe-S cluster assembly", ["iron_sulfur_biosynthesis", "iron_sulfur"]),
    (r"cobalamin|B12", ["cobalamin_binding", "cobalt_binding"]),
    (r"ferredoxin", ["iron_sulfur", "electron_transport"]),
    (r"flavin|FAD|FMN", ["fad_binding"]),
    (r"NAD|NADP|NADH", ["nad_binding"]),
    (r"PLP|pyridoxal", ["plp_binding"]),

    # Defense systems
    (r"CRISPR|Cas\d+", ["crispr_associated", "defense_system"]),
    (r"restriction|methylase.*DNA|DNA.*methylase", ["restriction_modification"]),
    (r"mRNA interferase", ["toxin_antitoxin", "rnase"]),
    (r"antitoxin", ["toxin_antitoxin", "antitoxin"]),

    # Cell division
    (r"\bFts[AZWLNQKX]\b", ["cell_division", "divisome"]),
    (r"\bPar[AB]\b", ["cell_division", "chromosome_partitioning"]),

    # Stress response
    (r"heat.*shock|cold.*shock", ["stress_response"]),
    (r"superoxide.*dismutase|catalase|peroxiredoxin", ["oxidative_stress"]),

    # Envelope / localization
    (r"SEC61|protein transport protein", ["transporter", "membrane"]),
    (r"membrane-associated protein|membrane protein|transmembrane protein", ["membrane"]),
    (r"MscS", ["membrane", "ion_channel"]),
    (r"multiple antibiotic resistance", ["antibiotic_resistance", "membrane"]),
]


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
        definition: KEGG definition text

    Returns:
        List of predicate IDs
    """
    definition = _resolve_ko_definition(ko_id, definition)
    predicates = set()

    # Direct mapping
    if ko_id in KEGG_TO_PREDICATES:
        predicates.update(KEGG_TO_PREDICATES[ko_id])

    # EC number mapping
    ec_numbers = parse_ec_numbers(definition)
    for ec in ec_numbers:
        predicates.update(get_predicates_for_ec(ec))

    # Pattern matching on definition
    for pattern, preds in KEGG_PATTERNS:
        if re.search(pattern, definition, re.IGNORECASE):
            predicates.update(preds)

    return sorted(predicates)


def _resolve_ko_definition(ko_id: str, definition: str = "") -> str:
    """Use local KO metadata when an annotation only stores score labels."""
    if not _is_informative_definition(ko_id, definition):
        return _load_ko_definitions().get(ko_id, definition or "")
    return definition


def _is_informative_definition(ko_id: str, definition: str = "") -> bool:
    """Return whether a KEGG definition can drive pattern/EC mapping."""
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
    "KEGG_PATTERNS",
    "KEGG_TO_PREDICATES",
    "get_predicates_for_ec",
    "get_predicates_for_kegg",
    "parse_ec_numbers",
]
