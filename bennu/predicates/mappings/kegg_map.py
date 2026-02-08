"""
KEGG ortholog to predicate mappings.

Maps KEGG KO accessions to semantic predicates using:
1. EC number classification (top-level enzyme classes)
2. Direct KO mappings for key pathways
3. Pattern-based matching on definitions
"""

from typing import Optional
import re


# ============================================================================
# EC NUMBER TO PREDICATE MAPPINGS
# ============================================================================
# Maps EC class prefixes to predicates

EC_TO_PREDICATES: dict[str, list[str]] = {
    # EC 1: Oxidoreductases
    "1": ["oxidoreductase"],
    "1.1": ["oxidoreductase", "dehydrogenase"],  # Acting on CH-OH
    "1.2": ["oxidoreductase"],  # Acting on aldehyde/oxo
    "1.3": ["oxidoreductase"],  # Acting on CH-CH
    "1.4": ["oxidoreductase", "aminotransferase"],  # Acting on CH-NH2
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
    "1.18": ["oxidoreductase", "iron_sulfur"],  # Acting on iron-sulfur
    "1.97": ["oxidoreductase"],  # Other oxidoreductases

    # EC 2: Transferases
    "2": ["transferase"],
    "2.1": ["transferase", "methyltransferase"],  # One-carbon groups
    "2.3": ["transferase", "acetyltransferase"],  # Acyltransferases
    "2.4": ["transferase", "glycosyltransferase"],  # Glycosyltransferases
    "2.5": ["transferase"],  # Alkyl/aryl transferases
    "2.6": ["transferase", "aminotransferase", "plp_binding"],  # N-group transferases
    "2.7": ["transferase", "kinase", "atp_binding"],  # Phosphotransferases
    "2.7.1": ["transferase", "kinase", "atp_binding"],  # Kinases (alcohol acceptors)
    "2.7.2": ["transferase", "kinase", "atp_binding"],  # Kinases (carboxyl acceptors)
    "2.7.4": ["transferase", "kinase", "atp_binding"],  # Phosphotransferases
    "2.7.7": ["transferase", "nucleotidyltransferase"],  # Nucleotidyltransferases
    "2.7.11": ["transferase", "serine_threonine_kinase", "kinase"],  # Ser/Thr kinases
    "2.7.12": ["transferase", "kinase"],  # Dual-specificity kinases
    "2.7.13": ["transferase", "sensor_kinase", "two_component"],  # Histidine kinases
    "2.8": ["transferase", "sulfurtransferase"],  # S-containing groups

    # EC 3: Hydrolases
    "3": ["hydrolase"],
    "3.1": ["hydrolase", "esterase"],  # Ester bonds
    "3.1.1": ["hydrolase", "esterase", "lipase"],  # Carboxylic ester
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
    "3.4.11": ["hydrolase", "protease", "aminopeptidase"],  # Aminopeptidases
    "3.4.21": ["hydrolase", "protease", "serine_protease"],  # Serine proteases
    "3.4.22": ["hydrolase", "protease", "cysteine_protease"],  # Cysteine proteases
    "3.4.23": ["hydrolase", "protease", "aspartic_protease"],  # Aspartic proteases
    "3.4.24": ["hydrolase", "protease", "metalloprotease", "metal_binding"],  # Metalloproteases
    "3.5": ["hydrolase", "amidase"],  # C-N bonds (non-peptide)
    "3.6": ["hydrolase", "atpase"],  # Acid anhydrides
    "3.6.1": ["hydrolase", "atpase", "atp_binding"],  # In phosphorus anhydrides
    "3.6.3": ["hydrolase", "atpase", "atp_binding", "transporter"],  # Catalyzing transmembrane movement
    "3.6.4": ["hydrolase", "atpase", "atp_binding", "helicase"],  # Acting on acid anhydrides

    # EC 4: Lyases
    "4": ["lyase"],
    "4.1": ["lyase", "decarboxylase"],  # C-C lyases
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
    "6.3.1": ["ligase", "synthetase", "ammonia_assimilation"],  # Acid-ammonia ligases
    "6.3.2": ["ligase", "synthetase", "peptide_synthesis"],  # Acid-amino acid ligases
    "6.4": ["ligase", "carboxylase", "biotin_binding"],  # C-C bonds
    "6.5": ["ligase", "ligase_dna", "dna_repair"],  # Phosphoric ester bonds
    "6.6": ["ligase"],  # N-metal bonds

    # EC 7: Translocases
    "7": ["translocase", "transporter", "membrane"],
    "7.1": ["translocase", "transporter", "electron_transport"],  # Linked to oxidoreduction
    "7.2": ["translocase", "transporter", "ion_transporter"],  # Linked to hydrolysis
    "7.3": ["translocase", "transporter"],  # Linked to decarboxylation
    "7.4": ["translocase", "transporter"],  # Linked to redox
    "7.5": ["translocase", "transporter", "atp_binding"],  # Linked to ATP hydrolysis
    "7.6": ["translocase", "transporter", "gtp_binding"],  # Linked to GTP hydrolysis
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
    "K00532": ["hydrogenase", "nife_hydrogenase", "hydrogen_metabolism"],  # hydB, hydrogenase large subunit
    "K00533": ["hydrogenase", "nife_hydrogenase", "hydrogen_metabolism"],  # hydA, hydrogenase small subunit
    "K00534": ["hydrogenase", "nife_hydrogenase", "hydrogen_metabolism"],  # hydG
    "K06281": ["hydrogenase", "nife_hydrogenase", "hydrogen_metabolism"],  # hoxH bidirectional
    "K06282": ["hydrogenase", "nife_hydrogenase", "hydrogen_metabolism"],  # hoxU bidirectional
    "K18008": ["hydrogenase", "fefe_hydrogenase", "hydrogen_metabolism"],  # hydA, FeFe hydrogenase

    # Group 1 - Uptake hydrogenases (H2-oxidizing)
    "K00436": ["hydrogenase", "nife_hydrogenase", "hydrogenase_uptake", "hydrogen_metabolism"],  # hyaA, uptake large
    "K05922": ["hydrogenase", "nife_hydrogenase", "hydrogenase_uptake", "hydrogen_metabolism"],  # hupL, uptake large
    "K05927": ["hydrogenase", "nife_hydrogenase", "hydrogenase_uptake", "hydrogen_metabolism"],  # hupS, uptake small

    # Group 3 - Bidirectional / F420-reducing hydrogenases
    "K00437": ["hydrogenase", "nife_hydrogenase", "hydrogenase_group3", "hydrogen_metabolism"],  # hoxU bidirectional
    "K00440": ["hydrogenase", "nife_hydrogenase", "hydrogenase_group3", "f420_dependent", "hydrogen_metabolism"],  # frhA
    "K00441": ["hydrogenase", "nife_hydrogenase", "hydrogenase_group3", "f420_dependent", "hydrogen_metabolism"],  # frhB
    "K00443": ["hydrogenase", "nife_hydrogenase", "hydrogenase_group3", "f420_dependent", "hydrogen_metabolism"],  # frhG

    # Group 4 - Energy-conserving / H2-evolving hydrogenases
    "K15826": ["hydrogenase", "nife_hydrogenase", "ech_hydrogenase", "hydrogenase_group4", "hydrogen_metabolism"],  # echA
    "K15827": ["hydrogenase", "nife_hydrogenase", "ech_hydrogenase", "hydrogenase_group4", "hydrogen_metabolism"],  # echB
    "K15828": ["hydrogenase", "nife_hydrogenase", "ech_hydrogenase", "hydrogenase_group4", "hydrogen_metabolism"],  # echC
    "K15829": ["hydrogenase", "nife_hydrogenase", "ech_hydrogenase", "hydrogenase_group4", "hydrogen_metabolism"],  # echD
    "K15830": ["hydrogenase", "nife_hydrogenase", "ech_hydrogenase", "hydrogenase_group4", "hydrogen_metabolism"],  # echE
    "K15831": ["hydrogenase", "nife_hydrogenase", "ech_hydrogenase", "hydrogenase_group4", "hydrogen_metabolism"],  # echF
    "K14126": ["hydrogenase", "nife_hydrogenase", "mbh_hydrogenase", "hydrogenase_group4", "hydrogen_metabolism", "membrane"],  # mbhJ
    "K14127": ["hydrogenase", "nife_hydrogenase", "mbh_hydrogenase", "hydrogenase_group4", "hydrogen_metabolism", "membrane"],  # mbhK
    "K14128": ["hydrogenase", "nife_hydrogenase", "mbh_hydrogenase", "hydrogenase_group4", "hydrogen_metabolism", "membrane"],  # mbhL
    "K14129": ["hydrogenase", "nife_hydrogenase", "mbh_hydrogenase", "hydrogenase_group4", "hydrogen_metabolism", "membrane"],  # mbhM
    "K14130": ["hydrogenase", "nife_hydrogenase", "mbh_hydrogenase", "hydrogenase_group4", "hydrogen_metabolism", "membrane"],  # mbhN

    # Hydrogenase maturation
    "K04651": ["hydrogenase_maturation", "nickel_binding"],  # hypA
    "K04652": ["hydrogenase_maturation"],  # hypB
    "K04653": ["hydrogenase_maturation"],  # hypC
    "K04654": ["hydrogenase_maturation"],  # hypD
    "K04655": ["hydrogenase_maturation"],  # hypE
    "K04656": ["hydrogenase_maturation"],  # hypF

    # -------------------------------------------------------------------------
    # NITROGEN FIXATION
    # -------------------------------------------------------------------------
    # NOTE: Only nifHDK (core nitrogenase subunits) should trigger nitrogen_fixation.
    # Accessory genes (nifBNE) are involved in cofactor biosynthesis but organisms
    # can have these genes without actually fixing nitrogen (genes repurposed).
    "K02586": ["nitrogenase", "nitrogen_fixation", "iron_sulfur"],  # nifH - CORE
    "K02591": ["nitrogenase", "nitrogen_fixation", "molybdenum_binding"],  # nifD - CORE
    "K02592": ["nitrogenase", "nitrogen_fixation"],  # nifK - CORE
    # Accessory genes - do NOT indicate nitrogen fixation without nifHDK
    "K02587": ["nitrogenase_maturation", "iron_sulfur"],  # nifB - FeMo-co biosynthesis
    "K02594": ["nitrogenase_maturation"],  # nifN - FeMo-co scaffold
    "K02593": ["nitrogenase_maturation"],  # nifE - FeMo-co scaffold

    # -------------------------------------------------------------------------
    # DENITRIFICATION
    # -------------------------------------------------------------------------
    "K00370": ["denitrification", "nitrate_reduction", "molybdenum_binding"],  # narG
    "K00371": ["denitrification", "nitrate_reduction", "iron_sulfur"],  # narH
    "K00374": ["denitrification", "nitrate_reduction"],  # narI
    "K02567": ["denitrification", "nitrite_reductase"],  # napA
    "K00368": ["denitrification", "nitrite_reductase", "heme_binding"],  # nirK
    "K15864": ["denitrification", "nitrite_reductase", "heme_binding"],  # nirS
    "K00376": ["denitrification", "nitrous_oxide_reductase", "copper_binding"],  # nosZ
    "K04561": ["denitrification", "nitric_oxide_reductase"],  # norB

    # -------------------------------------------------------------------------
    # SULFUR METABOLISM
    # -------------------------------------------------------------------------
    "K11180": ["sulfate_reduction", "sulfur_metabolism"],  # dsrA
    "K11181": ["sulfate_reduction", "sulfur_metabolism"],  # dsrB
    "K00394": ["sulfate_reduction", "sulfur_metabolism", "atp_binding"],  # aprA
    "K00395": ["sulfate_reduction", "sulfur_metabolism"],  # aprB
    "K00958": ["sulfur_assimilation", "sulfur_metabolism", "atp_binding"],  # sat
    "K00860": ["sulfur_assimilation", "sulfur_metabolism"],  # cysC
    "K17218": ["sulfur_oxidation", "sulfur_metabolism"],  # soxA
    "K17222": ["sulfur_oxidation", "sulfur_metabolism"],  # soxB
    "K17223": ["sulfur_oxidation", "sulfur_metabolism"],  # soxC
    "K17226": ["sulfur_oxidation", "sulfur_metabolism"],  # soxY
    "K17227": ["sulfur_oxidation", "sulfur_metabolism"],  # soxZ

    # -------------------------------------------------------------------------
    # METHANOGENESIS
    # -------------------------------------------------------------------------
    # NOTE: Only MCR (methyl-CoM reductase) subunits are DEFINITIVE for methanogenesis.
    # Other enzymes (HDR, Fwd, Mtr, Mer) are shared with non-methanogenic archaea.
    "K00399": ["methanogenesis", "mcr_complex", "one_carbon_metabolism", "nickel_binding"],  # mcrA - DEFINITIVE
    "K00401": ["methanogenesis", "mcr_complex", "one_carbon_metabolism"],  # mcrB - DEFINITIVE
    "K00402": ["methanogenesis", "mcr_complex", "one_carbon_metabolism"],  # mcrG - DEFINITIVE
    # HDR, Fwd, Mtr, Mer - shared with non-methanogens (e.g., Altiarchaeota)
    "K03388": ["archaeal_one_carbon", "one_carbon_metabolism", "iron_sulfur"],  # hdrA
    "K03389": ["archaeal_one_carbon", "one_carbon_metabolism"],  # hdrB
    "K03390": ["archaeal_one_carbon", "one_carbon_metabolism"],  # hdrC
    "K00200": ["archaeal_one_carbon", "one_carbon_metabolism", "wood_ljungdahl"],  # fwdA - formylmethanofuran dehydrogenase
    "K00201": ["archaeal_one_carbon", "one_carbon_metabolism", "wood_ljungdahl"],  # fwdB
    "K00577": ["archaeal_one_carbon", "one_carbon_metabolism", "cobalamin_binding"],  # mtrA
    "K00320": ["archaeal_one_carbon", "one_carbon_metabolism", "fad_binding"],  # mer

    # Methanotrophy
    "K10944": ["methanotrophy", "one_carbon_metabolism", "copper_binding"],  # pmoA
    "K10945": ["methanotrophy", "one_carbon_metabolism"],  # pmoB
    "K10946": ["methanotrophy", "one_carbon_metabolism"],  # pmoC
    "K16157": ["methanotrophy", "one_carbon_metabolism"],  # mmoX
    "K16158": ["methanotrophy", "one_carbon_metabolism"],  # mmoY
    "K16159": ["methanotrophy", "one_carbon_metabolism"],  # mmoZ

    # -------------------------------------------------------------------------
    # WOOD-LJUNGDAHL PATHWAY
    # -------------------------------------------------------------------------
    "K00192": ["wood_ljungdahl", "one_carbon_metabolism", "nickel_binding"],  # cdhC (CODH/ACS)
    "K00193": ["wood_ljungdahl", "one_carbon_metabolism"],  # cdhD
    "K00194": ["wood_ljungdahl", "one_carbon_metabolism"],  # cdhE
    "K00197": ["wood_ljungdahl", "one_carbon_metabolism", "iron_sulfur"],  # cdhA
    "K00198": ["wood_ljungdahl", "co_oxidation", "one_carbon_metabolism"],  # cooS (CODH)

    # -------------------------------------------------------------------------
    # PHOTOSYNTHESIS
    # -------------------------------------------------------------------------
    "K02703": ["photosynthesis", "photosystem_ii"],  # psbA (D1)
    "K02706": ["photosynthesis", "photosystem_ii"],  # psbD (D2)
    "K02689": ["photosynthesis", "photosystem_i"],  # psaA
    "K02690": ["photosynthesis", "photosystem_i"],  # psaB
    "K02634": ["photosynthesis", "electron_transport"],  # petA (cytochrome f)
    "K02635": ["photosynthesis", "electron_transport"],  # petB (cytochrome b6)
    "K02636": ["photosynthesis", "electron_transport", "iron_sulfur"],  # petC (Rieske)
    # NOTE: RuBisCO alone does NOT mean Calvin cycle - RuBisCO-like proteins (RLPs) exist.
    # Only assign rubisco and carbon_fixation, NOT calvin_cycle (requires PRK confirmation).
    "K01601": ["rubisco", "carbon_fixation"],  # rbcL (RuBisCO large) - NOT calvin_cycle without PRK
    "K01602": ["rubisco", "carbon_fixation"],  # rbcS (RuBisCO small) - NOT calvin_cycle without PRK

    # -------------------------------------------------------------------------
    # ATP SYNTHESIS
    # -------------------------------------------------------------------------
    "K02111": ["atp_synthesis", "energy_metabolism", "atp_binding"],  # atpA (F1 alpha)
    "K02112": ["atp_synthesis", "energy_metabolism", "atp_binding"],  # atpD (F1 beta)
    "K02113": ["atp_synthesis", "energy_metabolism"],  # atpG (F1 gamma)
    "K02109": ["atp_synthesis", "energy_metabolism", "membrane"],  # atpB (F0 a)
    "K02108": ["atp_synthesis", "energy_metabolism", "membrane"],  # atpE (F0 c)

    # -------------------------------------------------------------------------
    # CENTRAL CARBON METABOLISM
    # -------------------------------------------------------------------------
    # Glycolysis
    "K00844": ["glycolysis", "central_metabolism", "kinase"],  # HK (hexokinase)
    "K01810": ["glycolysis", "central_metabolism", "isomerase"],  # GPI
    "K00850": ["glycolysis", "central_metabolism", "kinase"],  # pfkA (PFK)
    "K01623": ["glycolysis", "central_metabolism", "lyase"],  # FBA (aldolase)
    "K00134": ["glycolysis", "central_metabolism", "dehydrogenase", "nad_binding"],  # GAPDH
    "K00927": ["glycolysis", "central_metabolism", "kinase"],  # PGK
    "K01689": ["glycolysis", "central_metabolism"],  # ENO (enolase)
    "K00873": ["glycolysis", "central_metabolism", "kinase"],  # pykF (pyruvate kinase)

    # TCA cycle
    "K01647": ["tca_cycle", "central_metabolism"],  # gltA (citrate synthase)
    "K01681": ["tca_cycle", "central_metabolism"],  # ACO (aconitase)
    "K00031": ["tca_cycle", "central_metabolism", "dehydrogenase", "nad_binding"],  # IDH
    "K00164": ["tca_cycle", "central_metabolism", "dehydrogenase"],  # OGDH (2-oxoglutarate DH)
    "K01902": ["tca_cycle", "central_metabolism", "ligase"],  # sucC (succinyl-CoA synthetase)
    "K00239": ["tca_cycle", "central_metabolism", "dehydrogenase", "fad_binding"],  # sdhA (SDH)
    "K01676": ["tca_cycle", "central_metabolism"],  # fumA (fumarase)
    "K00024": ["tca_cycle", "central_metabolism", "dehydrogenase", "nad_binding"],  # mdh (malate DH)

    # -------------------------------------------------------------------------
    # TRANSPORTERS
    # -------------------------------------------------------------------------
    "K02003": ["transporter", "abc_transporter", "abc_atpase", "atp_binding"],  # ABC.CD.A
    "K02004": ["transporter", "abc_transporter", "abc_permease"],  # ABC.CD.P
    "K02020": ["transporter", "abc_transporter", "peptide_transporter"],  # oppA
    "K02035": ["transporter", "abc_transporter", "phosphate_transporter"],  # pstS
    "K02036": ["transporter", "abc_transporter", "phosphate_transporter"],  # pstC
    "K02037": ["transporter", "abc_transporter", "phosphate_transporter"],  # pstA
    "K02038": ["transporter", "abc_transporter", "phosphate_transporter", "atp_binding"],  # pstB
    "K02040": ["transporter", "abc_transporter", "sulfate_transporter"],  # sbp
    "K02041": ["transporter", "abc_transporter", "sulfate_transporter"],  # cysP
    "K02046": ["transporter", "abc_transporter", "sugar_transporter"],  # mglA
    "K02055": ["transporter", "abc_transporter", "amino_acid_transporter"],  # livK
    "K03088": ["sigma_factor", "regulator"],  # rpoE (sigma-E)

    # -------------------------------------------------------------------------
    # TWO-COMPONENT SYSTEMS
    # -------------------------------------------------------------------------
    "K07638": ["two_component", "sensor_kinase", "regulator"],  # envZ
    "K07659": ["two_component", "response_regulator", "regulator"],  # ompR
    "K07648": ["two_component", "sensor_kinase", "regulator"],  # phoR
    "K07657": ["two_component", "response_regulator", "regulator"],  # phoB
    "K07678": ["two_component", "sensor_kinase", "regulator"],  # narX
    "K07684": ["two_component", "response_regulator", "regulator", "narl_family"],  # narL
    "K02483": ["two_component", "sensor_kinase", "regulator", "chemotaxis"],  # cheA
    "K03407": ["two_component", "response_regulator", "chemotaxis"],  # cheY

    # -------------------------------------------------------------------------
    # SECRETION SYSTEMS
    # -------------------------------------------------------------------------
    "K03070": ["secretion_system", "sec_pathway", "translocase"],  # secA
    "K03071": ["secretion_system", "sec_pathway", "membrane"],  # secY
    "K03072": ["secretion_system", "sec_pathway", "membrane"],  # secE
    "K03073": ["secretion_system", "sec_pathway", "membrane"],  # secG
    "K03076": ["secretion_system", "sec_pathway"],  # secB
    "K03116": ["secretion_system", "tat_pathway"],  # tatA
    "K03117": ["secretion_system", "tat_pathway"],  # tatB
    "K03118": ["secretion_system", "tat_pathway", "membrane"],  # tatC
    "K03205": ["secretion_system", "type_iii_secretion"],  # yscN
    "K03219": ["secretion_system", "type_iii_secretion"],  # yscV
    "K03194": ["secretion_system", "type_vi_secretion"],  # tssB
    "K03195": ["secretion_system", "type_vi_secretion"],  # tssC

    # -------------------------------------------------------------------------
    # MOBILE ELEMENTS / DEFENSE
    # -------------------------------------------------------------------------
    "K07481": ["transposase", "mobile_element"],  # IS transposase
    "K01356": ["integrase", "mobile_element", "recombinase"],  # int (phage integrase)
    "K03529": ["crispr_associated", "defense_system"],  # cas1
    "K09951": ["crispr_associated", "defense_system"],  # cas2
    "K07012": ["crispr_associated", "cas_nuclease", "defense_system"],  # cas3
    "K19086": ["crispr_associated", "cas_nuclease", "defense_system"],  # cas9
    "K19087": ["crispr_associated", "cas_nuclease", "defense_system"],  # cas12a

    # -------------------------------------------------------------------------
    # CELL DIVISION
    # -------------------------------------------------------------------------
    "K03531": ["cell_division", "ftsz", "gtp_binding"],  # ftsZ
    "K03587": ["cell_division", "divisome"],  # ftsA
    "K03589": ["cell_division", "divisome"],  # ftsK
    "K03590": ["cell_division", "divisome"],  # ftsL
    "K03591": ["cell_division", "divisome"],  # ftsN
    "K03592": ["cell_division", "divisome"],  # ftsQ
    "K03593": ["cell_division", "divisome"],  # ftsW
    "K03595": ["cell_division", "chromosome_partitioning"],  # parA
    "K03496": ["cell_division", "chromosome_partitioning"],  # parB

    # -------------------------------------------------------------------------
    # DNA REPLICATION
    # -------------------------------------------------------------------------
    "K02337": ["replication", "dna_polymerase"],  # dnaE (Pol III alpha)
    "K02338": ["replication", "dna_polymerase"],  # dnaN (Pol III beta)
    "K02340": ["replication", "dna_polymerase"],  # dnaQ (Pol III epsilon)
    "K02314": ["replication", "helicase", "atp_binding"],  # dnaB
    "K02316": ["replication", "primase"],  # dnaG
    "K02335": ["replication", "dna_polymerase"],  # dnaA
    "K02469": ["replication", "topoisomerase"],  # gyrA
    "K02470": ["replication", "topoisomerase"],  # gyrB

    # -------------------------------------------------------------------------
    # TRANSLATION
    # -------------------------------------------------------------------------
    "K02863": ["ribosomal_protein", "translation"],  # rpsA (S1)
    "K02886": ["ribosomal_protein", "translation"],  # rplA (L1)
    "K02906": ["ribosomal_protein", "translation"],  # rplL (L7/L12)
    "K02355": ["translation_factor", "translation", "gtp_binding"],  # EF-Tu
    "K02357": ["translation_factor", "translation", "gtp_binding"],  # EF-G
    "K02518": ["translation_factor", "translation"],  # IF-1
    "K02519": ["translation_factor", "translation", "gtp_binding"],  # IF-2
    "K02520": ["translation_factor", "translation"],  # IF-3
    "K01866": ["trna_synthetase", "translation"],  # tyrS
    "K01867": ["trna_synthetase", "translation"],  # leuS
    "K01868": ["trna_synthetase", "translation"],  # thrS

    # -------------------------------------------------------------------------
    # STRESS RESPONSE
    # -------------------------------------------------------------------------
    "K04043": ["heat_shock", "chaperone", "stress_response", "atp_binding"],  # dnaK (Hsp70)
    "K03686": ["heat_shock", "chaperone", "stress_response"],  # dnaJ
    "K04077": ["heat_shock", "chaperone", "stress_response"],  # groEL (Hsp60)
    "K04078": ["heat_shock", "chaperone", "stress_response"],  # groES
    "K04769": ["heat_shock", "stress_response"],  # htpG (Hsp90)
    "K03544": ["clp_protease", "chaperone", "atp_binding"],  # clpA
    "K03695": ["clp_protease", "chaperone", "atp_binding"],  # clpB
    "K01358": ["clp_protease", "protease"],  # clpP
    "K03671": ["oxidative_stress", "peroxiredoxin"],  # ahpC
    "K04564": ["oxidative_stress", "superoxide_dismutase", "manganese_binding"],  # sodA
    "K04565": ["oxidative_stress", "superoxide_dismutase", "iron_binding"],  # sodB
    "K03781": ["oxidative_stress", "catalase", "heme_binding"],  # katE
    "K03782": ["oxidative_stress", "catalase", "heme_binding"],  # katG

    # -------------------------------------------------------------------------
    # ANTIBIOTIC RESISTANCE
    # -------------------------------------------------------------------------
    "K01467": ["antibiotic_resistance", "beta_lactamase"],  # ampC
    "K18698": ["antibiotic_resistance", "beta_lactamase"],  # blaCTX-M
    "K17836": ["antibiotic_resistance", "aminoglycoside_resistance"],  # aac
    "K00984": ["antibiotic_resistance", "aminoglycoside_resistance"],  # strA
    "K18139": ["antibiotic_resistance", "efflux_pump", "multidrug_resistance"],  # acrA
    "K18140": ["antibiotic_resistance", "efflux_pump", "multidrug_resistance"],  # acrB

    # -------------------------------------------------------------------------
    # ADDITIONAL HIGH-ABUNDANCE KEGG ORTHOLOGS (from dataset analysis)
    # -------------------------------------------------------------------------
    # Glycosyltransferases - LPS/O-antigen biosynthesis
    "K20578": ["glycosyltransferase", "transferase", "carbohydrate_active"],  # glycosyltransferase EC:2.4.1.-
    "K23102": ["glycosyltransferase", "transferase", "lps_biosynthesis", "carbohydrate_active"],  # O18-antigen galactosyltransferase
    "K24654": ["glycosyltransferase", "transferase", "lps_biosynthesis", "carbohydrate_active"],  # O3-antigen galactosyltransferase
    "K24649": ["glycosyltransferase", "transferase", "lps_biosynthesis", "carbohydrate_active"],  # O111-antigen glucosyltransferase
    "K20438": ["glycosyltransferase", "transferase", "carbohydrate_active"],  # validoxylamine A glucosyltransferase
    "K23101": ["glycosyltransferase", "transferase", "lps_biosynthesis", "carbohydrate_active"],  # O18-antigen glucosyltransferase
    "K21586": ["glycosyltransferase", "transferase", "lps_biosynthesis", "carbohydrate_active"],  # O55-antigen galactosyltransferase
    "K23229": ["glycosyltransferase", "transferase", "lps_biosynthesis", "carbohydrate_active"],  # O55-antigen glycosyltransferase
    "K23105": ["glycosyltransferase", "transferase", "lps_biosynthesis", "carbohydrate_active"],  # O111-antigen colitosyltransferase
    "K23104": ["glycosyltransferase", "transferase", "lps_biosynthesis", "carbohydrate_active"],  # O111-antigen galactosyltransferase
    "K26073": ["glycosyltransferase", "transferase", "lps_biosynthesis", "carbohydrate_active"],  # O78-antigen glycosyltransferase
    "K23103": ["glycosyltransferase", "transferase", "lps_biosynthesis", "carbohydrate_active"],  # O18-antigen glycosyltransferase
    "K23130": ["glycosyltransferase", "transferase", "lps_biosynthesis", "carbohydrate_active"],  # O157-antigen glucosyltransferase
    "K23199": ["glycosyltransferase", "transferase", "lps_biosynthesis", "carbohydrate_active"],  # O7-antigen mannosyltransferase
    "K23131": ["glycosyltransferase", "transferase", "lps_biosynthesis", "carbohydrate_active"],  # O157-antigen glycosyltransferase
    "K21216": ["glycosyltransferase", "transferase", "carbohydrate_active"],  # glycosyltransferase
    "K21365": ["glycosyltransferase", "transferase", "lps_biosynthesis", "carbohydrate_active"],  # O86/O127-antigen GalNAc transferase
    "K23073": ["glycosyltransferase", "transferase", "lps_biosynthesis", "carbohydrate_active"],  # O16-antigen glucosyltransferase
    "K23250": ["glycosyltransferase", "transferase", "lps_biosynthesis", "carbohydrate_active"],  # O104-antigen galactosyltransferase
    "K23230": ["glycosyltransferase", "transferase", "lps_biosynthesis", "carbohydrate_active"],  # O55-antigen galactosyltransferase
    "K20573": ["glycosyltransferase", "transferase", "carbohydrate_active", "secondary_metabolism"],  # kanosaminyltransferase
    "K20589": ["glycosyltransferase", "transferase", "carbohydrate_active", "secondary_metabolism"],  # xylosyltransferase
    "K23251": ["glycosyltransferase", "transferase", "lps_biosynthesis", "carbohydrate_active"],  # O104-antigen galactosyltransferase

    # ABC transporters
    "K15558": ["transporter", "abc_transporter", "abc_atpase", "atp_binding"],  # phthalate transport ATP-binding
    "K20386": ["transporter", "abc_transporter", "abc_atpase", "atp_binding"],  # ABC subfamily B, bacterial

    # Selenate/chlorate reductases - anaerobic respiration
    "K27875": ["oxidoreductase", "anaerobic_respiration", "selenium_metabolism"],  # selenate reductase subunit B
    "K17051": ["oxidoreductase", "anaerobic_respiration", "selenium_metabolism"],  # selenate/chlorate reductase beta

    # Halogenases - secondary metabolism
    "K21256": ["oxidoreductase", "fad_binding", "secondary_metabolism"],  # flavin-dependent halogenase
    "K12714": ["oxidoreductase", "fad_binding", "secondary_metabolism"],  # clorobiocin halogenase

    # Methyltransferases
    "K16414": ["methyltransferase", "transferase", "sam_binding"],  # methyltransferase
    "K19889": ["methyltransferase", "transferase", "sam_binding", "secondary_metabolism"],  # O-methyltransferase
    "K27779": ["methyltransferase", "transferase", "sam_binding", "secondary_metabolism"],  # N-methyltransferase
    "K14374": ["methyltransferase", "transferase", "sam_binding", "secondary_metabolism"],  # avermectin O-methyltransferase
    "K12914": ["methyltransferase", "transferase", "sam_binding"],  # P-methyltransferase
    "K27844": ["methyltransferase", "transferase", "sam_binding", "secondary_metabolism"],  # ajmaline N-methyltransferase
    "K20592": ["methyltransferase", "transferase", "sam_binding"],  # SAM-dependent N-methyltransferase
    "K16403": ["methyltransferase", "transferase", "sam_binding"],  # O-methyltransferase
    "K12902": ["methyltransferase", "transferase", "sam_binding"],  # phosphonoacetaldehyde methylase
    "K21896": ["methyltransferase", "transferase", "sam_binding", "secondary_metabolism"],  # N-methyltransferase
    "K20594": ["methyltransferase", "transferase", "sam_binding", "secondary_metabolism"],  # gentamicin methyltransferase
    "K12633": ["methyltransferase", "transferase", "sam_binding"],  # N-methyl-transferase
    "K21541": ["methyltransferase", "transferase", "sam_binding", "terpene_synthesis"],  # squalene methyltransferase
    "K21542": ["methyltransferase", "transferase", "sam_binding", "terpene_synthesis"],  # botryococcene methyltransferase
    "K23153": ["methyltransferase", "transferase", "sam_binding", "secondary_metabolism"],  # nocamycin O-methyltransferase
    "K21458": ["methyltransferase", "transferase", "sam_binding", "secondary_metabolism"],  # pavine N-methyltransferase
    "K12705": ["methyltransferase", "transferase", "sam_binding", "secondary_metabolism"],  # novobiocin C8-methyltransferase
    "K00571": ["methyltransferase", "dna_methylase", "sam_binding"],  # site-specific DNA methyltransferase
    "K00558": ["methyltransferase", "dna_methylase", "sam_binding"],  # DNA (cytosine-5)-methyltransferase
    "K06223": ["methyltransferase", "dna_methylase", "sam_binding"],  # DNA adenine methylase

    # Acetyltransferases
    "K27979": ["acetyltransferase", "transferase", "secondary_metabolism"],  # altemicidin isoleucyltransferase
    "K20681": ["acetyltransferase", "transferase", "carbohydrate_active"],  # dTDP-sugar N-acetyltransferase
    "K24251": ["acetyltransferase", "transferase", "carbohydrate_active"],  # dTDP-D-Quip3N acetylase
    "K17939": ["acetyltransferase", "transferase", "carbohydrate_active"],  # GDP-perosamine N-acetyltransferase
    "K19313": ["acetyltransferase", "transferase", "secondary_metabolism"],  # aminocyclitol acetyltransferase
    "K20680": ["acetyltransferase", "transferase", "carbohydrate_active"],  # dTDP-sugar N-acetyltransferase
    "K26075": ["acetyltransferase", "transferase", "lps_biosynthesis"],  # O-acyltransferase

    # Dehydrogenases and reductases - carbohydrate metabolism
    "K13306": ["oxidoreductase", "dehydrogenase", "nad_binding", "carbohydrate_active"],  # dTDP-sugar reductase
    "K18836": ["oxidoreductase", "reductase"],  # reductase VcaE
    "K24122": ["oxidoreductase", "dehydrogenase", "terpene_synthesis"],  # limonene dehydrogenase
    "K15858": ["oxidoreductase", "reductase", "carbohydrate_active"],  # CDP-dideoxyhexulose reductase
    "K21271": ["oxidoreductase", "dehydrogenase", "secondary_metabolism"],  # aurachin B dehydrogenase
    "K19652": ["oxidoreductase", "dehydrogenase", "nad_binding", "carbohydrate_active"],  # dTDP-sugar dehydrogenase
    "K19859": ["oxidoreductase", "reductase", "carbohydrate_active"],  # dTDP-sugar reductase
    "K19180": ["oxidoreductase", "dehydrogenase", "nad_binding", "carbohydrate_active"],  # dTDP-sugar dehydrogenase
    "K13319": ["oxidoreductase", "reductase", "secondary_metabolism"],  # 4-ketoreductase
    "K19857": ["oxidoreductase", "reductase", "carbohydrate_active"],  # dTDP-sugar reductase
    "K22098": ["oxidoreductase", "dehydrogenase", "secondary_metabolism"],  # noscapine synthase
    "K20152": ["oxidoreductase", "reductase", "secondary_metabolism"],  # dehydrokanamycin reductase

    # Isomerases
    "K20930": ["isomerase", "lyase"],  # phosphinomethylmalate isomerase
    "K21212": ["lyase", "dehydratase", "carbohydrate_active"],  # NDP-hexose 2,3-dehydratase

    # Monooxygenases - cytochrome P450
    "K16415": ["oxidoreductase", "monooxygenase", "heme_binding"],  # cytochrome P450 monooxygenase
    "K15966": ["oxidoreductase", "monooxygenase"],  # monooxygenase
    "K13608": ["oxidoreductase", "monooxygenase", "secondary_metabolism"],  # senecionine N-oxygenase
    "K27593": ["oxidoreductase", "monooxygenase", "secondary_metabolism"],  # protoasukamycin monooxygenase
    "K00492": ["oxidoreductase", "monooxygenase"],  # trimethyluric acid monooxygenase
    "K20944": ["oxidoreductase", "monooxygenase"],  # resorcinol 4-hydroxylase

    # Aminotransferases
    "K21175": ["aminotransferase", "transferase", "plp_binding", "aromatic_aa_metabolism"],  # 2-amino-4-deoxychorismate synthase
    "K23397": ["aminotransferase", "transferase", "plp_binding"],  # hydroxydodecatetraenal aminotransferase
    "K20442": ["aminotransferase", "transferase", "plp_binding"],  # valolone aminotransferase
    "K20591": ["aminotransferase", "transferase", "plp_binding"],  # PLP-dependent aminotransferase
    "K26217": ["aminotransferase", "transferase", "plp_binding", "lps_biosynthesis"],  # UDP-L-arabinose aminotransferase
    "K21326": ["aminotransferase", "transferase", "plp_binding"],  # aminotransferase

    # Circadian clock and signaling
    "K08482": ["atpase", "atp_binding", "signaling"],  # KaiC circadian clock protein

    # Lyases
    "K20615": ["lyase", "amino_acid_biosynthesis"],  # capreomycidine synthase
    "K22004": ["lyase", "decarboxylase"],  # trans-aconitate decarboxylase

    # Hydrolases and deaminases
    "K05394": ["hydrolase", "dehalogenase"],  # atrazine chlorohydrolase
    "K21593": ["hydrolase", "amidase"],  # melamine deaminase
    "K24305": ["hydrolase", "carbohydrate_active"],  # UDP-sugar hydrolase
    "K21287": ["hydrolase", "amidase"],  # gamma-glutamylanilide hydrolase
    "K01471": ["hydrolase", "amidase"],  # 6-aminohexanoate cyclic dimer hydrolase
    "K03391": ["oxidoreductase", "monooxygenase"],  # pentachlorophenol monooxygenase

    # Synthases - secondary metabolism
    "K21934": ["transferase", "lyase"],  # 2-phosphonomethylmalate synthase
    "K12910": ["transferase"],  # 2-phosphinomethylmalate synthase

    # NRPS - nonribosomal peptide synthesis
    "K23646": ["nrps", "secondary_metabolism", "ligase"],  # ditryptophenaline NRPS
    "K15384": ["nrps", "polyketide_synthesis", "secondary_metabolism"],  # aspyridone synthetase hybrid
    "K16106": ["nrps", "secondary_metabolism"],  # bleomycin NRPS BlmVI
    "K16097": ["nrps", "secondary_metabolism", "antibiotic_biosynthesis"],  # gramicidin S synthase 2
    "K23888": ["nrps", "secondary_metabolism"],  # histidyltryptophyldiketopiperazine synthetase
    "K16108": ["nrps", "secondary_metabolism"],  # bleomycin NRPS BlmX
    "K16109": ["nrps", "secondary_metabolism"],  # bleomycin NRPS BlmIX
    "K23631": ["nrps", "secondary_metabolism"],  # aspergillic acid NRPS
    "K16112": ["nrps", "secondary_metabolism"],  # bleomycin NRPS BlmIV
    "K16096": ["nrps", "secondary_metabolism", "antibiotic_biosynthesis"],  # gramicidin S synthase 1
    "K16117": ["nrps", "secondary_metabolism", "antibiotic_biosynthesis"],  # pristinamycin I synthase 2
    "K15391": ["nrps", "polyketide_synthesis", "secondary_metabolism"],  # cyclopiazonic acid synthetase
    "K16116": ["nrps", "secondary_metabolism", "antibiotic_biosynthesis"],  # pristinamycin I synthetase 1
    "K23627": ["nrps", "secondary_metabolism"],  # cyclopeptine synthase
    "K12913": ["nrps", "secondary_metabolism", "antibiotic_biosynthesis"],  # phosphinothricin tripeptide synthetase PhsC
    "K12912": ["nrps", "secondary_metabolism", "antibiotic_biosynthesis"],  # phosphinothricin tripeptide synthetase phsB
    "K22581": ["nrps", "secondary_metabolism"],  # tenuazonic acid synthetase

    # PKS - polyketide synthesis
    "K16395": ["polyketide_synthesis", "secondary_metabolism"],  # epothilone synthetase B
    "K15393": ["oxidoreductase", "secondary_metabolism"],  # beta-cyclopiazonate dehydrogenase
    "K27305": ["polyketide_synthesis", "secondary_metabolism"],  # polyketide synthase
    "K14245": ["polyketide_synthesis", "secondary_metabolism"],  # PKS ketosynthase
    "K21107": ["fatty_acid_synthesis", "lipid_metabolism"],  # medium-chain fatty acid-ACP ligase

    # Archaeal flagella
    "K07332": ["flagellum", "atpase", "atp_binding"],  # archaeal flagellar protein FlaI
    "K07333": ["flagellum", "membrane"],  # archaeal flagellar protein FlaJ

    # Ligases and related
    "K12719": ["ligase", "secondary_metabolism"],  # L-proline-[carrier protein] ligase
    "K14417": ["ligase", "secondary_metabolism"],  # 4-chlorobenzoate-CoA ligase
    "K23108": ["ligase", "secondary_metabolism"],  # benzoate-CoA ligase
    "K20419": ["ligase", "secondary_metabolism"],  # gamma-glutamylanilide synthase

    # Ferredoxin components
    "K24245": ["iron_sulfur", "electron_transport", "metal_binding"],  # O-demethylase ferredoxin
    "K22362": ["iron_sulfur", "electron_transport", "reductase"],  # alkene monooxygenase ferredoxin reductase
    "K20809": ["iron_sulfur", "electron_transport"],  # halobenzoate dioxygenase electron transfer

    # Radical SAM enzymes
    "K20575": ["radical_sam", "iron_sulfur", "lyase"],  # radical SAM diol-dehydratase
    "K20593": ["radical_sam", "iron_sulfur", "cobalamin_binding", "methyltransferase"],  # cobalamin-dependent radical SAM methyltransferase

    # Miscellaneous enzymes
    "K17048": ["oxidoreductase", "iron_sulfur"],  # ethylbenzene hydroxylase beta subunit
    "K20439": ["phosphatase", "isomerase"],  # phosphohexomutase/phosphatase
    "K19632": ["oxidoreductase", "carbohydrate_active"],  # CDP-abequose synthase
    "K27986": ["transferase", "secondary_metabolism"],  # branched-chain oxoacid:malonyl-ACP transferase
    "K19884": ["oxidoreductase", "secondary_metabolism"],  # chloro-L-tryptophan oxidase
    "K13549": ["oxidoreductase", "dehydrogenase"],  # 2-deoxy-scyllo-inosamine dehydrogenase
    "K27193": ["esterase", "hydrolase"],  # carbendazim hydrolase
    "K20586": ["oxidoreductase"],  # oxidoreductase
    "K17055": ["lyase", "terpene_synthesis"],  # eugenol synthase
    # WARNING: K05914 (luciferase) is often mis-assigned to NRPS proteins by KofamScan
    # due to AMP-binding domain similarity. If a protein has both K05914 and NRPS annotations,
    # the bioluminescence assignment is likely a FALSE POSITIVE.
    "K05914": ["oxidoreductase", "bioluminescence"],  # photinus luciferin monooxygenase - USE WITH CAUTION
    "K15635": ["isomerase", "glycolysis"],  # phosphoglycerate mutase (independent)

    # Transcription factors
    "K03124": ["transcription", "dna_binding"],  # transcription initiation factor TFIIB
    "K03120": ["transcription", "dna_binding"],  # TATA-box-binding protein (TBP)

    # DNA replication and repair
    "K02684": ["primase", "replication"],  # DNA primase small subunit
    "K04799": ["nuclease", "dna_repair", "base_excision_repair"],  # flap endonuclease-1
    "K01520": ["hydrolase", "nucleotide_metabolism"],  # dUTP diphosphatase

    # CRISPR-Cas proteins (direct mapping for missing definitions)
    "K19091": ["crispr_associated", "cas_domain", "nuclease", "defense_system"],  # Cas6 endoribonuclease
    "K19144": ["crispr_associated", "cas_domain", "nuclease", "defense_system"],  # Csx3 (CRISPR-associated)

    # ATPases
    "K03924": ["atpase", "atp_binding", "chaperone"],  # MoxR-like ATPase

    # Transposases
    "K07496": ["transposase", "mobile_element"],  # putative transposase

    # Asparagine synthase
    "K01953": ["ligase", "amino_acid_biosynthesis", "atp_binding"],  # asparagine synthase (glutamine-hydrolysing)

    # Pyruvate formate lyase
    "K04069": ["radical_sam", "iron_sulfur", "glycolysis", "fermentation"],  # pyruvate formate lyase activating enzyme

    # Bromoperoxidase
    "K05918": ["oxidoreductase", "halogenase"],  # 4-hydroxybenzoate brominase

    # Uncharacterized proteins
    "K07133": ["hypothetical"],  # uncharacterized protein

    # Viral proteins (may appear in metagenomes)
    "K26370": ["rna_polymerase", "transcription", "phage_related"],  # Iridovirus RNA polymerase
    "K21664": ["phage_related", "dna_binding"],  # KSHV latency-associated nuclear antigen

    # Ubiquitin system
    "K15343": ["ubiquitin_ligase", "protein_modification"],  # E3 ubiquitin-protein ligase
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
    (r"ABC.*transporter|ATP-binding cassette", ["transporter", "abc_transporter", "atp_binding"]),
    (r"\bchannel\b", ["transporter", "ion_channel"]),
    (r"\bsymporter\b", ["transporter", "symporter"]),
    (r"\bantiporter\b", ["transporter", "antiporter"]),
    (r"efflux", ["transporter", "efflux_pump"]),

    # Regulators
    (r"transcription.*regulator|regulator.*transcription", ["regulator", "transcription_factor"]),
    (r"response.*regulator", ["response_regulator", "two_component"]),
    (r"sensor.*kinase|histidine.*kinase", ["sensor_kinase", "two_component"]),
    (r"sigma.*factor", ["sigma_factor"]),
    (r"transcription.*factor|TFIIB|TBP", ["transcription", "dna_binding"]),

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
    (r"\bchaperone\b", ["chaperone", "stress_response"]),

    # Mobile elements
    (r"\btransposase\b", ["transposase", "mobile_element"]),
    (r"\bintegrase\b", ["integrase", "mobile_element"]),
    (r"\bphage\b", ["phage_related"]),
    (r"recombinase", ["recombinase", "mobile_element"]),

    # Signaling
    (r"circadian|clock.*protein|KaiC", ["signaling"]),
    (r"ATPase", ["atpase", "atp_binding"]),
    (r"GTPase", ["gtpase", "gtp_binding"]),

    # Cofactor/coenzyme related
    (r"radical.*SAM", ["radical_sam", "iron_sulfur"]),
    (r"cobalamin|B12", ["cobalamin_binding", "cobalt_binding"]),
    (r"ferredoxin", ["iron_sulfur", "electron_transport"]),
    (r"flavin|FAD|FMN", ["fad_binding"]),
    (r"NAD|NADP|NADH", ["nad_binding"]),
    (r"PLP|pyridoxal", ["plp_binding"]),

    # Defense systems
    (r"CRISPR|Cas\d+", ["crispr_associated", "defense_system"]),
    (r"restriction|methylase.*DNA|DNA.*methylase", ["restriction_modification"]),

    # Cell division
    (r"\bFts[AZWLNQKX]\b", ["cell_division", "divisome"]),
    (r"\bPar[AB]\b", ["cell_division", "chromosome_partitioning"]),

    # Stress response
    (r"heat.*shock|cold.*shock", ["stress_response"]),
    (r"superoxide.*dismutase|catalase|peroxiredoxin", ["oxidative_stress"]),
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


__all__ = [
    "EC_TO_PREDICATES",
    "KEGG_TO_PREDICATES",
    "KEGG_PATTERNS",
    "parse_ec_numbers",
    "get_predicates_for_ec",
    "get_predicates_for_kegg",
]
