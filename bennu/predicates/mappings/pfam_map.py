"""
PFAM domain to predicate mappings.

Maps PFAM accessions to semantic predicates. Uses both:
1. Direct mappings for specific domains
2. Pattern-based matching on domain names/descriptions
"""

from typing import List, Tuple
from pathlib import Path
import re

# ============================================================================
# DIRECT PFAM TO PREDICATE MAPPINGS
# ============================================================================
# Format: "PFAM_ACC": ["predicate1", "predicate2", ...]

PFAM_TO_PREDICATES: dict[str, list[str]] = {
    # -------------------------------------------------------------------------
    # TRANSPORTERS
    # -------------------------------------------------------------------------
    # ABC transporters
    "PF00005": ["transporter", "abc_transporter", "abc_atpase", "atp_binding"],  # ABC_tran
    "PF00664": ["transporter", "abc_transporter", "abc_atpase", "atp_binding"],  # ABC_membrane
    "PF01061": ["transporter", "abc_transporter", "abc_atpase", "atp_binding"],  # ABC2_membrane
    "PF00528": ["transporter", "abc_transporter", "abc_substrate_binding", "periplasmic"],  # BPD_transp_1
    "PF01547": ["transporter", "abc_transporter", "abc_substrate_binding"],  # SBP_bac_1
    "PF01497": ["transporter", "abc_transporter", "abc_permease"],  # Peripla_BP_2
    "PF02470": ["transporter", "abc_transporter", "abc_permease"],  # PBP2_ABC
    "PF00950": ["transporter", "abc_transporter", "abc_atpase"],  # ABC_ATPase
    "PF06472": ["transporter", "abc_transporter", "abc_permease", "sugar_transporter"],  # ABC_tran_Xtn
    "PF02653": ["transporter", "abc_transporter", "abc_permease", "peptide_transporter"],  # BPD_transp_2

    # MFS transporters
    "PF07690": ["transporter", "mfs_transporter", "membrane"],  # MFS_1
    "PF00083": ["transporter", "mfs_transporter", "sugar_transporter", "membrane"],  # Sugar_tr
    "PF13347": ["transporter", "mfs_transporter", "membrane"],  # MFS_2
    "PF05977": ["transporter", "mfs_transporter", "membrane"],  # MFS_3

    # Ion channels and transporters
    "PF00520": ["transporter", "ion_channel", "membrane"],  # Ion_trans
    "PF07885": ["transporter", "ion_channel", "membrane"],  # Ion_trans_2
    "PF01699": ["transporter", "ion_transporter", "sodium_transporter", "membrane"],  # Na_H_Exchanger
    "PF00999": ["transporter", "ion_transporter", "sodium_transporter"],  # Na_H_antiport_1

    # 7TM receptors
    "PF00001": ["membrane"],  # 7tm_1 - 7 transmembrane receptor (rhodopsin family)
    "PF00002": ["membrane"],  # 7tm_2 - 7 transmembrane receptor (Secretin family)
    "PF00003": ["membrane"],  # 7tm_3 - 7 transmembrane sweet-taste receptor of 3 GCPR

    # Porins
    "PF00267": ["transporter", "porin", "outer_membrane", "beta_barrel"],  # Porin_1
    "PF02264": ["transporter", "porin", "outer_membrane", "beta_barrel"],  # Porin_2
    "PF13609": ["transporter", "porin", "outer_membrane", "beta_barrel"],  # Porin_4

    # Efflux pumps
    "PF00873": ["transporter", "efflux_pump", "membrane", "antibiotic_resistance"],  # ACR_tran
    "PF01554": ["transporter", "efflux_pump", "membrane"],  # MatE
    "PF02321": ["transporter", "efflux_pump", "membrane"],  # OEP

    # Amino acid transporters
    "PF01490": ["transporter", "amino_acid_transporter", "membrane"],  # AA_permease
    "PF00324": ["transporter", "amino_acid_transporter", "membrane"],  # AA_permease_2
    "PF13520": ["transporter", "amino_acid_transporter", "membrane"],  # AA_permease_N

    # Metal transporters
    "PF01566": ["transporter", "metal_transporter", "membrane"],  # Nramp
    "PF03547": ["transporter", "metal_transporter", "membrane"],  # MgtE
    "PF02535": ["transporter", "metal_transporter", "zinc_binding"],  # ZIP
    "PF01545": ["transporter", "metal_transporter", "copper_binding", "membrane"],  # Cation_efflux
    "PF00403": ["transporter", "metal_transporter", "heavy_metal_resistance"],  # HMA

    # Phosphate/sulfate transporters
    "PF01384": ["transporter", "phosphate_transporter", "membrane"],  # PHO4
    "PF00939": ["transporter", "sulfate_transporter", "membrane"],  # Na_sulph_symp
    "PF01740": ["transporter", "sulfate_transporter"],  # STAS

    # -------------------------------------------------------------------------
    # REGULATORS
    # -------------------------------------------------------------------------
    # Response regulators and two-component
    "PF00072": ["regulator", "response_regulator", "two_component"],  # Response_reg
    "PF00512": ["regulator", "sensor_kinase", "two_component", "kinase"],  # HisKA
    "PF02518": ["regulator", "sensor_kinase", "two_component", "atp_binding"],  # HATPase_c
    "PF00672": ["regulator", "sensor_kinase", "two_component"],  # HAMP
    "PF07730": ["regulator", "sensor_kinase", "two_component"],  # HisKA_3
    "PF13581": ["regulator", "response_regulator", "ompr_family"],  # HTH_27

    # LysR family
    "PF00126": ["regulator", "transcription_factor", "lysr_family", "dna_binding", "helix_turn_helix"],  # HTH_1
    "PF03466": ["regulator", "lysr_family"],  # LysR_substrate

    # TetR family
    "PF00440": ["regulator", "transcription_factor", "tetr_family", "dna_binding", "helix_turn_helix"],  # TetR_N
    "PF08361": ["regulator", "tetr_family"],  # TetR_C (C-terminal)

    # GntR family
    "PF00392": ["regulator", "transcription_factor", "gntr_family", "dna_binding", "helix_turn_helix"],  # GntR

    # AraC family
    "PF00165": ["regulator", "transcription_factor", "arac_family", "dna_binding", "helix_turn_helix"],  # HTH_AraC
    "PF12833": ["regulator", "arac_family"],  # HTH_18

    # LuxR family
    "PF00196": ["regulator", "transcription_factor", "luxr_family", "dna_binding", "helix_turn_helix"],  # GerE

    # MarR family
    "PF01047": ["regulator", "transcription_factor", "marr_family", "dna_binding"],  # MarR
    "PF12802": ["regulator", "marr_family"],  # MarR_2

    # LacI family
    "PF00356": ["regulator", "transcription_factor", "laci_family", "dna_binding"],  # LacI
    "PF13377": ["regulator", "laci_family"],  # Peripla_BP_3

    # IclR family
    "PF01614": ["regulator", "transcription_factor", "iclr_family", "dna_binding"],  # IclR

    # CRP/FNR family
    "PF00027": ["regulator", "transcription_factor", "crp_fnr_family", "dna_binding"],  # cNMP_binding
    "PF13545": ["regulator", "crp_fnr_family", "dna_binding", "helix_turn_helix"],  # Crp

    # Sigma factors
    "PF00140": ["sigma_factor", "regulator"],  # Sigma70_r1_2
    "PF04542": ["sigma_factor", "regulator"],  # Sigma70_r2
    "PF04539": ["sigma_factor", "regulator"],  # Sigma70_r3
    "PF04545": ["sigma_factor", "regulator"],  # Sigma70_r4
    "PF08281": ["sigma_factor", "regulator"],  # Sigma70_r4_2

    # -------------------------------------------------------------------------
    # SIGNALING DOMAINS
    # -------------------------------------------------------------------------
    "PF00989": ["signaling", "sensor_kinase"],  # PAS
    "PF08447": ["signaling", "sensor_kinase"],  # PAS_3
    "PF08448": ["signaling", "sensor_kinase"],  # PAS_4
    "PF13426": ["signaling", "sensor_kinase"],  # PAS_9
    "PF01590": ["signaling", "sensor_kinase"],  # GAF
    "PF13185": ["signaling", "sensor_kinase"],  # GAF_2
    "PF00990": ["signaling", "cyclic_dinucleotide", "diguanylate_cyclase"],  # GGDEF
    "PF00563": ["signaling", "cyclic_dinucleotide", "phosphodiesterase"],  # EAL
    "PF01966": ["signaling", "cyclic_dinucleotide", "phosphodiesterase"],  # HD
    "PF00571": ["regulatory", "cbs_domain", "nucleotide_binding"],  # CBS - regulatory domain
    "PF04185": ["signaling", "phosphodiesterase"],  # HD_assoc
    "PF01794": ["chemotaxis", "signaling"],  # Fer4_7

    # Chemotaxis
    "PF00015": ["chemotaxis", "signaling", "membrane"],  # MCPsignal

    # -------------------------------------------------------------------------
    # ENZYMES - OXIDOREDUCTASES
    # -------------------------------------------------------------------------
    # Dehydrogenases
    "PF00106": ["oxidoreductase", "dehydrogenase", "nad_binding", "cofactor_binding"],  # adh_short
    "PF00107": ["oxidoreductase", "dehydrogenase", "nad_binding", "metal_binding", "zinc_binding"],  # ADH_zinc_N
    "PF08240": ["oxidoreductase", "dehydrogenase", "nad_binding"],  # ADH_N
    "PF00465": ["oxidoreductase", "dehydrogenase", "iron_sulfur"],  # Fe-ADH
    "PF01073": ["oxidoreductase", "dehydrogenase", "fad_binding"],  # 3HCDH_N
    "PF02737": ["oxidoreductase", "dehydrogenase", "fad_binding"],  # 3HCDH
    "PF00389": ["oxidoreductase", "dehydrogenase", "fad_binding", "2fe2s"],  # 2-Hacid_dh_C
    "PF02826": ["oxidoreductase", "dehydrogenase", "fad_binding"],  # D-isomer_dh

    # Reductases
    "PF00175": ["oxidoreductase", "reductase", "nad_binding"],  # NAD_binding_1
    "PF13450": ["oxidoreductase", "reductase", "nad_binding"],  # NAD_binding_8
    "PF00724": ["oxidoreductase", "reductase", "nadp_binding"],  # NADP_Rossmann
    "PF01070": ["oxidoreductase", "reductase", "fmn_binding"],  # FMN_red
    "PF01512": ["oxidoreductase", "reductase", "fad_binding"],  # Flavin_Reduct

    # Oxidases and oxygenases
    "PF00067": ["oxidoreductase", "oxygenase", "monooxygenase", "heme_binding"],  # p450
    "PF01494": ["oxidoreductase", "oxidase", "fad_binding"],  # FAD_binding_3
    "PF01266": ["oxidoreductase", "oxidase", "fad_binding"],  # DAO
    "PF01593": ["oxidoreductase", "oxidase", "flavin_binding"],  # Amino_oxidase
    "PF04261": ["oxidoreductase", "dioxygenase", "iron_binding"],  # Dioxygenase_C
    "PF00903": ["oxidoreductase", "dioxygenase", "iron_binding"],  # Glyoxalase

    # Peroxidases and catalases
    "PF00141": ["oxidoreductase", "peroxidase", "heme_binding", "oxidative_stress"],  # Peroxidase
    "PF00199": ["oxidative_stress", "catalase", "heme_binding"],  # Catalase
    "PF06628": ["oxidative_stress", "catalase"],  # Catalase-rel

    # Thioredoxins and glutaredoxins
    "PF00085": ["oxidoreductase", "thioredoxin", "oxidative_stress"],  # Thioredoxin
    "PF00462": ["oxidoreductase", "glutaredoxin", "oxidative_stress"],  # Glutaredoxin
    "PF13098": ["oxidoreductase", "thioredoxin"],  # Thioredoxin_2

    # PQQ-dependent enzymes
    "PF01011": ["oxidoreductase", "pqq_binding", "dehydrogenase"],  # PQQ - glucose/alcohol dehydrogenase
    "PF13360": ["oxidoreductase", "pqq_binding"],  # PQQ_2 - PQQ-like repeat
    "PF13570": ["oxidoreductase", "pqq_binding"],  # PQQ_3 - PQQ enzyme repeat

    # -------------------------------------------------------------------------
    # HYDROGENASES AND RELATED
    # -------------------------------------------------------------------------
    # NOTE: Hydrogenase subtyping by PFAM is approximate. HydDB homology search
    # provides better classification. These are best-effort assignments:
    # - Group 1 [NiFe]: Uptake hydrogenases (H2 oxidation)
    # - Group 3 [NiFe]: Bidirectional, F420/NAD-coupled
    # - Group 4 [NiFe]: Energy-conserving, H2-evolving (Ech, Mbh)
    # - [FeFe]: Usually H2 production
    # - [NiFeSe]: Selenium-containing, often bidirectional

    # [NiFeSe] hydrogenases - selenium-containing, often bidirectional
    "PF00374": ["hydrogenase", "nife_hydrogenase", "nifese_hydrogenase", "hydrogen_metabolism", "nickel_binding", "selenium"],  # NiFeSe_Hases

    # [FeFe] hydrogenases - typically H2-producing
    "PF02256": ["hydrogenase", "fefe_hydrogenase", "hydrogen_metabolism", "iron_sulfur"],  # Fe_hyd_SSU
    "PF02906": ["hydrogenase", "fefe_hydrogenase", "hydrogen_metabolism"],  # Fe_hyd_lg_C

    # [NiFe] - generic (subtype determined by KEGG or HydDB)
    "PF14720": ["hydrogenase", "nife_hydrogenase", "hydrogen_metabolism"],  # NiFe_hyd_SSU_C

    # Group 4 - Energy-conserving hydrogenases (H2-evolving)
    "PF12234": ["hydrogenase", "nife_hydrogenase", "hydrogenase_group4", "hydrogen_metabolism"],  # Ech_Ehr_A
    "PF13244": ["hydrogenase", "nife_hydrogenase", "hydrogenase_group4", "mbh_hydrogenase", "membrane"],  # MbhD
    "PF20501": ["hydrogenase", "nife_hydrogenase", "hydrogenase_group4", "mbh_hydrogenase", "membrane"],  # MbhE

    # Group 3 - F420-reducing hydrogenases (bidirectional)
    "PF04422": ["hydrogenase", "nife_hydrogenase", "hydrogenase_group3", "f420_dependent"],  # FrhB_FdhB_N
    "PF04432": ["hydrogenase", "nife_hydrogenase", "hydrogenase_group3", "f420_dependent"],  # FrhB_FdhB_C

    # Hydrogenase maturation
    "PF01155": ["hydrogenase_maturation", "nickel_binding"],  # HypA
    "PF01750": ["hydrogenase_maturation", "protease"],  # HycI
    "PF01924": ["hydrogenase_maturation"],  # HypD

    # NOTE: PF00142 (Fer4_NifH) is defined later as generic iron_sulfur - it's NOT hydrogenase-specific

    # CO dehydrogenase
    "PF02552": ["co_oxidation", "oxidoreductase", "molybdenum_binding"],  # CO_deh_flav_C
    "PF03598": ["co_oxidation", "oxidoreductase"],  # CO_dehydrog
    "PF03599": ["co_oxidation", "oxidoreductase", "iron_sulfur"],  # CO_deh_flav_N
    "PF18537": ["co_oxidation", "wood_ljungdahl"],  # CODH_C

    # -------------------------------------------------------------------------
    # NITROGEN METABOLISM
    # -------------------------------------------------------------------------
    # NOTE: nitrogen_fixation predicate should ONLY come from KEGG nifHDK orthologs.
    # These PFAM domains are too broad - they're found in many non-nitrogenase proteins:
    # - PF00148 (Oxidored_molyb) - generic molybdopterin domain
    # - PF02579/PF00142 (Fer4_NifH) - Fe-S domains shared with hydrogenases
    # - PF01656 (CbiX) - cobalamin biosynthesis, NOT nitrogen fixation
    "PF00148": ["molybdenum_binding", "oxidoreductase"],  # Oxidored_molyb - generic molybdopterin
    "PF02579": ["iron_sulfur"],  # Fer4_NifH - Fe-S cluster, found in many proteins
    "PF00142": ["iron_sulfur"],  # Fer4_NifH variant - shared with hydrogenase
    "PF01656": ["iron_sulfur", "cobalamin_biosynthesis"],  # CbiX - cobalamin synthesis, NOT nif!
    "PF04055": ["radical_sam", "iron_sulfur", "enzyme"],  # Radical_SAM - versatile enzyme superfamily

    # Nitrate/nitrite reductases
    "PF00384": ["denitrification", "nitrate_reduction", "molybdenum_binding"],  # Molybdopterin
    "PF01568": ["denitrification", "nitrate_reduction"],  # Molydop_binding
    "PF02613": ["denitrification", "nitrite_reductase", "heme_binding"],  # Nitrite_red_cu
    "PF00115": ["respiration", "heme_binding"],  # COX1 (cytochrome oxidase)
    "PF00394": ["denitrification", "copper_binding"],  # Cu-oxidase

    # Ammonia metabolism
    "PF00120": ["nitrogen_metabolism", "ammonia_assimilation", "glutamine_synthetase"],  # Gln-synt_C
    "PF03951": ["nitrogen_metabolism", "ammonia_assimilation"],  # Gln-synt_N
    "PF00310": ["nitrogen_metabolism", "aminotransferase"],  # GATase_2

    # -------------------------------------------------------------------------
    # SULFUR METABOLISM
    # -------------------------------------------------------------------------
    "PF00890": ["sulfate_reduction", "sulfur_metabolism", "iron_sulfur"],  # DsrA (dissimilatory sulfite reductase)
    "PF04358": ["sulfate_reduction", "sulfur_metabolism"],  # DsrC
    "PF01087": ["sulfur_assimilation", "sulfur_metabolism"],  # Sir_G
    "PF12139": ["sulfur_metabolism", "iron_sulfur"],  # Fe4S4_Ni
    "PF02910": ["sulfur_metabolism", "thiosulfate"],  # Rhodanese

    # -------------------------------------------------------------------------
    # METHANOGENESIS
    # -------------------------------------------------------------------------
    # NOTE: methanogenesis predicate is STRICT - only MCR subunits (the definitive marker).
    # Other "methanogen" enzymes (F420, tetrahydromethanopterin) are found in non-methanogens
    # and should use archaeal_one_carbon instead.
    "PF02240": ["methanogenesis", "mcr_complex", "one_carbon_metabolism", "nickel_binding"],  # MCR_alpha - DEFINITIVE
    "PF02241": ["methanogenesis", "mcr_complex", "one_carbon_metabolism"],  # MCR_beta - DEFINITIVE
    "PF02249": ["methanogenesis", "mcr_complex", "one_carbon_metabolism"],  # MCR_gamma - DEFINITIVE
    "PF02007": ["archaeal_one_carbon", "one_carbon_metabolism"],  # MtrH - found in non-methanogens too
    "PF02512": ["cobalt_binding", "cobalamin_biosynthesis"],  # CobW_C - cobalamin, not methanogenesis-specific
    "PF00936": ["methanotrophy", "one_carbon_metabolism", "copper_binding"],  # MMO_alpha

    # -------------------------------------------------------------------------
    # PHOTOSYNTHESIS AND CARBON FIXATION
    # -------------------------------------------------------------------------
    # NOTE: RuBisCO alone does NOT mean Calvin cycle - RuBisCO-like proteins (RLPs)
    # are involved in nucleotide salvage. Calvin cycle requires RuBisCO + PRK.
    "PF00016": ["rubisco", "carbon_fixation"],  # RuBisCO_large - NOT calvin_cycle without PRK
    "PF00101": ["rubisco", "carbon_fixation"],  # RuBisCO_small
    "PF00485": ["prk", "calvin_cycle", "carbon_fixation", "kinase"],  # PRK - phosphoribulokinase, true Calvin marker
    "PF00223": ["photosynthesis", "photosystem_i", "iron_sulfur"],  # PsaA_PsaB
    "PF00421": ["photosynthesis", "photosystem_ii", "heme_binding"],  # PSII
    "PF00124": ["photosynthesis", "light_harvesting"],  # Photo_RC
    "PF00556": ["photosynthesis", "light_harvesting"],  # Antenna_bact
    "PF00032": ["respiration", "heme_binding"],  # Cytochrom_B_C
    "PF00033": ["heme_binding", "electron_transport"],  # Cytochrome_B - Cytochrome b/b6/petB
    "PF00034": ["heme_binding", "electron_transport"],  # Cytochrom_C - Cytochrome c
    "PF00042": ["heme_binding"],  # Globin

    # -------------------------------------------------------------------------
    # CENTRAL CARBON METABOLISM
    # -------------------------------------------------------------------------
    # Glycolysis
    "PF00044": ["glycolysis", "central_metabolism", "nad_binding"],  # Gp_dh_N (GAPDH)
    "PF02800": ["glycolysis", "central_metabolism"],  # Gp_dh_C
    "PF00113": ["glycolysis", "gluconeogenesis", "kinase"],  # Enolase_N
    "PF03952": ["glycolysis", "gluconeogenesis"],  # Enolase_C
    "PF00224": ["glycolysis", "kinase", "atp_binding"],  # PK (pyruvate kinase)
    "PF02887": ["glycolysis", "kinase", "atp_binding"],  # PK_C
    "PF00365": ["glycolysis", "isomerase"],  # PfkB (phosphofructokinase)
    "PF00162": ["glycolysis", "kinase"],  # Phosphofructokinase

    # TCA cycle
    "PF00549": ["tca_cycle", "central_metabolism", "ligase"],  # Ligase_CoA
    "PF02629": ["tca_cycle", "central_metabolism"],  # CoA_binding
    "PF00285": ["tca_cycle", "central_metabolism", "lyase"],  # Citrate_synt
    "PF02779": ["tca_cycle", "central_metabolism", "translocase"],  # 2-oxoacid_dh

    # Pentose phosphate
    "PF00479": ["pentose_phosphate", "central_metabolism"],  # G6PD_N
    "PF02781": ["pentose_phosphate", "central_metabolism"],  # G6PD_C
    "PF00378": ["pentose_phosphate", "central_metabolism"],  # ECH_1

    # -------------------------------------------------------------------------
    # ATP SYNTHESIS
    # -------------------------------------------------------------------------
    "PF00006": ["atp_synthesis", "energy_metabolism", "atp_binding"],  # ATP-synt_ab
    "PF02874": ["atp_synthesis", "energy_metabolism"],  # ATP-synt_ab_N
    "PF00213": ["atp_synthesis", "energy_metabolism"],  # OSCP
    "PF00137": ["atp_synthesis", "energy_metabolism", "membrane"],  # ATP-synt_C
    "PF05496": ["atp_synthesis", "energy_metabolism"],  # RuvB_N

    # -------------------------------------------------------------------------
    # AMINO ACID METABOLISM
    # -------------------------------------------------------------------------
    "PF00155": ["aminotransferase", "amino_acid_metabolism", "plp_binding"],  # Aminotran_1_2
    "PF00202": ["aminotransferase", "amino_acid_metabolism", "plp_binding"],  # Aminotran_3
    "PF01063": ["amino_acid_biosynthesis", "aromatic_aa_metabolism"],  # Aminotran_4
    "PF00291": ["amino_acid_biosynthesis", "aromatic_aa_metabolism"],  # DAHP_synth_1

    # -------------------------------------------------------------------------
    # LIPID METABOLISM
    # -------------------------------------------------------------------------
    "PF00109": ["fatty_acid_synthesis", "lipid_metabolism", "transferase"],  # Ketoacyl-synt
    "PF02801": ["fatty_acid_synthesis", "lipid_metabolism"],  # Ketoacyl-synt_C
    "PF00550": ["fatty_acid_synthesis", "lipid_metabolism"],  # PP-binding (ACP)
    "PF00698": ["fatty_acid_synthesis", "lipid_metabolism"],  # Acyl_transf_1
    "PF00173": ["fatty_acid_synthesis", "lipid_metabolism"],  # Cyt_b5

    "PF00441": ["fatty_acid_degradation", "lipid_metabolism", "beta_oxidation"],  # Acyl-CoA_dh_N
    "PF02770": ["fatty_acid_degradation", "lipid_metabolism", "beta_oxidation"],  # Acyl-CoA_dh_M
    "PF02771": ["fatty_acid_degradation", "lipid_metabolism", "beta_oxidation"],  # Acyl-CoA_dh_C
    "PF00378": ["fatty_acid_degradation", "lipid_metabolism"],  # ECH_1 (enoyl-CoA hydratase)

    # -------------------------------------------------------------------------
    # CELL WALL / PEPTIDOGLYCAN
    # -------------------------------------------------------------------------
    "PF00912": ["peptidoglycan", "murein_synthesis", "cell_wall", "transferase"],  # Transglycosyl
    "PF01476": ["peptidoglycan", "murein_synthesis", "cell_wall"],  # LysM
    "PF00768": ["peptidoglycan", "pbp", "cell_wall", "transpeptidase"],  # Transpeptidase
    "PF00905": ["peptidoglycan", "pbp", "cell_wall"],  # Penicillin_BP
    "PF07943": ["peptidoglycan", "lytic_transglycosylase", "cell_wall"],  # LT_GEWL
    "PF01464": ["peptidoglycan", "lytic_transglycosylase", "cell_wall", "hydrolase"],  # SLT

    # -------------------------------------------------------------------------
    # SURFACE STRUCTURES
    # -------------------------------------------------------------------------
    # Flagella
    "PF00669": ["flagellum", "flagellin", "cell_surface"],  # Flagellin_N
    "PF00700": ["flagellum", "flagellin", "cell_surface"],  # Flagellin_C
    "PF01514": ["flagellum", "flagellar_hook"],  # Flg_hook
    "PF02465": ["flagellum", "flagellar_motor", "membrane"],  # MotA_ExbB
    "PF01312": ["flagellum", "flagellar_motor"],  # FlhA

    # Pili
    "PF00114": ["pilus", "cell_surface"],  # Pilin
    "PF07963": ["pilus", "type_iv_pilus"],  # TadE
    "PF05157": ["pilus", "type_iv_pilus", "atpase"],  # T2SSE

    # Adhesins and cell surface proteins
    "PF05637": ["adhesin", "cell_surface", "repeat_domain"],  # Flg_new
    "PF07974": ["adhesin", "cell_surface", "repeat_domain"],  # EGF_2
    "PF13895": ["adhesin", "cell_surface", "repeat_domain"],  # Ig_3
    "PF01345": ["adhesin", "cell_surface", "repeat_domain"],  # DUF11 - archaeal/bacterial surface protein
    "PF00404": ["adhesin", "cell_surface", "dockerin"],  # Dockerin_1 - cellulosome/adhesion assembly
    "PF00963": ["adhesin", "cell_surface", "cohesin"],  # Cohesin - binds dockerin
    "PF00801": ["adhesin", "cell_surface", "repeat_domain"],  # PKD - Ig-like, cell adhesion
    "PF13385": ["adhesin", "cell_surface", "repeat_domain"],  # Laminin_G_3 - cell-matrix adhesion
    "PF02210": ["adhesin", "cell_surface", "repeat_domain"],  # Laminin_G_2
    "PF07705": ["adhesin", "cell_surface"],  # CARDB - cell adhesion related domain
    "PF05124": ["s_layer", "cell_surface"],  # S_layer_C - archaeal/bacterial S-layer
    "PF13620": ["cell_surface", "repeat_domain"],  # CarboxypepD_reg - surface regulatory domain
    "PF13517": ["cell_surface", "repeat_domain"],  # FG-GAP - cell surface repeat
    "PF00092": ["cell_surface", "repeat_domain"],  # VWA - von Willebrand factor, adhesion
    "PF05974": ["cell_surface", "repeat_domain"],  # TSP_3 - thrombospondin adhesin
    "PF02494": ["cell_surface", "repeat_domain"],  # HYR - cell adhesion domain
    "PF00354": ["cell_surface", "immune_related"],  # Pentaxin - pattern recognition
    "PF13229": ["cell_surface", "beta_helix", "repeat_domain"],  # Beta_helix - surface proteins

    # -------------------------------------------------------------------------
    # CHAPERONES
    # -------------------------------------------------------------------------
    "PF00012": ["chaperone", "hsp70", "stress_response", "atp_binding"],  # HSP70
    "PF00226": ["chaperone", "hsp60", "stress_response"],  # DnaJ
    "PF01321": ["chaperone", "stress_response"],  # CLP_N
    "PF00011": ["chaperone", "small_hsp", "stress_response"],  # HSP20
    "PF02518": ["chaperone", "stress_response", "atp_binding"],  # HATPase_c
    "PF00004": ["aaa_domain", "atp_binding", "atpase"],  # AAA
    "PF07724": ["aaa_domain", "atp_binding"],  # AAA_2
    "PF07728": ["aaa_domain", "atp_binding"],  # AAA_5
    "PF13173": ["aaa_domain", "atp_binding"],  # AAA_14
    "PF13304": ["aaa_domain", "atp_binding"],  # AAA_21
    "PF17862": ["aaa_domain", "atp_binding"],  # AAA_lid_3
    "PF02861": ["chaperone", "clp_protease"],  # Clp_N

    # -------------------------------------------------------------------------
    # DNA/RNA BINDING
    # -------------------------------------------------------------------------
    # DNA binding
    "PF01381": ["dna_binding", "helix_turn_helix"],  # HTH_3
    "PF12844": ["dna_binding", "helix_turn_helix"],  # HTH_19
    "PF13384": ["dna_binding", "helix_turn_helix"],  # HTH_23
    "PF00046": ["dna_binding", "helix_turn_helix"],  # Homeodomain
    "PF00010": ["dna_binding", "helix_loop_helix"],  # HLH
    "PF00096": ["dna_binding", "zinc_finger", "metal_binding"],  # zf-C2H2
    "PF13912": ["dna_binding", "zinc_finger"],  # zf-C2H2_6
    "PF00145": ["dna_binding"],  # DNA_methylase

    # RNA binding
    "PF00076": ["rna_binding"],  # RRM_1
    "PF00013": ["rna_binding"],  # KH_1
    "PF00035": ["rna_binding"],  # dsrm
    "PF01423": ["rna_binding"],  # LSM
    "PF00270": ["rna_binding", "helicase", "atp_binding"],  # DEAD

    # -------------------------------------------------------------------------
    # MOBILE ELEMENTS
    # -------------------------------------------------------------------------
    # Transposases
    "PF01526": ["transposase", "mobile_element"],  # DDE_1
    "PF01609": ["transposase", "mobile_element"],  # DDE_Tnp_IS1
    "PF01527": ["transposase", "mobile_element"],  # HTH_Tnp_Tc3_2
    "PF02371": ["transposase", "mobile_element"],  # Transposase_28
    "PF00665": ["integrase", "mobile_element", "recombinase"],  # rve
    "PF13683": ["integrase", "mobile_element"],  # rve_2
    "PF13358": ["integrase", "mobile_element", "dna_binding"],  # INT_C
    "PF00239": ["resolvase", "recombinase", "mobile_element", "dna_binding"],  # Resolvase

    # -------------------------------------------------------------------------
    # PHAGE-RELATED
    # -------------------------------------------------------------------------
    "PF03354": ["phage_related", "phage_terminase"],  # Terminase_1
    "PF04466": ["phage_related", "phage_terminase"],  # Terminase_3
    "PF05133": ["phage_related", "phage_portal"],  # Phage_portal
    "PF04860": ["phage_related", "phage_capsid"],  # Phage_capsid
    "PF05135": ["phage_related", "phage_tail"],  # Phage_tail_2
    "PF01520": ["phage_related", "holin", "membrane"],  # Phage_holin_1
    "PF04531": ["phage_related", "phage_lysin"],  # Phage_lysozyme

    # -------------------------------------------------------------------------
    # DEFENSE SYSTEMS
    # -------------------------------------------------------------------------
    # Restriction-modification
    "PF01420": ["restriction_modification", "methyltransferase_rm", "methyltransferase"],  # MethyltransfD12
    "PF00145": ["restriction_modification", "methyltransferase_rm", "dna_binding"],  # DNA_methylase
    "PF04851": ["restriction_modification", "restriction_enzyme", "nuclease"],  # ResIII
    "PF09019": ["restriction_modification", "restriction_enzyme"],  # HsdM_N

    # -------------------------------------------------------------------------
    # CRISPR-Cas COMPONENTS
    # -------------------------------------------------------------------------
    # NOTE: crispr_associated = KEGG-confirmed or CRISPR-specific domain
    #       cas_domain = has Cas-like domain (may be transposase or other)
    # System-level typing done by analyze_crispr_systems() operator

    # Unambiguous CRISPR markers - specific to CRISPR
    "PF09707": ["crispr_associated", "cas_domain", "defense_system"],  # Cas1 - adaptation module, CRISPR-specific
    "PF09481": ["crispr_associated", "cas_domain", "nuclease", "defense_system"],  # Cas2 - adaptation module, CRISPR-specific
    "PF09344": ["crispr_associated", "cas_domain", "defense_system"],  # Cas3_HD - Type I signature
    "PF01930": ["crispr_associated", "cas_domain", "defense_system"],  # Csh1 - CRISPR-specific

    # -------------------------------------------------------------------------
    # TOXIN-ANTITOXIN COMPONENTS
    # -------------------------------------------------------------------------
    # NOTE: These are DOMAIN-level tags. TA system calling requires pairing
    # analysis (toxin + antitoxin adjacent). Use detect_ta_systems() operator.
    # Many "toxin" domains have non-TA functions (PIN = tRNA processing, etc.)

    "PF04221": ["toxin_domain", "defense_system"],  # RelE - mRNA interferase (also ribosome rescue)
    "PF03693": ["antitoxin_domain", "defense_system"],  # YefM - antitoxin
    "PF01850": ["toxin_domain", "pin_domain", "ribonuclease"],  # PIN - also tRNA splicing
    "PF05973": ["toxin_domain", "defense_system"],  # ParE - DNA gyrase inhibitor
    "PF06414": ["toxin_domain", "defense_system"],  # Zeta_toxin
    "PF13560": ["toxin_domain", "defense_system"],  # HicA_toxin - mRNA interferase
    "PF07903": ["antitoxin_domain", "defense_system"],  # HicB_antitoxin
    "PF02604": ["antitoxin_domain", "defense_system"],  # AntA_antitoxin
    "PF15738": ["toxin_domain", "abortive_infection", "defense_system"],  # AbiE_toxin

    # -------------------------------------------------------------------------
    # SECRETION SYSTEM COMPONENTS
    # -------------------------------------------------------------------------
    # NOTE: These are COMPONENT-level tags. Locus-level system calling
    # (e.g., "this is a functional T6SS") requires clustering analysis via operator.
    # A single component doesn't mean the system is present/functional.

    "PF00437": ["secretion_component", "t1ss_component"],  # HlyD - Type I membrane fusion
    "PF00577": ["secretion_component", "sec_pathway", "signal_recognition"],  # Sec_GG
    "PF02416": ["secretion_component", "tat_pathway"],  # TatC
    "PF03798": ["secretion_component", "t3ss_component"],  # YscJ_FliF - T3SS basal body
    "PF05936": ["secretion_component", "t4ss_component", "conjugation"],  # VirB8
    "PF03743": ["secretion_component", "t6ss_component"],  # T6SS_IcmF

    # -------------------------------------------------------------------------
    # SECONDARY METABOLISM
    # -------------------------------------------------------------------------
    # PKS/NRPS
    "PF00668": ["secondary_metabolism", "polyketide_synthesis"],  # Condensation
    "PF00550": ["secondary_metabolism", "polyketide_synthesis", "nrps"],  # PP-binding
    "PF00501": ["secondary_metabolism", "nrps", "atp_binding"],  # AMP-binding

    # -------------------------------------------------------------------------
    # REPEAT DOMAINS
    # -------------------------------------------------------------------------
    "PF00515": ["repeat_domain", "tpr_repeat", "protein_binding"],  # TPR_1
    "PF07719": ["repeat_domain", "tpr_repeat"],  # TPR_2
    "PF13176": ["repeat_domain", "tpr_repeat"],  # TPR_7
    "PF13181": ["repeat_domain", "tpr_repeat"],  # TPR_8
    "PF13374": ["repeat_domain", "tpr_repeat"],  # TPR_10
    "PF13414": ["repeat_domain", "tpr_repeat"],  # TPR_11
    "PF13424": ["repeat_domain", "tpr_repeat"],  # TPR_12
    "PF13432": ["repeat_domain", "tpr_repeat"],  # TPR_16
    "PF13431": ["repeat_domain", "tpr_repeat"],  # TPR_17
    "PF14559": ["repeat_domain", "tpr_repeat"],  # TPR_19
    "PF02985": ["repeat_domain", "heat_repeat", "protein_binding"],  # HEAT
    "PF13646": ["repeat_domain", "heat_repeat"],  # HEAT_2
    "PF00400": ["repeat_domain", "wd40_repeat", "protein_binding"],  # WD40
    "PF00023": ["repeat_domain", "ankyrin_repeat", "protein_binding"],  # Ank
    "PF12796": ["repeat_domain", "ankyrin_repeat"],  # Ank_2
    "PF00560": ["repeat_domain", "lrr_repeat"],  # LRR_1
    "PF13855": ["repeat_domain", "lrr_repeat"],  # LRR_8
    "PF01344": ["repeat_domain", "kelch_repeat"],  # Kelch_1
    "PF08238": ["repeat_domain", "sel1_repeat"],  # Sel1

    # -------------------------------------------------------------------------
    # ANTIBIOTIC RESISTANCE
    # -------------------------------------------------------------------------
    "PF00144": ["antibiotic_resistance", "beta_lactamase", "hydrolase"],  # Beta-lactamase
    "PF00905": ["antibiotic_resistance", "pbp"],  # Penicillin_BP
    "PF01636": ["antibiotic_resistance", "aminoglycoside_resistance", "phosphotransferase"],  # APH
    "PF00903": ["antibiotic_resistance"],  # Glyoxalase
    "PF00873": ["antibiotic_resistance", "efflux_pump", "multidrug_resistance"],  # ACR_tran

    # -------------------------------------------------------------------------
    # HEAVY METAL RESISTANCE
    # -------------------------------------------------------------------------
    "PF00403": ["heavy_metal_resistance", "metal_transporter", "copper_binding"],  # HMA
    "PF03176": ["heavy_metal_resistance", "mercury_resistance"],  # MerT
    "PF00085": ["heavy_metal_resistance", "arsenic_resistance"],  # ArsC (thioredoxin fold)

    # -------------------------------------------------------------------------
    # METAL-BINDING (SPECIFIC METALS)
    # -------------------------------------------------------------------------
    # Nickel binding
    "PF00374": ["nickel_binding", "metal_binding"],  # NiFeSe_Hases (NiFe hydrogenase large subunit)
    "PF01155": ["nickel_binding", "metal_binding"],  # HypA (nickel incorporation)
    "PF00449": ["nickel_binding", "metal_binding", "urease"],  # Urease_alpha (urease active site)
    "PF01979": ["nickel_binding", "metal_binding"],  # Amidohydro_1 (urease family)
    "PF02492": ["nickel_binding", "metal_binding"],  # CobW (nickel/cobalt chaperone)
    "PF02240": ["nickel_binding", "metal_binding"],  # MCR_alpha (methyl-CoM reductase)
    "PF03063": ["nickel_binding", "metal_binding"],  # NikA_ABC (nickel ABC transporter)
    "PF04607": ["nickel_binding", "metal_binding"],  # RelA_SpoT (some have Ni coordination)

    # Molybdenum/tungsten binding (molybdopterin cofactor)
    "PF00384": ["molybdenum_binding", "metal_binding"],  # Molybdopterin (MoCo binding)
    "PF01568": ["molybdenum_binding", "metal_binding"],  # Molydop_binding (MoaD)
    "PF00994": ["molybdenum_binding", "metal_binding"],  # MoCF_biosynth (Mo cofactor carrier)
    "PF01315": ["molybdenum_binding", "metal_binding"],  # Ald_Xan_dh_C (aldehyde oxidase)
    "PF02738": ["molybdenum_binding", "metal_binding"],  # Ald_Xan_dh_C2 (xanthine dehydrogenase)
    "PF00174": ["molybdenum_binding", "metal_binding"],  # Oxidored_molyb (molybdopterin oxidoreductase)
    "PF01209": ["molybdenum_binding", "metal_binding"],  # Ubie_methyltran (ubiquinone biosynthesis)
    "PF02552": ["molybdenum_binding", "metal_binding"],  # CO_deh_flav_C (CO dehydrogenase)
    "PF03404": ["molybdenum_binding", "metal_binding"],  # Mo-co_dimer (Mo-co dimerization)
    "PF04879": ["molybdenum_binding", "metal_binding"],  # Molybdop_Fe4S4 (Mo-pt + Fe-S)

    # Copper binding
    "PF00127": ["copper_binding", "metal_binding"],  # Cu-oxidase (copper oxidase)
    "PF00394": ["copper_binding", "metal_binding"],  # Cu-oxidase_2 (multi-copper oxidase)
    "PF07731": ["copper_binding", "metal_binding"],  # Cu-oxidase_3 (laccase)
    "PF00924": ["copper_binding", "metal_binding"],  # MS_channel (mechanosensitive channel)
    "PF00936": ["copper_binding", "metal_binding"],  # BMC (methane monooxygenase)
    "PF02298": ["copper_binding", "metal_binding"],  # Cu_amine_oxid (copper amine oxidase)
    "PF01545": ["copper_binding", "metal_binding"],  # Cation_efflux (CDF family)
    "PF00403": ["copper_binding", "metal_binding"],  # HMA (heavy metal associated)
    "PF02978": ["copper_binding", "metal_binding"],  # Tyrosinase (di-copper center)
    "PF12706": ["copper_binding", "metal_binding"],  # Cupredoxin_1 (plastocyanin-like)
    "PF00355": ["copper_binding", "metal_binding"],  # Azurin (blue copper protein)
    "PF05048": ["copper_binding", "metal_binding", "metal_homeostasis"],  # NosD - copper insertion for NosZ
    "PF11617": ["copper_binding", "metal_binding"],  # Cu-binding_MopE - copper binding

    # Zinc binding (catalytic and structural)
    "PF00096": ["zinc_binding", "zinc_finger", "metal_binding"],  # zf-C2H2
    "PF00098": ["zinc_binding", "zinc_finger", "metal_binding"],  # zf-CCHC
    "PF13912": ["zinc_binding", "zinc_finger", "metal_binding"],  # zf-C2H2_6
    "PF00105": ["zinc_binding", "metal_binding"],  # zf-C4 (nuclear receptor)
    "PF01408": ["zinc_binding", "metal_binding"],  # GCN5L1 (zinc ribbon)
    "PF00107": ["zinc_binding", "metal_binding"],  # ADH_zinc_N (alcohol dehydrogenase)
    "PF00096": ["zinc_binding", "metal_binding"],  # zf-C2H2 (classical zinc finger)
    "PF00240": ["zinc_binding", "metal_binding"],  # ubiquitin (RING finger)
    "PF00097": ["zinc_binding", "metal_binding"],  # zf-C3HC4 (RING finger)
    "PF01909": ["zinc_binding", "metal_binding"],  # NTF2 (nuclear transport factor)
    "PF00136": ["zinc_binding", "metal_binding"],  # DNA_pol_B (DNA polymerase zinc finger)
    "PF01846": ["zinc_binding", "metal_binding"],  # FF (zinc-stabilized domain)
    "PF00653": ["zinc_binding", "metal_binding"],  # BIR (baculovirus inhibitor of apoptosis)
    "PF02622": ["zinc_binding", "metal_binding"],  # Zn_clus (zinc cluster)
    "PF00478": ["zinc_binding", "metal_binding"],  # IMPDH (IMP dehydrogenase)
    "PF01061": ["zinc_binding", "metal_binding"],  # ABC2_membrane (some have Zn site)

    # Cobalt/cobalamin (B12) binding
    "PF02310": ["cobalt_binding", "cobalamin_binding", "metal_binding"],  # B12-binding (adenosylcobalamin)
    "PF02607": ["cobalt_binding", "cobalamin_binding", "metal_binding"],  # B12-binding_2 (methylcobalamin)
    "PF02965": ["cobalt_binding", "cobalamin_binding", "metal_binding"],  # B12-bind_KAM (methylmalonyl-CoA mutase)
    "PF02512": ["cobalt_binding", "cobalamin_binding", "metal_binding"],  # CobW_C (cobalamin biosynthesis)
    "PF01903": ["cobalt_binding", "cobalamin_binding", "metal_binding"],  # CobD_Cbib (cobalamin biosynthesis)
    "PF02570": ["cobalt_binding", "metal_binding"],  # CobQ (cobyrinic acid synthase)
    "PF02492": ["cobalt_binding", "metal_binding"],  # CobW (cobalt chaperone)

    # Manganese binding
    "PF00080": ["manganese_binding", "metal_binding"],  # Sod_Fe_C (Mn/Fe superoxide dismutase)
    "PF02777": ["manganese_binding", "metal_binding"],  # Sod_Fe_N (Mn/Fe SOD N-terminal)
    "PF00149": ["manganese_binding", "metal_binding"],  # Metallophos (phosphatase, Mn-dependent)
    "PF01676": ["manganese_binding", "metal_binding"],  # Metalloenzyme (arginase-like)
    "PF00232": ["manganese_binding", "metal_binding"],  # Glyco_hydro_1 (some are Mn-dependent)

    # Calcium binding
    "PF00036": ["calcium_binding", "metal_binding"],  # EF-hand (EF-hand calcium-binding)
    "PF13499": ["calcium_binding", "metal_binding"],  # EF-hand_7 (EF-hand variant)
    "PF01023": ["calcium_binding", "metal_binding"],  # S-100 (calcium-binding protein)
    "PF00353": ["calcium_binding", "metal_binding"],  # RTX (RTX calcium-binding repeat)
    "PF01302": ["calcium_binding", "metal_binding"],  # CAP_GLY (calcium-binding repeat)
    "PF07699": ["calcium_binding", "metal_binding"],  # GCC2_GCC3 (calcium-binding)

    # Iron-sulfur clusters and ferredoxins
    "PF00037": ["iron_sulfur", "ferredoxin", "metal_binding", "electron_transport"],  # Fer4 (4Fe-4S binding)
    "PF13247": ["iron_sulfur", "ferredoxin", "metal_binding"],  # Fer4_7 (4Fe-4S ferredoxin)
    "PF12838": ["iron_sulfur", "ferredoxin", "metal_binding", "electron_transport"],  # Fer4_7 actual (4Fe-4S ferredoxin)
    "PF12800": ["iron_sulfur", "ferredoxin", "metal_binding", "electron_transport"],  # Fer4_4 (4Fe-4S cluster)
    "PF12837": ["iron_sulfur", "ferredoxin", "metal_binding", "electron_transport"],  # Fer4_6 (4Fe-4S cluster)
    "PF13237": ["iron_sulfur", "ferredoxin", "metal_binding", "electron_transport"],  # Fer4_10 (4Fe-4S dicluster)
    "PF13187": ["iron_sulfur", "ferredoxin", "metal_binding", "electron_transport"],  # Fer4_9 (4Fe-4S cluster)
    "PF13353": ["iron_sulfur", "ferredoxin", "metal_binding", "electron_transport"],  # Fer4_12 (4Fe-4S cluster)
    "PF12798": ["iron_sulfur", "ferredoxin", "metal_binding", "electron_transport"],  # Fer4_3 (4Fe-4S cluster)
    "PF12797": ["iron_sulfur", "ferredoxin", "metal_binding"],  # Fer4_2 (4Fe-4S ferredoxin-type)
    "PF00111": ["iron_sulfur", "ferredoxin", "metal_binding", "2fe2s"],  # Fer2 (2Fe-2S ferredoxin)
    "PF01077": ["iron_sulfur", "metal_binding"],  # ETF_alpha (electron transfer flavoprotein)
    "PF00210": ["iron_sulfur", "metal_binding"],  # Ferritin (iron storage, NOT ferredoxin)
    "PF13510": ["iron_sulfur", "ferredoxin", "metal_binding", "2fe2s"],  # Fer2_4 (2Fe-2S Rieske-type)
    "PF22117": ["iron_sulfur", "ferredoxin", "metal_binding", "electron_transport", "respiration"],  # Fer4_Nqo3 (Complex I ferredoxin)
    "PF13484": ["iron_sulfur", "ferredoxin", "metal_binding", "electron_transport"],  # Fer4_16 (4Fe-4S cluster)

    # Oxygen binding (non-heme iron)
    "PF01814": ["oxygen_binding", "iron_binding", "metal_binding"],  # Hemerythrin - O2 carrier/sensor
    "PF14720": ["oxygen_binding", "iron_binding"],  # Hemerythrin_2

    # -------------------------------------------------------------------------
    # GLYCOSYLTRANSFERASES
    # -------------------------------------------------------------------------
    "PF00535": ["glycosyltransferase", "transferase", "carbohydrate_active"],  # Glycos_transf_2
    "PF00534": ["glycosyltransferase", "transferase", "carbohydrate_active"],  # Glycos_transf_1
    "PF13692": ["glycosyltransferase", "transferase", "carbohydrate_active"],  # Glyco_trans_1_4
    "PF13439": ["glycosyltransferase", "transferase", "carbohydrate_active"],  # Glyco_transf_4
    "PF13641": ["glycosyltransferase", "transferase", "carbohydrate_active"],  # Glyco_tranf_2_3
    "PF13579": ["glycosyltransferase", "transferase", "carbohydrate_active"],  # Glyco_trans_4_4
    "PF10111": ["glycosyltransferase", "transferase", "carbohydrate_active"],  # Glyco_tranf_2_2

    # -------------------------------------------------------------------------
    # METHYLTRANSFERASES (non R-M)
    # -------------------------------------------------------------------------
    "PF08241": ["methyltransferase", "transferase", "sam_binding"],  # Methyltransf_11
    "PF13847": ["methyltransferase", "transferase", "sam_binding"],  # Methyltransf_31
    "PF13649": ["methyltransferase", "transferase", "sam_binding"],  # Methyltransf_25
    "PF13489": ["methyltransferase", "transferase", "sam_binding"],  # Methyltransf_23
    "PF08242": ["methyltransferase", "transferase", "sam_binding"],  # Methyltransf_12
    "PF01555": ["methyltransferase", "dna_methylase", "restriction_modification"],  # N6_N4_Mtase

    # -------------------------------------------------------------------------
    # DNA REPAIR
    # -------------------------------------------------------------------------
    "PF00154": ["dna_repair", "recombinational_repair", "atp_binding"],  # RecA
    "PF00270": ["dna_repair", "helicase", "atp_binding"],  # DEAD
    "PF06733": ["dna_repair", "base_excision_repair"],  # NUDIX
    "PF01419": ["dna_repair", "mismatch_repair"],  # 4HB_MCP_1
    "PF03167": ["dna_repair", "sos_response"],  # UmuC
    "PF18765": ["dna_repair", "base_excision_repair", "polymerase"],  # Polbeta - DNA polymerase beta

    # -------------------------------------------------------------------------
    # TRANSLATION
    # -------------------------------------------------------------------------
    # Ribosomal proteins
    "PF00252": ["ribosomal_protein", "translation"],  # Ribosomal_L16
    "PF00281": ["ribosomal_protein", "translation"],  # Ribosomal_L5
    "PF00297": ["ribosomal_protein", "translation"],  # Ribosomal_L6

    # tRNA synthetases
    "PF00152": ["trna_synthetase", "translation", "atp_binding", "amino_acid_metabolism"],  # tRNA-synt_2
    "PF01336": ["trna_synthetase", "translation", "amino_acid_metabolism"],  # OB_NTP_bind
    "PF01409": ["trna_synthetase", "translation"],  # tRNA-synt_2b

    # Translation factors
    "PF00009": ["translation_factor", "translation", "gtp_binding"],  # GTP_EFTU
    "PF03764": ["translation_factor", "translation"],  # EFG_IV
    "PF00679": ["translation_factor", "translation"],  # EFG_C

    # -------------------------------------------------------------------------
    # CELL DIVISION
    # -------------------------------------------------------------------------
    "PF00091": ["cell_division", "ftsz", "gtp_binding"],  # Tubulin
    "PF12327": ["cell_division", "ftsz"],  # FtsZ_C
    "PF01656": ["cell_division", "chromosome_partitioning"],  # CbiX (ParA-like)
    "PF03796": ["cell_division", "divisome"],  # DnaB_C

    # -------------------------------------------------------------------------
    # ADDITIONAL HIGH-ABUNDANCE DOMAINS (from dataset analysis)
    # -------------------------------------------------------------------------
    # AAA+ ATPase family domains
    "PF06745": ["atpase", "atp_binding", "signaling"],  # ATPase - KaiC circadian clock ATPase
    "PF13481": ["aaa_domain", "atp_binding", "atpase"],  # AAA_25
    "PF13614": ["aaa_domain", "atp_binding", "atpase"],  # AAA_31
    "PF13238": ["aaa_domain", "atp_binding", "atpase"],  # AAA_18
    "PF13175": ["aaa_domain", "atp_binding", "atpase"],  # AAA_15
    "PF13207": ["aaa_domain", "atp_binding", "atpase"],  # AAA_17
    "PF13401": ["aaa_domain", "atp_binding", "atpase"],  # AAA_22
    "PF01637": ["atpase", "atp_binding"],  # ATPase_2 - archaeal ATPase
    "PF17863": ["aaa_domain", "atp_binding"],  # AAA_lid_2

    # Epimerases and related NAD-dependent enzymes
    "PF01370": ["epimerase", "isomerase", "nad_binding", "carbohydrate_active"],  # Epimerase - NAD-dependent
    "PF16363": ["epimerase", "dehydrogenase", "nad_binding", "carbohydrate_active"],  # GDP_Man_Dehyd
    "PF04321": ["epimerase", "nad_binding", "carbohydrate_active"],  # RmlD_sub_bind - dTDP-sugar epimerase

    # HAD (haloacid dehalogenase) superfamily - phosphatases
    "PF00702": ["hydrolase", "phosphatase"],  # Hydrolase - HAD superfamily
    "PF13419": ["hydrolase", "phosphatase"],  # HAD_2
    "PF12710": ["hydrolase", "phosphatase"],  # HAD
    "PF08282": ["hydrolase", "phosphatase"],  # Hydrolase_3
    "PF13242": ["hydrolase", "phosphatase"],  # Hydrolase_like

    # tRNA synthetases - translation machinery
    "PF09334": ["trna_synthetase", "translation", "amino_acid_metabolism"],  # tRNA-synt_1g - Met-tRNA synthetase
    "PF00133": ["trna_synthetase", "translation", "amino_acid_metabolism"],  # tRNA-synt_1 - class I (I,L,M,V)
    "PF00587": ["trna_synthetase", "translation", "amino_acid_metabolism"],  # tRNA-synt_2b - class II core
    "PF03129": ["trna_synthetase", "translation", "rna_binding"],  # HGTP_anticodon - anticodon binding
    "PF08264": ["trna_synthetase", "translation", "rna_binding"],  # Anticodon_1 - anticodon binding

    # GTPases and ribosome-associated factors
    "PF01926": ["gtpase", "gtp_binding", "ribosomal_protein", "translation"],  # MMR_HSR1 - 50S ribosome-binding GTPase
    "PF03144": ["translation_factor", "translation", "gtp_binding"],  # GTP_EFTU_D2 - EF-Tu domain 2

    # Pyridine nucleotide-disulphide oxidoreductases
    "PF07992": ["oxidoreductase", "fad_binding", "nad_binding"],  # Pyr_redox_2
    "PF00070": ["oxidoreductase", "fad_binding", "nad_binding"],  # Pyr_redox
    "PF13738": ["oxidoreductase", "fad_binding", "nad_binding"],  # Pyr_redox_3

    # Metallo-beta-lactamase superfamily (diverse hydrolases)
    "PF00753": ["hydrolase", "metal_binding", "zinc_binding"],  # Lactamase_B - metallo-beta-lactamase fold

    # Helicase domains
    "PF00271": ["helicase", "atp_binding", "dna_repair"],  # Helicase_C - helicase C-terminal
    "PF13749": ["helicase", "atpase", "atp_binding"],  # HATPase_c_4 - DNA helicase recG

    # Transcription factors and DNA-binding
    "PF13412": ["dna_binding", "helix_turn_helix", "winged_helix", "regulator"],  # HTH_24 - winged HTH
    "PF01022": ["regulator", "transcription_factor", "dna_binding", "helix_turn_helix"],  # HTH_5 - ArsR family
    "PF00382": ["transcription", "dna_binding"],  # TFIIB - transcription factor TFIIB
    "PF00352": ["transcription", "dna_binding", "regulator"],  # TBP - TATA-binding protein
    "PF13404": ["regulator", "transcription_factor", "asnc_family", "dna_binding"],  # HTH_AsnC-type
    "PF01037": ["regulator", "asnc_family"],  # AsnC_trans_reg - Lrp/AsnC ligand binding

    # Purine biosynthesis - AIR synthase
    "PF00586": ["amino_acid_metabolism", "purine_metabolism", "ligase"],  # AIRS - AIR synthase N-term
    "PF02769": ["amino_acid_metabolism", "purine_metabolism"],  # AIRS_C - AIR synthase C-term

    # Polysaccharide biosynthesis
    "PF02719": ["glycosyltransferase", "transferase", "carbohydrate_active", "lps_biosynthesis"],  # Polysacc_synt_2
    "PF13524": ["glycosyltransferase", "transferase", "carbohydrate_active"],  # Glyco_trans_1_2
    "PF13632": ["glycosyltransferase", "transferase", "carbohydrate_active"],  # Glyco_trans_2_3
    "PF20706": ["glycosyltransferase", "transferase", "carbohydrate_active", "defense_system"],  # GT4-conflict

    # Glutamine amidotransferases
    "PF00117": ["glutamine_amidotransferase", "transferase", "amino_acid_metabolism"],  # GATase - class I
    "PF13522": ["glutamine_amidotransferase", "transferase"],  # GATase_6
    "PF13537": ["glutamine_amidotransferase", "transferase"],  # GATase_7

    # Aminotransferases
    "PF00266": ["aminotransferase", "transferase", "plp_binding", "amino_acid_metabolism"],  # Aminotran_5 - class V

    # Nucleotidyltransferases
    "PF00483": ["nucleotidyltransferase", "transferase"],  # NTP_transferase
    "PF12804": ["nucleotidyltransferase", "transferase", "cofactor_biosynthesis"],  # NTP_transf_3 - MobA-like

    # Secretion system components
    "PF18895": ["secretion_component", "t4ss_component", "pilus", "conjugation"],  # T4SS_pilin
    "PF00482": ["secretion_component", "t2ss_component"],  # T2SSF

    # Asparagine synthesis
    "PF00733": ["ligase", "amino_acid_biosynthesis", "asparagine_synthesis"],  # Asn_synthase

    # Amino acid kinases
    "PF00696": ["kinase", "amino_acid_metabolism", "amino_acid_biosynthesis", "atp_binding"],  # AA_kinase

    # Methyltransferases
    "PF05175": ["methyltransferase", "transferase", "sam_binding"],  # MTS - methyltransferase small domain
    "PF02384": ["methyltransferase", "dna_methylase", "restriction_modification"],  # N6_Mtase - N-6 DNA methylase
    "PF02475": ["methyltransferase", "transferase", "trna_modification"],  # TRM5-TYW2_MTfase
    "PF01135": ["methyltransferase", "transferase", "protein_modification"],  # PCMT - protein repair methyltransferase

    # Regulatory/ligand-binding domains
    "PF01842": ["regulatory", "amino_acid_metabolism"],  # ACT domain - allosteric regulator

    # NAD-binding domains (various)
    "PF07993": ["nad_binding", "cofactor_binding", "oxidoreductase"],  # NAD_binding_4
    "PF03446": ["nad_binding", "cofactor_binding", "pentose_phosphate"],  # NAD_binding_2 - 6-phosphogluconate DH

    # Phosphoribosyltransferases - nucleotide metabolism
    "PF00156": ["transferase", "nucleotide_metabolism", "purine_metabolism"],  # Pribosyltran

    # FAD-dependent oxidoreductases
    "PF12831": ["oxidoreductase", "fad_binding"],  # FAD_oxidored

    # Aminotransferase-like enzymes
    "PF01041": ["aminotransferase", "transferase", "plp_binding", "secondary_metabolism"],  # DegT_DnrJ_EryC1

    # Hexapeptide repeat transferases
    "PF00132": ["transferase", "acetyltransferase", "repeat_domain"],  # Hexapep
    "PF14602": ["transferase", "acetyltransferase", "repeat_domain"],  # Hexapep_2

    # Acetyltransferases (GNAT family)
    "PF00583": ["acetyltransferase", "transferase"],  # Acetyltransf_1 - GNAT family
    "PF13508": ["acetyltransferase", "transferase"],  # Acetyltransf_7

    # Histidine biosynthesis
    "PF00977": ["amino_acid_biosynthesis", "histidine_biosynthesis"],  # His_biosynth

    # Nitroreductases
    "PF00881": ["oxidoreductase", "reductase", "fmn_binding"],  # Nitroreductase

    # Iron transport
    "PF02421": ["transporter", "metal_transporter", "iron_binding", "gtp_binding"],  # FeoB_N - ferrous iron transport

    # GHMP kinases
    "PF00288": ["kinase", "transferase", "atp_binding"],  # GHMP_kinases_N
    "PF08544": ["kinase", "transferase"],  # GHMP_kinases_C

    # Aconitase - TCA cycle
    "PF00330": ["tca_cycle", "central_metabolism", "iron_sulfur", "lyase"],  # Aconitase

    # Arsenic resistance ATPase
    "PF02374": ["atpase", "atp_binding", "heavy_metal_resistance", "arsenic_resistance"],  # ArsA_ATPase

    # Additional ferredoxin/iron-sulfur domains
    "PF13183": ["iron_sulfur", "ferredoxin", "metal_binding", "electron_transport"],  # Fer4_8
    "PF14697": ["iron_sulfur", "ferredoxin", "metal_binding", "electron_transport"],  # Fer4_21

    # TPR repeats (additional)
    "PF13174": ["repeat_domain", "tpr_repeat", "protein_binding"],  # TPR_6

    # NAD synthase
    "PF02540": ["ligase", "nad_biosynthesis", "cofactor_biosynthesis"],  # NAD_synthase

    # SIS domain - sugar isomerase
    "PF01380": ["isomerase", "carbohydrate_active", "sugar_metabolism"],  # SIS domain

    # Protein mannosyltransferase
    "PF13231": ["glycosyltransferase", "transferase", "carbohydrate_active"],  # PMT_2

    # Helicase HerA
    "PF01935": ["helicase", "atp_binding", "dna_repair"],  # DUF87 - HerA central domain

    # Acetohydroxy acid isomeroreductase
    "PF07991": ["isomerase", "reductase", "nadp_binding", "branched_chain_aa"],  # KARI_N

    # Magnesium chelatase
    "PF01078": ["chelatase", "atp_binding", "magnesium_binding", "heme_biosynthesis"],  # Mg_chelatase

    # F420 oxidoreductase - archaeal coenzyme (NOT methanogen-specific!)
    # F420 is found in methanogens, Actinobacteria, and many other archaea
    "PF03807": ["oxidoreductase", "cofactor_binding", "f420_dependent", "archaeal_one_carbon"],  # F420_oxidored

    # HMGL-like domain - lyases
    "PF00682": ["lyase", "central_metabolism"],  # HMGL-like

    # Acylphosphatase
    "PF00708": ["hydrolase", "phosphatase", "energy_metabolism"],  # Acylphosphatase

    # Carbohydrate kinases
    "PF00294": ["kinase", "carbohydrate_active", "central_metabolism"],  # PfkB family

    # Homocitrate synthase
    "PF22617": ["transferase", "amino_acid_biosynthesis"],  # HCS_D2

    # ParA - chromosome partitioning
    "PF10609": ["cell_division", "chromosome_partitioning", "atp_binding"],  # ParA

    # OMPdecase - pyrimidine biosynthesis
    "PF00215": ["lyase", "pyrimidine_metabolism", "nucleotide_metabolism"],  # OMPdecase

    # Cytidylyltransferase
    "PF01467": ["transferase", "nucleotidyltransferase", "lipid_metabolism"],  # CTP_transf_like

    # Cysteine-rich domain
    "PF02754": ["metal_binding", "zinc_binding"],  # CCG - cysteine-rich

    # Peptidase C26
    "PF07722": ["protease", "hydrolase", "cysteine_protease"],  # Peptidase_C26

    # Flavoproteins
    "PF02441": ["oxidoreductase", "fad_binding", "electron_transport"],  # Flavoprotein

    # Helix-hairpin-helix DNA binding
    "PF14520": ["dna_binding", "dna_repair"],  # HHH_5

    # Thiamine biosynthesis
    "PF02568": ["thiamine_biosynthesis", "cofactor_biosynthesis", "sulfurtransferase"],  # ThiI

    # Metalloesterases
    "PF12850": ["hydrolase", "phosphatase", "metal_binding"],  # Metallophos_2

    # Carbamoyl-phosphate synthase
    "PF02786": ["ligase", "atp_binding", "amino_acid_metabolism", "pyrimidine_metabolism"],  # CPSase_L_D2

    # Band 7 / SPFH domain - membrane scaffolding
    "PF01145": ["membrane", "signaling", "stress_response"],  # Band_7

    # Schlafen domain
    "PF04326": ["dna_binding", "helicase"],  # SLFN_AlbA_2

    # RNA-binding metallo-hydrolase
    "PF07521": ["hydrolase", "rna_binding", "metal_binding", "zinc_binding"],  # RMMBL

    # Queuosine biosynthesis
    "PF06508": ["cofactor_biosynthesis", "trna_modification"],  # QueC

    # PhoU domain - phosphate regulation
    "PF01895": ["regulator", "phosphate_transporter", "metal_binding"],  # PhoU

    # SecD/SecF - protein secretion
    "PF02355": ["secretion_component", "sec_pathway", "membrane"],  # SecD_SecF_C

    # SMC proteins - chromosome maintenance
    "PF02463": ["cell_division", "chromosome_partitioning", "atp_binding", "dna_binding"],  # SMC_N

    # Cas12f1-like - CRISPR defense OR transposase
    # NOTE: Cas12f1-like_TNB is shared with IS605 TnpB transposases. TnpB is the evolutionary
    # ancestor of Cas12; both share this domain but serve different functions. Adding both
    # predicates makes the ambiguity explicit - check annotations to distinguish.
    "PF07282": ["cas_domain", "transposase", "nuclease"],  # Cas12f1-like_TNB - AMBIGUOUS: Cas12 OR TnpB

    # LVIVD repeat
    "PF08309": ["repeat_domain", "cell_surface"],  # LVIVD

    # UPF0020 - methyltransferase-like
    "PF01170": ["methyltransferase", "transferase"],  # UPF0020

    # PRC barrel - stress response
    "PF05239": ["stress_response", "membrane"],  # PRC-barrel

    # MazE antitoxin
    "PF04014": ["antitoxin_domain", "defense_system"],  # MazE_antitoxin

    # Amidohydrolase
    "PF07969": ["hydrolase", "amidase", "metal_binding"],  # Amidohydro_3

    # FtsX - cell division permease
    "PF02687": ["cell_division", "divisome", "membrane", "transporter"],  # FtsX

    # VWA domains (additional)
    "PF13519": ["cell_surface", "repeat_domain", "protein_binding"],  # VWA_2

    # PQQ domains - note: PQQ_3 may appear as "PQQ_3" not "PF13570" in some datasets
    # Already have PF13570, but some datasets use different accession formats

    # DUF domains (domains of unknown function) - assign based on any known associations
    "PF13635": ["hypothetical"],  # DUF4143 - unknown function
    "PF24346": ["hypothetical"],  # DUF7507 - unknown function

    # TruD - tRNA pseudouridine synthase
    "PF01142": ["trna_modification", "translation", "isomerase"],  # TruD - pseudouridine synthase

    # HEAT repeats (additional)
    "PF13646": ["repeat_domain", "heat_repeat", "protein_binding"],  # HEAT_2

    # ANAPC5 - cell cycle
    "PF12862": ["cell_division", "repeat_domain"],  # ANAPC5

    # -------------------------------------------------------------------------
    # ADDITIONAL HIGH-ABUNDANCE DOMAINS (second pass)
    # -------------------------------------------------------------------------
    # AAA ATPase variants
    "PF07726": ["aaa_domain", "atp_binding", "atpase"],  # AAA_3
    "PF13476": ["aaa_domain", "atp_binding", "atpase"],  # AAA_23
    "PF13191": ["aaa_domain", "atp_binding", "atpase"],  # AAA_16
    "PF13671": ["aaa_domain", "atp_binding", "atpase"],  # AAA_33

    # Ferredoxins (additional)
    "PF13534": ["iron_sulfur", "ferredoxin", "metal_binding", "electron_transport"],  # Fer4_17

    # PAC2 domain - signaling
    "PF09754": ["signaling"],  # PAC2

    # S1 domain - RNA binding
    "PF00575": ["rna_binding", "ribosomal_protein"],  # S1

    # TrkA - potassium transport
    "PF02254": ["transporter", "ion_transporter", "regulatory"],  # TrkA_N
    "PF02080": ["transporter", "ion_transporter"],  # TrkA_C

    # LeuA - isopropylmalate synthase
    "PF08502": ["transferase", "amino_acid_biosynthesis", "branched_chain_aa"],  # LeuA_dimer

    # GatB - glutamyl-tRNA amidotransferase
    "PF02934": ["trna_synthetase", "translation"],  # GatB_N
    "PF02637": ["trna_synthetase", "translation"],  # GatB_Yqey

    # FKBP - peptidyl-prolyl isomerase
    "PF00254": ["isomerase", "chaperone"],  # FKBP_C

    # SRP54 - signal recognition particle
    "PF02881": ["signal_recognition", "gtp_binding"],  # SRP54_N

    # PKD domain
    "PF18911": ["cell_surface", "repeat_domain"],  # PKD_4

    # Phosphomutase
    "PF10143": ["isomerase", "mutase", "carbohydrate_active"],  # PhosphMutase

    # Eco57I - restriction enzyme
    "PF07669": ["restriction_enzyme", "restriction_modification", "nuclease"],  # Eco57I

    # S4 domain - RNA binding
    "PF01479": ["rna_binding"],  # S4

    # GGR - glycosyl hydrolase related
    "PF22578": ["hydrolase", "carbohydrate_active"],  # GGR_cat

    # TraB - conjugation
    "PF01963": ["conjugation", "mobile_element"],  # TraB_PrgY_gumN

    # MS channel variants
    "PF21082": ["transporter", "ion_channel", "membrane"],  # MS_channel_3rd
    "PF21088": ["transporter", "ion_channel", "membrane"],  # MS_channel_1st

    # V4R - metal binding regulatory
    "PF02830": ["regulatory", "metal_binding"],  # V4R

    # TPR variants
    "PF13428": ["repeat_domain", "tpr_repeat", "protein_binding"],  # TPR_14

    # RHH DNA binding
    "PF01402": ["dna_binding", "ribbon_helix_helix", "regulator"],  # RHH_1

    # Thymidylate kinase
    "PF02223": ["kinase", "nucleotide_metabolism", "atp_binding"],  # Thymidylate_kin

    # NUDIX hydrolase
    "PF00293": ["hydrolase", "nucleotide_metabolism"],  # NUDIX

    # Complex I subunit
    "PF00346": ["respiration", "electron_transport", "iron_sulfur"],  # Complex1_49kDa
    "PF00329": ["respiration", "electron_transport"],  # Complex1_30kDa

    # Sua5/yciO/yrdC - translation
    "PF01300": ["translation", "trna_modification"],  # Sua5_yciO_yrdC

    # RibD - riboflavin biosynthesis
    "PF01872": ["cofactor_biosynthesis", "reductase"],  # RibD_C

    # EF-Tu domains
    "PF14578": ["translation_factor", "translation", "gtp_binding"],  # GTP_EFTU_D4
    "PF03143": ["translation_factor", "translation"],  # GTP_EFTU_D3

    # RtcB - RNA ligase
    "PF01139": ["ligase", "rna_processing"],  # RtcB

    # CDP-OH phosphotransferase
    "PF01066": ["transferase", "lipid_metabolism"],  # CDP-OH_P_transf

    # Lactamase B variants
    "PF16661": ["hydrolase", "metal_binding", "zinc_binding"],  # Lactamase_B_6
    "PF13483": ["hydrolase", "metal_binding"],  # Lactamase_B_3

    # Beta-CASP nuclease
    "PF10996": ["nuclease", "hydrolase"],  # Beta-Casp

    # Cpn60/TCP1 chaperonin
    "PF00118": ["chaperone", "hsp60", "stress_response"],  # Cpn60_TCP1

    # TSP3 bacterial
    "PF18884": ["cell_surface", "repeat_domain"],  # TSP3_bac

    # dCMP deaminase
    "PF00383": ["hydrolase", "nucleotide_metabolism"],  # dCMP_cyt_deam_1

    # Phosphodiesterase
    "PF01663": ["phosphodiesterase", "hydrolase", "signaling"],  # Phosphodiest

    # TFIIE alpha
    "PF02002": ["transcription", "dna_binding"],  # TFIIE_alpha

    # DeoC - aldolase
    "PF01791": ["lyase", "carbohydrate_active"],  # DeoC

    # PPDK N-terminal
    "PF01326": ["kinase", "central_metabolism"],  # PPDK_N

    # Thi4 - thiamine biosynthesis
    "PF01946": ["thiamine_biosynthesis", "cofactor_biosynthesis"],  # Thi4

    # tRNA-Thr editing domain
    "PF08915": ["trna_synthetase", "translation"],  # tRNA-Thr_ED

    # Cdc6 lid
    "PF22703": ["replication", "atp_binding"],  # Cdc6_lid

    # TSP_3 - thrombospondin
    "PF02412": ["cell_surface", "repeat_domain"],  # TSP_3

    # RNase HII
    "PF01351": ["rnase", "hydrolase", "dna_repair"],  # RNase_HII

    # Shikimate dehydrogenase
    "PF01488": ["oxidoreductase", "dehydrogenase", "aromatic_aa_metabolism"],  # Shikimate_DH

    # GDE C-terminal
    "PF06202": ["hydrolase", "carbohydrate_active"],  # GDE_C

    # EamA transporter
    "PF00892": ["transporter", "membrane"],  # EamA

    # PEP-utilizers
    "PF02896": ["transferase", "central_metabolism"],  # PEP-utilizers_C
    "PF00391": ["transferase", "central_metabolism"],  # PEP-utilizers

    # PDDEXK nuclease
    "PF12705": ["nuclease", "hydrolase"],  # PDDEXK_1

    # Formyl transferase
    "PF00551": ["transferase", "one_carbon_metabolism"],  # Formyl_trans_N

    # PAP2 phosphatase
    "PF01569": ["phosphatase", "hydrolase", "lipid_metabolism"],  # PAP2

    # Cupin 2
    "PF07883": ["hydrolase"],  # Cupin_2

    # CTP synthetase
    "PF06418": ["ligase", "nucleotide_metabolism", "atp_binding"],  # CTP_synth_N

    # PTPS
    "PF01242": ["lyase", "cofactor_biosynthesis"],  # PTPS

    # Mrr restriction enzyme
    "PF04471": ["restriction_enzyme", "restriction_modification", "nuclease"],  # Mrr_cat

    # TrmB - tRNA modification
    "PF01978": ["methyltransferase", "trna_modification"],  # TrmB

    # CDC48 N-terminal
    "PF02359": ["aaa_domain", "atp_binding", "chaperone"],  # CDC48_N
    "PF02933": ["aaa_domain", "atp_binding"],  # CDC48_2

    # tRNA intron endonuclease
    "PF01974": ["nuclease", "rna_processing"],  # tRNA_int_endo
    "PF02778": ["nuclease", "rna_processing"],  # tRNA_int_endo_N

    # GARS - glycinamide ribonucleotide synthetase
    "PF01071": ["ligase", "purine_metabolism", "atp_binding"],  # GARS_A
    "PF02844": ["ligase", "purine_metabolism"],  # GARS_N

    # MTHFR - methylenetetrahydrofolate reductase
    "PF02219": ["oxidoreductase", "one_carbon_metabolism", "fad_binding"],  # MTHFR

    # PIN domain variants
    "PF13470": ["toxin_domain", "pin_domain"],  # PIN_3
    "PF17146": ["toxin_domain", "pin_domain"],  # PIN_6
    "PF18477": ["toxin_domain", "pin_domain"],  # PIN_9
    "PF10130": ["toxin_domain", "pin_domain"],  # PIN_2

    # Polysaccharide deacetylase/synthesis
    "PF01522": ["hydrolase", "carbohydrate_active"],  # Polysacc_deac_1
    "PF13440": ["glycosyltransferase", "carbohydrate_active"],  # Polysacc_synt_3

    # SHMT - serine hydroxymethyltransferase
    "PF00464": ["transferase", "one_carbon_metabolism", "plp_binding"],  # SHMT

    # GrpE - nucleotide exchange factor
    "PF01025": ["chaperone", "heat_shock"],  # GrpE

    # PolC DNA polymerase
    "PF24846": ["dna_polymerase", "replication"],  # PolC_DP2_cat
    "PF24844": ["dna_polymerase", "replication"],  # PolC_DP2_central
    "PF03833": ["dna_polymerase", "replication"],  # PolC_DP2_N

    # TsaD - tRNA modification
    "PF00814": ["trna_modification", "translation"],  # TsaD

    # RNase H variant
    "PF13482": ["rnase", "hydrolase"],  # RNase_H_2

    # Orn/Arg decarboxylase
    "PF02784": ["lyase", "decarboxylase", "amino_acid_metabolism"],  # Orn_Arg_deC_N
    "PF00278": ["lyase", "decarboxylase", "amino_acid_metabolism"],  # Orn_DAP_Arg_deC

    # ThiC - thiamine biosynthesis
    "PF01964": ["radical_sam", "thiamine_biosynthesis"],  # ThiC_Rad_SAM

    # TehB - tellurite resistance
    "PF03848": ["methyltransferase", "heavy_metal_resistance"],  # TehB

    # MreB - cell shape
    "PF06723": ["cell_division", "atp_binding"],  # MreB_Mbl

    # Lon protease C-terminal
    "PF05362": ["protease", "atp_binding", "chaperone"],  # Lon_C

    # POR - pyruvate:ferredoxin oxidoreductase
    "PF01855": ["oxidoreductase", "iron_sulfur", "central_metabolism"],  # POR_N
    "PF01558": ["oxidoreductase", "iron_sulfur", "central_metabolism"],  # POR

    # TruB - pseudouridine synthase
    "PF01509": ["trna_modification", "isomerase"],  # TruB_N
    "PF16198": ["trna_modification", "isomerase"],  # TruB_C_2

    # Prenyltransferase
    "PF01255": ["transferase", "isoprenoid_biosynthesis"],  # Prenyltransf

    # DNA pol B exonuclease
    "PF03104": ["dna_polymerase", "nuclease", "replication"],  # DNA_pol_B_exo1

    # FecCD - iron transport
    "PF01032": ["transporter", "abc_transporter", "iron_binding"],  # FecCD

    # Glyco_hydro_57
    "PF03065": ["glycosidase", "hydrolase", "carbohydrate_active"],  # Glyco_hydro_57

    # DJ-1/PfpI - oxidative stress
    "PF01965": ["oxidative_stress", "protease"],  # DJ-1_PfpI

    # Histone deacetylase
    "PF00850": ["hydrolase", "regulator"],  # Hist_deacetyl

    # GATase variants
    "PF13230": ["glutamine_amidotransferase", "transferase"],  # GATase_4
    "PF07685": ["glutamine_amidotransferase", "transferase"],  # GATase_3
    "PF13507": ["glutamine_amidotransferase", "transferase"],  # GATase_5

    # Thiolase
    "PF00108": ["transferase", "fatty_acid_degradation", "lipid_metabolism"],  # Thiolase_N
    "PF02803": ["transferase", "fatty_acid_degradation"],  # Thiolase_C
    "PF22691": ["transferase", "fatty_acid_degradation"],  # Thiolase_C_1

    # TrbC - conjugation
    "PF04956": ["conjugation", "pilus"],  # TrbC

    # EF-1 alpha
    "PF22594": ["translation_factor", "translation", "gtp_binding"],  # GTP-eEF1A_C
    "PF00736": ["translation_factor", "translation"],  # EF1_GNE

    # SepF - cell division
    "PF04472": ["cell_division", "divisome"],  # SepF

    # DHQS - shikimate pathway
    "PF01959": ["lyase", "aromatic_aa_metabolism"],  # DHQS

    # SBP bacterial
    "PF13531": ["transporter", "abc_transporter", "periplasmic"],  # SBP_bac_11

    # SepSecS - selenocysteine synthesis
    "PF05889": ["trna_synthetase", "selenium_metabolism"],  # SepSecS

    # CorA magnesium transporter
    "PF01544": ["transporter", "metal_transporter", "magnesium_binding"],  # CorA

    # DNA pol3 delta
    "PF13177": ["dna_polymerase", "replication"],  # DNA_pol3_delta2

    # TatD DNase
    "PF01026": ["dnase", "hydrolase", "dna_repair"],  # TatD_DNase

    # Nre - nitrogen regulation
    "PF04894": ["regulator", "nitrogen_metabolism"],  # Nre_N
    "PF04895": ["regulator", "nitrogen_metabolism"],  # Nre_C

    # DHO dehydrogenase
    "PF01180": ["oxidoreductase", "pyrimidine_metabolism"],  # DHO_dh

    # DapB - diaminopimelate biosynthesis
    "PF01113": ["oxidoreductase", "amino_acid_biosynthesis"],  # DapB_N
    "PF05173": ["oxidoreductase", "amino_acid_biosynthesis"],  # DapB_C

    # AP endonuclease
    "PF01261": ["nuclease", "dna_repair", "base_excision_repair"],  # AP_endonuc_2

    # MCM lid
    "PF17855": ["helicase", "replication"],  # MCM_lid

    # SWI2/SNF2
    "PF18766": ["helicase", "atp_binding", "regulator"],  # SWI2_SNF2

    # Glycogen synthase
    "PF05693": ["glycosyltransferase", "carbohydrate_active"],  # Glycogen_syn

    # EFG domains
    "PF14492": ["translation_factor", "translation"],  # EFG_III

    # HMG-CoA reductase
    "PF00368": ["oxidoreductase", "isoprenoid_biosynthesis"],  # HMG-CoA_red

    # DHDPS - dihydrodipicolinate synthase
    "PF00701": ["lyase", "amino_acid_biosynthesis"],  # DHDPS

    # TGT - tRNA-guanine transglycosylase
    "PF01702": ["transferase", "trna_modification"],  # TGT
    "PF14810": ["transferase", "trna_modification"],  # TGT_C2

    # EIF-2 alpha
    "PF07541": ["translation_factor", "translation"],  # EIF_2_alpha

    # SBDS domains
    "PF01172": ["ribosomal_protein", "translation"],  # SBDS_N
    "PF09377": ["ribosomal_protein", "translation"],  # SBDS_domain_II
    "PF20268": ["ribosomal_protein", "translation"],  # SBDS_C

    # Amidase
    "PF01425": ["hydrolase", "amidase"],  # Amidase

    # SMC/ScpB
    "PF04079": ["cell_division", "chromosome_partitioning"],  # SMC_ScpB

    # HxlR regulator
    "PF01638": ["regulator", "one_carbon_metabolism"],  # HxlR

    # IGPS - imidazole glycerol phosphate synthase
    "PF00218": ["transferase", "histidine_biosynthesis"],  # IGPS

    # HisG - histidine biosynthesis
    "PF01634": ["transferase", "histidine_biosynthesis"],  # HisG

    # Fibrillarin - rRNA methyltransferase
    "PF01269": ["methyltransferase", "rrna_modification"],  # Fibrillarin

    # NDK - nucleoside diphosphate kinase
    "PF00334": ["kinase", "nucleotide_metabolism"],  # NDK

    # dUTPase
    "PF00692": ["hydrolase", "nucleotide_metabolism"],  # dUTPase

    # NMD3 - ribosome biogenesis
    "PF04981": ["ribosomal_protein", "translation"],  # NMD3

    # PTH2 - peptidyl-tRNA hydrolase
    "PF01981": ["hydrolase", "translation"],  # PTH2

    # Histidinol dehydrogenase
    "PF00815": ["oxidoreductase", "histidine_biosynthesis"],  # Histidinol_dh

    # GMP synthetase
    "PF00958": ["ligase", "purine_metabolism"],  # GMP_synt_C

    # PCNA
    "PF00705": ["replication", "dna_binding"],  # PCNA_N
    "PF02747": ["replication", "dna_binding"],  # PCNA_C

    # TauE - sulfite exporter
    "PF01925": ["transporter", "sulfur_metabolism"],  # TauE

    # KARI C-terminal
    "PF01450": ["isomerase", "reductase", "branched_chain_aa"],  # KARI_C

    # Methyltransferase RsmB/F
    "PF01189": ["methyltransferase", "rrna_modification"],  # Methyltr_RsmB-F
    "PF17125": ["methyltransferase", "rrna_modification"],  # Methyltr_RsmF_N

    # Topo-VIb
    "PF09239": ["topoisomerase", "replication"],  # Topo-VIb_trans

    # Spt5-NGN
    "PF03439": ["transcription"],  # Spt5-NGN

    # QRPTase - quinolinate phosphoribosyltransferase
    "PF01729": ["transferase", "nad_biosynthesis"],  # QRPTase_C
    "PF02749": ["transferase", "nad_biosynthesis"],  # QRPTase_N

    # Prefoldin
    "PF01920": ["chaperone"],  # Prefoldin_2
    "PF02996": ["chaperone"],  # Prefoldin

    # NikR C-terminal
    "PF08753": ["regulator", "nickel_binding"],  # NikR_C

    # Fe-ADH
    "PF13685": ["oxidoreductase", "dehydrogenase", "iron_binding"],  # Fe-ADH_2

    # ArgJ - acetylornithine aminotransferase
    "PF01960": ["aminotransferase", "amino_acid_biosynthesis"],  # ArgJ

    # eIF-1a
    "PF01176": ["translation_factor", "translation"],  # eIF-1a

    # eIF-6
    "PF01912": ["translation_factor", "translation"],  # eIF-6

    # Lycopene cyclase
    "PF05834": ["lyase", "isoprenoid_biosynthesis"],  # Lycopene_cycl

    # Asparaginase
    "PF00710": ["hydrolase", "amino_acid_metabolism"],  # Asparaginase
    "PF17763": ["hydrolase", "amino_acid_metabolism"],  # Asparaginase_C

    # Homoserine dehydrogenase
    "PF00742": ["oxidoreductase", "amino_acid_biosynthesis"],  # Homoserine_dh

    # PNP/UDP
    "PF01048": ["transferase", "nucleotide_metabolism"],  # PNP_UDP_1

    # Topoisomerase VI
    "PF20768": ["topoisomerase", "replication"],  # Topo_VI_alpha

    # Rib-5-P isomerase
    "PF06026": ["isomerase", "pentose_phosphate"],  # Rib_5-P_isom_A

    # DnaJ C-terminal
    "PF01556": ["chaperone", "heat_shock"],  # DnaJ_C

    # SUI1
    "PF01253": ["translation_factor", "translation"],  # SUI1

    # UbiD - ubiquinone biosynthesis
    "PF01977": ["lyase", "cofactor_biosynthesis"],  # UbiD
    "PF20696": ["lyase", "cofactor_biosynthesis"],  # UbiD_C
    "PF20695": ["lyase", "cofactor_biosynthesis"],  # UbiD_N

    # eIF-5/eIF-2B
    "PF01873": ["translation_factor", "translation"],  # eIF-5_eIF-2B

    # SDH C-terminal
    "PF18317": ["oxidoreductase", "tca_cycle"],  # SDH_C

    # DHquinase
    "PF01487": ["hydrolase", "aromatic_aa_metabolism"],  # DHquinase_I

    # MannoseP isomerase
    "PF01050": ["isomerase", "carbohydrate_active"],  # MannoseP_isomer

    # TctA - tricarboxylate transporter
    "PF01970": ["transporter", "membrane"],  # TctA

    # TIM - triosephosphate isomerase
    "PF00121": ["isomerase", "glycolysis"],  # TIM

    # zf-HYPF
    "PF07503": ["zinc_binding", "hydrogenase_maturation"],  # zf-HYPF

    # Nop - nucleolar protein
    "PF01798": ["rna_binding", "rrna_modification"],  # Nop

    # CM_2 - chorismate mutase
    "PF01817": ["isomerase", "aromatic_aa_metabolism"],  # CM_2

    # Trp_syntA - tryptophan synthase
    "PF00290": ["transferase", "aromatic_aa_metabolism"],  # Trp_syntA

    # FtsJ methyltransferase
    "PF01728": ["methyltransferase", "rrna_modification"],  # FtsJ

    # STT3 - oligosaccharyltransferase
    "PF02516": ["glycosyltransferase", "membrane"],  # STT3

    # TP6A N-terminal
    "PF04406": ["topoisomerase", "replication"],  # TP6A_N

    # MipZ - cell division
    "PF09140": ["cell_division", "atp_binding"],  # MipZ

    # Transketolase C-terminal
    "PF02780": ["transferase", "pentose_phosphate"],  # Transketolase_C

    # PEPCK GTP
    "PF00821": ["lyase", "gluconeogenesis"],  # PEPCK_GTP
    "PF17297": ["lyase", "gluconeogenesis"],  # PEPCK_N

    # SMC hinge
    "PF06470": ["cell_division", "chromosome_partitioning"],  # SMC_hinge

    # Shikimate_dh_N
    "PF08501": ["oxidoreductase", "aromatic_aa_metabolism"],  # Shikimate_dh_N

    # -------------------------------------------------------------------------
    # BATCH 1: HIGH-ABUNDANCE UNMAPPED DOMAINS (Altiarchaeota)
    # -------------------------------------------------------------------------
    # PQQ_3 - PQQ enzyme repeat (no PF prefix in some datasets)
    "PQQ_3": ["oxidoreductase", "pqq_binding"],  # PQQ enzyme repeat

    # ATP-binding proteins
    "PF21690": ["atp_binding", "hypothetical"],  # MJ1010-like_2nd - uncharacterized ATP-binding

    # Uncharacterized protein families (UPFs)
    "PF03683": ["hypothetical"],  # UPF0175
    "PF01894": ["hypothetical"],  # UPF0047
    "PF00919": ["hypothetical"],  # UPF0004
    "PF03602": ["hypothetical"],  # Cons_hypoth95

    # RNases and nucleases
    "PF01934": ["rnase", "hydrolase"],  # HepT-like - Ribonuclease HepT-like
    "PF06250": ["nuclease", "hydrolase"],  # YhcG_C - PDDEXK nuclease domain

    # DUF domains
    "PF17761": ["hypothetical"],  # DUF1016_N
    "PF01865": ["hypothetical"],  # PhoU_div (DUF47)

    # Signaling and regulatory domains
    "PF09335": ["signaling"],  # VTT_dom - VTT domain
    "PF01871": ["regulatory"],  # AMMECR1
    "PF22629": ["regulatory", "amino_acid_metabolism"],  # ACT_AHAS_ss - AHAS small subunit ACT domain

    # Membrane proteins
    "PF04307": ["membrane"],  # YdjM inner membrane protein
    "PF01956": ["membrane"],  # EMC3_TMCO1 - integral membrane protein

    # Phosphatases and phosphoesterases
    "PF02272": ["phosphatase", "hydrolase"],  # DHHA1 domain
    "PF02811": ["phosphatase", "hydrolase", "metal_binding"],  # PHP domain - polymerase/histidinol phosphatase

    # Deaminases
    "PF14437": ["hydrolase", "deaminase"],  # MafB19-deam - MafB19-like deaminase
    "PF22769": ["hydrolase", "nucleotide_metabolism"],  # DCD - dCTP deaminase-like

    # AAA+ ATPase domains
    "PF21960": ["aaa_domain", "atp_binding", "atpase"],  # RCF1-5-like_lid - AAA+ ATPase lid domain

    # Translation-related
    "PF01916": ["translation", "amino_acid_metabolism"],  # DS - Deoxyhypusine synthase (eIF5A modification)
    "PF07973": ["trna_synthetase", "translation"],  # tRNA_SAD - tRNA synthetase additional domain
    "PF21485": ["translation_factor", "translation"],  # IF5A-like_N - Translation initiation factor 5A-like

    # Kinases
    "PF08543": ["kinase", "thiamine_biosynthesis", "cofactor_biosynthesis"],  # Phos_pyr_kin - Phosphomethylpyrimidine kinase
    "PF13189": ["kinase", "nucleotide_metabolism"],  # Cytidylate_kin2 - Cytidylate kinase-like

    # Transferases
    "PF01174": ["glutamine_amidotransferase", "transferase"],  # SNO - glutamine amidotransferase
    "PF02005": ["methyltransferase", "trna_modification"],  # TRM - N2,N2-dimethylguanosine tRNA methyltransferase
    "PF02366": ["glycosyltransferase", "transferase"],  # PMT - Dolichyl-phosphate-mannose-protein mannosyltransferase
    "PF13793": ["transferase", "nucleotide_metabolism"],  # Pribosyltran_N - ribose phosphate pyrophosphokinase

    # Secretion pathway
    "PF07549": ["secretion_component", "sec_pathway"],  # Sec_GG - SecD/SecF GG Motif

    # Hydrogenase maturation
    "PF22521": ["hydrogenase_maturation", "transferase"],  # HypF_C_2 - Carbamoyltransferase Kae1-like

    # Cell division / sporulation
    "PF01944": ["cell_division"],  # SpoIIM - Stage II sporulation protein M

    # Mobile element related
    "PF02498": ["mobile_element"],  # Bro-N - BRO family (associated with baculovirus)

    # DNA replication
    "PF09845": ["replication"],  # OapC - Origin-associated protein OapC
    "PF17207": ["helicase", "replication", "rna_binding"],  # MCM_OB - MCM OB domain
    "PF21473": ["dna_binding", "replication"],  # Ssb-like_OB - Single-stranded DNA binding protein OB fold

    # Structural/repeat domains
    "PF01875": ["hypothetical"],  # Memo - Memo-like protein (function uncertain)
    "PF00932": ["cell_division"],  # LTD - Lamin Tail Domain
    "PF23477": ["zinc_binding", "metal_binding"],  # zf_Tbcl_2 - treble clef zinc finger
    "PF06777": ["hypothetical"],  # HBB - Helical and beta-bridge domain
    "PF05833": ["rna_binding", "translation"],  # NFACT_N - NFACT N-terminal
    "PF03259": ["signaling"],  # Robl_LC7 - Roadblock/LC7 domain (dynein regulation)

    # Oxidoreductases
    "PF14512": ["oxidoreductase", "reductase"],  # TM1586_NiRdase - putative nitroreductase
    "PF04127": ["oxidoreductase", "fad_binding", "cofactor_biosynthesis"],  # DFP - DNA/pantothenate metabolism flavoprotein

    # Amino acid metabolism
    "PF20979": ["ligase", "amino_acid_biosynthesis"],  # Arginosuc_syn_C - Arginosuccinate synthase C-terminal
    "PF08541": ["fatty_acid_synthesis", "lipid_metabolism"],  # ACP_syn_III_C - 3-Oxoacyl-ACP synthase III C terminal

    # Signaling - HD domain associated
    "PF19276": ["signaling", "phosphodiesterase"],  # HD_assoc_2 - HD associated region
    "PF18306": ["hypothetical"],  # LDcluster4 - SLOG cluster4 family

    # Other enzymes
    "PF03853": ["hydrolase", "carbohydrate_active"],  # YjeF_N - YjeF-related protein (sugar phosphatase)
    "PF03109": ["kinase", "atp_binding"],  # ABC1 - ABC1 atypical kinase-like domain

    # -------------------------------------------------------------------------
    # BATCH 2: MORE HIGH-ABUNDANCE UNMAPPED DOMAINS
    # -------------------------------------------------------------------------
    # Transcription and RNA processing
    "PF01849": ["dna_binding", "regulator"],  # NAC domain - plant-specific transcription factors
    "PF17214": ["rna_binding", "transcription_termination"],  # KH_TffA - transcription termination factor
    "PF06093": ["zinc_binding", "transcription"],  # Spt4 - Spt4/RpoE2 zinc finger

    # DNA replication and repair
    "PF08542": ["replication", "atp_binding"],  # Rep_fac_C - Replication factor C C-terminal
    "PF02732": ["nuclease", "dna_repair"],  # ERCC4 domain - XPF endonuclease
    "PF00752": ["nuclease", "dna_repair"],  # XPG_N - XPG N-terminal domain (NER)
    "PF00867": ["nuclease", "dna_repair"],  # XPG_I - XPG I-region (NER)
    "PF09079": ["replication", "dna_binding", "winged_helix"],  # Cdc6_C - CDC6 C terminal
    "PF08646": ["replication", "dna_binding"],  # Rep_fac-A_C - Replication factor-A C terminal

    # Phosphoesterases
    "PF21763": ["phosphatase", "hydrolase"],  # DHH_CID domain

    # Carbohydrate metabolism
    "PF06564": ["glycosyltransferase", "carbohydrate_active"],  # CBP_BcsQ - cellulose biosynthesis

    # Pyridoxal biosynthesis
    "PF01680": ["cofactor_biosynthesis", "pyridoxal_biosynthesis"],  # SOR_SNZ family

    # Translation machinery
    "PF03484": ["trna_synthetase", "translation"],  # B5 - tRNA synthetase B5 domain
    "PF19303": ["trna_synthetase", "translation", "rna_binding"],  # Anticodon_3 - anticodon binding
    "PF17777": ["ribosomal_protein", "translation"],  # RL10P_insert - 60S ribosomal protein L10P
    "PF01951": ["rna_processing", "ligase"],  # Archease - RNA ligase cofactor
    "PF09173": ["translation_factor", "translation", "gtp_binding"],  # eIF2_C - initiation factor eIF2 gamma
    "PF08704": ["methyltransferase", "trna_modification"],  # GCD14 - tRNA methyltransferase
    "PF21800": ["rna_binding", "translation"],  # KH_KRR1_2nd - KRR1 ribosome biogenesis
    "PF18195": ["trna_synthetase", "translation"],  # GatD_N - Glu-tRNA amidotransferase

    # Regulatory domains
    "PF05368": ["regulatory", "oxidoreductase"],  # NmrA-like - NAD(P)-binding regulatory

    # Kinases
    "PF07714": ["kinase", "serine_threonine_kinase", "tyrosine_kinase"],  # PK_Tyr_Ser-Thr

    # tRNA modification
    "PF22641": ["trna_modification", "translation"],  # TiaS_TCKD domain
    "PF01994": ["methyltransferase", "trna_modification"],  # Trm56 - tRNA 2'-O-methyltransferase
    "PF08608": ["trna_modification", "translation"],  # Wyosine_form - wyosine base formation
    "PF21133": ["nucleotidyltransferase", "translation"],  # CAA_C - CCA-adding enzyme C-terminal

    # Nucleotide metabolism
    "PF07831": ["transferase", "nucleotide_metabolism"],  # PYNP_C - pyrimidine nucleoside phosphorylase
    "PF10397": ["lyase", "purine_metabolism"],  # ADSL_C - Adenylosuccinate lyase C-terminus

    # Chaperones and protein folding
    "PF22199": ["isomerase", "chaperone"],  # FKBP26_IF - FK506-binding protein
    "PF00684": ["chaperone", "zinc_binding"],  # DnaJ_CXXCXGXG - DnaJ central domain

    # GTPases
    "PF08438": ["gtpase", "gtp_binding"],  # YGR210-like_G4 - Obg-like GTPase

    # RNA processing
    "PF01868": ["rnase", "rna_processing"],  # RNase_P-MRP_p29 - RNase P subunit
    "PF01900": ["rnase", "rna_processing"],  # RNase_P_Rpp14 - Rpp14/Pop5 family
    "PF05670": ["rna_binding", "translation"],  # NFACT-R_1 - NFACT RNA binding domain
    "PF10469": ["ligase", "rna_processing"],  # AKAP7_NLS - 2'5' RNA ligase-like

    # Purine biosynthesis
    "PF18072": ["transferase", "purine_metabolism"],  # FGAR-AT_linker - FGAM synthase linker

    # OB-fold domains
    "PF01796": ["dna_binding", "rna_binding"],  # OB_ChsH2_C - OB-fold domain

    # Oxidoreductases
    "PF17147": ["oxidoreductase", "iron_sulfur", "central_metabolism"],  # PFOR_II - pyruvate:ferredoxin oxidoreductase

    # Proteases and peptidases
    "PF00026": ["protease", "hydrolase"],  # Asp - Eukaryotic aspartyl protease
    "PF20436": ["aaa_domain", "protease", "atp_binding"],  # LonB_AAA-LID - archaeal LonB
    "PF02517": ["protease", "membrane"],  # Rce1-like - CAAX prenyl endopeptidase

    # Ligases
    "PF21948": ["ligase", "cofactor_binding"],  # LplA-B_cat - Lipoyl protein ligase

    # Zinc finger/metal binding
    "PF12172": ["zinc_binding", "metal_binding"],  # zf-ChsH2 - rubredoxin-like zinc ribbon

    # DNA binding
    "PF00633": ["dna_binding", "dna_repair"],  # HHH - Helix-hairpin-helix motif

    # Aldolases
    "PF03737": ["lyase", "aldolase"],  # RraA-like - Aldolase/RraA

    # NTPases
    "PF01725": ["hydrolase", "nucleotide_metabolism"],  # Ham1p_like - NTP pyrophosphatase

    # Membrane trafficking
    "PF04893": ["membrane", "vesicle_trafficking"],  # Yip1 domain

    # Amino acid metabolism
    "PF00800": ["lyase", "aromatic_aa_metabolism"],  # PDT - Prephenate dehydratase
    "PF01958": ["oxidoreductase", "amino_acid_metabolism"],  # Asp_DH_C - Aspartate dehydrogenase

    # LPS biosynthesis
    "PF06293": ["kinase", "lps_biosynthesis"],  # Kdo - Lipopolysaccharide kinase

    # Phosphatases
    "PF00459": ["phosphatase", "hydrolase"],  # Inositol_P - Inositol monophosphatase family

    # -------------------------------------------------------------------------
    # BATCH 3: MORE UNMAPPED DOMAINS
    # -------------------------------------------------------------------------
    # Cofactor biosynthesis
    "PF01874": ["transferase", "cofactor_biosynthesis"],  # CitG - ATP:dephospho-CoA triphosphoribosyl transferase
    "PF02445": ["lyase", "nad_biosynthesis", "cofactor_biosynthesis"],  # NadA - Quinolinate synthetase A
    "PF10120": ["transferase", "thiamine_biosynthesis"],  # ThiN - Thiamine-phosphate synthase
    "PF01973": ["kinase", "cofactor_biosynthesis"],  # MptE-like - 6-hydroxymethylpterin diphosphokinase

    # tRNA modification and translation
    "TYW2_N_arc": ["trna_modification", "translation"],  # TYW2_N_arc - tRNA wybutosine synthesis
    "PF02676": ["methyltransferase", "trna_modification"],  # TYW3 - tRNA wybutosine methyltransferase
    "PF05746": ["trna_synthetase", "translation", "rna_binding"],  # DALR_1 - anticodon binding
    "PF04414": ["hydrolase", "translation"],  # tRNA_deacylase - D-aminoacyl-tRNA deacylase
    "PF01207": ["oxidoreductase", "trna_modification"],  # Dus - Dihydrouridine synthase

    # Amino acid metabolism
    "PF01837": ["transferase", "sulfur_metabolism", "amino_acid_biosynthesis"],  # HcyBio - Homocysteine biosynthesis
    "PF02153": ["oxidoreductase", "aromatic_aa_metabolism"],  # PDH_N - Prephenate dehydrogenase
    "PF20463": ["oxidoreductase", "aromatic_aa_metabolism"],  # PDH_C - Prephenate dehydrogenase
    "PF04715": ["transferase", "aromatic_aa_metabolism"],  # Anth_synt_I_N - Anthranilate synthase
    "PF01053": ["transferase", "plp_binding", "sulfur_metabolism"],  # Cys_Met_Meta_PP - Cys/Met metabolism
    "PF01862": ["lyase", "decarboxylase", "amino_acid_metabolism"],  # PvlArgDC - pyruvoyl-dependent arginine decarboxylase
    "PF21570": ["oxidoreductase", "amino_acid_metabolism"],  # ArgZ-like_C_2nd - Arginine dihydrolase

    # DUF domains with some functional hints
    "PF01982": ["hypothetical"],  # CTP-dep_RFKase (DUF120)
    "PF08445": ["hypothetical"],  # FR47-like protein

    # RNA processing
    "PF01876": ["rnase", "rna_processing"],  # RNase_P_p30 - RNase P subunit
    "PF14382": ["nuclease", "rna_processing"],  # ECR1_N - Exosome complex exonuclease
    "PF01137": ["ligase", "rna_processing"],  # RTC - RNA 3'-terminal phosphate cyclase
    "PF05189": ["ligase", "rna_processing"],  # RTC_insert - RTC insert domain
    "PF22505": ["rnase", "rna_processing"],  # RNase_J_b_CASP - Ribonuclease J

    # Phospholipase and signaling
    "PF13091": ["hydrolase", "lipase", "signaling"],  # PLDc_2 - Phospholipase D

    # Stress response
    "PF00582": ["stress_response"],  # Usp - Universal stress protein family

    # DNA binding and transcription
    "PF00808": ["dna_binding", "regulator"],  # CBFD_NFYB_HMF - Histone-like transcription factor/archaeal histone
    "PF03551": ["regulator", "transcription_factor", "dna_binding"],  # PadR - PadR-like family
    "PF01325": ["regulator", "iron_binding", "dna_binding"],  # Fe_dep_repress - Iron dependent repressor N-terminal
    "PF02742": ["regulator", "metal_binding"],  # Fe_dep_repr_C - Iron dependent repressor C-terminal
    "PF03444": ["regulator", "dna_binding", "winged_helix"],  # HrcA_DNA-bdg - HrcA DNA-binding

    # Cell surface and repeat domains
    "PF01839": ["repeat_domain", "cell_surface"],  # FG-GAP repeat

    # Epigenetics/chromatin
    "PF23613": ["acetyltransferase", "trna_modification"],  # ELP3_N - Elongator complex acetyltransferase
    "PF08433": ["regulator", "transcription"],  # KTI12 - chromatin associated protein

    # DNA repair and replication
    "PF13298": ["ligase", "dna_repair"],  # LigD_N - DNA Ligase D phosphoesterase
    "PF09376": ["nuclease", "dna_repair"],  # NurA - NurA nuclease
    "PF01870": ["nuclease", "dna_repair", "recombinational_repair"],  # Hjc - Holliday junction resolvase
    "PF22175": ["dna_repair", "base_excision_repair"],  # Ogg-HhH - 8-oxoguanine DNA glycosylase
    "PF09846": ["replication"],  # OapB - Origin-associated protein

    # GTPases and signaling
    "PF16897": ["gtpase", "translation"],  # MMR_HSR1_Xtn - ribosome-binding GTPase
    "PF07088": ["atp_binding", "membrane"],  # GvpD_P-loop - gas vesicle protein

    # Tetrahydromethanopterin pathway - found in methanogens AND non-methanogens (e.g., Altiarchaeota)
    "PF01913": ["transferase", "archaeal_one_carbon", "one_carbon_metabolism"],  # FTR - Formylmethanofuran-tetrahydromethanopterin formyltransferase
    "PF02741": ["transferase", "archaeal_one_carbon"],  # FTR_C - FTR proximal lobe
    "PF09176": ["oxidoreductase", "archaeal_one_carbon", "one_carbon_metabolism"],  # Mpt_N - Methylene-tetrahydromethanopterin dehydrogenase

    # Oxidoreductases
    "PF07584": ["regulatory"],  # BatA - Aerotolerance regulator
    "PF21349": ["oxidative_stress", "iron_binding"],  # RUBY_RBDX - Rubrerythrin rubredoxin-like domain
    "PF02525": ["oxidoreductase", "fad_binding"],  # Flavodoxin_2 - Flavodoxin-like fold

    # tRNA modification additional
    "PF08489": ["trna_modification", "translation"],  # TiaS_FLD domain
    "PF02926": ["rna_binding", "trna_modification"],  # THUMP domain

    # Kinases
    "PF02110": ["kinase", "thiamine_biosynthesis"],  # HK - Hydroxyethylthiazole kinase

    # Ribosome biogenesis
    "PF08068": ["ribosomal_protein", "translation"],  # DKCLD - pseudouridine synthase
    "PF04135": ["rna_binding", "rrna_modification"],  # Nop10p - Nucleolar RNA-binding protein

    # Proteases and peptidases
    "PF14464": ["protease", "hydrolase"],  # Prok-JAB - JAB domain (deubiquitinating)

    # tRNA synthetase domains
    "PF03483": ["trna_synthetase", "translation"],  # B3_4 domain

    # Lipid metabolism
    "PF02353": ["methyltransferase", "lipid_metabolism"],  # CMAS - Mycolic acid cyclopropane synthetase

    # Amino acid biosynthesis
    "PF10369": ["regulatory", "branched_chain_aa"],  # ALS_ss_C - Acetolactate synthase small subunit
    "PF02006": ["ligase", "cofactor_biosynthesis"],  # PPS_PS - Phosphopantothenate synthetase
    "PF00475": ["lyase", "histidine_biosynthesis"],  # IGPD - Imidazoleglycerol-phosphate dehydratase

    # Iron-sulfur clusters
    "PF13459": ["iron_sulfur", "ferredoxin", "metal_binding", "electron_transport"],  # Fer4_15 - 4Fe-4S single cluster
    "PF12801": ["iron_sulfur", "ferredoxin", "metal_binding", "electron_transport"],  # Fer4_5 - 4Fe-4S binding
    "PF04060": ["iron_sulfur", "metal_binding"],  # FeS - putative Fe-S cluster

    # Isoprenoid biosynthesis
    "PF22700": ["lyase", "isoprenoid_biosynthesis"],  # MVD-like_N - Diphosphomevalonate decarboxylase

    # Transporters
    "PF03773": ["transporter", "membrane"],  # ArsP_1 - Predicted permease

    # Histidine biosynthesis
    "PF01502": ["hydrolase", "histidine_biosynthesis"],  # PRA-CH - Phosphoribosyl-AMP cyclohydrolase

    # Folate biosynthesis
    "PF02649": ["hydrolase", "cofactor_biosynthesis"],  # GCHY-1 - GTP cyclohydrolase (folate)

    # Chromosome partitioning
    "ParBc": ["cell_division", "chromosome_partitioning", "dna_binding"],  # ParBc - ParB/Spo0J DNA-binding

    # Carbohydrate metabolism
    "PF12439": ["hydrolase", "carbohydrate_active"],  # GDE_N - Glycogen debranching enzyme
    "PF00834": ["isomerase", "pentose_phosphate"],  # Ribul_P_3_epim - Ribulose-phosphate 3 epimerase
    "PF04412": ["tca_cycle", "iron_sulfur", "lyase"],  # AcnX - Aconitase X

    # One-carbon/folate metabolism
    "PF02289": ["hydrolase", "one_carbon_metabolism"],  # MCH - Cyclohydrolase

    # Pyrimidine metabolism
    "PF01948": ["regulatory", "pyrimidine_metabolism"],  # PyrI - Aspartate carbamoyltransferase regulatory
    "PF02748": ["regulatory", "pyrimidine_metabolism", "metal_binding"],  # PyrI_C - ATCase regulatory metal-binding
    "PF02511": ["transferase", "nucleotide_metabolism"],  # Thy1 - Thymidylate synthase complementing protein

    # Dehydratases
    "PF00920": ["lyase", "dehydratase", "branched_chain_aa"],  # ILVD_EDD_N - Dihydroxy-acid dehydratase
    "PF24877": ["lyase", "dehydratase", "branched_chain_aa"],  # ILV_EDD_C - Dihydroxy-acid dehydratase C-term

    # Phosphatases
    "PF14582": ["phosphatase", "hydrolase", "metal_binding"],  # Metallophos_3 - calcineurin superfamily

    # Replication
    "PF22090": ["replication"],  # Gins51_C - GINS complex subunit

    # Zinc finger domains
    "PF07754": ["zinc_binding", "metal_binding"],  # HVO_2753_ZBP - Small zinc finger protein

    # Secretion
    "PF03911": ["secretion_component", "sec_pathway", "membrane"],  # Sec61_beta - Sec61beta family
    "PF03703": ["signaling"],  # bPH_2 - Bacterial PH domain

    # Iron-sulfur cluster assembly
    "PF01458": ["iron_sulfur_biosynthesis", "iron_sulfur"],  # SUFBD_core - SUF system FeS assembly

    # Glycine cleavage
    "PF01597": ["amino_acid_degradation", "one_carbon_metabolism"],  # GCV_H - Glycine cleavage H-protein

    # Chromosome condensation
    "PF02616": ["cell_division", "chromosome_partitioning"],  # SMC_ScpA - Segregation and condensation ScpA

    # DNA-binding archaeal histone
    "PF01918": ["dna_binding", "regulator"],  # Alba - archaeal chromatin protein

    # Cell division/sporulation
    "PF07441": ["cell_division"],  # BofA - SigmaK-factor processing regulatory

    # ATP-grasp ligases
    "PF08443": ["ligase", "atp_binding"],  # RimK - RimK-like ATP-grasp domain

    # -------------------------------------------------------------------------
    # BATCH 4: MORE UNMAPPED DOMAINS (45-count range)
    # -------------------------------------------------------------------------
    # RNA processing and defense
    "PF01936": ["nuclease", "rna_processing"],  # NYN domain - PIN-like RNase

    # Phosphatases
    "PF00300": ["phosphatase", "hydrolase"],  # His_Phos_1 - Histidine phosphatase

    # Membrane proteins
    "PF04011": ["membrane"],  # LemA family
    "PF13768": ["cell_surface", "protein_binding"],  # VWA_3 - von Willebrand factor type A
    "PF13421": ["membrane", "signaling"],  # Band_7_1 - SPFH domain

    # Ligases
    "PF16177": ["ligase", "central_metabolism"],  # ACAS_N - Acetyl-CoA synthetase N-terminus

    # Proteases
    "PF06550": ["protease", "membrane"],  # SPP - Signal-peptide peptidase
    "PF19583": ["hydrolase", "beta_lactamase"],  # ODP - ODP family beta lactamase

    # Macro domain (ADP-ribose binding)
    "PF01661": ["regulatory", "signaling"],  # Macro domain

    # Translation release factors
    "PF18859": ["translation_factor", "translation"],  # acVLRF1 - release factor

    # Glycosyltransferases
    "PF22627": ["glycosyltransferase", "carbohydrate_active"],  # AglB_core-like - archaeal glycosylation

    # NTPases/pyrophosphatases
    "PF03819": ["hydrolase", "nucleotide_metabolism"],  # MazG - NTP pyrophosphohydrolase

    # Lyases (alliinase)
    "PF04864": ["lyase", "amino_acid_metabolism"],  # Alliinase_C

    # AAA domains
    "PF13555": ["aaa_domain", "atp_binding"],  # AAA_29 - P-loop containing AAA

    # Regulatory domains
    "PF02938": ["regulatory"],  # GAD domain

    # ChlI - magnesium chelatase
    "PF13541": ["chelatase", "atp_binding", "heme_biosynthesis"],  # ChlI - Mg-chelatase subunit

    # -------------------------------------------------------------------------
    # BATCH 5-6: MORE UNMAPPED DOMAINS
    # -------------------------------------------------------------------------
    # Conjugation/secretion
    "PF19044": ["conjugation", "atp_binding"],  # P-loop_TraG - TraG P-loop domain

    # Phosphatases
    "PF14378": ["phosphatase", "hydrolase"],  # PAP2_3 - PAP2 superfamily

    # Hydrogenases
    "PF02662": ["hydrogenase", "iron_sulfur"],  # FlpD - Methyl-viologen-reducing hydrogenase delta subunit

    # AAA domains
    "PF13245": ["aaa_domain", "atp_binding", "atpase"],  # AAA_19
    "PF12846": ["aaa_domain", "atp_binding", "atpase"],  # AAA_10
    "PF13500": ["aaa_domain", "atp_binding", "atpase"],  # AAA_26

    # Iron-sulfur cluster
    "PF04068": ["iron_sulfur", "ferredoxin", "metal_binding"],  # Fer4_RLI - Possible Fer4-like domain

    # Metal binding/stress
    "PF03091": ["metal_binding", "stress_response"],  # CutA1 - divalent ion tolerance

    # Cytochrome biogenesis
    "PF13386": ["heme_biosynthesis", "membrane"],  # DsbD_2 - Cytochrome C biogenesis

    # Amino acid metabolism
    "PF21571": ["oxidoreductase", "amino_acid_metabolism"],  # ArgZ-like_C_1st - Arginine dihydrolase

    # Histones/DNA binding
    "PF00125": ["dna_binding", "regulator"],  # Histone - Core histone

    # Membrane proteins
    "PF01914": ["membrane", "transporter"],  # MarC - MarC family integral membrane protein

    # Fatty acid synthesis
    "PF08659": ["fatty_acid_synthesis", "lipid_metabolism"],  # KR - ketoreductase domain

    # Histidine biosynthesis
    "PF08029": ["transferase", "histidine_biosynthesis"],  # HisG_C - HisG C-terminal

    # Phosphatases (PHP-related)
    "PF13263": ["phosphatase", "hydrolase"],  # PHP_C - PHP-associated

    # Helicase domains
    "PF13361": ["helicase", "dna_repair"],  # UvrD_C - UvrD-like helicase C-terminal

    # tRNA modification
    "PF21238": ["trna_modification", "isomerase"],  # Pus10_C - Pseudouridine synthase

    # DNA repair
    "PF12826": ["dna_binding", "dna_repair"],  # HHH_2 - Helix-hairpin-helix motif

    # Amino acid biosynthesis
    "PF14698": ["lyase", "amino_acid_biosynthesis"],  # ASL_C2 - Argininosuccinate lyase

    # Lipid metabolism (archaeal)
    "PF01864": ["transferase", "lipid_metabolism"],  # CarS-like - CDP-archaeol synthase

    # Hypothetical/uncharacterized
    "PF20582": ["hypothetical"],  # UPF0758_N

    # Cobalamin-independent methionine synthase
    "PF01717": ["transferase", "amino_acid_biosynthesis", "one_carbon_metabolism"],  # Meth_synt_2

    # Dehydrogenases
    "PF04455": ["oxidoreductase", "dehydrogenase", "amino_acid_metabolism"],  # Saccharop_dh_N - LOR/SDH

    # Oxidoreductases
    "PF22725": ["oxidoreductase"],  # GFO_IDH_MocA_C3
    "PF02894": ["oxidoreductase"],  # GFO_IDH_MocA_C

    # Signaling domains
    "PF08495": ["signaling"],  # FIST N domain

    # Lipoate protein ligase
    "PF10437": ["ligase", "cofactor_binding"],  # Lip_prot_lig_C - Bacterial lipoate protein ligase

    # Glycolysis/gluconeogenesis
    "PF06560": ["isomerase", "glycolysis"],  # GPI - Glucose-6-phosphate isomerase

    # Toxins
    "PF05707": ["toxin", "membrane"],  # Zot - Zonular occludens toxin

    # Aldolases
    "PF00596": ["lyase", "aldolase", "central_metabolism"],  # Aldolase_II - Class II Aldolase

    # RNA processing
    "PF04032": ["rnase", "rna_processing"],  # Rpr2 - RNase P Rpr2 subunit

    # Hypothetical
    "PF09383": ["hypothetical"],  # NIL domain

    # Isomerase domains
    "PF13580": ["isomerase", "carbohydrate_active"],  # SIS_2 - SIS domain

    # Kinases
    "PF02224": ["kinase", "nucleotide_metabolism"],  # Cytidylate_kin

    # Aromatic amino acid metabolism
    "PF00697": ["isomerase", "aromatic_aa_metabolism"],  # PRAI - PRA isomerase

    # Repeat domains
    "PF13576": ["repeat_domain"],  # Pentapeptide_3 - Pentapeptide repeats

    # Hydrolases
    "PF00884": ["hydrolase", "sulfur_metabolism"],  # Sulfatase

    # Transporters
    "PF01594": ["transporter", "signaling"],  # AI-2E_transport - AI-2E family transporter
    "PF00909": ["transporter", "nitrogen_metabolism"],  # Ammonium_transp - Ammonium Transporter
    "PF01769": ["transporter", "metal_transporter", "magnesium_binding"],  # MgtE - Divalent cation transporter
    "PF03083": ["transporter", "sugar_transporter"],  # MtN3_slv - Sugar efflux transporter

    # Transcription factors
    "PF09341": ["transcription", "regulator"],  # Pcc1 - Transcription factor Pcc1

    # Motif domains
    "PF01493": ["nucleotide_binding"],  # GXGXG motif

    # Chromosomal proteins
    "PF05854": ["dna_binding", "regulator"],  # MC1 - Non-histone chromosomal protein

    # Cobalt transport
    "PF02361": ["transporter", "cobalt_binding"],  # CbiQ - Cobalt transport protein

    # Iron-sulfur cluster assembly
    "PF01592": ["iron_sulfur_biosynthesis", "iron_sulfur"],  # NifU_N - NifU-like N terminal

    # Purine biosynthesis
    "PF02843": ["ligase", "purine_metabolism"],  # GARS_C
    "PF01808": ["transferase", "purine_metabolism"],  # AICARFT_IMPCHas - bienzyme

    # Nucleotidyltransferases
    "PF01148": ["nucleotidyltransferase", "lipid_metabolism"],  # CTP_transf_1

    # RNases
    "PF13691": ["rnase", "hydrolase"],  # Lactamase_B_4 - tRNase Z endonuclease

    # Ligases
    "PF03099": ["ligase", "cofactor_biosynthesis"],  # BPL_LplA_LipB - Biotin/lipoate ligase

    # Replication
    "PF14551": ["helicase", "replication"],  # MCM_N - MCM N-terminal domain

    # Cell surface
    "PF17210": ["cell_surface", "repeat_domain"],  # SdrD_B - SdrD B-like domain

    # Kinases
    "PF01202": ["kinase", "aromatic_aa_metabolism"],  # SKI - Shikimate kinase

    # One-carbon metabolism
    "PF00763": ["oxidoreductase", "one_carbon_metabolism"],  # THF_DHG_CYH - THF dehydrogenase/cyclohydrolase
    "PF02882": ["oxidoreductase", "one_carbon_metabolism", "nad_binding"],  # THF_DHG_CYH_C

    # DNA repair
    "PF10108": ["nuclease", "dna_repair", "replication"],  # DNA_pol_B_exo2 - exonuclease

    # Repeat domains
    "PF01436": ["repeat_domain"],  # NHL repeat

    # Cofactor biosynthesis
    "PF19288": ["lyase", "cofactor_biosynthesis"],  # CofH_C - F420 biosynthesis

    # Cell surface
    "PF22352": ["cell_surface", "repeat_domain"],  # K319L-like_PKD

    # Regulatory
    "PF13715": ["regulatory"],  # CarbopepD_reg_2

    # Transporters (cyclin M)
    "PF01595": ["transporter", "membrane", "metal_transporter"],  # CNNM - Cyclin M transmembrane

    # tRNA modification
    "PF01980": ["methyltransferase", "trna_modification"],  # TrmO_N - tRNA-methyltransferase

    # Translation quality control
    "PF09382": ["translation", "rna_binding"],  # RQC domain

    # Helicase domains
    "PF00570": ["helicase", "dna_repair"],  # HRDC domain

    # Isoprenoid biosynthesis
    "PF01884": ["transferase", "lipid_metabolism"],  # PcrB family

    # Metallopeptidases
    "PF19289": ["protease", "metal_binding"],  # PmbA_TldD_3rd
    "PF01523": ["protease", "metal_binding"],  # PmbA_TldD_1st
    "PF01863": ["protease", "metal_binding"],  # YgjP-like - metallopeptidase

    # Iron-sulfur cluster binding
    "PF13186": ["iron_sulfur", "radical_sam"],  # SPASM domain

    # Membrane repeats
    "PF04193": ["membrane", "repeat_domain"],  # PQ-loop repeat

    # Hydrogenase maturation
    "PF01455": ["hydrogenase_maturation"],  # HupF_HypC family

    # Translation factors
    "PF18854": ["translation_factor", "translation"],  # baeRF_family10 - release factor

    # Oxidoreductases
    "PF02852": ["oxidoreductase", "fad_binding"],  # Pyr_redox_dim

    # Restriction-modification
    "PF04313": ["restriction_enzyme", "restriction_modification"],  # HSDR_N - Type I restriction

    # NTPases
    "PF03266": ["hydrolase", "atp_binding"],  # NTPase_1

    # Hypothetical
    "PF03685": ["hypothetical"],  # UPF0147

    # Histidine biosynthesis
    "PF01503": ["hydrolase", "histidine_biosynthesis"],  # PRA-PH - Phosphoribosyl-ATP pyrophosphohydrolase

    # TCA cycle
    "PF01989": ["lyase", "tca_cycle"],  # AcnX_swivel_put - Aconitase X swivel

    # Purine biosynthesis
    "PF02700": ["transferase", "purine_metabolism"],  # PurS - FGAM synthase

    # RNA processing
    "PF11969": ["hydrolase", "rna_processing"],  # DcpS_C - mRNA decapping
    "PF08494": ["helicase", "rna_binding"],  # DEAD_assoc - DEAD/H associated

    # Mobile elements
    "PF03050": ["transposase", "mobile_element"],  # DDE_Tnp_IS66

    # Helicase domains
    "PF19306": ["helicase", "dna_binding", "winged_helix"],  # WH_Lhr - Large helicase-related

    # CRISPR-associated
    "PF09484": ["crispr_associated", "defense_system"],  # Cas_TM1802

    # Topoisomerase
    "PF00521": ["topoisomerase", "replication"],  # DNA_topoisoIV - DNA gyrase/topoisomerase IV

    # Stress response
    "PF06146": ["stress_response", "phosphate_transporter"],  # PsiE - Phosphate-starvation-inducible

    # Metal binding
    "PF03692": ["metal_binding", "zinc_binding"],  # CxxCxxCC - zinc/iron-chelating

    # Pterin biosynthesis
    "PF10131": ["cofactor_biosynthesis", "membrane"],  # PTPS_related

    # Regulatory (TOBE)
    "PF03459": ["regulatory", "transporter"],  # TOBE domain

    # Nucleotide metabolism
    "PF01230": ["hydrolase", "nucleotide_metabolism"],  # HIT domain

    # Aromatic metabolism
    "PF00857": ["hydrolase", "aromatic_aa_metabolism"],  # Isochorismatase domain

    # Cell cycle
    "PF12895": ["cell_division"],  # ANAPC3 - Anaphase-promoting complex

    # Amino acid degradation
    "PF00208": ["oxidoreductase", "dehydrogenase", "amino_acid_degradation"],  # ELFV_dehydrog

    # DNA repair
    "PF03215": ["replication", "dna_repair", "atp_binding"],  # Rad17 - P-loop domain

    # Antibiotic resistance
    "PF02673": ["transporter", "antibiotic_resistance"],  # BacA - Bacitracin resistance

    # Nucleotide metabolism
    "PF02867": ["oxidoreductase", "reductase", "nucleotide_metabolism"],  # Ribonuc_red_lgC - Ribonucleotide reductase

    # Binding domains
    "PF12849": ["binding"],  # PBP_like_2 - PBP superfamily domain
    "PF00017": ["protein_binding"],  # SH2 - SH2 domain
    "PF00018": ["protein_binding"],  # SH3_1 - SH3 domain

    # Defense systems
    "PF20720": ["defense_system", "atp_binding"],  # nSTAND3 - Novel STAND NTPase

    # DNA repair
    "PF04002": ["dna_repair", "hydrolase"],  # RadC - RadC-like JAB domain

    # -------------------------------------------------------------------------
    # BATCH 7-8: MORE UNMAPPED DOMAINS
    # -------------------------------------------------------------------------
    # Carbohydrate metabolism
    "PF00908": ["isomerase", "carbohydrate_active"],  # dTDP_sugar_isom - dTDP-4-dehydrorhamnose epimerase
    "PF01950": ["phosphatase", "gluconeogenesis"],  # FBPase_3 - Fructose-1,6-bisphosphatase
    "PF00316": ["phosphatase", "gluconeogenesis"],  # FBPase - Fructose-1-6-bisphosphatase N-terminal
    "PF18913": ["phosphatase", "gluconeogenesis"],  # FBPase_C - Fructose-1-6-bisphosphatase C-terminal
    "PF00923": ["transferase", "pentose_phosphate"],  # TAL_FSA - Transaldolase
    "PF00343": ["transferase", "carbohydrate_active"],  # Phosphorylase - Carbohydrate phosphorylase

    # Hypothetical/uncharacterized
    "PF21748": ["hypothetical"],  # UPF0150

    # Proteases
    "PF10263": ["protease", "dna_binding"],  # SprT-like family

    # Translation
    "PF18006": ["trna_synthetase", "translation", "selenium_metabolism"],  # SepRS_C - O-phosphoseryl-tRNA synthetase

    # Sulfur metabolism
    "PF01507": ["oxidoreductase", "sulfur_metabolism"],  # PAPS_reduct - Phosphoadenosine phosphosulfate reductase

    # Central metabolism
    "PF00390": ["oxidoreductase", "central_metabolism", "nad_binding"],  # malic - Malic enzyme
    "PF03949": ["oxidoreductase", "central_metabolism", "nad_binding"],  # Malic_M - Malic enzyme NAD binding

    # Cold shock
    "PF00313": ["cold_shock", "stress_response", "rna_binding"],  # CSD - Cold-shock DNA-binding

    # RNA binding
    "PF01985": ["rna_binding"],  # CRS1_YhbY - CRM domain

    # Transporters
    "PF03030": ["transporter", "energy_metabolism"],  # H_PPase - Inorganic H+ pyrophosphatase

    # Hydantoinase
    "PF01968": ["hydrolase", "amino_acid_metabolism"],  # Hydantoinase_A

    # Iron-sulfur proteins
    "PF06397": ["iron_sulfur", "metal_binding"],  # Desulfoferrod_N - Desulfoferrodoxin
    "PF01880": ["iron_sulfur", "oxidative_stress"],  # Desulfoferrodox - Desulfoferrodoxin

    # One-carbon metabolism
    "PF08714": ["one_carbon_metabolism", "oxidoreductase"],  # Fae - Formaldehyde-activating enzyme

    # Phosphatases
    "PF01937": ["phosphatase", "hydrolase"],  # ARMT1-like_dom - Damage-control phosphatase

    # Carbohydrate-binding
    "PF10633": ["carbohydrate_binding", "carbohydrate_active"],  # NPCBM_assoc

    # Cofactor biosynthesis
    "PF02558": ["oxidoreductase", "cofactor_biosynthesis"],  # ApbA - Ketopantoate reductase

    # Helicase/chromatin remodeling
    "PF00176": ["helicase", "atp_binding", "regulator"],  # SNF2-rel_dom

    # Hypothetical
    "PF22167": ["hypothetical"],  # PH0730-like_N

    # DNA replication
    "PF04042": ["dna_polymerase", "replication"],  # DNA_pol_E_B - DNA polymerase epsilon subunit B

    # tRNA modification
    "PF22082": ["trna_modification", "sulfur_metabolism"],  # TtuA_LIM_N - 2-thiouridine synthetase

    # Cobalamin metabolism
    "PF08267": ["transferase", "amino_acid_biosynthesis", "one_carbon_metabolism"],  # Meth_synt_1 - Cobalamin-independent synthase

    # Iron-sulfur cluster assembly
    "PF01106": ["iron_sulfur_biosynthesis", "iron_sulfur"],  # NifU

    # Toxin-antitoxin
    "PF02452": ["toxin_domain", "defense_system"],  # PemK_toxin - PemK/MazF-like toxin

    # Translation
    "PF17835": ["ribosomal_protein", "translation", "gtp_binding"],  # NOG1_N - NOG1 N-terminal

    # Molybdenum metabolism
    "PF02634": ["molybdenum_binding", "oxidoreductase"],  # FdhD-NarQ - formate dehydrogenase accessory

    # Amino acid metabolism
    "PF01168": ["isomerase", "amino_acid_metabolism", "plp_binding"],  # Ala_racemase_N - Alanine racemase

    # Hypothetical
    "PF03479": ["hypothetical"],  # PCC - Plants and Prokaryotes Conserved

    # DNA repair
    "PF02481": ["dna_repair", "recombinational_repair"],  # DNA_processg_A - DNA recombination-mediator

    # Nitrogen regulation
    "PF00543": ["regulatory", "nitrogen_metabolism"],  # P-II - Nitrogen regulatory protein

    # Cell surface
    "PF17957": ["cell_surface", "repeat_domain"],  # Big_7 - Bacterial Ig domain

    # Amino acid degradation
    "PF02812": ["oxidoreductase", "amino_acid_degradation"],  # ELFV_dehydrog_N - Glu/Leu/Phe/Val dehydrogenase

    # Secretion
    "PF12666": ["secretion_component", "t3ss_component"],  # PrgI - T3SS needle protein

    # Nucleases
    "PF00565": ["nuclease", "hydrolase"],  # SNase - Staphylococcal nuclease

    # Chaperones
    "PF19026": ["chaperone"],  # UBA_HYPK

    # Ribosome biogenesis
    "PF06858": ["ribosomal_protein", "gtp_binding", "translation"],  # NOG1 - Nucleolar GTP-binding protein

    # Antibiotic resistance
    "PF04892": ["antibiotic_resistance", "membrane"],  # VanZ - VanZ family

    # Metal transporters
    "PF01297": ["transporter", "abc_transporter", "zinc_binding"],  # ZnuA - Zinc-uptake component

    # Transketolase
    "PF00456": ["transferase", "pentose_phosphate", "thiamine_binding"],  # Transketolase_N
    "PF22613": ["transferase", "pentose_phosphate"],  # Transketolase_C_1

    # Beta propeller
    "PF09826": ["repeat_domain"],  # Beta_propel - Beta propeller domain

    # Formylmethanofuran dehydrogenase - in methanogenesis but also Wood-Ljungdahl
    "PF02663": ["molybdenum_binding", "archaeal_one_carbon", "wood_ljungdahl"],  # FmdE - Formylmethanofuran dehydrogenase

    # Cell surface
    "PF00028": ["adhesin", "cell_surface", "repeat_domain"],  # Cadherin - Cadherin domain
    "PF00039": ["cell_surface", "repeat_domain"],  # fn1 - Fibronectin type I domain
    "PF00040": ["cell_surface", "repeat_domain"],  # fn2 - Fibronectin type II domain
    "PF00041": ["cell_surface", "repeat_domain"],  # fn3 - Fibronectin type III domain
    "PF00047": ["cell_surface", "repeat_domain"],  # ig - Immunoglobulin domain

    # Phosphoesterases
    "PF02833": ["phosphatase", "hydrolase"],  # DHHA2 domain

    # Halogenases
    "PF04820": ["halogenase", "oxidoreductase", "fad_binding"],  # Trp_halogenase

    # Transposases
    "PF03400": ["transposase", "mobile_element"],  # DDE_Tnp_IS1

    # Transporters
    "PF03600": ["transporter", "tca_cycle"],  # CitMHS - Citrate transporter

    # Hypothetical
    "PF08734": ["hypothetical"],  # GYD domain
    "PF16285": ["hypothetical"],  # DUF4931_N

    # Replication
    "PF17244": ["replication", "rna_binding"],  # CDC24_OB3 - Cell division control protein OB domain

    # RNA binding
    "PF22891": ["rna_binding", "translation"],  # KH_PNO1_2nd

    # CRISPR-associated
    "PF06023": ["crispr_associated", "nuclease", "defense_system"],  # Csa1 - CRISPR exonuclease

    # tRNA modification
    "PF18093": ["methyltransferase", "trna_modification"],  # Trm5_N

    # Transporters
    "PF01758": ["transporter", "symporter", "membrane"],  # SBF - Sodium Bile acid symporter

    # Coenzyme metabolism
    "PF14574": ["kinase", "atp_binding"],  # RACo_C_ter - ASKHA domain
    "PF17651": ["binding"],  # Raco_middle

    # Translation
    "PF08207": ["translation_factor", "translation"],  # EFP_N - Elongation factor P

    # Immune/defense
    "PF20773": ["protease", "defense_system"],  # InhA-like_MAM - Immune inhibitor A

    # Metal transporters
    "PF01988": ["transporter", "metal_transporter", "iron_binding"],  # VIT1 family

    # Translation factors
    "PF01287": ["translation_factor", "translation"],  # eIF-5a - Elongation factor 5A

    # Nucleotide metabolism
    "PF01928": ["hydrolase", "nucleotide_metabolism"],  # CYTH domain

    # Membrane proteins
    "PF25145": ["membrane", "signaling"],  # NfeD1b_N

    # CO dehydrogenase
    "PF19436": ["co_oxidation", "oxidoreductase"],  # ACS_CODH_B_C - ACS/CODH beta subunit

    # Dehydrogenases
    "PF25137": ["oxidoreductase", "dehydrogenase", "iron_binding"],  # ADH_Fe_C - Fe-containing alcohol dehydrogenase

    # Nucleases
    "PF08378": ["nuclease", "hydrolase"],  # NERD - Nuclease-related domain

    # Membrane proteins
    "PF08308": ["membrane"],  # PEGA domain

    # Glycosylation
    "PF04756": ["glycosyltransferase", "membrane"],  # OST3_OST6 - oligosaccharyltransferase

    # Acetyltransferases
    "PF24894": ["acetyltransferase", "transferase"],  # Hexapep_GlmU

    # Metal binding
    "PF04434": ["zinc_binding", "metal_binding"],  # SWIM zinc finger

    # Membrane proteins
    "PF06271": ["membrane"],  # RDD family

    # Cell wall
    "PF03023": ["transporter", "peptidoglycan", "cell_wall"],  # MurJ - Lipid II flippase

    # NTPases
    "PF07693": ["atp_binding", "atpase"],  # KAP_NTPase - KAP family P-loop

    # Dehydrogenases
    "PF02866": ["oxidoreductase", "dehydrogenase"],  # Ldh_1_C - lactate/malate dehydrogenase
    "PF00056": ["oxidoreductase", "dehydrogenase", "nad_binding"],  # Ldh_1_N - lactate/malate dehydrogenase

    # Thioredoxin fold
    "PF01323": ["oxidoreductase", "thioredoxin"],  # DSBA - DSBA-like thioredoxin

    # Toxin-antitoxin
    "PF13638": ["toxin_domain", "pin_domain"],  # PIN_4

    # Molybdenum cofactor
    "PF06463": ["molybdenum_binding", "cofactor_biosynthesis"],  # Mob_synth_C

    # GTPases
    "PF00025": ["gtpase", "gtp_binding"],  # Arf - ADP-ribosylation factor

    # Hydrogenase maturation
    "PF17788": ["hydrogenase_maturation", "transferase"],  # HypF_C - HypF Kae1-like

    # Lipid metabolism
    "PF04191": ["methyltransferase", "phospholipid_metabolism"],  # PEMT - Phospholipid methyltransferase

    # tRNA modification
    "PF01938": ["trna_modification", "rna_binding"],  # TRAM domain

    # Dehydrogenases
    "PF01210": ["oxidoreductase", "dehydrogenase", "lipid_metabolism"],  # NAD_Gly3P_dh_N - Glycerol-3-phosphate dehydrogenase

    # Chemotaxis
    "PF04509": ["chemotaxis", "signaling"],  # CheC

    # Carbohydrate metabolism
    "PF02358": ["phosphatase", "carbohydrate_active"],  # Trehalose_PPase

    # Defense systems
    "PF13676": ["defense_system", "signaling"],  # TIR_2 - TIR domain
    "PF01582": ["defense_system", "signaling"],  # TIR domain

    # Metal transporters
    "PF01169": ["transporter", "metal_transporter", "calcium_binding"],  # GDT1 - Divalent cation/proton antiporter

    # Nickel metabolism
    "PF01969": ["nickel_binding", "hydrogenase_maturation"],  # Ni_insertion

    # Metal binding
    "PF01906": ["metal_binding", "heavy_metal_resistance"],  # YbjQ_1 - heavy-metal-binding

    # Binding domains
    "PF12727": ["binding"],  # PBP_like - PBP superfamily

    # Glutathione metabolism
    "PF02955": ["ligase", "atp_binding"],  # GSH-S_ATP - Prokaryotic glutathione synthetase

    # Molybdenum cofactor
    "PF01967": ["molybdenum_binding", "cofactor_biosynthesis"],  # MoaC

    # Nucleases
    "PF01541": ["nuclease", "hydrolase"],  # GIY-YIG - GIY-YIG catalytic domain

    # Protease inhibitors
    "PF00079": ["protease", "regulator"],  # Serpin - serine protease inhibitor

    # Transporters
    "PF02690": ["transporter", "phosphate_transporter"],  # Na_Pi_cotrans - Na+/Pi-cotransporter

    # DNA repair
    "PF10576": ["iron_sulfur", "dna_repair"],  # EndIII_4Fe-2S - Endonuclease III iron-sulfur

    # GTPases
    "PF00071": ["gtpase", "gtp_binding", "signaling"],  # Ras family

    # Molybdenum cofactor
    "PF03454": ["molybdenum_binding", "cofactor_biosynthesis"],  # MoeA_C
    "PF03453": ["molybdenum_binding", "cofactor_biosynthesis"],  # MoeA_N
    "PF03473": ["molybdenum_binding", "cofactor_biosynthesis"],  # MOSC domain

    # Fatty acid synthesis
    "PF13561": ["oxidoreductase", "fatty_acid_synthesis"],  # adh_short_C2 - Enoyl-ACP reductase

    # Translation
    "PF07581": ["translation", "ribosomal_protein"],  # Glug motif

    # Aldolases
    "PF01116": ["lyase", "aldolase", "glycolysis"],  # F_bP_aldolase - Fructose-bisphosphate aldolase class-II

    # Anticodon binding
    "PF03147": ["iron_sulfur", "trna_synthetase"],  # FDX-ACB - Ferredoxin-fold anticodon binding

    # Polysaccharide biosynthesis
    "PF04230": ["transferase", "carbohydrate_active"],  # PS_pyruv_trans - Polysaccharide pyruvyl transferase

    # Adhesins
    "PF05658": ["adhesin", "cell_surface"],  # YadA_head

    # Iron-sulfur clusters
    "PF13370": ["iron_sulfur", "ferredoxin", "metal_binding", "electron_transport"],  # Fer4_13
    "PF17179": ["iron_sulfur", "ferredoxin", "metal_binding", "electron_transport"],  # Fer4_22

    # Transcription
    "PF02969": ["transcription", "regulator"],  # TAF - TATA box binding protein associated factor

    # Hypothetical
    "PF10370": ["hypothetical"],  # Rv2993c-like_N

    # RNA repair
    "PF10127": ["dna_polymerase", "dna_repair"],  # RlaP - RNA repair DNA polymerase beta

    # NOTE: PF04422/PF04432 (FrhB) defined in HYDROGENASES section above

    # Secretion
    "PF05753": ["secretion_component", "sec_pathway"],  # TRAP_beta - Translocon-associated protein

    # Stress response
    "PF03168": ["stress_response"],  # LEA_2 - Late embryogenesis abundant protein

    # Virulence
    "PF13310": ["virulence"],  # Virulence_RhuM

    # Proteases
    "PF13975": ["protease", "mobile_element"],  # gag-asp_proteas - retroviral aspartyl protease

    # Kinases
    "PF04008": ["kinase", "nucleotide_metabolism"],  # Adenosine_kin - Adenosine kinase

    # Hypothetical
    "PF04289": ["hypothetical"],  # DUF447_N

    # Restriction-modification
    "PF12161": ["restriction_modification", "methyltransferase_rm"],  # HsdM_N

    # Exonucleases
    "PF00929": ["nuclease", "hydrolase"],  # RNase_T - Exonuclease

    # Defense/stress
    "PF13630": ["defense_system"],  # SdpI - SdpI protein family

    # DNA repair
    "PF08713": ["dna_repair", "base_excision_repair"],  # DNA_alkylation - DNA alkylation repair

    # AAA domains
    "PF13479": ["aaa_domain", "atp_binding", "atpase"],  # AAA_24

    # Decarboxylases
    "PF02627": ["lyase", "decarboxylase"],  # CMD - Carboxymuconolactone decarboxylase

    # CRISPR
    "PF03787": ["crispr_associated", "defense_system"],  # RAMPs - RAMP superfamily

    # Translation
    "PF03966": ["translation", "methyltransferase"],  # Trm112p

    # Repeat domains
    "PF00805": ["repeat_domain"],  # Pentapeptide - Pentapeptide repeats

    # Hypothetical with hints
    "PF01784": ["hypothetical"],  # DUF34_NIF3

    # Lactate racemase
    "PF21113": ["isomerase", "central_metabolism"],  # LarA_C - Lactate racemase
    "PF09861": ["isomerase", "central_metabolism"],  # Lar_N - Lactate racemase

    # Metal transporters
    "PF16916": ["transporter", "zinc_binding"],  # ZT_dimer - Zinc Transporter dimerization

    # CO dehydrogenase
    "PF10133": ["co_oxidation", "nickel_binding"],  # CooT - CO dehydrogenase accessory

    # Molybdenum cofactor
    "PF03205": ["molybdenum_binding", "cofactor_biosynthesis"],  # MobB

    # Hypothetical
    "PF09851": ["hypothetical"],  # SHOCT

    # Aldolases
    "PF08013": ["lyase", "aldolase", "carbohydrate_active"],  # GatZ_KbaZ-like - D-tagatose aldolase

    # Deaminases
    "PF13382": ["hydrolase", "deaminase"],  # Adenine_deam_C - Adenine deaminase

    # Cofactor biosynthesis (SAM)
    "PF20257": ["transferase", "cofactor_biosynthesis"],  # SAM_HAT_C - SAM hydroxide adenosyltransferase

    # Phosphoesterases
    "PF13286": ["phosphatase", "hydrolase"],  # HD_assoc - Phosphohydrolase-associated

    # Proteases
    "PF19290": ["protease", "metal_binding"],  # PmbA_TldD_2nd

    # Membrane/signaling
    "PF01957": ["membrane", "signaling"],  # NfeD

    # Cytochrome biogenesis
    "PF02683": ["heme_biosynthesis", "membrane"],  # DsbD_TM

    # Transcription
    "PF21715": ["regulator", "dna_binding"],  # CggR_N - CggR DNA binding

    # Transferases
    "PF02543": ["transferase", "amino_acid_metabolism"],  # Carbam_trans_N - Carbamoyltransferase
    "PF16861": ["transferase", "amino_acid_metabolism"],  # Carbam_trans_C

    # -------------------------------------------------------------------------
    # BATCH 9: ALTIARCHAEOTA HIGH-ABUNDANCE UNMAPPED DOMAINS
    # -------------------------------------------------------------------------
    # Transporters and membrane proteins
    "PF01891": ["transporter", "cobalt_binding", "membrane"],  # CbiM - Cobalt uptake transmembrane
    "PF04039": ["transporter", "antiporter", "membrane"],  # MnhB - Na+/H+ antiporter subunit
    "PF01899": ["transporter", "antiporter", "membrane"],  # MNHE - Na+/H+ antiporter subunit
    "PF03334": ["transporter", "antiporter", "membrane"],  # PhaG_MnhG_YufB - Na+/H+ antiporter subunit
    "PF00474": ["transporter", "symporter", "membrane"],  # SSF - Sodium:solute symporter
    "PF03471": ["transporter", "membrane"],  # CorC_HlyC - Transporter associated domain
    "PF07155": ["transporter", "cofactor_biosynthesis"],  # ECF-ribofla_trS - ECF riboflavin transporter
    "PF00654": ["transporter", "ion_channel", "membrane"],  # Voltage_CLC - Chloride channel
    "PF05978": ["transporter", "ion_channel", "membrane", "regulatory"],  # UNC-93 - Ion channel regulatory
    "PF01810": ["transporter", "amino_acid_transporter", "membrane"],  # LysE - Lysine translocator
    "PF00893": ["transporter", "efflux_pump", "multidrug_resistance", "membrane"],  # Multi_Drug_Res - SMR

    # Phosphatases
    "PF00481": ["phosphatase", "hydrolase", "regulatory"],  # PP2C - Protein phosphatase 2C
    "PF13672": ["phosphatase", "hydrolase"],  # PP2C_2 - Protein phosphatase 2C
    "PF01451": ["phosphatase", "hydrolase"],  # LMWPc - Low molecular weight phosphotyrosine phosphatase

    # Signaling and regulatory domains
    "PF13023": ["signaling", "phosphodiesterase"],  # HD_3 - HD domain (phosphodiesterase)
    "PF10442": ["signaling"],  # FIST_C - FIST C domain
    "PF04969": ["chaperone", "stress_response"],  # CS - CS domain (co-chaperone)
    "PF08402": ["regulatory", "transporter"],  # TOBE_2 - TOBE domain

    # Oxidoreductases and flavoproteins
    "PF00258": ["oxidoreductase", "fmn_binding", "electron_transport"],  # Flavodoxin_1 - Flavodoxin
    "PF12724": ["oxidoreductase", "fmn_binding", "electron_transport"],  # Flavodoxin_5 - Flavodoxin
    "PF12682": ["oxidoreductase", "fmn_binding", "electron_transport"],  # Flavodoxin_4 - Flavodoxin
    "PF03227": ["oxidoreductase", "reductase", "thioredoxin"],  # GILT - Lysosomal thiol reductase
    "PF13434": ["oxidoreductase", "oxygenase", "amino_acid_metabolism"],  # Lys_Orn_oxgnase - Lysine/ornithine oxygenase
    "PF02943": ["oxidoreductase", "iron_sulfur", "thioredoxin"],  # FeThRed_B - Ferredoxin thioredoxin reductase
    "PF00171": ["oxidoreductase", "dehydrogenase"],  # Aldedh - Aldehyde dehydrogenase
    "PF00676": ["oxidoreductase", "dehydrogenase", "central_metabolism"],  # E1_dh - Dehydrogenase E1 component
    "PF00317": ["oxidoreductase", "reductase", "nucleotide_metabolism"],  # Ribonuc_red_lgN - Ribonucleotide reductase

    # Iron-sulfur cluster proteins
    "PF04324": ["iron_sulfur", "ferredoxin", "2fe2s", "metal_binding"],  # Fer2_BFD - BFD-like [2Fe-2S] binding
    "PF13746": ["iron_sulfur", "ferredoxin", "metal_binding", "electron_transport"],  # Fer4_18 - 4Fe-4S dicluster
    "PF14691": ["iron_sulfur", "ferredoxin", "metal_binding", "oxidoreductase"],  # Fer4_20 - 4Fe-4S cluster

    # DNA repair and nucleases
    "PF01844": ["nuclease", "hydrolase", "dna_repair"],  # HNH - HNH endonuclease
    "PF14279": ["nuclease", "hydrolase"],  # HNH_5 - HNH endonuclease
    "PF02144": ["dna_repair", "replication"],  # Rad1 - Repair protein Rad1/Rec1/Rad17
    "PF02151": ["dna_repair", "nucleotide_excision_repair"],  # UVR - UvrB/uvrC motif
    "PF17760": ["dna_repair", "nucleotide_excision_repair"],  # UvrA_inter - UvrA interaction domain
    "PF17757": ["dna_repair", "nucleotide_excision_repair"],  # UvrB_inter - UvrB interaction domain
    "PF12344": ["dna_repair", "nucleotide_excision_repair"],  # UvrB - Ultra-violet resistance protein B
    "PF22920": ["rnase", "dna_repair"],  # UvrC_RNaseH - UvrC RNase H-like domain
    "PF08459": ["rnase", "nuclease", "dna_repair"],  # UvrC_RNaseH_dom - UvrC RNase H endonuclease
    "PF01939": ["nuclease", "dna_repair"],  # NucS_C - Endonuclease NucS C-terminal
    "PF21003": ["nuclease", "dna_repair"],  # NucS_N - Endonuclease NucS N-terminal

    # Transposases and mobile elements
    "PF12762": ["transposase", "mobile_element"],  # DDE_Tnp_IS1595 - Transposase domain
    "PF01797": ["transposase", "mobile_element"],  # Y1_Tnp - Transposase IS200 like
    "PF07592": ["transposase", "mobile_element"],  # DDE_Tnp_ISAZ013 - Transposase DDE domain
    "PF13546": ["transposase", "mobile_element", "nuclease"],  # DDE_5 - DDE superfamily endonuclease

    # Restriction-modification systems
    "PF12950": ["restriction_modification", "restriction_enzyme"],  # TaqI_C - TaqI-like specificity domain
    "PF22679": ["restriction_modification", "restriction_enzyme"],  # T1R_D3-like - Type I restriction enzyme
    "PF08463": ["restriction_modification", "restriction_enzyme"],  # EcoEI_R_C - EcoEI R protein C-terminal
    "PF15514": ["restriction_modification", "restriction_enzyme"],  # ThaI - Restriction endonuclease ThaI
    "PF03235": ["restriction_modification", "restriction_enzyme"],  # GmrSD_N - GmrSD restriction endonuclease
    "PF09491": ["restriction_modification", "restriction_enzyme"],  # RE_AlwI - AlwI restriction endonuclease
    "PF06300": ["restriction_modification", "restriction_enzyme"],  # Tsp45I - Tsp45I restriction enzyme

    # CRISPR-associated
    "PF01905": ["crispr_associated", "regulator", "defense_system"],  # DevR - CRISPR-associated DevR/Csa2

    # Repeat domains
    "PF07721": ["repeat_domain", "tpr_repeat", "protein_binding"],  # TPR_4 - Tetratricopeptide repeat
    "PF00415": ["repeat_domain", "regulatory"],  # RCC1 - Regulator of chromosome condensation repeat
    "PF13540": ["repeat_domain", "regulatory"],  # RCC1_2 - Regulator of chromosome condensation repeat

    # Isomerases and rotamases
    "PF00639": ["isomerase", "chaperone"],  # Rotamase - PPIC-type peptidyl-prolyl isomerase
    "PF13145": ["isomerase", "chaperone"],  # Rotamase_2 - PPIC-type PPIASE
    "PF13616": ["isomerase", "chaperone"],  # Rotamase_3 - PPIC-type PPIASE

    # Lyases and central metabolism
    "PF01293": ["lyase", "gluconeogenesis", "central_metabolism"],  # PEPCK_ATP - PEP carboxykinase
    "PF10415": ["lyase", "tca_cycle", "central_metabolism"],  # FumaraseC_C - Fumarase C C-terminus
    "PF18376": ["lyase", "decarboxylase", "isoprenoid_biosynthesis"],  # MDD_C - Mevalonate decarboxylase C-term

    # Transferases
    "PF00581": ["sulfurtransferase", "transferase", "sulfur_metabolism"],  # Rhodanese - Rhodanese-like
    "PF01206": ["sulfurtransferase", "transferase", "sulfur_metabolism"],  # TusA - Sulfurtransferase TusA
    "PF02277": ["transferase", "nucleotide_metabolism"],  # DBI_PRT - Phosphoribosyltransferase
    "PF02348": ["nucleotidyltransferase", "lipid_metabolism"],  # CTP_transf_3 - Cytidylyltransferase
    "PF04140": ["methyltransferase", "lipid_metabolism"],  # ICMT - Isoprenylcysteine carboxyl methyltransferase
    "PF05219": ["methyltransferase", "transferase"],  # DREV - DREV methyltransferase
    "PF04101": ["glycosyltransferase", "carbohydrate_active"],  # Glyco_tran_28_C - Glycosyltransferase 28

    # Kinases
    "PF01712": ["kinase", "nucleotide_metabolism"],  # dNK - Deoxynucleoside kinase
    "PF08645": ["kinase", "nuclease", "dna_repair"],  # PNK3P - Polynucleotide kinase 3 phosphatase

    # Hydrolases
    "PF00484": ["hydrolase", "lyase", "zinc_binding"],  # Pro_CA - Carbonic anhydrase
    "PF20393": ["hydrolase", "lyase"],  # Pro_CA_2 - Putative carbonic anhydrase
    "PF12917": ["hydrolase", "nucleotide_metabolism"],  # YfbR-like - 5'-deoxynucleotidase
    "PF02585": ["hydrolase", "carbohydrate_active"],  # PIG-L - GlcNAc-PI de-N-acetylase
    "PF16347": ["hydrolase", "sulfur_metabolism"],  # SGSH_C - N-sulphoglucosamine sulphohydrolase
    "PF02633": ["hydrolase", "amidase"],  # Creatininase - Creatinine amidohydrolase
    "PF00723": ["glycosidase", "hydrolase", "carbohydrate_active"],  # Glyco_hydro_15 - Glycosyl hydrolases 15

    # Stress response and phage shock
    "PF04024": ["stress_response", "membrane"],  # PspC - PspC domain (phage shock)
    "PF04012": ["stress_response", "membrane"],  # PspA_IM30 - PspA/IM30 family
    "PF22743": ["stress_response"],  # PspAA - PspA-Associated protein
    "PF03631": ["stress_response", "membrane"],  # Virul_fac_BrkB - Virulence factor BrkB

    # Replication and helicase
    "PF21120": ["helicase", "replication", "dna_binding"],  # MCM_WH_arc - Archaeal MCM winged-helix
    "PF05916": ["replication", "helicase"],  # Sld5 - GINS complex helical bundle domain
    "PF05872": ["helicase", "dna_repair"],  # HerA_C - Helicase HerA-like C-terminal
    "PF08519": ["replication", "dna_binding"],  # RFC1 - Replication factor RFC1 C terminal

    # AAA+ ATPases
    "PF13604": ["aaa_domain", "atp_binding", "atpase"],  # AAA_30
    "PF13514": ["aaa_domain", "atp_binding", "atpase"],  # AAA_27
    "PF13521": ["aaa_domain", "atp_binding", "atpase"],  # AAA_28
    "PF09336": ["aaa_domain", "atp_binding"],  # Vps4_C - Vps4 oligomerisation domain

    # Translation and tRNA-related
    "PF02403": ["trna_synthetase", "translation"],  # Seryl_tRNA_N - Seryl-tRNA synthetase N-terminal
    "PF22625": ["rna_binding", "rnase"],  # ECR1_N_2 - Exosome RNA binding protein
    "PF13636": ["rna_binding", "methyltransferase"],  # Methyltranf_PUA - RNA-binding PUA domain

    # Cell surface and adhesion
    "PF12245": ["cell_surface", "repeat_domain"],  # Big_3_2 - Bacterial Ig-like domain
    "PF13750": ["cell_surface", "repeat_domain"],  # Big_3_3 - Bacterial Ig-like domain
    "PF08487": ["cell_surface"],  # VIT - Vault protein inter-alpha-trypsin domain
    "PF01833": ["cell_surface", "adhesin"],  # TIG - IPT/TIG domain

    # Cofactor biosynthesis
    "PF06968": ["cofactor_biosynthesis", "thiamine_biosynthesis"],  # BATS - Biotin/Thiamin Synthesis
    "PF01243": ["oxidoreductase", "cofactor_biosynthesis", "pyridoxal_biosynthesis"],  # PNPOx_N - PNP oxidase
    "PF02581": ["thiamine_biosynthesis", "cofactor_biosynthesis"],  # TMP-TENI - Thiamine monophosphate synthase

    # Membrane and gas vesicle
    "PF00741": ["membrane", "cell_surface"],  # Gas_vesicle - Gas vesicle protein
    "PF03073": ["membrane", "signaling"],  # TspO_MBR - TspO/MBR family
    "PF11127": ["membrane"],  # YgaP-like_TM - Inner membrane protein

    # Ligases
    "PF15632": ["ligase", "atp_binding"],  # ATPgrasp_Ter - ATP-grasp in biosynthetic pathway
    "PF01268": ["ligase", "one_carbon_metabolism", "atp_binding"],  # FTHFS - Formate-tetrahydrofolate ligase

    # Archaeal regulators
    "PF08350": ["regulator", "archaeal_one_carbon"],  # FilR1_middle - archaeal C1 metabolism regulator
    # NOTE: PF13244/PF20501 (MbhD/MbhE) defined in HYDROGENASES section above

    # Protease inhibitors and proteases
    "PF05922": ["protease", "regulator"],  # Inhibitor_I9 - Peptidase inhibitor I9
    "PF13365": ["protease", "hydrolase"],  # Trypsin_2 - Trypsin-like peptidase

    # Transferases - galactose metabolism
    "PF02744": ["transferase", "carbohydrate_active"],  # GalP_UDP_tr_C - Galactose-1-phosphate uridyl transferase

    # Signal recognition
    "PF09439": ["signal_recognition", "gtp_binding"],  # SRPRB - SRP receptor beta subunit

    # Secretion
    "PF04021": ["secretion_component", "sec_pathway"],  # Class_IIIsignal - Class III signal peptide
    "PF02810": ["secretion_component", "sec_pathway"],  # SEC-C - SEC-C motif

    # RNA processing and modification
    "PF22640": ["isomerase", "carbohydrate_active"],  # ManC_GMP_beta-helix - MannoseP isomerase-like
    "PF04703": ["archaeal_one_carbon", "one_carbon_metabolism"],  # FaeA - FaeA-like protein (formaldehyde-activating)

    # DNA mismatch repair
    "PF00488": ["dna_repair", "mismatch_repair"],  # MutS_V - MutS domain V

    # MgtC family - virulence
    "PF02308": ["membrane", "stress_response"],  # MgtC - MgtC family

    # Deaminases
    "PF06559": ["hydrolase", "deaminase", "nucleotide_metabolism"],  # DCD_N - dCTP deaminase N-terminal
    "PF22569": ["hydrolase", "deaminase", "nucleotide_metabolism"],  # DCD_C - dCTP deaminase C-terminal

    # SAM-related
    "PF01887": ["transferase", "cofactor_biosynthesis"],  # SAM_HAT_N - SAM hydroxide adenosyltransferase N-term

    # Exosortase
    "PF09721": ["membrane", "cell_surface"],  # Exosortase_EpsH - Transmembrane exosortase

    # Glycine cleavage
    "PF02347": ["amino_acid_degradation", "one_carbon_metabolism"],  # GDC-P - Glycine cleavage system P-protein

    # DUF domains and hypothetical (high counts)
    "PF22665": ["hypothetical"],  # DUF6293_C - DUF6293 C-terminal winged helix
    "PF10543": ["hypothetical"],  # ORF6N - ORF6N domain
    "PF22357": ["hypothetical"],  # AF1548-like_C
    "PF20009": ["hypothetical"],  # GEVED - GEVED domain
    "PF09378": ["hypothetical"],  # HAS-barrel - HAS barrel domain
    "PF19810": ["hypothetical"],  # HFX_2341_N
    "PF05899": ["hydrolase"],  # Cupin_3 - EutQ-like cupin domain
    "PF01987": ["hypothetical"],  # AIM24 - Mitochondrial biogenesis AIM24
    "PF22689": ["transferase", "purine_metabolism"],  # FGAR-AT_PurM_N-like - FGAM synthase related
    "PF04457": ["rna_binding"],  # MJ1316 - RNA cyclic group end recognition
    "PF05401": ["methyltransferase", "transferase"],  # NodS - Nodulation protein S (N-methyltransferase)
    "PF24883": ["hypothetical"],  # NPHP3_N - Nephrocystin 3 N-terminal

    # Acetyltransferases (GST family)
    "PF13417": ["transferase", "glutathione"],  # GST_N_3 - Glutathione S-transferase N-terminal

    # Additional sulfur metabolism
    "PF14269": ["transferase", "sulfur_metabolism"],  # Arylsulfotran_2 - Arylsulfotransferase
    "PF05935": ["transferase", "sulfur_metabolism"],  # Arylsulfotrans - Arylsulfotransferase

    # Sigma factors
    "PF07638": ["sigma_factor", "regulator"],  # Sigma70_ECF - ECF sigma factor

    # Toxins
    "PF03495": ["toxin", "calcium_binding"],  # Binary_toxB - Binary toxin B Ca-binding
    "PF20126": ["toxin_domain"],  # TumE - Toxin TumE

    # Additional AAA domains
    "PF13087": ["aaa_domain", "atp_binding", "atpase"],  # AAA_12

    # Cell division
    "PF04824": ["cell_division", "chromosome_partitioning"],  # Rad21_Rec8 - Conserved region of Rad21/Rec8
    "PF15511": ["cell_division", "dna_binding"],  # CENP-T_C - Centromere kinetochore component

    # Lipid metabolism
    "PF00614": ["hydrolase", "lipase", "signaling"],  # PLDc - Phospholipase D
    "PF13396": ["nuclease", "lipase"],  # PLDc_N - Phospholipase_D-nuclease N-terminal
    "PF00781": ["kinase", "lipid_metabolism"],  # DAGK_cat - Diacylglycerol kinase

    # Various enzymes and domains
    "PF07931": ["kinase", "antibiotic_resistance"],  # CPT - Chloramphenicol phosphotransferase
    "PF22614": ["ion_channel", "regulatory", "calcium_binding"],  # Slo-like_RCK - Ca-activated K+ channel
    "PF07786": ["acetyltransferase", "transferase"],  # HGSNAT_cat - N-acetyltransferase
    "PF03060": ["oxidoreductase", "monooxygenase"],  # NMO - Nitronate monooxygenase
    "PF13593": ["transporter", "membrane"],  # SBF_like - SBF-like CPA transporter
    "PF22277": ["iron_binding", "metal_binding"],  # EncFtn-like - Encapsulin-like ferritin
    "PF12349": ["membrane", "regulatory"],  # Sterol-sensing - Sterol-sensing domain
    "PF02018": ["carbohydrate_binding", "carbohydrate_active"],  # CBM_4_9 - Carbohydrate binding module
    "PF05768": ["oxidoreductase", "glutaredoxin"],  # Glrx-like - Glutaredoxin-like domain
    "PF01161": ["binding", "lipid_metabolism"],  # PBP - Phosphatidylethanolamine-binding protein
    "PF00078": ["dna_polymerase", "mobile_element"],  # RVT_1 - Reverse transcriptase
    "PF02225": ["protease", "hydrolase"],  # PA - PA domain (protease associated)
    "PF05729": ["atpase", "atp_binding", "signaling"],  # NACHT - NACHT domain (NTPase)
    "PF08666": ["cell_surface", "adhesin"],  # SAF - SAF domain
    "PF04900": ["ribosomal_protein", "translation"],  # Fcf1 - Fcf1
    "PF03641": ["lyase", "decarboxylase", "amino_acid_metabolism"],  # Lysine_decarbox - Lysine decarboxylase
    "PF14681": ["transferase", "nucleotide_metabolism"],  # UPRTase - Uracil phosphoribosyltransferase
    "PF09136": ["glycosidase", "hydrolase", "carbohydrate_active"],  # Glucodextran_B - Glucodextranase
    "PF14667": ["glycosyltransferase", "carbohydrate_active"],  # Polysacc_synt_C - Polysaccharide biosynthesis
    "PF13289": ["hydrolase", "regulator"],  # SIR2_2 - SIR2-like domain (deacetylase)

    # LPS biosynthesis
    "PF06176": ["kinase", "lps_biosynthesis"],  # WaaY - LPS core biosynthesis protein

    # More amino acid metabolism
    "PF01262": ["oxidoreductase", "dehydrogenase", "amino_acid_metabolism"],  # AlaDh_PNT_C - Alanine dehydrogenase
    "PF04898": ["oxidoreductase", "amino_acid_metabolism"],  # Glu_syn_central - Glutamate synthase central

    # SbcC/RAD50 DNA repair
    "PF13558": ["dna_repair", "atp_binding"],  # SbcC_Walker_B - SbcC/RAD50-like Walker B

    # MrpF - membrane
    "PF04066": ["transporter", "membrane"],  # MrpF_PhaF - Multiple resistance and pH regulation

    # Mechanosensitive channel
    "PF01741": ["transporter", "ion_channel", "membrane"],  # MscL - Large-conductance mechanosensitive channel

    # Conjugation
    "PF13728": ["conjugation", "mobile_element"],  # TraF - F plasmid transfer operon

    # Decarboxylases
    "PF02261": ["lyase", "decarboxylase", "amino_acid_metabolism"],  # Asp_decarbox - Aspartate decarboxylase

    # Kelch repeat domains
    "PF13418": ["repeat_domain", "oxidoreductase"],  # Kelch_4 - Galactose oxidase central domain
    "PF13415": ["repeat_domain", "oxidoreductase"],  # Kelch_3 - Galactose oxidase central domain

    # Glycosyl hydrolases
    "PF04041": ["glycosidase", "hydrolase", "carbohydrate_active"],  # Glyco_hydro_130 - beta-1,4-mannooligosaccharide phosphorylase

    # Ornithine metabolism
    "PF02423": ["isomerase", "amino_acid_metabolism"],  # OCD_Mu_crystall - Ornithine cyclodeaminase

    # Signaling
    "PF13190": ["signaling"],  # PDGLE - PDGLE domain

    # SelR domain - selenium/methionine metabolism
    "PF01641": ["oxidoreductase", "selenium_metabolism"],  # SelR - SelR domain (methionine sulfoxide reductase)

    # TAT signal
    "PF10518": ["secretion_component", "tat_pathway"],  # TAT_signal - TAT pathway signal sequence

    # MgtE N-terminal
    "PF03448": ["transporter", "magnesium_binding", "regulatory"],  # MgtE_N - MgtE intracellular N domain

    # Defense system (Thoeris)
    "PF08937": ["defense_system", "signaling"],  # ThsB_TIR - Thoeris protein ThsB TIR-like

    # Antitoxin
    "PF13338": ["antitoxin_domain", "regulator"],  # AbiEi_4 - AbiEi antitoxin

    # NeuB - sialic acid biosynthesis
    "PF03102": ["transferase", "carbohydrate_active"],  # NeuB - NeuB family (sialic acid synthase)

    # WD40-like repeat
    "PF07676": ["repeat_domain", "wd40_repeat"],  # PD40 - WD40-like Beta Propeller Repeat

    # GTP cyclohydrolase
    "PF00925": ["hydrolase", "cofactor_biosynthesis"],  # GTP_cyclohydro2 - GTP cyclohydrolase II

    # Sulfotransferase
    "PF13469": ["transferase", "sulfur_metabolism"],  # Sulfotransfer_3 - Sulfotransferase family

    # ATP-sulfurylase
    "PF01747": ["transferase", "atp_binding", "sulfur_metabolism"],  # ATP-sulfurylase

    # Trypsin
    "PF00089": ["protease", "hydrolase"],  # Trypsin

    # Aldo/keto reductase
    "PF00248": ["oxidoreductase", "reductase"],  # Aldo_ket_red - Aldo/keto reductase family

    # PEP carboxylase
    "PF14010": ["lyase", "central_metabolism"],  # PEPcase_2 - Phosphoenolpyruvate carboxylase

    # ASCH domain
    "PF04266": ["rna_binding"],  # ASCH - ASCH domain (RNA-binding)

    # Aminoglycoside resistance
    "PF10706": ["transferase", "antibiotic_resistance"],  # Aminoglyc_resit - Aminoglycoside adenylyltransferase

    # HHH motif
    "PF12836": ["dna_binding", "dna_repair"],  # HHH_3 - Helix-hairpin-helix motif

    # N-acetylmannosamine epimerase
    "PF04131": ["isomerase", "epimerase", "carbohydrate_active"],  # NanE - N-acetylmannosamine-6-phosphate epimerase

    # Substrate binding proteins
    "PF13343": ["transporter", "abc_transporter", "abc_substrate_binding"],  # SBP_bac_6 - Bacterial SBP
    "PF13416": ["transporter", "abc_transporter", "abc_substrate_binding"],  # SBP_bac_8 - Bacterial SBP

    # Pyrroline-5-carboxylate reductase
    "PF14748": ["oxidoreductase", "reductase", "amino_acid_metabolism"],  # P5CR_dimer - P5C reductase

    # CheW chemotaxis
    "PF01584": ["chemotaxis", "signaling"],  # CheW - CheW-like domain

    # ERCC3/RAD25/XPB helicase
    "PF16203": ["helicase", "dna_repair", "nucleotide_excision_repair"],  # ERCC3_RAD25_C - XPB C-terminal helicase

    # CAP domain
    "PF00188": ["cell_surface"],  # CAP - Cysteine-rich secretory protein family

    # Additional lower-count but important domains
    "PF14377": ["protein_binding"],  # UBM - Ubiquitin binding region
    "PF04982": ["membrane"],  # TM_HPP - HPP transmembrane region
    "PF10080": ["iron_sulfur", "membrane"],  # FtrD-like - Membrane iron-sulfur protein
    "PF00050": ["protease", "regulator"],  # Kazal_1 - Kazal-type serine protease inhibitor
    "PF06819": ["protease", "hydrolase"],  # Arc_PepC - Archaeal Peptidase A24 C-terminal
    "PF21433": ["kinase", "isoprenoid_biosynthesis"],  # M3K_C - Mevalonate-3-kinase C-terminal
    "PF13378": ["lyase", "central_metabolism"],  # MR_MLE_C - Enolase C-terminal domain-like
    "PF01954": ["hypothetical"],  # AF2212-like
    "PF09084": ["thiamine_biosynthesis", "cofactor_biosynthesis"],  # NMT1 - NMT1/THI5 like
    "PF13376": ["defense_system"],  # OmdA - Bacteriocin-protection
    "PF14811": ["hypothetical"],  # TPD - Protein of unknown function TPD
    "PF20990": ["membrane", "hypothetical"],  # DUF2207_C - Predicted membrane protein
    "PF09516": ["restriction_modification", "restriction_enzyme"],  # RE_CfrBI - CfrBI restriction endonuclease
    "PF22633": ["carbohydrate_binding"],  # F5_F8_type_C_2 - NedA-like galactose-binding
    "PF22557": ["rna_binding", "dna_binding"],  # DuOB - Dual OB-containing domain
    "PF14470": ["signaling"],  # bPH_3 - Bacterial PH domain
    "PF00754": ["carbohydrate_binding"],  # F5_F8_type_C - F5/8 type C domain

    # -------------------------------------------------------------------------
    # BATCH 10: HIGH-ABUNDANCE UNMAPPED DOMAINS (Altiarchaeota dataset)
    # -------------------------------------------------------------------------
    # ABC transporter related
    "PF20275": ["transporter", "abc_transporter"],  # CTD10 - ABC-3C systems C-terminal domain

    # Sugar/carbohydrate metabolism
    "PF05523": ["transferase", "carbohydrate_active"],  # FdtA - WxcM-like C-terminal
    "PF10509": ["kinase", "carbohydrate_active"],  # GalKase_gal_bdg - Galactokinase galactose-binding
    "PF05116": ["phosphatase", "hydrolase", "carbohydrate_active"],  # S6PP - Sucrose-6F-phosphate phosphohydrolase
    "PF00342": ["isomerase", "glycolysis", "central_metabolism"],  # PGI - Phosphoglucose isomerase
    "PF01204": ["glycosidase", "hydrolase", "carbohydrate_active"],  # Trehalase
    "PF00722": ["glycosidase", "hydrolase", "carbohydrate_active"],  # Glyco_hydro_16
    "PF03537": ["glycosidase", "hydrolase", "carbohydrate_active"],  # Glyco_hydro_114
    "PF24517": ["carbohydrate_binding", "carbohydrate_active"],  # CBM96 - Carbohydrate-binding module 96
    "PF22422": ["glycosidase", "hydrolase", "carbohydrate_active"],  # MGH1-like_GH - Mannosylglycerate hydrolase

    # RNA/nucleotide processing
    "PF22600": ["nucleotidyltransferase", "rna_processing"],  # MTPAP-like_central - Poly(A) RNA polymerase
    "PF00075": ["rnase", "hydrolase", "dna_repair"],  # RNase_H
    "PF00849": ["trna_modification", "isomerase"],  # PseudoU_synth_2 - RNA pseudouridylate synthase
    "PF05635": ["rnase", "rna_processing"],  # 23S_rRNA_IVP - 23S rRNA-intervening sequence protein

    # Regulatory domains
    "PF06445": ["regulatory", "binding"],  # GyrI-like - Small molecule binding domain
    "PF11495": ["regulator", "transcription_factor", "dna_binding"],  # Regulator_TrmB - Archaeal transcriptional regulator
    "PF07085": ["regulatory"],  # DRTGG domain
    "PF13492": ["signaling", "regulatory"],  # GAF_3 - GAF domain
    "PF13411": ["regulator", "transcription_factor", "dna_binding", "helix_turn_helix"],  # MerR_1 - MerR HTH family
    "PF03965": ["regulator", "antibiotic_resistance"],  # Penicillinase_R - Penicillinase repressor
    "PF21086": ["regulatory", "amino_acid_metabolism"],  # ACT_PSP_2 - ACT domain

    # Mobile elements/phage
    "PF07929": ["phage_related", "hypothetical"],  # PRiA4_ORF3 - Plasmid pRiA4b ORF-3-like
    "PF01695": ["transposase", "mobile_element", "atp_binding"],  # IstB_IS21 - IstB-like ATP binding
    "PF01548": ["transposase", "mobile_element"],  # DEDD_Tnp_IS110 - Transposase
    "PF04986": ["transposase", "mobile_element"],  # Y2_Tnp - Putative transposase
    "PF13359": ["transposase", "mobile_element", "nuclease"],  # DDE_Tnp_4 - DDE superfamily endonuclease
    "PF22483": ["transposase", "mobile_element"],  # Mu-transpos_C_2 - Mu transposase C-terminal

    # Signaling domains
    "PF18171": ["signaling"],  # LSDAT_prok - SLOG in TRPM, prokaryote
    "PF02743": ["signaling", "sensor_kinase"],  # dCache_1 - Cache domain

    # Restriction-modification
    "PF22722": ["restriction_enzyme", "restriction_modification"],  # NA-iREase1 - NACHT-associated inactive Restriction Endonuclease
    "PF03852": ["nuclease", "dna_repair", "restriction_modification"],  # Vsr - DNA mismatch endonuclease
    "PF04556": ["restriction_enzyme", "restriction_modification"],  # DpnII restriction endonuclease
    "PF07510": ["restriction_enzyme", "restriction_modification"],  # GmrSD_C - GmrSD restriction endonuclease
    "PF09520": ["restriction_enzyme", "restriction_modification"],  # RE_TdeIII - Type II restriction endonuclease
    "PF09549": ["restriction_enzyme", "restriction_modification"],  # RE_Bpu10I - Bpu10I restriction endonuclease
    "PF09568": ["restriction_enzyme", "restriction_modification"],  # RE_MjaI - MjaI restriction endonuclease
    "PF06044": ["restriction_enzyme", "restriction_modification"],  # DpnI - Dam-replacing family
    "PF17726": ["restriction_modification", "dna_binding", "helix_turn_helix"],  # DpnI_C - Dam-replacing HTH domain
    "PF13588": ["restriction_enzyme", "restriction_modification"],  # HSDR_N_2 - Type I restriction enzyme
    "PF18643": ["restriction_enzyme", "restriction_modification"],  # RE_BsaWI - BsaWI restriction endonuclease

    # RNA binding and translation
    "PF01997": ["rna_binding", "translation"],  # Translin family
    "PF04410": ["rna_binding", "rrna_modification"],  # Gar1 - Gar1/Naf1 RNA binding region
    "PF09190": ["trna_synthetase", "translation", "rna_binding"],  # DALR_2 - DALR anticodon binding domain
    "PF19269": ["trna_synthetase", "translation", "rna_binding"],  # Anticodon_2 - Anticodon binding domain
    "PF13742": ["rna_binding", "translation"],  # tRNA_anti_2 - OB-fold nucleic acid binding
    "PF14306": ["rna_binding"],  # PUA_2 - PUA-like domain
    "PF09180": ["trna_synthetase", "translation"],  # ProRS-C_1 - Prolyl-tRNA synthetase C-terminal
    "PF21266": ["rna_binding", "rna_processing"],  # RRP4_S1 - RRP4 S1 domain
    "PF21189": ["rna_binding"],  # PHA02142 - OB-fold domain

    # Gas vesicle proteins
    "PF06386": ["membrane", "cell_surface"],  # GvpL_GvpF - Gas vesicle synthesis protein

    # Cell surface proteins
    "PF05762": ["cell_surface", "protein_binding"],  # VWA_CoxE - VWA domain containing CoxE-like
    "PF17963": ["cell_surface", "repeat_domain"],  # Big_9 - Bacterial Ig domain
    "PF19077": ["cell_surface", "repeat_domain"],  # Big_13 - Bacterial Ig-like domain
    "PF03640": ["cell_surface", "repeat_domain", "secreted"],  # Lipoprotein_15 - Secreted repeat
    "PF05345": ["cell_surface", "repeat_domain"],  # He_PIG - Putative Ig domain

    # Sulfur/methionine metabolism
    "PF00685": ["sulfurtransferase", "transferase", "sulfur_metabolism"],  # Sulfotransfer_1 - Sulfotransferase
    "PF01625": ["oxidoreductase", "oxidative_stress", "sulfur_metabolism"],  # PMSR - Peptide methionine sulfoxide reductase

    # Kinases
    "PF00406": ["kinase", "nucleotide_metabolism", "atp_binding"],  # ADK - Adenylate kinase
    "PF05191": ["kinase", "regulatory"],  # ADK_lid - Adenylate kinase active site lid

    # Proteases and peptidases
    "PF12770": ["protease", "hydrolase"],  # CHAT domain
    "PF00077": ["protease", "mobile_element"],  # RVP - Retroviral aspartyl protease
    "PF08126": ["protease", "regulatory"],  # Propeptide_C25 - Propeptide
    "PF03413": ["protease", "regulatory"],  # PepSY - Peptidase propeptide and YPEB domain

    # Amino acid metabolism
    "PF01571": ["transferase", "one_carbon_metabolism", "amino_acid_degradation"],  # GCV_T - Glycine cleavage aminomethyltransferase
    "PF08669": ["amino_acid_degradation", "one_carbon_metabolism"],  # GCV_T_C - Glycine cleavage T-protein C-terminal
    "PF21478": ["oxidoreductase", "amino_acid_degradation"],  # GcvP2_C - Glycine dehydrogenase C-terminal
    "PF02388": ["transferase", "peptidoglycan", "cell_wall"],  # FemAB - peptidoglycan synthesis
    "PF00962": ["hydrolase", "deaminase", "nucleotide_metabolism"],  # A_deaminase - Adenosine deaminase

    # Cofactor biosynthesis
    "PF00899": ["cofactor_biosynthesis", "thiamine_biosynthesis", "ligase"],  # ThiF - ThiF family
    "PF03232": ["cofactor_biosynthesis", "oxidoreductase"],  # COQ7 - Ubiquinone biosynthesis protein
    "PF02597": ["cofactor_biosynthesis", "thiamine_biosynthesis"],  # ThiS - ThiS family
    "PF05402": ["cofactor_biosynthesis", "pqq_binding"],  # PqqD - Coenzyme PQQ synthesis protein D
    "PF05165": ["hydrolase", "cofactor_biosynthesis"],  # GCH_III - GTP cyclohydrolase III
    "PF12900": ["oxidoreductase", "cofactor_biosynthesis", "pyridoxal_biosynthesis"],  # Pyridox_ox_2

    # Oxidoreductases
    "PF01995": ["regulatory", "nitrogen_metabolism"],  # NRD1_2 - NrpR regulatory domains
    "PF08534": ["oxidoreductase", "thioredoxin", "oxidative_stress"],  # Redoxin
    "PF00578": ["oxidoreductase", "peroxiredoxin", "oxidative_stress"],  # AhpC-TSA - peroxiredoxin family
    "PF10417": ["oxidoreductase", "peroxiredoxin", "oxidative_stress"],  # 1-cysPrx_C - 1-Cys peroxiredoxin C-terminal
    "PF00725": ["oxidoreductase", "dehydrogenase", "fatty_acid_degradation"],  # 3HCDH - 3-hydroxyacyl-CoA dehydrogenase
    "PF12641": ["oxidoreductase", "fmn_binding", "electron_transport"],  # Flavodoxin_3
    "PF16912": ["oxidoreductase", "dehydrogenase"],  # Glu_dehyd_C - Glucose dehydrogenase C-terminus
    "PF05222": ["oxidoreductase", "dehydrogenase", "amino_acid_metabolism"],  # AlaDh_PNT_N - Alanine dehydrogenase
    "PF16653": ["oxidoreductase", "dehydrogenase", "amino_acid_metabolism"],  # Sacchrp_dh_C

    # Hypothetical/uncharacterized
    "PF03699": ["hypothetical"],  # UPF0182 - Uncharacterised protein family
    "PF11028": ["glycosyltransferase", "membrane"],  # TMEM260-like - Protein O-mannosyl-transferase
    "PF05958": ["methyltransferase", "trna_modification"],  # tRNA_U5-meth_tr - tRNA uracil-5-methyltransferase
    "PF21941": ["hypothetical"],  # SMEK_N - SMEK domain
    "PF04139": ["cell_division", "dna_repair"],  # Rad9 - checkpoint protein
    "PF03741": ["transporter", "metal_transporter", "membrane"],  # TerC - Integral membrane protein TerC family
    "PF23743": ["repeat_domain", "wd40_repeat"],  # Beta-prop_BBS7 - BBS7 beta-propeller
    "PF08668": ["hydrolase", "phosphodiesterase"],  # HDOD domain
    "PF12544": ["isomerase", "radical_sam", "iron_sulfur"],  # LAM_C - Lysine-2,3-aminomutase
    "PF16881": ["radical_sam", "iron_sulfur", "cofactor_biosynthesis"],  # LIAS_N - lipoyl synthase N-terminal
    "PF13513": ["repeat_domain", "heat_repeat"],  # HEAT_EZ - HEAT-like repeat
    "PF04945": ["hypothetical"],  # YHS domain
    "PF04199": ["lyase"],  # Cyclase - Putative cyclase
    "PF03372": ["nuclease", "phosphatase", "hydrolase"],  # Exo_endo_phos - Endonuclease/Exonuclease/phosphatase
    "PF13280": ["regulatory", "defense_system"],  # WYL domain
    "PF15630": ["cell_division", "dna_binding"],  # CENP-S protein
    "PF16289": ["toxin_domain", "pin_domain"],  # PIN_12 domain
    "PF09250": ["dna_polymerase", "primase", "replication"],  # Prim-Pol - Bifunctional DNA primase/polymerase
    "PF05597": ["lipid_metabolism", "membrane"],  # Phasin - Poly(hydroxyalcanoate) granule associated protein
    "PF01057": ["mobile_element", "helicase", "atp_binding"],  # Parvo_NS1 - Parvovirus non-structural protein
    "Kelch_5": ["repeat_domain", "kelch_repeat"],  # Kelch_5 - Kelch repeat
    "PF22296": ["hypothetical"],  # bAvd-like
    "PF04464": ["glycosyltransferase", "carbohydrate_active", "cell_wall"],  # Glyphos_transf - glycerophosphotransferase
    "PF04151": ["protease", "regulatory", "secreted"],  # PPC - Bacterial pre-peptidase C-terminal domain

    # Defense systems
    "PF09407": ["antitoxin_domain", "defense_system"],  # AbiEi_1 - AbiEi antitoxin
    "PF05488": ["effector_domain", "t6ss_component"],  # PAAR_motif - T6SS effector tip

    # Chemotaxis
    "PF01739": ["methyltransferase", "chemotaxis"],  # CheR - CheR methyltransferase SAM binding
    "PF03975": ["chemotaxis", "signaling"],  # CheD - chemotactic sensory transduction
    "PF01339": ["hydrolase", "chemotaxis"],  # CheB_methylest - CheB methylesterase

    # Stress/osmotic response
    "PF02566": ["oxidative_stress", "stress_response"],  # OsmC - OsmC-like protein

    # SAM/methylation
    "PF02675": ["lyase", "decarboxylase", "sam_binding"],  # AdoMet_dc - S-adenosylmethionine decarboxylase

    # Metal binding
    "PF09360": ["iron_sulfur", "metal_binding"],  # zf-CDGSH - Iron-binding zinc finger
    "PF23475": ["zinc_binding", "metal_binding"],  # zf-Tbcl_FmdE - FmdE treble clef zinc finger

    # Transferases
    "PF13712": ["glycosyltransferase", "transferase", "carbohydrate_active"],  # Glyco_tranf_2_5
    "PF02397": ["glycosyltransferase", "transferase", "carbohydrate_active"],  # Bac_transf - Bacterial sugar transferase
    "PF01757": ["acetyltransferase", "transferase", "lipid_metabolism"],  # Acyl_transf_3
    "PF04188": ["glycosyltransferase", "membrane"],  # Mannosyl_trans2 - Mannosyltransferase (PIG-V)
    "PF04138": ["glycosyltransferase", "membrane"],  # GtrA_DPMS_TM - GtrA/DPMS transmembrane domain
    "PF14907": ["nucleotidyltransferase", "transferase"],  # NTP_transf_5

    # GST/glutathione
    "PF02798": ["transferase", "oxidative_stress"],  # GST_N - Glutathione S-transferase N-terminal
    "PF00043": ["transferase", "oxidative_stress"],  # GST_C - Glutathione S-transferase, C-terminal domain

    # Tubulin/cell structure
    "PF13809": ["cell_division", "gtp_binding"],  # Tubulin_2 - Tubulin like

    # Membrane proteins
    "PF08294": ["membrane", "transporter"],  # TIM21
    "PF04972": ["membrane", "signaling"],  # BON domain

    # Isomerases
    "PF04126": ["isomerase", "chaperone"],  # Cyclophil_like - Cyclophilin-like

    # DNA repair/replication
    "PF01612": ["nuclease", "dna_polymerase", "replication"],  # DNA_pol_A_exo1 - 3'-5' exonuclease
    "PF00476": ["dna_polymerase", "replication"],  # DNA_pol_A - DNA polymerase family A
    "PF02132": ["dna_repair", "zinc_binding"],  # RecR_ZnF - RecR Cys4-zinc finger
    "PF00817": ["dna_repair", "sos_response"],  # IMS - impB/mucB/samB family
    "PF08696": ["helicase", "nuclease", "replication"],  # Dna2 - DNA replication factor

    # Selenium metabolism
    "PF03841": ["trna_synthetase", "selenium_metabolism"],  # SelA - L-seryl-tRNA selenium transferase

    # Helicase
    "PF05970": ["helicase", "atp_binding", "replication"],  # PIF1-like helicase
    "PF22527": ["helicase", "atp_binding"],  # DEXQc_Suv3 - DEXQ-box helicase

    # Binding proteins
    "PF13458": ["transporter", "abc_transporter", "abc_substrate_binding"],  # Peripla_BP_6
    "PF13433": ["transporter", "abc_transporter", "abc_substrate_binding"],  # Peripla_BP_5

    # Cell wall related
    "PF01510": ["hydrolase", "peptidoglycan", "cell_wall"],  # Amidase_2 - N-acetylmuramoyl-L-alanine amidase

    # Repeat domains
    "PF01391": ["repeat_domain", "cell_surface"],  # Collagen - Collagen triple helix repeat
    "PF13599": ["repeat_domain"],  # Pentapeptide_4
    "PF13371": ["repeat_domain", "tpr_repeat", "protein_binding"],  # TPR_9
    "PF17874": ["repeat_domain", "tpr_repeat", "regulatory"],  # TPR_MalT

    # Calcium binding
    "PF11535": ["calcium_binding", "metal_binding"],  # Calci_bind_CcbP

    # GTPases
    "PF03308": ["gtpase", "gtp_binding", "cobalamin_binding"],  # MeaB - Methylmalonyl Co-A mutase-associated GTPase

    # Various
    "PF03481": ["trna_modification", "cofactor_biosynthesis"],  # Sua5_C - Threonylcarbamoyl-AMP synthase
    "PF03350": ["hypothetical"],  # UPF0114
    "PF00931": ["atp_binding", "signaling", "defense_system"],  # NB-ARC domain
    "PF07648": ["protease", "regulator"],  # Kazal_2 - serine protease inhibitor
    "PF04383": ["dna_binding", "mobile_element"],  # KilA-N domain
    "PF18046": ["isomerase", "chaperone"],  # FKBP26_C - FKBP26 C-terminal
    "PF01094": ["receptor", "signaling", "membrane"],  # ANF_receptor - Receptor family ligand binding
    "PF07228": ["phosphatase", "cell_division"],  # SpoIIE - Stage II sporulation protein E
    "PF12822": ["transporter", "abc_transporter"],  # ECF_trnsprt - ECF transporter
    "PF08889": ["transferase", "carbohydrate_active"],  # WbqC - WbqC-like protein
    "PF12973": ["oxidoreductase", "chromate_resistance"],  # Cupin_7 - ChrR Cupin-like domain
    "PF00190": ["hydrolase", "oxidoreductase"],  # Cupin_1 - Cupin domain
    "PF01239": ["transferase", "isoprenoid_biosynthesis"],  # PPTA - Protein prenyltransferase
    "PF02537": ["membrane", "stress_response"],  # CRCB - Camphor Resistance (CrcB)
    "PF03750": ["crispr_associated", "defense_system"],  # Csm2_III-A - CRISPR Type III-A
    "PF00480": ["kinase", "carbohydrate_active", "regulatory"],  # ROK family
    "PF04143": ["transporter", "sulfate_transporter"],  # Sulf_transp - Sulphur transport
    "PF10531": ["binding"],  # SLBB domain
    "PF15993": ["hypothetical"],  # Fuseless
    "PF03845": ["transporter", "membrane"],  # Spore_permease - Spore germination protein
    "PF07504": ["protease", "regulatory"],  # FTP - Fungalysin/Thermolysin Propeptide
    "PF04208": ["transferase", "archaeal_one_carbon", "one_carbon_metabolism"],  # MtrA - Tetrahydromethanopterin S-methyltransferase
    "PF20774": ["protease", "defense_system"],  # InhA-like_VEG - Immune inhibitor A-like metallopeptidase
    "PF10335": ["nucleotidyltransferase", "hypothetical"],  # DUF294_C
    "PF02900": ["dioxygenase", "oxidoreductase", "aromatic_aa_metabolism"],  # LigB - aromatic ring-opening dioxygenase
    "PF23441": ["oxidoreductase", "dehydrogenase", "nad_binding"],  # SDR - SDR-like rossmann domain
    "PF06728": ["membrane", "glycosyltransferase"],  # PIG-U - GPI transamidase subunit
    "PF07524": ["dna_binding", "regulatory"],  # Bromo_TP - Bromodomain associated
    "PF22912": ["zinc_binding", "dna_polymerase", "replication"],  # zf-DPOE - DNA polymerase-epsilon zinc finger
    "PF04932": ["glycosyltransferase", "carbohydrate_active", "lps_biosynthesis"],  # Wzy_C - O-Antigen ligase
    "PF14815": ["hydrolase", "nucleotide_metabolism"],  # NUDIX_4
    "PF20797": ["rnase", "hydrolase"],  # HepT-like_2
    "PF05552": ["transporter", "ion_channel", "membrane"],  # MS_channel_1st_1 - Mechanosensitive ion channel

    # -------------------------------------------------------------------------
    # BATCH 11: ADDITIONAL UNMAPPED DOMAINS (counts 5-7)
    # -------------------------------------------------------------------------
    # Helicases and DNA-related
    "PF13538": ["helicase", "dna_repair"],  # UvrD_C_2 - UvrD-like helicase C-terminal
    "PF13333": ["integrase", "mobile_element"],  # rve_2 - Integrase core domain

    # Reverse transcriptase/mobile elements
    "PF13456": ["dna_polymerase", "mobile_element"],  # RVT_3 - Reverse transcriptase

    # Zinc fingers
    "PF12773": ["zinc_binding", "metal_binding"],  # DZR - Zinc ribbon domain

    # Cofactor biosynthesis
    "PF02424": ["cofactor_biosynthesis", "fad_biosynthesis"],  # ApbE - thiamine/FAD biosynthesis

    # Dehydrogenases
    "PF16653": ["oxidoreductase", "dehydrogenase", "amino_acid_metabolism"],  # Sacchrp_dh_C - Saccharopine dehydrogenase

    # Gas vesicle
    "PF05121": ["membrane", "cell_surface"],  # GvpK - Gas vesicle protein

    # Nucleotidases
    "PF06941": ["hydrolase", "nucleotide_metabolism"],  # NT5C - 5'-nucleotidase

    # Phosphotransferases
    "PF01627": ["kinase", "two_component"],  # Hpt - Histidine phosphotransferase domain
    "PF01591": ["kinase", "carbohydrate_active"],  # 6PF2K - 6-phosphofructo-2-kinase

    # Isomerases
    "PF07385": ["isomerase", "carbohydrate_active"],  # Lyx_isomer - Sugar isomerase

    # Molybdenum metabolism
    "PF05161": ["cofactor_biosynthesis", "molybdenum_binding"],  # MOFRL - MoaD/ThiS family

    # Protein-protein interaction
    "PF07677": ["binding", "membrane"],  # A2M_recep - Alpha-2-macroglobulin receptor

    # Amino acid metabolism
    "PF00490": ["lyase", "heme_biosynthesis"],  # ALAD - Delta-aminolevulinic acid dehydratase

    # Lipases
    "PF01738": ["hydrolase", "esterase"],  # DLH - Dienelactone hydrolase

    # Ferritins
    "PF13668": ["iron_binding", "metal_binding", "oxidative_stress"],  # Ferritin_2

    # Archaeal domains
    "PF13020": ["oxidoreductase", "hypothetical"],  # NOV_C

    # Hydrogenase related
    "PF07449": ["hydrogenase_maturation"],  # HyaE - Hydrogenase expression/formation protein

    # Acetyl-CoA metabolism
    "PF13336": ["hydrolase", "tca_cycle"],  # AcetylCoA_hyd_C - Acetyl-CoA hydrolase C-terminal

    # DNA binding
    "PF00308": ["replication", "atp_binding", "dna_binding"],  # Bac_DnaA - Bacterial DnaA

    # AAA domains
    "PF13086": ["aaa_domain", "atp_binding", "atpase"],  # AAA_11

    # Citrate synthase
    "PF24948": ["tca_cycle", "central_metabolism"],  # Citrate_synth_N

    # RNA ligases
    "PF09511": ["ligase", "rna_processing"],  # RNA_lig_T4_1 - RNA ligase

    # Glycosyltransferases
    "PF21969": ["glycosyltransferase", "carbohydrate_active"],  # MGS_GT

    # Phosphatases
    "PF04029": ["phosphatase", "hydrolase"],  # 2-ph_phosp - 2-phosphoglycerate phosphatase

    # Acetyl-CoA synthetases
    "PF06094": ["ligase", "lipid_metabolism"],  # GGACT - Acetyl-CoA synthetase-like

    # iPGM (phosphoglycerate mutase)
    "PF06415": ["isomerase", "glycolysis"],  # iPGM_N - independent phosphoglycerate mutase

    # Glutamine synthetase
    "PF18318": ["ligase", "nitrogen_metabolism", "amino_acid_biosynthesis"],  # Gln-synt_C-ter

    # RHS repeats
    "PF05593": ["repeat_domain", "toxin"],  # RHS_repeat

    # Thiamin biosynthesis
    "PF13379": ["cofactor_biosynthesis", "thiamine_biosynthesis"],  # NMT1_2

    # Glutamine synthetase
    "PF12437": ["ligase", "nitrogen_metabolism"],  # GSIII_N - Glutamine synthetase III N-terminal

    # STAND NTPases
    "PF20703": ["atp_binding", "signaling", "defense_system"],  # nSTAND1

    # ABC transporters
    "PF04069": ["transporter", "abc_transporter", "abc_substrate_binding"],  # OpuAC - Glycine betaine binding

    # IMS domains
    "PF11799": ["dna_repair"],  # IMS_C

    # Glycosyltransferases
    "PF13704": ["glycosyltransferase", "transferase", "carbohydrate_active"],  # Glyco_tranf_2_4

    # Nucleases
    "PF13366": ["nuclease", "hydrolase"],  # PDDEXK_3

    # Poly(A) polymerase
    "PF01743": ["nucleotidyltransferase", "rna_processing"],  # PolyA_pol

    # Diguanylate cyclases
    "PF08859": ["signaling", "cyclic_dinucleotide", "diguanylate_cyclase"],  # DGC - diguanylate cyclase

    # Nucleases
    "PF01927": ["nuclease", "rna_processing"],  # Mut7-C

    # DNA repair
    "PF04423": ["dna_repair", "zinc_binding"],  # Rad50_zn_hook

    # Carbohydrate binding
    "PF06739": ["carbohydrate_binding"],  # SBBP

    # Adaptin
    "PF01602": ["membrane", "vesicle_trafficking"],  # Adaptin_N

    # CRISPR
    "PF17953": ["crispr_associated", "defense_system"],  # Csm4_C

    # DNA excision
    "PF03013": ["nuclease", "mobile_element"],  # Pyr_excise - Pyrimidine excision

    # Molybdenum cofactor
    "PF02391": ["cofactor_biosynthesis", "molybdenum_binding"],  # MoaE

    # EHD binding
    "PF10622": ["membrane", "signaling"],  # Ehbp

    # Cobalamin biosynthesis
    "PF01955": ["cofactor_biosynthesis", "cobalamin_biosynthesis"],  # CbiZ

    # Porphobilinogen
    "PF01379": ["lyase", "heme_biosynthesis"],  # Porphobil_deam - Porphobilinogen deaminase

    # Kelch repeats
    "PF07646": ["repeat_domain", "kelch_repeat"],  # Kelch_2

    # Formyl transferase
    "PF02911": ["transferase", "one_carbon_metabolism"],  # Formyl_trans_C

    # DrsE family
    "PF13686": ["membrane"],  # DrsE_2

    # Helicase
    "PF08706": ["helicase", "replication"],  # D5_N

    # Restriction enzymes
    "PF11463": ["restriction_enzyme", "restriction_modification"],  # R-HINP1I

    # Transposases
    "PF14319": ["transposase", "mobile_element"],  # Zn_Tnp_IS91

    # CRISPR
    "PF18211": ["crispr_associated", "defense_system"],  # Csm1_B

    # Cell surface
    "PF20251": ["cell_surface", "repeat_domain"],  # Big_14

    # Alcohol dehydrogenases
    "PF13602": ["oxidoreductase", "dehydrogenase", "zinc_binding"],  # ADH_zinc_N_2

    # DNA binding
    "PF17723": ["dna_binding", "regulator"],  # RHH_8 - Ribbon-helix-helix

    # Vutriplase
    "PF02592": ["protease", "hydrolase"],  # Vut_1

    # Restriction modification
    "PF12008": ["restriction_modification", "restriction_enzyme"],  # EcoR124_C

    # HMG-CoA synthase
    "PF01154": ["transferase", "lipid_metabolism", "isoprenoid_biosynthesis"],  # HMG_CoA_synt_N

    # Methylation
    "PF18755": ["methyltransferase", "transferase"],  # RAMA

    # -------------------------------------------------------------------------
    # BATCH 12: MORE UNMAPPED DOMAINS (counts 4+)
    # -------------------------------------------------------------------------
    # Methyltransferases
    "PF04989": ["methyltransferase", "transferase", "carbohydrate_active"],  # RMNT_CmcI - Rhamnosyl O-methyltransferase

    # Iron-sulfur proteins
    "PF00301": ["iron_sulfur", "metal_binding", "electron_transport"],  # Rubredoxin

    # E3 ligase
    "PF12483": ["ligase", "protein_modification"],  # GIDE - E3 Ubiquitin ligase

    # Photosynthesis
    "PF18933": ["photosynthesis"],  # PsbP_2 - PsbP-like protein
    "PF01789": ["photosynthesis"],  # PsbP - oxygen-evolving complex

    # Proteases
    "PF01327": ["protease", "hydrolase", "metal_binding"],  # Pep_deformylase - Polypeptide deformylase
    "PF22148": ["protease", "regulatory"],  # Fervidolysin_NPro-like - N-terminal prodomain

    # Secretion systems
    "PF16576": ["secretion_component", "t1ss_component", "membrane"],  # HlyD_D23 - membrane fusion protein
    "PF13437": ["secretion_component", "t1ss_component", "membrane"],  # HlyD_3 - HlyD family secretion protein
    "PF18962": ["secretion_component", "cell_surface"],  # Por_Secre_tail - Secretion C-terminal sorting

    # Pathogenicity/virulence
    "PF11731": ["virulence", "hypothetical"],  # Cdd1 - Pathogenicity locus

    # Ligases/transferases
    "PF25084": ["translation_factor", "translation"],  # LbH_EIF2B - EIF2B subunit epsilon
    "PF00432": ["transferase", "isoprenoid_biosynthesis"],  # Prenyltrans - Prenyltransferase

    # Antitoxins
    "PF13274": ["antitoxin_domain", "defense_system"],  # SocA_Panacea - Antitoxin SocA-like

    # Phosphatases/hydrolases
    "PF10110": ["hydrolase", "phosphodiesterase", "membrane"],  # GPDPase_memb - glycerophosphoryl diester phosphodiesterase
    "PF04608": ["phosphatase", "hydrolase", "lipid_metabolism"],  # PgpA - Phosphatidylglycerophosphatase A
    "PF19580": ["nuclease", "phosphatase", "hydrolase"],  # Exo_endo_phos_3

    # Aromatic metabolism
    "PF01947": ["lyase", "aromatic_aa_metabolism"],  # Rv2949c-like - Chorismate pyruvate-lyase
    "PF02567": ["lyase", "secondary_metabolism"],  # PhzC-PhzF - Phenazine biosynthesis

    # Cell surface
    "PF11867": ["restriction_enzyme", "restriction_modification"],  # T1RH-like_C - Type I restriction enzyme
    "PF00431": ["cell_surface", "repeat_domain"],  # CUB domain

    # Hypothetical
    "PF21827": ["hypothetical"],  # New_glue protein family
    "PF23343": ["replication", "mobile_element"],  # REP_ORF2-G2P - Replication-associated protein

    # ABC transporters
    "PF14524": ["transporter", "abc_transporter", "abc_permease"],  # Wzt_C - Wzt C-terminal domain
    "PF12399": ["transporter", "abc_transporter", "abc_atpase"],  # BCA_ABC_TP_C - Branched-chain amino acid ABC transporter

    # TPR repeats
    "PF07720": ["repeat_domain", "tpr_repeat", "protein_binding"],  # TPR_3

    # Phage related
    "PF04865": ["phage_related", "phage_baseplate"],  # Baseplate_J - Baseplate J-like protein
    "PF17289": ["phage_related", "phage_terminase", "nuclease"],  # Terminase_6C - Terminase RNaseH-like domain

    # Cofactor biosynthesis
    "PF01227": ["hydrolase", "cofactor_biosynthesis", "folate_biosynthesis"],  # GTP_cyclohydroI - GTP cyclohydrolase I
    "PF01983": ["transferase", "cofactor_biosynthesis"],  # CofC - Guanylyl transferase
    "PF02548": ["transferase", "cofactor_biosynthesis"],  # Pantoate_transf - Ketopantoate hydroxymethyltransferase
    "PF08546": ["oxidoreductase", "reductase", "cofactor_biosynthesis"],  # ApbA_C - Ketopantoate reductase

    # Cell surface/binding
    "PF13205": ["cell_surface", "repeat_domain"],  # Big_5 - Bacterial Ig-like domain
    "PF22544": ["cell_surface", "repeat_domain"],  # HYDIN_VesB_CFA65-like_Ig
    "PF07678": ["binding", "cell_surface"],  # TED_complement - A-macroglobulin TED domain
    "PF07691": ["carbohydrate_binding", "cell_surface"],  # PA14 domain

    # Signaling
    "PF17203": ["signaling", "sensor_kinase"],  # sCache_3_2 - Single cache domain
    "PF17202": ["signaling", "sensor_kinase"],  # sCache_3_3 - Single cache domain

    # DNA repair
    "PF21999": ["dna_polymerase", "dna_repair"],  # IMS_HHH_1 - DNA polymerase-iota thumb domain
    "PF11798": ["dna_binding", "dna_repair"],  # IMS_HHH motif
    "PF10107": ["nuclease", "dna_repair", "recombinational_repair"],  # Endonuc_Holl - Holliday junction resolvase

    # TIR domain
    "PF10137": ["defense_system", "signaling"],  # CAP12-PCTIR_TIR - TIR domain

    # Transporters
    "PF02659": ["transporter", "metal_transporter", "manganese_binding"],  # Mntp - Manganese efflux pump
    "PF11449": ["transporter", "heavy_metal_resistance", "membrane"],  # ArsP_2 - heavy-metal exporter

    # Membrane proteins
    "PF18917": ["membrane", "stress_response"],  # LiaI-LiaF-like_TM1

    # Translation
    "PF02686": ["trna_synthetase", "translation"],  # GatC - Glu-tRNAGln amidotransferase C subunit

    # Dehydrogenases
    "PF03971": ["oxidoreductase", "dehydrogenase", "tca_cycle"],  # IDH - Monomeric isocitrate dehydrogenase

    # Kelch repeats
    "PF24681": ["repeat_domain", "kelch_repeat"],  # Kelch_KLHDC2_KLHL20_DRC7

    # Iron-sulfur assembly
    "PF01883": ["iron_sulfur_biosynthesis", "iron_sulfur"],  # FeS_assembly_P

    # Chemotaxis
    "PF03705": ["chemotaxis", "methyltransferase"],  # CheR_N - CheR all-alpha domain

    # Repeat domains
    "PF01535": ["repeat_domain", "regulatory"],  # PPR repeat

    # Transposases
    "PF13751": ["transposase", "mobile_element", "nuclease"],  # DDE_Tnp_1_6
    "PF13737": ["transposase", "mobile_element", "nuclease"],  # DDE_Tnp_1_5

    # CRISPR
    "PF23400": ["crispr_associated", "defense_system"],  # CARF_Card1 - Card1 CARF domain

    # Heme biosynthesis
    "PF01208": ["lyase", "decarboxylase", "heme_biosynthesis"],  # URO-D - Uroporphyrinogen decarboxylase

    # Toxin-antitoxin
    "PF18478": ["toxin_domain", "pin_domain"],  # PIN_10

    # Cell surface
    "PF14686": ["carbohydrate_active", "polysaccharide_lyase"],  # fn3_3 - Polysaccharide lyase family 4

    # Regulatory
    "PF00325": ["regulator", "transcription_factor", "crp_fnr_family", "dna_binding"],  # Crp - bacterial regulatory proteins

    # Cysteine-rich
    "PF19114": ["hypothetical"],  # EsV_1_7_cys - cysteine-rich motif

    # Sulfur metabolism
    "PF02635": ["sulfurtransferase", "transferase", "sulfur_metabolism"],  # DsrE - DsrE/DsrF-like family

    # Restriction enzymes
    "PF09571": ["restriction_enzyme", "restriction_modification"],  # RE_XcyI

    # Amino acid biosynthesis
    "PF14821": ["lyase", "amino_acid_biosynthesis"],  # Thr_synth_N - Threonine synthase N terminus

    # Gluconolactonase
    "PF08450": ["hydrolase"],  # SGL - SMP-30/Gluconolactonase/LRE-like

    # Sugar isomerases
    "PF01182": ["isomerase", "carbohydrate_active"],  # Glucosamine_iso - Glucosamine-6-phosphate isomerase

    # Various
    "PF17942": ["hypothetical"],  # Morc6_S5 - ribosomal protein S5 domain-like
    "PF10544": ["hypothetical"],  # T5orf172 domain
    "PF10592": ["defense_system"],  # AIPR protein
    "PF02536": ["regulator", "transcription_termination"],  # mTERF - mitochondrial transcription termination factor
    "PF00436": ["dna_binding", "replication"],  # SSB - Single-strand binding protein
    "PF13784": ["transferase", "signaling"],  # Fic_N - Fic/DOC family N-terminal
    "PF14445": ["zinc_binding", "protein_modification"],  # Prok-RING_2 - Prokaryotic RING finger
    "PF00376": ["regulator", "transcription_factor", "dna_binding", "helix_turn_helix"],  # MerR family
    "PF17284": ["transferase", "amino_acid_metabolism"],  # Spermine_synt_N - Spermidine synthase
    "PF02666": ["lyase", "decarboxylase", "lipid_metabolism"],  # PS_Dcarbxylase - Phosphatidylserine decarboxylase
    "PF03601": ["hypothetical"],  # Cons_hypoth698
    "PF00496": ["transporter", "abc_transporter", "abc_substrate_binding"],  # SBP_bac_5 - solute-binding protein
    "PF03460": ["iron_sulfur", "denitrification", "sulfate_reduction"],  # NIR_SIR_ferr - Nitrite/Sulfite reductase ferredoxin

    # -------------------------------------------------------------------------
    # BATCH 13: REMAINING UNMAPPED DOMAINS (counts 3-4)
    # -------------------------------------------------------------------------
    # Restriction enzymes
    "PF14511": ["restriction_enzyme", "restriction_modification"],  # RE_EcoO109I
    "PF09233": ["restriction_enzyme", "restriction_modification"],  # Endonuc-EcoRV

    # One-carbon metabolism
    "PF01812": ["ligase", "one_carbon_metabolism", "atp_binding"],  # 5-FTHF_cyc-lig - 5-formyltetrahydrofolate cyclo-ligase

    # DNA binding
    "PF02178": ["dna_binding"],  # AT_hook

    # NTPases
    "PF12643": ["hydrolase", "nucleotide_metabolism"],  # MazG-like

    # Flagella
    "PF18998": ["flagellum", "flagellin"],  # Flg_new_2

    # Secretion
    "PF08700": ["secretion_component", "membrane"],  # VPS51_Exo84_N

    # Deaminases
    "PF14438": ["hydrolase", "deaminase"],  # SM-ATX

    # Archaeal domains
    "PF17791": ["hypothetical"],  # MG3

    # Toxins
    "PF14487": ["toxin", "transferase"],  # DarT - ADP-ribosyltransferase toxin

    # Chemotaxis
    "PF22673": ["chemotaxis", "signaling"],  # MCP-like_PDC_1

    # Phosphatases
    "PF02562": ["atpase", "atp_binding", "stress_response"],  # PhoH

    # DNA repair
    "PF01149": ["dna_repair", "base_excision_repair"],  # Fapy_DNA_glyco - Formamidopyrimidine-DNA glycosylase

    # Sulfotransferases
    "PF03567": ["sulfurtransferase", "transferase"],  # Sulfotransfer_2

    # Transcription
    "PF02291": ["transcription", "regulator"],  # TFIID-31kDa

    # Glycosyltransferases
    "PF01762": ["glycosyltransferase", "transferase", "carbohydrate_active"],  # Galactosyl_T

    # Membrane proteins
    "PF06813": ["membrane", "transporter"],  # Nodulin-like

    # Hypothetical
    "PF23262": ["hypothetical"],  # NFD4_C

    # Nucleases
    "PF00773": ["rnase", "hydrolase", "rna_processing"],  # RNB - Ribonuclease II/R

    # Transporters
    "PF00532": ["transporter", "abc_transporter", "abc_substrate_binding"],  # Peripla_BP_1

    # Binding domains
    "PF00207": ["binding", "protease"],  # A2M - Alpha-2-macroglobulin family

    # Glutamate metabolism
    "PF05201": ["oxidoreductase", "amino_acid_metabolism"],  # GlutR_N - Glutamyl-tRNA reductase

    # Gas vesicles
    "PF05800": ["membrane", "cell_surface"],  # GvpO - Gas vesicle protein
    "PF05120": ["membrane", "cell_surface"],  # GvpG - Gas vesicle protein

    # Glycoside hydrolases
    "PF14885": ["glycosidase", "hydrolase", "carbohydrate_active"],  # GHL15

    # Heparan sulfate
    "PF07940": ["hydrolase", "carbohydrate_active"],  # Hepar_II_III_C

    # Ligases
    "PF16575": ["kinase", "rna_processing"],  # CLP1_P

    # Kinases
    "PF00485": ["kinase", "carbohydrate_active"],  # PRK - Phosphoribulokinase

    # Potassium transport
    "PF02669": ["transporter", "ion_transporter", "membrane"],  # KdpC - Potassium-transporting ATPase

    # Isoprenoid biosynthesis
    "PF01222": ["oxidoreductase", "reductase", "isoprenoid_biosynthesis"],  # ERG4_ERG24

    # Hypothetical
    "PF22784": ["hypothetical"],  # PTP-SAK

    # Carotenoid biosynthesis
    "PF18916": ["lyase", "isoprenoid_biosynthesis"],  # Lycopene_cyc - Lycopene cyclase

    # Repeat domains
    "PF07661": ["repeat_domain", "membrane"],  # MORN_2

    # CRISPR
    "PF14526": ["crispr_associated", "defense_system"],  # Cass2

    # Carbohydrate binding
    "PF16841": ["carbohydrate_binding", "carbohydrate_active"],  # CBM60

    # Conjugation
    "PF11130": ["conjugation", "mobile_element"],  # TraC_F_IV

    # Nucleotide metabolism
    "PF05014": ["transferase", "nucleotide_metabolism"],  # Nuc_deoxyrib_tr - Nucleoside deoxyribosyltransferase

    # Sorting
    "PF15902": ["membrane", "vesicle_trafficking"],  # Sortilin-Vps10

    # tRNA editing
    "PF04073": ["hydrolase", "trna_modification"],  # tRNA_edit

    # DNA binding
    "PF06831": ["dna_binding", "dna_repair"],  # H2TH - Helix-2-turn-helix

    # Glycosyltransferases
    "PF03808": ["glycosyltransferase", "transferase", "carbohydrate_active"],  # Glyco_tran_WecG

    # Hydrolases
    "PF06439": ["hydrolase", "carbohydrate_active"],  # 3keto-disac_hyd

    # Phage
    "PF23899": ["phage_related", "phage_portal"],  # SU10_portal

    # Glyoxalase
    "PF12681": ["hydrolase", "oxidative_stress"],  # Glyoxalase_2

    # Bioluminescence
    "PF05893": ["oxidoreductase", "bioluminescence"],  # LuxC

    # Sugar metabolism
    "PF03214": ["glycosyltransferase", "carbohydrate_active"],  # RGP - Reversibly glycosylated polypeptide

    # Membrane proteins
    "PF01184": ["transporter", "membrane"],  # Gpr1_Fun34_YaaH

    # NAD metabolism
    "PF02233": ["transporter", "membrane"],  # PNTB - NAD(P) transhydrogenase

    # Bleomycin resistance
    "PF22677": ["hypothetical"],  # Ble-like_N

    # Methyltransferases
    "PF05711": ["methyltransferase", "transferase", "secondary_metabolism"],  # TylF

    # Archaeal
    "PF01835": ["hypothetical"],  # MG2

    # RNase
    "PF11977": ["rnase", "hydrolase"],  # RNase_Zc3h12a

    # Beta propeller
    "PF24981": ["repeat_domain", "wd40_repeat"],  # Beta-prop_ATRN-LZTR1

    # FAE
    "PF20434": ["esterase", "hydrolase", "carbohydrate_active"],  # BD-FAE - feruloyl esterase

    # Isoprenoid
    "PF01128": ["transferase", "isoprenoid_biosynthesis"],  # IspD - 2-C-methyl-D-erythritol 4-phosphate cytidylyltransferase

    # Hypothetical
    "PF22785": ["hypothetical"],  # Tc-R-P

    # Nucleotide metabolism
    "PF15891": ["transferase", "nucleotide_metabolism"],  # Nuc_deoxyri_tr2

    # A2M binding
    "PF07703": ["binding"],  # A2M_BRD - A-macroglobulin bait region domain

    # Iron uptake
    "PF04773": ["regulator", "iron_binding"],  # FecR - FecR protein

    # Esterases
    "PF05448": ["esterase", "hydrolase", "carbohydrate_active"],  # AXE1 - Acetyl xylan esterase

    # Restriction enzymes
    "PF19778": ["restriction_enzyme", "restriction_modification"],  # RE_endonuc

    # ATPases
    "PF21138": ["helicase", "atp_binding"],  # SMUBP-2_HCS1_1B

    # Protease inhibitors
    "PF00014": ["protease", "regulator"],  # Kunitz_BPTI - Kunitz/BPTI serine protease inhibitor

    # HD domain
    "PF24391": ["hydrolase", "phosphodiesterase"],  # HD-CE

    # Acetyltransferases
    "PF10686": ["acetyltransferase", "transferase"],  # YAcAr

    # Photosynthesis
    "PF14870": ["photosynthesis"],  # PSII_BNR - Photosystem II BNR

    # Pyrrolysine
    "PF21360": ["ligase", "amino_acid_biosynthesis"],  # PylC-like_N - Pyrrolysine biosynthesis

    # AIG2 family
    "PF13772": ["hydrolase"],  # AIG2_2

    # HPPK
    "PF01288": ["kinase", "cofactor_biosynthesis", "folate_biosynthesis"],  # HPPK - 2-amino-4-hydroxy-6-hydroxymethyldihydropteridine pyrophosphokinase

    # ICE elements
    "PF05315": ["mobile_element"],  # ICEA - Integrative conjugative element protein

    # tRNA modification
    "PF08351": ["acetyltransferase", "trna_modification"],  # TmcA_N

    # Carbohydrate binding
    "PF08305": ["carbohydrate_binding", "carbohydrate_active"],  # NPCBM - N-terminal domain of pectin methylesterase inhibitors

    # Decarboxylases
    "PF01276": ["lyase", "decarboxylase"],  # OKR_DC_1 - Orn/Lys/Arg decarboxylase

    # Heme biosynthesis
    "PF02602": ["transferase", "heme_biosynthesis"],  # HEM4 - Uroporphyrinogen-III synthase

    # Glycoside hydrolases
    "PF02449": ["glycosidase", "hydrolase", "carbohydrate_active"],  # Glyco_hydro_42 - Beta-galactosidase

    # ECH family
    "PF16113": ["lyase", "fatty_acid_degradation"],  # ECH_2 - Enoyl-CoA hydratase

    # Glycoside hydrolases
    "PF02435": ["glycosidase", "hydrolase", "carbohydrate_active"],  # Glyco_hydro_68 - Levansucrase

    # Poly(A) polymerase
    "PF12627": ["rna_binding", "nucleotidyltransferase"],  # PolyA_pol_RNAbd

    # Lyases
    "PF03328": ["lyase", "central_metabolism"],  # HpcH_HpaI - HpcH/HpaI aldolase

    # AAD
    "PF18785": ["hypothetical"],  # Inv-AAD

    # Primase
    "PF08708": ["primase", "replication"],  # PriCT_1

    # MZB
    "PF09369": ["hypothetical"],  # MZB

    # RNases
    "PF14622": ["rnase", "hydrolase"],  # Ribonucleas_3_3

    # DNA binding
    "PF12651": ["dna_binding", "regulator", "ribbon_helix_helix"],  # RHH_3

    # DNA polymerase
    "PF00712": ["dna_polymerase", "replication"],  # DNA_pol3_beta - DNA polymerase III beta subunit

    # Membrane
    "PF22570": ["membrane", "stress_response"],  # LiaF-TM

    # Restriction enzymes
    "PF12183": ["restriction_enzyme", "restriction_modification"],  # NotI

    # Hypothetical
    "PF21347": ["hypothetical"],  # DUF3108_like

    # Cell division
    "PF21980": ["cell_division"],  # MksE

    # Glycoside hydrolases
    "PF12899": ["glycosidase", "hydrolase", "carbohydrate_active"],  # Glyco_hydro_100

    # LPS biosynthesis
    "PF02706": ["membrane", "lps_biosynthesis", "regulatory"],  # Wzz - Chain length determinant protein

    # Archaeal flagella
    "PF04659": ["flagellum", "cell_surface"],  # Arch_fla_DE - Archaeal flagella assembly protein

    # Exonucleases
    "PF09810": ["nuclease", "hydrolase"],  # Exo5 - Exonuclease V

    # Transposases
    "PF13701": ["transposase", "mobile_element", "nuclease"],  # DDE_Tnp_1_4

    # Sulfatase
    "PF03781": ["hydrolase", "sulfur_metabolism"],  # FGE-sulfatase - Formylglycine-generating enzyme

    # FtsK/SpoIIIE
    "PF01580": ["cell_division", "dna_binding", "atp_binding"],  # FtsK_SpoIIIE - DNA translocase

    # DMP
    "PF14300": ["hypothetical"],  # DMP19

    # Condensin
    "PF12717": ["cell_division", "chromosome_partitioning"],  # Cnd1

    # Sporulation
    "PF09579": ["cell_division"],  # Spore_YtfJ

    # Viral
    "PF01693": ["phage_related", "hypothetical"],  # Cauli_VI

    # GIIM
    "PF08388": ["membrane", "signaling"],  # GIIM

    # Chemotaxis
    "PF05228": ["signaling", "chemotaxis"],  # CHASE4

    # Arginine deiminase
    "PF02274": ["hydrolase", "amino_acid_degradation"],  # ADI - Arginine deiminase

    # Zinc finger
    "PF06467": ["zinc_binding", "protein_modification"],  # zf-FCS

    # Potassium transport
    "PF02702": ["sensor_kinase", "two_component", "ion_transporter"],  # KdpD

    # PFO
    "PF12367": ["hydrolase", "toxin"],  # PFO_beta_C - Perfringolysin O

    # GST
    "PF13409": ["transferase", "oxidative_stress"],  # GST_N_2 - Glutathione S-transferase N-terminal

    # YhfC
    "PF10086": ["membrane"],  # YhfC

    # Archaeal
    "PF24271": ["hypothetical"],  # HVO_2833_C

    # MEDS
    "PF14417": ["metal_binding"],  # MEDS - Metal ion-dependent adhesion site

    # -------------------------------------------------------------------------
    # BATCH 14: FINAL UNMAPPED DOMAINS TO GET BELOW 500
    # -------------------------------------------------------------------------
    # Glycoside hydrolases
    "PF00251": ["glycosidase", "hydrolase", "carbohydrate_active"],  # Glyco_hydro_32N - GH32 N-terminal

    # Chaperones
    "PF18050": ["isomerase", "chaperone"],  # Cyclophil_like2 - Cyclophilin-like
    "PF00166": ["chaperone", "hsp60", "stress_response"],  # Cpn10 - Chaperonin 10 Kd subunit

    # Iron-sulfur proteins
    "PF06902": ["iron_sulfur", "ferredoxin", "metal_binding", "electron_transport"],  # Fer4_19 - Divergent 4Fe-4S

    # Signaling domains
    "PF17201": ["signaling", "sensor_kinase"],  # Cache_3-Cache_2 - Cache fusion domain
    "PF08269": ["signaling", "sensor_kinase"],  # dCache_2 - Cache domain

    # Flavin binding
    "PF07311": ["cofactor_binding", "fad_binding"],  # Dodecin - flavin binding protein

    # Transposases
    "PF13610": ["transposase", "mobile_element", "nuclease"],  # DDE_Tnp_IS240
    "PF04693": ["transposase", "mobile_element"],  # DDE_Tnp_2 - Archaeal putative transposase

    # Pilus biogenesis
    "PF23989": ["pilus", "type_iv_pilus", "atpase"],  # PilB3_C - PilB3-like C-terminal

    # DNA repair/helicase
    "PF09111": ["dna_binding", "helicase", "replication"],  # SLIDE domain
    "PF14490": ["helicase", "dna_repair", "dna_binding"],  # HHH_RecD2 - ATP-dependent RecD2 DNA helicase
    "PF22982": ["helicase", "atp_binding", "dna_binding"],  # WHD_HRQ1 - helicase HRQ1 winged helix

    # Electron transport
    "PF03116": ["electron_transport", "iron_sulfur", "respiration"],  # NQR2_RnfD_RnfE - NADH:quinone oxidoreductase

    # Cytochromes (heme-containing electron carriers)
    "PF13435": ["cytochrome", "electron_transport", "respiration", "heme_binding"],  # Cytochrome_C554
    "PF13442": ["cytochrome", "electron_transport", "respiration", "heme_binding"],  # Cytochrome_CBB3
    "PF14522": ["cytochrome", "electron_transport", "respiration", "heme_binding"],  # Cytochrome_C7
    "PF22085": ["cytochrome", "electron_transport", "respiration", "heme_binding"],  # NorB_cytochrome_c-like
    "PF00033": ["cytochrome", "electron_transport", "respiration", "heme_binding"],  # Cytochrome_B
    "PF11783": ["cytochrome", "electron_transport", "respiration", "heme_binding"],  # Cytochrome_cB

    # Translation
    "PF00472": ["translation_factor", "translation"],  # RF-1 domain - release factor

    # Cell surface
    "PF09134": ["cell_surface", "adhesin"],  # Invasin_D3 - Invasin domain 3

    # Methyltransferase
    "PF22458": ["methyltransferase", "iron_sulfur", "rrna_modification"],  # RsmF-B_ferredox

    # Secondary metabolism
    "PF10604": ["secondary_metabolism", "polyketide_synthesis"],  # Polyketide_cyc2
    "PF03364": ["secondary_metabolism", "polyketide_synthesis", "lipid_metabolism"],  # Polyketide_cyc

    # Virulence
    "PF03534": ["virulence", "toxin"],  # SpvB - Salmonella virulence

    # Phosphatases/transferases
    "PF03767": ["phosphatase", "hydrolase"],  # Acid_phosphat_B - HAD superfamily
    "PF01515": ["transferase", "central_metabolism"],  # PTA_PTB - Phosphate acetyl/butaryl transferase

    # DNA polymerase
    "PF02767": ["dna_polymerase", "replication", "sliding_clamp"],  # DNA_pol3_beta_2 - DNA polymerase III beta central domain

    # Membrane proteins
    "PF09586": ["membrane"],  # YfhO - Bacterial membrane protein
    "PF10639": ["membrane", "hypothetical"],  # TMEM234 - Putative transmembrane family

    # Transporters
    "PF00230": ["transporter", "membrane", "ion_channel"],  # MIP - Major intrinsic protein (aquaporin)
    "PF20293": ["transporter", "abc_transporter"],  # MC6 - ABC-3C Middle Component 6

    # Cell surface
    "PF16158": ["cell_surface", "repeat_domain"],  # N_BRCA1_IG - Ig-like domain

    # Metal binding
    "PF15616": ["metal_binding", "heavy_metal_resistance"],  # TerY_C - TerY-C metal binding

    # Sigma factor
    "PF14532": ["regulator", "sigma_factor", "atp_binding"],  # Sigma54_activ_2 - Sigma-54 interaction domain

    # -------------------------------------------------------------------------
    # NAME-ONLY HMMs (pipelines without PF accessions)
    # -------------------------------------------------------------------------
    "NosD": ["copper_binding", "metal_binding"],  # Periplasmic copper-binding protein
    "Beta_helix": ["beta_helix"],  # Right-handed beta helix region
    "TFIIB": ["transcription", "dna_binding"],  # Transcription factor IIB
    "TBP": ["transcription", "dna_binding", "regulator"],  # TATA-binding protein
    "HHH_5": ["dna_binding"],  # Helix-hairpin-helix domain
    "baeRF_family10": ["translation_factor", "translation"],  # Release factor family
    "acVLRF1": ["translation_factor", "translation"],  # Release factor family
    "YdjM": ["membrane"],  # Inner membrane protein
    "VWA_2": ["protein_binding"],  # von Willebrand factor type A domain
    "SPASM": ["radical_sam", "iron_sulfur"],  # Iron-sulfur cluster-binding domain
    "SMC_N": ["chromosome_partitioning"],  # SMC N-terminal domain
    "RimK": ["ligase", "protein_modification"],  # ATP-grasp ligase
    "PTH2": ["hydrolase", "translation"],  # Peptidyl-tRNA hydrolase
    "PLDc_2": ["lipase", "hydrolase"],  # PLD-like domain
    "MYG1_exonuc": ["nuclease", "rna_processing"],  # 3'-5' exonuclease
    "Lipocalin_2": ["beta_barrel"],  # Lipocalin-like domain
    "Lactamase_B_3": ["hydrolase", "metal_binding"],  # Metallo-beta-lactamase superfamily
    "Lactamase_B_2": ["hydrolase", "metal_binding"],  # Metallo-beta-lactamase superfamily
    "Lactamase_B": ["hydrolase", "metal_binding"],  # Metallo-beta-lactamase superfamily
    "HD": ["phosphatase", "hydrolase", "metal_binding"],  # HD phosphohydrolase
    "GmrSD_N": ["restriction_enzyme", "restriction_modification"],  # GmrSD restriction endonuclease
    "GmrSD_C": ["restriction_enzyme", "restriction_modification"],  # GmrSD restriction endonuclease
    "Glyco_tranf_2_3": ["glycosyltransferase", "transferase", "carbohydrate_active"],  # Glycosyltransferase family 2
    "GTP_EFTU": ["gtpase", "gtp_binding", "translation_factor", "translation"],  # EF-Tu GTP-binding
    "GSH-S_ATP": ["ligase"],  # Glutathione synthetase (ATP-grasp)
    "GIY-YIG": ["nuclease", "hydrolase"],  # GIY-YIG nuclease domain
    "GATase": ["glutamine_amidotransferase", "transferase"],  # Glutamine amidotransferase class I
    "FtsJ": ["methyltransferase", "rrna_modification"],  # rRNA methyltransferase
    "Flg_new_2": ["repeat_domain"],  # InlB B-repeat domain
    "FeThRed_B": ["oxidoreductase", "iron_sulfur", "thioredoxin"],  # Ferredoxin thioredoxin reductase
    "ERCC4": ["nuclease", "dna_repair"],  # ERCC4 nuclease domain
    "DevR": ["crispr_associated", "defense_system"],  # CRISPR-associated regulator
    "DSPc": ["phosphatase", "hydrolase"],  # Dual specificity phosphatase
    "DNA_processg_A": ["dna_repair", "recombinational_repair"],  # DprA recombination mediator
    "DNA_pol_B_exo1": ["dna_polymerase", "nuclease", "replication"],  # DNA polymerase B exonuclease
    "DNA_pol_B": ["dna_polymerase", "replication"],  # DNA polymerase family B
    "DEAD_2": ["rna_binding", "helicase", "atp_binding"],  # DEAD-box helicase
    "CxxCxxCC": ["zinc_binding", "metal_binding"],  # Metal-chelating motif
    "Csc2": ["crispr_associated", "cas_domain", "defense_system"],  # Cascade subunit Csc2
    "Cas_csx3": ["crispr_associated", "cas_domain", "defense_system"],  # Cas/Csx3
    "CDP-OH_P_transf": ["transferase", "lipid_metabolism", "phospholipid_metabolism"],  # CDP-alcohol phosphatidyltransferase
    "CBM96": ["carbohydrate_active", "carbohydrate_binding"],  # Carbohydrate-binding module 96
    "NikA-like": ["conjugation", "mobile_element"],  # Mobilization protein NikA
    "CBFD_NFYB_HMF": ["histone", "chromatin", "dna_binding"],  # Histone-like TF / archaeal histone
    "CAP12-PCTIR_TIR": ["defense_system"],  # TIR-domain anti-phage effector
    "Bromo_TP": ["dna_binding", "transcription"],  # Bromodomain-associated, DNA-binding
    "BRCT": ["dna_repair", "protein_binding"],  # BRCT phosphopeptide-binding domain
    "AAA_lid_8": ["aaa_domain", "atp_binding"],  # AAA lid domain
    "AAA_5": ["aaa_domain", "atp_binding"],  # AAA domain (dynein-related)
    "AAA_21": ["aaa_domain", "atp_binding"],  # AAA domain
    "AAA_19": ["aaa_domain", "atp_binding", "atpase"],  # AAA domain
    "AAA_15": ["aaa_domain", "atp_binding", "atpase"],  # AAA ATPase domain
    "Roc": ["gtpase", "gtp_binding"],  # Ras-like GTPase (Roc)
    "Ras": ["gtpase", "gtp_binding", "signaling"],  # Ras family GTPase
    "Arf": ["gtpase", "gtp_binding", "signaling"],  # ADP-ribosylation factor
    "Gtr1_RagA": ["gtpase", "gtp_binding"],  # RagA-like GTPase
    "MMR_HSR1": ["gtpase", "gtp_binding", "ribosomal_protein", "translation"],  # Ribosome-binding GTPase
    "DEAD": ["dna_repair", "helicase", "atp_binding"],  # DEAD/DEAH box helicase
    "AAA_lid_3": ["aaa_domain", "atp_binding"],  # AAA+ lid domain
    "AAA_23": ["aaa_domain", "atp_binding", "atpase"],  # AAA domain
    "RCF1-5-like_lid": ["aaa_domain", "atp_binding", "atpase"],  # AAA+ ATPase lid domain
    "Eco57I": ["restriction_enzyme", "restriction_modification", "nuclease"],  # Eco57I restriction enzyme
    "ResIII": ["restriction_modification", "restriction_enzyme", "nuclease"],  # Type III restriction enzyme
    "TaqI_C": ["restriction_modification", "restriction_enzyme"],  # TaqI-like specificity domain
    "XPG_N": ["nuclease", "dna_repair"],  # XPG N-terminal domain (NER)
    "XPG_I": ["nuclease", "dna_repair"],  # XPG I-region (NER)
    "UvrD_C": ["helicase", "dna_repair"],  # UvrD-like helicase C-terminal
    "UvrD_C_2": ["helicase", "dna_repair"],  # UvrD-like helicase C-terminal
    "SNF2-rel_dom": ["helicase", "atp_binding", "regulator"],  # SNF2-related helicase
    "MCM_lid": ["helicase", "replication"],  # MCM AAA-lid domain
    "MCM_OB": ["helicase", "replication", "rna_binding"],  # MCM OB domain
    "PCNA_N": ["replication", "sliding_clamp"],  # Sliding clamp (PCNA) N-terminal domain
    "Rep_fac_C": ["replication", "clamp_loader"],  # Replication factor C C-terminal domain
    "Rep_fac-A_C": ["replication", "clamp_loader"],  # Replication factor A C-terminal domain
    "Rad17": ["replication", "clamp_loader"],  # Rad17 clamp loader subunit
    "DNA_pol3_delta2": ["replication", "clamp_loader"],  # DNA polymerase III delta subunit
    "RNR_Alpha": ["oxidoreductase", "reductase", "nucleotide_metabolism"],  # Ribonucleotide reductase alpha
    "RNR-II_ins_dom": ["oxidoreductase", "reductase", "nucleotide_metabolism"],  # RNR insertion domain
    "RuvC": ["nuclease", "dna_repair", "recombinational_repair"],  # Holliday junction resolvase
    "UDG": ["dna_repair", "base_excision_repair"],  # Uracil-DNA glycosylase
    "HNH": ["nuclease", "hydrolase"],  # HNH endonuclease
    "HNH_2": ["nuclease", "hydrolase"],  # HNH endonuclease
    "HNH_5": ["nuclease", "hydrolase"],  # HNH endonuclease
    "NERD": ["nuclease", "hydrolase"],  # Nuclease-related domain
    "SWIM": ["zinc_binding"],  # SWIM zinc finger
    "Metallophos": ["phosphatase", "hydrolase", "metal_binding"],  # Metallophosphatase
    "Metallophos_2": ["phosphatase", "hydrolase", "metal_binding"],  # Metallophosphatase
    "dUTPase": ["hydrolase", "nucleotide_metabolism"],  # dUTP diphosphatase
    "dCMP_cyt_deam_1": ["deaminase", "hydrolase", "nucleotide_metabolism"],  # Cytidine deaminase
    "Thymidylate_kin": ["kinase", "transferase", "nucleotide_metabolism"],  # Thymidylate kinase
    "NUDIX": ["hydrolase", "nucleotide_metabolism"],  # Nudix hydrolase
    "MazG": ["hydrolase", "nucleotide_metabolism"],  # MazG pyrophosphohydrolase
    "MazG-like": ["hydrolase", "nucleotide_metabolism"],  # MazG-like hydrolase
    "NT5C": ["phosphatase", "nucleotide_metabolism"],  # 5' nucleotidase
    "tRNA_NucTran2_2": ["nucleotidyltransferase", "trna_modification"],  # tRNA nucleotidyltransferase
    "PolyA_pol": ["polymerase"],  # Poly(A) polymerase
    "PadR": ["regulator", "transcription_factor", "dna_binding"],  # PadR transcriptional regulator
    "Response_reg": ["regulator", "response_regulator", "two_component"],  # Response regulator receiver domain
    "Y1_Tnp": ["transposase", "mobile_element", "nuclease"],  # IS200-like transposase
    "Prok-JAB": ["deubiquitinase", "protease"],  # JAB deubiquitinase-like
    "JAB": ["deubiquitinase", "protease"],  # JAB deubiquitinase-like
    "PTP-SAK": ["phosphatase", "hydrolase"],  # Protein tyrosine phosphatase
    "Tc-R-P": ["phosphatase", "hydrolase"],  # DSP-PTPase phosphatase
    "PAP2": ["phosphatase", "hydrolase"],  # PAP2 family phosphatase
    "PAP2_3": ["phosphatase", "hydrolase"],  # PAP2 family phosphatase
    "ubiquitin": ["ubiquitin_like"],  # Ubiquitin family
    "ThiF": ["ubiquitin_activation", "protein_modification"],  # E1-like ubiquitin activation
    "Histone": ["histone"],  # Histone proteins
}


def _load_external_pfam_mappings() -> dict[str, list[str]]:
    """Load PFAM predicate mappings from a TSV file (data-driven extension)."""
    candidates = [
        Path("data/reference/pfam_predicate_map.tsv"),
        Path(__file__).parent.parent.parent.parent / "data/reference/pfam_predicate_map.tsv",
    ]
    path = next((p for p in candidates if p.exists()), None)
    if path is None:
        return {}

    mappings: dict[str, list[str]] = {}
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            key = parts[0].strip()
            preds = [p.strip() for p in parts[1].split(",") if p.strip()]
            if not key or not preds:
                continue
            mappings[key] = preds

    return mappings


def _merge_external_pfam_mappings(base: dict[str, list[str]]) -> dict[str, list[str]]:
    """Merge external mappings into the in-code PFAM mapping table."""
    merged = {k: list(v) for k, v in base.items()}
    external = _load_external_pfam_mappings()
    for key, preds in external.items():
        if key in merged:
            merged[key] = sorted(set(merged[key]) | set(preds))
        else:
            merged[key] = preds
    return merged


PFAM_TO_PREDICATES = _merge_external_pfam_mappings(PFAM_TO_PREDICATES)


# ============================================================================
# PATTERN-BASED MAPPINGS
# ============================================================================
# Format: (pattern_regex, [predicates])
# Applied to PFAM name or description when no direct mapping exists

PFAM_PATTERNS: list[tuple[str, list[str]]] = [
    # Transporters
    (r"\btransport\b", ["transporter"]),
    (r"ABC.*(transporter|permease|ATPase)", ["transporter", "abc_transporter"]),
    (r"MFS|Major.facilitat", ["transporter", "mfs_transporter"]),
    (r"\bchannel\b", ["transporter", "ion_channel"]),
    (r"\bporin\b", ["transporter", "porin", "outer_membrane"]),
    (r"\befflux\b", ["transporter", "efflux_pump"]),
    (r"\bsymport", ["transporter", "symporter"]),
    (r"\bantiport", ["transporter", "antiporter"]),
    (r"\bpermease\b", ["transporter", "membrane"]),
    (r"FeoB|ferrous.*iron.*transport", ["transporter", "metal_transporter", "iron_binding"]),

    # Enzymes
    (r"dehydrogenase", ["oxidoreductase", "dehydrogenase"]),
    (r"reductase", ["oxidoreductase", "reductase"]),
    (r"oxidase", ["oxidoreductase", "oxidase"]),
    (r"oxygenase", ["oxidoreductase", "oxygenase"]),
    (r"peroxidase", ["oxidoreductase", "peroxidase"]),
    (r"kinase", ["kinase", "transferase"]),
    (r"phosphatase", ["phosphatase", "hydrolase"]),
    (r"protease|peptidase", ["protease", "hydrolase"]),
    (r"nuclease|endonuclease|exonuclease", ["nuclease", "hydrolase"]),
    (r"lipase", ["lipase", "hydrolase"]),
    (r"esterase", ["esterase", "hydrolase"]),
    (r"hydrolase", ["hydrolase"]),
    (r"transferase", ["transferase"]),
    (r"lyase", ["lyase"]),
    (r"isomerase", ["isomerase"]),
    (r"ligase|synthetase", ["ligase"]),
    (r"synthase", ["synthase"]),
    (r"methyltransferase|methylase", ["methyltransferase", "transferase"]),
    (r"acetyltransferase", ["acetyltransferase", "transferase"]),
    (r"glycosyltransferase", ["glycosyltransferase", "transferase"]),

    # Regulators
    (r"transcription.*regulator|regulator.*transcription", ["regulator", "transcription_factor"]),
    (r"DNA.binding.*transcription|transcription.*factor", ["regulator", "transcription_factor", "dna_binding"]),
    (r"helix.turn.helix|HTH", ["dna_binding", "helix_turn_helix"]),
    (r"response.regulator", ["regulator", "response_regulator", "two_component"]),
    (r"sensor.*kinase|histidine.*kinase", ["sensor_kinase", "two_component"]),
    (r"LysR", ["regulator", "lysr_family"]),
    (r"TetR", ["regulator", "tetr_family"]),
    (r"GntR", ["regulator", "gntr_family"]),
    (r"AraC", ["regulator", "arac_family"]),
    (r"LuxR", ["regulator", "luxr_family"]),
    (r"MarR", ["regulator", "marr_family"]),
    (r"LacI", ["regulator", "laci_family"]),
    (r"sigma.*factor", ["sigma_factor", "regulator"]),

    # Signaling
    (r"GGDEF", ["signaling", "cyclic_dinucleotide", "diguanylate_cyclase"]),
    (r"EAL|c.di.GMP", ["signaling", "cyclic_dinucleotide"]),
    (r"PAS.*domain", ["signaling"]),
    (r"GAF.*domain", ["signaling"]),

    # Binding
    (r"DNA.bind|binds.*DNA", ["dna_binding"]),
    (r"RNA.bind|binds.*RNA", ["rna_binding"]),
    (r"ATP.bind|binds.*ATP|ATPase", ["atp_binding"]),
    (r"GTP.bind|binds.*GTP|GTPase", ["gtp_binding"]),
    (r"NAD.bind|NAD\+|NADH|NADP", ["nad_binding", "cofactor_binding"]),
    (r"FAD.bind|FAD", ["fad_binding", "cofactor_binding"]),
    (r"FMN", ["fmn_binding", "cofactor_binding"]),
    (r"PLP|pyridoxal", ["plp_binding", "cofactor_binding"]),
    (r"heme|haem|cytochrome", ["heme_binding", "iron_binding", "cytochrome"]),
    (r"iron.sulfur|Fe.S|4Fe.4S|2Fe.2S|ferredoxin", ["iron_sulfur", "metal_binding"]),
    (r"zinc.*finger|Zn.*finger", ["zinc_finger", "zinc_binding", "metal_binding"]),
    (r"nickel|urease|NiFe|Ni-Fe", ["nickel_binding", "metal_binding"]),
    (r"molybdo|MoCo|molybdopterin|tungsten", ["molybdenum_binding", "metal_binding"]),
    (r"copper|cupredoxin|plastocyanin|laccase|Cu-", ["copper_binding", "metal_binding"]),
    (r"cobalt|cobalamin|B12|adenosylcob|methylcob", ["cobalt_binding", "cobalamin_binding", "metal_binding"]),
    (r"manganese|Mn-|arginase", ["manganese_binding", "metal_binding"]),
    (r"calcium|EF.hand|calmodulin", ["calcium_binding", "metal_binding"]),
    (r"zinc.*bind|Zn.*bind|zinc.*site", ["zinc_binding", "metal_binding"]),
    (r"metal.*bind", ["metal_binding"]),

    # Membrane
    (r"membrane|transmembrane", ["membrane"]),
    (r"outer.*membrane", ["outer_membrane", "membrane"]),
    (r"signal.*peptide|secretion.*signal", ["secreted"]),
    (r"periplasmic", ["periplasmic"]),

    # Cell structures
    (r"flagell", ["flagellum"]),
    (r"pil[iu]|fimbr", ["pilus"]),
    (r"adhesin", ["adhesin", "cell_surface"]),
    (r"peptidoglycan|murein", ["peptidoglycan", "cell_wall"]),
    (r"capsul", ["capsule"]),
    (r"LPS|lipopolysaccharide|lipid.*A", ["lps_biosynthesis"]),

    # Mobile elements
    (r"transposase", ["transposase", "mobile_element"]),
    (r"integrase", ["integrase", "mobile_element"]),
    (r"recombinase", ["recombinase", "mobile_element"]),
    (r"resolvase", ["resolvase", "mobile_element"]),
    (r"phage|bacteriophage", ["phage_related"]),

    # Stress/defense
    (r"heat.*shock|HSP\d+", ["heat_shock", "stress_response", "chaperone"]),
    (r"cold.*shock|CSP", ["cold_shock", "stress_response"]),
    (r"chaperone|chaperonin", ["chaperone", "stress_response"]),
    (r"catalase", ["catalase", "oxidative_stress"]),
    (r"superoxide.*dismutase|SOD", ["superoxide_dismutase", "oxidative_stress"]),
    (r"thioredoxin", ["thioredoxin", "oxidative_stress"]),
    (r"glutaredoxin", ["glutaredoxin", "oxidative_stress"]),
    (r"beta.lactamase", ["beta_lactamase", "antibiotic_resistance"]),
    (r"antibiotic.*resistance", ["antibiotic_resistance"]),
    (r"multidrug|MDR", ["multidrug_resistance", "antibiotic_resistance"]),
    (r"arsenic|ars[ABC]", ["arsenic_resistance", "heavy_metal_resistance"]),
    (r"mercury|mer[ABCDE]", ["mercury_resistance", "heavy_metal_resistance"]),

    # Metabolism
    # NOTE: \b word boundary prevents "dehydrogenase" from matching "hydrogenase"
    (r"\bhydrogenase\b", ["hydrogenase", "hydrogen_metabolism"]),
    # NOTE: Be VERY careful with nitrogen fixation patterns. The NifH/frxC domain family includes
    # both nitrogenase reductase (NifH) AND ferredoxin:plastoquinone reductase (frxC, photosynthesis).
    # Also, molybdopterin domains have descriptions like "Nitrogenase component 1 type Oxidoreductase"
    # but are NOT nitrogenase. Only match "nitrogenase reductase" or "nitrogenase subunit" specifically.
    (r"nitrogenase\s+(reductase|subunit|NifH)", ["nitrogenase", "nitrogen_fixation"]),
    (r"nitrate.*reductase|nar[GHI]", ["nitrate_reduction", "denitrification"]),
    (r"nitrite.*reductase|nir[KS]", ["denitrification"]),
    # NOTE: "methano" pattern should NOT assign methanogenesis - many tetrahydromethanopterin
    # enzymes are found in non-methanogens. Only MCR subunits are definitive.
    (r"methano(?!genesis)", ["archaeal_one_carbon", "one_carbon_metabolism"]),  # tetrahydromethanopterin enzymes
    (r"methyl.*coenzyme.*M.*reductase|MCR[_-]?[ABG]", ["methanogenesis", "mcr_complex"]),  # ONLY MCR is definitive
    # NOTE: RuBisCO alone does NOT mean Calvin cycle - RuBisCO-like proteins (RLPs) exist.
    # Only assign rubisco and carbon_fixation, NOT calvin_cycle (requires PRK confirmation).
    (r"RuBisCO|ribulose.*bisphosphate.*carboxylase", ["rubisco", "carbon_fixation"]),
    (r"photosystem", ["photosynthesis"]),
    (r"glycolysis|glycolytic", ["glycolysis", "central_metabolism"]),
    (r"TCA|citric.*acid|Krebs", ["tca_cycle", "central_metabolism"]),
    (r"pentose.*phosphate", ["pentose_phosphate", "central_metabolism"]),
    (r"fatty.*acid.*synth", ["fatty_acid_synthesis", "lipid_metabolism"]),
    (r"beta.*oxidation|fatty.*acid.*degrad", ["fatty_acid_degradation", "lipid_metabolism"]),

    # CAZy-like patterns
    (r"cellulase", ["cellulase", "glycoside_hydrolase", "carbohydrate_active"]),
    (r"chitinase", ["chitinase", "glycoside_hydrolase", "carbohydrate_active"]),
    (r"amylase", ["amylase", "glycoside_hydrolase", "carbohydrate_active"]),
    (r"xylanase", ["xylanase", "glycoside_hydrolase", "carbohydrate_active"]),
    (r"pectinase|pectin.*lyase", ["pectinase", "carbohydrate_active"]),

    # Structural
    (r"TPR.*repeat", ["repeat_domain", "tpr_repeat"]),
    (r"WD.?40|beta.*propeller", ["repeat_domain", "wd40_repeat"]),
    (r"LRR|leucine.rich.*repeat", ["repeat_domain", "lrr_repeat"]),
    (r"ankyrin", ["repeat_domain", "ankyrin_repeat"]),
    (r"coiled.coil", ["coiled_coil"]),

    # Defense systems
    (r"restriction.*enzyme|restriction.*endonuclease", ["restriction_enzyme", "restriction_modification"]),
    # NOTE: CRISPR patterns should assign cas_domain (not crispr_associated) since many Cas-like
    # domains are shared with transposases. Only direct mappings for specific Cas1/Cas2/etc. get crispr_associated.
    (r"CRISPR|Cas\d+", ["cas_domain", "defense_system"]),
    # NOTE: toxin/antitoxin patterns should use domain-level predicates, not system-level
    (r"toxin.antitoxin|TA.*system", ["toxin_domain", "antitoxin_domain", "defense_system"]),

    # Information processing
    (r"DNA.*polymerase|pol[AB]|dna[EGI]", ["dna_polymerase", "replication"]),
    (r"helicase", ["helicase"]),
    (r"topoisomerase|gyr[AB]", ["topoisomerase", "replication"]),
    (r"primase", ["primase", "replication"]),
    (r"RNA.*polymerase|rpo[ABC]", ["rna_polymerase", "transcription"]),
    (r"ribosom", ["ribosomal_protein", "translation"]),
    (r"tRNA.*synthetase|aminoacyl.*tRNA|tRNA.synt", ["trna_synthetase", "translation"]),
    (r"elongation.*factor|EF.Tu|EF.G", ["translation_factor", "translation"]),

    # AAA+ ATPases
    (r"\bAAA\b", ["aaa_domain", "atp_binding", "atpase"]),
    (r"ATPase", ["atpase", "atp_binding"]),

    # Epimerases
    (r"\bepimerase\b", ["epimerase", "isomerase"]),
    (r"dehydratase|dTDP.*dehydratase", ["dehydratase", "lyase"]),

    # HAD superfamily
    (r"haloacid|HAD", ["hydrolase", "phosphatase"]),

    # Glycosyltransferases (broader patterns)
    (r"Glycos.*transf|Glyco.*trans", ["glycosyltransferase", "transferase", "carbohydrate_active"]),
    (r"polysaccharide.*synth", ["glycosyltransferase", "carbohydrate_active"]),
    (r"GT\d+|glycosyl.*transfer", ["glycosyltransferase", "transferase"]),

    # Aminotransferases (broader patterns)
    (r"Aminotran|aminotransfer", ["aminotransferase", "transferase", "plp_binding"]),

    # Ferredoxins/iron-sulfur (broader patterns)
    (r"\bFer4\b|\bFer2\b|ferredoxin", ["iron_sulfur", "ferredoxin", "metal_binding", "electron_transport"]),

    # TPR repeats
    (r"\bTPR\b|tetratricopeptide", ["repeat_domain", "tpr_repeat", "protein_binding"]),

    # HEAT repeats
    (r"\bHEAT\b", ["repeat_domain", "heat_repeat", "protein_binding"]),

    # Methyltransferases (broader)
    (r"Methyltransf|Mtase|MTase", ["methyltransferase", "transferase", "sam_binding"]),

    # Acetyltransferases
    (r"acetyltransf|GNAT", ["acetyltransferase", "transferase"]),

    # Radical SAM enzymes
    (r"Radical.*SAM|radical.SAM", ["radical_sam", "iron_sulfur", "enzyme"]),

    # PQQ-dependent enzymes
    (r"\bPQQ\b", ["oxidoreductase", "pqq_binding"]),

    # Cytochrome/heme binding
    (r"cytochrome|Cytochrom", ["cytochrome", "heme_binding", "electron_transport", "respiration"]),

    # Secretion systems
    (r"T2SS|Type.II.*secret", ["secretion_component", "t2ss_component"]),
    (r"T3SS|Type.III.*secret", ["secretion_component", "t3ss_component"]),
    (r"T4SS|Type.IV.*secret", ["secretion_component", "t4ss_component"]),
    (r"T6SS|Type.VI.*secret", ["secretion_component", "t6ss_component"]),
    (r"\bSec[ABDEFGY]\b|SecD|SecF", ["secretion_component", "sec_pathway"]),
    (r"\bTat[ABC]\b|twin.arginine", ["secretion_component", "tat_pathway"]),

    # Biosynthesis patterns
    (r"biosynth", ["biosynthesis"]),
    (r"NAD.*synth", ["nad_biosynthesis", "cofactor_biosynthesis"]),
    (r"thiamin|ThiI|TPP", ["thiamine_biosynthesis", "cofactor_biosynthesis"]),
    (r"cobalamin|B12|Cob[A-Z]", ["cobalamin_biosynthesis", "cofactor_biosynthesis"]),

    # Cell division
    (r"\bFts[AZWLNQKX]\b|divisome|division", ["cell_division", "divisome"]),
    (r"\bPar[AB]\b|partition", ["cell_division", "chromosome_partitioning"]),
    (r"\bSMC\b|structural.*maintenance", ["cell_division", "chromosome_partitioning"]),

    # Toxin-antitoxin (broader patterns)
    (r"\bantitoxin\b|MazE|HicB|RelB|ParD|VapB", ["antitoxin_domain", "defense_system"]),
    (r"\btoxin\b|MazF|HicA|RelE|ParE|VapC", ["toxin_domain", "defense_system"]),
    (r"\bPIN\b.*domain", ["toxin_domain", "pin_domain"]),

    # GTPases
    (r"\bGTPase\b|GTP.*binding|G.protein", ["gtpase", "gtp_binding"]),
    (r"50S.*ribosome.*GTPase|ribosome.*GTPase", ["gtpase", "ribosomal_protein", "translation"]),

    # Nitroreductases
    (r"nitroreductase", ["oxidoreductase", "reductase"]),

    # Flavoproteins
    (r"flavoprotein|flavin", ["oxidoreductase", "fad_binding"]),

    # NOTE: CRISPR pattern defined earlier in defense systems section

    # Amidohydrolases
    (r"amidohydrolase|Amidohydro", ["hydrolase", "amidase", "metal_binding"]),

    # Carbamoyl phosphate synthase
    (r"carbamoyl.*phosphate|CPSase", ["ligase", "amino_acid_metabolism"]),

    # Magnesium chelatase
    (r"chelatase", ["chelatase", "atp_binding"]),

    # F420 coenzyme - found in archaea AND Actinobacteria, NOT methanogen-specific
    (r"F420", ["oxidoreductase", "f420_dependent", "archaeal_one_carbon"]),

    # Acylphosphatase
    (r"acylphosphatase", ["hydrolase", "phosphatase"]),

    # ParA/ParB families
    (r"ParA.*NTPase|MinD", ["cell_division", "atp_binding"]),

    # Hemerythrin (oxygen binding)
    (r"hemerythrin|Hemerythrin", ["oxygen_binding", "iron_binding"]),

    # -------------------------------------------------------------------------
    # ADDITIONAL PATTERNS FOR COMPREHENSIVE COVERAGE
    # -------------------------------------------------------------------------
    # Intein and splicing
    (r"[Ii]ntein|splicing", ["mobile_element"]),

    # HEPN domain - RNA processing/defense
    (r"\bHEPN\b", ["nuclease", "defense_system"]),

    # KOW/KH domains - RNA binding
    (r"\bKOW\b|\bKH_\d", ["rna_binding", "translation"]),

    # ABC transporter related
    (r"MacB|ABC.*permease", ["transporter", "abc_transporter"]),

    # Adenylsuccinate synthetase
    (r"Adenyls|adenylsuccin", ["ligase", "purine_metabolism"]),

    # Topoisomerase/primase
    (r"Toprim|topoisom", ["topoisomerase", "replication"]),

    # ATP-grasp fold
    (r"ATP.grasp", ["ligase", "atp_binding"]),

    # tRNA modification
    (r"tRNA.*trans|tRNA_Me|tRNA.*modif", ["trna_modification", "translation"]),

    # Zinc ribbon transcription factors
    (r"Zn_Ribbon|zinc.*ribbon", ["zinc_binding", "dna_binding", "transcription"]),

    # PAC domain
    (r"\bPAC\b", ["signaling"]),

    # PrmA - ribosomal protein methyltransferase
    (r"\bPrmA\b", ["methyltransferase", "translation"]),

    # PUA domain - RNA binding
    (r"\bPUA\b", ["rna_binding"]),

    # Rad51/RecA family
    (r"Rad51|RecA", ["dna_repair", "recombinational_repair", "atp_binding"]),

    # UDPG dehydrogenase family
    (r"UDPG.*dh|UDP.*dehydrogen", ["oxidoreductase", "dehydrogenase", "carbohydrate_active"]),

    # MCM - DNA replication
    (r"\bMCM\b", ["helicase", "replication", "atp_binding"]),

    # S1 domain - RNA binding
    (r"\bS1\b.*domain|Ribosomal.*S1", ["rna_binding", "translation"]),

    # TrkA - potassium transport regulatory
    (r"\bTrkA\b|\bTrk[HK]\b", ["transporter", "ion_transporter", "regulatory"]),

    # D-Ala-D-Ala ligase
    (r"Dala.*lig|D-Ala", ["ligase", "peptidoglycan", "cell_wall"]),

    # AbiEii - abortive infection
    (r"\bAbi[A-Z]+\b", ["defense_system", "abortive_infection"]),

    # NRDD - anaerobic ribonucleotide reductase
    (r"\bNRDD\b", ["oxidoreductase", "reductase", "nucleotide_metabolism"]),

    # SRP54 - signal recognition particle
    (r"\bSRP\d*\b", ["signal_recognition", "gtp_binding", "secretion_component"]),

    # LeuA - leucine biosynthesis
    (r"\bLeuA\b", ["amino_acid_biosynthesis", "transferase"]),

    # Isocitrate dehydrogenase
    (r"Iso_dh|isocitrate.*dehydrogen", ["oxidoreductase", "dehydrogenase", "tca_cycle"]),

    # Proteasome
    (r"[Pp]roteasome", ["protease", "hydrolase"]),

    # AdoHcyase - SAM metabolism
    (r"AdoHcy|S-adenosyl.*homocysteine", ["hydrolase", "sam_binding"]),

    # HhH-GPD - DNA repair glycosylase
    (r"HhH.*GPD|glycosylase", ["dna_repair", "base_excision_repair"]),

    # UbiA - ubiquinone biosynthesis
    (r"\bUbiA\b|ubiquinone", ["transferase", "cofactor_biosynthesis"]),

    # Gate domain
    (r"\bGate\b", ["membrane"]),

    # LigT/phosphoesterase
    (r"LigT|PEase|phosphoesterase", ["hydrolase", "phosphatase"]),

    # bpMoxR - AAA ATPase
    (r"MoxR|bpMoxR", ["aaa_domain", "atp_binding", "chaperone"]),

    # SNAP - soluble NSF attachment protein
    (r"\bSNAP\b", ["membrane", "vesicle_trafficking"]),

    # ATP synthase subunits
    (r"ATP.synt|ATPsynthase", ["atp_synthesis", "energy_metabolism"]),

    # Semialdehyde dehydrogenase
    (r"[Ss]emialdhyde.*dh|semialdehyde.*dehydrogen", ["oxidoreductase", "dehydrogenase", "amino_acid_metabolism"]),

    # RNase PH
    (r"RNase.*PH|Rnase_PH", ["rnase", "hydrolase", "rna_processing"]),

    # Release factor
    (r"\beRF\d|release.*factor", ["translation_factor", "translation"]),

    # RrnaAD - rRNA adenine methylase
    (r"RrnaAD|rRNA.*methyl", ["methyltransferase", "rrna_modification"]),

    # RNA polymerase
    (r"RNA.*pol|Rpb\d", ["rna_polymerase", "transcription"]),

    # GatB - Glu-tRNA amidotransferase
    (r"\bGatB\b|Glu.*tRNA.*amidotrans", ["trna_synthetase", "translation"]),

    # IF-2/translation initiation
    (r"\bIF-\d|initiation.*factor", ["translation_factor", "translation"]),

    # TGS domain
    (r"\bTGS\b", ["regulatory", "gtp_binding"]),

    # FKBP - peptidyl-prolyl isomerase
    (r"\bFKBP\b", ["isomerase", "chaperone"]),

    # MetW - methionine biosynthesis
    (r"\bMet[WS]\b", ["amino_acid_biosynthesis"]),

    # Sigma54 activator
    (r"Sigma54.*activat", ["regulator", "sigma_factor", "atp_binding"]),

    # Aconitase
    (r"[Aa]conitase", ["tca_cycle", "iron_sulfur", "lyase"]),

    # PAS domain variants
    (r"\bPAS_\d", ["signaling"]),

    # ACT domain variants
    (r"\bACT_\d", ["regulatory", "amino_acid_metabolism"]),

    # PKD domain
    (r"\bPKD\b", ["cell_surface", "repeat_domain"]),

    # Transglutaminase
    (r"[Tt]ransglut", ["transferase"]),

    # Phosphomutase
    (r"[Pp]hospho[Mm]utase|PGM_PMM", ["isomerase", "mutase", "carbohydrate_active"]),

    # AIRC - purine biosynthesis
    (r"\bAIRC\b", ["lyase", "purine_metabolism"]),

    # Fic domain - AMPylation
    (r"\bFic\b", ["transferase"]),

    # LAGLIDADG - homing endonuclease
    (r"LAGLIDADG", ["nuclease", "mobile_element"]),

    # Rhomboid protease
    (r"[Rr]homboid", ["protease", "membrane"]),

    # GIDA - tRNA modification
    (r"\bGIDA\b", ["oxidoreductase", "trna_modification"]),

    # ThiG - thiamine biosynthesis
    (r"\bThiG\b", ["thiamine_biosynthesis", "cofactor_biosynthesis"]),

    # DHH domain - phosphoesterase
    (r"\bDHH\b", ["phosphatase", "hydrolase"]),

    # ATP-cone - allosteric ATP binding
    (r"ATP.cone", ["regulatory", "atp_binding"]),

    # Proton antiporter
    (r"[Pp]roton.*antipo", ["transporter", "antiporter", "membrane"]),

    # FeoA - iron transport
    (r"\bFeoA\b", ["transporter", "iron_binding"]),

    # Rubrerythrin - oxidative stress
    (r"[Rr]ubrerythrin", ["oxidative_stress", "iron_binding"]),

    # MGS - metal-binding
    (r"\bMGS\b", ["metal_binding"]),

    # OrfB/IS605 - insertion sequence
    (r"OrfB.*IS|IS\d+.*transpos", ["transposase", "mobile_element"]),

    # RIO kinase
    (r"\bRIO\d?\b", ["kinase", "translation"]),

    # GGR - glycosyl hydrolase
    (r"\bGGR\b", ["hydrolase", "carbohydrate_active"]),

    # S4 domain - RNA binding
    (r"\bS4\b.*domain", ["rna_binding"]),

    # TraB/conjugation
    (r"\bTraB\b|conjugat", ["conjugation", "mobile_element"]),

    # HEAT PBS - phycobilisome
    (r"HEAT.*PBS|phycobili", ["repeat_domain", "photosynthesis"]),

    # DUF patterns - mark as hypothetical but some have hints
    (r"\bDUF\d+\b", ["hypothetical"]),

    # General oxidoreductase patterns
    (r"oxidored|Oxidored", ["oxidoreductase"]),

    # General binding patterns
    (r"_bind\b|binding", ["binding"]),

    # Synthase/synthetase patterns
    (r"_synt\b|synth\b", ["synthase"]),

    # Hydrolase patterns
    (r"_hydro\b|hydrol", ["hydrolase"]),

    # Ligase patterns
    (r"_lig\b", ["ligase"]),

    # Domain with C-terminal/N-terminal
    (r"_[NC]\b", []),  # Skip these as they're domain fragments
]


def get_predicates_for_pfam(accession: str, name: str = "", description: str = "") -> list[str]:
    """
    Get predicates for a PFAM domain.

    First checks direct mappings, then applies pattern matching.

    Args:
        accession: PFAM accession (e.g., "PF00005")
        name: PFAM short name (e.g., "ABC_tran")
        description: PFAM description

    Returns:
        List of predicate IDs
    """
    predicates = set()

    # Direct mapping
    if accession in PFAM_TO_PREDICATES:
        predicates.update(PFAM_TO_PREDICATES[accession])

    # Pattern matching on name and description
    text = f"{name} {description}".lower()
    for pattern, preds in PFAM_PATTERNS:
        if re.search(pattern, text, re.IGNORECASE):
            predicates.update(preds)

    return sorted(predicates)


__all__ = [
    "PFAM_TO_PREDICATES",
    "PFAM_PATTERNS",
    "get_predicates_for_pfam",
]
