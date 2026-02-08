"""
Predicate vocabulary - comprehensive definitions of all predicates.

Predicates are organized by category:
- enzyme: Enzyme classes and subtypes (EC-based)
- transport: Transport and localization
- regulation: Regulatory and signaling functions
- metabolism: Metabolic pathway associations
- cazy: Carbohydrate-active enzyme categories
- structure: Structural/domain features
- binding: Binding and interaction domains
- envelope: Cell envelope and surface proteins
- mobile: Mobile genetic elements and defense
- stress: Stress response and resistance
- size: Size-based predicates
- annotation: Annotation status predicates
- composition: Sequence composition predicates
- info_processing: DNA/RNA/protein synthesis machinery
- division: Cell division and chromosome partitioning
- topology: Transmembrane topology predictions (requires pyTMHMM)
- viral: Viral protein functions (from VOGdb)
"""

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class PredicateVocab:
    """Vocabulary entry for a predicate."""
    predicate_id: str
    name: str
    description: str
    category: str
    parent: Optional[str] = None  # For hierarchical predicates


# ============================================================================
# ENZYME CLASSES (EC-based)
# ============================================================================
ENZYME_PREDICATES = [
    # Top-level EC classes
    PredicateVocab("oxidoreductase", "Oxidoreductase", "EC 1.x.x.x - catalyzes oxidation-reduction reactions", "enzyme"),
    PredicateVocab("transferase", "Transferase", "EC 2.x.x.x - transfers functional groups", "enzyme"),
    PredicateVocab("hydrolase", "Hydrolase", "EC 3.x.x.x - cleaves bonds using water", "enzyme"),
    PredicateVocab("lyase", "Lyase", "EC 4.x.x.x - cleaves bonds without hydrolysis or oxidation", "enzyme"),
    PredicateVocab("isomerase", "Isomerase", "EC 5.x.x.x - catalyzes isomerization", "enzyme"),
    PredicateVocab("ligase", "Ligase", "EC 6.x.x.x - joins molecules using ATP", "enzyme"),
    PredicateVocab("translocase", "Translocase", "EC 7.x.x.x - catalyzes movement across membranes", "enzyme"),

    # Oxidoreductase subtypes
    PredicateVocab("dehydrogenase", "Dehydrogenase", "NAD/NADP-dependent dehydrogenase", "enzyme", "oxidoreductase"),
    PredicateVocab("oxidase", "Oxidase", "Oxygen-dependent oxidase", "enzyme", "oxidoreductase"),
    PredicateVocab("reductase", "Reductase", "Reduction catalyst", "enzyme", "oxidoreductase"),
    PredicateVocab("peroxidase", "Peroxidase", "Peroxide-dependent oxidoreductase", "enzyme", "oxidoreductase"),
    PredicateVocab("oxygenase", "Oxygenase", "Oxygen incorporation enzyme", "enzyme", "oxidoreductase"),
    PredicateVocab("monooxygenase", "Monooxygenase", "Single oxygen atom incorporation", "enzyme", "oxidoreductase"),
    PredicateVocab("dioxygenase", "Dioxygenase", "Both oxygen atoms incorporated", "enzyme", "oxidoreductase"),

    # Transferase subtypes
    PredicateVocab("kinase", "Kinase", "Phosphotransferase (adds phosphate)", "enzyme", "transferase"),
    PredicateVocab("methyltransferase", "Methyltransferase", "Methyl group transferase", "enzyme", "transferase"),
    PredicateVocab("acetyltransferase", "Acetyltransferase", "Acetyl group transferase", "enzyme", "transferase"),
    PredicateVocab("glycosyltransferase", "Glycosyltransferase", "Sugar residue transferase", "enzyme", "transferase"),
    PredicateVocab("aminotransferase", "Aminotransferase", "Amino group transferase (transaminase)", "enzyme", "transferase"),
    PredicateVocab("phosphatase", "Phosphatase", "Phosphate-removing enzyme", "enzyme", "hydrolase"),

    # Hydrolase subtypes
    PredicateVocab("protease", "Protease", "Peptide bond hydrolase", "enzyme", "hydrolase"),
    PredicateVocab("peptidase", "Peptidase", "Peptide hydrolase", "enzyme", "hydrolase"),
    PredicateVocab("nuclease", "Nuclease", "Nucleic acid hydrolase", "enzyme", "hydrolase"),
    PredicateVocab("dnase", "DNase", "DNA hydrolase", "enzyme", "nuclease"),
    PredicateVocab("rnase", "RNase", "RNA hydrolase", "enzyme", "nuclease"),
    PredicateVocab("lipase", "Lipase", "Lipid ester hydrolase", "enzyme", "hydrolase"),
    PredicateVocab("esterase", "Esterase", "Ester bond hydrolase", "enzyme", "hydrolase"),
    PredicateVocab("glycosidase", "Glycosidase", "Glycosidic bond hydrolase", "enzyme", "hydrolase"),
    PredicateVocab("amidase", "Amidase", "Amide bond hydrolase", "enzyme", "hydrolase"),
    PredicateVocab("atpase", "ATPase", "ATP hydrolase", "enzyme", "hydrolase"),
    PredicateVocab("gtpase", "GTPase", "GTP hydrolase", "enzyme", "hydrolase"),

    # Ligase subtypes
    PredicateVocab("synthase", "Synthase", "Bond-forming enzyme (no ATP)", "enzyme"),
    PredicateVocab("synthetase", "Synthetase", "Bond-forming enzyme (ATP-dependent)", "enzyme", "ligase"),
    PredicateVocab("carboxylase", "Carboxylase", "CO2-adding enzyme", "enzyme", "ligase"),
]

# ============================================================================
# TRANSPORT & LOCALIZATION
# ============================================================================
TRANSPORT_PREDICATES = [
    # General transport
    PredicateVocab("transporter", "Transporter", "Any transport function", "transport"),
    PredicateVocab("membrane", "Membrane protein", "Membrane-associated protein", "transport"),
    PredicateVocab("transmembrane", "Transmembrane", "Spans the membrane", "transport", "membrane"),

    # ABC transporters
    PredicateVocab("abc_transporter", "ABC transporter", "ATP-binding cassette transporter", "transport", "transporter"),
    PredicateVocab("abc_atpase", "ABC ATPase", "ABC transporter ATP-binding component", "transport", "abc_transporter"),
    PredicateVocab("abc_permease", "ABC permease", "ABC transporter membrane component", "transport", "abc_transporter"),
    PredicateVocab("abc_substrate_binding", "ABC substrate binding", "ABC transporter substrate-binding component", "transport", "abc_transporter"),

    # Other transporter families
    PredicateVocab("mfs_transporter", "MFS transporter", "Major facilitator superfamily transporter", "transport", "transporter"),
    PredicateVocab("ion_channel", "Ion channel", "Ion channel protein", "transport", "transporter"),
    PredicateVocab("porin", "Porin", "Outer membrane porin", "transport", "transporter"),
    PredicateVocab("efflux_pump", "Efflux pump", "Drug/metabolite efflux system", "transport", "transporter"),
    PredicateVocab("symporter", "Symporter", "Co-transporter (same direction)", "transport", "transporter"),
    PredicateVocab("antiporter", "Antiporter", "Exchange transporter (opposite direction)", "transport", "transporter"),
    PredicateVocab("uniporter", "Uniporter", "Single-substrate transporter", "transport", "transporter"),

    # Substrate-specific transporters
    PredicateVocab("sugar_transporter", "Sugar transporter", "Carbohydrate transport", "transport", "transporter"),
    PredicateVocab("amino_acid_transporter", "Amino acid transporter", "Amino acid transport", "transport", "transporter"),
    PredicateVocab("peptide_transporter", "Peptide transporter", "Peptide/oligopeptide transport", "transport", "transporter"),
    PredicateVocab("ion_transporter", "Ion transporter", "Inorganic ion transport", "transport", "transporter"),
    PredicateVocab("metal_transporter", "Metal transporter", "Metal ion transport", "transport", "transporter"),
    PredicateVocab("phosphate_transporter", "Phosphate transporter", "Phosphate transport", "transport", "transporter"),
    PredicateVocab("sulfate_transporter", "Sulfate transporter", "Sulfate transport", "transport", "transporter"),
    PredicateVocab("nitrate_transporter", "Nitrate transporter", "Nitrate transport", "transport", "transporter"),

    # Localization
    PredicateVocab("secreted", "Secreted", "Contains secretion signal", "transport"),
    PredicateVocab("periplasmic", "Periplasmic", "Periplasmic localization", "transport"),
    PredicateVocab("cell_surface", "Cell surface", "Surface-exposed protein", "transport"),
    PredicateVocab("cytoplasmic", "Cytoplasmic", "Cytoplasmic protein", "transport"),
    PredicateVocab("outer_membrane", "Outer membrane", "Outer membrane protein", "transport", "membrane"),
    PredicateVocab("inner_membrane", "Inner membrane", "Inner/cytoplasmic membrane protein", "transport", "membrane"),
]

# ============================================================================
# REGULATION & SIGNALING
# ============================================================================
REGULATION_PREDICATES = [
    # General regulation
    PredicateVocab("regulator", "Regulator", "Transcriptional regulator", "regulation"),
    PredicateVocab("transcription_factor", "Transcription factor", "DNA-binding transcription factor", "regulation", "regulator"),
    PredicateVocab("repressor", "Repressor", "Transcriptional repressor", "regulation", "regulator"),
    PredicateVocab("activator", "Activator", "Transcriptional activator", "regulation", "regulator"),

    # Two-component systems
    PredicateVocab("two_component", "Two-component system", "Two-component signaling system", "regulation"),
    PredicateVocab("sensor_kinase", "Sensor kinase", "Sensor histidine kinase", "regulation", "two_component"),
    PredicateVocab("response_regulator", "Response regulator", "Response regulator protein", "regulation", "two_component"),

    # Signaling domains
    PredicateVocab("signaling", "Signaling", "Signal transduction protein", "regulation"),
    PredicateVocab("cyclic_dinucleotide", "Cyclic dinucleotide signaling", "c-di-GMP/c-di-AMP signaling", "regulation", "signaling"),
    PredicateVocab("diguanylate_cyclase", "Diguanylate cyclase", "c-di-GMP synthesis (GGDEF)", "regulation", "cyclic_dinucleotide"),
    PredicateVocab("phosphodiesterase", "Phosphodiesterase", "c-di-GMP hydrolysis (EAL/HD-GYP)", "regulation", "cyclic_dinucleotide"),
    PredicateVocab("serine_threonine_kinase", "Ser/Thr kinase", "Serine/threonine protein kinase", "regulation", "signaling"),
    PredicateVocab("tyrosine_kinase", "Tyrosine kinase", "Tyrosine protein kinase", "regulation", "signaling"),
    PredicateVocab("phosphoprotein_phosphatase", "Phosphoprotein phosphatase", "Protein dephosphorylation", "regulation", "signaling"),

    # Regulator families
    PredicateVocab("lysr_family", "LysR family", "LysR-type transcriptional regulator", "regulation", "transcription_factor"),
    PredicateVocab("tetr_family", "TetR family", "TetR-type transcriptional regulator", "regulation", "transcription_factor"),
    PredicateVocab("gntr_family", "GntR family", "GntR-type transcriptional regulator", "regulation", "transcription_factor"),
    PredicateVocab("arac_family", "AraC family", "AraC-type transcriptional regulator", "regulation", "transcription_factor"),
    PredicateVocab("luxr_family", "LuxR family", "LuxR-type transcriptional regulator", "regulation", "transcription_factor"),
    PredicateVocab("marr_family", "MarR family", "MarR-type transcriptional regulator", "regulation", "transcription_factor"),
    PredicateVocab("laci_family", "LacI family", "LacI-type transcriptional regulator", "regulation", "transcription_factor"),
    PredicateVocab("iclr_family", "IclR family", "IclR-type transcriptional regulator", "regulation", "transcription_factor"),
    PredicateVocab("asnc_family", "AsnC family", "AsnC/Lrp-type transcriptional regulator", "regulation", "transcription_factor"),
    PredicateVocab("crp_fnr_family", "CRP/FNR family", "CRP/FNR-type transcriptional regulator", "regulation", "transcription_factor"),
    PredicateVocab("ompr_family", "OmpR family", "OmpR-type response regulator", "regulation", "response_regulator"),
    PredicateVocab("narl_family", "NarL family", "NarL-type response regulator", "regulation", "response_regulator"),

    # Sigma factors
    PredicateVocab("sigma_factor", "Sigma factor", "RNA polymerase sigma factor", "regulation"),
    PredicateVocab("anti_sigma", "Anti-sigma factor", "Sigma factor antagonist", "regulation"),
]

# ============================================================================
# METABOLIC PATHWAYS
# ============================================================================
METABOLISM_PREDICATES = [
    # Central carbon metabolism
    PredicateVocab("central_metabolism", "Central metabolism", "Central carbon metabolism enzyme", "metabolism"),
    PredicateVocab("glycolysis", "Glycolysis", "Glycolytic pathway enzyme", "metabolism", "central_metabolism"),
    PredicateVocab("gluconeogenesis", "Gluconeogenesis", "Gluconeogenic enzyme", "metabolism", "central_metabolism"),
    PredicateVocab("tca_cycle", "TCA cycle", "Tricarboxylic acid cycle enzyme", "metabolism", "central_metabolism"),
    PredicateVocab("pentose_phosphate", "Pentose phosphate pathway", "PPP enzyme", "metabolism", "central_metabolism"),
    PredicateVocab("entner_doudoroff", "Entner-Doudoroff pathway", "ED pathway enzyme", "metabolism", "central_metabolism"),
    PredicateVocab("glyoxylate_cycle", "Glyoxylate cycle", "Glyoxylate bypass enzyme", "metabolism", "central_metabolism"),

    # Energy metabolism
    PredicateVocab("energy_metabolism", "Energy metabolism", "Energy generation/conservation", "metabolism"),
    PredicateVocab("respiration", "Respiration", "Respiratory chain component", "metabolism", "energy_metabolism"),
    PredicateVocab("aerobic_respiration", "Aerobic respiration", "Oxygen-dependent respiration", "metabolism", "respiration"),
    PredicateVocab("anaerobic_respiration", "Anaerobic respiration", "Anaerobic electron transport", "metabolism", "respiration"),
    PredicateVocab("atp_synthesis", "ATP synthesis", "ATP synthase component", "metabolism", "energy_metabolism"),
    PredicateVocab("fermentation", "Fermentation", "Fermentative metabolism", "metabolism", "energy_metabolism"),
    PredicateVocab("electron_transport", "Electron transport", "Electron transport chain", "metabolism", "energy_metabolism"),
    PredicateVocab("cytochrome", "Cytochrome", "Cytochrome protein (heme-containing electron carrier)", "metabolism", "electron_transport"),
    PredicateVocab("ferredoxin", "Ferredoxin", "Ferredoxin (iron-sulfur electron carrier)", "metabolism", "electron_transport"),

    # Amino acid metabolism
    PredicateVocab("amino_acid_metabolism", "Amino acid metabolism", "Amino acid biosynthesis/catabolism", "metabolism"),
    PredicateVocab("amino_acid_biosynthesis", "Amino acid biosynthesis", "Amino acid synthesis", "metabolism", "amino_acid_metabolism"),
    PredicateVocab("amino_acid_degradation", "Amino acid degradation", "Amino acid catabolism", "metabolism", "amino_acid_metabolism"),
    PredicateVocab("aromatic_aa_metabolism", "Aromatic AA metabolism", "Aromatic amino acid pathways", "metabolism", "amino_acid_metabolism"),
    PredicateVocab("branched_chain_aa", "Branched-chain AA metabolism", "BCAA biosynthesis/degradation", "metabolism", "amino_acid_metabolism"),

    # Lipid metabolism
    PredicateVocab("lipid_metabolism", "Lipid metabolism", "Fatty acid and lipid pathways", "metabolism"),
    PredicateVocab("fatty_acid_synthesis", "Fatty acid synthesis", "FAS pathway", "metabolism", "lipid_metabolism"),
    PredicateVocab("fatty_acid_degradation", "Fatty acid degradation", "Beta-oxidation", "metabolism", "lipid_metabolism"),
    PredicateVocab("phospholipid_metabolism", "Phospholipid metabolism", "Membrane lipid synthesis", "metabolism", "lipid_metabolism"),
    PredicateVocab("isoprenoid_biosynthesis", "Isoprenoid biosynthesis", "Terpenoid backbone synthesis", "metabolism", "lipid_metabolism"),

    # Nucleotide metabolism
    PredicateVocab("nucleotide_metabolism", "Nucleotide metabolism", "Purine/pyrimidine pathways", "metabolism"),
    PredicateVocab("purine_metabolism", "Purine metabolism", "Purine biosynthesis/salvage", "metabolism", "nucleotide_metabolism"),
    PredicateVocab("pyrimidine_metabolism", "Pyrimidine metabolism", "Pyrimidine biosynthesis/salvage", "metabolism", "nucleotide_metabolism"),

    # Cofactor biosynthesis
    PredicateVocab("cofactor_biosynthesis", "Cofactor biosynthesis", "Vitamin/cofactor synthesis", "metabolism"),
    PredicateVocab("nad_biosynthesis", "NAD biosynthesis", "NAD/NADP synthesis", "metabolism", "cofactor_biosynthesis"),
    PredicateVocab("fad_biosynthesis", "FAD biosynthesis", "FAD/FMN synthesis", "metabolism", "cofactor_biosynthesis"),
    PredicateVocab("cobalamin_biosynthesis", "Cobalamin biosynthesis", "Vitamin B12 synthesis", "metabolism", "cofactor_biosynthesis"),
    PredicateVocab("folate_biosynthesis", "Folate biosynthesis", "Folate/tetrahydrofolate synthesis", "metabolism", "cofactor_biosynthesis"),
    PredicateVocab("thiamine_biosynthesis", "Thiamine biosynthesis", "Vitamin B1 synthesis", "metabolism", "cofactor_biosynthesis"),
    PredicateVocab("pyridoxal_biosynthesis", "Pyridoxal biosynthesis", "Vitamin B6 synthesis", "metabolism", "cofactor_biosynthesis"),
    PredicateVocab("heme_biosynthesis", "Heme biosynthesis", "Heme/porphyrin synthesis", "metabolism", "cofactor_biosynthesis"),
    PredicateVocab("iron_sulfur_biosynthesis", "Fe-S cluster biosynthesis", "Iron-sulfur cluster assembly", "metabolism", "cofactor_biosynthesis"),

    # Nitrogen metabolism
    PredicateVocab("nitrogen_metabolism", "Nitrogen metabolism", "Nitrogen assimilation/dissimilation", "metabolism"),
    PredicateVocab("nitrogen_fixation", "Nitrogen fixation", "N2 fixation - requires nifHDK", "metabolism", "nitrogen_metabolism"),
    PredicateVocab("nitrogenase_maturation", "Nitrogenase maturation", "FeMo-co biosynthesis (nifBNE)", "metabolism", "nitrogen_metabolism"),
    PredicateVocab("nitrification", "Nitrification", "Ammonia to nitrate oxidation", "metabolism", "nitrogen_metabolism"),
    PredicateVocab("denitrification", "Denitrification", "Nitrate to N2 reduction", "metabolism", "nitrogen_metabolism"),
    PredicateVocab("ammonia_assimilation", "Ammonia assimilation", "Ammonia incorporation", "metabolism", "nitrogen_metabolism"),
    PredicateVocab("nitrate_reduction", "Nitrate reduction", "Nitrate reductase activity", "metabolism", "nitrogen_metabolism"),
    PredicateVocab("urea_metabolism", "Urea metabolism", "Urease/urea cycle", "metabolism", "nitrogen_metabolism"),

    # Sulfur metabolism
    PredicateVocab("sulfur_metabolism", "Sulfur metabolism", "Sulfur compound metabolism", "metabolism"),
    PredicateVocab("sulfate_reduction", "Sulfate reduction", "Dissimilatory sulfate reduction", "metabolism", "sulfur_metabolism"),
    PredicateVocab("sulfur_oxidation", "Sulfur oxidation", "Sulfur compound oxidation", "metabolism", "sulfur_metabolism"),
    PredicateVocab("sulfur_assimilation", "Sulfur assimilation", "Sulfate to cysteine", "metabolism", "sulfur_metabolism"),

    # One-carbon metabolism
    PredicateVocab("one_carbon_metabolism", "One-carbon metabolism", "C1 compound metabolism", "metabolism"),
    PredicateVocab("methanogenesis", "Methanogenesis", "Methane production - requires MCR complex", "metabolism", "one_carbon_metabolism"),
    PredicateVocab("mcr_complex", "MCR complex", "Methyl-CoM reductase - definitive methanogen marker", "metabolism", "methanogenesis"),
    PredicateVocab("archaeal_one_carbon", "Archaeal C1 metabolism", "Tetrahydromethanopterin pathway (non-methanogen)", "metabolism", "one_carbon_metabolism"),
    PredicateVocab("f420_dependent", "F420-dependent", "Coenzyme F420 utilizing enzyme", "metabolism", "one_carbon_metabolism"),
    PredicateVocab("methanotrophy", "Methanotrophy", "Methane oxidation", "metabolism", "one_carbon_metabolism"),
    PredicateVocab("methylotrophy", "Methylotrophy", "Methanol/methylamine oxidation", "metabolism", "one_carbon_metabolism"),
    PredicateVocab("wood_ljungdahl", "Wood-Ljungdahl pathway", "Reductive acetyl-CoA pathway", "metabolism", "one_carbon_metabolism"),
    PredicateVocab("co_oxidation", "CO oxidation", "Carbon monoxide oxidation", "metabolism", "one_carbon_metabolism"),

    # Carbon fixation
    PredicateVocab("carbon_fixation", "Carbon fixation", "CO2 fixation pathway", "metabolism"),
    PredicateVocab("rubisco", "RuBisCO", "Ribulose-1,5-bisphosphate carboxylase (does NOT imply Calvin cycle)", "metabolism", "carbon_fixation"),
    PredicateVocab("prk", "Phosphoribulokinase", "PRK - true Calvin cycle marker", "metabolism", "calvin_cycle"),
    PredicateVocab("calvin_cycle", "Calvin cycle", "RuBisCO-based carbon fixation (requires PRK)", "metabolism", "carbon_fixation"),
    PredicateVocab("reverse_tca", "Reverse TCA cycle", "rTCA carbon fixation", "metabolism", "carbon_fixation"),
    PredicateVocab("3hp_bicycle", "3-HP bicycle", "3-hydroxypropionate cycle", "metabolism", "carbon_fixation"),
    PredicateVocab("dicarboxylate_4hb", "DC/4-HB cycle", "Dicarboxylate/4-hydroxybutyrate", "metabolism", "carbon_fixation"),

    # Photosynthesis
    PredicateVocab("photosynthesis", "Photosynthesis", "Light-dependent energy capture", "metabolism"),
    PredicateVocab("photosystem_i", "Photosystem I", "PSI component", "metabolism", "photosynthesis"),
    PredicateVocab("photosystem_ii", "Photosystem II", "PSII component", "metabolism", "photosynthesis"),
    PredicateVocab("light_harvesting", "Light harvesting", "Antenna complex", "metabolism", "photosynthesis"),
    PredicateVocab("bacteriochlorophyll", "Bacteriochlorophyll", "Bacterial photosynthesis", "metabolism", "photosynthesis"),

    # Secondary metabolism
    PredicateVocab("secondary_metabolism", "Secondary metabolism", "Natural product biosynthesis", "metabolism"),
    PredicateVocab("polyketide_synthesis", "Polyketide synthesis", "PKS pathway", "metabolism", "secondary_metabolism"),
    PredicateVocab("nrps", "NRPS", "Non-ribosomal peptide synthesis", "metabolism", "secondary_metabolism"),
    PredicateVocab("terpene_synthesis", "Terpene synthesis", "Terpenoid biosynthesis", "metabolism", "secondary_metabolism"),
    PredicateVocab("siderophore_biosynthesis", "Siderophore biosynthesis", "Iron chelator synthesis", "metabolism", "secondary_metabolism"),
    PredicateVocab("antibiotic_biosynthesis", "Antibiotic biosynthesis", "Antimicrobial compound synthesis", "metabolism", "secondary_metabolism"),

    # Hydrogen metabolism
    PredicateVocab("hydrogen_metabolism", "Hydrogen metabolism", "H2 production/consumption", "metabolism"),
    PredicateVocab("hydrogenase", "Hydrogenase", "H2 oxidation/production enzyme", "metabolism", "hydrogen_metabolism"),
    PredicateVocab("hyddb_annotated", "HydDB annotated", "Annotated by HydDB HMMs", "metabolism", "hydrogenase"),

    # NiFe hydrogenases (Ni-Fe active site)
    PredicateVocab("nife_hydrogenase", "NiFe hydrogenase", "Nickel-iron hydrogenase", "metabolism", "hydrogenase"),
    PredicateVocab("nifese_hydrogenase", "NiFeSe hydrogenase", "Selenium-containing NiFe hydrogenase", "metabolism", "nife_hydrogenase"),

    # NiFe Group 1: Membrane-bound H2 uptake
    PredicateVocab("nife_group1", "NiFe Group 1", "Membrane-bound respiratory uptake hydrogenases", "metabolism", "nife_hydrogenase"),
    PredicateVocab("uptake_hydrogenase", "Uptake hydrogenase", "H2-oxidizing respiratory hydrogenase", "metabolism", "nife_group1"),

    # NiFe Group 2: Cytoplasmic H2 sensors and uptake
    PredicateVocab("nife_group2", "NiFe Group 2", "Cytoplasmic H2 sensors and cyanobacterial uptake", "metabolism", "nife_hydrogenase"),
    PredicateVocab("cytoplasmic_hydrogenase", "Cytoplasmic hydrogenase", "Soluble cytoplasmic hydrogenase", "metabolism", "hydrogenase"),
    PredicateVocab("h2_sensor", "H2 sensor", "Regulatory H2-sensing hydrogenase", "metabolism", "nife_group2"),

    # NiFe Group 3: Bidirectional/cofactor-coupled
    PredicateVocab("nife_group3", "NiFe Group 3", "Bidirectional cytoplasmic hydrogenases", "metabolism", "nife_hydrogenase"),
    PredicateVocab("bidirectional_hydrogenase", "Bidirectional hydrogenase", "Reversible H2 production/consumption", "metabolism", "hydrogenase"),
    PredicateVocab("f420_reducing", "F420-reducing", "Coenzyme F420-coupled hydrogenase", "metabolism", "nife_group3"),
    PredicateVocab("nadp_reducing", "NAD(P)-reducing", "NAD(P)-coupled hydrogenase", "metabolism", "nife_group3"),
    PredicateVocab("methyl_viologen_reducing", "MV-reducing", "Methyl viologen-reducing hydrogenase", "metabolism", "nife_group3"),

    # NiFe Group 4: Energy-conserving H2-evolving
    PredicateVocab("nife_group4", "NiFe Group 4", "Energy-conserving membrane-bound H2-evolving", "metabolism", "nife_hydrogenase"),
    PredicateVocab("h2_evolving", "H2-evolving", "Proton-reducing H2-producing hydrogenase", "metabolism", "hydrogenase"),
    PredicateVocab("energy_conserving_hydrogenase", "Energy-conserving", "Couples H2 evolution to ion gradient", "metabolism", "nife_group4"),
    PredicateVocab("formate_coupled", "Formate-coupled", "Formate hydrogen lyase complex", "metabolism", "nife_group4"),
    PredicateVocab("co_coupled", "CO-coupled", "CO dehydrogenase-coupled hydrogenase", "metabolism", "nife_group4"),
    PredicateVocab("mbh_hydrogenase", "Mbh hydrogenase", "Membrane-bound hydrogenase complex", "metabolism", "nife_group4"),
    PredicateVocab("ech_hydrogenase", "Ech hydrogenase", "Energy-conserving hydrogenase (Ech)", "metabolism", "nife_group4"),

    # FeFe hydrogenases (Fe-Fe active site, H-cluster)
    PredicateVocab("fefe_hydrogenase", "FeFe hydrogenase", "Iron-iron hydrogenase (H-cluster)", "metabolism", "hydrogenase"),

    # FeFe Group A: Monomeric fermentative
    PredicateVocab("fefe_groupA", "FeFe Group A", "Monomeric cytoplasmic fermentative hydrogenases", "metabolism", "fefe_hydrogenase"),
    PredicateVocab("monomeric_fefe", "Monomeric FeFe", "Single-subunit FeFe hydrogenase", "metabolism", "fefe_groupA"),
    PredicateVocab("fermentative_hydrogenase", "Fermentative hydrogenase", "H2-evolving during fermentation", "metabolism", "h2_evolving"),

    # FeFe Group B: Electron-bifurcating
    PredicateVocab("fefe_groupB", "FeFe Group B", "Electron-bifurcating hydrogenases", "metabolism", "fefe_hydrogenase"),
    PredicateVocab("bifurcating_hydrogenase", "Bifurcating hydrogenase", "Electron-bifurcating hydrogenase", "metabolism", "fefe_groupB"),

    # FeFe Group C: Sensory/regulatory
    PredicateVocab("fefe_groupC", "FeFe Group C", "Sensory/regulatory FeFe hydrogenases", "metabolism", "fefe_hydrogenase"),
    PredicateVocab("sensory_hydrogenase", "Sensory hydrogenase", "Regulatory H2-sensing FeFe", "metabolism", "fefe_groupC"),

    # Fe-only hydrogenases (Hmd)
    PredicateVocab("fe_only_hydrogenase", "Fe-only hydrogenase", "Iron-only Hmd hydrogenase (methanogens)", "metabolism", "hydrogenase"),
    PredicateVocab("methanogen_hydrogenase", "Methanogen hydrogenase", "Methanogen-associated Hmd", "metabolism", "fe_only_hydrogenase"),

    # Hydrogenase accessory/maturation
    PredicateVocab("hydrogenase_maturation", "Hydrogenase maturation", "Hydrogenase assembly factors (HypABCDEF)", "metabolism", "hydrogen_metabolism"),
    PredicateVocab("membrane_bound_hydrogenase", "Membrane-bound", "Hydrogenase with membrane subunits", "metabolism", "hydrogenase"),
]

# ============================================================================
# CARBOHYDRATE-ACTIVE ENZYMES (CAZy)
# ============================================================================
CAZY_PREDICATES = [
    # Main classes
    PredicateVocab("carbohydrate_active", "Carbohydrate-active", "Any CAZy enzyme", "cazy"),
    PredicateVocab("glycoside_hydrolase", "Glycoside hydrolase", "GH family enzyme", "cazy", "carbohydrate_active"),
    PredicateVocab("glycosyltransferase_cazy", "Glycosyltransferase (CAZy)", "GT family enzyme", "cazy", "carbohydrate_active"),
    PredicateVocab("polysaccharide_lyase", "Polysaccharide lyase", "PL family enzyme", "cazy", "carbohydrate_active"),
    PredicateVocab("carbohydrate_esterase", "Carbohydrate esterase", "CE family enzyme", "cazy", "carbohydrate_active"),
    PredicateVocab("carbohydrate_binding", "Carbohydrate-binding module", "CBM family", "cazy", "carbohydrate_active"),
    PredicateVocab("auxiliary_activity", "Auxiliary activity", "AA family enzyme", "cazy", "carbohydrate_active"),

    # Substrate-specific
    PredicateVocab("cellulase", "Cellulase", "Cellulose degradation", "cazy", "glycoside_hydrolase"),
    PredicateVocab("chitinase", "Chitinase", "Chitin degradation", "cazy", "glycoside_hydrolase"),
    PredicateVocab("amylase", "Amylase", "Starch degradation", "cazy", "glycoside_hydrolase"),
    PredicateVocab("xylanase", "Xylanase", "Xylan degradation", "cazy", "glycoside_hydrolase"),
    PredicateVocab("pectinase", "Pectinase", "Pectin degradation", "cazy"),
    PredicateVocab("mannanase", "Mannanase", "Mannan degradation", "cazy", "glycoside_hydrolase"),
    PredicateVocab("lytic_polysaccharide_monooxygenase", "LPMO", "Lytic polysaccharide monooxygenase", "cazy", "auxiliary_activity"),
]

# ============================================================================
# STRUCTURAL & PROTEIN FEATURES
# ============================================================================
STRUCTURE_PREDICATES = [
    # Domain architecture
    PredicateVocab("single_domain", "Single domain", "Single domain protein", "structure"),
    PredicateVocab("multi_domain", "Multi-domain", "3+ distinct domains", "structure"),

    # Repeat domains
    PredicateVocab("repeat_domain", "Repeat domain", "Contains repeat elements", "structure"),
    PredicateVocab("tpr_repeat", "TPR repeat", "Tetratricopeptide repeat", "structure", "repeat_domain"),
    PredicateVocab("wd40_repeat", "WD40 repeat", "WD40/beta-propeller repeat", "structure", "repeat_domain"),
    PredicateVocab("lrr_repeat", "LRR repeat", "Leucine-rich repeat", "structure", "repeat_domain"),
    PredicateVocab("ankyrin_repeat", "Ankyrin repeat", "Ankyrin repeat domain", "structure", "repeat_domain"),
    PredicateVocab("kelch_repeat", "Kelch repeat", "Kelch repeat domain", "structure", "repeat_domain"),
    PredicateVocab("sel1_repeat", "Sel1 repeat", "Sel1-like repeat", "structure", "repeat_domain"),

    # Structural features
    PredicateVocab("coiled_coil", "Coiled-coil", "Coiled-coil structure", "structure"),
    PredicateVocab("beta_barrel", "Beta barrel", "Beta-barrel fold", "structure"),
    PredicateVocab("alpha_helical", "Alpha helical", "Predominantly alpha-helical", "structure"),
    PredicateVocab("intrinsically_disordered", "Intrinsically disordered", "Low complexity/disordered", "structure"),
]

# ============================================================================
# BINDING & INTERACTION
# ============================================================================
BINDING_PREDICATES = [
    # Nucleic acid binding
    PredicateVocab("dna_binding", "DNA binding", "DNA-binding domain", "binding"),
    PredicateVocab("rna_binding", "RNA binding", "RNA-binding domain", "binding"),
    PredicateVocab("helix_turn_helix", "Helix-turn-helix", "HTH DNA-binding motif", "binding", "dna_binding"),
    PredicateVocab("zinc_finger", "Zinc finger", "Zinc finger DNA-binding", "binding", "dna_binding"),
    PredicateVocab("winged_helix", "Winged helix", "Winged helix-turn-helix", "binding", "dna_binding"),
    PredicateVocab("ribbon_helix_helix", "Ribbon-helix-helix", "RHH DNA-binding", "binding", "dna_binding"),

    # Nucleotide binding
    PredicateVocab("nucleotide_binding", "Nucleotide binding", "NTP/NDP binding", "binding"),
    PredicateVocab("atp_binding", "ATP binding", "ATP/ADP binding", "binding", "nucleotide_binding"),
    PredicateVocab("gtp_binding", "GTP binding", "GTP/GDP binding", "binding", "nucleotide_binding"),
    PredicateVocab("p_loop", "P-loop", "Walker A motif NTPase", "binding", "nucleotide_binding"),
    PredicateVocab("aaa_domain", "AAA+ domain", "AAA+ ATPase domain", "binding", "atp_binding"),

    # Cofactor binding
    PredicateVocab("cofactor_binding", "Cofactor binding", "Cofactor-binding domain", "binding"),
    PredicateVocab("nad_binding", "NAD binding", "NAD(P)-binding domain", "binding", "cofactor_binding"),
    PredicateVocab("fad_binding", "FAD binding", "FAD-binding domain", "binding", "cofactor_binding"),
    PredicateVocab("fmn_binding", "FMN binding", "FMN-binding domain", "binding", "cofactor_binding"),
    PredicateVocab("plp_binding", "PLP binding", "Pyridoxal phosphate binding", "binding", "cofactor_binding"),
    PredicateVocab("coenzyme_a_binding", "CoA binding", "Coenzyme A binding", "binding", "cofactor_binding"),
    PredicateVocab("biotin_binding", "Biotin binding", "Biotin-binding domain", "binding", "cofactor_binding"),
    PredicateVocab("thiamine_binding", "Thiamine binding", "Thiamine pyrophosphate binding", "binding", "cofactor_binding"),

    # Metal binding
    PredicateVocab("metal_binding", "Metal binding", "Metal ion binding site", "binding"),
    PredicateVocab("iron_binding", "Iron binding", "Iron/heme binding", "binding", "metal_binding"),
    PredicateVocab("iron_sulfur", "Iron-sulfur cluster", "Fe-S cluster binding", "binding", "metal_binding"),
    PredicateVocab("heme_binding", "Heme binding", "Heme/porphyrin binding", "binding", "iron_binding"),
    PredicateVocab("zinc_binding", "Zinc binding", "Zinc ion binding", "binding", "metal_binding"),
    PredicateVocab("copper_binding", "Copper binding", "Copper ion binding", "binding", "metal_binding"),
    PredicateVocab("manganese_binding", "Manganese binding", "Manganese ion binding", "binding", "metal_binding"),
    PredicateVocab("calcium_binding", "Calcium binding", "Calcium ion binding", "binding", "metal_binding"),
    PredicateVocab("magnesium_binding", "Magnesium binding", "Magnesium ion binding", "binding", "metal_binding"),
    PredicateVocab("molybdenum_binding", "Molybdenum binding", "Molybdopterin cofactor", "binding", "metal_binding"),
    PredicateVocab("nickel_binding", "Nickel binding", "Nickel ion binding", "binding", "metal_binding"),
    PredicateVocab("cobalt_binding", "Cobalt binding", "Cobalt/cobalamin binding", "binding", "metal_binding"),
    PredicateVocab("cobalamin_binding", "Cobalamin binding", "Vitamin B12 binding", "binding", "cobalt_binding"),

    # Protein-protein interaction
    PredicateVocab("protein_binding", "Protein binding", "Protein-protein interaction", "binding"),
    PredicateVocab("ubiquitin_like", "Ubiquitin-like", "Ubiquitin/UBL domain", "binding", "protein_binding"),
]

# ============================================================================
# CELL ENVELOPE & SURFACE
# ============================================================================
ENVELOPE_PREDICATES = [
    # Cell wall
    PredicateVocab("cell_wall", "Cell wall", "Cell wall synthesis/modification", "envelope"),
    PredicateVocab("peptidoglycan", "Peptidoglycan", "Peptidoglycan metabolism", "envelope", "cell_wall"),
    PredicateVocab("murein_synthesis", "Murein synthesis", "PG biosynthesis (Mur enzymes)", "envelope", "peptidoglycan"),
    PredicateVocab("pbp", "Penicillin-binding protein", "PG transpeptidase/transglycosylase", "envelope", "peptidoglycan"),
    PredicateVocab("lytic_transglycosylase", "Lytic transglycosylase", "PG lytic enzyme", "envelope", "peptidoglycan"),

    # LPS/lipid A
    PredicateVocab("lps_biosynthesis", "LPS biosynthesis", "Lipopolysaccharide synthesis", "envelope"),
    PredicateVocab("lipid_a_biosynthesis", "Lipid A biosynthesis", "Lipid A/endotoxin synthesis", "envelope", "lps_biosynthesis"),
    PredicateVocab("o_antigen", "O-antigen", "O-antigen biosynthesis", "envelope", "lps_biosynthesis"),

    # Capsule/EPS
    PredicateVocab("capsule", "Capsule", "Capsule biosynthesis", "envelope"),
    PredicateVocab("exopolysaccharide", "Exopolysaccharide", "EPS biosynthesis", "envelope"),

    # Surface structures
    PredicateVocab("adhesin", "Adhesin", "Surface adhesion protein", "envelope"),
    PredicateVocab("pilus", "Pilus", "Pilus/fimbria component", "envelope"),
    PredicateVocab("type_iv_pilus", "Type IV pilus", "Type IV pilus system", "envelope", "pilus"),
    PredicateVocab("flagellum", "Flagellum", "Flagellar component", "envelope"),
    PredicateVocab("flagellar_motor", "Flagellar motor", "Flagellar motor/switch", "envelope", "flagellum"),
    PredicateVocab("flagellar_hook", "Flagellar hook", "Flagellar hook protein", "envelope", "flagellum"),
    PredicateVocab("flagellin", "Flagellin", "Flagellar filament protein", "envelope", "flagellum"),
    PredicateVocab("chemotaxis", "Chemotaxis", "Chemotaxis signaling", "envelope"),

    # S-layer
    PredicateVocab("s_layer", "S-layer", "Surface layer protein", "envelope"),
]

# ============================================================================
# MOBILE GENETIC ELEMENTS & DEFENSE
# ============================================================================
MOBILE_PREDICATES = [
    # Mobile elements
    PredicateVocab("mobile_element", "Mobile element", "Mobile genetic element", "mobile"),
    PredicateVocab("transposase", "Transposase", "Transposable element enzyme", "mobile", "mobile_element"),
    PredicateVocab("integrase", "Integrase", "Site-specific integrase", "mobile", "mobile_element"),
    PredicateVocab("resolvase", "Resolvase", "Site-specific resolvase", "mobile", "mobile_element"),
    PredicateVocab("recombinase", "Recombinase", "Site-specific recombinase", "mobile", "mobile_element"),
    PredicateVocab("insertion_sequence", "Insertion sequence", "IS element protein", "mobile", "transposase"),

    # Phage-related
    PredicateVocab("phage_related", "Phage-related", "Phage-derived protein", "mobile"),
    PredicateVocab("phage_integrase", "Phage integrase", "Phage integration enzyme", "mobile", "phage_related"),
    PredicateVocab("phage_terminase", "Phage terminase", "DNA packaging enzyme", "mobile", "phage_related"),
    PredicateVocab("phage_portal", "Phage portal", "Capsid portal protein", "mobile", "phage_related"),
    PredicateVocab("phage_capsid", "Phage capsid", "Capsid structural protein", "mobile", "phage_related"),
    PredicateVocab("phage_tail", "Phage tail", "Tail structural protein", "mobile", "phage_related"),
    PredicateVocab("phage_baseplate", "Phage baseplate", "Baseplate component", "mobile", "phage_related"),
    PredicateVocab("phage_lysin", "Phage lysin", "Phage lytic enzyme", "mobile", "phage_related"),
    PredicateVocab("holin", "Holin", "Membrane pore for lysis", "mobile", "phage_related"),

    # Defense systems
    PredicateVocab("defense_system", "Defense system", "Anti-phage/MGE defense", "mobile"),
    PredicateVocab("restriction_modification", "Restriction-modification", "R-M system component", "mobile", "defense_system"),
    PredicateVocab("restriction_enzyme", "Restriction enzyme", "Restriction endonuclease", "mobile", "restriction_modification"),
    PredicateVocab("methyltransferase_rm", "DNA methyltransferase (R-M)", "R-M system methylase", "mobile", "restriction_modification"),
    PredicateVocab("cas_domain", "Cas-like domain", "Domain found in Cas proteins (may also be transposase)", "mobile", "defense_system"),
    PredicateVocab("crispr_associated", "CRISPR-associated", "Confirmed CRISPR-Cas system component", "mobile", "cas_domain"),
    PredicateVocab("cas_nuclease", "Cas nuclease", "CRISPR effector nuclease", "mobile", "crispr_associated"),

    # CRISPR-Cas Class 1 (multi-subunit effector)
    PredicateVocab("crispr_class1", "CRISPR Class 1", "Multi-subunit effector complex (Types I, III, IV)", "mobile", "crispr_associated"),

    # Type I - Cas3 signature (helicase-nuclease)
    PredicateVocab("crispr_type_i", "CRISPR Type I", "Cas3-based interference", "mobile", "crispr_class1"),
    PredicateVocab("crispr_type_i_a", "CRISPR Type I-A", "Type I-A (Cascade variant)", "mobile", "crispr_type_i"),
    PredicateVocab("crispr_type_i_b", "CRISPR Type I-B", "Type I-B (Cascade variant)", "mobile", "crispr_type_i"),
    PredicateVocab("crispr_type_i_c", "CRISPR Type I-C", "Type I-C (minimal Cascade)", "mobile", "crispr_type_i"),
    PredicateVocab("crispr_type_i_d", "CRISPR Type I-D", "Type I-D (Cyanobacteria)", "mobile", "crispr_type_i"),
    PredicateVocab("crispr_type_i_e", "CRISPR Type I-E", "Type I-E (E. coli Cascade)", "mobile", "crispr_type_i"),
    PredicateVocab("crispr_type_i_f", "CRISPR Type I-F", "Type I-F (Pseudomonas)", "mobile", "crispr_type_i"),
    PredicateVocab("crispr_type_i_g", "CRISPR Type I-G", "Type I-G variant", "mobile", "crispr_type_i"),

    # Type III - Cas10 signature (RNA-targeting)
    PredicateVocab("crispr_type_iii", "CRISPR Type III", "Cas10-based RNA interference", "mobile", "crispr_class1"),
    PredicateVocab("crispr_type_iii_a", "CRISPR Type III-A", "Type III-A Csm complex (DNA+RNA targeting)", "mobile", "crispr_type_iii"),
    PredicateVocab("crispr_type_iii_b", "CRISPR Type III-B", "Type III-B Cmr complex (RNA targeting)", "mobile", "crispr_type_iii"),
    PredicateVocab("crispr_type_iii_c", "CRISPR Type III-C", "Type III-C variant", "mobile", "crispr_type_iii"),
    PredicateVocab("crispr_type_iii_d", "CRISPR Type III-D", "Type III-D variant", "mobile", "crispr_type_iii"),

    # Type IV - Csf signature (minimal, no nuclease)
    PredicateVocab("crispr_type_iv", "CRISPR Type IV", "Csf-based system (no nuclease)", "mobile", "crispr_class1"),

    # CRISPR-Cas Class 2 (single effector)
    PredicateVocab("crispr_class2", "CRISPR Class 2", "Single-effector systems (Types II, V, VI)", "mobile", "crispr_associated"),

    # Type II - Cas9 signature
    PredicateVocab("crispr_type_ii", "CRISPR Type II", "Cas9-based DNA interference", "mobile", "crispr_class2"),
    PredicateVocab("crispr_type_ii_a", "CRISPR Type II-A", "Type II-A (canonical Cas9)", "mobile", "crispr_type_ii"),
    PredicateVocab("crispr_type_ii_b", "CRISPR Type II-B", "Type II-B variant", "mobile", "crispr_type_ii"),
    PredicateVocab("crispr_type_ii_c", "CRISPR Type II-C", "Type II-C (minimal)", "mobile", "crispr_type_ii"),

    # Type V - Cas12 signature
    PredicateVocab("crispr_type_v", "CRISPR Type V", "Cas12-based DNA interference", "mobile", "crispr_class2"),
    PredicateVocab("crispr_type_v_a", "CRISPR Type V-A", "Type V-A (Cas12a/Cpf1)", "mobile", "crispr_type_v"),
    PredicateVocab("crispr_type_v_b", "CRISPR Type V-B", "Type V-B (Cas12b/C2c1)", "mobile", "crispr_type_v"),

    # Type VI - Cas13 signature (RNA-targeting)
    PredicateVocab("crispr_type_vi", "CRISPR Type VI", "Cas13-based RNA interference", "mobile", "crispr_class2"),

    # Core adaptation module
    PredicateVocab("crispr_adaptation", "CRISPR adaptation", "Cas1/Cas2 adaptation module", "mobile", "crispr_associated"),

    # Accessory/regulatory
    PredicateVocab("crispr_accessory", "CRISPR accessory", "CRISPR-associated accessory protein", "mobile", "crispr_associated"),

    PredicateVocab("in_crispr_array", "In CRISPR array", "ORF overlapping CRISPR array (likely spurious)", "mobile"),
    PredicateVocab("toxin_domain", "Toxin-like domain", "Domain found in toxins (may have other functions)", "mobile", "defense_system"),
    PredicateVocab("antitoxin_domain", "Antitoxin-like domain", "Domain found in antitoxins", "mobile", "defense_system"),
    PredicateVocab("toxin_antitoxin", "Toxin-antitoxin system", "Confirmed paired TA system", "mobile", "defense_system"),
    PredicateVocab("toxin", "Toxin", "TA system toxin", "mobile", "toxin_antitoxin"),
    PredicateVocab("antitoxin", "Antitoxin", "TA system antitoxin", "mobile", "toxin_antitoxin"),
    PredicateVocab("abortive_infection", "Abortive infection", "Abi defense system", "mobile", "defense_system"),

    # Conjugation
    PredicateVocab("conjugation", "Conjugation", "Conjugative transfer", "mobile"),
    PredicateVocab("relaxase", "Relaxase", "Conjugative relaxase", "mobile", "conjugation"),
    PredicateVocab("type_iv_secretion", "Type IV secretion", "T4SS component", "mobile", "conjugation"),

    # Secretion systems - gene-level component tags
    PredicateVocab("secretion_system", "Secretion system", "Protein secretion apparatus", "mobile"),
    PredicateVocab("secretion_component", "Secretion component", "Component of protein secretion machinery", "mobile", "secretion_system"),
    PredicateVocab("t1ss_component", "T1SS component", "Type I secretion component", "mobile", "secretion_component"),
    PredicateVocab("t2ss_component", "T2SS component", "Type II secretion component", "mobile", "secretion_component"),
    PredicateVocab("t3ss_component", "T3SS component", "Type III secretion component", "mobile", "secretion_component"),
    PredicateVocab("t4ss_component", "T4SS component", "Type IV secretion component", "mobile", "secretion_component"),
    PredicateVocab("t5ss_component", "T5SS component", "Type V secretion/autotransporter component", "mobile", "secretion_component"),
    PredicateVocab("t6ss_component", "T6SS component", "Type VI secretion component", "mobile", "secretion_component"),
    PredicateVocab("effector_domain", "Effector domain", "Secretion system effector domain", "mobile", "secretion_system"),
    # Secretion systems - locus-level calls (from clustering analysis)
    PredicateVocab("type_i_secretion", "Type I secretion", "T1SS locus (confirmed)", "mobile", "secretion_system"),
    PredicateVocab("type_ii_secretion", "Type II secretion", "T2SS locus (confirmed)", "mobile", "secretion_system"),
    PredicateVocab("type_iii_secretion", "Type III secretion", "T3SS locus (confirmed)", "mobile", "secretion_system"),
    PredicateVocab("type_iv_secretion", "Type IV secretion", "T4SS locus (confirmed)", "mobile", "secretion_system"),
    PredicateVocab("type_v_secretion", "Type V secretion", "Autotransporter locus", "mobile", "secretion_system"),
    PredicateVocab("type_vi_secretion", "Type VI secretion", "T6SS locus (confirmed)", "mobile", "secretion_system"),
    PredicateVocab("sec_pathway", "Sec pathway", "General secretion pathway", "mobile", "secretion_system"),
    PredicateVocab("tat_pathway", "Tat pathway", "Twin-arginine translocation", "mobile", "secretion_system"),
]

# ============================================================================
# STRESS & ENVIRONMENT
# ============================================================================
STRESS_PREDICATES = [
    # General stress
    PredicateVocab("stress_response", "Stress response", "General stress response", "stress"),

    # Temperature
    PredicateVocab("heat_shock", "Heat shock", "Heat shock protein", "stress", "stress_response"),
    PredicateVocab("cold_shock", "Cold shock", "Cold shock protein", "stress", "stress_response"),

    # Oxidative stress
    PredicateVocab("oxidative_stress", "Oxidative stress", "ROS defense", "stress", "stress_response"),
    PredicateVocab("catalase", "Catalase", "H2O2 decomposition", "stress", "oxidative_stress"),
    PredicateVocab("superoxide_dismutase", "Superoxide dismutase", "Superoxide detoxification", "stress", "oxidative_stress"),
    PredicateVocab("peroxiredoxin", "Peroxiredoxin", "Peroxide reduction", "stress", "oxidative_stress"),
    PredicateVocab("thioredoxin", "Thioredoxin", "Thiol-disulfide oxidoreductase", "stress", "oxidative_stress"),
    PredicateVocab("glutaredoxin", "Glutaredoxin", "Glutathione-dependent redox", "stress", "oxidative_stress"),

    # Osmotic stress
    PredicateVocab("osmotic_stress", "Osmotic stress", "Osmotic adaptation", "stress", "stress_response"),
    PredicateVocab("osmoprotectant", "Osmoprotectant", "Compatible solute system", "stress", "osmotic_stress"),

    # Chaperones
    PredicateVocab("chaperone", "Chaperone", "Protein folding chaperone", "stress"),
    PredicateVocab("hsp70", "Hsp70", "DnaK/Hsp70 family", "stress", "chaperone"),
    PredicateVocab("hsp60", "Hsp60", "GroEL/Hsp60 chaperonin", "stress", "chaperone"),
    PredicateVocab("hsp90", "Hsp90", "HtpG/Hsp90 family", "stress", "chaperone"),
    PredicateVocab("small_hsp", "Small HSP", "sHsp/alpha-crystallin family", "stress", "chaperone"),
    PredicateVocab("clp_protease", "Clp protease", "Clp protease/chaperone", "stress", "chaperone"),
    PredicateVocab("lon_protease", "Lon protease", "Lon ATP-dependent protease", "stress", "chaperone"),

    # Resistance
    PredicateVocab("antibiotic_resistance", "Antibiotic resistance", "Antimicrobial resistance", "stress"),
    PredicateVocab("beta_lactamase", "Beta-lactamase", "Beta-lactam resistance", "stress", "antibiotic_resistance"),
    PredicateVocab("aminoglycoside_resistance", "Aminoglycoside resistance", "AG-modifying enzyme", "stress", "antibiotic_resistance"),
    PredicateVocab("multidrug_resistance", "Multidrug resistance", "MDR efflux/other", "stress", "antibiotic_resistance"),

    # Heavy metals
    PredicateVocab("heavy_metal_resistance", "Heavy metal resistance", "Metal detoxification", "stress"),
    PredicateVocab("arsenic_resistance", "Arsenic resistance", "Arsenic detox (ars)", "stress", "heavy_metal_resistance"),
    PredicateVocab("mercury_resistance", "Mercury resistance", "Mercury detox (mer)", "stress", "heavy_metal_resistance"),
    PredicateVocab("copper_resistance", "Copper resistance", "Copper homeostasis", "stress", "heavy_metal_resistance"),
    PredicateVocab("zinc_resistance", "Zinc resistance", "Zinc homeostasis", "stress", "heavy_metal_resistance"),

    # DNA repair
    PredicateVocab("dna_repair", "DNA repair", "DNA damage repair", "stress"),
    PredicateVocab("base_excision_repair", "Base excision repair", "BER pathway", "stress", "dna_repair"),
    PredicateVocab("nucleotide_excision_repair", "Nucleotide excision repair", "NER pathway", "stress", "dna_repair"),
    PredicateVocab("mismatch_repair", "Mismatch repair", "MMR pathway", "stress", "dna_repair"),
    PredicateVocab("recombinational_repair", "Recombinational repair", "HR/RecA pathway", "stress", "dna_repair"),
    PredicateVocab("sos_response", "SOS response", "DNA damage response", "stress", "dna_repair"),
]

# ============================================================================
# SIZE-BASED PREDICATES
# ============================================================================
SIZE_PREDICATES = [
    PredicateVocab("tiny", "Tiny", "Protein < 50 aa", "size"),
    PredicateVocab("small", "Small", "Protein 50-150 aa", "size"),
    PredicateVocab("medium", "Medium", "Protein 150-400 aa", "size"),
    PredicateVocab("large", "Large", "Protein 400-1000 aa", "size"),
    PredicateVocab("giant", "Giant", "Protein > 1000 aa", "size"),
    PredicateVocab("massive", "Massive", "Protein > 2000 aa", "size"),
]

# ============================================================================
# ANNOTATION STATUS
# ============================================================================
ANNOTATION_PREDICATES = [
    PredicateVocab("unannotated", "Unannotated", "No functional annotations", "annotation"),
    PredicateVocab("hypothetical", "Hypothetical", "Hypothetical/uncharacterized annotation", "annotation"),
    PredicateVocab("well_annotated", "Well annotated", "3+ confident annotations", "annotation"),
    PredicateVocab("confident_hit", "Confident hit", "E-value < 1e-10", "annotation"),
    PredicateVocab("weak_hit", "Weak hit", "Only weak annotations (e-value > 1e-5)", "annotation"),
    PredicateVocab("single_source", "Single source", "Annotations from only one source", "annotation"),
    PredicateVocab("multi_source", "Multi source", "Annotations from multiple sources", "annotation"),
    PredicateVocab("pfam_annotated", "PFAM annotated", "Has PFAM domain annotation", "annotation"),
    PredicateVocab("kegg_annotated", "KEGG annotated", "Has KEGG ortholog annotation", "annotation"),
    PredicateVocab("cazy_annotated", "CAZy annotated", "Has CAZy family annotation", "annotation"),
    PredicateVocab("vog_annotated", "VOGdb annotated", "Has VOGdb ortholog annotation", "annotation"),
    PredicateVocab("hyddb_annotated", "HydDB annotated", "Has HydDB hydrogenase classification", "annotation"),
]

# ============================================================================
# COMPOSITION-BASED
# ============================================================================
COMPOSITION_PREDICATES = [
    PredicateVocab("gc_outlier", "GC outlier", "GC content >2 std from genome mean", "composition"),
    PredicateVocab("gc_high", "High GC", "GC content > 60%", "composition"),
    PredicateVocab("gc_low", "Low GC", "GC content < 40%", "composition"),
    PredicateVocab("low_complexity", "Low complexity", "High proportion of LC regions", "composition"),
]

# ============================================================================
# INFORMATION PROCESSING
# ============================================================================
INFO_PROCESSING_PREDICATES = [
    # Replication
    PredicateVocab("replication", "DNA replication", "Replication machinery", "info_processing"),
    PredicateVocab("dna_polymerase", "DNA polymerase", "DNA synthesis", "info_processing", "replication"),
    PredicateVocab("helicase", "Helicase", "DNA/RNA unwinding", "info_processing"),
    PredicateVocab("primase", "Primase", "RNA primer synthesis", "info_processing", "replication"),
    PredicateVocab("topoisomerase", "Topoisomerase", "DNA topology control", "info_processing", "replication"),
    PredicateVocab("ligase_dna", "DNA ligase", "DNA strand joining", "info_processing", "replication"),
    PredicateVocab("sliding_clamp", "Sliding clamp", "Processivity clamp (PCNA/beta)", "info_processing", "replication"),
    PredicateVocab("clamp_loader", "Clamp loader", "Clamp loader complex subunit", "info_processing", "replication"),

    # Transcription
    PredicateVocab("transcription", "Transcription", "RNA synthesis machinery", "info_processing"),
    PredicateVocab("rna_polymerase", "RNA polymerase", "RNA synthesis", "info_processing", "transcription"),
    PredicateVocab("transcription_elongation", "Transcription elongation", "Elongation factors", "info_processing", "transcription"),
    PredicateVocab("transcription_termination", "Transcription termination", "Termination factors", "info_processing", "transcription"),
    PredicateVocab("chromatin", "Chromatin", "DNA packaging/chromatin proteins", "info_processing", "dna_binding"),
    PredicateVocab("histone", "Histone", "Histone protein", "info_processing", "chromatin"),

    # Translation
    PredicateVocab("translation", "Translation", "Protein synthesis machinery", "info_processing"),
    PredicateVocab("ribosomal_protein", "Ribosomal protein", "Ribosome component", "info_processing", "translation"),
    PredicateVocab("trna_synthetase", "tRNA synthetase", "Aminoacyl-tRNA synthesis", "info_processing", "translation"),
    PredicateVocab("translation_factor", "Translation factor", "IF/EF/RF factor", "info_processing", "translation"),
    PredicateVocab("rrna_modification", "rRNA modification", "rRNA processing/modification", "info_processing", "translation"),
    PredicateVocab("trna_modification", "tRNA modification", "tRNA processing/modification", "info_processing", "translation"),

    # RNA processing
    PredicateVocab("rna_processing", "RNA processing", "RNA modification/degradation", "info_processing"),
    PredicateVocab("rnase", "RNase", "RNA degradation", "info_processing", "rna_processing"),
    PredicateVocab("rna_helicase", "RNA helicase", "RNA unwinding", "info_processing", "rna_processing"),
]

# ============================================================================
# CELL DIVISION
# ============================================================================
DIVISION_PREDICATES = [
    PredicateVocab("cell_division", "Cell division", "Division/septation machinery", "division"),
    PredicateVocab("ftsz", "FtsZ", "FtsZ division ring", "division", "cell_division"),
    PredicateVocab("divisome", "Divisome", "Division complex component", "division", "cell_division"),
    PredicateVocab("chromosome_partitioning", "Chromosome partitioning", "ParA/ParB system", "division"),
]

# ============================================================================
# TOPOLOGY PREDICATES (from TM helix prediction)
# ============================================================================
TOPOLOGY_PREDICATES = [
    # Transmembrane topology
    PredicateVocab("transmembrane_predicted", "TM predicted", "Has predicted transmembrane helices", "topology"),
    PredicateVocab("single_pass_membrane", "Single-pass TM", "Single transmembrane helix", "topology", "transmembrane_predicted"),
    PredicateVocab("multi_pass_membrane", "Multi-pass TM", "2+ transmembrane helices", "topology", "transmembrane_predicted"),
    PredicateVocab("polytopic_membrane", "Polytopic TM", "4+ transmembrane helices", "topology", "multi_pass_membrane"),

    # Topology orientation (N-terminus location)
    PredicateVocab("n_in_topology", "N-in topology", "N-terminus cytoplasmic", "topology"),
    PredicateVocab("n_out_topology", "N-out topology", "N-terminus extracellular/periplasmic", "topology"),

    # Soluble prediction
    PredicateVocab("soluble_predicted", "Soluble predicted", "No predicted TM helices", "topology"),
]

# ============================================================================
# ADDITIONAL PREDICATES (expanded vocabulary for comprehensive coverage)
# ============================================================================
ADDITIONAL_PREDICATES = [
    # Additional enzyme subtypes
    PredicateVocab("nitrogenase", "Nitrogenase", "Nitrogen fixation enzyme", "metabolism", "nitrogen_fixation"),
    PredicateVocab("nitrite_reductase", "Nitrite reductase", "Nitrite to NO/ammonia reduction", "metabolism", "denitrification"),
    PredicateVocab("nitric_oxide_reductase", "Nitric oxide reductase", "NO to N2O reduction", "metabolism", "denitrification"),
    PredicateVocab("nitrous_oxide_reductase", "Nitrous oxide reductase", "N2O to N2 reduction", "metabolism", "denitrification"),
    PredicateVocab("glutamine_synthetase", "Glutamine synthetase", "Glutamine synthesis", "metabolism", "ammonia_assimilation"),
    PredicateVocab("glutamine_amidotransferase", "Glutamine amidotransferase", "Glutamine-dependent amidation", "enzyme", "transferase"),
    PredicateVocab("urease", "Urease", "Urea hydrolysis", "metabolism", "urea_metabolism"),

    # Lyase/isomerase subtypes
    PredicateVocab("decarboxylase", "Decarboxylase", "CO2-removing enzyme", "enzyme", "lyase"),
    PredicateVocab("dehydratase", "Dehydratase", "Water-removing enzyme", "enzyme", "lyase"),
    PredicateVocab("epimerase", "Epimerase", "Stereoisomer interconversion", "enzyme", "isomerase"),
    PredicateVocab("mutase", "Mutase", "Intramolecular group transfer", "enzyme", "isomerase"),
    PredicateVocab("racemase", "Racemase", "Enantiomer interconversion", "enzyme", "isomerase"),

    # Hydrolase subtypes
    PredicateVocab("dehalogenase", "Dehalogenase", "Halogen removal enzyme", "enzyme", "hydrolase"),
    PredicateVocab("cysteine_protease", "Cysteine protease", "Cys-dependent peptidase", "enzyme", "protease"),

    # Transferase subtypes
    PredicateVocab("nucleotidyltransferase", "Nucleotidyltransferase", "Nucleotide transfer enzyme", "enzyme", "transferase"),
    PredicateVocab("sulfurtransferase", "Sulfurtransferase", "Sulfur transfer enzyme", "enzyme", "transferase"),
    PredicateVocab("phosphotransferase", "Phosphotransferase", "Phosphate transfer enzyme", "enzyme", "kinase"),

    # Oxidoreductase subtypes
    PredicateVocab("halogenase", "Halogenase", "Halogen incorporation enzyme", "enzyme", "oxidoreductase"),

    # Cofactor binding - additional
    PredicateVocab("sam_binding", "SAM binding", "S-adenosylmethionine binding", "binding", "cofactor_binding"),
    PredicateVocab("pqq_binding", "PQQ binding", "Pyrroloquinoline quinone binding", "binding", "cofactor_binding"),
    PredicateVocab("nadp_binding", "NADP binding", "NADP-binding domain", "binding", "nad_binding"),
    PredicateVocab("flavin_binding", "Flavin binding", "General flavin cofactor binding", "binding", "cofactor_binding"),

    # Metal-binding - additional
    PredicateVocab("oxygen_binding", "Oxygen binding", "O2 carrier/sensor protein", "binding", "metal_binding"),

    # Iron-sulfur cluster subtypes
    PredicateVocab("2fe2s", "2Fe-2S cluster", "2Fe-2S iron-sulfur cluster", "binding", "iron_sulfur"),

    # Radical enzymes
    PredicateVocab("radical_sam", "Radical SAM", "Radical SAM enzyme superfamily", "enzyme"),

    # Structure - additional
    PredicateVocab("beta_helix", "Beta helix", "Beta-helical structure", "structure"),
    PredicateVocab("heat_repeat", "HEAT repeat", "HEAT repeat domain", "structure", "repeat_domain"),
    PredicateVocab("helix_loop_helix", "Helix-loop-helix", "HLH DNA-binding motif", "binding", "dna_binding"),
    PredicateVocab("cbs_domain", "CBS domain", "Cystathionine beta-synthase domain", "binding", "nucleotide_binding"),
    PredicateVocab("pin_domain", "PIN domain", "PilT N-terminus RNase domain", "enzyme", "nuclease"),
    PredicateVocab("dockerin", "Dockerin", "Cellulosome dockerin domain", "structure"),
    PredicateVocab("cohesin", "Cohesin domain", "Cellulosome cohesin domain", "structure"),

    # Transport - additional
    PredicateVocab("sodium_transporter", "Sodium transporter", "Na+ transport protein", "transport", "ion_transporter"),
    PredicateVocab("signal_recognition", "Signal recognition", "Signal recognition particle", "transport", "secretion_system"),

    # Metabolism - additional pathways
    PredicateVocab("selenium_metabolism", "Selenium metabolism", "Selenoprotein/selenate metabolism", "metabolism"),
    PredicateVocab("thiosulfate", "Thiosulfate metabolism", "Thiosulfate oxidation/reduction", "metabolism", "sulfur_metabolism"),
    PredicateVocab("sugar_metabolism", "Sugar metabolism", "General sugar metabolism", "metabolism"),
    PredicateVocab("histidine_biosynthesis", "Histidine biosynthesis", "His biosynthesis pathway", "metabolism", "amino_acid_biosynthesis"),
    PredicateVocab("asparagine_synthesis", "Asparagine synthesis", "Asn biosynthesis", "metabolism", "amino_acid_biosynthesis"),

    # Cell wall - additional
    PredicateVocab("transpeptidase", "Transpeptidase", "Peptidoglycan crosslinking", "envelope", "peptidoglycan"),

    # Chelatase
    PredicateVocab("chelatase", "Chelatase", "Metal insertion into porphyrin", "enzyme"),

    # Miscellaneous enzyme/function terms
    PredicateVocab("enzyme", "Enzyme", "General enzymatic activity", "enzyme"),
    PredicateVocab("regulatory", "Regulatory domain", "Regulatory/allosteric domain", "regulation"),
    PredicateVocab("biosynthesis", "Biosynthesis", "General biosynthetic function", "metabolism"),
    PredicateVocab("polymerase", "Polymerase", "Nucleic acid polymerase", "info_processing"),
    PredicateVocab("dna_methylase", "DNA methylase", "DNA methylation enzyme", "info_processing", "dna_binding"),
    PredicateVocab("protein_modification", "Protein modification", "Post-translational modification", "enzyme"),
    PredicateVocab("ubiquitin_activation", "Ubiquitin activation", "E1-like ubiquitin activation", "enzyme", "protein_modification"),
    PredicateVocab("ubiquitin_ligase", "Ubiquitin ligase", "E3 ubiquitin ligase", "enzyme", "protein_modification"),
    PredicateVocab("deubiquitinase", "Deubiquitinase", "Ubiquitin-specific protease", "enzyme", "protein_modification"),
    PredicateVocab("metal_homeostasis", "Metal homeostasis", "Metal ion homeostasis", "stress"),
    PredicateVocab("immune_related", "Immune related", "Immune system component", "stress"),
    PredicateVocab("bioluminescence", "Bioluminescence", "Light production", "metabolism"),
    PredicateVocab("beta_oxidation", "Beta oxidation", "Fatty acid beta-oxidation", "metabolism", "fatty_acid_degradation"),
    PredicateVocab("vesicle_trafficking", "Vesicle trafficking", "Membrane vesicle transport", "transport"),
    PredicateVocab("binding", "Binding domain", "General binding function", "binding"),
]


# ============================================================================
# VIRAL PROTEINS (VOGdb-derived)
# ============================================================================
# Predicates specific to viral proteins, derived from VOGdb annotations.
# Some overlap with mobile element predicates but these are more specific.
VIRAL_PREDICATES = [
    # VOGdb functional categories
    PredicateVocab("viral_replication", "Viral replication", "Involved in viral genome replication", "viral"),
    PredicateVocab("viral_structure", "Viral structure", "Virion structural protein", "viral"),
    PredicateVocab("host_beneficial", "Host-beneficial", "Viral protein beneficial to host", "viral"),
    PredicateVocab("virus_beneficial", "Virus-beneficial", "Viral protein beneficial to virus", "viral"),

    # Capsid/Head
    PredicateVocab("viral_capsid", "Viral capsid", "Capsid/head protein", "viral", "viral_structure"),
    PredicateVocab("viral_scaffold", "Viral scaffold", "Capsid assembly scaffold", "viral", "viral_structure"),

    # Portal
    PredicateVocab("viral_portal", "Viral portal", "Portal vertex protein", "viral", "viral_structure"),

    # Tail components
    PredicateVocab("viral_tail", "Viral tail", "Tail structural protein", "viral", "viral_structure"),
    PredicateVocab("viral_tail_fiber", "Tail fiber", "Tail fiber/spike protein", "viral", "viral_tail"),
    PredicateVocab("viral_tape_measure", "Tape measure", "Tail length determination protein", "viral", "viral_tail"),
    PredicateVocab("viral_tail_sheath", "Tail sheath", "Contractile tail sheath", "viral", "viral_tail"),
    PredicateVocab("viral_baseplate", "Baseplate", "Baseplate structural protein", "viral", "viral_tail"),

    # DNA packaging
    PredicateVocab("viral_terminase", "Terminase", "DNA packaging terminase", "viral"),
    PredicateVocab("viral_packaging", "DNA packaging", "Involved in genome packaging", "viral"),
    PredicateVocab("viral_maturation_protease", "Maturation protease", "Head maturation protease", "viral"),

    # Lysis
    PredicateVocab("lysis", "Lysis", "Involved in host lysis", "viral"),
    PredicateVocab("endolysin", "Endolysin", "Peptidoglycan hydrolase for lysis", "viral", "lysis"),
    PredicateVocab("lysin", "Lysin", "Lytic enzyme", "viral", "lysis"),
    PredicateVocab("lysozyme", "Lysozyme", "Muramidase activity", "viral", "lysis"),
    PredicateVocab("spanin", "Spanin", "Outer membrane disruption", "viral", "lysis"),
    PredicateVocab("lysis_inhibitor", "Lysis inhibitor", "Inhibits lysis timing", "viral", "lysis"),

    # Lysogeny
    PredicateVocab("lysogeny", "Lysogeny", "Involved in lysogenic cycle", "viral"),
    PredicateVocab("excisionase", "Excisionase", "Prophage excision", "viral", "lysogeny"),
    PredicateVocab("repressor", "Repressor", "Transcriptional repressor", "viral"),
    PredicateVocab("anti_repressor", "Anti-repressor", "Repressor antagonist", "viral"),

    # Host interaction / Anti-defense
    PredicateVocab("anti_defense", "Anti-defense", "Counters host defense systems", "viral"),
    PredicateVocab("anti_crispr", "Anti-CRISPR", "CRISPR-Cas inhibitor", "viral", "anti_defense"),
    PredicateVocab("anti_restriction", "Anti-restriction", "Restriction enzyme inhibitor", "viral", "anti_defense"),
    PredicateVocab("host_attachment", "Host attachment", "Host cell attachment", "viral"),
    PredicateVocab("receptor_binding", "Receptor binding", "Host receptor binding", "viral", "host_attachment"),
    PredicateVocab("dna_injection", "DNA injection", "Genome injection into host", "viral"),

    # Replication-associated
    PredicateVocab("dna_polymerase", "DNA polymerase", "DNA replication polymerase", "viral", "viral_replication"),
    PredicateVocab("rna_polymerase", "RNA polymerase", "RNA polymerase", "viral", "viral_replication"),
    PredicateVocab("helicase", "Helicase", "DNA/RNA unwinding", "viral", "viral_replication"),
    PredicateVocab("primase", "Primase", "Primer synthesis", "viral", "viral_replication"),
    PredicateVocab("ssb_protein", "SSB protein", "Single-strand DNA binding", "viral", "viral_replication"),

    # Nucleases
    PredicateVocab("endonuclease", "Endonuclease", "Internal DNA/RNA cleavage", "viral"),
    PredicateVocab("exonuclease", "Exonuclease", "Terminal DNA/RNA degradation", "viral"),
    PredicateVocab("homing_endonuclease", "Homing endonuclease", "Selfish genetic element", "viral"),

    # DNA modification
    PredicateVocab("dna_modification", "DNA modification", "Modifies DNA bases", "viral"),
    PredicateVocab("dam_methylase", "Dam methylase", "GATC adenine methylation", "viral", "dna_modification"),
    PredicateVocab("dcm_methylase", "Dcm methylase", "CCWGG cytosine methylation", "viral", "dna_modification"),

    # Regulatory
    PredicateVocab("transcription_factor", "Transcription factor", "Regulates transcription", "viral"),
    PredicateVocab("sigma_factor", "Sigma factor", "RNA polymerase sigma subunit", "viral", "transcription_factor"),
    PredicateVocab("anti_sigma", "Anti-sigma factor", "Sigma factor inhibitor", "viral", "transcription_factor"),

    # Other enzymes
    PredicateVocab("deaminase", "Deaminase", "Removes amino groups", "viral"),
]


# ============================================================================
# COMBINED VOCABULARY
# ============================================================================
ALL_PREDICATES = (
    ENZYME_PREDICATES +
    TRANSPORT_PREDICATES +
    REGULATION_PREDICATES +
    METABOLISM_PREDICATES +
    CAZY_PREDICATES +
    STRUCTURE_PREDICATES +
    BINDING_PREDICATES +
    ENVELOPE_PREDICATES +
    MOBILE_PREDICATES +
    STRESS_PREDICATES +
    SIZE_PREDICATES +
    ANNOTATION_PREDICATES +
    COMPOSITION_PREDICATES +
    INFO_PROCESSING_PREDICATES +
    DIVISION_PREDICATES +
    TOPOLOGY_PREDICATES +
    VIRAL_PREDICATES +
    ADDITIONAL_PREDICATES
)

# Build lookup dictionaries
PREDICATE_BY_ID = {p.predicate_id: p for p in ALL_PREDICATES}
PREDICATES_BY_CATEGORY = {}
for p in ALL_PREDICATES:
    if p.category not in PREDICATES_BY_CATEGORY:
        PREDICATES_BY_CATEGORY[p.category] = []
    PREDICATES_BY_CATEGORY[p.category].append(p)


def get_predicate(predicate_id: str) -> Optional[PredicateVocab]:
    """Get a predicate definition by ID."""
    return PREDICATE_BY_ID.get(predicate_id)


def list_predicates(category: Optional[str] = None) -> list[PredicateVocab]:
    """List all predicates, optionally filtered by category."""
    if category:
        return PREDICATES_BY_CATEGORY.get(category, [])
    return ALL_PREDICATES


def list_categories() -> list[str]:
    """List all predicate categories."""
    return sorted(PREDICATES_BY_CATEGORY.keys())


def get_hierarchy(predicate_id: str) -> list[str]:
    """Get full hierarchy for a predicate (from most specific to root)."""
    result = [predicate_id]
    current = PREDICATE_BY_ID.get(predicate_id)
    while current and current.parent:
        result.append(current.parent)
        current = PREDICATE_BY_ID.get(current.parent)
    return result


__all__ = [
    "PredicateVocab",
    "ALL_PREDICATES",
    "PREDICATE_BY_ID",
    "PREDICATES_BY_CATEGORY",
    "get_predicate",
    "list_predicates",
    "list_categories",
    "get_hierarchy",
]
