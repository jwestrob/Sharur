"""
PFAM domain to predicate mappings.

Maps PFAM accessions to semantic predicates. Uses both:
1. Direct mappings for specific domains
2. Pattern-based matching on domain names/descriptions
"""

from pathlib import Path
import re

# ============================================================================
# DIRECT PFAM TO PREDICATE MAPPINGS
# ============================================================================
# Format: "PFAM_ACC": ["predicate1", "predicate2", ...]

PFAM_MAPPING_FILE = Path(__file__).with_name("data") / "pfam_predicates.tsv"
EXTERNAL_PFAM_MAPPING_CANDIDATES = [
    Path("data/reference/pfam_predicate_map.tsv"),
    Path(__file__).parent.parent.parent.parent / "data/reference/pfam_predicate_map.tsv",
]


def _load_pfam_mapping_file(path: Path) -> dict[str, list[str]]:
    """Load PFAM predicate mappings from TSV with duplicate validation.

    Columns:
        1. PFAM accession or PFAM short name
        2. Comma-separated predicate IDs
        3. Optional note
        4. Optional merge policy. Duplicate keys require ``merge`` here.
    """
    mappings: dict[str, list[str]] = {}
    if not path.exists():
        raise FileNotFoundError(path)

    with path.open() as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            parts = raw_line.rstrip("\n").split("\t")
            if len(parts) < 2:
                raise ValueError(f"{path}:{line_number}: expected at least 2 TSV columns")

            key = parts[0].strip()
            preds = [pred.strip() for pred in parts[1].split(",") if pred.strip()]
            merge_policy = parts[3].strip().lower() if len(parts) >= 4 else ""
            if not key or not preds:
                raise ValueError(f"{path}:{line_number}: empty key or predicate list")

            if key in mappings and merge_policy != "merge":
                raise ValueError(
                    f"{path}:{line_number}: duplicate PFAM mapping for {key!r}; "
                    "set merge_policy=merge to make this explicit"
                )

            merged = set(mappings.get(key, [])) | set(preds)
            mappings[key] = sorted(merged)

    return mappings


def _load_packaged_pfam_mappings() -> dict[str, list[str]]:
    """Load curated PFAM mappings shipped with Sharur."""
    return _load_pfam_mapping_file(PFAM_MAPPING_FILE)


def _load_external_pfam_mappings() -> dict[str, list[str]]:
    """Load optional local PFAM predicate mapping extensions."""
    path = next((p for p in EXTERNAL_PFAM_MAPPING_CANDIDATES if p.exists()), None)
    if path is None:
        return {}
    return _load_pfam_mapping_file(path)


def _merge_pfam_mapping_layers() -> dict[str, list[str]]:
    """Merge packaged mappings with optional local extensions."""
    merged = _load_packaged_pfam_mappings()
    external = _load_external_pfam_mappings()
    for key, preds in external.items():
        merged[key] = sorted(set(merged.get(key, [])) | set(preds))
    return merged


PFAM_TO_PREDICATES = _merge_pfam_mapping_layers()


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
    (r"helix.turn.helix|helix.hairpin.helix|HTH|HhH", ["dna_binding", "helix_turn_helix"]),
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
