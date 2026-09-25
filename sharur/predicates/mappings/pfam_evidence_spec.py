"""Evidence definitions for Pfam-derived predicates.

For each predicate: ``go`` lists GO terms whose presence in a family's InterPro
GO annotation (including ancestors) supports it; ``text`` lists patterns that
must match the family's own Pfam name or description. Pattern conventions are
enforced by :func:`sharur.predicates.mappings.pfam_evidence.compile_evidence`:
lowercase alternatives are words (case-insensitive), alternatives containing
capitals are symbols (case-sensitive), every alternative starts at a word
boundary, and alternatives of four or fewer letters also end at one. A match
immediately followed by "-like" is a resemblance and never evidence.

A pattern states the predicate in Pfam's words or names the protein itself
(e.g. ``MutS``); domain-type tokens associated with a function by convention
(PIN, CBS, PAS, TIR, ...) are not evidence of that function.
"""
# Families whose Pfam name coincides with an unrelated protein's symbol.
# (family name, predicate) -> reason; these pairs never count as supported.
HOMONYMS = {
    ("Cas1_AcylT", "defense_system"): "fungal Cas1p capsule O-acetyltransferase, unrelated to CRISPR Cas1",
    ("Cas1_AcylT", "crispr_associated"): "fungal Cas1p capsule O-acetyltransferase, unrelated to CRISPR Cas1",
    ("Cas1_AcylT", "cas_domain"): "fungal Cas1p capsule O-acetyltransferase, unrelated to CRISPR Cas1",
}

E = {}
def d(pred, go=(), text=()):
    entry = E.setdefault(pred, {"go": [], "text": []})
    entry["go"] += [g for g in go if g not in entry["go"]]
    entry["text"] += [x for x in text if x not in entry["text"]]

# ---- Broad enzyme classes -------------------------------------------------
d("hydrolase", ["GO:0016787"], [r"hydrolase", r"\bpeptidase", r"protease", r"esterase", r"lipase", r"phosphatase", r"nuclease", r"glycosyl hydrolase", r"glycoside hydrolase", r"glyco_hydro", r"amidase", r"deacetylase", r"deaminase", r"lactamase", r"chitinase", r"lysozyme", r"diesterase", r"\bNUDIX\b", r"pyrophosphohydrolase", r"\bhydrolys", r"glucosidase|galactosidase|mannosidase|xylosidase|amylase|cellulase|glucanase|fucosidase|sialidase|hexosaminidase|arabinofuranosidase"])
d("transferase", ["GO:0016740"], [r"transferase", r"\bkinase", r"methylase", r"polymerase", r"transaminase", r"transglycosylase", r"\bMTase", r"Methyltransf", r"Glyco_trans", r"Glycos_trans", r"Acetyltransf", r"Acyl_transf", r"transpeptidase"])
d("oxidoreductase", ["GO:0016491"], [r"oxidoreductase", r"dehydrogenase", r"reductase", r"oxidase", r"oxygenase", r"peroxidase", r"catalase", r"hydroxylase", r"desaturase", r"thioredoxin", r"glutaredoxin", r"peroxiredoxin", r"dismutase", r"Oxidored", r"_dh\b", r"_DH\b", r"Pyr_redox"])
d("lyase", ["GO:0016829"], [r"lyase", r"aldolase", r"decarboxylase", r"dehydratase", r"hydratase", r"enolase", r"carboxy-?lyase"])
d("isomerase", ["GO:0016853"], [r"isomerase", r"epimerase", r"racemase", r"mutase", r"tautomerase", r"cyclase", r"topoisomerase", r"PPIase|peptidyl-prolyl"])
d("ligase", ["GO:0016874"], [r"ligase", r"synthetase", r"tRNA-synt"])
d("enzyme", ["GO:0003824"], [r"ase\b", r"enzyme"])
d("kinase", ["GO:0016301"], [r"kinase", r"Pkinase", r"PfkB", r"phosphotransferase"])
d("phosphatase", ["GO:0016791", "GO:0004721"], [r"phosphatase", r"\bPP2C", r"PAP2", r"phosphohydrolase", r"\bPTP", r"DSPc"])
d("methyltransferase", ["GO:0008168"], [r"methyl ?transferase", r"methylase", r"Methyltransf", r"\bMTase", r"MethylTransf", r"SAM-dependent", r"\bRsm", r"\bTrm"])
d("glycosyltransferase", ["GO:0016757"], [r"glyco?syl ?transferase", r"Glyco_trans", r"Glycos_trans", r"Glyco_tranf", r"transglycosylase", r"mannosyltransferase|glucosyltransferase|galactosyltransferase|rhamnosyltransferase"])
d("acetyltransferase", ["GO:0016407", "GO:0008080"], [r"acetyl ?transferase", r"Acetyltransf", r"GNAT", r"hexapeptide"])
d("nucleotidyltransferase", ["GO:0016779"], [r"nucleotidyl ?transferase", r"NTP_transf", r"adenylyltransferase", r"guanylyltransferase", r"uridylyltransferase", r"polymerase", r"cytidylyltransferase", r"PolyA_pol", r"NT_", r"MobA"])
d("phosphotransferase", ["GO:0016772"], [r"phosphotransferase", r"kinase"])
d("sulfurtransferase", ["GO:0016783"], [r"sulfurtransferase", r"rhodanese", r"\bThiI\b|ThiF|MoeB"])
d("aminotransferase", ["GO:0008483"], [r"aminotransferase", r"Aminotran", r"transaminase", r"amidotransferase", r"GATase"])
d("glutamine_amidotransferase", ["GO:0004359", "GO:0006541"], [r"glutamine amidotransferase", r"GATase", r"amidotransferase"])
d("synthase", [], [r"synthase"])
d("synthetase", [], [r"synthetase", r"tRNA-synt"])
d("dehydrogenase", [], [r"dehydrogenase", r"_dh\b", r"_DH\b", r"dehydrog", r"Ldh_1|Gp_dh|adh_|ADH_|Aldedh|IlvN|KARI|2-Hacid_dh|NAD_binding|Shikimate_DH|G6PD|6PGD|Sacchrp_dh|GFO_IDH|ELFV_dehydrog|Glu_dehyd|IDH|Iso_dh|DAO|FAD_binding|Sulfite_dh"])
d("reductase", [], [r"reductase", r"Oxidored", r"reductoisomerase", r"Pyr_redox"])
d("oxidase", ["GO:0016491"], [r"oxidase"])
d("oxygenase", ["GO:0004497", "GO:0051213"], [r"oxygenase", r"hydroxylase"])
d("monooxygenase", ["GO:0004497"], [r"monooxygenase", r"hydroxylase", r"p450"])
d("dioxygenase", ["GO:0051213"], [r"dioxygenase"])
d("peroxidase", ["GO:0004601"], [r"peroxidase", r"peroxiredoxin", r"AhpC", r"catalase"])
d("peroxiredoxin", ["GO:0051920"], [r"peroxiredoxin", r"AhpC", r"\bPrx", r"Redoxin"])
d("catalase", ["GO:0004096"], [r"catalase"])
d("superoxide_dismutase", ["GO:0004784"], [r"superoxide dismutase", r"Sod_"])
d("decarboxylase", ["GO:0016831"], [r"decarboxylase", r"carboxy-?lyase"])
d("dehydratase", ["GO:0016836"], [r"dehydratase", r"hydro-?lyase", r"hydratase"])
d("epimerase", ["GO:0016857", "GO:0016854"], [r"epimerase", r"racemase"])
d("racemase", ["GO:0036361", "GO:0016855", "GO:0016854"], [r"racemase", r"epimerase"])
d("mutase", ["GO:0016866"], [r"mutase", r"PGM", r"phosphoglucomutase|phosphomannomutase"])
d("carboxylase", ["GO:0016885", "GO:0004075", "GO:0008964", "GO:0016984"], [r"carboxylase", r"carboxyltransferase"])
d("deaminase", ["GO:0019239"], [r"deaminase", r"deiminase", r"dCMP_cyt_deam"])
d("amidase", ["GO:0016811", "GO:0016810"], [r"amidase", r"amidohydrolase", r"Amidohydro", r"deacylase", r"desuccinylase|aspartoacylase", r"acylase", r"Creatininase"])
d("esterase", ["GO:0016788", "GO:0052689"], [r"esterase", r"lipase", r"Abhydrolase", r"hydrolase_4|Hydrolase_4", r"thioesterase", r"diesterase", r"Esterase"])
d("lipase", ["GO:0004806", "GO:0016042", "GO:0004620"], [r"lipase", r"phospholipase"])
d("amylase", ["GO:0016160", "GO:0004556"], [r"amylase"])
d("cellulase", ["GO:0008810", "GO:0030245"], [r"cellulase", r"endoglucanase", r"cellulose"])
d("chitinase", ["GO:0004568", "GO:0006032"], [r"chitinase"])
d("lysozyme", ["GO:0003796"], [r"lysozyme", r"Glyco_hydro_25"])
d("glycosidase", ["GO:0016798", "GO:0004553"], [r"glycosi?d(ase|e hydrolase)", r"glycosyl hydrolase", r"Glyco_hydro", r"glucosidase|galactosidase|mannosidase|xylosidase|amylase|cellulase|glucanase|fucosidase|sialidase|hexosaminidase|arabinofuranosidase|chitinase|lysozyme|glucuronidase"])
d("phosphodiesterase", ["GO:0008081", "GO:0004114", "GO:0071111"], [r"phosphodiesterase", r"diesterase", r"\bEAL\b", r"HD-GYP", r"HDOD", r"PDEase", r"GDPD", r"Metallophos", r"DHH"])
d("endonuclease", ["GO:0004519"], [r"endonuclease", r"endonuc", r"\bHNH", r"GIY-YIG", r"restriction", r"\bRuvC", r"\bXPG"])
d("exonuclease", ["GO:0004527"], [r"exonuclease", r"exonuc", r"\bDnaQ|RNase_T|DEDDh|RecJ"])
d("rnase", ["GO:0004540"], [r"\bRNase", r"ribonuclease", r"RNase_PH|YbeY|Ribonuc_"])
d("dnase", ["GO:0004536"], [r"\bDNase", r"deoxyribonuclease", r"restriction endonuclease"])
d("helicase", ["GO:0004386"], [r"helicase", r"\bDEAD", r"Helicase_C", r"ResIII", r"UvrD", r"SNF2", r"\bRecQ", r"DnaB|MCM", r"\bRecG"])
d("topoisomerase", ["GO:0003916", "GO:0006265"], [r"topoisomerase", r"gyrase", r"Topo", r"Toprim", r"DNA_topoisoIV"])
d("protease", ["GO:0008233", "GO:0006508"], [r"protease", r"peptidase", r"proteinase", r"Peptidase_", r"\bClp", r"\bLon", r"signal peptidase", r"\bCAAX", r"Rhomboid", r"\bFtsH", r"\bHtrA|DegP|Trypsin"])
d("peptidase", ["GO:0008233"], [r"peptidase", r"protease", r"Peptidase_", r"proteinase", r"Rhomboid", r"\bJAB|MPN"])
d("cysteine_protease", ["GO:0008234"], [r"cysteine (protease|peptidase)", r"C1[0-9]?\b.*peptidase", r"Peptidase_C", r"caspase", r"\bUCH|OTU|Josephin|USP"])
d("deubiquitinase", ["GO:0101005", "GO:0004843", "GO:0016579"], [r"deubiquitin", r"ubiquitin.*(hydrolase|protease|specific)", r"\bUCH|OTU|JAB|MPN|USP"])
d("clp_protease", ["GO:0009368"], [r"\bClp"])
d("aaa_domain", [], [r"\bAAA\b|AAA_|AAA\+"])
d("p_loop", [], [r"P-?loop", r"NTPase"])
d("gtpase", ["GO:0003924"], [r"GTPase", r"GTP_EFTU", r"MMR_HSR1", r"\bRas\b|\bRoc\b|\bArf\b|Gtr1_RagA|\bEra\b|\bFeoB|Dynamin|Septin|Tubulin|FtsZ|Obg|EngA|IIGP"])
d("gtp_binding", ["GO:0005525", "GO:0003924", "GO:0019001"], [r"GTP", r"GTPase", r"\bRas\b|\bRoc\b|\bArf\b|Gtr1_RagA|MMR_HSR1|Tubulin|FtsZ|Septin|Dynamin|FeoB|Obg|EngA|IIGP|Gtr1"])
d("atp_binding", ["GO:0005524", "GO:0016887"], [r"\bATP", r"ATPase", r"\bAAA", r"ABC_tran", r"kinase", r"helicase", r"\bP-?loop", r"synthetase", r"ligase", r"\bDEAD", r"\bSMC", r"chaperon", r"HSP70|HSP90|Cpn60|GroEL|DnaK|HtpG|HATPase", r"ResIII|UvrD|SNF2|MCM|RecA|DnaB|FtsK|SecA|KaiC|ArsA|MipZ|ParA|CbiA|MinD|Mg_chelatase|Sigma54_activat"])
d("nucleotide_binding", ["GO:0000166"], [r"nucleotide[- ]binding", r"NTP", r"\bATP|GTP"])
d("amp_binding", [], [r"AMP-binding"])

# ---- Cofactors / metals ---------------------------------------------------
d("metal_binding", ["GO:0046872"], [r"metal", r"zinc|\bZn|iron|\bFe\b|\bFe-|copper|\bCu|nickel|\bNi|cobalt|manganese|magnesium|calcium|molybd|tungsten|heme|haem|cytochrome|ferritin|rubredoxin|zinc[- ]?finger|zf-|metallo|\bRING|ferric|siderophore|chelat|Fer[24]|Rieske|hemerythrin|cupin|cupredoxin|Cu_bind|Cu-oxidase|Radical_SAM|4Fe-4S|2Fe-2S|FeS|Fe-S|SufE|NifU|HypA|Ni_insertion|CbiX|B12", r"Sod_"])
d("zinc_binding", ["GO:0008270"], [r"zinc|\bZn|zf-|zinc[- ]?finger|\bRING|ADH_zinc|Zn_", r"ribbon"])
d("zinc_finger", ["GO:0008270"], [r"zinc[- ]?finger|zf-|\bRING\b|Zn_ribbon|Zn-ribbon|Zn_Ribbon"])
d("iron_binding", ["GO:0005506", "GO:0008198"], [r"iron|\bFe\b|\bFe-|ferric|ferrous|ferritin|heme|haem|hemerythrin|rubredoxin|siderophore|Fe-ADH|Fe_", r"FeoB|Dps|\bDPS\b"])
d("2fe2s", ["GO:0051537"], [r"2Fe-2S|Fer2|Rieske"])
d("heme_binding", ["GO:0020037", "GO:0046906"], [r"heme|haem|cytochrome|\bCytochrom|globin|catalase|peroxidase|p450|\bCyt_|Cytochrome_C|Cytochrom_C|Cytochrom_B|COX|HemS|Haem"])
d("copper_binding", ["GO:0005507"], [r"copper|\bCu\b|\bCu_|Cu-|cupredoxin|COX2|plastocyanin|azurin|CopC|CopD|CueO|NosD|NosL|multicopper"])
d("manganese_binding", ["GO:0030145"], [r"manganese|\bMn\b|\bMn-|Sod_Fe|superoxide dismutase|arginase|MntH"])
d("magnesium_binding", ["GO:0000287"], [r"magnesium|\bMg\b|\bMg-|Mg_chelatase|CorA|MgtE|enolase|RuBisCO"])
d("calcium_binding", ["GO:0005509"], [r"calcium|\bCa\b|\bCa-|\bCa2|EF-hand|EF_hand|cadherin|calx|HemolysinCabind|RTX|S_100|Calx-beta|GDT1|dockerin"])
d("molybdenum_binding", ["GO:0030151", "GO:0043546", "GO:0030151"], [r"molybd|\bMo\b|\bMo-|MoaA|MoaB|MoaC|MoaD|MoaE|MoeA|MoeB|MobA|MobB|MogA|Molybdop|MOSC|Mo-co|MoCF|tungsten|\bW\b-"])
d("molybdenum_cofactor", ["GO:0032324", "GO:0043546"], [r"molybd"])
d("fad_binding", ["GO:0050660", "GO:0071949"], [r"\bFAD|flavin adenine|FAD_binding|FAD_oxidored|flavoprotein"])
d("fmn_binding", ["GO:0010181"], [r"\bFMN|flavin mononucleotide|Flavodoxin|flavodoxin|Oxidored_FMN|FMN_red|Nitroreductase|Flavin_Reduct|FMN_dh|FMN_bind"])
d("flavin_binding", ["GO:0050660", "GO:0010181"], [r"flavin|\bFAD|\bFMN|flavodoxin|flavoprotein|Flavoprotein|ETF"])
d("pqq_binding", ["GO:0070968"], [r"\bPQQ|pyrroloquinoline"])
d("plp_binding", ["GO:0030170"], [r"pyridoxal|\bPLP|PALP|Pyridoxal"])
d("thiamine_binding", ["GO:0030976", "GO:0030975"], [r"thiamin|\bTPP|TPP_enzyme|thiamine pyrophosphate"])
d("biotin_binding", ["GO:0009374"], [r"biotin"])
d("sam_binding", ["GO:1904047"], [r"S-adenosyl|\bSAM[- ]binding|\bSAM[- ]dependent|Radical_SAM|radical SAM|SAM-dependent"])
d("radical_sam", ["GO:0051539"], [r"radical ?SAM|Radical_SAM|SPASM|\bRS_|DUF4953"])
d("coenzyme_a_binding", ["GO:0120225"], [r"\bCoA\b|coenzyme A"])
d("f420_dependent", [], [r"F420|F_420|coenzyme F420"])
d("oxygen_binding", ["GO:0019825"], [r"oxygen[- ]binding|globin|hemerythrin|hemocyanin|NiFe_hyd_SSU"])
d("cytochrome", [], [r"cytochrom|apocytochrom|Apocytochr|\bCyt_|Cytochrom_|COX|\bcyt\b"])
d("thioredoxin", ["GO:0015035"], [r"thioredoxin|Thioredoxin|\bTrx|Redoxin|DsbA|DsbB|DsbD|Glutaredoxin|glutaredoxin|\bTlpA|AhpC|SCO1|ResA"])
d("glutaredoxin", ["GO:0015038"], [r"glutaredoxin"])
d("redox", ["GO:0016491"], [r"redox|thioredoxin|oxidoreductase"])

# ---- Transport / membrane -------------------------------------------------
d("transporter", ["GO:0005215", "GO:0006810", "GO:0055085", "GO:0022857"], [r"transport|transporter|permease|channel|porin|symporter|antiporter|exporter|importer|efflux|carrier|uptake|\bABC\b|ABC_|ABC-|ABC2|SBP_bac|BPD_transp|Peripla_BP|MFS|Sugar_tr|AA_permease|Aa_trans|OPT|Na_H|Cation_ATPase|E1-E2|TonB|OmpA|pump|translocase|translocator|flippase|MatE|EamA|DMT|TrkH|TrkA|K_trans|Ion_trans|MscS|MscL|CorA|Mg_trans|MgtE|ZIP|Nramp|FTR1|FeoB|FeoA|CbiM|CbiQ|CbiN|NikA|OppC|oligo_HPY|DctM|DctQ|DctP|TRAP|SecY|SecE|SecG|SecA|TatA|TatC|Tat_|YidC|MacB|FtsX|HlyD|OEP|TolC|Porin|LamB|Secretin|GSPII|T2SS|PTS|EIIC|Na_Ca_ex|CPA|LysE|ArsB|ACR_tran|MMPL|Mem_trans|MotA|ExbB|TolQ|VIT|CDF|Cation_efflux|Bestrophin|Aquaporin|MIP|Ammonium_transp|Sulfate_transp|SLC|PhoU|Pi_transp|PHO|NCS|Xan_ur_permease|Nucleoside_tran|BCCT|SSS|SNF|Na_sulph|TauE|TSUP|DUF81|CitMHS|ArsA|Anion_ATPase|ABC_membrane|ABC_tran"])
d("abc_transporter", ["GO:0140359", "GO:0043190", "GO:0042626"], [r"\bABC\b|ABC_|ABC-|ABC2|ATP-binding cassette|SBP_bac|BPD_transp|Peripla_BP|binding-protein-dependent|OppC|oligo_HPY|NikA|MacB|FtsX|CbiQ|CbiM|CbiN|ECF|TOBE|MlaD|MlaE|NMT1|Phosphonate-bd|Lipoprotein_9|ABC_membrane|ABC_tran|SBP_bac_[0-9]|Peripla_BP_[0-9]|OpuAC|Mod_|ModA|PstS|PBP_like|Lipoprotein_17"])
d("abc_atpase", ["GO:0016887", "GO:0140359"], [r"ABC_tran|ABC transporter$|ABC.*ATP|ATP-binding cassette|ABC_ATPase"])
d("abc_permease", ["GO:0055085"], [r"ABC.*(membrane|permease|transmembrane)|BPD_transp|binding-protein-dependent|OppC|FtsX|MacB|CbiQ|ABC2_membrane|ABC_membrane|MlaE|permease"])
d("abc_substrate_binding", [], [r"solute-binding|substrate-binding|SBP_bac|Peripla_BP|periplasmic binding|extracellular solute|binding protein|NMT1|Lipoprotein_9|OpuAC|ModA|PstS|NikA"])
d("mfs_transporter", [], [r"\bMFS|major facilitator|Sugar_tr|MFS_"])
d("symporter", ["GO:0015293"], [r"symporter|Na\+.*(solute|coupled)|SSS|SNF|BCCT|Sodium|Na_sulph|DctA|SDF|Na_Ala_symp|NhaD|Glt_symporter|DAACS|SLC5|SLC6|Proton-dependent|OPT|PTR2|Sugar_tr|MFS"])
d("antiporter", ["GO:0015297"], [r"antiporter|exchanger|Na_H|Na\+/H\+|NhaA|NhaB|CPA|Mrp|Mnh|MNHE|MnhB|PhaG|Proton_antipo|Na_Ca_ex|MatE|GDT1|Sulfate_transp"])
d("ion_channel", ["GO:0005216", "GO:0015267"], [r"channel|\bIon_trans|MscS|MscL|Kch|CLC|Voltage_CLC|Bestrophin|Aquaporin|MIP|porin|Porin|Lig_chan|Mg_trans_NIPA|IRK|Glutamate receptor|ANF_receptor|Cation_efflux"])
d("ion_transporter", ["GO:0015075", "GO:0006811", "GO:0034220"], [r"\bion\b|cation|anion|sodium|potassium|\bNa\b|\bNa\+|\bK\+|proton|H\+|Na_|K_|CorA|Mg_|MgtE|ZIP|Nramp|Ion_trans|TrkH|TrkA|Kch|antiporter|symporter|Cation|Anion|chloride|Voltage_CLC|NhaA|Mrp|Mnh"])
d("metal_transporter", ["GO:0000041", "GO:0030001", "GO:0046873"], [r"metal|cobalt|nickel|zinc|\bZn|iron|ferric|ferrous|manganese|magnesium|copper|cadmium|mercury|Cbi[MNQO]|NikA|CorA|Mg_trans|MgtE|ZIP|Zip|Nramp|NRAMP|FTR1|FeoA|FeoB|MntH|CDF|Cation_efflux|Cation_ATPase|E1-E2|HMA|CopC|CopD|CusB|CzcD|ZupT|VIT|TonB|CbiM|CbiQ|siderophore|Heavy"])
d("sodium_transporter", ["GO:0015081", "GO:0006814"], [r"sodium|\bNa\b|\bNa\+|Na_|Mrp|Mnh|MNHE|NhaA|NQR|Rnf|OAD"])
d("phosphate_transporter", ["GO:0006817", "GO:0005315"], [r"phosphate transport|Pi_transp|PhoU|PHO4|PstS|PstA|PstC|Pst[ABC]|Phosphonate|phosphonate|PBP_like|Na_Pi|SPX"])
d("sulfate_transporter", ["GO:0008272", "GO:0015116"], [r"sulfate|sulphate|thiosulfate|Sulfate_transp|SulP|CysA|CysT|CysW|Sbp|TauE"])
d("sugar_transporter", ["GO:0008643", "GO:0051119"], [r"sugar|carbohydrate|glucose|maltose|ribose|xylose|arabinose|galactose|fructose|lactose|sucrose|Sugar_tr|SBP_bac_1|SBP_bac_8|Peripla_BP_1|Peripla_BP_4|MalF|MalG|PTS|EIIC|EIIB|EIIA|MFS_2|SWEET|TRAP"])
d("amino_acid_transporter", ["GO:0006865", "GO:0015171"], [r"amino[- ]acid.*(transport|permease)|AA_permease|Aa_trans|LysE|RhtB|BCCT|DAACS|Glt_symporter|SBP_bac_3|BPD_transp_2|Branched-chain amino acid|LIV|Aminoacid|AzlC|AzlD|ThrE|Trp_Tyr_perm|Gly_transporter"])
d("peptide_transporter", ["GO:0015833", "GO:0042938"], [r"peptide.*transport|oligopeptide|dipeptide|OppC|oligo_HPY|SBP_bac_5|OPT|PTR2|Peptide_transp|NikA|DppC"])
d("porin", ["GO:0015288"], [r"porin|Porin|LamB|OmpA|Omp85|OprB|OmpW|OmpH|Omp_"])
d("outer_membrane", ["GO:0019867", "GO:0009279"], [r"outer membrane|OmpA|Omp|porin|Porin|LamB|TonB_dep|TonB_dep_Rec|OEP|TolC|Secretin|BamA|Omp85|LptD|LptE|POTRA|Lipoprotein|Slam|Autotransporter|ShlB|FhaC|Usher|PapC|CsgG|Wza|PulD|PilQ|OMP|LolB|OprF|OmpH|YaeT|Bac_surface_Ag|LPS_assembly"])
d("tonb_dependent_receptor", ["GO:0015344", "GO:0038023", "GO:0015891"], [r"TonB|Plug"])
d("efflux_pump", ["GO:0042910", "GO:0015562", "GO:0046618"], [r"efflux|exporter|ACR_tran|AcrB|MatE|MATE|EmrE|SMR|Multi_Drug_Res|HlyD|TolC|OEP|MacB|MFS_1|NorM|CusA|CzcA|RND|TetR"])
d("multidrug_resistance", ["GO:0042910", "GO:0015562"], [r"multidrug|multi-drug|Multi_Drug|ACR_tran|AcrB|MatE|MATE|EmrE|SMR|NorM|drug"])
d("uniporter", [], [r"uniporter"])
d("periplasmic", ["GO:0042597"], [r"periplasm|Peripla"])
d("inner_membrane", ["GO:0005886"], [r"inner membrane|cytoplasmic membrane|plasma membrane"])
d("s_layer", ["GO:0030115"], [r"S-layer|S_layer|SLH|surface layer"])
d("adhesin", ["GO:0007155", "GO:0098609"], [r"adhesin|adhesion|Cadherin|cadherin|invasin|intimin|Hemagglutinin|haemagglutin|YadA|Hep_Hag|fibronectin|fn[123]|Cna|SdrG|SdrD|MucBP|Big_|Dockerin|Cohesin|PKD|FIVAR|Rib|collagen|Collagen|lectin|Lectin|fimbri|pilin"])
d("dockerin", [], [r"dockerin"])
d("cohesin", [], [r"cohesin"])
d("beta_barrel", ["GO:0015288", "GO:0019867"], [r"beta[- ]barrel|b-barrel|porin|Porin|Lipocalin|OmpA|Omp|TonB_dep|LamB|barrel"])
d("secretion_component", ["GO:0009306", "GO:0015031", "GO:0030254", "GO:0030255", "GO:0015628", "GO:0033103", "GO:0043952", "GO:0043953"], [r"secretion|secretory|secretin|Secretin|export|Sec[ABDEFGY]|SecD|SecF|SecA|SecY|Sec61|SEC-C|YajC|SRP|signal recognition|Tat[A-C]|TatA|TatC|TAT|twin-arginine|T[1-9]SS|GSP|T2SS|T3SS|T4SS|T6SS|type (I|II|III|IV|V|VI|VII|IX)|VirB|VirD|Trb|Tra[A-Z]|Hcp|VgrG|ImpA|ImpB|ImpC|ImpG|ImpJ|VasE|Vas|Tss|FliP|FliQ|FliR|FlhA|FlhB|YscN|InvA|HlyD|HlyB|TolC|OEP|PrgI|Por_Secre|T9SS|Usher|PapC|chaperone-usher|FhaC|ShlB|Autotransporter|Class_IIIsignal|signal peptide|Sortase|translocase|translocon|TRAP_beta|EMC|PilQ|PulD|pilus|Pilin|Flp"])
d("t1ss_component", ["GO:0030253"], [r"type I secretion|T1SS|HlyD|HlyB|TolC|OEP|RTX|HemolysinCabind|PrtD"])
d("t2ss_component", ["GO:0015628", "GO:0015627"], [r"type II secretion|T2SS|T2SSE|T2SSF|T2SSG|GSP|GspD|Secretin|PulD|general secretion"])
d("t3ss_component", ["GO:0030254", "GO:0030257"], [r"type III secretion|T3SS|YscN|YscC|YscV|InvA|SctN|SctV|SctC|PrgI|FliP|FliQ|FliR|FlhA|FlhB|EscN|EscV|HrcN|HrcV|Flagellar export"])
d("t4ss_component", ["GO:0030255", "GO:0043684"], [r"type IV secretion|T4SS|VirB|VirD4|Trb|TrbI|TraC|TraD|TraG|TraK|TraL|TraE|conjugal|conjugation|T4SS_pilin|TrwB"])
d("t6ss_component", ["GO:0033103", "GO:0033104"], [r"type VI secretion|T6SS|Hcp|VgrG|ImpA|ImpB|ImpC|ImpG|ImpJ|VasE|VasA|VasB|VasK|TssB|TssC|TssK|TssL|TssM|IcmF|PAAR|FHA_2|T6SS_"])
d("sec_pathway", ["GO:0043952", "GO:0065002", "GO:0006605", "GO:0031522"], [r"\bSec[ABDEFGY]|SecD|SecF|SecA|SecY|SecE|SecG|Sec61|SEC-C|Sec_GG|YajC|preprotein translocase|protein export|SRP|signal recognition|signal peptidase|Class_IIIsignal|TRAP_beta|translocon"])
d("tat_pathway", ["GO:0043953", "GO:0033281"], [r"\bTat|twin[- ]arginine|TAT_signal|MttA|Hcf106"])
d("secretion_system", ["GO:0009306"], [r"secretion"])
d("conjugation", ["GO:0000746", "GO:0009291"], [r"conjugat|Tra[A-Z]|Trb|T4SS|VirB|Mob|relaxase|TrwB|TraG|TrbI|T4SS_pilin|Tcp|TcpE|Tfp"])
d("relaxase", [], [r"relaxase|Mob"])

# ---- Motility / pili ------------------------------------------------------
d("flagellum", ["GO:0009288", "GO:0071973", "GO:0097588", "GO:0071978"], [r"flagell|Flg|Fli[A-Z]|Flh|Mot[AB]|FlaG|FlaI|FlaJ|FlaC|archaell|Arch_flagellin"])
d("flagellin", ["GO:0005198"], [r"flagellin|Flagellin|Arch_flagellin|FliC"])
d("flagellar_motor", ["GO:0009288", "GO:0071973"], [r"flagellar motor|MotA|MotB|FliG|FliM|FliN|FliY|OmpA"])
d("flagellar_hook", [], [r"hook|FlgE|FlgK|FlgL|Flg_hook"])
d("pilus", ["GO:0009289", "GO:0044096"], [r"pilus|pili\b|pilin|Pilin|fimbri|Pil[A-Z]|T4SS_pilin|Flp|Tad|CpaB|Usher|PapC|type IV pil|N_methyl|GspG|PulG"])
d("type_iv_pilus", ["GO:0044096", "GO:0043683"], [r"type IV pil|type 4 pil|Pil[A-Z]|N_methyl|T4P|Tfp|pilin|PilT|PilB|PilC|PilM|PilN|PilO|PilQ|T2SSE|T2SSF|GSPII"])
d("motility", ["GO:0048870", "GO:0071973", "GO:0097588"], [r"motil|flagell|archaell|pil"])
d("chemotaxis", ["GO:0006935"], [r"chemotaxis|Che[ABDRWYZ]|MCP|methyl-accepting|CheB|CheW|CheR|CheY|CheA|CheZ|CheC|CheD|HAMP"])

# ---- Nucleic acids / info processing -------------------------------------
d("rna_binding", ["GO:0003723", "GO:0019843", "GO:0000049", "GO:0003729"], [r"RNA[- ]binding|RNA_bind|RRM|KH|S1|\bPUA|THUMP|TRAM|\bS4|ribosom|Ribosomal|tRNA|rRNA|\bRNase|ribonuclease|pseudouridine|PseudoU|RNA_pol|RNA polymerase|RNA helicase|DEAD|Helicase_C|CSD|cold[- ]shock|Hfq|Sm|LSM|YTH|DUF55|KOW|NusA|NusB|NusG|Rho|anticodon|tRNA-synt|methyltransferase.*RNA|RNA.*methyltransferase|TruB|TruD|RluA|Trm|Nop|Fibrillarin|Brix|IF|EF|release factor|eRF|RF-1|Cas|CRISPR|RAMP|ProQ|CsrA|YhbY|PIN|Nob1|RNase_PH|Exosome|Csl4|Rrp|L7Ae|Ribosomal_L7Ae|Gar1|Nop10|Pop|Rpp|DUF11|tRNA_bind|RNA_binding|zf-RNPHF|Zn-ribbon.*RNA|YbeY"])
d("helix_turn_helix", ["GO:0003677", "GO:0043565"], [r"helix[- ]turn[- ]helix|HTH|wHTH|winged helix|Sigma70_r4|MarR|TetR|LysR|GntR|AraC|LacI|AsnC|ArsR|IclR|DeoR|MerR|Crp|LuxR|Trans_reg|Rrf2|PadR|HxlR|TrmB|Fe_dep_repress|DtxR|Phage_CI_repr|Cro|XRE|HTH_"])
d("winged_helix", ["GO:0003677"], [r"winged[- ]helix|wHTH|MarR|Trans_reg_C|OmpR|HTH_5|ArsR|PadR|HxlR|DtxR|Fe_dep_repress|MotA_activ|Crp|Rrf2|Penicillinase_R|TrmB|HTH_"])
d("ribbon_helix_helix", ["GO:0003677"], [r"ribbon[- ]helix[- ]helix|RHH|CopG|MetJ|Arc|ParD|RelB|Ribbon"])
d("helix_loop_helix", [], [r"helix[- ]loop[- ]helix|HLH"])
d("transcription_factor", ["GO:0003700", "GO:0140110", "GO:0006355", "GO:0001216", "GO:0001217"], [r"transcription(al)? (factor|regulator|activator|repressor)|regulatory protein|TFIIB|TBP|TATA|Tfb|TFIIS|TFE|sigma|Sigma|HTH|helix[- ]turn[- ]helix|MarR|TetR|LysR|GntR|AraC|LacI|AsnC|ArsR|IclR|DeoR|MerR|Crp|Fnr|LuxR|OmpR|Trans_reg|Rrf2|PadR|HxlR|TrmB|CBFD_NFYB|Histone-like transcription|AbrB|SpoVT|CopG|MetJ|Penicillinase_R|Fur|DtxR|NikR|ModE|LexA|CodY|Sigma54_activat|Response_reg|zf-|WhiB|CarD|DksA|GreA|NusA|NusG|Spt"])
d("activator", ["GO:0045893", "GO:0001216"], [r"activator|Sigma54_activat|enhancer-binding"])
d("repressor", ["GO:0045892", "GO:0001217"], [r"repressor|LexA|TetR|LacI|Penicillinase_R|BlaI|MecI|Fur|DtxR|NikR|ArgR|PurR|TrpR|CodY|Rrf2|AbrB|SpoVT|MarR|Phage_CI_repr|Cro|HipB"])
d("sigma_factor", ["GO:0016987", "GO:0006352"], [r"sigma|Sigma"])
d("anti_sigma", ["GO:0016989"], [r"anti-sigma|anti_sigma|RsbW|SpoIIAA|STAS|RseA|FecR|FlgM"])
d("arac_family", [], [r"AraC|HTH_18|HTH_AraC"])
d("asnc_family", [], [r"AsnC|Lrp"])
d("crp_fnr_family", [], [r"Crp|CRP|Fnr|FNR|cNMP|Cyclic nucleotide|cAMP|HTH_Crp"])
d("gntr_family", [], [r"GntR"])
d("iclr_family", [], [r"IclR"])
d("laci_family", [], [r"LacI|Peripla_BP_3"])
d("luxr_family", [], [r"LuxR|GerE|HTH_LUXR"])
d("lysr_family", [], [r"LysR"])
d("marr_family", [], [r"MarR"])
d("ompr_family", [], [r"OmpR|Trans_reg_C"])
d("tetr_family", [], [r"TetR"])
d("narl_family", [], [r"NarL|LuxR|GerE"])
d("response_regulator", ["GO:0000160", "GO:0000156"], [r"response regulator|Response_reg|receiver|CheY|OmpR|NarL|Trans_reg_C"])
d("sensor", [], [r"sensor|sensing|PAS|GAF|Cache"])
d("transcription", ["GO:0006351", "GO:0003899", "GO:0006355", "GO:0003700"], [r"transcription|RNA[- ]pol|RNA polymerase|sigma|Sigma|NusA|NusG|GreA|Rho|TFIIB|TBP|Tfb|TFIIS|Spt4|Spt5|RNA_pol|DksA|CarD|anti-termination|termination"])
d("transcription_termination", ["GO:0006353", "GO:0031564"], [r"terminat|Rho|NusA|NusB|NusG|anti-termination"])
d("transcription_elongation", ["GO:0006354"], [r"elongation|GreA|NusA|NusG|Spt|TFIIS"])
d("rna_polymerase", ["GO:0003899", "GO:0000428"], [r"RNA[- ]pol|RNA polymerase|RNA_pol|RpoE|Rpb|Rpo"])
d("dna_polymerase", ["GO:0003887", "GO:0034061"], [r"DNA[- ]pol|DNA polymerase|DNA_pol|PolY|IMS|Pol III|Pol_|PolB|PolD|Poly|DNA_primase|clamp|PCNA|Rep_fac|RFC|DnaQ|delta|epsilon|RNA_pol_Rpb"])
d("polymerase", ["GO:0034061", "GO:0003899", "GO:0016779"], [r"polymerase|Pol_|DNA_pol|RNA_pol|PolyA_pol|primase|Poly"])
d("primase", ["GO:0003896"], [r"primase|DnaG|Toprim|PriA|PriL|PriS|Pri[LS]|Prim_Pol"])
d("sliding_clamp", ["GO:0030337", "GO:0009360"], [r"clamp|PCNA|DNA_pol3_beta|beta subunit"])
d("clamp_loader", ["GO:0003689", "GO:0005663"], [r"clamp[- ]loader|clamp loading|Replication factor C|Rep_fac_C|RFC|DNA_pol3_delta|DNA_pol3_gamma|DNA_pol3_tau|HolA|HolB|Rad17|delta"])
d("ligase_dna", ["GO:0003909", "GO:0003910"], [r"DNA[- ]ligase|DNA_ligase|ligase.*DNA"])
d("dna_repair", ["GO:0006281", "GO:0000725", "GO:0006298", "GO:0006284", "GO:0006289", "GO:0006302"], [r"repair|Rad|RecA|RecF|RecN|RecO|RecR|RecJ|RecQ|RecG|RecB|RecC|RecD|RuvA|RuvB|RuvC|UvrA|UvrB|UvrC|UvrD|MutS|MutL|MutY|MutT|MutM|Mfd|UDG|uracil[- ]DNA|glycosylase|Fapy|HhH-GPD|Endonuclease_V|AP_endonuc|endonuclease (III|IV|V)|photolyase|DNA_photolyase|Photolyase|Mre11|Rad50|NurA|HerA|XPG|ERCC|Sbc|DinB|UmuC|IMS|PolY|LexA|SOS|Ku|LigD|DNA_ligase|Nuc-?ase|NERD|DEAD|ResIII|helicase|SMC_N|SNF2|ATP-dependent DNA|Rep_fac-A|DNA_pol_B_exo|Exonuc|exonuclease|RadC|RmuC|Lhr|Hef|YqgF|RadA|DNA_mis_repair|Tex|NUDIX|Nudix|MGMT|DNA_binding_1|Methyltransf_1N|AlkB|2OG-FeII|Mpg|Tag|AlkA|3mg|DNA_glycosylase|DnaQ|RNase_T"])
d("mismatch_repair", ["GO:0006298", "GO:0030983"], [r"mismatch|MutS|MutL|MutH|Vsr|DNA_mis_repair"])
d("base_excision_repair", ["GO:0006284", "GO:0019104"], [r"base[- ]excision|glycosylase|UDG|uracil[- ]DNA|Fapy|HhH-GPD|AP_endonuc|AP endonuclease|apurinic|MutY|MutM|Endonuclease_V|Nth|Nei|Tag|AlkA|3mg|Mpg|DNA_glycosylase|OGG|EndoIII|Exo_endo_phos|endonuclease IV|Xylose isomerase-like TIM barrel"])
d("nucleotide_excision_repair", ["GO:0006289", "GO:0009381"], [r"nucleotide[- ]excision|excinuclease|Uvr[ABC]|UvrD|XPG|XPD|XPB|ERCC|Rad2|Rad3|Mfd|DEAD_2|Helicase_C_2"])
d("recombinational_repair", ["GO:0000725", "GO:0000724", "GO:0006310"], [r"recombinat|RecA|RecF|RecN|RecO|RecR|RecJ|RecQ|RecG|RecB|RuvA|RuvB|RuvC|Rad51|Rad52|Rad50|Mre11|Sbc|RadA|Holliday|NurA|HerA|Lhr|Hef|double[- ]strand break|RmuC|YqgF"])
d("sos_response", ["GO:0009432"], [r"\bSOS|LexA|UmuC|UmuD|DinB|DinD|RecA|SulA|DinF|DinG|DinI|LEXA"])
d("recombinase", ["GO:0000150", "GO:0006310", "GO:0015074"], [r"recombinase|resolvase|Resolvase|integrase|Phage_integrase|Phage_int|Recombinase|RecA|Rad51|invertase|XerC|XerD|Tn3|Serine recombinase|Tyrosine recombinase|site-specific"])
d("resolvase", ["GO:0000150", "GO:0006310"], [r"resolvase|Resolvase|Holliday|RuvC|YqgF|Hjc|invertase"])
d("integrase", ["GO:0015074", "GO:0008907"], [r"integrase|Phage_integrase|Phage_int|rve|DNA integration"])
d("transposase", ["GO:0004803", "GO:0006313", "GO:0032196"], [r"transpos|Transposase|DDE|Tnp|IS[0-9]|Y1_Tnp|OrfB_|OrfB_IS605|OrfB_Zn|HTH_Tnp|rve|MULE|Tn7|Tn3|TnpB|Cas12f1-like"])
d("mobile_element", ["GO:0004803", "GO:0006313", "GO:0032196", "GO:0015074", "GO:0006310", "GO:0000150"], [r"transpos|Transposase|DDE|Tnp|\bIS[0-9]|Y1_Tnp|OrfB|rve|MULE|Tn[0-9]|integrase|Phage_int|resolvase|recombinase|invertase|excisionase|relaxase|Mob|mobiliz|conjugat|plasmid|phage|Phage|prophage|intron|reverse transcriptase|RVT|retron|Retron|homing|HNH|GIY-YIG|LAGLIDADG|intein|Intein|Hint|insertion|Group II|maturase|IS66|IS200|IS21|IS3|IS4|IS5|ISL3|IS30|IS110|IS256|IS481|IS630|IS982|IS1380|ISC|TnpB|Cas12f1-like|TniQ|TnsA|TnsB|TnsC"])
d("insertion_sequence", [], [r"\bIS[0-9]|insertion sequence|Tnp|Transposase"])
d("homing_endonuclease", [], [r"homing|LAGLIDADG|GIY-YIG|HNH|intein"])
d("excisionase", [], [r"excisionase|Xis|Excisionase"])
d("phage_integrase", ["GO:0015074"], [r"Phage_int|phage integrase|integrase"])
d("ssb_protein", ["GO:0003697"], [r"single[- ]strand(ed)?[- ]DNA[- ]binding|SSB|RPA|Rep_fac-A"])
d("chromatin", ["GO:0000785", "GO:0030527", "GO:0006325"], [r"histone|Histone|chromatin|nucleosome|HMG|Alba|Sul7|Cren7|MC1|CC1|CBFD_NFYB|HU|IHF|H-NS|Lsr2|Nucleoid|NAP|Bac_DNA_binding"])
d("histone", ["GO:0000786", "GO:0030527"], [r"histone|Histone|CBFD_NFYB|HMF|HMf"])
d("divisome", ["GO:1990586", "GO:0032153"], [r"divisome|Fts[A-Z]|FtsZ|ZipA|ZapA|SepF|Septum|septal"])
d("ftsz", ["GO:0043093"], [r"FtsZ|Tubulin/FtsZ|Tubulin"])
d("ribosomal_protein", ["GO:0003735", "GO:0005840"], [r"ribosomal protein|Ribosomal|ribosom.*subunit"])
d("trna_synthetase", ["GO:0004812", "GO:0043039", "GO:0006418"], [r"tRNA[- ]synthetase|tRNA-synt|aminoacyl[- ]tRNA|tRNA ligase|anticodon[- ]binding|Anticodon|tRNA_anti|tRNA-synt_|DALR|TGS|SelR|tRNA_edit|B3_4|B5|Val_tRNA|FDX-ACB|GAD|tRNA_bind"])
d("trna_modification", ["GO:0006400", "GO:0008175", "GO:0030488", "GO:0002098", "GO:0006400"], [r"tRNA", r"queuosine|wybutosine|archaeosine|thiouridine|pseudouridine|TruA|TruB|TruD|RluA|Trm|Tgt|TYW|Tad|TilS|MnmA|MnmE|MnmG|GidA|ThiI|Thg1|Kae1|YgjD|Sua5|YrdC|TsaB|TsaC|TsaD|TsaE|Mia|MiaA|MiaB|QueA|QueC|QueD|QueE|QueF|DUS|Dus|NSUN|PUA|THUMP|TRAM|Elp3|CTU|Ncs|DUF55|Pcc1|Cgi121|KEOPS|Bud32|RtcB|Rtc|CCA|tRNA_NucTran"])
d("rrna_modification", ["GO:0000154", "GO:0008649", "GO:0031167"], [r"rRNA|RsmA|RsmB|RsmC|RsmD|RsmE|RsmF|RsmG|RsmH|RsmI|RsmJ|RlmA|RlmB|RlmC|RlmD|RlmE|RlmF|RlmG|RlmH|RlmI|RlmJ|RlmK|RlmL|RlmM|RlmN|Fibrillarin|Nop|Nep1|KsgA|Dim1|ErmC|Erm|RrmJ|FtsJ|SpoU|Methyltr_RsmB|RluA|RluB|RluC|RluD|RsuA|TlyA|DUF|Gar1|Nop10|Cbf5|Pseudouridine"])
d("pseudouridine_synthase", ["GO:0009982", "GO:0001522"], [r"pseudouridine|PseudoU|TruA|TruB|TruD|RluA|RsuA|Pus|Cbf5"])
d("ribosome_biogenesis", ["GO:0042254"], [r"ribosome biogenesis|ribosome assembly|biogenesis|RIO|Rio|Nob1|Brix|Nop|Fibrillarin|KsgA|RimM|RbfA|Era|Obg|EngA|Der|YihA|YchF|RsgA|RimP|Rim|RlmN|MMR_HSR1|GTPase"])
d("cold_shock", ["GO:0009409"], [r"cold[- ]shock|CSD|Csp"])
d("heat_shock", ["GO:0009408"], [r"heat[- ]shock|HSP|Hsp|DnaK|DnaJ|GroEL|GroES|Cpn|HtpG|ClpB|IbpA|small heat|alpha crystallin|HSP20|HSP33|HSP90|HSP70"])
d("hsp60", ["GO:0006457", "GO:0051082"], [r"Cpn60|GroEL|chaperonin|TCP-1|HSP60|Hsp60"])
d("hsp70", ["GO:0006457"], [r"HSP70|Hsp70|DnaK"])
d("hsp90", ["GO:0006457"], [r"HSP90|Hsp90|HtpG"])
d("small_hsp", [], [r"HSP20|alpha crystallin|small heat|IbpA|sHSP"])
d("oxidative_stress", ["GO:0006979", "GO:0016209", "GO:0004601", "GO:0004784", "GO:0004096", "GO:0051920"], [r"oxidative|peroxid|catalase|superoxide|Sod_|Dps|DPS|ferritin|AhpC|AhpD|Redoxin|peroxiredoxin|Prx|thioredoxin reductase|glutathione peroxidase|GSHPx|MsrA|MsrB|SelR|methionine sulfoxide|OsmC|Ohr|OhrR|Rubrerythrin|rubrerythrin|Desulfoferrodoxin|Rubredoxin|superoxide reductase|SOR|NADH_oxidase|OxyR|SoxR|PerR|Fur|Bacterioferritin|YjgF|RidA|Glyoxalase|DUF.*oxid|Chloroperoxidase|Haloperoxidase|Catalase|Peroxidase|peroxide"])
d("osmotic_stress", ["GO:0006970"], [r"osmo|OsmC|OpuA|OpuC|BetT|ProP|ProU|Osmotic"])
d("thiosulfate", [], [r"thiosulfate|thiosulphate|rhodanese|Sox[ABCDXYZ]|TST"])

# ---- Metabolism -----------------------------------------------------------
d("amino_acid_degradation", ["GO:0009063", "GO:0006552", "GO:0006574", "GO:0006559", "GO:0006572", "GO:0019464", "GO:0006548", "GO:0006527", "GO:0006554", "GO:0006550", "GO:0006567", "GO:0006561"], [r"degradation|catabol|deiminase|ADI|arginase|Arginase|GCV|glycine cleavage|Gcv|ELFV_dehydrog|glutamate dehydrogenase|leucine dehydrogenase|phenylalanine dehydrogenase|valine dehydrogenase|AlaDh|alanine dehydrogenase|Proline dehydrogenase|Pro_dh|ProDH|Urocanase|HutU|HutH|HutI|HutG|Tryptophanase|Beta_elim_lyase|Asparaginase|Glutaminase|Ald_Xan|Thr_dehydrat|Ser_dehydrat|SDH|Aminotran|TPP_enzyme|2-oxoacid|OKR|Ornithine_cyclodeaminase|OCD_Mu_crystall"])
d("branched_chain_aa", ["GO:0009082", "GO:0009081", "GO:0009097", "GO:0009098", "GO:0009099", "GO:0009083"], [r"branched[- ]chain|isoleucine|leucine|valine|ILVD|KARI|IlvN|IlvB|IlvC|IlvD|IlvE|LeuA|LeuB|LeuC|LeuD|IPMS|isopropylmalate|acetohydroxy|acetolactate|ketol-acid|dihydroxy-acid|BCAT|LIV|BCCT"])
d("histidine_biosynthesis", ["GO:0000105"], [r"histidine|His[A-IZ]|HisG|IGPD|PRA-|PRAI|His_biosynth|HisIE|HisF|HisH|HisA|HisB|HisC|HisD|ATP-PRT|imidazoleglycerol"])
d("asparagine_synthesis", ["GO:0006529", "GO:0004066"], [r"asparagine synth|Asn_synthase|AsnA|AsnB"])
d("glutamine_synthetase", ["GO:0004356"], [r"glutamine synthetase|Gln-synt|GlnA"])
d("ammonia_assimilation", ["GO:0004356", "GO:0006542", "GO:0019676", "GO:0015930", "GO:0016639"], [r"glutamine synthetase|Gln-synt|GlnA|glutamate synthase|GltB|GltD|GOGAT|Glu_synthase|glutamate dehydrogenase|ELFV_dehydrog|ammonia|ammonium|AmtB|Ammonium_transp|GlnB|P-II|PII"])
d("nitrogen_metabolism", ["GO:0006807", "GO:0071941", "GO:0042128", "GO:0009399", "GO:0019333", "GO:0006542", "GO:0019740"], [r"nitrogen|nitrate|nitrite|nitric|nitrous|ammonia|ammonium|urea|urease|Urease|nitrogenase|Nitrogenase|Nif|nif|glutamine synthetase|Gln-synt|GlnB|P-II|Amt|Ammonium|Nar[GHIJ]|Nap[AB]|Nir[BDKS]|Nor[BC]|Nos[ZDLF]|Nrf|Nitr_red|Molybdopterin.*nitrate|cyanase|Cyanate|nitrilase|Nitrilase|CN_hydrolase|hydroxylamine|HAO|AmoA|AmoB|AmoC|Amo|PmoA|pmoA"])
d("nitrate_reduction", ["GO:0008940", "GO:0042128", "GO:0042126", "GO:0009325", "GO:0050464"], [r"nitrate reductase|nitrate reduction|Nitrate_red|Nar[GHIJ]|Nap[AB]|NasA|NasC|NR_"])
d("nitrite_reductase", ["GO:0098809", "GO:0050421", "GO:0008942", "GO:0042279"], [r"nitrite reductase|Nir[BDKS]|NirK|NirS|NrfA|Nitr_red|NIR_SIR|Cytochrom_C552"])
d("denitrification", ["GO:0019333", "GO:0050421", "GO:0016966", "GO:0050304"], [r"denitrif|nitric[- ]oxide reductase|nitrous[- ]oxide reductase|NorB|NorC|NosZ|NosD|NosL|NosF|NosY|NirK|NirS|Nir[SK]|Nor[BC]|Nos[ZDLFY]|NapA|NarG"])
d("nitrogen_fixation", ["GO:0009399", "GO:0016163"], [r"nitrogen fixation|nitrogenase|Nif|nif[A-Z]|Fer4_NifH|Oxidored_nitro|NifU|NifB|NifN|NifE|NifT|NifW|NifZ|FixA|FixB|Rnf"])
d("nitrogenase", ["GO:0016163"], [r"nitrogenase|Fer4_NifH|Oxidored_nitro|NifD|NifK|NifH|VnfD|AnfD"])
d("nitrogenase_maturation", [], [r"Nif[BENTUVWXYZ]|nitrogenase.*(cofactor|maturation|assembly|synthesis)|FeMo-?co|NifB|NifE|NifN|NifU|NifV"])
d("nitrification", ["GO:0019329"], [r"nitrif|ammonia monooxygenase|AmoA|AmoB|AmoC|hydroxylamine|HAO|nitrite oxidoreductase|NxrA|NxrB"])
d("ammonia_oxidation", ["GO:0019329"], [r"ammonia monooxygenase|AmoA|AmoB|AmoC|Amo|hydroxylamine|HAO|ammonia oxidation"])
d("urease", ["GO:0009039", "GO:0043419"], [r"urease|Ure[A-GJ]"])
d("urea_metabolism", ["GO:0019627", "GO:0043419"], [r"urea|urease|Ure[A-GJ]|arginase|Urea"])
d("sulfur_oxidation", ["GO:0019417", "GO:0070221", "GO:0019418", "GO:0070225"], [r"sulfur oxid|sulfide oxid|sulfite oxid|thiosulfate oxid|Sox[A-Z]|Sqr|SQR|sulfide:quinone|Fcc|FCSD|Sor|sulfur oxygenase|DsrE|TusA|Hdr-like"])
d("sulfur_assimilation", ["GO:0000103", "GO:0070814", "GO:0019344", "GO:0006790"], [r"sulfate assimilation|sulfur assimilation|Cys[ACDHIJNKE]|PAPS|APS kinase|sulfate adenylyltransferase|ATP-sulfurylase|Sat\b|CysK|cysteine synth"])
d("selenium_metabolism", ["GO:0001887", "GO:0004756", "GO:0016260", "GO:0001514"], [r"seleno|selenium|Sel[ABDX]|SelD|SelA|SelB|SelR|selenocysteine|selenophosphate|NiFeSe|YbbB|SelU"])
d("carbohydrate_active", ["GO:0005975", "GO:0016798", "GO:0016757", "GO:0030246", "GO:0016829", "GO:0000272"], [r"glyco|Glyco|glycosyl|carbohydrate|sugar|saccharide|glucan|chitin|cellul|xylan|pectin|pectate|amylase|amylo|starch|glycogen|mannan|galactan|fructan|levan|dextran|pullulan|CBM|cohesin|dockerin|Dockerin|Cohesin|polysaccharide|Polysacc|lyase.*(pectate|alginate|chondroitin|hyaluronate|heparin)|lysozyme|chitinase|cellulase|xylanase|mannanase|galactosidase|glucosidase|mannosidase|xylosidase|fucosidase|sialidase|glucuronidase|hexosaminidase|arabinofuranosidase|trehalase|invertase|sucrase|kinase.*(sugar|carbohydrate)|carbohydrate kinase|PfkB|FGGY|ROK|Aldolase|aldolase|epimerase|UDP|dTDP|GDP|CDP|NDP|nucleotide-sugar|nucleotide sugar|rhamnose|RmlA|RmlB|RmlC|RmlD|WecB|Wec|Wbp|Wza|Wzx|Wzy|Wzz|Epimerase|NAD_binding_4|GDP_Man|Polysacc_deac|NodB|deacetylase|carbohydrate esterase|CE_|Esterase_PHB|PGM_PMM|phosphoglucomutase|phosphomannomutase|Transglycosylase|transglycosylase|SLT|Transgly|Hexapep|Glyco_tranf|Glyco_trans|Glycos_trans|Glyco_hydro|Alpha-amylase|Cellulase|Chitin|CBM_|SBBP|F5_F8_type_C|RicinB|Laminin_G|PA14|Big|NPCBM|Gal_Lectin|Lectin|lectin|beta-propeller.*(sugar|carbohydrate)"])
d("carbohydrate_binding", ["GO:0030246", "GO:0030247", "GO:0008061", "GO:2001070", "GO:0030248"], [r"carbohydrate[- ]binding|sugar[- ]binding|CBM|chitin[- ]binding|Chitin_bind|cellulose[- ]binding|starch[- ]binding|lectin|Lectin|RicinB|F5_F8_type_C|Laminin_G|PA14|NPCBM|Gal_Lectin|glycan[- ]binding|polysaccharide[- ]binding|SBBP|Beta-propeller|beta-propeller|Dockerin|Cohesin|X2|Big_"])
d("carbohydrate_esterase", ["GO:0016788"], [r"carbohydrate esterase|CE_|acetyl xylan|pectin ?esterase|Polysacc_deac|NodB|deacetylase|Esterase"])
d("glycoside_hydrolase", ["GO:0004553", "GO:0016798"], [r"glycosi?d(ase|e hydrolase)|glycosyl hydrolase|Glyco_hydro|amylase|cellulase|chitinase|lysozyme|glucosidase|galactosidase|mannosidase|xylosidase|fucosidase|glucanase"])
d("polysaccharide_lyase", ["GO:0016837"], [r"polysaccharide lyase|pectate lyase|Pectate_lyase|alginate lyase|Alginate_lyase|chondroitinase|heparinase|hyaluronate lyase|PL_|Lyase_|lyase.*polysaccharide"])
d("xylanase", ["GO:0031176"], [r"xylan|xylanase"])
d("mannanase", ["GO:0016985"], [r"mannan|mannanase"])
d("pectinase", [], [r"pectin|pectate|Pectate"])
d("sugar_metabolism", ["GO:0005975", "GO:0044262", "GO:0019318"], [r"sugar|carbohydrate|glucose|fructose|galactose|mannose|xylose|arabinose|ribose|rhamnose|fucose|sucrose|maltose|lactose|trehalose|glycerol|hexose|pentose|PfkB|FGGY|ROK|Aldolase|aldolase|isomerase.*(sugar|ribose|xylose|glucose|mannose|galactose|arabinose|triose)|Epimerase|epimerase|GFO_IDH|DeoR|PTS|PGM_PMM|TIM|Tim|Glyco|UDP|dTDP|GDP|NDP|Polysacc"])
d("glycolysis", ["GO:0006096", "GO:0006094", "GO:0061621"], [r"glycoly|glyceraldehyde|Gp_dh|GAPDH|phosphoglycerate|PGK|PGAM|triose[- ]?phosphate|TIM|Enolase|enolase|pyruvate kinase|PK|hexokinase|glucokinase|phosphofructokinase|PFK|PfkB|fructose[- ]bisphosphate|F_bP|FBP|Aldolase|aldolase|glucose-6-phosphate isomerase|PGI|PEP|ADP-specific|AFOR|GAPOR|iPGM|apgM|Metalloenzyme|2,3-bisphosphoglycerate"])
d("gluconeogenesis", ["GO:0006094"], [r"gluconeogen|FBPase|fructose-1,6-bisphosphatase|F_bP|PEPCK|PEP carboxykinase|PEP_carboxykinase|PPDK|PEP-utilizers|pyruvate phosphate dikinase|phosphoenolpyruvate synthase|PPS|Pyruvate carboxylase|PC|Aldolase|aldolase|Gp_dh|glyceraldehyde|phosphoglycerate"])
d("pentose_phosphate", ["GO:0006098", "GO:0009052", "GO:0009051"], [r"pentose|G6PD|6PGD|6PGL|6-phosphogluconate|transketolase|Transketolase|transaldolase|Transaldolase|ribulose|Ribul|ribose-5-phosphate|Rib_5-P|RPE|RpiA|RpiB|LacAB_rpiB|glucose-6-phosphate dehydrogenase|NAD_binding_2|Sacchrp_dh_NADP|oxidative pentose"])
d("glyoxylate_cycle", ["GO:0006097", "GO:0004451", "GO:0004474"], [r"glyoxylate|isocitrate lyase|malate synthase|Isocitrate_lyase|ICL|Malate_synthase|MS_"])
d("fermentation", ["GO:0006113", "GO:0019660", "GO:0019664"], [r"ferment|lactate dehydrogenase|alcohol dehydrogenase|pyruvate formate[- ]lyase|PFL|butyrate|butanol|acetate kinase|AckA|phosphotransacetylase|Pta|ethanol"])
d("energy_metabolism", ["GO:0006091", "GO:0015980", "GO:0022900", "GO:0006119", "GO:0015986", "GO:0009060", "GO:0009061"], [r"energy|ATP synth|ATP-synt|respirat|electron (transfer|transport)|oxidative phosphorylation|NADH|quinone|cytochrom|ferredoxin|hydrogenase|formate dehydrogenase|proton|Complex1|COX|Rieske|H\+|Na\+"])
d("respiration", ["GO:0009060", "GO:0009061", "GO:0022900", "GO:0045333", "GO:0004129", "GO:0015990", "GO:0019646"], [r"respirat|cytochrome (c|b|bd|d|o|aa3|cbb3)|oxidase|COX|NADH.*(dehydrogenase|oxidoreductase)|NADH_dh|Complex1|quinone|quinol|ubiquinone|menaquinone|Rieske|succinate dehydrogenase|Sdh|fumarate reductase|Frd|nitrate reductase|Nar|Nap|DMSO|TMAO|Molybdopterin|Molydop|electron transfer|NDH|NQR|Rnf|Oxidored_q|Proton_antipo|terminal oxidase|Cyt_bd|CydA|CydB|Cytochrom_C|Cytochrome_C|COX[0-9]|Cytochrom_B|Cytochrome_B|Qor|DsrMK|HdrD|HdrE|Hdr"])
d("aerobic_respiration", ["GO:0009060", "GO:0004129", "GO:0019646"], [r"aerobic respiration|cytochrome c oxidase|COX|Cyt_bd|CydA|CydB|cytochrome bd|cytochrome o|cytochrome aa3|cbb3|CcoN|CcoO|CcoP|CcoQ|Heme-copper|HCO|quinol oxidase|terminal oxidase|alternative oxidase|AOX"])
d("anaerobic_respiration", ["GO:0009061"], [r"anaerobic respiration|nitrate reductase|fumarate reductase|DMSO reductase|TMAO reductase|tetrathionate reductase|Nar[GHIJ]|Nap[AB]|Frd|DmsA|DmsB|DmsC|TorA|TtrA|Psr|Phs|Sre|Srr|arsenate reductase|ArrA|selenate|chlorate|perchlorate|Molybdopterin"])
d("terminal_oxidase", ["GO:0004129", "GO:0016682", "GO:0009486"], [r"terminal oxidase|cytochrome (c|bd|o|aa3|cbb3|d) oxidase|COX|Cyt_bd|CydA|CydB|CcoN|CcoO|CcoP|Heme-copper|alternative oxidase|AOX|quinol oxidase|Cytochrom_C_oxidase|COX[0-9]"])
d("cytochrome_c_oxidase", ["GO:0004129", "GO:0045277"], [r"cytochrome c oxidase|COX|CcoN|CcoO|CcoP|Cytochrom_C_oxidase|Heme-copper"])
d("cytochrome_bd_oxidase", ["GO:0070069", "GO:0016682"], [r"cytochrome bd|Cyt_bd|CydA|CydB|CydX|bd-type"])
d("atp_synthesis", ["GO:0015986", "GO:0045259", "GO:0046933", "GO:0046961", "GO:0033178"], [r"ATP synth|ATP-synt|ATP synthase|OSCP|V-ATPase|V_ATPase|A-type|A-ATPase|vATP|ATP-synt_[A-Z]|H\+-ATPase|F0|F1|FOF1|Na\+-ATPase|ntp"])
d("photosynthesis", ["GO:0015979", "GO:0009521", "GO:0009522", "GO:0009523", "GO:0030076", "GO:0009765", "GO:0019684", "GO:0009767", "GO:0009772"], [r"photosynth|photosystem|Photosystem|PSI|PSII|Psa[A-Z]|Psb[A-Z]|PsbP|reaction cent|Photo_RC|light[- ]harvest|LHC|antenna|chlorophyll|bacteriochlorophyll|phycobili|phycocyanin|Phycobilisome|carotenoid|YCF|Ycf|thylakoid|plastocyanin|Cyt_b6|cytochrome b6|Cytochrom_B6|PsbO|PsbU|PsbV|NdhF|Chl|Bch|Pcb|IsiA|OCP|PSII_BNR"])
d("photosystem_i", ["GO:0009522", "GO:0015979"], [r"photosystem I\b|photosystem I |Photosystem I\b|PsaA|PsaB|PsaC|PsaD|PsaE|PsaF|PsaL|PSI\b|Psa[A-Z]"])
d("photosystem_ii", ["GO:0009523", "GO:0009654", "GO:0015979"], [r"photosystem II|Photosystem II|PSII|Psb[A-Z]|PsbP|PsbO|PsbU|PsbV|YCF48|D1|D2|CP43|CP47|oxygen-evolving"])
d("bacteriochlorophyll", ["GO:0030494"], [r"bacteriochlorophyll|Bch"])
d("carbon_fixation", ["GO:0015977", "GO:0016984", "GO:0019253", "GO:0043427"], [r"carbon fixation|CO2 fixation|RuBisCO|rubisco|ribulose[- ]bisphosphate carboxylase|phosphoribulokinase|PRK|Calvin|carbonic anhydrase|Pro_CA|Carb_anhydrase|PEPcase|PEP carboxylase|reductive|acetyl-CoA carboxylase|Wood-Ljungdahl|CODH|ACS|4-hydroxybutyryl|3-hydroxypropionate|citrate lyase|ATP-citrate|2-oxoglutarate:ferredoxin|pyruvate:ferredoxin|formate dehydrogenase"])
d("calvin_cycle", ["GO:0019253", "GO:0016984"], [r"Calvin|RuBisCO|rubisco|ribulose[- ]bisphosphate|phosphoribulokinase|PRK|sedoheptulose|SBPase|CP12"])
d("rubisco", ["GO:0016984"], [r"RuBisCO|rubisco|ribulose[- ]bisphosphate carboxylase"])
d("prk", ["GO:0008974"], [r"phosphoribulokinase|\bPRK\b"])
d("reverse_tca", [], [r"ATP-citrate lyase|citrate lyase|2-oxoglutarate:ferredoxin|OGFOR|KOR|fumarate reductase|reductive TCA|rTCA"])
d("3hp_bicycle", [], [r"3-hydroxypropionate|malonyl-CoA reductase|propionyl-CoA synthase"])
d("dicarboxylate_4hb", [], [r"4-hydroxybutyryl|dicarboxylate"])
d("wood_ljungdahl", ["GO:0043885", "GO:0019385", "GO:0006084", "GO:0043884", "GO:0018492", "GO:0006730"], [r"Wood[- ]Ljungdahl|carbon[- ]monoxide dehydrogenase|CO dehydrogenase|CODH|acetyl-CoA synthase|ACS|Cdh[A-E]|CdhC|CdhD|CO_dh|formyltetrahydrofolate|FTHFS|Methylenetetrahydrofolate|MTHFR|MetF|FolD|methenyl|formate[- ]tetrahydrofolate|THF|tetrahydrofolate|corrinoid|AcsE|AcsC|AcsD|CooS|CooC|CooT|CODH_A|Formylmethanofuran|Fmd|Fwd|FTR|Mtd|Mch|Mer|Mtr|methanopterin|H4MPT|CO_dh"])
d("co_oxidation", ["GO:0043885", "GO:0018492", "GO:0008805"], [r"carbon[- ]monoxide|CO dehydrogenase|CODH|CooS|CoxL|CoxM|CoxS|CO_dh|Cdh|CooT|CooC|Ald_Xan"])
d("archaeal_one_carbon", ["GO:0006730", "GO:0015948", "GO:0019386"], [r"methanopterin|H4MPT|methanofuran|Formylmethanofuran|Fmd|Fwd|FTR|Mtd|Mch|Mer|Mtr|MtrH|F420|coenzyme M|CoM|heterodisulfide|Hdr|MCR|methyl-coenzyme M|methanogen|FmdE|FTR_C|CdhD|CdhC|CO_dh|THMPT|Mtx|MtaA|MtbA|MttB|MtaB|MtmB|methyltransferase.*methano"])
d("methanotrophy", ["GO:0015947", "GO:0015049", "GO:0018662"], [r"methane monooxygenase|MMO|pmo|Pmo|mmo|methanotroph|methanol dehydrogenase|Mxa|Xox|AmoA|PmoA|particulate methane|soluble methane"])
d("methylotrophy", ["GO:0015945", "GO:0046170"], [r"methanol|methylamine|methylotroph|formaldehyde|Mxa|Xox|Mau|trimethylamine|DMSP|methylated"])
d("f420_reducing", [], [r"F420[- ]reducing|coenzyme F420 hydrogenase|FrhA|FrhB|FrhG|F420H2|Frh"])
d("hydrogen_metabolism", ["GO:0006007", "GO:0008901", "GO:0033748", "GO:0047806", "GO:0050454", "GO:0051912"], [r"hydrogenase|\bH2\b|hydrogen|NiFe|FeFe|Fe_hyd|NiFeSe|Hyd[A-Z]|Hup[A-Z]|Hox|Hyp[A-F]|HycI|HyaE|HupF|HypC|HypF|HypD|HypA|HypB|HypE|Ni_insertion|Frh|Mvh|Ech|Mbh|Eha|Ehb|Hyf|Hyc"])
d("hydrogenase", ["GO:0008901", "GO:0033748", "GO:0047806", "GO:0050454", "GO:0051912", "GO:0047985", "GO:0009375"], [r"hydrogenase"])
d("nife_hydrogenase", ["GO:0008901", "GO:0033748", "GO:0047806", "GO:0016151"], [r"NiFe|Ni-Fe|nickel[- ](dependent|iron) hydrogenase|NiFeSe|Nickel-dependent hydrogenase|NiFe_hyd|Ni_hydr|Complex1_49kDa.*hydrogenase|MbhD|MbhE|Mbh|Ech|Hyc|Hyf"])
d("nifese_hydrogenase", [], [r"NiFeSe"])
d("fefe_hydrogenase", ["GO:0008901"], [r"Fe[- ]only hydrogenase|FeFe|\[FeFe\]|Iron hydrogenase|iron[- ]only|Fe_hyd|iron hydrogenase|Iron only hydrogenase"])
d("fe_only_hydrogenase", ["GO:0047068"], [r"\bHmd|methylenetetrahydromethanopterin dehydrogenase|HMD"])
d("hydrogenase_maturation", ["GO:0051604", "GO:0016151", "GO:0065003"], [r"hydrogenase.*(maturation|expression|formation|assembly|protease)|Hyp[A-F]|HypA|HypB|HypC|HypD|HypE|HypF|HycI|HupF|HupG|HyaE|HyaF|HybE|HybG|HydE|HydF|HydG|zf-HYPF|Ni_insertion|nickel incorporation|NiFe_hyd_mat|HupF_HypC|HycI|Hydrogenase/urease nickel incorporation|HypF_C|Kae1-like|Carbamoyltransferase"])
d("mbh_hydrogenase", [], [r"\bMBH\b|Mbh|mbh[A-N]"])
# Ech is HydDB Group 4e; "energy-converting hydrogenase" also names Eha/Ehb (Groups 4h/4i).
d("ech_hydrogenase", [], [r"\bEch\b|ech hydrogenase"])
d("tetrapyrrole", [], [])
d("fad_biosynthesis", ["GO:0006747", "GO:0009231", "GO:0003919"], [r"FAD synth|FMN adenylyltransferase|riboflavin|Rib[A-H]|FAD_syn|Flavokinase|FAD synthetase"])
d("isoprenoid_biosynthesis", ["GO:0008299", "GO:0019288", "GO:0019287", "GO:0016114", "GO:0006720"], [r"isoprenoid|terpen|terpene|prenyl|polyprenyl|farnesyl|geranyl|squalene|hopene|carotenoid|Lycopene|MEP|DXP|DXS|DXR|IspC|IspD|IspE|IspF|IspG|IspH|MVA|mevalonate|HMG-CoA|IPP|IDI|isopentenyl|UbiA|MenA|polyprenyltransferase|trans-isoprenyl|Polyprenyl_synt|SQS_PSY|CrtB|CrtI|Amino_oxidase|Lycopene_cycl|ERG|Lanosterol|IspA|GGPS|HpnC|HpnD|YgbB|MECDP|LytB|GcpE|NHL|DXP_synthase|DXP_reductoisom"])
d("terpene_synthesis", ["GO:0016114", "GO:0010333"], [r"terpen|terpene|squalene|hopene|Terpene_synth|SQS_PSY|carotenoid"])
d("fatty_acid_synthesis", ["GO:0006633", "GO:0004315", "GO:0004316", "GO:0004318", "GO:0003989"], [r"fatty[- ]acid (biosynth|synth)|beta[- ]ketoacyl|ketoacyl|KAS|Fab[A-Z]|FabA|FabB|FabD|FabF|FabG|FabH|FabI|FabK|FabZ|enoyl[- ]ACP|Enoyl|ACP|acyl carrier|malonyl|Malonyl|Acyl_transf_1|ketoacyl-synt|Ketoacyl-synt|KR|adh_short|PP-binding|biotin carboxyl|carboxyltransferase|ACC|AccA|AccB|AccC|AccD|Biotin_carb|Biotin_lipoyl|CT|desaturase|FA_desaturase"])
d("fatty_acid_degradation", ["GO:0009062", "GO:0006635", "GO:0003995", "GO:0004300", "GO:0003857", "GO:0003988"], [r"beta[- ]oxidation|fatty[- ]acid (degradation|catabol|oxidation)|Acyl-CoA_dh|acyl[- ]CoA dehydrogenase|ECH|enoyl-CoA hydratase|3HCDH|hydroxyacyl-CoA dehydrogenase|thiolase|Thiolase|Fad[A-M]|FadA|FadB|FadD|FadE|FadH|FadJ|AMP-binding|acyl-CoA synthetase|ETF|Acyl-CoA oxidase"])
d("beta_oxidation", ["GO:0006635", "GO:0003995", "GO:0004300", "GO:0003857"], [r"beta[- ]oxidation|Acyl-CoA_dh|acyl[- ]CoA dehydrogenase|ECH|enoyl-CoA hydratase|3HCDH|hydroxyacyl-CoA|thiolase|Thiolase|Fad[ABEJ]"])
d("steroid_metabolism", [], [r"steroid|sterol|hopanoid|3Beta_HSD"])
d("polyketide_synthesis", ["GO:0030639", "GO:0034081"], [r"polyketide|PKS|ketoacyl-synt|Ketoacyl-synt|KR|Acyl_transf_1|PS-DH|ACP|PP-binding|SnoaL|Chalcone|CHS|type III PKS|ketosynthase|aromatase|cyclase"])
d("nrps", ["GO:0017000", "GO:0019184"], [r"nonribosomal|non-ribosomal|NRPS|Condensation|AMP-binding|PP-binding|Thioesterase|peptide synthetase"])
d("secondary_metabolism", ["GO:0044550", "GO:0019748", "GO:0017000"], [r"secondary metabol|antibiotic|polyketide|PKS|nonribosomal|NRPS|SnoaL|Condensation|AMP-binding|PP-binding|Thioesterase|siderophore|lantibiotic|Lant|LanB|LanC|LanM|thiopeptide|bacteriocin|microcin|YcaO|TfuA|halogenase|Trp_halogenase|cytochrome P450|p450|O-methyltransferase|glycosyltransferase|terpen|terpene|DegT|MMPL|Cyclase|cyclase|aromatase|Methyltransf_2|Methyltransf_3"])
d("antibiotic_biosynthesis", ["GO:0017000"], [r"antibiotic (biosynth|synth)|polyketide|NRPS|lantibiotic|bacteriocin|microcin"])
d("siderophore_biosynthesis", ["GO:0019290"], [r"siderophore|IucA|IucC|FhuF|enterobactin|Ent[A-F]|AMP-binding.*siderophore"])
d("nucleotide_metabolism", ["GO:0009117", "GO:0006163", "GO:0006220", "GO:0009165", "GO:0009116", "GO:0009262", "GO:0055086", "GO:0004748"], [r"nucleotide|nucleoside|purine|pyrimidine|ribonucleotide|deoxyribonucleotide|dNTP|NTP|nucleobase|adenine|guanine|cytosine|uracil|thymine|thymidylate|Thymidylat|ThyA|ThyX|dUTP|dUTPase|dCTP|dCMP|dUMP|uridine|cytidine|adenosine|guanosine|inosine|xanthine|hypoxanthine|IMP|GMP|AMP|UMP|CMP|TMP|Pur[A-UKT]|PurA|PurB|PurC|PurD|PurE|PurF|PurH|PurK|PurL|PurM|PurN|PurQ|PurS|PurT|Pyr[A-IR]|PyrB|PyrC|PyrD|PyrE|PyrF|PyrG|PyrH|PyrI|CTP|CPSase|Carbamoyl|OTCace|aspartate carbamoyltransferase|DHO|dihydroorotate|DHOD|orotate|OMP|Orotidine|UPRT|PRTase|Pribosyltran|phosphoribosyl|PRPP|Ribonuc_red|ribonucleotide reductase|RNR|NrdA|NrdB|NrdD|NrdJ|Nrd|ATP-cone|dNK|Thymidylate_kin|Thymidyl|adenylate kinase|ADK|Adenylat|Guanylate_kin|guanylate kinase|NDK|nucleoside diphosphate kinase|Cytidylate_kin|UMP_kinase|APRT|HGPRT|XPRT|NUDIX|MazG|HAM1|ITPase|dUTPase|DUTPase|dCMP_cyt_deam|Deoxycytidine|cytidine deaminase|CDA|adenosine deaminase|AMP deaminase|A_deaminase|Amidohydro|AIRC|SAICAR|AICAR|FGAM|GAR|IMPDH|GMP_synt|GMPS|AdSS|ADSL|NT5|5'-nucleotidase|NT5C|HD|YfbR|DGTP|dGTPase|SAMHD1|RNA_ligase|tRNA_NucTran"])
d("purine_metabolism", ["GO:0006163", "GO:0006164", "GO:0006144", "GO:0009113", "GO:0006188", "GO:0006177"], [r"purine|adenine|guanine|xanthine|hypoxanthine|inosine|adenosine|guanosine|IMP|GMP|AMP|Pur[A-UKT]|AIRC|SAICAR|AICAR|FGAM|GAR|IMPDH|GMP_synt|AdSS|ADSL|APRT|HGPRT|XPRT|adenylosuccinate|phosphoribosylformylglycinamidine|phosphoribosylaminoimidazole|amidophosphoribosyltransferase|A_deaminase|Adenosine deaminase|Amidohydro_1|Urate|urate|allantoin|Allantoinase|MoCF|Xanthine|Ald_Xan|MazG|NUDIX|Nudix|ITPase|HAM1|adenylate kinase|ADK|Adenylat|Guanylate_kin"])
d("pyrimidine_metabolism", ["GO:0006220", "GO:0006221", "GO:0006212", "GO:0009220"], [r"pyrimidine|cytosine|uracil|thymine|thymidine|uridine|cytidine|UMP|CMP|TMP|dUMP|dTMP|dUTP|CTP|Pyr[A-IR]|PyrB|PyrC|PyrD|PyrE|PyrF|PyrG|PyrH|PyrI|CPSase|Carbamoyl|aspartate carbamoyltransferase|OTCace|DHO|dihydroorotate|DHOD|DHODB|orotate|OMP|Orotidine|UPRT|Thymidylat|ThyA|ThyX|dUTPase|DUTPase|dCMP_cyt_deam|dCTP|cytidine deaminase|CDA|Cytidylate_kin|UMP_kinase|Thymidylate_kin|Thymidyl|uridine phosphorylase|thymidine phosphorylase|PNP_UDP|YjeF|Dihydroorotase|Amidohydro"])
d("nucleotide_binding", ["GO:0000166"], [r"nucleotide[- ]binding|\bNTP|\bATP|\bGTP|NAD|FAD|FMN|CoA|cyclic nucleotide|cNMP"])
d("cyclic_dinucleotide", ["GO:0035438", "GO:0052621", "GO:0071111", "GO:0061501", "GO:0106408"], [r"cyclic[- ]di|c-di-GMP|c-di-AMP|cdGMP|cGAMP|GGDEF|EAL|HD-GYP|PilZ|DAC|DisA|CdaA|GdpP|Diguanylate|diguanylate|cyclic nucleotide|CBASS|CD-NTase|SMODS|cGAS|TIR|Cap[0-9]|Pycsar|CAP12|NucC|Cyclase|CARF|SAVED|Card1|Csm6|Csx1"])
d("diguanylate_cyclase", ["GO:0052621"], [r"diguanylate cyclase|GGDEF"])
d("cdgmp_binding", ["GO:0035438"], [r"c-di-GMP|cyclic-di-GMP|cyclic di-GMP|PilZ|GGDEF|EAL|YcgR|BcsA|MshEN|CdgR"])
d("phosphoprotein_phosphatase", ["GO:0004721", "GO:0004722", "GO:0004725"], [r"protein phosphatase|phosphoprotein phosphatase|PP2C|SpoIIE|PTP|DSPc|tyrosine phosphatase|Ser/Thr phosphatase|serine/threonine phosphatase|Metallophos"])
d("serine_threonine_kinase", ["GO:0004674"], [r"serine/threonine|Ser/Thr|Pkinase|protein kinase|RIO|APH|Kinase-like|ABC1|Kdo|Rio1|RIO1|HipA"])
d("tyrosine_kinase", ["GO:0004713"], [r"tyrosine kinase|Pkinase_Tyr|PK_Tyr|Wzc|BY-kinase|CapB"])
d("protein_modification", ["GO:0036211", "GO:0006464", "GO:0043687", "GO:0018193"], [r"modification|kinase|phosphatase|methyltransferase|acetyltransferase|ADP-ribosyl|glycosyltransferase|lipoyl|biotin|ubiquitin|SUMO|sortase|Sortase|transglutaminase|Transglut|lipoprotein|Lgt|Lnt|LspA|prolyl|PPIase|disulfide|Dsb|arginine methyltransferase|PRMT|protein[- ]arginine|protein[- ]lysine|isoaspartyl|PCMT|PIMT|Fic|Doc|FIC|AMPylat|adenylyltransferase|YdiU|SelO|Protein-L-isoaspartate|TIG|ribosomal protein.*(acetyl|methyl)|RimI|RimJ|RimL|PrmA|PrmB|PrmC|ThiF|E1|ubiquitin-activating|MoeB|UBA|LipA|LipB|BirA|BPL|SelR|MsrA|MsrB"])
d("ubiquitin_like", ["GO:0031386"], [r"ubiquitin|Ubiquitin|SUMO|ThiS|MoaD|SAMP|Urm1|UBL|ubl|Ub-like|ubiquitin-like"])
d("ubiquitin_activation", ["GO:0004839", "GO:0008641"], [r"ubiquitin[- ]activating|E1|ThiF|MoeB|UBA|UBACT|E1-like"])
d("ubiquitin_ligase", ["GO:0061630", "GO:0004842"], [r"ubiquitin[- ](protein )?ligase|E3|RING|HECT|U-box|UBOX|NEL|SspH|IpaH"])
d("proteasome", ["GO:0000502", "GO:0005839", "GO:0010498"], [r"proteasome|Proteasome|PAN|ARC|Mpa|PrcA|PrcB|HslV|Pup|PafA|Dop"])
d("lon_protease", ["GO:0004176", "GO:0004252"], [r"\bLon\b|Lon_|LON|Lon protease"])

# ---- Defense / MGE / phage ------------------------------------------------
d("restriction_modification", ["GO:0009307", "GO:0009007", "GO:0008170", "GO:0003886", "GO:0015667", "GO:0009036"], [r"restriction|Restriction|modification methylase|RM\b|R-M|\bR\.|Mrr|HsdR|HsdS|HsdM|Hsd|EcoRI|EcoRII|EcoRV|Eco57I|EcoEI|EcoR124|DpnI|DpnII|BsaWI|TaqI|Tsp45I|AlwI|CfrBI|TdeIII|Bpu10I|MjaI|XcyI|HinP1I|NotI|ThaI|EcoO109I|GmrSD|Methylase_S|N6_N4_Mtase|N6_Mtase|DNA_methylase|MethylTransf|DNA methylase|DNA adenine methylase|Dam|Dcm|type I restriction|type II restriction|type III restriction|type IV restriction|McrA|McrB|McrBC|Mcr|Vsr|RE_|Endonuc-EcoRV|R-HINP1I|T1R|T1RH|ResIII|HSDR|Res_|DpnI_C|Mod"])
d("restriction_enzyme", ["GO:0009036", "GO:0015666"], [r"restriction (endonuclease|enzyme)|Restriction endonuclease|RE_|EcoRI|EcoRII|EcoRV|Endonuc-EcoRV|Eco57I|EcoEI|EcoR124|DpnII|DpnI|BsaWI|TaqI|Tsp45I|AlwI|CfrBI|TdeIII|Bpu10I|MjaI|XcyI|HinP1I|R-HINP1I|NotI|ThaI|EcoO109I|GmrSD|HsdR|HSDR|T1R|T1RH|ResIII|Mrr|McrA|McrB|type (I|II|III|IV) restriction|NA-iREase"])
d("methyltransferase_rm", ["GO:0009007", "GO:0009008", "GO:0003886", "GO:0009307"], [r"modification methylase|restriction[- ]modification|DNA (adenine |cytosine )?methyl|N6_N4_Mtase|N6_Mtase|DNA_methylase|Methylase_S|HsdM|HsdS|MethylTransf|Eco57I|TaqI|Dam|Dcm|M\.|DNA methylase"])
d("dna_methylase", ["GO:0009007", "GO:0009008", "GO:0003886", "GO:0006306"], [r"DNA (adenine |cytosine )?methyl|DNA methylase|DNA_methylase|N6_N4_Mtase|N6_Mtase|Dam|Dcm|modification methylase|C-5 cytosine|N-6 DNA|N-4|MethylTransf"])
d("dna_modification", ["GO:0006304", "GO:0006306"], [r"DNA modification|DNA methyl|Dnd|phosphorothioate|glucosyl|hydroxymethyl|DarT|DarG|ADP-ribosyl"])
d("crispr_associated", ["GO:0043571", "GO:0051607", "GO:0099048"], [r"CRISPR|Cas\d+(?![a-z])|Cas_|Csa[0-9]|Csc[0-9]|Cse[0-9]|Csm[0-9]|Csn[0-9]|Csx[0-9]|Cmr[0-9]|Cst[0-9]|Csd[0-9]|Csy|Csf|Csb|RAMP|DevR|CARF|Card1|Cas_Cas|CRISPR_|Cas12|Cas13|Cas9|Cse1|Csm1|Csm4|Csm2|CT1975|TM1802|Cas2CT1978|Cas4|DUF83"])
d("cas_domain", ["GO:0043571"], [r"CRISPR|Cas\d+(?![a-z])|Cas_|Csa|Csc|Cse|Csm|Csn|Csx|Cmr|Cst|Csd|RAMP|CT1975|TM1802|DUF83"])
d("cas_nuclease", ["GO:0004519", "GO:0043571"], [r"Cas\d+(?![a-z]).*(nuclease|endonuclease)|CRISPR.*(nuclease|endonuclease)|Cas3|Cas9|Cas12|Cas13|Csn1|Cas2|Cas1|Cas4"])
d("antitoxin_domain", ["GO:0110001", "GO:0097351"], [r"antitoxin|Antitoxin|antidote|Antidote|PhdYeFM|ParD|MazE|RelB|HicB|HigA|YefM|DinJ|VapB|MqsA|CcdA|AbiEi|SocA|Panacea|PaRep2a|ImmA|DUF.*antitox|Epsilon|HipB|GhoS|Kis"])
d("antitoxin", ["GO:0110001", "GO:0097351"], [r"antitoxin|Antitoxin|antidote"])
d("toxin", ["GO:0090729"], [r"toxin|Toxin|colicin|bacteriocin|hemolysin|Hemolysin|haemolysin|cytolysin|RTX|Zonular|Zot|SpvB|binary toxin|Binary_toxB|TumE|DarT|enterotoxin|exotoxin|Tox-"])
d("toxin_antitoxin", ["GO:0110001", "GO:0090729", "GO:0097351"], [r"toxin[- ]antitoxin|toxin|antitoxin"])
d("abortive_infection", ["GO:0051607", "GO:0099046"], [r"abortive|Abi[A-Z]?|AbiE|AbiG|AbiD|AbiF|AbiH|AbiJ|AbiK|AbiL|AbiP|AbiQ|AbiR|AbiT|AbiU|AbiV|AbiZ|Lit|PrrC|RexB"])
d("immune_related", ["GO:0006955", "GO:0045087"], [r"immun|immune|Immune"])
d("anti_defense", [], [r"anti-?CRISPR|Acr|anti-restriction|Ocr|ArdA|anti-defense"])
d("anti_crispr", [], [r"anti-?CRISPR|Acr[A-Z]"])
d("anti_restriction", [], [r"anti-restriction|Ocr|ArdA|ArdB|KlcA"])
d("phage_related", ["GO:0019028", "GO:0044423", "GO:0019069", "GO:0019058", "GO:0098003"], [r"phage|Phage|prophage|virus|viral|Virus|capsid|Capsid|portal|Portal|terminase|Terminase|tail|Tail|baseplate|Baseplate|holin|Holin|lysin|Lysin|endolysin|spanin|Rz|Gp[0-9]|gp[0-9]|head|Head|integrase|excisionase|repressor.*phage|CI|Cro|Mu|P2|P4|lambda|T4|T7|HK97|Siphovirus|Myovirus|Podovirus|Caudovirales|DUF.*phage|Phage_"])
d("phage_capsid", ["GO:0019028", "GO:0046798"], [r"capsid|Capsid|coat protein|major head|head protein|HK97|Phage_cap|Phage_capsid"])
d("phage_portal", ["GO:0046798", "GO:0019068"], [r"portal|Portal|Phage_portal"])
d("phage_terminase", ["GO:0019069", "GO:0004519"], [r"terminase|Terminase|Phage_term|DNA packaging"])
d("phage_tail", ["GO:0098003", "GO:0098015", "GO:0098024"], [r"tail|Tail|Phage_tail|sheath|tube|tape measure|Tape_meas"])
d("phage_baseplate", ["GO:0098025"], [r"baseplate|Baseplate|Phage_base|Gp5|gp5|GpJ|GpI|Phage_GPD"])
d("phage_lysin", ["GO:0003796", "GO:0019835"], [r"lysin|Lysin|endolysin|Endolysin|lysozyme|Lysozyme|amidase|Amidase_2|holin|Holin|spanin|Rz"])
d("holin", ["GO:0019835"], [r"holin|Holin|Phage_holin"])
d("viral_capsid", ["GO:0019028", "GO:0005198"], [r"capsid|Capsid|coat protein|Baculo_VP39|major capsid|viral structural|VP[0-9]"])
d("viral_structure", ["GO:0019028", "GO:0044423", "GO:0005198"], [r"capsid|Capsid|coat|virion|viral|virus|Baculo_VP|VP[0-9]|tail|Tail|portal|Portal|structural protein"])
d("receptor_binding", ["GO:0046789", "GO:0019062"], [r"receptor[- ]binding|receptor binding|RBP|A2M_recep|tail fiber|Tail_fiber|attachment"])
d("host_attachment", [], [r"attachment|adhesin|tail fiber"])

# ---- Structural / repeats / misc -----------------------------------------
d("tpr_repeat", [], [r"TPR|Tetratricopeptide|tetratricopeptide|TPR_|Sel1|SEL1|PPR|HAT|Suf"])
d("wd40_repeat", [], [r"WD40|WD\b|WD domain|WD_|WD-40|beta-propeller|YVTN|PQQ|NHL|RCC1"])
d("ankyrin_repeat", [], [r"ankyrin|Ank\b|Ank_|Ankyrin"])
d("lrr_repeat", [], [r"leucine[- ]rich|LRR"])
d("kelch_repeat", [], [r"kelch|Kelch"])
d("heat_repeat", [], [r"HEAT|Armadillo|Arm|HEAT_"])
d("sel1_repeat", [], [r"Sel1|SEL1"])
d("beta_helix", [], [r"beta[- ]helix|Beta_helix|PbH1|pectate lyase|Pectate_lyase|Hexapep|hexapeptide|LbH|pentapeptide"])
d("binding", ["GO:0005488"], [r"binding|bind|Bind"])
d("coiled_coil", [], [r"coiled[- ]coil|Coiled|CC\b"])
d("alpha_helical", [], [r"helical|helix|Helix"])
d("unannotated", [], [r"\bDUF[0-9]|unknown function|uncharacteri[sz]ed|hypothetical"])
d("enzyme", ["GO:0003824"], [r"ase\b|ase |enzyme|synthase|synthetase|kinase|reductase|dehydrogenase|transferase|hydrolase|isomerase|lyase|ligase|oxidase|mutase|epimerase|racemase|deaminase|phosphatase|esterase|peptidase|protease|nuclease|helicase|polymerase|cyclase|carboxylase|decarboxylase|aldolase|hydratase|dehydratase|Glyco_|Aminotran|Methyltransf|Acetyltransf|Radical_SAM|ATP-grasp|AMP-binding|Amidohydro|Abhydrolase|Hydrolase|HAD|NUDIX|TIM|Epimerase|adh_|ADH_|DAO|FAD_|Pyr_redox|Oxidored"])
d("murein_synthesis", ["GO:0009252"], [r"Mur[A-J]|Mur_ligase|peptidoglycan (biosynth|synth)|MraY|FtsW|RodA|SEDS|PBP|Transpeptidase|transglycosylase|Transgly|lipid II|Ddl|Dala_Dala|UDP-N-acetylmuramate|UDP-N-acetylglucosamine enolpyruvyl"])
d("pbp", ["GO:0008658", "GO:0009002", "GO:0008955"], [r"penicillin[- ]binding|PBP|Transpeptidase|transpeptidase|PASTA|Peptidase_S11|Peptidase_S13|carboxypeptidase|Beta-lactamase"])
d("transpeptidase", ["GO:0008658", "GO:0008955", "GO:0071972"], [r"transpeptidase|Transpeptidase|YkuD|L,D-transpeptidase|Ldt|Sortase|sortase|PBP"])
d("lytic_transglycosylase", ["GO:0008933", "GO:0016998"], [r"lytic transglycosylase|transglycosylase|Transglycosylase|SLT|MltA|MltB|MltC|MltD|MltE|MltF|MltG|Lytic_trans|3D"])
d("lps_biosynthesis", ["GO:0009103", "GO:0009244", "GO:0009245", "GO:0008653"], [r"lipopolysaccharide|LPS|lipid A|Lpx[A-M]|LpxA|LpxB|LpxC|LpxD|LpxH|LpxK|LpxL|LpxM|Kdo|KdsA|KdsB|KdsC|KdsD|WaaA|WaaC|WaaF|WaaG|WaaL|Waa|Rfa|O-antigen|O_antigen|Wzx|Wzy|Wzz|Wbp|Wbb|Rml|heptose|Hep|GmhA|GmhB|HldE|LptA|LptB|LptC|LptD|LptE|LptF|LptG|LPS_assembly|Lipid_A|glycosyltransferase.*(LPS|lipopolysaccharide)|LPG_synthase|lysylphosphatidylglycerol"])
d("lipid_a_biosynthesis", ["GO:0009245"], [r"lipid A|Lpx[A-M]|LpxA|LpxC|LpxD|LpxH|LpxK"])
d("o_antigen", ["GO:0009243"], [r"O-antigen|O_antigen|Wzx|Wzy|Wzz|Wbp|Wbb|Rfb"])
d("exopolysaccharide", ["GO:0000271", "GO:0045226"], [r"exopolysaccharide|EPS|capsul|Wza|Wzc|Wzb|polysaccharide (biosynth|export|synth)|Polysacc_synt|PolySac|Wzx|Wzy|Wzz"])
d("capsule", ["GO:0045227"], [r"capsul|Cap[A-D]|Kps|Wza|Wzc|Wzi"])
d("s_layer", ["GO:0030115"], [r"S-layer|S_layer|SLH|surface layer"])
d("vesicle_trafficking", ["GO:0016192", "GO:0006886", "GO:0006888", "GO:0048193"], [r"vesicle|Vesicle|trafficking|SNARE|Sec[0-9]|Vps|VPS|ESCRT|coatomer|COP|clathrin|Clathrin|adaptin|Adaptin|Rab|Ypt|Sar1|Arf|dynamin|Dynamin|exocyst|Exocyst|Exo[0-9]|Exo84|VPS51|TRAPP|Golgi|endosom|Endosom"])
d("large", [], [])
d("gtpase", ["GO:0003924"], [r"GTPase|GTP_EFTU|MMR_HSR1|Ras|Roc|Arf|Gtr1_RagA|Era|FeoB|Dynamin|Septin|Tubulin|FtsZ|Obg|EngA|IIGP|Rab|Ypt|Sar1|Gtr1|Rho|RIO|YihA|YchF|RsgA|HflX|TrmE|MnmE|Der|LepA|SelB|IF2|EF-Tu|EF-G|EF_G|GTP-binding|SRP54|FtsY|FlhF|MinD|MipZ|ParA|CbiA|Mrp|NUBPL|ApbC|YjeQ|Nog|Lsg|Rbg|Drg|Gem|Rag|Miro|Roc|COR|GBP|IRG|Mx"])
d("ras_family", [], [r"\bRas\b|Ras_|Roc|Rab|Rho|Ran|Arf|Sar|Gtr1_RagA|Miro"])
d("dehalogenase", ["GO:0019120"], [r"dehalogenase|HAD|haloacid|Haloacid|reductive dehalogenase|RdhA"])
d("halogenase", ["GO:0016491"], [r"halogenase|Trp_halogenase|haloperoxidase|chloroperoxidase|Chloroperoxidase|brominase|flavin-dependent halogen"])
d("bioluminescence", ["GO:0008218"], [r"luciferase|Luciferase|Lux[A-G]|bioluminescen|LuxC|LuxD|LuxE|LuxG"])
d("tautomerase", [], [r"tautomerase|Tautomerase"])
d("glyoxalase", ["GO:0004462", "GO:0004416"], [r"glyoxalase|Glyoxalase|lactoylglutathione"])
d("polyhydroxybutyrate", ["GO:0042619", "GO:0042621"], [r"polyhydroxy|PHB|Pha[A-Z]|PhaC|PhaZ|PhaR|PhaP|Esterase_PHB|poly\(3-hydroxy|poly-beta-hydroxy|PHA"])
d("dna_gyrase_inhibitor", [], [r"gyrase inhibitor|GyrI|SbmC|ParE|CcdB|Qnr|YacG|Pentapeptide"])
d("metal_homeostasis", ["GO:0055065", "GO:0006879", "GO:0055076", "GO:0030003", "GO:0046916"], [r"homeostasis|metal|iron|ferritin|Ferritin|Bacterioferritin|Dps|DPS|Fur|DtxR|NikR|ZntR|CueR|CsoR|ArsR|CopZ|CopY|CopA|HMA|metallochaperone|ferric uptake|siderophore|TonB|Fe_dep_repress|MerR|YciI"])
d("mercury_resistance", ["GO:0046689", "GO:0050787", "GO:0016152", "GO:0018836"], [r"mercur|Mer[A-TR]|MerA|MerB|MerC|MerD|MerE|MerF|MerP|MerR|MerT|organomercur"])
d("copper_resistance", ["GO:0046688"], [r"copper resistance|Cop[A-DZ]|Cus[ABCF]|CueO|Pco|CopC|CopD"])
d("zinc_resistance", [], [r"zinc resistance|Czc|ZntA|zinc efflux"])
d("antibiotic_resistance", ["GO:0046677", "GO:0008800", "GO:0017001", "GO:0042910", "GO:0071236"], [r"antibiotic|resistan|lactamase|Lactamase|Beta-lactamase|beta-lactamase|chloramphenicol|tetracycline|Tet[A-Z]|TetR|aminoglycoside|Aminoglyc|streptomycin|kanamycin|gentamicin|neomycin|vancomycin|Van[A-Z]|VanZ|bleomycin|Bleomycin|fosfomycin|Fos[A-C]|macrolide|Erm|Mph|Msr|lincosamide|Lnu|quinolone|Qnr|rifampin|Rif|sulfonamide|Sul[1-3]|trimethoprim|Dfr|bacitracin|BacA|polymyxin|colistin|Mcr|Arn|multidrug|Multi_Drug|efflux|Efflux|ACR_tran|AcrB|MatE|MATE|EmrE|SMR|NorM|penicillinase|Penicillinase|BlaI|BlaR|MecI|MecA|MecR|CPT|Chloramphenicol|APH|aminoglycoside phosphotransferase|ANT|AAC|Aminoglyc_resit|Aminoglycoside"])
d("beta_lactamase", ["GO:0008800", "GO:0030655", "GO:0017001"], [r"beta[- ]lactamase|Beta-lactamase|Lactamase|lactamase|penicillinase|carbapenemase|cephalosporinase|ODP|OXA|TEM|SHV|CTX|KPC|NDM|VIM|IMP|AmpC"])
d("aminoglycoside_resistance", ["GO:0046677"], [r"aminoglycoside|Aminoglyc|streptomycin|kanamycin|gentamicin|neomycin|tobramycin|amikacin|spectinomycin|APH.*aminoglycoside|AAC|ANT|Aminoglyc_resit"])
d("multidrug_resistance", ["GO:0042910", "GO:0015562", "GO:0046677"], [r"multidrug|multi-drug|Multi_Drug|ACR_tran|AcrB|AcrD|AcrF|MatE|MATE|EmrE|SMR|NorM|drug resistance|drug efflux|MdtK|MdfA|EmrB|EmrD|QacA|QacE"])
d("detox", [], [])
d("osmoprotectant", ["GO:0006970"], [r"osmoprotect|glycine betaine|betaine|choline|ectoine|proline betaine|OpuA|OpuC|ProU|BetA|BetB|BetT|EctA|EctB|EctC"])

# ---- Cell / physiology ----------------------------------------------------
d("respiration", ["GO:0009060", "GO:0009061", "GO:0022900", "GO:0045333"], [])  # merged with earlier
d("nitrate_transporter", ["GO:0015706", "GO:0015112"], [r"nitrate transport|NarK|NarU|NRT|Nitrate_transp|NasA"])
d("cytoplasmic", ["GO:0005737"], [r"cytoplasm|cytosol"])
d("ion_transporter", ["GO:0015075", "GO:0006811", "GO:0034220"], [])  # merged
d("synthase", [], [r"synthase"])
d("translocase", ["GO:0065002", "GO:0015450", "GO:0043952"], [r"translocase|translocator|TatC|SecY|SecA|YidC|Oxa1|flippase|MurJ|MviN|FtsW|RodA|scramblase|translocation"])
d("urea_metabolism", ["GO:0019627", "GO:0043419", "GO:0009039"], [r"urea|urease|Ure[A-GJ]|arginase|Arginase|Urea|allophanate|urea carboxylase|UAH"])
d("energy_metabolism", ["GO:0006091", "GO:0015980", "GO:0022900", "GO:0006119", "GO:0015986"], [])  # merged
d("cbs_domain", [], [r"\bCBS\b|CBS domain"])
d("chelatase", ["GO:0051002", "GO:0016851", "GO:0004325"], [r"chelatase"])
d("effector_domain", [], [r"effector"])
d("iron_sulfur_biosynthesis", ["GO:0016226"], [r"iron[- ]sulfur cluster (assembly|biosynth|insertion|scaffold)|Fe-S cluster (assembly|biogenesis)|\bSuf[A-E]\b|\bIsc[ASU]\b|NifU|Nfu|SufE|IscU|ApbC"])
d("nife_group3", [], [r"F420[- ]reducing|coenzyme F420 hydrogenase|NAD-reducing hydrogenase|methyl-viologen-reducing hydrogenase"])
d("nife_group4", [], [r"\bMBH\b|energy[- ]converting hydrogenase|formate hydrogenlyase"])
d("pin_domain", [], [r"\bPIN\b|PIN_|PIN domain|PIN-like"])

d("signaling", ["GO:0007165", "GO:0023052", "GO:0000160", "GO:0035556", "GO:0019932"], [r"signal(ing|l?ing| transduction)|sensor|sensing|second messenger|cyclic[- ]di|c-di-GMP|c-di-AMP|diguanylate|two[- ]component|histidine kinase|response regulator|chemotaxis|chemoreceptor|methyl-accepting|receptor"])

d("regulator", ["GO:0065007", "GO:0006355", "GO:0003700", "GO:0140110", "GO:0019222", "GO:0000160", "GO:0045892", "GO:0045893"], [r"regulat|repressor|activator|transcription(al)? factor|anti-sigma|sigma factor|anti-?terminat|MarR|TetR|LysR|GntR|AraC|LacI|AsnC|Lrp|ArsR|IclR|DeoR|MerR|Crp|Fnr|LuxR|OmpR|NarL|Rrf2|PadR|HxlR|TrmB|CodY|Fur\b|DtxR|NikR|ModE|LexA|BlaI|MecI|NrdR|ArgR|PurR|TrpR|PhoU|antitoxin|anti-toxin"])

d("regulatory", ["GO:0065007", "GO:0019222", "GO:0006355"], [r"regulat|repressor|activator|sensor|anti-sigma|inhibitor"])

d("sensor_kinase", ["GO:0000155", "GO:0004673"], [r"(sensor|histidine) kinase|sensor histidine|HisKA|His_kinase|CheA"])

d("two_component", ["GO:0000160", "GO:0000155"], [r"two[- ]component|response regulator|histidine kinase|HisKA|Response_reg|receiver domain|phosphotransfer"])

d("mcr_complex", ['GO:0050524', 'GO:0015948'], ['methyl[- ]coenzyme M reductase', 'MCR_|MCRA\\b'])

d("methanogenesis", ['GO:0015948', 'GO:0050524', 'GO:0019386', 'GO:0019385', 'GO:0019387'], ['methanogen|methyl[- ]coenzyme M|coenzyme M|heterodisulfide|methanopterin|methanofuran|formylmethanofuran|trimethylamine|methylamine', 'MCR_|MCRA\\b|Mtr[A-H]|MtrH|Fmd|Fwd|Mtd\\b|Mch\\b|FTR\\b|Hdr[A-E]|Frh[ABG]|Mvh[ADG]|Eha\\b|Ehb\\b|MtaA|MtaB|MttB|MtbA|MtmB|CdhC|CdhD|CO_dh'])

d("nickel_binding", ['GO:0016151'], ['nickel|urease|NiFe', 'Ni_|HypA|HypB|UreE|UreG|NikA|NikR|MCR_|CO_dh|CdhC|Ni-'])

d("iron_sulfur", ['GO:0051536', 'GO:0051539', 'GO:0051537', 'GO:0051538'], ['iron[- ]sulfur|iron[- ]sulphur|4fe-4s|2fe-2s|3fe-4s|radical sam|fe-s cluster', 'Fer[24]|Fe-?S\\b|FeS\\b|Rieske|Radical_SAM|NifU|SufE|IscU|DHODB_Fe-S|Fe-S_bind|Fer4_NifH|Molybdop_Fe4S4|SPASM'])

# Ferredoxin is the carrier protein; Fe-S binding domains (Fer2/Fer4, "4Fe-4S binding")
# and enzymes that use ferredoxin as a partner carry iron_sulfur / their activity instead.
d("ferredoxin", [], [r"ferredoxin(?![- ]?(oxidoreductase|reductase|hydrogenase|dependent|thioredoxin|type|:|nadp))"])

d("electron_transport", ['GO:0009055', 'GO:0022900'], ['electron (transfer|transport)|flavodoxin|cytochrom|rubredoxin|nadh.*(dehydrogenase|oxidoreductase)|quinone|quinol|4fe-4s|2fe-2s', 'Fer[24]|Rieske|ETF\\b|Cyt_|Cytochrom|Oxidored_q|Complex1|COX\\d|NQR|Rnf[A-G]|Hdr[A-E]|DsrC|DsrMK|Flavodoxin'])

d("cell_surface", ['GO:0009986', 'GO:0005618', 'GO:0030312', 'GO:0019867', 'GO:0007155', 'GO:0030115', 'GO:0009279'], ['surface|s-layer|cell wall|cell-wall|adhesin|adhesion|invasin|intimin|autotransporter|sortase|cellulosome|cohesin|dockerin|cadherin|haemagglutin|hemagglutinin|flagell|pilin|fimbri|outer membrane', 'SLH\\b|S_layer|LPXTG|Gram_pos_anchor|YadA|Hep_Hag|Por_Secre|T9SS|Choline_bind|CW_binding'])

d("secreted", ['GO:0005576'], ['secreted|extracellular|excreted'])

d("nuclease", ['GO:0004518'], ['nuclease|ribonuclease|deoxyribonuclease|endonuclease|exonuclease|restriction endonuclease|restriction enzyme|holliday junction resolvase|resolvase', 'RNase|DNase|Endonuc|Exonuc|HNH|GIY-YIG|RuvC|NERD|XPG|ERCC4|UvrC|RecJ|Mrr_cat|Cas2\\b|Cas4\\b|Cas_Cas4|Cas3|Cas9|Csn1|Cas12|Cas13|Vsr\\b'])

d("repeat_domain", [], ['repeat|tetratricopeptide|pentapeptide|hexapeptide|ankyrin|leucine[- ]rich|kelch|beta-propeller|beta propeller|propeller', 'TPR\\b|TPR_|WD40|Kelch|HEAT\\b|HEAT_|Ank\\b|Ank_|LRR|Sel1|PPR\\b|Pentapeptide|Hexapep|Arm\\b|RCC1|PD40|YVTN|NHL\\b'])

d("tca_cycle", ['GO:0006099', 'GO:0006101', 'GO:0006104', 'GO:0006105', 'GO:0006106', 'GO:0006107', 'GO:0006108'], ['citrate synth|citrate \\(si\\)|aconit|isocitrate dehydrogenase|2-oxoglutarate|oxoglutarate|alpha-ketoglutarate|succinyl-coa|succinate dehydrogenase|fumarase|fumarate hydratase|fumarate reductase|malate dehydrogenase|tricarboxylic', 'Citrate_synt|Aconitase|AcnX|IDH\\b|Iso_dh|OGDH|OGFOR|SucC|SucD|Succ_CoA|Ligase_CoA|SDH|Sdh[A-D]|Fumerase|Mdh\\b|Transket_pyr'])

d("central_metabolism", ['GO:0006096', 'GO:0006094', 'GO:0006099', 'GO:0006098', 'GO:0006090', 'GO:0006086', 'GO:0006097', 'GO:0006091'], ['glycoly|gluconeogen|citrate synth|aconit|isocitrate|oxoglutarate|succinyl-coa|fumarase|malate dehydrogenase|malic enzyme|oxaloacetate|pyruvate|phosphoenolpyruvate|enolase|glyceraldehyde|phosphoglycerate|triose[- ]?phosphate|fructose[- ]bisphosphate|phosphofructokinase|hexokinase|glucokinase|transketolase|transaldolase|glucose-6-phosphate dehydrogenase|6-phosphogluconate|acetate kinase|phosphotransacetylase|lactate', 'Gp_dh|GAPDH|PGK\\b|PGAM|TIM\\b|Enolase|PfkB|PFK\\b|F_bP|FBPase|Aldolase|Transketolase|Transaldolase|G6PD|6PGD|Citrate_synt|Aconitase|IDH\\b|Iso_dh|E1_dh|POR\\b|POR_N|PFOR|OGFOR|SucC|SucD|Succ_CoA|Ligase_CoA|Ldh_1|Mdh\\b|malic|PEPCK|PEPcase|PPDK|PEP-utilizers|PK\\b|PK_C|PTA_PTB|AckA|AFOR|Fumerase|Lar_N|LarA'])

d("aromatic_aa_metabolism", ['GO:0009072', 'GO:0009073', 'GO:0000162', 'GO:0009094', 'GO:0006571', 'GO:0006568', 'GO:0006570', 'GO:0006559', 'GO:0009423'], ['aromatic amino|tryptophan|tyrosine|phenylalanine|chorismate|prephenate|shikimate|anthranilate|dehydroquinate|dehydroquinase|3-dehydroquinate|dahp|epsp synthase|indole-3-glycerol', 'Trp_synt|PRAI|IGPS|PDT\\b|PDH_[NC]|CM_\\d|Chorismate|Shikimate|SKI\\b|DHQ|DHquinase|DAHP|EPSP|Anth_synt'])

d("amino_acid_metabolism", ['GO:0006520', 'GO:0008652', 'GO:0009063', 'GO:0006418', 'GO:0004812', 'GO:0008483', 'GO:0030170'], ['amino[- ]acid|glutamate|glutamine|aspartate|asparagin|alanine|arginine|argininosuccinate|arginase|ornithine|citrulline|cysteine|glycine|histidine|isoleucine|leucine|lysine|methionine|phenylalanine|proline|pyrroline|serine|threonine|tryptophan|tyrosine|valine|homoserine|homocysteine|diaminopimelate|dipicolinate|shikimate|chorismate|prephenate|anthranilate|aminotransferase|transaminase|glutamine amidotransferase', 'Aminotran|PALP|GATase|AlaDh|Ala_racemase|Asp_DH|Orn_|OCD_Mu|ELFV_dehydrog|Sacchrp_dh|KARI|IlvN|ILVD|LeuA|IPMS|IGPS|PRAI|Trp_synt|PDT\\b|DHDPS|DapB|Dap[A-F]\\b|DAP_epim|Homoserine|Thr_synth|Meth_synt|Met_synt|Cys_Met_Meta|CysK|SerA|P5CR|Gln-synt|Asn_synthase|Arginosuc|ArgZ|Arginase|ASL_C|OTCace|Carbam_trans|CPSase|AA_kinase|Semialdhyde_dh|GCV_|Gcv|GDC-P|ADI\\b|Asparaginase|GlutR|tRNA-synt'])


d("lipid_metabolism", ['GO:0006629', 'GO:0008610', 'GO:0016042', 'GO:0006631', 'GO:0006633', 'GO:0006635', 'GO:0006644', 'GO:0008654', 'GO:0046486'], ['lipid|fatty[- ]acid|acyl[- ]coa|acyl carrier|ketoacyl|enoyl|hydroxyacyl|thiolase|phospholipid|glycerophospholipid|diacylglycerol|cardiolipin|phospholipase|lipase|acyltransferase|glycerol[- ]3[- ]phosphate|sterol|hopanoid|isoprenoid|lipid a|lipopolysaccharide|cyclopropane|desaturase|polyketide|polyhydroxy|phosphatidyl.*(synthase|transferase|decarboxylase|kinase)|cdp-diacylglycerol|cdp-archaeol|archaeol', 'Acyl-CoA|ACP\\b|ACP_|PP-binding|Fab[A-Z]\\b|ECH_|ECH\\b|3HCDH|Thiolase|Acyl_transf|ketoacyl|Ketoacyl|PLDc|PAP2|CDP-OH_P_transf|Pls[BCXY]|Lpx[A-M]|DAGK|NAD_Gly3P|CMAS|FA_desaturase|Lipase|HMG_CoA_synt|CTP_transf_1|CarS'])

d("phospholipid_metabolism", ['GO:0006644', 'GO:0008654', 'GO:0046486', 'GO:0006650'], ['phospholipid|glycerophospho|cdp-diacylglycerol|cardiolipin|phospholipase|lysophospholipid|phosphatidyl.*(synthase|transferase|decarboxylase|kinase|methyltransferase)|diacylglycerol kinase|glycerol[- ]3[- ]phosphate|lysylphosphatidylglycerol', 'CDP-OH_P_transf|Pls[BCXY]|PgsA|PLDc|GDPD|NAD_Gly3P|DAGK|LPG_synthase|PEMT'])

d("one_carbon_metabolism", ['GO:0006730', 'GO:0035999', 'GO:0046653', 'GO:0015948', 'GO:0015947', 'GO:0006760'], ['one[- ]carbon|tetrahydrofolate|formyl|formate|formaldehyde|methylene|methenyl|methanopterin|methanofuran|coenzyme m|methanol|methylamine|glycine cleavage|serine hydroxymethyltransferase|thymidylate synth', 'THF_|FTHFS|5-FTHF|MTHFR|MetF|FolD|GCV_|Gcv|GDC-P|SHMT|Formyl_trans|FmdE|Fmd|Fwd|FTR\\b|Mtr[A-H]|MCR_|Fae\\b|FaeA|Thymidylat'])


d("sulfate_reduction", ['GO:0019420', 'GO:0018551', 'GO:0009973'], ['dissimilatory|sulfate reduc|adenylylsulfate reductase|aps reductase|atp-sulfurylase', 'Dsr[A-Z]\\b|DsrC|DsrMK|Apr[AB]|APS-reductase|ATP-sulfurylase|Qmo[A-C]|Sat\\b'])

d("sulfur_metabolism", ['GO:0006790', 'GO:0000103', 'GO:0019420', 'GO:0019418', 'GO:0070814', 'GO:0016783'], ['sulfur|sulphur|sulfate|sulphate|sulfite|sulphite|sulfide|sulphide|thiosulfate|sulfane|polysulfide|tetrathionate|adenylylsulfate|atp-sulfurylase|sulfurtransferase|rhodanese|cysteine desulfurase', 'PAPS|Cys[CDHIJNKE]\\b|Dsr[A-Z]\\b|DsrC|DsrE|Sox[A-Z]\\b|Sqr\\b|SQR|Fcc|TusA|Tus[A-E]\\b|IscS|SufS|NifS|Cys_desulf|Rhodanese|APS-reductase|ATP-sulfurylase|TauE'])

d("cell_wall", ['GO:0005618', 'GO:0009252', 'GO:0071555', 'GO:0008360', 'GO:0009273', 'GO:0042546'], ['cell wall|peptidoglycan|murein|penicillin[- ]binding|d-alanyl|d-ala|carboxypeptidase|muramoyl|lipid ii|teichoic|pseudomurein|diaminopimelate|undecaprenyl', 'Mur[A-J]\\b|MurJ|Mur_ligase|Dala_Dala|Ddl[AB]?\\b|PBP\\b|PBP_|PBP5|Transpeptidase|Peptidase_S11|Peptidase_S13|VanY|Transgly|MraY|UppP|BacA|Alr\\b|Ala_racemase|Amidase_2|Amidase_3|YkuD|RodA|FtsW'])

d("peptidoglycan", ['GO:0009252', 'GO:0000270', 'GO:0009254', 'GO:0042834', 'GO:0008955'], ['peptidoglycan|murein|penicillin[- ]binding|d-alanyl|carboxypeptidase|muramoyl|lipid ii|transpeptidase', 'Mur[A-J]\\b|MurJ|Mur_ligase|Dala_Dala|PBP\\b|PBP_|PBP5|Transpeptidase|Peptidase_S11|Peptidase_S13|Transgly|MraY|YkuD|Amidase_2|Amidase_3'])

d("membrane", ['GO:0016020', 'GO:0005886'], ['membrane|transmembrane|permease|transporter|channel|porin|symporter|antiporter|exporter|importer|efflux|flippase|translocon|translocase|integral', 'TM\\b|TM_|_TM\\b|TMEM|MFS|Mem_trans|OmpA|MotA|MotB|ExbB|ExbD|TolQ|TolR|MscS|MscL|MS_channel|CorA|MgtE|ZIP\\b|Zip\\b|Nramp|FTR1|FeoB|Cation_ATPase|E1-E2|ATP-synt_[ABC]\\b|ATP-synt_C|COX\\d|Cyt_bd|CcmB|CcdA|FtsX|MacB|Sec[DEFGY]\\b|SecD|SecY|Sec61|YidC|Oxa1|Tat[ABC]\\b|TatA|TatC|Rhomboid|DoxX|DedA|VKOR|EamA|MatE|Na_H_Exchanger|Proton_antipo|MNHE|MnhB|PhaG_MnhG|CDP-OH_P_transf|PAP2|LPG_synthase|Complex1_\\d+kDa|NQR|Rnf[A-G]\\b|PTS_EIIC'])

d("atpase", ['GO:0016887'], ['atpase|atp-ase|atp hydrolysis', 'AAA\\b|AAA_|AAA\\+|ABC_tran|ATPase|P-loop'])

d("arsenic_resistance", ['GO:0046685', 'GO:0015446', 'GO:0008794', 'GO:0015105', 'GO:0030612'], ['arsen', 'Ars[BCDHMR]\\b|ACR3|Acr3'])

d("heavy_metal_resistance", ['GO:0046686', 'GO:0010038', 'GO:0046690', 'GO:0046689', 'GO:0046685'], ['heavy[- ]metal|cadmium|mercur|arsen|chromate|tellur|silver|antimon|metal resistance|metal tolerance|divalent (ion|cation) tolerance', 'Cad[AC]\\b|Czc[A-D]\\b|Cop[A-DZ]\\b|Cus[ABCF]\\b|CueO|Pco[A-E]\\b|Mer[A-TR]\\b|Ars[BCDHMR]\\b|ChrA|ChrR|Ter[A-Z]\\b|TerY|Sil[A-P]\\b|ZntA|HMA\\b|CDF\\b|Cation_efflux'])

d("chaperone", ['GO:0006457', 'GO:0051082', 'GO:0044183', 'GO:0140662', 'GO:0003755'], ['chaperon|heat[- ]shock|chaperonin|prefoldin|trigger factor|peptidyl-prolyl|rotamase|ppiase|cyclophilin|foldase|protein folding', 'HSP\\d|Hsp\\d|HSP20|HSP70|HSP90|HSP33|DnaK|DnaJ|GroEL|GroES|Cpn\\d|TCP-1|HtpG|ClpB|Clp_N|IbpA|FKBP|Cyclophil|Pro_isomerase|Rotamase|SurA|SecB\\b|AHSA1'])

d("stress_response", ['GO:0006950', 'GO:0006979', 'GO:0009408', 'GO:0009409', 'GO:0006970'], ['stress|shock|universal stress|oxidative|peroxid|catalase|superoxide|osmotic|starvation|tolerance|resistance|chaperon|methionine sulfoxide|phage shock', 'HSP\\d|Hsp\\d|HSP20|HSP70|Usp\\b|UspA|Dps\\b|OsmC|PspA|PspC|Csp[A-E]?\\b|CSD\\b|DnaK|DnaJ|GroEL|Cpn\\d|ClpB|IbpA|MsrA|MsrB|SelR|Ohr|Sod_|RpoS|RpoH|AHSA1'])

d("dna_binding", ['GO:0003677', 'GO:0043565', 'GO:0003700'], ['dna[- ]binding|helix[- ]turn[- ]helix|winged helix|zinc[- ]?finger|histone|ribbon-helix-helix|transcription(al)? (factor|regulator|repressor|activator)|repressor|dna methylase|dna[- ]methyltransferase|recombinase|integrase|resolvase|single[- ]strand(ed)?[- ]dna|double[- ]strand(ed)?[- ]dna|dna polymerase|dna primase|topoisomerase|gyrase|dna ligase|dna glycosylase|nucleoid|chromosomal protein|chromatin', 'HTH|wHTH|zf-|RHH\\b|RHH_|HHH\\b|HHH_|HhH|MerR|TetR|LysR|LacI|GntR|AraC|MarR|ArsR|AsnC|Lrp\\b|IclR|DeoR|Crp\\b|Fnr\\b|LuxR|GerE|OmpR|Trans_reg|Sigma70|PadR|HxlR|Rrf2|CopG|MetJ|Arc\\b|AbrB|SpoVT|CBFD_NFYB|TBP\\b|TFIIB|Bac_DNA_binding|HU\\b|IHF|H-NS|Alba\\b|MC1\\b|SSB\\b|DnaA|DNA_methylase|N6_N4_Mtase|N6_Mtase|Methylase_S|HsdS|Phage_int|Resolvase|Myb|HLH\\b'])

d("protein_binding", ['GO:0005515'], ['protein[- ]protein interaction|protein[- ]binding|interacting|scaffold|adaptor'])

d("hypothetical", [], ['unknown function|uncharacteri[sz]ed|hypothetical|unknown', 'DUF\\d|UPF\\d'])

d("defense_system", ['GO:0051607', 'GO:0099046', 'GO:0006952', 'GO:0009307', 'GO:0043571'], ['defen[cs]e|anti[- ]?phage|antiviral|abortive infection|restriction|crispr|abortive|gasdermin|viperin|argonaute|retron|cbass|pycsar|immunity protein|bacteriocin immunity', 'CRISPR|Cas\\d+(?![a-z])|Cas_|Csa\\d|Csc\\d|Cse\\d|Csm\\d|Csn\\d|Csx\\d|Cmr\\d|Cst\\d|Csd\\d|RAMP|DevR|CARF|Csm[124]|Card1|Abi[A-Z]\\b|AbiE|RE_|Eco57I|EcoRI|EcoRII|EcoRV|EcoEI|EcoR124|DpnI|DpnII|BsaWI|TaqI|Tsp45I|AlwI|CfrBI|TdeIII|Bpu10I|MjaI|XcyI|HinP1I|R-HINP1I|NotI|ThaI|EcoO109I|GmrSD|HsdR|HSDR|T1R|ResIII|Mrr|McrBC|Methylase_S|DpnI_C|BREX|Brx[A-Z]|Pgl[XZ]|Dnd[A-E]|Gabija|Thoeris|ThsA|ThsB|Septu|PtuA|PtuB|Lamassu|Zorya|Hachiman|Kiwa|Wadjet|Jet[A-D]|Druantia|Shedu|Retron|CAP12|Pycsar|Avs\\d|SIR2|Sir2|GSDM|pAgo|PIWI|Mokosh|Ceres|Borvo|Nhi\\b|Shango|Menshen|Olokun|AIPR|Tiamat|PrrC|RloC|WYL|Stealth'])

d("toxin_domain", ['GO:0090729', 'GO:0004540'], ['toxin|colicin|bacteriocin|polymorphic toxin|ribosome inactivating', 'PemK|MazF|RelE|ParE|HigB|HicA|YafQ|YoeB|Doc\\b|VapC|Zeta|HipA|MqsR|GhoT|Hok\\b|SymE|CcdB|Kid\\b|Txe|YafO|Gp49|DUF891|RHS|Rhs|Ntox|Tox-|Zot\\b|SpvB|Binary_tox|TumE|DarT|CptA|RnlA|LsoA|GinA|Tae\\d|Tde\\d|Tge\\d|Tle\\d|Tse\\d'])

d("signal_recognition", ['GO:0048500', 'GO:0006614'], ['signal recognition|signal peptide|signal sequence', 'SRP\\b|SRP_|SRP\\d|Class_IIIsignal|TAT_signal'])

d("rna_processing", ['GO:0006396', 'GO:0016070', 'GO:0008033', 'GO:0006364'], ["rna processing|rna maturation|ribonuclease|splic|maturase|trna processing|rrna processing|exosome|polyadenylat|poly\\(a\\)|decapping|capping|rna ligase|rna 3'-terminal phosphate cyclase|rna cyclase|intron", 'RNase|Rnase|RNAse|Rrp\\d|RRP\\d|RNase_PH|PNPase|Nob1|RIO\\d|Fibrillarin|Nop\\d|Brix|Pop\\d|Rpp\\d|YbeY|RnpA|RNase_P|CCA\\b|CAA_C|RtcB|RtcA|RTC\\b|RTC_|NYN\\b|Mut7|DcpS|ECR1|CLP1|MTPAP|PolyA_pol'])

d("replication", ['GO:0006260', 'GO:0006261', 'GO:0006275', 'GO:0003896', 'GO:0003887'], ['replicat|dna polymerase|dna primase|primase|replication origin|origin|initiator|helicase loader|sliding clamp|clamp loader|single[- ]strand(ed)?[- ]dna[- ]binding|topoisomerase|gyrase', 'DNA_pol|PolC|Pol_|DnaA|DnaB|DnaC|DnaG|DnaQ|DnaD|DnaI|MCM\\b|MCM_|ORC\\d|Cdc6|GINS|Sld5|PCNA|RFC\\d|Rep_fac|Rad17|SSB\\b|Prim|PriA|PriL|PriS|Topo|Toprim|DNA_gyrase|HolA|HolB|Tus\\b|RepA|Rep_trans|Oap[A-C]'])

d("cell_division", ['GO:0051301', 'GO:0000917', 'GO:0007049', 'GO:0032153'], ['cell division|cell cycle|divisome|septum|septation|septal|cytokinesis|z-ring|anaphase', 'Fts[A-Z]\\b|FtsZ|FtsK|FtsX|FtsI|ZipA|Zap[AB]\\b|Min[CDE]\\b|SepF|MraZ|DivIC|DivIVA|Cdv[ABC]\\b|Septin|SpoIIIE|FtsK_SpoIIIE|MipZ|Tubulin/FtsZ'])

d("chromosome_partitioning", ['GO:0007059', 'GO:0051304', 'GO:0030261'], ['partition|segregation|condensation|condensin|cohesin complex', 'Par[ABM]\\b|ParA|ParB|Spo0J|Soj\\b|SMC\\b|SMC_|Scp[AB]\\b|ScpA|ScpB|Muk[BEF]\\b|KorB|CbiA|MinD|Xer[CD]\\b|Cnd1|MksE'])

d("light_harvesting", ['GO:0030076', 'GO:0009765', 'GO:0030089', 'GO:0016168'], ['light[- ]harvest|antenna|phycobili|phycocyanin|phycoerythrin', 'LHC\\b|PucA|PucB|PufA|PufB|IsiA'])

d("membrane", [], [r"intramembrane"])
d("dna_binding", [], [r"winged[- ]helix"])
d("winged_helix", [], [r"winged[- ]helix"])
d("nucleotide_metabolism", [], [r"nucleotidase|deoxy(cytidine|uridine|adenosine|guanosine|thymidine|ribonucle)", r"Ham1\b|HAM1"])
d("central_metabolism", [], [r"citrate lyase|acetyl-coa synthetase|acetyl-coenzyme a synthetase"])
d("isomerase", [], [r"PPIASE|PPIase"])
d("nuclease", [], [r"endonuc|exonuc", r"RNAse"])
d("rnase", [], [r"RNAse|Rnase"])
d("defense_system", [], [r"conflict system"])
d("biosynthesis", ["GO:0009058"], [r"biosynth"])

d("cofactor_biosynthesis", ['GO:0051188', 'GO:0009108', 'GO:0006783', 'GO:0009236', 'GO:0009228', 'GO:0009231', 'GO:0009435', 'GO:0006777', 'GO:0046656', 'GO:0009234', 'GO:0006744', 'GO:0015940', 'GO:0009102', 'GO:0042823', 'GO:0015937', 'GO:0032324', 'GO:2001118'], ['(biotin|thiamin\\w*|riboflavin|flavin|folate|pterin|molybdopterin|molybdenum cofactor|moco|cobalamin|cobyrinic|cobinamide|corrin|siroheme|heme|haem|porphyrin|lipoate|lipoic acid|pantothen\\w*|pantoate|coenzyme a|ubiquinone|menaquinone|nad\\+?|nicotinate|nicotinamide|f420|coenzyme m|methanopterin|pqq|pyrroloquinoline quinone|mycofactocin|queuosine|tetrahydrofolate|pyridoxal|pyridoxine|pyridoxamine|pyridoxal[- ]phosphate|plp)[\\w\\-/,() ]{0,24}?(biosynth|synthase|synthetase|synthesis|formation)', 'cofactor (biosynth|synth)|coenzyme (biosynth|synth)|gtp cyclohydrolase|dihydropteroate synth|dihydrofolate|uroporphyrinogen|protoporphyrinogen|coproporphyrinogen|precorrin|cobyrinic|adenosylcobinamide|quinolinate|dephospho-coa|phosphopantothenate', 'Bio[A-FW]\\b|Thi[A-HLMOSW]\\b|ThiI\\b|Thi4|TMP-TENI|Rib[A-H]\\b|RibD|Pdx[A-JST]\\b|PdxA|PNPOx|Fol[A-EKP]\\b|DHPS|DHFR|GTP_cyclohydro|Moa[A-E]\\b|Moe[AB]\\b|MoeA|Mob[AB]\\b|MobB|MoCF|MOFRL|Cob[A-Z]\\b|Cbi[A-Z]\\b|Hem[A-NY]\\b|Lip[AB]\\b|LIAS|Pan[B-E]\\b|PanE|ApbA|Coa[A-E]\\b|Ubi[A-JX]\\b|UbiD|Men[A-H]\\b|Cof[A-H]\\b|CofC|Fbi[AB]\\b|Com[A-E]\\b|Pqq[A-G]\\b|PqqD|Que[A-F]\\b|QRPTase|NAD_synthase|NadA|MptE|PTPS|HPPK|PPS_PS|CitG|SOR_SNZ|GCH_III|Phos_pyr_kin'])

d("thiamine_biosynthesis", ['GO:0009228', 'GO:0009229', 'GO:0036172'], ['thiamin\\w*[\\w\\-/,() ]{0,24}?(biosynth|synthase|synthetase|synthesis|formation)', 'thiazole|hydroxymethylpyrimidine|thiamine[- ]phosphate synthase|thiamine monophosphate synthase', 'Thi[A-HLMOSW]\\b|ThiI\\b|Thi4|TenA|TMP-TENI|ThiS\\b|ThiN\\b|Phos_pyr_kin|NMT1|BATS'])

d("cobalamin_biosynthesis", ['GO:0009236', 'GO:0009235'], ['(cobalamin|corrin|cobyrinic|cobinamide|precorrin|b12)[\\w\\-/,() ]{0,24}?(biosynth|synthase|synthetase|synthesis|formation)', 'precorrin|cobyrinic|adenosylcobinamide|cobaltochelatase|cobalamin[- ]5[- ]phosphate synthase', 'Cob[A-Z]\\b|Cbi[A-Z]\\b|CbiG|CbiZ|CobT|CobW|CobS|BtuR|EutT|PduO|SirB'])

d("folate_biosynthesis", ['GO:0046656', 'GO:0046654', 'GO:0009396'], ['(folate|pterin|tetrahydrofolate|dihydrofolate)[\\w\\-/,() ]{0,24}?(biosynth|synthase|synthetase|synthesis|formation)', 'dihydropteroate|dihydrofolate reductase|gtp cyclohydrolase|aminodeoxychorismate|hydroxymethyldihydropterin', 'Fol[A-EKP]\\b|DHPS|DHFR|GTP_cyclohydro|Pab[A-C]\\b|HPPK|DHNA'])

d("nad_biosynthesis", ['GO:0009435', 'GO:0019363', 'GO:0034628'], ['(nad|nicotinate|nicotinamide|quinolinate)[\\w\\-/,() ]{0,24}?(biosynth|synthase|synthetase|synthesis|formation)', 'nad\\+? synth|nad\\(\\+\\) synth|quinolinate|nicotinate[- ]nucleotide|nicotinamide[- ]nucleotide|nmn adenylyltransferase|namn adenylyltransferase|nad kinase', 'Nad[A-EKR]\\b|NAD_synthase|NAD_kinase|NAPRTase|QRPTase|Pnc[AB]\\b'])

d("pyridoxal_biosynthesis", ['GO:0042823', 'GO:0008615'], ['(pyridoxal|pyridoxine|pyridoxamine|plp|vitamin b6)[\\w\\-/,() ]{0,24}?(biosynth|synthase|synthetase|synthesis|formation)', "pyridoxamine 5'-phosphate oxidase|pyridoxine 5'-phosphate oxidase|pyridoxal kinase|pyridoxine kinase", 'Pdx[A-JST]\\b|PdxA|PNPOx|Pyridox_ox'])

d("heme_biosynthesis", ['GO:0006783', 'GO:0006779', 'GO:0033014'], ['(heme|haem|porphyrin|siroheme|tetrapyrrole)[\\w\\-/,() ]{0,24}?(biosynth|synthase|synthetase|synthesis|formation)', 'uroporphyrinogen|protoporphyrinogen|coproporphyrinogen|porphobilinogen|ferrochelatase|glutamyl-trna reductase|aminolevulinate', 'Hem[A-NY]\\b|ALAD|PBGD|UROD|CPOX|PPOX|GlutR|Ahb[A-D]\\b|Nir[DHJ]\\b'])

d("cobalamin_binding", ['GO:0031419'], ['cobalamin|corrinoid|vitamin b12|b12[- ](binding|dependent)', 'B12-binding|B12_binding|CbiX'])

d("cobalt_binding", ['GO:0050897'], ['cobalt|cobalamin|corrinoid|vitamin b12|b12[- ](binding|dependent)', 'Cbi[MNQOX]\\b|CbiX|CbiM|CbiQ|B12-binding'])

d("translation_factor", ['GO:0003743', 'GO:0003746', 'GO:0003747', 'GO:0008135', 'GO:0006413', 'GO:0006414', 'GO:0006415'], ['(?<!transcription )(?<!transcriptional )(?<!transcription-)(initiation|elongation|release|recycling) factor|translation factor|translation (initiation|elongation|termination)', 'eIF|eEF|eRF|IF[123]\\b|IF-[123]|EF-[GPT]|EF-Tu|EF_|GTP_EFTU|RRF\\b|RF-1|RF-3|RF[13]_|SelB|LepA|EFP|SUI1|IF-2B|Hpf|RaiA|EFG_|IF5A|acVLRF1|baeRF'])

d("translation", ['GO:0006412', 'GO:0003735', 'GO:0005840', 'GO:0043039', 'GO:0004812', 'GO:0003743', 'GO:0003746', 'GO:0003747', 'GO:0006414', 'GO:0006413', 'GO:0006415'], ['ribosomal protein|ribosomal subunit|translation|aminoacyl[- ]trna|trna[- ]synthetase|anticodon|(?<!transcription )(?<!transcriptional )(?<!transcription-)(initiation|elongation|release|recycling) factor', 'Ribosomal_|Ribosom_|RL\\d|tRNA-synt|tRNA_synt|eIF|eEF|eRF|IF[123]\\b|IF-[123]|EF-[GPT]|EF-Tu|GTP_EFTU|RRF\\b|RF-1|RF-3|SelB|LepA|EFP|SUI1|IF-2B|Hpf|RaiA|EFG_|IF5A|acVLRF1|baeRF|PTH2?\\b|Peptidyl_tRNA_hyd'])

d("amino_acid_biosynthesis", ['GO:0008652', 'GO:0009073', 'GO:0009082', 'GO:0009085', 'GO:0009089', 'GO:0006526', 'GO:0000105', 'GO:0009086', 'GO:0006571', 'GO:0000162', 'GO:0009094', 'GO:0006561', 'GO:0006564', 'GO:0009097', 'GO:0009098', 'GO:0009099', 'GO:0006535', 'GO:0006537', 'GO:0006542', 'GO:0019344', 'GO:0009067', 'GO:0019877', 'GO:0009423'], ['(amino acid|glutamate|glutamine|aspartate|asparagine|lysine|arginine|histidine|methionine|tryptophan|tyrosine|phenylalanine|leucine|isoleucine|valine|threonine|serine|cysteine|proline|glycine|homoserine|homocysteine|ornithine|diaminopimelate|dihydrodipicolinate|chorismate|shikimate)[\\w\\-/,() ]{0,24}?(biosynth|synthase|synthetase|synthesis|formation)', 'homoserine dehydrogenase|dihydrodipicolinate|argininosuccinate|ornithine carbamoyltransferase|glutamine synthetase|anthranilate synthase|tryptophan synthase|threonine synthase|cysteine synthase|methionine synthase|isopropylmalate|acetohydroxy|acetolactate synthase|dihydroxy-acid dehydratase|ketol-acid|diaminopimelate', 'DapB|Dap[A-F]\\b|DHDPS|LysX|LysW|Arg[A-HJ]\\b|OTCace|His[A-IZ]\\b|IGPD|Trp[A-E]\\b|Trp_synt|Anth_synt|Chorismate_synt|EPSP|DHQ_synth|SKI\\b|DAHP|PDT\\b|ILVD|KARI|IlvN|LeuA|IPMS|Homoserine|Thr_synth|Meth_synt|Met_synt|CysK|SerA|SerB|SerC|Pro[ABC]\\b|P5CR|Gln-synt|Asn_synthase|Arginosuc'])

d("nad_binding", ['GO:0051287', 'GO:0070403', 'GO:0050661'], ['(nad\\(p\\)h?|nad\\(p\\)|nadph?|nadh?)[- ](binding|dependent)|rossmann', 'NAD_binding|NADP_binding|NAD\\(P\\)-binding|NADP-binding|ADH_zinc|ADH_N|Ldh_1_N|adh_short|SDR\\b|2-Hacid_dh_C|3HCDH_N|Semialdhyde_dh|NAD_kinase'])

d("nadp_binding", ['GO:0050661'], ['(nadph?|nad\\(p\\)h?)[- ](binding|dependent)', 'NADP_binding|NADP-binding|NAD\\(P\\)-binding'])

d("cofactor_binding", ['GO:0048037', 'GO:0050662'], ['(nad\\(p\\)h?|nad\\(p\\)|nadph?|nadh?|fad|fmn|flavin|f420|pqq|biotin|lipoyl|lipoate|thiamin\\w*|tpp|pyridoxal[- ]phosphate|pyridoxal|plp|cobalamin|b12|heme|haem|s-adenosylmethionine|coenzyme a|coa|molybdopterin|pterin)[- ](binding|dependent)|cofactor[- ]binding|coenzyme[- ]binding', 'NAD_binding|NADP_binding|FAD_binding|FMN_bind|B12-binding|PQQ\\b|Biotin_lipoyl|Lipoyl|TPP_enzyme|F420_oxidored'])
d("iron_sulfur_biosynthesis", [], [r"SufBD|SufB\b|SufD\b|SufC\b|SufS\b|SufE\b"])
d("mobile_element", [], [r"Tni[ABQ]\b|Tns[A-E]\b"])

# ---- Cofactor classes implied by GO activity definitions --------------------
# GO's "NAD or NADP as acceptor" / "acting on NAD(P)H" oxidoreductase classes
# mirror ENZYME's NAD(P) sub-subclasses; SAM-dependent methyltransferase
# activity requires S-adenosylmethionine; "flavin as acceptor" classes mirror
# ENZYME's flavin sub-subclasses. IDs checked against go-basic.obo.
d("sam_binding", ["GO:0008757"])
d("nad_binding", ["GO:0016616", "GO:0016620", "GO:0016628", "GO:0016639", "GO:0016646", "GO:0016651",
                  "GO:0016680", "GO:0016696", "GO:0016723", "GO:0016726", "GO:0016731", "GO:0046857"])
d("flavin_binding", ["GO:0046997", "GO:0052890"])

# ---- Component-level equivalents of system predicates -----------------------
# Maps emit component-level predicates only; system-level ones (defense_system,
# toxin_antitoxin, abortive_infection, ...) come from system callers. A proposed
# system predicate resolves to its component equivalent (vocabulary
# COMPONENT_EQUIVALENT) and is judged on the specs below.
d("abi_domain", [], list(E["abortive_infection"]["text"]))
d("defense_component", list(E["defense_system"]["go"]) + ["GO:0110001", "GO:0097351"],
  list(E["defense_system"]["text"]) + [r"toxin[- ]antitoxin"])
d("pectinase", ["GO:0004650", "GO:0030570", "GO:0047490", "GO:0030599"])
