"""Tests for the predicate generator system."""

import pytest

from sharur.predicates.generator import (
    PredicateGenerator,
    AnnotationRecord,
    ProteinRecord,
    generate_predicates_for_proteins,
)
from sharur.predicates.vocabulary import (
    get_predicate,
    list_predicates,
    list_categories,
    get_hierarchy,
)
from sharur.predicates.mappings.pfam_map import get_predicates_for_pfam
from sharur.predicates.mappings.kegg_map import get_predicates_for_kegg, get_predicates_for_ec
from sharur.predicates.mappings.cazy_map import get_predicates_for_cazy
from sharur.predicates.mappings.vog_map import get_vog_predicates


class TestVocabulary:
    """Tests for predicate vocabulary."""

    def test_get_predicate_exists(self):
        """Should return predicate definition for existing ID."""
        pred = get_predicate("transporter")
        assert pred is not None
        assert pred.predicate_id == "transporter"
        assert pred.category == "transport"

    def test_get_predicate_not_exists(self):
        """Should return None for non-existent ID."""
        pred = get_predicate("nonexistent_predicate")
        assert pred is None

    def test_list_predicates_all(self):
        """Should return all predicates."""
        preds = list_predicates()
        assert len(preds) > 100  # We defined many predicates

    def test_list_predicates_by_category(self):
        """Should filter predicates by category."""
        transport_preds = list_predicates(category="transport")
        assert len(transport_preds) > 0
        assert all(p.category == "transport" for p in transport_preds)

    def test_list_categories(self):
        """Should return all categories."""
        categories = list_categories()
        assert "transport" in categories
        assert "enzyme" in categories
        assert "metabolism" in categories
        assert "cazy" in categories

    def test_get_hierarchy(self):
        """Should return predicate hierarchy."""
        hierarchy = get_hierarchy("abc_transporter")
        assert "abc_transporter" in hierarchy
        assert "transporter" in hierarchy  # Parent

    def test_get_hierarchy_for_leaf(self):
        """Should return single element for predicate without parent."""
        hierarchy = get_hierarchy("transporter")
        assert "transporter" in hierarchy


class TestPfamMapping:
    """Tests for PFAM mapping."""

    def test_direct_mapping(self):
        """Should map known PFAM domains."""
        preds = get_predicates_for_pfam("PF00005", "ABC_tran", "ABC transporter")
        assert "transporter" in preds
        assert "abc_transporter" in preds
        assert "atp_binding" in preds

    def test_response_regulator(self):
        """Should map response regulator domain."""
        preds = get_predicates_for_pfam("PF00072", "Response_reg", "Response regulator")
        assert "regulator" in preds
        assert "response_regulator" in preds
        assert "two_component" in preds

    def test_pattern_matching(self):
        """Should apply pattern matching for unmapped domains."""
        # Unknown accession but has transport in description
        preds = get_predicates_for_pfam("PF99999", "Unknown", "Some transport protein")
        assert "transporter" in preds

    def test_enzyme_domains(self):
        """Should map enzyme domains."""
        # Dehydrogenase
        preds = get_predicates_for_pfam("PF00106", "adh_short", "short chain dehydrogenase")
        assert "oxidoreductase" in preds
        assert "dehydrogenase" in preds
        assert "nad_binding" in preds

    def test_name_only_mapping(self):
        """Should map PFAM name-only HMMs."""
        preds = get_predicates_for_pfam("Ras", "Ras", "")
        assert "gtpase" in preds
        assert "gtp_binding" in preds

    def test_name_only_transcription(self):
        """Should map transcription factors by PFAM name."""
        preds = get_predicates_for_pfam("TFIIB", "TFIIB", "")
        assert "transcription" in preds
        assert "dna_binding" in preds

    def test_ubiquitin_like_mapping(self):
        """Should map ubiquitin-like domains by PFAM name."""
        preds = get_predicates_for_pfam("ubiquitin", "ubiquitin", "")
        assert "ubiquitin_like" in preds

    def test_histone_mapping(self):
        """Should map histone domains by PFAM name."""
        preds = get_predicates_for_pfam("Histone", "Histone", "")
        assert "histone" in preds

    def test_sliding_clamp_mapping(self):
        """Should map sliding clamp domains."""
        preds = get_predicates_for_pfam("PCNA_N", "PCNA_N", "")
        assert "sliding_clamp" in preds

    def test_clamp_loader_mapping(self):
        """Should map clamp loader domains."""
        preds = get_predicates_for_pfam("Rad17", "Rad17", "")
        assert "clamp_loader" in preds


class TestKeggMapping:
    """Tests for KEGG mapping."""

    def test_direct_mapping(self):
        """Should map known KEGG orthologs."""
        preds = get_predicates_for_kegg("K00532", "hydrogenase large subunit")
        assert "hydrogenase" in preds
        assert "hydrogen_metabolism" in preds

    def test_nitrogenase(self):
        """Should map nitrogenase."""
        preds = get_predicates_for_kegg("K02586", "nitrogenase iron protein NifH [EC:1.18.6.1]")
        assert "nitrogenase" in preds
        assert "nitrogen_fixation" in preds

    def test_ec_extraction(self):
        """Should extract EC numbers from definition."""
        preds = get_predicates_for_kegg("K00001", "alcohol dehydrogenase [EC:1.1.1.1]")
        assert "oxidoreductase" in preds
        assert "dehydrogenase" in preds

    def test_ec_to_predicates(self):
        """Should map EC numbers correctly."""
        # Oxidoreductase
        preds = get_predicates_for_ec("1.1.1.1")
        assert "oxidoreductase" in preds
        assert "dehydrogenase" in preds

        # Hydrolase
        preds = get_predicates_for_ec("3.2.1.4")
        assert "hydrolase" in preds
        assert "glycosidase" in preds

        # Kinase
        preds = get_predicates_for_ec("2.7.1.1")
        assert "transferase" in preds
        assert "kinase" in preds

    def test_direct_dna_methylase(self):
        """Should map DNA methylases by KO."""
        preds = get_predicates_for_kegg("K00558", "")
        assert "methyltransferase" in preds
        assert "dna_methylase" in preds

    def test_direct_crispr(self):
        """Should map CRISPR Cas proteins by KO."""
        preds = get_predicates_for_kegg("K19091", "")
        assert "crispr_associated" in preds
        assert "nuclease" in preds

    def test_direct_primase(self):
        """Should map primase by KO."""
        preds = get_predicates_for_kegg("K02684", "")
        assert "primase" in preds
        assert "replication" in preds

    def test_direct_ubiquitin_ligase(self):
        """Should map ubiquitin ligase by KO."""
        preds = get_predicates_for_kegg("K15343", "")
        assert "ubiquitin_ligase" in preds
        assert "protein_modification" in preds


class TestVogMapping:
    """Tests for VOGdb mapping."""

    def test_ribonucleotide_reductase(self):
        """Should map ribonucleotide reductase descriptions."""
        preds = get_vog_predicates("VOG00000", description="ribonucleotide reductase alpha")
        assert "nucleotide_metabolism" in preds
        assert "oxidoreductase" in preds

    def test_transcription_initiation_factor(self):
        """Should map transcription initiation factors."""
        preds = get_vog_predicates(
            "VOG00001",
            description="transcription initiation factor IIB-like protein",
        )
        assert "transcription" in preds
        assert "transcription_factor" in preds

    def test_leucine_rich_repeat(self):
        """Should map leucine rich repeat proteins."""
        preds = get_vog_predicates(
            "VOG00002",
            description="leucine rich repeat domain containing protein",
        )
        assert "lrr_repeat" in preds

    def test_transposase(self):
        """Should map mobile element proteins."""
        preds = get_vog_predicates("VOG00003", description="transposase")
        assert "transposase" in preds
        assert "mobile_element" in preds

    def test_ubiquitin_patterns(self):
        """Should map ubiquitin-related descriptions."""
        preds = get_vog_predicates("VOG00004", description="E3 ubiquitin ligase")
        assert "ubiquitin_ligase" in preds
        preds = get_vog_predicates("VOG00005", description="ubiquitin protease")
        assert "deubiquitinase" in preds

    def test_histone_patterns(self):
        """Should map histone descriptions."""
        preds = get_vog_predicates("VOG00006", description="histone-like protein")
        assert "histone" in preds


class TestCazyMapping:
    """Tests for CAZy mapping."""

    def test_glycoside_hydrolase(self):
        """Should map GH families."""
        preds = get_predicates_for_cazy("GH5")
        assert "carbohydrate_active" in preds
        assert "glycoside_hydrolase" in preds
        assert "cellulase" in preds

    def test_lpmo(self):
        """Should map LPMO families."""
        preds = get_predicates_for_cazy("AA9")
        assert "carbohydrate_active" in preds
        assert "auxiliary_activity" in preds
        assert "lytic_polysaccharide_monooxygenase" in preds
        assert "copper_binding" in preds

    def test_cbm(self):
        """Should map CBM families."""
        preds = get_predicates_for_cazy("CBM1")
        assert "carbohydrate_active" in preds
        assert "carbohydrate_binding" in preds

    def test_fallback_pattern(self):
        """Should use pattern fallback for unknown families."""
        preds = get_predicates_for_cazy("GH999")
        assert "glycoside_hydrolase" in preds
        assert "carbohydrate_active" in preds


class TestPredicateGenerator:
    """Tests for the predicate generator."""

    def test_basic_generation(self):
        """Should generate predicates from annotations."""
        gen = PredicateGenerator()
        protein = ProteinRecord(protein_id="test", sequence_length=500)
        annotations = [
            AnnotationRecord(
                source="pfam",
                accession="PF00005",
                name="ABC_tran",
                description="ABC transporter",
                evalue=1e-50,
            )
        ]
        preds = gen.generate_for_protein(protein, annotations)
        assert "transporter" in preds
        assert "abc_transporter" in preds
        assert "pfam_annotated" in preds
        assert "confident_hit" in preds

    def test_size_predicates(self):
        """Should add size predicates."""
        gen = PredicateGenerator()

        # Tiny
        protein = ProteinRecord(protein_id="test", sequence_length=30)
        preds = gen.generate_for_protein(protein, [])
        assert "tiny" in preds
        assert "unannotated" in preds

        # Giant
        protein = ProteinRecord(protein_id="test", sequence_length=1500)
        preds = gen.generate_for_protein(protein, [])
        assert "giant" in preds

        # Massive
        protein = ProteinRecord(protein_id="test", sequence_length=3000)
        preds = gen.generate_for_protein(protein, [])
        assert "giant" in preds
        assert "massive" in preds

    def test_unannotated(self):
        """Should mark unannotated proteins."""
        gen = PredicateGenerator()
        protein = ProteinRecord(protein_id="test", sequence_length=300)
        preds = gen.generate_for_protein(protein, [])
        assert "unannotated" in preds

    def test_hypothetical(self):
        """Should detect hypothetical annotations."""
        gen = PredicateGenerator()
        protein = ProteinRecord(protein_id="test", sequence_length=300)
        annotations = [
            AnnotationRecord(
                source="pfam",
                accession="PF12345",
                name="DUF999",
                description="Domain of unknown function",
                evalue=1e-20,
            )
        ]
        preds = gen.generate_for_protein(protein, annotations)
        assert "hypothetical" in preds

    def test_confidence_levels(self):
        """Should set confidence predicates."""
        gen = PredicateGenerator()
        protein = ProteinRecord(protein_id="test", sequence_length=300)

        # Confident hit
        annotations = [
            AnnotationRecord(
                source="pfam",
                accession="PF00005",
                evalue=1e-50,
            )
        ]
        preds = gen.generate_for_protein(protein, annotations)
        assert "confident_hit" in preds

        # Weak hit
        annotations = [
            AnnotationRecord(
                source="pfam",
                accession="PF00005",
                evalue=1e-3,
            )
        ]
        preds = gen.generate_for_protein(protein, annotations)
        assert "weak_hit" in preds

    def test_multi_domain(self):
        """Should detect multi-domain proteins."""
        gen = PredicateGenerator()
        protein = ProteinRecord(protein_id="test", sequence_length=800)
        annotations = [
            AnnotationRecord(source="pfam", accession="PF00001", evalue=1e-20),
            AnnotationRecord(source="pfam", accession="PF00002", evalue=1e-20),
            AnnotationRecord(source="pfam", accession="PF00003", evalue=1e-20),
        ]
        preds = gen.generate_for_protein(protein, annotations)
        assert "multi_domain" in preds
        assert "well_annotated" in preds

    def test_multi_source(self):
        """Should detect multi-source annotations."""
        gen = PredicateGenerator()
        protein = ProteinRecord(protein_id="test", sequence_length=500)
        annotations = [
            AnnotationRecord(source="pfam", accession="PF00005", evalue=1e-20),
            AnnotationRecord(source="kegg", accession="K00001", evalue=1e-20),
        ]
        preds = gen.generate_for_protein(protein, annotations)
        assert "multi_source" in preds
        assert "pfam_annotated" in preds
        assert "kegg_annotated" in preds

    def test_direct_access_predicates(self):
        """Should add direct access predicates when enabled."""
        gen = PredicateGenerator(include_direct_access=True)
        protein = ProteinRecord(protein_id="test", sequence_length=500)
        annotations = [
            AnnotationRecord(source="pfam", accession="PF00005", evalue=1e-20),
            AnnotationRecord(source="kegg", accession="K00001", evalue=1e-20),
            AnnotationRecord(source="cazy", accession="GH5", evalue=1e-20),
        ]
        preds = gen.generate_for_protein(protein, annotations)
        assert "pfam:PF00005" in preds
        assert "kegg:K00001" in preds
        assert "cazy:GH5" in preds

    def test_direct_access_disabled(self):
        """Should not add direct access predicates when disabled."""
        gen = PredicateGenerator(include_direct_access=False)
        protein = ProteinRecord(protein_id="test", sequence_length=500)
        annotations = [
            AnnotationRecord(source="pfam", accession="PF00005", evalue=1e-20),
        ]
        preds = gen.generate_for_protein(protein, annotations)
        assert "pfam:PF00005" not in preds
        assert "transporter" in preds  # Still gets semantic predicates

    def test_hierarchy_expansion(self):
        """Should expand hierarchy when enabled."""
        gen = PredicateGenerator(expand_hierarchy=True)
        protein = ProteinRecord(protein_id="test", sequence_length=500)
        annotations = [
            AnnotationRecord(source="pfam", accession="PF00005", evalue=1e-20),
        ]
        preds = gen.generate_for_protein(protein, annotations)
        # abc_transporter should bring in transporter
        assert "abc_transporter" in preds
        assert "transporter" in preds


class TestGenerateFromDatabase:
    """Tests for database-based predicate generation."""

    def test_generate_predicates_for_proteins(self, store):
        """Should generate predicates from database."""
        predicates = generate_predicates_for_proteins(store)
        assert len(predicates) > 0

        # Check a protein with annotations
        if "prot_001" in predicates:
            preds = predicates["prot_001"]
            assert "pfam_annotated" in preds
            # prot_001 has NiFe-hydrogenase annotation
            assert "hydrogenase" in preds or "confident_hit" in preds

    def test_generate_for_specific_proteins(self, store):
        """Should generate for specific protein subset."""
        predicates = generate_predicates_for_proteins(
            store, protein_ids=["prot_001", "prot_002"]
        )
        assert len(predicates) <= 2
        assert "prot_001" in predicates or len(predicates) == 0

    def test_unannotated_proteins_get_predicates(self, store):
        """Unannotated proteins should get size predicates."""
        predicates = generate_predicates_for_proteins(store)
        # prot_004 is giant and unannotated in test data
        if "prot_004" in predicates:
            preds = predicates["prot_004"]
            assert "giant" in preds


class TestTopologyCategory:
    """Tests for topology predicates in vocabulary."""

    def test_topology_category_exists(self):
        """Should have topology category."""
        categories = list_categories()
        assert "topology" in categories

    def test_topology_predicates(self):
        """Should have topology predicates."""
        preds = list_predicates(category="topology")
        pred_ids = [p.predicate_id for p in preds]
        assert "transmembrane_predicted" in pred_ids
        assert "single_pass_membrane" in pred_ids
        assert "multi_pass_membrane" in pred_ids
        assert "soluble_predicted" in pred_ids

    def test_topology_hierarchy(self):
        """Topology predicates should have correct hierarchy."""
        hierarchy = get_hierarchy("single_pass_membrane")
        assert "single_pass_membrane" in hierarchy
        assert "transmembrane_predicted" in hierarchy

        hierarchy = get_hierarchy("polytopic_membrane")
        assert "polytopic_membrane" in hierarchy
        assert "multi_pass_membrane" in hierarchy
        assert "transmembrane_predicted" in hierarchy


class TestMetalBindingMappings:
    """Tests for metal-binding PFAM mappings."""

    def test_nickel_binding(self):
        """Should map nickel-binding domains."""
        # NiFe hydrogenase
        preds = get_predicates_for_pfam("PF00374", "NiFeSe_Hases", "NiFe hydrogenase")
        assert "nickel_binding" in preds
        assert "metal_binding" in preds

        # Urease
        preds = get_predicates_for_pfam("PF00449", "Urease_alpha", "Urease alpha subunit")
        assert "nickel_binding" in preds

    def test_molybdenum_binding(self):
        """Should map molybdenum-binding domains."""
        preds = get_predicates_for_pfam("PF00384", "Molybdopterin", "Molybdopterin cofactor")
        assert "molybdenum_binding" in preds
        assert "metal_binding" in preds

    def test_copper_binding(self):
        """Should map copper-binding domains."""
        preds = get_predicates_for_pfam("PF00394", "Cu-oxidase_2", "Multi-copper oxidase")
        assert "copper_binding" in preds
        assert "metal_binding" in preds

    def test_cobalt_binding(self):
        """Should map cobalamin-binding domains."""
        preds = get_predicates_for_pfam("PF02310", "B12-binding", "B12 binding domain")
        assert "cobalt_binding" in preds
        assert "cobalamin_binding" in preds
        assert "metal_binding" in preds

    def test_zinc_binding(self):
        """Should map zinc-binding domains."""
        preds = get_predicates_for_pfam("PF00096", "zf-C2H2", "Zinc finger C2H2")
        assert "zinc_binding" in preds
        assert "zinc_finger" in preds
        assert "metal_binding" in preds

    def test_iron_sulfur(self):
        """Should map iron-sulfur cluster domains."""
        preds = get_predicates_for_pfam("PF00037", "Fer4", "4Fe-4S binding domain")
        assert "iron_sulfur" in preds
        assert "metal_binding" in preds

    def test_metal_binding_patterns(self):
        """Should match metal-binding patterns in descriptions."""
        # Nickel pattern
        preds = get_predicates_for_pfam("PF99999", "Unknown", "nickel-dependent enzyme")
        assert "nickel_binding" in preds

        # Copper pattern
        preds = get_predicates_for_pfam("PF99998", "Unknown", "copper oxidase domain")
        assert "copper_binding" in preds

        # Cobalt/B12 pattern
        preds = get_predicates_for_pfam("PF99997", "Unknown", "cobalamin binding protein")
        assert "cobalt_binding" in preds
        assert "cobalamin_binding" in preds


class TestTopologyModule:
    """Tests for topology prediction module."""

    def test_is_available(self):
        """Should check pyTMHMM availability."""
        from sharur.predicates.topology import is_available
        # Just check it returns a boolean (pyTMHMM may or may not be installed)
        result = is_available()
        assert isinstance(result, bool)

    def test_predict_topology_unavailable(self):
        """Should return None if pyTMHMM not available."""
        from sharur.predicates.topology import is_available, predict_topology
        if not is_available():
            result = predict_topology("MLALIFVLFFGLLASVLGG")
            assert result is None

    def test_topology_prediction_class(self):
        """Should have TopologyPrediction dataclass."""
        from sharur.predicates.topology import TopologyPrediction
        pred = TopologyPrediction(
            sequence_length=100,
            num_tm_helices=2,
            tm_helix_positions=[(10, 30), (50, 70)],
            topology="o" * 10 + "M" * 20 + "i" * 20 + "M" * 20 + "o" * 30,
            n_terminus_location="outside",
        )
        assert pred.is_transmembrane
        assert pred.is_multi_pass
        assert not pred.is_single_pass
        assert not pred.is_polytopic

    def test_get_topology_predicates(self):
        """Should generate predicates from topology prediction."""
        from sharur.predicates.topology import TopologyPrediction, get_topology_predicates

        # Multi-pass membrane protein
        pred = TopologyPrediction(
            sequence_length=300,
            num_tm_helices=4,
            tm_helix_positions=[(10, 30), (50, 70), (90, 110), (130, 150)],
            topology="o" * 300,
            n_terminus_location="outside",
        )
        preds = get_topology_predicates(pred)
        assert "transmembrane_predicted" in preds
        assert "multi_pass_membrane" in preds
        assert "polytopic_membrane" in preds
        assert "n_out_topology" in preds

        # Soluble protein
        pred = TopologyPrediction(
            sequence_length=200,
            num_tm_helices=0,
            tm_helix_positions=[],
            topology="i" * 200,
            n_terminus_location="inside",
        )
        preds = get_topology_predicates(pred)
        assert "soluble_predicted" in preds
        assert "transmembrane_predicted" not in preds

    def test_generator_without_sequence(self):
        """Generator should work without sequence (no topology predicates)."""
        gen = PredicateGenerator()
        protein = ProteinRecord(protein_id="test", sequence_length=300)
        preds = gen.generate_for_protein(protein, [])
        # Should have size predicates but not topology
        assert "medium" in preds
        assert "unannotated" in preds
        # No topology predicates without sequence
        assert "transmembrane_predicted" not in preds
        assert "soluble_predicted" not in preds
