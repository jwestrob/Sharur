"""Tests for V2 atom generator."""

import pytest

from sharur.predicates.generator import AnnotationRecord, ProteinRecord
from sharur.predicates_v2.generator import AtomGenerator
from sharur.predicates_v2.model import ClaimRelation, SemanticFacet
from sharur.predicates_v2.rules import clear_caches


@pytest.fixture(autouse=True)
def _reset_caches():
    """Clear rule caches before each test."""
    clear_caches()
    yield
    clear_caches()


class TestAtomGenerator:
    """Tests for AtomGenerator class."""

    def test_basic_pfam_annotation(self):
        """Should generate typed atoms from a PFAM annotation."""
        gen = AtomGenerator()
        protein = ProteinRecord(protein_id="test", sequence_length=500)
        annotations = [
            AnnotationRecord(
                source="pfam",
                accession="PF00005",
                name="ABC_tran",
                description="ABC transporter",
                evalue=1e-50,
                score=200.0,
            )
        ]
        atoms = gen.generate_atoms(protein, annotations)

        # Should have atoms for transporter, abc_transporter, etc.
        atom_ids = {a.atom_id for a in atoms}
        assert "transporter" in atom_ids
        assert "abc_transporter" in atom_ids
        assert "atp_binding" in atom_ids

        # Check that atoms have correct facets
        transporter_atoms = [a for a in atoms if a.atom_id == "transporter"]
        assert len(transporter_atoms) >= 1
        assert transporter_atoms[0].facet == SemanticFacet.role
        assert transporter_atoms[0].source_db == "pfam"
        assert transporter_atoms[0].evidence_evalue == 1e-50
        assert transporter_atoms[0].evidence_score == 200.0

    def test_kegg_annotation(self):
        """Should generate typed atoms from a KEGG annotation."""
        gen = AtomGenerator()
        protein = ProteinRecord(protein_id="test", sequence_length=500)
        annotations = [
            AnnotationRecord(
                source="kegg",
                accession="K00532",
                description="hydrogenase large subunit",
                evalue=1e-40,
                score=150.0,
            )
        ]
        atoms = gen.generate_atoms(protein, annotations)

        atom_ids = {a.atom_id for a in atoms}
        assert "hydrogenase" in atom_ids
        assert "hydrogen_metabolism" in atom_ids

        # KEGG default relation should be implies
        h_atoms = [a for a in atoms if a.atom_id == "hydrogenase"]
        assert len(h_atoms) >= 1
        assert h_atoms[0].relation == ClaimRelation.implies
        assert h_atoms[0].facet == SemanticFacet.activity

    def test_kofam_annotation_uses_kegg_mapping(self):
        """KOFAM K accessions should generate KEGG-derived atoms."""
        gen = AtomGenerator()
        protein = ProteinRecord(protein_id="test", sequence_length=500)
        annotations = [
            AnnotationRecord(
                source="kofam",
                accession="K00532",
                description="hydrogenase large subunit",
                evalue=1e-40,
                score=150.0,
            )
        ]
        atoms = gen.generate_atoms(protein, annotations)

        h_atoms = [a for a in atoms if a.atom_id == "hydrogenase"]
        assert h_atoms
        assert h_atoms[0].relation == ClaimRelation.implies
        assert h_atoms[0].source_db == "kofam"

    def test_size_atoms(self):
        """Should generate size atoms from protein properties."""
        gen = AtomGenerator()

        # Giant protein
        protein = ProteinRecord(protein_id="test", sequence_length=1500)
        atoms = gen.generate_atoms(protein, [])

        atom_ids = {a.atom_id for a in atoms}
        assert "giant" in atom_ids

        giant_atoms = [a for a in atoms if a.atom_id == "giant"]
        assert giant_atoms[0].facet == SemanticFacet.size_class
        assert giant_atoms[0].relation == ClaimRelation.implies
        assert giant_atoms[0].source_db == "_property"

    def test_massive_size(self):
        """Should generate both giant and massive for very large proteins."""
        gen = AtomGenerator()
        protein = ProteinRecord(protein_id="test", sequence_length=3000)
        atoms = gen.generate_atoms(protein, [])

        atom_ids = {a.atom_id for a in atoms}
        assert "giant" in atom_ids
        assert "massive" in atom_ids

    def test_unannotated_status(self):
        """Should mark unannotated proteins."""
        gen = AtomGenerator()
        protein = ProteinRecord(protein_id="test", sequence_length=300)
        atoms = gen.generate_atoms(protein, [])

        atom_ids = {a.atom_id for a in atoms}
        assert "unannotated" in atom_ids

        unannotated = [a for a in atoms if a.atom_id == "unannotated"]
        assert unannotated[0].facet == SemanticFacet.quality_flag
        assert unannotated[0].relation == ClaimRelation.implies

    def test_multi_domain_status(self):
        """Should detect multi-domain proteins."""
        gen = AtomGenerator()
        protein = ProteinRecord(protein_id="test", sequence_length=800)
        annotations = [
            AnnotationRecord(source="pfam", accession="PF00001", evalue=1e-20),
            AnnotationRecord(source="pfam", accession="PF00002", evalue=1e-20),
            AnnotationRecord(source="pfam", accession="PF00003", evalue=1e-20),
        ]
        atoms = gen.generate_atoms(protein, annotations)

        atom_ids = {a.atom_id for a in atoms}
        assert "multi_domain" in atom_ids
        assert "well_annotated" in atom_ids

    def test_unresolved_accession(self):
        """Should emit unresolved atom for unmapped accessions."""
        gen = AtomGenerator()
        protein = ProteinRecord(protein_id="test", sequence_length=300)
        # Use an accession that won't map to any predicate
        annotations = [
            AnnotationRecord(
                source="pfam",
                accession="PF_FAKE_UNMAPPED_12345",
                name="UnknownDomain",
                description="",
                evalue=1e-20,
            )
        ]
        atoms = gen.generate_atoms(protein, annotations)

        unresolved = [a for a in atoms if a.relation == ClaimRelation.unresolved]
        assert len(unresolved) >= 1
        assert unresolved[0].atom_id.startswith("unresolved:")
        assert unresolved[0].source_accession == "PF_FAKE_UNMAPPED_12345"

    def test_no_direct_access_atoms(self):
        """Should NOT include direct-access atoms (pfam:PF00005, etc.)."""
        gen = AtomGenerator()
        protein = ProteinRecord(protein_id="test", sequence_length=500)
        annotations = [
            AnnotationRecord(
                source="pfam", accession="PF00005", evalue=1e-20,
            )
        ]
        atoms = gen.generate_atoms(protein, annotations)

        # No direct access predicates should appear as atoms
        for atom in atoms:
            assert ":" not in atom.atom_id or atom.atom_id.startswith(("unresolved:", "_source_witness:"))

    def test_no_annotated_flag_atoms(self):
        """Should NOT include *_annotated atoms (handled as status)."""
        gen = AtomGenerator()
        protein = ProteinRecord(protein_id="test", sequence_length=500)
        annotations = [
            AnnotationRecord(
                source="pfam", accession="PF00005", evalue=1e-20,
            )
        ]
        atoms = gen.generate_atoms(protein, annotations)

        annotated_atoms = [a for a in atoms if a.atom_id.endswith("_annotated")]
        # pfam_annotated etc. should NOT appear as annotation-derived atoms
        # (they may appear as status atoms from _annotation_status_atoms)
        for a in annotated_atoms:
            assert a.source_db == "_computed"

    def test_evidence_passthrough(self):
        """Evidence values should be passed through from annotation."""
        gen = AtomGenerator()
        protein = ProteinRecord(protein_id="test", sequence_length=500)
        annotations = [
            AnnotationRecord(
                source="pfam",
                accession="PF00005",
                evalue=3.14e-42,
                score=187.3,
            )
        ]
        atoms = gen.generate_atoms(protein, annotations)

        # All annotation-derived atoms should carry the evidence values
        ann_atoms = [a for a in atoms if a.source_db == "pfam" and not a.atom_id.startswith("_source_witness:")]
        assert len(ann_atoms) > 0
        for atom in ann_atoms:
            assert atom.evidence_evalue == 3.14e-42
            assert atom.evidence_score == 187.3

    def test_cazy_annotation(self):
        """Should generate atoms from CAZy annotation."""
        gen = AtomGenerator()
        protein = ProteinRecord(protein_id="test", sequence_length=400)
        annotations = [
            AnnotationRecord(
                source="cazy", accession="GH5", evalue=1e-30,
            )
        ]
        atoms = gen.generate_atoms(protein, annotations)

        atom_ids = {a.atom_id for a in atoms}
        assert "carbohydrate_active" in atom_ids
        assert "glycoside_hydrolase" in atom_ids
        assert "cellulase" in atom_ids

        # CAZy default relation should be supports
        gh_atoms = [a for a in atoms if a.atom_id == "cellulase"]
        assert gh_atoms[0].relation == ClaimRelation.supports

    def test_hyddb_annotation(self):
        """Should generate atoms from HydDB annotation."""
        gen = AtomGenerator()
        protein = ProteinRecord(protein_id="test", sequence_length=600)
        annotations = [
            AnnotationRecord(
                source="hyddb", accession="NiFe", evalue=1e-50,
            )
        ]
        atoms = gen.generate_atoms(protein, annotations)

        atom_ids = {a.atom_id for a in atoms}
        assert "nife_hydrogenase" in atom_ids
        assert "hydrogenase" in atom_ids

        # HydDB default relation should be implies
        h_atoms = [a for a in atoms if a.atom_id == "nife_hydrogenase"]
        assert h_atoms[0].relation == ClaimRelation.implies

    def test_multi_source_annotations(self):
        """Should handle annotations from multiple sources."""
        gen = AtomGenerator()
        protein = ProteinRecord(protein_id="test", sequence_length=500)
        annotations = [
            AnnotationRecord(source="pfam", accession="PF00005", evalue=1e-20),
            AnnotationRecord(source="kegg", accession="K00001",
                             description="alcohol dehydrogenase [EC:1.1.1.1]",
                             evalue=1e-30),
        ]
        atoms = gen.generate_atoms(protein, annotations)

        # Should have atoms from both sources
        sources = {a.source_db for a in atoms if a.source_db not in ("_property", "_computed")}
        assert "pfam" in sources
        assert "kegg" in sources

    def test_hierarchy_expansion(self):
        """Should expand predicate hierarchy."""
        gen = AtomGenerator(expand_hierarchy=True)
        protein = ProteinRecord(protein_id="test", sequence_length=500)
        annotations = [
            AnnotationRecord(
                source="pfam", accession="PF00005", evalue=1e-20,
            )
        ]
        atoms = gen.generate_atoms(protein, annotations)

        atom_ids = {a.atom_id for a in atoms}
        # abc_transporter should bring in transporter via hierarchy
        assert "abc_transporter" in atom_ids
        assert "transporter" in atom_ids

    def test_gc_outlier_atom(self):
        """Should generate GC outlier atom when appropriate."""
        gen = AtomGenerator()
        protein = ProteinRecord(
            protein_id="test",
            sequence_length=300,
            gc_content=0.7,
            contig_gc_mean=0.4,
            contig_gc_std=0.05,
        )
        atoms = gen.generate_atoms(protein, [])

        atom_ids = {a.atom_id for a in atoms}
        assert "gc_outlier" in atom_ids
        assert "gc_high" in atom_ids

    def test_defensefinder_system_annotation(self):
        """Should generate atoms from system-validated DefenseFinder."""
        gen = AtomGenerator()
        protein = ProteinRecord(protein_id="test", sequence_length=400)
        annotations = [
            AnnotationRecord(
                source="defensefinder_system",
                accession="RM_type_II",
                name="RM/Type_II",
                evalue=1e-20,
            )
        ]
        atoms = gen.generate_atoms(protein, annotations)

        atom_ids = {a.atom_id for a in atoms}
        assert "defense_system" in atom_ids
        assert "restriction_modification" in atom_ids
        assert "rm_type_ii" in atom_ids

        # System-validated should be implies
        ds_atoms = [a for a in atoms if a.atom_id == "defense_system"]
        assert ds_atoms[0].relation == ClaimRelation.implies

    @pytest.mark.parametrize(
        ("system_name", "expected_atom"),
        [
            ("Eleos/Eleos", "defense_eleos"),
            ("HEC-09/HEC-09", "defense_hec"),
            ("Shedu/Shedu", "defense_shedu"),
            ("pAgo_LongA/pAgo_LongA", "defense_pago"),
            ("AbiE/AbiE", "abi_e"),
            ("PD-Lambda-5/PD-Lambda-5", "defense_pd_lambda"),
        ],
    )
    def test_defensefinder_system_subtype_atoms(self, system_name, expected_atom):
        """Validated DefenseFinder system names should emit subtype atoms."""
        gen = AtomGenerator()
        protein = ProteinRecord(protein_id="test", sequence_length=400)
        annotations = [
            AnnotationRecord(
                source="defensefinder_system",
                accession="sys_def_1",
                name=system_name,
                evalue=1e-20,
            )
        ]
        atoms = gen.generate_atoms(protein, annotations)

        atom_ids = {a.atom_id for a in atoms}
        assert "defense_system" in atom_ids
        assert expected_atom in atom_ids

    def test_txsscan_system_annotation_uses_vocabulary_ids(self):
        """System-validated TXSScan should emit lower-case V2 atom IDs."""
        gen = AtomGenerator()
        protein = ProteinRecord(protein_id="test", sequence_length=400)
        annotations = [
            AnnotationRecord(
                source="txsscan_system",
                accession="sys_sec_1",
                name="T5aSS/T5aSS",
            )
        ]
        atoms = gen.generate_atoms(protein, annotations)

        atom_ids = {a.atom_id for a in atoms}
        assert "secretion_system" in atom_ids
        assert "type_v_secretion" in atom_ids
        assert "type_V_secretion" not in atom_ids
