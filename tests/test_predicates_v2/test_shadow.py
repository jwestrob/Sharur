"""Tests for V2 shadow diff."""

import pytest

from sharur.predicates.generator import AnnotationRecord, ProteinRecord
from sharur.predicates_v2.shadow import (
    format_shadow_diff_jsonl,
    format_shadow_diff_report,
    shadow_diff,
)
from sharur.predicates_v2.rules import clear_caches
from sharur.predicates_v2.composites import clear_composites_cache


@pytest.fixture(autouse=True)
def _reset_caches():
    """Clear caches before each test."""
    clear_caches()
    clear_composites_cache()
    yield
    clear_caches()
    clear_composites_cache()


class TestShadowDiff:
    """Tests for shadow_diff function."""

    def test_empty_input(self):
        """Should handle empty protein list."""
        result = shadow_diff([], {})
        assert result["summary"]["total_proteins"] == 0
        assert result["per_protein"] == []

    def test_single_unannotated_protein(self):
        """Should process a single unannotated protein."""
        proteins = [ProteinRecord(protein_id="test", sequence_length=300)]
        result = shadow_diff(proteins, {})

        assert result["summary"]["total_proteins"] == 1
        # Both V1 and V2 should produce similar predicates for unannotated
        # (medium, unannotated)

    def test_annotated_protein(self):
        """Should process a protein with annotations."""
        proteins = [ProteinRecord(protein_id="test", sequence_length=500)]
        annotations = {
            "test": [
                AnnotationRecord(
                    source="pfam",
                    accession="PF00005",
                    name="ABC_tran",
                    description="ABC transporter",
                    evalue=1e-50,
                    score=200.0,
                )
            ]
        }
        result = shadow_diff(proteins, annotations)

        assert result["summary"]["total_proteins"] == 1
        # Should have some predicate counts
        assert result["summary"]["v1_unique_predicates"] > 0
        assert result["summary"]["v2_unique_predicates"] > 0

    def test_diff_structure(self):
        """Per-protein diff should have expected keys."""
        proteins = [ProteinRecord(protein_id="test", sequence_length=500)]
        annotations = {
            "test": [
                AnnotationRecord(
                    source="pfam",
                    accession="PF00005",
                    evalue=1e-50,
                )
            ]
        }
        result = shadow_diff(proteins, annotations)

        # If there are diffs, check structure
        for diff_entry in result["per_protein"]:
            assert "protein_id" in diff_entry
            assert "v1_only" in diff_entry
            assert "v2_only" in diff_entry
            assert "common_count" in diff_entry
            assert "v2_composites" in diff_entry
            assert "v2_unresolved_count" in diff_entry

    def test_summary_structure(self):
        """Summary should have expected keys."""
        proteins = [ProteinRecord(protein_id="test", sequence_length=300)]
        result = shadow_diff(proteins, {})

        summary = result["summary"]
        assert "total_proteins" in summary
        assert "match_count" in summary
        assert "diff_count" in summary
        assert "match_rate_pct" in summary
        assert "v1_unique_predicates" in summary
        assert "v2_unique_predicates" in summary
        assert "v2_total_atoms" in summary
        assert "v2_unresolved_atoms" in summary
        assert "v2_unique_unresolved" in summary
        assert "top_20_atoms" in summary

    def test_multiple_proteins(self):
        """Should handle multiple proteins."""
        proteins = [
            ProteinRecord(protein_id="prot_001", sequence_length=300),
            ProteinRecord(protein_id="prot_002", sequence_length=1500),
            ProteinRecord(protein_id="prot_003", sequence_length=50),
        ]
        result = shadow_diff(proteins, {})

        assert result["summary"]["total_proteins"] == 3
        # All should have size predicates
        assert result["summary"]["match_count"] + result["summary"]["diff_count"] == 3

    def test_unresolved_accessions_counted(self):
        """Unresolved accessions should be counted in summary."""
        proteins = [ProteinRecord(protein_id="test", sequence_length=300)]
        annotations = {
            "test": [
                AnnotationRecord(
                    source="pfam",
                    accession="PF_FAKE_UNMAPPED_99999",
                    name="NoMatch",
                    description="",
                    evalue=1e-20,
                )
            ]
        }
        result = shadow_diff(proteins, annotations)
        assert result["summary"]["v2_unresolved_atoms"] >= 1

    def test_top_atoms_present(self):
        """Should have top atom frequency data."""
        proteins = [
            ProteinRecord(protein_id=f"prot_{i}", sequence_length=300)
            for i in range(5)
        ]
        result = shadow_diff(proteins, {})

        # Should have at least some atoms in top list
        assert len(result["summary"]["top_20_atoms"]) > 0

    def test_top_atom_structure(self):
        """Top atom entries should have expected fields."""
        proteins = [ProteinRecord(protein_id="test", sequence_length=300)]
        result = shadow_diff(proteins, {})

        for atom_info in result["summary"]["top_20_atoms"]:
            assert "atom_id" in atom_info
            assert "total_count" in atom_info
            assert "facets" in atom_info
            assert "relations" in atom_info


class TestFormatShadowDiffReport:
    """Tests for report formatting."""

    def test_report_is_markdown(self):
        """Report should be valid markdown."""
        proteins = [ProteinRecord(protein_id="test", sequence_length=300)]
        diff = shadow_diff(proteins, {})
        report = format_shadow_diff_report(diff)

        assert report.startswith("# Predicate V2 Shadow Diff Report")
        assert "## Summary" in report
        assert "## Top 20 Atoms" in report

    def test_report_contains_stats(self):
        """Report should contain summary statistics."""
        proteins = [ProteinRecord(protein_id="test", sequence_length=300)]
        diff = shadow_diff(proteins, {})
        report = format_shadow_diff_report(diff)

        assert "Total proteins: 1" in report
        assert "Exact match:" in report


class TestFormatShadowDiffJsonl:
    """Tests for JSONL formatting."""

    def test_empty_diff(self):
        """Should handle no diffs."""
        diff = {"per_protein": [], "summary": {}}
        jsonl = format_shadow_diff_jsonl(diff)
        assert jsonl == ""

    def test_jsonl_format(self):
        """Each line should be valid JSON."""
        diff = {
            "per_protein": [
                {
                    "protein_id": "test",
                    "v1_only": ["foo"],
                    "v2_only": ["bar"],
                    "common_count": 5,
                    "v2_composites": [],
                    "v2_unresolved_count": 0,
                }
            ],
            "summary": {},
        }
        jsonl = format_shadow_diff_jsonl(diff)
        lines = jsonl.strip().split("\n")
        assert len(lines) == 1

        import json
        parsed = json.loads(lines[0])
        assert parsed["protein_id"] == "test"
        assert parsed["v1_only"] == ["foo"]
