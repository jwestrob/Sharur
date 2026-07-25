"""First-class case inspection and context-cohort comparisons."""

from __future__ import annotations

import pytest

from sharur.operators import Sharur
from sharur.storage.duckdb_store import DuckDBStore


def test_inspect_separates_structured_calls_from_raw_profiles(case_database):
    case = Sharur(case_database, read_only=True).inspect(
        "target_plus",
        entity_type="system",
        upstream_orfs=2,
        downstream_orfs=2,
    )

    assert case.record.entity.source_table == "defense_systems"
    assert case.record.entity.name == "Target"
    assert case.record.components == ["fg_plus_p3", "fg_plus_p4"]
    assert {call.call_id for call in case.record.named_calls} == {"target_plus"}

    component = next(
        protein for protein in case.record.proteins if protein.protein_id == "fg_plus_p3"
    )
    assert {annotation.source for annotation in component.annotations} == {"defensefinder"}
    assert component.named_calls[0].call_type == "Target"
    assert "NAMED `Target`" in case.to_markdown()
    assert "OBSERVED `defensefinder:raw_profile`" in case.to_markdown()


def test_asymmetric_context_uses_normalized_biological_strand(case_database):
    case = Sharur(case_database, read_only=True).inspect(
        "target_minus",
        entity_type="system",
        window=9,
        upstream_orfs=2,
        downstream_orfs=1,
    )

    context = case.record.context_window
    assert context is not None
    assert context.default_orfs == 9
    assert context.upstream_orfs == 2
    assert context.downstream_orfs == 1
    assert context.orientation == "-"
    assert context.requested_min_gene_index == 2
    assert context.requested_max_gene_index == 6
    assert context.realized_upstream_orfs == 2
    assert context.realized_downstream_orfs == 1

    by_id = {protein.protein_id: protein for protein in case.record.proteins}
    assert by_id["fg_minus_p6"].relative_orf == -2
    assert by_id["fg_minus_p6"].region_role == "upstream"
    assert by_id["fg_minus_p2"].relative_orf == 1
    assert by_id["fg_minus_p2"].region_role == "downstream"
    assert by_id["fg_minus_p3"].strand == "-"


def test_compare_context_reports_exact_counts_and_deduplicates_replicons(
    case_database,
):
    case = Sharur(case_database, read_only=True).inspect(
        "target_plus",
        entity_type="system",
    )

    by_entity = case.compare_context(
        features=["pfam:PFTEST"],
        window=2,
        min_components=2,
        deduplicate_by="entity",
    )
    assert (
        by_entity.composite.foreground_positive,
        by_entity.composite.foreground_total,
    ) == (2, 2)
    assert (
        by_entity.composite.background_positive,
        by_entity.composite.background_total,
    ) == (2, 3)
    assert by_entity.composite.fisher_p == pytest.approx(0.6)
    assert any("1 background" in limitation for limitation in by_entity.limitations)

    by_replicon = case.compare_context(
        features=["pfam:PFTEST"],
        window=2,
        min_components=2,
        deduplicate_by="replicon",
    )
    assert (
        by_replicon.composite.background_positive,
        by_replicon.composite.background_total,
    ) == (1, 2)
    assert by_replicon.composite.fisher_p == pytest.approx(0.5)
    assert "2/2" in by_replicon.to_markdown()
    assert by_replicon.recipe["upstream_orfs"] == 2
    assert by_replicon.recipe["downstream_orfs"] == 2


def test_feature_distance_and_asymmetric_comparison_windows_are_independent(
    case_database,
):
    case = Sharur(case_database, read_only=True).inspect(
        "target_plus",
        entity_type="system",
    )
    comparison = case.compare_context(
        features=[
            {
                "kind": "annotation",
                "source": "pfam",
                "value": "PFTEST",
                "max_orfs": 1,
            }
        ],
        window=7,
        upstream_orfs=2,
        downstream_orfs=1,
        min_components=2,
        deduplicate_by="entity",
    )

    assert comparison.context_window.default_orfs == 7
    assert comparison.context_window.upstream_orfs == 2
    assert comparison.context_window.downstream_orfs == 1
    assert comparison.composite.foreground_positive == 0
    assert comparison.composite.foreground_total == 2


def test_projection_views_do_not_make_locus_resolution_ambiguous(case_database):
    store = DuckDBStore(case_database)
    store.conn.execute(
        """
        INSERT INTO loci(
            locus_id, locus_type, contig_id, start, end_coord, confidence
        ) VALUES ('test_bgc', 'bgc', 'fg_plus_contig', 901, 1470, 0.9)
        """
    )
    store.conn.executemany(
        """
        INSERT INTO locus_proteins(locus_id, protein_id, position)
        VALUES ('test_bgc', ?, ?)
        """,
        [("fg_plus_p3", 1), ("fg_plus_p4", 2)],
    )
    store.close()

    case = Sharur(case_database, read_only=True).inspect(
        "test_bgc",
        entity_type="locus",
        window=1,
    )

    assert case.record.entity.source_table == "loci"
    assert case.record.entity.name == "bgc"
    assert case.record.components == ["fg_plus_p3", "fg_plus_p4"]


def test_compare_context_rejects_ambiguous_or_invalid_python_arguments(
    case_database,
):
    case = Sharur(case_database, read_only=True).inspect(
        "target_plus",
        entity_type="system",
    )

    with pytest.raises(ValueError, match="feature IDs must be unique"):
        case.compare_context(
            features=[
                {
                    "feature_id": "duplicate",
                    "kind": "annotation",
                    "source": "pfam",
                    "value": "PFTEST",
                },
                {
                    "feature_id": "duplicate",
                    "kind": "annotation",
                    "source": "pfam",
                    "value": "PFTEST",
                    "max_orfs": 1,
                },
            ],
            window=2,
        )
    with pytest.raises(ValueError, match="combine"):
        case.compare_context(  # type: ignore[arg-type]
            features=["pfam:PFTEST"],
            combine="neither",
        )
    with pytest.raises(ValueError, match="overlap"):
        case.compare_context(
            features=["pfam:PFTEST"],
            foreground_ids=["target_plus"],
            background_ids=["target_plus"],
        )


class TestCallerProfileEvidenceLevel:
    """Raw caller profile hits must be distinguishable from observed domains.

    Only the assembled projection (`defensefinder_system`) sits behind the
    provenance boundary; the caller's raw per-profile hits stay in
    observed_annotations with the same shape as a PFAM domain. In the
    Dormibacteria pilot 32,853 proteins carried a DefenseFinder profile hit
    while only 1,879 sat inside an assembled system -- 94.3% were fragments the
    caller explicitly declined to call. Unmarked, a reading agent re-reports
    them as discoveries.
    """

    @staticmethod
    def _mark(observed, named_calls, projection_sources=("defensefinder_system",)):
        from sharur.operators.cases import _mark_caller_profile_evidence

        _mark_caller_profile_evidence(observed, named_calls, list(projection_sources))
        return observed

    def test_unassembled_profile_hit_is_marked(self):
        observed = {"p1": [{"source": "defensefinder", "accession": "Gabija__GajB_2"}]}
        out = self._mark(observed, {})
        assert out["p1"][0]["evidence_level"] == "caller_profile_unassembled"

    def test_assembled_profile_hit_is_marked(self):
        observed = {"p1": [{"source": "defensefinder", "accession": "Wadjet__JetA_II"}]}
        named = {"p1": [{"system_source": "defensefinder_system", "call_type": "Wadjet_II"}]}
        out = self._mark(observed, named)
        assert out["p1"][0]["evidence_level"] == "caller_profile_assembled"

    def test_pfam_domains_are_left_alone(self):
        observed = {"p1": [{"source": "pfam", "accession": "PF00485"}]}
        out = self._mark(observed, {})
        assert "evidence_level" not in out["p1"][0]

    def test_no_projection_sources_is_a_noop(self):
        observed = {"p1": [{"source": "defensefinder", "accession": "Gabija__GajB_2"}]}
        out = self._mark(observed, {}, projection_sources=())
        assert "evidence_level" not in out["p1"][0]

    def test_mixed_protein_marks_only_the_caller_hit(self):
        observed = {
            "p1": [
                {"source": "pfam", "accession": "PF00004"},
                {"source": "defensefinder", "accession": "Pycsar__AG_cyclase"},
            ]
        }
        out = self._mark(observed, {})
        assert "evidence_level" not in out["p1"][0]
        assert out["p1"][1]["evidence_level"] == "caller_profile_unassembled"
