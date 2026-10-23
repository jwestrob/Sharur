"""DefenseFinder hit selection: MacSyFinder/DefenseFinder thresholds, FAM dedup, Stage 07 load."""

from __future__ import annotations

from importlib.machinery import SourceFileLoader
from pathlib import Path

import duckdb
import pytest

from sharur.defensefinder import selection as sel
from tests.defensefinder_fixtures import VERSION, hit, hits_frame, write_model_set


KB_PATH = Path(__file__).resolve().parents[1] / "src" / "ingest" / "07_build_knowledge_base.py"


@pytest.fixture
def model_set(tmp_path):
    return sel.load_model_set(write_model_set(tmp_path / "models"))


def _kept(model_set, *rows, **params):
    selected, stats, _ = sel.select_hits(hits_frame(*rows), model_set, sel.SelectionParams(**params))
    return list(zip(selected["sequence_id"], selected["hmm_name"], strict=True)), stats


class TestDefaults:
    def test_thresholds_follow_the_installed_macsyfinder_and_defensefinder_defaults(self):
        # macsylib 1.0.4 MacsyDefaults; DefenseFinder 2.0.1 overrides coverage to 0.4.
        assert sel.MACSYFINDER_I_EVALUE_SEL == 1e-3
        assert sel.MACSYFINDER_COVERAGE_PROFILE == 0.5
        assert sel.MACSYFINDER_E_VALUE_SEARCH == 0.1
        assert sel.DEFENSEFINDER_COVERAGE_PROFILE == 0.4
        params = sel.SelectionParams()
        assert (params.use_ga, params.i_evalue_sel, params.coverage_profile) == (True, 1e-3, 0.4)


class TestSelection:
    def test_short_protein_hit_at_moderate_evalue_with_full_coverage_is_kept(self, model_set):
        # E 1e-12 failed the former blanket 1e-15 cutoff; coverage 0.8 of a 200-position profile.
        kept, stats = _kept(model_set, hit("p1", "ToyB_internal", evalue=1e-12, hmm_from=21, hmm_to=180))
        assert kept == [("p1", "ToyB_internal")]
        assert stats.selected_rows == 1

    def test_low_profile_coverage_is_dropped(self, model_set):
        kept, stats = _kept(
            model_set,
            hit("p1", "ToyB_internal", evalue=1e-30, hmm_from=1, hmm_to=60),   # 0.30, no GA
            hit("p2", "Toy__ToyA", score=90, evalue=1e-30, hmm_from=1, hmm_to=30),  # 0.30, GA passed
        )
        assert kept == []
        assert stats.failed_coverage == 2

    def test_coverage_default_is_defensefinders_0_4(self, model_set):
        row = hit("p1", "ToyB_internal", evalue=1e-20, hmm_from=1, hmm_to=90)  # 90/200 = 0.45
        assert _kept(model_set, row)[0] == [("p1", "ToyB_internal")]
        assert _kept(model_set, row, coverage_profile=sel.MACSYFINDER_COVERAGE_PROFILE)[0] == []

    def test_ga_governs_profiles_that_carry_one(self, model_set):
        kept, stats = _kept(
            model_set,
            # below GA (25 bits) despite a strong E-value: dropped
            hit("p_below", "Toy__ToyA", score=24.0, evalue=1e-12),
            # above GA with full-sequence E 0.5: GA replaces the E search cutoff (0.1)
            hit("p_ga", "Toy__ToyA", score=26.0, evalue=0.5, i_evalue=1e-4),
            # the same E-value on a profile without GA fails the E search cutoff
            hit("p_noga", "ToyB_internal", score=26.0, evalue=0.5, i_evalue=1e-4,
                hmm_from=1, hmm_to=180),
        )
        assert kept == [("p_ga", "Toy__ToyA")]
        assert stats.failed_threshold == 2

    def test_ga_needs_both_sequence_and_domain_scores(self, model_set):
        kept, _ = _kept(model_set, hit("p1", "Toy__ToyA", score=40.0, dom_score=20.0))
        assert kept == []

    def test_i_evalue_selection_applies_on_top_of_ga(self, model_set):
        # MacSyFinder report.py filters every reported domain by i-Evalue, GA or not.
        kept, stats = _kept(model_set, hit("p1", "Toy__ToyA", score=30.0, evalue=1e-5, i_evalue=0.01))
        assert kept == []
        assert stats.failed_i_evalue == 1

    def test_profiles_resolve_by_internal_name_and_unknown_ones_are_counted(self, model_set):
        kept, stats = _kept(
            model_set,
            hit("p1", "ToyB_internal", hmm_from=1, hmm_to=180),
            hit("p2", "Not_In_Model_Set"),
        )
        assert kept == [("p1", "ToyB_internal")]
        assert stats.unknown_profile == 1
        assert stats.unknown_profile_names == ["Not_In_Model_Set"]

    def test_missing_alignment_coordinates_are_not_selected(self, model_set):
        row = hit("p1", "Toy__ToyA")
        row["hmm_from"] = row["hmm_to"] = None
        kept, stats = _kept(model_set, row)
        assert kept == []
        assert stats.missing_alignment == 1


    def test_antidefense_only_profiles_are_out_of_scope(self, model_set):
        # Anti__AntiX is referenced only by an ADF model; DefenseFinder searches it with -a only.
        kept, stats = _kept(model_set, hit("p1", "Anti__AntiX", score=90.0), hit("p2", "Solo__SoloA"))
        assert kept == [("p2", "Solo__SoloA")]
        assert stats.out_of_scope_profile == 1
        everything, _ = _kept(model_set, hit("p1", "Anti__AntiX", score=90.0),
                              subfamilies=("DefenseFinder", "RM", "ADF"))
        assert everything == [("p1", "Anti__AntiX")]

    def test_profile_scope_matches_the_callers_default_subfamilies(self):
        from sharur.colocation import DEFAULT_MODEL_SUBFAMILIES

        assert sel.DEFAULT_SUBFAMILIES == DEFAULT_MODEL_SUBFAMILIES["defensefinder"]


class TestAnnotationRows:
    def test_fam_variant_dedup_still_applies(self, model_set):
        rows, stats = sel.build_annotations(
            hits_frame(
                hit("p1", "RM_Type_X__MTases_FAM_1", score=50.0, hmm_from=1, hmm_to=140),
                hit("p1", "RM_Type_X__MTases_FAM_2", score=60.0, hmm_from=1, hmm_to=140),
                hit("p2", "RM_Type_X__MTases_FAM_1", score=55.0, hmm_from=1, hmm_to=140),
                hit("p1", "Solo__SoloA", score=70.0),
            ),
            model_set,
        )
        assert stats.selected_rows == 4
        assert stats.fam_dedup_removed == 1
        assert sorted(zip(rows.protein_id, rows.accession, strict=True)) == [
            ("p1", "RM_Type_X__MTases_FAM_2"), ("p1", "Solo__SoloA"), ("p2", "RM_Type_X__MTases_FAM_1"),
        ]

    def test_rows_follow_stage07_column_conventions(self, model_set):
        rows, _ = sel.build_annotations(
            hits_frame(hit("p1", "Solo__SoloA", score=70.0, dom_score=65.0, evalue=1e-25,
                           i_evalue=1e-24, ali=(11, 99))),
            model_set,
        )
        row = rows.iloc[0]
        assert (row.source, row.accession, row["name"], row.description) == (
            "defensefinder", "Solo__SoloA", "Solo__SoloA", "")
        # full-sequence E-value, domain score (MacSyFinder's hit score), alignment span
        assert (row.evalue, row.score, row.start_aa, row.end_aa) == (1e-25, 65.0, 11, 99)

    def test_model_version_and_load_are_recorded(self, model_set, tmp_path):
        assert model_set.version == VERSION
        assert model_set.models_dir == (tmp_path / "models").resolve()
        _rows, stats = sel.build_annotations(hits_frame(hit("p1", "Solo__SoloA")), model_set)
        conn = duckdb.connect()
        sel.record_source_version(conn, model_set, sel.SelectionParams(), stats, "hits.tsv", "abc")
        record = conn.execute(
            "SELECT source, model_set, model_version, profile_count, raw_rows, loaded_rows, "
            "parameters FROM annotation_source_versions"
        ).fetchone()
        assert record[:6] == ("defensefinder", "defense-finder-models", VERSION, 6, 1, 1)
        assert '"coverage_profile": 0.4' in record[6]
        assert '"subfamilies": ["DefenseFinder", "RM"]' in record[6]

    def test_model_set_resolution_prefers_explicit_then_environment(self, tmp_path, monkeypatch):
        family = write_model_set(tmp_path / "env_models", version="1.2.3")
        monkeypatch.setenv(sel.MODELS_ENV, str(family))
        assert sel.find_model_set().version == "1.2.3"
        explicit = write_model_set(tmp_path / "explicit", version="4.5.6")
        assert sel.find_model_set(explicit / "profiles").version == "4.5.6"


def _load_kb_module():
    return SourceFileLoader("kb_build_df", str(KB_PATH)).load_module()


def test_stage07_loads_defensefinder_hits_with_macsyfinder_selection(tmp_path):
    module = _load_kb_module()
    builder_cls, outputs_cls = module.KnowledgeBaseBuilder, module.PipelineOutputs
    # The blanket e-value cutoff no longer applies to DefenseFinder.
    assert "defensefinder" not in module.DEFAULT_EVALUE_THRESHOLDS
    data = tmp_path / "data"
    outputs = outputs_cls(**{
        name: data / name for name in (
            "stage00_dir", "stage01_dir", "stage02_dir", "stage03_dir", "stage04_dir",
            "stage05a_dir", "stage05b_dir", "stage05c_dir", "stage06_dir",
        )
    })
    hits_dir = outputs.stage04_dir / "defensefinder_results"
    hits_dir.mkdir(parents=True)
    hits_frame(
        hit("p1", "Solo__SoloA", evalue=1e-12),                         # kept (former cutoff dropped it)
        hit("p2", "Toy__ToyA", score=20.0),                             # below GA
        hit("p3", "ToyB_internal", hmm_from=1, hmm_to=40),              # coverage 0.2
        hit("p4", "RM_Type_X__MTases_FAM_1", score=50.0, hmm_from=1, hmm_to=140),
        hit("p4", "RM_Type_X__MTases_FAM_2", score=45.0, hmm_from=1, hmm_to=140),
    ).to_csv(hits_dir / "DefenseFinder_hits_df.tsv", sep="\t", index=False)
    family = write_model_set(tmp_path / "models")

    builder = builder_cls(outputs, tmp_path / "sharur.duckdb", force=True, threads=1,
                          defensefinder_models=family)
    builder._init_db()
    try:
        builder.conn.execute("INSERT INTO bins (bin_id) VALUES ('b1')")
        builder.conn.execute("INSERT INTO contigs (contig_id, bin_id, length) VALUES ('c1', 'b1', 9000)")
        for i in range(1, 5):
            builder.conn.execute(
                "INSERT INTO proteins (protein_id, contig_id, bin_id, start, end_coord, strand, "
                "gene_index) VALUES (?, 'c1', 'b1', ?, ?, '+', ?)",
                [f"p{i}", i * 1000, i * 1000 + 600, i],
            )
        builder._load_annotations()
        rows = builder.conn.execute(
            "SELECT protein_id, accession FROM annotations WHERE source = 'defensefinder' ORDER BY 1"
        ).fetchall()
        record = builder.conn.execute(
            "SELECT model_version, raw_rows, selected_rows, loaded_rows FROM annotation_source_versions"
        ).fetchone()
    finally:
        builder._release_db()
    assert rows == [("p1", "Solo__SoloA"), ("p4", "RM_Type_X__MTases_FAM_1")]
    assert record == (VERSION, 5, 3, 2)
