"""Optional assembly/host-evidence sidecar behavior."""

from __future__ import annotations

import gzip

import pytest

from sharur.assembly_evidence import (
    AssemblyEvidenceStore,
    _tetranucleotide_counts,
    compute_composition_evidence,
    import_contig_evidence,
)
from sharur.capabilities import CapabilityState, build_capability_brief
from sharur.operators import Sharur


def test_tetranucleotide_signature_is_strand_invariant():
    sequence = "AACCGTTAAGCTCCGAT"
    reverse_complement = sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]

    assert _tetranucleotide_counts(sequence) == _tetranucleotide_counts(reverse_complement)


def test_tabular_import_is_discovered_by_case_inspection(case_database, tmp_path):
    input_path = tmp_path / "coverage.tsv"
    input_path.write_text(
        "bin_id\tcontig_id\tcoverage_mean\tmapped_reads\tcaller_note\n"
        "fg_plus\tfg_plus_contig\t12.5\t400\twell_supported\n",
        encoding="utf-8",
    )
    sidecar = case_database.parent / "assembly_evidence.duckdb"

    result = import_contig_evidence(
        input_path,
        sidecar,
        source="coverage_pipeline_v1",
        validate_db_path=case_database,
    )
    assert result["rows_written"] == 1
    assert result["validation"]["state"] == "validated"

    case = Sharur(case_database, read_only=True).inspect(
        "target_plus",
        entity_type="system",
    )
    evidence = case.record.assembly_evidence
    assert case.record.assembly_evidence_state == "available"
    assert evidence is not None
    assert evidence.coverage_mean == pytest.approx(12.5)
    assert evidence.coverage_ratio_to_bin_median == pytest.approx(1.0)
    assert evidence.mapped_reads == 400
    assert evidence.metadata["caller_note"] == "well_supported"
    assert evidence.metadata["import_source"] == "coverage_pipeline_v1"


def test_composition_is_opt_in_scalar_only_and_preserves_existing_fields(
    case_database,
    tmp_path,
):
    input_path = tmp_path / "coverage.tsv"
    input_path.write_text(
        "bin_id\tcontig_id\tcoverage_mean\tproper_pair_fraction\n"
        "fg_plus\tfg_plus_contig\t9.0\t0.92\n",
        encoding="utf-8",
    )
    sidecar = case_database.parent / "assembly_evidence.duckdb"
    import_contig_evidence(
        input_path,
        sidecar,
        source="mapping",
        validate_db_path=case_database,
        hash_input=False,
    )

    fasta = tmp_path / "fg_plus.fna.gz"
    with gzip.open(fasta, mode="wt", encoding="utf-8") as handle:
        handle.write(
            ">fg_plus_contig\n" + "AAAACCCCGGGGTTTT" * 80 + "\n"
            ">unbinned_comparator\n" + "ATATATATATATATAT" * 80 + "\n"
        )
    result = compute_composition_evidence(
        {"fg_plus": fasta},
        sidecar,
        validate_db_path=None,
    )
    assert result["assemblies"][0]["vectors_persisted"] is False

    with AssemblyEvidenceStore(sidecar, read_only=True) as store:
        evidence = store.get("fg_plus", "fg_plus_contig")
        columns = {row[0] for row in store.conn.execute("DESCRIBE contig_evidence").fetchall()}
    assert evidence is not None
    assert evidence.coverage_mean == pytest.approx(9.0)
    assert evidence.proper_pair_fraction == pytest.approx(0.92)
    assert evidence.gc_zscore is not None
    assert evidence.tetranucleotide_distance is not None
    assert evidence.tetranucleotide_percentile is not None
    assert evidence.metadata["import_source"] == "mapping"
    assert evidence.metadata["vectors_persisted"] is False
    assert len(evidence.metadata["composition_content_sha256"]) == 64
    assert len(result["assemblies"][0]["content_sha256"]) == 64
    assert not any("vector" in column for column in columns)


def test_preflight_reports_optional_sidecar_without_affecting_required_state(
    case_database,
):
    absent = build_capability_brief(
        case_database,
        include_execution=False,
    ).get("assembly_evidence")
    assert absent.state == CapabilityState.unavailable
    assert absent.required is False
    assert absent.evidence["composition_computation"] == "not_run"

    sidecar = case_database.parent / "assembly_evidence.duckdb"
    with AssemblyEvidenceStore(sidecar) as store:
        store.upsert(
            [
                {
                    "bin_id": "fg_plus",
                    "contig_id": "fg_plus_contig",
                    "coverage_mean": 4.0,
                }
            ]
        )
    available = build_capability_brief(
        case_database,
        include_execution=False,
    ).get("assembly_evidence")
    assert available.state == CapabilityState.available
    assert available.evidence["rows"] == 1
    assert available.evidence["coverage_rows"] == 1
    assert available.evidence["vectors_persisted_by_sharur"] is False


def test_import_validation_fails_closed_on_unknown_replicon(
    case_database,
    tmp_path,
):
    input_path = tmp_path / "bad.tsv"
    input_path.write_text(
        "bin_id\tcontig_id\tcoverage_mean\nmissing\tmissing_contig\t5\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="do not match"):
        import_contig_evidence(
            input_path,
            tmp_path / "evidence.duckdb",
            validate_db_path=case_database,
        )

    assert not (tmp_path / "evidence.duckdb").exists()


def test_composition_leaves_undefined_metrics_null_and_never_uses_core_as_sidecar(
    case_database,
    tmp_path,
):
    fasta = tmp_path / "ambiguous.fna"
    fasta.write_text(
        ">only_ns\nNNNNNNNN\n>canonical\nAAAACCCC\n",
        encoding="utf-8",
    )
    sidecar = tmp_path / "assembly_evidence.duckdb"
    compute_composition_evidence(
        {"test_bin": fasta},
        sidecar,
        validate_db_path=None,
    )
    with AssemblyEvidenceStore(sidecar, read_only=True) as store:
        ambiguous = store.get("test_bin", "only_ns")

    assert ambiguous is not None
    assert ambiguous.gc_zscore is None
    assert ambiguous.tetranucleotide_distance is None

    with pytest.raises(ValueError, match="separate sidecar"):
        compute_composition_evidence(
            {"test_bin": fasta},
            case_database,
            validate_db_path=case_database,
        )
