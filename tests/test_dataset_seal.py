"""Dataset seal generation, verification, and CLI contracts."""

from __future__ import annotations

import hashlib
import json
import shutil

import pytest
from typer.testing import CliRunner

from sharur.cli import app
from sharur.dataset_seal import (
    DatasetSealError,
    build_dataset_seal,
    verify_dataset_seal,
    write_dataset_seal,
)
from sharur.ops.ledger import RunLedger
from sharur.storage.duckdb_store import DuckDBStore


def _build_dataset(root):
    db_path = root / "sharur.duckdb"
    store = DuckDBStore(db_path)
    store.execute(
        "INSERT INTO bins(bin_id, n_contigs, total_length) VALUES ('bin1', 1, 300)"
    )
    store.execute(
        "INSERT INTO contigs(contig_id, bin_id, length) VALUES ('contig1', 'bin1', 300)"
    )
    store.execute(
        """
        INSERT INTO proteins(
            protein_id, contig_id, bin_id, start, end_coord, strand, sequence,
            sequence_length
        ) VALUES ('p1', 'contig1', 'bin1', 1, 300, '+', 'MA', 2)
        """
    )
    store.execute(
        """
        INSERT INTO annotations(
            annotation_id, protein_id, source, accession, name
        ) VALUES (1, 'p1', 'observed_domains', 'T0001', 'Observed test domain')
        """
    )
    store.execute(
        """
        INSERT INTO semantic_state(
            protein_id, activities, roles, architecture, localization, topology,
            size_class, quality_flags, composite_predicates, unresolved_count
        ) VALUES ('p1', [], [], [], [], '{}', 'standard', [], [], 0)
        """
    )
    store.execute(
        """
        INSERT INTO semantic_terms(
            protein_id, term_id, term_kind, facet, relation, source_db,
            source_accession
        ) VALUES ('p1', 'standard', 'state', 'size', 'observed', '', '')
        """
    )
    store.execute(
        "INSERT INTO protein_predicates(protein_id, predicates) VALUES ('p1', [])"
    )
    store.execute(
        """
        CREATE TABLE custom_caller_output (
            system_id VARCHAR,
            system_type VARCHAR
        )
        """
    )
    store.execute(
        "INSERT INTO custom_caller_output VALUES ('call-1', 'caller-emitted-name')"
    )
    store.close()

    stage00 = root / "stage00_prepared"
    stage00.mkdir()
    assembly = stage00 / "assembly.fna"
    assembly.write_text(">contig1\nACGTACGT\n")
    assembly_sha = hashlib.sha256(assembly.read_bytes()).hexdigest()
    report = {
        "schema_version": 1,
        "status": "passed",
        "summary": {
            "files_scanned": 1,
            "records_scanned": 1,
            "total_bases": 8,
            "errors": 0,
            "warnings": 0,
        },
        "files": [
            {
                "filename": assembly.name,
                "checksums": {"sha256": assembly_sha},
            }
        ],
    }
    (stage00 / "input_integrity.json").write_text(json.dumps(report))
    manifest = {
        "version": "0.2.0",
        "stage": "stage00_prepare_inputs",
        "integrity_status": "passed",
        "source_set_sha256": "source-set-test",
        "summary": report["summary"],
        "genomes": [
            {
                "filename": assembly.name,
                "genome_id": "assembly",
                "file_size": assembly.stat().st_size,
                "sha256": assembly_sha,
                "format_valid": True,
                "sequence_count": 1,
                "total_length": 8,
                "output_path": str(assembly),
            }
        ],
    }
    (stage00 / "processing_manifest.json").write_text(json.dumps(manifest))

    survey = root / "survey"
    survey.mkdir()
    (survey / "findings.jsonl").write_text(
        '{"id":"F1","summary":"verified fixture","verification":[]}\n'
    )
    return db_path, assembly_sha


def test_seal_records_inputs_live_callers_counts_and_canonical_artifacts(tmp_path):
    db_path, assembly_sha = _build_dataset(tmp_path)

    seal = build_dataset_seal(db_path, max_hash_bytes=0)

    assert seal["schema_version"] == 1
    assert seal["seal_strength"] == "structural"
    assert len(seal["dataset_id"]) == 64
    assert seal["identity"]["inputs"]["assemblies"][0]["sha256"] == assembly_sha
    counts = {
        item["table"]: item["rows"]
        for item in seal["identity"]["database"]["tables"]
    }
    assert counts["proteins"] == 1
    callers = seal["identity"]["annotation_contract"]["structured_callers"]
    custom = next(item for item in callers if item["table"] == "custom_caller_output")
    assert custom["emitted_types"] == [
        {"value": "caller-emitted-name", "rows": 1}
    ]
    artifacts = {
        item["path"]: item for item in seal["identity"]["canonical_artifacts"]
    }
    assert "survey/findings.jsonl" in artifacts
    assert artifacts["sharur.duckdb"]["digest"]["scope"] == "sampled"
    capability_ids = {
        item["capability_id"]
        for item in seal["provenance"]["capability_brief"]["capabilities"]
    }
    assert not any(item.startswith("execution_") for item in capability_ids)


def test_seal_verification_detects_database_content_drift(tmp_path):
    db_path, _assembly_sha = _build_dataset(tmp_path)
    seal_path = tmp_path / "dataset.seal.json"
    write_dataset_seal(build_dataset_seal(db_path), seal_path)

    unchanged = verify_dataset_seal(seal_path)
    assert unchanged.valid
    assert unchanged.changed_sections == ()

    store = DuckDBStore(db_path)
    store.execute(
        "UPDATE annotations SET name = 'Changed observation' WHERE annotation_id = 1"
    )
    store.close()

    changed = verify_dataset_seal(seal_path)
    assert not changed.valid
    assert "canonical_artifacts" in changed.changed_sections
    assert "sharur.duckdb" in changed.changed_artifacts


def test_seal_verifies_after_dataset_directory_is_copied(tmp_path):
    original = tmp_path / "original"
    original.mkdir()
    db_path, _assembly_sha = _build_dataset(original)
    seal_path = original / "dataset.seal.json"
    write_dataset_seal(build_dataset_seal(db_path), seal_path)

    copied = tmp_path / "copied"
    shutil.copytree(original, copied)

    verification = verify_dataset_seal(copied / "dataset.seal.json")
    assert verification.valid
    assert verification.database_path == str(copied / "sharur.duckdb")


def test_full_seal_hashes_every_discovered_canonical_artifact(tmp_path):
    db_path, _assembly_sha = _build_dataset(tmp_path)

    seal = build_dataset_seal(db_path, hash_large_files=True, max_hash_bytes=0)

    assert seal["seal_strength"] == "content"
    assert {
        artifact["digest"]["scope"]
        for artifact in seal["identity"]["canonical_artifacts"]
    } == {"full"}


def test_operational_ledger_drift_does_not_change_scientific_identity(tmp_path):
    db_path, _assembly_sha = _build_dataset(tmp_path)
    seal_path = tmp_path / "dataset.seal.json"
    write_dataset_seal(build_dataset_seal(db_path), seal_path)

    ledger = RunLedger(tmp_path / "sharur_ops.db")
    ledger.create_run("analysis", tmp_path, created_by="test")
    ledger.close()

    verification = verify_dataset_seal(seal_path)
    assert verification.valid
    assert verification.provenance_changed


def test_seal_refuses_active_ingest_run(tmp_path):
    db_path, _assembly_sha = _build_dataset(tmp_path)
    ledger = RunLedger(tmp_path / "sharur_ops.db")
    ledger.create_run("ingest", tmp_path, created_by="test")
    ledger.close()

    with pytest.raises(DatasetSealError, match="an ingest run is active"):
        build_dataset_seal(db_path)


def test_seal_cli_refuses_overwrite_and_verify_exits_nonzero_on_drift(tmp_path):
    db_path, _assembly_sha = _build_dataset(tmp_path)
    seal_path = tmp_path / "custom-seal.json"
    runner = CliRunner()

    created = runner.invoke(
        app,
        ["seal", "--db", str(db_path), "--output", str(seal_path)],
    )
    assert created.exit_code == 0, created.output
    assert seal_path.is_file()

    refused = runner.invoke(
        app,
        ["seal", "--db", str(db_path), "--output", str(seal_path)],
    )
    assert refused.exit_code == 1
    assert "already exists" in refused.output

    verified = runner.invoke(
        app,
        ["verify-seal", str(seal_path), "--db", str(db_path), "--format", "json"],
    )
    assert verified.exit_code == 0, verified.output
    assert json.loads(verified.output)["valid"] is True

    findings = tmp_path / "survey" / "findings.jsonl"
    findings.write_text(findings.read_text() + '{"id":"F2"}\n')
    drifted = runner.invoke(
        app,
        ["verify-seal", str(seal_path), "--db", str(db_path), "--format", "json"],
    )
    assert drifted.exit_code == 1
    assert json.loads(drifted.output)["valid"] is False
