"""Hermetic Stage-00 FASTA integrity contract tests."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
STAGE00 = REPO_ROOT / "src" / "ingest" / "00_prepare_inputs.py"


def _run_stage00(input_dir: Path, output_dir: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        [
            sys.executable,
            str(STAGE00),
            "--input-dir",
            str(input_dir),
            "--output-dir",
            str(output_dir),
        ],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )


def test_stage00_accepts_valid_set_and_writes_both_contracts(tmp_path):
    input_dir = tmp_path / "input"
    output_dir = tmp_path / "stage00"
    input_dir.mkdir()
    (input_dir / "bin-a.fna").write_text(">bin_a_contig_1\nACGTNN\n")
    (input_dir / "bin-b.fa").write_text(">bin_b_contig_1\nGGRYTT\n")

    result = _run_stage00(input_dir, output_dir)

    assert result.returncode == 0, result.stdout + result.stderr
    integrity = json.loads((output_dir / "input_integrity.json").read_text())
    manifest = json.loads((output_dir / "processing_manifest.json").read_text())
    assert integrity["status"] == "passed"
    assert integrity["summary"] == {
        "errors": 0,
        "files_scanned": 2,
        "records_scanned": 2,
        "total_bases": 12,
        "warnings": 0,
    }
    assert manifest["integrity_status"] == "passed"
    assert len(manifest["source_set_sha256"]) == 64
    assert all(len(genome["sha256"]) == 64 for genome in manifest["genomes"])
    assert (output_dir / "bin-a.fna").is_symlink()
    assert (output_dir / "bin-b.fa").is_symlink()


def test_stage00_rejects_terminal_header_without_sequence(tmp_path):
    input_dir = tmp_path / "input"
    output_dir = tmp_path / "stage00"
    input_dir.mkdir()
    (input_dir / "truncated.fna").write_text(">complete\nACGT\n>terminal_header")

    result = _run_stage00(input_dir, output_dir)

    assert result.returncode == 1
    integrity = json.loads((output_dir / "input_integrity.json").read_text())
    assert integrity["status"] == "failed"
    assert "empty_sequence" in {issue["code"] for issue in integrity["issues"]}
    assert not (output_dir / "processing_manifest.json").exists()
    assert not (output_dir / "truncated.fna").exists()


def test_stage00_rejects_record_ids_repeated_across_assemblies(tmp_path):
    input_dir = tmp_path / "input"
    output_dir = tmp_path / "stage00"
    input_dir.mkdir()
    (input_dir / "first.fna").write_text(">shared_contig\nACGT\n")
    (input_dir / "second.fna").write_text(">shared_contig\nTGCA\n")

    result = _run_stage00(input_dir, output_dir)

    assert result.returncode == 1
    integrity = json.loads((output_dir / "input_integrity.json").read_text())
    issue = next(
        item
        for item in integrity["issues"]
        if item["code"] == "duplicate_record_id_across_files"
    )
    assert issue["record_id"] == "shared_contig"
    assert issue["related_file"] == "first.fna"


def test_stage00_rejects_filename_normalization_collisions(tmp_path):
    input_dir = tmp_path / "input"
    output_dir = tmp_path / "stage00"
    input_dir.mkdir()
    (input_dir / "same-bin.fna").write_text(">contig_a\nACGT\n")
    (input_dir / "same_bin.fa").write_text(">contig_b\nTGCA\n")

    result = _run_stage00(input_dir, output_dir)

    assert result.returncode == 1
    integrity = json.loads((output_dir / "input_integrity.json").read_text())
    assert "duplicate_genome_id" in {
        issue["code"] for issue in integrity["issues"]
    }
