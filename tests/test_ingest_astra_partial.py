"""Stage 04 keeps an interrupted search's hits out of the files stage 07 loads."""

import importlib.util
import os
import stat
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location(
    "stage04", Path(__file__).resolve().parents[1] / "src" / "ingest" / "04_astra_scan.py")
stage04 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(stage04)

FAKE_CLI = """#!/bin/sh
# Writes one genome's hits to tmp_results/, as a per-genome search does, then exits with $FAKE_EXIT.
out=""
while [ $# -gt 0 ]; do [ "$1" = "--outdir" ] && out="$2"; shift; done
mkdir -p "$out/tmp_results"
printf 'sequence_id\\thmm_name\\tbitscore\\n' > "$out/tmp_results/g1_results.tsv"
printf 'g1_1\\tDolB\\t40.0\\n' >> "$out/tmp_results/g1_results.tsv"
exit "${FAKE_EXIT:-0}"
"""


@pytest.fixture
def search_env(tmp_path, monkeypatch):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    cli = bin_dir / "fakesearch"
    cli.write_text(FAKE_CLI)
    cli.chmod(cli.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    proteins = tmp_path / "all_protein_symlinks"
    proteins.mkdir()
    for g in ("g1", "g2", "g3"):
        (proteins / f"{g}.faa").write_text(">x\nM\n")
    return proteins, tmp_path / "out", ("fakesearch", str(cli))


@pytest.mark.parametrize("exit_code, status, hits_name", [
    ("0", "success", "DefenseFinder_hits_df.tsv"),
    ("1", "partial", "DefenseFinder_hits_df.partial.tsv"),
])
def test_only_a_finished_search_writes_the_loadable_hits_file(search_env, monkeypatch, exit_code, status, hits_name):
    proteins, out, cli = search_env
    monkeypatch.setenv("FAKE_EXIT", exit_code)

    result = stage04.run_single_astra_scan("DefenseFinder", proteins, out, threads=1, cli=cli)

    results_dir = out / "defensefinder_results"
    assert result["execution_status"] == status and Path(result["hits_file"]).name == hits_name
    assert (results_dir / "DefenseFinder_hits_df.tsv").exists() == (status == "success")
    if status == "partial":
        assert "after 1 of 3 genomes" in result["error_message"]
