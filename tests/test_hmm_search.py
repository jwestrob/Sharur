"""Stage 04's HMM-search CLI: Aksha first, legacy Astra fallback, registry adoption."""

from __future__ import annotations

import importlib.util
import json
import subprocess
from pathlib import Path

import pytest

from sharur import diagnostics, hmm_search

STAGE04 = Path(__file__).resolve().parents[1] / "src" / "ingest" / "04_astra_scan.py"


def _stage04():
    spec = importlib.util.spec_from_file_location("stage04_scan", STAGE04)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _registry(path: Path, entries: dict[str, str | None]) -> Path:
    """A registry with ``{name: installation_dir}``; None marks a catalog-only entry."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"db_path": "$HOME/.config/Astra", "db_urls": [
        {"name": n, "installed": d is not None, "installation_dir": d or "", "molecule_type": "protein"}
        for n, d in entries.items()]}))
    return path


@pytest.fixture
def registries(tmp_path, monkeypatch):
    aksha = tmp_path / "aksha" / "hmm_databases.json"
    astra = tmp_path / "astra" / "hmm_databases.json"
    monkeypatch.setattr(hmm_search, "registry_path", lambda cli: aksha if cli == "aksha" else astra)
    return aksha, astra


def test_resolver_prefers_aksha_then_falls_back_to_astra():
    both = {"aksha": "/bin/aksha", "astra": "/bin/astra"}
    assert hmm_search.resolve_cli(both.get, env={}) == ("aksha", "/bin/aksha")
    assert hmm_search.resolve_cli({"astra": "/bin/astra"}.get, env={}) == ("astra", "/bin/astra")
    assert hmm_search.resolve_cli({}.get, env={}) is None
    forced = {hmm_search.ENV_OVERRIDE: "Astra"}
    assert hmm_search.resolve_cli(both.get, env=forced) == ("astra", "/bin/astra")


def test_commands_carry_the_per_database_threshold_flags():
    def cmd(db, cutoffs=True):
        return hmm_search.build_search_command("aksha", db, Path("p"), Path("o"), 8, cutoffs)

    assert cmd("PFAM")[:2] == ["aksha", "search"]
    assert cmd("PFAM")[2:10] == ["--prot_in", "p", "--installed_hmms", "PFAM", "--outdir", "o", "--threads", "8"]
    assert cmd("PFAM")[10:] == ["--cut_ga"]
    assert cmd("DefenseFinder")[10:] == ["--cut_ga", "--write_macsyfinder"]
    assert cmd("TXSScan")[10:] == ["--write_macsyfinder"]
    assert cmd("KOFAM")[10:] == ["--cut_ga", "--cascade"]
    assert cmd("KOFAM", cutoffs=False)[10:] == []
    assert cmd("VOGdb")[10:] == []


def test_registration_problem_points_astra_only_databases_at_adoption(tmp_path, registries):
    aksha, astra = registries
    kofam = tmp_path / "KOFAM"
    kofam.mkdir()
    _registry(aksha, {"PFAM": str(tmp_path), "KOFAM": None, "VOGdb": None})
    _registry(astra, {"KOFAM": str(kofam)})

    assert hmm_search.registration_problem("aksha", "PFAM") is None
    assert "sharur adopt-astra-hmms" in hmm_search.registration_problem("aksha", "KOFAM")
    assert "aksha initialize --hmms VOGdb" in hmm_search.registration_problem("aksha", "VOGdb")
    assert hmm_search.registration_problem("astra", "VOGdb") is None


def test_adoption_registers_existing_directories_and_keeps_a_backup(tmp_path):
    kofam, txss = tmp_path / "KOFAM", tmp_path / "TXSScan"
    kofam.mkdir()
    aksha = _registry(tmp_path / "aksha.json", {"PFAM": "/elsewhere/PFAM", "KOFAM": None, "TXSScan": None})
    astra = _registry(tmp_path / "astra.json", {"PFAM": "/old/PFAM", "KOFAM": str(kofam), "TXSScan": str(txss)})

    assert hmm_search.adopt_astra_registry(aksha, astra, dry_run=True) == [("KOFAM", str(kofam))]
    assert not (tmp_path / "aksha.json.bak").exists()

    assert hmm_search.adopt_astra_registry(aksha, astra) == [("KOFAM", str(kofam))]
    assert hmm_search.installed_databases(aksha) == {"PFAM": "/elsewhere/PFAM", "KOFAM": str(kofam)}
    assert "KOFAM" not in hmm_search.installed_databases(tmp_path / "aksha.json.bak")
    assert hmm_search.adopt_astra_registry(aksha, astra) == []


def test_stage04_runs_aksha_and_reads_its_hits_table(tmp_path, registries, monkeypatch):
    aksha, astra = registries
    _registry(aksha, {"PFAM": str(tmp_path)})
    stage = _stage04()
    calls = []

    def fake_run(cmd, **kwargs):
        calls.append(cmd)
        out = Path(cmd[cmd.index("--outdir") + 1])
        (out / "PFAM_hits_df.tsv").write_text("sequence_id\thmm_name\np1\tA\np1\tB\np2\tA\n")
        return subprocess.CompletedProcess(cmd, 0, "", "")

    monkeypatch.setattr(stage.subprocess, "run", fake_run)
    result = stage.run_single_astra_scan("PFAM", tmp_path / "faa", tmp_path / "out", 4,
                                         cli=("aksha", "/bin/aksha"))

    assert calls[0][:2] == ["aksha", "search"] and "--cut_ga" in calls[0]
    assert result["execution_status"] == "success" and result["cli"] == "aksha"
    assert (result["total_hits"], result["unique_proteins"], result["unique_domains"]) == (3, 2, 2)


def test_stage04_reports_an_unregistered_database_without_searching(tmp_path, registries, monkeypatch):
    aksha, astra = registries
    _registry(aksha, {"KOFAM": None})
    stage = _stage04()
    monkeypatch.setattr(stage.subprocess, "run", lambda *a, **k: pytest.fail("searched"))

    result = stage.run_single_astra_scan("KOFAM", tmp_path, tmp_path / "out", 4, cli=("aksha", "/bin/aksha"))

    assert result["execution_status"] == "failed"
    assert "aksha initialize --hmms KOFAM" in result["error_message"]


def _labels(checks):
    return {c.label: c for c in checks}


def test_doctor_flags_missing_runtime_and_astra_only_databases(tmp_path, registries, monkeypatch):
    aksha, astra = registries
    pfam = tmp_path / "PFAM"
    pfam.mkdir()
    _registry(aksha, {"PFAM": None})
    _registry(astra, {"PFAM": str(pfam)})
    monkeypatch.setattr(hmm_search, "resolve_cli", lambda: ("aksha", "/bin/aksha"))
    monkeypatch.setattr(diagnostics, "find_spec", lambda name: None)

    checks = _labels(diagnostics.check_hmm_search())

    assert checks["aksha"].status == diagnostics.OK
    assert checks["aksha runtime"].status == diagnostics.MISSING
    assert checks["aksha registry"].status == diagnostics.MISSING
    assert "sharur adopt-astra-hmms" in checks["aksha registry"].detail


def test_doctor_accepts_legacy_astra_with_a_warning(tmp_path, registries, monkeypatch):
    aksha, astra = registries
    _registry(astra, {"PFAM": str(tmp_path)})
    monkeypatch.setattr(hmm_search, "resolve_cli", lambda: ("astra", "/bin/astra"))

    checks = diagnostics.check_hmm_search()

    assert [c.status for c in checks] == [diagnostics.WARN, diagnostics.OK]
    assert not diagnostics.has_core_failure(checks)
