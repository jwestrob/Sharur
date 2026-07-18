"""Tests for install diagnostics used by ``sharur doctor``."""

from types import SimpleNamespace

from sharur import diagnostics


def test_ingest_entrypoint_detects_stale_editable_install(monkeypatch):
    monkeypatch.setattr(diagnostics.shutil, "which", lambda _name: None)
    monkeypatch.setattr(
        diagnostics,
        "distribution",
        lambda _name: SimpleNamespace(version="0.1.0", entry_points=[]),
    )

    check = diagnostics.check_ingest_entrypoint()

    assert check.status == diagnostics.MISSING
    assert check.core is True
    assert "pip install -e" in check.detail


def test_ingest_entrypoint_passes_when_metadata_and_path_agree(monkeypatch):
    entrypoint = SimpleNamespace(
        name="sharur-ingest",
        group="console_scripts",
    )
    monkeypatch.setattr(
        diagnostics.shutil,
        "which",
        lambda _name: "/tmp/bin/sharur-ingest",
    )
    monkeypatch.setattr(
        diagnostics,
        "distribution",
        lambda _name: SimpleNamespace(
            version="0.1.0",
            entry_points=[entrypoint],
        ),
    )

    check = diagnostics.check_ingest_entrypoint()

    assert check.status == diagnostics.OK
    assert check.core is True


def test_empty_core_reference_directory_is_a_core_failure(tmp_path):
    check = diagnostics._check_dir_db(
        "Astra HMMs",
        tmp_path,
        "annotation HMM databases",
        core=True,
    )

    assert check.status == diagnostics.MISSING
    assert diagnostics.has_core_failure([check]) is True


def test_minced_is_core_because_stage05c_runs_by_default():
    spec = next(spec for spec in diagnostics.TOOLS if spec.name == "minced")

    assert spec.core is True
