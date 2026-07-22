"""Regression tests for local Foldseek discovery and result contracts."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

import sharur.operators.foldseek as foldseek_operator
from sharur.operators.base import OperatorContext


def _write_executable(path: Path, body: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"#!/bin/sh\n{body}\n")
    path.chmod(0o755)
    return path


def test_find_foldseek_binary_prefers_valid_environment_override(
    tmp_path,
    monkeypatch,
):
    override = _write_executable(
        tmp_path / "override" / "foldseek",
        'test "$1" = version || exit 2\necho override-version',
    )
    path_binary = _write_executable(
        tmp_path / "path" / "foldseek",
        'test "$1" = version || exit 2\necho path-version',
    )
    monkeypatch.setenv("FOLDSEEK_BINARY", str(override))
    monkeypatch.setenv("PATH", str(path_binary.parent))
    monkeypatch.setenv("PYENV_ROOT", str(tmp_path / "empty-pyenv"))

    assert foldseek_operator._find_foldseek_binary() == str(override.resolve())


def test_find_foldseek_binary_skips_broken_pyenv_shim(
    tmp_path,
    monkeypatch,
):
    shim = _write_executable(
        tmp_path / "shims" / "foldseek",
        'echo "pyenv: foldseek: command not found" >&2\nexit 127',
    )
    pyenv_root = tmp_path / "pyenv"
    real_binary = _write_executable(
        pyenv_root / "versions" / "miniconda" / "envs" / "phold" / "bin" / "foldseek",
        'test "$1" = version || exit 2\necho 10.0',
    )
    monkeypatch.delenv("FOLDSEEK_BINARY", raising=False)
    monkeypatch.delenv("FOLDSEEK_BIN", raising=False)
    monkeypatch.setenv("PATH", str(shim.parent))
    monkeypatch.setenv("PYENV_ROOT", str(pyenv_root))

    assert foldseek_operator._find_foldseek_binary() == str(real_binary.resolve())


def test_find_local_database_resolves_symlinked_prefix(
    tmp_path,
    monkeypatch,
):
    real_prefix = tmp_path / "databases" / "pdb"
    real_prefix.parent.mkdir()
    real_prefix.write_bytes(b"database")
    Path(f"{real_prefix}.dbtype").write_bytes(b"\x00\x00\x00\x00")

    local_dir = tmp_path / "local"
    local_dir.mkdir()
    (local_dir / "pdb100").symlink_to(real_prefix)
    monkeypatch.setattr(foldseek_operator, "LOCAL_DB_DIR", local_dir)

    assert foldseek_operator._find_local_database("pdb100") == real_prefix.resolve()


def test_find_local_database_rejects_prefix_without_dbtype(
    tmp_path,
    monkeypatch,
):
    local_dir = tmp_path / "local"
    prefix = local_dir / "pdb100" / "pdb100"
    prefix.parent.mkdir(parents=True)
    prefix.write_bytes(b"incomplete database")
    monkeypatch.setattr(foldseek_operator, "LOCAL_DB_DIR", local_dir)

    assert foldseek_operator._find_local_database("pdb100") is None


def test_search_foldseek_local_parses_identity_and_coverage(
    tmp_path,
    monkeypatch,
):
    pdb_path = tmp_path / "query.pdb"
    pdb_path.write_text("END\n")
    database_path = tmp_path / "pdb"
    database_path.write_bytes(b"database")

    monkeypatch.setattr(
        foldseek_operator,
        "_find_foldseek_binary",
        lambda: "/opt/foldseek/bin/foldseek",
    )
    monkeypatch.setattr(
        foldseek_operator,
        "_find_local_database",
        lambda database: database_path,
    )

    def fake_run(cmd, *, capture_output, check, text):
        assert capture_output is True
        assert check is False
        assert text is True
        columns = cmd[cmd.index("--format-output") + 1].split(",")
        assert {"fident", "qcov", "tcov"}.issubset(columns)
        row = {
            "query": "query",
            "target": "1abc_A",
            "evalue": "1e-20",
            "bits": "123.4",
            "fident": "0.42",
            "qcov": "0.80",
            "tcov": "0.60",
            "qstart": "2",
            "qend": "101",
            "tstart": "4",
            "tend": "103",
            "taxname": "Test taxon",
            "theader": "Test protein",
        }
        Path(cmd[4]).write_text("\t".join(row[column] for column in columns) + "\n")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(foldseek_operator.subprocess, "run", fake_run)

    hits = foldseek_operator.search_foldseek_local(
        str(pdb_path),
        database="pdb100",
        top_k=5,
    )

    assert hits is not None
    assert len(hits) == 1
    assert hits[0]["seq_identity"] == pytest.approx(0.42)
    assert hits[0]["qcov"] == pytest.approx(0.80)
    assert hits[0]["tcov"] == pytest.approx(0.60)
    assert hits[0]["qstart"] == 2
    assert hits[0]["description"] == "Test protein"


def test_search_foldseek_returns_formatted_data_and_raw_hits(
    tmp_path,
    monkeypatch,
):
    pdb_path = tmp_path / "query.pdb"
    pdb_path.write_text("END\n")
    hit = {
        "database": "pdb100",
        "target": "1abc_A",
        "evalue": 1e-20,
        "score": 123.4,
        "seq_identity": 0.42,
        "qcov": 0.80,
        "tcov": 0.60,
        "description": "Test protein",
        "taxon": "Test taxon",
        "qstart": 2,
        "qend": 101,
        "tstart": 4,
        "tend": 103,
    }
    monkeypatch.setattr(
        foldseek_operator,
        "search_foldseek_local",
        lambda pdb_path, database, top_k: [hit],
    )

    result = foldseek_operator.search_foldseek(
        str(pdb_path),
        databases=["pdb100"],
        top_k=5,
    )

    assert isinstance(result.data, str)
    assert "1abc_A" in result.data
    assert result.raw == [hit]
    assert result.records == [hit]
    assert result.meta.rows == 1


def test_search_foldseek_returns_explicit_empty_payload(
    tmp_path,
    monkeypatch,
):
    pdb_path = tmp_path / "query.pdb"
    pdb_path.write_text("END\n")
    monkeypatch.setattr(
        foldseek_operator,
        "search_foldseek_local",
        lambda pdb_path, database, top_k: [],
    )

    result = foldseek_operator.search_foldseek(
        str(pdb_path),
        databases=["pdb100"],
    )

    assert result.data == "No structural homologs found."
    assert result.raw == []
    assert result.records == []
    assert result.meta.rows == 0
    assert result.status == "empty"


def test_search_foldseek_for_protein_tags_raw_hits(
    tmp_path,
    monkeypatch,
):
    pdb_path = tmp_path / "query.pdb"
    pdb_path.write_text("END\n")
    hit = {"target": "1abc_A", "evalue": 1e-20}
    with OperatorContext("search_foldseek", {}) as context:
        search_result = context.make_result(
            data="formatted hits",
            rows=1,
            raw=[hit],
        )
    monkeypatch.setattr(
        foldseek_operator,
        "search_foldseek",
        lambda pdb_path, databases, top_k: search_result,
    )

    result = foldseek_operator.search_foldseek_for_protein(
        None,
        "protein_1",
        pdb_path=str(pdb_path),
    )

    assert result.raw[0]["query_protein_id"] == "protein_1"


def test_list_foldseek_databases_returns_formatted_data_and_raw_records(
    monkeypatch,
):
    databases = [
        {"name": "pdb100", "description": "Protein Data Bank"},
        {"name": "afdb50", "description": "AlphaFold DB clustered at 50%"},
    ]

    class Response:
        def json(self):
            return databases

    monkeypatch.setattr(
        foldseek_operator.requests,
        "get",
        lambda url, timeout: Response(),
    )

    result = foldseek_operator.list_foldseek_databases()

    assert isinstance(result.data, str)
    assert "pdb100" in result.data
    assert result.raw == databases
    assert result.records == databases
