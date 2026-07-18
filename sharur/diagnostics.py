"""Install verification for Sharur — the backing logic for ``sharur doctor``.

This module centralizes the list of external tools and reference databases the
pipeline depends on. The tool set here mirrors what the ingest pipeline shells
out to (see the ``--skip-*`` flags in ``sharur/ingest_cli.py`` and the stage
scripts under ``src/ingest/``); the reference-DB locations mirror the hardcoded
paths in ``sharur/operators/foldseek.py`` and ``sharur/colocation.py``.

Kept deliberately dependency-light (stdlib only) so it stays import-safe and can
run even in a partially-provisioned environment.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from importlib.metadata import PackageNotFoundError, distribution
from pathlib import Path


# Result status values.
OK = "ok"
WARN = "warn"
MISSING = "missing"

_VERSION_TIMEOUT_S = 10


@dataclass
class Check:
    """The outcome of a single diagnostic check."""

    label: str
    status: str  # OK / WARN / MISSING
    detail: str
    core: bool
    purpose: str = ""


@dataclass
class ToolSpec:
    """An external binary the pipeline expects on ``$PATH``."""

    name: str
    binaries: tuple[str, ...]  # candidate names to look up via shutil.which
    version_args: tuple[str, ...]
    core: bool
    purpose: str


# Core = required for the default ingest/annotation path. Optional tools power
# opt-in QC, BGC, CAZy, structure, and validated-system extensions.
TOOLS: tuple[ToolSpec, ...] = (
    ToolSpec("prodigal", ("prodigal",), ("-v",), True, "gene calling (stage 03)"),
    ToolSpec("diamond", ("diamond",), ("version",), True, "protein alignment"),
    ToolSpec("hmmsearch", ("hmmsearch",), ("-h",), True, "HMM search (HMMER)"),
    ToolSpec("astra", ("astra",), ("--version",), True, "HMM annotation (stage 04)"),
    ToolSpec("quast", ("quast.py", "quast"), ("--version",), False, "assembly QC (stage 01)"),
    ToolSpec("dfast_qc", ("dfast_qc",), ("--version",), False, "assembly QC (stage 02)"),
    ToolSpec("minced", ("minced",), ("--version",), True, "CRISPR arrays (stage 05c)"),
    ToolSpec("gecco", ("gecco",), ("--version",), False, "BGC detection"),
    ToolSpec("run_dbcan", ("run_dbcan",), ("--version",), False, "CAZyme annotation"),
    ToolSpec("foldseek", ("foldseek",), ("version",), False, "structural homology search"),
    ToolSpec(
        "defense-finder",
        ("defense-finder", "defense_finder"),
        ("--version",),
        False,
        "validated defense systems",
    ),
)

# Reference-database locations (mirror the paths hardcoded elsewhere in the
# codebase; see module docstring). These are provisioned out-of-band.
ASTRA_DB_DIR = Path.home() / ".config" / "Astra"
FOLDSEEK_DB_DIR = Path.home() / ".foldseek"
MACSY_MODEL_DIRS = (
    Path.home() / ".macsyfinder" / "models",
    Path.home() / ".mdmlab" / "macsyfinder" / "models",
)


def _probe_version(binary: str, args: tuple[str, ...]) -> str:
    """Best-effort version string for ``binary``. Never raises."""
    try:
        proc = subprocess.run(  # noqa: S603 - trusted local tool invocation
            [binary, *args],
            capture_output=True,
            text=True,
            timeout=_VERSION_TIMEOUT_S,
        )
    except (OSError, subprocess.SubprocessError):
        return "present"
    output = f"{proc.stdout}\n{proc.stderr}"
    for line in output.splitlines():
        line = line.strip()
        if line and re.search(r"\d", line):
            return line[:48]
    return "present"


def check_tool(spec: ToolSpec) -> Check:
    """Resolve one tool: found → OK (+version); absent → MISSING (core) / WARN."""
    for binary in spec.binaries:
        resolved = shutil.which(binary)
        if resolved:
            return Check(
                label=spec.name,
                status=OK,
                detail=_probe_version(binary, spec.version_args),
                core=spec.core,
                purpose=spec.purpose,
            )
    return Check(
        label=spec.name,
        status=MISSING if spec.core else WARN,
        detail="not found on PATH",
        core=spec.core,
        purpose=spec.purpose,
    )


def _check_dir_db(
    label: str,
    path: Path,
    purpose: str,
    *,
    core: bool = False,
) -> Check:
    """Report a reference-DB directory as OK (with contents) or WARN if absent."""
    if path.is_dir():
        entries = sorted(p.name for p in path.iterdir() if not p.name.startswith("."))
        if entries:
            listed = ", ".join(entries[:6]) + (" …" if len(entries) > 6 else "")
            return Check(label, OK, f"{path} ({listed})", core, purpose)
        return Check(
            label,
            MISSING if core else WARN,
            f"{path} (empty)",
            core,
            purpose,
        )
    return Check(
        label,
        MISSING if core else WARN,
        f"not found at {path}",
        core,
        purpose,
    )


def check_ingest_entrypoint() -> Check:
    """Detect stale editable installs that omit the primary ingest command."""
    executable = shutil.which("sharur-ingest")
    try:
        installed = distribution("sharur")
        entrypoints = {
            entry.name
            for entry in installed.entry_points
            if entry.group == "console_scripts"
        }
        version = installed.version
    except PackageNotFoundError:
        entrypoints = set()
        version = "not installed"

    if executable and "sharur-ingest" in entrypoints:
        return Check(
            "sharur-ingest",
            OK,
            f"{executable} (package {version})",
            True,
            "primary staged-ingest CLI",
        )

    details = []
    if not executable:
        details.append("executable not on PATH")
    if "sharur-ingest" not in entrypoints:
        details.append(f"package {version} has no console entrypoint")
    details.append('repair: pip install -e ".[dev]"')
    return Check(
        "sharur-ingest",
        MISSING,
        "; ".join(details),
        True,
        "primary staged-ingest CLI",
    )


def check_reference_dbs() -> list[Check]:
    checks = [
        _check_dir_db(
            "Astra HMMs",
            ASTRA_DB_DIR,
            "annotation HMM databases",
            core=True,
        ),
        _check_dir_db("Foldseek DBs", FOLDSEEK_DB_DIR, "structure search databases"),
    ]
    macsy = next((d for d in MACSY_MODEL_DIRS if d.is_dir()), None)
    if macsy is not None:
        checks.append(_check_dir_db("DefenseFinder models", macsy, "MacSyFinder model definitions"))
    else:
        checks.append(
            Check(
                "DefenseFinder models",
                WARN,
                f"not found at {MACSY_MODEL_DIRS[0]}",
                False,
                "MacSyFinder model definitions",
            )
        )
    return checks


def check_api_keys() -> list[Check]:
    if os.environ.get("ESM_API_KEY"):
        detail, status = "set", OK
    else:
        detail, status = "unset (structure prediction disabled)", WARN
    return [Check("ESM_API_KEY", status, detail, False, "ESM3 structure prediction API")]


def run_all_checks() -> list[Check]:
    """Run every diagnostic and return the flat list of results."""
    results: list[Check] = [check_ingest_entrypoint()]
    results.extend(check_tool(spec) for spec in TOOLS)
    results.extend(check_reference_dbs())
    results.extend(check_api_keys())
    return results


def has_core_failure(checks: list[Check]) -> bool:
    """True if any *core* component is missing (used by ``doctor --strict``)."""
    return any(c.core and c.status == MISSING for c in checks)
