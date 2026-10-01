"""The HMM-search CLI behind stage 04: Aksha, with its predecessor Astra as a fallback.

Aksha (github.com/jwestrob/aksha) is Astra renamed and extended; ``aksha search``
accepts every flag stage 04 passes (``--prot_in``, ``--installed_hmms``,
``--outdir``, ``--threads``, ``--cut_ga``, ``--cascade``, ``--write_macsyfinder``)
and writes the same ``{DATABASE}_hits_df.tsv`` per database.

The two tools share database *storage* (the packaged catalog's ``db_path`` is
``$HOME/.config/Astra`` for both) and registry schema, but each keeps its own
registry of installed databases. Aksha starts with an empty registry, so a
database installed through Astra is on disk and still unknown to Aksha until
its record is adopted (:func:`adopt_astra_registry`) or it is reinstalled with
``aksha initialize --hmms NAME``.

Stdlib only: ``aksha.initialize`` imports the compiled runtime, so the config
location is computed here with the same rules platformdirs applies.
"""

from __future__ import annotations

import json
import os
import shutil
import sys
from pathlib import Path
from typing import Callable

# Overrides the CLI choice: "aksha" or "astra".
ENV_OVERRIDE = "SHARUR_HMM_SEARCH"
CLIS = ("aksha", "astra")

LEGACY_DB_DIR = Path.home() / ".config" / "Astra"
REGISTRY_FILE = "hmm_databases.json"

# Databases whose profiles all carry GA thresholds.
GA_DATABASES = {"PFAM", "HYDDB", "DEFENSEFINDER"}
# Databases whose hits feed MacSyFinder --previous-run co-location validation.
MACSYFINDER_DATABASES = {"DEFENSEFINDER", "TXSSCAN"}


def aksha_config_dir() -> Path:
    """Aksha's config directory (``platformdirs.user_config_dir("Aksha")``)."""
    if sys.platform == "darwin":
        return Path.home() / "Library" / "Application Support" / "Aksha"
    if sys.platform.startswith("win"):
        return Path(os.environ.get("APPDATA") or Path.home() / "AppData" / "Roaming") / "Aksha"
    xdg = os.environ.get("XDG_CONFIG_HOME", "").strip()
    return (Path(xdg) if xdg else Path.home() / ".config") / "Aksha"


def registry_path(cli: str) -> Path:
    """The installed-database registry ``cli`` reads."""
    return (aksha_config_dir() if cli == "aksha" else LEGACY_DB_DIR) / REGISTRY_FILE


def installed_databases(path: Path) -> dict[str, str] | None:
    """``{name: installation_dir}`` for databases marked installed, or None if unreadable."""
    try:
        catalog = json.loads(Path(path).read_text())
    except (OSError, ValueError):
        return None
    return {
        entry["name"]: entry.get("installation_dir") or ""
        for entry in catalog.get("db_urls", [])
        if isinstance(entry, dict) and entry.get("name") and entry.get("installed")
    }


def resolve_cli(
    which: Callable[[str], str | None] = shutil.which,
    env: dict[str, str] | None = None,
) -> tuple[str, str] | None:
    """(name, path) of the HMM-search CLI: Aksha first, then Astra, or the override."""
    env = os.environ if env is None else env
    forced = (env.get(ENV_OVERRIDE) or "").strip().lower()
    for name in ((forced,) if forced in CLIS else CLIS):
        path = which(name)
        if path:
            return name, path
    return None


def build_search_command(
    cli: str,
    database: str,
    prot_in: Path,
    outdir: Path,
    threads: int,
    use_cutoffs: bool = True,
) -> list[str]:
    """The ``search`` invocation for one database; identical flags for Aksha and Astra."""
    cmd = [
        cli, "search",
        "--prot_in", str(prot_in),
        "--installed_hmms", database,
        "--outdir", str(outdir),
        "--threads", str(threads),
    ]
    db = database.upper()
    # PFAM, HydDB and DefenseFinder carry GA thresholds on every profile.
    # DefenseFinder's GA values are permissive, so the load-time e-value
    # filter in stage 07 stays as a safety net.
    if use_cutoffs and db in GA_DATABASES:
        cmd.append("--cut_ga")
    # MacSyFinder-compatible hmmsearch output for --previous-run co-location validation.
    if db in MACSYFINDER_DATABASES:
        cmd.append("--write_macsyfinder")
    # KOFAM: --cascade applies per-profile adaptive thresholds; --cut_ga alone
    # would apply one global threshold.
    if use_cutoffs and db == "KOFAM":
        cmd.extend(["--cut_ga", "--cascade"])
    # VOGdb and CANT-HYD have no GA thresholds; stage 07 filters them at load
    # time via DEFAULT_EVALUE_THRESHOLDS.
    return cmd


def registration_problem(cli: str, database: str) -> str | None:
    """Why ``cli`` cannot find ``database`` in its registry, with the fix; None when it can.

    An unreadable registry yields None so the CLI reports its own error.
    """
    if cli != "aksha":
        return None
    installed = installed_databases(registry_path("aksha"))
    if installed is None or installed.get(database):
        return None
    legacy = installed_databases(registry_path("astra")) or {}
    if legacy.get(database) and Path(legacy[database]).is_dir():
        return (
            f"{database} is installed at {legacy[database]} and registered only in Astra's "
            f"registry; adopt it into Aksha's with `sharur adopt-astra-hmms` "
            f"(or reinstall: aksha initialize --hmms {database})"
        )
    return f"{database} is not installed for Aksha; install it with: aksha initialize --hmms {database}"


def adopt_astra_registry(
    aksha_registry: Path | None = None,
    astra_registry: Path | None = None,
    dry_run: bool = False,
) -> list[tuple[str, str]]:
    """Mark Astra-installed databases installed in Aksha's registry.

    Copies ``installed`` and ``installation_dir`` for each database Astra lists
    as installed whose directory exists and that Aksha's catalog names but
    lists as uninstalled. Files stay where they are. Returns the adopted
    ``(name, installation_dir)`` pairs.
    """
    aksha_registry = Path(aksha_registry or registry_path("aksha"))
    astra_registry = Path(astra_registry or registry_path("astra"))
    legacy = installed_databases(astra_registry)
    if legacy is None:
        raise FileNotFoundError(f"no readable Astra registry at {astra_registry}")
    try:
        catalog = json.loads(aksha_registry.read_text())
    except (OSError, ValueError) as exc:
        raise FileNotFoundError(
            f"no readable Aksha registry at {aksha_registry}; "
            "create it with `aksha initialize --show_available`"
        ) from exc
    adopted = []
    for entry in catalog.get("db_urls", []):
        source = legacy.get(entry.get("name"))
        if entry.get("installed") or not source or not Path(source).is_dir():
            continue
        entry["installed"] = True
        entry["installation_dir"] = source
        adopted.append((entry["name"], source))
    if adopted and not dry_run:
        backup = aksha_registry.with_name(aksha_registry.name + ".bak")
        shutil.copy2(aksha_registry, backup)
        tmp = aksha_registry.with_name(aksha_registry.name + ".tmp")
        tmp.write_text(json.dumps(catalog, indent=4))
        os.replace(tmp, aksha_registry)
    return adopted
