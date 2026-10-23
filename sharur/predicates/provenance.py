"""Which predicate maps produced a database's predicates.

Every predicate generation (full or subset) appends a row to
``predicate_provenance`` recording the maps and rules in force: the shipped
Pfam and CAZy maps, the VOG rules, the locally built KEGG map (``sharur setup-kegg``), the vocabulary,
the V2 configuration, and the semantic fingerprint. :func:`map_status`
compares a database's stamps with the maps installed now, so preflight can
report predicates that an older map produced.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from sharur import __version__


PACKAGE = Path(__file__).resolve().parents[1]
PFAM_MAP = PACKAGE / "predicates/mappings/data/pfam_predicates.tsv"
CAZY_MAP = PACKAGE / "predicates/mappings/data/cazy_predicates.tsv"
VOG_RULES = PACKAGE / "predicates/mappings/vog_map.py"
VOCABULARY = PACKAGE / "predicates/vocabulary.py"
# Source checkout keeps the V2 config at config/predicates_v2; wheels package it
# at sharur/predicates_v2/config.
V2_CONFIG = next((path for path in (PACKAGE.parent / "config/predicates_v2", PACKAGE / "predicates_v2/config")
                  if path.is_dir()), PACKAGE / "predicates_v2/config")
KEGG_RULES = (
    "kegg_brite_predicates.tsv",
    "kegg_module_predicates.tsv",
    "kegg_predicate_proposals.tsv",
    "kegg_predicate_proposal_patterns.tsv",
    "kegg_hyddb_snapshot.tsv",
    "kegg_swissprot_consensus.tsv",
)

PROVENANCE_COLUMNS = (
    "generation_id BIGINT PRIMARY KEY",
    "generated_at TIMESTAMP NOT NULL",
    "scope VARCHAR NOT NULL",  # full | subset
    "protein_count BIGINT NOT NULL",
    "semantic_fingerprint VARCHAR",
    "pfam_map_sha256 VARCHAR NOT NULL",
    "pfam_sources VARCHAR",
    "kegg_map_sha256 VARCHAR",  # NULL: no local KEGG build
    "kegg_release VARCHAR",
    "kegg_rules_sha256 VARCHAR NOT NULL",
    "cazy_map_sha256 VARCHAR NOT NULL",
    "vog_rules_sha256 VARCHAR NOT NULL",
    "vocabulary_sha256 VARCHAR NOT NULL",
    "v2_config_sha256 VARCHAR NOT NULL",
    "sharur_version VARCHAR NOT NULL",
    "git_commit VARCHAR",
    "reconciled_from BIGINT",          # scope=reconciled: the base generation it brings current
    "reconciled_accessions VARCHAR",   # scope=reconciled: JSON list of accessions whose mapping changed
)
CREATE_SQL = f"CREATE TABLE IF NOT EXISTS predicate_provenance ({', '.join(PROVENANCE_COLUMNS)})"
_LATER_COLUMNS = ("reconciled_from BIGINT", "reconciled_accessions VARCHAR")

# Fields compared by map_status; everything else is context.
COMPARED = ("pfam_map_sha256", "kegg_map_sha256", "kegg_rules_sha256", "cazy_map_sha256", "vog_rules_sha256",
            "vocabulary_sha256", "v2_config_sha256")
# Generations that describe the whole database: a full run, or a verified reconcile of one (see reconcile()).
BASE_SCOPES = ("full", "reconciled")
# Map fields a reconcile can bring current. Predicate generation reads KEGG only through the built map
# (kegg_predicates.tsv); the KEGG rule files act only through that build, so a KO-level diff of the built
# map covers both. Map-derived predicates are evaluated per protein; locus atoms (nitrogenase_loci) and
# validated-system atoms use fixed KO sets and caller tables, not the map. Regenerating the proteins that
# carry a changed KO is therefore complete. Any other changed field needs a full regeneration.
RECONCILABLE = ("kegg_map_sha256", "kegg_rules_sha256")


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _sha256_files(paths) -> str:
    digest = hashlib.sha256()
    for path in paths:
        if path.is_file():
            digest.update(f"{path.name}\0".encode())
            digest.update(path.read_bytes())
            digest.update(b"\0")
    return digest.hexdigest()


def _pfam_sources(path: Path) -> str:
    """The source lines of the Pfam map header (Pfam, GO, ENZYME, Swiss-Prot releases)."""
    lines = []
    with open(path) as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            text = line[1:].strip()
            if text.split(":", 1)[0] in ("Pfam", "GO", "ENZYME", "Swiss-Prot"):
                lines.append(text)
    return " | ".join(lines)


def _git_commit() -> str | None:
    try:
        out = subprocess.run(["git", "-C", str(PACKAGE.parent), "rev-parse", "HEAD"],
                             capture_output=True, text=True, timeout=5, check=False)
    except (OSError, subprocess.SubprocessError):
        return None
    return out.stdout.strip() or None


def current_maps() -> dict[str, Any]:
    """The maps and rules predicate generation would use right now."""
    from sharur.predicates.mappings.kegg_map import KEGG_MAPPING_FILE

    kegg_sha = kegg_release = None
    if KEGG_MAPPING_FILE is not None:
        kegg_sha = _sha256_file(KEGG_MAPPING_FILE)
        prov = KEGG_MAPPING_FILE.with_name("provenance.json")
        if prov.exists():
            kegg_release = json.loads(prov.read_text()).get("kegg_release")
    data = PFAM_MAP.parent
    return {
        "pfam_map_sha256": _sha256_file(PFAM_MAP),
        "pfam_sources": _pfam_sources(PFAM_MAP),
        "kegg_map_sha256": kegg_sha,
        "kegg_release": kegg_release,
        "kegg_rules_sha256": _sha256_files(data / name for name in KEGG_RULES),
        "cazy_map_sha256": _sha256_file(CAZY_MAP),
        "vog_rules_sha256": _sha256_file(VOG_RULES),
        "vocabulary_sha256": _sha256_file(VOCABULARY),
        "v2_config_sha256": _sha256_files(sorted(V2_CONFIG.glob("*.yaml"))) if V2_CONFIG.is_dir() else "",
        "sharur_version": __version__,
        "git_commit": _git_commit(),
    }


def ensure_table(store) -> None:
    store.execute(CREATE_SQL)
    for column in _LATER_COLUMNS:  # tables created before these columns existed
        store.execute(f"ALTER TABLE predicate_provenance ADD COLUMN IF NOT EXISTS {column}")


def record_generation(store, *, scope: str, protein_count: int, semantic_fingerprint: str | None = None,
                      reconciled_from: int | None = None, reconciled_accessions: list[str] | None = None) -> dict:
    """Append one provenance row for a predicate generation and return it."""
    ensure_table(store)
    row = {
        "generated_at": datetime.now(timezone.utc).replace(tzinfo=None),
        "scope": scope,
        "protein_count": protein_count,
        "semantic_fingerprint": semantic_fingerprint or None,
        **current_maps(),
        "reconciled_from": reconciled_from,
        "reconciled_accessions": json.dumps(reconciled_accessions) if reconciled_accessions is not None else None,
    }
    next_id = store.execute("SELECT COALESCE(MAX(generation_id), 0) + 1 FROM predicate_provenance")[0][0]
    columns = ["generation_id", *row]
    store.conn.execute(
        f"INSERT INTO predicate_provenance ({', '.join(columns)}) VALUES ({', '.join('?' * len(columns))})",
        [next_id, *row.values()],
    )
    return {"generation_id": next_id, **row}


def latest_stamps(store) -> dict[str, Any] | None:
    """The latest base generation (full, or a verified reconcile) and any subset generations after it, or None."""
    tables = {r[0] for r in store.execute(
        "SELECT table_name FROM information_schema.tables WHERE table_name = 'predicate_provenance'")}
    if not tables:
        return None
    present = {r[0] for r in store.execute(
        "SELECT column_name FROM information_schema.columns WHERE table_name = 'predicate_provenance'")}
    columns = [c.split()[0] for c in PROVENANCE_COLUMNS]
    selected = [c for c in columns if c in present]
    rows = []
    for r in store.execute(f"SELECT {', '.join(selected)} FROM predicate_provenance ORDER BY generation_id"):
        row = dict.fromkeys(columns)
        row.update(zip(selected, r, strict=True))
        rows.append(row)
    base = [r for r in rows if r["scope"] in BASE_SCOPES]
    if not base:
        return {"full": None, "subsets": rows}
    last = base[-1]
    return {"full": last, "subsets": [r for r in rows if r["generation_id"] > last["generation_id"]]}


@dataclass(frozen=True)
class MapStatus:
    state: str  # current | stale | unstamped
    changed: tuple[str, ...] = ()
    stamp: dict[str, Any] | None = None
    installed: dict[str, Any] | None = None
    mixed_subsets: int = 0  # later subset generations that used other maps than the full one


def map_status(store) -> MapStatus:
    """Compare the database's latest full-generation stamp with the installed maps."""
    installed = current_maps()
    stamps = latest_stamps(store)
    if not stamps or stamps["full"] is None:
        return MapStatus("unstamped", installed=installed)
    stamp = stamps["full"]
    changed = tuple(field for field in COMPARED if stamp.get(field) != installed.get(field))
    mixed = sum(1 for sub in stamps["subsets"] if any(sub.get(f) != stamp.get(f) for f in COMPARED))
    return MapStatus("stale" if changed else "current", changed, stamp, installed, mixed)


def kegg_map_changes(old_map: Path, new_map: Path) -> list[str]:
    """KOs whose predicate set differs between two built KEGG maps (added, removed or changed)."""
    from sharur.predicates.mappings.kegg_map import _load_kegg_mapping_file

    old, _ = _load_kegg_mapping_file(Path(old_map))
    new, _ = _load_kegg_mapping_file(Path(new_map))
    return sorted(ko for ko in old.keys() | new.keys() if old.get(ko, []) != new.get(ko, []))


class ReconcileError(RuntimeError):
    """The database cannot be brought current by a scoped regeneration."""


def reconcile(store, old_kegg_map: Path, regenerate) -> dict[str, Any]:
    """Bring a database's predicates current after a KEGG-map-only change, without a full regeneration.

    Verifies that the database's latest base stamp differs from the installed maps only in KEGG fields and
    that ``old_kegg_map`` is the map that stamp recorded (by sha256). Then regenerates exactly the proteins
    carrying a KO whose predicate set changed (``regenerate(protein_ids)``, normally
    ``generate_v2(store, protein_ids, update_legacy_predicates=True)``) and records a ``reconciled``
    generation, which map_status treats as current. Raises ReconcileError instead of guessing.
    """
    from sharur.predicates.mappings.kegg_map import KEGG_MAPPING_FILE

    stamps = latest_stamps(store)
    if not stamps or stamps["full"] is None:
        raise ReconcileError("no full generation stamp: run a full `sharur compute-predicates` first")
    base = stamps["full"]
    installed = current_maps()
    changed = [f for f in COMPARED if base.get(f) != installed.get(f)]
    if not changed:
        return {"state": "current", "base_generation": base["generation_id"], "proteins": 0, "accessions": []}
    blocking = [f for f in changed if f not in RECONCILABLE]
    if blocking:
        raise ReconcileError(f"{', '.join(blocking)} changed; only a full regeneration covers that")
    if KEGG_MAPPING_FILE is None:
        raise ReconcileError("no local KEGG build installed")
    old_sha = _sha256_file(Path(old_kegg_map))
    if old_sha != base.get("kegg_map_sha256"):
        raise ReconcileError(f"{old_kegg_map} (sha256 {old_sha[:16]}) is not the KEGG map the database was "
                             f"stamped with ({str(base.get('kegg_map_sha256'))[:16]})")
    kos = kegg_map_changes(Path(old_kegg_map), KEGG_MAPPING_FILE)
    proteins: list[str] = []
    if kos:
        proteins = [r[0] for r in store.execute(
            "SELECT DISTINCT protein_id FROM annotations WHERE source = 'kegg' AND accession IN "
            "(SELECT UNNEST(?::VARCHAR[])) ORDER BY 1", [kos])]
    if proteins:
        regenerate(proteins)
    row = record_generation(store, scope="reconciled", protein_count=len(proteins),
                            semantic_fingerprint=base.get("semantic_fingerprint"),
                            reconciled_from=base["generation_id"], reconciled_accessions=kos)
    return {"state": "reconciled", "base_generation": base["generation_id"], "generation_id": row["generation_id"],
            "changed_fields": changed, "accessions": kos, "proteins": len(proteins)}
