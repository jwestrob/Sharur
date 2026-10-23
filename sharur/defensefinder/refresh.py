"""Staged, validated reload of DefenseFinder evidence from a per-domain hits table.

The reload never edits the database file it is given. It copies the database to a
staging file, then in the copy:

1. replaces every ``defensefinder`` annotation row with the hits selected from the
   table by MacSyFinder/DefenseFinder rules (:mod:`sharur.defensefinder.selection`),
   and records the model set and parameters in ``annotation_source_versions``;
2. reruns the co-location caller (``sharur.colocation``), which replaces the
   ``defense_systems`` table, ``defensefinder_system`` annotations and their
   ``system_proteins`` rows;
3. regenerates V2 semantics and the V1 compatibility cache for every protein whose
   DefenseFinder evidence existed before or exists after;
4. validates that everything outside that scope is unchanged.

A dry run reports the difference and discards the copy. ``publish`` atomically
replaces the database (the staging file must live on the same filesystem), keeps
the prior file as a hard-linked (or copied) backup beside it, and re-seals when a
dataset seal exists.
"""

from __future__ import annotations

import json
import logging
import os
import re
import shutil
import time
from collections import Counter
from dataclasses import asdict, dataclass, field
from pathlib import Path

import duckdb
import pandas as pd

from sharur.defensefinder.selection import (
    SOURCE,
    SOURCE_VERSIONS_TABLE,
    ModelSet,
    SelectionParams,
    SelectionStats,
    annotation_rows,
    dedup_fam_variants,
    file_sha256,
    load_model_set,
    read_hits,
    record_source_version,
    select_hits,
)


SYSTEM_SOURCE = "defensefinder_system"
SYSTEM_TABLE = "defense_systems"
_PROTEIN_SCOPED = ("semantic_atoms", "semantic_state", "semantic_terms", "protein_predicates")
_REPLACED = (SYSTEM_TABLE, SOURCE_VERSIONS_TABLE)
_APPENDED = "predicate_provenance"


logger = logging.getLogger(__name__)


class ReloadValidationError(RuntimeError):
    """The staged database failed validation; the source database was left unchanged."""


@dataclass
class Snapshot:
    """DefenseFinder state of one database."""

    annotation_rows: int = 0
    annotation_proteins: int = 0
    system_annotation_rows: int = 0
    systems: int = 0
    genomes_with_system: int = 0
    systems_by_type: dict[str, int] = field(default_factory=dict)
    genomes_by_type: dict[str, int] = field(default_factory=dict)
    systems_by_family: dict[str, int] = field(default_factory=dict)
    genomes_by_family: dict[str, int] = field(default_factory=dict)
    proteins_by_accession: dict[str, int] = field(default_factory=dict)
    proteins_by_profile_family: dict[str, int] = field(default_factory=dict)


@dataclass
class ReloadReport:
    database: Path
    dry_run: bool
    hits_path: Path
    hits_sha256: str = ""
    model_family: str = ""
    model_version: str | None = None
    profile_count: int = 0
    params: SelectionParams = field(default_factory=SelectionParams)
    stats: SelectionStats = field(default_factory=SelectionStats)
    sensitivity: dict[str, int] = field(default_factory=dict)
    before: Snapshot = field(default_factory=Snapshot)
    after: Snapshot = field(default_factory=Snapshot)
    previous_model_version: str | None = None
    affected_proteins: int = 0
    unresolved_replicons: int = 0
    checks: dict[str, str] = field(default_factory=dict)
    new_tables: list[str] = field(default_factory=list)
    backup: Path | None = None
    resealed: bool = False
    elapsed_s: float = 0.0

    # --- derived views ------------------------------------------------- #
    @staticmethod
    def _changes(key: str, before: dict[str, int], after: dict[str, int],
                 genomes_before: dict[str, int], genomes_after: dict[str, int]) -> pd.DataFrame:
        keys = sorted(set(before) | set(after))
        frame = pd.DataFrame({
            key: keys,
            "systems_before": [before.get(k, 0) for k in keys],
            "systems_after": [after.get(k, 0) for k in keys],
            "genomes_before": [genomes_before.get(k, 0) for k in keys],
            "genomes_after": [genomes_after.get(k, 0) for k in keys],
        })
        frame["delta"] = frame["systems_after"] - frame["systems_before"]
        order = frame["delta"].abs().sort_values(ascending=False, kind="stable").index
        return frame.loc[order].reset_index(drop=True)

    def type_changes(self) -> pd.DataFrame:
        """Systems per model (system_type). Model renames between releases appear here."""
        b, a = self.before, self.after
        return self._changes("system_type", b.systems_by_type, a.systems_by_type,
                             b.genomes_by_type, a.genomes_by_type)

    def family_changes(self) -> pd.DataFrame:
        """Systems per model family (definitions folder), comparable across releases."""
        b, a = self.before, self.after
        return self._changes("family", b.systems_by_family, a.systems_by_family,
                             b.genomes_by_family, a.genomes_by_family)

    def to_text(self) -> str:
        b, a = self.before, self.after
        families = self.family_changes()
        gained = families[(families.systems_before == 0) & (families.systems_after > 0)]
        lost = families[(families.systems_before > 0) & (families.systems_after == 0)]
        types = self.type_changes()
        new_types = int(((types.systems_before == 0) & (types.systems_after > 0)).sum())
        gone_types = int(((types.systems_before > 0) & (types.systems_after == 0)).sum())
        fam_before, fam_after = b.proteins_by_profile_family, a.proteins_by_profile_family
        fam_gained = sorted(set(fam_after) - set(fam_before))
        fam_lost = sorted(set(fam_before) - set(fam_after))
        s = self.stats
        lines = [
            f"DefenseFinder reload {'DRY RUN' if self.dry_run else 'PUBLISHED'}: {self.database}",
            f"Model set: {self.model_family} {self.model_version or '(unversioned)'} "
            f"({self.profile_count:,} profiles); before-state families from "
            f"{self.previous_model_version or 'the same set'}",
            f"Hits: {self.hits_path} (sha256 {self.hits_sha256[:12]})",
            f"Selection: {json.dumps(asdict(self.params), sort_keys=True)}",
            f"  raw domain rows {s.raw_rows:,}; failed GA/search threshold {s.failed_threshold:,}; "
            f"failed i-Evalue {s.failed_i_evalue:,}; failed coverage {s.failed_coverage:,}; "
            f"unknown profile {s.unknown_profile:,}; profile outside {'/'.join(self.params.subfamilies)} "
            f"models {s.out_of_scope_profile:,}; no alignment {s.missing_alignment:,}",
            f"  selected {s.selected_rows:,}; FAM-variant dedup removed {s.fam_dedup_removed:,}; "
            f"loaded {s.loaded_rows:,}",
            f"Selection sensitivity (raw rows kept): {self.sensitivity}",
            f"defensefinder annotation rows: {b.annotation_rows:,} -> {a.annotation_rows:,} "
            f"(proteins {b.annotation_proteins:,} -> {a.annotation_proteins:,})",
            f"Systems: {b.systems:,} -> {a.systems:,}; genomes with any system "
            f"{b.genomes_with_system:,} -> {a.genomes_with_system:,}; system families "
            f"{len(b.systems_by_family)} -> {len(a.systems_by_family)}; models (system_type) "
            f"{len(b.systems_by_type)} -> {len(a.systems_by_type)} ({new_types} newly called, "
            f"{gone_types} no longer called, including release renames; per-model table in TSV)",
            f"defensefinder_system annotation rows: {b.system_annotation_rows:,} -> "
            f"{a.system_annotation_rows:,}",
            f"System families newly called ({len(gained)}): "
            + ", ".join(f"{r.family} ({r.systems_after} systems, {r.genomes_after} genomes)"
                        for r in gained.itertuples()),
            f"System families no longer called ({len(lost)}): "
            + ", ".join(f"{r.family} ({r.systems_before} systems, {r.genomes_before} genomes)"
                        for r in lost.itertuples()),
            "Largest family changes (systems): "
            + ", ".join(f"{r.family} {r.systems_before}->{r.systems_after}"
                        for r in families.head(12).itertuples()),
            f"Profile families with hits, newly present ({len(fam_gained)}): {', '.join(fam_gained)}",
            f"Profile families with hits, no longer present ({len(fam_lost)}): {', '.join(fam_lost)}",
            f"Proteins with regenerated semantics: {self.affected_proteins:,}",
        ]
        if self.unresolved_replicons:
            lines.append(f"Replicons skipped by the conflict-solver budget: {self.unresolved_replicons}")
        lines.append("Checks:")
        lines.extend(f"  {name}: {result}" for name, result in self.checks.items())
        if self.new_tables:
            lines.append(f"New tables: {self.new_tables}")
        if self.backup:
            lines.append(f"Backup: {self.backup}")
        if not self.dry_run:
            lines.append(f"Resealed: {self.resealed}")
        lines.append(f"Elapsed: {self.elapsed_s:.1f}s")
        return "\n".join(lines)

    def write(self, report_dir: Path, label: str) -> list[Path]:
        """Write the summary, per-type and per-family tables, and a JSON record."""
        report_dir.mkdir(parents=True, exist_ok=True)
        written = []
        summary = report_dir / f"{label}_reload_summary.txt"
        summary.write_text(self.to_text() + "\n")
        written.append(summary)

        types = report_dir / f"{label}_systems_by_type.tsv"
        self.type_changes().to_csv(types, sep="\t", index=False)
        written.append(types)
        families = report_dir / f"{label}_systems_by_family.tsv"
        self.family_changes().to_csv(families, sep="\t", index=False)
        written.append(families)

        b, a = self.before, self.after
        keys = sorted(set(b.proteins_by_profile_family) | set(a.proteins_by_profile_family))
        profile_families = pd.DataFrame({
            "profile_family": keys,
            "proteins_before": [b.proteins_by_profile_family.get(k, 0) for k in keys],
            "proteins_after": [a.proteins_by_profile_family.get(k, 0) for k in keys],
        })
        profiles = report_dir / f"{label}_profile_families.tsv"
        profile_families.to_csv(profiles, sep="\t", index=False)
        written.append(profiles)

        record = report_dir / f"{label}_reload.json"
        payload = {
            "database": str(self.database), "dry_run": self.dry_run,
            "hits_path": str(self.hits_path), "hits_sha256": self.hits_sha256,
            "model_family": self.model_family, "model_version": self.model_version,
            "previous_model_version": self.previous_model_version,
            "profile_count": self.profile_count, "params": asdict(self.params),
            "selection": self.stats.to_dict(), "sensitivity": self.sensitivity,
            "before": asdict(self.before), "after": asdict(self.after),
            "affected_proteins": self.affected_proteins,
            "unresolved_replicons": self.unresolved_replicons,
            "checks": self.checks, "new_tables": self.new_tables,
            "backup": str(self.backup) if self.backup else None,
            "resealed": self.resealed, "elapsed_s": round(self.elapsed_s, 1),
        }
        record.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n")
        written.append(record)
        return written


# --------------------------------------------------------------------------- #
# Snapshots and fingerprints
# --------------------------------------------------------------------------- #
def _tables(conn) -> list[str]:
    return [r[0] for r in conn.execute(
        "SELECT table_name FROM information_schema.tables "
        "WHERE table_schema = 'main' AND table_type = 'BASE TABLE'"
    ).fetchall()]


# Type IV restriction: defense-finder-models 2.x RM_Type_IV / RM_Type_IV_1 (profiles
# RM_Type_IV__*) became McrBC, GmrSD_RM_Type_IV and Other_Type_IV in 3.x.
_TYPE_IV_SUCCESSORS = {"McrBC", "GmrSD_RM_Type_IV", "Other_Type_IV"}


def model_families(model_set: ModelSet | None) -> dict[str, str]:
    """Model name -> family (its definitions folder), comparable across releases.

    RM models are keyed by restriction type: 2.x keeps them all in ``RM/RM/``
    (``RM_Type_I``, ``RM_Type_IIG_2``, ``RM_Type_IV_1``), 3.x splits them into
    ``RM/RM_Type_I/`` .. ``RM/RM_Type_III/`` with chemistry-specific model names
    (``N6_RM_Type_I``) and moves Type IV into McrBC, GmrSD_RM_Type_IV and Other_Type_IV.
    """
    if model_set is None or model_set.root is None:
        return {}
    families: dict[str, str] = {}
    definitions = model_set.root / "definitions"
    for xml in sorted(definitions.rglob("*.xml")):
        relative = xml.relative_to(definitions).parts
        family = xml.parent.name if len(relative) > 2 else xml.stem
        if relative[0] == "RM" and family == "RM":
            family = re.sub(r"_\d+$", "", xml.stem)
        if family in _TYPE_IV_SUCCESSORS:
            family = "RM_Type_IV"
        families[xml.stem] = family
    return families


def align_families(previous: dict[str, str], current: dict[str, str]) -> dict[str, str]:
    """Previous-release model -> family, named as the current release names it.

    A previous family takes the current family most of its shared model names carry
    (2.x All_UG/UG13a -> 3.x DRT-like_UG/UG13a), else a case-insensitive match
    (GasderMIN -> Gasdermin), else keeps its own name.
    """
    votes: dict[str, Counter[str]] = {}
    for model, family in previous.items():
        if model in current:
            votes.setdefault(family, Counter())[current[model]] += 1
    by_lower = {family.lower(): family for family in set(current.values())}
    rename = {}
    for family in set(previous.values()):
        if family in votes:
            rename[family] = votes[family].most_common(1)[0][0]
        else:
            rename[family] = by_lower.get(family.lower(), family)
    return {model: rename[family] for model, family in previous.items()}


def profile_family_of(name_to_stem: dict[str, str]):
    """Profile NAME -> family: the profile file stem's prefix before ``__``."""
    def family(accession: str) -> str:
        return name_to_stem.get(accession, accession).split("__", 1)[0]
    return family


def _snapshot(conn, profile_family, system_family: dict[str, str]) -> Snapshot:
    snap = Snapshot()
    snap.annotation_rows, snap.annotation_proteins = conn.execute(
        "SELECT COUNT(*), COUNT(DISTINCT protein_id) FROM annotations WHERE source = ?", [SOURCE]
    ).fetchone()
    snap.system_annotation_rows = conn.execute(
        "SELECT COUNT(*) FROM annotations WHERE source = ?", [SYSTEM_SOURCE]
    ).fetchone()[0]
    pairs = conn.execute(
        "SELECT DISTINCT protein_id, accession FROM annotations WHERE source = ?", [SOURCE]
    ).fetchdf()
    if not pairs.empty:
        snap.proteins_by_accession = pairs.groupby("accession")["protein_id"].nunique().to_dict()
        pairs["family"] = pairs["accession"].map(profile_family)
        snap.proteins_by_profile_family = (
            pairs.groupby("family")["protein_id"].nunique().to_dict()
        )
    if SYSTEM_TABLE in _tables(conn):
        systems = conn.execute(f"SELECT system_type, genome_id FROM {SYSTEM_TABLE}").fetchdf()
        snap.systems = len(systems)
        snap.genomes_with_system = int(systems["genome_id"].nunique())
        if not systems.empty:
            snap.systems_by_type = systems.groupby("system_type").size().to_dict()
            snap.genomes_by_type = systems.groupby("system_type")["genome_id"].nunique().to_dict()
            systems["family"] = systems["system_type"].map(lambda t: system_family.get(t, t))
            snap.systems_by_family = systems.groupby("family").size().to_dict()
            snap.genomes_by_family = systems.groupby("family")["genome_id"].nunique().to_dict()
    for key in ("systems_by_type", "genomes_by_type", "systems_by_family", "genomes_by_family",
                "proteins_by_accession", "proteins_by_profile_family"):
        setattr(snap, key, {str(k): int(v) for k, v in getattr(snap, key).items()})
    return snap


def _fingerprint(conn, table: str, where: str = "", columns: list[str] | None = None) -> tuple[int, int]:
    """Row count and an order-free hash of a table. With ``columns``, only those columns are hashed, so a
    column added later (e.g. by a schema migration) does not change the fingerprint of unchanged rows."""
    row = ("struct_pack(" + ", ".join(f'"{c}" := t."{c}"' for c in columns) + ")") if columns else "t"
    count, digest = conn.execute(
        f'SELECT COUNT(*), COALESCE(SUM(hash({row})::HUGEINT), 0) FROM "{table}" t {where}'
    ).fetchone()
    return int(count), int(digest)


def _columns(conn, table: str) -> list[str]:
    return [r[0] for r in conn.execute(
        "SELECT column_name FROM information_schema.columns WHERE table_schema = 'main' AND table_name = ? "
        "ORDER BY ordinal_position", [table]).fetchall()]


def _max_generation(conn) -> int:
    if _APPENDED not in _tables(conn):
        return 0
    return int(conn.execute(f"SELECT COALESCE(MAX(generation_id), 0) FROM {_APPENDED}").fetchone()[0])


def _fingerprints(conn, provenance_max: int,
                  provenance_columns: list[str] | None = None) -> dict[str, tuple[int, int]]:
    """Fingerprints of everything the reload must leave unchanged.

    ``provenance_columns``: the predicate_provenance columns present before the reload. Regeneration appends
    a stamp and may add columns to that table (sharur.predicates.provenance.ensure_table); the stamps that
    existed before must keep their values in those original columns."""
    prints = {}
    for table in _tables(conn):
        if table in _REPLACED:
            continue
        if table == _APPENDED:
            prints[table] = _fingerprint(conn, table, f"WHERE generation_id <= {provenance_max}",
                                         columns=provenance_columns)
        elif table == "annotations":
            prints[table] = _fingerprint(
                conn, table, f"WHERE source NOT IN ('{SOURCE}', '{SYSTEM_SOURCE}')"
            )
        elif table == "system_proteins":
            prints[table] = _fingerprint(conn, table, f"WHERE system_source <> '{SYSTEM_SOURCE}'")
        elif table in _PROTEIN_SCOPED:
            prints[table] = _fingerprint(
                conn, table, "WHERE protein_id NOT IN (SELECT protein_id FROM _reload_affected)"
            )
        elif table.startswith("v2_generation_"):
            continue
        else:
            prints[table] = _fingerprint(conn, table)
    return prints


def _set_affected(conn, affected: list[str]) -> None:
    conn.execute("CREATE OR REPLACE TEMP TABLE _reload_affected (protein_id VARCHAR)")
    if affected:
        conn.register("_affected_df", pd.DataFrame({"protein_id": affected}))
        conn.execute("INSERT INTO _reload_affected SELECT protein_id FROM _affected_df")
        conn.unregister("_affected_df")


def _db_state(path: Path) -> tuple[int, int]:
    stat = path.stat()
    return stat.st_size, stat.st_mtime_ns


def selection_sensitivity(gates: pd.DataFrame, hits: pd.DataFrame,
                          params: SelectionParams) -> dict[str, int]:
    """Raw rows each candidate selection rule keeps, for comparison with the applied rule."""
    if gates.empty:
        return {}
    known = gates["known_profile"]
    cov = gates["profile_coverage"]
    i_evalue = pd.to_numeric(hits["i_evalue"], errors="coerce")
    full_evalue = pd.to_numeric(hits["evalue"], errors="coerce")
    ie_pass = i_evalue <= params.i_evalue_sel
    thr = gates["pass_threshold"]
    return {
        "applied": int(gates["selected"].sum()),
        "ga_ie_cov0.5_macsyfinder_default": int((known & thr & ie_pass & (cov >= 0.5)).sum()),
        "ga_ie_cov0.4_defensefinder_default": int((known & thr & ie_pass & (cov >= 0.4)).sum()),
        "ga_only_where_ga_else_ie_cov0.5": int((
            (gates["has_ga"] & gates["pass_ga_only"])
            | (known & ~gates["has_ga"] & ie_pass & (cov >= 0.5))
        ).sum()),
        "ie_cov0.4_without_ga": int((known & ie_pass & (cov >= 0.4)).sum()),
        "legacy_full_evalue_le_1e-15": int((known & (full_evalue <= 1e-15)).sum()),
    }


# --------------------------------------------------------------------------- #
# Reload
# --------------------------------------------------------------------------- #
def reload_defensefinder(
    db_path: str | Path,
    hits_path: str | Path,
    models: str | Path,
    *,
    dry_run: bool = True,
    work_dir: str | Path | None = None,
    params: SelectionParams | None = None,
    threads: int = 2,
    memory_limit: str = "4GB",
    reseal: bool = True,
    report_dir: str | Path | None = None,
    label: str | None = None,
    regenerate_semantics: bool = True,
    previous_models: str | Path | None = None,
    _before_publish=None,
) -> ReloadReport:
    """Rebuild DefenseFinder annotations and systems in a staged copy; publish when asked."""
    started = time.monotonic()
    db = Path(db_path).expanduser().resolve()
    hits_file = Path(hits_path).expanduser().resolve()
    params = params or SelectionParams()
    model_set = load_model_set(models)
    if model_set.models_dir is None:
        raise ReloadValidationError(
            f"{models} has no definitions/ beside its profiles; system calling needs the full "
            "defense-finder-models directory"
        )
    previous_set = load_model_set(previous_models) if previous_models is not None else None
    report = ReloadReport(
        database=db, dry_run=dry_run, hits_path=hits_file, params=params,
        model_family=model_set.family, model_version=model_set.version,
        profile_count=len(model_set.profiles),
        previous_model_version=previous_set.version if previous_set else None,
    )
    # Profile NAMEs -> file stems from both releases; the current one wins on conflict.
    name_to_stem = {p.name: p.stem for p in (previous_set.profiles if previous_set else [])}
    name_to_stem.update({p.name: p.stem for p in model_set.profiles})
    profile_family = profile_family_of(name_to_stem)
    families_after = model_families(model_set)
    families_before = {**families_after,
                       **align_families(model_families(previous_set), families_after)}
    config = {"threads": threads, "memory_limit": memory_limit}

    from sharur.storage.migrations import MIGRATIONS, get_current_version

    wal = db.with_name(db.name + ".wal")
    latest = max(v for v, _, _ in MIGRATIONS)
    if dry_run:
        # A dry run only reads the source: an outstanding WAL would make a file copy stale.
        if wal.exists():
            raise ReloadValidationError(f"{wal} exists; checkpoint {db} before a dry run")
        probe = duckdb.connect(str(db), read_only=True, config=config)
    else:
        # Read-write takes DuckDB's exclusive lock, so this fails while any other
        # process holds the database; CHECKPOINT folds a pending WAL in.
        probe = duckdb.connect(str(db), config=config)
        probe.execute("CHECKPOINT")
    try:
        current = get_current_version(probe)
    finally:
        probe.close()
    if current < latest:
        raise ReloadValidationError(
            f"{db} is at schema version {current}; run `sharur migrate --db {db}` "
            f"(latest {latest}) and re-seal before reloading DefenseFinder"
        )
    if wal.exists():
        raise ReloadValidationError(f"{wal} persists after CHECKPOINT; resolve it before reloading")

    stage_dir = Path(work_dir).expanduser().resolve() if work_dir else db.parent
    stage_dir.mkdir(parents=True, exist_ok=True)
    if not dry_run and os.stat(stage_dir).st_dev != os.stat(db.parent).st_dev:
        raise ReloadValidationError(
            f"--work-dir {stage_dir} is on another filesystem; publishing needs an atomic rename"
        )
    if shutil.disk_usage(stage_dir).free < db.stat().st_size * 1.3:
        raise ReloadValidationError(f"Insufficient space in {stage_dir} for a staged copy of {db}")

    # Read and select hits before copying, so a malformed table costs nothing.
    report.hits_sha256 = file_sha256(hits_file)
    hits = read_hits(hits_file)
    selected, stats, gates = select_hits(hits, model_set, params)
    rows = annotation_rows(selected)
    deduped = dedup_fam_variants(rows).reset_index(drop=True)
    stats.fam_dedup_removed = len(rows) - len(deduped)
    stats.loaded_rows = len(deduped)
    report.stats = stats
    report.sensitivity = selection_sensitivity(gates, hits, params)
    del hits, gates, selected, rows

    source_state = _db_state(db)
    staging = stage_dir / f"{db.name}.defensefinder-reload.staging"
    spill = staging.with_name(staging.name + ".tmp")
    for stale in (staging, staging.with_name(staging.name + ".wal")):
        stale.unlink(missing_ok=True)
    shutil.rmtree(spill, ignore_errors=True)
    logger.info("Copying %s to %s", db, staging)
    shutil.copy2(db, staging)
    if _db_state(db) != source_state:
        staging.unlink(missing_ok=True)
        raise ReloadValidationError(f"{db} changed while it was being copied")

    try:
        def check(name: str, ok: bool, detail: str) -> None:
            report.checks[name] = ("PASS " if ok else "FAIL ") + detail
            if not ok:
                raise ReloadValidationError(f"{name}: {detail}")

        conn = duckdb.connect(str(staging), config=config)
        try:
            conn.register("_df_rows", deduped)
            unknown = conn.execute(
                "SELECT COUNT(*), COUNT(DISTINCT r.protein_id) FROM _df_rows r "
                "ANTI JOIN proteins p USING (protein_id)"
            ).fetchone()
            check("hit_proteins_in_database", unknown[0] == 0,
                  f"{unknown[1]:,} selected-hit protein IDs absent from proteins")
            check("profiles_resolved", stats.unknown_profile == 0,
                  f"{stats.unknown_profile:,} hit rows name profiles absent from the model set "
                  f"{stats.unknown_profile_names[:5]}")

            logger.info("Snapshot before")
            report.before = _snapshot(conn, profile_family, families_before)
            before_members = {r[0] for r in conn.execute(
                "SELECT DISTINCT protein_id FROM annotations WHERE source IN (?, ?)",
                [SOURCE, SYSTEM_SOURCE],
            ).fetchall()}
            tables_before = set(_tables(conn))
            provenance_max = _max_generation(conn)

            logger.info("Replacing %s defensefinder rows", f"{report.before.annotation_rows:,}")
            conn.execute("BEGIN TRANSACTION")
            conn.execute("DELETE FROM annotations WHERE source = ?", [SOURCE])
            next_id = conn.execute(
                "SELECT COALESCE(MAX(annotation_id), 0) FROM annotations"
            ).fetchone()[0]
            conn.execute(
                """
                INSERT INTO annotations (annotation_id, protein_id, source, accession, name,
                                         description, evalue, score, start_aa, end_aa)
                SELECT ? + ROW_NUMBER() OVER (), protein_id, source, accession, name,
                       description, evalue, score, start_aa, end_aa
                FROM _df_rows
                """,
                [int(next_id)],
            )
            record_source_version(conn, model_set, params, stats, hits_file, report.hits_sha256)
            conn.execute("COMMIT")
            conn.unregister("_df_rows")
            after_raw = {r[0] for r in conn.execute(
                "SELECT DISTINCT protein_id FROM annotations WHERE source = ?", [SOURCE]
            ).fetchall()}
        finally:
            conn.close()

        from sharur.colocation import integrate_defense_results, validate_systems

        logger.info("Calling systems")
        systems_df, genes_df = validate_systems(
            db_path=staging, models_dir=model_set.models_dir, model_family=model_set.family,
            source=SOURCE, verbose=False, connect_config=config,
        )
        report.unresolved_replicons = len(systems_df.attrs.get("unresolved_replicons", []))

        # Every protein whose DefenseFinder-derived evidence can change: former and new raw
        # hits, former system members (annotations, member table, summary table, semantic
        # rows) and new members. Fingerprints of everything outside this scope are taken
        # before the caller writes; only the raw rows have changed so far.
        conn = duckdb.connect(str(staging), config=config)
        try:
            scope = set(before_members) | after_raw
            scope.update(str(p) for p in genes_df.get("protein_id", pd.Series(dtype=str)))
            for table, column in (("semantic_atoms", "source_db"), ("semantic_terms", "source_db"),
                                  ("system_proteins", "system_source")):
                if table in tables_before:
                    scope.update(r[0] for r in conn.execute(
                        f"SELECT DISTINCT protein_id FROM {table} WHERE {column} IN (?, ?)",
                        [SOURCE, SYSTEM_SOURCE],
                    ).fetchall())
            if SYSTEM_TABLE in tables_before:
                for (text,) in conn.execute(
                    f"SELECT protein_ids FROM {SYSTEM_TABLE} WHERE protein_ids IS NOT NULL"
                ).fetchall():
                    scope.update(cleaned for value in re.split(r"[,;\s]+", str(text).strip("[](){}"))
                                 if (cleaned := value.strip().strip("'\"")))
            affected = sorted(scope)
            _set_affected(conn, affected)
            logger.info("Fingerprinting rows outside %s affected proteins", f"{len(affected):,}")
            provenance_columns = _columns(conn, _APPENDED) if _APPENDED in _tables(conn) else None
            prints_before = _fingerprints(conn, provenance_max, provenance_columns)
        finally:
            conn.close()

        logger.info("Integrating %s systems", f"{len(systems_df):,}")
        integrated = integrate_defense_results(staging, systems_df, genes_df, connect_config=config)
        outside = set(integrated) - scope
        check("caller_scope_covered", not outside,
              f"{len(outside)} proteins invalidated by the caller outside the reload scope")
        report.affected_proteins = len(affected)

        if regenerate_semantics and affected:
            from sharur.predicates_v2.persistence import generate_and_persist_v2
            from sharur.storage.duckdb_store import DuckDBStore

            logger.info("Regenerating semantics for %s proteins", f"{len(affected):,}")
            store = DuckDBStore(staging, threads=threads, memory_limit=memory_limit)
            try:
                generate_and_persist_v2(store, protein_ids=affected, update_legacy_predicates=True,
                                        return_states=False, predict_topology=False)
            finally:
                store.close()

        conn = duckdb.connect(str(staging), config=config)
        try:
            for table in set(_tables(conn)) - tables_before:
                if table.startswith("v2_generation_") and not conn.execute(
                    f'SELECT COUNT(*) FROM "{table}"'
                ).fetchone()[0]:
                    conn.execute(f'DROP TABLE "{table}"')
            _set_affected(conn, affected)
            logger.info("Validating")
            prints_after = _fingerprints(conn, provenance_max, provenance_columns)
            for table, before in sorted(prints_before.items()):
                scope = ("outside regenerated proteins" if table in _PROTEIN_SCOPED
                         else "stamps before this reload" if table == _APPENDED
                         else "non-DefenseFinder sources" if table == "annotations"
                         else "other caller sources" if table == "system_proteins"
                         else "all rows")
                check(f"unchanged:{table}", prints_after.get(table) == before,
                      f"{before[0]:,} rows, {scope}")
            report.new_tables = sorted(set(_tables(conn)) - tables_before)

            report.after = _snapshot(conn, profile_family, families_after)
            check("annotation_rows_loaded", report.after.annotation_rows == stats.loaded_rows,
                  f"{report.after.annotation_rows:,} defensefinder rows; {stats.loaded_rows:,} selected")
            check("systems_persisted", report.after.systems == len(systems_df),
                  f"{report.after.systems:,} {SYSTEM_TABLE} rows; caller returned {len(systems_df):,}")
            dup = conn.execute("SELECT COUNT(*) - COUNT(DISTINCT annotation_id) FROM annotations").fetchone()[0]
            check("annotation_id_unique", dup == 0, f"{dup} duplicated annotation_id values")
            members = conn.execute(
                "SELECT COUNT(*) FROM system_proteins sp ANTI JOIN annotations a "
                "ON a.protein_id = sp.protein_id AND a.source = ? WHERE sp.system_source = ?",
                [SYSTEM_SOURCE, SYSTEM_SOURCE],
            ).fetchone()[0]
            check("system_members_annotated", members == 0,
                  f"{members} system_proteins rows without a {SYSTEM_SOURCE} annotation")
            if regenerate_semantics and affected:
                stamps = conn.execute(
                    f"SELECT scope FROM {_APPENDED} WHERE generation_id > {provenance_max}"
                ).fetchall()
                check("provenance_stamped", len(stamps) == 1 and stamps[0][0] == "subset",
                      f"{len(stamps)} new stamp(s): {stamps}")
                missing = conn.execute(
                    "SELECT COUNT(*) FROM _reload_affected a "
                    "LEFT JOIN semantic_state s USING (protein_id) WHERE s.protein_id IS NULL"
                ).fetchone()[0]
                check("regenerated_state_present", missing == 0,
                      f"{len(affected):,} regenerated proteins, {missing} missing state")
            conn.execute("CHECKPOINT")
        finally:
            conn.close()

        if report_dir is not None:
            report.elapsed_s = time.monotonic() - started
            report.write(Path(report_dir), label or db.parent.name)
        if dry_run:
            return report

        if _db_state(db) != source_state:
            raise ReloadValidationError(f"{db} changed during the reload; staged copy discarded")
        if _before_publish is not None:
            _before_publish(staging)

        backup = db.with_name(f"{db.name}.pre-defensefinder-reload-{time.strftime('%Y%m%dT%H%M%S')}")
        try:
            os.link(db, backup)
        except OSError:
            shutil.copy2(db, backup)
        os.replace(staging, db)
        report.backup = backup

        seal = db.parent / "dataset.seal.json"
        if reseal and seal.exists():
            from sharur.dataset_seal import seal_dataset

            seal_dataset(db, output_path=seal, force=True)
            report.resealed = True
        if report_dir is not None:
            report.elapsed_s = time.monotonic() - started
            report.write(Path(report_dir), label or db.parent.name)
        return report
    finally:
        for leftover in (staging, staging.with_name(staging.name + ".wal")):
            leftover.unlink(missing_ok=True)
        shutil.rmtree(spill, ignore_errors=True)
        report.elapsed_s = time.monotonic() - started
