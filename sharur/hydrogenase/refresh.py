"""Staged, validated replacement of hydrogenase-derived evidence.

The refresh never writes to the production database file. It copies the
database, reclassifies and regenerates V2 (plus the V1 compatibility cache) for
every protein whose hydrogenase-derived labels may change, validates the copy,
and only then swaps it in with an atomic rename. The prior file is kept as a
hard-linked (or copied) backup beside the database. Any failure before the
rename leaves production untouched.

Withdrawal is provenance-safe because regeneration rebuilds each affected
protein from all of its annotations: a label that another caller (KEGG, Pfam)
supports independently survives, while a label only ``hyddb_subgroup`` emitted
disappears with that source's rows.
"""

from __future__ import annotations

import csv
import os
import shutil
import time
from dataclasses import dataclass, field
from pathlib import Path

import duckdb

from sharur.hydrogenase.classifier import (
    CURATION_FLAG,
    SearchFn,
    classify,
    find_reference,
    run_diamond,
    summarize,
    write_classifications,
)
from sharur.hydrogenase.subgroups import RETIRED_TERMS


# Tables the refresh may change, and how far. Everything else must be identical.
_PROTEIN_SCOPED = ("semantic_atoms", "semantic_state", "semantic_terms", "protein_predicates")
_REPLACED = ("hydrogenase_classifications",)


class RefreshValidationError(RuntimeError):
    """The staged database failed validation; production was left unchanged."""


@dataclass
class Transition:
    protein_id: str
    before: tuple[str, ...]
    after: tuple[str, ...]
    outcome: str
    reference_label: str | None
    curation_reason: str

    @property
    def curation_before(self) -> bool:
        return CURATION_FLAG in self.before

    @property
    def curation_after(self) -> bool:
        return CURATION_FLAG in self.after


@dataclass
class RefreshReport:
    database: Path
    dry_run: bool
    reference_release: str = ""
    reference_sha256: str = ""
    proteins_classified: int = 0
    affected_proteins: int = 0
    outcomes: dict[str, int] = field(default_factory=dict)
    curation: dict[str, int] = field(default_factory=dict)
    transitions: list[Transition] = field(default_factory=list)
    checks: dict[str, str] = field(default_factory=dict)
    new_tables: list[str] = field(default_factory=list)
    backup: Path | None = None
    resealed: bool = False
    elapsed_s: float = 0.0

    @property
    def changed(self) -> list[Transition]:
        return [t for t in self.transitions if t.before != t.after]

    def curation_transitions(self) -> dict[str, int]:
        counts: dict[str, int] = {}
        for t in self.transitions:
            key = f"{'flagged' if t.curation_before else 'unflagged'}->{'flagged' if t.curation_after else 'unflagged'}"
            counts[key] = counts.get(key, 0) + 1
        return counts

    def label_changes(self) -> dict[str, dict[str, int]]:
        removed: dict[str, int] = {}
        added: dict[str, int] = {}
        for t in self.changed:
            for label in set(t.before) - set(t.after):
                removed[label] = removed.get(label, 0) + 1
            for label in set(t.after) - set(t.before):
                added[label] = added.get(label, 0) + 1
        return {"removed": removed, "added": added}

    def write_transitions(self, path: Path) -> None:
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["protein_id", "outcome", "reference_label", "curation_before",
                             "curation_after", "curation_reason", "labels_removed",
                             "labels_added", "labels_before", "labels_after"])
            for t in self.transitions:
                writer.writerow([
                    t.protein_id, t.outcome, t.reference_label or "", int(t.curation_before),
                    int(t.curation_after), t.curation_reason,
                    ";".join(sorted(set(t.before) - set(t.after))),
                    ";".join(sorted(set(t.after) - set(t.before))),
                    ";".join(sorted(t.before)), ";".join(sorted(t.after)),
                ])

    def to_text(self) -> str:
        lines = [
            f"Hydrogenase refresh {'DRY RUN' if self.dry_run else 'PUBLISHED'}: {self.database}",
            f"Reference: HydDB {self.reference_release} (sha256 {self.reference_sha256[:12]})",
            f"Proteins classified: {self.proteins_classified:,}; regenerated: {self.affected_proteins:,}",
            f"Outcomes: {self.outcomes}",
            f"Curation status: {self.curation}",
            f"Curation flag transitions: {self.curation_transitions()}",
            f"Proteins with changed derived labels: {len(self.changed):,}",
        ]
        changes = self.label_changes()
        lines.append(f"Labels removed (proteins): {dict(sorted(changes['removed'].items()))}")
        lines.append(f"Labels added (proteins): {dict(sorted(changes['added'].items()))}")
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


def _tables(conn) -> list[str]:
    return [r[0] for r in conn.execute(
        "SELECT table_name FROM information_schema.tables WHERE table_schema = 'main' AND table_type = 'BASE TABLE'"
    ).fetchall()]


def _fingerprint(conn, table: str, where: str = "") -> tuple[int, int]:
    count, digest = conn.execute(
        f'SELECT COUNT(*), COALESCE(SUM(hash(t)::HUGEINT), 0) FROM "{table}" t {where}'
    ).fetchone()
    return int(count), int(digest)


def _fingerprints(conn) -> dict[str, tuple[int, int]]:
    prints = {}
    for table in _tables(conn):
        if table in _REPLACED:
            continue
        if table == "annotations":
            prints[table] = _fingerprint(conn, table, "WHERE source <> 'hyddb_subgroup'")
            prints["annotations[hyddb]"] = _fingerprint(conn, table, "WHERE LOWER(source) = 'hyddb'")
        elif table in _PROTEIN_SCOPED:
            prints[table] = _fingerprint(
                conn, table, "WHERE protein_id NOT IN (SELECT protein_id FROM _refresh_affected)"
            )
        else:
            prints[table] = _fingerprint(conn, table)
    return prints


def _derived_labels(conn) -> dict[str, tuple[str, ...]]:
    labels: dict[str, set[str]] = {}
    for pid, acc in conn.execute(
        "SELECT protein_id, accession FROM annotations WHERE source = 'hyddb_subgroup'"
    ).fetchall():
        labels.setdefault(pid, set()).add(acc)
    return {pid: tuple(sorted(v)) for pid, v in labels.items()}


def _db_state(path: Path) -> tuple[int, int]:
    stat = path.stat()
    return stat.st_size, stat.st_mtime_ns


def refresh_hydrogenases(
    db_path: str | Path,
    *,
    threads: int = 4,
    dry_run: bool = True,
    reference_dir: Path | None = None,
    search: SearchFn = run_diamond,
    transitions_path: Path | None = None,
    reseal: bool = True,
    _before_publish=None,
) -> RefreshReport:
    """Reclassify hydrogenases and replace derived evidence via a validated staged copy."""
    started = time.monotonic()
    db = Path(db_path).expanduser().resolve()
    report = RefreshReport(database=db, dry_run=dry_run)
    reference = find_reference(reference_dir)
    report.reference_release, report.reference_sha256 = reference.release, reference.sha256

    wal = db.with_name(db.name + ".wal")
    # Opening read-write takes DuckDB's exclusive lock, so this fails while any
    # other process holds the database; CHECKPOINT folds a pending WAL in.
    from sharur.storage.migrations import MIGRATIONS, get_current_version

    probe = duckdb.connect(str(db))
    try:
        probe.execute("CHECKPOINT")
        current = get_current_version(probe)
    finally:
        probe.close()
    latest = max(v for v, _, _ in MIGRATIONS)
    if current < latest:
        raise RefreshValidationError(
            f"{db} is at schema version {current}; run `sharur migrate --db {db}` "
            f"(latest {latest}) and re-seal before refreshing hydrogenases"
        )
    if wal.exists():
        raise RefreshValidationError(f"{wal} persists after CHECKPOINT; resolve it before refreshing")

    free = shutil.disk_usage(db.parent).free
    if free < db.stat().st_size * 1.2:
        raise RefreshValidationError(f"Insufficient space beside {db} for a staged copy")

    source_state = _db_state(db)
    staging = db.with_name(f"{db.name}.hydrogenase-refresh.staging")
    for stale in (staging, staging.with_name(staging.name + ".wal")):
        stale.unlink(missing_ok=True)
    shutil.copy2(db, staging)

    try:
        conn = duckdb.connect(str(staging))
        try:
            before_labels = _derived_labels(conn)
            results = classify(conn, reference, threads=threads, search=search)
            after_pids = {item.protein_id for item in results if item.derived_labels}
            affected = sorted(set(before_labels) | after_pids)
            conn.execute("CREATE TEMP TABLE _refresh_affected (protein_id VARCHAR)")
            if affected:
                conn.executemany("INSERT INTO _refresh_affected VALUES (?)", [(p,) for p in affected])
            tables_before = set(_tables(conn))
            prints_before = _fingerprints(conn)
            totals_before = {
                table: conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
                for table in ("semantic_state", "protein_predicates") if table in tables_before
            }
            write_classifications(conn, results, reference)
        finally:
            conn.close()

        from sharur.predicates_v2.persistence import generate_and_persist_v2
        from sharur.storage.duckdb_store import DuckDBStore

        if affected:
            store = DuckDBStore(staging)
            try:
                generate_and_persist_v2(store, protein_ids=affected, update_legacy_predicates=True,
                                        return_states=False, predict_topology=False)
            finally:
                store.close()

        conn = duckdb.connect(str(staging))
        try:
            # V2 subset generation creates empty scratch tables; leave the schema as found.
            for table in set(_tables(conn)) - tables_before:
                if table.startswith("v2_generation_") and not conn.execute(
                    f'SELECT COUNT(*) FROM "{table}"'
                ).fetchone()[0]:
                    conn.execute(f'DROP TABLE "{table}"')
            conn.execute("CREATE TEMP TABLE _refresh_affected (protein_id VARCHAR)")
            if affected:
                conn.executemany("INSERT INTO _refresh_affected VALUES (?)", [(p,) for p in affected])
            _validate(conn, report, results, affected, prints_before, tables_before, totals_before)
            after_labels = _derived_labels(conn)
            conn.execute("CHECKPOINT")
        finally:
            conn.close()

        by_pid = {item.protein_id: item for item in results}
        for pid in affected:
            item = by_pid.get(pid)
            report.transitions.append(Transition(
                protein_id=pid,
                before=before_labels.get(pid, ()),
                after=after_labels.get(pid, ()),
                outcome=item.outcome if item else "no_hyddb_hit",
                reference_label=item.hit.label if item and item.hit else None,
                curation_reason=item.curation_reason if item else "",
            ))
        summary = summarize(results)
        report.outcomes, report.curation = summary.outcomes, summary.curation
        report.proteins_classified, report.affected_proteins = len(results), len(affected)
        if transitions_path:
            report.write_transitions(Path(transitions_path))

        if dry_run:
            return report

        if _db_state(db) != source_state:
            raise RefreshValidationError(f"{db} changed during the refresh; staged copy discarded")
        if _before_publish is not None:
            _before_publish(staging)

        backup = db.with_name(f"{db.name}.pre-hydrogenase-refresh-{time.strftime('%Y%m%dT%H%M%S')}")
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
        return report
    finally:
        for leftover in (staging, staging.with_name(staging.name + ".wal")):
            leftover.unlink(missing_ok=True)
        report.elapsed_s = time.monotonic() - started


def _validate(conn, report, results, affected, prints_before, tables_before, totals_before) -> None:
    def check(name: str, ok: bool, detail: str) -> None:
        report.checks[name] = ("PASS " if ok else "FAIL ") + detail
        if not ok:
            raise RefreshValidationError(f"{name}: {detail}")

    prints_after = _fingerprints(conn)
    for table, before in sorted(prints_before.items()):
        after = prints_after.get(table)
        scope = ("outside regenerated proteins" if table in _PROTEIN_SCOPED
                 else "raw HydDB HMM rows" if table == "annotations[hyddb]"
                 else "excluding hyddb_subgroup rows" if table == "annotations" else "all rows")
        check(f"unchanged:{table}", after == before, f"{before[0]:,} rows, {scope}")

    report.new_tables = sorted(set(_tables(conn)) - tables_before - set(_REPLACED))

    hyd_proteins = conn.execute(
        "SELECT COUNT(DISTINCT protein_id) FROM annotations WHERE LOWER(source) = 'hyddb'"
    ).fetchone()[0]
    rows, distinct = conn.execute(
        "SELECT COUNT(*), COUNT(DISTINCT protein_id) FROM hydrogenase_classifications"
    ).fetchone()
    check("one_classification_per_protein", rows == distinct == hyd_proteins == len(results),
          f"{rows:,} rows / {distinct:,} proteins / {hyd_proteins:,} HydDB-hit proteins")

    for table, before in totals_before.items():
        n = conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
        check(f"row_count:{table}", n == before, f"{n:,} rows after, {before:,} before")

    retired = ",".join(f"'{t}'" for t in sorted(RETIRED_TERMS))
    stale = conn.execute(f"""
        SELECT (SELECT COUNT(*) FROM annotations WHERE source = 'hyddb_subgroup' AND accession IN ({retired}))
             + (SELECT COUNT(*) FROM semantic_atoms WHERE source_db = 'hyddb_subgroup' AND atom_id IN ({retired}))
    """).fetchone()[0]
    check("retired_labels_withdrawn", stale == 0, f"{stale} hydrogenase-derived retired labels remain")

    missing = conn.execute("""
        SELECT COUNT(*) FROM _refresh_affected a
        LEFT JOIN semantic_state s USING (protein_id) WHERE s.protein_id IS NULL
    """).fetchone()[0]
    check("regenerated_state_present", missing == 0, f"{len(affected):,} regenerated proteins, {missing} missing state")
