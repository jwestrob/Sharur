"""Protein families in DuckDB: export, staging build, import, summaries, activation, assignment.

Several clustering runs live side by side; each run is one partition of the
proteins. A Sharur build clusters twice (subfamilies over all proteins, families
over the subfamily representatives) and so records two runs, the subfamily run
naming the family run as its ``parent_run``. An external MMseqs2 clustering can
be registered as a further run (:func:`import_clustering`).

Tables (created by this module, never by Stage 07):

``protein_family_runs``
    One row per clustering run: ``clustering_run``, ``level`` (``family``,
    ``subfamily`` or the label of an imported run), ``parent_run``, ``build_id``
    (shared by the two levels of one Sharur build), ``origin`` (``sharur`` or
    ``import``), ``is_default`` (the run behind the per-level views), tool and
    version, ``parameters``, the clustered proteins' count and sequence digest,
    and the family count.
``protein_families``
    One row per protein and run: ``protein_id``, ``family_id``,
    ``is_representative``, ``clustering_run``, ``bin_id``, ``assignment``
    (``cluster``, ``import``, or for proteins placed later ``search``,
    ``new_subfamily``, ``new_family``) and ``batch_id``.
``family_stats``
    One row per family and run: ``n_members``, ``n_genomes``,
    ``rep_protein_id``, ``rep_length``, ``parent_family_id`` (subfamily runs),
    ``n_subfamilies`` (family runs with a subfamily run), ``n_orders`` (GTDB
    orders in ``bins.taxonomy``), ``median_length``, ``pfam_fraction`` (members
    with any Pfam hit), ``dominant_pfam_architecture`` and its share of members,
    ``ko_fraction``, ``dominant_ko`` and its share, ``pfam_dark`` (no member has
    a Pfam hit).
``protein_family_assignments``
    One row per incremental assignment batch.

Views: ``genome_family_presence`` (``clustering_run, genome_id, feature_id,
n_proteins``) for every run; ``genome_protein_families`` and
``genome_protein_subfamilies`` (``genome_id, feature_id, n_proteins``) for the
default family and subfamily runs; ``protein_family_levels`` (``protein_id,
bin_id, family_id, subfamily_id``) for the default build.

Identifiers. A Sharur build's ``build_id`` is ``pfr-`` plus the first ten hex
digits of a SHA-256 over both levels' criteria, the MMseqs2 version, the input
sequence digest and the resulting partition; its runs are ``<build_id>-family``
and ``<build_id>-subfamily``. Family ids are ``fam:<tag>:<representative>`` and
subfamily ids ``sub:<tag>:<representative>``, ``<tag>`` being those ten hex
digits, so rerunning on the same proteins with the same parameters and getting
the same partition reproduces every id. An imported run keeps the name given at
import, and its family ids are ``fam:<tag>:<representative>`` with ``<tag>``
hashed from the name, parameters, digest and partition.

The sequence digest is the SHA-256 over ``protein_id<TAB>sequence<LF>`` lines in
``protein_id`` order, with trailing ``*`` removed from sequences.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import re
import shutil
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Iterable, Sequence

import duckdb

from sharur.protein_families import mmseqs
from sharur.protein_families.mmseqs import LevelParams

RUNS = "protein_family_runs"
MEMBERS = "protein_families"
STATS = "family_stats"
BATCHES = "protein_family_assignments"
RESCUES = "protein_family_rescues"
PRESENCE = "genome_family_presence"
GENOME_FAMILIES = "genome_protein_families"
GENOME_SUBFAMILIES = "genome_protein_subfamilies"
LEVELS_VIEW = "protein_family_levels"
GENOME_FAMILIES_CLUSTERED = "genome_protein_families_clustered"
GENOME_SUBFAMILIES_CLUSTERED = "genome_protein_subfamilies_clustered"
TABLES = (RUNS, MEMBERS, STATS, BATCHES, RESCUES)
VIEWS = (PRESENCE, GENOME_FAMILIES, GENOME_SUBFAMILIES, LEVELS_VIEW, GENOME_FAMILIES_CLUSTERED,
         GENOME_SUBFAMILIES_CLUSTERED)
SCHEMA_VERSION = 3
# A build whose share of assigned proteins exceeds this is due for a full re-cluster.
RECLUSTER_ASSIGNED_FRACTION = 0.10
# Catalog prefix of the dataset database attached as the source of proteins (staging-file work).
SRC = "__pf_src."
_RUN_NAME = re.compile(r"[A-Za-z0-9][\w.-]*")

DDL = {
    RUNS: f"""CREATE TABLE IF NOT EXISTS {RUNS} (
    clustering_run VARCHAR PRIMARY KEY,
    level VARCHAR NOT NULL,
    parent_run VARCHAR,
    build_id VARCHAR,
    origin VARCHAR NOT NULL,
    is_default BOOLEAN NOT NULL,
    created_at TIMESTAMP NOT NULL,
    tool VARCHAR,
    tool_version VARCHAR,
    parameters JSON,
    input_proteins BIGINT NOT NULL,
    input_digest VARCHAR NOT NULL,
    families BIGINT,
    source_db VARCHAR,
    schema_version INTEGER NOT NULL,
    details JSON
)""",
    MEMBERS: f"""CREATE TABLE IF NOT EXISTS {MEMBERS} (
    protein_id VARCHAR NOT NULL,
    family_id VARCHAR NOT NULL,
    is_representative BOOLEAN NOT NULL,
    clustering_run VARCHAR NOT NULL,
    bin_id VARCHAR NOT NULL,
    assignment VARCHAR NOT NULL,
    batch_id VARCHAR,
    PRIMARY KEY (clustering_run, protein_id)
)""",
    STATS: f"""CREATE TABLE IF NOT EXISTS {STATS} (
    family_id VARCHAR NOT NULL,
    clustering_run VARCHAR NOT NULL,
    n_members BIGINT NOT NULL,
    n_genomes BIGINT NOT NULL,
    rep_protein_id VARCHAR NOT NULL,
    rep_length INTEGER,
    parent_family_id VARCHAR,
    n_subfamilies BIGINT,
    n_orders BIGINT NOT NULL,
    median_length DOUBLE,
    pfam_fraction DOUBLE NOT NULL,
    dominant_pfam_architecture VARCHAR,
    dominant_pfam_architecture_fraction DOUBLE,
    ko_fraction DOUBLE NOT NULL,
    dominant_ko VARCHAR,
    dominant_ko_fraction DOUBLE,
    pfam_dark BOOLEAN NOT NULL,
    n_rescued BIGINT,
    PRIMARY KEY (clustering_run, family_id)
)""",
    BATCHES: f"""CREATE TABLE IF NOT EXISTS {BATCHES} (
    batch_id VARCHAR PRIMARY KEY,
    build_id VARCHAR NOT NULL,
    family_run VARCHAR NOT NULL,
    subfamily_run VARCHAR NOT NULL,
    created_at TIMESTAMP NOT NULL,
    tool_version VARCHAR,
    parameters JSON,
    input_proteins BIGINT NOT NULL,
    input_digest VARCHAR NOT NULL,
    details JSON
)""",
    RESCUES: f"""CREATE TABLE IF NOT EXISTS {RESCUES} (
    rescue_id VARCHAR NOT NULL,
    family_run VARCHAR NOT NULL,
    subfamily_run VARCHAR,
    base_family_run VARCHAR NOT NULL,
    protein_id VARCHAR NOT NULL,
    bin_id VARCHAR NOT NULL,
    partial VARCHAR,
    original_family_id VARCHAR NOT NULL,
    original_subfamily_id VARCHAR,
    family_id VARCHAR NOT NULL,
    subfamily_id VARCHAR,
    joined_subfamily BOOLEAN NOT NULL,
    hit_protein_id VARCHAR NOT NULL,
    identity DOUBLE NOT NULL,
    query_coverage DOUBLE NOT NULL,
    target_coverage DOUBLE,
    bits DOUBLE,
    evalue DOUBLE,
    query_length INTEGER,
    target_length INTEGER,
    PRIMARY KEY (family_run, protein_id)
)""",
}
# Columns added after a table first shipped; `ensure_schema` adds them to older files.
ADDED_COLUMNS = {STATS: [("n_rescued", "BIGINT")]}
SCHEMA = ";\n".join(DDL.values()) + ";"

# Secondary ART indexes on protein_families cost about 110 bytes per row each on disk (over 1.6 GB for two
# runs of 3.6 million proteins) and buy little: DuckDB scans the columns in a fraction of a second. Files
# written before 2026-11-12 carry two; `_finalize` drops them.
LEGACY_INDEXES = ("idx_protein_families_protein", "idx_protein_families_family")

_BASE_RUN = "CASE WHEN r.origin = 'rescue' THEN json_extract_string(r.details, '$.base_run') ELSE r.clustering_run END"

VIEW_SQL = f"""
CREATE OR REPLACE VIEW {PRESENCE} AS
SELECT clustering_run, bin_id AS genome_id, family_id AS feature_id, COUNT(*) AS n_proteins,
       COUNT(*) FILTER (WHERE assignment = 'rescue') AS n_rescued
FROM {MEMBERS} GROUP BY clustering_run, bin_id, family_id;
CREATE OR REPLACE VIEW {GENOME_FAMILIES} AS
SELECT f.bin_id AS genome_id, f.family_id AS feature_id, COUNT(*) AS n_proteins,
       COUNT(*) FILTER (WHERE f.assignment = 'rescue') AS n_rescued
FROM {MEMBERS} f JOIN {RUNS} r USING (clustering_run)
WHERE r.level = 'family' AND r.is_default GROUP BY f.bin_id, f.family_id;
CREATE OR REPLACE VIEW {GENOME_SUBFAMILIES} AS
SELECT f.bin_id AS genome_id, f.family_id AS feature_id, COUNT(*) AS n_proteins,
       COUNT(*) FILTER (WHERE f.assignment = 'rescue') AS n_rescued
FROM {MEMBERS} f JOIN {RUNS} r USING (clustering_run)
WHERE r.level = 'subfamily' AND r.is_default GROUP BY f.bin_id, f.family_id;
CREATE OR REPLACE VIEW {GENOME_FAMILIES_CLUSTERED} AS
SELECT f.bin_id AS genome_id, f.family_id AS feature_id, COUNT(*) AS n_proteins
FROM {MEMBERS} f
JOIN (SELECT {_BASE_RUN} AS clustering_run FROM {RUNS} r WHERE r.level = 'family' AND r.is_default) b
  USING (clustering_run)
GROUP BY f.bin_id, f.family_id;
CREATE OR REPLACE VIEW {GENOME_SUBFAMILIES_CLUSTERED} AS
SELECT f.bin_id AS genome_id, f.family_id AS feature_id, COUNT(*) AS n_proteins
FROM {MEMBERS} f
JOIN (SELECT {_BASE_RUN} AS clustering_run FROM {RUNS} r WHERE r.level = 'subfamily' AND r.is_default) b
  USING (clustering_run)
GROUP BY f.bin_id, f.family_id;
CREATE OR REPLACE VIEW {LEVELS_VIEW} AS
SELECT f.protein_id, f.bin_id, f.family_id, s.family_id AS subfamily_id,
       f.clustering_run AS family_run, s.clustering_run AS subfamily_run, f.assignment = 'rescue' AS rescued
FROM {MEMBERS} f
JOIN {RUNS} rf ON rf.clustering_run = f.clustering_run AND rf.level = 'family' AND rf.is_default
LEFT JOIN {RUNS} rs ON rs.parent_run = rf.clustering_run
LEFT JOIN {MEMBERS} s ON s.clustering_run = rs.clustering_run AND s.protein_id = f.protein_id;
"""


class FamilyError(RuntimeError):
    """A protein-family precondition failed; nothing was written."""


# --------------------------------------------------------------------------- #
# Small helpers
# --------------------------------------------------------------------------- #


def _q(path: Path | str) -> str:
    return "'" + str(path).replace("'", "''") + "'"


def clean_sequence(sequence: str | None) -> str:
    return (sequence or "").rstrip("*")


def build_tag(build_id: str) -> str:
    return build_id.split("-", 1)[1]


def family_id(tag: str, representative: str) -> str:
    return f"fam:{tag}:{representative}"


def subfamily_id(tag: str, representative: str) -> str:
    return f"sub:{tag}:{representative}"


def representative_of(feature_id: str) -> str:
    """The representative protein id inside a family or subfamily id."""
    return feature_id.split(":", 2)[2]


def has_families(con: duckdb.DuckDBPyConnection, catalog: str | None = None) -> bool:
    rows = con.execute(
        "SELECT table_name FROM information_schema.tables WHERE table_schema = 'main' "
        + ("AND table_catalog = ? " if catalog else "AND table_catalog = current_database() ")
        + "AND table_name IN (?, ?, ?)",
        [*([catalog] if catalog else []), MEMBERS, RUNS, STATS]).fetchall()
    return len(rows) == 3


def connect(path: Path | str, *, read_only: bool = False, threads: int = 4,
            memory_limit: str = "8GB") -> duckdb.DuckDBPyConnection:
    con = duckdb.connect(str(path), read_only=read_only)
    con.execute(f"SET threads = {int(threads)}")
    con.execute(f"SET memory_limit = '{memory_limit}'")
    con.execute("SET preserve_insertion_order = false")
    return con


def _sha(payload: Any) -> str:
    return hashlib.sha256(json.dumps(payload, sort_keys=True, default=str).encode()).hexdigest()


def _now() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


# --------------------------------------------------------------------------- #
# Export and digests
# --------------------------------------------------------------------------- #


def _protein_rows(con, prefix: str, ids_table: str | None, batch: int = 100_000):
    join = f"JOIN {ids_table} ids USING (protein_id)" if ids_table else ""
    cur = con.execute(f"SELECT p.protein_id, p.sequence FROM {prefix}proteins p {join} ORDER BY p.protein_id")
    while True:
        rows = cur.fetchmany(batch)
        if not rows:
            return
        yield from rows


def export_proteins(con, fasta: Path | None, *, prefix: str = "", ids_table: str | None = None) -> dict[str, Any]:
    """Write proteins (optionally those in ``ids_table``) as FASTA in ``protein_id`` order.

    Returns the protein count, the sequence digest and the count of empty records
    skipped. With ``fasta=None`` only the digest is computed. Sequences are never
    printed.
    """
    digest = hashlib.sha256()
    written = empty = 0
    handle = open(fasta, "w") if fasta else None
    try:
        for protein_id, sequence in _protein_rows(con, prefix, ids_table):
            if any(ch.isspace() for ch in protein_id):
                raise FamilyError(f"protein id {protein_id!r} contains whitespace; MMseqs2 would truncate it")
            seq = clean_sequence(sequence)
            if not seq:
                empty += 1
                continue
            digest.update(f"{protein_id}\t{seq}\n".encode())
            if handle:
                handle.write(f">{protein_id}\n{seq}\n")
            written += 1
    finally:
        if handle:
            handle.close()
    return {"proteins": written, "skipped_empty": empty, "digest": digest.hexdigest()}


def _digest_of(con, sql: str, prefix: str) -> dict[str, Any]:
    """Sequence digest of the proteins that ``sql`` (one protein_id column) selects."""
    con.execute(f"CREATE OR REPLACE TEMP TABLE _pf_ids AS {sql}")
    try:
        return export_proteins(con, None, prefix=prefix, ids_table="_pf_ids")
    finally:
        con.execute("DROP TABLE IF EXISTS _pf_ids")


# --------------------------------------------------------------------------- #
# Summaries
# --------------------------------------------------------------------------- #


def _architectures(con, prefix: str) -> list[tuple[str, str]]:
    """``(protein_id, architecture)`` for proteins in ``_pf_scope`` with Pfam hits (``sharur architecture`` rules)."""
    from sharur.architecture import Domain, compact, resolve  # noqa: PLC0415

    cur = con.execute(
        f"""SELECT a.protein_id, COALESCE(NULLIF(a.name, ''), a.accession), a.accession, a.start_aa,
                   a.end_aa, a.evalue
            FROM {prefix}annotations a JOIN _pf_scope USING (protein_id)
            WHERE lower(a.source) = 'pfam' ORDER BY a.protein_id""")
    out: list[tuple[str, str]] = []
    current: str | None = None
    hits: list[Domain] = []

    def flush() -> None:
        if current is None:
            return
        names = [d.name for d in resolve(hits)] or sorted({h.name for h in hits})
        out.append((current, compact(names)))

    while True:
        rows = cur.fetchmany(200_000)
        if not rows:
            break
        for pid, name, accession, start, end, evalue in rows:
            if pid != current:
                flush()
                current, hits = pid, []
            hits.append(Domain(name, accession, start, end, evalue, "pfam"))
    flush()
    return out


def build_stats(con, runs: Sequence[str], *, prefix: str = "") -> dict[str, int]:
    """(Re)build ``family_stats`` rows of ``runs`` from members and the dataset tables under ``prefix``."""
    import pandas as pd  # noqa: PLC0415

    runs = list(runs)
    if not runs:
        return {}
    con.execute(DDL[STATS])
    ensure_schema(con)
    con.execute("CREATE OR REPLACE TEMP TABLE _pf_runs AS SELECT unnest(?::VARCHAR[]) AS clustering_run", [runs])
    con.execute(f"""CREATE OR REPLACE TEMP TABLE _pf_scope AS
                    SELECT DISTINCT protein_id FROM {MEMBERS} JOIN _pf_runs USING (clustering_run)""")
    con.execute(f"""CREATE OR REPLACE TEMP TABLE _pf_prot AS
                    SELECT p.protein_id, p.sequence_length,
                           NULLIF(regexp_extract(b.taxonomy, 'o__([^;]+)', 1), '') AS order_name
                    FROM {prefix}proteins p JOIN _pf_scope USING (protein_id)
                    LEFT JOIN {prefix}bins b ON b.bin_id = p.bin_id""")
    con.register("_pf_arch_df", pd.DataFrame(_architectures(con, prefix), columns=["protein_id", "architecture"]))
    con.execute("CREATE OR REPLACE TEMP TABLE _pf_arch AS SELECT * FROM _pf_arch_df")
    con.unregister("_pf_arch_df")
    con.execute(f"""CREATE OR REPLACE TEMP TABLE _pf_ko AS
                    SELECT DISTINCT a.protein_id, a.accession AS ko
                    FROM {prefix}annotations a JOIN _pf_scope USING (protein_id)
                    WHERE lower(a.source) IN ('kofam', 'kegg') AND regexp_matches(a.accession, '^K[0-9]{{5}}$')""")
    info = {r[0]: (r[1], r[2]) for r in con.execute(
        f"SELECT clustering_run, parent_run, (SELECT min(c.clustering_run) FROM {RUNS} c "
        f"WHERE c.parent_run = r.clustering_run) FROM {RUNS} r").fetchall()}
    counts = {}
    for run in runs:
        parent, child = info.get(run, (None, None))
        con.execute(f"DELETE FROM {STATS} WHERE clustering_run = ?", [run])
        con.execute(f"""CREATE OR REPLACE TEMP TABLE _pf_facts AS
            SELECT f.protein_id, f.family_id, f.bin_id, f.is_representative, f.assignment, p.sequence_length,
                   p.order_name
            FROM {MEMBERS} f JOIN _pf_prot p USING (protein_id) WHERE f.clustering_run = ?""", [run])
        child_sql = (f"""(SELECT f.family_id, COUNT(DISTINCT c.family_id) AS n FROM _pf_facts f
                          JOIN {MEMBERS} c ON c.protein_id = f.protein_id AND c.clustering_run = '{child}'
                          GROUP BY 1)""" if child else "(SELECT NULL::VARCHAR AS family_id, NULL::BIGINT AS n WHERE false)")
        parent_sql = (f"""(SELECT f.family_id, p.family_id AS parent_family_id FROM _pf_facts f
                           JOIN {MEMBERS} p ON p.protein_id = f.protein_id AND p.clustering_run = '{parent}'
                           WHERE f.is_representative)""" if parent
                      else "(SELECT NULL::VARCHAR AS family_id, NULL::VARCHAR AS parent_family_id WHERE false)")
        con.execute(f"""
            INSERT INTO {STATS} (family_id, clustering_run, n_members, n_genomes, rep_protein_id, rep_length,
                                 parent_family_id, n_subfamilies, n_orders, median_length, pfam_fraction,
                                 dominant_pfam_architecture, dominant_pfam_architecture_fraction, ko_fraction,
                                 dominant_ko, dominant_ko_fraction, pfam_dark, n_rescued)
            WITH base AS (
                SELECT f.family_id, COUNT(*) AS n, COUNT(DISTINCT f.bin_id) AS n_genomes,
                       COUNT(DISTINCT f.order_name) AS n_orders, median(f.sequence_length) AS median_length,
                       COUNT(a.protein_id) AS n_pfam,
                       COUNT(*) FILTER (WHERE f.assignment = 'rescue') AS n_rescued,
                       min(f.protein_id) FILTER (WHERE f.is_representative) AS rep,
                       min(f.sequence_length) FILTER (WHERE f.is_representative) AS rep_length
                FROM _pf_facts f LEFT JOIN _pf_arch a USING (protein_id) GROUP BY f.family_id),
            arch AS (
                SELECT family_id, architecture, cnt FROM (
                    SELECT f.family_id, a.architecture, COUNT(*) AS cnt,
                           row_number() OVER (PARTITION BY f.family_id ORDER BY COUNT(*) DESC, a.architecture) AS r
                    FROM _pf_facts f JOIN _pf_arch a USING (protein_id) GROUP BY f.family_id, a.architecture)
                WHERE r = 1),
            ko_any AS (
                SELECT f.family_id, COUNT(DISTINCT k.protein_id) AS n_ko
                FROM _pf_facts f JOIN _pf_ko k USING (protein_id) GROUP BY f.family_id),
            ko AS (
                SELECT family_id, ko, cnt FROM (
                    SELECT f.family_id, k.ko, COUNT(DISTINCT k.protein_id) AS cnt,
                           row_number() OVER (PARTITION BY f.family_id
                                              ORDER BY COUNT(DISTINCT k.protein_id) DESC, k.ko) AS r
                    FROM _pf_facts f JOIN _pf_ko k USING (protein_id) GROUP BY f.family_id, k.ko)
                WHERE r = 1),
            child AS {child_sql},
            parent AS {parent_sql}
            SELECT b.family_id, ?, b.n, b.n_genomes, b.rep, b.rep_length, parent.parent_family_id, child.n,
                   b.n_orders, b.median_length, b.n_pfam / b.n, arch.architecture, arch.cnt / b.n,
                   COALESCE(ko_any.n_ko, 0) / b.n, ko.ko, ko.cnt / b.n, b.n_pfam = 0, b.n_rescued
            FROM base b LEFT JOIN arch USING (family_id) LEFT JOIN ko_any USING (family_id)
            LEFT JOIN ko USING (family_id) LEFT JOIN child USING (family_id) LEFT JOIN parent USING (family_id)""",
                    [run])
        counts[run] = con.execute(f"SELECT COUNT(*) FROM {STATS} WHERE clustering_run = ?", [run]).fetchone()[0]
        missing_rep = con.execute(f"SELECT COUNT(*) FROM {STATS} WHERE clustering_run = ? AND rep_protein_id IS NULL",
                                  [run]).fetchone()[0]
        if missing_rep:
            raise FamilyError(f"run {run}: {missing_rep:,} families have no representative row")
        con.execute(f"UPDATE {RUNS} SET families = ? WHERE clustering_run = ?", [counts[run], run])
    for temp in ("_pf_runs", "_pf_scope", "_pf_prot", "_pf_arch", "_pf_ko", "_pf_facts"):
        con.execute(f"DROP TABLE IF EXISTS {temp}")
    return counts


def ensure_schema(con, catalog: str | None = None) -> None:
    """Add columns introduced after a table first shipped (files written by older versions)."""
    prefix = f"{catalog}." if catalog else ""
    for table, columns in ADDED_COLUMNS.items():
        if _columns(con, table, catalog):
            for name, kind in columns:
                con.execute(f"ALTER TABLE {prefix}{table} ADD COLUMN IF NOT EXISTS {name} {kind}")


def _columns(con, table: str, catalog: str | None = None) -> list[str]:
    return [r[0] for r in con.execute(
        "SELECT column_name FROM information_schema.columns WHERE table_schema = 'main' AND table_name = ? "
        + ("AND table_catalog = ? " if catalog else "AND table_catalog = current_database() ")
        + "ORDER BY ordinal_position", [table, *([catalog] if catalog else [])]).fetchall()]


def _copy_rows(con, table: str, source: str, where: str) -> None:
    """INSERT the rows of ``source``.``table`` matching ``where`` by column name (older files lack new columns)."""
    target = _columns(con, table)
    available = set(_columns(con, table, source))
    cols = ", ".join(c for c in target if c in available)
    con.execute(f"INSERT INTO {table} ({cols}) SELECT {cols} FROM {source}.{table} t WHERE {where}")


def _finalize(con) -> None:
    for view in VIEWS:
        con.execute(f"DROP VIEW IF EXISTS {view}")
    for index in LEGACY_INDEXES:
        con.execute(f"DROP INDEX IF EXISTS {index}")
    con.execute(VIEW_SQL)


def _table_fingerprint(con, table: str, catalog: str | None = None) -> tuple[int, int]:
    cols = _columns(con, table, catalog)
    packed = "struct_pack(" + ", ".join(f'"{c}" := t."{c}"' for c in cols) + ")"
    source = f"{catalog}.{table}" if catalog else table
    n, h = con.execute(f"SELECT COUNT(*), COALESCE(SUM(hash({packed})::HUGEINT), 0) FROM {source} t").fetchone()
    return int(n), int(h)


def compact_family_file(path: Path | str, *, threads: int = 4, memory_limit: str = "8GB") -> dict[str, Any]:
    """Rewrite a staging file without free blocks or legacy indexes; replace it atomically.

    Every family table is copied whole and checked by row count and content fingerprint.
    """
    import os  # noqa: PLC0415

    path = Path(path).expanduser().resolve()
    tmp = path.with_name(path.name + ".compact.tmp")
    for stale in (tmp, Path(str(tmp) + ".wal")):
        stale.unlink(missing_ok=True)
    before = path.stat().st_size
    con = connect(tmp, threads=threads, memory_limit=memory_limit)
    try:
        con.execute(f"ATTACH {_q(path)} AS __pf_old (READ_ONLY)")
        present = [t for t in TABLES if _columns(con, t, "__pf_old")]
        con.execute(SCHEMA)
        checks = {}
        for table in present:
            _copy_rows(con, table, "__pf_old", "true")
            same_cols = set(_columns(con, table)) == set(_columns(con, table, "__pf_old"))
            new, old = _table_fingerprint(con, table), _table_fingerprint(con, table, "__pf_old")
            if new[0] != old[0] or (same_cols and new != old):
                raise FamilyError(f"{table}: copy differs from the original ({old} -> {new})")
            checks[table] = new[0]
        _finalize(con)
        con.execute("DETACH __pf_old")
        con.execute("CHECKPOINT")
    except Exception:
        con.close()
        tmp.unlink(missing_ok=True)
        Path(str(tmp) + ".wal").unlink(missing_ok=True)
        raise
    con.close()
    os.replace(tmp, path)
    return {"rows": checks, "bytes_before": before, "bytes_after": path.stat().st_size}


# --------------------------------------------------------------------------- #
# Staging build (a Sharur hierarchical clustering)
# --------------------------------------------------------------------------- #


def _read_tsv(con, table: str, path: Path) -> None:
    cols = "{'rep': 'VARCHAR', 'member': 'VARCHAR'}"
    con.execute(f"""CREATE OR REPLACE TEMP TABLE {table} AS SELECT * FROM read_csv(
        {_q(path)}, delim = '\t', header = false, columns = {cols}, quote = '', escape = '')""")


def _check_partition(con, table: str, label: str) -> None:
    n, distinct, reps, rep_rows = con.execute(f"""
        SELECT COUNT(*), COUNT(DISTINCT member), COUNT(DISTINCT rep),
               COUNT(*) FILTER (WHERE rep = member) FROM {table}""").fetchone()
    if n != distinct:
        raise FamilyError(f"{label}: {n - distinct:,} proteins appear in more than one cluster")
    if reps != rep_rows:
        raise FamilyError(f"{label}: {reps - rep_rows:,} representatives are missing from their own cluster")


def _staging_tmp(staging: Path) -> Path:
    tmp = staging.with_name(staging.name + ".tmp")
    for stale in (tmp, Path(str(tmp) + ".wal")):
        stale.unlink(missing_ok=True)
    return tmp


def build_staging(
    db: Path | str,
    clustering_dir: Path | str,
    staging: Path | str,
    *,
    export_manifest: dict | None = None,
    threads: int = 4,
    memory_limit: str = "8GB",
    verify_digest: bool = True,
    log: Callable[[str], None] = print,
) -> dict[str, Any]:
    """Write a staging DuckDB with a Sharur build (subfamily and family runs) for the proteins of ``db``.

    ``clustering_dir`` holds ``clustering.json``, ``level0.tsv.gz`` (subfamilies) and
    ``level1.tsv.gz`` (families of subfamily representatives) from
    :func:`sharur.protein_families.mmseqs.cluster_levels`; ``export.json`` beside them
    (or ``export_manifest``) records the clustered proteins' digest. ``db`` is
    opened read-only. The staging file is replaced.
    """
    clustering_dir = Path(clustering_dir)
    staging = Path(staging)
    manifest = json.loads((clustering_dir / "clustering.json").read_text())
    if len(manifest["levels"]) != 2:
        raise FamilyError(f"{clustering_dir}: expected 2 levels (subfamily, family), found {len(manifest['levels'])}")
    export_manifest = export_manifest or json.loads((clustering_dir / "export.json").read_text())
    levels = [LevelParams(**level["params"]) for level in manifest["levels"]]
    db = Path(db).expanduser().resolve()
    staging.parent.mkdir(parents=True, exist_ok=True)
    tmp = _staging_tmp(staging)
    con = connect(tmp, threads=threads, memory_limit=memory_limit)
    try:
        con.execute(f"ATTACH {_q(db)} AS __pf_src (READ_ONLY)")
        started = time.time()
        if verify_digest:
            live = export_proteins(con, None, prefix=SRC)
            if live["digest"] != export_manifest["digest"]:
                raise FamilyError(
                    f"{db}: protein sequences differ from the clustered export "
                    f"({live['proteins']:,} live vs {export_manifest['proteins']:,} exported proteins); "
                    "re-export and re-cluster")
            log(f"digest verified over {live['proteins']:,} proteins ({time.time() - started:.0f} s)")
        _read_tsv(con, "l0", clustering_dir / "level0.tsv.gz")
        _read_tsv(con, "l1", clustering_dir / "level1.tsv.gz")
        _check_partition(con, "l0", "subfamily level")
        _check_partition(con, "l1", "family level")
        n0, n_sub, n1, reps_unplaced, unknown = con.execute("""
            SELECT (SELECT COUNT(*) FROM l0), (SELECT COUNT(DISTINCT rep) FROM l0), (SELECT COUNT(*) FROM l1),
                   (SELECT COUNT(*) FROM (SELECT DISTINCT rep FROM l0) r ANTI JOIN l1 ON l1.member = r.rep),
                   (SELECT COUNT(*) FROM l0 ANTI JOIN __pf_src.proteins p ON p.protein_id = l0.member)""").fetchone()
        if n0 != export_manifest["proteins"]:
            raise FamilyError(f"clustering holds {n0:,} proteins, the export {export_manifest['proteins']:,}")
        if n_sub != n1 or reps_unplaced:
            raise FamilyError(f"{reps_unplaced:,} subfamily representatives lack a family")
        if unknown:
            raise FamilyError(f"{unknown:,} clustered proteins are absent from {db}")
        partition = con.execute("""
            SELECT sha256(string_agg(l0.member || chr(9) || l0.rep || chr(9) || l1.rep, chr(10)
                                     ORDER BY l0.member))
            FROM l0 JOIN l1 ON l1.member = l0.rep""").fetchone()[0]
        build_id = "pfr-" + _sha({"levels": [lv.to_dict() for lv in levels], "tool": manifest["tool"],
                                  "tool_version": manifest["tool_version"], "input_digest": export_manifest["digest"],
                                  "partition_digest": partition})[:10]
        tag = build_tag(build_id)
        family_run, subfamily_run = f"{build_id}-family", f"{build_id}-subfamily"
        con.execute(SCHEMA)
        common = {"origin": "sharur", "build_id": build_id, "is_default": True, "tool": manifest["tool"],
                  "tool_version": manifest["tool_version"], "input_proteins": n0,
                  "input_digest": export_manifest["digest"], "source_db": str(db)}
        details = {"partition_digest": partition, "threads": manifest.get("threads"),
                   "seconds": manifest.get("seconds"), "skipped_empty": export_manifest.get("skipped_empty", 0)}
        _insert_run(con, clustering_run=family_run, level="family", parent_run=None,
                    parameters={**levels[1].to_dict(), "clusters": "subfamily representatives"},
                    details={**details, "level_seconds": manifest["levels"][1].get("seconds")}, **common)
        _insert_run(con, clustering_run=subfamily_run, level="subfamily", parent_run=family_run,
                    parameters={**levels[0].to_dict(), "clusters": "all proteins"},
                    details={**details, "level_seconds": manifest["levels"][0].get("seconds")}, **common)
        con.execute(f"""
            INSERT INTO {MEMBERS}
            SELECT l0.member, 'fam:{tag}:' || l1.rep, l0.member = l1.rep, '{family_run}', p.bin_id, 'cluster', NULL
            FROM l0 JOIN l1 ON l1.member = l0.rep JOIN __pf_src.proteins p ON p.protein_id = l0.member""")
        con.execute(f"""
            INSERT INTO {MEMBERS}
            SELECT l0.member, 'sub:{tag}:' || l0.rep, l0.member = l0.rep, '{subfamily_run}', p.bin_id, 'cluster', NULL
            FROM l0 JOIN __pf_src.proteins p ON p.protein_id = l0.member""")
        log(f"build {build_id}: {n0:,} proteins, {n_sub:,} subfamilies")
        counts = build_stats(con, [family_run, subfamily_run], prefix=SRC)
        _finalize(con)
        con.execute("DETACH __pf_src")
        con.execute("CHECKPOINT")
    except Exception:
        con.close()
        tmp.unlink(missing_ok=True)
        Path(str(tmp) + ".wal").unlink(missing_ok=True)
        raise
    con.close()
    tmp.replace(staging)
    return {"build_id": build_id, "family_run": family_run, "subfamily_run": subfamily_run, "proteins": n0,
            "subfamilies": counts[subfamily_run], "families": counts[family_run], "staging": str(staging)}


def _insert_run(con, *, clustering_run: str, level: str, parent_run: str | None, build_id: str | None,
                origin: str, is_default: bool, tool: str | None, tool_version: str | None, parameters: Any,
                input_proteins: int, input_digest: str, source_db: str | None, details: dict) -> None:
    con.execute(f"INSERT INTO {RUNS} VALUES (?, ?, ?, ?, ?, ?, ?::TIMESTAMP, ?, ?, ?, ?, ?, NULL, ?, ?, ?)",
                [clustering_run, level, parent_run, build_id, origin, is_default, _now(), tool, tool_version,
                 json.dumps(parameters), input_proteins, input_digest, source_db, SCHEMA_VERSION,
                 json.dumps(details)])


# --------------------------------------------------------------------------- #
# Import of an external clustering
# --------------------------------------------------------------------------- #


def import_clustering(
    db: Path | str,
    tsv: Path | str,
    staging: Path | str,
    *,
    name: str,
    level: str = "imported",
    parameters: Any = None,
    tool: str = "mmseqs2",
    tool_version: str | None = None,
    set_default: bool = False,
    threads: int = 4,
    memory_limit: str = "8GB",
    log: Callable[[str], None] = print,
) -> dict[str, Any]:
    """Register an externally produced clustering as run ``name`` in ``staging``.

    ``tsv`` is a two-column ``representative<TAB>member`` table, as ``mmseqs
    createtsv`` writes it (plain or gzipped). Every member must be a protein of
    ``db`` (opened read-only) and appear once; each representative must be a
    member of its own cluster. ``parameters`` (the command line or a mapping) is
    stored verbatim. The staging file is created when absent and a run of the same
    name in it is replaced; other runs are kept. With ``set_default`` the run
    becomes the default of its ``level``.
    """
    if not _RUN_NAME.fullmatch(name) or name.startswith("pfr-"):
        raise FamilyError(f"run name {name!r}: use letters, digits, '.', '_' or '-' (the pfr- prefix is reserved)")
    db = Path(db).expanduser().resolve()
    staging = Path(staging)
    staging.parent.mkdir(parents=True, exist_ok=True)
    con = connect(staging, threads=threads, memory_limit=memory_limit)
    try:
        con.execute(f"ATTACH {_q(db)} AS __pf_src (READ_ONLY)")
        _read_tsv(con, "ext", Path(tsv))
        _check_partition(con, "ext", str(tsv))
        n, unknown, n_clusters = con.execute("""
            SELECT COUNT(*), (SELECT COUNT(*) FROM ext ANTI JOIN __pf_src.proteins p ON p.protein_id = ext.member),
                   COUNT(DISTINCT rep) FROM ext""").fetchone()
        if unknown:
            raise FamilyError(f"{unknown:,} of {n:,} clustered proteins are absent from {db}")
        not_covered = con.execute("SELECT COUNT(*) FROM __pf_src.proteins p ANTI JOIN ext ON ext.member = p.protein_id"
                                  ).fetchone()[0]
        digest = _digest_of(con, "SELECT member AS protein_id FROM ext", SRC)
        partition = con.execute("SELECT sha256(string_agg(member || chr(9) || rep, chr(10) ORDER BY member)) "
                                "FROM ext").fetchone()[0]
        tag = _sha({"name": name, "parameters": parameters, "digest": digest["digest"], "partition": partition})[:10]
        con.execute("BEGIN TRANSACTION")
        try:
            con.execute(SCHEMA)
            for table in (MEMBERS, STATS, RUNS):
                con.execute(f"DELETE FROM {table} WHERE clustering_run = ?", [name])
            if set_default:
                con.execute(f"UPDATE {RUNS} SET is_default = false WHERE level = ?", [level])
            _insert_run(con, clustering_run=name, level=level, parent_run=None, build_id=None, origin="import",
                        is_default=set_default, tool=tool, tool_version=tool_version, parameters=parameters,
                        input_proteins=digest["proteins"], input_digest=digest["digest"], source_db=str(db),
                        details={"tag": tag, "partition_digest": partition, "source_tsv": str(Path(tsv).resolve()),
                                 "proteins_not_covered": not_covered})
            con.execute(f"""
                INSERT INTO {MEMBERS}
                SELECT ext.member, 'fam:{tag}:' || ext.rep, ext.member = ext.rep, ?, p.bin_id, 'import', NULL
                FROM ext JOIN __pf_src.proteins p ON p.protein_id = ext.member""", [name])
            counts = build_stats(con, [name], prefix=SRC)
            _finalize(con)
            con.execute("COMMIT")
        except Exception:
            con.execute("ROLLBACK")
            raise
        con.execute("DETACH __pf_src")
        con.execute("CHECKPOINT")
    finally:
        con.close()
    log(f"run {name}: {n:,} proteins in {n_clusters:,} clusters; {not_covered:,} proteins of the database not covered")
    return {"clustering_run": name, "level": level, "proteins": n, "families": counts[name],
            "proteins_not_covered": not_covered, "tag": tag, "staging": str(staging)}


# --------------------------------------------------------------------------- #
# Activation (copy staged runs into a dataset database)
# --------------------------------------------------------------------------- #


def _check_digests(con, alias: str, runs: list[str]) -> None:
    verified: set[tuple[str, str]] = set()  # the levels of one build share their clustered proteins
    for run, n, recorded, build in con.execute(
            f"SELECT clustering_run, input_proteins, input_digest, build_id FROM {alias}.{RUNS} "
            "WHERE list_contains(?::VARCHAR[], clustering_run) ORDER BY clustering_run", [runs]).fetchall():
        if build and (build, recorded) in verified:
            continue
        verified.add((build, recorded))
        live = _digest_of(con, f"SELECT protein_id FROM {alias}.{MEMBERS} WHERE clustering_run = '{run}' "
                               "AND batch_id IS NULL", "")
        if live["digest"] != recorded or live["proteins"] != n:
            raise FamilyError(f"run {run}: sequences of its {n:,} proteins changed since clustering; re-cluster")
    for batch, n, recorded, family_run in con.execute(
            f"SELECT batch_id, input_proteins, input_digest, family_run FROM {alias}.{BATCHES} "
            "WHERE list_contains(?::VARCHAR[], family_run)", [runs]).fetchall():
        live = _digest_of(con, f"SELECT protein_id FROM {alias}.{MEMBERS} WHERE clustering_run = '{family_run}' "
                               f"AND batch_id = '{batch}'", "")
        if live["digest"] != recorded or live["proteins"] != n:
            raise FamilyError(f"assignment batch {batch}: sequences changed since assignment")


def load_staging(con: duckdb.DuckDBPyConnection, staging: Path | str, *, runs: Sequence[str] | None = None,
                 verify_digest: bool = True, log: Callable[[str], None] = print) -> dict[str, Any]:
    """Copy runs of ``staging`` (all, or ``runs``) into the database open on ``con``.

    Runs of the same name are replaced; other runs in the database are kept. A
    loaded default run takes over its level's default. Every staged protein must
    exist in ``proteins`` with the same ``bin_id`` and, with ``verify_digest``, the
    sequences of each run's and assignment batch's proteins must match the digest
    recorded when they were clustered, imported or assigned.
    """
    con.execute(f"ATTACH {_q(Path(staging).expanduser().resolve())} AS __pf_stg (READ_ONLY)")
    try:
        if not has_families(con, catalog="__pf_stg"):
            raise FamilyError(f"{staging} holds no protein-family tables")
        staged = [r[0] for r in con.execute(f"SELECT clustering_run FROM __pf_stg.{RUNS} ORDER BY 1").fetchall()]
        selected = list(runs) if runs else staged
        unknown = sorted(set(selected) - set(staged))
        if unknown:
            raise FamilyError(f"{staging} holds no run {', '.join(unknown)}")
        con.execute("CREATE OR REPLACE TEMP TABLE _pf_load AS SELECT unnest(?::VARCHAR[]) AS clustering_run",
                    [selected])
        missing, moved = con.execute(f"""
            SELECT COUNT(*) FILTER (WHERE p.protein_id IS NULL),
                   COUNT(*) FILTER (WHERE p.protein_id IS NOT NULL AND p.bin_id IS DISTINCT FROM f.bin_id)
            FROM __pf_stg.{MEMBERS} f JOIN _pf_load USING (clustering_run)
            LEFT JOIN proteins p USING (protein_id)""").fetchone()
        if missing or moved:
            raise FamilyError(f"{missing:,} staged proteins are absent from the database and {moved:,} sit in "
                              "another genome; the staging file belongs to another protein set")
        if verify_digest:
            _check_digests(con, "__pf_stg", selected)
            log("sequence digests verified")
        con.execute("BEGIN TRANSACTION")
        try:
            for view in VIEWS:
                con.execute(f"DROP VIEW IF EXISTS {view}")
            con.execute(SCHEMA)
            ensure_schema(con)
            defaults = [r[0] for r in con.execute(
                f"SELECT DISTINCT level FROM __pf_stg.{RUNS} JOIN _pf_load USING (clustering_run) WHERE is_default"
            ).fetchall()]
            if defaults:
                con.execute(f"UPDATE {RUNS} SET is_default = false WHERE list_contains(?::VARCHAR[], level)",
                            [defaults])
            loaded = "(SELECT clustering_run FROM _pf_load)"
            for table in (MEMBERS, STATS, RUNS):
                con.execute(f"DELETE FROM {table} WHERE clustering_run IN {loaded}")
                _copy_rows(con, table, "__pf_stg", f"t.clustering_run IN {loaded}")
            for table in (BATCHES, RESCUES):
                con.execute(f"DELETE FROM {table} WHERE family_run IN {loaded}")
                if _columns(con, table, "__pf_stg"):
                    _copy_rows(con, table, "__pf_stg", f"t.family_run IN {loaded}")
            _finalize(con)
            con.execute("COMMIT")
        except Exception:
            con.execute("ROLLBACK")
            raise
        con.execute("DROP TABLE IF EXISTS _pf_load")
    finally:
        con.execute("DETACH __pf_stg")
    return status(con)


def _table_exists(con, name: str) -> bool:
    return bool(con.execute("SELECT 1 FROM information_schema.tables WHERE table_schema = 'main' AND table_name = ? "
                            "AND table_catalog = current_database()", [name]).fetchone())


def base_run_of(con, run: str) -> str:
    """The clustering run a rescue run was derived from (the run itself otherwise)."""
    origin, details = con.execute(f"SELECT origin, details FROM {RUNS} WHERE clustering_run = ?", [run]).fetchone()
    return json.loads(details)["base_run"] if origin == "rescue" else run


def rescue_runs_of(con, family_run: str) -> tuple[str | None, str | None]:
    """The rescue runs (family, subfamily) derived from a base family run, if any."""
    family = con.execute(f"SELECT min(clustering_run) FROM {RUNS} WHERE origin = 'rescue' AND "
                         "json_extract_string(details, '$.base_run') = ?", [family_run]).fetchone()[0]
    sub = (con.execute(f"SELECT min(clustering_run) FROM {RUNS} WHERE parent_run = ?", [family]).fetchone()[0]
           if family else None)
    return family, sub


def default_runs(con) -> tuple[str | None, str | None]:
    """The default family run and its subfamily run (``None`` when absent)."""
    family = con.execute(f"SELECT min(clustering_run) FROM {RUNS} WHERE level = 'family' AND is_default"
                         ).fetchone()[0]
    sub = (con.execute(f"SELECT min(clustering_run) FROM {RUNS} WHERE parent_run = ?", [family]).fetchone()[0]
           if family else None)
    return family, sub


def status(con) -> dict[str, Any]:
    """Runs and coverage of the family layer of the database open on ``con``."""
    if not has_families(con):
        return {"present": False}
    runs = []
    for row in con.execute(f"""
            SELECT r.clustering_run, r.level, r.parent_run, r.build_id, r.origin, r.is_default, r.families,
                   r.input_proteins, (SELECT COUNT(*) FROM {MEMBERS} m WHERE m.clustering_run = r.clustering_run),
                   r.parameters, r.created_at
            FROM {RUNS} r ORDER BY r.created_at, r.clustering_run""").fetchall():
        runs.append(dict(zip(["clustering_run", "level", "parent_run", "build_id", "origin", "is_default",
                              "families", "clustered_proteins", "proteins", "parameters", "created_at"],
                             [*row[:9], json.loads(row[9]) if row[9] else None, str(row[10])])))
    family_run, subfamily_run = default_runs(con)
    has_proteins = _table_exists(con, "proteins")
    out: dict[str, Any] = {"present": True, "runs": runs, "default_family_run": family_run,
                           "default_subfamily_run": subfamily_run,
                           "proteins": con.execute("SELECT COUNT(*) FROM proteins").fetchone()[0] if has_proteins
                           else None}
    if family_run is None:
        out.update({"proteins_with_family": 0, "proteins_without_family": out["proteins"], "by_assignment": {},
                    "families": 0, "subfamilies": 0, "assigned_fraction": 0.0, "recluster_recommended": False})
        return out
    by_assignment = dict(con.execute(f"SELECT assignment, COUNT(*) FROM {MEMBERS} WHERE clustering_run = ? "
                                     "GROUP BY 1 ORDER BY 1", [family_run]).fetchall())
    total = sum(by_assignment.values())
    assigned = total - sum(by_assignment.get(k, 0) for k in ("cluster", "import", "rescue"))
    out.update({
        "proteins_with_family": total,
        "proteins_without_family": (con.execute(
            f"SELECT COUNT(*) FROM proteins p ANTI JOIN (SELECT protein_id FROM {MEMBERS} "
            "WHERE clustering_run = ?) f USING (protein_id)", [family_run]).fetchone()[0] if has_proteins else None),
        "by_assignment": by_assignment,
        "families": next(r["families"] for r in runs if r["clustering_run"] == family_run),
        "subfamilies": next((r["families"] for r in runs if r["clustering_run"] == subfamily_run), 0),
        "assigned_fraction": assigned / total if total else 0.0,
    })
    out["recluster_recommended"] = out["assigned_fraction"] > RECLUSTER_ASSIGNED_FRACTION
    base = base_run_of(con, family_run)
    rescue_family, _ = rescue_runs_of(con, base)
    if rescue_family:
        details = json.loads(con.execute(f"SELECT details FROM {RUNS} WHERE clustering_run = ?",
                                         [rescue_family]).fetchone()[0])
        out["rescue"] = {"family_run": rescue_family, "base_family_run": base,
                         "included_by_default": rescue_family == family_run,
                         **{k: details.get(k) for k in ("rescue_id", "candidates", "rescued", "joined_subfamily")}}
    return out


# --------------------------------------------------------------------------- #
# Incremental assignment
# --------------------------------------------------------------------------- #


@dataclass
class AssignResult:
    batch_id: str | None
    new_proteins: int
    by_assignment: dict[str, int] = field(default_factory=dict)
    new_subfamilies: int = 0
    new_families: int = 0
    rescue: dict[str, Any] | None = None

    def to_dict(self) -> dict[str, Any]:
        return self.__dict__.copy()


def _fasta_for(con, ids: Iterable[str], fasta: Path) -> int:
    import pandas as pd  # noqa: PLC0415

    con.register("_pf_ids_df", pd.DataFrame({"protein_id": sorted(set(ids))}))
    con.execute("CREATE OR REPLACE TEMP TABLE _pf_sel AS SELECT * FROM _pf_ids_df")
    con.unregister("_pf_ids_df")
    try:
        return export_proteins(con, fasta, ids_table="_pf_sel")["proteins"]
    finally:
        con.execute("DROP TABLE IF EXISTS _pf_sel")


def _level_params(parameters: str) -> LevelParams:
    data = json.loads(parameters)
    return LevelParams(**{k: v for k, v in data.items() if k in LevelParams.__dataclass_fields__})


def assign(con: duckdb.DuckDBPyConnection, work_dir: Path | str, *, family_run: str | None = None,
           threads: int | None = None, binary: str | None = None,
           log: Callable[[str], None] = print) -> AssignResult:
    """Give every protein without a row in a Sharur build a family and subfamily, without re-clustering.

    The build is the default family run and its subfamily run, or ``family_run``.

    1. Search the new proteins against the subfamily representatives with the
       subfamily criteria; a hit joins that subfamily (``search``).
    2. Cluster the remaining proteins with the subfamily criteria; search each new
       subfamily's representative against the family representatives with the
       family criteria; a hit joins that family (``new_subfamily``).
    3. Cluster the remaining new subfamily representatives with the family
       criteria into new families (``new_family``).

    Existing families are never merged or split. Statistics of both runs are rebuilt.
    """
    if not has_families(con):
        raise FamilyError("the database holds no protein-family tables; activate a staging file first")
    family_run = family_run or default_runs(con)[0]
    if family_run is None:
        raise FamilyError("no default family run; name one with family_run")
    family_run = base_run_of(con, family_run)
    subfamily_run = con.execute(f"SELECT min(clustering_run) FROM {RUNS} WHERE parent_run = ?",
                                [family_run]).fetchone()[0]
    rows = {r[0]: r[1:] for r in con.execute(
        f"SELECT clustering_run, origin, build_id, parameters FROM {RUNS} "
        "WHERE list_contains(?::VARCHAR[], clustering_run)", [[family_run, subfamily_run]]).fetchall()}
    if family_run not in rows or subfamily_run not in rows or rows[family_run][0] != "sharur":
        raise FamilyError("assignment needs a Sharur build (a family run with its subfamily run)")
    build_id = rows[family_run][1]
    tag = build_tag(build_id)
    fam_level, sub_level = _level_params(rows[family_run][2]), _level_params(rows[subfamily_run][2])
    new = [r[0] for r in con.execute(
        f"SELECT protein_id FROM proteins ANTI JOIN (SELECT protein_id FROM {MEMBERS} WHERE clustering_run = ?) f "
        "USING (protein_id) WHERE rtrim(COALESCE(sequence, ''), '*') <> '' ORDER BY protein_id",
        [family_run]).fetchall()]
    if not new:
        return AssignResult(None, 0)
    missing_reps = con.execute(f"""
        SELECT COUNT(*) FROM (SELECT protein_id FROM {MEMBERS} WHERE is_representative
                              AND list_contains(?::VARCHAR[], clustering_run)) r
        ANTI JOIN proteins USING (protein_id)""", [[family_run, subfamily_run]]).fetchone()[0]
    if missing_reps:
        raise FamilyError(f"{missing_reps:,} representatives are gone from proteins; a full re-cluster is required")
    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)
    new_fasta, sub_fasta, fam_fasta = work / "new.faa", work / "subfamily_reps.faa", work / "family_reps.faa"
    _fasta_for(con, new, new_fasta)
    sub_of_rep = dict(con.execute(f"SELECT protein_id, family_id FROM {MEMBERS} WHERE clustering_run = ? "
                                  "AND is_representative", [subfamily_run]).fetchall())
    fam_of_rep = dict(con.execute(f"SELECT protein_id, family_id FROM {MEMBERS} WHERE clustering_run = ? "
                                  "AND is_representative", [family_run]).fetchall())
    fam_of_sub = dict(con.execute(f"SELECT s.family_id, f.family_id FROM {MEMBERS} s JOIN {MEMBERS} f "
                                  "ON f.protein_id = s.protein_id AND f.clustering_run = ? "
                                  "WHERE s.clustering_run = ? AND s.is_representative",
                                  [family_run, subfamily_run]).fetchall())
    _fasta_for(con, sub_of_rep, sub_fasta)
    _fasta_for(con, fam_of_rep, fam_fasta)
    placed: dict[str, tuple[str, str, str]] = {}  # protein -> (subfamily_id, family_id, assignment)
    new_reps: set[str] = set()
    hits = mmseqs.search_best(new_fasta, sub_fasta, work / "search_subfamilies", sub_level,
                              threads=threads, binary=binary)
    for query, hit in hits.items():
        sid = sub_of_rep[hit[0]]
        placed[query] = (sid, fam_of_sub[sid], "search")
    log(f"{len(hits):,} of {len(new):,} new proteins joined existing subfamilies")
    rest = [p for p in new if p not in placed]
    new_sub = new_fam = 0
    if rest:
        _fasta_for(con, rest, work / "unplaced.faa")
        mmseqs.cluster_levels(work / "new_subfamilies", [sub_level], fasta=work / "unplaced.faa",
                              threads=threads, binary=binary)
        members_of: dict[str, list[str]] = {}
        for rep, member in mmseqs.read_level(work / "new_subfamilies" / "level0.tsv.gz"):
            members_of.setdefault(rep, []).append(member)
        new_sub = len(members_of)
        new_reps.update(members_of)
        _fasta_for(con, members_of, work / "new_subfamily_reps.faa")
        fam_hits = mmseqs.search_best(work / "new_subfamily_reps.faa", fam_fasta, work / "search_families",
                                      fam_level, threads=threads, binary=binary)
        family_for_rep = {rep: fam_of_rep[hit[0]] for rep, hit in fam_hits.items()}
        orphans = [rep for rep in members_of if rep not in family_for_rep]
        if orphans:
            _fasta_for(con, orphans, work / "orphan_reps.faa")
            mmseqs.cluster_levels(work / "new_families", [fam_level], fasta=work / "orphan_reps.faa",
                                  threads=threads, binary=binary)
            fam_reps = set()
            for rep, member in mmseqs.read_level(work / "new_families" / "level0.tsv.gz"):
                family_for_rep[member] = family_id(tag, rep)
                fam_reps.add(rep)
            new_fam = len(fam_reps)
            new_reps.update(fam_reps)
        for rep, members in members_of.items():
            label = "new_subfamily" if rep in fam_hits else "new_family"
            for member in members:
                placed[member] = (subfamily_id(tag, rep), family_for_rep[rep], label)
    if len(placed) != len(new):
        raise FamilyError(f"assignment placed {len(placed):,} of {len(new):,} proteins")

    import pandas as pd  # noqa: PLC0415

    frame = pd.DataFrame([(p, *placed[p]) for p in new], columns=["protein_id", "subfamily_id", "family_id",
                                                                  "assignment"])
    frame["is_sub_rep"] = [p in new_reps and placed[p][0] == subfamily_id(tag, p) for p in new]
    frame["is_fam_rep"] = [placed[p][1] == family_id(tag, p) for p in new]
    con.register("_pf_new_df", frame)
    con.execute("CREATE OR REPLACE TEMP TABLE _pf_new AS SELECT * FROM _pf_new_df")
    con.unregister("_pf_new_df")
    digest = export_proteins(con, None, ids_table="_pf_new")
    batch_id = "pfa-" + _sha({"build": build_id, "digest": digest["digest"],
                              "rows": hashlib.sha256(frame.sort_values("protein_id").to_csv(index=False).encode())
                              .hexdigest()})[:10]
    by_assignment = {k: int(v) for k, v in frame["assignment"].value_counts().items()}
    con.execute("BEGIN TRANSACTION")
    try:
        for run, id_col, rep_col in ((family_run, "family_id", "is_fam_rep"),
                                     (subfamily_run, "subfamily_id", "is_sub_rep")):
            con.execute(f"""
                INSERT INTO {MEMBERS}
                SELECT n.protein_id, n.{id_col}, n.{rep_col}, ?, p.bin_id, n.assignment, ?
                FROM _pf_new n JOIN proteins p USING (protein_id)""", [run, batch_id])
        build_stats(con, [family_run, subfamily_run])
        con.execute(f"INSERT INTO {BATCHES} VALUES (?, ?, ?, ?, ?::TIMESTAMP, ?, ?, ?, ?, ?)",
                    [batch_id, build_id, family_run, subfamily_run, _now(),
                     mmseqs.mmseqs_version(mmseqs.resolve_binary(binary)),
                     json.dumps({"subfamily": sub_level.to_dict(), "family": fam_level.to_dict(),
                                 "search_max_seqs": mmseqs.SEARCH_MAX_SEQS}),
                     digest["proteins"], digest["digest"],
                     json.dumps({"by_assignment": by_assignment, "new_subfamilies": new_sub,
                                 "new_families": new_fam})])
        _finalize(con)
        con.execute("COMMIT")
    except Exception:
        con.execute("ROLLBACK")
        raise
    con.execute("DROP TABLE IF EXISTS _pf_new")
    shutil.rmtree(work, ignore_errors=True)
    result = AssignResult(batch_id, len(new), by_assignment, new_sub, new_fam)
    rescue_family, _ = rescue_runs_of(con, family_run)
    if rescue_family:  # the rescue layer covers the build's proteins: recompute it with its own parameters
        row = con.execute(f"SELECT parameters, is_default FROM {RUNS} WHERE clustering_run = ?",
                          [rescue_family]).fetchone()
        params = RescueParams(**json.loads(row[0])["rescue"]["params"])
        result.rescue = rescue(con, work / "rescue", params=params, family_run=family_run, set_default=row[1],
                               threads=threads, binary=binary, log=log).to_dict()
        shutil.rmtree(work, ignore_errors=True)
    return result


# --------------------------------------------------------------------------- #
# Rescue of truncated gene calls (a layer of derived runs over a build)
# --------------------------------------------------------------------------- #

TRUNCATED_FLAGS = ("10", "01", "11")
# Whether Stage 09 and `sharur families run` add the rescue layer, and whether an added layer becomes the
# default (views include rescued members), unless asked otherwise. False: the layer failed its pre-set
# validation criteria (docs/guides/protein-families.md), so it is opt-in and, when added, excluded.
RESCUE_BY_DEFAULT = False


@dataclass(frozen=True)
class RescueParams:
    """Acceptance rule for placing a truncated singleton into an existing family."""

    min_query_coverage: float = 0.8   # share of the CANDIDATE the alignment must cover
    min_seq_id: float | None = None   # None: the family run's identity floor
    min_family_size: int = 2          # members of the hit's family (the candidate excluded)
    sensitivity: float = mmseqs.SENSITIVITY
    evalue: float = mmseqs.EVALUE
    max_seqs: int = mmseqs.SEARCH_MAX_SEQS

    def __post_init__(self) -> None:
        if not 0 < self.min_query_coverage <= 1:
            raise ValueError(f"min_query_coverage must be in (0, 1], got {self.min_query_coverage}")
        if self.min_seq_id is not None and not 0 <= self.min_seq_id <= 1:
            raise ValueError(f"min_seq_id must be in [0, 1], got {self.min_seq_id}")
        if self.min_family_size < 2:
            raise ValueError(f"min_family_size must be at least 2, got {self.min_family_size}")

    def to_dict(self) -> dict[str, Any]:
        return {k: getattr(self, k) for k in self.__dataclass_fields__}


@dataclass
class RescueResult:
    rescue_id: str | None
    family_run: str | None
    subfamily_run: str | None
    candidates: int
    rescued: int
    joined_subfamily: int
    included_by_default: bool

    def to_dict(self) -> dict[str, Any]:
        return self.__dict__.copy()


def _rescue_scope(con, family_run: str, subfamily_run: str, prefix: str, min_family_size: int) -> int:
    """Temp tables ``_pf_cand`` (truncated singletons) and ``_pf_targ`` (members of families of
    >= ``min_family_size``) for one base build; returns the candidate count."""
    flags = "', '".join(TRUNCATED_FLAGS)
    has_partial = "partial" in _columns(con, "proteins", prefix.rstrip(".") or None)
    partial = "p.partial" if has_partial else "NULL"
    con.execute(f"""
        CREATE OR REPLACE TEMP TABLE _pf_cand AS
        SELECT f.protein_id, f.bin_id, {partial} AS partial, f.family_id AS original_family_id,
               s.family_id AS original_subfamily_id
        FROM {MEMBERS} f
        JOIN {STATS} st ON st.clustering_run = f.clustering_run AND st.family_id = f.family_id
        JOIN {prefix}proteins p USING (protein_id)
        LEFT JOIN {MEMBERS} s ON s.clustering_run = ? AND s.protein_id = f.protein_id
        WHERE f.clustering_run = ? AND st.n_members = 1 AND {partial} IN ('{flags}')
          AND rtrim(COALESCE(p.sequence, ''), '*') <> ''""", [subfamily_run, family_run])
    con.execute(f"""
        CREATE OR REPLACE TEMP TABLE _pf_targ AS
        SELECT f.protein_id, f.family_id, s.family_id AS subfamily_id, st.n_members
        FROM {MEMBERS} f
        JOIN {STATS} st ON st.clustering_run = f.clustering_run AND st.family_id = f.family_id
        LEFT JOIN {MEMBERS} s ON s.clustering_run = ? AND s.protein_id = f.protein_id
        WHERE f.clustering_run = ? AND st.n_members >= ?""", [subfamily_run, family_run, min_family_size])
    return con.execute("SELECT COUNT(*) FROM _pf_cand").fetchone()[0]


def rescue_inputs(con, out_dir: Path | str, *, prefix: str = "", family_run: str | None = None,
                  params: RescueParams | None = None) -> dict[str, Any]:
    """Write ``candidates.faa`` and ``targets.faa`` for a rescue search run elsewhere (e.g. a cluster node)."""
    params = params or RescueParams()
    family_run, subfamily_run, _ = _base_build(con, family_run)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    n = _rescue_scope(con, family_run, subfamily_run, prefix, params.min_family_size)
    cand = export_proteins(con, out_dir / "candidates.faa", prefix=prefix, ids_table="_pf_cand")
    con.execute("CREATE OR REPLACE TEMP TABLE _pf_targ_ids AS SELECT protein_id FROM _pf_targ")
    targ = export_proteins(con, out_dir / "targets.faa", prefix=prefix, ids_table="_pf_targ_ids")
    con.execute("DROP TABLE IF EXISTS _pf_targ_ids")
    return {"family_run": family_run, "candidates": n, "candidate_digest": cand["digest"],
            "targets": targ["proteins"], "target_digest": targ["digest"]}


def _base_build(con, family_run: str | None) -> tuple[str, str, str]:
    family_run = family_run or default_runs(con)[0]
    if family_run is None:
        raise FamilyError("no default family run")
    family_run = base_run_of(con, family_run)
    origin, build_id = con.execute(f"SELECT origin, build_id FROM {RUNS} WHERE clustering_run = ?",
                                   [family_run]).fetchone()
    subfamily_run = con.execute(f"SELECT min(clustering_run) FROM {RUNS} WHERE parent_run = ?",
                                [family_run]).fetchone()[0]
    if origin != "sharur" or subfamily_run is None:
        raise FamilyError("the rescue needs a Sharur build (a family run with its subfamily run)")
    return family_run, subfamily_run, build_id


def drop_rescue(con, *, family_run: str | None = None) -> list[str]:
    """Remove the rescue layer of a build; its base runs become the defaults again if the layer was."""
    base = base_run_of(con, family_run or default_runs(con)[0])
    rescue_family, rescue_sub = rescue_runs_of(con, base)
    runs = [r for r in (rescue_family, rescue_sub) if r]
    if not runs:
        return []
    was_default = con.execute(f"SELECT is_default FROM {RUNS} WHERE clustering_run = ?", [rescue_family]).fetchone()[0]
    base_sub = con.execute(f"SELECT min(clustering_run) FROM {RUNS} WHERE parent_run = ?", [base]).fetchone()[0]
    con.execute("BEGIN TRANSACTION")
    try:
        for table in (MEMBERS, STATS, RUNS):
            con.execute(f"DELETE FROM {table} WHERE list_contains(?::VARCHAR[], clustering_run)", [runs])
        if _columns(con, RESCUES):
            con.execute(f"DELETE FROM {RESCUES} WHERE family_run = ?", [rescue_family])
        if was_default:
            con.execute(f"UPDATE {RUNS} SET is_default = true WHERE list_contains(?::VARCHAR[], clustering_run)",
                        [[base, base_sub]])
        _finalize(con)
        con.execute("COMMIT")
    except Exception:
        con.execute("ROLLBACK")
        raise
    return runs


def rescue(
    con: duckdb.DuckDBPyConnection,
    work_dir: Path | str,
    *,
    prefix: str = "",
    params: RescueParams | None = None,
    family_run: str | None = None,
    set_default: bool = RESCUE_BY_DEFAULT,
    hits: dict[str, list[mmseqs.Hit]] | None = None,
    tool_version: str | None = None,
    threads: int | None = None,
    binary: str | None = None,
    log: Callable[[str], None] = print,
) -> RescueResult:
    """Place truncated singletons into existing families as a separate layer over a Sharur build.

    Candidates are proteins whose gene call runs off a contig end (Prodigal partial
    flag ``10``, ``01`` or ``11``) and whose family in the build has one member. Each is
    searched against every member of families with at least ``min_family_size``
    members; the best hit by bit score is accepted when the alignment covers at least
    ``min_query_coverage`` of the candidate and identity reaches ``min_seq_id`` (the
    family run's floor by default). An accepted candidate joins the hit's family, and
    the hit's subfamily when identity also reaches the subfamily floor.

    The build's runs stay unchanged. The layer is a pair of derived runs,
    ``<run>-rescued``, holding the build's membership with the rescued proteins moved
    (``assignment = 'rescue'``), plus one ``protein_family_rescues`` row per rescued
    protein with the hit, identity and coverages. With ``set_default`` the derived runs
    become the defaults, so the per-level views include rescued members (counted in
    ``n_rescued``). An earlier layer of the same build is replaced. ``prefix`` names the
    catalog holding ``proteins`` (``SRC`` for a staging file); ``hits`` passes search
    results computed elsewhere (:func:`sharur.protein_families.mmseqs.read_hits` over the
    files of :func:`rescue_inputs`), with ``tool_version`` naming the MMseqs2 that made them.
    """
    import pandas as pd  # noqa: PLC0415

    params = params or RescueParams()
    ensure_schema(con)
    con.execute(DDL[RESCUES])
    family_run, subfamily_run, build_id = _base_build(con, family_run)
    run_params = {r[0]: json.loads(r[1]) for r in con.execute(
        f"SELECT clustering_run, parameters FROM {RUNS} WHERE list_contains(?::VARCHAR[], clustering_run)",
        [[family_run, subfamily_run]]).fetchall()}
    min_id = params.min_seq_id if params.min_seq_id is not None else run_params[family_run]["min_seq_id"]
    sub_floor = run_params[subfamily_run]["min_seq_id"]
    n_cand = _rescue_scope(con, family_run, subfamily_run, prefix, params.min_family_size)
    log(f"{n_cand:,} truncated singletons in {family_run}")
    n_targ = con.execute("SELECT COUNT(*) FROM _pf_targ").fetchone()[0]
    searched = hits is None and n_cand > 0 and n_targ > 0
    if searched:
        work = Path(work_dir)
        work.mkdir(parents=True, exist_ok=True)
        export_proteins(con, work / "candidates.faa", prefix=prefix, ids_table="_pf_cand")
        con.execute("CREATE OR REPLACE TEMP TABLE _pf_targ_ids AS SELECT protein_id FROM _pf_targ")
        export_proteins(con, work / "targets.faa", prefix=prefix, ids_table="_pf_targ_ids")
        con.execute("DROP TABLE IF EXISTS _pf_targ_ids")
        level = LevelParams(min_seq_id=min_id, coverage=params.min_query_coverage, cov_mode=2,
                            sensitivity=params.sensitivity, evalue=params.evalue)
        hits = mmseqs.search_hits(work / "candidates.faa", work / "targets.faa", work / "search", level, top=5,
                                  threads=threads, binary=binary, max_seqs=params.max_seqs)
        shutil.rmtree(work, ignore_errors=True)
    hits = hits or {}
    candidates = {r[0] for r in con.execute("SELECT protein_id FROM _pf_cand").fetchall()}
    targets = {r[0]: (r[1], r[2]) for r in con.execute(
        "SELECT t.protein_id, t.family_id, t.subfamily_id FROM _pf_targ t "
        "JOIN (SELECT DISTINCT target FROM (SELECT unnest(?::VARCHAR[]) AS target)) h ON h.target = t.protein_id",
        [sorted({h.target for hs in hits.values() for h in hs})]).fetchall()}
    rows = []
    for query in sorted(candidates & set(hits)):
        best = next((h for h in sorted(hits[query], key=mmseqs._rank)
                     if h.target in targets and h.target != query and h.fident >= min_id
                     and h.qcov >= params.min_query_coverage), None)
        if best is None:
            continue
        fam, sub = targets[best.target]
        joined_sub = best.fident >= sub_floor and sub is not None
        rows.append((query, fam, sub if joined_sub else None, joined_sub, best.target, best.fident, best.qcov,
                     best.tcov, best.bits, best.evalue, best.qlen, best.tlen))
    accepted = pd.DataFrame(rows, columns=["protein_id", "family_id", "subfamily_id", "joined_subfamily",
                                           "hit_protein_id", "identity", "query_coverage", "target_coverage",
                                           "bits", "evalue", "query_length", "target_length"])
    con.execute("""CREATE OR REPLACE TEMP TABLE _pf_acc_raw (
        protein_id VARCHAR, family_id VARCHAR, subfamily_id VARCHAR, joined_subfamily BOOLEAN,
        hit_protein_id VARCHAR, identity DOUBLE, query_coverage DOUBLE, target_coverage DOUBLE, bits DOUBLE,
        evalue DOUBLE, query_length INTEGER, target_length INTEGER)""")
    if len(accepted):
        con.register("_pf_acc_df", accepted)
        con.execute("INSERT INTO _pf_acc_raw SELECT * FROM _pf_acc_df")
        con.unregister("_pf_acc_df")
    con.execute("CREATE OR REPLACE TEMP TABLE _pf_acc AS SELECT a.*, c.bin_id, c.partial, c.original_family_id, "
                "c.original_subfamily_id FROM _pf_acc_raw a JOIN _pf_cand c USING (protein_id)")
    con.execute("DROP TABLE _pf_acc_raw")
    rescue_id = "pfx-" + _sha({"base": family_run, "params": params.to_dict(), "min_seq_id": min_id,
                               "rows": accepted.sort_values("protein_id")[["protein_id", "family_id",
                                                                           "hit_protein_id"]].to_csv(index=False)
                               if len(accepted) else ""})[:10]
    rescue_family, rescue_sub = f"{family_run}-rescued", f"{subfamily_run}-rescued"
    drop_rescue(con, family_run=family_run)
    con.execute("BEGIN TRANSACTION")
    try:
        con.execute(f"""
            INSERT INTO {MEMBERS} (protein_id, family_id, is_representative, clustering_run, bin_id, assignment, batch_id)
            SELECT f.protein_id, COALESCE(a.family_id, f.family_id),
                   CASE WHEN a.protein_id IS NULL THEN f.is_representative ELSE false END, ?, f.bin_id,
                   CASE WHEN a.protein_id IS NULL THEN f.assignment ELSE 'rescue' END, f.batch_id
            FROM {MEMBERS} f LEFT JOIN _pf_acc a USING (protein_id) WHERE f.clustering_run = ?""",
                    [rescue_family, family_run])
        con.execute(f"""
            INSERT INTO {MEMBERS} (protein_id, family_id, is_representative, clustering_run, bin_id, assignment, batch_id)
            SELECT f.protein_id, CASE WHEN a.joined_subfamily THEN a.subfamily_id ELSE f.family_id END,
                   CASE WHEN a.joined_subfamily THEN false ELSE f.is_representative END, ?, f.bin_id,
                   CASE WHEN a.protein_id IS NULL THEN f.assignment ELSE 'rescue' END, f.batch_id
            FROM {MEMBERS} f LEFT JOIN _pf_acc a USING (protein_id) WHERE f.clustering_run = ?""",
                    [rescue_sub, subfamily_run])
        con.execute(f"""
            INSERT INTO {RESCUES} (rescue_id, family_run, subfamily_run, base_family_run, protein_id, bin_id, partial,
                                   original_family_id, original_subfamily_id, family_id, subfamily_id,
                                   joined_subfamily, hit_protein_id, identity, query_coverage, target_coverage,
                                   bits, evalue, query_length, target_length)
            SELECT ?, ?, ?, ?, protein_id, bin_id, partial, original_family_id, original_subfamily_id, family_id,
                   subfamily_id, joined_subfamily, hit_protein_id, identity, query_coverage, target_coverage,
                   bits, evalue, query_length, target_length
            FROM _pf_acc""", [rescue_id, rescue_family, rescue_sub, family_run])
        digest = _digest_of(con, f"SELECT protein_id FROM {MEMBERS} WHERE clustering_run = '{rescue_family}' "
                                 "AND batch_id IS NULL", prefix)
        n_sub = int(accepted["joined_subfamily"].sum()) if len(accepted) else 0
        details = {"rescue_id": rescue_id, "candidates": n_cand, "rescued": len(accepted),
                   "joined_subfamily": n_sub}
        if tool_version is None and searched:
            tool_version = mmseqs.mmseqs_version(mmseqs.resolve_binary(binary))
        rescue_meta = {"params": params.to_dict(), "min_seq_id_applied": min_id, "subfamily_floor": sub_floor,
                       "candidates": "Prodigal partial " + "/".join(TRUNCATED_FLAGS) + ", family singleton"}
        if set_default:
            con.execute(f"UPDATE {RUNS} SET is_default = false WHERE level IN ('family', 'subfamily')")
        for run, base, level, parent in ((rescue_family, family_run, "family", None),
                                         (rescue_sub, subfamily_run, "subfamily", rescue_family)):
            _insert_run(con, clustering_run=run, level=level, parent_run=parent, build_id=build_id,
                        origin="rescue", is_default=set_default, tool="mmseqs2", tool_version=tool_version,
                        parameters={**run_params[base], "rescue": rescue_meta},
                        input_proteins=digest["proteins"], input_digest=digest["digest"], source_db=None,
                        details={**details, "base_run": base})
        build_stats(con, [rescue_family, rescue_sub], prefix=prefix)
        _finalize(con)
        con.execute("COMMIT")
    except Exception:
        con.execute("ROLLBACK")
        raise
    for temp in ("_pf_cand", "_pf_targ", "_pf_acc"):
        con.execute(f"DROP TABLE IF EXISTS {temp}")
    log(f"rescued {len(accepted):,} of {n_cand:,} truncated singletons ({n_sub:,} also joined a subfamily)")
    return RescueResult(rescue_id, rescue_family, rescue_sub, n_cand, len(accepted), n_sub, set_default)


def base_fingerprint(con, family_run: str | None = None) -> dict[str, list[int]]:
    """Order-free fingerprint (rows, summed row hashes) of a build's base runs, excluding the default flag.

    A rescue layer leaves it unchanged; :func:`rescue_staging` checks that.
    """
    family_run = base_run_of(con, family_run or default_runs(con)[0])
    sub = con.execute(f"SELECT min(clustering_run) FROM {RUNS} WHERE parent_run = ?", [family_run]).fetchone()[0]
    runs = [family_run, sub]
    out = {}
    # columns a rescue may add (ADDED_COLUMNS) or flip (is_default) stay out of the fingerprint
    for table, exclude in ((MEMBERS, ()), (STATS, ("n_rescued",)), (RUNS, ("is_default",))):
        cols = [c for c in _columns(con, table) if c not in exclude]
        packed = "struct_pack(" + ", ".join(f'"{c}" := t."{c}"' for c in cols) + ")"
        n, h = con.execute(f"SELECT COUNT(*), COALESCE(SUM(hash({packed})::HUGEINT), 0) FROM {table} t "
                           "WHERE list_contains(?::VARCHAR[], t.clustering_run)", [runs]).fetchone()
        out[table] = [int(n), int(h)]
    return out


def rescue_staging(
    db: Path | str,
    staging: Path | str,
    *,
    params: RescueParams | None = None,
    set_default: bool = RESCUE_BY_DEFAULT,
    drop: bool = False,
    hits: dict[str, list[mmseqs.Hit]] | None = None,
    tool_version: str | None = None,
    threads: int = 4,
    memory_limit: str = "8GB",
    mmseqs_threads: int | None = None,
    binary: str | None = None,
    compact: bool = True,
    log: Callable[[str], None] = print,
) -> dict[str, Any]:
    """Add (or with ``drop`` remove) the rescue layer of a staging file, replacing the file atomically.

    The layer is written into a copy-on-write clone of ``staging`` (proteins read from
    ``db``, read-only), the build's base runs are checked unchanged by fingerprint, and
    the clone replaces ``staging`` with one rename, so readers see the old or the new
    file whole. With ``compact`` the result is then rewritten without free blocks
    (:func:`compact_family_file`), also replaced atomically.
    """
    import os  # noqa: PLC0415

    from sharur.staging import clone_file  # noqa: PLC0415

    staging = Path(staging).expanduser().resolve()
    tmp = staging.with_name(staging.name + ".rescue.tmp")
    for stale in (tmp, Path(str(tmp) + ".wal")):
        stale.unlink(missing_ok=True)
    method = clone_file(staging, tmp)
    con = connect(tmp, threads=threads, memory_limit=memory_limit)
    try:
        con.execute(f"ATTACH {_q(Path(db).expanduser().resolve())} AS __pf_src (READ_ONLY)")
        before = base_fingerprint(con)
        if drop:
            out: dict[str, Any] = {"dropped": drop_rescue(con)}
        else:
            out = rescue(con, staging.parent / (staging.name + ".rescue-work"), prefix=SRC, params=params,
                         set_default=set_default, hits=hits, tool_version=tool_version, threads=mmseqs_threads,
                         binary=binary, log=log).to_dict()
        after = base_fingerprint(con)
        if before != after:
            raise FamilyError(f"the build's base runs changed during the rescue: {before} -> {after}")
        con.execute("DETACH __pf_src")
        con.execute("CHECKPOINT")
    except Exception:
        con.close()
        tmp.unlink(missing_ok=True)
        Path(str(tmp) + ".wal").unlink(missing_ok=True)
        raise
    con.close()
    os.replace(tmp, staging)
    out.update({"staging": str(staging), "clone": method, "base_fingerprint": before})
    if compact:
        out["compact"] = compact_family_file(staging, threads=threads, memory_limit=memory_limit)
    return out


# --------------------------------------------------------------------------- #
# Pipeline
# --------------------------------------------------------------------------- #


def write_export(db: Path | str, out_dir: Path | str, *, threads: int = 4, memory_limit: str = "8GB") -> dict:
    """Export the proteins of ``db`` to ``out_dir/proteins.faa`` and ``out_dir/export.json``."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    con = connect(db, read_only=True, threads=threads, memory_limit=memory_limit)
    try:
        result = export_proteins(con, out_dir / "proteins.faa")
    finally:
        con.close()
    result.update({"source_db": str(Path(db).expanduser().resolve()), "fasta": "proteins.faa",
                   "created_at": time.strftime("%Y-%m-%dT%H:%M:%S%z")})
    (out_dir / "export.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def gzip_file(path: Path) -> Path:
    target = path.with_name(path.name + ".gz")
    with open(path, "rb") as src, gzip.open(target, "wb", compresslevel=6) as dst:
        shutil.copyfileobj(src, dst, length=1 << 20)
    path.unlink()
    return target


def default_duckdb_memory_limit(fraction: float = 0.25) -> str:
    """DuckDB memory limit at ``fraction`` of the allocation (at least 2 GB)."""
    total = mmseqs.allocated_memory_bytes()
    return f"{max(2, int((total or 8 * 1024 ** 3) * fraction / 1024 ** 3))}GB"


def run_pipeline(
    db: Path | str,
    out_dir: Path | str,
    *,
    levels: list[LevelParams] | None = None,
    threads: int | None = None,
    duckdb_threads: int | None = None,
    memory_limit: str | None = None,
    split_memory_limit: str | None = None,
    keep_fasta: bool = False,
    binary: str | None = None,
    rescue_params: RescueParams | None = None,
    rescue_include: bool = RESCUE_BY_DEFAULT,
    log: Callable[[str], None] = print,
) -> dict[str, Any]:
    """Export the proteins of ``db``, cluster both levels, and write ``out_dir/protein_families.duckdb``.

    With ``rescue_params`` the truncated-gene rescue layer is added (:func:`rescue`), and made
    the default when ``rescue_include``.

    ``db`` is only read. ``out_dir`` receives ``export.json``, ``clustering/``
    (``level0.tsv.gz`` subfamilies, ``level1.tsv.gz`` families, ``clustering.json``)
    and the staging DuckDB; the exported FASTA is removed unless ``keep_fasta``.
    """
    db = Path(db).expanduser().resolve()
    out_dir = Path(out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    levels = levels or mmseqs.default_levels()
    threads = threads or mmseqs.allocated_cpus()
    duckdb_threads = duckdb_threads or threads
    memory_limit = memory_limit or default_duckdb_memory_limit()
    started = time.time()
    export = write_export(db, out_dir, threads=duckdb_threads, memory_limit=memory_limit)
    log(f"exported {export['proteins']:,} proteins (digest {export['digest'][:12]})")
    if not export["proteins"]:
        raise FamilyError(f"{db} holds no protein sequences")
    clustering = out_dir / "clustering"
    shutil.rmtree(clustering, ignore_errors=True)
    manifest = mmseqs.cluster_levels(clustering, levels, fasta=out_dir / "proteins.faa", threads=threads,
                                     binary=binary, split_memory_limit=split_memory_limit)
    shutil.copyfile(out_dir / "export.json", clustering / "export.json")
    result = build_staging(db, clustering, out_dir / "protein_families.duckdb", export_manifest=export,
                           threads=duckdb_threads, memory_limit=memory_limit, log=log)
    if not keep_fasta:
        (out_dir / "proteins.faa").unlink(missing_ok=True)
    if rescue_params is not None:
        result["rescue"] = rescue_staging(db, out_dir / "protein_families.duckdb", params=rescue_params,
                                          set_default=rescue_include, threads=duckdb_threads, memory_limit=memory_limit,
                                          mmseqs_threads=threads, binary=binary, log=log)
    result.update({"seconds": round(time.time() - started, 1), "clustering_seconds": manifest["seconds"],
                   "mmseqs_version": manifest["tool_version"], "threads": threads,
                   "levels": [lv.to_dict() for lv in levels]})
    return result
