"""Remote-homology layer in DuckDB: input export, result loading, staging, activation.

Work directory (``export``)::

    WORK/export.json                 selection rule, clustering run, counts, sequence digest
    WORK/families.tsv                one row per selected family (family_index, family_key, ...)
    WORK/fasta/<bucket>/<key>.faa    up to ``max_sequences`` unique member sequences, representative first
    WORK/out/<bucket>/<key>/         runner.py output: .hhr(.gz) files, hits.tsv, done.json
    WORK/cdd/chunks/chunk_N.faa      unique sequences of every protein, named by MD5 (``export-cdd``)
    WORK/cdd/out/chunk_N.tsv         runner.py cdd output

Staging DuckDB (``load-families``, ``load-cdd``)::

    family_annotation_runs     one row per family-tier search: clustering run, selection rule,
                               databases with names and versions, parameters, tool versions, counts
    family_annotation_status   one row per selected family: ok / error / pending, timings
    family_annotations         one row per hit at probability >= 50% (the first 20 for UniRef30)
    cdd_hits                   CDD hits in the ``annotations`` shape plus ``sequence_md5``
    cdd_runs                   the CDD search: database name and version, E-value, counts

Activation copies these into a dataset database through the staged-write protocol: family
rows into ``family_annotations`` (and the run and status tables), CDD rows into
``annotations`` with ``source = 'cdd'`` (with a row in ``annotation_source_versions``), and
creates the ``family_resemblance`` view, which applies the "resembles" rule.
"""

from __future__ import annotations

import hashlib
import json
import time
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

import duckdb

from sharur.remote_homology import (
    CDD_MAX_EVALUE,
    RULE_MIN_ALIGNED_COLS,
    RULE_MIN_PROBABILITY,
    RULE_SOURCES,
    runner,
)


RUNS = "family_annotation_runs"
STATUS = "family_annotation_status"
HITS = "family_annotations"
CDD_HITS = "cdd_hits"
CDD_RUNS = "cdd_runs"
RESEMBLANCE_VIEW = "family_resemblance"
SOURCE_VERSIONS = "annotation_source_versions"
CDD_SOURCE = "cdd"
TABLES = (RUNS, STATUS, HITS)
# Tables an activation can change (diffed before publish).
ACTIVATION_TABLES = [RUNS, STATUS, HITS, "annotations", SOURCE_VERSIONS, "protein_families", "family_stats",
                     "protein_family_runs", "protein_family_assignments"]
SCHEMA_VERSION = 1
MAX_SEQUENCES = 250
MIN_GENOMES = 3

DDL = {
    RUNS: f"""CREATE TABLE IF NOT EXISTS {RUNS} (
    search_run VARCHAR PRIMARY KEY,
    clustering_run VARCHAR NOT NULL,
    created_at TIMESTAMP NOT NULL,
    selection JSON NOT NULL,
    msa_kinds VARCHAR[] NOT NULL,
    targets VARCHAR[] NOT NULL,
    databases JSON NOT NULL,
    parameters JSON NOT NULL,
    tool_versions JSON,
    families_selected BIGINT NOT NULL,
    families_ok BIGINT NOT NULL,
    families_error BIGINT NOT NULL,
    families_pending BIGINT NOT NULL,
    export_digest VARCHAR NOT NULL,
    config_digest VARCHAR NOT NULL,
    source_db VARCHAR,
    schema_version INTEGER NOT NULL,
    details JSON
)""",
    STATUS: f"""CREATE TABLE IF NOT EXISTS {STATUS} (
    search_run VARCHAR NOT NULL,
    clustering_run VARCHAR NOT NULL,
    family_id VARCHAR NOT NULL,
    status VARCHAR NOT NULL,
    n_members BIGINT,
    n_genomes BIGINT,
    n_unique BIGINT,
    n_sequences BIGINT,
    rep_protein_id VARCHAR,
    rep_length INTEGER,
    n_hits BIGINT,
    verified BOOLEAN,
    seconds DOUBLE,
    steps JSON,
    error VARCHAR,
    PRIMARY KEY (clustering_run, family_id)
)""",
    HITS: f"""CREATE TABLE IF NOT EXISTS {HITS} (
    family_id VARCHAR NOT NULL,
    clustering_run VARCHAR NOT NULL,
    search_run VARCHAR NOT NULL,
    source VARCHAR NOT NULL,
    db_name VARCHAR NOT NULL,
    db_version VARCHAR NOT NULL,
    msa_kind VARCHAR NOT NULL,
    hit_rank INTEGER NOT NULL,
    hit_accession VARCHAR NOT NULL,
    hit_name VARCHAR,
    hit_description VARCHAR,
    probability DOUBLE NOT NULL,
    evalue DOUBLE,
    score DOUBLE,
    aligned_cols INTEGER NOT NULL,
    identities DOUBLE,
    q_start INTEGER,
    q_end INTEGER,
    q_length INTEGER,
    t_start INTEGER,
    t_end INTEGER,
    t_length INTEGER,
    unknown_function BOOLEAN NOT NULL,
    hit_label VARCHAR,
    PRIMARY KEY (clustering_run, family_id, source, msa_kind, hit_rank)
)""",
    CDD_HITS: f"""CREATE TABLE IF NOT EXISTS {CDD_HITS} (
    protein_id VARCHAR NOT NULL,
    source VARCHAR NOT NULL,
    accession VARCHAR NOT NULL,
    name VARCHAR,
    description TEXT,
    evalue DOUBLE,
    score DOUBLE,
    start_aa INTEGER,
    end_aa INTEGER,
    sequence_md5 VARCHAR NOT NULL,
    pssm_id VARCHAR,
    model_length INTEGER,
    model_start INTEGER,
    model_end INTEGER
)""",
    CDD_RUNS: f"""CREATE TABLE IF NOT EXISTS {CDD_RUNS} (
    cdd_run VARCHAR PRIMARY KEY,
    created_at TIMESTAMP NOT NULL,
    db_name VARCHAR NOT NULL,
    db_version VARCHAR NOT NULL,
    db_prefix VARCHAR,
    max_evalue DOUBLE NOT NULL,
    tool_version VARCHAR,
    unique_sequences BIGINT NOT NULL,
    proteins BIGINT NOT NULL,
    hit_rows BIGINT NOT NULL,
    proteins_with_hits BIGINT NOT NULL,
    export_digest VARCHAR NOT NULL,
    source_db VARCHAR,
    details JSON
)""",
}

INDEXES = f"""
CREATE INDEX IF NOT EXISTS idx_family_annotations_family ON {HITS}(family_id);
CREATE INDEX IF NOT EXISTS idx_family_annotations_accession ON {HITS}(hit_accession);
"""


def resemblance_view_sql(min_probability: float = RULE_MIN_PROBABILITY,
                         min_cols: int = RULE_MIN_ALIGNED_COLS) -> str:
    """The ``family_resemblance`` view: family hits that meet the "resembles" rule."""
    sources = ", ".join(f"'{s}'" for s in RULE_SOURCES)
    return f"""CREATE OR REPLACE VIEW {RESEMBLANCE_VIEW} AS
SELECT family_id, clustering_run, search_run, source, db_name, db_version, msa_kind, hit_rank, hit_accession,
       hit_name, hit_description, probability, evalue, aligned_cols, q_start, q_end, q_length, t_start, t_end,
       t_length, unknown_function,
       CASE WHEN unknown_function THEN 'unknown_function' ELSE 'characterized' END AS hit_class,
       'resembles' AS evidence
FROM {HITS}
WHERE source IN ({sources}) AND probability >= {float(min_probability)} AND aligned_cols >= {int(min_cols)}"""


class RemoteHomologyError(RuntimeError):
    """A remote-homology precondition failed; nothing was written."""


def _q(path: Path | str) -> str:
    return "'" + str(path).replace("'", "''") + "'"


def connect(path: Path | str, *, read_only: bool = False, threads: int = 4,
            memory_limit: str = "8GB") -> duckdb.DuckDBPyConnection:
    con = duckdb.connect(str(path), read_only=read_only)
    con.execute(f"SET threads = {int(threads)}")
    con.execute(f"SET memory_limit = '{memory_limit}'")
    con.execute("SET preserve_insertion_order = false")
    return con


def _now() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def _sha(payload: Any) -> str:
    return hashlib.sha256(json.dumps(payload, sort_keys=True, default=str).encode()).hexdigest()


def _has_tables(con, names: Sequence[str], catalog: str | None = None) -> bool:
    rows = con.execute(
        "SELECT DISTINCT table_name FROM information_schema.tables WHERE table_schema = 'main' "
        + ("AND table_catalog = ? " if catalog else "")
        + "AND list_contains(?::VARCHAR[], table_name)", [*([catalog] if catalog else []), list(names)]).fetchall()
    return len(rows) == len(set(names))


# --------------------------------------------------------------------------- #
# Selection and export (family tier)
# --------------------------------------------------------------------------- #
FAMILY_TABLES = ("protein_families", "family_stats", "protein_family_runs")


def _family_catalog(con, db: Path, families_db: Path | None) -> str:
    """Attach the protein-family tables (from ``families_db`` or the dataset) and return their prefix."""
    if families_db is not None:
        con.execute(f"ATTACH {_q(Path(families_db).expanduser().resolve())} AS fam (READ_ONLY)")
        if not _has_tables(con, FAMILY_TABLES, "fam"):
            raise RemoteHomologyError(f"{families_db} holds no protein-family tables")
        return "fam."
    if not _has_tables(con, FAMILY_TABLES, "src"):
        raise RemoteHomologyError(f"{db} holds no protein-family tables; pass a families staging file "
                                  "(sharur families run/import) with --families")
    return "src."


def _resolve_run(con, fam: str, run: str | None) -> str:
    runs = {r[0]: r[1:] for r in con.execute(
        f"SELECT clustering_run, level, is_default FROM {fam}protein_family_runs").fetchall()}
    if run is None:
        defaults = sorted(r for r, (level, is_default) in runs.items() if level == "family" and is_default)
        if defaults:
            return defaults[0]
        if len(runs) == 1:
            return next(iter(runs))
        raise RemoteHomologyError(f"no default family run among {', '.join(sorted(runs))}; name one with --run")
    if run not in runs:
        raise RemoteHomologyError(f"no clustering run {run!r} (runs: {', '.join(sorted(runs))})")
    return run


def select_families(con, fam: str, run: str, *, min_genomes: int = MIN_GENOMES, require_pfam_dark: bool = True,
                    family_ids: Sequence[str] | None = None) -> None:
    """Create TEMP table ``_rh_sel``: the families of ``run`` in scope, ordered by ``family_id``.

    Scope: families found in ``min_genomes`` or more genomes whose members carry no Pfam hit
    (``family_stats.pfam_dark``; Sharur stores Pfam at gathering thresholds), or the families
    named in ``family_ids`` whatever their statistics.
    """
    where = ["s.clustering_run = ?"]
    params: list[Any] = [run]
    if family_ids:
        con.execute("CREATE OR REPLACE TEMP TABLE _rh_ids AS SELECT unnest(?::VARCHAR[]) AS family_id",
                    [list(family_ids)])
        where.append("s.family_id IN (SELECT family_id FROM _rh_ids)")
    else:
        where.append("s.n_genomes >= ?")
        params.append(int(min_genomes))
        if require_pfam_dark:
            where.append("s.pfam_dark")
    con.execute(f"""CREATE OR REPLACE TEMP TABLE _rh_sel AS
        SELECT row_number() OVER (ORDER BY s.family_id) - 1 AS family_index, s.family_id, s.n_members,
               s.n_genomes, s.rep_protein_id, s.rep_length
        FROM {fam}family_stats s WHERE {' AND '.join(where)}""", params)
    if family_ids:
        found = con.execute("SELECT COUNT(*) FROM _rh_sel").fetchone()[0]
        if found != len(set(family_ids)):
            raise RemoteHomologyError(f"{len(set(family_ids)) - found} named families are not in run {run}")


def family_key(index: int) -> str:
    return f"f{index:07d}"


def export_families(
    db: Path | str,
    out_dir: Path | str,
    *,
    families_db: Path | str | None = None,
    run: str | None = None,
    min_genomes: int = MIN_GENOMES,
    require_pfam_dark: bool = True,
    max_sequences: int = MAX_SEQUENCES,
    family_ids: Sequence[str] | None = None,
    threads: int = 4,
    memory_limit: str = "8GB",
    log: Callable[[str], None] = print,
) -> dict[str, Any]:
    """Write the family-tier inputs of ``run`` to ``out_dir`` (``db`` and ``families_db`` are only read).

    Each selected family gets one FASTA of its unique member sequences: the representative's
    first, then up to ``max_sequences - 1`` others in MD5 order of their sequence (a fixed,
    content-defined sample). Headers are protein IDs; a sequence shared by several members is
    written once under the lowest of their IDs (the representative's ID for its own sequence).
    """
    db = Path(db).expanduser().resolve()
    out = Path(out_dir).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    started = time.time()
    con = connect(":memory:", threads=threads, memory_limit=memory_limit)
    try:
        con.execute(f"ATTACH {_q(db)} AS src (READ_ONLY)")
        fam = _family_catalog(con, db, Path(families_db) if families_db else None)
        run = _resolve_run(con, fam, run)
        select_families(con, fam, run, min_genomes=min_genomes, require_pfam_dark=require_pfam_dark,
                        family_ids=family_ids)
        n_sel = con.execute("SELECT COUNT(*) FROM _rh_sel").fetchone()[0]
        if not n_sel:
            raise RemoteHomologyError(f"run {run}: no family meets the selection")
        con.execute(f"""CREATE OR REPLACE TEMP TABLE _rh_uniq AS
            WITH m AS (
                SELECT s.family_index, f.protein_id, f.is_representative, rtrim(p.sequence, '*') AS seq
                FROM _rh_sel s JOIN {fam}protein_families f ON f.family_id = s.family_id AND f.clustering_run = ?
                JOIN src.proteins p ON p.protein_id = f.protein_id
                WHERE rtrim(COALESCE(p.sequence, ''), '*') <> ''),
            u AS (
                SELECT family_index, seq, bool_or(is_representative) AS has_rep,
                       COALESCE(min(protein_id) FILTER (WHERE is_representative), min(protein_id)) AS protein_id
                FROM m GROUP BY family_index, seq)
            SELECT family_index, protein_id, seq,
                   row_number() OVER (PARTITION BY family_index ORDER BY has_rep DESC, md5(seq), protein_id) AS k,
                   count(*) OVER (PARTITION BY family_index) AS n_unique
            FROM u""", [run])
        missing_rep = con.execute("""
            SELECT COUNT(*) FROM _rh_sel s
            LEFT JOIN _rh_uniq u ON u.family_index = s.family_index AND u.k = 1 AND u.protein_id = s.rep_protein_id
            WHERE u.protein_id IS NULL""").fetchone()[0]
        if missing_rep:
            raise RemoteHomologyError(f"{missing_rep} families lack their representative's sequence in {db}")
        digest = hashlib.sha256()
        rows = []
        cur = con.execute("""SELECT s.family_index, s.family_id, s.n_members, s.n_genomes, s.rep_protein_id,
                                    s.rep_length, u.n_unique, u.protein_id, u.seq, u.k
                             FROM _rh_sel s JOIN _rh_uniq u USING (family_index)
                             WHERE u.k <= ? ORDER BY s.family_index, u.k""", [int(max_sequences)])
        current, handle, record = None, None, None
        while True:
            batch = cur.fetchmany(50_000)
            if not batch:
                break
            for index, fid, n_members, n_genomes, rep, _rep_length, n_unique, pid, seq, _k in batch:
                if index != current:
                    if handle:
                        handle.close()
                        rows.append(record)
                    current = index
                    key = family_key(index)
                    rel = Path("fasta") / f"{index // 1000:03d}" / f"{key}.faa"
                    (out / rel).parent.mkdir(parents=True, exist_ok=True)
                    handle = open(out / rel, "w")
                    # rep_length is the residue count written (proteins.sequence_length may count a stop)
                    record = {"family_index": index, "family_key": key, "family_id": fid, "n_members": n_members,
                              "n_genomes": n_genomes, "n_unique": n_unique, "n_sequences": 0,
                              "rep_protein_id": rep, "rep_length": len(seq), "fasta": str(rel)}
                text = f">{pid}\n{seq}\n"
                handle.write(text)
                digest.update(text.encode())
                record["n_sequences"] += 1
        if handle:
            handle.close()
            rows.append(record)
        with open(out / "families.tsv", "w") as tsv:
            tsv.write("\t".join(runner.FAMILY_TSV_COLUMNS) + "\n")
            for r in rows:
                tsv.write("\t".join(str(r[c]) for c in runner.FAMILY_TSV_COLUMNS) + "\n")
        run_row = con.execute(f"SELECT level, input_digest, parameters FROM {fam}protein_family_runs "
                              "WHERE clustering_run = ?", [run]).fetchone()
        manifest = {
            "kind": "family_tier", "clustering_run": run, "level": run_row[0], "clustering_input_digest": run_row[1],
            "selection": {"min_genomes": min_genomes, "require_pfam_dark": require_pfam_dark,
                          "max_sequences": max_sequences, "family_ids": sorted(family_ids) if family_ids else None},
            "families": len(rows), "sequences": sum(r["n_sequences"] for r in rows),
            "member_proteins": sum(r["n_members"] for r in rows), "export_digest": digest.hexdigest(),
            "source_db": str(db), "families_db": str(Path(families_db).resolve()) if families_db else None,
            "created_at": time.strftime("%Y-%m-%dT%H:%M:%S"), "seconds": round(time.time() - started, 1)}
        (out / "export.json").write_text(json.dumps(manifest, indent=2) + "\n")
    finally:
        con.close()
    log(f"exported {manifest['families']:,} families ({manifest['sequences']:,} sequences) of run {run}")
    return manifest


# --------------------------------------------------------------------------- #
# Export (CDD)
# --------------------------------------------------------------------------- #
def export_cdd(db: Path | str, out_dir: Path | str, *, chunk_size: int = 2000, threads: int = 4,
               memory_limit: str = "8GB", log: Callable[[str], None] = print) -> dict[str, Any]:
    """Write every unique protein sequence of ``db`` (MD5 headers, MD5 order) in chunks for RPS-BLAST."""
    db = Path(db).expanduser().resolve()
    out = Path(out_dir).expanduser().resolve() / "cdd"
    chunks = out / "chunks"
    chunks.mkdir(parents=True, exist_ok=True)
    for stale in chunks.glob("chunk_*.faa"):
        stale.unlink()
    started = time.time()
    con = connect(db, read_only=True, threads=threads, memory_limit=memory_limit)
    try:
        n_proteins = con.execute("SELECT COUNT(*) FROM proteins WHERE rtrim(COALESCE(sequence, ''), '*') <> ''"
                                 ).fetchone()[0]
        cur = con.execute("""SELECT md5(seq) AS h, any_value(seq) FROM (
                                 SELECT rtrim(sequence, '*') AS seq FROM proteins
                                 WHERE rtrim(COALESCE(sequence, ''), '*') <> '')
                             GROUP BY h ORDER BY h""")
        digest = hashlib.sha256()
        n = n_chunks = 0
        handle = None
        while True:
            batch = cur.fetchmany(50_000)
            if not batch:
                break
            for h, seq in batch:
                if n % chunk_size == 0:
                    if handle:
                        handle.close()
                    handle = open(chunks / f"chunk_{n_chunks:05d}.faa", "w")
                    n_chunks += 1
                text = f">{h}\n{seq}\n"
                handle.write(text)
                digest.update(text.encode())
                n += 1
        if handle:
            handle.close()
    finally:
        con.close()
    manifest = {"kind": "cdd", "proteins": n_proteins, "unique_sequences": n, "chunks": n_chunks,
                "chunk_size": chunk_size, "export_digest": digest.hexdigest(), "source_db": str(db),
                "created_at": time.strftime("%Y-%m-%dT%H:%M:%S"), "seconds": round(time.time() - started, 1)}
    (out / "export.json").write_text(json.dumps(manifest, indent=2) + "\n")
    log(f"exported {n:,} unique sequences of {n_proteins:,} proteins in {n_chunks:,} chunks")
    return manifest


# --------------------------------------------------------------------------- #
# Loading search results into a staging DuckDB
# --------------------------------------------------------------------------- #
_SOURCE_DB_KEY = {"hhsearch_pfam": "pfam", "hhsearch_pdb": "pdb", "hhblits_uniref": "uniref"}


def _family_status(work: Path, family: dict, digest: str) -> dict:
    d = runner.family_dir(work, family)
    for name in ("done.json", "error.json"):
        path = d / name
        if path.exists():
            try:
                record = json.loads(path.read_text())
            except json.JSONDecodeError:
                continue
            if name == "done.json" and record.get("config_digest") != digest:
                return {"status": "stale", "record": record}
            return {"status": record.get("status", "error"), "record": record}
    return {"status": "pending", "record": {}}


def load_families(
    work: Path | str,
    staging: Path | str,
    config: Path | str,
    *,
    allow_partial: bool = False,
    threads: int = 4,
    memory_limit: str = "8GB",
    log: Callable[[str], None] = print,
) -> dict[str, Any]:
    """Load the family-tier results of ``work`` into ``staging`` (created when absent).

    Every selected family gets a ``family_annotation_status`` row (``ok``, ``error`` or
    ``pending``); hits of ``ok`` families go to ``family_annotations`` with the name and
    version of their database from ``config``. Rows of the same clustering run in the staging
    file are replaced. Unless ``allow_partial``, every family must be ``ok``.
    """
    work = Path(work).expanduser().resolve()
    cfg = runner.load_config(config)
    digest = runner.family_digest(cfg)
    manifest = json.loads((work / "export.json").read_text())
    families = runner.read_families(work)
    run = manifest["clustering_run"]
    search_run = "fhs-" + _sha({"clustering_run": run, "config": digest,
                                "export": manifest["export_digest"]})[:10]
    status_rows, counts = [], {"ok": 0, "error": 0, "pending": 0, "stale": 0}
    versions_seen: dict[str, Any] = {}
    combined = work / "hits_combined.tsv"
    n_hits = 0
    with open(combined, "w") as out:
        out.write("family_id\t" + "\t".join(runner.hhr.HIT_COLUMNS) + "\n")
        for fam in families:
            st = _family_status(work, fam, digest)
            rec = st["record"]
            status = st["status"]
            if status == "ok":
                hits_path = runner.family_dir(work, fam) / "hits.tsv"
                if not hits_path.exists():
                    status = "pending"
                else:
                    with open(hits_path) as handle:
                        next(handle)
                        for line in handle:
                            out.write(f"{fam['family_id']}\t{line}")
                            n_hits += 1
            counts[status if status in counts else "error"] += 1
            if status == "stale":
                status = "pending"
            status_rows.append((search_run, run, fam["family_id"], status, fam["n_members"], fam["n_genomes"],
                                fam["n_unique"], fam["n_sequences"], fam["rep_protein_id"], fam["rep_length"],
                                rec.get("n_hits") if status == "ok" else None,
                                rec.get("verified") if status == "ok" else None, rec.get("seconds"),
                                json.dumps(rec.get("steps")) if rec.get("steps") else None,
                                rec.get("error") or rec.get("warning")))
    pending = counts["pending"] + counts["stale"]
    if (pending or counts["error"]) and not allow_partial:
        combined.unlink(missing_ok=True)
        raise RemoteHomologyError(f"{counts['ok']:,} of {len(families):,} families are done ({counts['error']:,} "
                                  f"errors, {pending:,} pending); finish the run or pass allow_partial")
    for summary in sorted((work / "logs").glob("families_shard*.json")):
        try:
            versions_seen = json.loads(summary.read_text()).get("versions") or versions_seen
        except json.JSONDecodeError:
            continue
    db_cols = {}
    for source, key in _SOURCE_DB_KEY.items():
        db = cfg["databases"].get(key)
        if db:
            db_cols[source] = (db["name"], db["version"])
    databases = {k: {f: v.get(f) for f in ("prefix", "name", "version", "kind", "titles", "built_from")}
                 for k, v in cfg["databases"].items() if k in (*cfg["family_targets"], "uniref")}
    staging = Path(staging).expanduser().resolve()
    staging.parent.mkdir(parents=True, exist_ok=True)
    con = connect(staging, threads=threads, memory_limit=memory_limit)
    try:
        con.execute("BEGIN TRANSACTION")
        for ddl in (DDL[RUNS], DDL[STATUS], DDL[HITS]):
            con.execute(ddl)
        for table in (HITS, STATUS, RUNS):
            con.execute(f"DELETE FROM {table} WHERE clustering_run = ?", [run])
        import pandas as pd  # noqa: PLC0415

        status_df = pd.DataFrame(status_rows, columns=[
            "search_run", "clustering_run", "family_id", "status", "n_members", "n_genomes", "n_unique",
            "n_sequences", "rep_protein_id", "rep_length", "n_hits", "verified", "seconds", "steps", "error"])
        status_df = status_df.astype({"n_hits": "Int64", "verified": "boolean", "seconds": "Float64"})
        con.register("_rh_status_df", status_df)
        con.execute(f"INSERT INTO {STATUS} SELECT * FROM _rh_status_df")
        con.unregister("_rh_status_df")
        con.execute("CREATE OR REPLACE TEMP TABLE _rh_db AS SELECT * FROM (VALUES "
                    + ", ".join("(?, ?, ?)" for _ in db_cols) + ") t(source, db_name, db_version)",
                    [x for s, (n, v) in db_cols.items() for x in (s, n, v)])
        con.execute(f"""
            INSERT INTO {HITS}
            SELECT h.family_id, ?, ?, h.source, d.db_name, d.db_version, h.msa_kind, h.hit_rank, h.hit_accession,
                   NULLIF(h.hit_name, ''), NULLIF(h.hit_description, ''), h.probability, h.evalue, h.score,
                   h.aligned_cols, h.identities, h.q_start, h.q_end, s.rep_length, h.t_start, h.t_end, h.t_length,
                   h.unknown_function, h.hit_label
            FROM read_csv({_q(combined)}, delim = '\t', header = true, quote = '', escape = '',
                          all_varchar = false, columns = {{
                              'family_id': 'VARCHAR', 'source': 'VARCHAR', 'msa_kind': 'VARCHAR',
                              'hit_rank': 'INTEGER', 'hit_accession': 'VARCHAR', 'hit_name': 'VARCHAR',
                              'hit_description': 'VARCHAR', 'probability': 'DOUBLE', 'evalue': 'DOUBLE',
                              'score': 'DOUBLE', 'aligned_cols': 'INTEGER', 'identities': 'DOUBLE',
                              'q_start': 'INTEGER', 'q_end': 'INTEGER', 't_start': 'INTEGER', 't_end': 'INTEGER',
                              't_length': 'INTEGER', 'unknown_function': 'BOOLEAN', 'hit_label': 'VARCHAR'}}) h
            JOIN _rh_db d USING (source)
            JOIN {STATUS} s ON s.family_id = h.family_id AND s.clustering_run = ?""", [run, search_run, run])
        loaded = con.execute(f"SELECT COUNT(*) FROM {HITS} WHERE clustering_run = ?", [run]).fetchone()[0]
        if loaded != n_hits:
            raise RemoteHomologyError(f"{n_hits:,} hit rows read but {loaded:,} loaded (a source lacks a database "
                                      "entry in the configuration)")
        con.execute(f"INSERT INTO {RUNS} VALUES (?, ?, ?::TIMESTAMP, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                    [search_run, run, _now(), json.dumps(manifest["selection"]), cfg["msa_kinds"],
                     cfg["family_targets"], json.dumps(databases), json.dumps(cfg["search"]),
                     json.dumps(versions_seen), len(families), counts["ok"], counts["error"], pending,
                     manifest["export_digest"], digest, manifest.get("source_db"), SCHEMA_VERSION,
                     json.dumps({"clustering_input_digest": manifest.get("clustering_input_digest"),
                                 "families_db": manifest.get("families_db"), "work_dir": str(work),
                                 "config_path": cfg["_path"], "sequences": manifest.get("sequences"),
                                 "member_proteins": manifest.get("member_proteins")})])
        con.execute(INDEXES)
        con.execute(resemblance_view_sql())
        con.execute("COMMIT")
    except Exception:
        con.execute("ROLLBACK")
        raise
    finally:
        con.close()
        combined.unlink(missing_ok=True)
    result = {"search_run": search_run, "clustering_run": run, "families": len(families), **counts,
              "hit_rows": n_hits, "staging": str(staging)}
    log(f"{search_run}: {counts['ok']:,} of {len(families):,} families loaded, {n_hits:,} hit rows")
    return result


_CDD_COLUMNS = ", ".join(f"'{c}': '{t}'" for c, t in zip(
    runner.CDD_COLUMNS, ("VARCHAR", "VARCHAR", "INTEGER", "INTEGER", "INTEGER", "INTEGER", "INTEGER", "INTEGER",
                         "DOUBLE", "DOUBLE", "VARCHAR")))


def load_cdd(
    work: Path | str,
    staging: Path | str,
    db: Path | str,
    config: Path | str,
    *,
    max_evalue: float = CDD_MAX_EVALUE,
    allow_partial: bool = False,
    threads: int = 4,
    memory_limit: str = "8GB",
    log: Callable[[str], None] = print,
) -> dict[str, Any]:
    """Load the RPS-BLAST results of ``work/cdd`` into ``staging`` as ``cdd_hits`` (protein level).

    Hits at ``max_evalue`` or below are expanded from unique sequences (MD5) to every protein of
    ``db`` (opened read-only) with that sequence. CDD titles ``"cd00009, AAA, description"``
    give accession, name and description. ``cdd_hits`` in the staging file is replaced.
    """
    work = Path(work).expanduser().resolve()
    cfg = runner.load_config(config)
    cdd = runner.database(cfg, "cdd")
    manifest = json.loads((work / "cdd" / "export.json").read_text())
    chunks = sorted((work / "cdd" / "chunks").glob("chunk_*.faa"))
    out = work / "cdd" / "out"
    done = [c for c in chunks if (out / f"{c.stem}.done").exists() and (out / f"{c.stem}.tsv").exists()]
    if len(done) != len(chunks) and not allow_partial:
        raise RemoteHomologyError(f"{len(done):,} of {len(chunks):,} CDD chunks are done; finish the run first")
    if not done:
        raise RemoteHomologyError("no finished CDD chunk")
    records = [json.loads((out / f"{c.stem}.done").read_text()) for c in done]
    evalues = {r.get("evalue") for r in records}
    if None in evalues:
        raise RemoteHomologyError("CDD chunk records lack the E-value of their search")
    if min(float(e) for e in evalues) < max_evalue:
        raise RemoteHomologyError(f"the search kept hits at E <= {min(evalues)}, stricter than {max_evalue}; "
                                  "rerun it with the configured cdd_evalue")
    db = Path(db).expanduser().resolve()
    staging = Path(staging).expanduser().resolve()
    staging.parent.mkdir(parents=True, exist_ok=True)
    # chunks without a single hit leave empty tables, which the CSV reader cannot sniff
    files = [str(out / f"{c.stem}.tsv") for c in done if (out / f"{c.stem}.tsv").stat().st_size > 0]
    con = connect(staging, threads=threads, memory_limit=memory_limit)
    try:
        con.execute(f"ATTACH {_q(db)} AS src (READ_ONLY)")
        if files:
            con.execute(f"""CREATE OR REPLACE TEMP TABLE _rh_raw AS
                SELECT * FROM read_csv(?, delim = '\t', header = false, quote = '', escape = '',
                                       columns = {{{_CDD_COLUMNS}}})""", [files])
        else:
            con.execute("CREATE OR REPLACE TEMP TABLE _rh_raw (" + ", ".join(
                f"{c} {t}" for c, t in zip(runner.CDD_COLUMNS, ("VARCHAR", "VARCHAR", "INTEGER", "INTEGER", "INTEGER",
                                                               "INTEGER", "INTEGER", "INTEGER", "DOUBLE", "DOUBLE",
                                                               "VARCHAR"))) + ")")
        raw_rows = con.execute("SELECT COUNT(*) FROM _rh_raw").fetchone()[0]
        con.execute("""CREATE OR REPLACE TEMP TABLE _rh_md5 AS
            SELECT protein_id, md5(rtrim(sequence, '*')) AS h FROM src.proteins
            WHERE rtrim(COALESCE(sequence, ''), '*') <> ''""")
        unknown = con.execute("SELECT COUNT(DISTINCT qseqid) FROM _rh_raw ANTI JOIN _rh_md5 m ON m.h = qseqid"
                              ).fetchone()[0]
        if unknown:
            raise RemoteHomologyError(f"{unknown:,} searched sequences match no protein of {db}; the export "
                                      "belongs to another protein set")
        con.execute("BEGIN TRANSACTION")
        con.execute(DDL[CDD_HITS])
        con.execute(DDL[CDD_RUNS])
        con.execute(f"DELETE FROM {CDD_HITS}")
        con.execute(f"DELETE FROM {CDD_RUNS}")
        con.execute(f"""
            INSERT INTO {CDD_HITS}
            SELECT m.protein_id, '{CDD_SOURCE}',
                   COALESCE(NULLIF(trim(split_part(r.stitle, ',', 1)), ''), split_part(r.sseqid, '|', 3)),
                   NULLIF(trim(split_part(r.stitle, ',', 2)), ''),
                   NULLIF(trim(regexp_extract(r.stitle, '^[^,]*,[^,]*,(.*)$', 1)), ''),
                   r.evalue, r.bitscore, r.qstart, r.qend, r.qseqid, split_part(r.sseqid, '|', 3),
                   r.slen, r.sstart, r.send
            FROM _rh_raw r JOIN _rh_md5 m ON m.h = r.qseqid
            WHERE r.evalue <= ?""", [float(max_evalue)])
        hit_rows, proteins_with = con.execute(
            f"SELECT COUNT(*), COUNT(DISTINCT protein_id) FROM {CDD_HITS}").fetchone()
        n_proteins = con.execute("SELECT COUNT(*) FROM _rh_md5").fetchone()[0]
        versions = {r.get("rpsblast") for r in
                    (json.loads(p.read_text()).get("versions", {}) for p in (work / "logs").glob("cdd_shard*.json"))
                    if r.get("rpsblast")}
        cdd_run = "cdd-" + _sha({"db": [cdd["name"], cdd["version"]], "evalue": max_evalue,
                                 "export": manifest["export_digest"]})[:10]
        con.execute(f"INSERT INTO {CDD_RUNS} VALUES (?, ?::TIMESTAMP, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                    [cdd_run, _now(), cdd["name"], cdd["version"], cdd["prefix"], float(max_evalue),
                     "; ".join(sorted(versions)) or None, manifest["unique_sequences"], n_proteins, hit_rows,
                     proteins_with, manifest["export_digest"], str(db),
                     json.dumps({"chunks": len(chunks), "chunks_done": len(done), "raw_rows": raw_rows,
                                 "work_dir": str(work), "config_path": cfg["_path"],
                                 "search_evalue": sorted(e for e in evalues if e is not None)})])
        con.execute("COMMIT")
    except Exception:
        try:
            con.execute("ROLLBACK")
        except duckdb.Error:
            pass
        raise
    finally:
        con.close()
    result = {"cdd_run": cdd_run, "chunks": len(chunks), "chunks_done": len(done), "raw_rows": raw_rows,
              "hit_rows": hit_rows, "proteins_with_hits": proteins_with, "proteins": n_proteins,
              "staging": str(staging)}
    log(f"{cdd_run}: {hit_rows:,} CDD rows on {proteins_with:,} of {n_proteins:,} proteins")
    return result


# --------------------------------------------------------------------------- #
# Activation (copy staged rows into a dataset database)
# --------------------------------------------------------------------------- #
def _attach_staging(con, path: Path, alias: str) -> set[str]:
    con.execute(f"ATTACH {_q(path)} AS {alias} (READ_ONLY)")
    return {r[0] for r in con.execute(
        "SELECT table_name FROM information_schema.tables WHERE table_catalog = ? AND table_schema = 'main'",
        [alias]).fetchall()}


def load_staging(con: duckdb.DuckDBPyConnection, staging: Sequence[Path | str], *,
                 families_staging: Path | str | None = None, verify_digest: bool = True,
                 log: Callable[[str], None] = print) -> dict[str, Any]:
    """Copy the remote-homology rows of ``staging`` file(s) into the database open on ``con``.

    Family tier: the clustering run of every staged search must exist in the database (or in
    ``families_staging``, whose runs are loaded first) with the clustering input digest the
    search recorded, and every staged family must exist in it; rows of the same clustering run
    and source are replaced. CDD: every staged protein must exist with the sequence the hit was
    computed for (MD5); ``annotations`` rows with ``source = 'cdd'`` are replaced and the
    database name and version are recorded in ``annotation_source_versions``.
    """
    from sharur.defensefinder.selection import _SOURCE_VERSIONS_SQL  # noqa: PLC0415

    summary: dict[str, Any] = {"family_tier": None, "cdd": None}
    if families_staging is not None:
        from sharur.protein_families import store as pf  # noqa: PLC0415

        needed: set[str] = set()
        for path in staging:
            con.execute(f"ATTACH {_q(Path(path).expanduser().resolve())} AS _rh_peek (READ_ONLY)")
            try:
                if _has_tables(con, [RUNS], "_rh_peek"):
                    needed |= {r[0] for r in con.execute(f"SELECT clustering_run FROM _rh_peek.{RUNS}").fetchall()}
            finally:
                con.execute("DETACH _rh_peek")
        pf.load_staging(con, families_staging, runs=sorted(needed) or None, verify_digest=verify_digest, log=log)
        summary["families_loaded"] = sorted(needed)
    for i, path in enumerate(staging):
        path = Path(path).expanduser().resolve()
        alias = f"rh{i}"
        tables = _attach_staging(con, path, alias)
        try:
            if {RUNS, STATUS, HITS} <= tables:
                summary["family_tier"] = _load_family_tier(con, alias, path)
                log(f"family tier: {summary['family_tier']}")
            if CDD_HITS in tables:
                summary["cdd"] = _load_cdd_rows(con, alias, path, verify_digest, _SOURCE_VERSIONS_SQL)
                log(f"cdd: {summary['cdd']}")
            if not ({RUNS, STATUS, HITS} <= tables or CDD_HITS in tables):
                raise RemoteHomologyError(f"{path} holds no remote-homology tables")
        finally:
            con.execute(f"DETACH {alias}")
    # the SQL surface's resembles macros read these tables; create them empty when a layer is absent
    from sharur.protein_families.store import has_families  # noqa: PLC0415

    con.execute(_SOURCE_VERSIONS_SQL)
    if has_families(con) and not _has_tables(con, [HITS]):
        for ddl in (DDL[RUNS], DDL[STATUS], DDL[HITS]):
            con.execute(ddl)
        con.execute(resemblance_view_sql())
    summary["status"] = status(con)
    return summary


def _load_family_tier(con, alias: str, path: Path) -> dict[str, Any]:
    from sharur.protein_families.store import has_families  # noqa: PLC0415

    if not has_families(con):
        raise RemoteHomologyError("the database holds no protein-family tables; pass the families staging "
                                  "file (--families-staging) or activate the clustering run first")
    runs = con.execute(f"SELECT search_run, clustering_run, details FROM {alias}.{RUNS}").fetchall()
    for search_run, run, details in runs:
        row = con.execute("SELECT input_digest FROM protein_family_runs WHERE clustering_run = ?", [run]).fetchone()
        if row is None:
            raise RemoteHomologyError(f"search {search_run} ran on clustering run {run}, which the database lacks")
        recorded = json.loads(details or "{}").get("clustering_input_digest")
        if recorded and recorded != row[0]:
            raise RemoteHomologyError(f"clustering run {run} in the database differs from the one searched "
                                      "(input digest); re-export and re-run the search")
    unknown = con.execute(f"""SELECT COUNT(*) FROM {alias}.{STATUS} s
                              ANTI JOIN family_stats f ON f.family_id = s.family_id
                                                      AND f.clustering_run = s.clustering_run""").fetchone()[0]
    if unknown:
        raise RemoteHomologyError(f"{unknown:,} searched families are absent from family_stats")
    con.execute("BEGIN TRANSACTION")
    try:
        con.execute(f"DROP VIEW IF EXISTS {RESEMBLANCE_VIEW}")
        for ddl in (DDL[RUNS], DDL[STATUS], DDL[HITS]):
            con.execute(ddl)
        con.execute(f"""CREATE OR REPLACE TEMP TABLE _rh_keys AS
                        SELECT DISTINCT clustering_run, source FROM {alias}.{HITS}""")
        con.execute(f"DELETE FROM {HITS} WHERE (clustering_run, source) IN (SELECT clustering_run, source FROM _rh_keys)")
        staged_runs = [r[1] for r in runs]
        con.execute(f"DELETE FROM {STATUS} WHERE list_contains(?::VARCHAR[], clustering_run)", [staged_runs])
        con.execute(f"DELETE FROM {RUNS} WHERE list_contains(?::VARCHAR[], clustering_run)", [staged_runs])
        for table in (RUNS, STATUS, HITS):
            con.execute(f"INSERT INTO {table} SELECT * FROM {alias}.{table}")
        con.execute(INDEXES)
        con.execute(resemblance_view_sql())
        con.execute("COMMIT")
    except Exception:
        con.execute("ROLLBACK")
        raise
    n_hits, n_fam = con.execute(f"SELECT COUNT(*), COUNT(DISTINCT family_id) FROM {alias}.{HITS}").fetchone()
    return {"staging": str(path), "searches": [r[0] for r in runs], "clustering_runs": sorted(set(staged_runs)),
            "hit_rows": n_hits, "families_with_hits": n_fam}


def _load_cdd_rows(con, alias: str, path: Path, verify_digest: bool, source_versions_sql: str) -> dict[str, Any]:
    run = con.execute(f"SELECT * FROM {alias}.{CDD_RUNS}").fetchdf()
    if len(run) != 1:
        raise RemoteHomologyError(f"{path}: expected one CDD run record, found {len(run)}")
    run = run.iloc[0]
    missing = con.execute(f"SELECT COUNT(DISTINCT protein_id) FROM {alias}.{CDD_HITS} ANTI JOIN proteins USING (protein_id)"
                          ).fetchone()[0]
    if missing:
        raise RemoteHomologyError(f"{missing:,} proteins with CDD hits are absent from the database")
    if verify_digest:
        changed = con.execute(f"""
            SELECT COUNT(*) FROM (SELECT DISTINCT protein_id, sequence_md5 FROM {alias}.{CDD_HITS}) h
            JOIN proteins p USING (protein_id) WHERE md5(rtrim(p.sequence, '*')) <> h.sequence_md5""").fetchone()[0]
        if changed:
            raise RemoteHomologyError(f"{changed:,} proteins changed sequence since the CDD search; rerun it")
    con.execute("BEGIN TRANSACTION")
    try:
        con.execute(source_versions_sql)
        removed = con.execute("SELECT COUNT(*) FROM annotations WHERE source = ?", [CDD_SOURCE]).fetchone()[0]
        con.execute("DELETE FROM annotations WHERE source = ?", [CDD_SOURCE])
        base = con.execute("SELECT COALESCE(MAX(annotation_id), 0) FROM annotations").fetchone()[0]
        con.execute(f"""
            INSERT INTO annotations (annotation_id, protein_id, source, accession, name, description, evalue,
                                     score, start_aa, end_aa)
            SELECT ? + row_number() OVER (ORDER BY protein_id, start_aa, accession, evalue), protein_id, source,
                   accession, name, description, evalue, score, start_aa, end_aa
            FROM {alias}.{CDD_HITS}""", [int(base)])
        loaded = con.execute("SELECT COUNT(*) FROM annotations WHERE source = ?", [CDD_SOURCE]).fetchone()[0]
        params = {"cdd_run": run["cdd_run"], "max_evalue": float(run["max_evalue"]), "db_prefix": run["db_prefix"],
                  "tool_version": run["tool_version"], "unique_sequences": int(run["unique_sequences"]),
                  "export_digest": run["export_digest"], "evidence": "resembles"}
        con.execute(f"DELETE FROM {SOURCE_VERSIONS} WHERE source = ?", [CDD_SOURCE])
        con.execute(f"""INSERT INTO {SOURCE_VERSIONS} (source, model_set, model_version, profiles_path, selection_rule,
                                    parameters, hits_path, raw_rows, selected_rows, loaded_rows, loaded_at)
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?::TIMESTAMP)""",
                    [CDD_SOURCE, run["db_name"], run["db_version"], run["db_prefix"],
                     f"rpsblast E <= {float(run['max_evalue']):g}, every HSP, unique sequences expanded to proteins",
                     json.dumps(params, sort_keys=True), str(path),
                     int(json.loads(run["details"] or "{}").get("raw_rows") or 0), int(run["hit_rows"]), loaded,
                     _now()])
        con.execute("COMMIT")
    except Exception:
        con.execute("ROLLBACK")
        raise
    return {"staging": str(path), "rows_removed": removed, "rows_loaded": loaded,
            "proteins_with_hits": int(run["proteins_with_hits"]), "db": f"{run['db_name']} {run['db_version']}"}


def status(con) -> dict[str, Any]:
    """Searches, coverage and rule counts of the remote-homology layer of the database on ``con``."""
    out: dict[str, Any] = {"family_tier": None, "cdd": None}
    if _has_tables(con, [RUNS, STATUS, HITS]):
        searches = []
        for row in con.execute(f"""SELECT search_run, clustering_run, created_at, msa_kinds, targets, databases,
                                          families_selected, families_ok, families_error, families_pending
                                   FROM {RUNS} ORDER BY created_at""").fetchall():
            searches.append(dict(zip(["search_run", "clustering_run", "created_at", "msa_kinds", "targets",
                                      "databases", "families_selected", "families_ok", "families_error",
                                      "families_pending"], [*row[:5], json.loads(row[5]), *row[6:]])))
        by_run: dict[str, dict[str, int]] = {}
        if _has_tables(con, [RESEMBLANCE_VIEW]):
            for run, characterized, unknown_only in con.execute(f"""
                    SELECT clustering_run, COUNT(*) FILTER (WHERE has_char),
                           COUNT(*) FILTER (WHERE NOT has_char)
                    FROM (SELECT clustering_run, family_id, bool_or(NOT unknown_function) AS has_char
                          FROM {RESEMBLANCE_VIEW} GROUP BY ALL) GROUP BY 1 ORDER BY 1""").fetchall():
                by_run[run] = {"characterized": characterized, "unknown_function_only": unknown_only}
        out["family_tier"] = {
            "searches": searches,
            "families_resembling_characterized": sum(v["characterized"] for v in by_run.values()),
            "families_resembling_unknown_function_only": sum(v["unknown_function_only"] for v in by_run.values()),
            "families_resembling_any": sum(v["characterized"] + v["unknown_function_only"] for v in by_run.values()),
            "by_clustering_run": by_run}
    if _has_tables(con, ["annotations"]):
        n, proteins = con.execute("SELECT COUNT(*), COUNT(DISTINCT protein_id) FROM annotations WHERE source = ?",
                                  [CDD_SOURCE]).fetchone()
        if n:
            version = None
            if _has_tables(con, [SOURCE_VERSIONS]):
                row = con.execute(f"SELECT model_set, model_version FROM {SOURCE_VERSIONS} WHERE source = ?",
                                  [CDD_SOURCE]).fetchone()
                version = f"{row[0]} {row[1]}" if row else None
            out["cdd"] = {"rows": n, "proteins": proteins, "database": version}
    return out
