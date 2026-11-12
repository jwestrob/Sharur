"""Remote-homology layer in DuckDB: input export, result loading, staging, activation.

Work directory (``export``)::

    WORK/export.json                 selection rule, clustering run, counts, sequence digest
    WORK/families.tsv                one row per selected family (family_index, family_key, ...)
    WORK/fasta/<bucket>/<key>.faa    up to ``max_sequences`` unique member sequences, representative first
    WORK/out/<bucket>/<key>/         runner.py output: .hhr(.gz) files, hits.tsv, done.json
    WORK/cdd/chunks/chunk_N.faa      unique sequences of every protein, named by MD5 (``export-cdd``)
    WORK/cdd/out/chunk_N.tsv         runner.py cdd output (raw RPS-BLAST tables; stay on the cluster)
    WORK/cdd/compact/                hits.parquet, models.parquet, compact.json (``compact-cdd``)

Staging DuckDB (``load-families``, ``load-cdd``)::

    family_annotation_runs     one row per family-tier search: clustering run, selection rule,
                               databases with names and versions, parameters, tool versions, counts
    family_annotation_status   one row per selected family: ok / error / pending, timings
    family_annotations         one row per hit at probability >= 50% (the first 20 for UniRef30)
    cdd_hits                   CDD hits per unique sequence (MD5): model id, E-value, bits, ranges
    cdd_models                 the CDD models hit: accession, short name, description, length
    cdd_runs                   the CDD search: database name and version, E-value, counts

Activation copies these into a dataset database through the staged-write protocol: family
rows into ``family_annotations`` (and the run and status tables), CDD hits expanded from
sequences to every protein with that sequence into ``annotations`` with ``source = 'cdd'``
(after checking that the database's unique sequences are the ones searched; with a row in
``annotation_source_versions``), and
creates the ``family_resemblance`` view, which applies the "resembles" rule.
"""

from __future__ import annotations

import hashlib
import json
import time
from collections.abc import Callable, Iterator, Sequence
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
CDD_MODELS = "cdd_models"
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
    hit_id BIGINT NOT NULL,
    sequence_md5 VARCHAR NOT NULL,
    pssm_id INTEGER NOT NULL,
    evalue FLOAT NOT NULL,
    score FLOAT,
    start_aa INTEGER,
    end_aa INTEGER,
    model_start INTEGER,
    model_end INTEGER,
    best_in_region BOOLEAN NOT NULL
)""",
    CDD_MODELS: f"""CREATE TABLE IF NOT EXISTS {CDD_MODELS} (
    pssm_id INTEGER PRIMARY KEY,
    accession VARCHAR NOT NULL,
    name VARCHAR,
    description VARCHAR,
    model_length INTEGER
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
    sequences_with_hits BIGINT NOT NULL,
    best_rows BIGINT NOT NULL,
    region_rule VARCHAR NOT NULL,
    max_overlap DOUBLE NOT NULL,
    export_digest VARCHAR NOT NULL,
    source_db VARCHAR,
    details JSON
)""",
}

# family_annotations keeps its primary key only; scans by family or accession take a fraction of a second,
# and secondary indexes would cost about 100 bytes per row. Files written with them lose them here.
INDEXES = f"""
DROP INDEX IF EXISTS idx_family_annotations_family;
DROP INDEX IF EXISTS idx_family_annotations_accession;
"""
STAGING_TABLES = (RUNS, STATUS, HITS, CDD_HITS, CDD_MODELS, CDD_RUNS)


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


# Reserved attachment names (a database file named like an alias would collide with it).
_SRC, _FAM, _STG, _PEEK, _FAMFILE, _OLD = "__rh_src", "__rh_fam", "__rh_stg", "__rh_peek", "__rh_famfile", "__rh_old"


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
        con.execute(f"ATTACH {_q(Path(families_db).expanduser().resolve())} AS {_FAM} (READ_ONLY)")
        if not _has_tables(con, FAMILY_TABLES, _FAM):
            raise RemoteHomologyError(f"{families_db} holds no protein-family tables")
        return f"{_FAM}."
    if not _has_tables(con, FAMILY_TABLES, _SRC):
        raise RemoteHomologyError(f"{db} holds no protein-family tables; pass a families staging file "
                                  "(sharur families run/import) with --families")
    return f"{_SRC}."


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
        con.execute(f"ATTACH {_q(db)} AS {_SRC} (READ_ONLY)")
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
                JOIN {_SRC}.proteins p ON p.protein_id = f.protein_id
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
_HAS_SEQ = "rtrim(COALESCE(sequence, ''), '*') <> ''"


def unique_sequences(con, prefix: str = "") -> Iterator[tuple[str, str]]:
    """``(md5, sequence)`` of every distinct protein sequence (trailing '*' removed), in MD5 order."""
    cur = con.execute(f"""SELECT md5(seq) AS h, any_value(seq) FROM (
                              SELECT rtrim(sequence, '*') AS seq FROM {prefix}proteins WHERE {_HAS_SEQ})
                          GROUP BY h ORDER BY h""")
    while True:
        batch = cur.fetchmany(50_000)
        if not batch:
            return
        yield from batch


def sequence_set_digest(con, prefix: str = "") -> tuple[int, str]:
    """Count and SHA-256 of the unique sequences as :func:`export_cdd` writes them (its ``export_digest``)."""
    digest = hashlib.sha256()
    n = 0
    for h, seq in unique_sequences(con, prefix):
        digest.update(f">{h}\n{seq}\n".encode())
        n += 1
    return n, digest.hexdigest()


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
    digest = hashlib.sha256()
    n = n_chunks = 0
    handle = None
    try:
        n_proteins = con.execute(f"SELECT COUNT(*) FROM proteins WHERE {_HAS_SEQ}").fetchone()[0]
        for h, seq in unique_sequences(con):
            if n % chunk_size == 0:
                if handle:
                    handle.close()
                handle = open(chunks / f"chunk_{n_chunks:05d}.faa", "w")  # noqa: SIM115 - one file per chunk
                n_chunks += 1
            text = f">{h}\n{seq}\n"
            handle.write(text)
            digest.update(text.encode())
            n += 1
    finally:
        if handle:
            handle.close()
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
_SOURCE_DB_KEY = {"hhsearch_pfam": "pfam", "hhsearch_pdb": "pdb", "hhblits_pdb": "pdb", "hhblits_uniref": "uniref"}


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


# Which CDD hits activation writes to ``annotations``.
CDD_ROWS_BEST = "best_per_region"
CDD_ROWS_ALL = "all"
CDD_ROW_CHOICES = (CDD_ROWS_BEST, CDD_ROWS_ALL)
# Same overlap measure and default as sharur.architecture.resolve (tests keep the defaults equal);
# defined here so this module imports without the rest of Sharur on cluster nodes.
MAX_OVERLAP = 0.5
REGION_RULE = ("greedy by bit score (ties: lower E-value, then lower model id): a hit is kept when its query "
               "range shares at most max_overlap of the shorter of the two ranges with every hit kept before "
               "it (the overlap measure of sharur.architecture.resolve)")


def best_per_region(hits: Sequence[tuple[float, float, int, int, int]], max_overlap: float = MAX_OVERLAP
                    ) -> list[bool]:
    """Which hits of one sequence are the best model of their region.

    ``hits`` are ``(bit score, E-value, model id, start, end)``. Hits are taken by bit score
    (ties: lower E-value, then lower model id); a hit is kept when, against every hit kept
    before it, the residues both query ranges share are at most ``max_overlap`` times the
    length of the shorter range. The best-scoring hit is always kept, so a sequence with any
    hit keeps at least one.
    """
    order = sorted(range(len(hits)), key=lambda i: (-(hits[i][0] or 0.0), hits[i][1], hits[i][2]))
    kept: list[tuple[int, int]] = []
    flags = [False] * len(hits)
    for i in order:
        start, end = hits[i][3], hits[i][4]
        length = end - start + 1
        if all(min(end, e) - max(start, s) + 1 <= max_overlap * min(length, e - s + 1) for s, e in kept):
            kept.append((start, end))
            flags[i] = True
    return flags


def _mark_best_per_region(con, max_overlap: float, log: Callable[[str], None]) -> int:
    """Set ``cdd_hits.best_in_region`` with :func:`best_per_region`, sequence by sequence; return the kept count."""
    started = time.time()
    cur = con.execute(f"""SELECT hit_id, sequence_md5, score, evalue, pssm_id, start_aa, end_aa FROM {CDD_HITS}
                          ORDER BY sequence_md5, hit_id""")
    keep: list[int] = []
    current, ids, hits = None, [], []

    def flush() -> None:
        for hit_id, flag in zip(ids, best_per_region(hits, max_overlap), strict=True):
            if flag:
                keep.append(hit_id)

    while True:
        batch = cur.fetchmany(500_000)
        if not batch:
            break
        for hit_id, md5, score, evalue, pssm, start, end in batch:
            if md5 != current:
                if ids:
                    flush()
                current, ids, hits = md5, [], []
            ids.append(hit_id)
            hits.append((score, evalue, pssm, start, end))
    if ids:
        flush()
    import pandas as pd  # noqa: PLC0415

    con.register("_rh_keep_df", pd.DataFrame({"hit_id": keep}, dtype="int64"))
    con.execute(f"UPDATE {CDD_HITS} SET best_in_region = hit_id IN (SELECT hit_id FROM _rh_keep_df)")
    con.unregister("_rh_keep_df")
    log(f"best model per region: {len(keep):,} hits kept (max overlap {max_overlap:g}, "
        f"{time.time() - started:.0f} s)")
    return len(keep)


_CDD_TYPES = ("VARCHAR", "VARCHAR", "INTEGER", "INTEGER", "INTEGER", "INTEGER", "INTEGER", "INTEGER", "DOUBLE",
              "DOUBLE", "VARCHAR")
_CDD_COLUMNS = ", ".join(f"'{c}': '{t}'" for c, t in zip(runner.CDD_COLUMNS, _CDD_TYPES, strict=True))


def compact_cdd(work: Path | str, *, max_evalue: float = CDD_MAX_EVALUE, allow_partial: bool = False,
                threads: int = 4, memory_limit: str = "8GB", log: Callable[[str], None] = print) -> dict[str, Any]:
    """Reduce the raw RPS-BLAST tables of ``work/cdd/out`` to two Parquet files in ``work/cdd/compact``.

    ``hits.parquet``: one row per hit at ``max_evalue`` or below and unique sequence (MD5),
    sorted by MD5, with the model's numeric id, E-value, bits and both ranges.
    ``models.parquet``: one row per CDD model hit (id, accession, short name, description,
    length) from titles ``"cd00009, AAA, description"``. ``compact.json`` records counts and
    the search's E-value. Run it beside the raw tables (on the cluster); only the compact
    files travel. Unless ``allow_partial``, every chunk must be done.
    """
    work = Path(work).expanduser().resolve()
    manifest = json.loads((work / "cdd" / "export.json").read_text())
    out = work / "cdd" / "out"
    done = sorted(p.with_suffix(".tsv") for p in out.glob("chunk_*.done") if p.with_suffix(".tsv").exists())
    n_chunks = manifest["chunks"]
    if len(done) != n_chunks and not allow_partial:
        raise RemoteHomologyError(f"{len(done):,} of {n_chunks:,} CDD chunks are done; finish the run first")
    if not done:
        raise RemoteHomologyError("no finished CDD chunk")
    records = [json.loads(p.with_suffix(".done").read_text()) for p in done]
    evalues = {r.get("evalue") for r in records}
    if None in evalues:
        raise RemoteHomologyError("CDD chunk records lack the E-value of their search")
    if min(float(e) for e in evalues) < max_evalue:
        raise RemoteHomologyError(f"the search kept hits at E <= {min(evalues)}, stricter than {max_evalue}; "
                                  "rerun it with the configured cdd_evalue")
    files = [str(p) for p in done if p.stat().st_size > 0]  # a chunk without hits leaves an empty table
    compact = work / "cdd" / "compact"
    compact.mkdir(parents=True, exist_ok=True)
    con = connect(":memory:", threads=threads, memory_limit=memory_limit)
    try:
        if files:
            listing = "[" + ", ".join(_q(f) for f in files) + "]"
            con.execute(f"""CREATE VIEW raw AS SELECT * FROM read_csv({listing}, delim = '\t', header = false,
                            quote = '', escape = '', columns = {{{_CDD_COLUMNS}}})""")
        else:
            con.execute("CREATE TABLE raw (" + ", ".join(
                f"{c} {t}" for c, t in zip(runner.CDD_COLUMNS, _CDD_TYPES, strict=True)) + ")")
        raw_rows = con.execute("SELECT COUNT(*) FROM raw").fetchone()[0]
        con.execute(f"""COPY (
            SELECT qseqid AS sequence_md5, CAST(split_part(sseqid, '|', 3) AS INTEGER) AS pssm_id,
                   CAST(evalue AS FLOAT) AS evalue, CAST(bitscore AS FLOAT) AS score, qstart AS start_aa,
                   qend AS end_aa, sstart AS model_start, send AS model_end
            FROM raw WHERE evalue <= {float(max_evalue)}
            ORDER BY sequence_md5, start_aa, pssm_id)
            TO {_q(compact / "hits.parquet")} (FORMAT parquet, COMPRESSION zstd)""")
        con.execute(f"""COPY (
            SELECT CAST(split_part(sseqid, '|', 3) AS INTEGER) AS pssm_id,
                   COALESCE(NULLIF(trim(split_part(any_value(stitle), ',', 1)), ''), split_part(sseqid, '|', 3))
                       AS accession,
                   NULLIF(trim(split_part(any_value(stitle), ',', 2)), '') AS name,
                   NULLIF(trim(regexp_extract(any_value(stitle), '^[^,]*,[^,]*,(.*)$', 1)), '') AS description,
                   any_value(slen) AS model_length
            FROM raw WHERE evalue <= {float(max_evalue)} GROUP BY sseqid ORDER BY pssm_id)
            TO {_q(compact / "models.parquet")} (FORMAT parquet, COMPRESSION zstd)""")
        rows, seqs = con.execute(f"SELECT COUNT(*), COUNT(DISTINCT qseqid) FROM raw WHERE evalue <= {float(max_evalue)}"
                                 ).fetchone()
    finally:
        con.close()
    versions = sorted({v for v in (json.loads(p.read_text()).get("versions", {}).get("rpsblast")
                                   for p in (work / "logs").glob("cdd_shard*.json")) if v})
    info = {"chunks": n_chunks, "chunks_done": len(done), "raw_rows": raw_rows, "rows": rows,
            "sequences_with_hits": seqs, "max_evalue": max_evalue, "search_evalue": sorted(float(e) for e in evalues),
            "rpsblast": versions, "export_digest": manifest["export_digest"],
            "databases": sorted({json.dumps(r.get("database"), sort_keys=True) for r in records}),
            "created_at": time.strftime("%Y-%m-%dT%H:%M:%S")}
    (compact / "compact.json").write_text(json.dumps(info, indent=2) + "\n")
    log(f"compacted {raw_rows:,} RPS-BLAST rows to {rows:,} hits at E <= {max_evalue:g} on {seqs:,} sequences")
    return info


def load_cdd(
    work: Path | str,
    staging: Path | str,
    config: Path | str,
    *,
    max_evalue: float = CDD_MAX_EVALUE,
    max_overlap: float = MAX_OVERLAP,
    allow_partial: bool = False,
    threads: int = 4,
    memory_limit: str = "8GB",
    log: Callable[[str], None] = print,
) -> dict[str, Any]:
    """Load the CDD hits of ``work/cdd`` into ``staging``: ``cdd_hits`` per unique sequence, ``cdd_models``.

    Every hit at ``max_evalue`` or below is kept in the staging file; ``best_in_region`` marks the
    best model of each region of each sequence (:func:`best_per_region` with ``max_overlap``),
    the rows activation writes to ``annotations`` by default.

    Uses ``work/cdd/compact`` (from :func:`compact_cdd`), compacting the raw tables first when
    they are present and the compact files are not. Hits stay keyed by sequence MD5; activation
    expands them to every protein of the target database with that sequence. The CDD tables of
    the staging file are replaced.
    """
    work = Path(work).expanduser().resolve()
    cfg = runner.load_config(config)
    cdd = runner.database(cfg, "cdd")
    manifest = json.loads((work / "cdd" / "export.json").read_text())
    compact = work / "cdd" / "compact"
    if not (compact / "compact.json").exists():
        compact_cdd(work, max_evalue=max_evalue, allow_partial=allow_partial, threads=threads,
                    memory_limit=memory_limit, log=log)
    info = json.loads((compact / "compact.json").read_text())
    if info["export_digest"] != manifest["export_digest"]:
        raise RemoteHomologyError("the compact CDD files belong to another export; rerun compact-cdd")
    if info["chunks_done"] != info["chunks"] and not allow_partial:
        raise RemoteHomologyError(f"{info['chunks_done']:,} of {info['chunks']:,} CDD chunks were compacted")
    if float(info["max_evalue"]) != float(max_evalue):
        raise RemoteHomologyError(f"compact files hold hits at E <= {info['max_evalue']}, not {max_evalue}")
    staging = Path(staging).expanduser().resolve()
    staging.parent.mkdir(parents=True, exist_ok=True)
    cdd_run = "cdd-" + _sha({"db": [cdd["name"], cdd["version"]], "evalue": max_evalue,
                             "export": manifest["export_digest"], "max_overlap": max_overlap})[:10]
    con = connect(staging, threads=threads, memory_limit=memory_limit)
    try:
        con.execute("BEGIN TRANSACTION")
        for table in (CDD_HITS, CDD_MODELS, CDD_RUNS):
            con.execute(f"DROP TABLE IF EXISTS {table}")
            con.execute(DDL[table])
        con.execute(f"""INSERT INTO {CDD_HITS}
            SELECT row_number() OVER (ORDER BY sequence_md5, start_aa, pssm_id, end_aa, score DESC) AS hit_id,
                   sequence_md5, pssm_id, evalue, score, start_aa, end_aa, model_start, model_end, false
            FROM read_parquet({_q(compact / 'hits.parquet')})""")
        con.execute(f"INSERT INTO {CDD_MODELS} SELECT * FROM read_parquet({_q(compact / 'models.parquet')})")
        rows, seqs = con.execute(f"SELECT COUNT(*), COUNT(DISTINCT sequence_md5) FROM {CDD_HITS}").fetchone()
        orphans = con.execute(f"SELECT COUNT(*) FROM (SELECT DISTINCT pssm_id FROM {CDD_HITS}) "
                              f"ANTI JOIN {CDD_MODELS} USING (pssm_id)").fetchone()[0]
        if rows != info["rows"] or orphans:
            raise RemoteHomologyError(f"compact files are inconsistent ({rows:,} rows vs {info['rows']:,} recorded, "
                                      f"{orphans} models missing)")
        best = _mark_best_per_region(con, max_overlap, log)
        con.execute(f"INSERT INTO {CDD_RUNS} VALUES (?, ?::TIMESTAMP, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                    [cdd_run, _now(), cdd["name"], cdd["version"], cdd["prefix"], float(max_evalue),
                     "; ".join(info.get("rpsblast") or []) or None, manifest["unique_sequences"],
                     manifest["proteins"], rows, seqs, best, REGION_RULE, float(max_overlap),
                     manifest["export_digest"], manifest.get("source_db"),
                     json.dumps({**{k: info[k] for k in ("chunks", "chunks_done", "raw_rows", "search_evalue")},
                                 "work_dir": str(work), "config_path": cfg["_path"]})])
        con.execute("COMMIT")
    except Exception:
        try:
            con.execute("ROLLBACK")
        except duckdb.Error:
            pass
        raise
    finally:
        con.close()
    result = {"cdd_run": cdd_run, "chunks": info["chunks"], "chunks_done": info["chunks_done"],
              "raw_rows": info["raw_rows"], "hit_rows": rows, "best_in_region_rows": best,
              "max_overlap": max_overlap, "sequences_with_hits": seqs,
              "unique_sequences": manifest["unique_sequences"], "proteins": manifest["proteins"],
              "staging": str(staging)}
    log(f"{cdd_run}: {rows:,} CDD hits on {seqs:,} of {manifest['unique_sequences']:,} unique sequences")
    return result


# --------------------------------------------------------------------------- #
# Activation (copy staged rows into a dataset database)
# --------------------------------------------------------------------------- #
def _attach_staging(con, path: Path, alias: str) -> set[str]:
    con.execute(f"ATTACH {_q(path)} AS {alias} (READ_ONLY)")
    return {r[0] for r in con.execute(
        "SELECT table_name FROM information_schema.tables WHERE table_catalog = ? AND table_schema = 'main'",
        [alias]).fetchall()}


FAMILY_RUNS_SEARCHED = "searched"
FAMILY_RUNS_ALL = "all"


def _family_runs_in(con, path: Path) -> list[str]:
    con.execute(f"ATTACH {_q(path)} AS {_FAMFILE} (READ_ONLY)")
    try:
        if not _has_tables(con, ["protein_family_runs"], _FAMFILE):
            raise RemoteHomologyError(f"{path} holds no protein-family tables")
        return [r[0] for r in con.execute(f"SELECT clustering_run FROM {_FAMFILE}.protein_family_runs ORDER BY 1"
                                          ).fetchall()]
    finally:
        con.execute(f"DETACH {_FAMFILE}")


def load_staging(con: duckdb.DuckDBPyConnection, staging: Sequence[Path | str], *,
                 families_staging: Path | str | Sequence[Path | str] | None = None,
                 family_runs: str = FAMILY_RUNS_SEARCHED, verify_digest: bool = True,
                 cdd_rows: str = CDD_ROWS_BEST, log: Callable[[str], None] = print) -> dict[str, Any]:
    """Copy the remote-homology rows of ``staging`` file(s) into the database open on ``con``.

    Protein-family staging files (``families_staging``, one or several, loaded first in the order
    given; a later file's run replaces a same-named one) contribute the searched clustering runs
    (``family_runs = 'searched'``) or every run they hold (``'all'``), so one activation can load
    a dataset's own build, imported clusterings and the remote-homology rows together.

    Family tier: the clustering run of every staged search must exist in the database (or in
    ``families_staging``) with the clustering input digest the
    search recorded, and every staged family must exist in it; rows of the same clustering run
    and source are replaced. CDD: every staged protein must exist with the sequence the hit was
    computed for (MD5); ``annotations`` rows with ``source = 'cdd'`` are replaced and the
    database name and version are recorded in ``annotation_source_versions``. ``cdd_rows``
    chooses the CDD rows written: ``best_per_region`` (default; the best model of each region,
    ``cdd_hits.best_in_region``) or ``all``. The full per-sequence table stays in the staging
    file either way.
    """
    from sharur.defensefinder.selection import _SOURCE_VERSIONS_SQL  # noqa: PLC0415

    summary: dict[str, Any] = {"family_tier": None, "cdd": None}
    if family_runs not in (FAMILY_RUNS_SEARCHED, FAMILY_RUNS_ALL):
        raise RemoteHomologyError(f"family_runs must be '{FAMILY_RUNS_SEARCHED}' or '{FAMILY_RUNS_ALL}'")
    if isinstance(families_staging, (str, Path)):
        families_staging = [families_staging]
    families_staging = [Path(f).expanduser().resolve() for f in (families_staging or [])]
    if families_staging:
        from sharur.protein_families import store as pf  # noqa: PLC0415

        needed: set[str] = set()
        for path in staging:
            con.execute(f"ATTACH {_q(Path(path).expanduser().resolve())} AS {_PEEK} (READ_ONLY)")
            try:
                if _has_tables(con, [RUNS], _PEEK):
                    needed |= {r[0] for r in con.execute(f"SELECT clustering_run FROM {_PEEK}.{RUNS}").fetchall()}
            finally:
                con.execute(f"DETACH {_PEEK}")
        loaded: list[dict[str, Any]] = []
        for fam_file in families_staging:
            available = _family_runs_in(con, fam_file)
            runs = available if family_runs == FAMILY_RUNS_ALL else sorted(set(available) & needed)
            if not runs:
                continue
            pf.load_staging(con, fam_file, runs=runs, verify_digest=verify_digest, log=log)
            loaded.append({"file": str(fam_file), "runs": runs})
        found = {r for item in loaded for r in item["runs"]}
        if family_runs == FAMILY_RUNS_SEARCHED and needed - found:
            log(f"searched runs taken from the database itself: {', '.join(sorted(needed - found))}")
        summary["families_loaded"] = loaded
    for i, path in enumerate(staging):
        path = Path(path).expanduser().resolve()
        alias = f"{_STG}{i}"
        tables = _attach_staging(con, path, alias)
        try:
            if {RUNS, STATUS, HITS} <= tables:
                summary["family_tier"] = _load_family_tier(con, alias, path)
                log(f"family tier: {summary['family_tier']}")
            if {CDD_HITS, CDD_MODELS, CDD_RUNS} <= tables:
                summary["cdd"] = _load_cdd_rows(con, alias, path, verify_digest, _SOURCE_VERSIONS_SQL, log,
                                                cdd_rows=cdd_rows)
                log(f"cdd: {summary['cdd']}")
            if not ({RUNS, STATUS, HITS} <= tables or {CDD_HITS, CDD_MODELS, CDD_RUNS} <= tables):
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


def _load_cdd_rows(con, alias: str, path: Path, verify_digest: bool, source_versions_sql: str,
                   log: Callable[[str], None], *, cdd_rows: str = CDD_ROWS_BEST) -> dict[str, Any]:
    if cdd_rows not in CDD_ROW_CHOICES:
        raise RemoteHomologyError(f"cdd_rows must be one of {', '.join(CDD_ROW_CHOICES)}, not {cdd_rows!r}")
    run = con.execute(f"SELECT * FROM {alias}.{CDD_RUNS}").fetchdf()
    if len(run) != 1:
        raise RemoteHomologyError(f"{path}: expected one CDD run record, found {len(run)}")
    run = run.iloc[0]
    if verify_digest:
        n, digest = sequence_set_digest(con)
        if digest != run["export_digest"]:
            raise RemoteHomologyError(
                f"the database's {n:,} unique protein sequences differ from the {int(run['unique_sequences']):,} "
                "that were searched (proteins were added, removed or changed sequence); rerun the CDD search")
        log(f"unique-sequence digest verified ({n:,} sequences)")
    con.execute("BEGIN TRANSACTION")
    try:
        con.execute(source_versions_sql)
        removed = con.execute("SELECT COUNT(*) FROM annotations WHERE source = ?", [CDD_SOURCE]).fetchone()[0]
        con.execute("DELETE FROM annotations WHERE source = ?", [CDD_SOURCE])
        base = con.execute("SELECT COALESCE(MAX(annotation_id), 0) FROM annotations").fetchone()[0]
        con.execute(f"""CREATE OR REPLACE TEMP TABLE _rh_pmd5 AS
            SELECT protein_id, md5(rtrim(sequence, '*')) AS sequence_md5 FROM proteins WHERE {_HAS_SEQ}""")
        all_rows = con.execute(f"""SELECT COALESCE(SUM(h.n), 0) FROM _rh_pmd5 p JOIN
            (SELECT sequence_md5, COUNT(*) AS n FROM {alias}.{CDD_HITS} GROUP BY 1) h USING (sequence_md5)""").fetchone()[0]
        where = "WHERE h.best_in_region" if cdd_rows == CDD_ROWS_BEST else ""
        con.execute(f"""
            INSERT INTO annotations (annotation_id, protein_id, source, accession, name, description, evalue,
                                     score, start_aa, end_aa)
            SELECT ? + row_number() OVER (ORDER BY p.protein_id, h.start_aa, m.accession, h.evalue), p.protein_id,
                   '{CDD_SOURCE}', m.accession, m.name, m.description, h.evalue, h.score, h.start_aa, h.end_aa
            FROM {alias}.{CDD_HITS} h JOIN {alias}.{CDD_MODELS} m USING (pssm_id) JOIN _rh_pmd5 p USING (sequence_md5)
            {where}""", [int(base)])
        loaded, proteins = con.execute("SELECT COUNT(*), COUNT(DISTINCT protein_id) FROM annotations "
                                       "WHERE source = ?", [CDD_SOURCE]).fetchone()
        rule = (f"rpsblast E <= {float(run['max_evalue']):g}; "
                + (f"best model per region (max overlap {float(run['max_overlap']):g}); "
                   if cdd_rows == CDD_ROWS_BEST else "every hit; ")
                + "unique sequences expanded to proteins")
        params = {"cdd_run": run["cdd_run"], "max_evalue": float(run["max_evalue"]), "db_prefix": run["db_prefix"],
                  "tool_version": run["tool_version"], "unique_sequences": int(run["unique_sequences"]),
                  "rows_stored": cdd_rows, "region_rule": run["region_rule"], "max_overlap": float(run["max_overlap"]),
                  "sequence_hits": int(run["hit_rows"]), "sequence_best_hits": int(run["best_rows"]),
                  "protein_rows_all": int(all_rows), "protein_rows_stored": int(loaded),
                  "proteins_with_hits": int(proteins), "full_table": f"cdd_hits and cdd_models in {path}",
                  "export_digest": run["export_digest"], "evidence": "resembles"}
        con.execute(f"DELETE FROM {SOURCE_VERSIONS} WHERE source = ?", [CDD_SOURCE])
        con.execute(f"""INSERT INTO {SOURCE_VERSIONS} (source, model_set, model_version, profiles_path, selection_rule,
                                    parameters, hits_path, raw_rows, selected_rows, loaded_rows, loaded_at)
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?::TIMESTAMP)""",
                    [CDD_SOURCE, run["db_name"], run["db_version"], run["db_prefix"], rule,
                     json.dumps(params, sort_keys=True), str(path), int(all_rows), int(loaded), int(loaded), _now()])
        con.execute("DROP TABLE IF EXISTS _rh_pmd5")
        con.execute("COMMIT")
    except Exception:
        con.execute("ROLLBACK")
        raise
    return {"staging": str(path), "rows_removed": removed, "rows_loaded": loaded, "rows_all": int(all_rows),
            "rows_stored": cdd_rows, "proteins_with_hits": proteins, "db": f"{run['db_name']} {run['db_version']}"}


def _fingerprint(con, table: str, alias: str | None, columns: Sequence[str]) -> tuple[int, int]:
    name = f"{alias}.{table}" if alias else table
    packed = "struct_pack(" + ", ".join(f'"{c}" := t."{c}"' for c in columns) + ")"
    return con.execute(f"SELECT COUNT(*), COALESCE(SUM(hash({packed})::HUGEINT), 0) FROM {name} t").fetchone()


def compact_staging_file(path: Path | str, *, threads: int = 4, memory_limit: str = "8GB") -> dict[str, Any]:
    """Rewrite a remote-homology staging file without free blocks or legacy indexes; replace it atomically.

    Every remote-homology table is copied whole (columns matched by name) and checked by row count
    and a content fingerprint; the original is replaced with one rename.
    """
    import os  # noqa: PLC0415

    path = Path(path).expanduser().resolve()
    tmp = path.with_name(path.name + ".compact.tmp")
    for stale in (tmp, Path(str(tmp) + ".wal")):
        stale.unlink(missing_ok=True)
    before = path.stat().st_size
    con = connect(tmp, threads=threads, memory_limit=memory_limit)
    try:
        con.execute(f"ATTACH {_q(path)} AS {_OLD} (READ_ONLY)")
        present = [t for t in STAGING_TABLES if _has_tables(con, [t], _OLD)]
        if not present:
            raise RemoteHomologyError(f"{path} holds no remote-homology tables")
        rows = {}
        for table in present:
            con.execute(DDL[table])
            con.execute(f"INSERT INTO {table} BY NAME SELECT * FROM {_OLD}.{table}")
            cols = [r[0] for r in con.execute(
                "SELECT column_name FROM information_schema.columns WHERE table_catalog = current_database() "
                "AND table_name = ? ORDER BY ordinal_position", [table]).fetchall()]
            new, old = _fingerprint(con, table, None, cols), _fingerprint(con, table, _OLD, cols)
            if new != old:
                raise RemoteHomologyError(f"{table}: copy differs from the original ({old} -> {new})")
            rows[table] = int(new[0])
        if HITS in present:
            con.execute(resemblance_view_sql())
        con.execute(f"DETACH {_OLD}")
        con.execute("CHECKPOINT")
    except Exception:
        con.close()
        tmp.unlink(missing_ok=True)
        Path(str(tmp) + ".wal").unlink(missing_ok=True)
        raise
    con.close()
    os.replace(tmp, path)
    return {"rows": rows, "bytes_before": before, "bytes_after": path.stat().st_size}


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
            version, params = None, {}
            if _has_tables(con, [SOURCE_VERSIONS]):
                row = con.execute(f"SELECT model_set, model_version, parameters FROM {SOURCE_VERSIONS} "
                                  "WHERE source = ?", [CDD_SOURCE]).fetchone()
                if row:
                    version = f"{row[0]} {row[1]}"
                    params = json.loads(row[2] or "{}")
            out["cdd"] = {"rows": n, "proteins": proteins, "database": version,
                          "rows_stored": params.get("rows_stored"), "rows_all_hits": params.get("protein_rows_all"),
                          "max_overlap": params.get("max_overlap"), "sequence_hits": params.get("sequence_hits"),
                          "full_table": params.get("full_table")}
    if _has_tables(con, [CDD_RUNS]):  # a staging file: the CDD search waiting for activation
        row = con.execute(f"""SELECT cdd_run, db_name, db_version, max_evalue, unique_sequences, proteins, hit_rows,
                                     sequences_with_hits, best_rows, max_overlap FROM {CDD_RUNS}""").fetchone()
        if row:
            out["cdd_staged"] = dict(zip(("cdd_run", "db_name", "db_version", "max_evalue", "unique_sequences",
                                          "proteins", "sequence_hits", "sequences_with_hits", "sequence_best_hits",
                                          "max_overlap"), row))
    return out
