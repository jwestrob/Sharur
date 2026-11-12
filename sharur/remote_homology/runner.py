#!/usr/bin/env python3
"""Run the remote-homology searches on a prepared work directory (standard library only).

The work directory comes from ``sharur remote-homology export`` (family tier) or
``sharur remote-homology export-cdd`` (CDD). This script needs Python 3.9+, HH-suite 3, MAFFT
and BLAST+, never DuckDB or Sharur, so it runs on any cluster node that sees the work
directory and the search databases::

    python3 runner.py check    --config CONFIG.json
    python3 runner.py families --work WORK --config CONFIG.json [--shard I --n-shards N]
    python3 runner.py cdd      --work WORK --config CONFIG.json [--shard I --n-shards N]

Workers. ``--threads`` (default ``$SLURM_CPUS_ON_NODE`` or the host CPU count) is split into
``threads // --threads-per-worker`` worker processes (default 4 threads each), so a 24-CPU
job runs six searches at once and every allocated CPU is used.

Family tier, per family (``out/<bucket>/<family_key>/``):

1. ``aln.fa``: ``mafft --auto`` of the family FASTA (representative first; copied when the
   family has one unique sequence).
2. ``family_only``: ``hhsearch -i aln.fa -M first`` against each target database.
3. ``enriched`` (when listed in ``msa_kinds``): ``hhblits -i aln.fa -M first`` against UniRef30
   (``-n 2 -e 0.001``) writes ``enr.a3m`` and ``uniref.hhr``; ``hhsearch -i enr.a3m`` against
   each target database.
4. Every ``.hhr`` is parsed into ``hits.tsv`` (hits at probability >= 50%, the first 20 for
   UniRef30). The query length of every result must equal the representative's length.
5. ``done.json`` records status, step times, tool versions and the configuration digest; the
   ``.hhr`` files are gzipped and the alignments removed (``--keep-intermediates`` keeps them).

A family with ``done.json`` for the current configuration digest is skipped, so a rerun
continues where a job stopped. Steps whose output exists are skipped as well. A lock
directory per family keeps two jobs off one family.

CDD, per chunk of unique sequences (``cdd/chunks/chunk_NNNNN.faa``): ``rpsblast -evalue 0.01``
writes ``cdd/out/chunk_NNNNN.tsv`` and a ``.done`` marker.

Configuration (JSON; site paths live here, never in code)::

    {"hhsuite_bin": "/path/to/hhsuite/bin", "mafft": "mafft", "rpsblast": "rpsblast",
     "databases": {
       "pfam":   {"prefix": "/path/to/pfam_hh/pfam", "name": "Pfam", "version": "39.0", "kind": "pfam"},
       "pdb":    {"prefix": "/path/to/pdb100", "name": "PDB100", "version": "YYYY-MM-DD", "kind": "pdb",
                  "titles": "/path/to/pdb_titles.tsv"},
       "uniref": {"prefix": "/path/to/UniRef30_YYYY_MM", "name": "UniRef30", "version": "YYYY_MM",
                  "kind": "uniref"},
       "cdd":    {"prefix": "/path/to/cdd/Cdd", "name": "CDD", "version": "3.xx"}},
     "family_targets": ["pfam", "pdb"], "msa_kinds": ["enriched", "family_only"],
     "search": {"hhblits_iterations": 2, "hhblits_evalue": 0.001, "max_hits": 100,
                "min_probability": 50, "uniref_max_hits": 20, "cdd_evalue": 0.01}}
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import shutil
import socket
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path


try:  # imported as sharur.remote_homology.runner
    from sharur.remote_homology import hhr
except ImportError:  # run as a script from its own directory
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import hhr  # type: ignore[no-redef]

SEARCH_DEFAULTS = {
    "hhblits_iterations": 2,
    "hhblits_evalue": 0.001,
    "max_hits": 100,
    "min_probability": hhr.STORE_MIN_PROBABILITY,
    "uniref_max_hits": 20,
    "cdd_evalue": 0.01,
}
DEFAULT_TARGETS = ("pfam", "pdb")
DEFAULT_MSA_KINDS = (hhr.MSA_FAMILY_ONLY,)
TARGET_SOURCE = {"pfam": hhr.SOURCE_PFAM, "pdb": hhr.SOURCE_PDB}
CDD_OUTFMT = "6 qseqid sseqid qlen qstart qend slen sstart send evalue bitscore stitle"
CDD_COLUMNS = ("qseqid", "sseqid", "qlen", "qstart", "qend", "slen", "sstart", "send", "evalue", "bitscore",
               "stitle")
FAMILY_TSV_COLUMNS = ("family_index", "family_key", "family_id", "n_members", "n_genomes", "n_unique",
                      "n_sequences", "rep_protein_id", "rep_length", "fasta")


class ConfigError(RuntimeError):
    pass


# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
def load_config(path: str | Path) -> dict:
    cfg = json.loads(Path(path).read_text())
    cfg.setdefault("hhsuite_bin", "")
    cfg.setdefault("mafft", "mafft")
    cfg.setdefault("rpsblast", "rpsblast")
    cfg.setdefault("databases", {})
    cfg["search"] = {**SEARCH_DEFAULTS, **cfg.get("search", {})}
    cfg["family_targets"] = list(cfg.get("family_targets") or DEFAULT_TARGETS)
    cfg["msa_kinds"] = list(cfg.get("msa_kinds") or DEFAULT_MSA_KINDS)
    for kind in cfg["msa_kinds"]:
        if kind not in hhr.MSA_KINDS:
            raise ConfigError(f"msa_kinds: unknown kind {kind!r} (use {', '.join(hhr.MSA_KINDS)})")
    for target in cfg["family_targets"]:
        if target not in TARGET_SOURCE:
            raise ConfigError(f"family_targets: unknown target {target!r} (use {', '.join(TARGET_SOURCE)})")
    cfg["_path"] = str(Path(path).resolve())
    return cfg


def family_digest(cfg: dict) -> str:
    """Digest of everything that changes family-tier results (databases, kinds, search parameters)."""
    dbs = {k: {f: v.get(f) for f in ("prefix", "name", "version", "kind")} for k, v in cfg["databases"].items()
           if k in (*cfg["family_targets"], *(["uniref"] if hhr.MSA_ENRICHED in cfg["msa_kinds"] else []))}
    payload = {"databases": dbs, "targets": sorted(cfg["family_targets"]), "msa_kinds": sorted(cfg["msa_kinds"]),
               "search": {k: cfg["search"][k] for k in sorted(cfg["search"]) if k != "cdd_evalue"}}
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()[:16]


def tool(cfg: dict, name: str) -> str:
    if name in ("hhblits", "hhsearch", "hhmake", "ffindex_build"):
        return os.path.join(cfg["hhsuite_bin"], name) if cfg["hhsuite_bin"] else name
    return cfg.get(name, name)


def database(cfg: dict, key: str) -> dict:
    db = cfg["databases"].get(key)
    if not db or not db.get("prefix"):
        raise ConfigError(f"configuration names no '{key}' database (databases.{key}.prefix)")
    for field in ("name", "version"):
        if not db.get(field):
            raise ConfigError(f"databases.{key}.{field} is required: every stored row carries it")
    return db


def _db_present(key: str, prefix: str) -> bool:
    if key == "cdd":
        return any(Path(prefix + ext).exists() for ext in (".pal", ".rps", ".00.rps"))
    return any(Path(f"{prefix}_{part}.ffindex").exists() for part in ("hhm", "a3m"))


def check_database(cfg: dict, key: str) -> dict:
    """The configured database ``key``, after checking that its files exist."""
    db = database(cfg, key)
    if not _db_present(key, db["prefix"]):
        raise ConfigError(f"database {key}: nothing found at prefix {db['prefix']}")
    return db


def tool_version(command: list[str]) -> str | None:
    try:
        out = subprocess.run(command, check=False, capture_output=True, text=True)
    except OSError:
        return None
    text = (out.stdout or "") + (out.stderr or "")
    for line in text.splitlines():
        line = line.strip()
        if line:
            return line[:200]
    return None


def versions(cfg: dict) -> dict:
    return {"hhblits": tool_version([tool(cfg, "hhblits"), "-h"]),
            "hhsearch": tool_version([tool(cfg, "hhsearch"), "-h"]),
            "mafft": tool_version([tool(cfg, "mafft"), "--version"]),
            "rpsblast": tool_version([tool(cfg, "rpsblast"), "-version"])}


def check(cfg: dict) -> dict:
    """Tools on PATH and database files present; raises ConfigError on a missing piece."""
    problems = []
    for name in ("hhblits", "hhsearch", "mafft", "rpsblast"):
        exe = tool(cfg, name)
        if shutil.which(exe) is None:
            problems.append(f"{name}: {exe} not found")
    for key, db in cfg["databases"].items():
        prefix = db.get("prefix", "")
        if not _db_present(key, prefix):
            problems.append(f"database {key}: nothing found at prefix {prefix}")
        if db.get("titles") and not Path(db["titles"]).exists():
            problems.append(f"database {key}: titles file {db['titles']} missing")
    report = {"versions": versions(cfg), "problems": problems, "family_digest": family_digest(cfg)}
    if problems:
        raise ConfigError("; ".join(problems))
    return report


# --------------------------------------------------------------------------- #
# Small helpers
# --------------------------------------------------------------------------- #
def allocated_cpus() -> int:
    value = os.environ.get("SLURM_CPUS_ON_NODE") or os.environ.get("SLURM_CPUS_PER_TASK")
    if value:
        return int(value)
    try:
        return len(os.sched_getaffinity(0))  # type: ignore[attr-defined]
    except AttributeError:
        return os.cpu_count() or 1


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header, chunks = None, []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def run(command: list[str], **kwargs) -> float:
    started = time.time()
    proc = subprocess.run(command, capture_output=True, text=True, check=False, **kwargs)
    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout or "").strip().splitlines()[-3:]
        raise RuntimeError(f"{Path(command[0]).name} exited {proc.returncode}: {' | '.join(tail)[:400]}")
    return time.time() - started


def write_json(path: Path, payload: dict) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(payload, indent=1, sort_keys=True, default=str) + "\n")
    os.replace(tmp, path)


def gzip_in_place(path: Path) -> Path:
    target = path.with_name(path.name + ".gz")
    with open(path, "rb") as src, gzip.open(target, "wb", compresslevel=6) as dst:
        shutil.copyfileobj(src, dst)
    path.unlink()
    return target


def read_hhr(path: Path) -> hhr.HHRResult:
    if path.exists():
        return hhr.parse_hhr(path)
    gz = path.with_name(path.name + ".gz")
    with gzip.open(gz, "rt", errors="replace") as handle:
        return hhr.parse_hhr_text(handle.read())


def hhr_present(path: Path) -> bool:
    return path.exists() or path.with_name(path.name + ".gz").exists()


class FamilyLock:
    """``mkdir`` lock per family; a lock older than ``stale_hours`` is taken over."""

    def __init__(self, directory: Path, stale_hours: float):
        self.path = directory / ".lock"
        self.stale = stale_hours * 3600

    def acquire(self) -> bool:
        try:
            self.path.mkdir()
        except FileExistsError:
            try:
                age = time.time() - self.path.stat().st_mtime
            except FileNotFoundError:
                return self.acquire()
            if age < self.stale:
                return False
            shutil.rmtree(self.path, ignore_errors=True)
            try:
                self.path.mkdir()
            except FileExistsError:
                return False
        (self.path / "owner").write_text(json.dumps({"host": socket.gethostname(), "pid": os.getpid(),
                                                     "job": os.environ.get("SLURM_JOB_ID"), "time": time.time()}))
        return True

    def release(self) -> None:
        shutil.rmtree(self.path, ignore_errors=True)


# --------------------------------------------------------------------------- #
# Family tier
# --------------------------------------------------------------------------- #
def read_families(work: Path) -> list[dict]:
    rows = []
    with open(work / "families.tsv") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            values = line.rstrip("\n").split("\t")
            row = dict(zip(header, values))
            for key in ("family_index", "n_members", "n_genomes", "n_unique", "n_sequences", "rep_length"):
                row[key] = int(row[key])
            rows.append(row)
    return rows


def family_dir(work: Path, family: dict) -> Path:
    return work / "out" / f"{family['family_index'] // 1000:03d}" / family["family_key"]


def _hhsearch(cfg: dict, inp: Path, db: dict, out: Path, threads: int, *, master_first: bool) -> float:
    s = cfg["search"]
    command = [tool(cfg, "hhsearch"), "-i", str(inp), *(["-M", "first"] if master_first else []),
               "-d", db["prefix"], "-o", str(out) + ".tmp", "-cpu", str(threads), "-v", "0",
               "-Z", str(s["max_hits"]), "-B", str(s["max_hits"])]
    seconds = run(command)
    os.replace(str(out) + ".tmp", out)
    return seconds


def run_family(job: tuple) -> dict:
    """Search one family; returns its ``done.json`` payload (or a skip/lock/error record)."""
    work, family, cfg, threads, force, keep, stale_hours = job
    work = Path(work)
    d = family_dir(work, family)
    d.mkdir(parents=True, exist_ok=True)
    digest = family_digest(cfg)
    done_path = d / "done.json"
    if not force and done_path.exists():
        try:
            done = json.loads(done_path.read_text())
        except json.JSONDecodeError:
            done = {}
        if done.get("config_digest") == digest and done.get("status") == "ok":
            return {"family_key": family["family_key"], "status": "skipped"}
    lock = FamilyLock(d, stale_hours)
    if not lock.acquire():
        return {"family_key": family["family_key"], "status": "locked"}
    started = time.time()
    steps: dict[str, float] = {}
    record = {"family_key": family["family_key"], "family_id": family["family_id"],
              "family_index": family["family_index"], "config_digest": digest, "host": socket.gethostname(),
              "job": os.environ.get("SLURM_JOB_ID"), "threads": threads, "n_sequences": family["n_sequences"],
              "rep_length": family["rep_length"], "msa_kinds": cfg["msa_kinds"],
              "targets": cfg["family_targets"], "started": time.strftime("%Y-%m-%dT%H:%M:%S")}
    try:
        if force:
            for old in list(d.glob("*.hhr*")) + [d / "hits.tsv", d / "aln.fa", d / "enr.a3m", done_path]:
                old.unlink(missing_ok=True)
        fasta = work / family["fasta"]
        aln = d / "aln.fa"
        outputs: list[tuple[str, str, str, Path]] = []  # (source, msa_kind, kind, hhr path)
        for kind in cfg["msa_kinds"]:
            for target in cfg["family_targets"]:
                outputs.append((TARGET_SOURCE[target], kind, cfg["databases"][target].get("kind", target),
                                d / f"{target}_{kind}.hhr"))
        if hhr.MSA_ENRICHED in cfg["msa_kinds"]:
            outputs.append((hhr.SOURCE_UNIREF, hhr.MSA_ENRICHED, "uniref", d / "uniref.hhr"))
        need_aln = any(not hhr_present(p) for *_, p in outputs)
        enr = d / "enr.a3m"
        need_enr = hhr.MSA_ENRICHED in cfg["msa_kinds"] and any(
            not hhr_present(p) for _, k, _, p in outputs if k == hhr.MSA_ENRICHED)
        if (need_aln or need_enr) and not aln.exists():
            if family["n_sequences"] <= 1:
                shutil.copyfile(fasta, aln)
                steps["mafft"] = 0.0
            else:
                with open(str(aln) + ".tmp", "w") as handle:
                    steps["mafft"] = _mafft(cfg, fasta, handle)
                os.replace(str(aln) + ".tmp", aln)
        if need_enr and (not enr.exists() or not hhr_present(d / "uniref.hhr")):
            uniref = database(cfg, "uniref")
            s = cfg["search"]
            steps["hhblits"] = run([tool(cfg, "hhblits"), "-i", str(aln), "-M", "first", "-d", uniref["prefix"],
                                    "-n", str(s["hhblits_iterations"]), "-e", str(s["hhblits_evalue"]),
                                    "-cpu", str(threads), "-o", str(d / "uniref.hhr") + ".tmp",
                                    "-oa3m", str(enr) + ".tmp", "-v", "0"])
            os.replace(str(enr) + ".tmp", enr)
            os.replace(str(d / "uniref.hhr") + ".tmp", d / "uniref.hhr")
            (d / "uniref.hhr.gz").unlink(missing_ok=True)
        for source, kind, _dbkind, path in outputs:
            if source == hhr.SOURCE_UNIREF or hhr_present(path):
                continue
            target = path.name.split("_", 1)[0]
            db = database(cfg, target)
            if kind == hhr.MSA_ENRICHED:
                steps[f"hhsearch_{target}_{kind}"] = _hhsearch(cfg, enr, db, path, threads, master_first=False)
            else:
                steps[f"hhsearch_{target}_{kind}"] = _hhsearch(cfg, aln, db, path, threads, master_first=True)
        # parse and verify
        s = cfg["search"]
        rows: list[tuple] = []
        match_columns: dict[str, int | None] = {}
        for source, kind, dbkind, path in outputs:
            result = read_hhr(path)
            match_columns[path.name] = result.match_columns
            titles = _titles(cfg, dbkind)
            rows.extend(hhr.hit_rows(result, source=source, msa_kind=kind, kind=dbkind,
                                     min_probability=s["min_probability"],
                                     max_hits=s["uniref_max_hits"] if source == hhr.SOURCE_UNIREF else None,
                                     titles=titles))
        mismatched = {k: v for k, v in match_columns.items() if v != family["rep_length"]}
        tmp = d / "hits.tsv.tmp"
        with open(tmp, "w") as handle:
            handle.write("\t".join(hhr.HIT_COLUMNS) + "\n")
            for row in rows:
                handle.write("\t".join(_cell(v) for v in row) + "\n")
        os.replace(tmp, d / "hits.tsv")
        record.update(status="ok", n_hits=len(rows), match_columns=match_columns, verified=not mismatched,
                      steps=steps, seconds=round(time.time() - started, 2))
        if mismatched:
            record["warning"] = f"query length differs from the representative ({family['rep_length']}): {mismatched}"
        if not keep and not mismatched:
            for _source, _kind, _dbkind, path in outputs:
                if path.exists():
                    gzip_in_place(path)
            aln.unlink(missing_ok=True)
            enr.unlink(missing_ok=True)
        record["finished"] = time.strftime("%Y-%m-%dT%H:%M:%S")
        write_json(done_path, record)
        return record
    except Exception as exc:
        record.update(status="error", error=f"{type(exc).__name__}: {exc}"[:600], steps=steps,
                      seconds=round(time.time() - started, 2))
        write_json(d / "error.json", record)
        return record
    finally:
        lock.release()


def _mafft(cfg: dict, fasta: Path, handle) -> float:
    started = time.time()
    proc = subprocess.run([tool(cfg, "mafft"), "--auto", "--thread", "1", "--quiet", str(fasta)], check=False, stdout=handle,
                          stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"mafft exited {proc.returncode}: {(proc.stderr or '').strip()[-300:]}")
    return time.time() - started


_TITLE_CACHE: dict[str, dict[str, str]] = {}


def _titles(cfg: dict, dbkind: str) -> dict[str, str] | None:
    if dbkind != "pdb":
        return None
    for db in cfg["databases"].values():
        if db.get("kind") == "pdb" and db.get("titles"):
            path = db["titles"]
            if path not in _TITLE_CACHE:
                _TITLE_CACHE[path] = hhr.read_titles(path)
            return _TITLE_CACHE[path]
    return None


def _cell(value) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value).replace("\t", " ").replace("\n", " ")


def select_shard(items: list, shard: int, n_shards: int) -> list:
    if not 0 <= shard < n_shards:
        raise ValueError(f"--shard {shard} outside 0..{n_shards - 1}")
    return [item for i, item in enumerate(items) if i % n_shards == shard]


def cmd_families(ns: argparse.Namespace) -> int:
    cfg = load_config(ns.config)
    work = Path(ns.work).resolve()
    families = read_families(work)
    # biggest first, so the slowest families start early and the shards stay balanced
    families.sort(key=lambda f: (-f["n_sequences"], -f["rep_length"], f["family_index"]))
    if ns.only:
        wanted = set(Path(ns.only).read_text().split())
        families = [f for f in families if f["family_key"] in wanted or f["family_id"] in wanted]
    mine = select_shard(families, ns.shard, ns.n_shards)
    if ns.limit:
        mine = mine[: ns.limit]
    threads = ns.threads or allocated_cpus()
    per = max(1, min(ns.threads_per_worker, threads))
    workers = max(1, threads // per)
    if hhr.MSA_ENRICHED in cfg["msa_kinds"]:
        database(cfg, "uniref")
    for target in cfg["family_targets"]:
        database(cfg, target)
    (work / "logs").mkdir(exist_ok=True)
    started = time.time()
    print(f"families shard {ns.shard}/{ns.n_shards}: {len(mine)} families, {workers} workers x {per} threads, "
          f"msa_kinds {','.join(cfg['msa_kinds'])}, targets {','.join(cfg['family_targets'])}, "
          f"config digest {family_digest(cfg)}", flush=True)
    counts: dict[str, int] = {}
    jobs = [(str(work), f, cfg, per, ns.force, ns.keep_intermediates, ns.stale_lock_hours) for f in mine]
    with ProcessPoolExecutor(workers) as pool:
        for record in pool.map(run_family, jobs):
            counts[record["status"]] = counts.get(record["status"], 0) + 1
            if record["status"] in ("ok", "error"):
                print(f"{record['family_key']} {record['status']} {record.get('seconds', 0):.1f}s"
                      + (f" {record['error']}" if record["status"] == "error" else ""), flush=True)
    elapsed = time.time() - started
    summary = {"shard": ns.shard, "n_shards": ns.n_shards, "families": len(mine), "counts": counts,
               "workers": workers, "threads_per_worker": per, "cpus": threads, "seconds": round(elapsed, 1),
               "host": socket.gethostname(), "job": os.environ.get("SLURM_JOB_ID"),
               "config_digest": family_digest(cfg), "versions": versions(cfg),
               "finished": time.strftime("%Y-%m-%dT%H:%M:%S")}
    write_json(work / "logs" / f"families_shard{ns.shard:04d}_of{ns.n_shards:04d}_{os.getpid()}.json", summary)
    print(json.dumps(summary), flush=True)
    return 1 if counts.get("error") else 0


# --------------------------------------------------------------------------- #
# CDD
# --------------------------------------------------------------------------- #
def run_cdd_chunk(job: tuple) -> dict:
    chunk, out_dir, cfg, threads = job
    chunk, out_dir = Path(chunk), Path(out_dir)
    out = out_dir / (chunk.stem + ".tsv")
    done = out_dir / (chunk.stem + ".done")
    if done.exists() and out.exists():
        return {"chunk": chunk.stem, "status": "skipped"}
    db = database(cfg, "cdd")
    started = time.time()
    try:
        run([tool(cfg, "rpsblast"), "-query", str(chunk), "-db", db["prefix"], "-evalue",
             str(cfg["search"]["cdd_evalue"]), "-num_threads", str(threads), "-outfmt", CDD_OUTFMT,
             "-out", str(out) + ".tmp"])
        os.replace(str(out) + ".tmp", out)
        with open(chunk) as handle:
            n_seq = sum(1 for line in handle if line.startswith(">"))
        record = {"chunk": chunk.stem, "status": "ok", "sequences": n_seq, "seconds": round(time.time() - started, 2),
                  "threads": threads, "host": socket.gethostname(), "job": os.environ.get("SLURM_JOB_ID"),
                  "database": {k: db.get(k) for k in ("prefix", "name", "version")},
                  "evalue": cfg["search"]["cdd_evalue"]}
        write_json(done, record)
        return record
    except Exception as exc:
        return {"chunk": chunk.stem, "status": "error", "error": f"{type(exc).__name__}: {exc}"[:600]}


def cmd_cdd(ns: argparse.Namespace) -> int:
    cfg = load_config(ns.config)
    work = Path(ns.work).resolve()
    chunks = sorted((work / "cdd" / "chunks").glob("chunk_*.faa"))
    if not chunks:
        raise SystemExit(f"{work}/cdd/chunks holds no chunk_*.faa (run `sharur remote-homology export-cdd`)")
    mine = select_shard(chunks, ns.shard, ns.n_shards)
    if ns.limit:
        mine = mine[: ns.limit]
    out_dir = work / "cdd" / "out"
    out_dir.mkdir(parents=True, exist_ok=True)
    threads = ns.threads or allocated_cpus()
    per = max(1, min(ns.threads_per_worker, threads))
    workers = max(1, threads // per)
    database(cfg, "cdd")
    started = time.time()
    print(f"cdd shard {ns.shard}/{ns.n_shards}: {len(mine)} chunks, {workers} workers x {per} threads", flush=True)
    counts: dict[str, int] = {}
    with ProcessPoolExecutor(workers) as pool:
        for record in pool.map(run_cdd_chunk, [(str(c), str(out_dir), cfg, per) for c in mine]):
            counts[record["status"]] = counts.get(record["status"], 0) + 1
            print(json.dumps(record), flush=True)
    summary = {"shard": ns.shard, "n_shards": ns.n_shards, "chunks": len(mine), "counts": counts,
               "workers": workers, "threads_per_worker": per, "cpus": threads,
               "seconds": round(time.time() - started, 1), "host": socket.gethostname(),
               "job": os.environ.get("SLURM_JOB_ID"), "versions": {"rpsblast": versions(cfg)["rpsblast"]}}
    (work / "logs").mkdir(exist_ok=True)
    write_json(work / "logs" / f"cdd_shard{ns.shard:04d}_of{ns.n_shards:04d}_{os.getpid()}.json", summary)
    print(json.dumps(summary), flush=True)
    return 1 if counts.get("error") else 0


# --------------------------------------------------------------------------- #
# Command line
# --------------------------------------------------------------------------- #
def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("check", help="tools and databases of a configuration")
    p.add_argument("--config", required=True)
    for name in ("families", "cdd"):
        p = sub.add_parser(name)
        p.add_argument("--work", required=True, help="work directory from `sharur remote-homology export`")
        p.add_argument("--config", required=True)
        p.add_argument("--shard", type=int, default=0)
        p.add_argument("--n-shards", type=int, default=1)
        p.add_argument("--threads", type=int, default=None, help="default: $SLURM_CPUS_ON_NODE or all CPUs")
        p.add_argument("--threads-per-worker", type=int, default=4)
        p.add_argument("--limit", type=int, default=None, help="at most this many items of the shard (pilots)")
        if name == "families":
            p.add_argument("--only", default=None, help="file of family keys or ids to run (others skipped)")
            p.add_argument("--force", action="store_true", help="redo families that are already done")
            p.add_argument("--keep-intermediates", action="store_true",
                           help="keep aln.fa and enr.a3m and leave .hhr files uncompressed")
            p.add_argument("--stale-lock-hours", type=float, default=6.0)
    ns = parser.parse_args(argv)
    try:
        if ns.command == "check":
            print(json.dumps(check(load_config(ns.config)), indent=2))
            return 0
        return cmd_families(ns) if ns.command == "families" else cmd_cdd(ns)
    except ConfigError as exc:
        print(f"configuration: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
