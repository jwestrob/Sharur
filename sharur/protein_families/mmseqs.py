#!/usr/bin/env python3
"""MMseqs2 commands behind Sharur protein families (standard library only).

This module runs on any host with Python 3.10+ and an ``mmseqs`` binary, so a
cluster job can execute it as a plain script without a Sharur installation::

    python3 sharur/protein_families/mmseqs.py cluster --fasta proteins.faa \\
        --out-dir clustering --levels 0.5,0.3 --threads "$SLURM_CPUS_ON_NODE"

Clustering is hierarchical. Level 0 clusters every input protein at the first
identity threshold; each further level clusters the representatives of the
level before it (the UniRef-style cascade), so every level-0 cluster lies
inside exactly one level-1 cluster. Sharur uses two levels: subfamilies
(``--levels SUBFAMILY,FAMILY``) and families. A single level (``--levels 0.3``)
is a direct clustering of all proteins.

Every level runs ``mmseqs cluster`` with the same alignment criteria: minimum
sequence identity over the alignment (``--alignment-mode 3``, ``--seq-id-mode
0``), coverage of both sequences (``-c 0.8 --cov-mode 0``), sensitivity 6,
50 prefilter hits per query, greedy set cover (``--cluster-mode 0``) and
``--cluster-reassign 1``. Outputs per level: ``level<i>.tsv.gz`` with
``representative<TAB>member`` rows (members include the representative), plus
``clustering.json`` with the exact commands, parameters, MMseqs2 version,
counts and timings. Temporary MMseqs2 files are removed when a level finishes.

``search`` assigns query proteins to the best representative that meets one
level's criteria; :mod:`sharur.protein_families.store` uses it to place the
proteins of newly added genomes into existing families.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import subprocess
import sys
import time
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path


# Defaults: families at 30% identity and subfamilies at a tighter identity, both
# with 80% coverage of both sequences. docs/guides/protein-families.md describes
# how the subfamily threshold was chosen.
FAMILY_MIN_SEQ_ID = 0.3
SUBFAMILY_MIN_SEQ_ID = 0.5
COVERAGE = 0.8
COV_MODE = 0
SENSITIVITY = 6.0
MAX_SEQS = 50
CLUSTER_MODE = 0
CLUSTER_REASSIGN = 1
EVALUE = 1e-3
# Assignment searches representatives only, which are far more diverse than the
# full set, so the prefilter keeps more candidates per query than clustering.
SEARCH_MAX_SEQS = 300


@dataclass(frozen=True)
class LevelParams:
    """Alignment and clustering criteria of one level."""

    min_seq_id: float
    coverage: float = COVERAGE
    cov_mode: int = COV_MODE
    sensitivity: float = SENSITIVITY
    max_seqs: int = MAX_SEQS
    cluster_mode: int = CLUSTER_MODE
    cluster_reassign: int = CLUSTER_REASSIGN
    evalue: float = EVALUE

    def criteria_args(self) -> list[str]:
        return ["--min-seq-id", _num(self.min_seq_id), "-c", _num(self.coverage),
                "--cov-mode", str(self.cov_mode), "--alignment-mode", "3", "--seq-id-mode", "0",
                "-e", _num(self.evalue), "-s", _num(self.sensitivity)]

    def cluster_args(self) -> list[str]:
        return [*self.criteria_args(), "--max-seqs", str(self.max_seqs),
                "--cluster-mode", str(self.cluster_mode),
                "--cluster-reassign", str(self.cluster_reassign)]

    def search_args(self, max_seqs: int = SEARCH_MAX_SEQS) -> list[str]:
        return [*self.criteria_args(), "--max-seqs", str(max_seqs)]

    def to_dict(self) -> dict:
        return asdict(self)


def default_levels(subfamily: float = SUBFAMILY_MIN_SEQ_ID, family: float = FAMILY_MIN_SEQ_ID,
                   **overrides) -> list[LevelParams]:
    """Subfamily then family level, sharing every criterion except identity."""
    if not 0 < family < subfamily <= 1:
        raise ValueError(f"need 0 < family identity ({family}) < subfamily identity ({subfamily}) <= 1")
    base = LevelParams(min_seq_id=subfamily, **overrides)
    return [base, replace(base, min_seq_id=family)]


def _num(value: float) -> str:
    return f"{value:g}"


# --------------------------------------------------------------------------- #
# Environment
# --------------------------------------------------------------------------- #


def allocated_cpus() -> int:
    """CPUs this process may use: the Slurm allocation, else the host count."""
    return int(os.environ.get("SLURM_CPUS_ON_NODE") or os.cpu_count() or 1)


def allocated_memory_bytes() -> int | None:
    """Memory of the Slurm allocation, else physical memory, else None."""
    per_node = os.environ.get("SLURM_MEM_PER_NODE")
    if per_node:
        return int(per_node) * 1024 ** 2
    per_cpu = os.environ.get("SLURM_MEM_PER_CPU")
    if per_cpu:
        return int(per_cpu) * 1024 ** 2 * allocated_cpus()
    try:
        return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    except (ValueError, OSError, AttributeError):
        return None


def default_split_memory_limit(fraction: float = 0.7) -> str | None:
    """MMseqs2 ``--split-memory-limit`` at ``fraction`` of the allocation, in whole GB."""
    total = allocated_memory_bytes()
    if not total:
        return None
    return f"{max(1, int(total * fraction / 1024 ** 3))}G"


def mmseqs_version(binary: str = "mmseqs") -> str:
    out = subprocess.run([binary, "version"], check=True, capture_output=True, text=True)
    return out.stdout.strip()


def resolve_binary(binary: str | None = None) -> str:
    candidate = binary or os.environ.get("SHARUR_MMSEQS") or "mmseqs"
    found = shutil.which(candidate)
    if not found:
        raise FileNotFoundError(
            f"MMseqs2 executable {candidate!r} not found; install MMseqs2 or set SHARUR_MMSEQS")
    return found


# --------------------------------------------------------------------------- #
# Running
# --------------------------------------------------------------------------- #


@dataclass
class _Runner:
    binary: str
    log_path: Path | None = None
    commands: list[list[str]] = field(default_factory=list)

    def __call__(self, *args: str) -> None:
        command = [self.binary, *args]
        self.commands.append(command)
        print("+ " + " ".join(command), flush=True)
        if self.log_path is None:
            subprocess.run(command, check=True)
            return
        with open(self.log_path, "a") as log:
            log.write("+ " + " ".join(command) + "\n")
            log.flush()
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
        if result.returncode != 0:
            raise RuntimeError(f"{' '.join(command[:2])} exited {result.returncode}; see {self.log_path}")


def _gzip_tsv(src: Path, dst: Path) -> int:
    rows = 0
    with open(src, "rb") as handle, gzip.open(dst, "wb", compresslevel=6) as out:
        for line in handle:
            out.write(line)
            rows += 1
    src.unlink()
    return rows


def _remove_db(prefix: Path) -> None:
    """Remove an MMseqs2 database and its companion files."""
    for path in prefix.parent.glob(prefix.name + "*"):
        if path.name == prefix.name or path.name[len(prefix.name)] in "._":
            if path.is_dir():
                shutil.rmtree(path)
            else:
                path.unlink()


def cluster_levels(
    out_dir: Path,
    levels: list[LevelParams],
    *,
    fasta: Path | None = None,
    seqdb: Path | None = None,
    threads: int | None = None,
    binary: str | None = None,
    split_memory_limit: str | None = None,
    keep_databases: bool = False,
) -> dict:
    """Cluster ``fasta`` (or a prebuilt ``seqdb``) level by level; return the manifest.

    Level ``i`` clusters the representatives of level ``i - 1``. Each level writes
    ``level<i>.tsv.gz`` (``representative<TAB>member``) into ``out_dir``.
    """
    if not levels:
        raise ValueError("at least one level is required")
    for upper, lower in zip(levels, levels[1:]):
        if lower.min_seq_id >= upper.min_seq_id:
            raise ValueError("each level needs a lower identity threshold than the level before it")
    if (fasta is None) == (seqdb is None):
        raise ValueError("give exactly one of fasta or seqdb")
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    binary = resolve_binary(binary)
    threads = threads or allocated_cpus()
    split = split_memory_limit or default_split_memory_limit()
    run = _Runner(binary, out_dir / "mmseqs.log")
    started = time.time()
    manifest: dict = {
        "tool": "mmseqs2",
        "tool_version": mmseqs_version(binary),
        "threads": threads,
        "split_memory_limit": split,
        "input": str(fasta) if fasta else str(seqdb),
        "levels": [],
    }
    db_dir = out_dir / "db"
    db_dir.mkdir(exist_ok=True)
    if seqdb is None:
        seqdb = db_dir / "seqDB"
        t0 = time.time()
        run("createdb", str(fasta), str(seqdb), "--dbtype", "1", "--threads", str(threads))
        manifest["createdb_seconds"] = round(time.time() - t0, 1)
    current = Path(seqdb)
    created: list[Path] = []
    for index, level in enumerate(levels):
        t0 = time.time()
        clu = db_dir / f"level{index}_clu"
        tmp = out_dir / f"tmp_level{index}"
        _remove_db(clu)
        args = ["cluster", str(current), str(clu), str(tmp), *level.cluster_args(),
                "--threads", str(threads), "--remove-tmp-files", "1"]
        if split:
            args += ["--split-memory-limit", split]
        run(*args)
        tsv = out_dir / f"level{index}.tsv"
        run("createtsv", str(current), str(current), str(clu), str(tsv), "--threads", str(threads))
        rows = _gzip_tsv(tsv, out_dir / f"level{index}.tsv.gz")
        shutil.rmtree(tmp, ignore_errors=True)
        n_clusters = _count_index_entries(clu)
        record = {"level": index, "params": level.to_dict(), "input_sequences": _count_index_entries(current),
                  "clusters": n_clusters, "members": rows, "seconds": round(time.time() - t0, 1),
                  "tsv": f"level{index}.tsv.gz"}
        manifest["levels"].append(record)
        print(f"level {index} (min-seq-id {level.min_seq_id:g}): {record['input_sequences']:,} sequences -> "
              f"{n_clusters:,} clusters in {record['seconds']:.0f} s", flush=True)
        created.append(clu)
        if index + 1 < len(levels):
            reps = db_dir / f"level{index}_reps"
            _remove_db(reps)
            run("createsubdb", str(clu), str(current), str(reps))
            run("createsubdb", str(clu), str(current) + "_h", str(reps) + "_h")
            created.append(reps)
            current = reps
    manifest["commands"] = [c[1:] for c in run.commands]
    manifest["seconds"] = round(time.time() - started, 1)
    if not keep_databases:
        shutil.rmtree(db_dir, ignore_errors=True)
    (out_dir / "clustering.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def _count_index_entries(prefix: Path) -> int:
    index = Path(str(prefix) + ".index")
    with open(index, "rb") as handle:
        return sum(1 for _ in handle)


def search_best(
    query_fasta: Path,
    target_fasta: Path,
    out_dir: Path,
    level: LevelParams,
    *,
    threads: int | None = None,
    binary: str | None = None,
    max_seqs: int = SEARCH_MAX_SEQS,
    split_memory_limit: str | None = None,
) -> dict[str, tuple[str, float, float, float, float]]:
    """Best target per query among hits that meet ``level``'s identity and coverage.

    Returns ``{query: (target, fident, qcov, tcov, bits)}``; the best hit has the
    highest bit score (ties: higher identity, then target ID).
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    binary = resolve_binary(binary)
    threads = threads or allocated_cpus()
    split = split_memory_limit or default_split_memory_limit()
    run = _Runner(binary, out_dir / "mmseqs.log")
    qdb, tdb, res, tmp = (out_dir / n for n in ("queryDB", "targetDB", "resultDB", "tmp"))
    for db in (qdb, tdb, res):
        _remove_db(db)
    run("createdb", str(query_fasta), str(qdb), "--dbtype", "1", "--threads", str(threads))
    run("createdb", str(target_fasta), str(tdb), "--dbtype", "1", "--threads", str(threads))
    args = ["search", str(qdb), str(tdb), str(res), str(tmp), *level.search_args(max_seqs),
            "--threads", str(threads), "--remove-tmp-files", "1"]
    if split:
        args += ["--split-memory-limit", split]
    run(*args)
    tsv = out_dir / "hits.tsv"
    run("convertalis", str(qdb), str(tdb), str(res), str(tsv),
        "--format-output", "query,target,fident,qcov,tcov,bits", "--threads", str(threads))
    best: dict[str, tuple[str, float, float, float, float]] = {}
    with open(tsv) as handle:
        for line in handle:
            query, target, fident, qcov, tcov, bits = line.rstrip("\n").split("\t")
            hit = (target, float(fident), float(qcov), float(tcov), float(bits))
            if hit[1] < level.min_seq_id or _covered(level, hit[2], hit[3]) < level.coverage:
                continue
            kept = best.get(query)
            if kept is None or (hit[4], hit[1], kept[0]) > (kept[4], kept[1], hit[0]):
                best[query] = hit
    for db in (qdb, tdb, res):
        _remove_db(db)
    shutil.rmtree(tmp, ignore_errors=True)
    tsv.unlink(missing_ok=True)
    return best


def _covered(level: LevelParams, qcov: float, tcov: float) -> float:
    """Coverage that ``--cov-mode`` tests: 0 both sequences, 1 the target, 2 the query."""
    return {0: min(qcov, tcov), 1: tcov, 2: qcov}.get(level.cov_mode, min(qcov, tcov))


def write_hits(best: dict[str, tuple[str, float, float, float, float]], path: Path) -> None:
    with open(path, "w") as out:
        out.write("query\ttarget\tfident\tqcov\ttcov\tbits\n")
        for query in sorted(best):
            target, fident, qcov, tcov, bits = best[query]
            out.write(f"{query}\t{target}\t{fident:g}\t{qcov:g}\t{tcov:g}\t{bits:g}\n")


def read_level(path: Path):
    """Yield ``(representative, member)`` rows of a ``level<i>.tsv.gz`` file."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            rep, member = line.rstrip("\n").split("\t")
            yield rep, member


# --------------------------------------------------------------------------- #
# Command line
# --------------------------------------------------------------------------- #


def _parse_levels(text: str) -> list[float]:
    values = [float(v) for v in text.split(",") if v.strip()]
    if not values:
        raise argparse.ArgumentTypeError("give at least one identity threshold")
    return values


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    cl = sub.add_parser("cluster", help="Hierarchical MMseqs2 clustering of a protein FASTA.")
    src = cl.add_mutually_exclusive_group(required=True)
    src.add_argument("--fasta", type=Path)
    src.add_argument("--seqdb", type=Path, help="Prebuilt MMseqs2 sequence database (reused, left in place).")
    cl.add_argument("--out-dir", type=Path, required=True)
    cl.add_argument("--levels", type=_parse_levels,
                    default=[SUBFAMILY_MIN_SEQ_ID, FAMILY_MIN_SEQ_ID],
                    help="Comma-separated identity thresholds, strictly decreasing "
                         f"(default {SUBFAMILY_MIN_SEQ_ID:g},{FAMILY_MIN_SEQ_ID:g}: subfamily, family).")
    cl.add_argument("--coverage", type=float, default=COVERAGE)
    cl.add_argument("--cov-mode", type=int, default=COV_MODE)
    cl.add_argument("--sensitivity", type=float, default=SENSITIVITY)
    cl.add_argument("--max-seqs", type=int, default=MAX_SEQS)
    cl.add_argument("--cluster-mode", type=int, default=CLUSTER_MODE)
    cl.add_argument("--cluster-reassign", type=int, default=CLUSTER_REASSIGN)
    cl.add_argument("--threads", type=int, default=None, help="Default: $SLURM_CPUS_ON_NODE or all CPUs.")
    cl.add_argument("--split-memory-limit", default=None,
                    help="MMseqs2 --split-memory-limit (default: 70%% of the allocation).")
    cl.add_argument("--mmseqs", default=None, help="MMseqs2 executable (default: $SHARUR_MMSEQS or mmseqs).")
    cl.add_argument("--keep-databases", action="store_true")
    se = sub.add_parser("search", help="Best representative per query that meets one level's criteria.")
    se.add_argument("--query", type=Path, required=True, help="Query protein FASTA.")
    se.add_argument("--target", type=Path, required=True, help="Representative protein FASTA.")
    se.add_argument("--out", type=Path, required=True, help="TSV of best hits (query, target, fident, qcov, tcov, bits).")
    se.add_argument("--work-dir", type=Path, default=None, help="Scratch (default: beside --out).")
    se.add_argument("--min-seq-id", type=float, default=FAMILY_MIN_SEQ_ID)
    se.add_argument("--coverage", type=float, default=COVERAGE)
    se.add_argument("--cov-mode", type=int, default=COV_MODE, help="0 both sequences, 1 target, 2 query.")
    se.add_argument("--sensitivity", type=float, default=SENSITIVITY)
    se.add_argument("--max-seqs", type=int, default=SEARCH_MAX_SEQS)
    se.add_argument("--threads", type=int, default=None)
    se.add_argument("--split-memory-limit", default=None)
    se.add_argument("--mmseqs", default=None)
    args = parser.parse_args(argv)
    if args.command == "search":
        level = LevelParams(min_seq_id=args.min_seq_id, coverage=args.coverage, cov_mode=args.cov_mode,
                            sensitivity=args.sensitivity)
        work = args.work_dir or args.out.with_name(args.out.name + ".work")
        best = search_best(args.query, args.target, work, level, threads=args.threads, binary=args.mmseqs,
                           max_seqs=args.max_seqs, split_memory_limit=args.split_memory_limit)
        write_hits(best, args.out)
        shutil.rmtree(work, ignore_errors=True)
        print(f"{len(best):,} queries with a hit -> {args.out}")
    if args.command == "cluster":
        levels = [LevelParams(min_seq_id=v, coverage=args.coverage, cov_mode=args.cov_mode,
                              sensitivity=args.sensitivity, max_seqs=args.max_seqs,
                              cluster_mode=args.cluster_mode, cluster_reassign=args.cluster_reassign)
                  for v in args.levels]
        manifest = cluster_levels(args.out_dir, levels, fasta=args.fasta, seqdb=args.seqdb,
                                  threads=args.threads, binary=args.mmseqs,
                                  split_memory_limit=args.split_memory_limit,
                                  keep_databases=args.keep_databases)
        print(json.dumps({k: v for k, v in manifest.items() if k != "commands"}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
