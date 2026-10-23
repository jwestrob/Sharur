#!/usr/bin/env python3
"""
CAZyme Classification Pipeline (dbCAN 3-tool consensus)

Recapitulates the standard dbCAN annotation pipeline:
  Tool 1: DIAMOND blastp vs CAZy.dmnd        (e-value ≤ 1e-102)
  Tool 2: Astra/HMMER vs dbCAN.hmm           (i-evalue ≤ 1e-15, HMM coverage ≥ 0.35)
  Tool 3: Astra/HMMER vs dbCAN-sub.hmm       (i-evalue ≤ 1e-15, HMM coverage ≥ 0.35)
  Consensus: retain genes called by ≥ 2 of 3 tools

HMMER searches are delegated to Aksha (aksha search --hmm_in), which handles
parallelization and PyHMMER optimization. Only DIAMOND is run directly, in its
default sensitivity mode with -k 1, as run_dbcan runs it.

If dbCAN-sub.hmm is not available, falls back to 2-tool mode (DIAMOND + dbCAN HMM)
and requires both tools to agree.

Usage:
    python scripts/classify_cazymes.py --db data/my_dataset/sharur.duckdb --threads 8

    # Compute on a cluster against a read-only database, load elsewhere:
    python scripts/classify_cazymes.py --db DB --threads N --results-only --out cazy.tsv
    python scripts/classify_cazymes.py --db DB --load-results cazy.tsv

    # A protein FASTA (file or directory of .faa) in place of a database:
    python scripts/classify_cazymes.py --proteins all.faa --threads N --results-only --out cazy.tsv

Sharded runs (Slurm arrays; scripts/submit_cazy_shards.sh submits all three):
    --stage search --shard-index i --n-shards N --work-dir W
        dbCAN.hmm profile part i against every protein; shard 0 also runs DIAMOND on
        every protein.
    --stage sub --shard-index i --n-shards N --work-dir W     (after every search shard)
        dbCAN-sub profile part i against the candidates from all search shards.
    --stage merge --n-shards N --work-dir W --out cazy.tsv    (after every sub shard)
        Consensus over the combined hits.

    Neither search splits the proteins. HMMER i-evalues scale with the number of
    target sequences in a search (Z), so the HMM searches split the profiles and every
    part searches the single-process target set. DIAMOND's top target for a protein in
    default mode depends on the other proteins in its run, so DIAMOND runs once on
    every protein, as in the single-process run. The merged consensus table equals
    the single-process table.

References:
    Zheng J, et al. (2023) dbCAN3: automated carbohydrate-active enzyme and
    substrate annotation. Nucleic Acids Res 51(W1):W115-W121.
"""

import argparse
import heapq
import json
import mmap
import os
import shutil
import subprocess
import tempfile
import re
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from collections import defaultdict

import time as _time

import duckdb
import numpy as np
import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[0].parent))
from sharur.predicates.mappings.cazy_map import get_predicates_for_cazy, get_cazy_class


def _log(msg: str) -> None:
    """Print with flush so output appears immediately under redirection."""
    print(msg, flush=True)

# ── dbCAN standard thresholds ──────────────────────────────────────────────────
# run_dbcan and the dbCAN3 server defaults: DIAMOND E < 1e-102; HMMER (dbCAN and
# dbCAN-sub) E < 1e-15 with HMM coverage > 0.35. A 1e-18 DIAMOND cut passed 11% of
# the SR-VP archaeal proteome (224,612 / 1,976,915) against 1.35% for dbCAN.hmm.
DIAMOND_EVALUE = 1e-102
HMMER_IEVALUE  = 1e-15
HMMER_COVERAGE = 0.35
MIN_TOOLS      = 2

CAZY_PREFIX_RE = re.compile(r'^(GH|GT|PL|CE|CBM|AA)\d+')


# ── Database discovery ─────────────────────────────────────────────────────────

def get_dbcan_path(override: str | Path | None = None) -> Path:
    """Find the dbCAN database directory (must contain CAZy.dmnd + dbCAN.hmm).

    ``override`` or ``$SHARUR_DBCAN_DIR`` names it explicitly; otherwise the
    standard locations are searched.
    """
    explicit = override or os.environ.get("SHARUR_DBCAN_DIR")
    if explicit:
        path = Path(explicit)
        if (path / "CAZy.dmnd").exists() and (path / "dbCAN.hmm").exists():
            return path
        raise FileNotFoundError(f"dbCAN directory {path} lacks CAZy.dmnd + dbCAN.hmm")
    candidates = [
        Path(__file__).parent.parent / "data/dbcan_db",
        Path(__file__).parent.parent / "data/reference/dbcan",
        Path("data/dbcan_db"),
        Path("data/reference/dbcan"),
        Path.home() / ".sharur/dbcan",
    ]
    for path in candidates:
        if (path / "CAZy.dmnd").exists() and (path / "dbCAN.hmm").exists():
            return path
    raise FileNotFoundError(
        "dbCAN databases not found. Need at least CAZy.dmnd + dbCAN.hmm.\n"
        "Download: aws s3 cp s3://dbcan/db_v5-2_9-13-2025/ data/dbcan_db/ "
        "--no-sign-request --recursive"
    )


def dbcan_hmm_files(dbcan: Path) -> tuple[Path, Path, bool]:
    """(dbCAN.hmm, dbCAN-sub HMM, whether the sub-family file is usable)."""
    # dbCAN releases name the sub-family profiles dbCAN_sub.hmm; older docs used dbCAN-sub.hmm.
    sub = next((dbcan / n for n in ("dbCAN-sub.hmm", "dbCAN_sub.hmm") if (dbcan / n).exists()),
               dbcan / "dbCAN-sub.hmm")
    return dbcan / "dbCAN.hmm", sub, sub.exists() and sub.stat().st_size > 0


# ── Protein extraction ─────────────────────────────────────────────────────────

def count_fasta_records(prot_dir: str) -> int:
    """Header lines across the directory's .faa files."""
    n = 0
    for faa in sorted(Path(prot_dir).glob("*.faa")):
        with open(faa, "rb") as fh:
            n += sum(1 for line in fh if line.startswith(b">"))
    return n


def resolve_proteins(db_path: str | None, proteins: str | None = None,
                     verbose: bool = True) -> tuple:
    """(prot_dir, n_proteins, is_temp) from an explicit FASTA, or else from the database.

    ``proteins`` names a protein FASTA file or a directory of .faa files and needs
    no database.
    """
    if proteins is None:
        return find_or_extract_proteins(db_path, verbose)
    path = Path(proteins).resolve()
    if path.is_dir():
        prot_dir, is_temp = str(path), False
    elif path.is_file():
        prot_dir, is_temp = tempfile.mkdtemp(prefix="cazyme_prot_"), True
        (Path(prot_dir) / "proteins.faa").symlink_to(path)
    else:
        raise FileNotFoundError(f"protein input {path} does not exist")
    n = count_fasta_records(prot_dir)
    _log(f"[cazyme] Using protein input {path} ({n:,} records)")
    return prot_dir, n, is_temp


def find_or_extract_proteins(db_path: str, verbose: bool = True) -> tuple:
    """Find existing source .faa or extract proteins from database.

    Astra expects --prot_in <directory> containing .faa files.

    Returns (prot_dir, n_proteins, is_temp).
    - prot_dir: Path to directory containing .faa file(s)
    - n_proteins: total protein count
    - is_temp: whether prot_dir is a temp directory that should be cleaned up
    """
    db_dir = Path(db_path).parent

    # Prefer existing FAA files over extracting from DB
    # Check stage03 all_protein_symlinks first (fastest — no extraction needed)
    stage03 = db_dir / "stage03_prodigal" / "genomes" / "all_protein_symlinks"
    if stage03.exists():
        faas = list(stage03.glob("*.faa"))
        if faas:
            db = duckdb.connect(db_path, read_only=True)
            n = db.execute("SELECT COUNT(*) FROM proteins").fetchone()[0]
            db.close()
            _log(f"[cazyme] Using stage03 FAA directory ({len(faas)} files, {n:,} proteins)")
            return str(stage03), n, False

    # Check for single consolidated FAA
    source_faa = db_dir / "source" / "proteins_db_ids.faa"
    if source_faa.exists():
        db = duckdb.connect(db_path, read_only=True)
        n = db.execute("SELECT COUNT(*) FROM proteins").fetchone()[0]
        db.close()
        _log(f"[cazyme] Using existing protein file: {source_faa}")
        tmp_dir = tempfile.mkdtemp(prefix="cazyme_prot_")
        (Path(tmp_dir) / "proteins.faa").symlink_to(source_faa.resolve())
        return tmp_dir, n, True

    # Fall back to extracting from DB
    _log(f"[cazyme] Checking DB for protein sequences...")
    db = duckdb.connect(db_path, read_only=True)
    n = db.execute(
        "SELECT COUNT(*) FROM proteins WHERE sequence IS NOT NULL AND sequence != ''"
    ).fetchone()[0]
    db.close()
    _log(f"[cazyme] {n:,} proteins with sequences in DB")

    if n == 0:
        return None, 0, False

    # Extract from database
    _log(f"[cazyme] Extracting {n:,} protein sequences from database...")

    tmp_dir = tempfile.mkdtemp(prefix="cazyme_prot_")
    fasta_path = Path(tmp_dir) / "proteins.faa"

    # One streamed scan: paging with LIMIT/OFFSET and no ORDER BY can skip or
    # repeat rows under a parallel scan.
    db = duckdb.connect(db_path, read_only=True, config={"threads": 4})
    written = 0
    cursor = db.execute(
        "SELECT protein_id, sequence FROM proteins "
        "WHERE sequence IS NOT NULL AND sequence != '' ORDER BY protein_id"
    )
    with open(fasta_path, 'w') as f:
        while rows := cursor.fetchmany(100_000):
            for pid, seq in rows:
                f.write(f">{pid}\n{seq}\n")
                written += 1
            if written % 500_000 < 100_000:
                _log(f"[cazyme]   extracted {written:,} sequences...")
    db.close()

    _log(f"[cazyme] Wrote {written:,} sequences")
    return tmp_dir, n, True


def _subset_fasta(prot_dir: str, protein_ids: set,
                  verbose: bool = True) -> str:
    """Extract a subset of proteins from prot_dir into a new temp directory.

    Returns path to new temp directory containing the subset .faa.
    Caller is responsible for cleanup.
    """
    faa_files = sorted(Path(prot_dir).glob("*.faa"))
    if not faa_files:
        raise FileNotFoundError(f"No .faa files in {prot_dir}")

    sub_dir = tempfile.mkdtemp(prefix="cazyme_sub_")
    sub_faa = Path(sub_dir) / "subset.faa"
    written = 0
    with open(sub_faa, 'w') as fout:
        for faa in faa_files:  # every file: a stage03 directory holds one FAA per genome
            with open(faa) as fin:
                write_this = False
                for line in fin:
                    if line.startswith('>'):
                        pid = line[1:].split()[0]
                        write_this = pid in protein_ids
                        if write_this:
                            fout.write(line)
                            written += 1
                    elif write_this:
                        fout.write(line)

    if verbose:
        print(f"  Subset for dbCAN-sub: {written:,} / {len(protein_ids):,} "
              f"candidate proteins")
    return sub_dir


# ── Tool 1: DIAMOND ────────────────────────────────────────────────────────────

def parse_cazy_families(subject_id: str) -> list[str]:
    """Extract CAZy family names from a dbCAN DIAMOND subject ID.

    Format: ACCESSION|FAMILY[_SUBFAM] or ACCESSION|FAM1|FAM2 (multi-domain).
    Returns base families (e.g., ["GH13"] or ["CE4", "GT2"]).
    """
    families = []
    for part in subject_id.split("|")[1:]:
        part = part.strip()
        if part and CAZY_PREFIX_RE.match(part):
            families.append(re.split(r'_', part)[0])
    return families


DIAMOND_CACHE = "cazyme_diamond_hits.tsv"


def load_diamond_cache(db_path: str,
                       verbose: bool = True) -> dict[str, dict[str, tuple]] | None:
    """Load cached DIAMOND results from cazyme_diamond_hits.tsv beside the database.

    The cache has its own file: the consensus table (cazyme_classification.tsv)
    merges HMM families and must never be read back as DIAMOND evidence.

    Returns {protein_id: {family_class: (evalue, bitscore)}} or None if no cache.
    """
    return _read_diamond_cache(Path(db_path).parent / DIAMOND_CACHE, verbose)


def _read_diamond_cache(tsv: Path, verbose: bool = True) -> dict[str, dict[str, tuple]] | None:
    if not tsv.exists():
        return None

    try:
        df = pd.read_csv(tsv, sep='\t')
    except Exception:
        return None

    if df.empty or 'protein_id' not in df.columns:
        return None

    if verbose:
        print(f"Loading cached DIAMOND results from {tsv.name} "
              f"({df['protein_id'].nunique():,} proteins)")

    hits: dict[str, dict[str, tuple]] = defaultdict(dict)
    for _, row in df.iterrows():
        pid = row['protein_id']
        fam = row.get('family_class', '')
        ev = float(row.get('evalue', 1.0))
        bs = float(row.get('bitscore', 0.0))
        if ev > DIAMOND_EVALUE:  # a cache written under a looser cut is refiltered here
            continue
        if fam and CAZY_PREFIX_RE.match(fam):
            if fam not in hits[pid] or ev < hits[pid][fam][0]:
                hits[pid][fam] = (ev, bs)

    if verbose:
        print(f"  DIAMOND (cached): {len(hits):,} proteins")
    return dict(hits) if hits else None


def _save_diamond_cache(hits: dict, data_dir: Path,
                        verbose: bool = True) -> None:
    """Save DIAMOND hits to cazyme_diamond_hits.tsv for future reuse."""
    rows = []
    for pid, fam_dict in hits.items():
        for fam, (ev, bs) in fam_dict.items():
            rows.append({
                'protein_id': pid,
                'family_class': fam,
                'cazy_class': get_cazy_class(fam) or '',
                'evalue': ev,
                'bitscore': bs,
            })
    if rows:
        df = pd.DataFrame(rows).sort_values(['protein_id', 'family_class'], kind='stable')
        out = data_dir / DIAMOND_CACHE
        df.to_csv(out, sep='\t', index=False)
        if verbose:
            print(f"  Saved DIAMOND cache to {out.name}")


def _add_diamond_line(hits: dict, line: str) -> bool:
    """Fold one outfmt-6 line (qseqid sseqid pident evalue bitscore) into ``hits``.

    The first line with a family's lowest E-value wins, so a protein's result
    depends on its own lines alone.
    """
    fields = line.rstrip('\n').split('\t')
    if len(fields) < 5:
        return False
    qid, sid = fields[0], fields[1]
    ev, bs = float(fields[3]), float(fields[4])
    for fam in parse_cazy_families(sid):
        if fam not in hits[qid] or ev < hits[qid][fam][0]:
            hits[qid][fam] = (ev, bs)
    return True


def parse_diamond_tsv(path: Path) -> dict[str, dict[str, tuple]]:
    """{protein_id: {family_class: (evalue, bitscore)}} from saved DIAMOND outfmt-6 output."""
    hits: dict[str, dict[str, tuple]] = defaultdict(dict)
    with open(path) as fh:
        for line in fh:
            _add_diamond_line(hits, line)
    return dict(hits)


def run_diamond(prot_dir: str, dbcan_path: Path, threads: int = 4,
                evalue: float = DIAMOND_EVALUE,
                verbose: bool = True,
                raw_out: Path | None = None) -> dict[str, dict[str, tuple]]:
    """DIAMOND blastp against CAZy.dmnd in default mode, top target only (-k 1).

    ``raw_out`` keeps DIAMOND's outfmt-6 lines.

    Returns {protein_id: {family_class: (evalue, bitscore)}}.
    """
    dmnd = dbcan_path / "CAZy.dmnd"
    # Find .faa file(s) in prot_dir
    faa_files = sorted(Path(prot_dir).glob("*.faa"))
    if not faa_files:
        _log("Warning: No .faa files found for DIAMOND")
        return {}

    # Concatenate FAA files, filtering empty sequences and stop codons
    _log(f"[diamond] Concatenating {len(faa_files)} FAA files (filtering empty seqs + stop codons)...")
    t_cat = _time.time()
    # A private directory: the parent of prot_dir is shared ($TMPDIR, or a dataset's
    # stage03 tree), and two concurrent runs writing one fixed path corrupt each other's query.
    work_dir = Path(tempfile.mkdtemp(prefix="cazyme_diamond_"))
    combined = work_dir / "diamond_combined.faa"
    n_written = 0
    n_skipped = 0

    def flush(header, seq_lines):
        nonlocal n_written, n_skipped
        if header is None:
            return
        if seq_lines:
            fout.write(header)
            fout.writelines(seq_lines)
            n_written += 1
        else:
            n_skipped += 1

    with open(combined, "w") as fout:
        for faa in faa_files:
            with open(faa) as fin:
                header = None
                seq_lines: list[str] = []
                for line in fin:
                    if line.startswith(">"):
                        flush(header, seq_lines)
                        header = line
                        seq_lines = []
                    else:
                        cleaned = line.rstrip("\n").replace("*", "")
                        if cleaned:
                            seq_lines.append(cleaned + "\n")
                flush(header, seq_lines)
    fasta_path = str(combined)
    _log(f"[diamond] Concatenated in {_time.time() - t_cat:.1f}s ({n_written:,} seqs, {n_skipped} empty skipped)")

    _log(f"[diamond] Running DIAMOND blastp against {dmnd.name} ({Path(fasta_path).stat().st_size / 1e6:.0f} MB query)...")

    cmd = [
        "diamond", "blastp",
        "--db", str(dmnd), "--query", fasta_path,
        "--outfmt", "6", "qseqid", "sseqid", "pident", "evalue", "bitscore",
        "--max-target-seqs", "1",
        "--evalue", str(evalue),
        "--threads", str(threads),
        "--tmpdir", str(work_dir),
    ]

    # stderr goes to a file: an unread stderr pipe can fill and stall DIAMOND
    err_handle = tempfile.TemporaryFile(mode="w+")
    try:
        proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=err_handle,
            text=True, bufsize=1,
        )
    except FileNotFoundError:
        shutil.rmtree(work_dir, ignore_errors=True)
        _log("Warning: diamond not found on PATH")
        return {}

    hits: dict[str, dict[str, tuple]] = defaultdict(dict)
    n_lines = 0
    raw = open(raw_out, "w") if raw_out is not None else None
    try:
        for line in proc.stdout:
            if raw is not None:
                raw.write(line)
            if not _add_diamond_line(hits, line):
                continue
            n_lines += 1
            if verbose and n_lines % 500_000 == 0:
                _log(f"  [diamond] Parsed {n_lines:,} hits, {len(hits):,} proteins so far...")

        proc.wait()
    except BaseException:
        proc.kill()
        raise
    finally:
        if raw is not None:
            raw.close()

    err_handle.seek(0)
    stderr_out = err_handle.read()
    err_handle.close()
    shutil.rmtree(work_dir, ignore_errors=True)
    if proc.returncode != 0:
        # fatal: an empty DIAMOND set would silently turn the 3-tool consensus into HMM-only
        raise RuntimeError(f"DIAMOND failed (rc={proc.returncode}): {stderr_out[-500:]}")

    if verbose and stderr_out:
        lines = stderr_out.strip().split('\n')
        for line in [l for l in lines[:-3] if l.startswith("Algorithm")] + lines[-3:]:
            _log(f"  {line}")

    if verbose:
        _log(f"  DIAMOND: {len(hits):,} proteins with CAZy hits ({n_lines:,} hit lines)")
    return dict(hits)


# ── Tools 2 & 3: Astra/HMMER ──────────────────────────────────────────────────

def hmm_lengths(hmm_path: Path) -> dict[str, int]:
    """Profile NAME -> LENG from a HMMER3 text file (for the dbCAN coverage filter)."""
    lengths: dict[str, int] = {}
    name = None
    with open(hmm_path, "rb") as fh, mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ) as mm:
        for m in re.finditer(rb"^(NAME|LENG)[ \t]+(\S+)", mm, re.M):
            if m.group(1) == b"NAME":
                name = m.group(2).decode()
            elif name is not None:
                lengths[name] = int(m.group(2))
                name = None
    return lengths


def hmm_profile_order(hmm_path: Path) -> dict[str, int]:
    """Profile NAME -> position in a HMMER3 text file; names must be unique."""
    with open(hmm_path, "rb") as fh, mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ) as mm:
        names = [m.group(1).decode() for m in re.finditer(rb"^NAME\s+(\S+)", mm, re.M)]
    order = {name: k for k, name in enumerate(names)}
    if len(order) != len(names):
        raise ValueError(f"{hmm_path}: profile names repeat, so shard rows cannot be ordered")
    return order


def write_hmm_part(hmm_path: Path, part: int, n_parts: int, out: Path) -> int:
    """Write profiles part, part + n_parts, part + 2 n_parts, ... of ``hmm_path`` to ``out``.

    Striding spreads each family's consecutive subfamily profiles over every part,
    which balances the search cost. Returns the number of profiles written.
    """
    with open(hmm_path, "rb") as fh, mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ) as mm:
        starts = [m.start() for m in re.finditer(rb"^HMMER3", mm, re.M)]
        n_names = sum(1 for _ in re.finditer(rb"^NAME\s", mm, re.M))
        if not starts or len(starts) != n_names:
            raise ValueError(f"{hmm_path}: {len(starts)} HMMER3 records but {n_names} NAME lines")
        ends = starts[1:] + [len(mm)]
        with open(out, "wb") as fout:
            for k in range(part, len(starts), n_parts):
                fout.write(mm[starts[k]:ends[k]])
    return len(range(part, len(starts), n_parts))


def parse_hmm_family(hmm_name: str) -> str | None:
    """Extract CAZy family_class from an HMM profile name.

    dbCAN.hmm:     "GH13.hmm"                                 → "GH13"
    dbCAN-sub.hmm: "PL1_e18.hmm|PL1_7:82|4.2.2.2:31|PL1:1"   → "PL1"
    """
    name = hmm_name.split('.hmm')[0] if '.hmm' in hmm_name else hmm_name
    name = name.split('|')[0]
    m = CAZY_PREFIX_RE.match(name)
    return m.group(0) if m else None


def _hmm_search_command() -> str:
    """Aksha is the production HMM-search CLI; legacy astra forwards to it."""
    return "aksha" if shutil.which("aksha") else "astra"


def run_astra_hmm(prot_dir: str, label: str,
                  installed_hmm: str | None = None,
                  hmm_path: Path | None = None,
                  threads: int = 4, ievalue: float = HMMER_IEVALUE,
                  verbose: bool = True,
                  keep_tsv: Path | None = None) -> dict[str, dict[str, tuple]]:
    """Run Astra hmmsearch against an HMM database.

    Uses --installed_hmms for Astra-installed databases (e.g., dbCAN),
    or --hmm_in for custom HMM files (e.g., dbCAN-sub.hmm).

    Returns {protein_id: {family_class: (i_evalue, dom_bitscore)}}.
    """
    if verbose:
        name = installed_hmm or (hmm_path.name if hmm_path else "?")
        print(f"Running Astra search against {name}...")

    out_dir = tempfile.mkdtemp(prefix=f"cazyme_{label}_")

    cmd = [_hmm_search_command(), "search", "--prot_in", prot_dir,
           "--outdir", out_dir, "--threads", str(threads)]
    if installed_hmm:
        cmd.extend(["--installed_hmms", installed_hmm])
    elif hmm_path:
        cmd.extend(["--hmm_in", str(hmm_path)])
    else:
        return {}

    result = subprocess.run(cmd, capture_output=True, text=True)

    # fatal: a missing tool would silently shrink the consensus
    if result.returncode != 0:
        shutil.rmtree(out_dir, ignore_errors=True)
        raise RuntimeError(f"HMM search ({label}) failed (rc={result.returncode}): {result.stderr[-500:]}")

    # Find output TSV (name depends on --installed_hmms vs --hmm_in)
    tsvs = list(Path(out_dir).glob("*_hits_df.tsv"))
    if not tsvs:
        shutil.rmtree(out_dir, ignore_errors=True)
        raise RuntimeError(f"HMM search ({label}) wrote no *_hits_df.tsv")

    lengths = hmm_lengths(hmm_path) if hmm_path else None
    hits = _parse_astra_tsv(tsvs[0], ievalue, lengths)
    if keep_tsv is not None:  # keep the raw per-domain table for audit and reuse
        shutil.copy2(tsvs[0], keep_tsv)
    shutil.rmtree(out_dir, ignore_errors=True)

    if verbose:
        print(f"  {label}: {len(hits):,} proteins with hits")
    return hits


def load_astra_tsv_cache(tsv_path: Path, ievalue: float = HMMER_IEVALUE,
                         verbose: bool = True,
                         lengths: dict[str, int] | None = None) -> dict[str, dict[str, tuple]] | None:
    """Load pre-existing Astra hits TSV. Returns None if file doesn't exist."""
    if not tsv_path.exists():
        return None
    if verbose:
        print(f"Loading cached Astra results from {tsv_path.name}...")
    return _parse_astra_tsv(tsv_path, ievalue, lengths)


ASTRA_COLUMNS = {'sequence_id', 'hmm_name', 'i_evalue', 'dom_bitscore', 'hmm_from', 'hmm_to'}


def _astra_best_rows(path: Path,
                     ievalue_threshold: float,
                     lengths: dict[str, int] | None = None,
                     min_coverage: float = HMMER_COVERAGE,
                     order: dict[str, int] | None = None) -> pd.DataFrame | None:
    """Filtered Astra rows reduced to one per (protein, family): the lowest i-evalue,
    the earliest row on a tie.

    Columns: sequence_id, family, i_evalue, dom_bitscore, prof, row. ``row`` is the
    row's position in ``path``; ``prof`` is its profile's position in ``order`` (0
    without one), so (prof, row) orders rows from several profile-part tables as one
    single-process table. Returns None when no row passes.
    """
    try:
        # identifiers stay strings whatever a file's rows look like
        df = pd.read_csv(path, sep='\t', usecols=lambda c: c in ASTRA_COLUMNS,
                         dtype={'sequence_id': str, 'hmm_name': str})
    except Exception as e:
        print(f"Warning: Failed to parse Astra TSV: {e}")
        return None
    if df.empty:
        return None
    df['row'] = np.arange(len(df))

    # Filter by i-evalue
    df = df[df['i_evalue'] < ievalue_threshold]
    if df.empty:
        return None
    if lengths and {'hmm_from', 'hmm_to'} <= set(df.columns):
        leng = df['hmm_name'].map(lengths)
        if leng.isna().any():
            raise ValueError(f"{path.name}: {int(leng.isna().sum())} hits name profiles absent "
                             "from the HMM file's NAME/LENG records")
        before = len(df)
        df = df[(df['hmm_to'] - df['hmm_from']) / leng >= min_coverage]
        _log(f"  [{path.name}] HMM coverage >= {min_coverage}: {len(df):,} of {before:,} domain hits")
        if df.empty:
            return None
    else:
        _log(f"  [{path.name}] no HMM coordinates or profile lengths: i-evalue filter only")

    # Extract family from hmm_name
    families = {name: parse_hmm_family(name) for name in df['hmm_name'].unique()}
    df = df.assign(family=df['hmm_name'].map(families)).dropna(subset=['family'])
    if df.empty:
        return None
    df['prof'] = df['hmm_name'].map(order) if order is not None else 0

    # Per protein x family, the first row (in table order) with the lowest i-evalue
    df = (df.sort_values(['i_evalue', 'prof', 'row'], kind='stable')
            .drop_duplicates(['sequence_id', 'family'], keep='first'))
    return df[['sequence_id', 'family', 'i_evalue', 'dom_bitscore', 'prof', 'row']]


def _hits_from_rows(frames: list[pd.DataFrame | None]) -> dict[str, dict[str, tuple]]:
    """{protein_id: {family_class: (i_evalue, dom_bitscore)}} from ``_astra_best_rows`` frames.

    Several frames (profile parts of one search) reduce together, ties going to the
    earlier (prof, row), as within one table."""
    frames = [f for f in frames if f is not None]
    if not frames:
        return {}
    df = frames[0] if len(frames) == 1 else (
        pd.concat(frames, ignore_index=True)
          .sort_values(['i_evalue', 'prof', 'row'], kind='stable')
          .drop_duplicates(['sequence_id', 'family'], keep='first'))
    result: dict[str, dict[str, tuple]] = {}
    for pid, fam, ev, sc in zip(df['sequence_id'].tolist(), df['family'].tolist(),
                                df['i_evalue'].tolist(), df['dom_bitscore'].tolist()):
        result.setdefault(pid, {})[fam] = (ev, sc)
    return result


def _parse_astra_tsv(path: Path,
                     ievalue_threshold: float,
                     lengths: dict[str, int] | None = None,
                     min_coverage: float = HMMER_COVERAGE) -> dict[str, dict[str, tuple]]:
    """Parse Astra hits TSV and filter by i-evalue and, given profile lengths, HMM coverage.

    Astra columns: sequence_id, hmm_name, bitscore, evalue, c_evalue,
                   i_evalue, env_from, env_to, dom_bitscore[, ali_from, ali_to,
                   hmm_from, hmm_to]

    Coverage is (hmm_to - hmm_from) / profile length, as run_dbcan computes it.
    A table without HMM coordinates, or a run without profile lengths, is
    filtered on i-evalue only and says so.

    Returns {protein_id: {family_class: (best_i_evalue, best_dom_bitscore)}}, the
    first row with the lowest i-evalue for each protein and family.
    """
    return _hits_from_rows([_astra_best_rows(path, ievalue_threshold, lengths, min_coverage)])


# ── Consensus ──────────────────────────────────────────────────────────────────

def apply_consensus(diamond: dict, hmm: dict, sub: dict,
                    verbose: bool = True) -> pd.DataFrame:
    """Apply dbCAN consensus: keep proteins called by ≥2 of 3 tools.

    For each consensus protein × family, stores the best e-value and score
    from whichever tool provided them.

    Returns DataFrame with columns: protein_id, family_class, cazy_class,
    evalue, score, n_tools, tools.
    """
    all_pids = set(diamond) | set(hmm) | set(sub)

    rows = []
    for pid in all_pids:
        tools = []
        if pid in diamond:
            tools.append('diamond')
        if pid in hmm:
            tools.append('dbcan_hmm')
        if pid in sub:
            tools.append('dbcan_sub')

        if len(tools) < MIN_TOOLS:
            continue

        # Merge all families from agreeing tools, keep best e-value per family
        families: dict[str, tuple] = {}
        for src in (hmm, sub, diamond):  # priority order
            if pid not in src:
                continue
            for fam, (ev, sc) in src[pid].items():
                if fam not in families or ev < families[fam][0]:
                    families[fam] = (ev, sc)

        tools_str = ','.join(sorted(tools))
        for fam, (ev, sc) in families.items():
            rows.append({
                'protein_id': pid,
                'family_class': fam,
                'cazy_class': get_cazy_class(fam) or '',
                'evalue': ev,
                'score': sc,
                'n_tools': len(tools),
                'tools': tools_str,
            })

    # sorted rows: the table is identical however the hits were gathered
    df = (pd.DataFrame(rows).sort_values(['protein_id', 'family_class'], kind='stable',
                                         ignore_index=True)
          if rows else pd.DataFrame())

    if verbose:
        n_d, n_h, n_s = len(diamond), len(hmm), len(sub)
        n_cons = df['protein_id'].nunique() if not df.empty else 0
        print(f"\n  Per-tool protein counts:")
        print(f"    DIAMOND:    {n_d:>8,}")
        print(f"    dbCAN HMM:  {n_h:>8,}")
        print(f"    dbCAN-sub:  {n_s:>8,}")
        print(f"  Consensus (≥{MIN_TOOLS} tools): {n_cons:,} proteins")

    return df


# ── Main pipeline ──────────────────────────────────────────────────────────────

def _artifact_dir(db_path: str | None, results_out: str | Path | None) -> Path:
    """Where the consensus table and the per-tool side files go: beside the database,
    or beside ``results_out`` for a run on a protein FASTA alone."""
    if db_path:
        return Path(db_path).parent
    if results_out:
        return Path(results_out).resolve().parent
    raise ValueError("a run on --proteins without --db needs --out")


def _log_dbcan(dbcan: Path) -> tuple[Path, Path, bool]:
    dbcan_hmm, dbcan_sub, has_sub = dbcan_hmm_files(dbcan)
    _log(f"dbCAN databases: {dbcan}")
    _log(f"  CAZy.dmnd:     yes")
    _log(f"  dbCAN.hmm:     yes ({dbcan_hmm.stat().st_size / 1e6:.0f} MB)")
    _log(f"  dbCAN-sub.hmm: {'yes' if has_sub else 'NO — 2-tool mode'}"
         f"{f' ({dbcan_sub.stat().st_size / 1e6:.0f} MB)' if has_sub else ''}")
    return dbcan_hmm, dbcan_sub, has_sub


def _write_results(results_df: pd.DataFrame, out: Path, db_path: str | None, write_db: bool,
                   update_predicates: bool, n_proteins: int, verbose: bool) -> None:
    results_df.to_csv(out, sep='\t', index=False)
    _log(f"[cazyme] Consensus table: {out}")
    if write_db:
        load_cazyme_results(db_path, results_df, update_predicates=update_predicates,
                            n_proteins=n_proteins, verbose=verbose)


def classify_cazymes(
    db_path: str | None,
    threads: int = 4,
    update_predicates: bool = True,
    verbose: bool = True,
    dbcan_dir: str | Path | None = None,
    write_db: bool = True,
    results_out: str | Path | None = None,
    proteins: str | None = None,
) -> pd.DataFrame:
    """Full dbCAN 3-tool consensus CAZyme classification.

    Runs DIAMOND + Astra/HMMER (dbCAN.hmm) + Astra/HMMER (dbCAN-sub.hmm),
    applies ≥2-tool consensus filter, loads annotations as source='cazy',
    and updates predicates.

    Args:
        db_path: Path to sharur.duckdb (optional with ``proteins`` and ``write_db=False``)
        threads: CPU threads for DIAMOND and Astra
        update_predicates: Whether to update protein_predicates
        verbose: Print progress
        dbcan_dir: dbCAN release directory (default: $SHARUR_DBCAN_DIR, then
            the standard locations)
        write_db: False computes the consensus without opening the database
            for writing (load it later with ``load_cazyme_results``)
        results_out: Where to write the consensus table (default:
            cazyme_classification.tsv beside the database)
        proteins: Protein FASTA file or directory of .faa files to classify
            (default: discovered beside, or extracted from, the database)

    Returns:
        DataFrame with consensus classification results
    """
    t_total = _time.time()
    _log(f"[cazyme] Starting CAZyme classification pipeline")

    # Find databases
    try:
        dbcan = get_dbcan_path(dbcan_dir)
    except FileNotFoundError as e:
        _log(f"Error: {e}")
        return pd.DataFrame()

    dbcan_hmm, dbcan_sub, has_sub = _log_dbcan(dbcan) if verbose else dbcan_hmm_files(dbcan)
    data_dir = _artifact_dir(db_path, results_out)

    # Prepare proteins (explicit FASTA, existing .faa beside the database, or extract from it)
    _log(f"[cazyme] Finding protein sequences...")
    prot_dir, n_proteins, is_temp = resolve_proteins(db_path, proteins, verbose)
    if prot_dir is None:
        if verbose:
            _log("No protein sequences found")
        return pd.DataFrame()
    _log(f"[cazyme] {n_proteins:,} proteins ready ({_time.time() - t_total:.1f}s)")

    try:
        # Tool 1: DIAMOND (use cached results if available)
        t1 = _time.time()
        _log(f"[cazyme] Tool 1: DIAMOND search...")
        diamond_hits = _read_diamond_cache(data_dir / DIAMOND_CACHE, verbose=verbose)
        if diamond_hits is None:
            diamond_hits = run_diamond(prot_dir, dbcan, threads=threads, verbose=verbose)
            _save_diamond_cache(diamond_hits, data_dir, verbose=verbose)
        _log(f"[cazyme] Tool 1 done: {len(diamond_hits):,} hits ({_time.time() - t1:.1f}s)")

        # Tool 2: HMMER vs dbCAN.hmm from the same release as CAZy.dmnd and
        # dbCAN_sub.hmm (the HMM-search registry's "dbCAN" entry can be an
        # older release, so it is the fallback only)
        t2 = _time.time()
        _log(f"[cazyme] Tool 2: dbCAN HMM search...")
        hmm_cache = data_dir / "annotations" / "dbCAN_hits_df.tsv"
        hmm_hits = load_astra_tsv_cache(hmm_cache, verbose=verbose, lengths=hmm_lengths(dbcan_hmm))
        if hmm_hits is None:
            hmm_hits = run_astra_hmm(
                prot_dir, "dbCAN_HMM",
                hmm_path=dbcan_hmm,
                threads=threads, verbose=verbose,
                keep_tsv=data_dir / "cazyme_dbcan_hmm_hits.tsv",
            )
        _log(f"[cazyme] Tool 2 done: {len(hmm_hits):,} hits ({_time.time() - t2:.1f}s)")

        # Tool 3: Astra/HMMER vs dbCAN-sub.hmm — only on proteins already
        # hit by ≥1 other tool.  Consensus needs ≥2, so proteins seen by
        # neither DIAMOND nor dbCAN HMM can never pass even with a sub hit.
        sub_hits: dict = {}
        if has_sub:
            t3 = _time.time()
            candidate_pids = set(diamond_hits) | set(hmm_hits)
            _log(f"[cazyme] Tool 3: dbCAN-sub HMM on {len(candidate_pids):,} candidates...")
            if candidate_pids:
                sub_dir = _subset_fasta(prot_dir, candidate_pids, verbose=verbose)
                try:
                    sub_hits = run_astra_hmm(
                        sub_dir, "dbCAN_sub",
                        hmm_path=dbcan_sub,
                        threads=threads, verbose=verbose,
                        keep_tsv=data_dir / "cazyme_dbcan_sub_hits.tsv",
                    )
                finally:
                    shutil.rmtree(sub_dir, ignore_errors=True)
            _log(f"[cazyme] Tool 3 done: {len(sub_hits):,} hits ({_time.time() - t3:.1f}s)")

        # Consensus filter
        _log(f"[cazyme] Applying {MIN_TOOLS}-tool consensus filter...")
        results_df = apply_consensus(diamond_hits, hmm_hits, sub_hits, verbose=verbose)
    finally:
        if is_temp:
            shutil.rmtree(prot_dir, ignore_errors=True)

    if results_df.empty:
        if verbose:
            print("No CAZymes passed consensus filter")
        return results_df

    out = Path(results_out) if results_out else data_dir / "cazyme_classification.tsv"
    _write_results(results_df, out, db_path, write_db, update_predicates, n_proteins, verbose)
    _log(f"[cazyme] Total {_time.time() - t_total:.1f}s")
    return results_df


# ── Sharded runs ───────────────────────────────────────────────────────────────
#
# search shard i: dbCAN.hmm profile part i against every protein; shard 0 also runs
#   DIAMOND on every protein. DIAMOND's top target for a protein (default mode, -k 1)
#   depends on the other proteins in its run, so its query stays whole.
# sub shard i (after every search shard): dbCAN-sub profile part i against the
#   candidate subset built from all search shards.
# merge (after every sub shard): rebuild each tool's single-process table and apply
#   the consensus.
#
# HMMER reports i-evalue = P-value x Z, with Z the number of target sequences in the
# search. Splitting profiles keeps each search's target set, and so its Z, equal to
# the single-process run's: the whole proteome for dbCAN.hmm and the whole candidate
# subset for dbCAN-sub. Part i holds every N-th profile from profile i. Aksha writes
# each profile's rows together, in profile order, so the merge restores the
# single-process table row for row by a stable sort on profile position, and ties
# between rows resolve the same way too.

HMM_TSV_HEADER = ("sequence_id\thmm_name\tbitscore\tevalue\tc_evalue\ti_evalue\t"
                  "env_from\tenv_to\tdom_bitscore\tali_from\tali_to\thmm_from\thmm_to\n")


def _shard_dir(work_dir: Path, stage: str, index: int) -> Path:
    return Path(work_dir) / f"{stage}_{index:04d}"


def _proteome_fingerprint(prot_dir: str) -> dict:
    faas = sorted(Path(prot_dir).glob("*.faa"))
    return {"records": count_fasta_records(prot_dir), "bytes": sum(f.stat().st_size for f in faas)}


def _mark_done(shard_dir: Path, info: dict) -> None:
    tmp = shard_dir / "DONE.json.tmp"
    tmp.write_text(json.dumps(info, indent=1) + "\n")
    tmp.replace(shard_dir / "DONE.json")


def finished_shards(work_dir: Path, stage: str, n_shards: int) -> list[dict]:
    """Every shard's completion record for ``stage``; raises if one is missing or disagrees."""
    infos = []
    for i in range(n_shards):
        done = _shard_dir(work_dir, stage, i) / "DONE.json"
        if not done.exists():
            raise RuntimeError(f"{stage} shard {i} of {n_shards} has not finished ({done} missing)")
        infos.append(json.loads(done.read_text()))
    for info in infos:
        for key in ("n_shards", "proteins", "dbcan_dir"):
            if info[key] != infos[0][key]:
                raise RuntimeError(f"{stage} shards disagree on {key}: {infos[0][key]} vs {info[key]}")
    if infos[0]["n_shards"] != n_shards:
        raise RuntimeError(f"{stage} shards were run with {infos[0]['n_shards']} shards, not {n_shards}")
    return infos


def _hmm_part_search(target_dir: str, hmm_path: Path, part: int, n_parts: int, label: str,
                     threads: int, out_tsv: Path, verbose: bool) -> None:
    """Search profile part ``part`` of ``hmm_path`` against every sequence in ``target_dir``."""
    tmp = Path(tempfile.mkdtemp(prefix=f"cazyme_{label}_part_"))
    try:
        part_hmm = tmp / f"{hmm_path.stem}_part{part}.hmm"
        if write_hmm_part(hmm_path, part, n_parts, part_hmm) == 0:
            out_tsv.write_text(HMM_TSV_HEADER)
            return
        run_astra_hmm(target_dir, f"{label}_part{part}", hmm_path=part_hmm,
                      threads=threads, verbose=verbose, keep_tsv=out_tsv)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def _parse_shard_file(kind: str, path: str, lengths: dict | None, order: dict | None):
    """DIAMOND hits (a dict) for ``diamond``; ``_astra_best_rows`` for an HMM part table."""
    if kind == "diamond":
        return parse_diamond_tsv(Path(path))
    return _astra_best_rows(Path(path), HMMER_IEVALUE, lengths, order=order)


def _parse_files(jobs: list[tuple], threads: int) -> list:
    """Parse (kind, path, profile lengths, profile order) jobs, one file per worker process."""
    if threads > 1 and len(jobs) > 1:
        with ProcessPoolExecutor(max_workers=min(threads, len(jobs))) as pool:
            return list(pool.map(_parse_shard_file, *zip(*[(k, str(p), l, o) for k, p, l, o in jobs])))
    return [_parse_shard_file(kind, str(path), lengths, order) for kind, path, lengths, order in jobs]


def _fold_in_order(parts: list[dict]) -> dict[str, dict[str, tuple]]:
    """Fold the search shards' DIAMOND hits (shard 0 holds them all)."""
    out: dict[str, dict[str, tuple]] = {}
    for hits in parts:
        for pid, fams in hits.items():
            cur = out.setdefault(pid, {})
            for fam, val in fams.items():
                if fam not in cur or val[0] < cur[fam][0]:
                    cur[fam] = val
    return out


def combine_hmm_tables(parts: list[Path], hmm_path: Path, out: Path) -> Path:
    """Write the profile-part tables as one table in the single-process row order.

    Each part lists its profiles' rows together, in profile-file order, as Aksha
    writes them; a k-way merge on profile position restores the single-process
    table. A part out of profile order stops the merge.
    """
    order = hmm_profile_order(hmm_path)

    def keyed_rows(part: Path):
        last = -1
        with open(part) as fh:
            if fh.readline() != HMM_TSV_HEADER:
                raise ValueError(f"{part}: unexpected header")
            for line in fh:
                if not line.strip():
                    continue
                name = line.split("\t", 2)[1]
                if name not in order:
                    raise ValueError(f"{part}: profile {name} is absent from {hmm_path}")
                key = order[name]
                if key < last:
                    raise ValueError(f"{part}: rows are out of profile order at {name}")
                last = key
                yield key, line if line.endswith("\n") else line + "\n"

    with open(out, "w") as fout:
        fout.write(HMM_TSV_HEADER)
        for _, line in heapq.merge(*(keyed_rows(p) for p in parts), key=lambda kv: kv[0]):
            fout.write(line)
    return out


def shard_candidates(work_dir: Path, n_shards: int, dbcan: Path, threads: int = 1) -> tuple[set, int, int]:
    """Proteins with a DIAMOND or dbCAN.hmm hit in any search shard (the dbCAN-sub subset),
    with the DIAMOND and dbCAN.hmm protein counts."""
    dbcan_hmm, _, _ = dbcan_hmm_files(dbcan)
    hmm_len = hmm_lengths(dbcan_hmm)
    jobs = [("diamond", _shard_dir(work_dir, "search", i) / "diamond.tsv", None, None) for i in range(n_shards)]
    jobs += [("hmm", _shard_dir(work_dir, "search", i) / "dbcan_hmm.tsv", hmm_len, None) for i in range(n_shards)]
    parsed = _parse_files(jobs, threads)
    diamond = set().union(*parsed[:n_shards])
    hmm = set().union(*(set(f['sequence_id']) for f in parsed[n_shards:] if f is not None))
    return diamond | hmm, len(diamond), len(hmm)


def shard_search(prot_dir: str, dbcan: Path, index: int, n_shards: int, work_dir: Path,
                 threads: int = 4, verbose: bool = True) -> dict:
    """Search shard ``index``: dbCAN.hmm profile part on every protein; shard 0 adds DIAMOND."""
    t0 = _time.time()
    out = _shard_dir(work_dir, "search", index)
    out.mkdir(parents=True, exist_ok=True)
    dbcan_hmm, _, _ = dbcan_hmm_files(dbcan)
    info = {"stage": "search", "index": index, "n_shards": n_shards,
            "proteins": _proteome_fingerprint(prot_dir), "dbcan_dir": str(dbcan.resolve())}

    t1 = _time.time()
    if index == 0:
        _log(f"[cazyme] search shard 0/{n_shards}: DIAMOND on every protein...")
        run_diamond(prot_dir, dbcan, threads=threads, verbose=verbose, raw_out=out / "diamond.tsv")
    else:
        (out / "diamond.tsv").write_text("")
    t_diamond = _time.time() - t1

    _log(f"[cazyme] search shard {index}/{n_shards}: dbCAN.hmm profile part {index} on every protein...")
    t2 = _time.time()
    _hmm_part_search(prot_dir, dbcan_hmm, index, n_shards, "dbCAN_HMM", threads,
                     out / "dbcan_hmm.tsv", verbose)
    t_hmm = _time.time() - t2

    info["seconds"] = {"diamond": round(t_diamond, 1), "dbcan_hmm": round(t_hmm, 1),
                       "total": round(_time.time() - t0, 1)}
    _mark_done(out, info)
    _log(f"[cazyme] search shard {index}/{n_shards} done: {info['seconds']}")
    return info


def shard_sub(prot_dir: str, dbcan: Path, index: int, n_shards: int, work_dir: Path,
              threads: int = 4, verbose: bool = True) -> dict:
    """Sub shard ``index``: dbCAN-sub profile part on the candidates from every search shard."""
    t0 = _time.time()
    searched = finished_shards(work_dir, "search", n_shards)
    fingerprint = _proteome_fingerprint(prot_dir)
    if fingerprint != searched[0]["proteins"]:
        raise RuntimeError(f"protein input {fingerprint} differs from the search shards' "
                           f"{searched[0]['proteins']}")
    _, dbcan_sub, has_sub = dbcan_hmm_files(dbcan)
    out = _shard_dir(work_dir, "sub", index)
    out.mkdir(parents=True, exist_ok=True)

    candidates, n_diamond, n_hmm = shard_candidates(work_dir, n_shards, dbcan, threads)
    _log(f"[cazyme] sub shard {index}/{n_shards}: {len(candidates):,} candidates "
         f"({n_diamond:,} DIAMOND, {n_hmm:,} dbCAN.hmm)")
    if has_sub and candidates:
        sub_dir = _subset_fasta(prot_dir, candidates, verbose=verbose)
        try:
            _hmm_part_search(sub_dir, dbcan_sub, index, n_shards, "dbCAN_sub", threads,
                             out / "dbcan_sub.tsv", verbose)
        finally:
            shutil.rmtree(sub_dir, ignore_errors=True)

    info = {"stage": "sub", "index": index, "n_shards": n_shards, "proteins": fingerprint,
            "dbcan_dir": str(dbcan.resolve()), "candidates": len(candidates), "has_sub": has_sub,
            "seconds": {"total": round(_time.time() - t0, 1)}}
    _mark_done(out, info)
    _log(f"[cazyme] sub shard {index}/{n_shards} done: {info['seconds']}")
    return info


def merge_shards(work_dir: Path, n_shards: int, dbcan: Path, data_dir: Path,
                 threads: int = 1, verbose: bool = True) -> tuple[pd.DataFrame, int]:
    """Consensus over every shard's hits; writes the single-process run's side files
    (DIAMOND cache, raw HMM tables) to ``data_dir``. Returns (table, protein records)."""
    searched = finished_shards(work_dir, "search", n_shards)
    subbed = finished_shards(work_dir, "sub", n_shards)
    if subbed[0]["proteins"] != searched[0]["proteins"]:
        raise RuntimeError("sub shards searched a different protein input than the search shards")
    dbcan_hmm, dbcan_sub, _ = dbcan_hmm_files(dbcan)

    hmm_parts = [_shard_dir(work_dir, "search", i) / "dbcan_hmm.tsv" for i in range(n_shards)]
    sub_parts = [_shard_dir(work_dir, "sub", i) / "dbcan_sub.tsv" for i in range(n_shards)]
    has_sub_tables = all(p.exists() for p in sub_parts)  # dbCAN-sub ran (a release with it, and candidates)
    jobs = [("diamond", _shard_dir(work_dir, "search", i) / "diamond.tsv", None, None) for i in range(n_shards)]
    jobs += [("hmm", p, hmm_lengths(dbcan_hmm), hmm_profile_order(dbcan_hmm)) for p in hmm_parts]
    if has_sub_tables:
        sub_len, sub_order = hmm_lengths(dbcan_sub), hmm_profile_order(dbcan_sub)
        jobs += [("sub", p, sub_len, sub_order) for p in sub_parts]
    parsed = _parse_files(jobs, threads)
    diamond_hits = _fold_in_order(parsed[:n_shards])
    hmm_hits = _hits_from_rows(parsed[n_shards:2 * n_shards])
    sub_hits = _hits_from_rows(parsed[2 * n_shards:]) if has_sub_tables else {}
    _log(f"[cazyme] merged {n_shards} shards: {len(diamond_hits):,} DIAMOND, {len(hmm_hits):,} "
         f"dbCAN.hmm, {len(sub_hits):,} dbCAN-sub proteins ({subbed[0]['candidates']:,} candidates)")

    _save_diamond_cache(diamond_hits, data_dir, verbose=verbose)
    combine_hmm_tables(hmm_parts, dbcan_hmm, data_dir / "cazyme_dbcan_hmm_hits.tsv")
    if has_sub_tables:
        combine_hmm_tables(sub_parts, dbcan_sub, data_dir / "cazyme_dbcan_sub_hits.tsv")

    _log(f"[cazyme] Applying {MIN_TOOLS}-tool consensus filter...")
    return apply_consensus(diamond_hits, hmm_hits, sub_hits, verbose=verbose), searched[0]["proteins"]["records"]


def run_stage(stage: str, db_path: str | None, proteins: str | None, work_dir: str | Path,
              n_shards: int, shard_index: int | None = None, threads: int = 4,
              dbcan_dir: str | Path | None = None, results_out: str | Path | None = None,
              write_db: bool = False, update_predicates: bool = True,
              verbose: bool = True) -> pd.DataFrame | dict:
    """One step of a sharded run: ``search`` or ``sub`` for shard ``shard_index``, or ``merge``."""
    t0 = _time.time()
    dbcan = get_dbcan_path(dbcan_dir)
    work_dir = Path(work_dir)
    if stage == "merge":
        if verbose:
            _log_dbcan(dbcan)
        data_dir = _artifact_dir(db_path, results_out)
        results_df, n_proteins = merge_shards(work_dir, n_shards, dbcan, data_dir, threads, verbose)
        if results_df.empty:
            _log("No CAZymes passed consensus filter")
            return results_df
        out = Path(results_out) if results_out else data_dir / "cazyme_classification.tsv"
        _write_results(results_df, out, db_path, write_db, update_predicates, n_proteins, verbose)
        _log(f"[cazyme] merge done ({_time.time() - t0:.1f}s)")
        return results_df

    if shard_index is None or not 0 <= shard_index < n_shards:
        raise ValueError(f"--shard-index must lie in [0, {n_shards})")
    prot_dir, _, is_temp = resolve_proteins(db_path, proteins, verbose)
    if prot_dir is None:
        raise RuntimeError("no protein sequences found")
    try:
        step = shard_search if stage == "search" else shard_sub
        return step(prot_dir, dbcan, shard_index, n_shards, work_dir, threads, verbose)
    finally:
        if is_temp:
            shutil.rmtree(prot_dir, ignore_errors=True)


def load_cazyme_results(
    db_path: str,
    results_df: pd.DataFrame,
    update_predicates: bool = True,
    n_proteins: int | None = None,
    verbose: bool = True,
) -> dict[str, int]:
    """Replace source='cazy' annotations with a consensus table and update predicates.

    ``results_df`` is the table ``classify_cazymes`` writes (protein_id,
    family_class, cazy_class, evalue, score, n_tools, tools). Rows naming
    proteins absent from the database are dropped and counted.
    """
    _log(f"[cazyme] Opening DB for annotation loading...")
    db = duckdb.connect(db_path)
    known = set(r[0] for r in db.execute(
        "SELECT protein_id FROM proteins WHERE protein_id IN (SELECT UNNEST(?::VARCHAR[]))",
        [sorted(results_df['protein_id'].astype(str).unique())]).fetchall())
    unknown = int((~results_df['protein_id'].astype(str).isin(known)).sum())
    if unknown:
        _log(f"[cazyme] Dropping {unknown:,} rows for proteins absent from the database")
        results_df = results_df[results_df['protein_id'].astype(str).isin(known)].copy()
    if n_proteins is None:
        n_proteins = db.execute("SELECT COUNT(*) FROM proteins").fetchone()[0]

    # Clear existing CAZy annotations
    existing = db.execute(
        "SELECT COUNT(*) FROM annotations WHERE source = 'cazy'"
    ).fetchone()[0]
    if existing > 0:
        _log(f"[cazyme] Clearing {existing:,} existing CAZy annotations...")
        db.execute("DELETE FROM annotations WHERE source = 'cazy'")

    _log(f"[cazyme] Loading {len(results_df):,} CAZy annotations...")

    next_id = db.execute(
        "SELECT COALESCE(MAX(annotation_id), 0) FROM annotations"
    ).fetchone()[0]

    ann = results_df[['protein_id', 'family_class', 'evalue', 'score', 'tools']].copy()
    ann['annotation_id'] = range(next_id + 1, next_id + 1 + len(ann))
    ann['source'] = 'cazy'
    ann['accession'] = ann['family_class']
    ann['name'] = ann['family_class']
    ann['description'] = ann['family_class']
    ann['start_aa'] = None
    ann['end_aa'] = None

    insert_df = ann[['annotation_id', 'protein_id', 'source', 'accession',
                      'name', 'description', 'evalue', 'score',
                      'start_aa', 'end_aa']]
    db.register("tmp_cazy", insert_df)
    db.execute("""
        INSERT INTO annotations (
            annotation_id, protein_id, source, accession, name, description,
            evalue, score, start_aa, end_aa
        ) SELECT * FROM tmp_cazy
    """)
    db.unregister("tmp_cazy")

    if verbose:
        print(f"Loaded {len(insert_df):,} annotations")

    # ── Update predicates ──────────────────────────────────────────────────
    if update_predicates:
        t_pred = _time.time()
        _log(f"[cazyme] Updating predicates...")

        # Collect all new predicates per protein first (avoid row-by-row UPDATE)
        protein_families = (
            results_df.groupby('protein_id')['family_class'].apply(list).to_dict()
        )

        pid_new_preds: dict[str, set] = {}
        for pid, families in protein_families.items():
            new_preds = set()
            for fam in families:
                new_preds.update(get_predicates_for_cazy(fam))
                new_preds.add(f"cazy:{fam}")
            if new_preds:
                pid_new_preds[pid] = new_preds

        # Fetch current predicates for all affected proteins in one query
        if pid_new_preds:
            pids_list = list(pid_new_preds.keys())
            current_rows = db.execute(
                "SELECT protein_id, predicates FROM protein_predicates "
                "WHERE protein_id IN (SELECT UNNEST(?::VARCHAR[]))",
                [pids_list],
            ).fetchall()
            current_map = {r[0]: set(r[1]) for r in current_rows}

            updated = 0
            for pid, new_preds in pid_new_preds.items():
                existing = current_map.get(pid, set())
                # Remove old cazy: predicates, add new ones
                merged = {p for p in existing if not p.startswith("cazy:")} | new_preds
                if merged != existing:
                    db.execute(
                        "UPDATE protein_predicates SET predicates = ?, "
                        "updated_at = CURRENT_TIMESTAMP WHERE protein_id = ?",
                        [list(merged), pid],
                    )
                    updated += 1

            db.commit()
            _log(f"[cazyme] Updated predicates for {updated:,} proteins ({_time.time() - t_pred:.1f}s)")
        else:
            _log(f"[cazyme] No predicates to update")

    # ── Summary ────────────────────────────────────────────────────────────
    if verbose:
        n_prot = results_df['protein_id'].nunique()
        pct = 100 * n_prot / n_proteins if n_proteins > 0 else 0
        print(f"\n=== CAZyme CLASSIFICATION SUMMARY ===")
        print(f"Consensus proteins: {n_prot:,} / {n_proteins:,} ({pct:.1f}% of proteome)")
        if 'cazy_class' in results_df.columns:
            print(f"\nBy CAZy class:")
            cls_counts = results_df.groupby('cazy_class')['protein_id'].nunique()
            for cls, cnt in cls_counts.sort_values(ascending=False).items():
                if cls:
                    print(f"  {cls}: {cnt:,}")
        print(f"\nTop 20 families:")
        fam_counts = results_df.groupby('family_class')['protein_id'].nunique()
        for fam, cnt in fam_counts.sort_values(ascending=False).head(20).items():
            print(f"  {fam}: {cnt:,}")

    db.close()
    return {"annotations": len(insert_df), "proteins": int(results_df['protein_id'].nunique()),
            "dropped_unknown": unknown}


# ── CLI ────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Classify CAZymes using dbCAN 3-tool consensus",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--db", default=None,
        help="Path to sharur.duckdb (needed to load results; optional with --proteins)",
    )
    parser.add_argument(
        "--proteins", default=None,
        help="Protein FASTA file or directory of .faa files (default: beside or from --db)",
    )
    parser.add_argument(
        "--threads", type=int, default=4,
        help="CPU threads for DIAMOND and Astra (default: 4)",
    )
    parser.add_argument(
        "--no-update", action="store_true",
        help="Skip predicate updates",
    )
    parser.add_argument(
        "--dbcan-dir", default=None,
        help="dbCAN release directory (default: $SHARUR_DBCAN_DIR, then standard locations)",
    )
    parser.add_argument(
        "--results-only", action="store_true",
        help="Compute the consensus table without writing to the database",
    )
    parser.add_argument(
        "--out", default=None,
        help="Consensus table path (default: cazyme_classification.tsv beside the database)",
    )
    parser.add_argument(
        "--load-results", default=None, metavar="TSV",
        help="Load a precomputed consensus table into the database and exit",
    )
    parser.add_argument(
        "--quiet", action="store_true",
        help="Suppress progress output",
    )
    shards = parser.add_argument_group("sharded runs (one Slurm array task per shard)")
    shards.add_argument(
        "--stage", choices=("search", "sub", "merge"), default=None,
        help="search: dbCAN.hmm profile part i on every protein (shard 0 adds DIAMOND); "
             "sub: dbCAN-sub profile part i on the candidates (after every search shard); "
             "merge: consensus over all shards (after every sub shard)",
    )
    shards.add_argument("--n-shards", type=int, default=None, help="Number of shards N")
    shards.add_argument("--shard-index", type=int, default=None,
                        help="Shard i in [0, N) for --stage search/sub (e.g. $SLURM_ARRAY_TASK_ID)")
    shards.add_argument("--work-dir", default=None, help="Directory holding the shard outputs")

    args = parser.parse_args()
    if args.db is None and (args.proteins is None or args.load_results):
        parser.error("--db is required unless --proteins names the input")
    if args.db is None and not args.results_only and (args.stage in (None, "merge")):
        parser.error("loading results needs --db; add --results-only to only write the table")

    if args.stage:
        if args.n_shards is None or args.n_shards < 1 or args.work_dir is None:
            parser.error("--stage needs --n-shards N (N >= 1) and --work-dir")
        result = run_stage(
            args.stage, args.db, args.proteins, args.work_dir, args.n_shards,
            shard_index=args.shard_index, threads=args.threads, dbcan_dir=args.dbcan_dir,
            results_out=args.out, write_db=not args.results_only,
            update_predicates=not args.no_update, verbose=not args.quiet,
        )
        return 0 if args.stage != "merge" or not result.empty else 1

    if args.load_results:
        stats = load_cazyme_results(args.db, pd.read_csv(args.load_results, sep='\t'),
                                    update_predicates=not args.no_update, verbose=not args.quiet)
        _log(f"[cazyme] Loaded {stats}")
        return 0

    results = classify_cazymes(
        db_path=args.db,
        threads=args.threads,
        update_predicates=not args.no_update,
        verbose=not args.quiet,
        dbcan_dir=args.dbcan_dir,
        write_db=not args.results_only,
        results_out=args.out,
        proteins=args.proteins,
    )

    return 0 if not results.empty else 1


if __name__ == "__main__":
    exit(main())
