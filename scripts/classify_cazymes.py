#!/usr/bin/env python3
"""
CAZyme Classification Pipeline (dbCAN 3-tool consensus)

Recapitulates the standard dbCAN annotation pipeline:
  Tool 1: DIAMOND blastp vs CAZy.dmnd        (e-value ≤ 1e-102)
  Tool 2: Astra/HMMER vs dbCAN.hmm           (i-evalue < 1e-15 at run_dbcan Z, HMM coverage > 0.35)
  Tool 3: Astra/HMMER vs dbCAN-sub.hmm       (i-evalue < 1e-15 at run_dbcan Z, HMM coverage > 0.35)
  Consensus: retain genes called by ≥ 2 of 3 tools, with the families run_dbcan 5.2.9
             recommends for them (HMM domains after its overlap filter and grouping)

HMMER searches are delegated to Aksha (aksha search --hmm_in), which handles
parallelization and PyHMMER optimization. Only DIAMOND is run directly, in its
default sensitivity mode with -k 1, as run_dbcan runs it. The HMM i-evalue rescaling
matches run_dbcan 5.2.9 and later, which fix HMMER's Z at the profile count; 5.2.8
and earlier use pyhmmer's default Z, the number of proteins searched.

The DIAMOND cache and the raw dbCAN.hmm table are reused only when their provenance
stamp (protein content hash, reference-file content hash, search settings) matches.

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
import hashlib
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
# run_dbcan (bcb-unl/run_dbcan, dbcan/annotation/pyhmmer_search.py) searches every input
# protein with HMMER Z fixed at the HMM file's profile count (875 dbCAN.hmm, 25,618
# dbCAN-sub in the 2025-09-13 release) and keeps domains with i-evalue < 1e-15 and
# (hmm_to - hmm_from + 1) / LENG > 0.35; dbcan_ievalue_scale converts Aksha's Z.
DIAMOND_EVALUE = 1e-102
HMMER_IEVALUE  = 1e-15
HMMER_COVERAGE = 0.35
MIN_TOOLS      = 2

CAZY_PREFIX_RE = re.compile(r'^(GH|GT|PL|CE|CBM|AA)\d+')

# run_dbcan 5.2.9 (dbcan/constants/OverviewGenerator_constants.py:48, CAZY_ID_PATTERN): the
# prefix that marks a CAZy token in a CAZy.dmnd subject ID. It also matches accessions such
# as AAO30670.1, so a token becomes a family only when it is a whole family or subfamily name.
RUN_DBCAN_CAZY_ID_RE = re.compile(r'^(GH|GT|CBM|AA|CE|PL)')
CAZY_FAMILY_TOKEN_RE = re.compile(r'^(GH|GT|PL|CE|CBM|AA)\d+(_\d+)?$')

# run_dbcan 5.2.9 overlap thresholds: per-tool domain filter
# (dbcan/constants/process_utils_constants.py:3) and hmm/sub grouping
# (dbcan/constants/OverviewGenerator_constants.py:44).
RUN_DBCAN_OVERLAP_RATIO = 0.5
RUN_DBCAN_GROUP_OVERLAP = 0.5
RUN_DBCAN_GT2_PREFIX, RUN_DBCAN_GT2_NAME = "GT2_", "GT2.hmm"  # pyhmmer_search.py renames GT2 profiles

# Settings a DIAMOND cache must have been produced with to be reused (see diamond_stamp).
DIAMOND_SETTINGS = {
    "program": "diamond blastp",
    "sensitivity": "default",
    "max_target_seqs": 1,
    "evalue": DIAMOND_EVALUE,
    "subject_parser": "run_dbcan 5.2.9 OverviewGenerator.extract_cazy_id",
}


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
    """(dbCAN.hmm, dbCAN-sub HMM, whether the sub-family file is usable).

    dbCAN releases name the sub-family profiles dbCAN_sub.hmm; older docs used dbCAN-sub.hmm.
    A directory holding both names is ambiguous (a release downloaded over an older one)
    and stops the run.
    """
    present = [dbcan / n for n in ("dbCAN-sub.hmm", "dbCAN_sub.hmm") if (dbcan / n).exists()]
    if len(present) > 1:
        raise ValueError(f"{dbcan} holds both dbCAN-sub.hmm and dbCAN_sub.hmm; keep the file of the "
                         "release that CAZy.dmnd and dbCAN.hmm come from and remove the other")
    sub = present[0] if present else dbcan / "dbCAN-sub.hmm"
    return dbcan / "dbCAN.hmm", sub, sub.exists() and sub.stat().st_size > 0


# ── Protein extraction ─────────────────────────────────────────────────────────

def count_fasta_records(prot_dir: str) -> int:
    """Header lines across the directory's .faa files."""
    n = 0
    for faa in sorted(Path(prot_dir).glob("*.faa")):
        with open(faa, "rb") as fh:
            n += sum(1 for line in fh if line.startswith(b">"))
    return n


# ── Provenance stamps for reusable caches ──────────────────────────────────────
#
# A cache (the DIAMOND hits, the dbCAN.hmm table) is reused only when its stamp equals
# the stamp of the run about to use it: the protein input's content hash, the reference
# file's size and content hash, and the search settings. A cache without a stamp, or
# with any other stamp, is ignored and the tool reruns.

_BLOCK = 1 << 24
_DIGESTS: dict[tuple, dict] = {}


def file_digest(path: Path) -> dict:
    """{'name', 'bytes', 'blake2b'} of a reference file: its full-content BLAKE2b-128 hash
    (about 3.5 s for the 2.2 GB CAZy.dmnd on local SSD), memoized per size and mtime."""
    path = Path(path)
    st = path.stat()
    key = (str(path.resolve()), st.st_size, st.st_mtime_ns)
    if key not in _DIGESTS:
        h = hashlib.blake2b(digest_size=16)
        with open(path, "rb") as fh:
            while chunk := fh.read(_BLOCK):
                h.update(chunk)
        _DIGESTS[key] = {"name": path.name, "bytes": st.st_size, "blake2b": h.hexdigest()}
    return _DIGESTS[key]


def protein_fingerprint(prot_dir: str) -> dict:
    """Record count, byte count and BLAKE2b-128 content hash of the directory's .faa files,
    read in file-name order."""
    h = hashlib.blake2b(digest_size=16)
    records = size = 0
    for faa in sorted(Path(prot_dir).glob("*.faa")):
        with open(faa, "rb") as fh:
            for line in fh:
                h.update(line)
                size += len(line)
                records += line.startswith(b">")
    return {"records": records, "bytes": size, "blake2b": h.hexdigest()}


def _stamp_path(table: Path) -> Path:
    return Path(table).with_name(Path(table).name + ".stamp.json")


def _normalized(stamp: dict) -> dict:
    return json.loads(json.dumps(stamp, sort_keys=True))


def write_stamp(table: Path, stamp: dict) -> None:
    tmp = _stamp_path(table).with_suffix(".tmp")
    tmp.write_text(json.dumps(stamp, indent=1, sort_keys=True) + "\n")
    tmp.replace(_stamp_path(table))


def clear_stamp(table: Path) -> None:
    """Drop a table's stamp before rewriting it, so an interrupted write leaves no stamp."""
    _stamp_path(table).unlink(missing_ok=True)


def stamp_matches(table: Path, stamp: dict, label: str) -> bool:
    """True when ``table`` exists with a stamp equal to ``stamp``; logs why a table is ignored."""
    table = Path(table)
    if not table.exists():
        return False
    sp = _stamp_path(table)
    if not sp.exists():
        _log(f"[cazyme] Ignoring {label} {table}: it carries no provenance stamp")
        return False
    try:
        found = json.loads(sp.read_text())
    except (OSError, ValueError) as e:
        _log(f"[cazyme] Ignoring {label} {table}: unreadable stamp ({e})")
        return False
    expected = _normalized(stamp)
    if found != expected:
        differs = sorted(k for k in set(found) | set(expected) if found.get(k) != expected.get(k))
        _log(f"[cazyme] Ignoring {label} {table}: stamp differs in {', '.join(differs)}")
        return False
    return True


def diamond_stamp(proteins: dict, dbcan: Path) -> dict:
    """What a DIAMOND cache must match: the protein input, CAZy.dmnd and DIAMOND_SETTINGS."""
    return {"kind": "diamond", "proteins": proteins, "reference": file_digest(dbcan / "CAZy.dmnd"),
            "settings": DIAMOND_SETTINGS}


def hmm_table_stamp(proteins: dict, hmm_path: Path) -> dict:
    """What a raw dbCAN.hmm table must match: the protein input and the profile file."""
    return {"kind": "dbcan_hmm", "proteins": proteins, "reference": file_digest(hmm_path),
            "settings": {"search": "aksha search --hmm_in, every protein (Z = target count)"}}


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
    """Base CAZy families in a CAZy.dmnd subject ID, chosen as run_dbcan chooses them.

    Mirrors run_dbcan 5.2.9 OverviewGenerator.extract_cazy_id (dbcan/IO/OverviewGenerator.py:157-198):
    when a token names a FASTA file (JGI subjects such as ``GH3|355645|<genome>.aa.fasta``),
    the ID is cut at its first all-digit token; the CAZy tokens are then those from the
    first one matching the CAZy prefix onward. The family can be the first token: 359,744
    of the 4,098,879 subjects in the db_v5-2_9-13-2025 CAZy.dmnd are JGI subjects. A token
    counts as a family when it is a whole family or subfamily name (``GH13``, ``GH13_32``),
    which leaves out accessions that merely start with a prefix (``AAO30670.1``).

    Returns base families in subject order, each once (``["CBM48", "GH13"]``).
    """
    parts = [p.strip() for p in subject_id.split("|")]
    if any("fasta" in p.lower() for p in parts):
        for idx, tok in enumerate(parts):
            if tok.isdigit():
                parts = parts[:idx]
                break
    first = next((i for i, t in enumerate(parts) if RUN_DBCAN_CAZY_ID_RE.match(t)), None)
    if first is None:
        return []
    families: list[str] = []
    for tok in parts[first:]:
        if RUN_DBCAN_CAZY_ID_RE.match(tok) and CAZY_FAMILY_TOKEN_RE.match(tok):
            fam = tok.split("_")[0]
            if fam not in families:
                families.append(fam)
    return families


DIAMOND_CACHE = "cazyme_diamond_hits.tsv"
HMM_TABLE = "cazyme_dbcan_hmm_hits.tsv"
SUB_TABLE = "cazyme_dbcan_sub_hits.tsv"


def side_files(db_path: str | None, results_out: str | Path | None) -> dict[str, Path]:
    """Per-tool side files of one run: {'diamond', 'hmm', 'sub'} -> path.

    With ``results_out`` they sit beside it under its stem (``cazy_arc.tsv`` ->
    ``cazy_arc.cazyme_diamond_hits.tsv``), so runs with different outputs never share a
    file; without it, beside the database under the default names.
    """
    if results_out:
        out = Path(results_out).resolve()
        base, prefix = out.parent, f"{out.stem}."
    elif db_path:
        base, prefix = Path(db_path).resolve().parent, ""
    else:
        raise ValueError("a run on --proteins without --db needs --out")
    return {"diamond": base / f"{prefix}{DIAMOND_CACHE}", "hmm": base / f"{prefix}{HMM_TABLE}",
            "sub": base / f"{prefix}{SUB_TABLE}"}


def load_diamond_cache(db_path: str, stamp: dict,
                       verbose: bool = True) -> dict[str, dict[str, tuple]] | None:
    """Load cazyme_diamond_hits.tsv beside the database when its stamp equals ``stamp``.

    The cache has its own file: the consensus table (cazyme_classification.tsv)
    merges HMM families and must never be read back as DIAMOND evidence.

    Returns {protein_id: {family_class: (evalue, bitscore)}} or None if no usable cache.
    """
    return _read_diamond_cache(side_files(db_path, None)["diamond"], stamp, verbose)


def _read_diamond_cache(tsv: Path, stamp: dict,
                        verbose: bool = True) -> dict[str, dict[str, tuple]] | None:
    if not stamp_matches(tsv, stamp, "DIAMOND cache"):
        return None

    try:
        df = pd.read_csv(tsv, sep='\t', dtype={'protein_id': str, 'family_class': str})
    except Exception as e:
        _log(f"[cazyme] Ignoring DIAMOND cache {tsv}: unreadable ({e})")
        return None

    if 'protein_id' not in df.columns:
        _log(f"[cazyme] Ignoring DIAMOND cache {tsv}: no protein_id column")
        return None

    if verbose:
        print(f"Loading cached DIAMOND results from {tsv.name} "
              f"({df['protein_id'].nunique():,} proteins)")

    hits: dict[str, dict[str, tuple]] = defaultdict(dict)
    for pid, fam, ev, bs in zip(df['protein_id'], df['family_class'], df['evalue'], df['bitscore']):
        ev, bs = float(ev), float(bs)
        if ev > DIAMOND_EVALUE:
            continue
        if isinstance(fam, str) and CAZY_PREFIX_RE.match(fam):
            if fam not in hits[pid] or ev < hits[pid][fam][0]:
                hits[pid][fam] = (ev, bs)

    if verbose:
        print(f"  DIAMOND (cached): {len(hits):,} proteins")
    return dict(hits)


def _save_diamond_cache(hits: dict, out: Path, stamp: dict,
                        verbose: bool = True) -> None:
    """Save DIAMOND hits to ``out`` with their provenance stamp, for reuse by a matching run.

    The stamp is removed before the table is rewritten and written after it, so an
    interrupted write leaves a table that no run reuses.
    """
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
    clear_stamp(out)
    df = (pd.DataFrame(rows).sort_values(['protein_id', 'family_class'], kind='stable')
          if rows else pd.DataFrame(columns=['protein_id', 'family_class', 'cazy_class', 'evalue', 'bitscore']))
    df.to_csv(out, sep='\t', index=False)
    write_stamp(out, stamp)
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
    # Every failure raises: an empty DIAMOND set would silently turn the 3-tool consensus
    # into HMM-only, and a sharded search shard would record success without diamond.tsv.
    if shutil.which("diamond") is None:
        raise RuntimeError("diamond is not on PATH; Tool 1 of the consensus cannot run")
    faa_files = sorted(Path(prot_dir).glob("*.faa"))
    if not faa_files:
        raise FileNotFoundError(f"no .faa files in {prot_dir} for DIAMOND")

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
    except FileNotFoundError as e:
        shutil.rmtree(work_dir, ignore_errors=True)
        err_handle.close()
        raise RuntimeError("diamond could not be started; Tool 1 of the consensus cannot run") from e

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


def dbcan_ievalue_scale(hmm_path: Path, n_targets: int, n_profiles: int | None = None) -> float:
    """Factor that moves Aksha i-evalues onto run_dbcan's scale.

    i-evalue = P-value x Z. Aksha (pyhmmer's default) sets Z to the number of target
    sequences searched; run_dbcan fixes Z at the number of profiles in the HMM file
    (hmmscan's convention, kept by its pyhmmer port), so its E-values do not depend
    on how many proteins are searched. Multiplying by n_profiles / n_targets gives
    run_dbcan's i-evalue for every domain row.
    """
    if n_targets <= 0:
        raise ValueError(f"{hmm_path.name}: no target sequences to scale E-values by")
    if n_profiles is None:
        n_profiles = len(hmm_lengths(hmm_path))
    if n_profiles <= 0:
        raise ValueError(f"{hmm_path.name} has no NAME/LENG profile records")
    return n_profiles / n_targets


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
                  keep_tsv: Path | None = None,
                  ievalue_scale: float = 1.0,
                  order: dict[str, int] | None = None,
                  domains: bool = False):
    """Run Astra hmmsearch against an HMM database.

    Uses --installed_hmms for Astra-installed databases (e.g., dbCAN),
    or --hmm_in for custom HMM files (e.g., dbCAN-sub.hmm).

    ``ievalue_scale`` multiplies each i-evalue before the cut (``dbcan_ievalue_scale``
    puts it on run_dbcan's scale); the kept raw table holds Aksha's own values.

    Returns {protein_id: {family_class: (i_evalue, dom_bitscore)}}; with ``domains``,
    (that dict, the passing domain rows from ``_astra_passing_rows``).
    """
    if verbose:
        name = installed_hmm or (hmm_path.name if hmm_path else "?")
        print(f"Running Astra search against {name}...")

    if installed_hmm:
        database = ["--installed_hmms", installed_hmm]
    elif hmm_path:
        database = ["--hmm_in", str(hmm_path)]
    else:
        raise ValueError(f"HMM search ({label}) names no HMM database")
    out_dir = tempfile.mkdtemp(prefix=f"cazyme_{label}_")
    cmd = [_hmm_search_command(), "search", "--prot_in", prot_dir,
           "--outdir", out_dir, "--threads", str(threads), *database]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        # fatal: a missing tool would silently shrink the consensus
        if result.returncode != 0:
            raise RuntimeError(f"HMM search ({label}) failed (rc={result.returncode}): {result.stderr[-500:]}")
        # Find output TSV (name depends on --installed_hmms vs --hmm_in)
        tsvs = list(Path(out_dir).glob("*_hits_df.tsv"))
        if not tsvs:
            raise RuntimeError(f"HMM search ({label}) wrote no *_hits_df.tsv")
        lengths = hmm_lengths(hmm_path) if hmm_path else None
        rows = _astra_passing_rows(tsvs[0], ievalue, lengths, order=order, ievalue_scale=ievalue_scale)
        if keep_tsv is not None:  # keep the raw per-domain table for audit and reuse
            clear_stamp(keep_tsv)
            shutil.copy2(tsvs[0], keep_tsv)
    finally:
        shutil.rmtree(out_dir, ignore_errors=True)

    hits = _hits_from_rows([_best_per_family(rows)])
    if verbose:
        print(f"  {label}: {len(hits):,} proteins with hits")
    return (hits, rows) if domains else hits


ASTRA_COLUMNS = {'sequence_id', 'hmm_name', 'i_evalue', 'dom_bitscore', 'hmm_from', 'hmm_to',
                 'ali_from', 'ali_to'}
DOMAIN_COLUMNS = ['sequence_id', 'hmm_name', 'family', 'i_evalue', 'dom_bitscore',
                  'ali_from', 'ali_to', 'prof', 'row']


def _astra_passing_rows(path: Path,
                        ievalue_threshold: float,
                        lengths: dict[str, int] | None = None,
                        min_coverage: float = HMMER_COVERAGE,
                        order: dict[str, int] | None = None,
                        ievalue_scale: float = 1.0) -> pd.DataFrame | None:
    """Every Astra domain row that passes the dbCAN cuts, in table order.

    i-evalues are multiplied by ``ievalue_scale`` before the cut and are reported
    on that scale. Coverage is (hmm_to - hmm_from + 1) / LENG > ``min_coverage``,
    as run_dbcan computes it. Rows of profiles outside the CAZy classes (cohesin,
    dockerin, SLH) stay, with no family, because run_dbcan's overlap filter and
    family recommendation see them.

    Columns: DOMAIN_COLUMNS. ``row`` is the row's position in ``path``; ``prof`` is its
    profile's position in ``order`` (0 without one), so (prof, row) orders rows from
    several profile-part tables as one single-process table. ``ali_from``/``ali_to``
    (the alignment's target coordinates, run_dbcan's Target From/To) are NaN in a table
    without them. Returns None when no row passes.
    """
    try:
        # identifiers stay strings whatever a file's rows look like
        df = pd.read_csv(path, sep='\t', usecols=lambda c: c in ASTRA_COLUMNS,
                         dtype={'sequence_id': str, 'hmm_name': str})
    except pd.errors.EmptyDataError:
        return None
    if df.empty:
        return None
    df['row'] = np.arange(len(df))

    # Filter by i-evalue, on run_dbcan's Z when a scale is given
    if ievalue_scale != 1.0:
        df['i_evalue'] = df['i_evalue'] * ievalue_scale
    df = df[df['i_evalue'] < ievalue_threshold]
    if df.empty:
        return None
    if lengths and {'hmm_from', 'hmm_to'} <= set(df.columns):
        leng = df['hmm_name'].map(lengths)
        if leng.isna().any():
            raise ValueError(f"{path.name}: {int(leng.isna().sum())} hits name profiles absent "
                             "from the HMM file's NAME/LENG records")
        before = len(df)
        df = df[(df['hmm_to'] - df['hmm_from'] + 1) / leng > min_coverage]
        _log(f"  [{path.name}] HMM coverage > {min_coverage}: {len(df):,} of {before:,} domain hits")
        if df.empty:
            return None
    else:
        _log(f"  [{path.name}] no HMM coordinates or profile lengths: i-evalue filter only")

    families = {name: parse_hmm_family(name) for name in df['hmm_name'].unique()}
    df = df.assign(family=df['hmm_name'].map(families))
    df['prof'] = df['hmm_name'].map(order) if order is not None else 0
    for col in ('ali_from', 'ali_to'):
        if col not in df.columns:
            df[col] = np.nan
    return df[DOMAIN_COLUMNS].reset_index(drop=True)


def _best_per_family(df: pd.DataFrame | None) -> pd.DataFrame | None:
    """One row per (protein, CAZy family): the lowest i-evalue, the earliest (prof, row) on a tie."""
    if df is None:
        return None
    df = df.dropna(subset=['family'])
    if df.empty:
        return None
    df = (df.sort_values(['i_evalue', 'prof', 'row'], kind='stable')
            .drop_duplicates(['sequence_id', 'family'], keep='first'))
    return df[['sequence_id', 'family', 'i_evalue', 'dom_bitscore', 'prof', 'row']]


def _astra_best_rows(path: Path,
                     ievalue_threshold: float,
                     lengths: dict[str, int] | None = None,
                     min_coverage: float = HMMER_COVERAGE,
                     order: dict[str, int] | None = None,
                     ievalue_scale: float = 1.0) -> pd.DataFrame | None:
    """Filtered Astra rows reduced to one per (protein, CAZy family) (``_best_per_family``)."""
    return _best_per_family(_astra_passing_rows(path, ievalue_threshold, lengths, min_coverage,
                                                order, ievalue_scale))


def _hits_from_rows(frames: list[pd.DataFrame | None]) -> dict[str, dict[str, tuple]]:
    """{protein_id: {family_class: (i_evalue, dom_bitscore)}} from ``_best_per_family`` frames.

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


def _domain_table(frames: list[pd.DataFrame | None]) -> pd.DataFrame | None:
    """The passing domain rows of a search's profile parts as one table."""
    frames = [f for f in frames if f is not None]
    if not frames:
        return None
    return frames[0] if len(frames) == 1 else pd.concat(frames, ignore_index=True)


def _parse_astra_tsv(path: Path,
                     ievalue_threshold: float,
                     lengths: dict[str, int] | None = None,
                     min_coverage: float = HMMER_COVERAGE,
                     ievalue_scale: float = 1.0) -> dict[str, dict[str, tuple]]:
    """Parse Astra hits TSV and filter by i-evalue and, given profile lengths, HMM coverage.

    Astra columns: sequence_id, hmm_name, bitscore, evalue, c_evalue,
                   i_evalue, env_from, env_to, dom_bitscore[, ali_from, ali_to,
                   hmm_from, hmm_to]

    Coverage is (hmm_to - hmm_from + 1) / profile length, as run_dbcan computes it.
    A table without HMM coordinates, or a run without profile lengths, is
    filtered on i-evalue only and says so.

    Returns {protein_id: {family_class: (best_i_evalue, best_dom_bitscore)}}, the
    first row with the lowest i-evalue for each protein and family.
    """
    return _hits_from_rows([_astra_best_rows(path, ievalue_threshold, lengths, min_coverage,
                                             ievalue_scale=ievalue_scale)])


# ── run_dbcan family recommendation ────────────────────────────────────────────
#
# run_dbcan 5.2.9 reports, for a protein that ≥2 tools hit, the families of its HMM
# domains: each tool's domain table is first thinned of strongly overlapping domains
# (dbcan/process/process_utils.py:102-172, filter_overlaps), then, when both HMM tools
# hit, dbCAN.hmm and dbCAN-sub domains are grouped by overlap and one domain per group
# is kept (dbcan/IO/OverviewGenerator.py:200-221 and 337-409); with one HMM tool, every
# domain of that tool is kept; DIAMOND's families are used only when neither HMM tool
# hit (OverviewGenerator.py:290-333), which a ≥2-tool protein reaches only in a
# DIAMOND-only configuration. The functions below follow that code line for line.
#
# Domain entries are tuples (name, target_from, target_to, tool, i_evalue, score), tool
# being 'hmm' or 'sub' as in run_dbcan. One documented difference: run_dbcan orders a
# protein's domains by (from, to) with pandas' default quicksort, which leaves the order
# of domains with identical coordinates unspecified; here those ties go to profile-file
# order, so single-process and sharded runs agree.

def run_dbcan_domain_name(hmm_name: str, tool: str) -> str:
    """The name run_dbcan 5.2.9 shows for a dbCAN.hmm ('hmm') or dbCAN-sub ('sub') domain.

    GT2 profiles are renamed GT2.hmm at search time (annotation/pyhmmer_search.py:148-149);
    dbCAN-sub names keep their '.hmm' segments without the suffix
    (process/process_dbcan_sub.py:205-208); dbCAN.hmm names lose '.hmm'
    (IO/OverviewGenerator.py:137-141).
    """
    name = RUN_DBCAN_GT2_NAME if RUN_DBCAN_GT2_PREFIX in hmm_name else hmm_name
    if tool == "sub":
        parts = [p.split(".hmm")[0] for p in name.split("|") if p.endswith(".hmm")]
        return ";".join(parts) if parts else "-"
    return name.split(".hmm")[0] if ".hmm" in name else name


def run_dbcan_filter_overlaps(entries: list[tuple]) -> list[tuple]:
    """process_utils.py:130-161: one protein's domains, sorted by (from, to), thinned of
    domains that overlap the last kept one by more than half of either length; of two
    such domains the lower i-evalue stays, the earlier on a tie."""
    keep: list[tuple] = []
    for cur in entries:
        if not keep:
            keep.append(cur)
            continue
        last = keep[-1]
        overlap = min(last[2], cur[2]) - max(last[1], cur[1])
        if overlap > 0:
            r_last = overlap / max(1, last[2] - last[1])
            r_cur = overlap / max(1, cur[2] - cur[1])
            if r_last > RUN_DBCAN_OVERLAP_RATIO or r_cur > RUN_DBCAN_OVERLAP_RATIO:
                if last[4] > cur[4]:
                    keep[-1] = cur
            else:
                keep.append(cur)
        else:
            keep.append(cur)
    return keep


def _run_dbcan_overlap(s1, e1, s2, e2) -> bool:
    """OverviewGenerator.calculate_overlap (lines 200-206)."""
    overlap = max(0, min(e1, e2) - max(s1, s2) + 1)
    return overlap / min(e1 - s1 + 1, e2 - s2 + 1) > RUN_DBCAN_GROUP_OVERLAP


def _run_dbcan_select_best(group: list[tuple]) -> tuple:
    """OverviewGenerator.select_best_result (lines 208-221)."""
    for entry in group:
        if "_" in entry[0] and entry[3] == "hmm":
            return entry
    subs = [e for e in group if e[3] == "sub"]
    if subs:
        return subs[0]
    hmms = [e for e in group if e[3] == "hmm"]
    if hmms:
        return hmms[0]
    return group[0]


def run_dbcan_graph_grouping(all_results: list[tuple]) -> list[tuple]:
    """OverviewGenerator.graph_based_grouping (lines 337-409)."""
    if not all_results:
        return []
    hmm_results = sorted([r for r in all_results if r[3] == "hmm"], key=lambda x: x[1])
    sub_results = sorted([r for r in all_results if r[3] == "sub"], key=lambda x: x[1])

    if len(sub_results) == 1 and len(hmm_results) > 1:
        sub = sub_results[0]
        if all(_run_dbcan_overlap(sub[1], sub[2], h[1], h[2]) for h in hmm_results):
            return [_run_dbcan_select_best(sub_results + hmm_results)]

    if len(hmm_results) == 1 and len(sub_results) > 1:
        hmm = hmm_results[0]
        if "_" in hmm[0] and len({s[0] for s in sub_results}) == 1:
            if all(_run_dbcan_overlap(hmm[1], hmm[2], s[1], s[2]) for s in sub_results):
                return [_run_dbcan_select_best([hmm] + sub_results)]

    processed_subs: set[int] = set()
    sub_groups: list[list[tuple]] = []
    for i, sub1 in enumerate(sub_results):
        if i in processed_subs:
            continue
        group = [sub1]
        processed_subs.add(i)
        changed = True
        while changed:
            changed = False
            for j, sub2 in enumerate(sub_results):
                if j in processed_subs:
                    continue
                if any(_run_dbcan_overlap(g[1], g[2], sub2[1], sub2[2]) for g in group):
                    group.append(sub2)
                    processed_subs.add(j)
                    changed = True
        sub_groups.append(group)

    processed_hmms: set[int] = set()
    all_groups: list[list[tuple]] = []
    for sub_group in sub_groups:
        group = sub_group.copy()
        for i, hmm in enumerate(hmm_results):
            if any(_run_dbcan_overlap(s[1], s[2], hmm[1], hmm[2]) for s in sub_group):
                group.append(hmm)
                processed_hmms.add(i)
        all_groups.append(group)

    for i, hmm in enumerate(hmm_results):
        if i in processed_hmms:
            continue
        group = [hmm]
        processed_hmms.add(i)
        for j, hmm2 in enumerate(hmm_results):
            if j in processed_hmms or j == i:
                continue
            if _run_dbcan_overlap(hmm[1], hmm[2], hmm2[1], hmm2[2]):
                group.append(hmm2)
                processed_hmms.add(j)
        all_groups.append(group)

    return [_run_dbcan_select_best(g) for g in all_groups if g]


def run_dbcan_recommendation(hmm_entries: list[tuple], sub_entries: list[tuple],
                             diamond_families: dict[str, tuple] | None = None) -> list[tuple]:
    """The domains run_dbcan recommends for a ≥2-tool protein (OverviewGenerator.py:290-333).

    ``hmm_entries`` and ``sub_entries`` are the protein's overlap-filtered domains in table
    order. DIAMOND entries, used only when neither HMM tool hit, are
    (family, None, None, 'diamond', evalue, bitscore).
    """
    if hmm_entries and sub_entries:
        grouped = sorted(run_dbcan_graph_grouping(list(hmm_entries) + list(sub_entries)),
                         key=lambda x: x[1])
        seen: dict[tuple, tuple] = {}
        for entry in grouped:  # dict.fromkeys over "name(from-to)" strings
            seen.setdefault((entry[0], entry[1], entry[2]), entry)
        return list(seen.values())
    if hmm_entries:
        return list(hmm_entries)
    if sub_entries:
        return list(sub_entries)
    return [(fam, None, None, "diamond", ev, bs) for fam, (ev, bs) in (diamond_families or {}).items()]


def _filtered_entries(domains: pd.DataFrame | None, pids: set, tool: str) -> dict[str, list[tuple]]:
    """{protein: overlap-filtered domain entries in table order} for ``pids``."""
    if domains is None or not pids:
        return {}
    df = domains[domains['sequence_id'].isin(pids)]
    if df.empty:
        return {}
    if df[['ali_from', 'ali_to']].isna().any().any():
        raise ValueError(f"{tool} domain rows lack ali_from/ali_to, which run_dbcan's family "
                         "recommendation needs; rerun the HMM search with Aksha")
    df = df.sort_values(['sequence_id', 'ali_from', 'ali_to', 'prof', 'row'], kind='stable')
    names = {n: run_dbcan_domain_name(n, tool) for n in df['hmm_name'].unique()}
    out: dict[str, list[tuple]] = {}
    for pid, name, frm, to, ev, sc in zip(df['sequence_id'].tolist(), df['hmm_name'].map(names).tolist(),
                                          df['ali_from'].astype(int).tolist(), df['ali_to'].astype(int).tolist(),
                                          df['i_evalue'].tolist(), df['dom_bitscore'].tolist()):
        out.setdefault(pid, []).append((name, frm, to, tool, ev, sc))
    return {pid: run_dbcan_filter_overlaps(entries) for pid, entries in out.items()}


# ── Consensus ──────────────────────────────────────────────────────────────────

TOOL_LABELS = {"diamond": "diamond", "hmm": "dbcan_hmm", "sub": "dbcan_sub"}


def apply_consensus(diamond: dict, hmm: dict, sub: dict,
                    hmm_domains: pd.DataFrame | None, sub_domains: pd.DataFrame | None,
                    verbose: bool = True) -> pd.DataFrame:
    """Apply dbCAN consensus: keep proteins called by ≥2 of 3 tools, with run_dbcan's families.

    ``diamond``, ``hmm`` and ``sub`` ({protein: {family: (evalue, score)}}) decide which
    tools call a protein. ``hmm_domains`` and ``sub_domains`` (``_astra_passing_rows``
    tables) give the domains run_dbcan's recommendation chooses among
    (``run_dbcan_recommendation``); each recommended CAZy family becomes one row, with the
    i-evalue and score of its best recommended domain.

    Returns DataFrame with columns: protein_id, family_class, cazy_class, evalue, score,
    n_tools, tools (the tools calling the protein), family_tools (the tools whose hits
    include the family).
    """
    all_pids = set(diamond) | set(hmm) | set(sub)
    consensus: dict[str, list[str]] = {}
    for pid in all_pids:
        tools = [TOOL_LABELS[key] for key, src in (("diamond", diamond), ("hmm", hmm), ("sub", sub))
                 if pid in src]
        if len(tools) >= MIN_TOOLS:
            consensus[pid] = sorted(tools)

    hmm_entries = _filtered_entries(hmm_domains, set(consensus), "hmm")
    sub_entries = _filtered_entries(sub_domains, set(consensus), "sub")

    rows = []
    no_family = 0
    for pid, tools in consensus.items():
        chosen = run_dbcan_recommendation(hmm_entries.get(pid, []), sub_entries.get(pid, []),
                                          diamond.get(pid))
        families: dict[str, tuple] = {}
        for entry in chosen:
            fam = entry[0] if entry[3] == "diamond" else parse_hmm_family(entry[0])
            if fam is None:
                continue
            if fam not in families or entry[4] < families[fam][0]:
                families[fam] = (entry[4], entry[5])
        if not families:
            no_family += 1
        tools_str = ','.join(tools)
        for fam, (ev, sc) in families.items():
            support = [TOOL_LABELS[key] for key, src in (("diamond", diamond), ("hmm", hmm), ("sub", sub))
                       if fam in src.get(pid, {})]
            rows.append({
                'protein_id': pid,
                'family_class': fam,
                'cazy_class': get_cazy_class(fam) or '',
                'evalue': ev,
                'score': sc,
                'n_tools': len(tools),
                'tools': tools_str,
                'family_tools': ','.join(sorted(support)),
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
        print(f"  Consensus (≥{MIN_TOOLS} tools): {len(consensus):,} proteins, "
              f"{n_cons:,} with a recommended CAZy family ({no_family:,} without)")

    return df


# ── Main pipeline ──────────────────────────────────────────────────────────────

def _consensus_out(db_path: str | None, results_out: str | Path | None) -> Path:
    """The consensus table: ``results_out``, or cazyme_classification.tsv beside the database."""
    if results_out:
        return Path(results_out)
    if db_path:
        return Path(db_path).parent / "cazyme_classification.tsv"
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
    out = _consensus_out(db_path, results_out)
    side = side_files(db_path, results_out)

    # Prepare proteins (explicit FASTA, existing .faa beside the database, or extract from it)
    _log(f"[cazyme] Finding protein sequences...")
    prot_dir, n_proteins, is_temp = resolve_proteins(db_path, proteins, verbose)
    if prot_dir is None:
        if verbose:
            _log("No protein sequences found")
        return pd.DataFrame()
    _log(f"[cazyme] {n_proteins:,} proteins ready ({_time.time() - t_total:.1f}s)")

    try:
        fingerprint = protein_fingerprint(prot_dir)

        # Tool 1: DIAMOND, from a cache only when its stamp matches this input and release
        t1 = _time.time()
        _log(f"[cazyme] Tool 1: DIAMOND search...")
        d_stamp = diamond_stamp(fingerprint, dbcan)
        diamond_hits = _read_diamond_cache(side["diamond"], d_stamp, verbose=verbose)
        if diamond_hits is None:
            diamond_hits = run_diamond(prot_dir, dbcan, threads=threads, verbose=verbose)
            _save_diamond_cache(diamond_hits, side["diamond"], d_stamp, verbose=verbose)
        _log(f"[cazyme] Tool 1 done: {len(diamond_hits):,} hits ({_time.time() - t1:.1f}s)")

        # Tool 2: HMMER vs the release's own dbCAN.hmm, every protein. A raw table is
        # reused only with a matching stamp; Stage 04's registry table
        # (annotations/dbCAN_hits_df.tsv) carries none and is never reused.
        t2 = _time.time()
        _log(f"[cazyme] Tool 2: dbCAN HMM search...")
        hmm_scale = dbcan_ievalue_scale(dbcan_hmm, fingerprint["records"])
        hmm_order = hmm_profile_order(dbcan_hmm)
        h_stamp = hmm_table_stamp(fingerprint, dbcan_hmm)
        cached = [side["hmm"]] + ([Path(db_path).parent / "annotations" / "dbCAN_hits_df.tsv"] if db_path else [])
        reuse = next((t for t in cached if stamp_matches(t, h_stamp, "dbCAN.hmm table")), None)
        if reuse is not None:
            _log(f"[cazyme] Reusing dbCAN.hmm table {reuse}")
            hmm_domains = _astra_passing_rows(reuse, HMMER_IEVALUE, hmm_lengths(dbcan_hmm),
                                              order=hmm_order, ievalue_scale=hmm_scale)
            hmm_hits = _hits_from_rows([_best_per_family(hmm_domains)])
        else:
            hmm_hits, hmm_domains = run_astra_hmm(
                prot_dir, "dbCAN_HMM",
                hmm_path=dbcan_hmm,
                threads=threads, verbose=verbose,
                keep_tsv=side["hmm"],
                ievalue_scale=hmm_scale, order=hmm_order, domains=True,
            )
            write_stamp(side["hmm"], h_stamp)
        _log(f"[cazyme] Tool 2 done: {len(hmm_hits):,} hits ({_time.time() - t2:.1f}s)")

        # Tool 3: Astra/HMMER vs dbCAN-sub.hmm — only on proteins already
        # hit by ≥1 other tool.  Consensus needs ≥2, so proteins seen by
        # neither DIAMOND nor dbCAN HMM can never pass even with a sub hit.
        sub_hits: dict = {}
        sub_domains = None
        if has_sub:
            t3 = _time.time()
            candidate_pids = set(diamond_hits) | set(hmm_hits)
            _log(f"[cazyme] Tool 3: dbCAN-sub HMM on {len(candidate_pids):,} candidates...")
            if candidate_pids:
                sub_dir = _subset_fasta(prot_dir, candidate_pids, verbose=verbose)
                try:
                    sub_hits, sub_domains = run_astra_hmm(
                        sub_dir, "dbCAN_sub",
                        hmm_path=dbcan_sub,
                        threads=threads, verbose=verbose,
                        keep_tsv=side["sub"],
                        ievalue_scale=dbcan_ievalue_scale(dbcan_sub, count_fasta_records(sub_dir)),
                        order=hmm_profile_order(dbcan_sub), domains=True,
                    )
                finally:
                    shutil.rmtree(sub_dir, ignore_errors=True)
            _log(f"[cazyme] Tool 3 done: {len(sub_hits):,} hits ({_time.time() - t3:.1f}s)")

        # Consensus filter
        _log(f"[cazyme] Applying {MIN_TOOLS}-tool consensus filter...")
        results_df = apply_consensus(diamond_hits, hmm_hits, sub_hits, hmm_domains, sub_domains,
                                     verbose=verbose)
    finally:
        if is_temp:
            shutil.rmtree(prot_dir, ignore_errors=True)

    if results_df.empty:
        if verbose:
            print("No CAZymes passed consensus filter")
        return results_df

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
# HMMER reports i-evalue = P-value x Z. Aksha sets Z to the number of target sequences
# in the search; the parse multiplies by n_profiles / n_targets (dbcan_ievalue_scale)
# to reach run_dbcan's fixed Z = number of profiles in the HMM file. Splitting profiles
# keeps each search's target set, and so its Z, equal to the single-process run's: the
# whole proteome for dbCAN.hmm and the whole candidate subset for dbCAN-sub. Part i holds every N-th profile from profile i. Aksha writes
# each profile's rows together, in profile order, so the merge restores the
# single-process table row for row by a stable sort on profile position, and ties
# between rows resolve the same way too.

HMM_TSV_HEADER = ("sequence_id\thmm_name\tbitscore\tevalue\tc_evalue\ti_evalue\t"
                  "env_from\tenv_to\tdom_bitscore\tali_from\tali_to\thmm_from\thmm_to\n")


def _shard_dir(work_dir: Path, stage: str, index: int) -> Path:
    return Path(work_dir) / f"{stage}_{index:04d}"


_proteome_fingerprint = protein_fingerprint


def candidate_set_hash(candidates: set) -> str:
    """BLAKE2b-128 of the sorted candidate protein IDs: the dbCAN-sub target set."""
    h = hashlib.blake2b(digest_size=16)
    for pid in sorted(candidates):
        h.update(pid.encode())
        h.update(b"\n")
    return h.hexdigest()


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


def _parse_shard_file(kind: str, path: str, lengths: dict | None, order: dict | None,
                      ievalue_scale: float = 1.0, best: bool = False):
    """DIAMOND hits (a dict) for ``diamond``; for an HMM part table its passing domain rows
    (``_astra_passing_rows``), or with ``best`` one row per protein and family."""
    if kind == "diamond":
        return parse_diamond_tsv(Path(path))
    rows = _astra_passing_rows(Path(path), HMMER_IEVALUE, lengths, order=order, ievalue_scale=ievalue_scale)
    return _best_per_family(rows) if best else rows


def _parse_files(jobs: list[tuple], threads: int) -> list:
    """Parse (kind, path, profile lengths, profile order, i-evalue scale, best) jobs, one
    file per worker process."""
    if threads > 1 and len(jobs) > 1:
        with ProcessPoolExecutor(max_workers=min(threads, len(jobs))) as pool:
            return list(pool.map(_parse_shard_file, *zip(*[(k, str(p), l, o, z, b) for k, p, l, o, z, b in jobs])))
    return [_parse_shard_file(kind, str(path), lengths, order, scale, best)
            for kind, path, lengths, order, scale, best in jobs]


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


def shard_candidates(work_dir: Path, n_shards: int, dbcan: Path, n_records: int,
                     threads: int = 1) -> tuple[set, int, int]:
    """Proteins with a DIAMOND or dbCAN.hmm hit in any search shard (the dbCAN-sub subset),
    with the DIAMOND and dbCAN.hmm protein counts. ``n_records`` is the searched proteome size."""
    dbcan_hmm, _, _ = dbcan_hmm_files(dbcan)
    hmm_len = hmm_lengths(dbcan_hmm)
    scale = dbcan_ievalue_scale(dbcan_hmm, n_records, len(hmm_len))
    jobs = [("diamond", _shard_dir(work_dir, "search", i) / "diamond.tsv", None, None, 1.0, True)
            for i in range(n_shards)]
    jobs += [("hmm", _shard_dir(work_dir, "search", i) / "dbcan_hmm.tsv", hmm_len, None, scale, True)
             for i in range(n_shards)]
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
    (out / "DONE.json").unlink(missing_ok=True)  # a failed rerun must not leave an earlier record
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
    for stale in ("DONE.json", "dbcan_sub.tsv"):  # a failed rerun must not leave earlier output
        (out / stale).unlink(missing_ok=True)

    candidates, n_diamond, n_hmm = shard_candidates(work_dir, n_shards, dbcan, fingerprint["records"], threads)
    _log(f"[cazyme] sub shard {index}/{n_shards}: {len(candidates):,} candidates "
         f"({n_diamond:,} DIAMOND, {n_hmm:,} dbCAN.hmm)")
    sub_records = 0
    if has_sub and candidates:
        sub_dir = _subset_fasta(prot_dir, candidates, verbose=verbose)
        sub_records = count_fasta_records(sub_dir)  # Aksha's Z for this search
        try:
            _hmm_part_search(sub_dir, dbcan_sub, index, n_shards, "dbCAN_sub", threads,
                             out / "dbcan_sub.tsv", verbose)
        finally:
            shutil.rmtree(sub_dir, ignore_errors=True)

    info = {"stage": "sub", "index": index, "n_shards": n_shards, "proteins": fingerprint,
            "dbcan_dir": str(dbcan.resolve()), "candidates": len(candidates),
            "candidates_hash": candidate_set_hash(candidates),
            "sub_records": sub_records, "has_sub": has_sub,
            "seconds": {"total": round(_time.time() - t0, 1)}}
    _mark_done(out, info)
    _log(f"[cazyme] sub shard {index}/{n_shards} done: {info['seconds']}")
    return info


def merge_shards(work_dir: Path, n_shards: int, dbcan: Path, side: dict[str, Path],
                 threads: int = 1, verbose: bool = True) -> tuple[pd.DataFrame, int]:
    """Consensus over every shard's hits; writes the single-process run's side files
    (stamped DIAMOND cache and dbCAN.hmm table, raw dbCAN-sub table) to ``side``
    (``side_files``). Returns (table, protein records).

    The sub shards must agree on their candidate set, and that set must equal the one the
    current search shards give; every sub shard that searched (a release with dbCAN-sub
    and at least one candidate) must have left its table. Anything else stops the merge.
    """
    searched = finished_shards(work_dir, "search", n_shards)
    subbed = finished_shards(work_dir, "sub", n_shards)
    if subbed[0]["proteins"] != searched[0]["proteins"]:
        raise RuntimeError("sub shards searched a different protein input than the search shards")
    for key in ("candidates_hash", "sub_records", "has_sub"):
        if any(key not in info for info in subbed):
            raise RuntimeError(f"sub shards lack {key} (written by an older version); rerun the sub stage")
        if len({json.dumps(info[key]) for info in subbed}) > 1:
            raise RuntimeError(f"sub shards disagree on {key}; rerun the sub stage")
    dbcan_hmm, dbcan_sub, has_sub = dbcan_hmm_files(dbcan)
    if has_sub != subbed[0]["has_sub"]:
        raise RuntimeError(f"{dbcan} {'now has' if has_sub else 'no longer has'} a dbCAN-sub file; "
                           "rerun the sub stage")

    hmm_parts = [_shard_dir(work_dir, "search", i) / "dbcan_hmm.tsv" for i in range(n_shards)]
    sub_parts = [_shard_dir(work_dir, "sub", i) / "dbcan_sub.tsv" for i in range(n_shards)]
    hmm_len = hmm_lengths(dbcan_hmm)
    hmm_scale = dbcan_ievalue_scale(dbcan_hmm, searched[0]["proteins"]["records"], len(hmm_len))
    hmm_order = hmm_profile_order(dbcan_hmm)
    jobs = [("diamond", _shard_dir(work_dir, "search", i) / "diamond.tsv", None, None, 1.0, False)
            for i in range(n_shards)]
    jobs += [("hmm", p, hmm_len, hmm_order, hmm_scale, False) for p in hmm_parts]
    parsed = _parse_files(jobs, threads)
    diamond_hits = _fold_in_order(parsed[:n_shards])
    hmm_domains = _domain_table(parsed[n_shards:])
    hmm_hits = _hits_from_rows([_best_per_family(f) for f in parsed[n_shards:]])

    # the dbCAN-sub shards must have searched the candidates these search shards give
    candidates = set(diamond_hits) | set(hmm_hits)
    if candidate_set_hash(candidates) != subbed[0]["candidates_hash"]:
        raise RuntimeError(f"the search shards give {len(candidates):,} dbCAN-sub candidates, but the "
                           f"sub shards searched a different set ({subbed[0]['candidates']:,}); "
                           "rerun the sub stage")
    sub_hits: dict = {}
    sub_domains = None
    searched_sub = subbed[0]["has_sub"] and subbed[0]["candidates"] > 0
    if searched_sub:
        missing = [str(p) for p in sub_parts if not p.exists()]
        if missing:
            raise RuntimeError(f"{len(missing)} dbCAN-sub shard table(s) missing although the shards "
                               f"report a search: {', '.join(missing)}")
        sub_len, sub_order = hmm_lengths(dbcan_sub), hmm_profile_order(dbcan_sub)
        sub_scale = dbcan_ievalue_scale(dbcan_sub, subbed[0]["sub_records"], len(sub_len))
        sub_parsed = _parse_files([("sub", p, sub_len, sub_order, sub_scale, False) for p in sub_parts],
                                  threads)
        sub_domains = _domain_table(sub_parsed)
        sub_hits = _hits_from_rows([_best_per_family(f) for f in sub_parsed])
    _log(f"[cazyme] merged {n_shards} shards: {len(diamond_hits):,} DIAMOND, {len(hmm_hits):,} "
         f"dbCAN.hmm, {len(sub_hits):,} dbCAN-sub proteins ({subbed[0]['candidates']:,} candidates)")

    proteins = searched[0]["proteins"]
    _save_diamond_cache(diamond_hits, side["diamond"], diamond_stamp(proteins, dbcan), verbose=verbose)
    clear_stamp(side["hmm"])
    combine_hmm_tables(hmm_parts, dbcan_hmm, side["hmm"])
    write_stamp(side["hmm"], hmm_table_stamp(proteins, dbcan_hmm))
    if searched_sub:
        combine_hmm_tables(sub_parts, dbcan_sub, side["sub"])

    _log(f"[cazyme] Applying {MIN_TOOLS}-tool consensus filter...")
    table = apply_consensus(diamond_hits, hmm_hits, sub_hits, hmm_domains, sub_domains, verbose=verbose)
    return table, proteins["records"]


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
        out = _consensus_out(db_path, results_out)
        results_df, n_proteins = merge_shards(work_dir, n_shards, dbcan, side_files(db_path, results_out),
                                              threads, verbose)
        if results_df.empty:
            _log("No CAZymes passed consensus filter")
            return results_df
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


def _cazy_predicates(gen, families) -> set[str]:
    """What the predicate generator derives from CAZy annotations of ``families``
    (sharur/predicates/generator.py, source 'cazy', then its hierarchy expansion)."""
    preds = {"cazy_annotated"} if families else set()
    for fam in families:
        preds.update(get_predicates_for_cazy(fam))
        preds.add(f"cazy:{fam}")
    return gen._expand_hierarchy(preds)


def _rewrite_cazy_predicates(db, results_df: pd.DataFrame) -> dict[str, int]:
    """Rewrite the CAZy-derived predicates of every protein that carries a cazy: predicate
    or is in ``results_df``.

    A protein's earlier CAZy-derived predicates (those of its cazy: families, expanded
    through the predicate hierarchy) are removed, except predicates its other annotation
    sources also give (``hydrolase`` from a PFAM glycosyl hydrolase domain stays); then
    the predicates of its new families are added. Runs inside the caller's transaction.
    """
    from sharur.predicates.generator import AnnotationRecord, PredicateGenerator
    gen = PredicateGenerator(predict_topology=False)

    new_families = results_df.groupby('protein_id')['family_class'].apply(lambda s: sorted(set(s))).to_dict()
    current = {pid: set(preds or []) for pid, preds in db.execute(
        "SELECT protein_id, predicates FROM protein_predicates "
        "WHERE len(list_filter(predicates, x -> starts_with(x, 'cazy:'))) > 0 "
        "   OR protein_id IN (SELECT UNNEST(?::VARCHAR[]))", [sorted(new_families)]).fetchall()}
    affected = sorted(current)

    # predicates the proteins' other annotation sources give, which a CAZy rewrite keeps
    other: dict[str, set] = defaultdict(set)
    for pid, source, acc, name, desc, ev, sc in db.execute(
            "SELECT protein_id, source, accession, name, description, evalue, score FROM annotations "
            "WHERE source != 'cazy' AND protein_id IN (SELECT UNNEST(?::VARCHAR[]))", [affected]).fetchall():
        other[pid] |= gen._predicates_from_annotation(
            AnnotationRecord(source=source, accession=acc, name=name, description=desc, evalue=ev, score=sc))

    ids, lists = [], []
    stripped = 0
    for pid, existing in current.items():
        old_fams = {p.split(":", 1)[1] for p in existing if p.startswith("cazy:")}
        kept = gen._expand_hierarchy(other[pid]) if pid in other else set()
        merged = (existing - (_cazy_predicates(gen, old_fams) - kept)) | _cazy_predicates(gen, new_families.get(pid, []))
        if merged != existing:
            ids.append(pid)
            lists.append(sorted(merged))
            stripped += pid not in new_families
    if ids:
        # lists travel as unit-separator-joined strings: no pyarrow needed for a list column
        db.register("cazy_predicate_update", pd.DataFrame(
            {"protein_id": ids, "joined": ["\x1f".join(p) for p in lists]}))
        db.execute("UPDATE protein_predicates SET predicates = CASE WHEN u.joined = '' THEN []::VARCHAR[] "
                   "ELSE string_split(u.joined, chr(31)) END, updated_at = CURRENT_TIMESTAMP "
                   "FROM cazy_predicate_update u WHERE protein_predicates.protein_id = u.protein_id")
        db.unregister("cazy_predicate_update")
    no_row = len(set(new_families) - set(current))
    if no_row:
        _log(f"[cazyme] {no_row:,} consensus proteins have no protein_predicates row; left without predicates")
    return {"updated": len(ids), "stripped_only": stripped}


def load_cazyme_results(
    db_path: str,
    results_df: pd.DataFrame,
    update_predicates: bool = True,
    n_proteins: int | None = None,
    verbose: bool = True,
) -> dict[str, int]:
    """Replace source='cazy' annotations with a consensus table and rewrite CAZy predicates.

    ``results_df`` is the table ``classify_cazymes`` writes (protein_id,
    family_class, cazy_class, evalue, score, n_tools, tools[, family_tools]). Rows naming
    proteins absent from the database are dropped and counted. The whole replacement runs
    in one transaction: a failure leaves the earlier annotations and predicates in place.
    """
    _log(f"[cazyme] Opening DB for annotation loading...")
    results_df = results_df.assign(protein_id=results_df['protein_id'].astype(str))
    db = duckdb.connect(db_path)
    try:
        db.execute("BEGIN TRANSACTION")
        known = set(r[0] for r in db.execute(
            "SELECT protein_id FROM proteins WHERE protein_id IN (SELECT UNNEST(?::VARCHAR[]))",
            [sorted(results_df['protein_id'].unique())]).fetchall())
        unknown = int((~results_df['protein_id'].isin(known)).sum())
        if unknown:
            _log(f"[cazyme] Dropping {unknown:,} rows for proteins absent from the database")
            results_df = results_df[results_df['protein_id'].isin(known)].copy()
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

        ann = results_df[['protein_id', 'family_class', 'evalue', 'score']].copy()
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

        if update_predicates:
            t_pred = _time.time()
            _log(f"[cazyme] Updating predicates...")
            pred_stats = _rewrite_cazy_predicates(db, results_df)
            _log(f"[cazyme] Updated predicates for {pred_stats['updated']:,} proteins, "
                 f"{pred_stats['stripped_only']:,} of them no longer CAZymes "
                 f"({_time.time() - t_pred:.1f}s)")
        db.execute("COMMIT")
    except BaseException:
        try:
            db.execute("ROLLBACK")
        except duckdb.Error:
            pass
        raise
    finally:
        db.close()

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
        stats = load_cazyme_results(args.db, pd.read_csv(args.load_results, sep='\t',
                                                         dtype={'protein_id': str}),
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
