#!/usr/bin/env python3
"""
CAZyme Classification Pipeline (dbCAN 3-tool consensus)

Recapitulates the standard dbCAN annotation pipeline:
  Tool 1: DIAMOND blastp vs CAZy.dmnd        (e-value ≤ 1e-102)
  Tool 2: Astra/HMMER vs dbCAN.hmm           (i-evalue ≤ 1e-15, HMM coverage ≥ 0.35)
  Tool 3: Astra/HMMER vs dbCAN-sub.hmm       (i-evalue ≤ 1e-15, HMM coverage ≥ 0.35)
  Consensus: retain genes called by ≥ 2 of 3 tools

HMMER searches are delegated to Astra (astra search --hmm_in), which handles
parallelization and PyHMMER optimization. Only DIAMOND is run directly.

If dbCAN-sub.hmm is not available, falls back to 2-tool mode (DIAMOND + dbCAN HMM)
and requires both tools to agree.

Usage:
    python scripts/classify_cazymes.py --db data/my_dataset/sharur.duckdb --threads 8

    # Compute on a cluster against a read-only database, load elsewhere:
    python scripts/classify_cazymes.py --db DB --threads N --results-only --out cazy.tsv
    python scripts/classify_cazymes.py --db DB --load-results cazy.tsv

References:
    Zheng J, et al. (2023) dbCAN3: automated carbohydrate-active enzyme and
    substrate annotation. Nucleic Acids Res 51(W1):W115-W121.
"""

import argparse
import os
import shutil
import subprocess
import tempfile
import re
from pathlib import Path
from collections import defaultdict

import time as _time

import duckdb
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


# ── Protein extraction ─────────────────────────────────────────────────────────

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
    """Load cached DIAMOND results from cazyme_diamond_hits.tsv if present.

    The cache has its own file: the consensus table (cazyme_classification.tsv)
    merges HMM families and must never be read back as DIAMOND evidence.

    Returns {protein_id: {family_class: (evalue, bitscore)}} or None if no cache.
    """
    tsv = Path(db_path).parent / DIAMOND_CACHE
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
        df = pd.DataFrame(rows)
        out = data_dir / DIAMOND_CACHE
        df.to_csv(out, sep='\t', index=False)
        if verbose:
            print(f"  Saved DIAMOND cache to {out.name}")


def run_diamond(prot_dir: str, dbcan_path: Path, threads: int = 4,
                evalue: float = DIAMOND_EVALUE,
                verbose: bool = True) -> dict[str, dict[str, tuple]]:
    """DIAMOND blastp against CAZy.dmnd.

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
    with open(combined, "w") as fout:
        for faa in faa_files:
            with open(faa) as fin:
                header = None
                seq_lines: list[str] = []
                for line in fin:
                    if line.startswith(">"):
                        # Write previous record if it had sequence
                        if header is not None:
                            if seq_lines:
                                fout.write(header)
                                for sl in seq_lines:
                                    fout.write(sl)
                                n_written += 1
                            else:
                                n_skipped += 1
                        header = line
                        seq_lines = []
                    else:
                        cleaned = line.rstrip("\n").replace("*", "")
                        if cleaned:
                            seq_lines.append(cleaned + "\n")
                # Write last record
                if header is not None:
                    if seq_lines:
                        fout.write(header)
                        for sl in seq_lines:
                            fout.write(sl)
                        n_written += 1
                    else:
                        n_skipped += 1
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
        "--sensitive",
    ]

    # stderr goes to a file: an unread stderr pipe can fill and stall DIAMOND
    err_handle = tempfile.TemporaryFile(mode="w+")
    try:
        proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=err_handle,
            text=True, bufsize=1,
        )
    except FileNotFoundError:
        _log("Warning: diamond not found on PATH")
        return {}

    hits: dict[str, dict[str, tuple]] = defaultdict(dict)
    n_lines = 0
    try:
        for line in proc.stdout:
            line = line.rstrip('\n')
            if not line:
                continue
            fields = line.split('\t')
            if len(fields) < 5:
                continue
            qid, sid = fields[0], fields[1]
            ev, bs = float(fields[3]), float(fields[4])
            for fam in parse_cazy_families(sid):
                if fam not in hits[qid] or ev < hits[qid][fam][0]:
                    hits[qid][fam] = (ev, bs)
            n_lines += 1
            if verbose and n_lines % 500_000 == 0:
                _log(f"  [diamond] Parsed {n_lines:,} hits, {len(hits):,} proteins so far...")

        proc.wait()
    except BaseException:
        proc.kill()
        raise

    err_handle.seek(0)
    stderr_out = err_handle.read()
    err_handle.close()
    shutil.rmtree(work_dir, ignore_errors=True)
    if proc.returncode != 0:
        # fatal: an empty DIAMOND set would silently turn the 3-tool consensus into HMM-only
        raise RuntimeError(f"DIAMOND failed (rc={proc.returncode}): {stderr_out[-500:]}")

    if verbose and stderr_out:
        for line in stderr_out.strip().split('\n')[-3:]:
            _log(f"  {line}")

    if verbose:
        _log(f"  DIAMOND: {len(hits):,} proteins with CAZy hits ({n_lines:,} hit lines)")
    return dict(hits)


# ── Tools 2 & 3: Astra/HMMER ──────────────────────────────────────────────────

def hmm_lengths(hmm_path: Path) -> dict[str, int]:
    """Profile NAME -> LENG from a HMMER3 text file (for the dbCAN coverage filter)."""
    lengths: dict[str, int] = {}
    name = None
    with open(hmm_path) as fh:
        for line in fh:
            if line.startswith("NAME "):
                name = line.split(None, 1)[1].strip()
            elif line.startswith("LENG ") and name is not None:
                lengths[name] = int(line.split()[1])
                name = None
    return lengths


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

    Returns {protein_id: {family_class: (best_i_evalue, best_dom_bitscore)}}.
    """
    try:
        df = pd.read_csv(path, sep='\t')
    except Exception as e:
        print(f"Warning: Failed to parse Astra TSV: {e}")
        return {}

    if df.empty:
        return {}

    # Filter by i-evalue
    df = df[df['i_evalue'] < ievalue_threshold].copy()
    if df.empty:
        return {}
    if lengths and {'hmm_from', 'hmm_to'} <= set(df.columns):
        leng = df['hmm_name'].map(lengths)
        if leng.isna().any():
            raise ValueError(f"{path.name}: {int(leng.isna().sum())} hits name profiles absent "
                             "from the HMM file's NAME/LENG records")
        before = len(df)
        df = df[(df['hmm_to'] - df['hmm_from']) / leng >= min_coverage].copy()
        _log(f"  [{path.name}] HMM coverage >= {min_coverage}: {len(df):,} of {before:,} domain hits")
        if df.empty:
            return {}
    else:
        _log(f"  [{path.name}] no HMM coordinates or profile lengths: i-evalue filter only")

    # Extract family from hmm_name
    df['family'] = df['hmm_name'].apply(parse_hmm_family)
    df = df.dropna(subset=['family'])
    if df.empty:
        return {}

    # Group by protein × family, keep best i-evalue
    result: dict[str, dict[str, tuple]] = {}
    for _, row in df.iterrows():
        pid = row['sequence_id']
        fam = row['family']
        ev = row['i_evalue']
        sc = row['dom_bitscore']

        if pid not in result:
            result[pid] = {}
        if fam not in result[pid] or ev < result[pid][fam][0]:
            result[pid][fam] = (ev, sc)

    return result


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

    df = pd.DataFrame(rows) if rows else pd.DataFrame()

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

def classify_cazymes(
    db_path: str,
    threads: int = 4,
    update_predicates: bool = True,
    verbose: bool = True,
    dbcan_dir: str | Path | None = None,
    write_db: bool = True,
    results_out: str | Path | None = None,
) -> pd.DataFrame:
    """Full dbCAN 3-tool consensus CAZyme classification.

    Runs DIAMOND + Astra/HMMER (dbCAN.hmm) + Astra/HMMER (dbCAN-sub.hmm),
    applies ≥2-tool consensus filter, loads annotations as source='cazy',
    and updates predicates.

    Args:
        db_path: Path to sharur.duckdb
        threads: CPU threads for DIAMOND and Astra
        update_predicates: Whether to update protein_predicates
        verbose: Print progress
        dbcan_dir: dbCAN release directory (default: $SHARUR_DBCAN_DIR, then
            the standard locations)
        write_db: False computes the consensus without opening the database
            for writing (load it later with ``load_cazyme_results``)
        results_out: Where to write the consensus table (default:
            cazyme_classification.tsv beside the database)

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

    dbcan_hmm = dbcan / "dbCAN.hmm"
    # dbCAN releases name the sub-family profiles dbCAN_sub.hmm; older docs used dbCAN-sub.hmm.
    dbcan_sub = next((dbcan / n for n in ("dbCAN-sub.hmm", "dbCAN_sub.hmm") if (dbcan / n).exists()),
                     dbcan / "dbCAN-sub.hmm")
    has_sub = dbcan_sub.exists() and dbcan_sub.stat().st_size > 0

    if verbose:
        _log(f"dbCAN databases: {dbcan}")
        _log(f"  CAZy.dmnd:     yes")
        _log(f"  dbCAN.hmm:     yes ({dbcan_hmm.stat().st_size / 1e6:.0f} MB)")
        _log(f"  dbCAN-sub.hmm: {'yes' if has_sub else 'NO — 2-tool mode'}"
             f"{f' ({dbcan_sub.stat().st_size / 1e6:.0f} MB)' if has_sub else ''}")

    # Prepare proteins (find existing .faa or extract from DB)
    _log(f"[cazyme] Finding protein sequences...")
    prot_dir, n_proteins, is_temp = find_or_extract_proteins(db_path, verbose)
    if prot_dir is None:
        if verbose:
            _log("No protein sequences found")
        return pd.DataFrame()
    _log(f"[cazyme] {n_proteins:,} proteins ready ({_time.time() - t_total:.1f}s)")

    data_dir = Path(db_path).parent
    source_dir = data_dir / "source"

    try:
        # Tool 1: DIAMOND (use cached results if available)
        t1 = _time.time()
        _log(f"[cazyme] Tool 1: DIAMOND search...")
        diamond_hits = load_diamond_cache(db_path, verbose=verbose)
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

    out = Path(results_out) if results_out else Path(db_path).parent / "cazyme_classification.tsv"
    results_df.to_csv(out, sep='\t', index=False)
    _log(f"[cazyme] Consensus table: {out}")
    if write_db:
        load_cazyme_results(db_path, results_df, update_predicates=update_predicates,
                            n_proteins=n_proteins, verbose=verbose)
    _log(f"[cazyme] Total {_time.time() - t_total:.1f}s")
    return results_df


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
        "--db", required=True,
        help="Path to sharur.duckdb",
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

    args = parser.parse_args()

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
    )

    return 0 if not results.empty else 1


if __name__ == "__main__":
    exit(main())
