# dbCAN 3-Tool Consensus CAZyme Annotation Pipeline

## Overview

The standard dbCAN annotation pipeline uses a **3-tool consensus** approach to minimize false positives in carbohydrate-active enzyme (CAZyme) classification. Only proteins identified by at least 2 of 3 independent methods are retained.

**Reference:** Zheng J, et al. (2023) dbCAN3: automated carbohydrate-active enzyme and substrate annotation. *Nucleic Acids Res* 51(W1):W115-W121.

### The Three Tools

| Step | Tool | Database | Threshold | Purpose |
|------|------|----------|-----------|---------|
| 1 | DIAMOND blastp (default mode, `-k 1`) | `CAZy.dmnd` | e-value <= 1e-102 | Sequence homology to known CAZymes |
| 2 | HMMER (via Aksha) | `dbCAN.hmm` | i-evalue <= 1e-15, HMM coverage >= 0.35 | HMM profile matching to CAZy family models |
| 3 | HMMER (via Aksha) | `dbCAN-sub.hmm` | i-evalue <= 1e-15, HMM coverage >= 0.35 | Substrate-specific HMM sub-family models |

These are the run_dbcan settings, including DIAMOND's default sensitivity mode.

### Consensus Rule

- A protein must be called by **>= 2 of 3 tools** to pass consensus.
- If `dbCAN-sub.hmm` is unavailable, falls back to 2-tool mode (DIAMOND + dbCAN HMM), requiring **both** to agree.
- Step 3 runs only on the **union** of hit sets from steps 1 and 2 (optimization: proteins unseen by either tool can never reach 2/3 consensus even with a sub hit).

### CAZy Family Classes

| Prefix | Class | Description |
|--------|-------|-------------|
| GH | Glycoside Hydrolases | Break glycosidic bonds |
| GT | Glycosyltransferases | Form glycosidic bonds |
| PL | Polysaccharide Lyases | Cleave polysaccharides via beta-elimination |
| CE | Carbohydrate Esterases | Remove ester modifications |
| AA | Auxiliary Activities | Redox enzymes (lytic polysaccharide monooxygenases, etc.) |
| CBM | Carbohydrate-Binding Modules | Non-catalytic binding domains |

---

## Required Database Files

All database files should reside in `data/dbcan_db/`.

| File | Size | Required? | Notes |
|------|------|-----------|-------|
| `CAZy.dmnd` | ~2 GB | Yes | DIAMOND database for step 1 |
| `dbCAN.hmm` (+ `hmmpress` index) | ~130 MB | Yes | HMM profiles for step 2 |
| `dbCAN-sub.hmm` (+ `hmmpress` index) | ~2.3 GB | Recommended | Substrate HMM profiles for step 3; the release names it `dbCAN_sub.hmm`, and both names are accepted |
| `fam-substrate-mapping.tsv` | ~100 KB | Optional | Maps CAZy families to substrates |

### How to Download / Update

dbCAN distributes its databases from an AWS S3 bucket (`s3://dbcan/`, public,
region us-west-2); the release current in September 2026 is `db_v5-2_9-13-2025`.
`run_dbcan database --db_dir data/dbcan_db --aws_s3` downloads the same files.

```bash
mkdir -p data/dbcan_db && cd data/dbcan_db/
BASE=https://dbcan.s3.us-west-2.amazonaws.com/db_v5-2_9-13-2025

curl -fLO "$BASE/CAZy.dmnd"                    # step 1 (prebuilt DIAMOND database)
curl -fLO "$BASE/dbCAN.hmm" && hmmpress dbCAN.hmm   # step 2
curl -fLO "$BASE/dbCAN_sub.hmm" && hmmpress dbCAN_sub.hmm   # step 3, ~2.3 GB
curl -fLO "$BASE/fam-substrate-mapping.tsv"    # substrate mapping

# or everything in the release (includes PUL, TCDB, TF and peptidase databases):
aws s3 cp s3://dbcan/db_v5-2_9-13-2025/ . --no-sign-request --recursive
```

### Aksha Integration

Steps 2 and 3 pass the release's own `dbCAN.hmm` and `dbCAN_sub.hmm` to `aksha search --hmm_in`,
so all three tools come from one dbCAN release. The `"dbCAN"` entry in Aksha's HMM registry can
hold an older release and is left unused.

---

## Running the Pipeline

### Script Location

**Primary script:** `scripts/classify_cazymes.py`

This is the canonical, up-to-date implementation that handles all 3 tools, consensus filtering, database loading, and predicate updates.

**Legacy script:** `src/ingest/dbcan_cazyme.py` (DIAMOND-only, do NOT use for new datasets)

### Usage

```bash
# Standard run (3-tool consensus)
python scripts/classify_cazymes.py --db data/my_dataset/sharur.duckdb --threads 12

# Skip predicate updates (annotation-only)
python scripts/classify_cazymes.py --db data/my_dataset/sharur.duckdb --threads 12 --no-update

# Quiet mode
python scripts/classify_cazymes.py --db data/my_dataset/sharur.duckdb --threads 12 --quiet

# Cluster compute against a read-only database, loaded later
python scripts/classify_cazymes.py --db DB --threads 24 --dbcan-dir DIR --results-only --out cazy.tsv
python scripts/classify_cazymes.py --db DB --load-results cazy.tsv

# A protein FASTA (file or directory of .faa) needs no database
python scripts/classify_cazymes.py --proteins all.faa --threads 24 --results-only --out cazy.tsv
```

### What the Script Does

1. **Locates databases** in `--dbcan-dir`, `$SHARUR_DBCAN_DIR`, `data/dbcan_db/`, `data/reference/dbcan/`, or `~/.sharur/dbcan/`
2. **Finds proteins:** `--proteins`, else `stage03_prodigal/genomes/all_protein_symlinks/*.faa` or `source/proteins_db_ids.faa` beside the database, else extraction from the database
3. **Tool 1 (DIAMOND):** Runs `diamond blastp` in default mode with `-k 1` against `CAZy.dmnd` at e-value <= 1e-102 (the run_dbcan settings), or loads cached results from `cazyme_diamond_hits.tsv` if present (refiltered to the current cut)
4. **Tool 2 (Aksha/HMMER):** Runs `aksha search --hmm_in dbCAN.hmm` on every protein, or loads cached results from `annotations/dbCAN_hits_df.tsv`
5. **Tool 3 (Aksha/HMMER):** Subsets proteins to union of step 1+2 hits, runs `aksha search --hmm_in dbCAN_sub.hmm` against the subset
6. **Consensus:** Keeps proteins called by >= 2 tools; records which tools agreed
7. **Saves results** to `cazyme_classification.tsv` (or `--out`), one row per protein and family, sorted by protein and family
8. **Loads annotations** into DuckDB as `source='cazy'` (clears any existing CAZy annotations first), unless `--results-only`
9. **Updates predicates** (adds `cazy:GH13`, `carbohydrate_active`, etc.)

### Caching Behavior

The script keeps each tool's evidence beside the database (beside `--out` for a `--proteins` run):

| File | Contents | Reused by |
|------|----------|-----------|
| `cazyme_diamond_hits.tsv` | DIAMOND hits (protein_id, family_class, cazy_class, evalue, bitscore) | Step 1 of the next run |
| `annotations/dbCAN_hits_df.tsv` | Stage 04 dbCAN HMM hits, copied there by Stage 07 | Step 2 |
| `cazyme_dbcan_hmm_hits.tsv` | Raw Aksha table for step 2 | Audit |
| `cazyme_dbcan_sub_hits.tsv` | Raw Aksha table for step 3 | Audit |

The consensus table merges HMM families and is never read back as DIAMOND evidence.

### Integration with Stage 07

The `07_build_knowledge_base.py` ingest pipeline calls `classify_cazymes.py` automatically as part of the build process:

1. During `_load_annotations()`, raw dbCAN HMM TSV files are **skipped** (the line `if "dbcan" in tsv.as_posix().lower(): continue`)
2. During `_classify_cazymes()`, the full 3-tool consensus script is invoked
3. If stage04 Astra results include a `dbCAN_hits_df.tsv`, it is copied to `annotations/` for caching

This ensures that only consensus-filtered results enter the database as `source='cazy'`.

---

## Sharded Runs (Slurm arrays)

A large proteome splits into N shards that run as two Slurm arrays and one merge job.
`scripts/submit_cazy_shards.sh` writes the three job scripts to `WORK_DIR/slurm/` and submits
them chained with `afterok`:

```bash
scripts/submit_cazy_shards.sh -n 16 -w WORK_DIR -m 16 -e setup/env.sh \
    -s "-p standard -c 24" -- \
    --proteins all.faa --dbcan-dir "$SHARUR_DBCAN_DIR" --results-only --out cazy.tsv
```

| Job | Task i does | Starts after |
|-----|-------------|--------------|
| `cazy_search` (array 0..N-1) | `dbCAN.hmm` profile part i on every protein; task 0 also runs DIAMOND on every protein | submission |
| `cazy_sub` (array 0..N-1) | `dbCAN_sub.hmm` profile part i on the candidate subset built from every search shard | all search tasks |
| `cazy_merge` | Rebuilds each tool's single-process hits and tables, applies the consensus, writes `--out` and the side files above | all sub tasks |

Each step is also a plain command (`--stage search|sub --shard-index i --n-shards N --work-dir W`,
then `--stage merge --n-shards N --work-dir W --out cazy.tsv`), so a shard can be rerun alone.
Every shard searches afresh; the caches above serve single-process runs.

**Why no step splits the proteins.** HMMER reports i-evalue = P-value x Z, where Z is the number
of target sequences in the search (pyhmmer's default, which Aksha uses). A search over 1/N of the
proteins therefore reports i-evalues N times smaller and passes extra hits at the 1e-15 cut. Each
HMM part searches the same targets as the single-process run (every protein for `dbCAN.hmm`, the
whole candidate subset for `dbCAN_sub.hmm`), so Z, every i-evalue and every call are unchanged.
Part i holds every N-th profile from profile i, which spreads each family's consecutive subfamily
profiles over all parts and balances the work. Aksha writes each profile's rows together in
profile order; the merge reduces each part in parallel and breaks ties by profile position and row,
which reproduces the single-process table and its calls exactly.

DIAMOND chooses its seed-search algorithm from the query size: double-indexed for the 1,976,915
SR-VP archaeal proteins, query-indexed for a sixteenth of them. The top target (`-k 1`) then
differs for 22 of 19,152 proteins, while two runs on the whole set agree byte for byte. DIAMOND
therefore runs once on every protein, in search task 0; in default mode that takes under 3 minutes
for 2 million proteins. The candidate subset depends on every search shard, which is why dbCAN-sub
runs as a second array.

`tests/test_cazyme_sharding.py` checks that the merged table equals the single-process table
byte for byte, and that a protein-split HMM search changes calls.

Measured on biotite (24 CPUs per job) for the SR-VP archaeal set, 1,976,915 proteins and 38,883
dbCAN-sub candidates, with all 16 shards running at once:

| Stage | 16 shards | Single process |
|-------|-----------|----------------|
| search (DIAMOND + dbCAN.hmm) | 5.5 min; task 0 with DIAMOND 329 s, other tasks median 99 s | 167 s + 811 s |
| dbCAN-sub | 7.6 min; tasks 387-450 s | 105 min |
| merge | 3.5 min | |
| total | 16.6 min | about 2 h |

dbCAN-sub dominates and scales with the candidate count, so the shard count is best set by the
number of jobs that can run at once.

---

## Output Format

### DuckDB Annotations (source='cazy')

```sql
SELECT protein_id, accession, name, evalue, score
FROM annotations
WHERE source = 'cazy'
LIMIT 5;
```

Each row represents one protein-family pair. A protein with multiple CAZy domains (e.g., GH5 + CBM3) will have multiple rows.

| Column | Content |
|--------|---------|
| `source` | `'cazy'` |
| `accession` | CAZy family (e.g., `GH13`, `GT2`, `CBM48`) |
| `name` | Same as accession |
| `evalue` | Best e-value from whichever tool provided it |
| `score` | Best score from whichever tool provided it |

### cazyme_classification.tsv

Detailed results file saved alongside the database:

| Column | Description |
|--------|-------------|
| `protein_id` | Protein identifier |
| `family_class` | CAZy family (e.g., `GH13`) |
| `cazy_class` | Class abbreviation (e.g., `GH`, `GT`) |
| `evalue` | Best e-value across tools |
| `score` | Best score across tools |
| `n_tools` | Number of tools that detected this protein (2 or 3) |
| `tools` | Comma-separated list of agreeing tools |

### Predicates

Each CAZyme protein receives:
- `cazy:{family}` (e.g., `cazy:GH13`, `cazy:GT2`) -- specific family tag
- Class-level predicates from `sharur/predicates/mappings/cazy_map.py` (e.g., `carbohydrate_active`, `glycosyl_transferase`, `glycoside_hydrolase`)

---

## Troubleshooting

### "dbCAN databases not found"

Ensure `data/dbcan_db/` contains both `CAZy.dmnd` and `dbCAN.hmm`. The script checks multiple locations but requires at least these two files.

### dbCAN-sub not running (2-tool fallback)

If `dbCAN-sub.hmm` is missing or zero-size, the script falls back to 2-tool mode (DIAMOND + dbCAN HMM only). This requires **both** tools to agree, which is more conservative but may miss some valid CAZymes.

To enable 3-tool mode:
```bash
# Download dbCAN-sub.hmm (see above)
# Then hmmpress it:
cd data/dbcan_db/
hmmpress dbCAN-sub.hmm
```

### Long runtimes

The script sets no time limit on DIAMOND or Aksha, and a failed DIAMOND or HMM search stops the
run. For proteomes of a million proteins or more, use a sharded run (above).

### Existing annotations

Running the script clears all existing `source='cazy'` annotations before loading new ones. This is safe to re-run.

---

