# dbCAN 3-Tool Consensus CAZyme Annotation Pipeline

## Overview

The standard dbCAN annotation pipeline uses a **3-tool consensus** approach to minimize false positives in carbohydrate-active enzyme (CAZyme) classification. Only proteins identified by at least 2 of 3 independent methods are retained.

**Reference:** Zheng J, et al. (2023) dbCAN3: automated carbohydrate-active enzyme and substrate annotation. *Nucleic Acids Res* 51(W1):W115-W121.

### The Three Tools

| Step | Tool | Database | Threshold | Purpose |
|------|------|----------|-----------|---------|
| 1 | DIAMOND blastp | `CAZy.dmnd` | e-value <= 1e-18 | Sequence homology to known CAZymes |
| 2 | HMMER (via Astra) | `dbCAN.hmm` | i-evalue <= 1e-15 | HMM profile matching to CAZy family models |
| 3 | HMMER (via Astra) | `dbCAN-sub.hmm` | i-evalue <= 1e-15 | Substrate-specific HMM sub-family models |

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

| File | Size | Status | Required? | Notes |
|------|------|--------|-----------|-------|
| `CAZy.dmnd` | ~2.1 GB | **Present** | Yes | DIAMOND database for step 1 |
| `dbCAN.hmm` | ~130 MB | **Present** (hmmpress'd) | Yes | HMM profiles for step 2 |
| `dbCAN.hmm.h3{f,i,m,p}` | ~140 MB total | **Present** | Yes (auto-generated) | hmmpress index files |
| `dbCAN-sub.hmm` | ~2.3 GB | **Present** (NOT hmmpress'd) | Recommended | Substrate HMM profiles for step 3 |
| `dbCAN-sub.hmm.h3{f,i,m,p}` | — | **Missing** | Needed for step 3 | Must run `hmmpress dbCAN-sub.hmm` |
| `fam-substrate-mapping.tsv` | ~92 KB | **Present** | Optional | Maps CAZy families to substrates |

### How to Download / Update

dbCAN databases are hosted at: https://bcb.unl.edu/dbCAN2/download/Databases/

To download the latest version (V13 as of early 2026):

```bash
cd data/dbcan_db/

# dbCAN HMM profiles (step 2)
wget https://bcb.unl.edu/dbCAN2/download/Databases/V13/dbCAN-HMMdb-V13.txt
mv dbCAN-HMMdb-V13.txt dbCAN.hmm
hmmpress dbCAN.hmm

# DIAMOND database (step 1)
wget https://bcb.unl.edu/dbCAN2/download/Databases/V13/CAZyDB.fa
diamond makedb --in CAZyDB.fa --db CAZy.dmnd

# dbCAN-sub HMM profiles (step 3) — large file, ~2.3 GB
wget https://bcb.unl.edu/dbCAN2/download/Databases/V13/dbCAN_sub.hmm
mv dbCAN_sub.hmm dbCAN-sub.hmm
hmmpress dbCAN-sub.hmm

# Substrate mapping
wget https://bcb.unl.edu/dbCAN2/download/Databases/V13/fam-substrate-mapping-08252022.tsv
mv fam-substrate-mapping-08252022.tsv fam-substrate-mapping.tsv
```

Alternatively, using AWS (if still available):
```bash
aws s3 cp s3://dbcan/db_v5-2_9-13-2025/ data/dbcan_db/ --no-sign-request --recursive
```

### Astra Integration

The `dbCAN.hmm` file is registered in Astra (`~/.config/Astra/hmm_databases.json`) as `"dbCAN"` and symlinked from `~/.config/Astra/dbCAN/dbCAN.hmm`. This allows step 2 to use `--installed_hmms dbCAN`.

Step 3 (`dbCAN-sub.hmm`) is NOT registered in Astra. The classify script passes it via `--hmm_in` instead.

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
```

### What the Script Does

1. **Locates databases** in `data/dbcan_db/`, `data/reference/dbcan/`, or `~/.sharur/dbcan/`
2. **Finds or extracts proteins** from the DuckDB database (prefers `source/proteins_db_ids.faa` if available)
3. **Tool 1 (DIAMOND):** Runs `diamond blastp --sensitive` against `CAZy.dmnd` with e-value <= 1e-18, or loads cached results from `cazyme_classification.tsv` if present
4. **Tool 2 (Astra/HMMER):** Runs `astra search --installed_hmms dbCAN` against `dbCAN.hmm`, or loads cached results from `annotations/dbCAN_hits_df.tsv`
5. **Tool 3 (Astra/HMMER):** Subsets proteins to union of step 1+2 hits, runs `astra search --hmm_in dbCAN-sub.hmm` against the subset
6. **Consensus:** Keeps proteins called by >= 2 tools; records which tools agreed
7. **Loads annotations** into DuckDB as `source='cazy'` (clears any existing CAZy annotations first)
8. **Updates predicates** (adds `cazy:GH13`, `carbohydrate_active`, etc.)
9. **Saves results** to `cazyme_classification.tsv`

### Caching Behavior

The script caches intermediate results to avoid re-running expensive searches:

| Cache File | Contents | Used By |
|------------|----------|---------|
| `cazyme_classification.tsv` | DIAMOND hits (protein_id, family_class, evalue, bitscore) | Step 1 |
| `annotations/dbCAN_hits_df.tsv` | Astra dbCAN HMM hits | Step 2 |

Step 3 (dbCAN-sub) is NOT cached because it runs only on the intersection subset, which is fast relative to steps 1-2.

### Integration with Stage 07

The `07_build_knowledge_base.py` ingest pipeline calls `classify_cazymes.py` automatically as part of the build process:

1. During `_load_annotations()`, raw dbCAN HMM TSV files are **skipped** (the line `if "dbcan" in tsv.as_posix().lower(): continue`)
2. During `_classify_cazymes()`, the full 3-tool consensus script is invoked
3. If stage04 Astra results include a `dbCAN_hits_df.tsv`, it is copied to `annotations/` for caching

This ensures that only consensus-filtered results enter the database as `source='cazy'`.

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

### Astra search timeout

For very large datasets (>1M proteins), Astra HMM searches can take many hours. The script allows 12 hours per Astra search and 2 hours for DIAMOND. If these limits are exceeded, the affected tool returns empty results and is excluded from consensus.

### Existing annotations

Running the script clears all existing `source='cazy'` annotations before loading new ones. This is safe to re-run.

---

## Omnitrophota Dataset Status

### Current State (as of 2026-02-12)

| Item | Status | Details |
|------|--------|---------|
| `cazyme_classification.tsv` | **Present (DIAMOND-only)** | 591,215 rows, old format (no `n_tools`/`tools` columns) |
| `annotations/dbCAN_hits_df.tsv` | **Missing** | No cached Astra dbCAN HMM results |
| DuckDB `source='cazy'` | **Unknown** | May contain DIAMOND-only results or be empty |
| Protein file | **Present** | `source/proteins_combined.faa` (996 MB, ~2.9M proteins) |

### What Happened

The omni_production dataset was annotated with DIAMOND-only CAZyme classification (the `cazyme_classification.tsv` file has columns `protein_id, family_class, cazy_class, evalue, bitscore` -- no consensus columns). This predates the 3-tool consensus script.

### Steps to Run Full 3-Tool Pipeline

1. **Prepare dbCAN-sub.hmm** (currently NOT hmmpress'd):
   ```bash
   cd data/dbcan_db/
   hmmpress dbCAN-sub.hmm
   # This will take several minutes and produce ~5-7 GB of index files
   ```

2. **Run Astra dbCAN HMM search** (step 2 -- the longest step):
   ```bash
   # Option A: Let classify_cazymes.py handle it (recommended)
   # Option B: Pre-run Astra and cache the output
   mkdir -p data/omni_production/annotations/
   astra search --installed_hmms dbCAN \
       --prot_in data/omni_production/source/astra_input/ \
       --outdir data/omni_production/annotations/ \
       --threads 12
   # This produces annotations/dbCAN_hits_df.tsv which the script will cache-load
   ```

3. **Run the full consensus pipeline:**
   ```bash
   python scripts/classify_cazymes.py \
       --db data/omni_production/sharur.duckdb \
       --threads 12
   ```

   The script will:
   - Load cached DIAMOND results from existing `cazyme_classification.tsv`
   - Load cached Astra dbCAN HMM results from `annotations/dbCAN_hits_df.tsv` (if pre-run), or run Astra
   - Subset proteins and run dbCAN-sub HMMER
   - Apply consensus (>= 2/3 tools)
   - Overwrite `cazyme_classification.tsv` with consensus results
   - Load consensus annotations into DuckDB as `source='cazy'`
   - Update predicates

### Runtime Estimates (2.9M proteins)

| Step | Estimated Time | Notes |
|------|---------------|-------|
| hmmpress dbCAN-sub.hmm | 5-15 min | One-time setup |
| DIAMOND (step 1) | **Cached** | Already completed |
| Astra dbCAN HMM (step 2) | 2-8 hours | ~800 HMM profiles vs 2.9M proteins |
| Astra dbCAN-sub (step 3) | 1-4 hours | Subset only; ~25k profiles but fewer proteins |
| Consensus + DB load | < 5 min | Fast |
| **Total** | **3-12 hours** | Depends on thread count and system load |

### Important: DIAMOND Cache Compatibility

The existing `cazyme_classification.tsv` uses column name `bitscore` (old format). The consensus script's `load_diamond_cache()` function reads from this column via `row.get('bitscore', 0.0)`, so it will load correctly. After running, the file will be overwritten with consensus results using the new column names (`score`, `n_tools`, `tools`).

---

## Comparison: Hinthialibacterota (Reference Run)

The Hinthialibacterota dataset (41 genomes, 184K proteins) was successfully run with the full 3-tool pipeline:

- **Input:** 184,689 Astra dbCAN HMM hits (pre-filtered)
- **Output:** 20,945 consensus CAZyme annotations
- **Format:** Full consensus columns (`n_tools`, `tools` present)
- **All three tools** contributed: entries show `dbcan_hmm,dbcan_sub,diamond` combinations
- **Runtime:** Minutes (small dataset)

This confirms the pipeline works end-to-end. The key difference for omni_production is scale (16x more proteins).
