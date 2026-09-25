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

