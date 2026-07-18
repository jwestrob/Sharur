# Sharur Ingestion Pipeline

## Quick Start: Primary CLI

Use `sharur-ingest` as the default way to run the Sharur ingest pipeline. Do not add manual Astra database overrides unless you are intentionally doing a specialized run.

```bash
sharur-ingest \
  --input-dir /path/to/genome_fastas \
  --data-dir data/my_dataset \
  --output data/my_dataset/sharur.duckdb \
  --profile auto
```

If `sharur-ingest` is missing, refresh the editable install with
`pip install -e ".[embeddings]"`.

The default plan runs the core stages and skips optional QUAST, DFAST, GECCO, and the
deprecated legacy dbCAN helper. Opt in with `--with-quast`, `--with-dfast`,
`--with-gecco`, or `--with-legacy-dbcan`. The separate Stage 07 dbCAN three-tool
consensus classifier is opt-in with `--enable-cazymes`. `--dry-run` prints the exact
stage order and paths without creating the dataset directory.

The primary CLI builds a dependency DAG and writes its run, stage attempts, signatures,
heartbeats, resources, commands, and output snapshots to `data/my_dataset/sharur_ops.db`.
`--resume` is the default. Reuse requires an identical stage signature plus live, non-empty
declared outputs that match the successful attempt's recorded snapshots; missing or modified
outputs invalidate reuse. Signatures cover concrete commands, stage script content/stat,
resource requests, dependency signatures, and artifact snapshots. Directory snapshots use a
recursive metadata fingerprint, while Stage 00's checksum manifest remains the authoritative
content inventory. If an upstream stage executes rather than being reused, every downstream
dependent stage executes too.

Resource profiles:

| Flag | Execution |
|---|---|
| `--profile auto` | MPS when PyTorch reports a usable backend; bounded local CPU otherwise |
| `--profile local` | All stages local; Stage 06 CPU |
| `--profile mps` | Stage 06 MPS under a cross-process exclusive lock; other stages local |
| `--profile slurm` | Generate `slurm/*.sbatch`, local wrappers, logs, and dependency-linked `submit.sh` |

The SLURM profile gives Astra a 72-hour allocation and Stage 06 one CUDA GPU. MinCED remains
single-threaded and local. Bundle generation does not submit jobs unless
`--submit-slurm` is supplied.

## Manual Stage-by-Stage Pipeline

Use the individual stage scripts below when you need direct stage control, debugging, or partial reruns. Do not modify flags casually. Do not add `--databases` unless you explicitly want a non-default database set. Do not skip Stage 05c.

```bash
DATASET=my_dataset
INPUT=/path/to/genome/fastas

# Create the standard Sharur directory structure
mkdir -p data/$DATASET/{source,annotations,embeddings,structures,exploration,figures,reports,survey}

# Stage 00: Validate and organize input assemblies
python src/ingest/00_prepare_inputs.py -i $INPUT -o data/$DATASET/stage00_prepared --copy

# Stage 03: Gene calling (default --mode single is correct for most inputs)
python src/ingest/03_prodigal.py -i data/$DATASET/stage00_prepared -o data/$DATASET/stage03_prodigal --max-workers 8

# Stage 04: Functional annotation (defaults: PFAM, KOFAM, HydDB, DefenseFinder, dbCAN)
python src/ingest/04_astra_scan.py -i data/$DATASET/stage03_prodigal -o data/$DATASET/stage04_astra -t 12

# Stage 05c: CRISPR array detection (runs on nucleotide assemblies, NOT proteins)
python src/ingest/minced_crispr.py -i data/$DATASET/stage00_prepared -o data/$DATASET/stage05c_crispr

# Stage 07: Build knowledge base (auto-discovers all stage subdirectories)
python src/ingest/07_build_knowledge_base.py -d data/$DATASET -o data/$DATASET/sharur.duckdb
```

That is the manual stage sequence behind `sharur-ingest`. Stage 07 auto-discovers all `stage*` directories under `data/$DATASET/` and loads proteins, annotations, CRISPR arrays, and predicates.

```bash
# Stage 06: ESM2 embeddings (standard — required for ELSA synteny discovery)
# GPU recommended (MPS/CUDA); CPU works but 10-20x slower
python src/ingest/06_esm2_embeddings.py data/$DATASET/stage03_prodigal data/$DATASET/embeddings/
```

Stage 06 streams the Stage 03 FASTAs into the canonical `protein_embeddings.h5` (HDF5:
`protein_ids` string array + `embeddings` float32 matrix), writes an embedding manifest, and
then builds a persistent FAISS index, stable SQLite row-to-protein map, and atomic
active-generation manifest in a Torch-free process. The H5 remains the ELSA input and
canonical embedding artifact.

**ELSA synteny + genome browser (after Stage 06 + 07):**
```bash
# Run ELSA with Sharur annotations (skips redundant Astra PFAM scan)
elsa synteny \
    --db data/$DATASET/sharur.duckdb \
    --embeddings data/$DATASET/embeddings/protein_embeddings.h5 \
    --annotations-db data/$DATASET/sharur.duckdb \
    --store data/$DATASET/synteny/store \
    -o data/$DATASET/synteny/

# Populate genome browser DB with all annotation sources from Sharur
elsa browser data/$DATASET/synteny/ --store data/$DATASET/synteny/store \
    --annotations-db data/$DATASET/sharur.duckdb
```

The `--annotations-db` flag loads PFAM domains into the browser `genes` table and all annotation sources (KEGG, CAZy, DefenseFinder, etc.) into the browser `annotations_multi` table, eliminating the need for a separate Astra PFAM run.

---

> **WARNING: DO NOT pass `--databases` to Stage 04.**
> The defaults (`PFAM, KOFAM, HydDB, DefenseFinder, dbCAN`) are correct and cover all
> standard annotation sources. Passing explicit `--databases` flags causes errors because
> Typer `List[str]` options require repeated flags (`-d PFAM -d KOFAM -d HydDB`), not
> space-separated values. If you think you need to override databases, you almost certainly
> do not. See the Stage 04 section below for details.

> **WARNING: DO NOT use CRISPRCasFinder.**
> CRISPRCasFinder is installed in Astra as an HMM database but it is NOT a standard pipeline
> tool. It detects Cas protein domains (which DefenseFinder already handles better). For
> CRISPR *array* detection (the repeat-spacer structures), use MinCED (Stage 05c).

> **WARNING: Stage 05c (MinCED) IS a standard pipeline step.**
> It detects CRISPR repeat-spacer arrays in the nucleotide assemblies. Stage 07 loads these
> arrays into the `loci` table. Skipping this stage means CRISPR array coordinates are
> missing from the knowledge base.

---

## Pipeline Architecture

```
Stage 00  Input Preparation       (validate + organize genome FASTAs)
   |
Stage 03  Gene Calling            (Prodigal: genomes -> proteins + GFF)
   |
   +---> Stage 04  Annotation    (Astra/PyHMMer: proteins -> HMM hits)
   |
   +---> Stage 05c CRISPR Arrays (MinCED: genomes -> CRISPR array coords)
   |
Stage 07  Knowledge Base          (consolidate everything -> sharur.duckdb)
   |
Stage 06  Embeddings              (ESM2: proteins -> vector embeddings for ELSA)
   |
Stage 06i Persistent index        (CPU FAISS + stable protein-ID map)
```

**Primary CLI path:** `sharur-ingest`
**Standard pipeline stages:** 00, 03, 04, 05c, 07, 06, 06i
**Post-pipeline (standard):** 06 (embeddings) then 06i (persistent similarity index)
**Optional/deprecated:** 01 (QUAST QC), 02 (DFAST QC), 05a (GECCO BGC), 05b (dbCAN legacy)

---

## Stage Numbering

| Stage | Script | Status | Purpose |
|-------|--------|--------|---------|
| 00 | `00_prepare_inputs.py` | **Standard** | Validate and organize input genomes |
| 01 | `01_run_quast.py` | Optional | Assembly QC metrics (N50, etc.) |
| 02 | `02_dfast_qc.py` | Optional | DFAST-based quality control |
| 03 | `03_prodigal.py` | **Standard** | Protein-coding gene prediction |
| 04 | `04_astra_scan.py` | **Standard** | HMM-based functional annotation |
| 05a | `gecco_bgc.py` | Optional | Biosynthetic gene cluster detection |
| 05b | `dbcan_cazyme.py` | Deprecated | Legacy dbCAN (replaced by Stage 07 built-in) |
| 05c | `minced_crispr.py` | **Standard** | CRISPR repeat-spacer array detection |
| 06 | `06_esm2_embeddings.py` | **Standard** | Protein vector embeddings (ELSA input) |
| 06i | `vector_index_runner.py` | **Internal** | Persistent FAISS/ID sidecars in a Torch-free process |
| 07 | `07_build_knowledge_base.py` | **Standard** | Build DuckDB knowledge base |

Stages 01, 02, 05a, and 05b exist as scripts but are not required for a standard pipeline run. Stage 07 gracefully handles their absence.

---

## Stage Reference

### Stage 00: Prepare Inputs (`00_prepare_inputs.py`)

**What it does:** Validates FASTA format, checks for duplicate sequence IDs, generates checksums, and creates an organized output directory with symlinks (or copies) of input genomes plus a `processing_manifest.json`.

**Required inputs:** Directory of genome assembly FASTA files (`.fna`, `.fa`, `.fasta`).

**Outputs:**
- `stage00_prepared/` directory containing symlinked (or copied) genome files
- `stage00_prepared/processing_manifest.json` (consumed by Stage 03)

**Minimal invocation:**
```bash
python src/ingest/00_prepare_inputs.py -i /path/to/genomes -o data/DATASET/stage00_prepared
```

**Important options:**

| Flag | Default | Description |
|------|---------|-------------|
| `-i`, `--input-dir` | `data/raw` | Directory containing input genome FASTAs |
| `-o`, `--output-dir` | `data/stage00_prepared` | Output directory |
| `--copy` / `--symlink` | `--symlink` | Copy files instead of symlinking (use `--copy` if genomes are on a different filesystem) |
| `--force` | off | Overwrite existing output directory |
| `-e`, `--extensions` | `.fasta .fa .fna` | File extensions to search for |

**Common errors:**
- "Output directory already exists" -- add `--force` or delete the existing directory.
- Symlink failures -- use `--copy` if input files are on a different mount or filesystem.

---

### Stage 03: Gene Calling (`03_prodigal.py`)

**What it does:** Runs Prodigal gene prediction on each genome, producing protein FASTA (`.faa`) and nucleotide gene sequences (`.genes.fna`). Processes genomes in parallel. Creates a `genomes/all_protein_symlinks/` directory containing symlinks to all `.faa` files (this is what Stage 04 reads).

**Required inputs:** Stage 00 output directory (must contain `processing_manifest.json`).

**Outputs:**
- `stage03_prodigal/genomes/{genome_id}/{genome_id}.faa` -- protein sequences per genome
- `stage03_prodigal/genomes/{genome_id}/{genome_id}.genes.fna` -- nucleotide gene sequences
- `stage03_prodigal/genomes/all_protein_symlinks/` -- all `.faa` files symlinked here
- `stage03_prodigal/processing_manifest.json`

**Minimal invocation:**
```bash
python src/ingest/03_prodigal.py -i data/DATASET/stage00_prepared -o data/DATASET/stage03_prodigal --max-workers 8
```

**Important options:**

| Flag | Default | Description |
|------|---------|-------------|
| `-i`, `--input-dir` | `data/stage00_prepared` | Stage 00 output directory |
| `-o`, `--output-dir` | `data/stage03_prodigal` | Output directory |
| `-m`, `--mode` | `single` | Prodigal mode: `single` trains on each genome (better gene calls). Only use `meta` for very short contigs (<100 kb total) where Prodigal can't train. |
| `-w`, `--max-workers` | CPU count - 1 | Parallel workers |
| `-g`, `--genetic-code` | `11` | Genetic code table (11 = bacteria/archaea) |
| `--force` | off | Overwrite existing output |

**Common errors:**
- Prodigal segfault on short contigs or malformed input -- use Pyrodigal as a fallback (`pip install pyrodigal`).
- Timeout (>5 min per genome) -- check for extremely large genome files.
- "Input manifest not found" -- run Stage 00 first.

---

### Stage 04: Functional Annotation (`04_astra_scan.py`)

**What it does:** Runs Astra HMM searches against multiple domain databases. Processes databases sequentially (one at a time) to avoid resource contention. Each database produces a `{DATABASE}_hits_df.tsv` results file.

**Required inputs:** Stage 03 output directory (specifically the `genomes/all_protein_symlinks/` subdirectory).

**Outputs:**
- `stage04_astra/{database}_results/{DATABASE}_hits_df.tsv` -- hit tables per database
- `stage04_astra/processing_manifest.json`
- `stage04_astra/summary_stats.json`

**Minimal invocation:**
```bash
python src/ingest/04_astra_scan.py -i data/DATASET/stage03_prodigal -o data/DATASET/stage04_astra -t 12
```

That is all you need. The five default databases cover all standard annotation:

| Database | What it annotates | Cutoff strategy |
|----------|-------------------|-----------------|
| PFAM | Protein domain families | `--cut_ga` (gathering thresholds) |
| KOFAM | KEGG orthologs | `--cut_ga --cascade` (per-profile adaptive) |
| HydDB | Hydrogenase subgroups | `--cut_ga` |
| DefenseFinder | Defense systems (CRISPR Cas proteins, RM, CBASS, etc.) | `--cut_ga` (permissive; Stage 07 applies e-value 1e-15 filter) |
| dbCAN | CAZyme domain HMMs | No GA; Stage 07 applies e-value filter |

**Non-default databases (opt-in via `-d`):**

| Database | What it annotates | Cutoff strategy |
|----------|-------------------|-----------------|
| TXSScan | Secretion systems (T1SS–T9SS, T4P, Tad, Flagellum, MSH, ComM) | No GA; Stage 07 applies e-value 1e-10 filter |
| VOGdb | Viral ortholog groups | No GA; requires e-value 1e-15 filter |
| CANT-HYD | Hydrocarbon degradation markers | No GA; e-value 1e-10 filter |

DefenseFinder and TXSScan automatically produce `--write_macsyfinder` output for downstream MacSyFinder co-localization validation in Stage 07.

**Important options:**

| Flag | Default | Description |
|------|---------|-------------|
| `-i`, `--input-dir` | `data/stage03_prodigal` | Stage 03 output directory |
| `-o`, `--output-dir` | `data/stage04_astra` | Output directory |
| `-t`, `--threads` | `8` | Threads per database scan |
| `--force` | off | Overwrite existing output |

> **DO NOT pass `--databases` / `-d`.**
> The default list `["PFAM", "KOFAM", "HydDB", "DefenseFinder", "dbCAN"]` is correct.
> If you must override (rare), Typer requires **repeated flags**: `-d PFAM -d KOFAM -d HydDB`.
> Passing `-d "PFAM KOFAM HydDB"` or `--databases PFAM KOFAM HydDB` will fail or produce
> wrong results. But again: do not override the defaults.

> **DO NOT add VOGdb, CRISPRCasFinder, or CANT-HYD to the default database list.**
> - VOGdb has no GA thresholds (0/48,439 profiles) and produces noisy results that require careful filtering.
> - CRISPRCasFinder in Astra detects Cas protein domains, not CRISPR arrays. DefenseFinder is better for Cas proteins, and MinCED (Stage 05c) handles array detection.
> - CANT-HYD is specialized and only needed for specific hydrocarbon degradation analysis.
> These databases can be added for specialized analyses, but they are not standard pipeline databases.

> **TXSScan (secretion systems) is available but NOT a default database.**
> Add `-d TXSScan` alongside the defaults to enable secretion system detection (T1SS–T9SS,
> T4P, Tad, Flagellum, MSH, ComM). Requires repeated `-d` flags for all databases:
> ```bash
> python src/ingest/04_astra_scan.py -i ... -o ... -t 12 \
>     -d PFAM -d KOFAM -d HydDB -d DefenseFinder -d dbCAN -d TXSScan
> ```
> Stage 07 automatically validates TXSScan hits via MacSyFinder co-localization when
> `txsscan_results/macsyfinder_compat/` exists.

**Common errors:**
- "Protein symlink directory not found" -- run Stage 03 first.
- Astra not found -- ensure `astra` is on PATH (`which astra`). Installed via pyenv shim at `~/astra/`.
- KOFAM takes hours on large datasets (>100k proteins) -- this is normal, do not kill the process.
- HMM database not found -- check `~/.config/Astra/` for installed databases.

---

### Stage 05c: CRISPR Array Detection (`minced_crispr.py`)

**What it does:** Runs MinCED on each genome FASTA to detect CRISPR repeat-spacer arrays. Produces per-genome JSON files that Stage 07 loads into the `loci` table with `locus_type='crispr'`.

**This is a standard pipeline step.** CRISPR array detection is separate from Cas protein annotation (which DefenseFinder handles in Stage 04). MinCED finds the DNA repeat-spacer structures; DefenseFinder finds the protein machinery.

**Required inputs:** Stage 00 output directory (nucleotide genome FASTAs, not protein files).

**Outputs:**
- `stage05c_crispr/{genome_stem}_crispr.gff` -- GFF annotation per genome
- `stage05c_crispr/{genome_stem}_crispr.txt` -- MinCED text output
- `stage05c_crispr/{genome_stem}_crispr_arrays.json` -- parsed arrays (consumed by Stage 07)

**Minimal invocation:**
```bash
python src/ingest/minced_crispr.py -i data/DATASET/stage00_prepared -o data/DATASET/stage05c_crispr
```

**Important options:**

| Flag | Default | Description |
|------|---------|-------------|
| `-i`, `--input-dir` | (required) | Directory containing genome FASTAs (`.fna`, `.fa`, `.fasta`) |
| `-o`, `--output-dir` | (required) | Output directory |
| `--force` | off | Overwrite existing output |

**Common errors:**
- "minced executable not found on PATH" -- install MinCED or ensure the wrapper script is on PATH. Known location: `~/.pyenv/versions/miniconda3-latest/bin/minced`.
- Timeout (30 min per genome) -- check for extremely large genome files.
- No FASTA files found -- ensure input directory contains `.fna`, `.fa`, or `.fasta` files (not `.faa` protein files).

---

### Stage 07: Build Knowledge Base (`07_build_knowledge_base.py`)

**What it does:** Consolidates all upstream outputs into a single DuckDB database. Auto-discovers stage directories under the data directory. Loads proteins, annotations, CRISPR arrays, BGC loci (if present), classifies hydrogenases (if HydDB annotations exist), optionally runs dbCAN three-tool consensus with `--enable-cazymes`, validates defense/secretion systems via the co-location engine, then generates the V2 predicate tables and V2-derived legacy compatibility predicates.

**Required inputs:** Data directory containing at minimum `stage03_prodigal/` and `stage04_astra/`. Loads additional data from any other stage directories that exist.

**Stage directory auto-discovery:** Stage 07 looks for these subdirectories under `--data-dir`:

| Subdirectory | Source stage | What it loads |
|---|---|---|
| `stage00_prepared/` | Stage 00 | (not directly loaded; consumed by Stage 03) |
| `stage02_dfast_qc/` | Stage 02 | Bin metadata (completeness, contamination, taxonomy) |
| `stage03_prodigal/` | Stage 03 | Proteins, contigs, bin records |
| `stage04_astra/` | Stage 04 | All `*_hits_df.tsv` annotation files (except dbCAN raw, handled by CAZyme pipeline) |
| `stage05a_gecco/` | GECCO | BGC loci (if present) |
| `stage05b_dbcan/` | Legacy dbCAN | Legacy CAZyme JSON (backward compat) |
| `stage05c_crispr/` | Stage 05c | CRISPR array loci from `*_crispr_arrays.json` |
| `embeddings/` | Stage 06 | Canonical embeddings and count metadata |

Stage 07 also recognizes the legacy `stage06_embeddings/` layout when reading older datasets.

**Outputs:**
- `sharur.duckdb` -- complete knowledge base with tables: `bins`, `contigs`, `proteins`, `annotations`, `loci`, `locus_proteins`, `semantic_atoms`, `semantic_state`, `protein_predicates`, `feature_store`
- `reports/predicates_v2_review_queue.tsv` -- unresolved annotation accessions ranked for curation

**Minimal invocation:**
```bash
python src/ingest/07_build_knowledge_base.py -d data/DATASET -o data/DATASET/sharur.duckdb
```

**Important options:**

| Flag | Default | Description |
|------|---------|-------------|
| `-d`, `--data-dir` | `data` | Parent directory containing `stage*` subdirectories |
| `-o`, `--output` | `data/sharur.duckdb` | Output DuckDB path |
| `--force` | off | Overwrite existing database |
| `--enable-cazymes` | off | Run the slower dbCAN three-tool consensus classifier |

**What Stage 07 does automatically (no user intervention needed):**
1. Loads bins from Stage 02 manifest (or infers bins from Stage 03 FAA files)
2. Parses Prodigal FAA headers for protein coordinates, contigs, gene indices
3. Loads all annotation TSV files from Stage 04, applying e-value thresholds per source
4. Loads CRISPR arrays from Stage 05c JSON files into `loci` table
5. Loads BGC loci from Stage 05a (if present)
6. Computes length z-scores per bin
7. Runs HydDB subgroup classification (if HydDB annotations exist)
8. Optionally runs dbCAN three-tool consensus CAZyme classification when `--enable-cazymes` is set
9. Validates DefenseFinder hits via co-location rules → `defense_systems` table + `defensefinder_system` annotations
10. Validates TXSScan hits via co-location rules → `secretion_systems` table + `txsscan_system` annotations
11. Generates V2 semantic atoms/states for all proteins
12. Materializes `protein_predicates` from V2 for legacy query compatibility
13. Emits CRISPR-overlap quality flags in V2 and the compatibility table
14. Creates indexes

**E-value thresholds applied at load time:**

| Source | Threshold | Notes |
|--------|-----------|-------|
| pfam | 1e-5 | Safety net; already clean via `--cut_ga` |
| kegg | 1e-5 | Safety net; already clean via `--cut_ga --cascade` |
| hyddb | 1e-5 | Safety net; already clean via `--cut_ga` |
| defensefinder | 1e-15 | Strict; superfamily HMMs are permissive at GA threshold |
| txsscan | 1e-10 | No GA thresholds; secretion system HMMs |
| vogdb | 1e-15 | No GA thresholds; e-value filter required |
| cant_hyd | 1e-10 | No `--cut_ga` |
| cazy | 1e-5 | dbCAN has own thresholds |

**Common errors:**
- "FileExistsError: sharur.duckdb exists" -- add `--force` to overwrite.
- Missing `classify_hydrogenases` or `classify_cazymes` -- these are non-fatal; classification steps are skipped with a warning.
- Missing reference files (`data/reference/pfam_id_desc.tsv`, `data/reference/ko_list`) -- annotations will lack human-readable names but pipeline still works.

---

## Optional Stages

### Stage 01: QUAST Quality Assessment (`01_run_quast.py`)

Assembly QC metrics (N50, L50, contig counts). Not required for the standard pipeline but useful for assessing assembly quality before analysis.

### Stage 02: DFAST QC (`02_dfast_qc.py`)

Quality control using DFAST. Produces bin metadata (completeness, contamination, taxonomy) that Stage 07 can load.

### Stage 05a: BGC Detection (`gecco_bgc.py`)

Biosynthetic gene cluster detection using GECCO. Stage 07 loads results from `stage05a_gecco/combined_bgc_data.json` if present.

### Stage 05b: Legacy dbCAN (`dbcan_cazyme.py`)

**Deprecated.** For new datasets, use Stage 04 dbCAN annotations or opt into Stage 07's
three-tool consensus with `--enable-cazymes`; do not run this helper separately.

### Stage 06: ESM2 Embeddings (`06_esm2_embeddings.py`)

Generates 320-dimensional protein embeddings using ESM2. Requires PyTorch and Transformers.
GPU recommended. It streams `.faa` records and extendable H5 batches, so neither the
proteome nor the complete embedding matrix is accumulated in RAM. Pooling excludes padding
and special tokens; invalid or zero model output fails the stage. The manifest reports the
number of proteins truncated at the model residue limit. This stage can run any time after
Stage 03 and does not need the DuckDB database.

```bash
python src/ingest/06_esm2_embeddings.py \
  data/DATASET/stage03_prodigal \
  data/DATASET/embeddings/ \
  --device mps
```

Use `--device cpu|mps|cuda` to make placement explicit, `--skip-index` to defer sidecar
construction, and `--index-threads N` to bound the CPU FAISS build. The standard DAG uses
`--skip-index` and records the Torch-free CPU build as a separate `06i` attempt, so index
retries reuse the completed H5. Direct Stage 06 use launches the same builder in an isolated
child process. For an existing H5:

```bash
sharur build-vector-index --embeddings data/DATASET/embeddings/protein_embeddings.h5
```

---

## Output Database Schema

After a successful pipeline run, `sharur.duckdb` contains:

| Table | Key Columns | Description |
|-------|-------------|-------------|
| `bins` | `bin_id`, `completeness`, `contamination`, `taxonomy` | Genome/MAG metadata |
| `contigs` | `contig_id`, `bin_id`, `length`, `gc_content` | Contig records |
| `proteins` | `protein_id`, `contig_id`, `bin_id`, `sequence`, `sequence_length`, `gene_index` | All predicted proteins |
| `annotations` | `protein_id`, `source`, `accession`, `name`, `description`, `evalue`, `score` | HMM annotation hits |
| `loci` | `locus_id`, `locus_type`, `contig_id`, `start`, `end_coord` | CRISPR arrays, BGCs |
| `locus_proteins` | `locus_id`, `protein_id`, `position` | Protein membership in loci |
| `semantic_atoms` | `protein_id`, `atom_id`, `facet`, `relation`, `source_db`, `source_accession` | V2 evidence-backed semantic claims |
| `semantic_state` | `protein_id`, `activities`, `roles`, `architecture`, `localization`, `quality_flags`, `composite_predicates` | V2 resolved per-protein state |
| `protein_predicates` | `protein_id`, `predicates` (list) | V2-derived legacy compatibility predicate tags |
| `defense_systems` | `system_id`, `genome_id`, `system_type`, `protein_ids` | MacSyFinder-validated defense systems (if DefenseFinder ran) |
| `secretion_systems` | `system_id`, `genome_id`, `system_type`, `protein_ids` | MacSyFinder-validated secretion systems (if TXSScan ran) |
| `feature_store` | `protein_id`, `metric_name`, `metric_value` | Computed metrics (length z-score) |

**Query conventions (for downstream agents):**
- Use `name` column for domain names, NOT `annotation_id`
- Use `score` column for bitscores, NOT `bitscore`
- Always `COUNT(DISTINCT protein_id)` for protein counts -- repeat domains inflate `COUNT(*)`

---

## Troubleshooting

### Prodigal segfault
Prodigal can segfault on very short contigs or malformed FASTA. Install Pyrodigal as a fallback:
```bash
pip install pyrodigal
```

### MinCED not found
MinCED is a Java-based tool. Known wrapper script location:
```bash
~/.pyenv/versions/miniconda3-latest/bin/minced
```
Ensure this is on your PATH, or install via conda: `conda install -c bioconda minced`.

### Astra not found
Astra is installed at `~/astra/` and accessed via pyenv shim. Verify: `which astra`. Installed HMM databases live at `~/.config/Astra/`.

### Typer List[str] options
If you must override a `List[str]` option (like `--databases`), Typer requires **repeated flags**:
```bash
# CORRECT (but you should not need to do this)
python src/ingest/04_astra_scan.py -d PFAM -d KOFAM -d HydDB ...

# WRONG -- will fail or produce unexpected results
python src/ingest/04_astra_scan.py --databases "PFAM KOFAM HydDB" ...
python src/ingest/04_astra_scan.py --databases PFAM KOFAM HydDB ...
```

### KOFAM is slow
KOFAM on large datasets (>100k proteins) can take several hours. This is expected behavior due to the large number of KO profiles. Do not kill the process.

### "ModuleNotFoundError: No module named 'sharur'"
Install Sharur in development mode from the repo root:
```bash
pip install -e .
```

### DuckDB version mismatch
```bash
pip install --upgrade duckdb
```

### "CUDA out of memory" (Stage 06)
Reduce batch size in the embeddings script, or use CPU (10-20x slower but works).

---

## Dependencies

### Core Pipeline (Stages 00, 03, 04, 05c, 07)
- Python 3.8+
- DuckDB, Pandas, Biopython
- Typer, Rich (CLI framework)
- [Astra](https://github.com/Dreycey/Astra) (HMM search)
- Prodigal (gene calling)
- MinCED (CRISPR detection)

### Embeddings (Stage 06)
- PyTorch, Transformers (HuggingFace)
- GPU recommended

### Astra Installed Databases (`~/.config/Astra/`)
PFAM, KOFAM, VOGdb, HydDB, DefenseFinder, CRISPRCasFinder, CANT-HYD, dbCAN, padloc, TXSScan

Note: Not all installed databases are used in the standard pipeline. Only PFAM, KOFAM, HydDB, DefenseFinder, and dbCAN are defaults. TXSScan is available as a non-default opt-in.
