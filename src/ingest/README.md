# Sharur Ingestion Pipeline

## Quick Start: End-to-End Pipeline

Copy-paste this. Do not modify flags. Do not add `--databases`. Do not skip Stage 05c.

```bash
DATASET=my_dataset
INPUT=/path/to/genome/fastas

# Create the standard Sharur directory structure
mkdir -p data/$DATASET/{source,annotations,embeddings,structures,exploration,figures,reports,survey}

# Stage 00: Validate and organize input assemblies
python src/ingest/00_prepare_inputs.py -i $INPUT -o data/$DATASET/stage00_prepared --copy

# Stage 03: Gene calling (--mode meta for MAGs, --mode single for isolates)
python src/ingest/03_prodigal.py -i data/$DATASET/stage00_prepared -o data/$DATASET/stage03_prodigal --mode meta --max-workers 8

# Stage 04: Functional annotation (defaults: PFAM, KOFAM, HydDB, DefenseFinder, dbCAN)
python src/ingest/04_astra_scan.py -i data/$DATASET/stage03_prodigal -o data/$DATASET/stage04_astra -t 12

# Stage 05c: CRISPR array detection (runs on nucleotide assemblies, NOT proteins)
python src/ingest/minced_crispr.py -i data/$DATASET/stage00_prepared -o data/$DATASET/stage05c_crispr

# Stage 07: Build knowledge base (auto-discovers all stage subdirectories)
python src/ingest/07_build_knowledge_base.py -d data/$DATASET -o data/$DATASET/sharur.duckdb
```

That is the complete standard pipeline. Stage 07 auto-discovers all `stage*` directories under `data/$DATASET/` and loads proteins, annotations, CRISPR arrays, and predicates.

Optional post-pipeline step:
```bash
# Stage 06: ESM2 embeddings (GPU recommended; CPU works but 10-20x slower)
python src/ingest/06_esm2_embeddings.py data/$DATASET/stage03_prodigal data/$DATASET/embeddings/
```

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
Stage 06  Embeddings (optional)   (ESM2: proteins -> vector embeddings)
```

**Standard pipeline stages:** 00, 03, 04, 05c, 07
**Optional/post-pipeline:** 06 (embeddings)
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
| 06 | `06_esm2_embeddings.py` | Optional | Protein vector embeddings |
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
python src/ingest/03_prodigal.py -i data/DATASET/stage00_prepared -o data/DATASET/stage03_prodigal --mode meta
```

**Important options:**

| Flag | Default | Description |
|------|---------|-------------|
| `-i`, `--input-dir` | `data/stage00_prepared` | Stage 00 output directory |
| `-o`, `--output-dir` | `data/stage03_prodigal` | Output directory |
| `-m`, `--mode` | `single` | Prodigal mode: `meta` for MAGs, `single` for isolates |
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

**What it does:** Consolidates all upstream outputs into a single DuckDB database. Auto-discovers stage directories under the data directory. Loads proteins, annotations, CRISPR arrays, BGC loci (if present), generates predicates, classifies hydrogenases (if HydDB annotations exist), runs dbCAN 3-tool consensus CAZyme classification, and validates defense/secretion systems via MacSyFinder co-localization (if macsyfinder_compat output exists).

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
| `stage06_embeddings/` | Stage 06 | Embedding count metadata |

**Outputs:**
- `sharur.duckdb` -- complete knowledge base with tables: `bins`, `contigs`, `proteins`, `annotations`, `loci`, `locus_proteins`, `protein_predicates`, `feature_store`

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

**What Stage 07 does automatically (no user intervention needed):**
1. Loads bins from Stage 02 manifest (or infers bins from Stage 03 FAA files)
2. Parses Prodigal FAA headers for protein coordinates, contigs, gene indices
3. Loads all annotation TSV files from Stage 04, applying e-value thresholds per source
4. Loads CRISPR arrays from Stage 05c JSON files into `loci` table
5. Loads BGC loci from Stage 05a (if present)
6. Computes length z-scores per bin
7. Generates semantic predicates for all proteins (parallelized)
8. Flags proteins overlapping CRISPR arrays
9. Runs HydDB subgroup classification (if HydDB annotations exist)
10. Runs dbCAN 3-tool consensus CAZyme classification (DIAMOND + dbCAN.hmm + dbCAN-sub.hmm)
11. Validates DefenseFinder hits via MacSyFinder co-localization (if `defensefinder_results/macsyfinder_compat/` exists) → `defense_systems` table + `defensefinder_system` annotations
12. Validates TXSScan hits via MacSyFinder co-localization (if `txsscan_results/macsyfinder_compat/` exists) → `secretion_systems` table + `txsscan_system` annotations
13. Creates indexes

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

**Deprecated.** Stage 07 now runs dbCAN 3-tool consensus classification internally using `scripts/classify_cazymes.py`. Do not run this separately for new datasets.

### Stage 06: ESM2 Embeddings (`06_esm2_embeddings.py`)

Generates 320-dimensional protein embeddings using ESM2. Requires PyTorch, Transformers, and LanceDB. GPU recommended. Can run any time after Stage 03 (reads .faa files from the Prodigal output directory, does NOT need the DuckDB database).

```bash
python src/ingest/06_esm2_embeddings.py data/DATASET/stage03_prodigal data/DATASET/embeddings/
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
| `protein_predicates` | `protein_id`, `predicates` (list) | Semantic predicate tags |
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
- PyTorch, Transformers (HuggingFace), LanceDB
- GPU recommended

### Astra Installed Databases (`~/.config/Astra/`)
PFAM, KOFAM, VOGdb, HydDB, DefenseFinder, CRISPRCasFinder, CANT-HYD, dbCAN, padloc, TXSScan

Note: Not all installed databases are used in the standard pipeline. Only PFAM, KOFAM, HydDB, DefenseFinder, and dbCAN are defaults. TXSScan is available as a non-default opt-in.
