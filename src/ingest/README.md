# Sharur Ingestion Pipeline

This directory contains the staged ingestion pipeline for converting raw genomic data into Sharur-ready databases.

---

## Pipeline Overview

The ingestion pipeline is organized into numbered stages (00-07). Each stage is independent and can be run individually or as part of a full pipeline.

```
Raw Data → Stage 00-04 → Annotations → Stage 06-07 → Sharur Database
```

---

## Pipeline Stages

### Stage 00: Prepare Inputs (`00_prepare_inputs.py`)
**Purpose:** Standardize input files and create working directory structure

**Inputs:**
- Genome assemblies (`.fna`, `.fasta`)
- Existing gene calls (optional)
- Metadata files

**Outputs:**
- Standardized directory structure
- Validated input files

**Usage:**
```bash
python src/ingest/00_prepare_inputs.py \
  --input /path/to/raw/data \
  --output data/my_dataset_production
```

---

### Stage 01: Quality Assessment (`01_run_quast.py`)
**Purpose:** Run QUAST for assembly quality metrics

**Inputs:** Genome assemblies

**Outputs:**
- Assembly statistics (N50, L50, contigs, etc.)
- Quality report in `quast_results/`

**Usage:**
```bash
python src/ingest/01_run_quast.py \
  --assemblies data/my_dataset/source/*.fna \
  --output data/my_dataset/quast_results
```

---

### Stage 02: DFAST QC (`02_dfast_qc.py`)
**Purpose:** Quality control using DFAST annotation tool

**Inputs:** Genome assemblies

**Outputs:**
- Gene predictions
- Preliminary annotations
- QC metrics

---

### Stage 03: Gene Calling (`03_prodigal.py`)
**Purpose:** Predict protein-coding genes using Prodigal

**Inputs:** Genome assemblies (`.fna`)

**Outputs:**
- Protein sequences (`.faa`)
- Gene coordinates (`.gff`)
- Nucleotide sequences (`.ffn`)

**Usage:**
```bash
python src/ingest/03_prodigal.py \
  --input data/my_dataset/source/genomes.fna \
  --output data/my_dataset/source/proteins.faa \
  --mode meta  # or 'single' for isolate genomes
```

---

### Stage 04: Astra Annotation Scan (`04_astra_scan.py`)
**Purpose:** Run Astra HMM searches for comprehensive annotation

**HMM Databases:**
- PFAM (protein domains)
- KOFAM (KEGG orthologs)
- HydDB (hydrogenases)
- VOGdb (viral genes)
- DefenseFinder (defense systems)
- CRISPRCasFinder (CRISPR-Cas)

**Usage:**
```bash
python src/ingest/04_astra_scan.py \
  --proteins data/my_dataset/source/proteins.faa \
  --output data/my_dataset/annotations \
  --databases PFAM KOFAM HydDB \
  --threads 12
```

**Note:** This stage is a wrapper around Astra CLI. You can also run Astra directly (see QUICKSTART.md).

---

### Stage 05: ⚠️ RESERVED / DEPRECATED

**Status:** No stage 05 currently exists. The numbering gap is intentional to allow for future insertion of additional preprocessing steps without renumbering.

**Historical note:** Stage 05 was originally planned for CAZyme annotation via dbCAN but was moved to optional module (`dbcan_cazyme.py`).

---

### Stage 06: ESM2 Embeddings (`06_esm2_embeddings.py`)
**Purpose:** Generate protein embeddings for semantic similarity search

**Model:** `facebook/esm2_t6_8M_UR50D` (320-dimensional embeddings)

**Inputs:**
- Sharur database (`sharur.duckdb`)
- Protein sequences from database

**Outputs:**
- LanceDB vector store in `embeddings/`
- Embeddings for all proteins in database

**Usage:**
```bash
python src/ingest/06_esm2_embeddings.py \
  data/my_dataset/sharur.duckdb \
  data/my_dataset/embeddings/
```

**Requirements:**
- PyTorch
- Transformers (HuggingFace)
- LanceDB
- GPU recommended (CPU works but 10-20× slower)

**Time:** ~5-15 minutes for 10k proteins (GPU)

---

### Stage 07: Build Knowledge Base (`07_build_knowledge_base.py`)
**Purpose:** Consolidate all annotations into Sharur database and generate predicates

**Inputs:**
- Protein FASTA
- Annotation TSV files (PFAM, KEGG, HydDB, VOGdb, etc.)

**Outputs:**
- Complete `sharur.duckdb` with:
  - Proteins table
  - Annotations table
  - Predicates table
  - Ontology mappings

**Features:**
- Automatic predicate generation from annotations
- Hydrogenase subgroup classification (if HydDB annotations present)
- Metadata extraction (contig, gene index, size, etc.)
- Quality filtering (e-value thresholds)

**Usage:**
```bash
python src/ingest/07_build_knowledge_base.py \
  --proteins data/my_dataset/source/proteins.faa \
  --pfam data/my_dataset/annotations/pfam_results/PFAM_hits_df.tsv \
  --kegg data/my_dataset/annotations/kofam_results/KOFAM_hits_df.tsv \
  --hyddb data/my_dataset/annotations/hyddb_results/HydDB_hits_df.tsv \
  --output data/my_dataset/sharur.duckdb
```

---

## Optional Modules

These scripts provide specialized annotation but are not part of the core pipeline:

### CAZyme Annotation (`dbcan_cazyme.py`)
**Purpose:** Identify carbohydrate-active enzymes using dbCAN

**Databases:** CAZy families (GH, GT, PL, CE, AA, CBM)

**Usage:**
```bash
python src/ingest/dbcan_cazyme.py \
  --proteins data/my_dataset/source/proteins.faa \
  --output data/my_dataset/annotations/cazy_results
```

---

### BGC Detection (`gecco_bgc.py`)
**Purpose:** Detect biosynthetic gene clusters using GECCO

**Usage:**
```bash
python src/ingest/gecco_bgc.py \
  --proteins data/my_dataset/source/proteins.faa \
  --output data/my_dataset/annotations/bgc_results
```

---

### CRISPR Array Detection (`minced_crispr.py`)
**Purpose:** Identify CRISPR arrays using MinCED

**Inputs:** Genome nucleotide sequences (`.fna`)

**Usage:**
```bash
python src/ingest/minced_crispr.py \
  --genomes data/my_dataset/source/genomes.fna \
  --output data/my_dataset/annotations/crispr_arrays
```

---

## Recommended Workflow

### For New Datasets (Start from scratch):
```bash
# Full pipeline
python src/ingest/00_prepare_inputs.py ...   # Setup
python src/ingest/03_prodigal.py ...         # Gene calling
python src/ingest/04_astra_scan.py ...       # Annotations
python src/ingest/07_build_knowledge_base.py ...  # Database
python src/ingest/06_esm2_embeddings.py ...  # Embeddings
```

### For Existing Protein FASTAs (Quick start):
See `QUICKSTART.md` for streamlined workflow starting from proteins.

---

## Dependencies

### Core Pipeline
- Python 3.8+
- DuckDB
- Pandas
- Biopython
- [Astra](https://github.com/Dreycey/Astra)

### Embeddings (Stage 06)
- PyTorch
- Transformers
- LanceDB

### Optional Tools
- Prodigal (gene calling)
- QUAST (QC)
- dbCAN (CAZymes)
- GECCO (BGCs)
- MinCED (CRISPR arrays)

---

## Output Database Schema

After ingestion, `sharur.duckdb` contains:

### Tables

**`proteins`**
- `protein_id` (primary key)
- `sequence` (amino acid sequence)
- `length` (integer)
- `contig_id` (optional)
- `gene_index` (optional)

**`annotations`**
- `protein_id` (foreign key)
- `source` (PFAM, KEGG, etc.)
- `annotation_id` (PF00001, K00001, etc.)
- `evalue` (float)
- `score` (float)
- `description` (text, optional)

**`predicates`**
- `protein_id` (foreign key)
- `predicate` (e.g., "giant", "transporter", "crispr_associated")
- `source` (derived from annotation or direct)

---

## Troubleshooting

### "ModuleNotFoundError: No module named 'bennu'"
Install Sharur in development mode: `pip install -e .`

### "Astra command not found"
Install Astra and ensure it's in PATH: `which astra`

### "DuckDB version mismatch"
Update DuckDB: `pip install --upgrade duckdb`

### "CUDA out of memory" (Stage 06)
Reduce batch size in embeddings script or use CPU (slower)

---

## Contributing

When adding new pipeline stages:
1. Use next available number (08, 09, etc.) or insert at 05
2. Follow naming convention: `NN_description.py`
3. Include `--help` argparse documentation
4. Update this README
5. Add to `QUICKSTART.md` if core workflow

---

**Last Updated:** 2026-02-06
