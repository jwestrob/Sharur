# Data Organization for Sharur Analyses

This document describes the standard organization for genome analysis projects.

The standard ingest workflow starts from genome assembly FASTAs and should usually be run through `sharur-ingest`. The staged scripts in `src/ingest/` remain the manual/debugging path. Protein-only ingestion is a fallback path for datasets that do not have assemblies.

## Directory Structure

```
data/
├── archives/                           # Completed analyses (read-only)
│   ├── {organism}_{YYYY-MM-DD}/       # Archived with completion date
│   └── ...
├── {organism}_production/              # Active analysis
│   ├── sharur.duckdb                   # Main database
│   ├── sharur_ops.db                  # Run ledger + coordination/task state
│   ├── manifest.json                  # Analysis state and session continuity
│   │
│   ├── source/                        # Input files
│   │   ├── *.fna / *.fa / *.fasta     # Source genome assemblies (standard path)
│   │   └── *.faa                      # Pre-called proteins (protein-only fallback path)
│   │
│   ├── stage00_prepared/              # Validated + organized assemblies
│   ├── stage03_prodigal/              # Gene calls and protein FASTAs
│   ├── stage04_astra/                 # Annotation stage outputs
│   ├── stage05c_crispr/               # CRISPR array stage outputs
│   ├── stage05a_gecco/                # Optional BGC stage outputs
│   ├── stage06_embeddings/            # Legacy embedding location (if present)
│   │
│   ├── annotations/                   # Consolidated/derived annotation exports
│   │
│   ├── embeddings/                    # ESM2 protein embeddings
│   │   ├── embedding_manifest.json
│   │   ├── protein_embeddings.h5      # Canonical ESM2 embeddings
│   │   ├── protein_embeddings.index.json   # Atomic active-generation manifest
│   │   ├── protein_embeddings.*.faiss      # Persistent mmap-ready generation
│   │   └── protein_embeddings.*.ids.sqlite # Stable row ↔ protein-ID map
│   │
│   ├── slurm/                         # Generated ingest bundle/logs (SLURM profile)
│   │
│   ├── structures/                    # ESM3 structure predictions
│   │   ├── *.pdb                      # Predicted PDB files
│   │   └── foldseek_results.tsv       # Structural homology searches
│   │
│   ├── synteny/                       # ELSA synteny results + FAISS store
│   │   ├── store/                     # FAISS synteny vector store + metadata
│   │   └── results/                   # Anchors, blocks, clusters (loaded into DuckDB)
│   │
│   ├── exploration/                   # Exploration phase outputs
│   │   ├── findings.jsonl             # Documented loci
│   │   ├── state.json                 # Exploration progress
│   │   └── figures/                   # Phase-specific figures
│   │
│   ├── survey/                        # Survey phase outputs (if run separately)
│   │   ├── findings.jsonl
│   │   ├── state.json
│   │   └── figures/
│   │
│   ├── reports/                       # All generated reports
│   │   ├── survey.pdf
│   │   ├── exploration.pdf
│   │   ├── final.pdf
│   │   └── *.md                       # Markdown summaries
│   │
│   └── figures/                       # Top-level figures (cross-phase)
│       └── *.png
│
├── dbcan_db/                          # Reference databases (shared)
├── antismash_hmms/                    # antiSMASH/BGC reference HMMs (shared)
├── reference/                         # Reference data (shared)
└── example/                           # Test data (version controlled)
```

> **Archival note:** `data/archives/` below is the in-repo convention. In practice,
> completed datasets are usually moved to **external cold storage** (e.g. an external
> drive) to reclaim space rather than kept under `data/`. Use whichever fits, but keep
> the database + `exploration/` + `reports/` + `figures/` together so the analysis stays
> reproducible.

## New Genome Checklist

### 1. Pre-Analysis Setup

```bash
# Archive previous analysis if complete
PREV_ORGANISM="thorarchaeota"
ARCHIVE_DATE=$(date +%Y-%m-%d)
mkdir -p data/archives/${PREV_ORGANISM}_${ARCHIVE_DATE}
mv data/${PREV_ORGANISM}_production/* data/archives/${PREV_ORGANISM}_${ARCHIVE_DATE}/

# Create new production directory
NEW_ORGANISM="heimdall_megavirus"  # lowercase, underscores
mkdir -p data/${NEW_ORGANISM}_production/{source,annotations,embeddings,structures,exploration,survey,reports,figures}
```

### 2. Required Inputs

| Input | Required | Format | Notes |
|-------|----------|--------|-------|
| Genome FASTA | Yes | `.fna`, `.fa`, `.fasta` | Standard Sharur ingest starts here |
| Protein FASTA | Optional | `.faa` or `.faa.gz` | Fallback path only when assemblies are unavailable |
| Astra databases | Yes | Installed locally | Standard: PFAM, KOFAM, HydDB, DefenseFinder, dbCAN |
| MinCED | Yes | Executable | Required for standard CRISPR array ingest |
| Optional Astra DBs | Optional | Installed locally | TXSScan, VOGdb, CANT-HYD |

### 3. Ingestion Workflow

```bash
mkdir -p data/${NEW_ORGANISM}_production/{source,annotations,embeddings,structures,exploration,survey,reports,figures}

sharur-ingest \
  --input-dir /path/to/genomes \
  --data-dir data/${NEW_ORGANISM}_production \
  --output data/${NEW_ORGANISM}_production/sharur.duckdb \
  --profile auto
```

Resume is ledger-backed and enabled by default. Use `--no-resume` or `--force` for an
intentional rerun. For manual stage-by-stage execution and resource-profile details, see
`QUICKSTART.md` and `src/ingest/README.md`.

If you only have pre-called proteins, you can bootstrap a database with:

```bash
python scripts/ingest_protein_fasta.py \
  /path/to/proteins.faa.gz \
  --output data/${NEW_ORGANISM}_production/sharur.duckdb
```

That fallback path does not replace the standard stage pipeline, and it does not run Astra, MinCED, or Stage 07 for you.

### 4. Standard Naming Conventions

| Item | Convention | Example |
|------|------------|---------|
| Production directory | `{organism}_production` | `heimdall_megavirus_production` |
| Archive directory | `{organism}_{YYYY-MM-DD}` | `thorarchaeota_2026-02-03` |
| Source file | `source/{sample}.fna` | `source/Heimdall_Megavirus.fna` |
| Stage outputs | `stageXX_*` | `stage04_astra/` |
| Annotations | `annotations/{db}.tsv` | `annotations/pfam.tsv` |
| Embeddings | `embeddings/` | Canonical H5 + persistent FAISS/SQLite sidecars |
| Reports | `reports/{type}.pdf` | `reports/final.pdf` |
| Figures | `figures/{locus_name}.png` | `figures/escrt_locus.png` |

### 5. What Gets Archived

**Always archive:**
- `sharur.duckdb` - Main database
- `sharur_ops.db` - Run/stage history, task leases, and coordination state
- `manifest.json` - Analysis state
- `exploration/` - All findings and state
- `reports/` - Synthesis documents
- `figures/` - Key visualizations

**Archive if valuable:**
- `structures/` - ESM3 PDBs (can regenerate but slow)
- `annotations/` - Processed annotation TSVs

**Don't archive (regenerable):**
- `embeddings/` - Can regenerate from sequences
- Raw HMM outputs (if kept separately)

### 6. Gitignore Verification

The `.gitignore` should cover:
```
data/*                    # All data directories
!data/.gitkeep           # Keep structure
!data/example/           # Example data OK
*.duckdb                 # Databases
*.fasta, *.faa, *.fna    # Sequence files
*.pt, *.pth, *.h5        # Model weights
*.lance/                 # Vector stores
```

**Before committing, verify:**
```bash
git status --porcelain | grep "^?" | head -20
# Should NOT show any large data files
```

## Analysis Phases

### Phase 1: Ingestion & Annotation
- Run the staged assembly-based ingest pipeline
- Build `sharur.duckdb` via Stage 07
- Compute ESM2 embeddings

### Phase 2: Initial Exploration
- Run `/explore` or genome browser agents
- Document interesting loci to `findings.jsonl`
- Generate neighborhood figures for key regions

### Phase 3: Structure Prediction (if warranted)
- Select candidates (unannotated, giant, interesting)
- Run ESM3 prediction
- Search Foldseek for remote homologs
- Research PDB hits for functions

### Phase 4: Synthesis & Reporting
- Run report generator script
- Review and refine
- Archive when complete

## Organism-Specific Notes

### Single Genome (e.g., Thorarchaeota)
- No genome profiling phase
- Focus on locus-by-locus analysis
- Structure prediction valuable for novel proteins

### Multi-Genome Dataset (e.g., Altiarchaeota)
- Run genome profiling first
- Identify core vs. variable features
- Consider ecotype/subclade analysis
- Spawn parallel agents for genome deep-dives

### Viral Genomes (e.g., Megavirus)
- Typically smaller protein count
- Focus on novel genes, host interactions
- Compare to known giant virus families
- Look for host-derived genes

## Migration from Legacy Structure

Datasets created before 2026-02-04 may have the old structure:
- `proteins/` → `source/`
- `pfam_results/`, `kofam_results/`, `vogdb_results/` → `annotations/`
- `stage06_embeddings/` → `embeddings/`
- `*.pdf`, `*.md` reports at root → `reports/`
- `foldseek_results.tsv` at root → `structures/`

The Sharur code supports both structures for backwards compatibility.
