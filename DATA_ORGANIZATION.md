# Data Organization for Bennu Analyses

This document describes the standard organization for genome analysis projects.

## Directory Structure

```
data/
├── archives/                           # Completed analyses (read-only)
│   ├── {organism}_{YYYY-MM-DD}/       # Archived with completion date
│   └── ...
├── {organism}_production/              # Active analysis
│   ├── bennu.duckdb                   # Main database
│   ├── manifest.json                  # Analysis state and session continuity
│   │
│   ├── source/                        # Input files
│   │   └── proteins.faa               # Source protein FASTA
│   │
│   ├── annotations/                   # Annotation results (consolidated)
│   │   ├── pfam.tsv                   # PFAM domain hits
│   │   ├── kegg.tsv                   # KEGG KOfam results
│   │   ├── vogdb.tsv                  # VOGdb hits (viral)
│   │   └── cazy.tsv                   # CAZy hits (if applicable)
│   │
│   ├── embeddings/                    # ESM2 protein embeddings
│   │   ├── embedding_manifest.json
│   │   ├── lancedb/                   # LanceDB vector store
│   │   └── protein_embeddings.h5      # Raw embeddings (optional)
│   │
│   ├── structures/                    # ESM3 structure predictions
│   │   ├── *.pdb                      # Predicted PDB files
│   │   └── foldseek_results.tsv       # Structural homology searches
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
├── reference/                         # Reference data (shared)
└── example/                           # Test data (version controlled)
```

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
| Protein FASTA | Yes | `.faa` or `.faa.gz` | Prodigal-called or from annotation |
| PFAM annotations | Recommended | TSV/domtblout | From hmmsearch or Astra |
| KEGG annotations | Recommended | TSV | From KOfamScan or Astra |
| HydDB annotations | Recommended | TSV | Hydrogenase classification (Astra --installed_hmms HydDB --cut_ga) |
| Genome FASTA | Optional | `.fna` | For contig context |
| GFF annotations | Optional | `.gff3` | Gene coordinates |

### 3. Ingestion Workflow

```bash
# Step 1: Ingest proteins
python scripts/ingest_protein_fasta.py \
  --fasta /path/to/proteins.faa.gz \
  --output data/${NEW_ORGANISM}_production/bennu.duckdb

# Step 2: Add annotations (if available)
# Via Astra pipeline or direct import

# Step 3: Generate predicates
python -c "
from bennu.operators import Bennu
b = Bennu('data/${NEW_ORGANISM}_production/bennu.duckdb')
b.regenerate_predicates()
"

# Step 4: Generate embeddings (ESM2)
python src/ingest/06_esm2_embeddings.py \
  data/${NEW_ORGANISM}_production/bennu.duckdb \
  data/${NEW_ORGANISM}_production/embeddings/
```

### 4. Standard Naming Conventions

| Item | Convention | Example |
|------|------------|---------|
| Production directory | `{organism}_production` | `heimdall_megavirus_production` |
| Archive directory | `{organism}_{YYYY-MM-DD}` | `thorarchaeota_2026-02-03` |
| Source file | `source/proteins.faa` | `source/Heimdall_Megavirus.faa` |
| Annotations | `annotations/{db}.tsv` | `annotations/pfam.tsv` |
| Embeddings | `embeddings/` | Standard LanceDB structure |
| Reports | `reports/{type}.pdf` | `reports/final.pdf` |
| Figures | `figures/{locus_name}.png` | `figures/escrt_locus.png` |

### 5. What Gets Archived

**Always archive:**
- `bennu.duckdb` - Main database
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
- Ingest protein FASTA
- Run/import PFAM and KEGG annotations
- Generate predicates
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

The Bennu code supports both structures for backwards compatibility.
