# Sharur Quickstart: New Dataset Ingestion

**Goal:** Go from raw protein FASTA to exploration-ready Sharur database in ~30 minutes.

---

## Prerequisites

### Software
- Python 3.8+ with Sharur installed (`pip install -e .`)
- [Astra](https://github.com/Dreycey/Astra) annotation tool
- HMM databases installed via Astra:
  - `PFAM` (Pfam-A domains)
  - `KOFAM` (KEGG orthologs)
  - `HydDB` (hydrogenases)
  - `VOGdb` (virus orthologous groups, if viral)
  - `DefenseFinder` (optional, defense systems)

### Input Data
- **Protein FASTA** (`.faa` or `.faa.gz`)
  - From Prodigal, DFAST, or other gene caller
  - Must have unique protein IDs

---

## Step 0: Setup

```bash
# Set your dataset name (lowercase, underscores)
export DATASET="my_dataset_production"

# Create directory structure
mkdir -p data/${DATASET}/{source,annotations,embeddings,structures,exploration,figures,reports}

# Copy your protein FASTA
cp /path/to/your/proteins.faa.gz data/${DATASET}/source/
```

---

## Step 1: Run Astra Annotations

**This is the longest step (~10-20 minutes depending on dataset size).**

```bash
# Prepare: Astra needs a directory of FASTA files, not a single file
mkdir -p data/${DATASET}/source/proteins_dir
gunzip -c data/${DATASET}/source/proteins.faa.gz > data/${DATASET}/source/proteins_dir/proteins.faa

# Run PFAM (domain annotations)
astra search \
  --installed_hmms PFAM \
  --prot_in data/${DATASET}/source/proteins_dir \
  --outdir data/${DATASET}/annotations/pfam_results \
  --cut_ga \
  --threads 12

# Run KOFAM (KEGG orthologs)
astra search \
  --installed_hmms KOFAM \
  --prot_in data/${DATASET}/source/proteins_dir \
  --outdir data/${DATASET}/annotations/kofam_results \
  --cut_ga \
  --threads 12

# Run HydDB (hydrogenases, if prokaryotic)
astra search \
  --installed_hmms HydDB \
  --prot_in data/${DATASET}/source/proteins_dir \
  --outdir data/${DATASET}/annotations/hyddb_results \
  --cut_ga \
  --threads 12

# Run VOGdb (if viral dataset)
astra search \
  --installed_hmms VOGdb \
  --prot_in data/${DATASET}/source/proteins_dir \
  --outdir data/${DATASET}/annotations/vogdb_results \
  --cut_ga \
  --threads 12

# Cleanup temp directory
rm -r data/${DATASET}/source/proteins_dir
```

**Output:** Annotation TSV files in `annotations/*/`

---

## Step 2: Ingest into Sharur Database

```bash
# Import proteins
python scripts/ingest_protein_fasta.py \
  --fasta data/${DATASET}/source/proteins.faa.gz \
  --output data/${DATASET}/sharur.duckdb

# Import PFAM annotations
python -c "
import duckdb
import pandas as pd

db = duckdb.connect('data/${DATASET}/sharur.duckdb')
pfam = pd.read_csv('data/${DATASET}/annotations/pfam_results/PFAM_hits_df.tsv', sep='\t')
pfam.columns = ['protein_id', 'accession', 'evalue', 'score', 'bias']  # Adjust if needed
# Note: annotations table schema is (annotation_id, protein_id, source, accession, name, description, evalue, score, start_aa, end_aa)
# Use Stage 07 (07_build_knowledge_base.py) for proper loading instead of manual INSERT
db.close()
"

# Import KEGG annotations
python -c "
import duckdb
import pandas as pd

db = duckdb.connect('data/${DATASET}/sharur.duckdb')
kegg = pd.read_csv('data/${DATASET}/annotations/kofam_results/KOFAM_hits_df.tsv', sep='\t')
kegg.columns = ['protein_id', 'accession', 'evalue', 'score', 'bias']
# Note: annotations table schema is (annotation_id, protein_id, source, accession, name, description, evalue, score, start_aa, end_aa)
# Use Stage 07 (07_build_knowledge_base.py) for proper loading instead of manual INSERT
db.close()
"

# Repeat for HydDB, VOGdb, etc.
```

**Alternative:** Use `src/ingest/07_build_knowledge_base.py` if you have all annotations ready.

---

## Step 3: Generate Predicates

```bash
python -c "
from sharur.operators import Sharur
b = Sharur('data/${DATASET}/sharur.duckdb')
print('Generating predicates from annotations...')
# Use the regenerate_predicates script instead:
# python scripts/regenerate_predicates.py data/${DATASET}/sharur.duckdb
from sharur.predicates.vocabulary import ALL_PREDICATES
print(f'Vocabulary contains {len(ALL_PREDICATES)} defined predicates')
"
```

**Output:** Predicates stored in database, queryable via `b.search_by_predicates()`

---

## Step 4: Generate ESM2 Embeddings

```bash
python src/ingest/06_esm2_embeddings.py \
  data/${DATASET}/stage03_prodigal \
  data/${DATASET}/embeddings/
```

**Note:** Requires GPU for speed (CPU works but slow). Uses `facebook/esm2_t6_8M_UR50D` model.

**Output:** FAISS vector index in `embeddings/` for similarity search

---

## Step 5: Verify Database

```bash
python -c "
from sharur.operators import Sharur
b = Sharur('data/${DATASET}/sharur.duckdb')

total = b.store.execute('SELECT COUNT(*) FROM proteins')[0][0]
print(f'Total proteins: {total}')
annotated = b.store.execute(\"SELECT COUNT(DISTINCT protein_id) FROM annotations\")[0][0]
print(f'Annotated: {annotated} proteins')
unannotated = total - annotated
print(f'Unannotated: {unannotated} proteins')
from sharur.predicates.vocabulary import ALL_PREDICATES
print(f'Defined predicates: {len(ALL_PREDICATES)}')

# Test similarity search
print(f'Embeddings loaded: {b.vector_store is not None}')
"
```

**Expected output:**
```
Total proteins: 10234
Annotated: 7891 proteins
Unannotated: 2343 proteins
Defined predicates: 1523
Embeddings loaded: True
```

---

## Step 6: Start Exploring

### Option A: Claude Code Skills

If using Claude Code:
```bash
/explore --focus metabolism
/survey  # Comprehensive systematic survey
```

### Option B: Programmatic Analysis

```python
from sharur.operators import Sharur

b = Sharur('data/my_dataset_production/sharur.duckdb')

# Find interesting proteins
giants = b.search_by_predicates(has=["giant", "unannotated"])
defense = b.search_by_predicates(has=["crispr_associated"])
transporters = b.search_by_predicates(has=["transporter", "membrane_protein"])

# Examine neighborhoods (use ._raw to get underlying list from SharurResult)
giant_proteins = giants._raw
for pid in giant_proteins[:5]:
    b.visualize_neighborhood(pid, window=10,
                            output_path=f"figures/{pid}_neighborhood.png")

# Find similar proteins
similar = b.find_similar(giant_proteins[0], k=20)
```

### Option C: Generate Report

```bash
# Copy and adapt canonical report generator
cp scripts/generate_viral_genome_report.py scripts/generate_${DATASET}_report.py

# Edit ONLY these in the top section:
# - DB_PATH = "data/${DATASET}/sharur.duckdb"
# - DATA_DIR = Path("data/${DATASET}")
# - Update organism name in headers/titles
# DO NOT change OUTPUT_PDF - it uses standard "COMPREHENSIVE_REPORT.pdf"

# Generate report (outputs to data/${DATASET}/COMPREHENSIVE_REPORT.pdf)
python scripts/generate_${DATASET}_report.py
```

---

## Troubleshooting

### "No annotations found"
- Check Astra output TSV files exist and have data
- Verify column names match expected format
- Use `--cut_ga` flag with Astra for trusted cutoffs

### "Predicates not generating"
- Ensure annotations are in database: `SELECT COUNT(*) FROM annotations`
- Check accession format (should be like `PF00001`, `K00001`)
- Regenerate predicates: `python scripts/regenerate_predicates.py data/${DATASET}/sharur.duckdb`

### "Embeddings taking forever"
- Use GPU if available
- Reduce batch size if OOM errors
- ESM2 model will auto-download first time (~200MB)

### "Vector store not loading"
- Check `embeddings/` directory has `protein_embeddings.h5`
- Verify protein_id column exists: `db.open_table('protein_embeddings')`
- Regenerate if corrupted

---

## Next Steps

Once your database is ready:

1. **Exploration:** Use `/explore` or `/survey` skills (Claude Code)
2. **Structure prediction:** Select interesting unknowns, run ESM3
3. **Foldseek searches:** Find remote homologs via structural similarity
4. **Comparative analysis:** Load multiple datasets, use `/compare` skill
5. **Generate report:** Use canonical report generator for publication

See `AGENTS.md` for full agent workflows and canonical tools.

---

## Time Estimates

| Step | Time | Notes |
|------|------|-------|
| Astra annotations | 10-20 min | Depends on dataset size, CPU cores |
| Database ingestion | 1-2 min | Fast |
| Predicate generation | 30 sec | Depends on annotation count |
| ESM2 embeddings | 5-15 min | GPU: fast, CPU: slow |
| **Total** | **~30 min** | For typical dataset (10k proteins) |

---

## File Checklist

After completion, you should have:
```
data/my_dataset_production/
├── sharur.duckdb              # Core database (10-100 MB)
├── source/
│   └── proteins.faa.gz       # Original input
├── annotations/
│   ├── pfam_results/PFAM_hits_df.tsv
│   ├── kofam_results/KOFAM_hits_df.tsv
│   └── hyddb_results/HydDB_hits_df.tsv
├── embeddings/               # FAISS vector index (100-500 MB)
├── exploration/              # Created during analysis
├── figures/                  # Created during analysis
└── reports/                  # Created during analysis
```

**You're ready to explore!** 🚀
