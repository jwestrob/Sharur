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
pfam.columns = ['protein_id', 'annotation_id', 'evalue', 'score', 'bias']  # Adjust if needed
db.execute('INSERT INTO annotations SELECT * FROM pfam')
db.close()
"

# Import KEGG annotations
python -c "
import duckdb
import pandas as pd

db = duckdb.connect('data/${DATASET}/sharur.duckdb')
kegg = pd.read_csv('data/${DATASET}/annotations/kofam_results/KOFAM_hits_df.tsv', sep='\t')
kegg.columns = ['protein_id', 'annotation_id', 'evalue', 'score', 'bias']
db.execute('INSERT INTO annotations SELECT * FROM kegg')
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
b.regenerate_predicates()  # Uses PFAM, KEGG, HydDB, VOGdb
print(f'Generated {len(b.list_predicates())} unique predicates')
"
```

**Output:** Predicates stored in database, queryable via `b.search_by_predicates()`

---

## Step 4: Generate ESM2 Embeddings

```bash
python src/ingest/06_esm2_embeddings.py \
  data/${DATASET}/sharur.duckdb \
  data/${DATASET}/embeddings/
```

**Note:** Requires GPU for speed (CPU works but slow). Uses `facebook/esm2_t6_8M_UR50D` model.

**Output:** LanceDB vector store in `embeddings/` for similarity search

---

## Step 5: Verify Database

```bash
python -c "
from sharur.operators import Sharur
b = Sharur('data/${DATASET}/sharur.duckdb')

print(f'Total proteins: {b.total_proteins()}')
print(f'Annotated: {len(b.search(\"confident_hit\"))} proteins')
print(f'Unannotated: {len(b.search(\"unannotated\"))} proteins')
print(f'Available predicates: {len(b.list_predicates())}')

# Test similarity search
print(f'Embeddings loaded: {b.vector_store is not None}')
"
```

**Expected output:**
```
Total proteins: 10234
Annotated: 7891 proteins
Unannotated: 2343 proteins
Available predicates: 1523
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
giants = b.search("giant AND unannotated")
defense = b.search("crispr_associated OR restriction_modification")
transporters = b.search("transporter AND membrane_protein")

# Examine neighborhoods
for pid in giants[:5]:
    b.visualize_neighborhood(pid, window=10,
                            output_path=f"figures/{pid}_neighborhood.png")

# Find similar proteins
similar = b.find_similar(giants[0], k=20)
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
- Check annotation_id format (should be like `PF00001`, `K00001`)
- Run with verbose: `b.regenerate_predicates(verbose=True)`

### "Embeddings taking forever"
- Use GPU if available
- Reduce batch size if OOM errors
- ESM2 model will auto-download first time (~200MB)

### "Vector store not loading"
- Check `embeddings/` directory has LanceDB files
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
├── embeddings/               # LanceDB vector store (100-500 MB)
├── exploration/              # Created during analysis
├── figures/                  # Created during analysis
└── reports/                  # Created during analysis
```

**You're ready to explore!** 🚀
