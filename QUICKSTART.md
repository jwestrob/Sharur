# Sharur Quickstart: New Dataset Ingestion

**Goal:** Go from genome assembly FASTAs to an exploration-ready Sharur dataset.

## Canonical Path

The standard ingest workflow starts from nucleotide assemblies (`.fna`, `.fa`, `.fasta`) and should be run through `sharur-ingest`.

- Use `sharur-ingest` as the default interface for new dataset ingestion.
- Use the staged scripts in `src/ingest/` only when you need manual stage control, debugging, or rerunning one stage.

- Do not call `astra search` directly. Use `src/ingest/04_astra_scan.py`.
- Do not manually insert annotation rows into DuckDB. Use `src/ingest/07_build_knowledge_base.py`.
- Do not skip `src/ingest/minced_crispr.py` in the standard pipeline.

If you only have pre-called proteins and no assemblies, see [Alternative: Protein-Only Ingest](#alternative-protein-only-ingest). That path is a special case, not the default workflow.

## Prerequisites

### Software
- Python 3.10+ with Sharur installed: `pip install -e ".[dev]"`
- Prodigal
- [Astra](https://github.com/Dreycey/Astra)
- MinCED
- Optional but recommended GPU access for Stage 06 embeddings

If `sharur-ingest` is not available after install, refresh the editable install:
`pip install -e ".[dev]"`

### Astra databases
- Standard: `PFAM`, `KOFAM`, `HydDB`, `DefenseFinder`, `dbCAN`
- Optional: `TXSScan`, `VOGdb`, `CANT-HYD`

### Input data
- Genome assembly FASTAs in one directory
- Extensions supported by Stage 00: `.fna`, `.fa`, `.fasta`

## Standard Pipeline

```bash
sharur-ingest \
  --input-dir /path/to/genome_fastas \
  --data-dir data/my_dataset \
  --output data/my_dataset/sharur.duckdb \
  --force
```

This is the primary interface for running the standard pipeline. It orchestrates the staged workflow, including the standard stages 00, 03, 04, 05c, 07, and 06, while exposing skip flags for optional stages when needed.

## Manual Stage-by-Stage Pipeline

Use this only if you need direct control over individual stages or you are debugging a specific stage failure.

```bash
DATASET=my_dataset
INPUT=/path/to/genome/fastas

mkdir -p data/$DATASET/{source,annotations,embeddings,structures,exploration,figures,reports,survey}

# Stage 00: validate and organize assemblies
python src/ingest/00_prepare_inputs.py \
  -i $INPUT \
  -o data/$DATASET/stage00_prepared \
  --copy

# Stage 03: gene calling
python src/ingest/03_prodigal.py \
  -i data/$DATASET/stage00_prepared \
  -o data/$DATASET/stage03_prodigal \
  --max-workers 8

# Stage 04: annotation
python src/ingest/04_astra_scan.py \
  -i data/$DATASET/stage03_prodigal \
  -o data/$DATASET/stage04_astra \
  -t 12

# Stage 05c: CRISPR arrays
python src/ingest/minced_crispr.py \
  -i data/$DATASET/stage00_prepared \
  -o data/$DATASET/stage05c_crispr

# Stage 07: build DuckDB knowledge base + V2 predicates
python src/ingest/07_build_knowledge_base.py \
  -d data/$DATASET \
  -o data/$DATASET/sharur.duckdb

# Stage 06: embeddings (run after Stage 07)
python src/ingest/06_esm2_embeddings.py \
  data/$DATASET/stage03_prodigal \
  data/$DATASET/embeddings/
```

This is the manual reference sequence behind `sharur-ingest`. `04_astra_scan.py` supplies the correct per-database flags, and `07_build_knowledge_base.py` is the loader that consolidates proteins, annotations, loci, validated systems, and V2 predicates.

## What Each Stage Does

- `00_prepare_inputs.py`: validates assemblies and creates a processing manifest
- `03_prodigal.py`: produces per-genome protein FASTAs and the `all_protein_symlinks/` directory used by Stage 04
- `04_astra_scan.py`: runs PFAM, KOFAM, HydDB, DefenseFinder, and dbCAN with Sharur's expected settings
- `minced_crispr.py`: finds CRISPR repeat-spacer arrays from nucleotide assemblies
- `07_build_knowledge_base.py`: loads stage outputs into `sharur.duckdb`, integrates supported validation steps, writes `semantic_atoms`/`semantic_state`, and materializes legacy-compatible predicates
- `06_esm2_embeddings.py`: produces `embeddings/protein_embeddings.h5` for similarity search and ELSA

## Verify the Dataset

```bash
python - <<'PY'
from sharur.operators import Sharur

b = Sharur("data/my_dataset/sharur.duckdb")
total = b.store.execute("SELECT COUNT(*) FROM proteins")[0][0]
annotated = b.store.execute("SELECT COUNT(DISTINCT protein_id) FROM annotations")[0][0]
v2_states = b.store.execute("SELECT COUNT(*) FROM semantic_state")[0][0]
flat_predicates = b.store.execute("SELECT COUNT(*) FROM protein_predicates")[0][0]

print(f"Total proteins: {total}")
print(f"Annotated proteins: {annotated}")
print(f"V2 semantic states: {v2_states}")
print(f"Compatibility predicate rows: {flat_predicates}")
print(f"Embeddings loaded: {b.session.vector_store is not None}")
PY
```

Sanity checks:
- `proteins` should be non-zero
- `annotations` should be non-zero after Stage 04 + Stage 07
- `semantic_state` should match the protein count after Stage 07
- `protein_predicates` should match the protein count after Stage 07
- `Embeddings loaded: True` after Stage 06

## Start Exploring

### Claude Code

```bash
/survey
/explore --focus metabolism
```

### Python API

```python
from sharur.operators import Sharur

b = Sharur("data/my_dataset/sharur.duckdb")

giants = b.search_by_predicates(has=["giant", "unannotated"])
defense = b.search_by_predicates(has=["crispr_associated"])
similar = b.find_similar(giants._raw[0], k=20)
```

## Alternative: Protein-Only Ingest

If assemblies are unavailable and you only have pre-called proteins, you can bootstrap a database with:

```bash
python scripts/ingest_protein_fasta.py \
  /path/to/proteins.faa.gz \
  --output data/$DATASET/sharur.duckdb
```

Use this only when the standard assembly-based pipeline is impossible. It:
- loads bins, contigs, and proteins
- does not run Prodigal, Astra, or MinCED
- does not replace Stage 04 or Stage 07
- leaves annotation loading and downstream validation to you

## Common Mistakes

- Reconstructing the standard pipeline manually when `sharur-ingest` would do the job
- Running raw `astra search` commands instead of `src/ingest/04_astra_scan.py`
- Passing `--databases` as a space-separated list instead of repeated `-d` flags when overriding defaults
- Skipping `minced_crispr.py` and then expecting CRISPR array loci in DuckDB
- Treating legacy Stage `05b` dbCAN as the standard CAZyme path; standard ingest uses Stage 04 + Stage 07

## More Detail

- Full manual stage reference: `src/ingest/README.md`
- Tool-specific details: `docs/tools_reference.md`
- Dataset layout and archival conventions: `DATA_ORGANIZATION.md`
