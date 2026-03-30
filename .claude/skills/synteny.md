# /synteny — ELSA Synteny Discovery

Run ELSA embedding-based synteny detection on a Sharur dataset. Finds conserved gene neighborhoods across genomes by chaining cosine-similar protein pairs into collinear blocks.

## Prerequisites

1. **Embeddings must exist.** Stage 06 (`protein_embeddings.h5`) is required. If missing, run:
   ```bash
   python src/ingest/06_esm2_embeddings.py data/DATASET/stage03_prodigal data/DATASET/embeddings/
   ```

2. **ELSA must be installed** in the active Python environment:
   ```bash
   pip install -e /path/to/ELSA
   pip install faiss-cpu h5py
   ```

3. **macOS OMP fix** — set before running:
   ```bash
   export KMP_DUPLICATE_LIB_OK=TRUE
   ```

## Workflow

### Step 1: Run ELSA synteny

```bash
elsa synteny \
    --db data/DATASET/sharur.duckdb \
    --embeddings data/DATASET/embeddings/protein_embeddings.h5 \
    --store data/DATASET/synteny/store \
    -o data/DATASET/synteny/
```

Use `--store` to persist the FAISS index. Subsequent runs can reload with just `--store` (no `--db`/`--embeddings`).

### Step 2: Interpret results

Two output files in the output directory:

**`micro_chain_blocks.csv`** (~1.6M rows for DPANN) — one row per syntenic block (collinear gene run between two contigs):
- `block_id`, `cluster_id` — identifiers
- `query_genome`, `target_genome` — the two genomes
- `n_anchors` — number of anchor gene pairs (2 = gene pair, 20+ = large operon)
- `chain_score` — sum of cosine similarities
- `orientation` — `"same"` or `"inverted"`
- `query_anchor_genes`, `target_anchor_genes` — JSON arrays of **Sharur protein_ids** (e.g., `["DATDYP010000003.1_26", "DATDYP010000003.1_29"]`)

**`micro_chain_clusters.csv`** (~68k rows for DPANN, 107 MB) — clusters of overlapping blocks:
- `cluster_id`, `size` — cluster identifier and block count
- `genome_support` — number of genomes sharing this syntenic region
- `genes_json` — JSON dict: `{genome_id: [gene_id, ...]}` where gene_ids use **`contig_id:gene_index`** format (e.g., `"DATDYP010000003.1:25"`)

#### CRITICAL: Gene ID formats differ between files

| File | Gene ID format | Example | How to map to Sharur |
|------|---------------|---------|---------------------|
| `micro_chain_blocks.csv` | Sharur `protein_id` | `DATDYP010000003.1_26` | Direct lookup in `proteins` table |
| `micro_chain_clusters.csv` | `contig_id:gene_index` | `DATDYP010000003.1:25` | `WHERE contig_id = 'DATDYP010000003.1' AND gene_index = 25` |

These are the SAME protein (gene_index=25 → protein_id ending in _26 for Prodigal-named proteins; for accession-style IDs like `HZX44842.1`, use the contig+gene_index lookup). Always verify mappings against the `proteins` table.

### Step 3: Querying ELSA results

#### Forward query: Find the most conserved syntenic regions

```python
import pandas as pd, json
from sharur.operators import Sharur

b = Sharur("data/DATASET/sharur.duckdb")
clusters = pd.read_csv("data/DATASET/synteny/results/micro_chain_clusters.csv")
top = clusters.nlargest(10, "genome_support")

for _, row in top.iterrows():
    genes_by_genome = json.loads(row["genes_json"])
    print(f"cluster {row.cluster_id}: {row.size} blocks, {row.genome_support} genomes")
```

#### Reverse query: Find which ELSA cluster(s) contain a specific protein

**This is the most common agent query. You MUST run this — never invent cluster IDs.**

```python
import duckdb, json, csv

db = duckdb.connect("data/DATASET/sharur.duckdb", read_only=True)

# 1. Convert protein_id → ELSA cluster gene_id format (contig:gene_index)
protein_id = "YOUR_PROTEIN_ID"
r = db.execute(f"""
    SELECT contig_id, gene_index, bin_id,
           contig_id || ':' || gene_index AS elsa_gene_id
    FROM proteins WHERE protein_id = '{protein_id}'
""").fetchone()
elsa_gene_id = r[3]  # e.g. "DATDYP010000003.1:25"

# 2. Search the clusters CSV for this gene_id
#    genes_json is a JSON dict, so the gene_id appears as a quoted string
matches = []
with open("data/DATASET/synteny/results/micro_chain_clusters.csv") as f:
    reader = csv.DictReader(f)
    for row in reader:
        if elsa_gene_id in row["genes_json"]:
            matches.append({
                "cluster_id": int(row["cluster_id"]),
                "size": int(row["size"]),
                "genome_support": int(row["genome_support"]),
            })
print(f"Protein {protein_id} found in {len(matches)} ELSA clusters:")
for m in matches:
    print(f"  cluster {m['cluster_id']}: {m['size']} blocks, {m['genome_support']} genomes")
```

#### Batch reverse query: Find clusters for a set of proteins (e.g., all LANC_like)

```python
import duckdb, json, csv

db = duckdb.connect("data/DATASET/sharur.duckdb", read_only=True)

# 1. Get all protein_ids with the domain of interest + their ELSA gene_ids
targets = db.execute("""
    SELECT p.protein_id, p.contig_id || ':' || p.gene_index AS elsa_gene_id, p.bin_id
    FROM annotations a
    JOIN proteins p ON a.protein_id = p.protein_id
    WHERE a.source = 'pfam' AND a.name = 'LANC_like'
""").fetchdf()
target_ids = set(targets.elsa_gene_id)

# 2. Scan clusters for any of these gene_ids
hits = []
with open("data/DATASET/synteny/results/micro_chain_clusters.csv") as f:
    reader = csv.DictReader(f)
    for row in reader:
        genes_json = row["genes_json"]
        found = [gid for gid in target_ids if gid in genes_json]
        if found:
            genes_by_genome = json.loads(genes_json)
            hits.append({
                "cluster_id": int(row["cluster_id"]),
                "size": int(row["size"]),
                "genome_support": int(row["genome_support"]),
                "matched_gene_ids": found,
                "genomes": list(genes_by_genome.keys()),
            })

for h in hits:
    print(f"cluster {h['cluster_id']}: {h['genome_support']} genomes, "
          f"matched {len(h['matched_gene_ids'])} target genes")
    print(f"  genomes: {h['genomes'][:5]}")
```

#### Searching blocks for pairwise synteny between specific genomes

```python
import pandas as pd, json

blocks = pd.read_csv("data/DATASET/synteny/results/micro_chain_blocks.csv")

# Find all syntenic blocks between two genomes
pair = blocks[
    (blocks.query_genome == "GCA_021801225.1") &
    (blocks.target_genome == "GCA_027330505.1")
]
# Or blocks containing a specific protein_id (blocks use protein_id directly)
protein_id = "YOUR_PROTEIN_ID"
mask = blocks.query_anchor_genes.str.contains(protein_id, na=False) | \
       blocks.target_anchor_genes.str.contains(protein_id, na=False)
protein_blocks = blocks[mask]
```

### Step 4: Citing ELSA results

**NEVER cite a cluster_id you didn't obtain from an actual query.** ELSA cluster IDs are arbitrary integers — they look plausible at any value, which makes fabrication undetectable without verification.

When reporting ELSA synteny evidence in findings:

```jsonl
{
  "id": "E202",
  "title": "Lanthipeptide BGCs in Micrarchaeota",
  "elsa_evidence": {
    "cluster_id": 170014,
    "genome_support": 3,
    "size": 3,
    "genomes": ["GCA_021801225.1", "GCA_027330505.1", "GCA_027354795.1"],
    "query": "batch reverse query for LANC_like proteins"
  }
}
```

**Required fields when citing ELSA:**
- `cluster_id` — the actual cluster_id from the CSV
- `genome_support` — from the CSV, not estimated
- `genomes` — list of genomes in the cluster (parsed from genes_json)
- `query` — brief description of how you found this cluster

**If you cannot run the query** (CSV not available, too large, etc.), say so explicitly:
"ELSA synteny validation not performed" — do NOT invent cluster IDs.

## Key Parameters

| Parameter | Default | When to change |
|-----------|---------|----------------|
| `--similarity-threshold` | 0.85 | Lower to 0.7 for divergent species (different phyla) |
| `--max-gap` | 2 | Increase for fragmented assemblies or loci with many hypotheticals |
| `--min-chain-size` | 2 | Increase to 3+ to filter noise, reduce for sparse datasets |
| `--min-genome-support` | 2 | Increase to require broader conservation |

## Python API

For programmatic access within agents:

```python
from pathlib import Path
from elsa.adapter import load_proteins_from_duckdb, load_embeddings_h5, build_genes_dataframe
from elsa.store import SyntenyStore
from elsa.analyze.pipeline import run_chain_pipeline, ChainConfig

proteins = load_proteins_from_duckdb("data/DATASET/sharur.duckdb")
embeddings = load_embeddings_h5("data/DATASET/embeddings/protein_embeddings.h5")
genes_df = build_genes_dataframe(proteins, embeddings, normalize=True)

store = SyntenyStore.create(Path("data/DATASET/synteny/store"), genes_df)
# Or reload: store = SyntenyStore.load(Path("data/DATASET/synteny/store"))

summary = run_chain_pipeline(
    output_dir=Path("data/DATASET/synteny/"),
    config=ChainConfig(),
    genes_df=store.get_genes_df(),
    prebuilt_index=store.get_index_tuple(),
)
print(f"blocks={summary.num_blocks}, clusters={summary.num_clusters}")
```

## Adding new genomes to an existing store

```bash
elsa synteny \
    --store data/DATASET/synteny/store \
    --add-db data/NEW_GENOMES/sharur.duckdb \
    --add-embeddings data/NEW_GENOMES/embeddings/protein_embeddings.h5 \
    -o data/DATASET/synteny/results_combined/
```

Deduplicates by `protein_id` — already-indexed proteins are skipped. The FAISS index is rebuilt on add (~2s for 100k vectors).

## Output Directory Convention

```
data/DATASET/
├── synteny/
│   ├── store/                  # Persistent FAISS store
│   │   ├── index.faiss
│   │   ├── metadata.parquet
│   │   ├── embeddings.npy
│   │   └── config.json
│   ├── micro_chain_blocks.csv  # Syntenic blocks
│   └── micro_chain_clusters.csv # Block clusters
```
