# Compare Skill

Comparative genomics between genomes or genome groups.

**CRITICAL: You are a leaf agent. DO NOT spawn sub-agents or use the Task tool.**

**CONCURRENCY: DuckDB does not support concurrent writes. Only ONE agent should access a database at a time. The coordinator must run DB-accessing skills sequentially, not in parallel.**

> **Mandatory:** Follow the shared validation protocols in `_validation_protocols.md`.
> Verify accession names before reporting. Use COUNT(DISTINCT protein_id) for protein
> counts. Apply Context-First protocol for annotations averaging >10 hits/genome.

---

## Usage

```
/compare GENOME_A GENOME_B           # Pairwise comparison
/compare --all                       # All-vs-all comparison matrix
/compare --groups group1.txt group2.txt  # Group comparison
```

---

## Prompt

You are performing comparative genomic analysis between genomes in a metagenomic dataset.

### Pairwise Comparison

```python
from bennu.operators import Bennu
import pandas as pd

b = Bennu("data/DATASET/bennu.duckdb")

genome_a = "GENOME_A"
genome_b = "GENOME_B"

# Basic stats
stats = b.store.execute(f"""
    SELECT
        bin_id,
        COUNT(*) as n_proteins,
        AVG(sequence_length) as avg_length,
        MAX(sequence_length) as max_length
    FROM proteins
    WHERE bin_id IN ('{genome_a}', '{genome_b}')
    GROUP BY bin_id
""")
print("Basic Statistics:")
for row in stats:
    print(f"  {row[0]}: {row[1]} proteins, avg {row[2]:.0f} aa, max {row[3]} aa")
```

### Feature Comparison

```python
# Compare predicate profiles
KEY_PREDICATES = [
    'hydrogenase', 'crispr_associated', 'transposase', 'adhesin',
    'toxin_antitoxin', 'restriction_modification', 's_layer',
    'wood_ljungdahl', 'electron_transport', 'unannotated', 'giant'
]

comparison = []
for pred in KEY_PREDICATES:
    counts = b.store.execute(f"""
        SELECT p.bin_id, COUNT(*) as n
        FROM proteins p
        JOIN protein_predicates pp ON p.protein_id = pp.protein_id
        WHERE p.bin_id IN ('{genome_a}', '{genome_b}')
          AND '{pred}' = ANY(pp.predicates)
        GROUP BY p.bin_id
    """)
    count_dict = {row[0]: row[1] for row in counts}
    comparison.append({
        'predicate': pred,
        genome_a: count_dict.get(genome_a, 0),
        genome_b: count_dict.get(genome_b, 0),
        'diff': count_dict.get(genome_a, 0) - count_dict.get(genome_b, 0)
    })

df = pd.DataFrame(comparison)
print("\nFeature Comparison:")
print(df.to_string(index=False))
```

### Domain Profile Comparison

```python
# Top domains in each genome
def get_top_domains(genome, n=20):
    return b.store.execute(f"""
        SELECT a.name, COUNT(*) as n
        FROM annotations a
        JOIN proteins p ON a.protein_id = p.protein_id
        WHERE p.bin_id = '{genome}' AND a.source = 'pfam'
        GROUP BY a.name
        ORDER BY n DESC
        LIMIT {n}
    """)

domains_a = set(row[0] for row in get_top_domains(genome_a, 50))
domains_b = set(row[0] for row in get_top_domains(genome_b, 50))

shared = domains_a & domains_b
unique_a = domains_a - domains_b
unique_b = domains_b - domains_a

print(f"\nDomain overlap:")
print(f"  Shared: {len(shared)}")
print(f"  Unique to {genome_a}: {len(unique_a)}")
print(f"  Unique to {genome_b}: {len(unique_b)}")

if unique_a:
    print(f"\n  {genome_a} specific: {', '.join(list(unique_a)[:10])}")
if unique_b:
    print(f"\n  {genome_b} specific: {', '.join(list(unique_b)[:10])}")
```

### Ortholog Detection (via embedding similarity)

```python
# Find potential orthologs using ESM2 embeddings
# Proteins from A that are most similar to proteins in B

orthologs = b.store.execute(f"""
    WITH proteins_a AS (
        SELECT protein_id FROM proteins WHERE bin_id = '{genome_a}'
    ),
    proteins_b AS (
        SELECT protein_id FROM proteins WHERE bin_id = '{genome_b}'
    )
    -- This requires vector similarity search
    -- Implementation depends on LanceDB setup
""")

# Alternative: use annotation-based ortholog detection
shared_kos = b.store.execute(f"""
    SELECT a1.accession as ko,
           a1.protein_id as prot_a,
           a2.protein_id as prot_b
    FROM annotations a1
    JOIN annotations a2 ON a1.accession = a2.accession
    JOIN proteins p1 ON a1.protein_id = p1.protein_id
    JOIN proteins p2 ON a2.protein_id = p2.protein_id
    WHERE a1.source = 'kegg' AND a2.source = 'kegg'
      AND p1.bin_id = '{genome_a}'
      AND p2.bin_id = '{genome_b}'
""")
print(f"\nShared KEGG orthologs: {len(shared_kos)}")
```

### All-vs-All Comparison Matrix

```python
# Get all genomes
genomes = [row[0] for row in b.store.execute("SELECT DISTINCT bin_id FROM proteins")]

# Build feature matrix
import numpy as np

features = ['hydrogenase', 'crispr_associated', 'transposase', 'adhesin',
            'wood_ljungdahl', 'sulfur_metabolism', 'nitrogen_fixation']

matrix = np.zeros((len(genomes), len(features)))

for i, genome in enumerate(genomes):
    for j, feat in enumerate(features):
        count = b.store.execute(f"""
            SELECT COUNT(*) FROM proteins p
            JOIN protein_predicates pp ON p.protein_id = pp.protein_id
            WHERE p.bin_id = '{genome}' AND '{feat}' = ANY(pp.predicates)
        """)[0][0]
        # Normalize by genome size
        genome_size = b.store.execute(f"SELECT COUNT(*) FROM proteins WHERE bin_id = '{genome}'")[0][0]
        matrix[i, j] = count * 1000 / genome_size if genome_size > 0 else 0

# Cluster genomes by feature profile
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist

Z = linkage(matrix, method='ward')
```

### Visualization

```python
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

FIGURES_DIR = Path("data/DATASET/exploration/figures")
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Pairwise feature comparison bar chart
fig, ax = plt.subplots(figsize=(10, 6))
x = range(len(KEY_PREDICATES))
width = 0.35

ax.bar([i - width/2 for i in x], [c[genome_a] for c in comparison], width, label=genome_a)
ax.bar([i + width/2 for i in x], [c[genome_b] for c in comparison], width, label=genome_b)

ax.set_xticks(x)
ax.set_xticklabels(KEY_PREDICATES, rotation=45, ha='right')
ax.legend()
ax.set_ylabel('Count')
ax.set_title(f'Feature Comparison: {genome_a} vs {genome_b}')
plt.tight_layout()
plt.savefig(FIGURES_DIR / f"compare_{genome_a}_{genome_b}.png", dpi=150)

# All-vs-all heatmap
fig, ax = plt.subplots(figsize=(12, 10))
sns.clustermap(
    pd.DataFrame(matrix, index=genomes, columns=features),
    cmap='viridis',
    figsize=(12, 10)
)
plt.savefig(FIGURES_DIR / "genome_feature_clustering.png", dpi=150, bbox_inches='tight')
```

### Output Format

```markdown
## Comparative Analysis: {genome_a} vs {genome_b}

### Basic Statistics
| Metric | {genome_a} | {genome_b} |
|--------|------------|------------|
| Proteins | X | Y |
| Avg length | X aa | Y aa |
| Max length | X aa | Y aa |

### Feature Comparison
[Table of predicate counts with differences highlighted]

### Unique Features
**{genome_a} has but {genome_b} lacks:**
- [List distinctive features]

**{genome_b} has but {genome_a} lacks:**
- [List distinctive features]

### Shared Features
- [Core features present in both]

### Domain Analysis
- Shared top domains: X
- Unique to {genome_a}: Y
- Unique to {genome_b}: Z

### Biological Interpretation
[What do the differences suggest about lifestyle/ecology?]

### Figures
- compare_{genome_a}_{genome_b}.png: Feature comparison
```
