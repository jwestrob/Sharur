# Pathway Skill

Check completeness of metabolic pathways using KEGG annotations.

**CRITICAL: You are a leaf agent. DO NOT spawn sub-agents or use the Task tool.**

**CONCURRENCY: DuckDB does not support concurrent writes. Only ONE agent should access a database at a time. The coordinator must run DB-accessing skills sequentially, not in parallel.**

> **Mandatory:** Follow the shared validation protocols in `_validation_protocols.md`.
> Verify accession names before reporting. Use COUNT(DISTINCT protein_id) for protein
> counts. Apply Context-First protocol for annotations averaging >10 hits/genome.

---

## Usage

```
/pathway                              # Analyze all pathways
/pathway Wood-Ljungdahl               # Check specific pathway
/pathway --genome GCA_003598175.1     # Per-genome analysis
/pathway --compare                    # Cross-genome comparison
```

---

## Prompt

You are analyzing metabolic pathway completeness in a metagenomic dataset using KEGG ortholog (KO) annotations.

### Quick Analysis

```python
from bennu.operators import Bennu
from scripts.pathway_completeness import analyze_pathways, PATHWAYS, get_ko_presence, calculate_pathway_completeness

b = Bennu("data/DATASET/bennu.duckdb")

# Dataset-wide pathway completeness
results = analyze_pathways("data/DATASET/bennu.duckdb")
```

### Available Pathways

| Pathway | Key Marker | Biological Significance |
|---------|------------|------------------------|
| **Wood-Ljungdahl** | CODH/ACS complex | Anaerobic CO2 fixation, acetogenesis |
| **Methanogenesis** | MCR (K00399-402) | Methane production (archaea only) |
| **Calvin Cycle** | RuBisCO + PRK | Aerobic CO2 fixation |
| **TCA Cycle** | Complete cycle | Central carbon metabolism |
| **Glycolysis** | Full pathway | Sugar catabolism |
| **Nitrogen Fixation** | NifHDK | N2 → NH3 |
| **Sulfate Reduction** | DsrAB | Anaerobic respiration |
| **NiFe Hydrogenase** | Large + small subunit | H2 metabolism |
| **CRISPR-Cas (Type I)** | Cas3 + Cascade | Adaptive immunity |
| **CRISPR-Cas (Type III)** | Cas10/Csm1 | Adaptive immunity |

### Detailed Pathway Analysis

```python
# Get KOs present in dataset
all_kos = get_ko_presence(b)

# Check specific pathway
pathway_name = "Wood-Ljungdahl (Acetyl-CoA pathway)"
steps = PATHWAYS[pathway_name]
completeness, completed, missing = calculate_pathway_completeness(all_kos, steps)

print(f"{pathway_name}: {completeness:.0f}% complete")
print(f"Present: {', '.join(completed)}")
print(f"Missing: {', '.join(missing)}")
```

### Per-Genome Analysis

```python
# Get all genomes
genomes = [row[0] for row in b.store.execute("SELECT DISTINCT bin_id FROM proteins")]

# Build genome × pathway matrix
import pandas as pd

data = []
for genome in genomes:
    row = {'genome': genome}
    genome_kos = get_ko_presence(b, genome)

    for pathway_name, steps in PATHWAYS.items():
        completeness, _, _ = calculate_pathway_completeness(genome_kos, steps)
        row[pathway_name] = completeness

    data.append(row)

df = pd.DataFrame(data).set_index('genome')
print(df.to_string())
```

### Finding Pathway Genes

```python
# Find all proteins contributing to a pathway
pathway_kos = []
for step_name, kos in PATHWAYS["Wood-Ljungdahl (Acetyl-CoA pathway)"]:
    pathway_kos.extend(kos)

proteins = b.store.execute(f"""
    SELECT p.protein_id, p.bin_id, a.accession, a.name
    FROM proteins p
    JOIN annotations a ON p.protein_id = a.protein_id
    WHERE a.source = 'kegg' AND a.accession IN ({','.join(f"'{ko}'" for ko in pathway_kos)})
    ORDER BY p.bin_id, a.accession
""")

# Group by step
from collections import defaultdict
by_step = defaultdict(list)
for prot_id, genome, ko, name in proteins:
    for step_name, step_kos in PATHWAYS["Wood-Ljungdahl (Acetyl-CoA pathway)"]:
        if ko in step_kos:
            by_step[step_name].append((prot_id, genome, ko))

for step, prots in by_step.items():
    print(f"\n{step}:")
    for prot_id, genome, ko in prots[:3]:
        print(f"  {prot_id} ({genome}) - {ko}")
```

### Visualization

```python
import matplotlib.pyplot as plt
import seaborn as sns

# Heatmap of pathway completeness
fig, ax = plt.subplots(figsize=(14, 8))
sns.heatmap(
    df,
    cmap='RdYlGn',
    vmin=0, vmax=100,
    annot=True, fmt='.0f',
    ax=ax
)
ax.set_title('Pathway Completeness by Genome (%)')
plt.tight_layout()
plt.savefig("data/DATASET/exploration/figures/pathway_completeness.png", dpi=150)
```

### Syntrophic Metabolism Focus

For syntrophs (like Hinthialibacterota), key pathways to check:

**Energy conservation:**
- Group 4 hydrogenases (Mbh/Ech) - H2-evolving, energy-conserving
- Formate dehydrogenase - Alternative electron carrier
- Rnf complex - Ferredoxin:NAD+ oxidoreductase
- Electron-confurcating enzymes

**Carbon metabolism:**
- Beta-oxidation (fatty acid degradation)
- Benzoyl-CoA pathway (aromatic degradation)
- Wood-Ljungdahl (may run in reverse for CO2 reduction)

```python
# Custom syntrophic pathway check
SYNTROPHIC_MARKERS = {
    "Mbh hydrogenase": ["K18016", "K18017", "K18023"],  # Group 4
    "Ech hydrogenase": ["K15833", "K15834"],
    "Formate dehydrogenase": ["K00122", "K00123", "K00124", "K00125"],
    "Rnf complex": ["K03617", "K03618", "K03619", "K03620", "K03612", "K03613"],
    "Electron-confurcating FeFe": ["K17998", "K17999", "K18000"],
}

for marker, kos in SYNTROPHIC_MARKERS.items():
    present = [ko for ko in kos if ko in all_kos]
    print(f"{marker}: {len(present)}/{len(kos)} KOs present")
```

### Biological Interpretation

**Complete pathway (>90%)**
- Organism likely performs this function
- Check for gene clustering (operon structure)

**Partial pathway (50-90%)**
- May be functional with alternative enzymes
- Check for non-canonical enzymes
- May be degraded/non-functional

**Minimal pathway (<50%)**
- Unlikely to be functional
- Genes may have other roles
- Could be recent acquisition (incomplete transfer)

### Output

Report includes:
1. Dataset-wide completeness summary
2. Per-genome completeness matrix
3. Missing enzyme identification
4. Heatmap visualization
5. Gene lists for each pathway step
