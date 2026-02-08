# Defense Systems Skill

Analyze defense systems: CRISPR-Cas, restriction-modification, toxin-antitoxin, and other anti-phage mechanisms.

**CRITICAL: You are a leaf agent. DO NOT spawn sub-agents or use the Task tool.**

**CONCURRENCY: DuckDB does not support concurrent writes. Only ONE agent should access a database at a time. The coordinator must run DB-accessing skills sequentially, not in parallel.**

> **Mandatory:** Follow the shared validation protocols in `_validation_protocols.md`.
> Verify accession names before reporting. Use COUNT(DISTINCT protein_id) for protein
> counts. Apply Context-First protocol for annotations averaging >10 hits/genome.

---

## Usage

```
/defense                    # Full defense inventory
/defense --crispr           # CRISPR-Cas focus
/defense --genome GENOME_ID # Single genome analysis
```

---

## Prompt

You are analyzing defense systems in a metagenomic dataset. Provide a comprehensive inventory with proper classification and validation.

### Step 1: Run Built-in CRISPR Analysis

```python
from bennu.operators import Bennu
b = Bennu("data/DATASET/bennu.duckdb")

# Comprehensive CRISPR-Cas analysis
print(b.analyze_crispr_systems())
```

This returns:
- Total arrays and their types
- Complete vs incomplete systems
- Orphan arrays (no nearby cas genes)
- Arrays at contig edges (fragmentation artifacts)

### Step 2: Defense System Inventory

```python
# Count all defense-related predicates
defense_predicates = [
    'crispr_associated', 'cas_domain',
    'restriction_modification', 'restriction_enzyme', 'methyltransferase',
    'toxin_antitoxin', 'toxin_domain', 'antitoxin_domain',
    'defense_system', 'anti_crispr', 'anti_restriction',
    'abortive_infection', 'retron',
]

for pred in defense_predicates:
    count = b.store.execute(f"""
        SELECT COUNT(*) FROM protein_predicates
        WHERE '{pred}' = ANY(predicates)
    """)[0][0]
    if count > 0:
        print(f"{pred}: {count}")
```

### Step 2b: DefenseFinder / PADLOC Integration

Check for system-level defense annotations from DefenseFinder and PADLOC:

```python
# DefenseFinder annotations
df_systems = b.store.execute("""
    SELECT a.accession, a.name, a.description,
           COUNT(DISTINCT a.protein_id) as n_proteins,
           COUNT(DISTINCT p.bin_id) as n_genomes
    FROM annotations a
    JOIN proteins p ON a.protein_id = p.protein_id
    WHERE a.source = 'defensefinder'
    GROUP BY a.accession, a.name, a.description
    ORDER BY n_proteins DESC
""").fetchall()

if df_systems:
    print("DefenseFinder Systems:")
    for acc, name, desc, n_prot, n_gen in df_systems:
        print(f"  {name} ({acc}): {n_prot} proteins in {n_gen}/{n_genomes} genomes")
else:
    print("No DefenseFinder annotations found")

# PADLOC annotations (same pattern)
padloc_systems = b.store.execute("""
    SELECT a.accession, a.name, a.description,
           COUNT(DISTINCT a.protein_id) as n_proteins,
           COUNT(DISTINCT p.bin_id) as n_genomes
    FROM annotations a
    JOIN proteins p ON a.protein_id = p.protein_id
    WHERE a.source = 'padloc'
    GROUP BY a.accession, a.name, a.description
    ORDER BY n_proteins DESC
""").fetchall()

if padloc_systems:
    print("\nPADLOC Systems:")
    for acc, name, desc, n_prot, n_gen in padloc_systems:
        print(f"  {name} ({acc}): {n_prot} proteins in {n_gen}/{n_genomes} genomes")
```

### Step 2c: Cross-Reference Defense Detections

Compare DefenseFinder/PADLOC system-level calls with predicate-based detection to find discrepancies:

```python
# Predicate-based defense protein count
pred_defense = b.store.execute("""
    SELECT COUNT(DISTINCT protein_id) FROM protein_predicates
    WHERE 'defense_system' = ANY(predicates)
""").fetchone()[0]

# DefenseFinder protein count
df_defense = b.store.execute("""
    SELECT COUNT(DISTINCT protein_id) FROM annotations
    WHERE source = 'defensefinder'
""").fetchone()[0]

print(f"Defense proteins (predicates): {pred_defense}")
print(f"Defense proteins (DefenseFinder): {df_defense}")

# Find DefenseFinder proteins NOT in predicate defense set
if df_defense > 0:
    missing_from_preds = b.store.execute("""
        SELECT a.protein_id, a.name, a.accession
        FROM annotations a
        WHERE a.source = 'defensefinder'
          AND a.protein_id NOT IN (
            SELECT protein_id FROM protein_predicates
            WHERE 'defense_system' = ANY(predicates)
          )
        LIMIT 10
    """).fetchall()
    if missing_from_preds:
        print(f"\nDefenseFinder hits NOT in defense predicates ({len(missing_from_preds)} shown):")
        for pid, name, acc in missing_from_preds:
            print(f"  {pid}: {name} ({acc})")
```

### Step 3: CRISPR-Cas Classification

**Type I** (signature: Cas3 helicase-nuclease):
```python
# Find Cas3 proteins
cas3 = b.search_by_predicates(has=["cas3"])
print(f"Type I systems (Cas3): {len(cas3.data)} proteins")

# Check subtypes by Cas8 variants
for protein_id in cas3.data[:5]:
    print(b.get_neighborhood(protein_id, window=10))
```

**Type III** (signature: Cas10/Csm1):
```python
cas10 = b.store.execute("""
    SELECT p.protein_id, p.bin_id
    FROM proteins p
    JOIN annotations a ON p.protein_id = a.protein_id
    WHERE a.name ILIKE '%cas10%' OR a.name ILIKE '%csm1%' OR a.name ILIKE '%cmr2%'
""")
print(f"Type III systems (Cas10): {len(cas10)} proteins")
```

**Type V** (signature: Cas12):
```python
# Cas12 variants - watch for TnpB confusion!
# Use Context-First protocol: check neighborhood to distinguish real CRISPR loci
# from orphan transposase-like domains
cas12_candidates = b.store.execute("""
    SELECT p.protein_id, a.name
    FROM proteins p
    JOIN annotations a ON p.protein_id = a.protein_id
    WHERE a.name ILIKE '%cas12%'
""")

# Validate via neighborhood - true Cas12 should be near CRISPR arrays and other cas genes
for pid, domain in cas12_candidates:
    nbr = b.get_neighborhood(pid, window=8, all_annotations=True)
    # Check for CRISPR array, other cas genes, or repeat-spacer structures nearby
    # Isolated Cas12-like proteins without CRISPR context are likely TnpB transposases
    print(f"  {pid}: {domain} — check neighborhood for CRISPR array context")
```

### Step 4: Restriction-Modification Systems

```python
# R-M system components
rm_proteins = b.store.execute("""
    SELECT p.protein_id, p.bin_id, a.name
    FROM proteins p
    JOIN annotations a ON p.protein_id = a.protein_id
    WHERE a.name ILIKE '%restriction%'
       OR a.name ILIKE '%methyltransferase%'
       OR a.name ILIKE '%methylase%'
    ORDER BY p.bin_id, p.gene_index
""")

# Group by genome and check for paired R-M
from collections import defaultdict
by_genome = defaultdict(lambda: {'restriction': [], 'methyltransferase': []})
for pid, genome, name in rm_proteins:
    if 'restriction' in name.lower():
        by_genome[genome]['restriction'].append(pid)
    else:
        by_genome[genome]['methyltransferase'].append(pid)

for genome, counts in by_genome.items():
    print(f"{genome}: {len(counts['restriction'])} restriction, {len(counts['methyltransferase'])} MTases")
```

### Step 5: Toxin-Antitoxin Systems

```python
# TA system analysis
ta_proteins = b.search_by_predicates(has=["toxin_antitoxin"])
print(f"Confirmed TA pairs: {len(ta_proteins.data)}")

# Check for orphan toxins (dangerous!)
orphan_toxins = b.store.execute("""
    WITH toxins AS (
        SELECT p.protein_id, p.contig_id, p.gene_index
        FROM proteins p
        JOIN protein_predicates pp ON p.protein_id = pp.protein_id
        WHERE 'toxin_domain' = ANY(pp.predicates)
    ),
    antitoxins AS (
        SELECT p.contig_id, p.gene_index
        FROM proteins p
        JOIN protein_predicates pp ON p.protein_id = pp.protein_id
        WHERE 'antitoxin_domain' = ANY(pp.predicates)
    )
    SELECT t.protein_id
    FROM toxins t
    LEFT JOIN antitoxins a ON t.contig_id = a.contig_id
        AND ABS(t.gene_index - a.gene_index) <= 3
    WHERE a.gene_index IS NULL
""")
print(f"Orphan toxins (no nearby antitoxin): {len(orphan_toxins)}")
```

### Step 6: Novel/Rare Defense Systems

```python
# Check for recently discovered systems
novel_defense = [
    ('retron', 'Retron defense'),
    ('abortive_infection', 'Abortive infection'),
    ('thoeris', 'Thoeris'),
    ('gabija', 'Gabija'),
    ('lamassu', 'Lamassu'),
    ('septu', 'Septu'),
    ('hachiman', 'Hachiman'),
]

for pred, name in novel_defense:
    count = b.store.execute(f"""
        SELECT COUNT(*) FROM protein_predicates
        WHERE '{pred}' = ANY(predicates)
    """)[0][0]
    if count > 0:
        print(f"{name}: {count} proteins")
```

### Step 7: Defense Islands

Defense systems often cluster. Look for defense islands:

```python
# Find contigs with multiple defense genes
defense_rich = b.store.execute("""
    SELECT p.contig_id, COUNT(DISTINCT p.protein_id) as n_defense
    FROM proteins p
    JOIN protein_predicates pp ON p.protein_id = pp.protein_id
    WHERE 'crispr_associated' = ANY(pp.predicates)
       OR 'restriction_modification' = ANY(pp.predicates)
       OR 'toxin_antitoxin' = ANY(pp.predicates)
       OR 'defense_system' = ANY(pp.predicates)
    GROUP BY p.contig_id
    HAVING COUNT(DISTINCT p.protein_id) >= 5
    ORDER BY n_defense DESC
""")

print(f"Defense islands (≥5 defense genes): {len(defense_rich)}")
for contig, n in defense_rich[:5]:
    print(f"  {contig}: {n} defense genes")
```

### Step 8: Visualize Key Loci

```python
from pathlib import Path

FIGURES_DIR = Path("data/DATASET/exploration/figures")
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Visualize top CRISPR locus
if cas3.data:
    b.visualize_neighborhood(
        cas3.data[0],
        window=15,
        output_path=str(FIGURES_DIR / "crispr_type_i_locus.png")
    )
```

### Output Format

```markdown
## Defense Systems Inventory

### Summary
- CRISPR-Cas systems: X complete, Y incomplete, Z orphan arrays
- R-M systems: X complete pairs
- TA systems: X confirmed pairs, Y orphan toxins
- Other defense: [list]

### CRISPR-Cas Details
[Per-type breakdown with validation notes]

### Defense Islands
[List defense-rich contigs with gene inventories]

### Figures
- crispr_type_i_locus.png: Representative Type I system
- defense_island_contig_X.png: Dense defense region

### Biological Interpretation
[What does the defense repertoire suggest about viral pressure?]
```

### Step 9: Ecological Interpretation

Go beyond inventories — interpret the defense repertoire ecologically.

#### Defense Load Calculation

```python
# Calculate defense investment per genome
defense_load = b.store.execute("""
    SELECT p.bin_id,
           COUNT(DISTINCT CASE WHEN 'defense_system' = ANY(pp.predicates)
                               OR 'crispr_associated' = ANY(pp.predicates)
                               OR 'restriction_modification' = ANY(pp.predicates)
                               OR 'toxin_antitoxin' = ANY(pp.predicates)
                 THEN p.protein_id END) as n_defense,
           COUNT(DISTINCT p.protein_id) as n_total
    FROM proteins p
    LEFT JOIN protein_predicates pp ON p.protein_id = pp.protein_id
    GROUP BY p.bin_id
    ORDER BY n_defense DESC
""").fetchall()

for genome, n_defense, n_total in defense_load:
    pct = n_defense / n_total * 100 if n_total > 0 else 0
    print(f"  {genome}: {n_defense}/{n_total} ({pct:.1f}%)")
```

#### Interpretation Guide

| Defense Load | Interpretation |
|-------------|----------------|
| >10% of proteins | Heavy phage pressure — active arms race |
| 5-10% | Moderate — typical for free-living prokaryotes |
| <5% | Low — stable niche, reduced viral pressure, or obligate lifestyle |
| 0% | Check assembly quality before concluding |

#### Defense Diversity vs Depth

- **Many system types, few copies each** → diverse viral threats, broad defense
- **Few system types, many copies** → specialized defense against specific phage families
- **CRISPR-dominant** → adaptive immunity primary strategy, high spacer diversity expected
- **Innate-dominant** → constitutive defense, less phage-responsive

#### Co-occurrence Patterns

```python
# Which defense systems co-occur in the same genomes?
# Build a genome x defense-type matrix
defense_types = ['crispr_associated', 'restriction_modification',
                 'toxin_antitoxin', 'abortive_infection']

for i, dt1 in enumerate(defense_types):
    for dt2 in defense_types[i+1:]:
        both = b.store.execute(f"""
            SELECT COUNT(DISTINCT p.bin_id) FROM proteins p
            JOIN protein_predicates pp ON p.protein_id = pp.protein_id
            WHERE '{dt1}' = ANY(pp.predicates)
              AND p.bin_id IN (
                SELECT DISTINCT p2.bin_id FROM proteins p2
                JOIN protein_predicates pp2 ON p2.protein_id = pp2.protein_id
                WHERE '{dt2}' = ANY(pp2.predicates)
              )
        """).fetchone()[0]
        print(f"  {dt1} + {dt2}: {both} genomes")
```
