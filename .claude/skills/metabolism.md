# Metabolism Skill

Analyze energy and carbon metabolism: electron transport, hydrogenases, carbon fixation, fermentation.

**CRITICAL: You are a leaf agent. DO NOT spawn sub-agents or use the Task tool.**

**CONCURRENCY: DuckDB does not support concurrent writes. Only ONE agent should access a database at a time. The coordinator must run DB-accessing skills sequentially, not in parallel.**

> **Mandatory:** Follow the shared validation protocols in `_validation_protocols.md`.
> Verify accession names before reporting. Use COUNT(DISTINCT protein_id) for protein
> counts. Apply Context-First protocol for annotations averaging >10 hits/genome.

---

## Usage

```
/metabolism                     # Full metabolic overview
/metabolism --energy            # Focus on energy conservation
/metabolism --carbon            # Focus on carbon metabolism
/metabolism --genome GENOME_ID  # Single genome analysis
```

---

## Prompt

You are analyzing metabolic capabilities in a metagenomic dataset. Focus on energy conservation, carbon metabolism, and lifestyle inferences.

### Step 1: Overview

```python
from sharur.operators import Sharur
b = Sharur("data/DATASET/sharur.duckdb")

# Key metabolic predicates
metabolic_predicates = [
    # Energy
    'hydrogenase', 'hydrogenase_group3', 'hydrogenase_group4', 'mbh_hydrogenase',
    'nife_hydrogenase', 'fefe_hydrogenase', 'electron_transport', 'atp_synthesis',
    'cytochrome', 'ferredoxin', 'flavodoxin',

    # Carbon fixation
    'wood_ljungdahl', 'carbon_fixation', 'rubisco', 'calvin_cycle',
    'reverse_tca', 'reductive_acetyl_coa',

    # Central carbon
    'glycolysis', 'gluconeogenesis', 'tca_cycle', 'pentose_phosphate',

    # Fermentation
    'fermentation', 'alcohol_dehydrogenase', 'lactate_dehydrogenase',
    'acetate_kinase', 'pyruvate_formate_lyase',

    # Respiration
    'aerobic_respiration', 'anaerobic_respiration',
    'sulfate_reduction', 'nitrate_reduction', 'methanogenesis',

    # Other
    'nitrogen_fixation', 'sulfur_metabolism', 'oxidoreductase',
]

print("Metabolic Predicate Inventory:")
for pred in metabolic_predicates:
    count = b.store.execute(f"""
        SELECT COUNT(*) FROM protein_predicates
        WHERE '{pred}' = ANY(predicates)
    """)[0][0]
    if count > 0:
        print(f"  {pred}: {count}")
```

### Step 2: Hydrogenase Analysis

Hydrogenases are key for understanding energy metabolism, especially in anaerobes:

```python
# Hydrogenase classification
hydrogenase_types = {
    'Group 1 (uptake)': ['hydrogenase_uptake', 'K00436'],
    'Group 2 (sensory)': ['hydrogenase_sensory'],
    'Group 3 (bidirectional/F420)': ['hydrogenase_group3', 'K00440', 'K00441'],
    'Group 4 (energy-conserving)': ['hydrogenase_group4', 'mbh_hydrogenase', 'ech_hydrogenase'],
    'FeFe-hydrogenase': ['fefe_hydrogenase', 'K00532'],
}

print("\nHydrogenase Classification:")
for group, markers in hydrogenase_types.items():
    # Check predicates
    pred_counts = []
    for marker in markers:
        if marker.startswith('K'):
            count = b.store.execute(f"""
                SELECT COUNT(DISTINCT protein_id) FROM annotations
                WHERE accession = '{marker}'
            """)[0][0]
        else:
            count = b.store.execute(f"""
                SELECT COUNT(*) FROM protein_predicates
                WHERE '{marker}' = ANY(predicates)
            """)[0][0]
        if count > 0:
            pred_counts.append(f"{marker}:{count}")

    if pred_counts:
        print(f"  {group}: {', '.join(pred_counts)}")
```

### Step 3: Carbon Fixation Pathways

```python
from scripts.pathway_completeness import PATHWAYS, get_ko_presence, calculate_pathway_completeness

all_kos = get_ko_presence(b)

carbon_pathways = [
    "Wood-Ljungdahl (Acetyl-CoA pathway)",
    "Calvin Cycle (CO2 fixation)",
    "TCA Cycle",  # Reverse TCA is a carbon fixation pathway
]

print("\nCarbon Fixation/Central Metabolism:")
for pathway_name in carbon_pathways:
    if pathway_name in PATHWAYS:
        completeness, completed, missing = calculate_pathway_completeness(all_kos, PATHWAYS[pathway_name])
        status = "✓" if completeness > 80 else "◐" if completeness > 50 else "○"
        print(f"  {status} {pathway_name}: {completeness:.0f}%")
        if missing and completeness < 100:
            print(f"      Missing: {', '.join(missing[:3])}")
```

### Step 4: Electron Transport Chain

```python
# ETC components
etc_components = [
    ('Complex I (NADH dehydrogenase)', 'K00330', 'nuo'),
    ('Complex II (Succinate DH)', 'K00234', 'sdh'),
    ('Complex III (Cytochrome bc1)', 'K00412', 'cyt'),
    ('Complex IV (Cytochrome oxidase)', 'K02274', 'cox'),
    ('ATP synthase', 'K02117', 'atp'),
    ('Alternative oxidase', 'K17893', None),
]

print("\nElectron Transport Chain:")
for name, ko, pfam_hint in etc_components:
    count = b.store.execute(f"""
        SELECT COUNT(DISTINCT protein_id) FROM annotations
        WHERE accession = '{ko}'
    """)[0][0]
    status = "✓" if count > 0 else "○"
    print(f"  {status} {name}: {count} proteins")
```

### Step 5: Terminal Electron Acceptors

```python
# What can this organism respire?
tea_markers = {
    'Oxygen (aerobic)': ['K02274', 'K02275'],  # Cytochrome c oxidase
    'Nitrate': ['K00370', 'K00371', 'K00374'],  # Nitrate reductase
    'Sulfate': ['K00394', 'K00395'],  # APS reductase
    'Sulfur': ['K17218', 'K17219'],  # Sulfur reductase
    'Fumarate': ['K00244', 'K00245', 'K00246'],  # Fumarate reductase
    'DMSO/TMAO': ['K07306', 'K07307', 'K07308'],
    'Iron': ['K22620'],  # Ferric reductase
}

print("\nTerminal Electron Acceptors:")
for tea, kos in tea_markers.items():
    present = sum(1 for ko in kos if b.store.execute(
        f"SELECT COUNT(*) FROM annotations WHERE accession = '{ko}'"
    )[0][0] > 0)
    if present > 0:
        print(f"  ✓ {tea}: {present}/{len(kos)} markers")
```

### Step 6: Syntrophic Metabolism Markers

For syntrophs and other obligate partners:

```python
# Syntrophic lifestyle markers
syntrophic_markers = {
    'H2-evolving hydrogenase (Mbh)': ['mbh_hydrogenase', 'hydrogenase_group4'],
    'Formate dehydrogenase': ['K00122', 'K00123', 'K00124'],
    'Electron-confurcating enzymes': ['K00534'],  # NfnAB
    'Rnf complex': ['K03617', 'K03618', 'K03619'],
    'Electron transfer flavoprotein': ['K00249', 'K00250'],
}

print("\nSyntrophic Metabolism Markers:")
for marker_name, markers in syntrophic_markers.items():
    count = 0
    for m in markers:
        if m.startswith('K'):
            count += b.store.execute(
                f"SELECT COUNT(DISTINCT protein_id) FROM annotations WHERE accession = '{m}'"
            )[0][0]
        else:
            count += b.store.execute(
                f"SELECT COUNT(*) FROM protein_predicates WHERE '{m}' = ANY(predicates)"
            )[0][0]
    if count > 0:
        print(f"  ✓ {marker_name}: {count}")
```

### Step 7: Per-Genome Metabolic Profile

```python
import pandas as pd

genomes = [row[0] for row in b.store.execute("SELECT DISTINCT bin_id FROM proteins")]

key_features = ['hydrogenase', 'wood_ljungdahl', 'electron_transport',
                'sulfate_reduction', 'nitrogen_fixation', 'crispr_associated']

profiles = []
for genome in genomes:
    row = {'genome': genome}
    for feat in key_features:
        count = b.store.execute(f"""
            SELECT COUNT(*) FROM proteins p
            JOIN protein_predicates pp ON p.protein_id = pp.protein_id
            WHERE p.bin_id = '{genome}' AND '{feat}' = ANY(pp.predicates)
        """)[0][0]
        row[feat] = count
    profiles.append(row)

df = pd.DataFrame(profiles).set_index('genome')
print("\nPer-Genome Metabolic Features:")
print(df.to_string())
```

### Step 8: Visualize Key Metabolic Loci

```python
from pathlib import Path

FIGURES_DIR = Path("data/DATASET/exploration/figures")
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Find and visualize hydrogenase operons
hydrogenase_proteins = b.search_by_predicates(has=["hydrogenase"], limit=10)
if hydrogenase_proteins.data:
    b.visualize_neighborhood(
        hydrogenase_proteins.data[0],
        window=15,
        output_path=str(FIGURES_DIR / "hydrogenase_operon.png")
    )
```

### Output Format

```markdown
## Metabolic Analysis

### Energy Conservation
- **Primary strategy:** [aerobic/anaerobic/fermentative/syntrophic]
- **Hydrogenases:** [Group types present and their roles]
- **Electron transport:** [Complete/partial chain, complexes present]
- **Terminal acceptors:** [O2/NO3/SO4/etc.]

### Carbon Metabolism
- **Autotrophy:** [WL/Calvin/rTCA pathway presence]
- **Heterotrophy:** [Substrates that can be utilized]
- **Central metabolism:** [Glycolysis/TCA completeness]

### Lifestyle Inference
Based on metabolic gene content:
- [Hypothesis about organism lifestyle]
- [Environmental niche prediction]

### Notable Findings
- [Unusual metabolic features]
- [Missing expected pathways]
- [Interesting gene clusters]

### Figures
- hydrogenase_operon.png: Representative hydrogenase gene cluster
- pathway_completeness.png: Heatmap of pathway presence
```

### Hypothesis Tracking & Provenance

Register metabolic lifestyle hypotheses and log the analytical steps that led to them:

```python
# Log key analytical steps with provenance chaining
e1 = b.log_provenance("ETC complex inventory", "Complex I: 41, II: 38, III: 3, IV: 2, V: 41")
e2 = b.log_provenance("Terminal acceptor scan", "Only 2/41 genomes have cytochrome oxidase", parent_ids=[e1.entry_id])
e3 = b.log_provenance("Hydrogenase subtyping", "39/41 have Group 4 Mbh energy-conserving", parent_ids=[e1.entry_id])

# Propose hypothesis based on accumulated evidence
h = b.propose_hypothesis("Lineage is obligately syntrophic based on incomplete electron transport chain")

# Link evidence to hypothesis
b.add_evidence(h.hypothesis_id, "ETC analysis", "Missing Complex III and IV in 38/41 genomes", True, 0.9)
b.add_evidence(h.hypothesis_id, "Hydrogenase profile", "Universal Group 4 Mbh = H2-evolving", True, 0.85)

# Review all hypotheses
print(b.hypothesis_summary())
print(b.list_hypotheses())

# Render provenance DAG for the paper
b.render_provenance(title="Metabolic Analysis", output_path="figures/metabolism_provenance.mermaid")
```

Hypotheses persist across sessions — `b.resume()` shows active hypotheses automatically.
