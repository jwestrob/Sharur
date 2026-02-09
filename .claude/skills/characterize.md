# Characterize Skill

Deep characterization pipeline for unknown/ambiguous proteins. Combines multiple analysis methods to generate functional hypotheses.

**CRITICAL: You are a leaf agent. DO NOT spawn sub-agents or use the Task tool.**

**CONCURRENCY: DuckDB does not support concurrent writes. Only ONE agent should access a database at a time. The coordinator must run DB-accessing skills sequentially, not in parallel.**

> **Mandatory:** Follow the shared validation protocols in `_validation_protocols.md`.
> Verify accession names before reporting. Use COUNT(DISTINCT protein_id) for protein
> counts. Apply Context-First protocol for annotations averaging >10 hits/genome.

---

## Usage

```
/characterize protein_12345
/characterize --query "giant AND unannotated" --limit 10
/characterize --proteins prot_1,prot_2,prot_3
```

---

## Pipeline Overview

```
┌─────────────────┐
│  Input Protein  │
└────────┬────────┘
         │
    ┌────▼────┐
    │ Context │ ← Genomic neighborhood, co-occurring domains
    └────┬────┘
         │
    ┌────▼────┐
    │Structure│ ← ESM3 prediction (if <1024 aa)
    └────┬────┘
         │
    ┌────▼────┐
    │Foldseek │ ← Search PDB, AlphaFold DB
    └────┬────┘
         │
    ┌────▼────┐
    │Similar  │ ← ESM2 embedding similarity
    └────┬────┘
         │
    ┌────▼─────┐
    │Literature│ ← Search for domain/homolog info
    └────┬─────┘
         │
    ┌────▼─────┐
    │Hypothesis│ ← Synthesize evidence into prediction
    └──────────┘
```

---

## Prompt

You are characterizing unknown or poorly annotated proteins from a metagenomic dataset. Your goal is to generate functional hypotheses with supporting evidence.

### Step 1: Gather Context

```python
from sharur.operators import Sharur
b = Sharur("data/DATASET/sharur.duckdb")

# Get protein details
protein_id = "TARGET_PROTEIN"
print(b.get_protein(protein_id, verbosity=2))

# Get genomic neighborhood
print(b.get_neighborhood(protein_id, window=15))

# Check for domain hits (even weak ones)
domains = b.store.execute(f"""
    SELECT source, name, accession, evalue, description
    FROM annotations
    WHERE protein_id = '{protein_id}'
    ORDER BY evalue
""")

# Find similar proteins by embedding
similar = b.find_similar(protein_id, k=20)
```

### Step 2: Predict Structure (if appropriate)

**Criteria for structure prediction:**
- Protein is ≤1024 aa (ESM3 open model limit)
- Protein is genuinely unannotated or ambiguous
- Structure might reveal function (not just repeats)

```python
import os
from pathlib import Path

# Check for ESM API key
if not os.environ.get("ESM_API_KEY"):
    print("⚠️  ESM_API_KEY not set - skipping structure prediction")
else:
    structures_dir = Path("data/DATASET/structures")
    structures_dir.mkdir(exist_ok=True)

    pdb_path = structures_dir / f"{protein_id}.pdb"

    if not pdb_path.exists():
        result = b.predict_structure(protein_id, output_path=str(pdb_path))
        print(f"Structure predicted: pLDDT={result._raw.get('plddt_mean', 'N/A'):.2f}")
    else:
        print(f"Structure already exists: {pdb_path}")
```

### Step 3: Search Structural Homologs

```python
# Search Foldseek if structure exists
if pdb_path.exists():
    hits = b.search_foldseek(str(pdb_path), databases=["pdb100", "afdb50"], top_k=10)

    if hits.data:
        print("Structural homologs found:")
        for hit in hits.data[:5]:
            print(f"  {hit['target']}: E={hit['evalue']:.2e}, {hit.get('description', 'N/A')[:60]}")
```

### Step 4: Research Hits

For each significant Foldseek hit, look up the actual function:

```python
# For PDB hits - look up on RCSB
# WebFetch("https://www.rcsb.org/structure/XXXX", "What is this protein's function?")

# For AlphaFold hits - check UniProt
# WebFetch("https://www.uniprot.org/uniprotkb/ACCESSION", "What is the protein function?")
```

### Step 5: Synthesize Hypothesis

Combine all evidence into a functional hypothesis:

```markdown
## Characterization: {protein_id}

### Basic Info
- Length: XXX aa
- Genome: GENOME_ID
- Location: CONTIG:START-END (strand)

### Domain Evidence
- [List any domain hits, even weak ones]

### Genomic Context
- Neighboring genes suggest: [operon context]
- Co-occurs with: [related functions]

### Structural Evidence
- pLDDT: X.XX (confidence)
- Best Foldseek hit: TARGET (E=X.XX)
- Hit function: [researched function]

### Similar Proteins
- N similar proteins in dataset
- Most are annotated as: [common annotation]

### Hypothesis
**Predicted function:** [Your prediction]

**Confidence:** High/Medium/Low

**Supporting evidence:**
1. [Evidence point 1]
2. [Evidence point 2]
3. [Evidence point 3]

**Alternative interpretations:**
- [What else could it be?]

**Suggested validation:**
- [How to confirm this prediction]
```

### Step 6: Log Finding

```python
import json
from datetime import datetime
from pathlib import Path

EXPLORE_DIR = Path("data/DATASET/exploration")
EXPLORE_DIR.mkdir(exist_ok=True)

finding = {
    "timestamp": datetime.now().isoformat(),
    "category": "protein_characterization",
    "title": f"Characterized {protein_id}: [predicted function]",
    "description": "[Full hypothesis text]",
    "proteins": [protein_id],
    "evidence": {
        "length": XXX,
        "domains": [...],
        "foldseek_hits": [...],
        "similar_proteins": [...],
        "confidence": "medium",
    },
    "provenance": {
        "query": "b.search_foldseek('structures/protein_id.pdb', databases=['pdb100', 'afdb50'])",
        "raw_result": [{"target": "3GW6", "evalue": 1.2e-15, "description": "Tail fiber protein"}],
        "accession_verified": "N/A (structure-based, not accession-based)",
        "interpretation": "Structural homology to phage tail fiber (TM-score 0.72)"
    },
    "priority": "high"
}

with open(EXPLORE_DIR / "findings.jsonl", "a") as f:
    f.write(json.dumps(finding) + "\n")

# Log the analytical provenance chain
e1 = b.log_provenance(f"Domain search for {protein_id}", f"{len(domains)} domain hits")
e2 = b.log_provenance(f"Foldseek search for {protein_id}", "Top hit: 3GW6 tail fiber (E=1.2e-15)", parent_ids=[e1.entry_id])
e3 = b.log_provenance(f"Neighborhood analysis for {protein_id}", "Co-located with phage structural genes", parent_ids=[e1.entry_id])

# Register as a persistent hypothesis if the characterization yields a functional prediction
h = b.propose_hypothesis(f"{protein_id} functions as [predicted function]")
b.add_evidence(h.hypothesis_id, "Foldseek structural homology", "TM-score 0.72 to tail fiber", True, 0.7)
b.add_evidence(h.hypothesis_id, "Genomic context", "Adjacent to phage tail genes", True, 0.6)

# Review hypothesis state
print(b.hypothesis_summary())

# Render provenance DAG
b.render_provenance(title=f"Characterization: {protein_id}", output_path=f"figures/{protein_id}_provenance.mermaid")
```

---

## Batch Mode (Structural Analysis)

For batch structural analysis, first select candidates from prior exploration/survey findings, then run structure prediction and Foldseek in sequence.

### Read Prior Findings

```python
from pathlib import Path
import json

# Read exploration/survey findings to select candidates
findings_path = Path("data/YOUR_DATASET/exploration/findings.jsonl")

candidates = []
with open(findings_path) as f:
    for line in f:
        finding = json.loads(line)
        if finding.get('category') in ['novel_cluster', 'giant_protein', 'novel_feature',
                                         'giant_unannotated', 'unknown_cluster']:
            proteins = finding.get('proteins', [])
            candidates.extend(proteins)

# Deduplicate
candidates = list(dict.fromkeys(candidates))
print(f"Candidates for structure prediction: {len(candidates)}")
```

### Batch Structure Prediction

```python
structures_dir = Path("data/YOUR_DATASET/structures")
structures_dir.mkdir(exist_ok=True)

MAX_PREDICTIONS = 30
for i, protein_id in enumerate(candidates[:MAX_PREDICTIONS]):
    print(f"[{i+1}/{min(len(candidates), MAX_PREDICTIONS)}] Predicting {protein_id}...")

    output_path = structures_dir / f"{protein_id}.pdb"
    if output_path.exists():
        print(f"  Skipping (already exists)")
        continue

    result = b.predict_structure(protein_id, output_path=str(output_path))
    if result.rows > 0:
        data = result.data
        print(f"  pLDDT: {data.get('plddt', 'N/A'):.2f}")
    else:
        print(f"  Failed: {result.data}")
```

### Batch Foldseek Search

```python
import pandas as pd

results = []
for pdb_file in structures_dir.glob("*.pdb"):
    protein_id = pdb_file.stem
    print(f"Searching {protein_id}...")

    hits = b.search_foldseek(
        str(pdb_file),
        databases=["afdb50", "pdb100", "afdb-swissprot"],
        top_k=5
    )

    if hits.data and isinstance(hits.data, list):
        for hit in hits.data[:3]:
            results.append({
                "protein_id": protein_id,
                "target": hit.get("target"),
                "database": hit.get("database"),
                "evalue": hit.get("evalue"),
                "seq_identity": hit.get("seq_identity"),
                "description": hit.get("description", "")
            })
            print(f"  {hit.get('database')}: {hit.get('target')} (E={hit.get('evalue'):.2e})")
    else:
        print(f"  No significant hits")

# Save results
df = pd.DataFrame(results)
df.to_csv("data/YOUR_DATASET/foldseek_results.tsv", sep="\t", index=False)
print(f"\nStructural homologs found for {len(df['protein_id'].unique())} proteins")
```

### Results Documentation

Output files:
- `structures/*.pdb` — ESM3 predicted structures
- `foldseek_results.tsv` — All Foldseek hits in tabular format
- `exploration/findings.jsonl` — Documented structural insights

Create a summary document (`STRUCTURE_PREDICTION_FINDINGS.md`) containing:
1. Executive summary of key discoveries
2. Table of all predicted structures (Gene, Length, pLDDT, Status)
3. Detailed findings for each protein with Foldseek hits
4. Summary table (Gene, Function Prediction, Confidence, Key Evidence)

### Figure Legend Guidance

When including structure-based findings in reports:

```python
{
    "title": "Gene 584 is a putative RdRp",
    "evidence": {
        "foldseek_hits": [
            {"target": "AF-A0A6C0LPY5", "evalue": "8.6e-06", "identity": "14.3%"}
        ],
        "figure_legend": {
            "title": "RNA-dependent RNA Polymerase (Gene 584)",
            "caption": "Structure prediction for Gene 584 (620 aa) showing RdRp-like fold. "
                       "Foldseek identified multiple RdRp hits from viral metagenomes."
        }
    }
}
```

---

## Confidence Levels

**High confidence** requires 2+ of:
- Strong Foldseek hit (E < 1e-5) with known function
- Domain hit consistent with Foldseek
- Genomic context supports function
- Similar proteins have consistent annotation

**Medium confidence** requires 1 of:
- Single strong Foldseek hit
- Weak domain hit + supporting context
- Consistent annotation in similar proteins

**Low confidence:**
- Only weak evidence
- Conflicting signals
- Novel fold with no homologs

---

## Output

Results are saved to:
- `structures/{protein_id}.pdb` - Predicted structure
- `exploration/findings.jsonl` - Characterization finding
- Manifest is auto-updated with structure prediction
