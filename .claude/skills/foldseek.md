# Structure Prediction & Foldseek Skill

Structure-based functional annotation for proteins where sequence annotation fails. Combines ESM3 structure prediction with Foldseek structural homology search.

**When to use:**
- Unknown/ambiguous proteins with zero or few PFAM/KEGG hits
- Giant proteins (>1000 aa) where PFAM bitscore cutoffs miss real domains
- Proteins where fold might reveal function (DUFs, orphan proteins)
- Validating weak sequence-based annotations with structural evidence

**CONCURRENCY:** DuckDB reads may run in parallel with `read_only=True`. Structure/model
work must also obey the project-wide one-MPS-process rule.

> **Mandatory:** Follow the shared validation protocols in `_validation_protocols.md`.

---

## Usage

```
/foldseek protein_12345                       # Predict + search single protein
/foldseek --batch prot_1,prot_2,prot_3       # Batch predict + search
/foldseek --giant protein_67890              # Giant protein recovery workflow
/foldseek --search structures/existing.pdb   # Search existing PDB only
```

---

## Setup

```python
from sharur.operators import Sharur
b = Sharur("data/DATASET/sharur.duckdb", read_only=True)

# Requires ESM_API_KEY for structure prediction
# export ESM_API_KEY=your_key
```

---

## Structure Prediction (ESM3)

### For database proteins

```python
result = b.predict_structure("protein_id", output_path="structures/protein_id.pdb")
# Auto-cleans sequence (removes stop codons, non-standard AAs)
# Auto-updates manifest
```

### For arbitrary sequences (no DB lookup)

```python
from sharur.operators.structure import predict_structure_from_sequence
result = predict_structure_from_sequence(
    sequence="MKVL...",
    output_path="structures/custom.pdb",
    name="my_protein"
)
```

### Batch prediction

```python
result = b.batch_predict_structures(
    protein_ids=["id1", "id2", "id3"],
    output_dir="structures/",
    max_length=1024
)
# Filenames: {protein_id}.pdb
# Skips proteins > max_length (reports them as skipped)
```

### ESM3 limits and workarounds

- **Max 1024 aa** (ESM3 open model hard limit)
- Proteins >1024 aa return an error, not a truncated structure
- **Do NOT skip giant proteins.** Export FASTA for external folding:

```python
# For >1024 aa proteins, export for external tools
rows = b.store.execute("""
    SELECT protein_id, sequence, sequence_length
    FROM proteins WHERE protein_id = ?
""", [protein_id])

if rows[0][2] > 1024:
    # Option 1: ESMFold server (up to 2048 aa, sometimes longer)
    # https://esmatlas.com/resources?action=fold
    # Option 2: Local ESMFold on H200/A100 (no length limit with enough VRAM)
    # Option 3: AlphaFold2 via ColabFold
    with open(f"structures/{protein_id}.fasta", "w") as f:
        f.write(f">{protein_id}\n{rows[0][1]}\n")
    print(f"Exported {rows[0][2]} aa FASTA for external folding")
```

### pLDDT interpretation

| pLDDT | Meaning | Action |
|-------|---------|--------|
| >70 | Confident fold | Foldseek results reliable |
| 50-70 | Uncertain loops/regions | Foldseek core hits may still be valid |
| <50 | Likely disordered | Structure may be biologically relevant (IDR) or prediction failure |

**pTM** (predicted TM-score) >0.5 means the overall fold topology is likely correct.

---

## Foldseek Search

### Standard search

```python
hits = b.search_foldseek(
    "structures/protein.pdb",
    databases=["afdb50", "pdb100", "afdb-swissprot"],
    top_k=10
)
# Returns SharurResult:
#   hits.data is formatted Markdown
#   hits.raw / hits.records contain hit dictionaries
# Local hits also include qcov and tcov.
```

### Convenience: search for a database protein

```python
hits = b.search_foldseek_for_protein("protein_id")
# Looks for existing PDB in common locations
# Returns empty if no PDB exists (does NOT auto-predict)
```

### Format results for reporting

```python
from sharur.operators.foldseek import format_foldseek_hits
print(hits.data)
# Or reformat the structured payload with a different query label:
print(format_foldseek_hits(hits.raw, query_name="my_protein"))
```

### Local vs web search

- Local binary auto-detected and execution-tested; `FOLDSEEK_BINARY` overrides discovery
- Local databases may use nested prefixes or resolved symlink prefixes under `~/.foldseek/`
- Local is faster and has no rate limits
- Falls back to web API for databases not installed locally
- `b.list_foldseek_databases()` shows available web databases

### Score interpretation

| TM-score | Meaning | Function transfer? |
|----------|---------|-------------------|
| >0.7 | High confidence homology | Yes, with literature validation |
| 0.5-0.7 | Similar fold | Possible, verify with other evidence |
| 0.3-0.5 | Same fold superfamily | Fold only; function likely diverged |
| <0.3 | Not meaningful | Do not use |

**E-value** complements TM-score:
- <1e-10: Strong hit regardless of TM-score
- 1e-10 to 1e-3: Moderate; trust if TM-score also >0.5
- >1e-3: Weak; treat with extreme skepticism

**When TM-score and E-value disagree**, the more conservative interpretation wins. A low E-value with low TM-score may indicate a partial domain match rather than full-protein homology.

---

## Hit Interpretation (CRITICAL)

### NEVER guess from target names

Foldseek target strings like `AF-Q8ZZM4-F1` or `5fms_A` tell you nothing about function. The description field from Foldseek is often truncated or ambiguous. You MUST look up every hit you plan to report.

### For PDB hits (experimentally determined structures)

```python
# Parse PDB ID: first 4 characters of target
# e.g., "5fms_A" -> PDB ID = "5fms"
WebFetch("https://www.rcsb.org/structure/5fms",
         "What is this protein's function, organism, and biological context?")
```

### For AlphaFold/UniProt hits

```python
# Parse UniProt ID: strip "AF-" prefix and "-F1" suffix
# e.g., "AF-Q8ZZM4-F1" -> UniProt ID = "Q8ZZM4"
# e.g., "AF-A0A1P8B9K7-F1" -> UniProt ID = "A0A1P8B9K7"
WebFetch("https://www.uniprot.org/uniprotkb/Q8ZZM4",
         "What is this protein's function, gene name, and GO annotations?")
```

### Batch hit interpretation: dispatch `/literature`

For any analysis where you have more than 2-3 hits to interpret, dispatch the literature agent with Protocol B rather than doing individual WebFetch calls:

```
/literature foldseek --hits structures/foldseek_results.tsv
```

The literature skill (Protocol B) handles ID parsing, database lookup, TM-score interpretation, and cross-validation across multiple hits. This is the preferred workflow for publication-quality results.

### Cross-validation

Multiple Foldseek hits pointing to the same function = higher confidence. Hits pointing to different functions = the target is probably a structural fold shared across functional families (e.g., TIM barrel, Rossmann fold). In that case, report the fold, not a specific function.

---

## Giant Protein Annotation Recovery

Standard PFAM bitscore cutoffs (`--cut_ga`) are LENGTH-BIASED. Giant proteins (>1000 aa) routinely show zero hits because the gathering thresholds were calibrated on small proteins. E-values are NOT length-biased.

### When to apply

- Any protein >1000 aa with 0 PFAM hits
- Any protein >2000 aa with <3 PFAM hits
- Before reporting ANY giant protein as "unannotated"

### Recovery workflow

```bash
# 1. Extract the protein sequence
python3 -c "
from sharur.operators import Sharur
b = Sharur('data/DATASET/sharur.duckdb', read_only=True)
rows = b.store.execute('''
    SELECT protein_id, sequence FROM proteins WHERE protein_id = ?
''', ['PROTEIN_ID'])
row = rows[0]
with open('/tmp/giant.faa', 'w') as f:
    f.write(f'>{row[0]}\n{row[1]}\n')
"

# 2. Run hmmsearch with relaxed E-value threshold
hmmsearch --domE 1e-5 --noali --domtblout /tmp/giant_domains.tsv \
    ~/.config/Astra/PFAM/Pfam-A.hmm /tmp/giant.faa

# 3. Parse results — look for repeated domains spanning the sequence
```

### E-value interpretation for giant proteins

| E-value | Confidence | Notes |
|---------|-----------|-------|
| <1e-10 | High | Real domain, report confidently |
| 1e-10 to 1e-5 | Moderate | Meaningful for giants; check for repeat patterns |
| 1e-5 to 1e-3 | Weak | Only report with structural or contextual support |

### Common giant protein architectures

| Architecture | Domains | Typical function |
|-------------|---------|-----------------|
| TPR/HEAT/WD40 repeats | PF00515, PF02985, PF00400 | Scaffold/signaling platforms |
| Big_13/Big_3_3/Big_8 | PF17956, PF17957, PF17958 | Bacterial Ig-like adhesins |
| Beta_helix (RHS/YD-repeat) | PF05593 | Toxin delivery / autotransporters |
| ANK repeats | PF00023 | Protein-protein interaction |
| Cadherin/FN3 | PF00028/PF00041 | Cell-cell adhesion |
| VCBS/FG-GAP | PF13517/PF01839 | Repeat-rich surface proteins |

### Worked example: TPR solenoid recovery

A 3572 aa protein in the susan_genomes dataset had zero standard PFAM hits. Running `hmmsearch --domE 1e-5` revealed 18+ TPR repeats spanning the entire sequence, identifying it as a solenoid scaffold protein. Without the E-value recovery, this protein would have been reported as "unannotated."

---

## Complete Workflow: Unknown Protein Characterization

```python
from sharur.operators import Sharur
b = Sharur("data/DATASET/sharur.duckdb", read_only=True)

protein_id = "TARGET_PROTEIN"

# 1. Check what annotations already exist
annots = b.store.execute("""
    SELECT source, accession, name, score, evalue
    FROM annotations WHERE protein_id = ?
    ORDER BY evalue
""", [protein_id])

# 2. Check protein size
info = b.store.execute("""
    SELECT sequence_length, bin_id FROM proteins WHERE protein_id = ?
""", [protein_id])
length = info[0][0]

# 3. Get genomic neighborhood (context-first)
print(b.get_neighborhood(protein_id, window=8, all_annotations=True))

# 4. Predict structure (if <=1024 aa)
if length <= 1024:
    result = b.predict_structure(protein_id,
        output_path=f"structures/{protein_id}.pdb")
    # Check pLDDT before proceeding to Foldseek
    plddt = result.raw.get("plddt_mean")
    if plddt and plddt < 50:
        print(f"WARNING: Low pLDDT ({plddt:.1f}) — structure may be unreliable")
else:
    # Export for external folding
    seq = b.store.execute("SELECT sequence FROM proteins WHERE protein_id = ?",
        [protein_id])[0][0]
    with open(f"structures/{protein_id}.fasta", "w") as f:
        f.write(f">{protein_id}\n{seq}\n")
    print(f"Exported {length} aa for external folding (exceeds 1024 aa ESM3 limit)")

# 5. Run Foldseek (if PDB exists)
pdb_path = f"structures/{protein_id}.pdb"
hits = b.search_foldseek(pdb_path,
    databases=["afdb50", "pdb100", "afdb-swissprot"], top_k=10)

# 6. MANDATORY: Look up every hit you plan to report
# See "Hit Interpretation" section above
# For >2-3 hits, dispatch /literature agent with Protocol B

# 7. If giant (>1000 aa) and few/no sequence hits, run E-value recovery
# See "Giant Protein Annotation Recovery" section above

# 8. Log finding with structural provenance
```

---

## Output Specification

### Files produced

| File | Content |
|------|---------|
| `structures/{protein_id}.pdb` | Predicted structure |
| `structures/{protein_id}.fasta` | Exported FASTA for >1024 aa proteins |
| `structures/foldseek_results.tsv` | Tabular Foldseek results (optional, for batch) |
| `exploration/findings.jsonl` | Append findings with structural evidence |

### Finding schema for structural results

```json
{
  "id": "E015",
  "title": "DUF1015 is a SerK-family phosphotransferase (Foldseek TM=0.82)",
  "category": "structural_characterization",
  "description": "Foldseek search of the predicted ESM3 structure for protein X (DUF1015, 312 aa) returned PDB 5X0K (SerK phosphotransferase, E=1e-8, TM-score 0.82). The kinase fold and active-site conservation are consistent with phosphotransferase function. Genomic context shows co-localization with sugar transport genes.",
  "evidence": "ESM3 pLDDT=78.3; Foldseek top hit PDB 5X0K (TM=0.82, E=1e-8); 3/5 catalytic residues conserved",
  "n_genomes": 15,
  "provenance": {
    "query": "b.predict_structure('protein_X') -> b.search_foldseek('structures/protein_X.pdb')",
    "raw_result": "PDB 5X0K: SerK kinase, Thermotoga maritima, 1.8A resolution",
    "accession_verified": "WebFetch('https://www.rcsb.org/structure/5X0K', ...) confirmed phosphotransferase",
    "interpretation": "Structural homology (TM>0.7) to characterized kinase supports phosphotransferase function"
  },
  "figures": ["figures/protein_X_neighborhood.png"],
  "phase": "exploration"
}
```

**Every structural claim MUST include:**
1. pLDDT of the predicted structure
2. TM-score and E-value of the Foldseek hit
3. Provenance showing the hit was looked up (WebFetch or /literature), not guessed
4. Genomic context (neighborhood check)

---

## Common Pitfalls

1. **Guessing PDB function from name.** Target names like `5fms_A` or `AF-Q8ZZM4-F1` are identifiers, not functional descriptions. Always look them up.

2. **Skipping giant proteins.** "No PFAM hits" for a >1000 aa protein almost always means the bitscore cutoff failed, not that the protein has no domains. Run E-value recovery.

3. **Reporting low-confidence structures.** If pLDDT < 50, the structure is likely wrong. Don't run Foldseek on garbage structures.

4. **Conflating fold with function.** A TIM barrel (TM=0.6) could be any of dozens of enzyme families. Report the fold, not a specific function, unless TM > 0.7 AND the active site is conserved.

5. **Ignoring partial matches.** Foldseek may match only a subdomain (qstart-qend covers 30% of query). This means the protein has additional domains not captured by the hit. Report the match region explicitly.

6. **Not checking neighborhood.** A Foldseek hit tells you the fold. The genomic neighborhood tells you the function. A protein structurally similar to a toxin that sits in a biosynthetic operon is probably a biosynthetic enzyme, not a toxin.

---

## When to Dispatch Sub-Agents

| Situation | Dispatch |
|-----------|----------|
| >2-3 Foldseek hits need functional lookup | `/literature` (Protocol B) |
| Hit is a DUF that might be reclassified | `/literature` (Protocol A) |
| Hit is in a defense locus | `/defense` for system-level context |
| Protein is in a prophage region | `/prophage` for viral context |
| Giant protein needs full characterization | `/characterize` with structural evidence |

Return validated functional assignments to the calling agent. Include TM-score, E-value, provenance URL, and confidence level for every assignment.
