# Atlas Skill

Exhaustive bottom-up genome-by-genome reading of a metagenomic dataset. Dispatches subagents per genome batch (3-5 genomes each) to read through every annotation source, sample neighborhoods, and flag items for specialist follow-up.

**Atlas vs Survey vs Explore:**
- **Survey**: Top-down pass by functional category. Fast. Gets you started.
- **Atlas**: Bottom-up reading, genome by genome. Subagent per batch. Hours. Catches what survey misses.
- **Explore**: Hypothesis-driven follow-up on what survey/atlas found.

Run atlas AFTER survey to fill coverage gaps, or INSTEAD of survey when you want exhaustive per-genome inventories from the start.

**AGENT ARCHITECTURE:**
- Atlas is a **coordinator agent** that dispatches genome-batch subagents
- Genome-batch subagents are **leaf agents** that read and write findings
- Subagents CAN dispatch their own sub-subagents for specialist tasks (literature, foldseek)
- **Small datasets (<20 genomes):** subagents run in PARALLEL (read-only on DB, each writes own file)
- **Large datasets (20+):** batch subagents run SEQUENTIALLY (concurrent append contention)
- Literature/web search sub-subagents can always run in parallel (no DB writes)

**CONCURRENCY:** Subagents are read-only on DuckDB. For <20 genomes, dispatch all in parallel — each writes to its own `atlas/inventory_{genome}.jsonl`. Coordinator merges after. For 20+ genomes, batch and run sequentially with shared `atlas/inventories.jsonl`.

> **Mandatory:** Follow the shared validation protocols in `_validation_protocols.md`.

---

## Output Files

| File | Contents |
|------|----------|
| `atlas/inventories.jsonl` | One JSON line per genome: annotation census, flags, summary |
| `survey/findings.jsonl` | Notable findings (standard schema, `"phase": "atlas"`) |
| `atlas/ATLAS_SUMMARY.md` | Dataset-wide overview after all batches complete |
| `atlas/flags_collected.json` | Aggregated flags for specialist dispatch |

---

## Coordinator Workflow

```
1. Query DB for all genomes, get protein counts, sort by size
2. Split into batches of 3-5 genomes
3. For each batch (SEQUENTIALLY):
   a. Craft subagent prompt with genome list + boilerplate
   b. Dispatch subagent
   c. Verify subagent wrote inventory entries + findings
4. After all batches complete:
   a. Read all inventories, collect flags by type
   b. Dispatch specialist agents based on flags
   c. Write atlas/ATLAS_SUMMARY.md
```

### Step 1: Get Genomes

```python
from sharur.operators import Sharur
import json, os

b = Sharur("data/DATASET/sharur.duckdb", read_only=True)

genomes = b.store.execute("""
    SELECT bin_id,
           COUNT(*) as n_proteins,
           COUNT(DISTINCT contig_id) as n_contigs,
           MAX(sequence_length) as max_protein_length,
           ROUND(AVG(sequence_length), 0) as avg_protein_length
    FROM proteins
    GROUP BY bin_id
    ORDER BY n_proteins DESC
""")

print(f"Total genomes: {len(genomes)}")
for g in genomes[:5]:
    print(f"  {g[0]}: {g[1]} proteins, {g[2]} contigs")
```

### Step 2: Choose Dispatch Strategy

```python
if len(genomes) <= 20:
    # Small dataset: 1 agent per genome, run ALL in parallel
    # Each writes to its own file — no contention
    batches = [[g[0]] for g in genomes]
    parallel = True
else:
    # Large dataset: batch and run sequentially
    batch_size = min(5, max(3, len(genomes) // 10))
    batches = []
    for i in range(0, len(genomes), batch_size):
        batch = genomes[i:i + batch_size]
        batches.append([g[0] for g in batch])
    parallel = False

print(f"Strategy: {'parallel' if parallel else 'sequential'}, {len(batches)} {'agents' if parallel else 'batches'}")
```

### Step 3: Dispatch Subagents

For small datasets (parallel): dispatch ALL agents in a single message with `run_in_background=True`.
For large datasets (sequential): dispatch one batch at a time, wait for completion.

### Step 4: Collect Flags and Dispatch Specialists

```python
import json
from collections import defaultdict

# Read all inventories
inventories = []
with open("data/DATASET/atlas/inventories.jsonl") as f:
    for line in f:
        if line.strip():
            inventories.append(json.loads(line))

# Collect flags
flags = defaultdict(list)
for inv in inventories:
    for flag in inv.get("flags", []):
        # Extract flag type (e.g., "hyddb_unvalidated" from "hyddb_unvalidated:protein_id")
        flag_type = flag.split(":")[0] if ":" in flag else flag
        flags[flag_type].append({"genome": inv["genome"], "flag": flag})

# Dispatch map
dispatch_map = {
    "hyddb_unvalidated":       "/hydrogenase — neighborhood curation",
    "giant_unannotated":       "/characterize or /foldseek — structural homology",
    "high_prevalence_defense":  "/defense — with superfamily awareness warning",
    "unknown_clusters":         "/literature — functional ambiguity resolution",
    "ambiguous_annotation":     "/literature — domain vs function clarification",
    "novel_operon":             "/explore — locus characterization",
    "prophage_candidate":       "/prophage — viral element validation",
}

for flag_type, items in flags.items():
    specialist = dispatch_map.get(flag_type, "manual review")
    print(f"  {flag_type}: {len(items)} items -> {specialist}")
```

### Step 5: Write Summary

After specialists complete, write `atlas/ATLAS_SUMMARY.md` covering:
- Total genomes read, annotation coverage statistics
- Per-genome summaries (one paragraph each)
- Flags dispatched and resolution status
- Patterns visible only at per-genome resolution that survey missed

---

## Subagent Prompt Template

This is the prompt the coordinator sends to each genome-batch subagent. Copy it verbatim and fill in the bracketed values.

```
You are reading through {N} genomes exhaustively for the Atlas skill.

DB: data/{DATASET}/sharur.duckdb
Import: from sharur.operators import Sharur; b = Sharur("data/{DATASET}/sharur.duckdb", read_only=True)
Genomes in this batch: {genome_list}
Output directory: data/{DATASET}/atlas/
Draft findings file: data/{DATASET}/atlas/findings_{agent_id}.jsonl

## Database Column Reference
- Annotations table: 'name' (not annotation_id), 'score' (not bitscore)
- Proteins table: 'sequence_length' (not 'length')
- b.store.execute() returns a list — do NOT call .fetchall() or .fetchone()
- Always COUNT(DISTINCT protein_id) for protein counts — repeat domains inflate COUNT(*)
- MAG caveat: "Not detected" not "absent"

## Domain Documentation (READ ON DEMAND)
When you encounter a domain-specific situation, look up the relevant protocol doc:
Docs path: /Users/jacob/Documents/Obsidian Vault/sharur-docs/
Key docs:
  - hydrogenase-classification.md — HydDB curation, Complex I FP detection, neighborhood KEGG KOs
  - defense-system-validation.md — superfamily FP filtering, prevalence sanity checks
  - giant-protein-recovery.md — E-value recovery for >1000 aa, ESM3/Foldseek workflow
  - context-first-protocol.md — co-annotation validation, claim escalation ladder
  - mag-quality-interpretation.md — fragmentation checks, absence claim language
Use the Read tool to load any doc when you need it. Don't guess — look it up.

## Your Task

For EACH genome in your batch, perform the following steps and write one inventory
entry to the agent's unique inventory spool. Log notable findings to the agent's
unique draft findings spool with phase="atlas"; the coordinator performs a strict merge.

### Step 1: Annotation Source Census

```python
census = b.store.execute("""
    SELECT a.source, COUNT(DISTINCT a.protein_id) as n_proteins
    FROM annotations a
    JOIN proteins p ON a.protein_id = p.protein_id
    WHERE p.bin_id = '{genome}'
    GROUP BY a.source
    ORDER BY n_proteins DESC
""")

total_proteins = b.store.execute("""
    SELECT COUNT(*) FROM proteins WHERE bin_id = '{genome}'
""")[0][0]
```

Record: source -> protein count for every annotation source present.

### Step 2: Top Annotations

```python
top_annots = b.store.execute("""
    SELECT a.source, a.accession, a.name,
           COUNT(DISTINCT a.protein_id) as n
    FROM annotations a
    JOIN proteins p ON a.protein_id = p.protein_id
    WHERE p.bin_id = '{genome}'
    GROUP BY a.source, a.accession, a.name
    ORDER BY n DESC
    LIMIT 30
""")
```

### Step 3: Context-First Check (CRITICAL)

For any annotation with >10 hits in this genome, run a co-annotation check:

```python
co_annots = b.store.execute("""
    SELECT a2.source, a2.accession, a2.name,
           COUNT(DISTINCT a1.protein_id) as n
    FROM annotations a1
    JOIN annotations a2 ON a1.protein_id = a2.protein_id
    JOIN proteins p ON a1.protein_id = p.protein_id
    WHERE p.bin_id = '{genome}'
      AND a1.accession = '{high_hit_accession}'
      AND a2.accession != '{high_hit_accession}'
    GROUP BY a2.source, a2.accession, a2.name
    ORDER BY n DESC LIMIT 10
""")
```

Verdict for each:
- If top co-annotations are from a DIFFERENT enzyme family -> "superfamily_fp"
- If co-annotations support the claimed function -> "validated"
- If ambiguous -> "ambiguous_annotation" (flag for literature)

### Step 4: Unannotated Proteins

```python
unannotated = b.store.execute("""
    SELECT p.protein_id, p.sequence_length
    FROM proteins p
    LEFT JOIN annotations a ON p.protein_id = a.protein_id
    WHERE p.bin_id = '{genome}'
      AND a.protein_id IS NULL
    ORDER BY p.sequence_length DESC
""")
```

Report: total unannotated count, list of giant unannotated (>1000 aa).

### Step 5: Giant Proteins

```python
giants = b.store.execute("""
    SELECT p.protein_id, p.sequence_length,
           COUNT(DISTINCT a.accession) as n_annotations
    FROM proteins p
    LEFT JOIN annotations a ON p.protein_id = a.protein_id
    WHERE p.bin_id = '{genome}' AND p.sequence_length > 1000
    GROUP BY p.protein_id, p.sequence_length
    ORDER BY p.sequence_length DESC
""")
```

Flag unannotated giants (>1000 aa, 0 annotations) as "giant_unannotated".
Note: giant proteins (>1000 aa) with zero PFAM hits may need E-value recovery
(PFAM bitscore cutoffs are length-biased).

### Step 6: Sample Neighborhoods (5-8 per genome)

Pick proteins to sample from:
- Largest unannotated protein
- Protein with most unusual annotation (rare accession)
- A HydDB hit (if any) — check for Complex I false positives
- A DefenseFinder hit — check for superfamily inflation
- Cluster of 3+ consecutive unannotated proteins
- Any protein adjacent to a contig edge (fragmentation check)

```python
nbr = b.get_neighborhood(protein_id, window=8, all_annotations=True)
```

Record what you see. Flag anything needing specialist follow-up.

### Step 7: Write Inventory Entry

Append one JSON line per genome to atlas/inventories.jsonl:

```python
import json, os

inventory = {
    "genome": bin_id,
    "n_proteins": total_proteins,
    "n_contigs": n_contigs,
    "annotation_coverage": {source: count for source, count in census},
    "unannotated_count": len(unannotated),
    "unannotated_pct": round(len(unannotated) / total_proteins * 100, 1),
    "giant_proteins": len(giants),
    "giant_unannotated": [g[0] for g in giants if g[2] == 0],
    "neighborhoods_sampled": N_sampled,
    "flags": collected_flags,  # list of strings
    "context_first_checks": context_checks,  # list of dicts
    "summary": one_sentence_summary
}

os.makedirs("data/{DATASET}/atlas", exist_ok=True)
with open("data/{DATASET}/atlas/inventories.jsonl", "a") as f:
    f.write(json.dumps(inventory) + "\n")
```

### Step 8: Log Notable Findings

For anything genuinely notable (not routine), write to the unique draft spool:

```python
from sharur.core.analysis_record_io import append_finding_record

finding = {
    "id": "atlas-{genome}-NNN",  # e.g., atlas-mb104-001
    "title": "Self-contained title with all qualifiers",
    "category": "appropriate_category",
    "description": "Prose paragraph with biological interpretation.",
    "evidence": {"genome": bin_id, ...},
    "verification": [
        {"claim": "...", "query": "...", "expected": "..."},
    ],
    "n_genomes": 1,
    "provenance": {
        "query": "the SQL that produced this",
        "raw_result": "literal output",
        "interpretation": "what it means"
    },
    "figures": [],
    "related_findings": [],
    "phase": "atlas"
}

append_finding_record(
    "data/{DATASET}/atlas/findings_{agent_id}.jsonl",
    finding,
    phase="atlas",
    strict=False,  # draft only; canonical merge must be strict
)
```

What is "notable" — worth a finding entry:
- A genome that completely lacks a function present in >80% of others
- An unusual annotation (rare accession, <5% prevalence dataset-wide)
- A giant unannotated protein >2000 aa
- An operon-scale cluster of unannotated genes (5+ consecutive)
- A defense island (3+ defense systems co-located)
- A HydDB hit that needs neighborhood curation
- Any context-first check that overturns a functional assumption

What is NOT worth a finding entry:
- Routine annotation statistics (those go in the inventory)
- Expected annotations for the organism type
- Single unannotated proteins under 1000 aa
- Annotations already captured by survey findings

## Ambiguity Checks Without Priming

Do not give subagents a named list of expected false positives. Require them to inspect
the live schema for whatever curated callers exist, separate raw observed domains from
named caller output, and trigger co-annotation/neighborhood/specialist validation when
the evidence class is broad or unexpectedly prevalent. Record the check and verdict
without prescribing the answer in advance.

## Flag Types

Flags are strings collected from all subagents and used to dispatch specialists:

| Flag | Meaning | Specialist |
|------|---------|------------|
| `hyddb_unvalidated` | HydDB hit without PF00374 corroboration | /hydrogenase |
| `giant_unannotated` | Protein >1000 aa with zero annotation hits | /characterize or /foldseek |
| `giant_unannotated_Nk` | Protein >N000 aa unannotated (e.g., `giant_unannotated_5k`) | /foldseek (priority) |
| `high_prevalence_defense` | DefenseFinder hit at >50% genome prevalence | /defense with superfamily warning |
| `unknown_clusters` | 5+ consecutive unannotated proteins on one contig | /literature or /explore |
| `ambiguous_annotation` | Context-first check returned "ambiguous" | /literature |
| `novel_operon` | Conserved gene cluster not matching known systems | /explore |
| `prophage_candidate` | VOGdb/phage markers clustered on a contig | /prophage |
| `no_annotation_source_X` | Genome has zero hits from expected source (e.g., KEGG) | Check pipeline completeness |
| `contig_edge_giant` | Giant protein at contig edge — likely fragmented | Note in inventory, deprioritize |

Format: `flag_type` or `flag_type:protein_id` or `flag_type:details`

## Finding ID Convention

Atlas findings use the prefix `atlas-{genome_id}-NNN`:
- `atlas-mb104-001`, `atlas-mb104-002`, etc.
- This avoids collision with survey (`survey-NNN`) and explore (`ENNN`) IDs
- Cross-genome patterns found during coordinator summary use `atlas-cross-NNN`

## Coordinator: Writing ATLAS_SUMMARY.md

After all subagents complete and specialists have been dispatched, write a summary:

```markdown
# Atlas Summary: {DATASET}

## Overview
- **Genomes read:** N
- **Total proteins:** N
- **Mean annotation coverage:** N% (range: N%-N%)
- **Mean unannotated:** N% (range: N%-N%)

## Annotation Source Coverage
| Source | Genomes with hits | Mean proteins/genome | Notes |
...

## Per-Genome Summaries
(One paragraph per genome, sorted by size. Include: protein count, dominant
annotation sources, notable features, flags raised.)

## Flags Dispatched
| Flag type | Count | Specialist | Status |
...

## Patterns Not Visible in Survey
(Things that only become apparent when reading genome-by-genome:
  - Genomes with unusual annotation profiles
  - Consistent gaps across a clade
  - Annotation sources that fail for specific genomes
  - Co-occurrence patterns between rare features)

## Recommendations for Explore
(Specific hypotheses or loci that merit follow-up investigation.)
```

## Scaling Rules

| Dataset size | Agents | Parallelism | Output pattern |
|---|---|---|---|
| <20 genomes | 1 per genome | **Parallel** (read-only on DB) | `atlas/inventory_{genome}.jsonl` per agent, merge after |
| 20-100 genomes | 3-5 per batch | Sequential batches | `atlas/inventories.jsonl` (append) |
| 100-500 genomes | 5-10 per batch | Sequential batches | `atlas/inventories.jsonl` (append) |
| 500+ genomes | Sample + batch | Sequential batches | Sample representatives, not all |

**Why parallel works for small datasets:** Atlas subagents are read-only on DuckDB — they
only query, never write to the database. Each agent writes to its own JSONL file, avoiding
append contention. The coordinator merges all per-genome files after completion.

**Parallel output pattern:**
- Each agent writes: `atlas/inventory_{short_name}.jsonl` (single JSON line)
- Each agent appends findings to its own: `atlas/findings_{short_name}.jsonl`
- Coordinator merges: `cat atlas/inventory_*.jsonl > atlas/inventories.jsonl`
- Coordinator reads every draft finding, resolves duplicate IDs or conflicting claims,
  and appends accepted records to `survey/findings.jsonl` with
  `append_finding_record(..., strict=True)`. Never use `cat` for canonical findings.

**Performance:**
- Each subagent typically takes 2-5 minutes per genome
- If a subagent fails, re-dispatch it individually after the others complete

## Relationship to Other Phases

Atlas findings feed into the standard pipeline:

```
Survey  ──┐
           ├──> Explore ──> Deepen ──> Review ──> Manuscript
Atlas   ──┘
```

Atlas and survey findings both live in `survey/findings.jsonl` (distinguished by
`"phase": "atlas"` vs `"phase": "survey"`). The explore, deepen, and manuscript
phases consume them identically. The report manifest groups findings by category,
not by phase.
