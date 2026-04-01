# Findings Specification

**Audience:** All agents producing analytical outputs (survey, explore, characterize, defense, etc.)

**Purpose:** Every discovery must be recorded as a structured, verifiable finding in `findings.jsonl`. Prose reports (`.md`) are for human consumption; findings are for machines, auditors, and downstream agents.

---

## Location

Each analysis phase writes to its own findings file within the dataset directory:

```
data/{dataset}/survey/findings.jsonl
data/{dataset}/exploration/findings.jsonl
```

The coordinator (when active) may also maintain a unified `data/{dataset}/findings.jsonl` that merges across phases.

---

## Schema

Each line in `findings.jsonl` is a self-contained JSON object with these fields:

### Required Fields

| Field | Type | Description |
|-------|------|-------------|
| `id` | string | Stable ID. Convention: `{phase}-{NNN}` (e.g., `survey-001`, `explore-014`). Auto-increment within file. |
| `title` | string | Self-contained one-line summary with all qualifiers. Must be meaningful without reading the description. |
| `category` | string | Functional area: `energy_metabolism`, `defense_systems`, `cell_surface`, `novel_proteins`, `metal_resistance`, `carbon_metabolism`, `crispr`, `secondary_metabolism`, `general` |
| `description` | string | Prose paragraph with biological interpretation. Not a restatement of the title — explain *why this matters*. |
| `evidence` | object | Quantitative evidence. Structure depends on finding type (see below). |
| `verification` | array | **MANDATORY.** List of `{claim, query, expected}` triples. Every specific number in the title, description, or evidence must have a verification query. See Verification section. |
| `phase` | string | Which analysis phase produced this: `survey`, `exploration`, `characterization`, `defense`, `metabolism`, etc. |

### Recommended Fields

| Field | Type | Description |
|-------|------|-------------|
| `confidence` | float | 0.0–1.0. Agent's self-assessed confidence in the finding. |
| `novelty` | int | 0=routine, 1=interesting, 2=novel, 3=potentially_significant |
| `n_genomes` | int | Number of genomes where finding applies. For unbinned metagenomes, use `null`. |
| `protein_ids` | array | Key protein IDs supporting the finding. Include enough to reproduce, not all. |
| `contigs` | array | Key contig IDs. |
| `provenance` | object | `{source_report, query_used, accession_verified, agent_id}` |
| `figures` | array | Paths to any figures generated for this finding. |
| `related_findings` | array | IDs of related findings in any phase. |

---

## Evidence Payloads

Structure the `evidence` field based on what was found:

```json
// Gene/protein-level finding
{
  "protein_id": "IL_826_68",
  "contig_id": "CoronaMine_..._IL_826",
  "accession": "PF13437",
  "name": "MtrB_PioB",
  "score": 559.5,
  "sequence_length": 698
}

// Neighborhood/cassette finding
{
  "genes": [
    {"protein_id": "...", "annotation": "MtrB_PioB", "gene_index": 67},
    {"protein_id": "...", "annotation": "Paired_CXXCH_1", "gene_index": 68},
    {"protein_id": "...", "annotation": "MtrC/MtrF", "gene_index": 69}
  ],
  "contig_id": "...",
  "span_bp": 12000
}

// Statistical/community-level finding
{
  "metric": "fold_enrichment",
  "value": 32.4,
  "n_observed": 68,
  "n_expected": 2.1,
  "system_types_associated": 20,
  "total_in_metagenome": 312
}

// Defense system finding
{
  "system_type": "RM_Type_II",
  "n_systems": 599,
  "n_genes": 1842,
  "fp_reduction": "85.2%",
  "source_table": "defense_systems"
}
```

---

## Verification

**Every specific number must have a verification query.** This is the single most important field. Without it, findings are unauditable.

```json
{
  "verification": [
    {
      "claim": "62 MtrB_PioB proteins in the metagenome",
      "query": "SELECT COUNT(DISTINCT protein_id) FROM annotations WHERE name = 'MtrB_PioB'",
      "expected": 62
    },
    {
      "claim": "36 of 62 MtrB_PioB have Paired_CXXCH_1 within 2 genes",
      "query": "Python: see mtr_pathway_colocation.md methods section for full neighborhood scan",
      "expected": 36
    },
    {
      "claim": "DUF6088 appears 63 times in the metagenome",
      "query": "SELECT COUNT(DISTINCT protein_id) FROM annotations WHERE accession LIKE '%DUF6088%'",
      "expected": 63
    }
  ]
}
```

Queries can be:
- **SQL** against `sharur.duckdb` (preferred — most reproducible)
- **Shell commands** (awk, grep — for file-based checks)
- **Python one-liners** (for complex logic)
- **References** to methods in the prose report (acceptable when the full query is complex, but include the report path)

**If you cannot write a verification query for a number, do not include that number in the finding.**

**Decompositions require separate verification.** "68 DUF1016_N near defense, across 20 system types" is two claims needing two queries — one for the count, one for the system type diversity.

---

## ID Convention

```
survey-001, survey-002, ...      # Survey phase findings
explore-001, explore-002, ...    # Exploration phase findings
char-001, char-002, ...          # Characterization phase findings
defense-001, defense-002, ...    # Defense-specific findings
```

Auto-increment within each file. When reading an existing `findings.jsonl`, count existing entries to determine the next ID.

---

## Complete Example

```json
{
  "id": "explore-003",
  "title": "DUF6088 forms a conserved genomic module with AbiEii abortive infection toxin",
  "category": "defense_systems",
  "description": "DUF6088, a ~200 aa protein of unknown function, co-occurs with AbiEii toxin domains in 63.5% of instances (40/63 proteins). This association spans 6 different defense system types, suggesting DUF6088 is a genuine accessory component of abortive infection systems rather than a coincidental neighbor. It may function as a regulatory or modulatory subunit controlling AbiEii toxin activation.",
  "evidence": {
    "duf_accession": "DUF6088",
    "total_in_metagenome": 63,
    "near_defense": 40,
    "fraction_near_abiEii": 0.635,
    "defense_system_types_associated": 6,
    "example_loci": ["IL_777", "IL_22959", "IL_123011"]
  },
  "verification": [
    {
      "claim": "63 DUF6088 proteins in the metagenome",
      "query": "SELECT COUNT(DISTINCT protein_id) FROM annotations WHERE name LIKE '%DUF6088%'",
      "expected": 63
    },
    {
      "claim": "40 of 63 have AbiEii within 2 genes",
      "query": "See duf_defense_candidates.md methods: neighborhood scan with window=2 for AbiEii co-occurrence",
      "expected": 40
    }
  ],
  "confidence": 0.85,
  "novelty": 3,
  "n_genomes": null,
  "protein_ids": [
    "CoronaMine_BoilerAditFilter_100nm_15_12_2024_IL_777_6",
    "CoronaMine_BoilerAditFilter_100nm_15_12_2024_IL_123011_2",
    "CoronaMine_BoilerAditFilter_100nm_15_12_2024_IL_22959_1"
  ],
  "phase": "exploration",
  "provenance": {
    "source_report": "exploration/duf_defense_candidates.md",
    "agent_id": "explore_duf_defense"
  },
  "related_findings": ["survey-005"]
}
```

---

## When to Write Findings

Write a finding when you discover something that is:
- **Quantitative** — has specific counts, enrichments, or measurements
- **Verifiable** — you can write SQL or a command to reproduce the numbers
- **Interpretable** — it means something biologically, not just "X exists"

Do NOT write findings for:
- Routine annotation counts without interpretation ("1,276,251 PFAM annotations")
- Negative results ("no CRISPR arrays found") unless the absence is itself surprising
- Methodological notes ("used window=10 for neighborhood analysis")

---

## Integration with Skills

Every skill that produces analytical outputs should:
1. **Read this spec** before writing findings
2. **Write `findings.jsonl`** alongside the prose `.md` report
3. **Include verification queries** for every quantitative claim
4. **Use the ID convention** for the current phase

Reference this document in skill prompts:
```
Read docs/findings_spec.md for the required findings output format.
Write findings to {output_dir}/findings.jsonl alongside your prose report.
```
