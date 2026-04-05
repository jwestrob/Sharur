# Findings Specification

**Audience:** All agents producing analytical outputs (`survey`, `exploration`, `characterize`, `defense`, etc.)

**Purpose:** Every scientific discovery must be recorded as a structured, verifiable finding in `findings.jsonl`. Prose reports (`.md`) are for humans. Findings are the canonical scientific record for downstream agents, reports, and audits.

---

## Canonical Ownership

These files are the scientific source of truth:

```text
data/{dataset}/survey/findings.jsonl
data/{dataset}/exploration/findings.jsonl
data/{dataset}/exploration/hypotheses.json
```

Boundary rules:

- `survey/findings.jsonl` and `exploration/findings.jsonl` are the canonical findings archive.
- `exploration/hypotheses.json` is the canonical local hypothesis store.
- `manifest.json` is a derived summary/cache only. It may be regenerated from canonical findings and other dataset outputs.
- `sharur_ops.db` is coordination-only state for multi-agent runs. It is not the long-term scientific archive.

---

## Lean Agent Write Contract

Most agents should write only these fields:

| Field | Type | Requirement | Notes |
|-------|------|-------------|-------|
| `title` | string | Required | One-line summary with qualifiers included. |
| `category` | string | Required | Functional area such as `energy_metabolism`, `defense_systems`, `crispr`, `general`. |
| `description` | string | Required | Interpretation and significance, not just a restatement of the title. |
| `evidence` | object | Required | Structured supporting data. See Evidence Payloads. |
| `verification` | array | Required | List of `{claim, query, expected}` records. See Verification. |
| `protein_ids` | array | Optional | Key supporting proteins. Include enough to reproduce. |
| `contigs` | array | Optional | Key supporting contigs. |
| `falsification` | string | Required if novelty >= 2 | "This finding would be wrong if ___." Must be testable. See Falsification. |

Agents do **not** need to hand-author these in normal workflows:

- `id`
- `phase`
- `schema_version`
- `provenance`
- `related_findings`
- `novelty`
- `confidence`

The system may fill or derive those fields during normalization/storage.

---

## Canonical Stored Finding Shape

Each line in `findings.jsonl` should normalize to this stored shape:

### Core Fields

| Field | Type | Description |
|-------|------|-------------|
| `id` | string | Stable finding ID. Canonical convention: `{phase}-{NNN}` such as `survey-001` or `exploration-014`. |
| `schema_version` | string | Finding schema version. Current canonical value: `"2.0"`. |
| `phase` | string | Producing phase, such as `survey`, `exploration`, `characterization`, `defense`, `metabolism`. |
| `title` | string | Self-contained one-line summary. |
| `category` | string | Functional category. |
| `description` | string | Biological interpretation. |
| `evidence` | object | Structured evidence payload. |
| `verification` | array | Mandatory verification records. |

### Optional / Derived Fields

| Field | Type | Description |
|-------|------|-------------|
| `protein_ids` | array | Canonical top-level protein references for downstream tooling. |
| `contigs` | array | Canonical top-level contig references. |
| `falsification` | string | **Required for novelty >= 2.** States the specific condition that would break this finding. See Falsification. |
| `provenance` | object | Normalized provenance metadata. Preferred keys include `source_report`, `query_used`, `accession_verified`, `agent_id`. |
| `related_findings` | array | Related finding IDs across phases. |
| `novelty` | int | 0=routine, 1=interesting, 2=novel, 3=potentially_significant. |
| `confidence` | float | 0.0-1.0 self-assessed confidence. |
| `n_genomes` | int or null | Number of genomes where the finding applies. |
| `figures` | array | Figure paths associated with the finding. |

### Canonical Example

```json
{
  "id": "exploration-003",
  "schema_version": "2.0",
  "phase": "exploration",
  "title": "DUF6088 forms a conserved genomic module with AbiEii abortive infection toxin",
  "category": "defense_systems",
  "description": "DUF6088, a ~200 aa protein of unknown function, co-occurs with AbiEii toxin domains in 63.5% of instances (40/63 proteins). This association spans 6 different defense system types, suggesting DUF6088 is a genuine accessory component of abortive infection systems rather than a coincidental neighbor.",
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
      "query": "See exploration/duf_defense_candidates.md methods: neighborhood scan with window=2 for AbiEii co-occurrence",
      "expected": 40
    }
  ],
  "protein_ids": [
    "CoronaMine_BoilerAditFilter_100nm_15_12_2024_IL_777_6",
    "CoronaMine_BoilerAditFilter_100nm_15_12_2024_IL_123011_2",
    "CoronaMine_BoilerAditFilter_100nm_15_12_2024_IL_22959_1"
  ],
  "contigs": ["IL_777", "IL_22959", "IL_123011"],
  "provenance": {
    "source_report": "exploration/duf_defense_candidates.md",
    "query_used": "Neighborhood scan with window=2 for AbiEii co-occurrence",
    "agent_id": "explore_duf_defense"
  },
  "related_findings": ["survey-005"],
  "novelty": 3,
  "confidence": 0.85,
  "falsification": "Wrong if the AbiEii PFAM domain cross-reacts with non-defense proteins (e.g., SanaT toxins) making the co-occurrence an annotation artifact rather than a genuine functional association."
}
```

### Minimal Authoring Example

This is valid agent-facing input as long as the normalizer can infer the phase from the destination file:

```json
{
  "title": "DUF6088 forms a conserved genomic module with AbiEii abortive infection toxin",
  "category": "defense_systems",
  "description": "DUF6088 co-occurs with AbiEii toxin domains in a substantial fraction of neighborhoods, suggesting it is a genuine accessory component rather than a coincidental neighbor.",
  "evidence": {
    "duf_accession": "DUF6088",
    "total_in_metagenome": 63,
    "near_defense": 40
  },
  "verification": [
    {
      "claim": "63 DUF6088 proteins in the metagenome",
      "query": "SELECT COUNT(DISTINCT protein_id) FROM annotations WHERE name LIKE '%DUF6088%'",
      "expected": 63
    }
  ],
  "protein_ids": [
    "CoronaMine_BoilerAditFilter_100nm_15_12_2024_IL_777_6"
  ]
}
```

---

## Compatibility and Normalization Policy

Only the compatibility/normalization layer may interpret legacy top-level fields.

### Accepted Legacy Inputs

| Legacy Input | Normalized Behavior |
|--------------|---------------------|
| Top-level `genes` | Legacy-only compatibility input. Normalize to canonical `protein_ids`. |
| Top-level `priority` | Legacy-only compatibility input. Normalize to a derived summary signal such as "key finding". |
| `provenance.query` | Normalize to `provenance.query_used`. |
| Legacy IDs like `E001` or `D001` | Preserve on read if already present. Do not generate new legacy IDs. |

### Important Distinction

- `evidence["genes"]` is **valid canonical evidence** and should remain valid.
- Top-level `finding["genes"]` is **legacy compatibility input only**.
- Legacy string `evidence` values should be migrated to an object wrapper that preserves the original text, for example with keys such as `source_format`, `legacy_text`, and `statements`. New findings should not emit free-text evidence strings.

The normalizer may use `evidence["genes"]` to derive canonical `protein_ids`, but downstream code should read `protein_ids`, not top-level legacy `genes`.

### Validation Behavior

- `verification` is mandatory.
- Missing or malformed verification records must surface as validation issues.
- Do not silently invent verification, provenance, or IDs.

---

## Evidence Payloads

Structure the `evidence` field based on what was found.

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

- SQL against `sharur.duckdb` (preferred)
- Shell commands
- Python one-liners
- References to a concrete methods section when the logic is too large to inline

**If you cannot write a verification query for a number, do not include that number in the finding.**

**Decompositions require separate verification.** "68 DUF1016_N near defense, across 20 system types" is two claims needing two queries, not one.

---

## Falsification

**Required for novelty >= 2.** Before committing a finding, state the specific condition that would make it wrong. This forces you to articulate the failure mode and ideally test it before writing the finding.

The `falsification` field is a plain string answering: "This finding would be wrong if ___."

**Examples:**

```json
{
  "title": "SoxYZ carrier is part of the LanM operon",
  "falsification": "Wrong if SoxYZ and LanM are on different strands, >5 genes apart, or the association disappears when restricted to same-contig same-strand genes within 500bp intergenic distance."
}
```

```json
{
  "title": "YcaO + TfuA co-locate with LanM in 13 genomes",
  "falsification": "Wrong if these 13 genomes are all from 2-3 closely related species sharing whole-genome synteny, making this genome conservation rather than functional co-selection."
}
```

```json
{
  "title": "DUF6088 is a novel defense accessory component",
  "falsification": "Wrong if DUF6088 is already annotated as part of a known defense system in DefenseFinder or if the AbiEii PFAM domain cross-reacts with non-defense proteins."
}
```

**The falsification must be testable.** "This might not be real" is not a falsification. "This would be wrong if the 40/63 co-occurrence drops below 10/63 when restricted to same-contig same-strand pairs" is testable.

**When novelty >= 2, agents should test their own falsification before committing the finding.** If the test breaks the finding, revise or discard it. If it survives, include both the falsification statement and the result of the test.

---

## Locus Validation

Before claiming genes form a "conserved locus" or "operon," verify all of:

1. **Same contig** — genes share the same `contig_id`
2. **Co-oriented** — genes are on the same strand
3. **Short intergenic distance** — adjacent genes are <500bp apart (operonic spacing)
4. **Not just genome-level synteny** — the association should hold across phylogenetically diverse genomes, not just closely related strains sharing whole-chromosome gene order

A query pattern for checking:

```sql
SELECT p1.protein_id, p2.protein_id,
       p1.strand, p2.strand,
       p2.start - p1.end_coord AS intergenic_bp
FROM proteins p1
JOIN proteins p2 ON p1.contig_id = p2.contig_id
  AND p1.bin_id = p2.bin_id
  AND p2.gene_index = p1.gene_index + 1
WHERE p1.protein_id = ?
```

Claims about "conserved neighborhoods" that fail these checks should be downgraded from "operon" to "genomic proximity" and the finding description should reflect the weaker association.

---

## ID Convention

Canonical IDs use the full phase slug:

```text
survey-001, survey-002, ...
exploration-001, exploration-002, ...
characterization-001, characterization-002, ...
defense-001, defense-002, ...
metabolism-001, metabolism-002, ...
```

When generating a new ID, count existing entries in the phase file and emit the next `{phase}-{NNN}` value.

Legacy IDs such as `E001` or `D001` may still appear in older datasets. Preserve them on read; do not emit new ones.

---

## When to Write Findings

Write a finding when the result is:

- Quantitative
- Verifiable
- Interpretable

Do **not** write findings for:

- Routine annotation counts without interpretation
- Negative results unless the absence itself is the finding
- Methodological notes that are not scientific claims

---

## Integration with Skills

Every analytical skill should:

1. Read this spec before writing findings
2. Write `findings.jsonl` alongside the prose `.md` report
3. Include verification queries for every quantitative claim
4. Use canonical slug-style IDs for new findings

Reference snippet for skill prompts:

```text
Read docs/findings_spec.md for the required findings format.
Write findings to {output_dir}/findings.jsonl alongside your prose report.
```
