# Agent Knowledge Base (CLAUDE.md / AGENTS.md)

**Audience:** All AI agents working on Sharur (Claude Code, Codex, Antigravity, etc.)

**CLAUDE.md scope:** Project-level rules and navigation only. **NEVER write dataset-specific content here.** Dataset context belongs in each dataset's directory.

## Project Overview

Sharur is an agent-driven metagenomic exploration system. It's a data plane that makes large metagenomic datasets navigable by AI agents.

## Reference Docs (load on demand)

Detailed guides live in `docs/` and `.claude/skills/`. **Read the relevant doc before starting a task** — don't work from memory of what might be in them.

| When you are... | Read this |
|-----------------|-----------|
| Running a full 5-phase analysis pipeline | `docs/analysis_workflow.md` |
| Writing manuscripts or compiling reports | `docs/manuscript_guide.md` |
| Using tools (Astra, ELSA, ESM3, Foldseek, V2 atoms) | `docs/tools_reference.md` |
| Making biological claims or interpreting annotations | `docs/biological_interpretation.md` |
| Dispatching or working as a subagent | `docs/subagent_guide.md` |
| Working with V2 predicate system | `docs/predicates_v2.md` |
| Validating hydrogenases | `.claude/skills/hydrogenase.md` |
| Querying ELSA synteny results | `.claude/skills/synteny.md` |

**Other key docs:**

| Document | Purpose |
|----------|---------|
| `QUICKSTART.md` | **NEW DATASET INGESTION (START HERE)** |
| `DATA_ORGANIZATION.md` | Data directory structure, archival procedures |
| `src/ingest/README.md` | Ingestion pipeline stages (00-07) |
| `.claude/skills/_validation_protocols.md` | Shared validation protocols for all analysis skills |

## Skills (Claude Code)

| Skill | Purpose |
|-------|---------|
| `explore.md` | Curiosity-driven discovery, locus exploration |
| `survey.md` | Systematic comprehensive survey |
| `characterize.md` | Protein/locus characterization + structural analysis |
| `defense.md` | Defense system identification |
| `prophage.md` | Prophage & viral element detection |
| `metabolism.md` | Metabolic pathway reconstruction |
| `compare.md` | Cross-genome comparative analysis |
| `literature.md` | Literature/database research for functional ambiguity |
| `reviewer_2.md` | Adversarial manuscript claim verification |
| `atlas.md` | Exhaustive genome-by-genome reading |
| `foldseek.md` | ESM3 structure prediction + Foldseek search |
| `hydrogenase.md` | NiFe/FeFe hydrogenase validation |
| `synteny.md` | ELSA-powered synteny discovery |

---

## Core Rules (always apply)

### Data Integrity

**Never write structured data from memory.** When producing JSON, JSONL, TSV, or any output that references findings, annotations, or database entries, ALWAYS read the source data first. Write a script to process actual files — do not reconstruct from conversation history.

### Reproducible Findings

**Every specific number in a finding — in the title, description, evidence, or any other text field — must have a verification query.** This applies to counts, accessions, IDs, taxonomies, sizes, breakdowns — any concrete datum. The write-up phase is where fabrication happens; treat it as a verification step, not a summarization step.

**This includes decompositions.** If a finding says "9 genomes carry X, including 3 from phylum A and 4 from phylum B", that's three verification queries — one for the total AND one for each breakdown. Agents reliably verify headline totals but fabricate the per-genome, per-phylum, and per-category breakdowns in prose. If you can't verify a breakdown number, don't include it.

Findings must include a `verification` field — a list of claim/query/expected triples:

```jsonl
{
  "id": "E200",
  "verification": [
    {
      "claim": "9 genomes carry LANC_like",
      "query": "SELECT COUNT(DISTINCT p.bin_id) FROM annotations a JOIN proteins p ON a.protein_id = p.protein_id WHERE a.source = 'pfam' AND a.name = 'LANC_like'",
      "expected": 9
    },
    {
      "claim": "3 of those are Nanobdellota",
      "query": "SELECT COUNT(DISTINCT p.bin_id) FROM annotations a JOIN proteins p ON a.protein_id = p.protein_id WHERE a.source = 'pfam' AND a.name = 'LANC_like' AND p.taxonomy LIKE '%Nanobdellota%'",
      "expected": 3
    },
    {
      "claim": "ELSA cluster 170014 has genome_support=3",
      "query": "awk -F, '$1==170014 {print $3}' data/dpann/synteny/results/micro_chain_clusters.csv",
      "expected": "3"
    }
  ]
}
```

Queries can be SQL, shell commands, or Python one-liners — whatever reproduces the result. A verification agent can re-run them to confirm or deny. **If you cannot write a verification query for a number, do not write that number.**

### Database Queries

```python
# Column names: 'name' (not annotation_id), 'score' (not bitscore)
# Always COUNT(DISTINCT protein_id) — repeat domains inflate COUNT(*)
# NEVER use correlated subqueries — rewrite as JOINs
# Prefer Sharur operators over raw DuckDB:
result = b.search_by_predicates(has=["unannotated", "giant"]); proteins = result._raw
b.get_neighborhood(protein_id, window=10)
b.get_neighborhood(protein_id, window=5, all_annotations=True)
```

### Scientific Rigor

**Forbidden language:** "confirms/proves/demonstrates" (unless truly definitive), "unprecedented/first ever/groundbreaking", "paradigm-shifting/revolutionary"

**Required language:** "suggests/indicates/supports", "consistent with/compatible with", "to our knowledge", "provides evidence for"

**Common errors:** Domain presence ≠ function proof. MAG absence ≠ biological absence. Single marker ≠ pathway presence. "First in analysis" ≠ "first ever."

### MAG Interpretation

**Absence of evidence ≠ evidence of absence.** MAGs are inherently incomplete. Say "not detected in this MAG (N contigs)" — NOT "genome lacks X." Before claiming "A has X but B doesn't", verify B isn't just more fragmented.

### Manuscript Citations

**NEVER include literature citations from training memory.** Use `[CITE: topic]` placeholders and run the `/literature` agent before finalizing. See `docs/manuscript_guide.md`.

### Check for Functional Detail

**Don't stop at generic predicates — drill into subgroup-level detail.**

- **Hydrogenases:** `hydrogenase` → check `nife_group1`–`nife_group4`, `fefe_groupA`–`fefe_groupC`. Group 3 vs 4 reveals uptake vs evolution.
- **CRISPR:** `cas_domain` → check `type_i_crispr`/`type_ii_crispr`/`type_iii_crispr`, effectors, `loci` table.
- **Defense:** `defense_system` → check DefenseFinder system annotations for specific types. **CRITICAL: NEVER report raw DefenseFinder HMM hits as defense systems.** Raw hits have an 80%+ false positive rate. Only the `defense_systems` table (from MacSyFinder validation via `validate_defense_systems.py`) contains genuine system calls.
- **CAZy:** `carbohydrate_active` → check `cazy:GH5`, `cazy:GT2`, etc. GH families reveal substrate specificity.

---

## Key Files

| File | Purpose |
|------|---------|
| `sharur/predicates/vocabulary.py` | All predicate definitions |
| `sharur/predicates/generator.py` | Main predicate computation |
| `sharur/predicates/mappings/` | PFAM/KEGG/CAZy/VOGdb → predicate mappings |
| `sharur/predicates_v2/` | V2 semantic atom system |
| `config/predicates_v2/` | V2 YAML config: facets, relations, composites |
| `sharur/operators/structure.py` | ESM3 structure prediction |
| `sharur/operators/foldseek.py` | Foldseek structural homology search |
| `sharur/operators/manifest.py` | Analysis manifest for session continuity |
| `sharur/core/hypothesis_registry.py` | Persistent hypothesis store |

## Pipeline Notes

- **MinCED is single-threaded.** Run on the login node, not SLURM. A single-threaded job doesn't justify a cluster allocation.
- **KOFAM is slow.** On large datasets (>100k proteins) KOFAM can take many hours. This is normal. Submit as a SLURM job with generous walltime (48h+).

## Testing

```bash
python -m pytest tests/test_operators/test_predicate_generator.py -v
python -c "from sharur.predicates.vocabulary import ALL_PREDICATES, list_categories; print(f'Total: {len(ALL_PREDICATES)}'); print(f'Categories: {list_categories()}')"
```

## Standard Directory Structure

```
data/{dataset_name}/
├── sharur.duckdb                # Core database
├── manifest.json               # Analysis state
├── source/                     # Input files (.faa)
├── annotations/                # Annotation results
├── embeddings/                 # ESM2 embeddings (H5 + FAISS)
├── structures/                 # ESM3 PDBs + Foldseek results
├── synteny/                    # ELSA synteny results + FAISS store
├── exploration/                # Exploration outputs
├── survey/                     # Survey outputs
├── reports/                    # Generated reports
└── figures/                    # All figures
```

## TODO

- [ ] **Bake superfamily validation deeper into `/survey`** — known traps (HydDB, Mokosh/BREX, Radical_SAM), co-annotation checks, automatic specialist dispatch
- [ ] **Co-location engine (`sharur/colocation.py`) regression test** — run on Susan genomes and DPANN, compare against MacSyFinder output. Spec: `COLOCATION_ENGINE_SPEC.md`
