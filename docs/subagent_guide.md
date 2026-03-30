# Subagent Guide

**Load this doc when:** Dispatching analysis subagents, or when working as a subagent on a large dataset.

## Context Injection (CRITICAL)

**Subagents receive ONLY the prompt parameter.** They do NOT automatically get CLAUDE.md, skill specs (`.claude/skills/*.md`), conversation history, or memory files. The quality of subagent work depends entirely on what context the orchestrator includes in the prompt.

**When dispatching analysis subagents, include in the prompt:**
1. **Database path** and dataset context (genome count, phylum, environment)
2. **Domain documentation path** — tell subagents to read reference docs on demand
3. **Relevant skill spec content** — read the appropriate `.claude/skills/*.md` and include key sections
4. **Superfamily awareness rule** — co-annotation + neighborhood checks before functional claims
5. **Scientific rigor guidelines** — MAG interpretation ("not detected" not "absent"), forbidden language
6. **Database query patterns** — column names, COUNT(DISTINCT), no correlated subqueries
7. **Output specification** — what files to produce, what format, completion criteria

**Minimum boilerplate for any analysis subagent prompt:**
```
DB: data/DATASET/sharur.duckdb
Import: from sharur.operators import Sharur; b = Sharur("data/DATASET/sharur.duckdb")
Columns: 'name' (not annotation_id), 'score' (not bitscore)
Counts: always COUNT(DISTINCT protein_id) — repeat domains inflate COUNT(*)
MAG caveat: "Not detected" not "absent"
Findings: JSONL with id, title, category, description, evidence, n_genomes, provenance
Verification: EVERY specific number in a finding (totals AND breakdowns) needs a verification query.
  If you write "9 genomes have X, 3 from phylum A", that's TWO queries — one for 9, one for 3.
  If you can't write a query for a number, don't write the number.

## Reference Docs (READ ON DEMAND)
When you encounter a domain-specific situation, read the relevant doc before acting:
  - docs/biological_interpretation.md — context-first protocol, MAG caveats, claim escalation, scientific rigor
  - docs/tools_reference.md — Astra, ELSA, ESM3, Foldseek, V2 predicates
  - docs/manuscript_guide.md — findings schema, claim tracking, figure requirements
  - .claude/skills/hydrogenase.md — HydDB curation, Complex I FP detection, neighborhood validation
  - .claude/skills/synteny.md — ELSA query patterns, gene ID formats, cluster citing rules
Use the Read tool to load any doc when you need its protocol. Don't guess — look it up.
```

**For survey/explore subagents**, read the full skill spec and paste the protocol sections, finding schema, and validation rules into the prompt.

## Sub-Agent Protocol
- Sub-agents CAN spawn further sub-agents for specialist tasks (literature, foldseek, hydrogenase curation)
- Each sub-agent produces a discrete output (markdown report, JSONL findings)
- Parent agent synthesizes outputs from all sub-agents
- **DuckDB write locks** — run subagents SEQUENTIALLY if they write to the database

## Practical Tips
- **Parallel genome browser agents** work well (quarters or groups)
- **JSONL for findings** — easy to append, merge, and process
- **Check database schema** before writing queries (`DESCRIBE table_name`)
- **Research external data** via WebFetch — don't guess PDB functions
- **Don't create new visualization code** when operators exist

## Large Dataset Performance (>50k proteins)

**Rules:**
1. **Never** `b.search_proteins(query="")` on large datasets — use SQL counts or specific predicates
2. **Always** check result size before iterating (use `result._raw`)
3. **Always** limit iteration (e.g., `for pid in proteins[:10]`)
4. **Aggregate** in SQL, not Python loops
5. **If query takes >5 seconds**, stop and refine
6. **NEVER use correlated subqueries** — rewrite as JOINs

```python
# BAD
all_proteins = b.search_proteins(query="")
for pid in all_proteins._raw: ...

# GOOD
result = b.search_by_predicates(has=["hydrogenase", "membrane_protein"])
proteins = result._raw
for pid in proteins[:10]: ...

# GOOD — SQL aggregation
stats = b.store.execute("""
    SELECT COUNT(*) as total, AVG(sequence_length) as avg_size
    FROM proteins
""")[0]
```

## DuckDB Query Patterns for Large Datasets

**CRITICAL: Correlated subqueries destroy performance on >1M row tables.** Never put a subquery inside WHERE that references the outer row. Instead:

1. **Materialize the seed set** in a CTE with `contig_id` and `gene_index`
2. **JOIN** to find neighbors (same contig, gene_index ± window)
3. **JOIN** to annotations/predicates for neighbor features
4. **Use GROUP BY** for per-genome cross-tabs, not Python loops over genomes
5. **Use CASE WHEN** inside COUNT(DISTINCT ...) for multi-marker pivots
