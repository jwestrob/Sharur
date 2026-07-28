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
Import: from sharur.operators import Sharur; b = Sharur("data/DATASET/sharur.duckdb", read_only=True)
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
- Local read-only DuckDB sessions may run concurrently with explicit process
  budgets. Large coordinated campaigns use the shared query service described
  below. Run database writers sequentially. Use separate output files or the
  locked findings writer for report records.

### Passing confounds to a sub-agent

Quote a confound with its trigger condition. Never paraphrase one into a general
rule.

Every row in the confound table is keyed on an observable — the accession,
source or pattern that puts you at risk. Strip that key and the rule applies
wherever a word matches, which produces two failures at once: the check fires on
data it cannot inform, and the sub-agent reports having "passed" it, which reads
as rigour while being noise. A rule that cannot fail on the data in front of it
is worse than no rule, because it launders an unexamined claim.

```
BAD   "Only credit hydrogenase with NiFeSe_Hases."
GOOD  "hyddb group 4 / energy_conserving / ech_hydrogenase: the group-4
       catalytic subunit is homologous to Complex I NuoD, so most hits are
       Complex I. Require NiFeSe_Hases on the same protein to credit a group-4
       hydrogenase. This does not apply to other groups."
```

The general form: **when you observe X, the trap is Y, so require Z — and state
where it does not apply.** If a sub-agent reports clearing a confound, check
that the trigger was present at all.

## Coordinated Worker Protocol

Use the Ops HTTP service for distributed agents. One service process owns
`sharur_ops.db`; workers use per-agent credentials and communicate through
`SharurOps`.

```python
from sharur.ops.client import SharurOps

ops = SharurOps(
    "http://ops-host:8811",
    agent_id="worker-07",
    api_token=worker_token,
)

task = ops.claim_next_task(campaign_id=campaign_id)
if task is not None:
    # heartbeat_task(), complete_task(), and fail_task() automatically use the
    # attempt token cached from this claim.
    ops.heartbeat_task(task["id"])
```

The task is logical coordination. Its executor can be a persistent process,
local launcher, Slurm job, or array element. Keep scheduler packing and array
construction in the executor layer.

For a large shared dataset, use the same agent token with Sharur Query:

```python
from sharur.query import SharurQuery

query = SharurQuery(
    "http://query-host:8812",
    api_token=worker_token,
)
hits = query.search_predicates(
    has=["hydrogenase", "membrane_protein"],
    limit=100,
)
```

The service supplies one shared DuckDB cache, thread/memory/spill ceilings,
weighted admission, execution timeouts, cancellation, and bounded results.
Findings and task state continue through `SharurOps`. Include the query URL in
worker prompts when this mode is active.

Register each direct-mode worker with generic capacity and capabilities. Match
its DuckDB process budget to that registration:

```python
from sharur.operators import Sharur

b = Sharur(
    "data/DATASET/sharur.duckdb",
    read_only=True,
    duckdb_threads=4,
    duckdb_memory_limit="8GB",
    duckdb_temp_directory="/local-scratch/worker-07",
)
```

Task attempts have finite leases and opaque fencing tokens. A recovered or
replacement attempt has exclusive terminal-write authority, including when a
replacement process reuses the same agent ID.

Long logical tasks can persist a retry-visible checkpoint under that same
attempt fence:

```python
ops.put_task_checkpoint(
    task["id"],
    "progress",
    cursor=opaque_cursor,
    payload={"completed_units": completed_units},
)
```

A replacement attempt reads the prior value with
`get_task_checkpoint()`. Checkpoint writes retain the same token, attempt,
identity, status, and server-clock lease checks as completion. Batch
checkpoint writes at scientifically safe boundaries to control SQLite and
HTTP mutation volume.

Atlas campaigns use the specialized plan and packet contract in
`.claude/skills/atlas.md`: one task per genome, adaptive sequence-free packets
containing data from exactly one bin, packet-cursor checkpoints, a mandatory
zero-model-call packet census, and exact frame/contig/protein coverage
manifests. Dataset size changes the number of dynamically claimed tasks while
preserving exhaustive traversal.

### Hierarchical review workers

Read `docs/review_workflow.md` before serving a review profile. Register one
or more semantic capabilities, for example:

```python
registration = operator.register_agent(
    "independent-reviewer-openai-01",
    role="worker",
    capabilities=["profile:independent_openai"],
    max_concurrent_tasks=1,
    capacity_cpu_slots=1,
)
```

Claimed `scientific_review` tasks provide a frozen target, review tier,
execution profile, resolved provider/model/effort, policy hash, rubric,
blindness flags, and source-review manifest. The submitting identity must own
the task. Its review record must reproduce the task's target and execution
contract exactly.

Independent profiles receive empty source-review manifests. Adjudication
profiles receive the two frozen independent reviews. This data flow supplies
hierarchy across Codex and Claude executors while keeping model-session
creation in the executor layer.

Every promotion review appends executable verification results before the
controller advances it. Use `sharur-review status` for queue and funnel
metrics and `sharur-review trace` for a bounded provenance reconstruction.

See [`query_service.md`](query_service.md) for launch commands, endpoint
contracts, resource arithmetic, and access-path selection.

## Practical Tips

- **Parallel genome browser agents** work well (quarters or groups)
- **JSONL for findings** — easy to append, merge, and process
- **Check database schema** before writing queries (`DESCRIBE table_name`)
- **Research external data** via WebFetch — don't guess PDB functions
- **Don't create new visualization code** when operators exist

## Large Dataset Performance (>50k proteins)

**Rules:**

1. **Never** `b.search_proteins(query="")` on large datasets — use SQL counts or specific predicates
2. **Always** check result size before iterating (use `result.records`)
3. **Always** limit iteration (e.g., `for pid in proteins[:10]`)
4. **Aggregate** in SQL, not Python loops
5. **If query takes >5 seconds**, stop and refine
6. **NEVER use correlated subqueries** — rewrite as JOINs

```python
# BAD
all_proteins = b.search_proteins(query="")
for protein in all_proteins.records: ...

# GOOD
result = b.search_by_predicates(has=["hydrogenase", "membrane_protein"])
proteins = result.records
for protein in proteins[:10]: ...

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
