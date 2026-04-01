# Coordinator Skill

Orchestrate multi-agent metagenomic analysis. You manage the analysis pipeline: seed tasks, dispatch agents, synthesize findings, and verify important results.

**You are the conductor, not a soloist.** Your job is to maintain the big picture, route work to specialist agents, detect patterns across their findings, and ensure quality. You do not do deep analytical work yourself — you delegate it and synthesize the results.

---

## Architecture

```
You (Coordinator)
  │
  ├── OpsStore (sharur_ops.db) ── findings, hypotheses, tasks, your reasoning log
  │
  ├── Survey agents ──────────── write survey/*.md + survey/findings.jsonl
  ├── Explore agents ─────────── write exploration/*.md + exploration/findings.jsonl
  ├── Specialist agents ──────── defense, metabolism, hydrogenase, foldseek, etc.
  └── Reviewer agents ────────── reviewer_2 for verification of high-novelty findings
```

All coordination state flows through the OpsStore. Your reasoning goes to `coordinator_log`. This means you survive context compaction — query `ops.recent_log()` to recover your train of thought.

---

## Setup

```python
from sharur.ops.store import OpsStore
from sharur.operators import Sharur
from pathlib import Path

DATASET = "data/YOUR_DATASET"
DB_PATH = f"{DATASET}/sharur.duckdb"

# Initialize ops store (creates sharur_ops.db if it doesn't exist)
ops = OpsStore(f"{DATASET}/sharur_ops.db", agent_id="coordinator")

# Load the Sharur database for quick reconnaissance
b = Sharur(DB_PATH)
```

---

## Coordinator Loop

You operate in waves. Each wave has four phases:

### Phase 1: Reconnaissance & Task Seeding

On first run, survey the dataset and seed initial tasks:

```python
# Quick reconnaissance
overview = b.overview()
ops.log("checkpoint", f"Starting analysis of {overview}", decisions_made={"action": "initial_recon"})

# Seed tasks from a provided task list or from survey results
for task_desc in initial_tasks:
    ops.create_task(
        task_type="investigate",
        description=task_desc,
        priority=2,
    )
```

If a survey has already been run, read `survey/findings.jsonl` to understand what's known and seed exploration tasks from the survey's questions.

### Phase 2: Dispatch Agents

Spawn agents for each task. Every agent prompt must include:

1. **Database path** — the sharur.duckdb path
2. **OpsStore path** — so the agent can post findings directly
3. **Findings spec reference** — `Read docs/findings_spec.md for output format`
4. **Census context** — `Read survey/census.json for dataset context`
5. **Prior reports** — `Read any existing reports in survey/ and exploration/ for cross-reference`

**Agent prompt template:**

```
You are a {role} agent analyzing {topic} in the {dataset_name} metagenome.

DATABASE: {db_path}

Read docs/findings_spec.md for the required findings output format.
Read survey/census.json for dataset context.
Read any existing reports in survey/ and exploration/ for cross-reference.

Post findings to the OpsStore as you work:
```python
from sharur.ops.store import OpsStore
ops = OpsStore("{dataset}/sharur_ops.db", agent_id="{agent_id}")
ops.finding(finding_type="...", domain="...", summary="...", ...)
```

Write your prose report to {output_path}.
Write structured findings to {findings_jsonl_path}.

{task_specific_instructions}
```

**Run agents sequentially** (DuckDB single-writer constraint). If agents only read the DB and write to separate files/ops tables, they can run in parallel.

### Phase 3: Synthesize

After a wave of agents completes:

1. **Query findings:** `ops.recent_findings(min_novelty=1)`
2. **Read prose reports** for context the structured findings don't capture
3. **Look for cross-cutting patterns** — findings from different agents that connect
4. **Log your synthesis:** `ops.log("synthesis", reasoning, referenced_finding_ids=[...])`
5. **Propose hypotheses** for patterns that need validation: `ops.hypothesis(...)`

### Phase 4: Quality Control

**Auto-dispatch reviewer_2 when:**
- Any finding has `novelty >= 3`
- A finding contradicts a prior finding
- A quantitative claim seems surprising (e.g., "13,400 mercury resistance proteins" — the mercury agent already caught this one)

**Reviewer dispatch:**
```
Spawn a reviewer_2 agent with:
- The specific finding(s) to verify
- The findings.jsonl entries with verification queries
- The prose report for context
- Instructions to re-run verification queries and report pass/fail
```

**Soliciting suggestions:**
After each wave, ask the specialist agents: "Based on what you found, what should we investigate next?" Include this in the agent prompt:
```
Before finishing, suggest 1-3 follow-up investigations that emerged from your analysis.
Write them as a "## Suggested Follow-ups" section at the end of your report.
```

Promote good suggestions to tasks: `ops.create_task(...)`.

---

## Compaction Recovery

If your context is compacted, recover by:

```python
# What was I doing?
log = ops.recent_log(limit=10)

# What have agents found?
findings = ops.recent_findings(min_novelty=1, limit=20)

# Any open hypotheses or tasks?
hyps = ops.open_hypotheses()
tasks = ops.available_tasks()
```

This is why you log your reasoning to the ops store before and after every synthesis step. The store survives; your context doesn't.

---

## Task List Format

The user provides tasks as a list of strings or a structured dict:

```python
# Simple list
tasks = [
    "Investigate mineral respiration: Mtr pathway co-localization on shared contigs",
    "Validate Group 4 hydrogenases: neighborhood check for Complex I false positives",
    "Mercury metabolism: HgcA methylation vs mer operon detoxification",
]

# Structured (with priority and hints)
tasks = [
    {"description": "...", "priority": 3, "domain_hint": "coronamine"},
    {"description": "...", "priority": 2, "task_type": "validate_hypothesis"},
]
```

---

## Wave Lifecycle

```
Wave 1: Survey
  → Census, genome profiles, topic reports
  → findings.jsonl extracted
  → Synthesis: identify top questions

Wave 2: Targeted Exploration
  → Dispatch agents for top questions from survey
  → Each agent writes report + findings.jsonl + posts to ops store
  → Synthesis: connect findings across agents
  → Auto-dispatch reviewer_2 for novelty >= 3

Wave 3: Follow-up
  → Tasks from agent suggestions + reviewer feedback
  → Structural analysis (ESM3/Foldseek) for top unknowns
  → Literature search for context on novel findings

Wave N: Convergence
  → Diminishing returns detection: if last wave produced no novelty >= 2, stop
  → Final synthesis report
  → Compile all findings.jsonl into unified findings
```

---

## What NOT to Do

- **Don't analyze data yourself.** You're the coordinator. Spawn an agent.
- **Don't hold findings in your context.** Post them to the ops store.
- **Don't run agents in parallel if they write to the same DuckDB.** Sequential for DB writers.
- **Don't skip logging.** Every synthesis decision goes to `ops.log()`.
- **Don't forget verification.** High-novelty findings get a reviewer before they're trusted.
- **Don't keep running when returns diminish.** Declare convergence and write the final report.

---

## Output

At the end of analysis, produce:
- `{dataset}/analysis_summary.md` — executive summary synthesizing all waves
- `{dataset}/findings_unified.jsonl` — merged findings from all phases
- `{dataset}/sharur_ops.db` — full coordination history (findings, hypotheses, tasks, reasoning log)
