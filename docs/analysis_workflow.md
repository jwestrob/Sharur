# Analysis Workflow Reference

**Load this doc when:** Running a full 5-phase analysis pipeline, orchestrating subagents, or working on the Deepen/Review phases.

## Orchestration Principles

**Every task is a separate subagent with curated context.**

The orchestrator (the main conversation or a dedicated coordinator agent) does NOT do analysis
itself. It reads state, decides what to do next, crafts a focused prompt for each subagent,
dispatches, reads outputs, and repeats.

**Two agent levels:**

- **Coordinator agents** (orchestrator, Deepen coordinator): dispatch subagents, run loops,
  make phase-transition decisions. Can be the main conversation or a background agent.
- **Task agents** (survey topic, explore hypothesis, literature search, section writer):
  leaf agents that produce outputs and return. Do not dispatch further agents.

Each task agent:

- Gets **curated context**: only the findings, data, and instructions relevant to its task
- Has a **clear output specification**: what files to produce, what format
- Has **completion criteria**: findings must include provenance, locus-level findings must include figures
- Runs with **clean context**: no prior conversation history bleed

**Orchestrator responsibilities:**
1. Read current state (findings, figures, reports, gaps)
2. Decide what phase/task comes next
3. Craft a prompt with curated context for each subagent
4. Dispatch subagents (parallel when independent, sequential when dependent)
5. Read subagent outputs, verify completion criteria
6. Repeat until phase is complete, then advance to next phase

## Five-Phase Architecture

**CRITICAL: Always include independent review at the end of analysis.**

```
Orchestrator (reads state, crafts prompts, dispatches subagents)
│
├── 1. Survey (systematic coverage)
│   ├── Dispatch parallel topic subagents: metabolic, defense, surface, etc.
│   │   Each subagent → survey/findings.jsonl entries + summary figures
│   └── Orchestrator reads all survey findings, plans Explore
│
├── 2. Explore (hypothesis-driven discovery)
│   ├── Orchestrator reads survey → generates hypotheses
│   ├── Dispatch hypothesis-testing subagents (parallel)
│   │   Each subagent → exploration/findings.jsonl entries + locus figures
│   └── Orchestrator reads results, decides if more exploration needed
│
├── 3. Deepen (successive rounds — see below)
│   ├── Round 1: Read all findings → identify gaps → dispatch task agents
│   │   Gaps may need data analysis, literature research, or both (parallel)
│   ├── Read new findings → check for remaining gaps
│   ├── Round N: Dispatch more task agents if needed
│   └── Create report manifest + compile COMPREHENSIVE_REPORT.md
│
├── 4. Review (independent validation)
│   ├── Run: verify_claims.py → CLAIM_VERIFICATION.jsonl (mechanical)
│   ├── Dispatch /reviewer_2 agent → review/correction_queue.md (interpretive)
│   └── Fix all critical/meaningful corrections before proceeding
│
└── 5. Manuscript (section-by-section — see docs/manuscript_guide.md)
    ├── Dispatch section-writing subagents (sequential by topic)
    ├── Orchestrator assembles MANUSCRIPT.md from sections
    ├── Run: validate_provenance.py → PROVENANCE_AUDIT.md
    ├── Dispatch literature agent — this is NOT optional
    └── Apply citation corrections, re-render PDF
```

## Deepen Phase (Successive Rounds)

**Goal**: Turn preliminary findings into well-supported, publication-ready results. Unlike Explore
(which follows its own curiosity), Deepen is the orchestrator making deliberate decisions about
what needs more evidence — and running as many rounds as needed until gaps are filled.

**Workflow (iterative):**

**Round 1:**

1. **Read all findings** from Survey and Explore (`findings.jsonl`, `genome_profiles.tsv`)
2. **Gap analysis** — identify:
   - Findings with weak evidence (single annotation, no neighborhood context, no figure)
   - Unknown proteins that deserve characterization
   - Metabolic claims missing pathway completeness checks
   - Defense systems needing detailed inventory
   - Novel findings lacking literature context
   - Findings missing figures (locus-level findings MUST have a neighborhood diagram)
3. **Dispatch targeted task agents** — one per gap, each with curated context.
   Gaps may require data analysis, literature research, or both:

```python
# Data analysis task agent
Task(subagent_type="general-purpose",
     prompt="""Deepen finding survey-001 (hydrogenase survey).
     Gap: No neighborhood validation of NiFe Group 4 calls.
     DB: data/DATASET/sharur.duckdb
     Task: Pick 3-5 Group 4 hits, run get_neighborhood(pid, window=8),
     check for Hyf/Hyc operon context. Generate neighborhood figures.
     Output: Append findings to exploration/findings.jsonl with id='D001'.
     Save figures to figures/.""",
     run_in_background=True)

# Literature task agent (dispatched in parallel with data agents)
Task(subagent_type="general-purpose",
     prompt="""Literature search for finding E001 (DsrAB in 9/41 genomes).
     Questions: Is DsrAB known in other candidate phyla? What is the typical
     operon structure? Are reverse-DsrAB variants reported in non-sulfate-reducers?
     Output: Append to exploration/findings.jsonl with literature provenance.""",
     run_in_background=True)
```

4. **Read task agent outputs** — verify each produced findings + figures

**Round 2+ (if needed):**

5. **Re-read all findings** including new Deepen findings
6. **Check for remaining gaps** — did the subagents actually fill the gaps?
   - If a subagent's findings raise new questions → dispatch another round
   - If gaps remain unfilled (e.g., structure prediction failed) → note as limitation
   - If all major gaps filled → move on
7. **Stop criterion:** No more gaps that would weaken a manuscript claim, or diminishing returns

**After final round:**

8. **Create report manifest** — organize all findings into thematic sections
9. **Build report manifest PROGRAMMATICALLY** — read `findings.jsonl` files and group by category:
```python
python3 -c "
import json
from pathlib import Path
from collections import defaultdict
findings = []
for sub in ['survey', 'exploration']:
    fp = Path('data/DATASET/') / sub / 'findings.jsonl'
    if fp.exists():
        findings.extend(json.loads(l) for l in fp.read_text().splitlines() if l.strip())
by_cat = defaultdict(list)
for f in findings:
    by_cat[f.get('category','uncategorized')].append(f.get('id'))
for cat, ids in sorted(by_cat.items()):
    print(f'{cat}: {ids}')
"
```
**NEVER write report_manifest.json from memory.** Always read the actual findings files first.

10. **Generate comprehensive report:**
```bash
python scripts/compile_comprehensive_report.py --dataset data/DATASET_NAME/
```

**Key principles:**
- The orchestrator decides what's worth investing in — not every finding needs deepening
- Focus on: findings likely to be headline results, claims reviewers will challenge, unknowns that could change the story
- **Successive rounds** prevent single-pass gaps
- Each round's subagents are independent and can run in parallel

## Review Agent Responsibilities

**Goal**: Independent validation that strengthens confidence in findings and catches gaps before publication.

1. **Validate protein counts with direct queries**
   ```sql
   SELECT COUNT(DISTINCT a.protein_id) as n_proteins,
          COUNT(DISTINCT p.bin_id) as n_genomes
   FROM annotations a
   JOIN proteins p ON a.protein_id = p.protein_id
   WHERE a.accession = 'PF12345'
   ```

2. **Cross-check PFAM functions using predicates**
   ```python
   radical_sam = b.search_by_predicates(has=['radical_sam'])
   ```

3. **Verify statistical tests** — re-run independently, check effect sizes not just p-values

4. **Assess assembly quality for absence claims** — high-quality (<50 contigs) vs fragmented (>200 contigs)

5. **Perform missing analyses** — HydDB, CheckM2, pathway completeness as needed

6. **Publication readiness assessment** — specific blocking issues + actionable fixes
