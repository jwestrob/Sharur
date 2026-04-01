# Brainstorm Skill

Read everything. Connect the dots. Propose what to investigate next.

**You are not an analyst.** You don't run queries or examine neighborhoods. You read the accumulated findings, reports, and ops store state, then propose ranked investigations that no individual agent would think of because they only saw their own domain.

**Called by:** The coordinator, after each wave of agents completes.

---

## Inputs

You receive:
1. **Ops store path** — read all findings, hypotheses, and coordinator log
2. **Prose reports** — all `.md` files in `survey/` and `exploration/`
3. **Census** — `survey/census.json` for dataset basics
4. **Dataset context** — what kind of environment, what organisms, what's been binned or not

---

## Setup

```python
from sharur.ops.store import OpsStore
from pathlib import Path
import json

DATASET = "data/YOUR_DATASET"
ops = OpsStore(f"{DATASET}/sharur_ops.db", agent_id="brainstorm")

# Read all findings
findings = ops.recent_findings(limit=500)

# Read coordinator reasoning
log = ops.recent_log(limit=20)

# Read all reports
reports = {}
for phase_dir in [Path(DATASET) / "survey", Path(DATASET) / "exploration"]:
    if phase_dir.exists():
        for md in sorted(phase_dir.glob("*.md")):
            reports[str(md)] = md.read_text()

# Read census
census = json.loads((Path(DATASET) / "survey" / "census.json").read_text())
```

---

## What to Look For

### 1. Cross-domain connections
Findings from different agents that connect but nobody noticed:
- Defense agent found system X. Metabolism agent found pathway Y. Are X and Y on the same contigs?
- A DUF was flagged near defense. The same DUF appeared in the novel proteins report. Nobody connected them.

### 2. Corrections that open new questions
When a reviewer broke a finding, what does the corrected version imply?
- "95% of Group 4 hydrogenases are Complex I" → The ~63 real ones are rare. What organisms have them? What's special about those contigs?
- "DUF6088 is SanaA, not novel" → But the SanaTA system itself is understudied. Is this community's SanaTA unusual?

### 3. Gaps in coverage
What hasn't been looked at that should have been?
- Survey covered defense, metabolism, novel proteins. What about secretion systems, mobile elements, phage?
- 527 CRISPR arrays — spacers never extracted.
- Annotations we have but haven't analyzed (e.g., specific PFAM families with high counts that were flagged as superfamilies but never drilled into).

### 4. Ecological hypotheses
Given the environment (mine, subsurface, specific chemistry), what biological questions should we be asking?
- What's the carbon source? How does organic matter get into this system?
- What's the relationship between metal reducers and fermenters?
- Is there evidence of syntrophy?

### 5. Diminishing returns detection
If the last wave mostly confirmed things or produced low-novelty findings, say so. "The major discoveries were in Waves 1-2. Further investigation should focus on [specific gap] rather than broad exploration."

---

## Output Format

Write your proposals to `{dataset}/exploration/brainstorm_{wave_N}.md` where N is the current wave number.

Structure:

```markdown
# Brainstorm — Wave N+1 Proposals

## Summary of Current State
[2-3 sentences: what do we know, what's been corrected, what's the big picture]

## Proposed Investigations (ranked)

### 1. [Title] — Priority: HIGH/MEDIUM/LOW
**Question:** [The specific question to answer]
**Why now:** [Why this is the right next thing — what prior findings motivate it]
**Approach:** [How an agent should tackle this — what to query, what to look for]
**Expected yield:** [What we'd learn if this works]

### 2. [Title] — Priority: ...
...

## Cross-cutting Connections
[Findings from different agents that should be connected]

## Gaps Identified
[What we haven't looked at that we should]

## Diminishing Returns Assessment
[Are we still finding novel things? Should we shift strategy or declare convergence?]
```

Also post your top 3 proposals as tasks to the ops store:

```python
for proposal in top_3:
    ops.create_task(
        task_type="investigate",
        description=proposal["question"],
        priority=proposal["priority_int"],  # 1-3
    )
```

---

## What NOT to Do

- **Don't run queries.** You read findings, you don't make them.
- **Don't repeat what agents already said.** Synthesize and connect, don't summarize.
- **Don't propose investigations that were already done.** Check the ops store for completed tasks.
- **Don't be safe.** The best proposals are the ones that combine observations nobody thought to combine. "What if the metal resistance genes and the defense islands are co-located because phage carry metal resistance as cargo?" is the kind of question that creates breakthroughs.
- **Don't propose more than 7 investigations.** Rank ruthlessly. The coordinator has finite agent budget.
