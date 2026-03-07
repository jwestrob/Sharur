# Agent Knowledge Base (CLAUDE.md / AGENTS.md)

**Audience:** All AI agents working on Sharur (Claude Code, Codex, Antigravity, etc.)

This file provides shared context, canonical tools, and best practices for agent-driven metagenomic analysis.

**Quick Start:**
- **Data organization:** `DATA_ORGANIZATION.md`
- **Manuscript workflow:** Write `MANUSCRIPT.md` → render with pandoc (see § Manuscript Rendering)
- **Skills (Claude Code):** `.claude/skills/*.md`
- **Validation protocols:** `.claude/skills/_validation_protocols.md`

## Project Overview

Sharur is an agent-driven metagenomic exploration system. It's a data plane that makes large metagenomic datasets navigable by AI agents.

**CLAUDE.md scope:** This file contains project-level tools, patterns, and best practices only. **NEVER write dataset-specific content here** — no analysis results, dataset stats, findings, missing data notes, or per-dataset TODOs. Dataset-specific context belongs in each dataset's directory (manifest.json, ANALYSIS_PLAN.md, etc.).

## Key Documentation

| Document | Purpose |
|----------|---------|
| `QUICKSTART.md` | **NEW DATASET INGESTION (START HERE)** |
| `DATA_ORGANIZATION.md` | Data directory structure, archival procedures |
| `src/ingest/README.md` | Ingestion pipeline stages (00-07) documentation |
| `.claude/skills/_validation_protocols.md` | Shared validation protocols for all analysis skills |

## Skills & Workflows (Claude Code)

**Location:** `.claude/skills/*.md`

| Skill | Purpose | Key Outputs |
|-------|---------|-------------|
| `explore.md` | Curiosity-driven discovery, locus exploration | findings.jsonl, neighborhood figures |
| `survey.md` | Systematic comprehensive survey | genome_profiles.tsv, comparative analysis |
| `characterize.md` | Single protein/locus characterization + batch structural analysis | Detailed functional analysis, structure predictions |
| `defense.md` | Defense system identification | CRISPR, RM, CBASS inventories |
| `prophage.md` | Prophage & viral element detection | Prophage loci, misbinned phage contigs, loci table entries |
| `metabolism.md` | Metabolic pathway reconstruction | Pathway completeness, gaps |
| `compare.md` | Cross-genome comparative analysis | Synteny, orthology |
| `literature.md` | Systematic literature/database research for functional ambiguity | Domain characterization, Foldseek hit interpretation, pathway context |
| `reviewer_2.md` | Adversarial manuscript claim verification (Phase 4) | CLAIM_VERIFICATION.jsonl, correction_queue.md |
| `atlas.md` | Exhaustive bottom-up genome-by-genome reading (subagent per batch) | atlas/inventories.jsonl, flags for specialist dispatch |
| `foldseek.md` | Structure prediction (ESM3) + Foldseek structural homology search | PDB structures, foldseek_results.tsv |
| `hydrogenase.md` | NiFe/FeFe hydrogenase validation with neighborhood curation | Validated hydrogenase inventory |

**Critical Guidelines in Skills:**
- Manuscripts: write `MANUSCRIPT.md`, render with pandoc (see § Manuscript Rendering)
- Scientific rigor: No hyperbole, evidence-calibrated claims
- Data organization: Standard directory structure for report compatibility
- MAG interpretation: "Not detected" not "absent"

**For non-Claude Code agents:** Reference the workflows in these files, but execute using your native tooling.

---

## Standard Analysis Workflow

**CRITICAL: Always include independent review at the end of analysis.**

### Orchestration Principles

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

### Five-Phase Architecture

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
├── 3. Deepen (successive rounds — see § Deepen Phase)
│   ├── Round 1: Read all findings → identify gaps → dispatch task agents
│   │   Gaps may need data analysis, literature research, or both (parallel)
│   ├── Read new findings → check for remaining gaps
│   ├── Round N: Dispatch more task agents if needed
│   └── Create report manifest + compile COMPREHENSIVE_REPORT.md
│
├── 4. Review (independent validation — see § /reviewer_2 Skill)
│   ├── Run: verify_claims.py → CLAIM_VERIFICATION.jsonl (mechanical)
│   ├── Dispatch /reviewer_2 agent → review/correction_queue.md (interpretive)
│   └── Fix all critical/meaningful corrections before proceeding
│
└── 5. Manuscript (section-by-section — see § Manuscript Writing)
    ├── Dispatch section-writing subagents (sequential by topic)
    │   Each subagent → one Results section + its MANUSCRIPT_CLAIMS.jsonl entries
    ├── Orchestrator assembles MANUSCRIPT.md from sections
    ├── Run: validate_provenance.py → PROVENANCE_AUDIT.md
    ├── Fix any provenance gaps
    ├── Run: verify_claims.py + /reviewer_2 → correction_queue.md
    ├── Fix corrections
    ├── Dispatch literature agent — this is NOT optional
    └── Apply citation corrections, re-render PDF
```

### Manuscript Citation Policy (CRITICAL — BLOCKING REQUIREMENT)

**DO NOT include literature citations in manuscript drafts from training memory.**
**DO NOT finalize or deliver a manuscript PDF without running the literature agent.**

The literature agent is a **mandatory** step in manuscript production, not an optional polish.
A manuscript without literature-agent-verified citations is incomplete and must not be
presented as finished work.

When drafting a manuscript (MANUSCRIPT.md or equivalent):
1. Write all data-derived claims normally (these come from database queries)
2. For any claim that requires a literature citation, insert a placeholder: `[CITE: topic]`
   - Example: `"...obligate syntrophs with hydrogenotrophic partners [CITE: syntrophic metabolism review]"`
3. **Run the literature agent** (`/literature manuscript`) on the draft — **this step is not optional**:
   - Replace all `[CITE: ...]` placeholders with verified citations
   - Find additional references the draft missed
   - Verify any comparative claims ("largest known", "first reported", etc.)
   - Check for missing key references (the user's own prior work, landmark studies in the field)
   - Output: `literature_citations.jsonl` (structured provenance) + `citation_report.md` (human-readable audit)
4. Apply corrections from the citation report to the manuscript
5. Only after the literature agent completes and corrections are applied should the final PDF be rendered

**Rationale:** LLM training data citations are unreliable -- wrong years, wrong authors,
confabulated papers, misattributed findings. In the Omnitrophota manuscript, training-memory
citations produced a wrong journal, wrong author list, wrong title (Moreira et al. 2021),
a factual error (GGDEF:EAL "1:1" baseline), and missed the first author's own preprint.
The literature agent performs real-time web searches and records provenance (URL, quote,
verification status) for every citation, creating an auditable trail.

### Manuscript Writing (Section-by-Section Subagents)

**The manuscript is written by separate subagents, one per Results section.** Each subagent
produces both the section text AND its claims. The orchestrator assembles the full manuscript
from the sections.

**Orchestrator workflow:**

1. Read all findings + report manifest to determine section topics
2. For each Results section, dispatch a subagent:

```python
Task(subagent_type="general-purpose",
     prompt="""Write Results section: "Energy Metabolism and Respiratory Chain"
     Source findings: survey-001, survey-005, E002, D001
     DB: data/DATASET/sharur.duckdb (for verification queries)

     Output TWO files:
     1. sections/energy_metabolism.md — the section text with [CITE: ...] placeholders
     2. sections/energy_metabolism_claims.jsonl — one claim per factual assertion

     Rules for claims:
     - Every sentence with a number, percentage, or count → quantitative claim
     - Every "X suggests Y" or "consistent with Z" → interpretive claim
     - Every "more than", "absent in", "unlike" → comparative claim
     - Use log_claim() with claim_type, source_findings, section name
     - claim_id format: MC001, MC002, ... (sequential within this section)
     """)
```

3. After all section subagents complete, assemble:
   - Concatenate section `.md` files into `MANUSCRIPT.md` (add title, abstract, intro, discussion, methods)
   - Concatenate section `_claims.jsonl` files into `MANUSCRIPT_CLAIMS.jsonl` (renumber claim_ids globally)
   - Write Introduction and Discussion as additional subagent tasks

4. Run `python scripts/validate_provenance.py --dataset data/DATASET_NAME/`
   - Fix any provenance gaps or schema drift before proceeding
5. Dispatch literature agent
6. Apply citation corrections, re-render PDF

**What gets a claim entry:**
- Every quantitative assertion (counts, percentages, sizes)
- Every comparative claim ("larger than", "more than any known", "absent in")
- Every interpretive conclusion ("consistent with syntrophy", "suggests predation")
- Methodological claims that affect interpretation ("validated by PacBio", "false positive rate of 44%")

**What does NOT need a claim entry:**
- Background/introduction sentences sourced from literature (those are tracked by citations)
- Figure captions (the figure reference in the claim entry covers this)
- Methods section boilerplate (database versions, tool parameters)

**When to run `validate_provenance.py`:**
- **NOT during Phase 4 (Review)** — MANUSCRIPT.md and MANUSCRIPT_CLAIMS.jsonl don't exist yet
- **During Phase 5**, after assembling MANUSCRIPT.md + MANUSCRIPT_CLAIMS.jsonl but before the literature agent

### Manuscript Changelog (REQUIRED)

Every dataset with a manuscript must have a `MANUSCRIPT_CHANGELOG.md` alongside `MANUSCRIPT.md`.
When editing the manuscript, **always** append an entry to the changelog documenting:

1. **What changed** — specific before/after for each edit
2. **Why it changed** — the evidence, author feedback, or review finding that motivated it
3. **Lessons learned** — if the edit reveals a systemic issue (e.g., inflated claims, missing prior work)

This applies to all edits: citation corrections, factual fixes, interpretive reframing, figure updates,
and language tightening. The changelog creates a traceable editorial history that the author can review
and that future agents can consult to understand why the manuscript says what it says.

**Format:** Number revisions sequentially (Rev 0, Rev 1, ...) with date and summary table.
See `data/omni_production/MANUSCRIPT_CHANGELOG.md` for the reference example.

Use `b.propose_hypothesis()` and `b.add_evidence()` to track analytical reasoning across sessions. Hypotheses persist in `exploration/hypotheses.json` and appear in `b.resume()` output. Use `b.render_provenance()` to generate Mermaid DAG figures for publications.

### Deepen Phase (Successive Rounds)

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

8. **Create report manifest** — organize all findings into thematic sections:
```python
manifest = {
    "dataset": "DATASET_NAME",
    "sections": [
        {
            "title": "Energy Metabolism and Hydrogenases",
            "finding_ids": ["survey-001", "survey-005", "D001"],
            "figures": ["figures/hydrogenase_neighborhood.png"],
            "narrative": "exploration/hypothesis_hydrogenase.md"
        },
        # ... one section per major topic
    ],
    "excluded_findings": [
        {"id": "survey-009", "reason": "Data quality issue, not biological"}
    ]
}
```

9. **Build report manifest PROGRAMMATICALLY** — read `findings.jsonl` files and group by category:
```python
# Read actual findings, group by category, build sections
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
# Review output, curate into report_manifest.json sections
for cat, ids in sorted(by_cat.items()):
    print(f'{cat}: {ids}')
"
```
**NEVER write report_manifest.json from memory.** Always read the actual findings files first.
Methods QC findings (false positive corrections, annotation artifacts) go in `excluded_findings`.

10. **Generate comprehensive report:**
```bash
python scripts/compile_comprehensive_report.py --dataset data/DATASET_NAME/
```

**Key principles:**
- The orchestrator decides what's worth investing in — not every finding needs deepening
- Focus on: findings likely to be headline results, claims reviewers will challenge, unknowns that could change the story
- **Successive rounds** prevent single-pass gaps: if Round 1 reveals that "DsrAB is present" but doesn't check operon context, Round 2 dispatches a neighborhood subagent
- Each round's subagents are independent and can run in parallel

### Review Agent Responsibilities

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

### Proper PFAM Function Verification

**Use predicates and names, not accession memory.**

```python
# GOOD: Search by semantic predicate
radical_sam_proteins = b.search_by_predicates(has=['radical_sam'])

# GOOD: Validate function against database name
result = b.store.execute("""
    SELECT accession, name, COUNT(*)
    FROM annotations WHERE accession = 'PF04055'
    GROUP BY accession, name
""")[0]
print(f"{result[0]} = {result[1]}")  # PF04055 = Radical_SAM

# AVOID: Relying on accession-to-function memory
```

**For pathway claims**, verify multiple components and confirm PFAM names match expected enzymes.

### Context-First Analysis Protocol

**The domain tells you the fold. The neighbors tell you the function.**

**Superfamily awareness applies to ALL HMM-based annotation sources** (PFAM, KEGG, DefenseFinder, VOGdb, CAZy). When any accession averages >10 hits per genome, it likely describes a protein fold, not a specific function. Even below that threshold, if a claimed function appears in >50% of genomes, ask whether that prevalence makes biological sense.

Before reporting any HMM-based functional claim, run two checks:
(1) **Co-annotation**: what other domains do these proteins carry? The additional domains often reveal the actual function.
(2) **Genome context**: pick 3-5 examples and examine the ±8 gene neighborhood for preserved synteny, diagnostic domain fusions, or pathway-specific neighbors that support the claimed function.

```python
# 1. Co-annotation check
co_annots = b.store.execute("""
    SELECT a2.source, a2.accession, a2.name, COUNT(DISTINCT a1.protein_id) as n
    FROM annotations a1
    JOIN annotations a2 ON a1.protein_id = a2.protein_id
    WHERE a1.accession = 'K23108' AND a2.accession != 'K23108'
    GROUP BY a2.source, a2.accession, a2.name
    ORDER BY n DESC LIMIT 10
""")

# 2. Neighborhood with ALL annotation sources
result = b.get_neighborhood(protein_id, window=5, all_annotations=True)

# 3. KEGG REST API for pathway context
# curl -s https://rest.kegg.jp/link/module/ko:K23108
```

**Claim escalation:** "Contains domain X" → "Functions as Y" (co-annotations + neighborhood) → "Encodes pathway Z" (multiple markers + co-localization) → "Phylum performs W" (all above + conservation).

---

## Canonical Tools

### Report Tiers

Not every analysis needs the full publication pipeline. Match the output to the audience:

| Tier | Audience | Requires | Does NOT require |
|------|----------|----------|------------------|
| **Collaborator report** | Colleague who owns the data | manifest → COMPREHENSIVE_REPORT → authored MANUSCRIPT.md + figures → PDF | MANUSCRIPT_CLAIMS.jsonl, verify_claims.py, /reviewer_2, /literature |
| **Publication manuscript** | Journal submission | Everything above PLUS claims tracking, provenance audit, adversarial review, literature agent | — |

**Both tiers require:**
1. `reports/report_manifest.json` built programmatically from findings (see § Deepen Phase step 9)
2. `COMPREHENSIVE_REPORT.md` generated by `compile_comprehensive_report.py`
3. COMPREHENSIVE_REPORT used as scaffold for authored `MANUSCRIPT.md`
4. Figures generated and embedded

**The collaborator report is the most common output.** Do not skip the manifest + compilation steps just because it's not a journal submission.

### Standard Figures

Every report/manuscript should include at minimum:

| Figure | Content | When generated |
|--------|---------|----------------|
| Genome overview | Heatmap: genomes × key metrics (size, proteins, annotation rate, defense, CRISPR) | After survey |
| Metabolic matrix | Presence/absence of key pathway markers across genomes | After survey |
| Defense distribution | Defense system types × genomes | After survey + CRISPR analysis |
| Neighborhood diagrams | Gene arrow diagrams for notable loci (1 per major locus-level finding) | During explore/deepen |

**Conventions:**
- Output directory: `figures/`
- Naming: `figure{N}_{descriptive_name}.png` (e.g., `figure1_genome_overview.png`)
- Format: PNG, 300 dpi, tight layout
- Style: clean professional palette, readable at printed page size
- Generation: separate subagent with DB access, after survey/explore phases complete
- Reference in MANUSCRIPT.md with relative paths: `![Caption](figures/filename.png)`

### Manuscript Pipeline

The full provenance chain from structured findings to rendered manuscript:

```
findings.jsonl (survey/ + exploration/)     ← Phases 1-3
  ↓
reports/report_manifest.json                ← Phase 3 (Deepen): organize findings into sections
  ↓
COMPREHENSIVE_REPORT.md                     ← compile_comprehensive_report.py
  ↓
sections/*.md + sections/*_claims.jsonl     ← Phase 5: section-writing subagents (one per topic)
  ↓  orchestrator assembles
MANUSCRIPT.md + MANUSCRIPT_CLAIMS.jsonl     ← concatenated + renumbered
  ↓
PROVENANCE_AUDIT.md                         ← validate_provenance.py (full chain audit)
  ↓  fix any gaps
CLAIM_VERIFICATION.jsonl                    ← verify_claims.py (number re-derivation)
  ↓
review/correction_queue.md                  ← /reviewer_2 agent (adversarial review)
  ↓  fix corrections
literature_citations.jsonl                  ← /literature agent: verified citations
  ↓
MANUSCRIPT.pdf                              ← pandoc
```

#### Rendering

```bash
export PATH="/Library/TeX/texbin:$PATH"
cd data/DATASET_NAME/
pandoc MANUSCRIPT.md -o MANUSCRIPT.pdf --pdf-engine=xelatex -V geometry:margin=1in
```

**Requirements:** pandoc + xelatex. Figure paths relative to the Markdown file directory.

#### File Structure

| File | Purpose |
|------|---------|
| `MANUSCRIPT.md` | Complete Markdown manuscript |
| `MANUSCRIPT_CLAIMS.jsonl` | Maps each claim → source findings + citations |
| `MANUSCRIPT_CHANGELOG.md` | Revision history (see § Manuscript Changelog) |
| `COMPREHENSIVE_REPORT.md` | Auto-generated from findings + manifest |
| `PROVENANCE_AUDIT.md` | Audit output from validation script |
| `CLAIM_VERIFICATION.jsonl` | Per-claim verification records from verify_claims.py |
| `REVIEW_REPORT.md` | Human-readable claim verification summary |
| `review/correction_queue.md` | Prioritized corrections from /reviewer_2 |
| `review/interpretive_review.md` | Non-quantitative claim assessment |
| `reports/report_manifest.json` | Finding → section organization |
| `literature_citations.jsonl` | Verified citations with finding links |
| `figures/` | All referenced figures as PNG files |
| `MANUSCRIPT.pdf` | Rendered output |

#### Findings Schema (Extended)

`findings.jsonl` entries support these optional fields for provenance linking:

```jsonl
{
  "id": "survey-001",
  "title": "...",
  "category": "energy_metabolism",
  "description": "...",
  "evidence": "...",
  "n_genomes": 1366,
  "provenance": { "query": "...", "raw_result": "...", "accession_verified": "...", "interpretation": "..." },
  "figures": ["figures/hydrogenase_distribution.png"],
  "related_findings": ["survey-005", "E001"],
  "phase": "survey"
}
```

Existing findings without `figures`, `related_findings`, or `phase` remain valid.

#### Figure Completion Criteria

**Locus-level findings MUST include a neighborhood figure.** If a finding describes a specific
gene, operon, or genomic region, the subagent producing it must also generate a visualization.
This is a completion criterion, not a suggestion.

| Finding type | Figure requirement |
|---|---|
| Specific locus/operon (e.g., "DsrAB operon in 9 genomes") | Neighborhood diagram for 1-3 representative examples |
| Distribution claim (e.g., "CoxI in 25/41 genomes") | Optional summary plot (heatmap, bar chart) |
| Pathway completeness (e.g., "WL pathway in 30 genomes") | Optional pathway diagram with gene presence/absence |
| Interpretive/comparative (e.g., "ecotypes are non-overlapping") | Recommended comparison figure |

**When figures are generated:**
- **Survey subagents:** Summary/distribution figures for their topic area
- **Explore subagents:** Neighborhood diagrams for every locus-level finding
- **Deepen subagents:** Figures for any finding that was flagged as "missing figure" in gap analysis
- **Manuscript subagents:** No new figures — reference existing ones from `figures/`

#### Finding ID Convention

Every finding must have a stable `id` for downstream referencing by claims and the report manifest.

| Phase | Prefix | Example | Auto-generated by |
|-------|--------|---------|-------------------|
| Survey | `survey-` | `survey-001`, `survey-015` | `log_finding()` in survey.md |
| Exploration | `E` | `E001`, `E016` | `log_finding()` in explore.md |
| Deepen | `D` | `D001`, `D010` | `log_finding()` in explore.md (with `finding_id="D001"`) |

IDs auto-increment within each phase file. The `log_finding()` helpers in survey.md and
explore.md generate IDs automatically. For Deepen findings appended to `exploration/findings.jsonl`,
pass an explicit `finding_id="DNNN"` to distinguish them from exploration findings.

**A finding without an `id` cannot be referenced by claims or the report manifest.** The
compilation script shows `[???]` for missing IDs.

#### Report Manifest (`reports/report_manifest.json`)

See § Deepen Phase for schema and example.

#### Manuscript Claims (`MANUSCRIPT_CLAIMS.jsonl`)

Maps each factual claim in the manuscript to source findings and citations:

```jsonl
{
  "claim_id": "MC001",
  "section": "Abstract",
  "claim_text": "Energy-conserving Group 4 NiFe hydrogenases are present in 74.6% of genomes",
  "source_findings": ["survey-001"],
  "source_citations": [],
  "figures": [],
  "claim_type": "quantitative"
}
```

**Required:** `claim_id`, `claim_text`, `source_findings`, `claim_type`
**Recommended:** `section` (manuscript section where the claim appears)
**Optional:** `source_citations`, `figures`, `review_status`
**claim_type:** `quantitative` | `interpretive` | `comparative` | `methodological`

#### Literature Citations (Enhanced)

`literature_citations.jsonl` entries can link to findings they support:

```jsonl
{
  "citation_id": "touchon_2016",
  "supports_findings": ["D008", "D009"],
  "... existing fields ..."
}
```

#### Scripts

- **Compile report:** `python scripts/compile_comprehensive_report.py --dataset data/DATASET_NAME/` → `COMPREHENSIVE_REPORT.md`
- **Validate provenance:** `python scripts/validate_provenance.py --dataset data/DATASET_NAME/` → `PROVENANCE_AUDIT.md`
- **Verify claims:** `python scripts/verify_claims.py --dataset data/DATASET_NAME/ --auto-extract` → `CLAIM_VERIFICATION.jsonl` + `REVIEW_REPORT.md`
- **Validate defense systems:** `python scripts/validate_defense_systems.py data/DATASET_NAME/` → `defense_systems` table + `defensefinder_system` annotations
- **Validate secretion systems:** `python scripts/validate_secretion_systems.py data/DATASET_NAME/` → `secretion_systems` table + `txsscan_system` annotations
- Legacy report generators (`generate_paper_report.py`, etc.) are hardcoded for specific datasets — do NOT use for new manuscripts.

### Hydrogenase Classification
**Script:** `scripts/classify_hydrogenases.py`
**Requires:** HydDB HMMs via Astra, DIAMOND database (`data/reference/hyddb/HydDB_all.dmnd`)
**Output:** Subgroup-level classification (NiFe Group 1-4, FeFe A-C)
**Note:** Pipeline classifies all HydDB hits and assigns subgroup predicates. Hits lacking PFAM corroboration are tagged `hyddb_needs_curation` for agent neighborhood curation — see skill specs.

### CAZyme Classification (dbCAN 3-tool consensus)
**Script:** `scripts/classify_cazymes.py`
**Requires:** `data/dbcan_db/` with CAZy.dmnd, dbCAN.hmm, dbCAN-sub.hmm; dbCAN installed in Astra
**Output:** Consensus CAZyme annotations (`source='cazy'`, family-level: GH, GT, PL, CE, CBM, AA)
**Method:** DIAMOND (e-value ≤1e-18) + Astra/HMMER vs dbCAN.hmm (i-evalue ≤1e-15) + Astra/HMMER vs dbCAN-sub.hmm (i-evalue ≤1e-15), consensus ≥2 tools. Falls back to 2-tool if dbCAN-sub absent.
**Caching:** Reuses DIAMOND results from `cazyme_classification.tsv` and Astra HMM results from `annotations/dbCAN_hits_df.tsv`. dbCAN-sub runs only on proteins hit by ≥1 other tool.
**Pipeline:** Integrated into stage 07 (`_classify_cazymes()`) after predicates + hydrogenase classification. Stage 04 runs Astra dbCAN automatically.

### Embedding Visualization
**Script:** `scripts/visualize_embeddings.py`
**Requires:** `umap-learn`, `plotly` (interactive) or `matplotlib` (static)
**Usage:** `python scripts/visualize_embeddings.py --db data/DATASET/sharur.duckdb --output figures/umap.html --color-by genome`
**Color modes:** `genome` (by bin_id), `predicate` (highlight specific predicate), `annotation` (by PFAM/KEGG name)

### Local Foldseek
**Binary:** Auto-detected via `which foldseek`
**Database path:** `~/.foldseek/{db_name}/{db_name}` (e.g. `~/.foldseek/pdb100/pdb100`)
**Behavior:** `search_foldseek()` tries local binary first (`prefer_local=True` by default), falls back to web API for databases not installed locally. Local search is faster and has no rate limits.

---

## Check for Functional Detail (IMPORTANT)

**Don't stop at generic predicates — drill into subgroup-level detail when available.**

### Hydrogenases
If you see `hydrogenase` or `hydrogen_metabolism` predicates, **check for subgroup predicates**:
- `nife_group1`–`nife_group4`, `mbh_hydrogenase`, `ech_hydrogenase`
- `fefe_groupA`, `fefe_groupB`, `fefe_groupC`
- Group 3 vs Group 4 reveals uptake vs evolution, respiratory vs fermentative

### CRISPR Systems
`crispr_associated` or `cas_domain` → check `type_i_crispr`, `type_ii_crispr`, `type_iii_crispr`, effectors (`cas3`, `cas9`, `cas10`), and `loci` table for CRISPR arrays.

### Defense Systems
`defense_system` → check DefenseFinder source annotations for specific system types (RM, CBASS, BREX, DISARM, etc.)

### CAZy Enzymes
`carbohydrate_active` → check `cazy:GH5`, `cazy:GT2`, etc. GH families reveal substrate specificity.

**Why this matters:** "41 genomes have hydrogenases" tells you nothing. "5 have Group 4 energy-conserving, 6 have Group 3 F420-reducing" reveals metabolic diversity.

---

## Key Files

| File | Purpose |
|------|---------|
| `PREDICATE_PLAN.md` | Predicate system design and status |
| `sharur/predicates/vocabulary.py` | All predicate definitions |
| `sharur/predicates/generator.py` | Main predicate computation |
| `sharur/predicates/mappings/` | PFAM/KEGG/CAZy/VOGdb → predicate mappings |
| `sharur/predicates/topology.py` | TM helix prediction wrapper (pyTMHMM) |
| `sharur/operators/structure.py` | ESM3 structure prediction |
| `sharur/operators/foldseek.py` | Foldseek structural homology search |
| `sharur/operators/manifest.py` | Analysis manifest for session continuity |
| `sharur/core/hypothesis_registry.py` | Persistent cross-session hypothesis store |
| `sharur/core/provenance_renderer.py` | Mermaid DAG renderer for provenance figures |
| `sharur/reports/template.py` | PDF report generation with themes |

## Analysis Manifest System

Each dataset has a `manifest.json` for session continuity. Key APIs: `b.resume()` (status overview), `b.manifest.log_session(phase, note)`, `b.manifest.save()`. Structures and figures are auto-tracked. Migration: `python scripts/migrate_to_manifest.py data/my_dataset/`

## Quick Reference: Hypothesis Tracking & Provenance

```python
b = Sharur("data/YOUR_DATASET/sharur.duckdb")
h = b.propose_hypothesis("Group 4 NiFe hydrogenases are energy-conserving")
b.add_evidence(h.hypothesis_id, "NiFe Group 4 survey", "12/41 genomes", True, 0.8)
print(b.hypothesis_summary())  # One-line overview
# Provenance: e1 = b.log_provenance("query", "result"); e2 = b.log_provenance("next", "result", parent_ids=[e1.entry_id])
# DAG figure: b.render_provenance(title="Analysis DAG", output_path="figures/provenance.mermaid")
```

Hypotheses persist in `exploration/hypotheses.json` across sessions. `b.resume()` shows active hypotheses.

## Quick Reference: Structure Prediction & Foldseek

**Reference code:** `sharur/operators/structure.py` (implementation), `sharur/operators/__init__.py` (Sharur class integration, lines 541-615)

**Available methods:**

### 1. For database proteins (standard workflow)
```python
b = Sharur("data/DATASET/sharur.duckdb")
result = b.predict_structure("protein_id", output_path="structures/protein.pdb")  # requires ESM_API_KEY
# Auto-updates manifest with structure prediction
# Max length: 1024 aa (ESM3 API limit)
```

### 2. For arbitrary sequences (no database lookup)
```python
from sharur.operators.structure import predict_structure_from_sequence
result = predict_structure_from_sequence(
    sequence="MKVL...",  # raw AA sequence
    output_path="output.pdb",
    name="my_protein"  # optional identifier
)
# Use for external sequences, synthetic constructs, or protein fragments
```

### 3. Batch processing
```python
result = b.batch_predict_structures(
    protein_ids=["id1", "id2", "id3"],
    output_dir="structures/",
    max_length=1024  # skip proteins longer than this
)
# Auto-generates filenames: {protein_id}.pdb
```

**Foldseek structural homology search:**
```python
hits = b.search_foldseek("structures/protein.pdb", databases=["afdb50", "pdb100"], top_k=10)
hits = b.search_foldseek_for_protein("protein_id")  # convenience: uses existing PDB if available
b.list_foldseek_databases()  # afdb50, afdb-swissprot, pdb100
```

**TM-score:** >0.5 = similar fold; >0.7 = high confidence homology. pdb100 = experimental; afdb50 = AlphaFold predictions.

**Implementation details:** See `sharur/operators/structure.py` lines 34-318 for full ESM3 API integration, error handling, and manifest tracking.

## Astra Annotation Pipeline

Astra manages pre-installed HMM databases for annotation searches.

**Location:** `~/astra/` (source), installed via pyenv shim

**Installed databases** (`~/.config/Astra/`): PFAM, KOFAM, VOGdb, HydDB, DefenseFinder, CRISPRCasFinder, CANT-HYD, dbCAN, padloc, TXSScan

```bash
astra search --installed_hmms DefenseFinder --threads 12 \
    --prot_in <directory_with_faa_files> \
    --outdir <output_directory> \
    --cut_ga
```

**Notes:**
- `--prot_in` expects a **directory** containing `.faa` files, not a single file
- Output: `<outdir>/<database>_hits_df.tsv` (tab-separated)
- `--cut_ga` uses gathering thresholds (recommended)
- For single files: `mkdir source/ && cp proteins.faa source/`
- `--write_macsyfinder` produces MacSyFinder-compatible hmmsearch output for co-localization validation (auto-added for DefenseFinder and TXSScan in stage 04)

### Secretion System Identification (TXSScan)
**Script:** `scripts/validate_secretion_systems.py`
**Requires:** TXSScan HMMs via Astra, TXSScan models in MacSyFinder (`macsyfinder --install-models TXSScan`)
**Pipeline:** Non-default — add `-d TXSScan` to stage 04. Stage 07 auto-validates via MacSyFinder `--previous-run` when `txsscan_results/macsyfinder_compat/` exists.
**Output:** `secretion_systems` table + `txsscan_system` annotations in DuckDB. Predicates: `secretion_system`, `type_I_secretion`–`type_IX_secretion`, `type_IV_pilus`, `tad_pilus`, `flagellar`, `msh_pilus`, `competence`.
**Systems detected:** T1SS–T9SS, T4P, Tad, Flagellum, MSH, ComM (280 HMM profiles).

## Testing

```bash
python -m pytest tests/test_operators/test_predicate_generator.py -v
python -c "from sharur.predicates.vocabulary import ALL_PREDICATES, list_categories; print(f'Total: {len(ALL_PREDICATES)}'); print(f'Categories: {list_categories()}')"
python -c "from sharur.predicates.mappings.pfam_map import PFAM_TO_PREDICATES, PFAM_PATTERNS; print(f'Direct: {len(PFAM_TO_PREDICATES)}, Patterns: {len(PFAM_PATTERNS)}')"
```

---

## Predicate System Design Principles

- **Gene-level vs Locus-level**: Components get gene-level tags; system calls require clustering
- **Tiered predicates**: Generic domain → specific system (e.g., `cas_domain` → `crispr_associated`)
- **Biological accuracy**: Ferredoxins are electron carriers, not substrate indicators
- **Pathway completeness**: Single markers don't prove pathway function
- **Carbon fixation**: RuBisCO without PRK is NOT Calvin cycle (likely nucleotide salvage). `calvin_cycle` requires PRK.
- **Methanogenesis**: Only MCR complex triggers; H4MPT enzymes alone are insufficient
- **Topology**: pyTMHMM integration (optional: `pip install sharur[topology]`). TMbed planned for signal peptide prediction.

### PFAM Mapping Scaling

To avoid growing `pfam_map.py` indefinitely:
- Extension file: `data/reference/pfam_predicate_map.tsv` (PFAM_ID_OR_NAME \t pred1,pred2,... \t optional_note)
- `pfam_map.py` loads and merges this at import time
- Bulk expansions go in the TSV, curated core mappings stay in Python

---

## Code Standards & Best Practices

### Data Integrity Rule

**Never write structured data from memory.** When producing JSON, JSONL, TSV, or any output that
references findings, database entries, annotations, or other generated artifacts, ALWAYS read the
source data first. This applies to:
- `report_manifest.json` — read `findings.jsonl` to get actual finding IDs and categories
- `MANUSCRIPT_CLAIMS.jsonl` — verify numbers against the database before writing claims
- Any summary table — query the database, don't reconstruct from conversation history

If you need to group, filter, or organize existing data, write a script that reads and processes
the actual files. Do not type structured output from what you remember seeing earlier.

### Manuscript & Report Writing

**Standard output:** Each dataset produces `MANUSCRIPT.md` → `MANUSCRIPT.pdf` via pandoc (see § Manuscript Rendering under Canonical Tools).

**Do NOT:**
- Create new Python PDF generators — write Markdown and use pandoc
- Use hardcoded report scripts (`generate_paper_report.py`, etc.) for new datasets
- Create versioned or dated PDF filenames — always `MANUSCRIPT.pdf`

**Writing guidelines:**
- Avoid sensationalized framing — no "Mystery - Resolved" sections
- TnpB classification — report as transposases unless clearly within a CRISPR/Cas locus
- Use neutral, factual language
- Figures go in `figures/` subdirectory, referenced with relative paths in Markdown
- Number figures sequentially; when inserting a new figure, renumber all subsequent figures and their cross-references

### Visualization

**When modifying existing plotting code, make targeted edits only.** Read the script first, understand the layout, change only what was requested.

**Use Sharur's visualization operators:**
```python
b.visualize_neighborhood(protein_id, window=12, output_path="output.png")
```

**Multi-source locus diagrams:**
```bash
python scripts/plot_locus_multisource.py \
    --db data/dataset/sharur.duckdb \
    --protein PROTEIN_ID \
    --window 12 \
    --output figures/locus.png \
    --title "Custom Title"
```

Features: Multi-source annotation priority (Foldseek > DefenseFinder > PADLOC > PFAM/KEGG/VOGdb), clean label boxes, gene numbers below track, absolute genome coordinates, CRISPR array detection.

**Custom implementations** should use `dna_features_viewer` with `annotate_inline=False`, gene numbers below the track, absolute coordinates on x-axis. Color by functional category, not annotation source. For ambiguous annotations (Cas12f vs TnpB), use honest labels.

### Database Queries

```python
# Use 'name' column for domain names, not 'annotation_id'
# Use 'score' column, not 'bitscore'
# Always COUNT(DISTINCT protein_id) for protein counts — repeat domains inflate COUNT(*)

# Prefer Sharur operators over raw DuckDB:
result = b.search_by_predicates(has=["unannotated", "giant"]); proteins = result._raw
b.get_neighborhood(protein_id, window=10)
b.get_neighborhood(protein_id, window=5, all_annotations=True)
```

### External Data Lookup

**Research PDB hits via WebFetch** — don't guess protein functions:
```python
WebFetch("https://www.rcsb.org/structure/5fms", "What is the protein function?")
```

---

## Biological Interpretation Guidelines

### MAG Interpretation

**Cardinal Rule: Absence of evidence ≠ evidence of absence**

MAGs are inherently incomplete. A missing gene does NOT mean the organism lacks it.

| Contigs | Fragmentation | How to Interpret Absence |
|---------|---------------|--------------------------|
| <50 | Low | Reasonably reliable |
| 50-200 | Moderate | Include caveats |
| >200 | High | "Not detected" only |
| genes/contig <5 | Very high | Many genes likely missing |

**Language:** "No hydrogenases were detected in this MAG (N contigs)" — NOT "Genome X lacks hydrogenases."

**Comparative claims:** Before saying "A has X but B doesn't", verify B isn't just more fragmented.

### Giant Protein Annotation Recovery

**Standard PFAM bitscore cutoffs are LENGTH-BIASED.** Giant proteins (>1000 aa) often show zero hits. E-values are NOT length-biased — use them instead.

```bash
# For proteins >1000 aa with 0 standard hits
hmmsearch --domE 1e-5 ~/.config/Astra/PFAM/Pfam-A.hmm giant_protein.faa
```

**When to apply:** >1000 aa with 0 hits, >2000 aa with <3 hits, before reporting any giant protein as "unannotated."

**E-value interpretation:** ≤1e-10 high confidence; 1e-10 to 1e-5 moderate (meaningful for giants); 1e-5 to 1e-3 weak (check for repeat patterns).

**Common giant protein architectures:** Big_13/Big_3_3/Big_8 (adhesins), TPR/HEAT/WD40 (scaffolds), Beta_helix (autotransporters), ANK (signaling), Cadherin/FN3 (adhesion).

### Hydrogenase Classification

**Primary source: HydDB HMMs** — use HydDB over PFAM for hydrogenase typing.

#### Pipeline + Agent Curation

**Pipeline (`classify_hydrogenases.py`):** Classifies all HydDB hits via DIAMOND, assigns subgroup predicates (Groups 1-4, FeFe A-C). Hits with PFAM corroboration (PF00374 for NiFe, PF02906/PF02256 for FeFe) are high-confidence. Hits without are tagged `hyddb_needs_curation`.

**Agent curation (neighborhood-based):**
For `hyddb_needs_curation` hits, check ±8 gene neighborhood:

| Evidence | KEGG KOs | Verdict |
|----------|----------|---------|
| Hyf (hydrogenase-4) | K12136-K12145 | Real Group 4f |
| Hyc (formate hydrogenlyase) | K15828-K15833 | Real Group 4 |
| Maturation (HypA-F, HycI) | K04651-K04656, K03605 | Real hydrogenase |
| Complex I (nuoA-N) | K00330-K00343 | False positive |
| No evidence either way | — | Reject conservatively |

#### Key PFAM Domains

| Domain | Accession | Meaning |
|--------|-----------|---------|
| NiFeSe_Hases | PF00374 | NiFe large subunit (Groups 1-3 only) |
| Fe_hyd_lg_C | PF02906 | FeFe hydrogenase large subunit |
| Ni_hydr_CYTB | PF01292 | Cytochrome b — NOT a hydrogenase |
| Complex1_49kDa | PF00346 | NADH dehydrogenase — FALSE POSITIVE without neighborhood |

#### Subgroups

| NiFe Group | Function | Key Predicates |
|------------|----------|----------------|
| 1a-1l | Respiratory uptake | `nife_group1`, `uptake_hydrogenase` |
| 2a-2e | Cytoplasmic H2 sensors | `nife_group2`, `h2_sensor` |
| 3a-3d | Bidirectional, cofactor-coupled | `nife_group3`, `bidirectional_hydrogenase` |
| 4a-4i | Energy-conserving, H2-evolving | `nife_group4`, `mbh_hydrogenase`, `ech_hydrogenase` |

| FeFe Group | Function | Key Predicates |
|------------|----------|----------------|
| A1-A4 | Monomeric fermentative | `fefe_groupA` |
| B | Electron-bifurcating | `fefe_groupB`, `bifurcating_hydrogenase` |
| C1-C3 | Sensory/regulatory | `fefe_groupC` |

### Cytochrome Validation

**Always validate respiratory system claims, especially "lacks cytochromes."**

Three detection methods:
1. **Predicates:** `b.search_by_predicates(has=["cytochrome"])`
2. **Annotations:** `WHERE LOWER(name) LIKE '%cytochrome%'`
3. **Sequence motifs:** CxxCH heme-binding motif (`re.search(r'C..CH', seq)`)

If all three find nothing AND genome is high-quality (<50 contigs), absence is plausible. Otherwise, state uncertainty.

---

## Scientific Rigor

**Under-promise, over-deliver.** Rigorous, conservative science is more impressive than hyperbolic claims.

**Forbidden language:** "confirms/proves/demonstrates" (unless truly definitive), "unprecedented/first ever/groundbreaking", "paradigm-shifting/revolutionary", "Nature/Science-tier discovery"

**Required language:** "suggests/indicates/supports", "consistent with/compatible with", "to our knowledge" (after verification), "provides evidence for"

| Evidence | Appropriate Language |
|----------|---------------------|
| Sequence annotation | "annotated as", "contains domain" |
| Structural homology | "structural similarity suggests" |
| Genomic context | "co-located with", "may indicate" |
| Literature | "similar to", "consistent with" |
| Experimental | "demonstrates", "confirms" (OK!) |

**Common errors:** Domain presence ≠ function proof. MAG absence ≠ biological absence. Single marker ≠ pathway presence. Transposase ≠ mobile element proof (could be Cas12f). "First in analysis" ≠ "first ever."

---

## Subagent Strategies

### Sub-Agent Context Injection (CRITICAL)

**Subagents receive ONLY the prompt parameter.** They do NOT automatically get CLAUDE.md, skill specs (`.claude/skills/*.md`), conversation history, or memory files. The quality of subagent work depends entirely on what context the orchestrator includes in the prompt.

**When dispatching analysis subagents, include in the prompt:**
1. **Database path** and dataset context (genome count, phylum, environment)
2. **Domain documentation path** — tell subagents about the Obsidian docs vault so they can look up protocols on demand (see below)
3. **Relevant skill spec content** — read the appropriate `.claude/skills/*.md` and include key sections: finding schema, `log_finding()` helper, validation protocols, output specifications
4. **Superfamily awareness rule** — always include the Context-First Analysis Protocol (co-annotation + neighborhood checks before functional claims)
5. **Scientific rigor guidelines** — MAG interpretation ("not detected" not "absent"), claim escalation, forbidden language
6. **Database query patterns** — column names (`name` not `annotation_id`, `score` not `bitscore`), `COUNT(DISTINCT protein_id)`, no correlated subqueries
7. **Output specification** — what files to produce, what format, completion criteria

**Minimum boilerplate for any analysis subagent prompt:**
```
DB: data/DATASET/sharur.duckdb
Import: from sharur.operators import Sharur; b = Sharur("data/DATASET/sharur.duckdb")
Columns: 'name' (not annotation_id), 'score' (not bitscore)
Counts: always COUNT(DISTINCT protein_id) — repeat domains inflate COUNT(*)
MAG caveat: "Not detected" not "absent"
Findings: JSONL with id, title, category, description, evidence, n_genomes, provenance

## Domain Documentation (READ ON DEMAND)
When you encounter a domain-specific situation, look up the relevant protocol doc before acting:
Docs path: /Users/jacob/Documents/Obsidian Vault/sharur-docs/
Key docs:
  - hydrogenase-classification.md — HydDB curation, Complex I FP detection, neighborhood validation
  - defense-system-validation.md — superfamily FP filtering, prevalence sanity checks
  - giant-protein-recovery.md — E-value recovery for >1000 aa, ESM3/Foldseek workflow
  - context-first-protocol.md — co-annotation validation, claim escalation ladder
  - mag-quality-interpretation.md — fragmentation checks, absence claim language
  - query-patterns.md — DuckDB query templates, JOIN patterns, performance tips
  - neighborhood.md — get_neighborhood() usage, what to look for
  - predicates-overview.md — predicate system, search_by_predicates()
Use the Read tool to load any doc when you need its protocol. Don't guess — look it up.
```

**For survey/explore subagents**, read the full skill spec and paste the protocol sections, finding schema, and validation rules into the prompt. These specs contain ~1000+ lines of battle-tested methodology.

### Sub-Agent Protocol
- Sub-agents CAN spawn further sub-agents (recursive dispatch) for specialist tasks (literature, foldseek, hydrogenase curation)
- Each sub-agent produces a discrete output (markdown report, JSONL findings)
- Parent agent synthesizes outputs from all sub-agents
- **DuckDB write locks** — run subagents SEQUENTIALLY if they write to the database

### Practical Tips
- **Parallel genome browser agents** work well (quarters or groups)
- **JSONL for findings** — easy to append, merge, and process
- **Check database schema** before writing queries (`DESCRIBE table_name`)
- **Research external data** via WebFetch — don't guess PDB functions
- **Don't create new visualization code** when operators exist

---

## Large Dataset Performance (>50k proteins)

**Rules:**
1. **Never** `b.search_proteins(query="")` on large datasets — use `b.store.execute("SELECT COUNT(*) FROM proteins")[0][0]` or specific predicates
2. **Always** check result size before iterating (use `result._raw` to get the underlying list)
3. **Always** limit iteration (e.g., `for pid in proteins[:10]` after `proteins = result._raw`)
4. **Combine** specific predicates (AND/OR) to narrow results
5. **Aggregate** in SQL, not Python loops
6. **Limit** to 20-30 visualizations per analysis
7. **If query takes >5 seconds**, stop and refine
8. **NEVER use correlated subqueries** — rewrite as JOINs (see below)

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
    SELECT COUNT(*) as total, AVG(sequence_length) as avg_size, MAX(sequence_length) as max_size
    FROM proteins
""")[0]
```

### DuckDB Query Patterns for Large Datasets

**CRITICAL: Correlated subqueries destroy performance on >1M row tables.** Never put a subquery inside WHERE that references the outer row. Instead:

1. **Materialize the seed set** in a CTE with `contig_id` and `gene_index`
2. **JOIN** to find neighbors (same contig, gene_index ± window)
3. **JOIN** to annotations/predicates for neighbor features
4. **Use GROUP BY** for per-genome cross-tabs, not Python loops over genomes
5. **Use CASE WHEN** inside COUNT(DISTINCT ...) for multi-marker pivots

---

## Standard Directory Structure

```
data/{dataset_name}/
├── sharur.duckdb                # Core database
├── manifest.json               # Analysis state
├── source/                     # Input files (.faa)
├── annotations/                # Annotation results (pfam.tsv, kegg.tsv, etc.)
├── embeddings/                 # ESM2 embeddings (LanceDB)
├── structures/                 # ESM3 PDBs + Foldseek results
├── exploration/                # Exploration outputs
├── survey/                     # Survey outputs
├── reports/                    # Generated reports
└── figures/                    # Top-level figures
```

## TODO

- [ ] **Bake superfamily validation deeper into `/survey`** — known traps (HydDB, Mokosh/BREX, Radical_SAM), co-annotation checks, automatic specialist dispatch
