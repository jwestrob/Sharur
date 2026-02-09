# Agent Knowledge Base (CLAUDE.md / AGENTS.md)

**Audience:** All AI agents working on Sharur (Claude Code, Codex, Antigravity, etc.)

This file provides shared context, canonical tools, and best practices for agent-driven metagenomic analysis.

**Quick Start:**
- **Data organization:** `DATA_ORGANIZATION.md`
- **Canonical report generator:** `scripts/generate_viral_genome_report.py`
- **Skills (Claude Code):** `.claude/skills/*.md`
- **Validation protocols:** `.claude/skills/_validation_protocols.md`

## Project Overview

Sharur is an agent-driven metagenomic exploration system. It's a data plane that makes large metagenomic datasets navigable by AI agents.

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
| `metabolism.md` | Metabolic pathway reconstruction | Pathway completeness, gaps |
| `compare.md` | Cross-genome comparative analysis | Synteny, orthology |
| `literature.md` | Systematic literature/database research for functional ambiguity | Domain characterization, Foldseek hit interpretation, pathway context |

**Critical Guidelines in Skills:**
- Use canonical report generator: `scripts/generate_viral_genome_report.py`
- Scientific rigor: No hyperbole, evidence-calibrated claims
- Data organization: Standard directory structure for report compatibility
- MAG interpretation: "Not detected" not "absent"

**For non-Claude Code agents:** Reference the workflows in these files, but execute using your native tooling.

---

## Standard Analysis Workflow

**CRITICAL: Always include independent review at the end of analysis.**

### Five-Phase Architecture

```
Coordinator
│
├── 1. Survey (systematic coverage)
│   └── Spawns topic subagents: metabolic, defense, surface, etc.
│
├── 2. Explore (hypothesis-driven discovery)
│   └── Reads survey → generates hypotheses → spawns testing subagents
│
├── 3. Deepen (targeted follow-up on findings)
│   └── Coordinator reads all findings, identifies gaps, dispatches specialists
│
├── 4. Review (independent validation)
│   └── Validates claims, catches errors, assesses publication readiness
│
└── 5. Literature & Manuscript (MANDATORY before any manuscript is finalized)
    ├── Draft manuscript with [CITE: topic] placeholders (NO training-memory citations)
    ├── Run literature agent (`/literature manuscript`) — this is NOT optional
    ├── Literature agent verifies/finds citations with provenance
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

### Deepen Phase (Coordinator-Driven)

**Goal**: Turn preliminary findings into well-supported, publication-ready results. Unlike Explore (which follows its own curiosity), Deepen is the coordinator making deliberate decisions about what needs more evidence.

**Workflow:**

1. **Read all findings** from Survey and Explore (`findings.jsonl`, `genome_profiles.tsv`, `exploration_state.json`)
2. **Gap analysis** — identify:
   - Findings with weak evidence (single annotation, no neighborhood context)
   - Unknown proteins that deserve characterization
   - Metabolic claims missing pathway completeness checks
   - Defense systems needing detailed inventory
   - Novel findings lacking literature context
3. **Dispatch targeted subagents** based on gaps:

```python
# Characterize key unknowns — structure + foldseek + literature
Task(subagent_type="general-purpose",
     prompt="Characterize protein X: predict structure, run foldseek, search literature...",
     run_in_background=True)

# Verify metabolic pathway claims
Task(subagent_type="general-purpose",
     prompt="Check completeness of Wood-Ljungdahl pathway across all genomes...",
     run_in_background=True)

# Literature search on novel findings
Task(subagent_type="general-purpose",
     prompt="Search literature for precedents of [finding]...",
     run_in_background=True)

# Defense system deep-dive
Task(subagent_type="general-purpose",
     prompt="Complete defense system inventory: CRISPR typing, RM specificity...",
     run_in_background=True)
```

4. **Synthesize** — integrate subagent results back into findings, update hypotheses, write narrative

**Key principle:** The coordinator decides what's worth investing in. Not every finding needs deepening. Deepen is the "what would a reviewer ask?" step — focus effort on:
- Findings likely to be headline results
- Claims that reviewers will challenge
- Unknowns that could change the story

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
""").fetchone()
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
""").fetchall()

# 2. Neighborhood with ALL annotation sources
result = b.get_neighborhood(protein_id, window=5, all_annotations=True)

# 3. KEGG REST API for pathway context
# curl -s https://rest.kegg.jp/link/module/ko:K23108
```

**Claim escalation:** "Contains domain X" → "Functions as Y" (co-annotations + neighborhood) → "Encodes pathway Z" (multiple markers + co-localization) → "Phylum performs W" (all above + conservation).

---

## Canonical Tools

### Report Generation
**Script:** `scripts/generate_viral_genome_report.py` (37KB, proven generator)
**Use for:** Comprehensive viral genome reports with locus diagrams
**Do NOT:** Create new report generators from scratch

### Hydrogenase Classification
**Script:** `scripts/classify_hydrogenases.py`
**Requires:** HydDB HMMs via Astra, DIAMOND database (`data/reference/hyddb/HydDB_all.dmnd`)
**Output:** Subgroup-level classification (NiFe Group 1-4, FeFe A-C)
**Note:** Pipeline classifies all HydDB hits and assigns subgroup predicates. Hits lacking PFAM corroboration are tagged `hyddb_needs_curation` for agent neighborhood curation — see skill specs.

### Embedding Visualization
**Script:** `scripts/visualize_embeddings.py`
**Requires:** `umap-learn`, `plotly` (interactive) or `matplotlib` (static)
**Usage:** `python scripts/visualize_embeddings.py --db data/DATASET/sharur.duckdb --output figures/umap.html --color-by genome`
**Color modes:** `genome` (by bin_id), `predicate` (highlight specific predicate), `annotation` (by PFAM/KEGG name)

### Local Foldseek
**Binary:** Auto-detected via `which foldseek`
**Database path:** `~/.foldseek/{db_name}/{db_name}` (e.g. `~/.foldseek/pdb100/pdb100`)
**Behavior:** `search_foldseek()` tries local binary first (`prefer_local=True` by default), falls back to web API for databases not installed locally. Local search is faster and has no rate limits.

### Visualization
**Sharur operators:** `b.visualize_neighborhood()`, `b.visualize_domains()`
**Multi-source locus script:** `scripts/plot_locus_multisource.py`
**Do NOT:** Create matplotlib code from scratch when operators exist

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

Each dataset has a `manifest.json` for session continuity:

```python
from sharur.operators import Sharur
b = Sharur("data/YOUR_DATASET/sharur.duckdb")

print(b.resume())  # Status, findings, structures, hypotheses, recent activity

# Auto-tracking: structures and figures are automatically recorded
b.predict_structure("protein_id")
b.visualize_neighborhood("protein_id", title="My Figure", legend="Caption text")

# Manual updates
b.manifest.log_session("exploration", "Completed CRISPR locus analysis")
b.manifest.save()
```

**Migration:** `python scripts/migrate_to_manifest.py data/my_dataset/`

**Report from manifest:**
```python
from sharur.reports import generate_report_from_manifest, SharurReport
generate_report_from_manifest("data/my_dataset/manifest.json", "output.pdf", theme="viral")
```

## Quick Reference: Hypothesis Tracking & Provenance

```python
b = Sharur("data/YOUR_DATASET/sharur.duckdb")

# Propose a hypothesis (persists to exploration/hypotheses.json)
h = b.propose_hypothesis("Group 4 NiFe hydrogenases are energy-conserving")

# Add evidence after analysis
b.add_evidence(h.hypothesis_id, "NiFe Group 4 survey", "12/41 genomes", True, 0.8)
b.add_evidence(h.hypothesis_id, "Neighborhood check", "Hyf operon present", True, 0.9)

# Check state
print(b.hypothesis_summary())    # One-line-per-hypothesis overview
b.list_hypotheses()              # Full Hypothesis objects

# Explicit provenance logging with parent chaining
e1 = b.log_provenance("Count hydrogenases", "42 found")
e2 = b.log_provenance("Check neighborhoods", "12 with Hyf", parent_ids=[e1.entry_id])

# Render provenance DAG as Mermaid diagram
mermaid = b.render_provenance(title="Analysis DAG", output_path="figures/provenance.mermaid")
```

- `b.resume()` automatically shows active hypotheses
- Hypotheses persist across sessions and subagent runs
- `add_evidence()` accepts UUID or string hypothesis_id
- `log_provenance()` accepts UUID or string parent_ids

## Quick Reference: Structure Prediction & Foldseek

```python
b = Sharur("data/YOUR_DATASET/sharur.duckdb")

# Predict structure (requires ESM_API_KEY env var)
result = b.predict_structure("protein_id", output_path="structures/protein.pdb")

# Search structure against databases
b.list_foldseek_databases()  # afdb50, afdb-swissprot, pdb100
hits = b.search_foldseek("structures/protein.pdb", databases=["afdb50", "pdb100"], top_k=10)

# Convenience: search for protein (uses existing PDB if available)
hits = b.search_foldseek_for_protein("protein_id")
```

**Interpreting results:**
- TM-score > 0.5: Similar fold, likely related function
- TM-score > 0.7: High confidence homology
- pdb100 hits: Known structures with experimental validation
- afdb50 hits: AlphaFold predictions, check UniProt

### All-Atom Folding with Ligands (ESM3 Forge API)

**Status:** Not yet implemented in Sharur operators. The Forge API supports all-atom folding with ligands via `/api/v1/fold_all_atom` (proteins, DNA, RNA, ligands, covalent bonds). See API docs for `ProteinInput`, `LigandInput`, and `covalent_bonds` parameters.

**Alternatives:** AlphaFold3, RoseTTAFold All-Atom, AlphaFill (post-hoc enrichment).

## Astra Annotation Pipeline

Astra manages pre-installed HMM databases for annotation searches.

**Location:** `~/astra/` (source), installed via pyenv shim

**Installed databases** (`~/.config/Astra/`): PFAM, KOFAM, VOGdb, HydDB, DefenseFinder, CRISPRCasFinder, CANT-HYD

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

### Report Generation

**Two report scripts — use the right one:**

| Script | Purpose | When to use |
|--------|---------|-------------|
| `scripts/generate_exploration_report.py` | Auto-generated survey/exploration report | After `/survey` or `/explore` completes |
| `scripts/generate_paper_report.py` | Hand-curated paper-style report | Publication-quality output with narrative sections |

**Never create ad-hoc PDF generation** — always adapt existing report generators.

**Report filename convention:** Use ONE consistent filename per dataset: `data/DATASET_NAME/COMPREHENSIVE_REPORT.pdf`. Do not create versioned or dated filenames.

**Report writing guidelines:**
- Avoid sensationalized framing — no "Mystery - Resolved" sections
- TnpB classification — report as transposases unless clearly within a CRISPR/Cas locus
- Use neutral, factual language

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
b.search("unannotated AND giant")
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

### Sub-Agent Protocol
- Sub-agents are **leaf agents** — they should NOT spawn further agents
- Provide full context in the prompt (don't assume they know prior conversation)
- Each sub-agent produces a discrete output (markdown report, JSONL findings)
- Parent agent synthesizes outputs from all sub-agents

### Background Tasks
```python
result = Task(subagent_type="general-purpose", prompt="...", run_in_background=True)
# Returns immediately with task_id; use TaskOutput(task_id) to check
```

### Practical Tips
- **Parallel genome browser agents** work well (quarters or groups)
- **JSONL for findings** — easy to append, merge, and process
- **DuckDB write locks** — run subagents SEQUENTIALLY if they write to the database
- **Check database schema** before writing queries (`DESCRIBE table_name`)
- **Research external data** via WebFetch — don't guess PDB functions
- **Don't create new visualization code** when operators exist
- **Don't create ad-hoc report generators** — adapt existing ones

---

## Large Dataset Performance (>50k proteins)

**Rules:**
1. **Never** `b.search("")` on large datasets — use `b.total_proteins()` or specific predicates
2. **Always** check `len(result)` before iterating
3. **Always** use `[:10]` or `[:100]` limits when iterating
4. **Combine** specific predicates (AND/OR) to narrow results
5. **Aggregate** in SQL, not Python loops
6. **Limit** to 20-30 visualizations per analysis
7. **If query takes >5 seconds**, stop and refine
8. **NEVER use correlated subqueries** — rewrite as JOINs (see below)

```python
# BAD
all_proteins = b.search("")
for pid in all_proteins: ...

# GOOD
targets = b.search("hydrogenase AND membrane_protein")
for pid in targets[:10]: ...

# GOOD — SQL aggregation
stats = b.store.execute("""
    SELECT COUNT(*) as total, AVG(length) as avg_size, MAX(length) as max_size
    FROM proteins
""").fetchone()
```

### DuckDB Query Patterns for Large Datasets

**CRITICAL: Correlated subqueries destroy performance on >1M row tables.**

DuckDB cannot optimize nested `SELECT` inside `WHERE EXISTS` or `WHERE x = (SELECT ...)` when the outer query is large. These cause per-row subquery execution, eating memory and swap.

```sql
-- BAD: Correlated subquery — O(n * m), causes swap thrashing on 2.9M proteins
SELECT su.protein_id,
    EXISTS (
        SELECT 1 FROM proteins p2
        WHERE p2.contig_id = (SELECT contig_id FROM proteins WHERE protein_id = su.protein_id)
          AND ABS(p2.gene_index - (SELECT gene_index FROM proteins WHERE protein_id = su.protein_id)) BETWEEN 1 AND 3
    )
FROM sample_unann su

-- GOOD: Materialize context first, then JOIN
WITH sample AS (
    SELECT pp.protein_id, p.contig_id, p.gene_index
    FROM protein_predicates pp
    JOIN proteins p ON pp.protein_id = p.protein_id
    WHERE list_contains(pp.predicates, 'unannotated')
    ORDER BY RANDOM() LIMIT 500
)
SELECT DISTINCT s.protein_id
FROM sample s
JOIN proteins p2 ON s.contig_id = p2.contig_id
    AND ABS(p2.gene_index - s.gene_index) BETWEEN 1 AND 3
JOIN annotations a ON p2.protein_id = a.protein_id AND a.source = 'pfam'
```

**General pattern for neighborhood queries:**
1. **Materialize the seed set** in a CTE with `contig_id` and `gene_index`
2. **JOIN** to find neighbors (same contig, gene_index ± window)
3. **JOIN** to annotations/predicates for neighbor features
4. **Never** put a subquery inside WHERE that references the outer row

**For per-genome cross-tabs:**
```sql
-- BAD: Python loop over genomes
for genome in genomes:
    count = db.execute(f"SELECT COUNT(*) FROM ... WHERE bin_id = '{genome}'")

-- GOOD: Single GROUP BY query
SELECT p.bin_id,
    COUNT(DISTINCT CASE WHEN a.accession = 'K02274' THEN p.protein_id END) as coxI,
    COUNT(DISTINCT CASE WHEN a.accession = 'K02275' THEN p.protein_id END) as coxII,
    COUNT(DISTINCT CASE WHEN list_contains(pp.predicates, 'nife_group4') THEN p.protein_id END) as g4_hyd
FROM proteins p
LEFT JOIN annotations a ON p.protein_id = a.protein_id
LEFT JOIN protein_predicates pp ON p.protein_id = pp.protein_id
GROUP BY p.bin_id
```

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

## Active Datasets

### Omnitrophota (2026-02-08, In Progress)

**Database:** `data/omni_production/bennu.duckdb` (5.3 GB)
**Scale:** 1,831 MAGs (~1,637 unique), 355,904 contigs, 2,921,111 proteins, 11,249,502 annotations

| Source | Hits |
|--------|------|
| PFAM | 4,880,659 |
| KEGG | 3,940,782 |
| DefenseFinder | 2,424,782 |
| HydDB | 3,279 |

**Missing data:** No taxonomy (all "unknown"), no CheckM2, no GC content, no embeddings, no VOGdb/CAZy.

**Emerging model:** Omnitrophota are **sessile, surface-specialist syntrophs** that invest heavily in polysaccharide decoration and use fermentative metabolism to produce acetate and H2 for partners.

**Key findings from Phase 1 (Survey):**
- Mega-proteins up to 85,804 aa — entire biosynthesis pathways in single ORFs
- Group 4 NiFe hydrogenases in 74.6% of genomes (H2-evolving, obligate syntrophy)
- 2.3% of proteome = glycosyltransferases (2-3x normal bacteria)
- Extreme c-di-GMP signaling density (15+ GGDEF/PilZ per genome)
- Non-motile (<4% flagella) but universal type IV pili + adhesins
- Rnf + acetate kinase energy strategy — no canonical respiratory chain
- 46,646 Radical SAM proteins — dual defense (Viperin) + biosynthetic roles
- DUF1015 in 85% of genomes — top characterization target
- 194 duplicate genome pairs detected
- 13,422 giant unknown proteins (largest unannotated: 39,880 aa)

**Analysis plan:** `data/omni_production/ANALYSIS_PLAN.md`
**Survey outputs:** `data/omni_production/survey/`
**Exploration outputs:** `data/omni_production/exploration/`

**Status:** Phase 1 (Survey) complete. Phase 2 (Explore) in progress. Phases 3-4 pending.

### Hinthialibacterota (2026-02-05, Complete)

**Database:** `data/hinthialibacterota_production/bennu.duckdb`
41 genomes, 184,136 proteins. MANUSCRIPT.md with 30 literature citations. PDF rendered.

### GJALLARVIRUS (2026-02-06, Complete)

**Database:** `data/heimdall_megavirus_production/bennu.duckdb`
473 kb genome, 588 proteins. Comprehensive report complete.

## Archives

| Dataset | Date | Location | Size |
|---------|------|----------|------|
| Altiarchaeota (63 genomes) | 2026-02-03 | `data/archives/altiarchaeota_2026-02-03/` | 789 MB |
| Thorarchaeota/Heimdall | 2026-02-03 | `data/archives/thorarchaeota_2026-02-03/` | ~50 MB |
| Heimdall Megavirus | 2026-02-04 | `data/archives/heimdall_megavirus_2026-02-04/` | 16 MB |
| BioFrame DSL (retired) | 2026-02-07 | `data/archives/bioframe_2026-02-07/` | — |

## TODO

- [ ] **Build `/atlas` skill** — Hierarchical annotation census with Context-First protocol baked in
- [x] **Sharur rename** — completed 2026-02-09
- [ ] **Omnitrophota: Run CheckM2 + GTDB-Tk** for quality and taxonomy
- [ ] **Omnitrophota: VOGdb annotations** pending Aksha refactor
