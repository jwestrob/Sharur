# Survey Skill

Comprehensive, systematic documentation of a metagenomic dataset. Not just counting - actual science with synthesis and interpretation.

For curiosity-driven exploration, use `/explore` instead (run AFTER survey completes).

**AGENT ARCHITECTURE:**
- You CAN spawn subagents to write specific reports (metabolic_reconstruction.md, defense_systems.md, etc.)
- Your subagents are LEAF AGENTS - they CANNOT spawn their own subagents
- Tell each subagent: "You are a leaf agent. DO NOT use the Task tool. Write your report and return."

**CONCURRENCY: DuckDB does not support concurrent writes. Run subagents SEQUENTIALLY, not in parallel.**

> **Mandatory:** Follow the shared validation protocols in `_validation_protocols.md`.
> This file contains authoritative inline versions of those protocols with additional
> survey-specific context. Verify accession names, use COUNT(DISTINCT protein_id),
> and apply Context-First protocol for annotations averaging >10 hits/genome.

---

## Mission

**Document what's there, comprehensively.** Go deep in each functional area. Synthesize patterns. Identify what's core vs variable. Note questions for follow-up. Build understanding, not just catalogs.

---

## Writing Style

**Write prose paragraphs, not bullet lists.** Your reports should read like sections of a scientific paper — flowing narrative that interprets findings in context, not enumerated data dumps.

**BAD (bullet list):**
```
- 39/41 genomes have NiFe Group 4 hydrogenases
- 2 genomes lack hydrogenases entirely
- Group 4 Mbh complexes are dominant
- Maturation genes (HypA-F) co-located in 35/39 genomes
- FeFe hydrogenases found in 3 genomes
```

**GOOD (prose):**
```
Hydrogen metabolism is nearly universal in Hinthialibacterota, with 39 of 41 genomes
encoding NiFe hydrogenases. Group 4 Mbh-type energy-conserving complexes dominate,
consistent with a syntrophic lifestyle where H2 production drives interspecies electron
transfer. Maturation genes (HypA-F) co-localize with catalytic subunits in 35 of these
genomes, suggesting functional operons rather than orphan genes. The two genomes lacking
hydrogenases (GCA_029558215 and GCA_029558220) are both highly fragmented (>500 contigs),
so absence likely reflects assembly gaps rather than genuine metabolic loss. Three genomes
additionally encode FeFe hydrogenases, which may serve as fermentative H2-evolving enzymes
under conditions where the Mbh complex is insufficient.
```

**When bullets ARE appropriate:** Genuinely list-like content such as enzyme inventories, feature counts, or method descriptions. But even then, **open with a prose paragraph** that interprets the list before presenting it.

**Example — list with interpretive lead-in:**
```
The defense repertoire is dominated by innate immunity systems, with CRISPR-Cas
playing a secondary role. The following systems were identified across the dataset:

- CRISPR-Cas Type I-C: 24/41 genomes (59%)
- BREX Type 1: 18/41 genomes (44%)
- Paris defense system: 15/41 genomes (37%)
- DISARM Type 2: 8/41 genomes (20%)
```

---

## Workflow with Subagents

**After initial reconnaissance**, spawn subagents to write detailed reports for major functional areas:

### Example Subagent Spawning Pattern

```python
# Example spawning pattern
Task(
    subagent_type="general-purpose",
    description="Metabolic reconstruction analysis",
    prompt="""You are a leaf agent analyzing energy metabolism and carbon flow for the survey.

**CRITICAL: You are a LEAF AGENT. DO NOT spawn sub-agents or use the Task tool.**

**WRITING STYLE: Write in flowing prose paragraphs. Minimize bullet points. Your report
should read like a section of a scientific manuscript — interpret findings, explain
biological significance, and connect observations into a coherent narrative. Open every
section with a prose paragraph before any lists or tables.**

**CONTEXT-FIRST ANALYSIS: When a PFAM domain or KO averages >10 hits per genome, it is a
superfamily-level annotation. Do NOT use the name as a functional claim. Run co-annotation
checks, verify neighborhood context with b.get_neighborhood(pid, window=5, all_annotations=True),
and only claim specific function if supported by pathway-specific markers. Check accession
names in the database before reporting.**

Dataset: data/hinthialibacterota_v3/sharur.duckdb

Your task: Write a comprehensive metabolic_reconstruction.md report covering:
- Energy metabolism (respiration, fermentation, electron carriers)
- Carbon catabolism (what they eat)
- Biosynthesis (what they can make)
- Metabolic strategy narrative

Output: Write to data/hinthialibacterota_v3/survey/metabolic_reconstruction.md

Use the Sharur operators to query. Synthesize, don't just enumerate. See /Users/jacob/Documents/Sandbox/Bennu2/bennu/.claude/skills/survey.md for detailed guidance on each area."""
)

# Wait for it to complete, then spawn next subagent for cell_surface_biology.md, etc.
```

### Subagents Should Write To:
- `survey/metabolic_reconstruction.md`
- `survey/cell_surface_biology.md`
- `survey/secondary_metabolism.md`
- `survey/defense_immunity.md`
- `survey/novel_proteins.md`
- `survey/scratchpad.md` (shared - append observations)

### You (Survey Agent) Should Write:
- `survey/summary.md` - Executive overview synthesizing all subagent reports
- `survey/genome_profiles.tsv` - Feature matrix (do this first, before spawning subagents)
- `survey/findings.jsonl` - High-level structured findings

**Run subagents SEQUENTIALLY** (DuckDB limitation). Each subagent reads the database and writes to its own report file.

---

## Census Phase (MANDATORY FIRST STEP)

**Before spawning any subagents, run the census.** This gives you and every subagent a factual inventory of what the database contains, flags superfamilies, and prevents the single most common analysis error: treating fold-level annotations as functional claims.

### Step 1: Annotation Source Census

```python
# What annotation sources are loaded and how much coverage do they provide?
source_census = b.store.execute("""
    SELECT source,
           COUNT(DISTINCT protein_id) as n_proteins,
           COUNT(*) as n_annotations,
           COUNT(DISTINCT accession) as n_accessions
    FROM annotations
    GROUP BY source
    ORDER BY n_proteins DESC
""").fetchall()

n_genomes = b.store.execute("SELECT COUNT(DISTINCT bin_id) FROM proteins").fetchone()[0]
n_proteins = b.store.execute("SELECT COUNT(*) FROM proteins").fetchone()[0]

print(f"Dataset: {n_genomes} genomes, {n_proteins} proteins")
print(f"\nAnnotation Sources:")
for source, n_prot, n_annot, n_acc in source_census:
    print(f"  {source}: {n_prot} proteins, {n_annot} annotations, {n_acc} accessions")
```

### Step 2: Top 20 Annotations per Source + Superfamily Flagging

```python
# For each source, get top 20 annotations and flag superfamilies
for source, _, _, _ in source_census:
    top_annots = b.store.execute("""
        SELECT accession, name,
               COUNT(DISTINCT protein_id) as n_proteins
        FROM annotations
        WHERE source = ?
        GROUP BY accession, name
        ORDER BY n_proteins DESC
        LIMIT 20
    """, [source]).fetchall()

    print(f"\n--- {source.upper()} Top 20 ---")
    for accession, name, n_prot in top_annots:
        hits_per_genome = n_prot / n_genomes
        flag = " ** SUPERFAMILY **" if hits_per_genome > 10 else ""
        print(f"  {accession} ({name}): {n_prot} proteins ({hits_per_genome:.1f}/genome){flag}")
```

### Step 3: Predicate Frequency Summary

```python
pred_counts = b.store.execute("""
    SELECT pred, COUNT(*) as n
    FROM protein_predicates,
         LATERAL unnest(predicates) AS t(pred)
    GROUP BY pred
    ORDER BY n DESC
""").fetchall()

print(f"\nPredicate Summary ({len(pred_counts)} distinct predicates):")
for pred, n in pred_counts[:30]:
    print(f"  {pred}: {n}")
```

### Step 4: Save Census Outputs

```python
import json
from pathlib import Path

survey_dir = Path("data/YOUR_DATASET/survey")
survey_dir.mkdir(exist_ok=True)

# Machine-readable census
census_data = {
    "n_genomes": n_genomes,
    "n_proteins": n_proteins,
    "sources": {s: {"n_proteins": np, "n_annotations": na, "n_accessions": nacc}
                for s, np, na, nacc in source_census},
    "superfamilies": [
        {"source": source, "accession": acc, "name": name, "n_proteins": n,
         "hits_per_genome": round(n / n_genomes, 1)}
        for source, _, _, _ in source_census
        for acc, name, n in b.store.execute("""
            SELECT accession, name, COUNT(DISTINCT protein_id)
            FROM annotations WHERE source = ?
            GROUP BY accession, name
            HAVING COUNT(DISTINCT protein_id) / ? > 10
        """, [source, n_genomes]).fetchall()
    ],
    "predicates": {pred: n for pred, n in pred_counts}
}

with open(survey_dir / "census.json", "w") as f:
    json.dump(census_data, f, indent=2)

# Human-readable census
census_lines = [
    f"# Dataset Census",
    f"\n**Genomes:** {n_genomes}",
    f"**Proteins:** {n_proteins}",
    f"\n## Annotation Sources\n",
]
for source, np, na, nacc in source_census:
    census_lines.append(f"- **{source}**: {np} proteins, {na} annotations, {nacc} accessions")

if census_data["superfamilies"]:
    census_lines.append(f"\n## Superfamily-Level Annotations (>10 hits/genome)\n")
    census_lines.append("| Source | Accession | Name | Proteins | Hits/Genome |")
    census_lines.append("|--------|-----------|------|----------|-------------|")
    for sf in census_data["superfamilies"]:
        census_lines.append(
            f"| {sf['source']} | {sf['accession']} | {sf['name']} | "
            f"{sf['n_proteins']} | {sf['hits_per_genome']} |"
        )
    census_lines.append("\n> These annotations reflect structural folds, not specific functions.")
    census_lines.append("> Do NOT use their names as functional claims without neighborhood validation.")

with open(survey_dir / "census.md", "w") as f:
    f.write("\n".join(census_lines))

print(f"Census saved to {survey_dir / 'census.json'} and {survey_dir / 'census.md'}")
```

### Step 5: Distribute Census to Subagents

When spawning topic subagents (metabolic, defense, etc.), include the census summary in their prompt context:

```python
# Read census for subagent context
census_text = (survey_dir / "census.md").read_text()

Task(
    subagent_type="general-purpose",
    description="Metabolic reconstruction analysis",
    prompt=f"""You are a leaf agent analyzing energy metabolism for the survey.

**CENSUS CONTEXT (from mandatory census phase):**
{census_text}

**Superfamilies flagged above are fold-level annotations — do NOT claim specific
enzymatic function from them without neighborhood validation.**

... [rest of subagent prompt]
"""
)
```

---

## Recommended Subagent Strategy

1. **Run census phase** (above — mandatory, do this first)
2. **Create genome_profiles.tsv** (feature matrix for all genomes)
3. **Spawn topic-specific subagents sequentially** with census context included in each prompt

---

## Coverage Areas

Your survey should address all major functional domains. Not every analysis applies to every dataset - use judgment about what's relevant and what's not. State what you're doing and why. State what you're skipping and why.

### Energy Metabolism & Respiration

**What do these organisms use for energy?** Map the electron flow from substrate to terminal acceptor.

- **Respiration components**: Cytochromes, terminal oxidases, quinones
- **Hydrogenases**: If detected, subtype them using HydDB classifications (Group 1-4, Mbh/Ech)
- **Respiratory complexes**: Complex I/II/III/IV, alternative complexes (Rnf, Nqr)
- **Fermentation**: PFOR, lactate dehydrogenase, butyrate/propionate pathways
- **Energy conservation**: ATP synthase, Rnf complex, ion gradients
- **Electron carriers**: Ferredoxins, flavoproteins, iron-sulfur proteins

**Synthesis**: Don't just list components - trace the electron flow. "Electrons from NADH → Complex I → Quinone → Cytochrome bc1 → Cytochrome c → Complex IV → O2" tells a story that counts alone don't.

### Carbon Metabolism & Biosynthesis

**What can they eat? What can they make?**

- **Carbon catabolism**:
  - Carbohydrate degradation (CAZymes - what substrates?)
  - Aromatic compound degradation (benzoyl-CoA pathway, gentisate, etc.)
  - Fatty acid beta-oxidation
  - Amino acid catabolism
- **Carbon fixation**:
  - Calvin cycle (RuBisCO + PRK for true Calvin, not just RuBisCO alone)
  - Wood-Ljungdahl pathway (check completeness: ACS, CODH, all 8 enzymes)
  - 3-hydroxypropionate, reverse TCA, etc.
- **Central metabolism**: Glycolysis, TCA cycle, pentose phosphate pathway
- **Biosynthesis**:
  - Amino acid synthesis (which ones can they make?)
  - Nucleotide synthesis (purine, pyrimidine pathways)
  - Cofactor synthesis (folate, thiamine, B12, etc.)

**Metabolic reconstruction**: For each genome (or genome group), create a narrative of carbon/energy flow. What comes in, how is it processed, what comes out?

### Secondary Metabolism & Biosynthetic Gene Clusters

**What specialized metabolites do they produce?**

- **BGC analysis**: GECCO results are in `annotations/gecco.tsv` (if available) - interpret the predictions
- **NRPS/PKS systems**: Giant multi-domain enzymes, domain architecture
- **Terpene synthases**: Isoprenoid biosynthesis
- **Ribosomally synthesized peptides**: Lanthipeptides, sactipeptides, etc.
- **Antibiotic resistance**: Efflux pumps, modification enzymes, target mutations

**Product prediction**: What might these BGCs produce? Antimicrobials (competition)? Siderophores (iron acquisition)? Signaling molecules (quorum sensing)?

### Cell Surface Biology & Interactions

**How do they interact with their environment and other cells?**

- **Surface structures**:
  - Pili (Type IV, Tad, etc.) - motility or attachment?
  - Flagella - motile or non-motile?
  - S-layers - structural glycoproteins
  - Fimbriae - adhesive organelles
- **Adhesins**: Giant surface proteins, Ig-like domains, lectins
- **Secretion systems**: Type I-VI (what do they secrete? toxins? effectors?)
- **Membrane architecture**: LPS biosynthesis, peptidoglycan modifications
- **Protein export**: Sec, Tat, twin-arginine translocation

**Ecological context**: Surface biology reflects lifestyle. Motile flagellated cells vs sessile piliated cells tell different stories.

### Defense & Immunity

**How do they defend against phages, mobile elements, and competitors?**

- **CRISPR-Cas systems**: Types (I, II, III, etc.), array counts, spacer diversity (if available)
- **Innate immunity**:
  - Restriction-modification systems
  - CBASS, BREX, Pycsar, Druantia
  - Defensive islands (Paris, Jallet, Wadjet, Gabija, etc.)
- **Toxin-antitoxin systems**: Abundance, types
- **Prophage elements**: VOG annotations, integrases, lysins
- **Abortive infection**: Systems that kill the infected cell to save the population

**Investment level**: Heavy defense investment suggests constant phage pressure or competitive environments. Minimal defense might indicate stable, low-threat niches.

### Novel & Unannotated Proteins

**What don't we understand yet?**

- **Giant proteins (>1000 aa)**: Systematic annotation attempt
  - Run relaxed PFAM searches (`--domE 1e-5`) for annotation recovery
  - Domain architecture analysis
  - Potential functions based on genomic context
- **High-abundance unknowns**: If 20% of proteins are "hypothetical", investigate the most abundant ones
- **Unusual domain architectures**: Fusions that shouldn't exist, exotic combinations
- **Orphan genes**: No homologs anywhere - genuinely novel or poor annotation?

**Scratchpad candidates**: Note proteins that would benefit from structure prediction (AlphaFold3) or structural homology search (Foldseek) for later follow-up.

### Specialized Analyses (When Relevant)

**Hydrogenase Subtyping** (if hydrogenases detected):
- Use HydDB classifications to assign Groups 1-4, Mbh/Ech
- Distinguish uptake vs bidirectional vs energy-conserving
- Check for maturation genes (HypA-F, HydE-G)

**Ecotype Clustering** (if multiple genomes with variation):
- Build genome × feature matrix
- Hierarchical clustering to identify metabolic subgroups
- Statistical enrichment testing per cluster
- Generate dendrogram + heatmap

**Pathway Completeness Matrix** (if metabolic genes detected):
- For key pathways (WL, Calvin, methanogenesis, etc.), check presence/absence of all components
- Report as fractions: "Wood-Ljungdahl 5/8 enzymes" vs "8/8 complete"
- Save as `pathway_completeness.json`

**CRISPR Array Analysis** (if CRISPR widespread):
- Spacer diversity, array structure
- Integration site conservation
- Evidence of recent infection events

---

## Interpretation Framework

For every major claim, provide:

### Pattern
What does the data show? (quantitative, with numbers)

### Mechanism
Why does this pattern exist biologically? What selective pressures or biochemical constraints explain it?

### Counterevidence
What would make you doubt this interpretation? What alternative explanations exist?

### Uncertainty
What don't you know? What would you need to confirm this?

### Conservation Context
If multi-genome: Is this universal (100%)? Core (>80%)? Variable (20-80%)? Rare (<20%)? This matters for interpretation.

---

## Accession Verification (MANDATORY)

**NEVER assume a PFAM/KEGG accession encodes a specific function from memory.** Always verify the `name` field in the database before reporting.

This rule exists because a false finding (PF04055 claimed as "benzoyl-CoA reductase", actually "Radical_SAM") propagated through survey → exploration → manuscript before being caught. One SQL query would have prevented it.

**Before reporting ANY PFAM/KEGG-based claim:**

```python
# Verify accession matches expected function
name = b.store.execute(
    "SELECT DISTINCT name FROM annotations WHERE accession = ?", [accession]
).fetchone()
print(f"{accession} = {name[0]}")  # Confirm this is what you think it is
```

**Include verification in your reports:**
```markdown
NiFe hydrogenases (PF00374, verified as "NiFeSe_Hases" in database): 108 proteins in 40/41 genomes.
```

**Red flags requiring extra verification:**
- Single accession with >500 proteins (likely a superfamily, not a specific enzyme)
- "Universal" presence (>90% of genomes) of a functional claim
- Accession you haven't personally verified in this session

## Counting Proteins vs Annotation Rows

**ALWAYS use `COUNT(DISTINCT protein_id)` when reporting protein counts.** Repeat domains (WD40, TPR, FG-GAP, Ig-like) produce multiple annotation rows per protein. `COUNT(*)` gives annotation hits, not protein counts — these can differ by 5-10x.

```python
# ✗ WRONG — counts annotation rows
b.store.execute("SELECT COUNT(*) FROM annotations WHERE accession = 'PF00400'")

# ✓ RIGHT — counts unique proteins
b.store.execute("SELECT COUNT(DISTINCT protein_id) FROM annotations WHERE accession = 'PF00400'")
```

## Hydrogenase Curation Protocol

The pipeline classifies all HydDB hits with DIAMOND subgroup assignments and tags them with subgroup predicates (`nife_group1` through `nife_group4`, `fefe_groupA` through `fefe_groupC`). Hits that lack PFAM corroboration are tagged `hyddb_needs_curation` — these include legitimate Group 4 NiFe hydrogenases (Hyf/Hyc/Mbh/Ech) AND Complex I false positives that share HMM similarity.

**When the survey finds hydrogenases, curate the flagged hits by neighborhood:**

```python
# Step 1: Get all hydrogenases and the subset needing curation
all_hyddb = b.search_by_predicates(has=["nife_group4"])  # or broader: ["nife_hydrogenase"]
needs_curation = b.search_by_predicates(has=["hyddb_needs_curation"])

# Step 2: For each flagged hit, check ±8 genes
for pid in needs_curation[:30]:  # cap to avoid runaway loops
    nbr = b.get_neighborhood(pid, window=8, all_annotations=True)
    # Inspect neighbor KOs and PFAMs, log verdict
```

**Hydrogenase evidence** (any → likely real):
- KEGG: K12136-K12145 (hyfA-J), K15828-K15833 (hycB-G)
- KEGG: K04651-K04656, K03605 (maturation: HypA-F, HycI)
- PFAM: PF00374 (NiFeSe_Hases) on a neighboring gene
- HydDB annotation on a neighboring gene

**Complex I evidence** (→ false positive):
- KEGG: K00330-K00343 (nuoA-N)
- PFAM: PF00346, PF00329 (Complex1_49kDa, Complex1_30kDa) on neighbors

**Report with curation context:** "N NiFe hydrogenases classified (X Group 1-3 corroborated by PF00374, Y Group 4 confirmed by neighborhood markers, Z rejected as Complex I false positives)."

---

## Context-First Analysis Protocol

**The domain tells you the fold. The neighbors tell you the function.**

Annotation databases assign labels based on structural fold similarity, not biological function. When a PFAM HMM or KOfam profile matches an entire enzyme superfamily, hit counts are high but the name is misleading. Always validate functional claims with genomic context and co-annotations.

### Superfamily Awareness Rule

When a domain or KO averages more than ~10 hits per genome, it is likely a **superfamily-level annotation** — the HMM/profile matches a structural fold shared across many unrelated enzyme families. Treat such hits as structural information, not functional evidence.

| Hits per genome | Interpretation | Action |
|----------------|----------------|--------|
| 1-3 | Specific enzyme or system | Name is likely informative |
| 3-10 | Multi-copy gene family | Name is informative but check paralogs |
| 10-50 | Broad enzyme class | Name reflects fold, not specific function |
| >50 | Superfamily | Name is purely structural; no pathway info |

**When the superfamily flag triggers:**
1. Do NOT use the domain/KO name as a functional claim
2. Run co-annotation analysis (see below)
3. Check genomic neighborhood of representative proteins
4. Only claim specific function if supported by context or pathway-specific markers

### Co-Annotation Validation

Before claiming a protein functions as enzyme Y, check what *other* annotations it carries:

```python
# Standard co-annotation check for any suspicious KO/domain
co_annots = b.store.execute("""
    SELECT a2.source, a2.accession, a2.name, COUNT(DISTINCT a1.protein_id) as n
    FROM annotations a1
    JOIN annotations a2 ON a1.protein_id = a2.protein_id
    WHERE a1.accession = 'K23108'   -- your claimed KO
      AND a2.accession != 'K23108'  -- exclude self
    GROUP BY a2.source, a2.accession, a2.name
    ORDER BY n DESC LIMIT 10
""").fetchall()
# If top co-annotations are from a DIFFERENT enzyme family → superfamily cross-hit
```

### Neighborhood Validation

For protein-level functional claims that imply a pathway or system, use `b.get_neighborhood()` with `all_annotations=True` to see the full annotation context of surrounding genes:

```python
# View all annotations for genes ±5 around a protein of interest
result = b.get_neighborhood(protein_id, window=5, all_annotations=True)
print(result)  # Shows every annotation source per gene (PFAM, KEGG, HydDB, DefenseFinder, VOGdb, CAZy)
```

**What to look for:**
- **Neighbors consistent with claimed function** → Claim supported
- **Neighbors suggest a different system** → Revise interpretation
- **No informative neighbors** → State uncertainty
- **Protein at contig edge** → Note potential fragmentation

### Claim Escalation Ladder

| Claim Level | Evidence Required |
|------------|-------------------|
| "Contains domain X" | Accession name verification only |
| "Functions as enzyme Y" | Name verification + co-annotation check |
| "Genome encodes pathway Z" | Multiple pathway-specific markers + co-localization |
| "Phylum universally performs W" | All above + conservation across genomes |

### KEGG Pathway Lookup

When you find a KEGG KO and want to understand its pathway context, use the KEGG REST API:

```python
import subprocess

# What pathways contain this KO?
result = subprocess.run(
    ["curl", "-s", "https://rest.kegg.jp/link/pathway/ko:K23108"],
    capture_output=True, text=True
)

# What KEGG modules (minimum functional units) use this KO?
result = subprocess.run(
    ["curl", "-s", "https://rest.kegg.jp/link/module/ko:K23108"],
    capture_output=True, text=True
)

# Get the full module definition
result = subprocess.run(
    ["curl", "-s", "https://rest.kegg.jp/get/md:M00551"],
    capture_output=True, text=True
)
```

**KEGG modules are more useful than pathways** for metabolic claims because they define minimum functional units.

**Caveats for divergent organisms:**
- KOfam profiles are built from characterized enzymes, mostly from model organisms
- Deeply branching lineages may have homologs too divergent for detection
- Report: "X/Y steps detected; undetected steps may reflect divergence rather than absence"

**Local reference:** `data/reference/ko_list` contains all 27,324 KO definitions. `scripts/pathway_completeness.py` defines 12 major pathways with component KOs.

---

## Statistical Validation

For comparative claims, consider statistical validation when appropriate:

- **Comparing genome groups**: Fisher exact test (categorical) or t-test (continuous)
- **Correlation between features**: Spearman correlation
- **Effect size**: Odds ratio, Cohen's d - is the difference meaningful?
- **Multiple testing correction**: If testing many hypotheses, adjust p-values (Bonferroni, FDR)

But don't force statistics where they don't make sense:
- 3 genomes? Descriptive only, no stats
- Single genome? Focus on completeness and capabilities
- Small differences? Report them honestly without claiming significance

Save statistical results in `statistical_results.json` for reproducibility.

---

## Output Structure

Generate **separate reports for major topics** (easier to navigate than one giant document):

### Required Reports

1. **`metabolic_reconstruction.md`** - Energy, carbon, biosynthesis
   - Electron transport chains
   - Carbon catabolism pathways
   - Biosynthetic capabilities
   - Overall metabolic strategy

2. **`genome_profiles.tsv`** - Feature matrix
   - One row per genome
   - Columns: genome_id, n_proteins, n_contigs, annotation_rate, then feature counts
   - Used for downstream statistical analysis

3. **`summary.md`** - Executive overview
   - Dataset characteristics
   - Major findings (top 5-10)
   - Key patterns across genomes
   - Questions for follow-up

### Optional Reports (Generate if Relevant)

4. **`cell_surface_biology.md`** - Pili, adhesins, secretion systems
5. **`secondary_metabolism.md`** - BGCs, NRPS/PKS, antibiotics
6. **`defense_immunity.md`** - CRISPR, RM, innate immunity
7. **`novel_proteins.md`** - Giants, unknowns, unusual domains
8. **`ecotype_analysis.md`** - Clustering results (if multiple genomes)
9. **`pathway_completeness.json`** - Structured data for metabolic pathways

### Working Documents

10. **`scratchpad.md`** - Your working notes
    - Observations that don't fit structured findings
    - Questions that arise during analysis
    - Candidates for follow-up (structure prediction, Foldseek, etc.)
    - Cross-references between analyses
    - "I noticed X correlates with Y, should investigate"

11. **`findings.jsonl`** - Structured findings (append-only)
    - Major discoveries in machine-readable format
    - Each line is a JSON object with: title, category, description, evidence, significance

### Hypothesis Tracking & Provenance

Log analytical steps and link findings to hypotheses for cross-session persistence:

```python
# Log key survey steps with provenance chaining
e1 = b.log_provenance("Census: annotation sources", f"{len(source_census)} sources loaded")
e2 = b.log_provenance("Hydrogenase inventory", f"{n_hydro} NiFe in {n_genomes_hydro} genomes", parent_ids=[e1.entry_id])
e3 = b.log_provenance("ETC completeness check", "Complex III/IV absent in 38/41", parent_ids=[e1.entry_id])

# Propose hypotheses discovered during the survey
h = b.propose_hypothesis("Lineage is obligately syntrophic based on incomplete ETC")

# Link evidence from survey analyses
b.add_evidence(h.hypothesis_id, "Survey: NiFe distribution", f"{n} genomes positive", True, 0.7)
b.add_evidence(h.hypothesis_id, "Survey: ETC gaps", "Missing Complex III and IV in 38/41", True, 0.9)

# Review all hypotheses
print(b.hypothesis_summary())
print(b.list_hypotheses())

# Render provenance DAG for the paper
b.render_provenance(title="Survey Analysis", output_path="figures/survey_provenance.mermaid")
```

Hypotheses persist across sessions and appear in `b.resume()`.

---

## Flexibility & Judgment

**Not every analysis applies to every dataset.** Use judgment:

- **1 genome?** Focus on metabolic completeness, capabilities, and potential niche
- **50 genomes?** Find patterns, outliers, conserved vs variable features
- **Hydrogenases present?** Subtype them (HydDB provides classifications)
- **GECCO results available?** Interpret the BGC predictions (don't re-run)
- **Giant unknowns?** Try relaxed PFAM (`--domE 1e-5`) for annotation recovery
- **Low annotation rate (<50%)?** Note this limitation in interpretation

State what you're doing and why. State what you're skipping and why.

Examples:
- "Skipping ecotype clustering because only 3 genomes - insufficient for meaningful clustering"
- "Not subtyping hydrogenases because only 2 genes detected across all genomes"
- "Skipping CRISPR array analysis because <10% of genomes have CRISPR-Cas"
- "Running pathway completeness matrix because strong metabolic annotation coverage (78%)"

---

## Thoroughness vs Mechanical Execution

**Go deep, don't just check boxes.**

Bad survey:
- "Defense systems: 500 CRISPR genes, 200 RM genes, 100 TA genes"
- "Hydrogenases: 40 detected"
- "CAZymes: 1000 detected"

Good survey:
- "Defense investment is exceptionally high (15% of proteins) with universal CRISPR-Cas Type I-C systems in all 41 genomes plus extensive innate immunity (Paris, Jallet, Wadjet). This suggests constant phage pressure in their natural environment, consistent with their syntrophic lifestyle requiring stable partnerships despite mobile genetic element threats."

- "Hydrogenases are nearly universal (40/41 genomes, 98%) with Group 4 Mbh-type energy-conserving complexes dominant. The single exception (GCA_029558215_1) is highly fragmented (725 contigs), suggesting hydrogenase genes were likely missed rather than genuinely absent. This universality indicates H2 cycling is core to Hinthialibacterota metabolism."

**Synthesis over enumeration.** Every number should answer "so what biologically?"

---

## MAG Quality Considerations

Always report genome fragmentation when making claims:

```python
# Check fragmentation before making absence claims
contig_stats = b.store.execute("""
    SELECT bin_id, COUNT(DISTINCT contig_id) as n_contigs,
           ROUND(COUNT(*)::FLOAT / COUNT(DISTINCT contig_id), 1) as genes_per_contig
    FROM proteins GROUP BY bin_id ORDER BY n_contigs DESC
""")
```

| Contigs | Interpretation |
|---------|----------------|
| <50 | High quality - absence claims reasonable |
| 50-200 | Moderate - include caveats |
| >200 | Fragmented - use "not detected" not "absent" |
| genes/contig <5 | Very fragmented - many genes likely missing |

**Language:**
- ✅ "Not detected in this assembly"
- ✅ "Not recovered in MAG"
- ❌ "Lacks" or "Absent" (unless high-quality genome)

Before saying "Genome A has X but Genome B doesn't", verify B isn't just more fragmented.

---

## Cytochrome Validation

**Always validate respiratory system claims, especially cytochrome absence.**

Cytochromes are ubiquitous in bacterial respiration. Claims of absence require validation:

1. **Predicate search**: `b.search_by_predicates(has=["cytochrome"])`
2. **Sequence motif search**: Look for CxxCH heme-binding sites in protein sequences
3. **Raw annotation search**: Query for "cytochrome" in annotation names

Use all three methods before claiming cytochrome absence. If all three find nothing AND genome is high-quality (<50 contigs), absence is plausible. Otherwise, state uncertainty.

---

## Tips for Success

1. **Start with reconnaissance** - understand the dataset before diving into analyses
2. **Let questions emerge** - don't force a predetermined structure
3. **Synthesize, don't enumerate** - every count needs biological interpretation
4. **Note limitations** - fragmentation, annotation coverage, missing data
5. **Use scratchpad liberally** - capture observations and questions as you work
6. **Cross-reference analyses** - "The high CAZyme load (227/1000 proteins) is consistent with the complete aromatic degradation pathway, suggesting these organisms specialize in complex organic matter breakdown"
7. **Generate figures sparingly** - survey doesn't need neighborhood diagrams (that's /explore), but feature distribution plots, heatmaps, dendrograms are useful

---

## Example Workflow

```python
from sharur.operators import Sharur
from pathlib import Path

b = Sharur("data/YOUR_DATASET/sharur.duckdb")
output_dir = Path("data/YOUR_DATASET/survey")
output_dir.mkdir(exist_ok=True)

# 1. Reconnaissance
overview = b.overview()
n_genomes = overview['n_genomes']
annotation_coverage = overview['annotation_rate']

# Open scratchpad
scratchpad = open(output_dir / "scratchpad.md", "w")
scratchpad.write(f"# Survey Scratchpad\n\n")
scratchpad.write(f"Dataset: {n_genomes} genomes, {annotation_coverage:.1%} annotated\n\n")

# 2. Build genome profiles matrix
genome_stats = b.store.execute("""
    SELECT bin_id,
           COUNT(*) as n_proteins,
           COUNT(DISTINCT contig_id) as n_contigs,
           AVG(sequence_length) as avg_length
    FROM proteins
    GROUP BY bin_id
""").df()

genome_stats.to_csv(output_dir / "genome_profiles.tsv", sep="\t", index=False)

# 3. For each functional area, query and synthesize
# Energy metabolism
hydrogenases = b.search_by_predicates(has=["nife_hydrogenase"])
cytochromes = b.search_by_predicates(has=["cytochrome"])
# ... analyze and write to metabolic_reconstruction.md

# 4. Note observations in scratchpad
scratchpad.write(f"## Observations\n")
scratchpad.write(f"- Hydrogenases in {len(set(h['bin_id'] for h in hydrogenases.data))}/{n_genomes} genomes\n")
scratchpad.write(f"- TODO: Check if the 2 genomes without hydrogenases are fragmented\n\n")

scratchpad.close()
```

---

## Remember

You are building **systematic documentation with synthesis and interpretation**, not a mechanical checklist. Be thorough. Be thoughtful. Note what's interesting, what's surprising, what needs follow-up. Use the scratchpad as your laboratory notebook.

**The goal is comprehensive understanding, not task completion.**
