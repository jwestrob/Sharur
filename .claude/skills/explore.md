# Explore Skill

Curiosity-driven exploration of Sharur metagenomic data to discover interesting loci, test hypotheses, and build mechanistic models.

**Run AFTER `/survey` completes** - use survey results to identify questions worth investigating.

**AGENT ARCHITECTURE:**
- You CAN spawn subagents to test specific hypotheses or investigate interesting loci
- Your subagents are LEAF AGENTS - they CANNOT spawn their own subagents
- Tell each subagent: "You are a leaf agent. DO NOT use the Task tool. Investigate this question and return findings."

**CONCURRENCY: DuckDB does not support concurrent writes. Run subagents SEQUENTIALLY, not in parallel.**

> **Mandatory:** Follow the shared validation protocols in `_validation_protocols.md`.
> This file contains authoritative inline versions of those protocols with additional
> exploration-specific context. Verify accession names, use COUNT(DISTINCT protein_id),
> and apply Context-First protocol for annotations averaging >10 hits/genome.

---

## Start with Survey Results

**Before exploring**, read what the survey found:

```python
from pathlib import Path

survey_dir = Path("data/YOUR_DATASET/survey")

# Read survey outputs
if survey_dir.exists():
    # Executive summary
    summary = (survey_dir / "summary.md").read_text()

    # Genome feature matrix
    genome_profiles = pd.read_csv(survey_dir / "genome_profiles.tsv", sep="\t")

    # Topic reports
    metabolic = (survey_dir / "metabolic_reconstruction.md").read_text() if (survey_dir / "metabolic_reconstruction.md").exists() else None

    # Scratchpad notes
    scratchpad = (survey_dir / "scratchpad.md").read_text() if (survey_dir / "scratchpad.md").exists() else None
```

**Use survey findings to identify:**
- Surprising patterns that need explanation
- Conflicting observations that need resolution
- Outlier genomes worth investigating
- Gaps in understanding flagged in scratchpad
- Hypotheses to test

---

## Data Specification

**All exploration data lives in a single directory relative to the database.**

```
{DB_DIR}/exploration/
├── findings.jsonl           # Append-only log of all discoveries
├── exploration_state.json   # Session state (what's been searched/browsed)
├── exploration_summary.md   # Human-readable findings summary
├── genome_profiles.tsv      # Per-genome feature counts (multi-genome only)
└── genome_comparison.md     # Cross-genome comparative analysis (multi-genome only)
```

### File Formats

#### findings.jsonl
Append-only JSONL. Each line is a JSON object with **mandatory provenance fields**:
```json
{
  "timestamp": "2026-02-02T10:30:00",
  "category": "defense",
  "title": "High CRISPR load in GCA_018260655",
  "description": "18 CRISPR-associated proteins...",
  "proteins": ["protein_id_1", "protein_id_2"],
  "evidence": {"n_genomes": 45, "conservation": "71%"},
  "location": "contig_id:start-end",
  "priority": "high",
  "provenance": {
    "query": "b.store.execute(\"SELECT COUNT(*) FROM annotations a JOIN proteins p ON a.protein_id = p.protein_id WHERE a.accession = 'PF01396' AND p.bin_id = 'GCA_018260655'\").fetchone()",
    "raw_result": [18],
    "accession_verified": "PF01396 = 'CRISPR_assoc' (name confirmed in database)",
    "interpretation": "18 CRISPR-associated proteins in genome GCA_018260655"
  }
}
```

#### exploration_state.json
Session state tracking:
```json
{
  "searches_completed": ["giant+unannotated", "crispr_search"],
  "contigs_browsed": ["contig_1", "contig_2"],
  "genomes_explored": ["GCA_001", "GCA_002"],
  "predicates_searched": ["adhesin+giant"],
  "last_updated": "2026-02-02T10:30:00",
  "total_findings": 125
}
```

#### genome_profiles.tsv
Tab-separated, one row per genome. Columns:
- `bin_id`, `n_proteins`, `avg_length`, `max_length`, `avg_gc`, `annotated`, `annotation_rate`
- Then one column per tracked predicate (e.g., `hydrogenase`, `crispr_associated`, etc.)

#### genome_comparison.md
Markdown report with:
- Dataset overview (genome count, protein count, predicate count)
- Feature distribution matrix (which predicates are in how many genomes)
- Core vs variable feature categories
- Outlier genomes (by size, defense load, novelty, etc.)
- Key biological insights

### Report Generation

The report generator (`scripts/generate_exploration_report.py`) reads these files:
- `findings.jsonl` → Finding sections by category
- `genome_profiles.tsv` → Per-genome statistics table
- `genome_comparison.md` → Cross-genome analysis section

Report output: `{DB_DIR}/altiarchaeota_exploration_report.pdf`

---

## Usage

```
/explore                     # General exploration
/explore --focus metabolism  # Focus on metabolic systems
/explore --focus phage       # Focus on prophage regions
/explore --focus defense     # Focus on defense systems (CRISPR, R-M, TA)
/explore --focus novel       # Focus on unannotated/hypothetical proteins
/explore --focus adhesion    # Focus on cell adhesion/surface proteins
/explore --contig <id>       # Browse specific contig window-by-window
/explore --resume            # Resume exploration from saved state
```

## Prompt

You are exploring a metagenomic database using Sharur. Your mission is to find scientifically interesting loci that warrant further investigation.

**CRITICAL: You must persist all findings to disk as you discover them.** Do not rely on context memory alone.

### Database Location
The Sharur database is at: `data/sharur.duckdb` (or specified path)

### Persistence Infrastructure

**ALWAYS initialize the exploration session first:**

```python
import sys
import json
import os
from datetime import datetime
from pathlib import Path

sys.path.insert(0, '.')
from sharur.operators import Sharur

# Initialize Sharur
DB_PATH = "data/altiarchaeota_production/sharur.duckdb"  # Adjust as needed
b = Sharur(DB_PATH)

# Set up persistence directory
EXPLORE_DIR = Path(DB_PATH).parent / "exploration"
EXPLORE_DIR.mkdir(exist_ok=True)

FINDINGS_FILE = EXPLORE_DIR / "findings.jsonl"
STATE_FILE = EXPLORE_DIR / "exploration_state.json"
SUMMARY_FILE = EXPLORE_DIR / "exploration_summary.md"

def log_finding(category: str, title: str, description: str,
                proteins: list = None, evidence: dict = None,
                location: str = None, priority: str = "medium",
                provenance: dict = None):
    """Append a finding to the persistent findings log.

    PROVENANCE IS MANDATORY for any finding that makes a functional claim.
    The provenance dict must contain:
      - query: The exact code/SQL that produced the evidence
      - raw_result: The literal output (count, list, accession name, etc.)
      - accession_verified: For PFAM/KEGG claims, the NAME field from the database
      - interpretation: The human-readable claim derived from the raw result

    Example:
        provenance={
            "query": "b.store.execute(\"SELECT name FROM annotations WHERE accession = 'PF04055' LIMIT 1\").fetchone()",
            "raw_result": ["Radical_SAM"],
            "accession_verified": "PF04055 = 'Radical_SAM' (NOT benzoyl-CoA reductase)",
            "interpretation": "1,790 Radical SAM superfamily proteins (generic, not pathway-specific)"
        }
    """
    finding = {
        "timestamp": datetime.now().isoformat(),
        "category": category,
        "title": title,
        "description": description,
        "proteins": proteins or [],
        "evidence": evidence or {},
        "location": location,
        "priority": priority,  # high, medium, low
        "provenance": provenance or {}
    }
    with open(FINDINGS_FILE, "a") as f:
        f.write(json.dumps(finding) + "\n")
    print(f"[LOGGED] {category}: {title}")
    return finding

def load_findings():
    """Load all findings from the log."""
    if not FINDINGS_FILE.exists():
        return []
    findings = []
    with open(FINDINGS_FILE) as f:
        for line in f:
            if line.strip():
                findings.append(json.loads(line))
    return findings

def save_state(state: dict):
    """Save exploration state (what's been searched, etc.)."""
    state["last_updated"] = datetime.now().isoformat()
    with open(STATE_FILE, "w") as f:
        json.dump(state, f, indent=2)

def load_state():
    """Load exploration state."""
    if not STATE_FILE.exists():
        return {
            "searches_completed": [],
            "contigs_browsed": [],
            "proteins_examined": [],
            "genomes_explored": [],
            "focus_areas_done": []
        }
    with open(STATE_FILE) as f:
        return json.load(f)

def update_summary():
    """Regenerate the summary markdown from all findings."""
    findings = load_findings()
    if not findings:
        return

    # Group by category
    by_category = {}
    for f in findings:
        cat = f.get("category", "uncategorized")
        if cat not in by_category:
            by_category[cat] = []
        by_category[cat].append(f)

    lines = [
        "# Exploration Findings Summary",
        f"\n**Total findings:** {len(findings)}",
        f"**Last updated:** {datetime.now().isoformat()}",
        ""
    ]

    for category, items in sorted(by_category.items()):
        lines.append(f"\n## {category.replace('_', ' ').title()} ({len(items)} findings)\n")
        for item in sorted(items, key=lambda x: x.get("priority", "medium") == "high", reverse=True):
            priority_marker = "**[HIGH]** " if item.get("priority") == "high" else ""
            lines.append(f"### {priority_marker}{item['title']}")
            lines.append(f"{item['description']}")
            if item.get("proteins"):
                lines.append(f"\n**Proteins:** {', '.join(item['proteins'][:5])}")
            if item.get("location"):
                lines.append(f"\n**Location:** {item['location']}")
            lines.append("")

    with open(SUMMARY_FILE, "w") as f:
        f.write("\n".join(lines))
    print(f"[SUMMARY] Updated {SUMMARY_FILE}")

# Initialize state
state = load_state()
```

### Writing Style

**Write prose paragraphs, not bullet lists.** Exploration findings, synthesis documents, and hypothesis reports should read like sections of a scientific paper.

**Finding descriptions** in `log_finding()` should be complete sentences or paragraphs, not telegraphic fragments:

```python
# BAD - telegraphic bullet-style description
log_finding(
    category="defense",
    title="CRISPR-Cas Type I-C in genome X",
    description="Type I-C system. 6 cas genes. 42 spacers. Adjacent to RM system.",
    ...
)

# GOOD - prose description with interpretation
log_finding(
    category="defense",
    title="CRISPR-Cas Type I-C with extensive spacer array in genome X",
    description="A complete Type I-C CRISPR-Cas system spans genes 145-152, encoding "
        "the full Cascade complex (Cas5, Cas7, Cas8c) plus Cas3 helicase-nuclease. "
        "The adjacent CRISPR array contains 42 spacers, among the largest in the dataset, "
        "suggesting sustained phage pressure. An RM system (Type II) flanks the locus, "
        "creating a defense island that occupies 3% of this contig.",
    ...
)
```

**Synthesis documents** (`exploration_summary.md`, hypothesis reports) should flow as narrative prose. Interpret findings in biological context rather than listing them.

### Figure Legends

**EVERY call to `b.visualize_neighborhood()` MUST include `title=` and `legend=` parameters.** These are recorded in the manifest and rendered in the PDF report.

A good legend contains: (1) what the locus shows, (2) key genes highlighted, (3) the biological observation.

```python
# BAD - no legend context
b.visualize_neighborhood(pid, window=12, output_path="figures/locus.png")

# GOOD - informative title and legend
b.visualize_neighborhood(
    center_protein,
    window=12,
    output_path=f"exploration/figures/{bin_id}_defense_locus.png",
    title="Type III CRISPR-Cas Locus with RAMP Proteins",
    legend="CRISPR-Cas locus in genome X showing Cas10-Cmr2 complex with 4 RAMP "
        "proteins (red). Adjacent to a 42-spacer CRISPR array, indicating active "
        "adaptive immunity against mobile genetic elements."
)
```

**Legend writing tips:**
- Name the system or feature shown (e.g., "Wood-Ljungdahl operon", "Type IV pilus cluster")
- Mention the key proteins by name or domain
- State the biological interpretation (what this arrangement tells us)
- Include conservation context if relevant ("present in 38/41 genomes")

---

### Exploration Protocol

**IMPORTANT: Follow this protocol to avoid re-reading data and ensure persistence.**

#### 1. Check State First
```python
state = load_state()
print(f"Previously searched: {state['searches_completed']}")
print(f"Contigs browsed: {len(state['contigs_browsed'])}")
```

#### 2. Log Findings Immediately
When you find something interesting, log it RIGHT AWAY before moving on:
```python
log_finding(
    category="prophage",
    title="Novel phage tail protein cluster",
    description="33kb prophage region with uncharacterized tail fiber (3140aa)",
    proteins=["WJJK01000031.1_9", "WJJK01000031.1_12"],
    evidence={"foldseek_hit": "tip_attachment_protein", "evalue": "5.7e-84"},
    location="WJJK01000031.1:1-33000",
    priority="high"
)
```

#### 3. Update State After Each Search
```python
# After completing a search
state["searches_completed"].append("giant+unannotated")
save_state(state)

# After browsing a contig
state["contigs_browsed"].append("WJJK01000031.1")
save_state(state)
```

#### 4. Update Summary Periodically
```python
update_summary()  # Regenerates exploration_summary.md
```

### Core Sharur Operations

```python
# Dataset overview (do this once, save key stats to state)
overview = b.overview()
print(overview.data)

# Search by predicates (combine with has/lacks)
results = b.search_by_predicates(has=["giant", "unannotated"], lacks=["transposase"], limit=50)
print(results.data)

# Genomic neighborhood (window = proteins on each side)
neighborhood = b.get_neighborhood("protein_id", window=10)
print(neighborhood.data)

# Protein details
protein = b.get_protein("protein_id", verbosity=2)
print(protein.data)

# Find similar proteins by embedding
similar = b.find_similar("protein_id", k=10)
print(similar.data)

# Structure prediction (requires ESM_API_KEY, max 1024 aa)
structure = b.predict_structure("protein_id")
print(structure.data)
# Returns: structure.raw["pdb_path"]

# Foldseek search
foldseek = b.search_foldseek("/path/to/structure.pdb")
print(foldseek.data)

# Visualization (always include title and legend)
viz = b.visualize_neighborhood(
    "protein_id", window=10,
    title="Locus Description",
    legend="Genomic context around protein_id showing key functional genes."
)
print(viz.raw["image_path"])  # Path to PNG
```

### Available Predicates - Full Vocabulary

Use `search_by_predicates(has=[...], lacks=[...])` with ANY combination of these predicates. The system is hierarchical - searching for a parent (e.g., `oxidoreductase`) includes all children (e.g., `dehydrogenase`, `reductase`).

**Size predicates:**
- `tiny` (<50aa), `small` (50-150), `medium` (150-400), `large` (400-1000), `giant` (>1000), `massive` (>2000)

**Annotation status:**
- `unannotated` - No database hits at all
- `hypothetical` - Best hit is DUF/hypothetical/uncharacterized
- `well_annotated` - 3+ confident annotations
- `confident_hit` - E-value < 1e-10
- `weak_hit` - Only weak annotations
- `multi_domain` - 3+ distinct Pfam domains
- `single_domain` - Only one domain

**Enzymes (EC-based hierarchy):**
- `oxidoreductase` → `dehydrogenase`, `oxidase`, `reductase`, `peroxidase`, `oxygenase`, `monooxygenase`, `dioxygenase`
- `transferase` → `kinase`, `methyltransferase`, `acetyltransferase`, `glycosyltransferase`, `aminotransferase`
- `hydrolase` → `protease`, `peptidase`, `nuclease`, `dnase`, `rnase`, `lipase`, `esterase`, `glycosidase`, `atpase`, `gtpase`, `phosphatase`
- `lyase`, `isomerase`, `ligase` → `synthase`, `synthetase`, `carboxylase`
- `translocase`

**Transport & localization:**
- `transporter` → `abc_transporter`, `mfs_transporter`, `ion_channel`, `porin`, `efflux_pump`, `symporter`, `antiporter`
- Substrate: `sugar_transporter`, `amino_acid_transporter`, `peptide_transporter`, `ion_transporter`, `metal_transporter`, `phosphate_transporter`
- Localization: `membrane`, `transmembrane`, `secreted`, `periplasmic`, `cell_surface`, `outer_membrane`, `inner_membrane`

**Regulation & signaling:**
- `regulator` → `transcription_factor`, `repressor`, `activator`
- `two_component` → `sensor_kinase`, `response_regulator`
- `signaling` → `cyclic_dinucleotide`, `diguanylate_cyclase`, `phosphodiesterase`, `serine_threonine_kinase`
- Families: `lysr_family`, `tetr_family`, `gntr_family`, `arac_family`, `marr_family`, `laci_family`

**Metabolism - Energy:**
- `energy_metabolism` → `respiration`, `aerobic_respiration`, `anaerobic_respiration`, `atp_synthesis`, `fermentation`, `electron_transport`
- `hydrogen_metabolism` → `hydrogenase`, `nife_hydrogenase`, `fefe_hydrogenase`, `hydrogenase_maturation`

**Metabolism - Carbon:**
- `central_metabolism` → `glycolysis`, `gluconeogenesis`, `tca_cycle`, `pentose_phosphate`, `glyoxylate_cycle`
- `carbon_fixation` → `calvin_cycle`, `reverse_tca`, `wood_ljungdahl`, `3hp_bicycle`
- `one_carbon_metabolism` → `methanogenesis`, `methanotrophy`, `methylotrophy`, `co_oxidation`

**Metabolism - Nitrogen/Sulfur:**
- `nitrogen_metabolism` → `nitrogen_fixation`, `nitrification`, `denitrification`, `ammonia_assimilation`, `nitrate_reduction`
- `sulfur_metabolism` → `sulfate_reduction`, `sulfur_oxidation`, `sulfur_assimilation`

**Metabolism - Biosynthesis:**
- `amino_acid_metabolism`, `lipid_metabolism`, `nucleotide_metabolism`
- `cofactor_biosynthesis` → `nad_biosynthesis`, `fad_biosynthesis`, `cobalamin_biosynthesis`, `folate_biosynthesis`, `heme_biosynthesis`, `iron_sulfur_biosynthesis`
- `secondary_metabolism` → `polyketide_synthesis`, `nrps`, `terpene_synthesis`, `siderophore_biosynthesis`

**CAZy (carbohydrate-active):**
- `carbohydrate_active` → `glycoside_hydrolase`, `glycosyltransferase_cazy`, `polysaccharide_lyase`, `carbohydrate_esterase`, `carbohydrate_binding`
- Specific: `cellulase`, `chitinase`, `amylase`, `xylanase`, `pectinase`

**Binding domains:**
- `dna_binding` → `helix_turn_helix`, `zinc_finger`, `winged_helix`
- `rna_binding`
- `nucleotide_binding` → `atp_binding`, `gtp_binding`, `p_loop`, `aaa_domain`
- `cofactor_binding` → `nad_binding`, `fad_binding`, `plp_binding`, `coenzyme_a_binding`
- `metal_binding` → `iron_binding`, `iron_sulfur`, `heme_binding`, `zinc_binding`, `copper_binding`, `nickel_binding`, `molybdenum_binding`

**Cell envelope & surface:**
- `cell_wall` → `peptidoglycan`, `murein_synthesis`, `pbp`
- `lps_biosynthesis`, `capsule`, `exopolysaccharide`
- `adhesin`, `pilus` → `type_iv_pilus`
- `flagellum` → `flagellar_motor`, `flagellin`
- `chemotaxis`, `s_layer`

**Mobile elements & defense:**
- `mobile_element` → `transposase`, `integrase`, `resolvase`, `recombinase`, `insertion_sequence`
- `phage_related` → `phage_integrase`, `phage_terminase`, `phage_portal`, `phage_capsid`, `phage_tail`, `phage_lysin`, `holin`
- `defense_system` → `restriction_modification`, `restriction_enzyme`, `crispr_associated`, `cas_nuclease`, `toxin_antitoxin`, `toxin`, `antitoxin`, `abortive_infection`
- `secretion_system` → `type_i_secretion` through `type_vi_secretion`, `sec_pathway`, `tat_pathway`
- `conjugation` → `relaxase`, `type_iv_secretion`

**Stress & resistance:**
- `stress_response` → `heat_shock`, `cold_shock`
- `oxidative_stress` → `catalase`, `superoxide_dismutase`, `peroxiredoxin`, `thioredoxin`, `glutaredoxin`
- `chaperone` → `hsp70`, `hsp60`, `hsp90`, `small_hsp`, `clp_protease`, `lon_protease`
- `antibiotic_resistance` → `beta_lactamase`, `aminoglycoside_resistance`, `multidrug_resistance`
- `heavy_metal_resistance` → `arsenic_resistance`, `mercury_resistance`, `copper_resistance`
- `dna_repair` → `base_excision_repair`, `nucleotide_excision_repair`, `mismatch_repair`, `recombinational_repair`, `sos_response`

**Structural features:**
- `repeat_domain` → `tpr_repeat`, `wd40_repeat`, `lrr_repeat`, `ankyrin_repeat`, `kelch_repeat`
- `coiled_coil`, `beta_barrel`, `intrinsically_disordered`

**Information processing:**
- `replication` → `dna_polymerase`, `helicase`, `primase`, `topoisomerase`, `ligase_dna`
- `transcription` → `rna_polymerase`, `transcription_elongation`, `transcription_termination`
- `translation` → `ribosomal_protein`, `trna_synthetase`, `translation_factor`

**Composition:**
- `gc_outlier` - GC >2 std from genome mean (possible HGT)
- `gc_high` (>60%), `gc_low` (<40%)

**Topology (predicted):**
- `transmembrane_predicted`, `single_pass_membrane`, `multi_pass_membrane`, `polytopic_membrane`
- `soluble_predicted`

### Exploration Strategies Using Predicates

**Discovery-driven combinations (find the unexpected):**
```python
# Defense and mobile elements
b.search_by_predicates(has=["gc_outlier", "defense_system"])
b.search_by_predicates(has=["gc_outlier", "mobile_element"])

# Complex unknowns - multiple domains but poorly understood
b.search_by_predicates(has=["multi_domain"], lacks=["confident_hit"])

# Membrane proteins of unknown function
b.search_by_predicates(has=["transmembrane_predicted", "unannotated"])

# Surface biology
b.search_by_predicates(has=["adhesin"])
b.search_by_predicates(has=["s_layer"])
b.search_by_predicates(has=["type_iv_pilus"])

# Energy metabolism
b.search_by_predicates(has=["hydrogenase"])
b.search_by_predicates(has=["electron_transport", "iron_sulfur"])

# Defense systems diversity
b.search_by_predicates(has=["crispr_associated"])
b.search_by_predicates(has=["toxin_antitoxin"])
b.search_by_predicates(has=["restriction_modification"])

# Unusual enzyme combinations
b.search_by_predicates(has=["oxidoreductase", "giant"])
b.search_by_predicates(has=["glycosyltransferase", "multi_domain"])
```

**Avoid premature conclusions** - use predicates to FIND proteins, then examine context before interpreting.

---

## Hypothesis-Driven Exploration with Subagents

After completing Phase 0-3 (see below) and building an initial understanding, **generate testable hypotheses** and spawn subagents to investigate them.

### Workflow

1. **Complete systematic browsing first** (Phases 0-3)
2. **Identify puzzles and patterns** that need deeper investigation
3. **Formulate specific hypotheses** to test
4. **Spawn focused subagents** to investigate each hypothesis

### Example Hypothesis Generation (from Survey Results)

Survey found:
- "15 genomes have dockerin but no cohesin"
- "2 genomes lack hydrogenases but have high cytochrome counts"
- "Giant 1878aa protein with both dockerin and MtrC domains"

**Hypotheses to test:**
- H1: Orphan dockerin genomes are syntrophic partners that dock to other organisms
- H2: Cytochrome-rich genomes use DIET instead of H2 for electron transfer
- H3: Dockerin-cytochrome fusion enables partner-specific electron conduits

### Spawning Hypothesis-Testing Subagents

```python
# Example: Test if orphan dockerin genomes have distinct metabolic features
Task(
    subagent_type="general-purpose",
    description="Test orphan dockerin hypothesis",
    prompt="""You are a leaf agent testing a specific hypothesis.

**CRITICAL: You are a LEAF AGENT. DO NOT spawn sub-agents or use the Task tool.**

Dataset: data/hinthialibacterota_v3/sharur.duckdb

Hypothesis: "Orphan dockerin genomes (dockerin but no cohesin) are obligate syntrophic partners with distinct metabolic features compared to self-assembly genomes."

Your task:
1. Identify orphan vs self-assembly genomes
2. Compare their metabolic features (cytochromes, transporters, genome size)
3. Test for statistical differences (Fisher exact, t-tests)
4. Generate neighborhood figures for representative dockerin-cytochrome proteins
5. Document findings with evidence

Write results to: data/hinthialibacterota_v3/exploration/hypothesis_orphan_dockerins.md

Include:
- Statistical tests with p-values
- At least 2 gene neighborhood figures (save to exploration/figures/)
- Clear conclusion: hypothesis supported/rejected/needs modification"""
)

# Wait for completion, review results, spawn next hypothesis test
```

### Hypothesis Tracking API

Use the persistent hypothesis registry to track reasoning across sessions:

```python
# Propose a hypothesis (persists to exploration/hypotheses.json)
h = b.propose_hypothesis("Group 4 NiFe hydrogenases are energy-conserving in this lineage")

# Do analysis...
result = b.search_by_predicates(has=["nife_group4"])

# Link evidence to hypothesis
b.add_evidence(
    h.hypothesis_id,
    query="Search for NiFe Group 4 across all genomes",
    result_summary=f"Found in {len(result.data)} genomes",
    supports=True,
    confidence=0.8,
)

# Check hypothesis state
print(b.hypothesis_summary())

# Generate provenance figure for the paper
mermaid = b.render_provenance(
    title="Hydrogenase Analysis",
    output_path="figures/hydrogenase_provenance.mermaid",
)
```

Hypotheses persist across sessions -- `b.resume()` shows active hypotheses automatically. Use `b.list_hypotheses()` to get full hypothesis objects.

For explicit provenance chaining (building a DAG of analytical steps):

```python
e1 = b.log_provenance("Count hydrogenases", "42 found")
e2 = b.log_provenance("Check neighborhoods", "12 with Hyf operon", parent_ids=[e1.entry_id])
e3 = b.log_provenance("Statistical test", "p=0.003", parent_ids=[e1.entry_id, e2.entry_id])

# Render the full DAG including hypotheses
mermaid = b.render_provenance(title="Analysis", output_path="figures/provenance.mermaid")
```

### What Makes a Good Hypothesis for Subagent Investigation?

✅ **Good (specific, testable)**:
- "Genomes lacking hydrogenases compensate with higher cytochrome investment"
- "Giant unannotated proteins in synteny with cohesins are scaffolding components"
- "GC-outlier regions contain recently acquired defense systems"

❌ **Bad (too broad, unfocused)**:
- "Characterize defense systems"
- "Understand metabolism"
- "Find interesting proteins"

### Subagent Outputs Should Include:

- **Markdown report**: hypothesis_[topic].md with clear conclusion
- **Figures**: Gene neighborhoods, statistical plots (save to exploration/figures/)
- **Evidence**: Protein IDs, statistical tests, comparative data
- **Conclusion**: Hypothesis supported/rejected/needs modification/unclear

### Integration Back to Main Exploration

After subagents complete:
1. **Review their findings** and integrate into your synthesis
2. **Generate additional hypotheses** if needed
3. **Document what changed your mind** if initial intuition was wrong
4. **Build final explanatory model** incorporating all evidence

---

### Exploration Workflow (Four Phases)

**All phases are required for a complete exploration.**

**Critical mindset**: While following these phases, you are building an **explanatory model**, not just collecting data:
- Form hypotheses about what organizes the dataset (metabolic strategies? defense tradeoffs? ecological specialization?)
- Test competing explanations - when you see a pattern, consider multiple mechanisms
- Look for counterevidence - find exceptions to your model and understand why
- Document what changed your mind - if initial intuition was wrong, explain the pivot
- Build toward a narrative - your synthesis should explain dominant patterns and meaningful exceptions

The phases below provide the **workflow** for discovering loci and gathering evidence. Your **thinking** should be focused on interpretation and model-building throughout.

#### Phase 0: Dataset Assessment & Genome Profiling (do this FIRST)

Before diving into specific searches, understand the dataset structure:

```python
# Check how many genomes we have
n_genomes = b.store.execute("SELECT COUNT(DISTINCT bin_id) FROM proteins")[0][0]
print(f"Dataset contains {n_genomes} genomes")

# If multiple genomes, generate comparative profiles
if n_genomes > 1 and 'genome_profiles_generated' not in state.get('searches_completed', []):
    print("Multiple genomes detected - generating comparative profiles...")
    # See "Genome Profiling Workflow" section below
```

**For multi-genome datasets (>1 genome), ALWAYS generate genome profiles first.** This gives you:
- Bird's-eye view of functional variation across genomes
- Identifies outlier genomes to investigate
- Shows which features are core vs variable
- Prevents over-generalizing from a single genome's features

Skip to Phase 1 only for single-genome datasets.

#### Phase 1: Overview & Quick Searches (~20% of time)

1. **Check state** - Don't repeat work already done
2. **Run overview** - Understand dataset size and composition
3. **Quick predicate searches** - Find known interesting categories:
   - `has=['crispr_associated']` - Defense systems
   - `has=['hydrogenase']` - Energy metabolism
   - `has=['adhesin']` or `has=['s_layer']` - Surface biology
   - `has=['mobile_element']` - Genomic plasticity
   - `has=['multi_domain', 'unannotated']` - Complex unknowns
4. **Log findings immediately** - Don't accumulate in memory
5. **Note which contigs have interesting proteins** - For Phase 2

#### Phase 2: Genome-by-Genome Browsing (~50% of time) - THIS IS THE CORE

**This is where you discover biology.** Don't just check predicates - actually READ the genomes.

**BATCH APPROACH:** Read whole contigs via SQL, not individual API calls. This is faster and gives you context.

```python
# For each genome, browse its contigs in batches
genomes = b.store.execute("""
    SELECT bin_id, COUNT(*) as n_proteins, COUNT(DISTINCT contig_id) as n_contigs
    FROM proteins GROUP BY bin_id ORDER BY n_proteins DESC
""")

for bin_id, n_proteins, n_contigs in genomes:
    if bin_id in state.get('genomes_browsed', []):
        continue

    print(f"\n{'='*60}")
    print(f"BROWSING GENOME: {bin_id}")
    print(f"  {n_proteins} proteins across {n_contigs} contigs")
    print('='*60)

    # Get ALL proteins for this genome in one query (batch read)
    genome_proteins = b.store.execute("""
        SELECT p.protein_id, p.contig_id, p.gene_index, p.sequence_length,
               p.strand, pp.predicates
        FROM proteins p
        LEFT JOIN protein_predicates pp ON p.protein_id = pp.protein_id
        WHERE p.bin_id = ?
        ORDER BY p.contig_id, p.gene_index
    """, [bin_id])

    # Get annotations for this genome in one query (batch read)
    genome_annotations = b.store.execute("""
        SELECT a.protein_id, a.source, a.accession, a.name, a.evalue
        FROM annotations a
        JOIN proteins p ON a.protein_id = p.protein_id
        WHERE p.bin_id = ?
        ORDER BY a.protein_id, a.evalue
    """, [bin_id])

    # Group annotations by protein
    annot_by_protein = {}
    for row in genome_annotations:
        pid = row[0]
        if pid not in annot_by_protein:
            annot_by_protein[pid] = []
        annot_by_protein[pid].append(row[1:])  # (source, accession, name, evalue)

    # Now browse the genome contig by contig
    current_contig = None
    contig_proteins = []
    interesting_loci = []  # Collect for this genome

    for row in genome_proteins:
        pid, contig, gene_idx, length, strand, predicates = row
        annotations = annot_by_protein.get(pid, [])

        if contig != current_contig:
            # Process previous contig if exists
            if contig_proteins:
                loci = analyze_contig(current_contig, contig_proteins)
                interesting_loci.extend(loci)
            current_contig = contig
            contig_proteins = []

        contig_proteins.append({
            'protein_id': pid,
            'gene_index': gene_idx,
            'length': length,
            'strand': strand,
            'predicates': predicates or [],
            'annotations': annotations
        })

    # Don't forget last contig
    if contig_proteins:
        loci = analyze_contig(current_contig, contig_proteins)
        interesting_loci.extend(loci)

    # === PER-GENOME SYNTHESIS ===
    synthesize_genome_findings(bin_id, interesting_loci)

    state.setdefault('genomes_browsed', []).append(bin_id)
    save_state(state)


def analyze_contig(contig_id, proteins):
    """Scan a contig for interesting loci. Returns list of loci."""
    loci = []

    # Look for clusters of interest
    for i, p in enumerate(proteins):
        # Unannotated clusters (3+ in a row)
        if 'unannotated' in p['predicates']:
            cluster = get_unannotated_cluster(proteins, i)
            if len(cluster) >= 3:
                loci.append(('unannotated_cluster', cluster))

        # Defense islands
        if any(pred in p['predicates'] for pred in ['crispr_associated', 'restriction_modification', 'toxin_antitoxin']):
            island = get_defense_island(proteins, i)
            if len(island) >= 2:
                loci.append(('defense_island', island))

        # Giant proteins always interesting
        if 'giant' in p['predicates']:
            context = proteins[max(0,i-5):min(len(proteins),i+6)]
            loci.append(('giant_protein', context))

        # Prophage regions (integrase + structural)
        if 'phage_related' in p['predicates'] or 'integrase' in p['predicates']:
            region = get_prophage_region(proteins, i)
            if len(region) >= 3:
                loci.append(('prophage_region', region))

    return loci


def synthesize_genome_findings(bin_id, loci):
    """Synthesize findings for one genome and log them."""
    if not loci:
        return

    # Group by type
    by_type = {}
    for locus_type, proteins in loci:
        by_type.setdefault(locus_type, []).append(proteins)

    # Log summary for this genome
    summary_parts = []
    for locus_type, instances in by_type.items():
        summary_parts.append(f"{len(instances)} {locus_type.replace('_', ' ')}(s)")

    log_finding(
        category="genome_summary",
        title=f"Genome {bin_id}: {', '.join(summary_parts)}",
        description=f"Browsed {bin_id}, found {len(loci)} interesting loci",
        evidence={"genome": bin_id, "loci_counts": {k: len(v) for k,v in by_type.items()}},
        priority="medium"
    )

    # Log individual high-value loci
    for locus_type, proteins in loci:
        if locus_type == 'giant_protein' or len(proteins) >= 5:
            log_finding(
                category=locus_type,
                title=f"{locus_type.replace('_', ' ').title()} in {bin_id}",
                description=describe_locus(locus_type, proteins),
                proteins=[p['protein_id'] for p in proteins[:10]],
                location=f"{proteins[0]['protein_id'].rsplit('_',1)[0]}",
                priority="high" if locus_type == 'giant_protein' else "medium"
            )

            # Generate neighborhood visualization for high-value loci
            # Pick a representative protein (usually central or most interesting)
            if locus_type == 'giant_protein':
                center_protein = proteins[len(proteins)//2]['protein_id']
            else:
                # For clusters, pick the most annotated or central protein
                center_protein = proteins[0]['protein_id']

            figure_path = f"exploration/figures/{bin_id}_{locus_type}_{len(by_type[locus_type])}_neighborhood.png"
            EXPLORE_DIR.parent / "exploration" / "figures").mkdir(parents=True, exist_ok=True)

            viz_result = b.visualize_neighborhood(
                center_protein,
                window=12,  # Show ±12 genes for context
                output_path=str(EXPLORE_DIR.parent / figure_path),
                title=f"{locus_type.replace('_', ' ').title()} in {bin_id}",
                legend=f"{describe_locus(locus_type, proteins)} Neighborhood of "
                    f"{center_protein} showing ±12 genes of genomic context."
            )

            # You don't need to generate a figure for EVERY locus, but do generate them for:
            # - Giant proteins (always interesting)
            # - Large clusters (≥5 genes)
            # - Representative examples of each locus type
            # - Unusual or unexpected arrangements
```

**What to look for while browsing:**
- **Unannotated clusters** - 3+ unannotated genes in a row = potential novel operon
- **Mixed clusters** - Known enzyme flanked by unknowns = pathway with novel components
- **Defense islands** - CRISPR, restriction enzymes, toxins clustered together
- **Prophage regions** - Integrases, structural genes, lysins
- **Surface operons** - Adhesins, S-layer, secretion systems
- **Metabolic operons** - Enzymes that work together in a pathway
- **Giant proteins** - Always examine context of proteins >1000aa

#### Phase 3: Cross-Genome Synthesis & Deep Dives (~30% of time)

After browsing all (or many) genomes, synthesize patterns ACROSS genomes:

```python
# Cross-genome synthesis
def cross_genome_synthesis():
    """Identify patterns that span multiple genomes."""

    # Which locus types appear in which genomes?
    findings = load_findings()
    genome_summaries = [f for f in findings if f['category'] == 'genome_summary']

    # Aggregate locus counts across genomes
    locus_totals = {}
    genomes_with_locus = {}
    for f in genome_summaries:
        genome = f['evidence'].get('genome')
        for locus_type, count in f['evidence'].get('loci_counts', {}).items():
            locus_totals[locus_type] = locus_totals.get(locus_type, 0) + count
            genomes_with_locus.setdefault(locus_type, set()).add(genome)

    # Log cross-genome patterns
    for locus_type, total in sorted(locus_totals.items(), key=lambda x: -x[1]):
        n_genomes = len(genomes_with_locus[locus_type])
        conservation = n_genomes / len(genome_summaries) * 100

        log_finding(
            category="cross_genome_pattern",
            title=f"{locus_type.replace('_', ' ').title()}: {total} instances across {n_genomes} genomes ({conservation:.0f}%)",
            description=f"Found {total} {locus_type} loci across {n_genomes} genomes",
            evidence={
                "locus_type": locus_type,
                "total_instances": total,
                "n_genomes": n_genomes,
                "conservation_pct": round(conservation, 1)
            },
            priority="high" if conservation > 50 else "medium"
        )

cross_genome_synthesis()
```

**Deep dives** on the most interesting cross-genome patterns:
- For highly conserved loci (>80% genomes): What's the core vs variable component?
- For rare loci (<20% genomes): Is this contamination, HGT, or niche-specific?

**Visualization is mandatory** - generate neighborhood diagrams for representative examples:

```python
# Example: Visualize a conserved defense locus
for genome in representative_genomes:
    # Find a representative protein from the locus
    defense_protein = b.search_by_predicates(
        has=["crispr_associated"],
        bin_id=genome
    ).data[0]['protein_id']

    # Generate neighborhood figure with title and legend
    b.visualize_neighborhood(
        defense_protein,
        window=15,
        output_path=f"exploration/figures/{genome}_defense_locus.png",
        title=f"Defense Locus in {genome}",
        legend=f"CRISPR-Cas and associated defense genes in {genome}. "
            "Red arrow marks the query protein. Color indicates functional category."
    )
```

Don't just generate 3-5 figures and stop - that's treating minimums as targets. Generate figures for:
- Each major locus type you document
- Representative examples from multiple genomes (show conservation/variation)
- Anything unusual or unexpected
- Key findings that support your synthesis

### Minimum Requirements for Complete Exploration

- [ ] Dataset overview recorded
- [ ] **Genome profiles generated** (if multi-genome dataset)
- [ ] **At least 10 genomes browsed gene-by-gene** (batch SQL, not individual API calls)
- [ ] **Per-genome synthesis** for each browsed genome
- [ ] **Cross-genome synthesis** identifying patterns across genomes
- [ ] At least 5 predicate searches for targeted discovery
- [ ] At least 30 findings logged (including genome summaries)
- [ ] **Conservation reported** for all major findings (N genomes, % conservation)
- [ ] **At least 3-5 loci visualized** (preferably more based on dataset size and complexity - don't treat minimums as targets)
- [ ] PDF report generated

### What "Thorough" Means

A thorough exploration answers these questions:

#### 1. What's there? (Discovery)
- Major protein families and domains
- Metabolic capabilities
- Defense systems
- Mobile elements
- Surface/adhesion systems

#### 2. How conserved is it? (Conservation)
For EVERY major finding, report:
- N genomes with the feature
- N proteins total
- % conservation
```python
# Template query for conservation
b.store.execute("""
    SELECT COUNT(DISTINCT bin_id) as genomes, COUNT(*) as proteins
    FROM proteins p
    JOIN annotations a ON p.protein_id = a.protein_id
    WHERE a.name LIKE '%{feature}%' OR a.accession = '{accession}'
""")
```

#### 3. What's the architecture? (Characterization)
For interesting systems (adhesion, defense, etc.):
- Full list of components
- Domain architectures of key proteins
- Genomic organization (operon structure)

#### 4. What's missing? (Negative findings)
Explicitly check for and report absence of:
- Key pathway genes (MCR, PRK, etc.)
- Expected systems for the environment
- Genes present in related organisms

#### 5. What's unique? (Novelty)
- Lineage-specific expansions (e.g., 244 DUF11 proteins)
- Unusual domain combinations
- Conserved proteins of unknown function

#### 6. What remains unknown? (Gaps)
End every exploration with:
- List of unanswered questions
- Suggested follow-up analyses
- Key hypotheses to test experimentally

### Advanced Analysis (When Time Permits)

#### Genome-level Statistics
```python
# Proteins per genome
b.store.execute("""
    SELECT bin_id, COUNT(*) as n_proteins,
           AVG(sequence_length) as avg_length
    FROM proteins
    GROUP BY bin_id
    ORDER BY n_proteins DESC
""")

# Annotation rate per genome
b.store.execute("""
    SELECT p.bin_id,
           COUNT(DISTINCT p.protein_id) as total,
           COUNT(DISTINCT a.protein_id) as annotated
    FROM proteins p
    LEFT JOIN annotations a ON p.protein_id = a.protein_id
    GROUP BY p.bin_id
""")
```

#### Domain Co-occurrence
What domains appear together?
```python
# Find proteins with both Dockerin AND Cohesin (unusual)
b.store.execute("""
    SELECT p.protein_id, p.sequence_length
    FROM proteins p
    WHERE EXISTS (SELECT 1 FROM annotations a WHERE a.protein_id = p.protein_id AND a.name LIKE '%Dockerin%')
      AND EXISTS (SELECT 1 FROM annotations a WHERE a.protein_id = p.protein_id AND a.name LIKE '%Cohesin%')
""")
```

#### Size Distribution by Function
Are certain functions associated with giant proteins?
```python
b.store.execute("""
    SELECT a.name,
           AVG(p.sequence_length) as avg_size,
           MAX(p.sequence_length) as max_size,
           COUNT(*) as count
    FROM proteins p
    JOIN annotations a ON p.protein_id = a.protein_id
    GROUP BY a.name
    HAVING COUNT(*) > 10
    ORDER BY avg_size DESC
    LIMIT 20
""")
```

#### Finding Protein Families via Similarity
Use embeddings to find related unannotated proteins:
```python
# Find proteins similar to an interesting one
similar = b.find_similar("interesting_protein_id", k=50)
# Are they all unannotated? Might be a novel family
```

#### Cross-Genome Presence/Absence
Which genomes have which systems?
```python
# Presence/absence matrix for a set of features
features = ['Dockerin', 'CRISPR', 'Hydrogenase', 'RuBisCO']
for feature in features:
    result = b.store.execute(f"""
        SELECT COUNT(DISTINCT bin_id) FROM proteins p
        JOIN annotations a ON p.protein_id = a.protein_id
        WHERE a.name LIKE '%{feature}%'
    """)
    print(f"{feature}: {result[0][0]} genomes")

### What Makes a Locus Interesting?

- **Defense system islands** (CRISPR, R-M, TA modules)
- **Metabolic operons** with uncharacterized accessory genes
- **Cell adhesion systems** (Dockerin, Cohesin, HYR, PKD domains)
- **Prophage regions** with uncharacterized structural proteins
- **Gene clusters** with mixed annotation (known + unknown together)
- **Unusual domain combinations** (DUFs fused to known domains)
- **Surface proteins** (DUF11, CarboxypepD_reg, secretion signals)
- **Conserved unknowns** - unannotated proteins present across many genomes

### Context-First Analysis Protocol

**The domain tells you the fold. The neighbors tell you the function.**

Annotation databases assign labels based on structural fold similarity, not biological function. When a PFAM HMM or KOfam profile is broad enough to match an entire enzyme superfamily, the hit count is high but the name is misleading. A BcrAD_BadFG domain next to Rieske and cytochrome is a radical enzyme activase. The same domain next to BcrBC catalytic subunits would be benzoyl-CoA reductase. The domain is the same; the context determines function.

#### Superfamily Awareness Rule

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

**Real examples of this trap:**

| Annotation | Hits/genome | What name says | What it actually is |
|-----------|-------------|---------------|---------------------|
| K00134 (GAPDH) | ~1 | Glyceraldehyde-3-P dehydrogenase | Correct — specific enzyme |
| PF01869 (BcrAD_BadFG) | ~1.7 | Benzoyl-CoA reductase | Fe-protein ATPase superfamily (misleading!) |
| K23108 (benzoate-CoA ligase) | ~11.5 | Benzoate-CoA ligase | Adenylation domain superfamily (wrong!) |
| PF04055 (Radical_SAM) | ~43.7 | Radical SAM | Huge superfamily — zero pathway info |

#### Co-Annotation Validation

Before claiming a protein functions as enzyme Y, check what *other* annotations it carries. If a "benzoate-CoA ligase" (K23108) is co-annotated with NRPS adenylation domains, it's an adenylation domain enzyme, not a benzoate-specific ligase.

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

#### Neighborhood Validation

For any protein-level functional claim that implies a pathway or system, use `b.get_neighborhood()` with `all_annotations=True` to see the full annotation context of surrounding genes:

```python
# View all annotations for genes ±5 around a protein of interest
result = b.get_neighborhood(protein_id, window=5, all_annotations=True)
print(result)  # Shows every annotation source per gene

# Programmatic access to neighbor annotations
for p in result._raw['proteins']:
    print(p['protein_id'], p['length_aa'], 'aa')
    for source, annots in p.get('annotations', {}).items():
        for a in annots:
            print(f"  {source}: {a['accession']} {a['name']}")
```

**What to look for:**
- **Neighbors consistent with claimed function** → Claim supported
- **Neighbors suggest a different system** → Revise interpretation
- **No informative neighbors** (isolated gene, unannotated neighbors) → State uncertainty
- **Protein at contig edge** → Note potential fragmentation

#### Claim Escalation Ladder

Different claims require different evidence levels:

| Claim Level | Evidence Required | Example |
|------------|-------------------|---------|
| "Protein contains domain X" | Accession name verification only | "Contains Radical_SAM domain" |
| "Protein functions as enzyme Y" | Name verification + co-annotation check | "Functions as 2-hydroxyacyl-CoA dehydratase activase (K23876)" |
| "Genome encodes pathway Z" | Multiple pathway-specific markers + co-localization check | "Encodes 7/10 glycolysis enzymes" |
| "Phylum universally performs W" | All above + conservation analysis across genomes | "41/41 genomes encode NiFe hydrogenases (validated by PF00374)" |

**For pathway claims specifically:**
- A single marker gene does NOT prove a pathway
- Divergent organisms may lack detectable homologs for some steps even when the function is present — report "X/Y steps detected" rather than "pathway incomplete"
- Check whether pathway components are co-localized (same contig, within ±10 genes) — scattered components deserve more skepticism
- Use the KEGG REST API to check what pathway a KO belongs to and what the other required KOs are (see KEGG Pathway Lookup below)

#### KEGG Pathway Lookup

When you find a KEGG KO and want to understand its pathway context, use the KEGG REST API:

```python
# What pathways contain this KO?
import subprocess
result = subprocess.run(
    ["curl", "-s", "https://rest.kegg.jp/link/pathway/ko:K23108"],
    capture_output=True, text=True
)
print(result.stdout)
# Returns: ko:K23108    path:map00362  (benzoate degradation)

# What are ALL the KOs in that pathway?
result = subprocess.run(
    ["curl", "-s", "https://rest.kegg.jp/link/ko/path:map00362"],
    capture_output=True, text=True
)
# Returns every KO in the benzoate degradation pathway
# → Check which ones are present in your dataset

# What KEGG modules (minimum functional units) use this KO?
result = subprocess.run(
    ["curl", "-s", "https://rest.kegg.jp/link/module/ko:K23108"],
    capture_output=True, text=True
)
# Returns: ko:K23108    md:M00551 (benzoate degradation module)

# Get the full module definition (what's required for function)
result = subprocess.run(
    ["curl", "-s", "https://rest.kegg.jp/get/md:M00551"],
    capture_output=True, text=True
)
# Shows all steps in the module and alternative enzymes per step
```

**KEGG modules are more useful than pathways for metabolic claims** because they define minimum functional units. A complete module = a testable assertion ("does this genome have the minimum gene set for this function?").

**Important caveats for divergent organisms:**
- KOfam profiles are built from characterized enzymes, mostly from model organisms
- Deeply branching lineages (Asgard archaea, candidate phyla) may have homologs too divergent for detection
- When reporting pathway completeness, state: "X/Y steps detected; undetected steps may be present but too divergent for KOfam detection"
- Missing steps in well-conserved pathways (e.g., glycolysis) are more likely genuine absences than in less-conserved systems

**Local reference:** `data/reference/ko_list` contains all 27,324 KO definitions with EC numbers. `scripts/pathway_completeness.py` defines 12 major pathways with component KOs.

---

### Scientific Validation Checklist

**IMPORTANT: Don't overclaim. Validate findings before logging them as confirmed.**

#### MANDATORY: Accession-to-Function Verification

**NEVER assume a PFAM/KEGG accession encodes a specific function from memory.** Always verify against the database `name` field. This rule exists because a major false finding (PF04055 claimed as benzoyl-CoA reductase, actually Radical_SAM) propagated through an entire analysis and into a manuscript before being caught.

**Before logging ANY finding that references a PFAM/KEGG accession:**

```python
# STEP 1: Verify the accession name in the database
name_check = b.store.execute("""
    SELECT DISTINCT name FROM annotations WHERE accession = 'PF04055'
""").fetchall()
print(f"PF04055 = {name_check}")  # → [('Radical_SAM',)]
# If the name doesn't match your assumption, STOP. Your claim is wrong.

# STEP 2: Include the verified name in your finding's provenance
provenance = {
    "query": "SELECT COUNT(*), name FROM annotations WHERE accession = 'PF04055' GROUP BY name",
    "raw_result": [(1790, "Radical_SAM")],
    "accession_verified": "PF04055 = 'Radical_SAM' — this is a generic superfamily, NOT benzoyl-CoA reductase",
    "interpretation": "1,790 Radical SAM proteins (normal for anaerobic bacteria, not pathway-specific)"
}
```

**Red flags that should trigger extra verification:**
- Accession you haven't seen before in this dataset
- Claim that a single PFAM domain proves a complete pathway
- "Universal" claims (present in >90% of genomes) — especially if the domain is a large superfamily
- Protein counts >500 from a single accession — large superfamilies are usually generic

**The cost of checking is 1 SQL query. The cost of not checking was a retracted finding.**

#### Counting Proteins vs Annotation Rows

**ALWAYS use `COUNT(DISTINCT protein_id)` when reporting protein counts.** Repeat domains (WD40, TPR, FG-GAP, Ig-like, etc.) produce multiple annotation rows per protein. `COUNT(*)` on the annotations table gives you annotation hits, not protein counts — and these can differ by 5-10x for repeat domains.

```python
# ✗ WRONG — counts annotation rows (864 for WD40 in giant proteins)
b.store.execute("SELECT COUNT(*) FROM annotations WHERE accession = 'PF00400'")

# ✓ RIGHT — counts unique proteins (352 for WD40)
b.store.execute("SELECT COUNT(DISTINCT protein_id) FROM annotations WHERE accession = 'PF00400'")
```

When reporting counts, always state what you're counting: "352 proteins with WD40 domains" or "864 WD40 domain instances across 352 proteins."

#### Finding Titles Must Be Self-Contained

When logging a finding, the title must include all qualifiers. A reader should understand the scope from the title alone. Qualifiers dropped during summarization have caused errors.

- **Bad**: "Largest protein (5,461 aa) has zero annotations" — was actually "largest *unannotated* protein"
- **Good**: "Largest unannotated protein (5,461 aa, JAJVIK010000008.1_48) — candidate adhesin"

#### Hydrogenase Curation

HydDB classifies proteins by structural similarity to known hydrogenases, but NADH dehydrogenase (Complex I) shares the same [NiFe] binding site fold and produces ~44% false positives for NiFe calls. The automated pipeline (`classify_hydrogenases.py`) applies a PF00374-based filter that validates Groups 1-3 but **systematically rejects all Group 4 NiFe hydrogenases** (Hyf/Hyc/Mbh/Ech) because they diverged too far to retain the NiFeSe_Hases domain.

**When you encounter hydrogenases, run neighborhood-based curation on unvalidated HydDB hits:**

```python
# Step 1: Get pipeline-validated and unvalidated counts
validated = b.search_by_predicates(has=["nife_hydrogenase"])   # PF00374-validated (Groups 1-3)
all_hyddb = b.search_by_predicates(has=["hyddb:NiFe"])         # All HydDB NiFe calls
unvalidated = [p for p in all_hyddb if p not in set(validated)]

# Step 2: For each unvalidated hit, check ±8 gene neighborhood
for pid in unvalidated:
    nbr = b.get_neighborhood(pid, window=8, all_annotations=True)
    # Hydrogenase evidence (→ rescue):
    #   KEGG: K12136-K12145 (hyfA-J), K15828-K15833 (hycB-G)
    #   KEGG: K04651-K04656 (HypA-F maturation), K03605 (HycI)
    #   PFAM: PF00374 on a neighboring gene, HydDB annotation on neighbor
    # Complex I evidence (→ reject):
    #   KEGG: K00330-K00343 (nuoA-N)
    #   PFAM: PF00346, PF00329 on neighbors
    # Ambiguous (both or neither): reject conservatively
```

**In reports, state the curation method:** "97 NiFe hydrogenases validated by PF00374 (Groups 1-3), plus 12 rescued by neighborhood curation (8 Group 4f Hyf/Hyc + 4 Group 4e operon subunits), excluding 66 Complex I false positives."

#### Pathway Validation
Don't assume a pathway from one gene - check for the complete pathway:

| If you find... | Check for... | Why |
|----------------|--------------|-----|
| RuBisCO | PRK nearby | Without PRK, it's likely nucleotide salvage, not Calvin cycle |
| MCR (methyl-CoM reductase) | - | This IS the definitive methanogenesis marker |
| CODH | ACS, methyltransferases | Complete Wood-Ljungdahl needs all components |
| Nitrogenase (nifH) | nifD, nifK clustered | Need the complete complex |
| Hydrogenase | Neighborhood curation | See "Hydrogenase Curation" above — PF00374 misses Group 4 |

#### Archaeal-Specific Considerations
- **RuBisCO in archaea** is usually Form II/III for AMP recycling, NOT CO2 fixation
- **No peptidoglycan** in most archaea - don't expect MurA-F genes
- **Archaellum ≠ flagellum** - different genes, different mechanism
- **Ether lipids** - look for GGGP synthase, not fatty acid synthesis
- **Histones** common in archaea - don't assume eukaryotic contamination

#### Operon/Clustering Expectations
Functional systems should be clustered. If genes are scattered, be skeptical:

| System | Expected clustering |
|--------|---------------------|
| CRISPR-Cas | cas genes + CRISPR array together |
| Ribosomal proteins | Large operons (S10, spc, alpha) |
| Secretion systems | T2SS/T4SS genes clustered |
| Hydrogenases | Structural + maturation together |
| Biosynthetic pathways | Usually operonic |

#### Negative Results Are Findings
Report what's ABSENT - it's often more informative than what's present:
- No MCR → NOT a methanogen (even if has other "methanogenesis" genes)
- No PRK → RuBisCO probably not for CO2 fixation
- No CRISPR array → cas genes may be orphaned/non-functional

#### Conservation Analysis (ALWAYS DO THIS)

**For every interesting finding, check how conserved it is:**

```python
# How many genomes have this domain/annotation?
conservation = b.store.execute("""
    SELECT COUNT(DISTINCT p.bin_id) as n_genomes,
           COUNT(*) as n_proteins
    FROM proteins p
    JOIN annotations a ON p.protein_id = a.protein_id
    WHERE a.name LIKE '%Dockerin%'
""")
print(f"Found in {conservation[0][0]} of 126 genomes ({conservation[0][1]} proteins)")

# Which genomes have it?
genomes_with = b.store.execute("""
    SELECT DISTINCT p.bin_id
    FROM proteins p
    JOIN annotations a ON p.protein_id = a.protein_id
    WHERE a.name LIKE '%Dockerin%'
    LIMIT 20
""")
```

**Interpretation guide:**
- Present in >80% of genomes → Core gene, essential function
- Present in 20-80% → Variable, possibly niche-specific
- Present in <20% → Rare, possibly HGT or contamination
- Present in 1 genome → Very suspicious, check carefully

**Always report conservation in findings:**
```python
log_finding(
    ...,
    evidence={
        "n_genomes": 45,
        "n_proteins": 127,
        "conservation": "36% of genomes"  # ALWAYS INCLUDE THIS
    }
)
```

#### Red Flags
- **Too good to be true** - Finding a "complete" eukaryotic pathway in archaea
- **Single gene evidence** - One enzyme doesn't make a pathway
- **Annotation disagreement** - PFAM says X, KEGG says Y
- **Unusual size** - Protein much larger/smaller than expected for the domain

#### Ecological Context Matters
Consider the organism's environment when interpreting findings:

- **Subsurface/anaerobic** → Expect: hydrogenases, CODH, anaerobic respiration. Don't expect: photosynthesis, aerobic pathways
- **Hot springs/thermophiles** → Expect: reverse gyrase, heat-stable enzymes
- **Hypersaline** → Expect: compatible solute synthesis, ion pumps
- **Host-associated** → Expect: reduced genomes, transport-heavy
- **Predatory** → Expect: adhesins, secretion systems, lytic enzymes

#### When to Log a Correction
If you realize a previous finding was wrong or overclaimed:
```python
log_finding(
    category="correction",
    title="CORRECTION: [Original claim] is incorrect",
    description="[Why it's wrong and what the correct interpretation is]",
    evidence={"original_claim": "...", "why_wrong": "...", "correct_interpretation": "..."},
    priority="high"
)
```
Corrections are valuable - they show rigorous science.

#### Metagenome-Specific Pitfalls
MAGs (metagenome-assembled genomes) have unique issues:

- **Contamination** - Bins may contain genes from multiple organisms. If a gene seems out of place (e.g., bacterial gene in archaeal MAG), check contig coverage/GC
- **Incompleteness** - Missing genes may be real absences OR assembly gaps. Don't overclaim "loss of function"
- **Fragmented contigs** - Operons may be split across contigs. A "missing" gene might be on another contig
- **Annotation propagation errors** - Databases contain errors that get propagated. Unusual annotations warrant skepticism
- **Chimeric contigs** - Contigs can join sequences from different organisms. Look for sudden GC/coverage shifts

#### Structure Prediction Caveats
When using ESM3/Foldseek:
- **pLDDT < 0.7** - Low confidence, interpret structure cautiously
- **pTM < 0.5** - Global fold may be unreliable
- **Foldseek E-value** - E < 1e-10 is strong; E > 1e-3 is weak
- **Sequence identity** - Low seq ID + good structural match = possible convergent evolution or distant homology

### Finding Categories

Use these categories when logging findings:
- `prophage` - Prophage/phage-related regions
- `adhesion` - Cell adhesion systems (Dockerin/Cohesin, HYR, etc.)
- `surface_protein` - Surface/secreted proteins
- `defense` - CRISPR, R-M, TA systems
- `metabolism` - Metabolic enzymes and pathways
- `giant_unannotated` - Large proteins with no annotation
- `novel_domain` - Unusual or novel domain architectures
- `mobile_element` - Transposons, integrases, mobile elements
- `unknown_cluster` - Clusters of unannotated genes

### Output Format

At the end of exploration, ensure:
1. All findings are logged to `findings.jsonl`
2. State is saved to `exploration_state.json`
3. Summary is updated in `exploration_summary.md`

Then generate the PDF report and provide a summary to the user:

```python
# Generate PDF report
import subprocess
result = subprocess.run(
    ["python3", "scripts/generate_exploration_report.py"],
    capture_output=True, text=True
)
print(result.stdout)
if result.returncode != 0:
    print(f"Warning: PDF generation failed: {result.stderr}")
```

```
## Exploration Session Complete

**New findings this session:** X
**Total findings:** Y
**Files updated:**
- exploration/findings.jsonl
- exploration/exploration_state.json
- exploration/exploration_summary.md
- exploration/genome_profiles.tsv (if multi-genome)
- exploration/genome_comparison.md (if multi-genome)

**PDF Report:** {DB_DIR}/altiarchaeota_exploration_report.pdf

**Top discoveries:**
1. [Brief description of most important finding]
2. [Second most important]
3. [Third most important]

Run `/explore --resume` to continue exploration.
```

---

## Focus Instructions

### Default (no focus)
Complete all three phases:
1. Quick predicate searches across categories (giant, defense, metabolism, etc.)
2. **Browse at least 15-20 contigs systematically** - this is required, not optional
3. Deep dive on top 5 most interesting findings

Aim for breadth in Phase 1, depth in Phase 2.

### --focus metabolism
Focus on metabolic systems: hydrogenases, oxidoreductases, electron carriers, unusual cofactors.

### --focus phage
Focus on prophage regions: phage structural genes, integrases, uncharacterized phage proteins.

### --focus defense
Focus on defense systems: CRISPR arrays, R-M systems, TA modules, abortive infection.

### --focus novel
Focus on uncharacterized proteins: giant unannotated, DUF-only proteins, hypotheticals in operons.

### --focus adhesion
Focus on cell adhesion/surface systems: Dockerin/Cohesin, HYR domains, PKD domains, surface layer proteins.

### --resume
Load state and continue from where you left off. Check what's been searched and explore new areas.

### --contig <id>
Browse a specific contig window-by-window. Log interesting regions as you go.

---

## Genome Profiling Workflow (Multi-Genome Datasets)

**Required for any dataset with >1 genome. Run this in Phase 0 before other exploration.**

This produces:
1. A profile for each genome (functional capabilities, notable features)
2. Comparative matrices showing feature distribution
3. Identification of genome clusters and outliers

**Step 1: Generate per-genome statistics**

```python
# Get list of all genomes
genomes = b.store.execute("""
    SELECT bin_id, COUNT(*) as n_proteins
    FROM proteins
    GROUP BY bin_id
    ORDER BY n_proteins DESC
""")
print(f"Analyzing {len(genomes)} genomes")

# For each genome, compute a functional profile
genome_profiles = {}

for bin_id, n_proteins in genomes:
    if bin_id in state.get('genomes_explored', []):
        continue  # Skip already profiled

    profile = {
        'bin_id': bin_id,
        'n_proteins': n_proteins,
    }

    # Basic stats
    stats = b.store.execute("""
        SELECT
            AVG(sequence_length) as avg_length,
            MAX(sequence_length) as max_length,
            AVG(gc_content) as avg_gc
        FROM proteins WHERE bin_id = ?
    """, [bin_id])[0]
    profile['avg_length'] = round(stats[0], 1)
    profile['max_length'] = stats[1]
    profile['avg_gc'] = round(stats[2], 3) if stats[2] else None

    # Annotation rate
    annot = b.store.execute("""
        SELECT
            COUNT(DISTINCT p.protein_id) as total,
            COUNT(DISTINCT a.protein_id) as annotated
        FROM proteins p
        LEFT JOIN annotations a ON p.protein_id = a.protein_id
        WHERE p.bin_id = ?
    """, [bin_id])[0]
    profile['annotated'] = annot[1]
    profile['annotation_rate'] = round(annot[1] / annot[0] * 100, 1) if annot[0] > 0 else 0

    genome_profiles[bin_id] = profile
```

**Step 2: Count key features per genome using predicates**

```python
# Key predicates to profile (functional categories)
KEY_PREDICATES = [
    # Energy metabolism
    'hydrogenase', 'electron_transport', 'atp_synthesis',
    # Carbon metabolism
    'carbon_fixation', 'methanogenesis', 'wood_ljungdahl',
    # Defense systems
    'crispr_associated', 'restriction_modification', 'toxin_antitoxin', 'defense_system',
    # Mobile elements
    'transposase', 'integrase', 'phage_related', 'mobile_element',
    # Surface biology
    'adhesin', 's_layer', 'secretion_system', 'flagellum',
    # Annotation status
    'unannotated', 'hypothetical', 'multi_domain',
    # Size classes
    'giant', 'massive',
]

for bin_id, profile in genome_profiles.items():
    # Get all predicates for this genome's proteins
    pred_counts = b.store.execute("""
        SELECT unnest(pp.predicates) as pred, COUNT(*) as cnt
        FROM proteins p
        JOIN protein_predicates pp ON p.protein_id = pp.protein_id
        WHERE p.bin_id = ?
        GROUP BY pred
    """, [bin_id])

    pred_dict = {row[0]: row[1] for row in pred_counts}

    for pred in KEY_PREDICATES:
        profile[pred] = pred_dict.get(pred, 0)
```

**Step 3: Build comparative matrix**

```python
import csv
from pathlib import Path

# Output genome profiles to TSV for comparative analysis
output_file = EXPLORE_DIR / "genome_profiles.tsv"

# Define column order
columns = ['bin_id', 'n_proteins', 'avg_length', 'max_length', 'avg_gc',
           'annotated', 'annotation_rate'] + KEY_PREDICATES

with open(output_file, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=columns, delimiter='\t', extrasaction='ignore')
    writer.writeheader()
    for profile in genome_profiles.values():
        writer.writerow(profile)

print(f"Wrote {len(genome_profiles)} genome profiles to {output_file}")
```

**Step 4: Identify patterns and outliers**

```python
# Find genomes with unusual features
findings = []

for bin_id, profile in genome_profiles.items():
    # High defense burden
    defense_total = profile.get('crispr_associated', 0) + profile.get('restriction_modification', 0) + profile.get('toxin_antitoxin', 0)
    if defense_total > 30:
        findings.append(f"{bin_id}: High defense burden ({defense_total} defense genes)")

    # Low annotation rate (novel biology?)
    if profile.get('annotation_rate', 100) < 50:
        findings.append(f"{bin_id}: Low annotation rate ({profile['annotation_rate']}%) - novel biology?")

    # No hydrogenases (unusual for this lineage)
    if profile.get('hydrogenase', 0) == 0:
        findings.append(f"{bin_id}: No hydrogenases detected - check lifestyle")

    # High mobile element load
    mobile = profile.get('transposase', 0) + profile.get('integrase', 0)
    if mobile > 20:
        findings.append(f"{bin_id}: High mobile element load ({mobile} elements)")

    # Missing expected features
    if profile.get('flagellum', 0) == 0:
        findings.append(f"{bin_id}: No flagellum/archaellum detected")

# Log outliers
for finding in findings:
    print(finding)
```

**Step 5: Generate comparative summary**

```python
# Compute feature distribution statistics
summary_lines = ["# Genome Comparative Analysis\n"]
summary_lines.append(f"**Genomes analyzed:** {len(genome_profiles)}\n")

# For each key predicate, show distribution
summary_lines.append("\n## Feature Distribution Across Genomes\n")
summary_lines.append("| Feature | Genomes with | Mean count | Max count | Min (non-zero) |")
summary_lines.append("|---------|-------------|------------|-----------|----------------|")

for pred in KEY_PREDICATES:
    values = [p.get(pred, 0) for p in genome_profiles.values()]
    n_with = sum(1 for v in values if v > 0)
    non_zero = [v for v in values if v > 0]
    mean_val = sum(values) / len(values) if values else 0
    max_val = max(values) if values else 0
    min_nz = min(non_zero) if non_zero else 0

    summary_lines.append(f"| {pred} | {n_with}/{len(values)} ({round(n_with/len(values)*100)}%) | {mean_val:.1f} | {max_val} | {min_nz} |")

# Genome clusters by capability
summary_lines.append("\n## Genome Groupings\n")

# Group by presence/absence of key systems
has_hydrogenase = [b for b, p in genome_profiles.items() if p.get('hydrogenase', 0) > 0]
has_crispr = [b for b, p in genome_profiles.items() if p.get('crispr_associated', 0) > 0]
has_adhesin = [b for b, p in genome_profiles.items() if p.get('adhesin', 0) > 0]

summary_lines.append(f"- **Hydrogenase-positive:** {len(has_hydrogenase)} genomes")
summary_lines.append(f"- **CRISPR-positive:** {len(has_crispr)} genomes")
summary_lines.append(f"- **Adhesin-positive:** {len(has_adhesin)} genomes")

# Write summary
summary_file = EXPLORE_DIR / "genome_comparison.md"
with open(summary_file, 'w') as f:
    f.write('\n'.join(summary_lines))

print(f"Wrote comparative summary to {summary_file}")
```

**Step 6: Update state and log completion**

```python
# Mark all genomes as explored
state['genomes_explored'] = list(genome_profiles.keys())
save_state(state)

# Log the analysis as a finding
log_finding(
    category="comparative_analysis",
    title=f"Genome-by-genome analysis of {len(genome_profiles)} genomes",
    description=f"Generated functional profiles and comparative matrix. See genome_profiles.tsv and genome_comparison.md",
    evidence={
        "n_genomes": len(genome_profiles),
        "output_files": ["genome_profiles.tsv", "genome_comparison.md"]
    },
    priority="high"
)

update_summary()
```

#### Output Files

The `--genomes` mode produces:
- `exploration/genome_profiles.tsv` - Tab-separated table with one row per genome, columns for each feature count
- `exploration/genome_comparison.md` - Human-readable summary with distribution statistics and genome groupings

#### What to Look For

When reviewing genome profiles:

1. **Outlier genomes** - Unusually high/low counts for key features may indicate:
   - Different ecological niche
   - Contamination or chimeric assembly
   - Lineage-specific adaptations

2. **Correlated features** - Features that co-occur across genomes suggest:
   - Functional linkage
   - Shared evolutionary history
   - Ecological constraints

3. **Variable vs core features** - Features present in:
   - >90% of genomes = Core, essential
   - 50-90% = Variable, possibly niche-specific
   - <50% = Accessory, may be recent acquisitions or losses
   - 1-2 genomes = Suspicious, check for contamination

4. **Missing expected features** - If a feature expected for the lineage is absent:
   - Assembly gaps or incompleteness?
   - Genuine loss?
   - Alternative pathway?
