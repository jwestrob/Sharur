# Defense Systems Skill

Analyze defense systems: CRISPR-Cas, restriction-modification, toxin-antitoxin, and other anti-phage mechanisms.

**CONCURRENCY: DuckDB does not support concurrent writes. Only ONE agent should access a database at a time. The coordinator must run DB-accessing skills sequentially, not in parallel.**

> **Mandatory:** Read `docs/findings_spec.md` for the structured findings format.
> Write `findings.jsonl` alongside your prose reports. Every quantitative claim needs a verification query.

> **Mandatory:** Follow the shared validation protocols in `_validation_protocols.md`.
> Verify accession names before reporting. Use COUNT(DISTINCT protein_id) for protein
> counts. Apply Context-First protocol for annotations averaging >10 hits/genome.

> **Literature dispatch:** When you encounter ambiguous annotations, unknown Foldseek hits,
> or need to make comparative claims ("first known", "largest"), dispatch a literature agent.
> Read `.claude/skills/literature.md` for protocols.

---

## Usage

```
/defense                    # Full defense inventory
/defense --crispr           # CRISPR-Cas focus
/defense --genome GENOME_ID # Single genome analysis
```

---

## Prompt

You are analyzing defense systems in a metagenomic dataset. Provide a comprehensive inventory with proper classification and validation.

### Step 1: Run Built-in CRISPR Analysis

```python
from sharur.operators import Sharur
b = Sharur("data/DATASET/sharur.duckdb")

# Comprehensive CRISPR-Cas analysis
print(b.analyze_crispr_systems())
```

This returns:
- Total arrays and their types
- Complete vs incomplete systems
- Orphan arrays (no nearby cas genes)
- Arrays at contig edges (fragmentation artifacts)

### Step 2: Defense System Inventory

```python
# Count all defense-related predicates
defense_predicates = [
    'crispr_associated', 'cas_domain',
    'restriction_modification', 'restriction_enzyme', 'methyltransferase',
    'toxin_antitoxin', 'toxin_domain', 'antitoxin_domain',
    'defense_system', 'anti_crispr', 'anti_restriction',
    'abortive_infection', 'retron',
]

for pred in defense_predicates:
    count = b.store.execute(f"""
        SELECT COUNT(*) FROM protein_predicates
        WHERE '{pred}' = ANY(predicates)
    """)[0][0]
    if count > 0:
        print(f"{pred}: {count}")
```

### Step 2b: System-Level Defense Validation (CRITICAL)

**Astra's DefenseFinder HMM hits have a ~72% false positive rate.** Most hits are superfamily
matches (kinases, Rossmann folds, helicases) that are NOT part of real defense systems.
Always prefer system-validated annotations (`source='defensefinder_system'`) over raw HMM
hits (`source='defensefinder'`).

```python
# PREFERRED: System-validated defense annotations (co-localization checked)
# These come from running defense-finder proper (MacSyFinder system detection)
# Run: python scripts/run_defensefinder_systems.py data/DATASET/
validated = b.store.execute("""
    SELECT a.name, COUNT(DISTINCT a.protein_id) as n_proteins,
           COUNT(DISTINCT p.bin_id) as n_genomes
    FROM annotations a
    JOIN proteins p ON a.protein_id = p.protein_id
    WHERE a.source = 'defensefinder_system'
    GROUP BY a.name
    ORDER BY n_proteins DESC
""")
print("System-validated defense proteins:")
for name, n_prot, n_gen in validated:
    print(f"  {name}: {n_prot} proteins in {n_gen} genomes")

# System-level summary from defense_systems table
systems = b.store.execute("""
    SELECT system_type, system_subtype, COUNT(*) as n_systems,
           SUM(genes_count) as total_genes
    FROM defense_systems
    GROUP BY system_type, system_subtype
    ORDER BY n_systems DESC
""")
print("\nDefense systems by type:")
for stype, subtype, n_sys, n_genes in systems:
    print(f"  {stype}/{subtype}: {n_sys} systems, {n_genes} genes")
```

**If `defensefinder_system` annotations don't exist yet**, run the system validation:
```bash
python scripts/run_defensefinder_systems.py data/DATASET/ --workers 8
```
This takes ~2 minutes per 50k proteins. Uses `--db-type ordered_replicon` for co-localization.

**Fallback: raw HMM hits** (use with extreme caution, expect ~72% FPs):
```python
# Raw DefenseFinder HMM annotations (NOT system-validated)
df_systems = b.store.execute("""
    SELECT a.accession, a.name, a.description,
           COUNT(DISTINCT a.protein_id) as n_proteins,
           COUNT(DISTINCT p.bin_id) as n_genomes
    FROM annotations a
    JOIN proteins p ON a.protein_id = p.protein_id
    WHERE a.source = 'defensefinder'
    GROUP BY a.accession, a.name, a.description
    ORDER BY n_proteins DESC
""")

if df_systems:
    print("DefenseFinder HMM Hits (WARNING: ~72% FP rate):")
    for acc, name, desc, n_prot, n_gen in df_systems:
        print(f"  {name} ({acc}): {n_prot} proteins in {n_gen} genomes")
else:
    print("No DefenseFinder annotations found")

# PADLOC annotations (same pattern)
padloc_systems = b.store.execute("""
    SELECT a.accession, a.name, a.description,
           COUNT(DISTINCT a.protein_id) as n_proteins,
           COUNT(DISTINCT p.bin_id) as n_genomes
    FROM annotations a
    JOIN proteins p ON a.protein_id = p.protein_id
    WHERE a.source = 'padloc'
    GROUP BY a.accession, a.name, a.description
    ORDER BY n_proteins DESC
""")

if padloc_systems:
    print("\nPADLOC Systems:")
    for acc, name, desc, n_prot, n_gen in padloc_systems:
        print(f"  {name} ({acc}): {n_prot} proteins in {n_gen}/{n_genomes} genomes")
```

### Step 2c: Cross-Reference Defense Detections

Compare DefenseFinder/PADLOC system-level calls with predicate-based detection to find discrepancies:

```python
# Predicate-based defense protein count
pred_defense = b.store.execute("""
    SELECT COUNT(DISTINCT protein_id) FROM protein_predicates
    WHERE 'defense_system' = ANY(predicates)
""")[0][0]

# DefenseFinder protein count
df_defense = b.store.execute("""
    SELECT COUNT(DISTINCT protein_id) FROM annotations
    WHERE source = 'defensefinder'
""")[0][0]

print(f"Defense proteins (predicates): {pred_defense}")
print(f"Defense proteins (DefenseFinder): {df_defense}")

# Find DefenseFinder proteins NOT in predicate defense set
if df_defense > 0:
    missing_from_preds = b.store.execute("""
        SELECT a.protein_id, a.name, a.accession
        FROM annotations a
        WHERE a.source = 'defensefinder'
          AND a.protein_id NOT IN (
            SELECT protein_id FROM protein_predicates
            WHERE 'defense_system' = ANY(predicates)
          )
        LIMIT 10
    """)
    if missing_from_preds:
        print(f"\nDefenseFinder hits NOT in defense predicates ({len(missing_from_preds)} shown):")
        for pid, name, acc in missing_from_preds:
            print(f"  {pid}: {name} ({acc})")
```

### Step 3: CRISPR-Cas Classification

**Type I** (signature: Cas3 helicase-nuclease):
```python
# Find Cas3 proteins
cas3 = b.search_by_predicates(has=["cas3"])
print(f"Type I systems (Cas3): {len(cas3.data)} proteins")

# Check subtypes by Cas8 variants
for protein_id in cas3.data[:5]:
    print(b.get_neighborhood(protein_id, window=10))
```

**Type III** (signature: Cas10/Csm1):
```python
cas10 = b.store.execute("""
    SELECT p.protein_id, p.bin_id
    FROM proteins p
    JOIN annotations a ON p.protein_id = a.protein_id
    WHERE a.name ILIKE '%cas10%' OR a.name ILIKE '%csm1%' OR a.name ILIKE '%cmr2%'
""")
print(f"Type III systems (Cas10): {len(cas10)} proteins")
```

**Type V** (signature: Cas12):

For **Cas12a/b/c/d/e** (the larger, well-characterized effectors): use the standard `b.get_neighborhood(...)` validation pattern from `_validation_protocols.md` §3 (Context-First Protocol). Real Cas12a-e effectors are 1000-1500 aa, near CRISPR arrays, with adaptation modules.

For **Cas12f / Cas12f1 — DO NOT REPORT as a Cas finding by default.** The PFAM `Cas12f*` model cross-matches TnpB and IS200/IS605 OrfB at high rates. See `_validation_protocols.md` §8 for the four conditions that must ALL hold before elevating a Cas12f hit to a finding. If any condition fails, treat as TnpB-family RuvC homolog, inventory at PFAM-completeness level only, and do not open an explore agent on it.

### Step 4: Restriction-Modification Systems

```python
# R-M system components
rm_proteins = b.store.execute("""
    SELECT p.protein_id, p.bin_id, a.name
    FROM proteins p
    JOIN annotations a ON p.protein_id = a.protein_id
    WHERE a.name ILIKE '%restriction%'
       OR a.name ILIKE '%methyltransferase%'
       OR a.name ILIKE '%methylase%'
    ORDER BY p.bin_id, p.gene_index
""")

# Group by genome and check for paired R-M
from collections import defaultdict
by_genome = defaultdict(lambda: {'restriction': [], 'methyltransferase': []})
for pid, genome, name in rm_proteins:
    if 'restriction' in name.lower():
        by_genome[genome]['restriction'].append(pid)
    else:
        by_genome[genome]['methyltransferase'].append(pid)

for genome, counts in by_genome.items():
    print(f"{genome}: {len(counts['restriction'])} restriction, {len(counts['methyltransferase'])} MTases")
```

### Step 5: Toxin-Antitoxin Systems

```python
# TA system analysis
ta_proteins = b.search_by_predicates(has=["toxin_antitoxin"])
print(f"Confirmed TA pairs: {len(ta_proteins.data)}")

# Check for orphan toxins (dangerous!)
orphan_toxins = b.store.execute("""
    WITH toxins AS (
        SELECT p.protein_id, p.contig_id, p.gene_index
        FROM proteins p
        JOIN protein_predicates pp ON p.protein_id = pp.protein_id
        WHERE 'toxin_domain' = ANY(pp.predicates)
    ),
    antitoxins AS (
        SELECT p.contig_id, p.gene_index
        FROM proteins p
        JOIN protein_predicates pp ON p.protein_id = pp.protein_id
        WHERE 'antitoxin_domain' = ANY(pp.predicates)
    )
    SELECT t.protein_id
    FROM toxins t
    LEFT JOIN antitoxins a ON t.contig_id = a.contig_id
        AND ABS(t.gene_index - a.gene_index) <= 3
    WHERE a.gene_index IS NULL
""")
print(f"Orphan toxins (no nearby antitoxin): {len(orphan_toxins)}")
```

### Step 6: Novel/Rare Defense Systems

```python
# Check for recently discovered systems
novel_defense = [
    ('retron', 'Retron defense'),
    ('abortive_infection', 'Abortive infection'),
    ('thoeris', 'Thoeris'),
    ('gabija', 'Gabija'),
    ('lamassu', 'Lamassu'),
    ('septu', 'Septu'),
    ('hachiman', 'Hachiman'),
]

for pred, name in novel_defense:
    count = b.store.execute(f"""
        SELECT COUNT(*) FROM protein_predicates
        WHERE '{pred}' = ANY(predicates)
    """)[0][0]
    if count > 0:
        print(f"{name}: {count} proteins")
```

### Step 7: Defense Islands

Defense systems often cluster. Look for defense islands:

```python
# Find contigs with multiple defense genes
defense_rich = b.store.execute("""
    SELECT p.contig_id, COUNT(DISTINCT p.protein_id) as n_defense
    FROM proteins p
    JOIN protein_predicates pp ON p.protein_id = pp.protein_id
    WHERE 'crispr_associated' = ANY(pp.predicates)
       OR 'restriction_modification' = ANY(pp.predicates)
       OR 'toxin_antitoxin' = ANY(pp.predicates)
       OR 'defense_system' = ANY(pp.predicates)
    GROUP BY p.contig_id
    HAVING COUNT(DISTINCT p.protein_id) >= 5
    ORDER BY n_defense DESC
""")

print(f"Defense islands (≥5 defense genes): {len(defense_rich)}")
for contig, n in defense_rich[:5]:
    print(f"  {contig}: {n} defense genes")
```

### Step 8: Visualize Key Loci

> **Read `.claude/skills/visualize.md` before generating figures.** For any figure going into
> a report, use `plot_locus_multisource.py` (multi-source annotations, publication styling).
> Use `b.visualize_neighborhood()` for quick exploratory checks only. Note: Mokosh Type I
> and BREX pglW are serine/threonine kinase superfamily false positives — exclude from
> defense heatmaps (see survey-014).

```python
from pathlib import Path

FIGURES_DIR = Path("data/DATASET/exploration/figures")
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# For publication figures, prefer the multi-source CLI script:
# python scripts/plot_locus_multisource.py --db data/DATASET/sharur.duckdb \
#     --protein PROTEIN_ID --window 15 --output figures/crispr_locus.png

# For quick exploration:
if cas3.data:
    b.visualize_neighborhood(
        cas3.data[0],
        window=15,
        output_path=str(FIGURES_DIR / "crispr_type_i_locus.png")
    )
```

### Output Format

```markdown
## Defense Systems Inventory

### Summary
- CRISPR-Cas systems: X complete, Y incomplete, Z orphan arrays
- R-M systems: X complete pairs
- TA systems: X confirmed pairs, Y orphan toxins
- Other defense: [list]

### CRISPR-Cas Details
[Per-type breakdown with validation notes]

### Defense Islands
[List defense-rich contigs with gene inventories]

### Figures
- crispr_type_i_locus.png: Representative Type I system
- defense_island_contig_X.png: Dense defense region

### Biological Interpretation
[What does the defense repertoire suggest about viral pressure?]
```

### Step 9: Ecological Interpretation

Go beyond inventories — interpret the defense repertoire ecologically.

#### Defense Load Calculation

```python
# Calculate defense investment per genome
defense_load = b.store.execute("""
    SELECT p.bin_id,
           COUNT(DISTINCT CASE WHEN 'defense_system' = ANY(pp.predicates)
                               OR 'crispr_associated' = ANY(pp.predicates)
                               OR 'restriction_modification' = ANY(pp.predicates)
                               OR 'toxin_antitoxin' = ANY(pp.predicates)
                 THEN p.protein_id END) as n_defense,
           COUNT(DISTINCT p.protein_id) as n_total
    FROM proteins p
    LEFT JOIN protein_predicates pp ON p.protein_id = pp.protein_id
    GROUP BY p.bin_id
    ORDER BY n_defense DESC
""")

for genome, n_defense, n_total in defense_load:
    pct = n_defense / n_total * 100 if n_total > 0 else 0
    print(f"  {genome}: {n_defense}/{n_total} ({pct:.1f}%)")
```

#### Interpretation Guide

| Defense Load | Interpretation |
|-------------|----------------|
| >10% of proteins | Heavy phage pressure — active arms race |
| 5-10% | Moderate — typical for free-living prokaryotes |
| <5% | Low — stable niche, reduced viral pressure, or obligate lifestyle |
| 0% | Check assembly quality before concluding |

#### Defense Diversity vs Depth

- **Many system types, few copies each** → diverse viral threats, broad defense
- **Few system types, many copies** → specialized defense against specific phage families
- **CRISPR-dominant** → adaptive immunity primary strategy, high spacer diversity expected
- **Innate-dominant** → constitutive defense, less phage-responsive

#### Co-occurrence Patterns

```python
# Which defense systems co-occur in the same genomes?
# Build a genome x defense-type matrix
defense_types = ['crispr_associated', 'restriction_modification',
                 'toxin_antitoxin', 'abortive_infection']

for i, dt1 in enumerate(defense_types):
    for dt2 in defense_types[i+1:]:
        both = b.store.execute(f"""
            SELECT COUNT(DISTINCT p.bin_id) FROM proteins p
            JOIN protein_predicates pp ON p.protein_id = pp.protein_id
            WHERE '{dt1}' = ANY(pp.predicates)
              AND p.bin_id IN (
                SELECT DISTINCT p2.bin_id FROM proteins p2
                JOIN protein_predicates pp2 ON p2.protein_id = pp2.protein_id
                WHERE '{dt2}' = ANY(pp2.predicates)
              )
        """)[0][0]
        print(f"  {dt1} + {dt2}: {both} genomes")
```

---

## Defense Island Detection

**Goal:** Cluster DefenseFinder-annotated genes into genomic loci ("defense islands") and write them to the `loci`/`locus_proteins` tables for downstream analysis and prophage co-localization.

Defense systems cluster in bacterial genomes — often near prophage insertion sites, forming "defense hotspots." Detecting these islands turns point annotations into spatial loci.

### Step 10: Tag Genes for Defense Island Detection

```python
import json, uuid, re
from collections import defaultdict

# Load all proteins with gene positions
proteins = b.store.execute("""
    SELECT protein_id, contig_id, gene_index, start, end_coord, strand, bin_id
    FROM proteins
    WHERE gene_index IS NOT NULL
    ORDER BY contig_id, gene_index
""")

# Build contig -> gene list
contig_genes = defaultdict(list)
for pid, cid, gi, s, e, strand, bid in proteins:
    contig_genes[cid].append({
        "protein_id": pid, "contig_id": cid, "gene_index": gi,
        "start": s, "end_coord": e, "strand": strand, "bin_id": bid,
    })

# Load DefenseFinder annotations
df_annots = b.store.execute("""
    SELECT protein_id, accession, name, description
    FROM annotations
    WHERE source = 'defensefinder'
""")

# Index by protein_id
protein_defense = defaultdict(list)
for pid, acc, name, desc in df_annots:
    protein_defense[pid].append((acc, name or '', desc or ''))

print(f"DefenseFinder annotations: {len(df_annots):,} on {len(protein_defense):,} proteins")

# Tag each gene
def tag_gene_defense(pid, annots):
    """Tag: 'marker' if has DefenseFinder hit, else 'negative'."""
    extras = {"marker_accessions": [], "system_types": set()}

    for acc, name, desc in annots:
        extras["marker_accessions"].append(acc)
        # Extract system type from accession (e.g., "CBASS__CdnC" → "CBASS")
        if "__" in acc:
            extras["system_types"].add(acc.split("__")[0])
        elif name:
            extras["system_types"].add(name)

    if extras["marker_accessions"]:
        return "marker", extras
    return "negative", extras

for cid, genes in contig_genes.items():
    genes.sort(key=lambda g: g["gene_index"])
    for gene in genes:
        annots = protein_defense.get(gene["protein_id"], [])
        tag, extras = tag_gene_defense(gene["protein_id"], annots)
        gene["tag"] = tag
        gene["extras"] = extras
```

### Step 11: Gap-Based Island Clustering

Same algorithm as prophage detection, tuned for defense islands.

```python
# Defense island parameters
MAX_GAP = 5        # Max gap between defense genes
MIN_MARKERS = 3    # Min DefenseFinder-annotated genes
MIN_GENES = 3      # Min total genes in island (markers only — no "support" concept here)

def detect_defense_islands(genes, max_gap=MAX_GAP, min_markers=MIN_MARKERS, min_genes=MIN_GENES):
    """Gap-based clustering of DefenseFinder markers on one contig."""
    islands = []
    current = []
    gap_count = 0

    def emit(cluster):
        marker_count = sum(1 for g in cluster if g["tag"] == "marker")
        if marker_count >= min_markers and len(cluster) >= min_genes:
            islands.append({
                "genes": list(cluster),
                "marker_count": marker_count,
            })

    for gene in genes:
        if gene["tag"] == "marker":
            current.append(gene)
            gap_count = 0
        else:
            if current:
                gap_count += 1
                if gap_count > max_gap:
                    while current and current[-1]["tag"] == "negative":
                        current.pop()
                    emit(current)
                    current = []
                    gap_count = 0
                else:
                    current.append(gene)

    if current:
        while current and current[-1]["tag"] == "negative":
            current.pop()
        emit(current)

    return islands

# Run across all contigs
all_islands = []
for cid, genes in contig_genes.items():
    islands = detect_defense_islands(genes)
    for island in islands:
        island["contig_id"] = cid
        island["bin_id"] = genes[0]["bin_id"]
    all_islands.extend(islands)

print(f"Defense islands detected: {len(all_islands)}")
```

### Step 12: Score and Classify Defense Islands

Islands with multiple distinct defense system types are the most biologically significant — they represent defense hotspots where systems accumulate, likely near prophage insertion sites.

```python
def score_defense_island(island):
    """Score a defense island by system diversity and size.

    Returns: (confidence, level, system_types)
    """
    system_types = set()
    for g in island["genes"]:
        system_types.update(g["extras"].get("system_types", set()))

    marker_count = island["marker_count"]
    n_systems = len(system_types)

    if n_systems >= 3 and marker_count >= 5:
        return 0.95, "hotspot", system_types      # Multi-system defense hotspot
    if n_systems >= 2 and marker_count >= 5:
        return 0.90, "high", system_types          # Two+ systems, well-supported
    if n_systems >= 2 and marker_count >= 3:
        return 0.80, "high", system_types          # Two systems
    if marker_count >= 5:
        return 0.75, "medium", system_types        # Single system, many components
    if marker_count >= 3:
        return 0.60, "medium", system_types        # Minimum island
    return 0.40, "low", system_types

for island in all_islands:
    confidence, level, system_types = score_defense_island(island)
    island["confidence"] = confidence
    island["confidence_level"] = level
    island["system_types"] = system_types

# Summary
hotspots = [i for i in all_islands if i["confidence_level"] == "hotspot"]
high = [i for i in all_islands if i["confidence_level"] == "high"]
medium = [i for i in all_islands if i["confidence_level"] == "medium"]

print(f"Defense hotspots (3+ system types): {len(hotspots)}")
print(f"High confidence (2+ system types): {len(high)}")
print(f"Medium confidence: {len(medium)}")
```

### Step 13: Write Defense Islands to Loci Tables

```python
def write_defense_loci(islands, locus_type="defense_island", clear_existing=True):
    """Write defense islands to loci/locus_proteins tables."""

    if clear_existing:
        existing = b.store.execute(
            "SELECT COUNT(*) FROM loci WHERE locus_type = ?", [locus_type]
        )[0][0]
        if existing > 0:
            print(f"Clearing {existing} existing {locus_type} loci")
            b.store.execute("""
                DELETE FROM locus_proteins
                WHERE locus_id IN (SELECT locus_id FROM loci WHERE locus_type = ?)
            """, [locus_type])
            b.store.execute("DELETE FROM loci WHERE locus_type = ?", [locus_type])

    n_links = 0
    for island in islands:
        genes = island["genes"]
        locus_id = f"{locus_type}_{uuid.uuid4().hex[:12]}"

        start = min(g["start"] for g in genes)
        end_coord = max(g["end_coord"] for g in genes)

        metadata = {
            "marker_count": island["marker_count"],
            "total_genes": len(genes),
            "system_types": sorted(island.get("system_types", set())),
            "n_system_types": len(island.get("system_types", set())),
            "confidence_level": island["confidence_level"],
            "marker_accessions": sorted(set(
                a for g in genes for a in g["extras"].get("marker_accessions", [])
            )),
        }

        b.store.execute("""
            INSERT INTO loci (locus_id, locus_type, contig_id, start, end_coord, confidence, metadata)
            VALUES (?, ?, ?, ?, ?, ?, ?)
        """, [locus_id, locus_type, island["contig_id"], start, end_coord,
              island["confidence"], json.dumps(metadata)])

        for pos, gene in enumerate(genes):
            b.store.execute("""
                INSERT INTO locus_proteins (locus_id, protein_id, position)
                VALUES (?, ?, ?)
            """, [locus_id, gene["protein_id"], pos])
            n_links += 1

    b.store.commit()
    print(f"Wrote {len(islands)} {locus_type} loci with {n_links:,} protein links")

write_defense_loci(all_islands, locus_type="defense_island")
```

**IMPORTANT:** Use `locus_type='defense_island'` (not `'island'`). The Omnitrophota dataset used `'island'` historically but `'defense_island'` is the correct convention going forward. Do NOT touch existing loci of other types (crispr, prophage, viral_contig).

### Step 14: Characterize Defense Islands

```python
def characterize_defense_islands():
    """Analyze defense island composition and distribution."""

    loci = b.store.execute("""
        SELECT l.locus_id, l.contig_id, l.confidence, l.metadata,
               COUNT(lp.protein_id) as n_genes
        FROM loci l
        JOIN locus_proteins lp ON l.locus_id = lp.locus_id
        WHERE l.locus_type = 'defense_island'
        GROUP BY l.locus_id, l.contig_id, l.confidence, l.metadata
    """)

    print(f"\n=== DEFENSE ISLAND CHARACTERIZATION ===")
    print(f"Total defense islands: {len(loci)}")

    # Size distribution
    sizes = [n for _, _, _, _, n in loci]
    if sizes:
        print(f"Size (genes): min={min(sizes)}, median={sorted(sizes)[len(sizes)//2]}, max={max(sizes)}")

    # System type composition across all islands
    system_counts = defaultdict(int)
    multi_system = 0
    for lid, cid, conf, meta_json, n_genes in loci:
        meta = json.loads(meta_json) if isinstance(meta_json, str) else meta_json
        sys_types = meta.get("system_types", [])
        if len(sys_types) >= 2:
            multi_system += 1
        for st in sys_types:
            system_counts[st] += 1

    print(f"\nMulti-system islands: {multi_system}/{len(loci)} ({multi_system/len(loci)*100:.0f}%)")
    print("\nSystem types across defense islands:")
    for st, n in sorted(system_counts.items(), key=lambda x: -x[1])[:15]:
        print(f"  {st}: {n} islands")

    # Genome distribution
    genome_dist = b.store.execute("""
        SELECT p.bin_id, COUNT(DISTINCT l.locus_id) as n_islands
        FROM loci l
        JOIN locus_proteins lp ON l.locus_id = lp.locus_id
        JOIN proteins p ON lp.protein_id = p.protein_id
        WHERE l.locus_type = 'defense_island'
        GROUP BY p.bin_id
        ORDER BY n_islands DESC
    """)

    n_genomes = b.store.execute("SELECT COUNT(DISTINCT bin_id) FROM proteins")[0][0]
    print(f"\nGenomes with defense islands: {len(genome_dist)}/{n_genomes} ({len(genome_dist)/n_genomes*100:.1f}%)")

    return loci, system_counts, genome_dist

island_loci, system_counts, genome_dist = characterize_defense_islands()
```

### Step 15: Prophage Co-localization Analysis

Defense systems often accumulate near prophage insertion sites. If prophage loci exist, check for spatial association.

```python
# Check if prophage loci exist
prophage_count = b.store.execute("""
    SELECT COUNT(*) FROM loci WHERE locus_type = 'prophage'
""")[0][0]

if prophage_count > 0:
    # Find defense islands within 20 genes of a prophage
    colocalized = b.store.execute("""
        WITH defense_bounds AS (
            SELECT l.locus_id as defense_id, p.contig_id, p.bin_id,
                   MIN(p.gene_index) as d_start, MAX(p.gene_index) as d_end
            FROM loci l
            JOIN locus_proteins lp ON l.locus_id = lp.locus_id
            JOIN proteins p ON lp.protein_id = p.protein_id
            WHERE l.locus_type = 'defense_island'
            GROUP BY l.locus_id, p.contig_id, p.bin_id
        ),
        prophage_bounds AS (
            SELECT l.locus_id as prophage_id, p.contig_id,
                   MIN(p.gene_index) as p_start, MAX(p.gene_index) as p_end
            FROM loci l
            JOIN locus_proteins lp ON l.locus_id = lp.locus_id
            JOIN proteins p ON lp.protein_id = p.protein_id
            WHERE l.locus_type = 'prophage'
            GROUP BY l.locus_id, p.contig_id
        )
        SELECT db.defense_id, pb.prophage_id, db.bin_id, db.contig_id
        FROM defense_bounds db
        JOIN prophage_bounds pb ON db.contig_id = pb.contig_id
            AND (ABS(db.d_start - pb.p_end) <= 20 OR ABS(pb.p_start - db.d_end) <= 20)
    """)

    n_defense_near_prophage = len(set(r[0] for r in colocalized))
    n_prophage_near_defense = len(set(r[1] for r in colocalized))
    total_defense = len(island_loci)

    print(f"\nProphage-defense co-localization:")
    print(f"  Defense islands near prophages: {n_defense_near_prophage}/{total_defense} "
          f"({n_defense_near_prophage/total_defense*100:.0f}%)")
    print(f"  Prophages near defense islands: {n_prophage_near_defense}/{prophage_count}")
else:
    print("\nNo prophage loci detected — run /prophage first for co-localization analysis")
```

### Step 16: Visualize Representative Defense Islands

```python
from pathlib import Path

FIGURES_DIR = Path("data/DATASET/exploration/figures")
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Visualize:
# 1. A defense hotspot (3+ system types)
# 2. A large multi-system island
# 3. A defense island adjacent to a prophage (if co-localized)

def visualize_defense_island(locus_id, label, description):
    """Visualize a defense island using the center protein."""
    center = b.store.execute("""
        SELECT lp.protein_id FROM locus_proteins lp
        JOIN proteins p ON lp.protein_id = p.protein_id
        WHERE lp.locus_id = ?
        ORDER BY lp.position
        LIMIT 1 OFFSET (SELECT COUNT(*)/2 FROM locus_proteins WHERE locus_id = ?)
    """, [locus_id, locus_id])

    if center:
        locus_size = b.store.execute(
            "SELECT COUNT(*) FROM locus_proteins WHERE locus_id = ?", [locus_id]
        )[0][0]
        window = max(locus_size // 2 + 5, 12)

        b.visualize_neighborhood(
            center[0],
            window=window,
            output_path=str(FIGURES_DIR / f"defense_island_{label}.png"),
            title=f"Defense Island: {label}",
            legend=description,
        )
        return str(FIGURES_DIR / f"defense_island_{label}.png")
    return None
```

---

### Hypothesis Tracking & Provenance

Log analytical steps and register defense-related hypotheses:

```python
# Log key analytical steps with provenance chaining
e1 = b.log_provenance("Defense inventory", f"{n_crispr} CRISPR, {n_rm} RM, {n_ta} TA systems")
e2 = b.log_provenance("Defense load calculation", f"Mean {mean_pct:.1f}% of proteome", parent_ids=[e1.entry_id])
e3 = b.log_provenance("Co-occurrence analysis", "CRISPR + RM co-occur in 80% of genomes", parent_ids=[e1.entry_id])

# Propose hypotheses about defense strategy
h = b.propose_hypothesis("Heavy defense investment indicates constant phage pressure in syntrophic niche")
b.add_evidence(h.hypothesis_id, "Defense load", f"{mean_pct:.1f}% mean defense proteome", True, 0.8)
b.add_evidence(h.hypothesis_id, "CRISPR universality", f"Type I-C in {n_crispr_genomes}/{n_genomes} genomes", True, 0.7)

# Review all hypotheses
print(b.hypothesis_summary())
print(b.list_hypotheses())

# Render provenance DAG
b.render_provenance(title="Defense Analysis", output_path="figures/defense_provenance.mermaid")
```

Hypotheses persist across sessions — `b.resume()` shows active hypotheses automatically.
