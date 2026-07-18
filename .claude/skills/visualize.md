# Visualize Skill

Generate publication-quality gene neighborhood diagrams and domain architecture figures from Sharur databases.

**This is a reference skill, not an agent skill.** It documents the two visualization tools available and when to use each. Consult this when generating figures during survey, exploration, deepening, or manuscript work.

---

## Tool 1: `plot_locus_multisource.py` (CLI -- publication figures)

**Script:** `scripts/plot_locus_multisource.py`

Multi-source annotation priority locus diagrams. Uses the full annotation stack to label genes with the most informative available annotation.

### CLI Usage

```bash
python scripts/plot_locus_multisource.py \
    --db data/DATASET/sharur.duckdb \
    --protein PROTEIN_ID \
    --window 10 \
    --output figures/figureN_locus_name.png \
    --title "Descriptive Locus Title" \
    --subtitle "Custom subtitle (optional)"
```

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `--db` | Yes | -- | Path to sharur.duckdb |
| `--protein` | Yes | -- | Center protein ID |
| `--window` | No | 10 | Genes on each side of center |
| `--output` | No | `locus_{protein_id}.png` | Output PNG path |
| `--title` | No | Auto-generated | Figure title |
| `--subtitle` | No | Annotation priority legend | Custom subtitle line |

### Annotation Priority Stack

Labels are assigned per-gene using the highest-priority source available:

1. **Foldseek** structural homology (PDB descriptions, cleaned automatically)
2. **DefenseFinder** system-level calls (loaded from `defensefinder_results/` or `annotations/`)
3. **PADLOC** defense calls (loaded from `padloc_results/` or `annotations/`)
4. **PFAM / KEGG / VOGdb** -- best by e-value from the `annotations` table
   - KEGG KO IDs are translated to human-readable descriptions via `data/reference/ko_list`

If no annotation source produces a hit, the gene is labeled `?` and colored gray.

### Color Scheme (by annotation source)

| Source | Color | Hex |
|--------|-------|-----|
| Foldseek | Red | `#FF6B6B` |
| DefenseFinder | Teal | `#4ECDC4` |
| PADLOC | Light teal | `#95E1D3` |
| PFAM | Blue | `#3498DB` |
| KEGG | Purple | `#9B59B6` |
| VOGdb | Orange | `#F39C12` |
| CAZy | Green | `#2ECC71` |
| HydDB | Red | `#E74C3C` |
| Unannotated | Light gray | `#BDC3C7` |

### Features

- Gene numbers below the track (1-indexed, `g1`, `g2`, ...) with thin gray leader lines
- Absolute genome coordinates on x-axis with comma formatting
- `annotate_inline=False` -- labels are placed externally to avoid overlap
- Query protein highlighted with a star marker and dark red color (`#C0392B`)
- CRISPR arrays detected automatically from the `loci` table (within 1kb of the gene window) and rendered as light-red boxes with repeat/spacer counts
- Dynamic legend showing only annotation sources used in the current figure
- Foldseek PDB descriptions are cleaned: experimental method prefixes, PDB IDs, organism names, and state descriptors are stripped to produce concise protein function labels
- Figure size: 16x6 inches at 150 dpi

### Best For

- Publication figures in manuscripts and reports
- Locus-level findings that need multi-source annotation context
- Neighborhood validation during Context-First protocol checks
- Defense island, prophage, or metabolic operon visualization

---

## Tool 2: `b.visualize_neighborhood()` (Python API -- quick exploration)

**Module:** `sharur/operators/visualization.py`

Lightweight neighborhood diagrams accessible directly from the Sharur API. Colors genes by functional keyword rather than annotation source.

### Python API

```python
from sharur.operators import Sharur
b = Sharur("data/DATASET/sharur.duckdb", read_only=True)

# Neighborhood diagram
result = b.visualize_neighborhood(
    "PROTEIN_ID",
    window=10,
    output_path="figures/locus.png",
    figure_width=14,      # default 14
    show_labels=True,     # default True
    show_crispr_arrays=True,  # default True
)
print(result.data)  # path to saved PNG

# Domain architecture diagram
result = b.visualize_domains("PROTEIN_ID", output_path="figures/domains.png")
```

### Color Scheme (by functional keyword)

| Keyword in annotation | Color | Category |
|----------------------|-------|----------|
| transporter, permease | Green | Transport |
| oxidoreductase, dehydrogenase | Blue | Redox |
| transferase | Purple | Transfer |
| hydrolase, peptidase | Orange | Hydrolysis |
| kinase | Pink | Phosphorylation |
| phage, prophage | Red | Phage |
| transposase, integrase | Brown | Mobile elements |
| crispr, cas | Cyan | Defense |
| unannotated / NO HITS | Gray | Unknown |
| hypothetical | Dark gray | Hypothetical |
| (default) | Blue-gray | Other |

### Features

- Query protein marked with `>>> LABEL <<<` and colored red
- Absolute genome coordinates on x-axis
- CRISPR arrays from `loci` table (light blue boxes)
- Auto-records to analysis manifest
- Annotations from PFAM/KEGG/VOGdb only (best by e-value); does NOT include Foldseek, DefenseFinder, or PADLOC
- Figure size: 14x5 inches at 150 dpi

### Best For

- Quick exploratory looks during analysis sessions
- In-session inspection when browsing genomes
- Rapid iteration when you need to check many neighborhoods fast

---

## When to Use Which Tool

| Situation | Tool | Reason |
|-----------|------|--------|
| Manuscript or report figure | `plot_locus_multisource.py` | Full annotation stack, publication styling |
| Neighborhood validation (Context-First) | `plot_locus_multisource.py` | Shows DefenseFinder/PADLOC/Foldseek context |
| Quick check during exploration | `b.visualize_neighborhood()` | Faster, no CLI overhead |
| Batch browsing many loci | `b.visualize_neighborhood()` | Scriptable in Python loops |
| Defense island characterization | `plot_locus_multisource.py` | Shows DefenseFinder + PADLOC annotations |
| Locus with Foldseek structural data | `plot_locus_multisource.py` | Only tool that renders Foldseek labels |

---

## Standard Conventions

### Window Sizes

| Locus type | Recommended window | Rationale |
|-----------|-------------------|-----------|
| Small operon (3-5 genes) | `window=8` | Enough flanking context without clutter |
| Medium locus (defense, metabolic cluster) | `window=10-12` | Standard for most figures |
| Large genomic island | `window=15-20` | Captures full island + flanking regions |
| Giant protein context | `window=12` | Shows neighbors that may reveal function |

### Figure Naming

```
figures/figureN_descriptive_name.png
```

Examples:
- `figures/figure3_crispr_type_ic_locus.png`
- `figures/figure7_dockerin_cohesin_neighborhood.png`
- `exploration/figures/GCA_018260655_defense_island.png`

### Manuscript Integration

Do NOT put "Figure N." in matplotlib titles -- pandoc figure captions handle numbering. In `MANUSCRIPT.md`:

```markdown
![Figure 3: Type I-C CRISPR-Cas locus in genome GCA_018260655 showing Cascade complex (Cas5/7/8c), Cas3 helicase-nuclease, and adjacent 42-spacer array.](figures/figure3_crispr_type_ic_locus.png)
```

When inserting a new figure, renumber all subsequent figures and their cross-references in the manuscript.

### Visual Style Rules

- `annotate_inline=False` -- always use external labels, never inline
- Gene numbers below the track, 1-indexed
- Absolute genome coordinates on x-axis with comma formatting
- Color by annotation source (multisource script) or functional category (operator)
- For ambiguous annotations (e.g., Cas12f vs TnpB), use honest labels -- do not resolve the ambiguity in the figure label

---

## System Figure Caveats

Build system heatmaps only from purpose-built caller output found in the live schema.
Raw per-profile/domain hits may be visualized as observed evidence, but must not be
relabeled as systems. If prevalence is surprising, inspect caller membership,
co-annotations, and validation status without priming the review with a named expected error.

---

## Finding Protein IDs for Visualization

Common patterns for locating anchor proteins to center a figure on:

```python
# By annotation accession (e.g., find a Cas3 protein)
pid = b.store.execute("""
    SELECT a.protein_id FROM annotations a
    JOIN proteins p ON a.protein_id = p.protein_id
    WHERE a.accession = 'PF01930' AND p.bin_id = 'GENOME_ID'
    ORDER BY a.evalue LIMIT 1
""")[0][0]

# By predicate (e.g., find a giant unannotated protein)
results = b.search_by_predicates(has=["giant", "unannotated"], limit=5)
pid = results.data[0]['protein_id']

# By gene name pattern
pid = b.store.execute("""
    SELECT a.protein_id FROM annotations a
    WHERE LOWER(a.name) LIKE '%cas3%'
    ORDER BY a.evalue LIMIT 1
""")[0][0]

# Center protein of a locus (e.g., defense island or prophage)
pid = b.store.execute("""
    SELECT lp.protein_id FROM locus_proteins lp
    WHERE lp.locus_id = 'LOCUS_ID'
    ORDER BY lp.position
    LIMIT 1 OFFSET (SELECT COUNT(*)/2 FROM locus_proteins WHERE locus_id = 'LOCUS_ID')
""")[0][0]

# By DefenseFinder system type
pid = b.store.execute("""
    SELECT protein_id FROM annotations
    WHERE source = 'defensefinder' AND accession LIKE 'CBASS%'
    LIMIT 1
""")[0][0]
```

---

## Domain Architecture Diagrams

For single-protein domain layout (not neighborhood context):

```python
result = b.visualize_domains("PROTEIN_ID", output_path="figures/domains_protein.png")
```

This shows PFAM domain annotations arranged along the protein length. Useful for characterizing multi-domain proteins and giant proteins with complex architecture. Note: domain coordinates are approximate (evenly spaced) since the annotations table does not store domain start/end positions.
