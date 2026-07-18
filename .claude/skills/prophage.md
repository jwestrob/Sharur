# Prophage & Viral Element Detection Skill

Identify integrated prophage regions, viral remnants, and potentially misbinned phage contigs within metagenomic assemblies. Outputs structured loci to the `loci`/`locus_proteins` tables and findings to `findings.jsonl`.

**CONCURRENCY:** Read-only detection may run in parallel, but this workflow persists
curated loci. Serialize its DuckDB write section and never overlap it with another writer.

> **Mandatory:** Follow the shared validation protocols in `_validation_protocols.md`.
> Use COUNT(DISTINCT protein_id) for protein counts. Verify annotation accessions before reporting.

> **Literature dispatch:** When you encounter ambiguous annotations, unknown Foldseek hits,
> or need to make comparative claims ("first known", "largest"), dispatch a literature agent.
> Read `.claude/skills/literature.md` for protocols.

---

## Usage

```
/prophage                          # Full prophage detection + characterization
/prophage --genome GENOME_ID       # Single genome analysis
/prophage --characterize           # Skip detection, characterize existing prophage loci
/prophage --misbinned              # Focus on identifying misbinned phage contigs
```

---

## Critical Background: Why VOGdb Density Is Misleading

**Do NOT use raw VOGdb hit density as a prophage signal.** This was the single biggest lesson from Omnitrophota prophage detection:

- 94% of VOGs are category Xu (uncharacterized) — they match housekeeping genes with viral orthologs
- Many core bacterial genes (helicases, primases, recombinases) have VOGdb hits
- Raw VOGdb density produces ~76% false positive prophage calls

**The solution is a three-layer filter:**
1. **Category markers**: VOGs with Xs (viral-specific) or Xr (viral replication) functional categories
2. **Keyword markers**: VOGs whose description matches phage-specific terms (terminase, capsid, portal, tail, baseplate, lysin, holin, integrase, etc.)
3. **Diagnostic gene triage**: Candidate regions must contain at least one diagnostic phage gene (terminase, capsid, portal, tail fiber, integrase, or lysis cassette)

**"Replication" is NOT diagnostic.** Helicases and primases are shared between phage and host — do not count them as phage-specific evidence.

---

## Prompt

You are detecting integrated viral elements (prophages) in a metagenomic dataset. Your goal is to find genuine prophage loci, distinguish them from viral remnants/artifacts, identify potentially misbinned phage contigs, and write results to the database's loci tables.

### Step 1: Reconnaissance

Before any detection, understand what you're working with:

```python
import sys, json, uuid
from pathlib import Path
from collections import defaultdict
from datetime import datetime

sys.path.insert(0, '.')
from sharur.operators import Sharur

DB_PATH = "data/DATASET/sharur.duckdb"  # Adjust as needed
b = Sharur(DB_PATH)
DB_DIR = Path(DB_PATH).parent

# Dataset scale
n_genomes = b.store.execute("SELECT COUNT(DISTINCT bin_id) FROM proteins")[0][0]
n_proteins = b.store.execute("SELECT COUNT(*) FROM proteins")[0][0]
n_contigs = b.store.execute("SELECT COUNT(DISTINCT contig_id) FROM proteins")[0][0]
print(f"Dataset: {n_genomes} genomes, {n_proteins:,} proteins, {n_contigs:,} contigs")

# Check available annotation sources
sources = b.store.execute("""
    SELECT source, COUNT(*) as n_hits, COUNT(DISTINCT protein_id) as n_proteins
    FROM annotations
    GROUP BY source
    ORDER BY n_hits DESC
""")
for src, n_hits, n_prot in sources:
    print(f"  {src}: {n_hits:,} hits on {n_prot:,} proteins")

# Check for VOGdb — REQUIRED for prophage detection
has_vogdb = any(src == 'vogdb' for src, _, _ in sources)
if not has_vogdb:
    print("ERROR: No VOGdb annotations found. Cannot run prophage detection.")
    print("Run Astra with VOGdb first, then load annotations.")
    # EXIT — cannot proceed without VOGdb

# Check for existing prophage loci
existing = b.store.execute("""
    SELECT locus_type, COUNT(*) FROM loci
    WHERE locus_type IN ('prophage', 'viral_contig')
    GROUP BY locus_type
""")
if existing:
    print(f"Existing loci: {existing}")
    print("Will clear and re-detect unless --characterize mode")
```

### Step 2: Build VOGdb Marker Sets

Classify VOGs into phage-specific markers vs generic viral orthologs. This requires the VOGdb annotations reference file.

```python
import re

# Locate vog.annotations.tsv — check common paths
VOG_ANNOT_PATHS = [
    DB_DIR / "vog.annotations.tsv",
    DB_DIR.parent / "reference" / "vogdb" / "vog.annotations.tsv",
    Path("data/reference/vogdb/vog.annotations.tsv"),
    Path("~/.config/Astra/VOGdb/vog.annotations.tsv").expanduser(),
]

vog_annot_path = None
for p in VOG_ANNOT_PATHS:
    if p.exists():
        vog_annot_path = p
        break

if vog_annot_path is None:
    print("WARNING: vog.annotations.tsv not found — falling back to keyword-only detection")
    # Can still detect prophages using annotation name/description fields in the DB

# Parse VOGdb functional categories and descriptions
PROPHAGE_KEYWORDS = re.compile(
    r"\b("
    r"terminase|capsid|head|portal|tail|baseplate|spike|"
    r"lysin|holin|endolysin|spanin|"
    r"integrase|recombinase|excisionase|"
    r"packaging|scaffolding|"
    r"tape.measure|neck|sheath|tube"
    r")\b",
    re.IGNORECASE,
)

category_markers = set()  # Xs or Xr functional category
keyword_markers = set()   # Description matches phage keywords

if vog_annot_path:
    with open(vog_annot_path) as fh:
        header = fh.readline()
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            vog_id, func_cat, desc = parts[0], parts[3], parts[4]

            if "Xs" in func_cat or "Xr" in func_cat:
                category_markers.add(vog_id)
            if PROPHAGE_KEYWORDS.search(desc):
                keyword_markers.add(vog_id)

    all_markers = category_markers | keyword_markers
    print(f"VOGdb markers: {len(category_markers)} by category (Xs/Xr), "
          f"{len(keyword_markers)} by keyword, {len(all_markers)} total unique")
else:
    all_markers = set()
```

### Step 3: Tag Genes and Detect Prophage Regions

For each contig, walk through genes in order. Tag each gene as marker (phage-specific VOG), support (any VOGdb hit), or negative. Cluster markers using gap-based algorithm.

```python
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

# Load VOGdb annotations + cross-source annotations for integrase/terminase keywords
annot_rows = b.store.execute("""
    SELECT protein_id, source, accession, name, description
    FROM annotations
    WHERE source IN ('vogdb', 'pfam', 'kegg', 'defensefinder')
""")

protein_annots = defaultdict(list)
for pid, src, acc, name, desc in annot_rows:
    protein_annots[pid].append((src, acc, name or '', desc or ''))

# Tag each gene
INTEGRASE_KW = re.compile(r"\b(integrase|recombinase|excisionase)\b", re.IGNORECASE)
TERMINASE_KW = re.compile(r"\b(terminase)\b", re.IGNORECASE)

# Diagnostic phage gene keywords — these are what distinguish real prophages from remnants
DIAGNOSTIC_KEYWORDS = re.compile(
    r"\b(terminase|capsid|portal|tail|baseplate|head.morphogenesis|"
    r"integrase|lysis|lysin|holin|endolysin)\b",
    re.IGNORECASE,
)
# NOTE: "replication" (helicase/primase) is deliberately EXCLUDED — too many bacterial orthologs

def tag_gene(pid, annots, marker_vogs):
    """Tag a gene for prophage detection. Returns (tag, extras)."""
    extras = {
        "has_integrase": False, "has_terminase": False,
        "has_diagnostic": False, "diagnostic_types": set(),
        "marker_accessions": [],
    }
    is_marker = False
    is_support = False

    for source, accession, name, desc in annots:
        text = f"{name} {desc}"

        if source != "vogdb":
            # Check non-VOGdb sources for integrase/terminase
            if INTEGRASE_KW.search(text):
                extras["has_integrase"] = True
                extras["has_diagnostic"] = True
                extras["diagnostic_types"].add("integrase")
            if TERMINASE_KW.search(text):
                extras["has_terminase"] = True
                extras["has_diagnostic"] = True
                extras["diagnostic_types"].add("terminase")
            continue

        # VOGdb annotation
        is_support = True
        if accession in marker_vogs:
            is_marker = True
            extras["marker_accessions"].append(accession)

        # Check VOGdb names for diagnostic keywords
        if DIAGNOSTIC_KEYWORDS.search(text):
            extras["has_diagnostic"] = True
            for kw in ["terminase", "capsid", "portal", "tail", "baseplate",
                       "integrase", "lysis", "lysin", "holin", "endolysin"]:
                if re.search(rf"\b{kw}\b", text, re.IGNORECASE):
                    extras["diagnostic_types"].add(kw)
            if TERMINASE_KW.search(text):
                extras["has_terminase"] = True
            if INTEGRASE_KW.search(text):
                extras["has_integrase"] = True

    if is_marker:
        return "marker", extras
    if is_support:
        return "support", extras
    return "negative", extras

# Tag all genes
for cid, genes in contig_genes.items():
    genes.sort(key=lambda g: g["gene_index"])
    for gene in genes:
        annots = protein_annots.get(gene["protein_id"], [])
        tag, extras = tag_gene(gene["protein_id"], annots, all_markers)
        gene["tag"] = tag
        gene["extras"] = extras
```

### Step 4: Gap-Based Island Detection

Walk each contig and cluster marker/support genes, allowing small gaps of negative genes. Then apply the diagnostic gene filter.

```python
# Detection parameters — adjust based on dataset
MAX_GAP = 5        # Max consecutive non-viral genes within a prophage
MIN_MARKERS = 3    # Min marker (phage-specific) genes
MIN_GENES = 5      # Min total positive (marker + support) genes

def detect_prophage_regions(genes, max_gap=MAX_GAP, min_markers=MIN_MARKERS, min_genes=MIN_GENES):
    """Gap-based clustering on sorted gene list from one contig."""
    islands = []
    current = []
    gap_count = 0

    def emit(cluster):
        marker_count = sum(1 for g in cluster if g["tag"] == "marker")
        positive_count = sum(1 for g in cluster if g["tag"] in ("marker", "support"))
        if marker_count >= min_markers and positive_count >= min_genes:
            islands.append({
                "genes": list(cluster),
                "marker_count": marker_count,
                "positive_count": positive_count,
            })

    for gene in genes:
        if gene["tag"] in ("marker", "support"):
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

# Run detection across all contigs
all_candidates = []
for cid, genes in contig_genes.items():
    islands = detect_prophage_regions(genes)
    for island in islands:
        island["contig_id"] = cid
        island["bin_id"] = genes[0]["bin_id"]
    all_candidates.extend(islands)

print(f"Candidate prophage regions before triage: {len(all_candidates)}")
```

### Step 5: Diagnostic Gene Triage (CRITICAL FILTER)

This is where most false positives are eliminated. A real prophage needs at least one diagnostic phage gene — not just generic viral orthologs.

```python
def triage_prophage(island):
    """Classify a candidate prophage region.

    Returns: (verdict, confidence, diagnostic_genes)
        verdict: 'real_prophage', 'degraded_prophage', 'viral_remnant', 'false_positive'
        confidence: 0.0-1.0
        diagnostic_genes: set of diagnostic gene types found
    """
    diagnostic_types = set()
    has_integrase = False
    has_terminase = False

    for g in island["genes"]:
        diagnostic_types.update(g["extras"].get("diagnostic_types", set()))
        if g["extras"].get("has_integrase"):
            has_integrase = True
        if g["extras"].get("has_terminase"):
            has_terminase = True

    marker_count = island["marker_count"]
    total_genes = len(island["genes"])

    # Classification logic
    structural_genes = diagnostic_types - {"integrase", "lysis", "lysin", "holin", "endolysin"}
    has_structural = len(structural_genes) > 0  # capsid, portal, tail, terminase, baseplate

    if has_terminase and has_structural and marker_count >= 5:
        return "real_prophage", 0.95, diagnostic_types
    if (has_integrase or has_terminase) and has_structural:
        return "real_prophage", 0.85, diagnostic_types
    if has_structural and marker_count >= 5:
        return "real_prophage", 0.80, diagnostic_types
    if has_integrase and marker_count >= 5:
        return "degraded_prophage", 0.70, diagnostic_types
    if len(diagnostic_types) >= 1 and marker_count >= 3:
        return "degraded_prophage", 0.60, diagnostic_types
    if marker_count >= 3 and total_genes >= 8:
        return "viral_remnant", 0.50, diagnostic_types
    return "false_positive", 0.30, diagnostic_types

# Apply triage
real_prophages = []
degraded_prophages = []
viral_remnants = []
false_positives = []

for island in all_candidates:
    verdict, confidence, diag = triage_prophage(island)
    island["verdict"] = verdict
    island["confidence"] = confidence
    island["diagnostic_types"] = diag

    if verdict == "real_prophage":
        real_prophages.append(island)
    elif verdict == "degraded_prophage":
        degraded_prophages.append(island)
    elif verdict == "viral_remnant":
        viral_remnants.append(island)
    else:
        false_positives.append(island)

print(f"\nTriage results:")
print(f"  Real prophages: {len(real_prophages)}")
print(f"  Degraded prophages: {len(degraded_prophages)}")
print(f"  Viral remnants: {len(viral_remnants)}")
print(f"  False positives (discarded): {len(false_positives)}")
```

### Step 6: Misbinned Phage Contig Detection

Some contigs may be entire phage genomes accidentally binned into bacterial MAGs. These show:
- Small contig size (typically <100 genes)
- Very high VOGdb marker density (>40% of genes are markers)
- Diagnostic phage genes present
- Unusual GC content compared to host genome

```python
def check_misbinned_phages():
    """Find contigs that might be entire phage genomes."""

    misbinned_candidates = []

    # Per-contig statistics
    contig_stats = b.store.execute("""
        SELECT p.contig_id, p.bin_id,
               COUNT(DISTINCT p.protein_id) as n_genes,
               MIN(p.start) as contig_start,
               MAX(p.end_coord) as contig_end
        FROM proteins p
        WHERE p.gene_index IS NOT NULL
        GROUP BY p.contig_id, p.bin_id
        HAVING COUNT(DISTINCT p.protein_id) BETWEEN 10 AND 200
    """)

    for cid, bid, n_genes, cstart, cend in contig_stats:
        genes = contig_genes.get(cid, [])
        if not genes:
            continue

        # Count markers and support genes
        n_markers = sum(1 for g in genes if g.get("tag") == "marker")
        n_support = sum(1 for g in genes if g.get("tag") == "support")
        n_positive = n_markers + n_support

        marker_fraction = n_markers / n_genes if n_genes else 0
        positive_fraction = n_positive / n_genes if n_genes else 0

        # Collect diagnostic genes
        diag_types = set()
        for g in genes:
            diag_types.update(g.get("extras", {}).get("diagnostic_types", set()))

        # Criteria for misbinned phage:
        # - >40% of genes are phage markers
        # - OR >60% of genes have any VOGdb hit AND diagnostic genes present
        # - Small-to-medium contig (10-200 genes)
        if (marker_fraction > 0.4 and len(diag_types) >= 2) or \
           (positive_fraction > 0.6 and len(diag_types) >= 3 and n_genes <= 150):
            misbinned_candidates.append({
                "contig_id": cid,
                "bin_id": bid,
                "n_genes": n_genes,
                "n_markers": n_markers,
                "n_positive": n_positive,
                "marker_fraction": round(marker_fraction, 3),
                "positive_fraction": round(positive_fraction, 3),
                "diagnostic_types": diag_types,
                "contig_length_bp": cend - cstart if cstart and cend else None,
            })

    return misbinned_candidates

misbinned = check_misbinned_phages()
print(f"\nPotentially misbinned phage contigs: {len(misbinned)}")
for m in misbinned[:10]:
    print(f"  {m['contig_id']} ({m['bin_id']}): {m['n_genes']} genes, "
          f"{m['marker_fraction']:.0%} markers, diag: {m['diagnostic_types']}")
```

### Step 7: Write Prophage Loci to Database

Write detected prophages to the `loci` and `locus_proteins` tables.

```python
def write_prophage_loci(prophages, locus_type="prophage", clear_existing=True):
    """Write prophage regions to loci/locus_proteins tables."""

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
    for island in prophages:
        genes = island["genes"]
        locus_id = f"{locus_type}_{uuid.uuid4().hex[:12]}"

        start = min(g["start"] for g in genes)
        end_coord = max(g["end_coord"] for g in genes)

        metadata = {
            "verdict": island["verdict"],
            "marker_count": island["marker_count"],
            "positive_count": island["positive_count"],
            "total_genes": len(genes),
            "diagnostic_types": sorted(island.get("diagnostic_types", set())),
            "has_integrase": any(g["extras"].get("has_integrase") for g in genes),
            "has_terminase": any(g["extras"].get("has_terminase") for g in genes),
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
    print(f"Wrote {len(prophages)} {locus_type} loci with {n_links:,} protein links")

# Write real + degraded prophages (not remnants or FPs)
all_prophages = real_prophages + degraded_prophages
write_prophage_loci(all_prophages, locus_type="prophage")

# Optionally write misbinned phages as a separate locus type
if misbinned:
    # Convert misbinned contig candidates to locus format
    misbinned_loci = []
    for m in misbinned:
        genes = contig_genes.get(m["contig_id"], [])
        if genes:
            misbinned_loci.append({
                "genes": genes,
                "contig_id": m["contig_id"],
                "verdict": "misbinned_phage",
                "confidence": 0.75,
                "marker_count": m["n_markers"],
                "positive_count": m["n_positive"],
                "diagnostic_types": m["diagnostic_types"],
            })
    if misbinned_loci:
        write_prophage_loci(misbinned_loci, locus_type="viral_contig", clear_existing=True)
```

### Step 8: Characterize Prophage Loci

Go beyond detection — characterize what kinds of prophages are present.

```python
def characterize_prophages():
    """Analyze prophage loci for biological interpretation."""

    # Load all prophage loci
    loci = b.store.execute("""
        SELECT l.locus_id, l.contig_id, l.confidence, l.metadata,
               COUNT(lp.protein_id) as n_genes
        FROM loci l
        JOIN locus_proteins lp ON l.locus_id = lp.locus_id
        WHERE l.locus_type = 'prophage'
        GROUP BY l.locus_id, l.contig_id, l.confidence, l.metadata
    """)

    print(f"\n=== PROPHAGE CHARACTERIZATION ===")
    print(f"Total prophage loci: {len(loci)}")

    # Size distribution
    sizes = [n for _, _, _, _, n in loci]
    print(f"Size (genes): min={min(sizes)}, median={sorted(sizes)[len(sizes)//2]}, "
          f"max={max(sizes)}, mean={sum(sizes)/len(sizes):.1f}")

    # Genome distribution
    genome_dist = b.store.execute("""
        SELECT p.bin_id, COUNT(DISTINCT l.locus_id) as n_prophages
        FROM loci l
        JOIN locus_proteins lp ON l.locus_id = lp.locus_id
        JOIN proteins p ON lp.protein_id = p.protein_id
        WHERE l.locus_type = 'prophage'
        GROUP BY p.bin_id
        ORDER BY n_prophages DESC
    """)

    genomes_with = len(genome_dist)
    pct = genomes_with / n_genomes * 100
    print(f"Genomes with prophages: {genomes_with}/{n_genomes} ({pct:.1f}%)")

    # Confidence distribution
    for lid, cid, conf, meta_json, n_genes in loci:
        meta = json.loads(meta_json) if isinstance(meta_json, str) else meta_json
        # Classify by completeness
        # ... analyze diagnostic gene composition

    # Diagnostic gene composition across all prophages
    all_diag = defaultdict(int)
    for lid, cid, conf, meta_json, n_genes in loci:
        meta = json.loads(meta_json) if isinstance(meta_json, str) else meta_json
        for dt in meta.get("diagnostic_types", []):
            all_diag[dt] += 1

    print("\nDiagnostic gene prevalence across prophage loci:")
    for dt, n in sorted(all_diag.items(), key=lambda x: -x[1]):
        print(f"  {dt}: {n}/{len(loci)} ({n/len(loci)*100:.0f}%)")

    return genome_dist, all_diag

genome_dist, diag_composition = characterize_prophages()
```

### Step 9: Co-occurrence with Defense Islands

Prophages and defense systems co-evolve. Check for spatial association.

```python
# Check if defense islands exist in the database
defense_count = b.store.execute("""
    SELECT COUNT(*) FROM loci WHERE locus_type = 'island'
""")[0][0]

if defense_count > 0:
    # Find prophages near defense islands (within 20 genes)
    colocalized = b.store.execute("""
        WITH prophage_bounds AS (
            SELECT l.locus_id as prophage_id, p.contig_id, p.bin_id,
                   MIN(p.gene_index) as p_start, MAX(p.gene_index) as p_end
            FROM loci l
            JOIN locus_proteins lp ON l.locus_id = lp.locus_id
            JOIN proteins p ON lp.protein_id = p.protein_id
            WHERE l.locus_type = 'prophage'
            GROUP BY l.locus_id, p.contig_id, p.bin_id
        ),
        defense_bounds AS (
            SELECT l.locus_id as defense_id, p.contig_id,
                   MIN(p.gene_index) as d_start, MAX(p.gene_index) as d_end
            FROM loci l
            JOIN locus_proteins lp ON l.locus_id = lp.locus_id
            JOIN proteins p ON lp.protein_id = p.protein_id
            WHERE l.locus_type = 'island'
            GROUP BY l.locus_id, p.contig_id
        )
        SELECT pb.prophage_id, db.defense_id, pb.bin_id,
               pb.contig_id, pb.p_start, pb.p_end, db.d_start, db.d_end
        FROM prophage_bounds pb
        JOIN defense_bounds db ON pb.contig_id = db.contig_id
            AND (ABS(pb.p_start - db.d_end) <= 20 OR ABS(db.d_start - pb.p_end) <= 20)
    """)

    print(f"\nProphage-defense island co-localization: {len(colocalized)} pairs")
    # This is biologically expected — defense systems often cluster near prophage insertion sites
```

### Step 10: Visualize Representative Loci

Generate neighborhood figures for the most interesting prophage loci.

```python
FIGURES_DIR = DB_DIR / "exploration" / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Pick representative loci to visualize:
# 1. Highest-confidence complete prophage
# 2. Largest prophage
# 3. A degraded prophage (for comparison)
# 4. A misbinned phage contig (if any)

# Get a representative protein from each locus for visualization
def visualize_prophage_locus(locus_id, label, description):
    """Visualize a prophage locus using the center protein."""
    center = b.store.execute("""
        SELECT lp.protein_id FROM locus_proteins lp
        JOIN proteins p ON lp.protein_id = p.protein_id
        WHERE lp.locus_id = ?
        ORDER BY lp.position
        LIMIT 1 OFFSET (SELECT COUNT(*)/2 FROM locus_proteins WHERE locus_id = ?)
    """, [locus_id, locus_id])

    if center:
        # Calculate window to cover the full locus
        locus_size = b.store.execute(
            "SELECT COUNT(*) FROM locus_proteins WHERE locus_id = ?", [locus_id]
        )[0][0]
        window = max(locus_size // 2 + 3, 12)  # Ensure full coverage + flanking

        b.visualize_neighborhood(
            center[0],
            window=window,
            output_path=str(FIGURES_DIR / f"prophage_{label}.png"),
            title=f"Prophage Locus: {label}",
            legend=description,
        )
        return str(FIGURES_DIR / f"prophage_{label}.png")
    return None

# Visualize top 3-5 representative loci
# (Select based on your detection results)
```

### Step 11: Log Findings

Write structured findings to findings.jsonl with full provenance.

```python
# Set up findings logging
EXPLORE_DIR = DB_DIR / "exploration"
EXPLORE_DIR.mkdir(exist_ok=True)
FINDINGS_FILE = EXPLORE_DIR / "findings.jsonl"

from sharur.core.analysis_record_io import append_finding_record

def log_finding(category, title, description, evidence, verification=None,
                finding_id=None, n_genomes=None, priority="medium",
                provenance=None, figures=None, related_findings=None):
    finding = {
        "timestamp": datetime.now().isoformat(),
        "category": category,
        "title": title,
        "description": description,
        "evidence": evidence,
        "verification": verification,
        "n_genomes": n_genomes,
        "priority": priority,
        "provenance": provenance or {},
        "figures": figures or [],
        "related_findings": related_findings or [],
        "phase": "exploration",
    }
    if finding_id:
        finding["id"] = finding_id
    result = append_finding_record(FINDINGS_FILE, finding, phase="exploration")
    print(f"[LOGGED] {result.finding['id']}: {title}")
    return result.finding

# verification=None exists only so the strict canonical validator can emit a
# useful error. Every real call must supply claim/query/expected records.

# Log the main prophage finding
genomes_with_prophage = len(genome_dist)
log_finding(
    finding_id="prophage-001",
    category="prophage",
    title=f"Prophage detection: {len(real_prophages)} complete + {len(degraded_prophages)} degraded prophages in {genomes_with_prophage}/{n_genomes} genomes ({genomes_with_prophage/n_genomes*100:.1f}%)",
    description=f"Gap-based clustering of VOGdb markers (category Xs/Xr + keyword matching, "
        f"{len(all_markers)} marker VOGs) identified {len(all_candidates)} candidate regions. "
        f"Diagnostic gene triage (requiring terminase, capsid, portal, tail, integrase, or lysis) "
        f"retained {len(real_prophages)} high-confidence and {len(degraded_prophages)} degraded "
        f"prophages, discarding {len(false_positives)} false positives and {len(viral_remnants)} "
        f"remnants. Prophages are present in {genomes_with_prophage} of {n_genomes} genomes.",
    evidence={
        "n_candidates": len(all_candidates),
        "n_real": len(real_prophages),
        "n_degraded": len(degraded_prophages),
        "n_remnants": len(viral_remnants),
        "n_false_positives": len(false_positives),
        "n_genomes_with_prophage": genomes_with_prophage,
        "pct_genomes": round(genomes_with_prophage / n_genomes * 100, 1),
        "diagnostic_composition": dict(diag_composition),
    },
    verification=verification_records,
    n_genomes=genomes_with_prophage,
    priority="high",
    provenance={
        "query": "VOGdb marker classification + gap-based clustering + diagnostic gene triage",
        "raw_result": f"{len(all_candidates)} candidates → {len(real_prophages)+len(degraded_prophages)} retained",
        "accession_verified": "VOGdb Xs/Xr categories + keyword matching from vog.annotations.tsv",
        "interpretation": f"{genomes_with_prophage}/{n_genomes} genomes harbor detectable prophages",
    },
)

# Log misbinned phage finding if any
if misbinned:
    log_finding(
        finding_id="prophage-002",
        category="prophage",
        title=f"Potentially misbinned phage contigs: {len(misbinned)} contigs with >40% phage marker density",
        description=f"Identified {len(misbinned)} contigs where phage-specific VOGdb markers "
            f"constitute >40% of all genes, suggesting these may be complete phage genomes "
            f"accidentally binned into bacterial MAGs rather than integrated prophages.",
        evidence={
            "n_misbinned": len(misbinned),
            "contigs": [m["contig_id"] for m in misbinned[:10]],
        },
        verification=misbinned_verification_records,
        n_genomes=len(set(m["bin_id"] for m in misbinned)),
        priority="high",
    )
```

---

## Output Format

```markdown
## Prophage & Viral Element Analysis

### Detection Summary
- **Method**: VOGdb marker classification (N category + N keyword markers) → gap-based clustering → diagnostic gene triage
- **Candidates identified**: N regions
- **After triage**: N real prophages, N degraded, N remnants (discarded), N false positives (discarded)
- **Genomes with prophages**: N/M (X%)

### Prophage Size Distribution
- Min: N genes | Median: N genes | Max: N genes
- N with terminase, N with integrase, N with both

### Diagnostic Gene Composition
| Gene Type | Loci with Gene | % of Prophages |
|-----------|---------------|----------------|
| terminase | N | X% |
| capsid | N | X% |
| portal | N | X% |
| tail | N | X% |
| integrase | N | X% |
| lysis/lysin/holin | N | X% |

### Genome Distribution
- N genomes with 1 prophage, N with 2, N with 3+
- Top genomes by prophage load: [list]

### Misbinned Phage Contigs
[If any detected — list contigs with marker fractions and diagnostic gene types]

### Defense Co-evolution
[If defense islands detected — co-localization statistics]

### Biological Interpretation
[Phage pressure, lifestyle implications, comparative context]

### Figures
- prophage_complete_example.png: Representative complete prophage with terminase + capsid
- prophage_degraded_example.png: Degraded prophage showing gene erosion
- [Additional locus figures]
```

---

## Adaptation Notes

This skill adapts to the dataset:

- **No VOGdb annotations?** → Cannot run. Report this and suggest running Astra with VOGdb.
- **No vog.annotations.tsv?** → Fall back to keyword-only marker detection using the `name` and `description` fields in the annotations table. Less sensitive but still functional.
- **Small dataset (<5 genomes)?** → Reduce MIN_MARKERS to 2, increase figure output per genome.
- **Large dataset (>500 genomes)?** → Use SQL aggregation for statistics, sample representative loci for figures (don't visualize all), batch the contig walk if memory is a concern.
- **Defense islands already detected?** → Run co-localization analysis (Step 9).
- **No loci tables?** → Report the error. The loci/locus_proteins tables are part of the standard Sharur schema.

---

## Findings Categories

Use these finding IDs:
- `prophage-NNN` — prophage detection and characterization findings
- Category: `prophage` for all prophage-related findings
- Category: `misbinned_phage` for potentially misbinned contigs
- Category: `phage_defense_coevolution` for defense co-occurrence findings
