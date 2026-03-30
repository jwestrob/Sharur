# Hydrogenase Skill

Validate, classify, and interpret [NiFe]- and [FeFe]-hydrogenases in metagenomic datasets with mandatory neighborhood curation for all Group 4 calls.

**CONCURRENCY: DuckDB does not support concurrent writes. Only ONE agent should access a database at a time. The coordinator must run DB-accessing skills sequentially, not in parallel.**

> **Mandatory:** Follow the shared validation protocols in `_validation_protocols.md`.
> Verify accession names before reporting. Use COUNT(DISTINCT protein_id) for protein
> counts. Apply Context-First protocol for annotations averaging >10 hits/genome.

> **Literature dispatch:** When you encounter ambiguous annotations, unknown Foldseek hits,
> or need to make comparative claims ("first known", "largest"), dispatch a literature agent.
> Read `.claude/skills/literature.md` for protocols.

---

## BLOCKING REQUIREMENT: Neighborhood Validation

**You CANNOT report hydrogenase results without completing the full validation protocol below.** HydDB classification alone produces massive false positive rates. In real datasets, 20 of 22 "Group 4" hydrogenase calls were actually Complex I (NADH:ubiquinone oxidoreductase) NuoD subunits.

**Every Group 4 HydDB hit MUST have neighborhood validation before it appears in any finding, report, or manuscript.** There are no exceptions. If you cannot run neighborhood validation (e.g., database unavailable), you MUST state "HydDB Group 4 calls are unvalidated and likely contain >90% false positives" rather than reporting them as hydrogenases.

Groups 1-3 with PF00374 corroboration are high-confidence and do not require neighborhood curation. Groups 1-3 WITHOUT PF00374 should still be spot-checked (3-5 examples).

---

## Why This Exists

The catalytic subunit HMM for Group 4 [NiFe]-hydrogenases hits NuoD (Complex I subunit D, K00333) because they are homologous. Group 4 encompasses Ech, Mbh, Eha, Hyf, Hyc, and other energy-converting hydrogenases -- ALL of which share the catalytic subunit fold with respiratory Complex I.

PF00374 (NiFeSe_Hases) is diagnostic for Groups 1-3 ONLY. Group 4 proteins will NEVER have PF00374 corroboration, making neighborhood validation mandatory for every Group 4 call.

Without operon validation, a "Group 4g Ech" call is indistinguishable from a Complex I NuoD subunit.

---

## Usage

```
/hydrogenase                     # Full hydrogenase inventory with validation
/hydrogenase --group4-only       # Focus on Group 4 curation
/hydrogenase --genome GENOME_ID  # Single genome analysis
```

---

## Prompt

You are validating hydrogenase classifications in a metagenomic dataset. Every Group 4 HydDB call must be neighborhood-validated before reporting. Follow the steps below in order.

### Step 1: Inventory All HydDB Hits

```python
from sharur.operators import Sharur
b = Sharur("data/DATASET/sharur.duckdb")

# Total HydDB NiFe and FeFe hits
nife_all = b.store.execute("""
    SELECT COUNT(DISTINCT protein_id) FROM annotations
    WHERE source = 'hyddb' AND LOWER(name) LIKE '%nife%'
""")[0][0]

fefe_all = b.store.execute("""
    SELECT COUNT(DISTINCT protein_id) FROM annotations
    WHERE source = 'hyddb' AND LOWER(name) LIKE '%fefe%'
""")[0][0]

print(f"HydDB NiFe hits (raw, unvalidated): {nife_all}")
print(f"HydDB FeFe hits (raw): {fefe_all}")

# Subgroup breakdown from predicates
for pred in ['nife_group1', 'nife_group2', 'nife_group3', 'nife_group4',
             'fefe_groupA', 'fefe_groupB', 'fefe_groupC']:
    count = b.store.execute(f"""
        SELECT COUNT(DISTINCT protein_id) FROM protein_predicates
        WHERE '{pred}' = ANY(predicates)
    """)[0][0]
    if count > 0:
        print(f"  {pred}: {count}")

# Flag: how many need curation?
needs_curation = b.store.execute("""
    SELECT COUNT(DISTINCT protein_id) FROM protein_predicates
    WHERE 'hyddb_needs_curation' = ANY(predicates)
""")[0][0]
print(f"\nFlagged for curation (hyddb_needs_curation): {needs_curation}")
```

### Step 2: Separate High-Confidence from Mandatory-Curation Hits

```python
# HIGH CONFIDENCE: Groups 1-3 with PF00374 corroboration
validated_g13 = b.store.execute("""
    SELECT DISTINCT a_hyd.protein_id
    FROM annotations a_hyd
    JOIN annotations a_pfam ON a_hyd.protein_id = a_pfam.protein_id
    WHERE a_hyd.source = 'hyddb'
      AND a_pfam.accession = 'PF00374'
""")
validated_pids = set(r[0] for r in validated_g13)
print(f"PF00374-corroborated (Groups 1-3, high confidence): {len(validated_pids)}")

# HIGH CONFIDENCE: FeFe with PF02906 corroboration
fefe_validated = b.store.execute("""
    SELECT DISTINCT a_hyd.protein_id
    FROM annotations a_hyd
    JOIN annotations a_pfam ON a_hyd.protein_id = a_pfam.protein_id
    WHERE a_hyd.source = 'hyddb'
      AND a_pfam.accession IN ('PF02906', 'PF02256')
""")
fefe_validated_pids = set(r[0] for r in fefe_validated)
print(f"PF02906/PF02256-corroborated FeFe (high confidence): {len(fefe_validated_pids)}")

# MANDATORY CURATION: All Group 4 hits (NONE will have PF00374)
group4_all = b.store.execute("""
    SELECT DISTINCT pp.protein_id
    FROM protein_predicates pp
    WHERE 'nife_group4' = ANY(pp.predicates)
""")
group4_pids = [r[0] for r in group4_all]
print(f"Group 4 hits requiring neighborhood validation: {len(group4_pids)}")

# Also check: any NiFe hits without PF00374 that are NOT Group 4
# These may be misclassified or divergent -- spot-check them
uncorroborated_g13 = b.store.execute("""
    SELECT DISTINCT a.protein_id
    FROM annotations a
    JOIN protein_predicates pp ON a.protein_id = pp.protein_id
    WHERE a.source = 'hyddb'
      AND LOWER(a.name) LIKE '%nife%'
      AND a.protein_id NOT IN (
          SELECT DISTINCT protein_id FROM annotations WHERE accession = 'PF00374'
      )
      AND NOT ('nife_group4' = ANY(pp.predicates))
""")
uncorroborated_g13_pids = [r[0] for r in uncorroborated_g13]
print(f"Groups 1-3 without PF00374 (spot-check 3-5): {len(uncorroborated_g13_pids)}")
```

### Step 3: Neighborhood Validation for ALL Group 4 Hits (MANDATORY)

This is the core of the skill. Every Group 4 hit gets classified as **confirmed**, **rejected**, or **ambiguous**.

```python
from collections import defaultdict

# Complex I rejection markers (KEGG KOs)
COMPLEX_I_KOS = {
    'K00330', 'K00331', 'K00332', 'K00333', 'K00334', 'K00335',
    'K00336', 'K00337', 'K00338', 'K00339', 'K00340', 'K00341',
    'K00342', 'K00343',  # nuoA through nuoN
}

# Complex I rejection markers (PFAM accessions)
COMPLEX_I_PFAMS = {
    'PF00329',   # Complex1_30kDa
    'PF00346',   # Complex1_49kDa
    'PF01058',   # Oxidored_q6 (NADH:ubiquinone oxidoreductase)
    'PF00361',   # Proton_antipo_M (proton antiporter membrane subunit)
    'PF01059',   # Oxidored_q5_N
    'PF00507',   # Oxidored_q4 (NADH-ubiquinone oxidoreductase 24kDa)
}

# Complex I rejection markers (annotation names -- case insensitive)
COMPLEX_I_NAMES = [
    'nadhd', 'nadh dehydrogenase', 'nadh-ubiquinone', 'nadh:ubiquinone',
    'complex i', 'nuo', 'oxidored_q',
]

# Hyf confirmation markers (Group 4f -- hydrogenase-4)
HYF_KOS = {
    'K12136', 'K12137', 'K12138', 'K12139', 'K12140',
    'K12141', 'K12142', 'K12143', 'K12144', 'K12145',  # hyfA-J
}

# Hyc confirmation markers (formate hydrogenlyase)
HYC_KOS = {
    'K15828', 'K15829', 'K15830', 'K15831', 'K15832', 'K15833',  # hycB-G
}

# Ech confirmation markers (PFAMs and names)
ECH_NAMES = [
    'echabcdef', 'ech hydrogenase', 'energy-converting hydrogenase',
    'eche', 'echf', 'echa', 'echb', 'echc', 'echd',
]

# Maturation factors (confirm ANY hydrogenase nearby)
MATURATION_KOS = {
    'K04651', 'K04652', 'K04653', 'K04654', 'K04655', 'K04656',  # HypA-F
    'K03605',  # HycI endopeptidase
}

# General hydrogenase markers
HYDROGENASE_PFAMS = {
    'PF00374',   # NiFeSe_Hases (unlikely on Group 4 but confirmatory on neighbor)
    'PF02906',   # Fe_hyd_lg_C
    'PF14720',   # Fer4_20 (4Fe-4S in hydrogenase maturation)
}


def validate_group4_neighborhood(protein_id, window=6):
    """Validate a single Group 4 HydDB hit by neighborhood analysis.

    Returns: (verdict, evidence_summary, details)
        verdict: 'confirmed', 'rejected', or 'ambiguous'
        evidence_summary: one-line rationale
        details: dict with all markers found
    """
    nbr = b.get_neighborhood(protein_id, window=window, all_annotations=True)

    complex_i_evidence = []
    hydrogenase_evidence = []

    for gene in nbr:
        pid = gene.get('protein_id', '')
        if pid == protein_id:
            continue  # skip the query protein itself

        annots = gene.get('annotations', [])
        for annot in annots:
            acc = annot.get('accession', '')
            name = annot.get('name', '')
            name_lower = name.lower() if name else ''

            # Check Complex I markers
            if acc in COMPLEX_I_KOS:
                complex_i_evidence.append(f"KEGG:{acc} ({name})")
            if acc in COMPLEX_I_PFAMS:
                complex_i_evidence.append(f"PFAM:{acc} ({name})")
            for ci_name in COMPLEX_I_NAMES:
                if ci_name in name_lower:
                    complex_i_evidence.append(f"name:{name}")
                    break

            # Check hydrogenase markers
            if acc in HYF_KOS:
                hydrogenase_evidence.append(f"Hyf KEGG:{acc} ({name})")
            if acc in HYC_KOS:
                hydrogenase_evidence.append(f"Hyc KEGG:{acc} ({name})")
            if acc in MATURATION_KOS:
                hydrogenase_evidence.append(f"maturation KEGG:{acc} ({name})")
            if acc in HYDROGENASE_PFAMS:
                hydrogenase_evidence.append(f"PFAM:{acc} ({name})")
            for ech_name in ECH_NAMES:
                if ech_name in name_lower:
                    hydrogenase_evidence.append(f"Ech name:{name}")
                    break

    # Decision logic
    n_ci = len(set(complex_i_evidence))  # deduplicate
    n_hyd = len(set(hydrogenase_evidence))

    details = {
        'protein_id': protein_id,
        'complex_i_markers': list(set(complex_i_evidence)),
        'hydrogenase_markers': list(set(hydrogenase_evidence)),
        'n_complex_i': n_ci,
        'n_hydrogenase': n_hyd,
    }

    if n_ci >= 2 and n_hyd == 0:
        return 'rejected', f"Complex I ({n_ci} markers, 0 hydrogenase)", details
    if n_ci >= 2 and n_hyd >= 1:
        # Both present -- ambiguous, but Complex I evidence is stronger
        return 'ambiguous', f"Mixed ({n_ci} Complex I, {n_hyd} hydrogenase)", details
    if n_hyd >= 1 and n_ci == 0:
        return 'confirmed', f"Hydrogenase operon ({n_hyd} markers, 0 Complex I)", details
    if n_ci == 1 and n_hyd == 0:
        return 'rejected', f"Weak Complex I (1 marker, 0 hydrogenase) -- conservative reject", details
    # Neither evidence
    return 'rejected', f"No operon context (0 markers) -- conservative reject", details


# Run validation on ALL Group 4 hits
results = {'confirmed': [], 'rejected': [], 'ambiguous': []}

for pid in group4_pids:
    verdict, summary, details = validate_group4_neighborhood(pid)
    details['verdict'] = verdict
    details['summary'] = summary
    results[verdict].append(details)

print(f"\n=== GROUP 4 VALIDATION RESULTS ===")
print(f"Confirmed (genuine hydrogenase operon): {len(results['confirmed'])}")
print(f"Rejected (Complex I false positive):    {len(results['rejected'])}")
print(f"Ambiguous (mixed evidence):             {len(results['ambiguous'])}")

if len(group4_pids) > 0:
    fp_rate = len(results['rejected']) / len(group4_pids) * 100
    print(f"False positive rate: {fp_rate:.0f}%")
```

### Step 4: Classify Confirmed Hits by Subgroup and Physiology

For confirmed Group 4 hits, determine the specific subcomplex (Ech, Hyf, Hyc, Mbh, Eha) and physiological role.

```python
def classify_group4_subtype(details):
    """Classify a confirmed Group 4 hit by subcomplex.

    Returns: (subtype, physiological_role)
    """
    markers = ' '.join(details.get('hydrogenase_markers', []))
    markers_lower = markers.lower()

    if any('hyf' in m.lower() for m in details['hydrogenase_markers']):
        return 'Group 4f (Hyf)', 'formate hydrogenlyase / fermentative H2 evolution'
    if any('hyc' in m.lower() for m in details['hydrogenase_markers']):
        return 'Group 4f (Hyc)', 'formate hydrogenlyase complex'
    if any('ech' in m.lower() for m in details['hydrogenase_markers']):
        return 'Group 4g (Ech)', 'energy-conserving, ferredoxin-dependent ion translocation'
    if any('mbh' in m.lower() for m in details['hydrogenase_markers']):
        return 'Group 4e (Mbh)', 'membrane-bound hydrogenase, archaeal energy conservation'
    if any('eha' in m.lower() for m in details['hydrogenase_markers']):
        return 'Group 4d (Eha)', 'energy-converting, methanogen-associated'
    if any('maturation' in m.lower() for m in details['hydrogenase_markers']):
        return 'Group 4 (unresolved)', 'confirmed by maturation factors, subgroup indeterminate'
    return 'Group 4 (unresolved)', 'confirmed by operon context, subgroup indeterminate'

for hit in results['confirmed']:
    subtype, physiology = classify_group4_subtype(hit)
    hit['subtype'] = subtype
    hit['physiology'] = physiology
    print(f"  {hit['protein_id']}: {subtype} -- {physiology}")
```

### Step 5: Spot-Check Groups 1-3 Without PF00374

```python
# Spot-check uncorroborated Groups 1-3 (pick up to 5)
import random
spot_check = uncorroborated_g13_pids[:5] if len(uncorroborated_g13_pids) <= 5 else random.sample(uncorroborated_g13_pids, 5)

for pid in spot_check:
    verdict, summary, details = validate_group4_neighborhood(pid, window=6)
    # For Groups 1-3, we are checking if they have any supporting context
    print(f"  {pid}: {verdict} -- {summary}")
    # If most come back 'rejected', note that some Group 1-3 calls may also be FPs
```

### Step 6: Determine Physiological Context for Confirmed Hits

**Ech is NOT primarily a "hydrogen metabolism" enzyme.** It is a membrane-bound, multisubunit [NiFe]-hydrogenase complex (EchA-F) that couples reversible H2/H+ chemistry to ion translocation. It is functionally an energy-coupling device -- an evolutionary intermediate in the Complex I superfamily.

Determine what the organism is doing with each confirmed hydrogenase:

```python
# For each confirmed Group 4, check genome-level metabolic context
for hit in results['confirmed']:
    pid = hit['protein_id']
    bin_id = b.store.execute("""
        SELECT bin_id FROM proteins WHERE protein_id = ?
    """, [pid])[0][0]

    # Check for methanogenesis markers
    mcr = b.store.execute("""
        SELECT COUNT(DISTINCT pp.protein_id) FROM protein_predicates pp
        JOIN proteins p ON pp.protein_id = p.protein_id
        WHERE p.bin_id = ? AND 'methanogenesis' = ANY(pp.predicates)
    """, [bin_id])[0][0]

    # Check for Wood-Ljungdahl markers
    wl = b.store.execute("""
        SELECT COUNT(DISTINCT pp.protein_id) FROM protein_predicates pp
        JOIN proteins p ON pp.protein_id = p.protein_id
        WHERE p.bin_id = ? AND 'wood_ljungdahl' = ANY(pp.predicates)
    """, [bin_id])[0][0]

    # Check for fermentation markers
    ferm = b.store.execute("""
        SELECT COUNT(DISTINCT pp.protein_id) FROM protein_predicates pp
        JOIN proteins p ON pp.protein_id = p.protein_id
        WHERE p.bin_id = ? AND 'fermentation' = ANY(pp.predicates)
    """, [bin_id])[0][0]

    context = []
    if mcr > 0: context.append("methanogenesis")
    if wl > 0: context.append("Wood-Ljungdahl")
    if ferm > 0: context.append("fermentation")

    hit['metabolic_context'] = context or ['unknown']
    hit['bin_id'] = bin_id
    print(f"  {pid} ({bin_id}): {', '.join(hit['metabolic_context'])}")
```

### Step 7: Generate Neighborhood Figures

> **Read `.claude/skills/visualize.md` before generating figures.** For any figure going into
> a report, use `plot_locus_multisource.py` (multi-source annotations, publication styling).
> Use `b.visualize_neighborhood()` for quick exploratory checks only.

```python
from pathlib import Path

FIGURES_DIR = Path("data/DATASET/figures")
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Generate figures for ALL confirmed Group 4 hits (mandatory -- these are the interesting ones)
for hit in results['confirmed']:
    pid = hit['protein_id']
    label = pid.replace('.', '_').replace(':', '_')
    # Use plot_locus_multisource.py for publication figures:
    # python scripts/plot_locus_multisource.py \
    #     --db data/DATASET/sharur.duckdb \
    #     --protein {pid} --window 10 \
    #     --output figures/hydrogenase_{label}.png \
    #     --title "Confirmed {hit['subtype']}" \
    #     --subtitle "{hit['physiology']}"

    # For quick exploration:
    b.visualize_neighborhood(
        pid, window=10,
        output_path=str(FIGURES_DIR / f"hydrogenase_{label}.png")
    )

# Also generate 1-2 representative REJECTED loci to show the Complex I context
for hit in results['rejected'][:2]:
    pid = hit['protein_id']
    label = pid.replace('.', '_').replace(':', '_')
    b.visualize_neighborhood(
        pid, window=10,
        output_path=str(FIGURES_DIR / f"complex_i_fp_{label}.png")
    )
```

### Step 8: Compile Final Inventory and Write Findings

```python
import json
from pathlib import Path

# Build the validated inventory
n_genomes = b.store.execute("SELECT COUNT(DISTINCT bin_id) FROM proteins")[0][0]

# Count genomes with confirmed hydrogenases
confirmed_pids = [h['protein_id'] for h in results['confirmed']]
if confirmed_pids:
    placeholders = ','.join(['?'] * len(confirmed_pids))
    g4_genomes = b.store.execute(f"""
        SELECT COUNT(DISTINCT bin_id) FROM proteins
        WHERE protein_id IN ({placeholders})
    """, confirmed_pids)[0][0]
else:
    g4_genomes = 0

g13_genomes = b.store.execute("""
    SELECT COUNT(DISTINCT p.bin_id)
    FROM proteins p
    JOIN protein_predicates pp ON p.protein_id = pp.protein_id
    WHERE 'nife_group1' = ANY(pp.predicates)
       OR 'nife_group2' = ANY(pp.predicates)
       OR 'nife_group3' = ANY(pp.predicates)
""")[0][0]

# Build summary
summary = {
    "groups_1_3": {
        "total_proteins": len(validated_pids),
        "n_genomes": g13_genomes,
        "validation": "PF00374 corroborated",
    },
    "group_4": {
        "raw_hyddb_calls": len(group4_pids),
        "confirmed": len(results['confirmed']),
        "rejected_complex_i": len(results['rejected']),
        "ambiguous": len(results['ambiguous']),
        "false_positive_rate": f"{fp_rate:.0f}%" if len(group4_pids) > 0 else "N/A",
        "confirmed_genomes": g4_genomes,
    },
    "fefe": {
        "total_proteins": len(fefe_validated_pids),
        "validation": "PF02906/PF02256 corroborated",
    },
}

print(json.dumps(summary, indent=2))
```

---

## Neighborhood Marker Lookup Tables

### Complex I Rejection Markers

If >=2 of these are found in the +-6 gene neighborhood, the Group 4 call is a Complex I false positive.

| Source | Identifiers | Description |
|--------|-------------|-------------|
| KEGG | K00330-K00343 | NuoA through NuoN (NADH:ubiquinone oxidoreductase subunits) |
| PFAM | PF00329 | Complex1_30kDa |
| PFAM | PF00346 | Complex1_49kDa (NuoD homolog -- the reason HydDB hits it) |
| PFAM | PF01058 | Oxidored_q6 (NADH:ubiquinone oxidoreductase) |
| PFAM | PF00361 | Proton_antipo_M (proton antiporter, membrane domain) |
| PFAM | PF01059 | Oxidored_q5_N |
| PFAM | PF00507 | Oxidored_q4 (24kDa subunit) |
| Name | *nadh dehydrogenase*, *nadh-ubiquinone*, *nuo* | Case-insensitive name matches |

### Hydrogenase Confirmation Markers

| Subgroup | Source | Identifiers | Description |
|----------|--------|-------------|-------------|
| Hyf (4f) | KEGG | K12136-K12145 | Hydrogenase-4 components (hyfA-J) |
| Hyc (4f) | KEGG | K15828-K15833 | Formate hydrogenlyase (hycB-G) |
| Ech (4g) | Names | echA-F | Energy-converting hydrogenase subunits |
| Any NiFe | KEGG | K04651-K04656, K03605 | Maturation factors HypA-F + HycI endopeptidase |
| Any NiFe | PFAM | PF00374 | NiFeSe_Hases on a NEIGHBORING gene (not the hit itself) |
| Any | PFAM | PF14720 | Fer4_20 (4Fe-4S cluster in hydrogenase maturation) |
| FeFe | PFAM | PF02906, PF02256 | FeFe hydrogenase large subunit domains |

### Ech (Group 4g) Operon Structure

The complete Ech complex has 6 subunits. The catalytic subunit (EchB) is the one HydDB hits.

| Subunit | Function | Homology |
|---------|----------|----------|
| EchA | Membrane antiporter-like | NuoL/M/N homolog |
| EchB | Catalytic [NiFe] subunit | NuoD homolog (this is why HydDB hits Complex I) |
| EchC | FeS electron transfer | NuoI homolog |
| EchD | Bridge subunit | NuoH homolog |
| EchE | Membrane antiporter-like | NuoL/M/N homolog |
| EchF | Membrane antiporter-like | NuoL/M/N homolog |

**If you see NuoA-N context instead of EchA-F context, it is Complex I. Reject.**

---

## Operon Validation Requirements by Subgroup

### Group 4g (Ech)

Requires echA-F operon. Subunits: EchA/E/F (membrane antiporter-like, homologous to NuoL/M/N), EchB (catalytic, homologous to NuoD), EchC (FeS electron transfer), EchD (bridge). If you see NuoA-N context, it is Complex I. Reject.

### Group 4f (Hyf/Hyc)

Look for hydrogenase-4 components (K12136-K12145 for Hyf; K15828-K15833 for Hyc). HydG maturation regulator nearby is confirmatory.

### Group 4e (Mbh)

Membrane-bound hydrogenase from Pyrococcus-like archaea. Rare in bacteria. If bacterial and in NuoA-N context, it is Complex I. Reject.

### Group 4b/4c/4d

Various energy-converting hydrogenases. All require operon validation. Same rejection criteria as above.

### Groups 1-3 (PF00374-corroborated)

PF00374 corroboration makes these high-confidence. Still check neighborhood for functional context:
- **Group 1**: Uptake hydrogenase -- expect respiratory chain proximity (cytochromes, quinone reductase)
- **Group 2**: Sensory/regulatory -- expect regulatory genes, chemotaxis context
- **Group 3**: Cofactor-coupled bidirectional -- expect NAD/F420/NADP-binding subunits

---

## Known False Positive Patterns

These are documented patterns that have caused incorrect hydrogenase claims in real analyses.

### Pattern 1: NuoD as "Group 4g Ech" (most common)

**Frequency:** 90%+ of Group 4 false positives.
**Mechanism:** HydDB catalytic subunit HMM matches Complex I subunit D (NuoD, K00333) because they share an ancestral fold.
**How to catch:** Check +-6 genes for nuoA-N (K00330-K00343), Complex1_30kDa (PF00329), Complex1_49kDa (PF00346).
**Typical appearance:** The "Group 4" hit sits in the middle of a 10-14 gene nuo operon.

### Pattern 2: Ni_hydr_CYTB (PF01292) as Hydrogenase

**What it is:** Cytochrome b subunit of NiFe hydrogenases OR respiratory complexes.
**Problem:** PF01292 alone does not indicate a hydrogenase. It is part of the same superfamily.
**Rule:** PF01292 without PF00374 or HydDB corroboration is NOT a hydrogenase.

### Pattern 3: Orphan Catalytic Subunits

**What it is:** A single gene annotated as Ech-like without any operon context.
**Problem:** Could be Mbh, Eha, or another Group 4 subcomplex with different physiology. Could also be a gene fragment or translocated domain.
**Rule:** Orphan catalytic subunits without operon context are classified as "ambiguous" and excluded from counts.

### Pattern 4: Proton_antipo_M (PF00361) Confusion

**What it is:** Proton antiporter membrane domain shared by Complex I (NuoL/M/N), Ech (EchA/E/F), and Mrp antiporters.
**Problem:** PF00361 alone is not diagnostic for either Complex I or Ech.
**Rule:** PF00361 requires additional context -- if with nuoB/C/D/E (K00331-K00334), it is Complex I. If with HydDB-annotated catalytic subunit AND maturation genes, it may be Ech.

---

## Language Guidance

### Correct Terminology

- DO say: "HydDB classifies this as Group 4 [NiFe]-hydrogenase, specifically in the Ech subgroup"
- DO say: "membrane-bound, energy-conserving, ferredoxin-dependent [NiFe]-hydrogenase"
- DO say: "Group 4 hydrogenases are part of the Complex I superfamily; their primary role is energy conservation via ion translocation"
- DO say: "N NiFe hydrogenases validated by PF00374 (Groups 1-3), plus M confirmed by neighborhood curation, excluding P Complex I false positives"

### Incorrect or Misleading Terminology

- DO NOT say: "a hydrogenase that makes H2" full stop -- this undersells the bioenergetic function
- DO NOT say: finding Ech subunits means the organism "uses hydrogen as an energy source" -- Ech in methanogens primarily produces H2 while pumping ions, or uses ion motive force for thermodynamically unfavorable electron transfer. Directionality is context-dependent.
- DO NOT say: "22 Group 4 hydrogenases detected" without validation -- this is reporting raw HydDB calls, which are likely 90%+ false positives
- DO NOT say: "lacks hydrogenases" for a MAG -- say "no hydrogenases were detected (N contigs)"

### What Ech Actually Is

Ech is a membrane-bound, multisubunit [NiFe]-hydrogenase complex (EchA-F) that couples reversible H2/H+ chemistry to ion translocation. It is functionally an **energy-coupling device**, not primarily a hydrogen-metabolizing enzyme. It is an evolutionary intermediate in the Complex I superfamily -- Ech, Mbh, and respiratory Complex I all share a common ancestor.

Individual Ech subunits (especially membrane ones like EchA/E/F with PF00361) might belong to Complex I instead. You need the full echABCDEF operon to call Ech.

Orphan catalytic subunits flagged as Ech might belong to Mbh, Eha, or other Group 4 subcomplexes with different physiologies. Be specific about what you have evidence for.

### Claim Escalation for Hydrogenases

| Claim Level | Evidence Required | Example |
|------------|-------------------|---------|
| "HydDB classifies N proteins as Group 4" | Raw HydDB output | Inventory step only |
| "N confirmed Group 4 hydrogenases" | Neighborhood validation (this protocol) | After Step 3 |
| "Ech energy-conserving hydrogenase" | echA-F operon structure validated | After Step 4 |
| "Organism uses Ech for ion translocation" | Operon + metabolic context (methanogen? syntrophy?) | After Step 6 |

---

## Output Specification

### Required Findings (append to findings.jsonl)

Each hydrogenase analysis MUST produce at minimum:

**Finding 1: Validated Hydrogenase Inventory**

```jsonl
{
  "id": "survey-NNN or E/D-NNN",
  "title": "Hydrogenase inventory: X Groups 1-3 (PF00374-validated), Y Group 4 (neighborhood-confirmed), Z FeFe in N genomes; P/Q Group 4 calls rejected as Complex I (FP rate: R%)",
  "category": "energy_metabolism",
  "description": "...",
  "evidence": "PF00374 corroboration for Groups 1-3; neighborhood curation with Complex I rejection markers (nuoA-N / PF00329 / PF00346 / PF01058 / PF00361) for Group 4; FP rate: P/Q (R%)",
  "n_genomes": N,
  "provenance": {
    "query": "neighborhood validation of N Group 4 HydDB calls",
    "raw_result": "X confirmed, Y rejected, Z ambiguous",
    "accession_verified": "PF00374 = NiFeSe_Hases (Groups 1-3); PF00346 = Complex1_49kDa (rejection marker)",
    "interpretation": "..."
  },
  "figures": ["figures/hydrogenase_*.png"],
  "phase": "survey"
}
```

**Finding 2 (if Group 4 confirmed): Specific Subgroup Finding**

One finding per confirmed subgroup type (Ech, Hyf, Hyc, Mbh) with:
- Subunit composition from operon analysis
- Physiological role in metabolic context
- Neighborhood figure showing the operon

### Required Report Language

The final report section MUST include:
1. Raw HydDB call count AND validated count (never just one)
2. False positive rate for Group 4
3. Validation method description
4. Subgroup-level breakdown
5. Metabolic context for confirmed Group 4 hits

**Example paragraph:**

> HydDB classified 145 proteins as [NiFe]-hydrogenases across 42 genomes. Of these, 97 belong to Groups 1-3 and are corroborated by the NiFeSe_Hases domain (PF00374): 45 Group 1 (respiratory uptake), 12 Group 2 (sensory), and 40 Group 3 (bidirectional, cofactor-coupled). The remaining 48 Group 4 calls were subjected to mandatory neighborhood validation: 6 were confirmed as genuine energy-converting hydrogenases (4 Hyf/hydrogenase-4, 2 Ech) based on conserved operon context, while 40 were rejected as Complex I (NADH:ubiquinone oxidoreductase) NuoD false positives (83% FP rate) and 2 were ambiguous. The confirmed Ech-type hydrogenases, present in 2 genomes, are membrane-bound, ferredoxin-dependent [NiFe]-hydrogenases that couple reversible H2/H+ chemistry to ion translocation for energy conservation.

---

## Hypothesis Tracking and Provenance

```python
# Log the validation as a provenance chain
e1 = b.log_provenance(
    "HydDB raw inventory",
    f"{nife_all} NiFe + {fefe_all} FeFe raw HydDB calls"
)
e2 = b.log_provenance(
    "Group 4 neighborhood validation",
    f"{len(results['confirmed'])} confirmed, {len(results['rejected'])} rejected ({fp_rate:.0f}% FP)",
    parent_ids=[e1.entry_id]
)

# Propose hypotheses based on results
if results['confirmed']:
    h = b.propose_hypothesis(
        "Group 4 hydrogenases serve energy-conserving function in syntrophic/fermentative niche"
    )
    b.add_evidence(
        h.hypothesis_id,
        "Neighborhood validation",
        f"{len(results['confirmed'])} confirmed Group 4 in metabolic context: {set(c.get('metabolic_context', ['unknown'])[0] for c in results['confirmed'])}",
        True, 0.7
    )

# Review all hypotheses
print(b.hypothesis_summary())

# Render provenance DAG
b.render_provenance(
    title="Hydrogenase Validation",
    output_path="figures/hydrogenase_provenance.mermaid"
)
```

Hypotheses persist across sessions -- `b.resume()` shows active hypotheses automatically.

---

## NiFe Subgroup Reference

| NiFe Group | Function | Key Predicates | PF00374? |
|------------|----------|----------------|----------|
| 1a-1l | Respiratory uptake | `nife_group1`, `uptake_hydrogenase` | Yes |
| 2a-2e | Cytoplasmic H2 sensors | `nife_group2`, `h2_sensor` | Yes |
| 3a-3d | Bidirectional, cofactor-coupled | `nife_group3`, `bidirectional_hydrogenase` | Yes |
| 4a-4i | Energy-conserving, H2-evolving | `nife_group4`, `mbh_hydrogenase`, `ech_hydrogenase` | **NEVER** |

| FeFe Group | Function | Key Predicates |
|------------|----------|----------------|
| A1-A4 | Monomeric fermentative | `fefe_groupA` |
| B | Electron-bifurcating | `fefe_groupB`, `bifurcating_hydrogenase` |
| C1-C3 | Sensory/regulatory | `fefe_groupC` |

## Key PFAM Domains

| Domain | Accession | Meaning | Diagnostic? |
|--------|-----------|---------|-------------|
| NiFeSe_Hases | PF00374 | NiFe large subunit (Groups 1-3 only) | Yes for Groups 1-3 |
| Fe_hyd_lg_C | PF02906 | FeFe hydrogenase large subunit | Yes for FeFe |
| Complex1_49kDa | PF00346 | NADH dehydrogenase subunit D | **Complex I marker** |
| Complex1_30kDa | PF00329 | NADH dehydrogenase 30kDa subunit | **Complex I marker** |
| Oxidored_q6 | PF01058 | NADH:ubiquinone oxidoreductase | **Complex I marker** |
| Proton_antipo_M | PF00361 | Proton antiporter membrane domain | Shared: Complex I AND Ech |
| Ni_hydr_CYTB | PF01292 | Cytochrome b -- NOT a hydrogenase | Not diagnostic alone |
| Fer4_20 | PF14720 | 4Fe-4S in hydrogenase maturation | Confirmatory |
