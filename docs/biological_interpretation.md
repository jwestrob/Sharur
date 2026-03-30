# Biological Interpretation Guide

**Load this doc when:** Making biological claims, interpreting annotations, validating functional assignments, or writing about MAG-derived findings.

## Context-First Analysis Protocol

**The domain tells you the fold. The neighbors tell you the function.**

**Superfamily awareness applies to ALL HMM-based annotation sources** (PFAM, KEGG, DefenseFinder, VOGdb, CAZy). When any accession averages >10 hits per genome, it likely describes a protein fold, not a specific function. Even below that threshold, if a claimed function appears in >50% of genomes, ask whether that prevalence makes biological sense.

Before reporting any HMM-based functional claim, run two checks:
(1) **Co-annotation**: what other domains do these proteins carry?
(2) **Genome context**: pick 3-5 examples and examine the ±8 gene neighborhood.

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

## PFAM Function Verification

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

## MAG Interpretation

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

## Giant Protein Annotation Recovery

**Standard PFAM bitscore cutoffs are LENGTH-BIASED.** Giant proteins (>1000 aa) often show zero hits. E-values are NOT length-biased — use them instead.

```bash
# For proteins >1000 aa with 0 standard hits
hmmsearch --domE 1e-5 ~/.config/Astra/PFAM/Pfam-A.hmm giant_protein.faa
```

**When to apply:** >1000 aa with 0 hits, >2000 aa with <3 hits, before reporting any giant protein as "unannotated."

**E-value interpretation:** ≤1e-10 high confidence; 1e-10 to 1e-5 moderate; 1e-5 to 1e-3 weak.

**Common giant protein architectures:** Big_13/Big_3_3/Big_8 (adhesins), TPR/HEAT/WD40 (scaffolds), Beta_helix (autotransporters), ANK (signaling), Cadherin/FN3 (adhesion).

## Hydrogenase Classification

For full classification protocol, neighborhood validation markers, and false positive patterns, see `.claude/skills/hydrogenase.md`.

## Cytochrome Validation

**Always validate respiratory system claims, especially "lacks cytochromes."**

Three detection methods:
1. **Predicates:** `b.search_by_predicates(has=["cytochrome"])`
2. **Annotations:** `WHERE LOWER(name) LIKE '%cytochrome%'`
3. **Sequence motifs:** CxxCH heme-binding motif (`re.search(r'C..CH', seq)`)

If all three find nothing AND genome is high-quality (<50 contigs), absence is plausible.

## Scientific Rigor

**Under-promise, over-deliver.** Rigorous, conservative science is more impressive than hyperbolic claims.

**Forbidden language:** "confirms/proves/demonstrates" (unless truly definitive), "unprecedented/first ever/groundbreaking", "paradigm-shifting/revolutionary"

**Required language:** "suggests/indicates/supports", "consistent with/compatible with", "to our knowledge" (after verification), "provides evidence for"

| Evidence | Appropriate Language |
|----------|---------------------|
| Sequence annotation | "annotated as", "contains domain" |
| Structural homology | "structural similarity suggests" |
| Genomic context | "co-located with", "may indicate" |
| Literature | "similar to", "consistent with" |
| Experimental | "demonstrates", "confirms" (OK!) |

**Common errors:** Domain presence ≠ function proof. MAG absence ≠ biological absence. Single marker ≠ pathway presence. Transposase ≠ mobile element proof (could be Cas12f). "First in analysis" ≠ "first ever."

## Predicate System Design Principles

- **Gene-level vs Locus-level**: Components get gene-level tags; system calls require clustering
- **Tiered predicates**: Generic domain → specific system (e.g., `cas_domain` → `crispr_associated`)
- **Biological accuracy**: Ferredoxins are electron carriers, not substrate indicators
- **Pathway completeness**: Single markers don't prove pathway function
- **Carbon fixation**: RuBisCO without PRK is NOT Calvin cycle. `calvin_cycle` requires PRK.
- **Methanogenesis**: Only MCR complex triggers; H4MPT enzymes alone are insufficient
- **Predicate V2:** Typed semantic atoms with composite DSL. See `docs/predicates_v2.md`

### PFAM Mapping Scaling

Extension file: `data/reference/pfam_predicate_map.tsv`. Bulk expansions go in the TSV, curated core mappings stay in Python (`sharur/predicates/mappings/pfam_map.py`).
