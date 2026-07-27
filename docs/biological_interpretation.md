# Biological Interpretation Guide

**Load this doc when:** Making biological claims, interpreting annotations, validating functional assignments, or writing about MAG-derived findings.

## Known Annotation Confounds — CHECK THIS FIRST

**Keyed by what you OBSERVE, not by what you set out to study.** You will usually
meet these by running a query and getting rows back, with no prior reason to
suspect a trap. Scan this table against your result set *before* interpreting it.
Every row below has produced a wrong conclusion in a real session.

| If you observe | The trap | Required check before any claim |
|---|---|---|
| `source='hyddb'` / `hyddb_subgroup` = `nife_group4`, `energy_conserving_hydrogenase`, `ech_hydrogenase`, `h2_evolving` | Group 4 catalytic-subunit HMM is homologous to respiratory **Complex I NuoD/NqoD**. Most hits are Complex I, not hydrogenase. | Co-annotation: count `Complex1_49kDa` (PF00346) + `Complex1_30kDa` (PF00329) vs `NiFeSe_Hases` on the *same* proteins. Full protocol: `.claude/skills/hydrogenase.md` |
| `hyddb_subgroup` = `hyddb_needs_curation` | The caller itself is flagging low confidence. Often ~99% of hyddb rows in a taxon. | Report subgroup as OBSERVED/UNVERIFIED. Never emit a named hydrogenase class from a flagged row. |
| `RuBisCO_large` / `RuBisCO_large_N` | Form IV RuBisCO-like proteins (RLP) do methionine salvage, not carbon fixation. Large subunit alone proves nothing. | Require **PRK** (PF00485 / KEGG K00855) in the same bin. `RuBisCO_small` further discriminates form I (L8S8) from form II/IV. |
| `Ald_Xan_dh_C` (and molybdopterin oxidoreductase families) | Contains CoxL (form I CO dehydrogenase) **but also** xanthine dehydrogenase, aldehyde oxidoreductase, nicotinate dehydrogenase. | Never call CO oxidation from the family alone. Require the CoxL active-site motif (AYXCSFR) plus `coxM`/`coxS` co-localization. |
| Any accession averaging **>10 hits/genome**, or a claimed function in **>50% of genomes** | You are looking at a fold/superfamily, not a function. | Co-annotation + ±8-gene neighborhood on 3–5 examples (protocol below). |
| A trait "absent" because a name/description `LIKE` pattern returned zero | Naming artifact, not biology. `%glycogen%` misses the entire **GlgE** pathway (`Glyco_transf_5`, `GlgB_N`, `GlgP_C`, `GlgE_dom_N_S`, `CBM_48`). | Probe by accession family and pathway members, never by one English name. Report zeros only after ≥3 orthogonal patterns. |
| Any absence in a MAG | Fragmentation, not biology. | Check contig count + completeness first. See MAG Interpretation below. |
| Genes on **different contigs** of the same bin described as one system / operon / locus / cluster | A bin is a statistical grouping, not a chromosome. Two contigs may be different organisms, strain variants, or a misbin. Co-membership in a bin is **not** evidence of co-localization. | Every gene in the call must share one `contig_id`. If they do not, it is co-occurrence *in a bin* and must be labelled as such. See "Co-localization requires one contig" below. |

**Escalation rule:** an HMM row is an *observation*. A functional name requires the
co-annotation and neighborhood checks. A pathway claim requires multiple markers
co-localized. See "Claim escalation" below.

**Adding to this table:** when a confound costs you a wrong conclusion, add the row
here — keyed on the observable that misled you — and put the deep protocol in the
relevant `.claude/skills/*.md`. A fact filed only under its topic will not be found
by someone who does not yet know the topic applies.

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

### Co-localization requires one contig

A bin is a statistical grouping of contigs, not a chromosome. Contigs in one bin
may belong to different organisms, to strain variants, or to a misbin. Physical
adjacency is observable only *within* a contig.

> **A system, operon, locus, cluster, cassette or island may only be called from
> genes that share one `contig_id`.** Genes on different contigs of the same bin
> are co-occurrence in a bin — a weaker and different claim — and must be written
> that way, never as a locus.

Components split across contigs are equally consistent with one fragmented
system, two partial systems, or one system plus contamination. Nothing in the
data distinguishes them, so the linkage must be reported as unresolved.

Consequences for how work is set up:

- **Any co-localization caller must be scoped per contig, not per bin.** Tools
  that accept a whole genome as one ordered replicon will assemble systems
  across contig boundaries, especially under models that permit a system to be
  built from separate loci. Give the tool the contig as the unit of adjacency.
- **Any agent reading a multi-contig packet must be told** that adjacency in the
  packet is not adjacency in the genome.
- **Filter cross-contig calls as an assertion, not as the remedy.** Dropping them
  after a caller has chosen its best solution cannot undo a spurious assembly
  having outcompeted a correct one.

The cost of the rule is that a real system split by assembly fragmentation goes
uncalled. That is the correct trade: an unobservable linkage must not be
asserted. Report component-level evidence separately from system-level calls.

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
3. **Sequence motifs:** compute the CxxCH heme-binding motif locally and expose
   only stable protein IDs and match counts to the model

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
