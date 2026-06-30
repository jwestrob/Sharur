# Shared Validation Protocols

**Required reading for ALL Sharur analysis skills.** These protocols prevent known error classes that have caused false findings and manuscript retractions.

---

## 1. Accession Verification

**NEVER assume a PFAM/KEGG accession encodes a specific function from memory.** Always verify against the database `name` field before reporting.

This rule exists because PF04055 was falsely claimed as "benzoyl-CoA reductase" (actually "Radical_SAM") — the error propagated through findings, synthesis, and into a manuscript before being caught. One SQL query would have prevented it.

**Before reporting ANY PFAM/KEGG-based functional claim:**

```python
# STEP 1: Verify the accession name in the database
name_check = b.store.execute("""
    SELECT DISTINCT name FROM annotations WHERE accession = 'PF04055'
""")
print(f"PF04055 = {name_check}")  # → [('Radical_SAM',)]
# If the name doesn't match your assumption, STOP. Your claim is wrong.
```

**Red flags requiring extra verification:**
- Accession you haven't verified in this session
- Single PFAM domain cited as proof of a complete pathway
- "Universal" claims (>90% of genomes) — especially if the domain is a large superfamily
- Protein counts >500 from a single accession — large superfamilies are usually generic

**The cost of checking is 1 SQL query. The cost of not checking was a retracted finding.**

---

## 2. COUNT(DISTINCT protein_id) Rule

**ALWAYS use `COUNT(DISTINCT protein_id)` when reporting protein counts.** Repeat domains (WD40, TPR, FG-GAP, Ig-like, etc.) produce multiple annotation rows per protein. `COUNT(*)` on the annotations table gives annotation hits, not protein counts — these can differ by 5-10x.

```python
# WRONG — counts annotation rows (864 for WD40 in giant proteins)
b.store.execute("SELECT COUNT(*) FROM annotations WHERE accession = 'PF00400'")

# RIGHT — counts unique proteins (352 for WD40)
b.store.execute("SELECT COUNT(DISTINCT protein_id) FROM annotations WHERE accession = 'PF00400'")
```

When reporting counts, always state what you're counting: "352 proteins with WD40 domains" or "864 WD40 domain instances across 352 proteins."

---

## 3. Context-First Protocol

**The domain tells you the fold. The neighbors tell you the function.**

Annotation databases assign labels based on structural fold similarity, not biological function. When a PFAM HMM or KOfam profile matches an entire enzyme superfamily, hit counts are high but the name is misleading.

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
co_annots = b.store.execute("""
    SELECT a2.source, a2.accession, a2.name, COUNT(DISTINCT a1.protein_id) as n
    FROM annotations a1
    JOIN annotations a2 ON a1.protein_id = a2.protein_id
    WHERE a1.accession = 'K23108'   -- your claimed KO
      AND a2.accession != 'K23108'  -- exclude self
    GROUP BY a2.source, a2.accession, a2.name
    ORDER BY n DESC LIMIT 10
""")
# If top co-annotations are from a DIFFERENT enzyme family → superfamily cross-hit
```

### Neighborhood Validation

For protein-level functional claims that imply a pathway or system:

```python
result = b.get_neighborhood(protein_id, window=5, all_annotations=True)
# Shows PFAM, KEGG, HydDB, DefenseFinder, VOGdb, CAZy per gene
```

**What to look for:**
- Neighbors consistent with claimed function → Claim supported
- Neighbors suggest a different system → Revise interpretation
- No informative neighbors → State uncertainty
- Protein at contig edge → Note potential fragmentation

### Reading V2 Semantic State (preferred for facet-level claims)

Stage 07 now writes typed semantic atoms to `semantic_atoms` and aggregated
`semantic_state` rows in `sharur.duckdb`. **When making a functional claim,
read the atoms — not the raw annotation row.**

```python
# What does V2 actually claim about this protein?
state = b.get_semantic_state("CoronaMine_..._IL_777_6")
# state.activities = {"hydrogenase": {"relation": "implies", ...}}
# state.roles      = {"defense_system": {"relation": "supports", ...}}

# Pull every atom for finer detail (which source DB, which accession)
atoms = b.get_atoms("CoronaMine_..._IL_777_6")
```

Atom relations are the leverage that closes the iconic-member trap:

- `implies` — strong evidence (HydDB classification, KEGG ortholog). Cite this.
- `supports` — moderate (PFAM domain in family). Use as one line of evidence.
- `flags`   — weak / superfamily / pre-system-validation. Do NOT report as function.
- `excludes` — counter-evidence (Complex I subunit excludes hydrogenase).
- `unresolved` — no curated mapping; the accession needs human triage.

**Rule:** never report `(activity, X)` based on a single `flags` atom. If the
only support for "this protein does X" is a flags-level atom, write
"fold consistent with X family — substrate undetermined."

The compat layer materializes `protein_predicates` from V2 so
`b.search_by_predicates(has=[...])` continues to work, but it loses the relation
metadata. For high-stakes claims (manuscript-bound findings, novelty >= 2),
go through the V2 facade.

### Batch Neighborhood Validation (preferred for >5 candidates)

When you have a list of hits to validate (defense system candidates, hydrogenase candidates, CAZy guild members), **do not loop calling `get_neighborhood` one at a time** — that's the open-coded pattern defense/prophage/hydrogenase used to share. Call `batch_context_validate` instead: one SQL pass for N hits, returns a typed verdict per hit.

```python
from sharur.operators import batch_context_validate, summarize_verdicts

# Example: validate Group 4 NiFe hydrogenase candidates as Mbh-type
candidates = [r[0] for r in b.store.execute(
    "SELECT DISTINCT protein_id FROM annotations WHERE name LIKE '%MbhD%'"
)]
verdicts = batch_context_validate(
    b.store, candidates,
    required_any=["Mbh", "MnhA", "MnhD"],          # any of these in window = real
    blocklist=["NuoD", "NuoH", "NuoL", "NuoM"],    # Complex I markers = FP
    window=8,
)
print(summarize_verdicts(verdicts))
# → {"total": 132, "counts": {"validated": 63, "rejected_no_context": 65, "rejected_superfamily_match": 4}, ...}

# Iterate validated set for follow-up
real = [v.protein_id for v in verdicts if v.verdict == "validated"]
```

Verdicts: `validated`, `rejected_no_context`, `rejected_superfamily_match`, `ambiguous`, `protein_not_found`.

**When to use which:**
- `required_any` — "at least one of these markers nearby" (most common). E.g., any Mbh subunit, any Cas accessory.
- `required_all` — "all of these markers nearby". E.g., a complete cassette where each subunit is mandatory.
- `blocklist` — superfamily / Complex I-style false-positive tokens. Always include when a known superfamily collision exists.

The function does NOT consult the `defense_systems` table. For DefenseFinder hits, query that table directly (per CLAUDE.md: 80%+ FP rate on raw HMM hits).

### Claim Escalation Ladder

| Claim Level | Evidence Required |
|------------|-------------------|
| "Contains domain X" | Accession name verification only |
| "Functions as enzyme Y" | Name verification + co-annotation check |
| "Genome encodes pathway Z" | Multiple pathway-specific markers + co-localization |
| "Phylum universally performs W" | All above + conservation across genomes |

---

## 4. KEGG REST API for Pathway Context

When you find a KEGG KO and want to understand its pathway context:

```python
import subprocess

# What KEGG modules (minimum functional units) use this KO?
result = subprocess.run(
    ["curl", "-s", "https://rest.kegg.jp/link/module/ko:K23108"],
    capture_output=True, text=True
)

# Get the full module definition (what's required for function)
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

---

## 5. MAG Quality and Absence Claims

MAGs are inherently incomplete. A missing gene does NOT mean the organism lacks it.

| Contigs | Fragmentation | How to Interpret Absence |
|---------|---------------|--------------------------|
| <50 | Low | Reasonably reliable |
| 50-200 | Moderate | Include caveats |
| >200 | High | "Not detected" only |
| genes/contig <5 | Very high | Many genes likely missing |

**Language:**
- "No hydrogenases were detected in this MAG (N contigs)" — NOT "Genome X lacks hydrogenases."
- Before saying "A has X but B doesn't", verify B isn't just more fragmented.

```python
# Check fragmentation before making absence claims
contig_stats = b.store.execute("""
    SELECT bin_id, COUNT(DISTINCT contig_id) as n_contigs,
           ROUND(COUNT(*)::FLOAT / COUNT(DISTINCT contig_id), 1) as genes_per_contig
    FROM proteins GROUP BY bin_id ORDER BY n_contigs DESC
""")
```

---

## 6. Hydrogenase Curation

HydDB classifies proteins by structural similarity to known hydrogenases, but NADH dehydrogenase (Complex I) shares the same [NiFe] binding site fold and produces ~44% false positives for NiFe calls.

The pipeline (`classify_hydrogenases.py`) uses a PF00374 filter that validates Groups 1-3 but **systematically rejects all Group 4 NiFe hydrogenases** (Hyf/Hyc/Mbh/Ech) because they diverged too far.

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
    # Complex I evidence (→ reject):
    #   KEGG: K00330-K00343 (nuoA-N)
    # Ambiguous (both or neither): reject conservatively
```

**In reports, state the curation method:** "97 NiFe hydrogenases validated by PF00374 (Groups 1-3), plus 12 rescued by neighborhood curation (8 Group 4f + 4 Group 4e), excluding 66 Complex I false positives."

---

## 7. Findings Provenance

Every finding that makes a functional claim requires a `provenance` dict with:
- `query`: The exact code/SQL that produced the evidence
- `raw_result`: The literal output (count, list, accession name, etc.)
- `accession_verified`: For PFAM/KEGG claims, the NAME field from the database
- `interpretation`: The human-readable claim derived from the raw result

```python
provenance = {
    "query": "SELECT COUNT(*), name FROM annotations WHERE accession = 'PF04055' GROUP BY name",
    "raw_result": [(1790, "Radical_SAM")],
    "accession_verified": "PF04055 = 'Radical_SAM' — generic superfamily, NOT benzoyl-CoA reductase",
    "interpretation": "1,790 Radical SAM proteins (normal for anaerobic bacteria, not pathway-specific)"
}
```

This creates an audit trail that reviewers can independently verify.

---

## 8. Cas12f / TnpB / OrfB — the non-Cas RuvC trap (HARD RULE)

**Do not report Cas12f as a finding.** Cas12f, TnpB, and IS200/IS605 OrfB share a RuvC-like nuclease fold. PFAM `Cas12f`/`Cas12f1` matches all three with high false-positive rates. In metagenomes — especially archaeal — the overwhelming majority of "Cas12f" PFAM hits are non-Cas TnpB or IS-OrfB transposases. This trap has produced spurious "novel CRISPR effector" findings on multiple datasets despite existing guidance to discount it.

**Default behavior**: treat any PFAM `Cas12f*` hit as a TnpB-family RuvC homolog, NOT as a CRISPR effector. Inventory at PFAM-completeness level only; do not survey, explore, or characterize.

**The only exception** — a Cas12f finding may be elevated only when ALL of these hold:
1. The protein is **not** carrying IS200/IS605 OrfB co-annotation on the same protein record
2. The protein is **within 5 kb of a MinCED-detected CRISPR array** in the same bin
3. The locus carries **either Cas4-Cas1-Cas2 adaptation OR Tn7 (TnsB/C/TniQ) CAST machinery** within 20 genes
4. The locus is **not within 10 genes of a TnpA/Y1_Tnp transposase** (that combination is the diagnostic TnpB pairing)

A single failing condition disqualifies the finding. Drop it.

**Do not** open a Cas12f explore agent unless a survey-phase finding meets all four conditions for at least one locus.

**Reviewer_2 must reject** any finding that elevates "Cas12f" without all four conditions evidenced in the verification queries.

This rule supersedes any softer "watch for TnpB confusion" comments elsewhere in the skill set (`defense.md`, `manuscript_guide.md`).

---

## 9. Confidence in catalytic function from single-HMM annotations

Catalytic substrate identity is typically not well-resolved by a single PFAM domain hit alone. PFAM HMMs frequently match multiple protein families that share a fold but catalyze different reactions on different substrates — for instance, the same RubisCO large-subunit HMM scores across form I, form II, form III, and RuBisCO-like proteins, which span substantively different substrate ranges. KOFAM HMMs exhibit a comparable pattern at a different scale of divergence: a top-scoring KO assignment can sit just above many near-tied KO models that together indicate a fold superfamily rather than a confident substrate-specific call. The Sharur database surfaces only the best hit per protein, which removes the multi-hit evidence that would otherwise signal this ambiguity. Substrate-level claims from a single PFAM or KEGG annotation generally benefit from supporting context — neighboring partner subunits, pathway co-localization, or other contig-resident evidence — before the call is treated as confident.

---

## 10. KEGG / KOFAM scope of confidence

KEGG annotations are most reliable in canonical model-organism contexts. For non-model lineages — uncultured archaea, environmental anaerobes, novel phyla, divergent or basal clades — KOFAM HMM coverage has uneven gaps that vary by pathway and lineage, with the gaps often invisible at the time of analysis. In these contexts, default to PFAM-first scoring, cross-reference with KEGG where the KO is available, and read KEGG absence as uninformative rather than as evidence of pathway absence. KEGG remains the canonical pathway atlas — narrowing its scope of high-confidence inference for non-model lineages preserves its strengths without overstating its coverage.

---

## 11. Figure construction

Figures begin with the argument they make, not the data they show. Write the argument in one sentence before opening a plotting tool. The visualization type follows from the argument. After rendering, `Read` the PNG and verify the argument lands at a glance; iterate the layout until it does. The minimum protocol is render → Read → critique → adjust, repeated until the figure communicates clearly at the size humans will read it.

---

## 12. Lineage biology cross-reference

Before shipping a quantitative claim about a named lineage, verify the claim is consistent with published biology for that lineage. Where the claim contradicts standard pathway biology for a well-studied lineage — for instance, a documented producer scoring zero for the biosynthesis pathway, or an obligate anaerobe carrying a pathway whose markers require oxygen — web-search to confirm or to identify the methodological gap behind the apparent contradiction. The literature-validation step is cheap and catches the highest-impact errors.

---

## Self-Contained Finding Titles

Finding titles must include ALL qualifiers. A reader should understand the scope from the title alone:

- **Bad**: "Largest protein (5,461 aa) has zero annotations"
- **Good**: "Largest unannotated protein (5,461 aa, JAJVIK010000008.1_48) — candidate adhesin"

Qualifiers dropped during summarization have caused factual errors in manuscripts.
