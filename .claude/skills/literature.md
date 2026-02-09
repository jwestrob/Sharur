# Literature Research Skill

Systematic literature and database research to resolve functional ambiguity, interpret structural hits, understand domains, and provide biological context for metagenomic findings.

**CRITICAL: You are a leaf agent. DO NOT spawn sub-agents or use the Task tool.**

**Note:** This skill uses WebSearch and WebFetch only — no database writes. Can run in parallel with other skills.

> **Mandatory:** Follow the shared validation protocols in `_validation_protocols.md`.
> In particular: verify PFAM/KEGG accession names against the database before claiming function,
> use COUNT(DISTINCT protein_id) for protein counts, and apply Context-First protocol.

---

## Usage

```
/literature PF12345                    # Domain/DUF research
/literature foldseek AF-Q8ZZM4-F1     # Research a Foldseek hit
/literature kegg K23108               # KEGG module context
/literature organism Hinthialibacterota  # Lineage context
/literature defense CBASS             # Defense system research
/literature family "giant adhesin"    # Protein family research
/literature manuscript data/omni_production/MANUSCRIPT.md  # Verify & find all citations for manuscript
```

---

## Research Protocols

### Protocol A: Domain / DUF Research

**When to use:** A PFAM domain is found and its function is unclear, especially DUFs (Domains of Unknown Function) or large superfamilies where the name doesn't convey specific function.

**Search strategy:**

1. **InterPro lookup** — authoritative source for domain classification and clan membership:
   ```
   WebFetch("https://www.ebi.ac.uk/interpro/entry/pfam/PF{XXXXX}/",
            "What is this domain's function, what clan does it belong to, and has it been reclassified?")
   ```

2. **WebSearch for characterization:**
   ```
   WebSearch("[domain name] characterized function")
   WebSearch("[domain name] crystal structure mechanism")
   WebSearch("DUF[number] renamed reclassified")  # DUFs specifically
   ```

3. **Check clan membership** — if the domain belongs to a clan (superfamily), the clan tells you the fold but NOT the specific function. Note this in your output.

4. **Check if DUF has been reclassified** — many DUFs have been characterized and renamed since their original designation. InterPro will show the current status.

**Output fields:**
- `domain_id`: PFAM accession (e.g., PF12345)
- `domain_name`: Name from InterPro (may differ from database)
- `clan`: Superfamily clan if any (e.g., CL0023 = P-loop NTPase)
- `function_summary`: What is known about the domain's function
- `reclassified`: Whether a DUF has been renamed/characterized
- `key_references`: PMIDs or DOIs of characterization papers
- `confidence`: high (experimentally characterized) / medium (computational prediction) / low (uncharacterized DUF)

---

### Protocol B: Foldseek Hit Research

**When to use:** Structure search (Foldseek) returns PDB or AlphaFold hits that need functional interpretation. This is the primary way to assign function to unannotated proteins.

**ID parsing — extract researchable accessions from Foldseek target strings:**

| Foldseek format | Type | How to extract |
|----------------|------|----------------|
| `AF-Q8ZZM4-F1` | AlphaFold/UniProt | Strip `AF-` prefix and `-F1` suffix → `Q8ZZM4` |
| `AF-A0A1P8B9K7-F1` | AlphaFold/UniProt (long) | Same pattern → `A0A1P8B9K7` |
| `5fms_A` | PDB | First 4 characters → `5fms`, chain = `A` |
| `3gw6_B` | PDB | First 4 characters → `3gw6`, chain = `B` |

**Search strategy:**

1. **For PDB hits** — experimentally determined structures with published papers:
   ```
   WebFetch("https://www.rcsb.org/structure/{PDB_ID}",
            "What is this protein's function, organism, and biological context? What paper describes it?")
   ```

2. **For AlphaFold/UniProt hits** — predicted structures linked to UniProt annotations:
   ```
   WebFetch("https://www.uniprot.org/uniprotkb/{UNIPROT_ID}",
            "What is this protein's function, organism, gene name, and GO annotations?")
   ```

3. **Interpret TM-scores:**

   | TM-score | Interpretation |
   |----------|---------------|
   | > 0.7 | High confidence structural homology — function transfer likely valid |
   | 0.5–0.7 | Similar fold — function may be conserved but verify with other evidence |
   | 0.3–0.5 | Weak similarity — same fold superfamily but function may differ |
   | < 0.3 | Not meaningful — do not use for function prediction |

4. **Cross-validate** — if multiple Foldseek hits point to the same function, confidence increases. If hits point to different functions, the domain may be a structural fold shared across functional families.

**Output fields:**
- `query_protein`: Protein ID that was searched
- `target_id`: Foldseek target string
- `parsed_accession`: Extracted PDB/UniProt ID
- `target_function`: Function from PDB/UniProt
- `target_organism`: Source organism
- `tm_score`: Structural similarity score
- `function_prediction`: Inferred function for query protein
- `evidence_quality`: high (PDB hit, TM > 0.7) / medium (AF hit or TM 0.5-0.7) / low (TM < 0.5)

---

### Protocol C: KEGG Module Research

**When to use:** A KEGG KO is found and you need to understand its pathway context — what metabolic module it belongs to, what the complete module requires, and how much of the module is present in the dataset.

**Search strategy:**

1. **Find modules containing this KO:**
   ```bash
   curl -s "https://rest.kegg.jp/link/module/ko:{KO_ID}"
   # Returns: ko:K23108\tmd:M00551
   ```

2. **Get module definition (required components):**
   ```bash
   curl -s "https://rest.kegg.jp/get/md:{MODULE_ID}"
   # Returns: module name, definition (required KOs), class, pathway
   ```

3. **Get KO description:**
   ```bash
   curl -s "https://rest.kegg.jp/get/ko:{KO_ID}"
   # Returns: KO name, definition, linked pathways, genes
   ```

4. **Cross-reference with dataset** — check which other KOs from the module are present:
   ```python
   # After getting required KOs from module definition
   module_kos = ['K00001', 'K00002', 'K00003']  # from KEGG
   present = b.store.execute(f"""
       SELECT DISTINCT accession FROM annotations
       WHERE accession IN ({','.join(f"'{k}'" for k in module_kos)})
   """).fetchall()
   ```

**Caveats for divergent organisms:**
- KOfam profiles are built from characterized enzymes in model organisms
- Deeply branching lineages (Asgard archaea, CPR) may have divergent homologs below detection
- Report: "X/Y module steps detected; undetected steps may reflect sequence divergence"

**Output fields:**
- `ko_id`: Query KO
- `ko_name`: KO function name
- `modules`: List of KEGG modules containing this KO
- `module_completeness`: For each module, fraction of required KOs present in dataset
- `biological_role`: Summary of metabolic role
- `confidence`: high (>75% module complete) / medium (50-75%) / low (<50% or single marker)

---

### Protocol D: Organism / Lineage Context

**When to use:** Analyzing an understudied lineage (candidate phylum, novel class) and you need context about what's known — expected metabolism, ecology, related organisms, published MAG studies.

**Search strategy:**

1. **General lineage search:**
   ```
   WebSearch("[lineage name] metabolism ecology genome")
   WebSearch("[lineage name] metagenome-assembled genome MAG")
   ```

2. **GTDB lineage check** — for taxonomic context:
   ```
   WebSearch("[lineage name] GTDB taxonomy")
   ```

3. **Published MAG studies:**
   ```
   WebSearch("[lineage name] metabolic reconstruction")
   WebSearch("[lineage name] genomic potential")
   ```

4. **Related organisms** — if the lineage is poorly studied, look at the closest well-studied relative:
   ```
   WebSearch("[parent clade] characterized representative metabolism")
   ```

5. **Key review papers** — for major lineages:
   ```
   WebSearch("[lineage name] review nature microbiology")
   ```

**Output fields:**
- `lineage`: Taxonomic name and rank
- `known_lifestyle`: Autotroph/heterotroph, aerobe/anaerobe, symbiont/free-living
- `expected_capabilities`: Metabolic pathways previously reported in this lineage
- `notable_features`: Unusual biology (giant proteins, minimal genomes, etc.)
- `related_organisms`: Closest well-studied relatives
- `key_references`: PMIDs/DOIs of foundational studies
- `confidence`: high (well-studied lineage) / medium (a few MAG studies) / low (candidate phylum, minimal data)

---

### Protocol E: Defense System Research

**When to use:** Novel or rare defense systems are identified by DefenseFinder, PADLOC, or manual curation, and you need to understand mechanism, distribution, and biological significance.

**Search strategy:**

1. **System-specific search:**
   ```
   WebSearch("[system name] anti-phage defense mechanism")
   WebSearch("[system name] defense system bacteria archaea")
   ```

2. **Key review papers** — the field moves fast, prioritize recent work:
   ```
   WebSearch("[system name] Doron Sorek defense")       # Sorek lab - major discovery group
   WebSearch("[system name] Millman defense island")     # Defense islands
   WebSearch("[system name] Gao anti-phage")             # Gao et al. surveys
   WebSearch("[system name] Tesson DefenseFinder 2024")  # DefenseFinder updates
   ```

3. **For CRISPR subtypes:**
   ```
   WebSearch("CRISPR type [subtype] mechanism effector")
   WebSearch("Cas[number] function structure")
   ```

4. **For ambiguous Cas/transposase calls** (e.g., TnpB vs Cas12f):
   ```
   WebSearch("TnpB Cas12f IS200 IS605 evolution")
   WebSearch("[protein name] CRISPR transposon discrimination")
   ```

**Output fields:**
- `system_name`: Defense system type
- `mechanism`: How it protects against phage (abortive infection, nucleic acid degradation, etc.)
- `components`: Required genes/proteins
- `known_distribution`: Which lineages have it
- `phage_targets`: What types of phage it defends against (if known)
- `key_references`: PMIDs/DOIs
- `confidence`: high (experimentally validated) / medium (computationally predicted, widespread) / low (recently proposed, few examples)

---

### Protocol F: Protein Family Research

**When to use:** Unannotated proteins cluster together by embedding similarity, or a group of proteins shares an annotation but the functional diversity within the family is unclear.

**Search strategy:**

1. **Research annotated cluster members** — if some proteins in the cluster have annotations, use those as entry points:
   ```
   WebSearch("[annotation name] protein family function diversity")
   WebSearch("[annotation name] subfamilies classification")
   ```

2. **For giant/unusual proteins:**
   ```
   WebSearch("[predicted domain] giant protein bacteria")
   WebSearch("[size range] amino acid protein function prokaryote")
   ```

3. **For repeat-domain proteins** (WD40, TPR, HEAT, Ig-like):
   ```
   WebSearch("[repeat type] repeat protein function bacteria archaea")
   WebSearch("[repeat type] scaffold signaling prokaryote")
   ```

4. **For proteins with only structural hits** (Foldseek but no PFAM):
   ```
   # Use Protocol B first, then:
   WebSearch("[structural hit function] homolog [lineage name]")
   ```

**Output fields:**
- `family_name`: Best name for the protein family/cluster
- `annotated_members`: Count and identity of annotated proteins in the cluster
- `functional_diversity`: Range of functions within the family
- `lineage_distribution`: Whether family is lineage-specific or widespread
- `key_references`: PMIDs/DOIs
- `confidence`: high (well-characterized family) / medium (some members characterized) / low (mostly uncharacterized)

---

### Protocol G: Manuscript Citation Research

**When to use:** A manuscript draft exists with claims that need literature citations. This is the primary workflow for turning a data-driven manuscript into a citable publication. Run this BEFORE finalizing any manuscript.

**Input:** Path to manuscript markdown file (e.g., `data/omni_production/MANUSCRIPT.md`)

**Workflow:**

1. **Read the manuscript** and extract every claim that requires a citation:
   - Explicit placeholders: `[CITE: topic]`
   - Comparative claims: "largest known", "first reported", "exceeds", "unprecedented"
   - Background claims: "X is known to...", "previous studies showed..."
   - Methodological references: tools, databases, algorithms cited by name
   - Existing citations that need verification

2. **For each claim, perform a targeted search:**
   ```
   WebSearch("[specific claim keywords] [organism/domain] [year range]")
   ```
   Search at least 2-3 queries per claim to ensure coverage. Prioritize:
   - Primary research papers over reviews (unless the claim is about general knowledge)
   - Recent papers (2020+) for fast-moving fields (defense systems, protein structure)
   - The original discovery paper, not just papers that cite it

3. **For each candidate citation, fetch and verify:**
   ```
   WebFetch("[paper URL or DOI]",
            "What are the key findings? Does this paper support the claim: '[exact claim from manuscript]'?
             Extract the specific sentence or data point that supports this claim.")
   ```

4. **Record structured provenance for every citation** (see output format below)

5. **Check for missing citations:**
   - Search for the study organism/phylum + key findings to find prior work
   - Search for the authors' own prior publications on this topic
   - Search for preprints that may not appear in standard searches

6. **Verify existing citations** (if manuscript already has references):
   - Confirm author names, year, journal, title
   - Confirm the cited claim actually appears in the paper
   - Flag any citation that cannot be verified

**Output:** Write to `{dataset_dir}/exploration/literature_citations.jsonl` with one entry per citation:

```json
{
  "citation_id": "greening_2016_hydrogenase",
  "manuscript_claim": "H2 is a widely utilized energy source for microbial growth (Greening et al., 2016)",
  "claim_location": "Introduction, paragraph 3",
  "reference": {
    "authors": "Greening, C., Biswas, A., Carere, C.R., ..., Morales, S.E.",
    "year": 2016,
    "title": "Genomic and metagenomic surveys of hydrogenase distribution indicate H2 is a widely utilised energy source for microbial growth and survival",
    "journal": "The ISME Journal",
    "volume": "10",
    "pages": "761-777",
    "doi": "10.1038/ismej.2015.153"
  },
  "verification": {
    "source_url": "https://doi.org/10.1038/ismej.2015.153",
    "supporting_quote": "Our results indicate that H2 metabolism is widespread and diversified among microorganisms, suggesting H2 is a widely utilized energy source.",
    "quote_location": "Abstract",
    "verified": true,
    "verification_method": "WebFetch of DOI landing page"
  },
  "justification": "Supports the claim that H2 metabolism is central to syntrophic partnerships. Provides the hydrogenase classification framework (Groups 1-4) used throughout our analysis.",
  "confidence": "high",
  "notes": ""
}
```

**Also output a summary report** to `{dataset_dir}/exploration/citation_report.md`:

```markdown
# Citation Verification Report

## Verified Citations (N)
| # | Citation | Claim | Status | Quote |
|---|----------|-------|--------|-------|
| 1 | Greening et al. 2016 | H2 widely utilized | VERIFIED | "Our results indicate..." |

## Unverified Citations (N)
| # | Citation | Issue |
|---|----------|-------|
| 1 | Smith et al. 2019 | Paper not found at DOI |

## Missing Citations (N)
| # | Manuscript Claim | Suggested Citation | Justification |
|---|-----------------|-------------------|---------------|

## Citation Placeholders Resolved (N)
| # | Placeholder | Resolved To | Supporting Quote |
|---|-------------|-------------|-----------------|
```

**Critical rules:**
- **NEVER fabricate a citation.** If you can't find a paper to support a claim, say so.
- **Always include a supporting quote** — the specific sentence or data point from the paper that backs the manuscript's claim. This allows the human author to verify independently.
- **Flag "from training" citations** — if the manuscript already has a citation and you can verify it, great. If you cannot find it via web search, flag it as `"verified": false, "notes": "Could not locate via web search; may be from LLM training memory"`.
- **Include DOIs** whenever possible — these are permanent, machine-resolvable identifiers.
- **Check for the authors' own prior work** — search for `[first author surname] [organism name]` to find directly relevant prior publications that should be cited.

---

## Structured Output Format

All literature findings should be formatted as findings.jsonl entries for integration with other analysis outputs:

```json
{
  "category": "literature",
  "title": "PDB 3GW6 is a phage T4 tail fiber protein (TM-score 0.82)",
  "description": "Foldseek hit 3gw6_A (TM-score 0.82) is the T4 long tail fiber adhesin gp37, responsible for host recognition. The query protein's structural homology suggests a tail fiber function, consistent with its location in a viral locus containing baseplate and tail tube genes.",
  "provenance": {
    "query": "WebFetch('https://www.rcsb.org/structure/3GW6', 'What is this protein?')",
    "raw_result": "T4 long tail fiber adhesin gp37, binds OmpC on E. coli surface",
    "accession_verified": "N/A (structure-based, not PFAM/KEGG)",
    "interpretation": "Structural homology to T4 tail fiber supports phage tail fiber function"
  },
  "confidence": "high",
  "evidence_type": "structural_homology"
}
```

### Confidence Assessment Guide

| Level | Criteria |
|-------|----------|
| **high** | Experimentally validated function, PDB hit with TM > 0.7, well-characterized domain, >75% pathway complete |
| **medium** | Computational prediction with multiple supporting lines, AlphaFold hit, 50-75% pathway complete, moderately studied lineage |
| **low** | Single weak hit, uncharacterized DUF, <50% pathway complete, candidate phylum with minimal data |

### Evidence Type Tags

Use these in the `evidence_type` field to categorize the basis of the finding:

- `structural_homology` — Foldseek/DALI structural comparison
- `sequence_annotation` — PFAM/KEGG/InterPro domain identification
- `pathway_reconstruction` — KEGG module completeness analysis
- `genomic_context` — Neighborhood/co-localization evidence
- `literature_review` — Published experimental or computational study
- `lineage_context` — Organism/clade-level biological context

---

## Tips

- **Don't guess — look it up.** If you're uncertain about a protein's function based on its annotation name, use the appropriate protocol above. A 30-second WebFetch is cheaper than a retracted finding.
- **Multiple evidence types strengthen claims.** A Foldseek hit (Protocol B) plus genomic context (from explore/survey) plus literature support (Protocol D/E) is much stronger than any one alone.
- **Note when evidence is indirect.** "Structural homology to X suggests function Y" is honest. "Functions as Y" without experimental evidence is overclaiming.
- **Check publication dates.** Defense system biology and protein family characterization move fast. A DUF characterized in 2024 may not appear in older reviews.
- **Report negative results.** "No characterized function found for DUF3683 as of [date]" is a valid and useful finding — it tells future agents not to re-search.
