# How Sharur Predicates Are Built

Sharur predicates are statements about proteins drawn from annotation sources.
Each predicate traces to one annotation hit and to a recorded reason that the
annotation supports it. This document describes those reasons: the sources, the
rules that combine them, the thresholds, and the places where the rules make
judgment calls. To see the chain for one protein, run
`sharur why PROTEIN PREDICATE --db ...`.

Counts below are for the maps built on 2026-09-25 (Pfam 38.2 clans with the
installed Pfam-A.hmm, Swiss-Prot 2026_03, KEGG REST ko/brite dated 2026-09-24).

## Principles

1. **Observed versus named.** A protein's domains and orthology hits are
   observations. A predicate names a function, cofactor, pathway role or
   system, and ships only when a recorded source states it.
2. **Evidence per pair.** Every (Pfam family | KO | CAZy family) → predicate
   pair carries one evidence string. The test suite re-derives every string
   from the snapshot of the sources it came from.
3. **Proposals are separate from evidence.** Curators and patterns propose
   pairs; a proposal ships only with evidence. Unsupported proposals are
   written to a report and stay out of the map.
4. **Components and systems.** Single genes and domains yield component-level
   predicates. System-level predicates (a validated defense system, a paired
   toxin-antitoxin system, a secretion-system locus) come only from
   purpose-built system callers.

## Maps and their evidence

### Pfam families (`sharur/predicates/mappings/data/pfam_predicates.tsv`, shipped)

30,024 pairs over 17,020 families, built by `scripts/build_pfam_predicate_map.py`.

| Evidence | Pairs | Meaning |
|---|---|---|
| `text:<match>` | 23,042 | The family's Pfam name or description states the predicate |
| `go:<GO id>` | 3,930 | The family's InterPro GO annotation (pfam2go), closed over is_a/part_of ancestors, reaches a GO term anchored to the predicate |
| `swissprot:...` | 2,837 | Reviewed proteins carrying the family agree (below) |
| `enzyme:<EC> (<name>)` | 215 | The description names an Expasy ENZYME enzyme whose EC class maps to the predicate |

GO anchors and text patterns for each predicate live in
`sharur/predicates/mappings/pfam_evidence_spec.py`. Pfam family names and
descriptions come from the installed Pfam-A.hmm, then the current
Pfam-A.clans.tsv for families the installed release lacks.

### KEGG orthologs (built locally by `sharur setup-kegg`)

52,180 pairs over 19,501 KOs in the current local build. KEGG is copyright
Kanehisa Laboratories and KEGG REST is for academic use, so Sharur ships its
rules and ID-keyed counts only; each user builds the map from KEGG data on their
own machine (`--inputs` accepts a licensed KEGG copy).

| Evidence | Pairs | Meaning |
|---|---|---|
| `ec:<EC>` | 24,795 | KEGG's KO definition lists the EC number |
| `brite:<hierarchy> <path>` | 12,162 | KEGG places the KO under a BRITE node mapped in `kegg_brite_predicates.tsv` |
| `text:<match>` | 11,170 | The KO's KEGG symbols or name state the predicate |
| `swissprot:k/n` | 2,667 | Reviewed proteins linked to the KO agree (shipped table `kegg_swissprot_consensus.tsv`) |
| `module:<M number>` | 1,321 | The KO is in a KEGG module mapped in `kegg_module_predicates.tsv`, or joined to other components by `+`/`-` in a module definition (`complex_subunit`) |
| `hyddb:k/n` | 65 | For KOs KEGG names as hydrogenases: HydDB reference labels captured by the KO's KOfam profile |

BRITE and module rules map component-level predicates only. Module rules cover
central carbon, fatty acid and cofactor pathways, respiratory and photosynthetic
complexes, nitrogenase and SOX; modules whose components also serve other
pathways stay out (for example Hdr, DsrAB, NarGHI, AmoA/PmoA).

### CAZy families (`cazy_predicates.tsv`, shipped)

1,151 pairs over 346 families, built by `scripts/build_cazy_predicate_map.py`.
`class:<GH|GT|PL|CE|CBM|AA>` (973 pairs) states CAZy's class definitions;
`swissprot:...` (178) states family-level activities on which reviewed
proteins agree. Carbohydrate-binding modules take class predicates only. The CE
class states hydrolase activity, since it includes de-N-acetylases. Families
absent from the map take their class predicates; subfamilies resolve to their
family.

### VOGdb

VOG functional categories (`category:Xr` etc.) and consensus-description
patterns. A pattern match counts only when it passes the text conventions below;
`vog_map.vog_evidence` returns the evidence per predicate.

### EC numbers

`EC_TO_PREDICATES` in `kegg_map.py` maps EC classes to predicates, using the
most specific entry. Cofactor and oxygenase entries follow ENZYME's class
definitions (`enzclass.txt`): NAD(P)+ acceptor classes give `nad_binding`,
flavin acceptor classes `flavin_binding`, and 1.13/1.14 sub-subclasses separate
mono- from dioxygenases. Substrate predicates (cellulase, chitinase, amylase,
xylanase, mannanase, pectinase, lysozyme, LPMO) come only from ECs whose ENZYME
name states the substrate. `ec_sam_dependent.tsv` lists the 379 ECs whose
ENZYME reaction consumes S-adenosyl-L-methionine; these add `sam_binding`.

## Text evidence conventions

Name and description matching (Pfam, KEGG, VOG) follows one set of rules
(`pfam_evidence.vetted_match`):

- Lowercase alternatives match whole words case-insensitively; alternatives with
  capitals are case-sensitive symbols. Every alternative needs a leading word
  boundary, and short alternatives also a trailing one.
- Enzyme-class suffixes may sit inside compound words, with exclusions
  (dehydrogenase is not a hydrogenase; oxidoreductase not a reductase).
- `-like` matches are rejected.
- Relations name something else: "X inhibitor", "X-activating" reject the text;
  "X-binding", "X-associated", "X-interacting" skip that match (except for
  `*_binding` / `*_associated` predicates); "regulator of X", "permease for X"
  skip it.
- `HOMONYMS` lists family/predicate pairs where the wording is a known
  homonym.

## Reviewed-protein (Swiss-Prot) consensus

Swiss-Prot is used only when maps are built. For each reviewed protein Sharur
collects predicates from:

- curator-assigned EC numbers;
- GO terms with experimental evidence codes (EXP, IDA, IPI, IMP, IGI, IEP and
  high-throughput equivalents);
- curated `COFACTOR` annotations (ChEBI → metal, Fe-S, heme, PLP, flavin, TPP,
  NAD(P), molybdopterin, cobalamin predicates);
- `CATALYTIC ACTIVITY` reactions that consume S-adenosyl-L-methionine;
- for Pfam consensus, the KEGG-evidenced predicates of the protein's KOs.

A family or KO gains a predicate when (thresholds in `pfam_evidence.py`):

- **Coverage:** ≥80% of ≥5 reviewed carriers agree, with a Wilson 95% lower
  bound ≥0.5.
- **Ortholog majority:** when carriers span two or more KOs, most KOs agree
  (each judged by its own carriers), so one heavily reviewed ortholog cannot
  speak for a family.
- **Attribution (Pfam, CAZy):** ≥80% of ≥3 single-domain carriers agree, or,
  with fewer of them, no co-occurring family that carries the predicate on its
  own evidence appears in ≥90% of the supporting proteins.

Excluded from consensus: fold, repeat and bookkeeping predicates, and
subcellular compartments (cytoplasmic, secreted, periplasmic, membrane and
related), which depend on the organism's architecture.

The evidence string records the counts, for example
`swissprot:305/315 KOs 34/35 single-domain 90/90`. Setting `SHARUR_SWISSPROT`,
`SHARUR_SWISSPROT_KEGG` and `SHARUR_GO_OBO` (or `make recount`) recounts every
such string from the raw files.

## From annotations to predicates

1. Each annotation hit maps through the relevant table to predicates.
2. Annotation sources other than system callers (`defensefinder_system`,
   `txsscan_system`) contribute component-level predicates only; a system-level
   predicate from such a source becomes its component equivalent
   (`vocabulary.COMPONENT_EQUIVALENT`).
3. Hierarchy expansion adds each predicate's `parent` (is-a) chain. No
   component predicate has a system-level ancestor. `part_of` records membership
   in a system and is never expanded.
4. V2 atoms record the relation of each claim by source: KEGG/KOfam, HydDB
   class and system callers `implies`; Pfam, CAZy and HydDB subgroups
   `supports`; VOG and raw profile flags `flags`.

## Hydrogenases

HydDB subgroup labels come from Sharur's nearest-reference DIAMOND assignment
against HydDB references, interpreted through HydDB Table 1
(`sharur/hydrogenase/subgroups.py`); putative and unresolved subgroups carry
structure only. Each assignment records a Pfam catalytic-domain check and KOfam
support: whether the protein's KO hits capture HydDB references of the assigned
subgroup (`hydrogenase_classifications.ko_support`). KOfam support corroborates
an assignment and leaves it unchanged.

## Provenance and verification

- Every predicate generation appends a `predicate_provenance` row with hashes
  of the Pfam, KEGG and CAZy maps, VOG rules, vocabulary and V2 config, the KEGG
  release dates, and the git commit. `sharur preflight` and `sharur describe`
  report whether a database's predicates match the installed maps.
- `make snapshots predicate-maps recount` downloads dated source snapshots,
  rebuilds every map in dependency order, and re-verifies it; rebuilding from
  the same snapshot reproduces the shipped maps byte for byte.

## Known soft spots

- **Text evidence is the largest Pfam tier.** Family names are curated, and the
  conventions reject the common failure modes, but a name states what the
  family's characterized members do.
- **Complex-level ECs.** Swiss-Prot and KEGG assign an enzyme complex's EC to
  each subunit, so some activity predicates describe the complex a subunit
  belongs to; `complex_subunit` marks module complex components.
- **Curated cofactors are mostly inferred.** Most Swiss-Prot COFACTOR records
  are curator inferences by similarity; few cite experiments.
- **Reviewed proteins favor model organisms.** Ortholog-majority and
  attribution rules limit the effect; for lineages such as DPANN and CPR, many
  proteins carry no functional predicate.
- **Family-level claims.** A predicate on a family applies to every carrier;
  polyspecific families keep only what their reviewed members share.

## Editing the maps

Propose pairs in the proposal files (`pfam_predicate_proposals.tsv`,
`kegg_predicate_proposals.tsv`, `cazy_predicate_proposals.tsv`) or add rules
(`kegg_brite_predicates.tsv`, `kegg_module_predicates.tsv`, GO anchors and text
patterns in `pfam_evidence_spec.py`), then rebuild. The dropped-proposal reports
list what lacks evidence.
