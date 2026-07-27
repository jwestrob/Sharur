# Co-Location Engine: MacSyFinder Reimplementation Spec

**Status:** Implemented — see `sharur/colocation.py` (wired into `src/ingest/07_build_knowledge_base.py`).
This document is retained as the original design rationale; the shipped implementation has
since evolved beyond it (see **Divergences in implementation** at the end).
**Date:** 2026-03-08 (design); implemented and in production since 2026-03-08.
**Motivation:** MacSyFinder is the bottleneck in defense/secretion system validation. On DPANN (3,564 genomes, 3.4M proteins), a single run takes 20-30 minutes even with `--previous-run` (skipping hmmsearch). It also has fragile failure modes: filename mismatches, hit description parsing crashes, opaque error messages. All the data it needs is already in DuckDB. We can do the co-localization in-process in seconds.

## What MacSyFinder Actually Does

MacSyFinder's core algorithm has 5 steps:

1. **Load models** — parse XML definitions for each system type (mandatory/accessory/neutral/forbidden genes, quorum thresholds, gap constraints)
2. **Load HMM hits** — read hmmsearch output files (one per HMM profile)
3. **Cluster** — group hits on the same contig within `inter_gene_max_space` genes of each other
4. **Validate** — check each cluster against the model: enough mandatory genes? enough total genes? any forbidden genes?
5. **Resolve conflicts** — when clusters overlap (share proteins), pick the best-scoring non-overlapping set

We already have step 2's data in `annotations` table (source='defensefinder' or 'txsscan'). Steps 1, 3-5 are pure logic on that data + protein positions from `proteins` table.

## XML Model Schema

Each system is defined by an XML file:

```xml
<model inter_gene_max_space="5"        <!-- max gene-index gap between consecutive system genes -->
       min_mandatory_genes_required="2" <!-- minimum mandatory genes present -->
       min_genes_required="3"           <!-- minimum total genes (mandatory + accessory) -->
       multi_loci="False"              <!-- may a system be built from several separate
                                            loci on ONE replicon? NOT "across contigs":
                                            a replicon is a single contig, and a system
                                            may never span two. -->
       vers="2.0">

  <gene name="System__GeneA" presence="mandatory"/>
  <gene name="System__GeneB" presence="mandatory">
    <exchangeables>                     <!-- any of these satisfies GeneB -->
      <gene name="System__GeneB_v2"/>
      <gene name="System__GeneB_v3"/>
    </exchangeables>
  </gene>
  <gene name="System__GeneC" presence="accessory" loner="True"/>
  <gene name="System__GeneD" presence="neutral"/>
  <gene name="Other__Conflict" presence="forbidden"/>
</model>
```

### Gene Presence Types

| Type | Semantics | Counts toward quorum? |
|------|-----------|----------------------|
| `mandatory` | Must be present | Yes — both `min_mandatory_genes_required` AND `min_genes_required` |
| `accessory` | Optional, adds confidence | Yes — toward `min_genes_required` only |
| `neutral` | Recorded if found, no effect | No |
| `forbidden` | Presence → system rejected | Negative — invalidates the cluster |

### Gene Attributes

| Attribute | Values | Meaning |
|-----------|--------|---------|
| `loner` | `1`/`True`/`False` | Gene can exist outside `inter_gene_max_space` of the main cluster |
| `multi_system` | `1`/`True`/`False` | Gene can be assigned to multiple system instances |
| `multi_model` | `1` | Gene can match multiple model definitions |

### Exchangeables

Homolog groups. If a gene has `<exchangeables>`, ANY hit from the parent OR any exchangeable child satisfies that gene slot. The parent's presence type applies to the entire group.

### Model Corpus

| Source | Models | Profiles | multi_loci usage | loner usage |
|--------|--------|----------|-----------------|-------------|
| DefenseFinder | 556 | 1,178 | **Never** | Never |
| TXSScan | 19 | 280 | 11/19 systems | 24 genes |

Key implication: DefenseFinder is the simpler case (single-locus only). TXSScan needs `multi_loci` + `loner` support.

## Data Sources (Already in DuckDB)

```sql
-- HMM hits: which proteins matched which HMM profiles
SELECT protein_id, accession, name, score
FROM annotations
WHERE source = 'defensefinder'  -- or 'txsscan'

-- Protein positions: contig + gene order
SELECT protein_id, contig_id, bin_id, gene_index, start, end_coord, strand
FROM proteins
```

The `accession` in annotations corresponds to the HMM profile filename stem (e.g., `Gabija__GajA`). The XML model `<gene name="Gabija__GajA">` matches this accession.

## Algorithm

### Phase 1: Parse Models

```python
@dataclass
class GeneSpec:
    name: str                    # HMM profile name (= annotations.accession)
    presence: str                # mandatory | accessory | neutral | forbidden
    loner: bool = False
    multi_system: bool = False
    exchangeables: list[str] = field(default_factory=list)  # alternative profile names

@dataclass
class SystemModel:
    name: str                    # e.g., "Gabija", "T3SS"
    family: str                  # e.g., "DefenseFinder", "TXSScan"
    inter_gene_max_space: int
    min_mandatory_genes_required: int
    min_genes_required: int
    multi_loci: bool
    genes: list[GeneSpec]

    # Derived at parse time:
    mandatory_genes: set[str]    # all profile names that can satisfy a mandatory slot
    accessory_genes: set[str]    # all profile names for accessory slots
    neutral_genes: set[str]
    forbidden_genes: set[str]
    # gene_name → set of exchangeable alternatives (including self)
    exchangeable_groups: dict[str, set[str]]
```

Parse from XML with `xml.etree.ElementTree`. Recursively resolve exchangeables into flat lookup tables:
- `profile_to_gene_slot`: maps every profile name (including exchangeables) to its parent gene spec
- `profile_to_presence`: maps every profile name to mandatory/accessory/neutral/forbidden

### Phase 2: Load Hits from DuckDB

Single query per annotation source:

```sql
SELECT a.protein_id, a.accession, a.score,
       p.contig_id, p.bin_id, p.gene_index
FROM annotations a
JOIN proteins p ON a.protein_id = p.protein_id
WHERE a.source = 'defensefinder'
ORDER BY p.contig_id, p.gene_index
```

Result: a DataFrame/dict keyed by `(contig_id, gene_index)` → list of `(protein_id, accession, score)`.

### Phase 3: Cluster Hits per Contig

For each contig, iterate through gene positions that have at least one relevant HMM hit. Build clusters using a greedy gap-based approach:

```python
def cluster_hits_on_contig(
    contig_hits: list[tuple[int, str, str, float]],  # (gene_index, protein_id, accession, score)
    max_gap: int,
) -> list[list[tuple]]:
    """Group hits within max_gap gene positions of each other."""
    if not contig_hits:
        return []

    sorted_hits = sorted(contig_hits, key=lambda x: x[0])
    clusters = [[sorted_hits[0]]]

    for hit in sorted_hits[1:]:
        # Can this hit extend the current cluster?
        if hit[0] - clusters[-1][-1][0] <= max_gap:
            clusters[-1].append(hit)
        else:
            clusters.append([hit])

    return clusters
```

**Important nuance:** Different models have different `inter_gene_max_space` values, so
gap enforcement must be **per-model** — clustering once at a global max gap would admit
genes too far apart for a tighter model (a cluster valid at gap=20 may contain genes 15
apart, invalid for a gap=5 model). The implemented approach clusters at `max_gap` per
model. Since DefenseFinder models are all gap ≤ 10 and there are only 556+19 = 575 models
total, this is fast. For each model:
1. Filter hits to only those matching the model's gene set (mandatory + accessory + neutral + forbidden profiles)
2. Cluster those hits on each contig with the model's `inter_gene_max_space`
3. Validate each cluster

### Phase 4: Validate Clusters Against Models

For each (model, cluster) pair:

```python
def validate_cluster(
    cluster: list[tuple],        # (gene_index, protein_id, accession, score)
    model: SystemModel,
) -> Optional[SystemHit]:
    """Check if a cluster satisfies the model's quorum and constraints."""

    # 1. Map hit accessions to gene slots via exchangeable_groups
    satisfied_mandatory = set()
    satisfied_accessory = set()
    satisfied_neutral = set()
    has_forbidden = False

    for gene_idx, pid, accession, score in cluster:
        gene_slot = model.profile_to_gene_slot.get(accession)
        if gene_slot is None:
            continue  # Hit doesn't belong to this model

        presence = gene_slot.presence
        slot_name = gene_slot.name  # canonical name (parent, not exchangeable)

        if presence == "forbidden":
            has_forbidden = True
            break
        elif presence == "mandatory":
            satisfied_mandatory.add(slot_name)
        elif presence == "accessory":
            satisfied_accessory.add(slot_name)
        elif presence == "neutral":
            satisfied_neutral.add(slot_name)

    # 2. Check forbidden
    if has_forbidden:
        return None

    # 3. Check quorum
    n_mandatory = len(satisfied_mandatory)
    n_total = n_mandatory + len(satisfied_accessory)  # neutral doesn't count

    if n_mandatory < model.min_mandatory_genes_required:
        return None
    if n_total < model.min_genes_required:
        return None

    # 4. Gap constraint is enforced during per-model clustering (Phase 3), so a
    #    cluster reaching validation already satisfies the model's inter_gene_max_space.
    #    The implemented _validate_cluster returns None on any quorum/forbidden failure.

    # 5. Build result
    return SystemHit(
        model_name=model.name,
        genes=[(pid, accession) for _, pid, accession, _ in cluster],
        score=sum(s for _, _, _, s in cluster),
        contig_id=...,  # from context
        mandatory_satisfied=n_mandatory,
        total_satisfied=n_total,
    )
```

### Phase 4b: Loner Genes (TXSScan only)

For models with `loner` genes, after validating the core cluster:
1. Find all hits for the model's loner-eligible profiles across the genome (any contig)
2. If a loner hit is NOT within `inter_gene_max_space` of any cluster gene → it's a loner candidate
3. Add loner genes to the system hit (they count toward quorum)
4. Re-check quorum with loners included

This only matters for TXSScan (24 loner genes across 11 models). DefenseFinder never uses loners.

### Phase 4c: Multi-Loci (TXSScan only)

For models with `multi_loci=True`:
1. Collect ALL valid clusters for this model across ALL contigs in the genome
2. Merge them into a single system hit
3. Re-check quorum on the merged set

Again, only TXSScan. DefenseFinder is always single-locus.

### Phase 5: Conflict Resolution

Multiple system models may claim the same protein. Resolve by:

1. **Score each system hit**: sum of HMM scores for all constituent genes
2. **Build conflict graph**: system hits that share any protein are in conflict
3. **Greedy selection**: sort by score descending, greedily select non-conflicting hits

```python
def resolve_conflicts(hits: list[SystemHit]) -> list[SystemHit]:
    """Select best non-overlapping system assignments."""
    hits_sorted = sorted(hits, key=lambda h: h.score, reverse=True)
    used_proteins = set()
    selected = []

    for hit in hits_sorted:
        hit_proteins = {pid for pid, _ in hit.genes}
        if hit_proteins & used_proteins:
            continue  # Conflict with already-selected system
        selected.append(hit)
        used_proteins |= hit_proteins

    return selected
```

## Output Schema

> **Note (implementation):** The columns below describe the engine's intermediate
> result rows. `_integrate_results` remaps these into the live `defense_systems` /
> `secretion_systems` tables, which use `system_type` / `system_subtype` (not `type` /
> `subtype`), add `activity`, `sys_beg`, `sys_end`, `created_at`, and **drop `score` /
> `hit_id`**. A separate `system_proteins` table (`system_id, protein_id, system_source,
> position, profile_name, score`) records per-gene membership.

### systems table (live: `defense_systems` / `secretion_systems`)
| Column | Type | Description |
|--------|------|-------------|
| system_id | TEXT | `{genome_id}_{type}_{N}` |
| genome_id | TEXT | bin_id from proteins table |
| system_type | TEXT | System name (e.g., "Gabija", "T3SS") |
| system_subtype | TEXT | Model name (may differ from type for subtypes) |
| activity | TEXT | System activity class (defense_systems only) |
| protein_ids | TEXT | Comma-separated protein IDs |
| genes_count | INT | Number of genes in system |
| profile_names | TEXT | Comma-separated HMM profile names |
| sys_beg / sys_end | TEXT | First / last protein IDs spanning the system |
| created_at | TIMESTAMP | Validation run timestamp |

### annotations updates
For each gene in a validated system:
```sql
INSERT INTO annotations (protein_id, source, accession, name, score)
VALUES (pid, 'defensefinder_system', system_type, gene_profile, hmm_score)
```

## Implementation Plan

### File: `sharur/colocation.py`

Core engine, no external dependencies beyond stdlib + DuckDB.

```
sharur/colocation.py
├── parse_models(models_dir: Path) -> dict[str, SystemModel]
├── load_hits(db: duckdb.Connection, source: str) -> pd.DataFrame
├── cluster_and_validate(hits_df, models, prot_to_genome) -> list[SystemHit]
│   ├── _cluster_contig(hits, max_gap) -> list[Cluster]
│   ├── _validate_cluster(cluster, model) -> Optional[SystemHit]
│   ├── _add_loner_genes(hit, model, genome_hits) -> SystemHit
│   └── _merge_multi_loci(hits, model) -> SystemHit
├── resolve_conflicts(hits: list[SystemHit]) -> list[SystemHit]
└── validate_systems(db_path, models_dir, source) -> (systems_df, genes_df)
    # Top-level entry point, drop-in replacement for current scripts
```

### Integration

Replace the MacSyFinder subprocess calls in:
- `scripts/validate_defense_systems.py` — `_run_macsyfinder_on_dir()` → `colocation.validate_systems()`
- `scripts/validate_secretion_systems.py` — same replacement
- `src/ingest/07_build_knowledge_base.py` — calls the above scripts

### Performance Estimate

Current MacSyFinder on DPANN: ~20-30 minutes (8 workers, `--previous-run`).

In-process DuckDB approach:
- Phase 2 (load hits): <1s (single indexed query)
- Phase 3+4 (cluster + validate): O(contigs × models) — ~50k contigs × 575 models. Each model only has hits on a fraction of contigs. Realistic: ~2-5s for DefenseFinder, ~1s for TXSScan.
- Phase 5 (conflicts): negligible

**Expected: <10 seconds total.** 100-200x speedup.

### Testing

1. **Unit tests**: Parse known XML models, verify gene specs
2. **Cluster tests**: Known hit positions → expected clusters
3. **Validation tests**: Known clusters + models → expected system hits
4. **Integration test**: Run on Susan genomes (14 MAGs), compare output to MacSyFinder results
5. **Regression test**: Run on DPANN, compare system counts to MacSyFinder output

### Risks / Edge Cases

1. **MacSyFinder's actual clustering may be subtler** — they may use a different gap algorithm (e.g., allowing internal gaps as long as the total span is reasonable). Need to verify against their source code or by comparing outputs.

2. **Score normalization** — MacSyFinder may weight mandatory vs accessory differently in scoring. For conflict resolution, raw HMM score sum should work, but may diverge from MacSyFinder's choices in edge cases.

3. **Gene index gaps** — if Prodigal numbers genes non-contiguously (e.g., some predicted genes filtered out), `gene_index` gaps could produce false cluster splits. Verify that `gene_index` is dense per-contig in our `proteins` table.

4. **`multi_system` semantics** — a protein with `multi_system=True` can be assigned to multiple system instances. This means conflict resolution shouldn't exclude it. Need a flag in the conflict resolver.

5. **Accession matching** — our `annotations.accession` must exactly match the XML `<gene name="...">` values. Verify this mapping for both DefenseFinder and TXSScan profiles.

## Out of Scope (v1)

- **Custom model authoring** — v1 only reads existing DefenseFinder/TXSScan XMLs
- **Quorum scoring nuances** — MacSyFinder has a "best solution" algorithm that considers alternative system assignments globally. Our greedy resolver is simpler but should produce equivalent results in >99% of cases.
- **PADLOC models** — different format, different tool. Separate effort if needed.

## Dependencies

- Python stdlib (`xml.etree.ElementTree`, `dataclasses`, `collections`)
- `duckdb` (already required)
- `pandas` (already required, for output compatibility)

No new dependencies.

## Divergences in implementation

The shipped `sharur/colocation.py` faithfully implements the 5-phase algorithm above, but
has evolved past this design doc in the following ways:

- **Integration path:** Stage 07 (`07_build_knowledge_base.py`) calls `validate_systems` +
  `integrate_defense_results` / `integrate_secretion_results` **directly**. The older
  `scripts/validate_defense_systems.py` / `validate_secretion_systems.py` (MacSyFinder path)
  still exist but are not what the ingest pipeline invokes.
- **Output tables:** see the Output Schema note above — real columns are `system_type` /
  `system_subtype` (+ `activity`, `sys_beg`/`sys_end`, `created_at`), no `score`/`hit_id`,
  plus a `system_proteins` membership table.
- **Conflict resolution** sorts by `(mandatory_fraction, primary_mandatory_hits, score)` and
  exempts `multi_system` proteins — richer than the pure score-descending sort in Phase 5.
- **Post-design bug fixes** not reflected above: HMM internal-name→filename-stem mapping
  (`_build_hmm_name_mapping`), primary-mandatory tracking, multi_loci merge-then-validate
  ordering, and a forbidden-proximity check.
- **`load_hits`** returns a list of hit records, not a `pd.DataFrame`.

**Still outstanding:** the Testing section's unit/regression suite is not yet built — there
are no `tests/` covering `colocation.py`. This overlaps the CLAUDE.md TODO ("Co-location
engine regression test — run on Susan genomes and DPANN, compare against MacSyFinder").

## Concordance with the reference implementation

The engine is only worth having if it agrees with MacSyFinder, and agreement is
measured, not assumed. `scripts/check_colocation_concordance.py` compares the two
on `(genome_id, system_type, frozenset(protein_ids))`.

Two rules make the comparison meaningful:

- **Compare identity, not counts.** Two implementations can report the same
  number of systems while disagreeing about which proteins are in them.
- **Both sides must consume the same hits.** If the engine's input annotations
  came from a different search than the reference run, a mismatch measures
  search settings rather than co-location logic. Profiles outside Pfam generally
  carry no GA/TC/NC cutoffs, so a scan left on defaults will not reproduce the
  reference thresholds; take them from the reference output and set them
  explicitly.

Validate against the single-locus model set first. It exercises clustering and
quorum without `multi_loci`, `loner` or per-gene spacing overrides, so a failure
there is unambiguous. Only then move to model sets using those features, where
the reference's own behaviour is subtle.

### Known divergences and where they come from

Two implementation choices account for every observed difference from the
reference. Both are defensible; neither is obviously wrong; both change results.

**1. When to collapse competing profile hits on one protein.**

The reference keeps one best hit per protein position, sorting by score
descending and taking the first. Python's sort is stable, so when two profiles
tie on score the winner is whichever was processed first — an ordering that
follows profile iteration, not evidence. This engine breaks the same tie
deterministically on e-value.

The consequence appears wherever a short, promiscuous protein hits several
related profiles within a fraction of a bit. Whichever profile wins decides
which system the protein can join, so a tie decided differently removes a gene
from one system and offers it to another.

Matching the reference exactly here means reproducing an arbitrary ordering.
The deterministic rule is the better one to keep; the divergence should be
documented rather than engineered away.

**2. One reference call is unexplained by the reference's own rule.**

The gap rule itself is reproduced exactly: the allowed gap for an adjacent pair
comes from each gene's own `inter_gene_max_space`, resolved so an exchangeable
inherits its parent's spacing, taking the **minimum** when both are defined —
not the maximum — and falling back to the model value only when neither defines
one. Verified against the reference implementation line for line.

Despite that, a system has been observed reported by the reference as a single
locus while spanning more intervening genes than any applicable override allows,
with no loner, `multi_loci` or `multi_system` declared and no intervening hit to
bridge the gap. Applying the shared rule to that case yields two clusters, which
is what this engine produces. The divergence is therefore not a defect in the
reimplementation; it is behaviour in the reference that its own clustering
function does not account for, and it belongs upstream rather than in a local
workaround.

Do not "fix" this by loosening clustering to match. That would trade a correct,
explainable rule for agreement with an unexplained one, and would silently widen
every other model's clusters too.
