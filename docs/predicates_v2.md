# Predicate V2: Typed Semantic Atom System

## Overview

V2 replaces V1's flat boolean predicates with **typed semantic atoms** -- structured claims about proteins that carry facet, relation, and evidence metadata. This enables:

- **Faceted reasoning**: know _what kind_ of claim an annotation makes (activity vs. role vs. architecture)
- **Evidence strength tracking**: distinguish `implies` (KEGG ortholog) from `supports` (PFAM domain) from `flags` (superfamily hit)
- **Conflict resolution**: `excludes` beats `supports`; `implies` + `excludes` = conflict surfaced for review
- **Composite predicates**: declarative YAML rules that combine atoms into higher-order conclusions
- **Unmapped accession review**: surfaces annotation accessions that produce no semantic mapping, prioritized by frequency

V2 is the normal predicate backend for new Stage 07 knowledge-base builds.
Stage 07 writes `semantic_atoms`, `semantic_state`, and `semantic_terms`, then
materializes the legacy `protein_predicates` compatibility table from V2 output
so downstream code (`search_by_predicates`, reports, older scripts) continues
to work unchanged. Shadow comparison remains available for regression testing
and backend changes.

## Data Model

### SemanticFacet

What kind of claim an atom makes. Seven facets:

| Facet | Meaning | Examples |
|-------|---------|---------|
| `activity` | What the protein catalyzes | hydrogenase, kinase, calvin_cycle |
| `role` | Functional role in the cell | defense_system, transporter, regulator |
| `architecture` | Structural features | multi_domain, beta_barrel, coiled_coil |
| `localization` | Where the protein resides | membrane, periplasmic, secreted |
| `topology` | TM helix count, signal peptide | transmembrane_predicted, polytopic_membrane |
| `size_class` | Size category | tiny, small, medium, large, giant, massive |
| `quality_flag` | Annotation quality | unannotated, well_annotated, gc_outlier |

### ClaimRelation

How confidently the source evidence supports the atom. Deterministic rules only -- never LLM-generated.

| Relation | Meaning | Typical sources |
|----------|---------|-----------------|
| `implies` | Strong evidence | KEGG ortholog, HydDB classification, computed properties |
| `supports` | Moderate evidence | PFAM domain consistent with function, CAZy hit |
| `flags` | Weak/ambiguous | VOGdb hit, DefenseFinder HMM (pre-system-validation), superfamily |
| `excludes` | Counter-evidence | Complex I subunit excludes hydrogenase classification |
| `unresolved` | No curated mapping | Annotation accession not in any mapping dict |

### SemanticAtom

A single typed claim about a protein, derived from one annotation.

```python
@dataclass
class SemanticAtom:
    protein_id: str
    atom_id: str              # e.g. "hydrogenase", "membrane", "giant"
    facet: SemanticFacet       # activity, role, architecture, ...
    relation: ClaimRelation    # implies, supports, flags, excludes, unresolved
    source_accession: str      # PF00374, K00001, _computed
    source_db: str             # pfam, kegg, cazy, hyddb, _property, ...
    evidence_evalue: float     # Pass-through from HMM/DIAMOND (optional)
    evidence_score: float      # Pass-through bitscore (optional)
```

One protein typically has many atoms (one per annotation-to-predicate mapping).

### SemanticState

Per-protein aggregated view after resolving all atoms. Conflicts are handled, duplicates removed, unresolved accessions collected.

```python
@dataclass
class SemanticState:
    protein_id: str
    activities: list[str]           # Resolved activity atoms
    roles: list[str]                # Resolved role atoms
    architecture: list[str]         # Resolved architecture atoms
    localization: list[str]         # Resolved localization atoms
    topology: dict[str, Any]        # {"tm_count": 4, "signal_peptide": False}
    size_class: str                 # "giant", "medium", etc.
    quality_flags: list[str]        # ["unannotated", "gc_outlier"]
    composite_predicates: list[str] # ["giant_unannotated", "novel_membrane_protein"]
    unresolved_accessions: list[str]
```

## Facade API (Recommended)

The Sharur facade (`b = Sharur(...)`) exposes V2 through high-level methods. **This is the recommended way for agents and analysis scripts to use V2.**

### Generate V2 atoms + state

```python
from sharur.operators import Sharur
b = Sharur("data/my_dataset/sharur.duckdb")

# Run the full V2 pipeline: generate atoms, aggregate, evaluate composites, persist.
# Full-dataset runs do not return all states by default; read them from DuckDB.
b.generate_v2()

# Subset runs return states in memory by default.
states = b.generate_v2(protein_ids=["protein_123"])

# With review queue output (unmapped accessions ranked by frequency)
b.generate_v2(output_review_queue="review_queue.tsv")

# Optional for ad hoc/manual runs: materialize V2-derived flat predicates into
# the legacy compatibility table.
b.generate_v2(update_legacy_predicates=True)

# Process-parallel full refresh. The parent process owns DuckDB; workers run
# only the CPU transform.
b.generate_v2(
    workers=64,
    chunk_size=100_000,
    update_legacy_predicates=True,
)
```

### Parallel generation and resume

Full-dataset generation splits each database chunk into bounded worker
microbatches. The default microbatch target is about two tasks per worker,
capped at 5,000 proteins. Worker results return in input order and one parent
process appends atoms, states, search terms, compatibility predicates, and the
checkpoint through a single DuckDB connection.

This design holds bulk data memory near one database chunk plus its worker
shards. Each worker also carries one Python interpreter and one copy of the
predicate rule state, so that smaller component scales with worker count.
Stage 07 passes its resolved `--threads` value directly to V2; under Slurm this
resolves from `SLURM_CPUS_ON_NODE`.

Full refreshes append into constraint-free `v2_generation_*` tables. This keeps
the dominant write path sequential and leaves the canonical ART indexes out of
the per-chunk transaction. After the last checkpoint, DuckDB promotes the
complete generation in one transaction and builds the six query indexes once.
Successful publication drops the scratch tables; interrupted runs retain them.

Each chunk updates `v2_generation_checkpoint` in the same transaction as its
scratch rows. Resume validates the semantic code/config fingerprint,
source-table signature, persisted state count, and ordered protein boundary
before continuing:

```bash
# Standard parallel Stage 07 build
python src/ingest/07_build_knowledge_base.py \
  --data-dir data/my_dataset \
  --output data/my_dataset/sharur.duckdb \
  --threads 64 \
  --force

# Continue an interrupted V2 phase from its latest committed chunk
python src/ingest/07_build_knowledge_base.py \
  --data-dir data/my_dataset \
  --output data/my_dataset/sharur.duckdb \
  --threads 64 \
  --resume-v2

# Start a fresh V2 generation while preserving proteins, annotations, loci,
# validated callers, and all other completed upstream Stage 07 products
python src/ingest/07_build_knowledge_base.py \
  --data-dir data/my_dataset \
  --output data/my_dataset/sharur.duckdb \
  --threads 64 \
  --restart-v2
```

The facade exposes the same recovery path with
`b.generate_v2(workers=64, resume=True, return_states=False)`. Resume applies
to a compatible interrupted generation. `--restart-v2` selects a new semantic
generation contract and reuses the completed upstream database tables. Subset
refreshes retain their targeted replacement behavior.

Review queue generation performs a set-oriented aggregation over persisted
`semantic_atoms` after the final chunk. Its Python memory footprint therefore
stays independent of the number of unresolved atoms, and the same aggregation
can be regenerated after a recovered run.

Validated defense and secretion systems are folded in automatically during
generation. V2 reads `defense_systems` and `secretion_systems` directly and
creates synthetic `defensefinder_system` / `txsscan_system` evidence per member
protein, so validated system atoms do not depend on duplicate rows existing in
the generic `annotations` table. Membership is normalized into
`system_proteins` during generation; the original delimited `protein_ids` fields
remain only for compatibility.

Validated DefenseFinder system rows also emit controlled subtype atoms for
high-yield system families and obvious grouped subfamilies. Examples include
`rm_type_i`, `rm_type_ii`, `defense_eleos`, `defense_hec`, `defense_shedu`,
`defense_pago`, `abi_e`, `defense_pd_lambda`, `defense_retron`, and
`defense_mokosh`.

### Query semantic state

```python
# Get the resolved semantic state for a protein
state = b.get_semantic_state("protein_123")
state.activities     # ["hydrogenase", "nife_hydrogenase"]
state.roles          # ["defense_system"]
state.architecture   # ["multi_domain"]
state.size_class     # "giant"
state.composite_predicates  # ["energy_conserving_hydrogenase"]

# Get raw atoms (which annotation produced each claim)
atoms = b.get_atoms("protein_123")
for a in atoms:
    print(f"{a.atom_id} ({a.facet.value}): {a.relation.value} via {a.source_db}:{a.source_accession}")
```

### Explain one protein

Use `explain()` when you need to audit why a protein matched a search term.
There is one facade method, not separate V1/V2 variants.

```python
explanation = b.explain("protein_123")
explanation["semantic_state"]       # resolved V2 state dict
explanation["atoms"]                # raw atom evidence
explanation["direct_access_terms"]  # e.g. ["pfam:PF00005"]
explanation["composite_terms"]      # matched YAML composites
explanation["validated_systems"]    # normalized system membership rows

# Why did a composite match?
for witness in explanation["composite_explanations"]["restriction_modification_validated"]:
    print(
        witness["atom_id"],
        witness["relation"],
        witness["source_db"],
        witness["source_accession"],
    )
```

`composite_explanations` maps each matched composite to the positive atom
witnesses that satisfied its rule. `none_of` clauses contribute no witnesses
because their successful state is absence.

### Search by V2 facets and atoms

**Default limit=50.** All search methods return at most 50 results by default to avoid flooding context. Pass `limit=N` for larger result sets, or use SQL aggregation (`COUNT(DISTINCT protein_id)`) for counts.

```python
# Search by facet: find all proteins with a specific type of claim
b.search_by_facet("activity")                                      # any activity
b.search_by_facet("activity", atom_ids=["hydrogenase"])            # specific activity
b.search_by_facet("role", atom_ids=["defense_system"], relation="implies")  # strong evidence only
b.search_by_facet("activity", limit=500)                           # override default limit

# Search by atom presence/absence (V2-native analog of search_by_predicates)
b.search_by_atoms(has=["giant", "unannotated"])         # giant + unannotated
b.search_by_atoms(has=["hydrogenase"], lacks=["membrane"])  # hydrogenase, not membrane
b.search_by_atoms(has=["defense_system_validated"])     # validated DefenseFinder system members
b.search_by_atoms(has=["secretion_system_validated"])   # validated TXSScan system members
b.search_by_atoms(has=["pfam:PF00005"])                 # direct accession key

# Source/relation-aware atom evidence search
b.search_atoms(
    atom_id="hydrogenase",
    relation="implies",
    source_db="hyddb",
    limit=500,
)
```

Use the validated composites for system-level biological claims. Bare
`defense_system` and `secretion_system` atoms can come from weaker component or
domain evidence such as PFAM/supports or raw profile flags; they are useful for
candidate discovery, not for reporting validated systems.

### Composite predicates (DSL)

```python
# List all composite predicate definitions
for comp in b.list_composites():
    print(f"{comp.name}: {comp.description}")

# Composite predicates are auto-evaluated during generate_v2() and stored
# in semantic_state.composite_predicates. Query them via get_semantic_state:
state = b.get_semantic_state("protein_123")
if "nife_hydrogenase_validated" in state.composite_predicates:
    print("Validated NiFe hydrogenase!")
```

### Review queue and shadow diff

```python
# Unmapped accessions: what annotation hits have no semantic mapping yet?
queue = b.v2_review_queue(limit=20)
for entry in queue:
    print(f"{entry['accession']} ({entry['source_db']}): {entry['n_proteins']} proteins")

# SQL-native filters for curation runs
b.v2_review_queue(
    limit=10000,
    source=["pfam", "kofam", "kegg"],
    min_proteins=10,
    exclude_raw_system_profiles=True,
    output_tsv="reports/predicates_v2_review_queue.tsv",
)

# V1 vs V2 comparison
diff = b.run_shadow_diff(output_report="shadow_report.md")
print(f"Match rate: {diff['summary']['match_rate_pct']}%")
```

Use the shadow gate helper when promoting V2 into release or ingestion
workflows:

```python
from sharur.predicates_v2 import evaluate_shadow_gate

gate = evaluate_shadow_gate(
    diff,
    min_match_rate_pct=95.0,
    max_v2_unresolved_atoms=1000,
)
assert gate["passed"], gate["failures"]
```

### Legacy compatibility

Normal Stage 07 builds write `semantic_atoms`, `semantic_state`, and a
V2-derived `protein_predicates` compatibility table. For manual
`b.generate_v2()` calls, pass `update_legacy_predicates=True` when you also
want to refresh that compatibility table. The conversion preserves source
flags, confidence flags, topology aliases, and direct accession keys
(`pfam:PF00005`, `kegg:K00001`, etc.).

`search_by_predicates()` treats a partial `protein_predicates` cache as stale.
If complete V2 state is present, it repairs the compatibility table from
`semantic_state`/`semantic_atoms` before searching. If neither complete V2 state
nor a complete compatibility table exists, it fails loudly instead of returning
truncated results.

```python
# Only p1/p2 V2 rows are replaced. Existing V2 rows for other proteins remain.
b.generate_v2(protein_ids=["p1", "p2"])

# Replace the legacy compatibility cache with V2-compatible output.
b.generate_v2(update_legacy_predicates=True)
b.search_by_predicates(has=["giant", "unannotated"])
```

For shell refreshes, use the V2-backed CLI:

```bash
sharur compute-predicates \
  --db data/my_dataset/sharur.duckdb \
  --chunk-size 100000 \
  --workers 64 \
  --review-queue data/my_dataset/reports/predicates_v2_review_queue.tsv

# Recovery after an interrupted full refresh
sharur compute-predicates \
  --db data/my_dataset/sharur.duckdb \
  --chunk-size 100000 \
  --workers 64 \
  --resume
```

## Library API (Advanced)

For custom pipelines or extending V2 internals, import directly from `sharur.predicates_v2`:

```python
from sharur.predicates_v2 import (
    # Data model
    ClaimRelation, SemanticAtom, SemanticFacet, SemanticState,
    # Pipeline components
    AtomGenerator, aggregate_atoms, evaluate_composites, load_composites,
    create_v2_tables, generate_and_persist_v2,
    materialize_semantic_terms_from_v2, materialize_system_proteins,
    # V1 compatibility
    semantic_state_to_predicates,
    # Review queue
    build_review_queue, format_review_queue_tsv,
    # Shadow diff
    shadow_diff, run_shadow_diff_on_store, evaluate_shadow_gate,
)
```

### Step-by-step pipeline

```python
from sharur.predicates.generator import AnnotationRecord, ProteinRecord
from sharur.predicates_v2 import (
    AtomGenerator, aggregate_atoms, evaluate_composites,
    semantic_state_to_predicates,
)

# 1. Create generator
gen = AtomGenerator(expand_hierarchy=True, predict_topology=False)

# 2. Build records
protein = ProteinRecord(protein_id="prot_001", sequence_length=450, gc_content=0.35)
annotations = [
    AnnotationRecord(source="pfam", accession="PF00374", name="NiFeSe_Hases",
                     description="NiFe hydrogenase", evalue=1e-50, score=180.0),
    AnnotationRecord(source="kegg", accession="K00436", name="hoxH",
                     description="NiFe hydrogenase", evalue=1e-80, score=250.0),
]

# 3. Generate atoms
atoms = gen.generate_atoms(protein, annotations)

# 4. Aggregate into state
state = aggregate_atoms(protein.protein_id, atoms)

# 5. Evaluate composite predicates
composites = evaluate_composites(atoms, topology=state.topology)
state.composite_predicates = composites

# 6. Convert back to V1 format if needed
flat_predicates = semantic_state_to_predicates(state, atoms=atoms)
```

## Composite Predicate DSL

Composite predicates are declared in YAML and evaluated against a protein's atom set. They combine multiple atoms into higher-order conclusions.

### Config file

`config/predicates_v2/composites.yaml`

### Operators

#### `has_atom`

Check if a protein has a specific atom. Supports optional filters:

```yaml
# Basic: protein must have atom "hydrogenase"
- has_atom: hydrogenase

# With relation filter: atom must have exactly this relation
- has_atom: defense_system
  relation: implies

# With relation_in filter: atom relation must be one of these
- has_atom: nife_hydrogenase
  relation_in: [implies, supports]

# With source_db filter: atom must come from this database
- has_atom: defense_system
  relation: implies
  source_db: defensefinder_system
```

#### `all_of`

All sub-conditions must be true:

```yaml
requires:
  all_of:
    - has_atom: giant
    - has_atom: unannotated
```

#### `any_of`

At least one sub-condition must be true:

```yaml
requires:
  any_of:
    - has_atom: nife_group4
    - has_atom: mbh_hydrogenase
    - has_atom: ech_hydrogenase
```

#### `none_of`

No sub-condition may be true (exclusion):

```yaml
requires:
  all_of:
    - has_atom: nife_hydrogenase
      relation_in: [implies, supports]
  none_of:
    - has_atom: complex_I_subunit
      relation: implies
```

When a condition dict contains multiple compound keys (`all_of` + `none_of`), ALL must be satisfied.

#### Property operators

Operate on the topology dict from SemanticState:

```yaml
# Exact match
- property_equals: {tm_count: 0}

# Greater than or equal
- property_gte: {tm_count: 4}

# Less than or equal
- property_lte: {tm_count: 1}
```

### Full composite example

```yaml
nife_hydrogenase_validated:
  description: "NiFe hydrogenase confirmed by PFAM corroboration"
  facet: activity
  requires:
    all_of:
      - has_atom: nife_hydrogenase
        relation_in: [implies, supports]
      - has_atom: nickel_binding
    none_of:
      - has_atom: complex_I_subunit
        relation: implies
```

### Bundled composites (20 definitions)

| Name | Purpose |
|------|---------|
| `nife_hydrogenase_validated` | NiFe hydrogenase with PFAM corroboration, no Complex I |
| `fefe_hydrogenase_validated` | FeFe hydrogenase with iron-sulfur cluster |
| `energy_conserving_hydrogenase` | Group 4 NiFe (Ech/Mbh) |
| `novel_membrane_protein` | Unannotated + TM helices |
| `giant_unannotated` | Giant + no annotation |
| `giant_multi_domain` | Giant + multi_domain |
| `defense_system_validated` | DefenseFinder system-level validation |
| `crispr_cas_supported` | CRISPR-Cas associated protein supported by curated component evidence |
| `restriction_modification_validated` | RM system confirmed by DefenseFinder system-level validation |
| `abc_transporter_complete` | ABC transporter + ATPase |
| `polytopic_transporter` | Transporter + polytopic membrane |
| `carbon_fixation_active` | RuBisCO + PRK (Calvin cycle) |
| `methanogenesis_confirmed` | MCR complex present |
| `nitrogen_fixation_confirmed` | Nitrogenase + nitrogen_fixation |
| `cellulosome_component` | Dockerin or cohesin domain |
| `secretion_system_validated` | TXSScan system-level validation |
| `gc_outlier_unannotated` | GC outlier + unannotated (possible HGT) |
| `well_characterized` | well_annotated + multi_source |

## Config Files

All config lives in `config/predicates_v2/`.

### `facet_assignments.yaml`

Maps each V1 predicate ID to a SemanticFacet. 592 entries, grouped by facet.

```yaml
# Format: predicate_id: facet
hydrogenase: activity
defense_system: role
membrane: localization
multi_domain: architecture
giant: size_class
unannotated: quality_flag
```

Predicates not in this file fall back to their V1 vocabulary category, then default to `role`.

### `relation_overrides.yaml`

Controls the ClaimRelation assigned to atoms.

```yaml
source_defaults:
  kegg: implies       # KEGG orthologs are strong evidence
  pfam: supports      # PFAM domains are moderate evidence
  vogdb: flags        # VOGdb hits are weak/ambiguous
  hyddb: implies      # HydDB is authoritative for hydrogenases
  defensefinder: flags        # Raw HMM hits (pre-validation)
  defensefinder_system: implies  # System-validated hits

accession_overrides:
  PF00374: implies    # NiFeSe_Hases -- definitive NiFe marker
  PF04055: flags      # Radical_SAM -- broad superfamily
```

Source defaults apply to all accessions from that source. Accession overrides take precedence.

### `composites.yaml`

See "Composite Predicate DSL" section above.

### `duf_curations.yaml`

Empty template for DUF structural assignment curations. Intended for recording Foldseek-based DUF function assignments that should modify atom facets or relations.

## Adding New Rules

### Adding a new composite predicate

1. Edit `config/predicates_v2/composites.yaml`
2. Add a new entry:
   ```yaml
   my_new_composite:
     description: "What this composite detects"
     facet: role  # or activity, quality_flag, etc.
     requires:
       all_of:
         - has_atom: some_atom
         - has_atom: another_atom
   ```
3. Run tests: `python -m pytest tests/test_predicates_v2/test_composites.py -v`

### Adding a facet override

If a predicate is assigned to the wrong facet (e.g., a localization predicate categorized as `role`):

1. Edit `config/predicates_v2/facet_assignments.yaml`
2. Change the facet: `my_predicate: localization`
3. Run tests: `python -m pytest tests/test_predicates_v2/test_rules.py -v`

### Adding a relation override

To promote or demote a specific accession:

1. Edit `config/predicates_v2/relation_overrides.yaml`
2. Add under `accession_overrides`: `PF99999: implies`
3. Run tests: `python -m pytest tests/test_predicates_v2/test_generator.py -v`

### Curating the review queue

1. Run V2 on a dataset:
   ```python
   generate_and_persist_v2(store, output_review_queue="review_queue.tsv")
   ```
2. Open `review_queue.tsv` -- accessions sorted by priority (n_proteins * n_genomes)
3. For each high-priority unmapped accession:
   - Determine the correct predicate(s) it should map to
   - Add the mapping to the appropriate V1 mapping file (`pfam_map.py`, `kegg_map.py`, etc.)
   - Or add an accession override in `relation_overrides.yaml`
4. Re-run V2 to verify the accession is no longer unresolved

## DuckDB Tables

Current V2 tables are part of schema version 5. The core V2 tables can also be
created explicitly via
`create_v2_tables(store)`:

Schema version 5 adds contig provenance to validated system calls. Because the
pre-v5 caller did not enforce replicon-local MacSyFinder semantics, migration
quarantines those legacy named calls and removes semantic state/search-cache
rows for their former members. Raw DefenseFinder/TXSScan observations remain.
Rerun the current caller and refresh V2 for the affected proteins before using
named system evidence from a migrated legacy dataset.

### `semantic_atoms`

Raw evidence layer. One row per (protein, atom, source_accession) triple.

| Column | Type | Description |
|--------|------|-------------|
| protein_id | VARCHAR | FK to proteins |
| atom_id | VARCHAR | Semantic tag |
| facet | VARCHAR | activity, role, etc. |
| relation | VARCHAR | implies, supports, etc. |
| source_accession | VARCHAR | PF00374, K00001, _computed |
| source_db | VARCHAR | pfam, kegg, _property |
| evidence_evalue | DOUBLE | Pass-through e-value |
| evidence_score | DOUBLE | Pass-through bitscore |

Primary key: `(protein_id, atom_id, source_accession)`

### `semantic_state`

Resolved per-protein view. One row per protein.

| Column | Type | Description |
|--------|------|-------------|
| protein_id | VARCHAR | PK |
| activities | VARCHAR[] | Resolved activity atoms |
| roles | VARCHAR[] | Resolved role atoms |
| architecture | VARCHAR[] | Resolved architecture atoms |
| localization | VARCHAR[] | Resolved localization atoms |
| topology | JSON | {"tm_count": N, "signal_peptide": bool} |
| size_class | VARCHAR | giant, medium, etc. |
| quality_flags | VARCHAR[] | Annotation quality indicators |
| composite_predicates | VARCHAR[] | Matched composite names |
| unresolved_count | INTEGER | Number of unmapped accessions |
| updated_at | TIMESTAMP | Last generation time |

### `semantic_terms`

Materialized search view. It stores:

- `term_kind='atom'`: semantic atom IDs, with facet/relation/source metadata
- `term_kind='direct_access'`: accession keys such as `pfam:PF00005`
- `term_kind='composite'`: YAML composite predicate names

`search_by_atoms()` uses this table when populated and falls back to
`semantic_atoms`/`semantic_state` only for older databases that have not been
backfilled.

### `v2_generation_checkpoint`

Operational state for a full semantic refresh. The row stores the semantic
code/config fingerprint, source-table signature, status, ordered
`last_protein_id`, processed count, and total count. Chunk products and the
ordered boundary commit atomically. A completed refresh retains the row as
durable provenance for the generation boundary.

During an active full refresh, `v2_generation_atoms`,
`v2_generation_state`, `v2_generation_terms`, and `v2_generation_legacy`
hold the committed scratch generation. Their rows share a transaction boundary
with the checkpoint. Canonical tables remain publication targets and receive
the complete generation atomically at the end.

### `system_proteins`

Normalized member table for validated systems:

| Column | Type | Description |
|--------|------|-------------|
| system_id | VARCHAR | Validated DefenseFinder/TXSScan system ID |
| protein_id | VARCHAR | Member protein |
| system_source | VARCHAR | `defensefinder_system` or `txsscan_system` |
| position | INTEGER | Order within the validated system |
| profile_name | VARCHAR | Matched system profile, when available |
| score | DOUBLE | System-level score/count proxy |

## Testing

```bash
# V2 tests only
python -m pytest tests/test_predicates_v2/ -v

# V2 + schema migration tests
python -m pytest tests/test_predicates_v2 tests/test_schema_migration.py -q --no-cov
```

---

## Appendix: Design rationale & history

> Lifted from the original `PREDICATE_V2_SPEC.md` build plan (now retired). Preserved here
> for the *why* behind V2. Note one premise of that spec is now obsolete: V2 was designed to
> run in **shadow mode** alongside V1 with switchover "deferred to a separate PR." That
> switchover has since happened — V2 is the normal Stage 07 backend and materializes the
> legacy `protein_predicates` compat table (commit `ebc141a`).

### Why V2 — the three V1 limitations it set out to fix

The flat V1 predicate system (≈547 predicates, ~8000 lines of vocabulary + mappings +
generator) worked but had three structural limitations:

1. **Flat booleans.** `hydrogenase`, `membrane`, and `giant` were all the same type — no
   way to ask "give me all activity-type predicates" or "what localization claims exist for
   this protein," and no way to distinguish a structural observation from a functional claim.
2. **Hard-coded composite logic.** Rules like "NiFe hydrogenase requires PF00374 AND not
   Complex I" lived in Python methods, so adding a validation rule meant editing code, not
   config.
3. **No curation feedback.** Every new dataset had hundreds of annotation accessions that hit
   no mapping and silently produced no predicates, with no mechanism to surface
   "PF12345 appears in 847 proteins and has no mapping."

### Rule-loading strategy (why V1 mappings were not rewritten as YAML)

V2 deliberately **imports the existing V1 mapping dicts directly** rather than porting ~6200
lines of Python mappings to YAML:

1. `rules.py` imports `PFAM_TO_PREDICATES`/`PFAM_PATTERNS`, `KEGG_TO_PREDICATES`, etc. from
   `sharur/predicates/mappings/`.
2. For each emitted V1 predicate, V2 looks up its facet in `facet_assignments.yaml` and its
   relation in `relation_overrides.yaml`. Defaults: facet from the V1 category→facet table
   (below); relation by source — `kegg`/`hyddb`/`defensefinder_system` → `implies`,
   `pfam`/`cazy` → `supports`, `vogdb`/`defensefinder` (raw HMM) → `flags`.
3. YAML overrides take precedence — this is how specific PFAMs get promoted from `supports`
   to `implies` (e.g. PF00374 NiFeSe_Hases → `implies hydrogenase`).
4. Genuinely new rules with no V1 equivalent go directly in YAML.

### V1 category → V2 facet default mapping

```yaml
enzyme: activity        metabolism: activity     cazy: activity
transport: role         regulation: role         binding: role
mobile: role            stress: role             info_processing: role
division: role          viral: role
structure: architecture
envelope: localization
size: size_class
annotation: quality_flag   composition: quality_flag
topology: topology
```

Predicates needing a manual override (e.g. `membrane`, V1 category `transport` but facet
`localization`) are pinned in `facet_assignments.yaml`.

### Original non-goals (V2 v1 scope)

- DSL query language for agents (was future work; `search_by_atoms`/composites now cover much of this)
- Neighborhood/context-based atoms — V2 is intrinsic (per-protein) only; locus-level semantics remain future work
- Automatic curation of the review queue — curation is manual (agents add YAML entries)
- ~~Production switchover~~ — **done** (see note above)
- Migration of V1 mapping files to YAML — V2 imports the Python dicts; YAML is for overrides/new rules only
