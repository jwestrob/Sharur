# Predicate V2: Typed Semantic Atom System

## Overview

V2 replaces V1's flat boolean predicates with **typed semantic atoms** -- structured claims about proteins that carry facet, relation, and provenance metadata. This enables:

- **Faceted reasoning**: know _what kind_ of claim an annotation makes (activity vs. role vs. architecture)
- **Evidence strength tracking**: distinguish `implies` (KEGG ortholog) from `supports` (PFAM domain) from `flags` (superfamily hit)
- **Conflict resolution**: `excludes` beats `supports`; `implies` + `excludes` = conflict surfaced for review
- **Composite predicates**: declarative YAML rules that combine atoms into higher-order conclusions
- **Unmapped accession review**: surfaces annotation accessions that produce no semantic mapping, prioritized by frequency

V2 runs in **shadow mode** alongside V1. The compatibility layer (`semantic_state_to_predicates`) converts V2 output back to flat predicate lists, so all downstream code (search_by_predicates, reports) continues to work unchanged.

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

# Run the full V2 pipeline: generate atoms, aggregate, evaluate composites, persist
states = b.generate_v2()

# With review queue output (unmapped accessions ranked by frequency)
states = b.generate_v2(output_review_queue="review_queue.tsv")
```

### Query semantic state

```python
# Get the resolved semantic state for a protein
state = b.get_semantic_state("protein_123")
state.activities     # ["hydrogenase", "nife_hydrogenase"]
state.roles          # ["defense_system"]
state.architecture   # ["multi_domain"]
state.size_class     # "giant"
state.composite_predicates  # ["energy_conserving_hydrogenase"]

# Get raw atoms (full provenance: which annotation produced each claim)
atoms = b.get_atoms("protein_123")
for a in atoms:
    print(f"{a.atom_id} ({a.facet.value}): {a.relation.value} via {a.source_db}:{a.source_accession}")
```

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
b.search_by_atoms(has=["defense_system"], limit=5000)   # get all defense proteins
```

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

# V1 vs V2 comparison
diff = b.run_shadow_diff(output_report="shadow_report.md")
print(f"Match rate: {diff['summary']['match_rate_pct']}%")
```

### V1 compatibility

V2 produces flat V1-compatible predicate lists via the compat layer, so `search_by_predicates()` continues to work unchanged:

```python
# V1 predicate search still works
b.search_by_predicates(has=["giant", "unannotated"])
b.search_by_predicates(has=["hydrogenase"], lacks=["hypothetical"])
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
    # V1 compatibility
    semantic_state_to_predicates,
    # Review queue
    build_review_queue, format_review_queue_tsv,
    # Shadow diff
    shadow_diff, run_shadow_diff_on_store,
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
| `crispr_validated` | CRISPR-Cas with validated components |
| `restriction_modification_validated` | RM system confirmed |
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

Maps each V1 predicate ID to a SemanticFacet. 536 entries, grouped by facet.

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

V2 creates two tables (via `create_v2_tables(store)`):

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

Resolved per-protein view.

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

## Testing

```bash
# V2 tests only (131 tests)
python -m pytest tests/test_predicates_v2/ -v

# All tests (excluding pre-existing infrastructure failures)
python -m pytest tests/ --ignore=tests/test_ingest_smoke.py --ignore=tests/test_ingest_resources.py -q
```
