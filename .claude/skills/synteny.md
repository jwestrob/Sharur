# /synteny — Exact ELSA Synteny Access

Run ELSA embedding-based synteny discovery and query its normalized,
run-scoped result sidecar from Sharur. ELSA evidence describes inferred
conservation of gene neighborhoods. Functional system and pathway names still
require an exact call from a live, purpose-built annotation resource.

## Canonical contract

Each Sharur dataset may contain:

```text
data/DATASET/
├── sharur.duckdb                 # Core proteins and annotations
├── dataset.seal.json             # Stable dataset identity
├── synteny.duckdb                # Canonical ELSA result sidecar
└── synteny/
    ├── store/
    │   ├── index.faiss
    │   ├── metadata.parquet
    │   ├── embeddings.npy
    │   └── config.json
    └── results/
        ├── micro_chain_blocks.csv
        └── micro_chain_clusters.csv
```

`synteny.duckdb` is the agent query surface. It carries:

- a deterministic `run_id`, dataset identity, parameters, input hashes, ELSA
  version, commit, mapping version, validation results, and active-run state;
- exact protein IDs and coordinate-sorted ELSA position indices;
- orientation-resolved anchor pairs;
- every block and every cluster, including singleton clusters;
- explicit cluster loci and many-to-many protein membership with `anchor` or
  `context` roles.

Source cluster integers are local labels. Cite them together with `run_id` and
the canonical `cluster_key`.

## Preflight

Inspect the live capability before using ELSA evidence:

```python
from sharur.operators import Sharur

b = Sharur("data/DATASET/sharur.duckdb", read_only=True)
capability = b.capabilities().get("elsa_synteny")
print(capability.state.value)
print(capability.evidence)
```

Proceed when `elsa_synteny` is `available`. A `stale` state identifies schema,
active-run, or dataset-seal drift. An `unavailable` state identifies a dataset
that still needs materialization.

Query methods enforce the dataset-seal comparison. A stale sidecar raises
`SyntenyDatasetMismatchError`, and `inspect()` records `synteny_state="stale"`
with zero attached memberships.

Historical analysis requires an explicit, documented opt-in:

```python
b = Sharur(
    "data/DATASET/sharur.duckdb",
    read_only=True,
    allow_stale_synteny=True,
)
```

Use this only after inspecting the recorded drift and cite the run as
historical. A current biological claim requires a refreshed embedding and
chaining run when the changed records intersect the claim.

## Run ELSA

Prerequisites:

1. Stage 06 `protein_embeddings.h5` exists.
2. ELSA is installed in the active environment.
3. macOS runs export `KMP_DUPLICATE_LIB_OK=TRUE`.
4. One MPS or other heavy ELSA compute process runs at a time.

```bash
elsa synteny \
    --db data/DATASET/sharur.duckdb \
    --embeddings data/DATASET/embeddings/protein_embeddings.h5 \
    --annotations-db data/DATASET/sharur.duckdb \
    --store data/DATASET/synteny/store \
    --result-db data/DATASET/synteny.duckdb \
    --run-label DATASET-production \
    --jobs "$(sysctl -n hw.logicalcpu)" \
    -o data/DATASET/synteny/results/
```

The CLI writes the legacy CSV interchange files and materializes the normalized
sidecar. `--store` preserves the exact gene ordering and embedding index used
by the run.

### Materialize an existing checkpoint

Use the original result directory and its exact store:

```bash
elsa materialize-results data/DATASET/synteny/results/ \
    --store data/DATASET/synteny/store \
    --result-db data/DATASET/synteny.duckdb \
    --dataset-seal data/DATASET/dataset.seal.json \
    --parameters-file data/DATASET/synteny/run_parameters.json \
    --run-label DATASET-production \
    --threads "$(sysctl -n hw.logicalcpu)"
```

The materializer validates unique protein and block IDs, coordinate-ordered
store rows, valid block intervals, anchor-array shape, orientation, exact
anchor resolution onto the declared block side, locus endpoints, and
block-to-cluster coverage before marking a run `ready`.

## Agent query API

### Exact reverse membership

```python
result = b.synteny_for_protein("PROTEIN_ID")
print(result)
memberships = result.raw
```

Each row includes `run_id`, exact `protein_id`, `member_role`, `cluster_key`,
source cluster ID, cluster size, genome support, and a resolved locus.
Membership is many-to-many: retain every returned cluster.

### Batch reverse membership

```python
protein_ids = [
    row[0]
    for row in b.store.execute(
        """
        SELECT DISTINCT protein_id
        FROM annotations
        WHERE source = ? AND accession = ?
        ORDER BY protein_id
        """,
        ["SOURCE", "ACCESSION"],
    )
]

page = b.synteny_for_proteins(protein_ids, limit=500, offset=0)
print(page.meta.total_rows, page.meta.truncated)
memberships = page.raw
```

Paginate until `meta.truncated` is false.

### Exact anchor evidence

```python
anchors = b.synteny_anchor_blocks("PROTEIN_ID", limit=100)
```

This surface returns exact block partners after ELSA orientation has been
resolved.

### Run-scoped cluster expansion

```python
membership = b.synteny_for_protein("PROTEIN_ID").raw[0]
cluster = b.get_synteny_cluster(
    membership["cluster_key"],
    run_id=membership["run_id"],
    member_limit=500,
)
print(cluster)
```

Numeric source IDs are accepted only inside an explicit run:

```python
cluster = b.get_synteny_cluster(
    SOURCE_CLUSTER_ID,
    run_id="elsa-RUN_ID",
)
```

### Case enrichment

```python
case = b.inspect("PROTEIN_ID", entity_type="protein")
print(case.record.synteny_state)
print(case.record.synteny_memberships)
```

`inspect()` attaches bounded, typed ELSA memberships as `inferred` evidence
while keeping observed domains and caller-emitted names in their own evidence
classes.

### Conserved-cluster inventory

```python
import duckdb

sidecar = duckdb.connect(
    "data/DATASET/synteny.duckdb",
    read_only=True,
)
top = sidecar.execute(
    """
    SELECT run_id, cluster_key, source_cluster_id, size, genome_support,
           locus_count, member_count
    FROM current_elsa_clusters
    ORDER BY genome_support DESC, size DESC, cluster_key
    LIMIT 25
    """
).fetchdf()
```

## Citation contract

An ELSA evidence record includes:

```json
{
  "evidence_level": "inferred",
  "run_id": "elsa-RUN_ID",
  "cluster_key": "cluster:SOURCE_ID",
  "source_cluster_id": 123,
  "protein_id": "EXACT_PROTEIN_ID",
  "member_role": "anchor",
  "block_count": 12,
  "genome_support": 8,
  "locus_key": "cluster:SOURCE_ID:locus:0",
  "query": "Sharur.synteny_for_protein exact membership"
}
```

Required identity fields are `run_id`, `cluster_key`, and the exact queried
protein or locus. Support counts come from the same run-scoped row. Cluster
expansion supplies the explicit genomes, loci, and members when the claim
depends on them.

If the capability is unavailable, record: “ELSA synteny validation pending.”

## Legacy CSV quarantine

The CSV files serve checkpoint recovery, interchange, and manual audit. They
are unsuitable as the primary agent lookup surface because they serialize
large JSON values and omit a run namespace.

Legacy handling rules:

- `genes_json` coordinates refer to ELSA’s coordinate-sorted store order.
  Sharur `proteins.gene_index` belongs to a separate ingest contract. Resolve
  legacy coordinates through the exact store `metadata.parquet`, then verify
  the protein ID.
- Substring, regex, `LIKE`, and `str.contains` membership scans can match
  neighboring identifiers. Exact relational membership is the accepted
  lookup.
- Scalar `protein_id -> cluster_id` maps discard valid memberships. Use a
  one-to-many collection.
- Raw query and target anchor arrays are independently sorted by genomic
  coordinate. For inverted blocks, query position `i` pairs with target
  position `n - 1 - i`.

For a raw-block audit, use the validated parser:

```python
from sharur.elsa_blocks import elsa_anchor_pairs_from_block

pairs = elsa_anchor_pairs_from_block(block_row)
```

It accepts both CSV and parquet anchor-column names and rejects malformed
array lengths or orientations.

## Adding genomes

```bash
elsa synteny \
    --store data/DATASET/synteny/store \
    --add-db data/NEW_GENOMES/sharur.duckdb \
    --add-embeddings data/NEW_GENOMES/embeddings/protein_embeddings.h5 \
    --result-db data/DATASET/synteny.duckdb \
    --run-label DATASET-plus-new-genomes \
    --jobs "$(sysctl -n hw.logicalcpu)" \
    -o data/DATASET/synteny/results_combined/
```

The resulting run receives a new deterministic namespace and becomes the
active run after validation.
