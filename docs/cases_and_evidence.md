# Cases, context comparisons, and optional assembly evidence

Sharur can resolve a protein, caller-emitted system, locus, contig, or bin into
a first-class biological case. A case is both:

- a typed, serializable evidence record; and
- a connected object that can compare contexts, render a locus, compile
  claim-checked findings, and export a replayable evidence bundle.

The implementation is deliberately strict about annotation provenance and
deliberately lazy about expensive data.

## Inspect a case

```python
from sharur.operators import Sharur

b = Sharur("data/DATASET/sharur.duckdb", read_only=True)
case = b.inspect(
    "CALL_OR_PROTEIN_ID",
    entity_type="system",
    window=10,          # symmetric default
    upstream_orfs=4,    # optional side-specific override
    downstream_orfs=8,
)

print(case.to_markdown())
record = case.evidence()
```

The CLI exposes the same operation:

```bash
sharur inspect CALL_ID \
  --type system \
  --upstream 4 \
  --downstream 8 \
  --db data/DATASET/sharur.duckdb
```

`window` defaults to ten ORFs on each side. `upstream_orfs` and
`downstream_orfs` override that default independently. For a co-oriented
anchor, upstream and downstream follow transcriptional orientation. Current
`+`/`-` and legacy `1`/`-1` strand encodings are normalized. Mixed or unknown
anchor orientation falls back to ascending coordinate order and is recorded
as a limitation.

IDs that occur in more than one entity class fail as ambiguous. Pass
`entity_type`, and when needed `source_table` or `bin_id`, rather than letting
the operator guess.

## Evidence boundaries

Case records keep these categories separate:

- `annotations`: per-domain or per-profile observations, labeled `observed`;
- `named_calls`: exact names emitted by a structured caller resource, labeled
  `caller_named`;
- predicates: derived semantic features;
- assembly evidence: optional contig-level measurements from a separate
  sidecar.

Raw system-profile rows are never promoted into named system calls. The
operator discovers the live schema, resolves system membership through
`system_proteins`, and gets a system name from the corresponding structured
caller table. Multi-replicon system membership fails closed for neighborhood
construction.

## Compare ORF contexts

```python
comparison = case.compare_context(
    features=[
        "pfam:PF00589",
        {
            "kind": "other_called_system",
            "max_orfs": 6,
        },
    ],
    window=6,
    upstream_orfs=4,
    downstream_orfs=6,
    combine="all",
    min_components=2,
    require_full_context=True,
    deduplicate_by="replicon",
)

print(comparison.to_markdown())
```

Compact feature forms are:

| Form | Meaning |
|---|---|
| `pfam:PF00589` | Exact observed accession from the named source |
| `annotation:PF00589` | Exact observed accession from any source |
| `name:integrase` | Case-insensitive substring in observed name/description |
| `predicate:mobile_element` | Active semantic predicate |
| `system:RM_Type_III` | Exact structured caller-emitted system type |
| `other_called_system` | Any structured call other than the anchor call |
| `locus:crispr` | Exact structured locus type |

An explicit feature mapping can additionally set:

- `side`: `upstream`, `downstream`, `anchor`, or `either`;
- `max_orfs`: maximum absolute relative ORF displacement, where an adjacent
  flanking gene has displacement one;
- `include_anchor_call`: include the anchor call itself for system features.

For system and locus cases, the default foreground is the same exact
caller-emitted type as the inspected case. The default background is every
other type emitted by that same structured resource. Protein cases require
explicit foreground and background IDs because a single protein does not
define a defensible automatic cohort.

Controls are explicit:

- `require_full_context=True` excludes edge-censored neighborhoods;
- `min_components` filters partial calls;
- `deduplicate_by` chooses `entity`, `replicon`, or `bin` as the counting unit;
- foreground-bearing units are excluded from the background by default;
- `taxonomy_filter` and `same_taxonomy_rank` constrain both groups.

Each feature and the requested `all`/`any` composite gets a 2×2 table, odds
ratio with a 95% Wald interval (Haldane correction for zero cells), exact
Fisher test, and Benjamini–Hochberg q-value. The output retains the complete
entity matrix, denominators, exclusions, parameters, and limitations. These
are locus/replicon observations, not automatically phylogenetically
independent prevalence estimates.

The comparator fetches only requested anchor windows and only annotation
classes needed by the requested features; it does not materialize complete
replicons.

CLI example:

```bash
sharur compare-context CALL_ID \
  --type system \
  --feature pfam:PF00589 \
  --feature other_called_system \
  --combine all \
  --upstream 4 \
  --downstream 6 \
  --min-components 2 \
  --deduplicate-by replicon \
  --db data/DATASET/sharur.duckdb \
  --format json
```

## Optional assembly and host-assignment evidence

Assembly evidence is not part of the canonical Sharur schema. It lives in a
small `assembly_evidence.duckdb` beside `sharur.duckdb`, or at an explicitly
configured path:

```python
b = Sharur(
    "data/DATASET/sharur.duckdb",
    read_only=True,
    assembly_evidence_path="/path/to/assembly_evidence.duckdb",
)
```

If the sidecar is absent, case inspection continues and reports the optional
capability as unavailable. Preflight never creates it:

```bash
sharur preflight --db data/DATASET/sharur.duckdb --format json
```

Import a TSV, CSV, or JSONL containing `bin_id`, `contig_id`, and any of the
recognized scalar fields:

```bash
sharur import-assembly-evidence contig_metrics.tsv \
  --db data/DATASET/sharur.duckdb \
  --source mapping_pipeline_v2
```

Recognized fields include coverage mean/SD/CV, coverage relative to the bin
median, mapped reads, proper-pair fraction, insert-size summaries, SNV
counts/density, assembly-graph degree/component, contig taxonomy and
congruence, GC z-score, and tetranucleotide distance/percentile. Extra input
columns are retained in JSON metadata. Re-imports merge non-null fields and
metadata, so adding composition evidence does not erase coverage or linkage
measurements.

Core-dataset validation and input hashing are enabled by default. They can be
disabled explicitly when the external evidence uses a different namespace or
when avoiding an extra sequential read matters.

### Composition is a separate opt-in operation

Nothing in ingestion, preflight, `inspect`, or `compare_context` scans
assemblies or computes k-mers. To request it:

```bash
sharur compute-composition-evidence \
  --assembly BIN_A=/path/to/BIN_A.fna \
  --assembly BIN_B=/path/to/BIN_B.fna \
  --db data/DATASET/sharur.duckdb
```

This command streams the supplied FASTAs (plain text or gzip-compressed),
canonicalizes reverse-complement
4-mers, compares each contig with a leave-one-contig-out bin aggregate, and
persists only GC z-scores plus scalar cosine-distance/percentile summaries.
The 4-mer vectors exist only in memory for the current assembly and are never
stored by Sharur. A normalized FASTA-content SHA-256 is accumulated during
that same scan, without a second disk read. A single-contig assembly has no
defensible leave-one-out 4-mer reference, so its tetranucleotide distance
remains null; a contig with no canonical A/C/G/T bases also has no GC z-score.

## Claim compiler

Cases can create a finding draft whose claims have explicit levels and exact
evidence references:

```python
draft = case.propose_finding(
    title="Qualified finding title",
    category="defense_systems",
    description="Interpretation with the claim boundary stated explicitly.",
    novelty=2,
    falsification="Wrong if the caller cohort or observed feature is invalid.",
    comparison=comparison,
    claims=[
        {
            "text": "The structured caller emitted the recorded system type.",
            "level": "caller_named",
            "evidence_refs": [f"call:{case.entity_id}"],
        },
        {
            "text": "The context association is compatible with a linked module.",
            "level": "inferred",
            "evidence_refs": ["comparison:composite"],
        },
    ],
)

report = draft.validate_claims()
finding = draft.compile()  # strict by default
```

The compiler blocks unsupported evidence references, named claims labeled as
raw observations, definitive language without experimental evidence,
quantitative claims without sufficient replayable outputs, high-novelty
findings without falsification, and priority language without an attached
literature audit. A draft containing only unverified ideas is not a canonical
finding; keep it in the hypothesis registry until at least one supporting
claim is replayable.

## Compact evidence bundles

```python
path = draft.publish_bundle(
    "data/DATASET/exploration/bundles/case-name",
    include_sequences=True,
)
```

Or directly:

```bash
sharur compare-context CALL_ID \
  --type system \
  --feature pfam:PF00589 \
  --bundle /path/to/bundle \
  --db data/DATASET/sharur.duckdb
```

A bundle contains typed case JSON, a human-readable summary, comparison
matrix and recipe, an executable comparison verifier, a claim-compiled
finding, verification SQL/replay instructions, anchor-component FASTA, and a
hash manifest. It never copies the dataset DuckDB, embeddings, FAISS indexes,
raw reads, assemblies, or k-mer vectors.
