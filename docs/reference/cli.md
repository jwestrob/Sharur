# Command-line reference


Generated from the installed commands by `scripts/gen_cli_docs.py`. Every command also prints this text with `--help`.


## `sharur`

Sharur - Metagenomic dataset exploration CLI

**Usage**:

```console
$ sharur [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--version`: Show the Sharur version and exit.
* `--help`: Show this message and exit.

**Commands**:

* `overview`: Show dataset overview with summary...
* `genomes`: List genomes (MAGs) with optional filtering.
* `proteins`: List proteins with optional filtering.
* `neighborhood`: Show genomic neighborhood around a protein.
* `inspect`: Resolve an entity into a...
* `compare-context`: Run an exact, reproducible...
* `import-assembly-evidence`: Import optional scalar contig evidence...
* `compute-composition-evidence`: Explicitly scan FASTAs for scalar GC/4-mer...
* `search`: Search proteins by predicates or annotations.
* `compute-predicates`: Compute and store V2 predicates for proteins.
* `predicates`: List available predicates for search.
* `preflight`: Emit one typed dataset/runtime capability...
* `seal`: Write a portable integrity seal for a...
* `migrate`: Apply pending additive schema/index...
* `backfill-contig-context`: Record assembly contig lengths and...
* `verify-seal`: Recompute a dataset seal and report...
* `build-vector-index`: Build mmap-ready FAISS sidecars and a...
* `doctor`: Verify external tools, reference...
* `architecture`: Find proteins whose ordered domains match...
* `describe`: What a dataset holds: annotation sources,...
* `card`: Summarize one protein: context,...
* `why`: Explain why a protein carries a predicate:...
* `modules`: KEGG module completeness per genome, or...
* `setup-kegg`: Fetch KEGG data and build the KO ->...

### `sharur overview`

Show dataset overview with summary statistics.

Displays genome/protein counts, annotation coverage,
taxonomy distribution, and predicate summary.

**Usage**:

```console
$ sharur overview [OPTIONS]
```

**Options**:

* `-d, --db TEXT`: Path to DuckDB database  [default: data/sharur.duckdb]
* `-f, --format [markdown|json|jsonl|tsv]`: Output format: markdown, json, jsonl, or tsv  [default: markdown]
* `--help`: Show this message and exit.

### `sharur genomes`

List genomes (MAGs) with optional filtering.

Examples:
    sharur genomes --taxonomy Archaea
    sharur genomes --min-comp 90 --max-contam 5

**Usage**:

```console
$ sharur genomes [OPTIONS]
```

**Options**:

* `-t, --taxonomy TEXT`: Filter by taxonomy substring
* `--min-comp FLOAT`: Minimum completeness %
* `--max-contam FLOAT`: Maximum contamination %
* `-n, --limit INTEGER`: Maximum results  [default: 20]
* `-d, --db TEXT`: Path to DuckDB database  [default: data/sharur.duckdb]
* `-f, --format [markdown|json|jsonl|tsv]`: [default: markdown]
* `--help`: Show this message and exit.

### `sharur proteins`

List proteins with optional filtering.

Examples:
    sharur proteins --genome bin_001
    sharur proteins --min-length 2000 --unannotated

**Usage**:

```console
$ sharur proteins [OPTIONS]
```

**Options**:

* `-g, --genome TEXT`: Filter by genome (bin_id)
* `-c, --contig TEXT`: Filter by contig
* `--min-length INTEGER`: Minimum length (aa)
* `--max-length INTEGER`: Maximum length (aa)
* `--annotated / --unannotated`: Filter by annotation status
* `-n, --limit INTEGER`: Maximum results  [default: 50]
* `-d, --db TEXT`: Path to DuckDB database  [default: data/sharur.duckdb]
* `-f, --format [markdown|json|jsonl|tsv]`: [default: markdown]
* `--help`: Show this message and exit.

### `sharur neighborhood`

Show genomic neighborhood around a protein.

Displays proteins in context with coordinates, annotations,
and predicates as an ASCII table.

Example:
    sharur neighborhood PROTEIN_ID --window 15

**Usage**:

```console
$ sharur neighborhood [OPTIONS] PROTEIN_ID
```

**Arguments**:

* `PROTEIN_ID`: Protein ID as anchor  [required]

**Options**:

* `-w, --window INTEGER`: Genes on each side  [default: 10]
* `-v, --verbose`: Show predicate details
* `-d, --db TEXT`: Path to DuckDB database  [default: data/sharur.duckdb]
* `-f, --format [markdown|json|jsonl|tsv]`: [default: markdown]
* `--help`: Show this message and exit.

### `sharur inspect`

Resolve an entity into a provenance-separated, strand-aware case.

**Usage**:

```console
$ sharur inspect [OPTIONS] ENTITY_ID
```

**Arguments**:

* `ENTITY_ID`: Protein, system, locus, contig, or bin ID  [required]

**Options**:

* `--type TEXT`: Disambiguate as protein, system, locus, contig, or bin.
* `--bin TEXT`: Required only when a contig label occurs in multiple bins.
* `--source-table TEXT`: Structured caller table when a system/locus ID is ambiguous.
* `-w, --window INTEGER`: Default ORFs on each side.  [default: 10]
* `--upstream INTEGER`: Biological upstream ORFs; overrides --window on that side.
* `--downstream INTEGER`: Biological downstream ORFs; overrides --window on that side.
* `--assembly-evidence PATH`: Optional assembly_evidence.duckdb sidecar.
* `--include-sequences`: Embed context sequences in JSON output (can make output large).
* `--plot PATH`: Render the resolved locus to this PNG/SVG path.
* `--bundle PATH`: Write a compact, replayable evidence-bundle directory.
* `--bundle-sequences / --no-bundle-sequences`: Include anchor-component FASTA in --bundle.  [default: bundle-sequences]
* `--overwrite`: Replace an existing --bundle directory.
* `-d, --db TEXT`: Path to DuckDB database  [default: data/sharur.duckdb]
* `-f, --format [markdown|json]`: Output format: markdown or json.  [default: markdown]
* `--help`: Show this message and exit.

### `sharur compare-context`

Run an exact, reproducible foreground/background ORF-context comparison.

**Usage**:

```console
$ sharur compare-context [OPTIONS] ENTITY_ID
```

**Arguments**:

* `ENTITY_ID`: System, locus, or protein case ID  [required]

**Options**:

* `--feature TEXT`: Repeatable feature: pfam:ACCESSION, name:TEXT, predicate:ID, system:TYPE, locus:TYPE, or other_called_system.
* `--type TEXT`
* `--source-table TEXT`
* `--bin TEXT`
* `-w, --window INTEGER`: Default ORFs on each side.  [default: 10]
* `--upstream INTEGER`: Biological upstream ORFs; overrides --window.
* `--downstream INTEGER`: Biological downstream ORFs; overrides --window.
* `--foreground-id TEXT`: Explicit foreground entity ID; repeat for protein/custom cohorts.
* `--background-id TEXT`: Explicit background entity ID; repeat for protein/custom cohorts.
* `--combine TEXT`: Combine features with all or any.  [default: all]
* `--min-components INTEGER`: Exclude caller-emitted systems/loci with fewer components.  [default: 1]
* `--require-full-context / --allow-edge-censored`: Exclude contig-edge-censored neighborhoods.  [default: require-full-context]
* `--deduplicate-by TEXT`: Independent unit: entity, replicon, or bin.  [default: replicon]
* `--exclude-foreground-units / --allow-foreground-overlap`: Keep foreground-bearing units out of the background.  [default: exclude-foreground-units]
* `--taxonomy TEXT`: Require this taxonomy substring in both groups.
* `--same-taxonomy-rank TEXT`: Match the case at domain/phylum/class/order/family/genus/species.
* `--alternative TEXT`: Fisher alternative: greater, less, or two-sided.  [default: greater]
* `--bundle PATH`: Write the case, comparison matrix, recipe, and verifier.
* `--overwrite`
* `-d, --db TEXT`: Path to DuckDB database  [default: data/sharur.duckdb]
* `-f, --format [markdown|json]`: [default: markdown]
* `--help`: Show this message and exit.

### `sharur import-assembly-evidence`

Import optional scalar contig evidence without modifying the core database.

**Usage**:

```console
$ sharur import-assembly-evidence [OPTIONS] INPUT_PATH
```

**Arguments**:

* `INPUT_PATH`: TSV, CSV, or JSONL with bin_id and contig_id columns.  [required]

**Options**:

* `-d, --db PATH`: Core DuckDB used for validation and sidecar discovery.  [default: data/sharur.duckdb]
* `--sidecar PATH`: Output sidecar (default: assembly_evidence.duckdb beside --db).
* `--source TEXT`: Provenance label for these measurements.
* `--validate / --no-validate`: Require every (bin_id, contig_id) to exist in the core dataset.  [default: validate]
* `--hash-input / --skip-input-hash`: SHA-256 the input for provenance (one extra sequential read).  [default: hash-input]
* `--help`: Show this message and exit.

### `sharur compute-composition-evidence`

Explicitly scan FASTAs for scalar GC/4-mer evidence.

This command is never run by ingestion, preflight, or case inspection.
It reads every supplied FASTA, keeps 4-mer vectors in memory only, and
persists scalar leave-one-contig-out distances.

**Usage**:

```console
$ sharur compute-composition-evidence [OPTIONS]
```

**Options**:

* `--assembly TEXT`: Repeatable BIN_ID=/path/to/assembly.fna mapping.
* `-d, --db PATH`: Core DuckDB used for validation and sidecar discovery.  [default: data/sharur.duckdb]
* `--sidecar PATH`: Output sidecar (default: assembly_evidence.duckdb beside --db).
* `--validate / --no-validate`: Require FASTA contigs to exist in the core dataset.  [default: validate]
* `--help`: Show this message and exit.

### `sharur search`

Search proteins by predicates or annotations.

Predicate search uses set logic:
- --has: Protein must have ALL specified predicates (AND)
- --lacks: Protein must have NONE of specified predicates

Examples:
    sharur search --has giant,unannotated
    sharur search --has confident_hit --lacks hypothetical
    sharur search --annotation "hydrogenase"
    sharur search --accession PF00142

**Usage**:

```console
$ sharur search [OPTIONS]
```

**Options**:

* `--has TEXT`: Predicates that must be true (comma-separated)
* `--lacks TEXT`: Predicates that must be false (comma-separated)
* `-a, --annotation TEXT`: Annotation pattern to match
* `--accession TEXT`: Exact accession (e.g., PF00142)
* `-t, --taxonomy TEXT`: Taxonomy filter
* `-n, --limit INTEGER`: Maximum results  [default: 50]
* `-d, --db TEXT`: Path to DuckDB database  [default: data/sharur.duckdb]
* `-f, --format [markdown|json|jsonl|tsv]`: [default: markdown]
* `--help`: Show this message and exit.

### `sharur compute-predicates`

Compute and store V2 predicates for proteins.

This writes semantic_atoms and semantic_state, then materializes the
V2-derived protein_predicates compatibility table used by legacy
search and reports.

Run this after loading annotations to enable predicate-based search.

Examples:
    sharur compute-predicates --db data/sharur.duckdb
    sharur compute-predicates --protein PROT_001 --db data/sharur.duckdb

**Usage**:

```console
$ sharur compute-predicates [OPTIONS]
```

**Options**:

* `-d, --db TEXT`: Path to DuckDB database  [default: data/sharur.duckdb]
* `-p, --protein TEXT`: Compute for specific protein only
* `--chunk-size INTEGER`: V2 generation batch size  [default: 100000]
* `-t, --workers INTEGER`: Transform processes (defaults to Slurm CPUs, then host CPUs)
* `--worker-batch-size INTEGER`: Optional proteins per transform task
* `--pipeline-depth INTEGER RANGE`: Bounded full-dataset chunks overlapped across read/transform/write  [default: 2; x>=1]
* `--resume`: Continue the latest full V2 checkpoint
* `--review-queue PATH`: Write unresolved-accession review TSV
* `--help`: Show this message and exit.

### `sharur predicates`

List available predicates for search.

Predicates are boolean properties computed over proteins.
Use with 'sharur search --has PREDICATE'.

Categories include:
- enzyme: Enzyme classes (oxidoreductase, hydrolase, etc.)
- transport: Transporters and membrane proteins
- regulation: Regulators and signaling
- metabolism: Metabolic pathway enzymes
- cazy: Carbohydrate-active enzymes
- binding: Binding domains
- envelope: Cell surface and envelope
- mobile: Mobile elements and defense
- stress: Stress response and resistance
- size: Size-based (tiny, giant, massive)
- annotation: Annotation status

Examples:
    sharur predicates
    sharur predicates --category enzyme
    sharur predicates --category transport --hierarchy

**Usage**:

```console
$ sharur predicates [OPTIONS]
```

**Options**:

* `-c, --category TEXT`: Filter by category
* `-h, --hierarchy`: Show parent predicates
* `--help`: Show this message and exit.

### `sharur preflight`

Emit one typed dataset/runtime capability brief without modifying data.

**Usage**:

```console
$ sharur preflight [OPTIONS]
```

**Options**:

* `-d, --db TEXT`: Path to DuckDB database  [default: data/sharur.duckdb]
* `--assembly-evidence PATH`: Optional non-default assembly-evidence sidecar.
* `--synteny PATH`: Optional non-default normalized ELSA sidecar.
* `-f, --format [markdown|json]`: Output format: markdown or json  [default: markdown]
* `--skip-tools`: Skip slower external binary/version probes.
* `--strict`: Exit non-zero unless every required dataset capability is available.
* `--help`: Show this message and exit.

### `sharur seal`

Write a portable integrity seal for a completed dataset.

**Usage**:

```console
$ sharur seal [OPTIONS]
```

**Options**:

* `-d, --db PATH`: Path to the dataset DuckDB.  [default: data/sharur.duckdb]
* `-o, --output PATH`: Seal path (default: DATASET/dataset.seal.json).
* `--full`: Fully hash every canonical artifact, including large DuckDB/H5/index files.
* `--include-tools`: Record slower external tool and reference-database probes as provenance.
* `--force`: Replace an existing seal atomically.
* `--help`: Show this message and exit.

### `sharur migrate`

Apply pending additive schema/index migrations in a maintenance window.

**Usage**:

```console
$ sharur migrate [OPTIONS]
```

**Options**:

* `-d, --db PATH`: Path to a writable Sharur DuckDB.  [default: data/sharur.duckdb]
* `--help`: Show this message and exit.

### `sharur backfill-contig-context`

Record assembly contig lengths and Prodigal truncation flags in an existing database.

Older databases store each contig's length as the end of its last gene and
omit Prodigal's partial-gene flags; contig-edge context needs both.
Writes the canonical database: run in a maintenance window, then reseal.

**Usage**:

```console
$ sharur backfill-contig-context [OPTIONS]
```

**Options**:

* `-d, --db PATH`: Path to a writable Sharur DuckDB.  [default: data/sharur.duckdb]
* `--assemblies PATH`: Directory of assembly FASTAs named BIN_ID.fna|fa|fasta[.gz] (default: the dataset's stage00_prepared manifest).
* `--proteins PATH`: Prodigal output directory with *.faa files (default: the dataset's stage03_prodigal).
* `--help`: Show this message and exit.

### `sharur verify-seal`

Recompute a dataset seal and report canonical identity drift.

**Usage**:

```console
$ sharur verify-seal [OPTIONS] SEAL_PATH
```

**Arguments**:

* `SEAL_PATH`: Path to dataset.seal.json.  [required]

**Options**:

* `-d, --db PATH`: Current DuckDB path (needed if the seal was moved separately).
* `-f, --format [markdown|json]`: Output format: markdown or json.  [default: markdown]
* `--help`: Show this message and exit.

### `sharur build-vector-index`

Build mmap-ready FAISS sidecars and a stable disk-backed protein-ID map.

**Usage**:

```console
$ sharur build-vector-index [OPTIONS]
```

**Options**:

* `-e, --embeddings PATH`: Canonical protein_embeddings.h5 path.
* `-d, --db PATH`: DuckDB path; discovers embeddings beside the dataset.
* `--force`: Rebuild an already-valid index.
* `--chunk-size INTEGER`: Number of vectors streamed per FAISS add batch.  [default: 50000]
* `--nprobe INTEGER`: IVF partitions probed per query.  [default: 32]
* `-t, --threads INTEGER`: FAISS CPU build threads (default: FAISS runtime default).
* `--help`: Show this message and exit.

### `sharur doctor`

Verify external tools, reference databases, and API keys are available.

Default is informational (always exits 0) so it is safe as a CI smoke check
and container HEALTHCHECK. Pass --strict to fail when a core component is
missing.

**Usage**:

```console
$ sharur doctor [OPTIONS]
```

**Options**:

* `--strict`: Exit non-zero if any core tool/database is missing.
* `--help`: Show this message and exit.

### `sharur architecture`

Find proteins whose ordered domains match a pattern.

Tokens: domain names or accessions, globs (Big_*), '.' for any domain,
quantifiers (? * + {m} {m,} {m,n}; put a space before ? and *), groups
with alternatives ( A | B ), and anchors ^ (N-terminus) and $ (C-terminus).

**Usage**:

```console
$ sharur architecture [OPTIONS] PATTERN
```

**Arguments**:

* `PATTERN`: Domain pattern, e.g. "Big_2{5,} . VWA" or "^ SP ( Cadherin | Big_2 )+"  [required]

**Options**:

* `-d, --db TEXT`: Path to DuckDB database  [default: data/sharur.duckdb]
* `-s, --source TEXT`: Annotation source(s); repeatable.  [default: pfam]
* `-b, --bin TEXT`: Restrict to genome(s); repeatable.
* `-n, --limit INTEGER`: Records to show (the total is always counted).  [default: 50]
* `--max-overlap FLOAT`: Overlap (fraction of the shorter domain) above which the weaker hit is dropped.  [default: 0.5]
* `-f, --format [markdown|json]`: markdown or json  [default: markdown]
* `--help`: Show this message and exit.

### `sharur describe`

What a dataset holds: annotation sources, curated callers, genome metadata, predicate state.

**Usage**:

```console
$ sharur describe [OPTIONS]
```

**Options**:

* `-d, --db TEXT`: Path to DuckDB database  [default: data/sharur.duckdb]
* `-f, --format [markdown|json]`: markdown or json  [default: markdown]
* `--help`: Show this message and exit.

### `sharur card`

Summarize one protein: context, annotations, evidence-backed predicates (no sequences).

**Usage**:

```console
$ sharur card [OPTIONS] PROTEIN_ID
```

**Arguments**:

* `PROTEIN_ID`: Protein ID  [required]

**Options**:

* `-d, --db TEXT`: Path to DuckDB database  [default: data/sharur.duckdb]
* `-w, --window INTEGER`: Neighbors on each side.  [default: 5]
* `-f, --format [markdown|json]`: markdown or json  [default: markdown]
* `--help`: Show this message and exit.

### `sharur why`

Explain why a protein carries a predicate: annotation hits, map evidence, expansion chain.

**Usage**:

```console
$ sharur why [OPTIONS] PROTEIN_ID PREDICATE
```

**Arguments**:

* `PROTEIN_ID`: Protein ID  [required]
* `PREDICATE`: Predicate ID, e.g. nad_binding  [required]

**Options**:

* `-d, --db TEXT`: Path to DuckDB database  [default: data/sharur.duckdb]
* `-f, --format [markdown|json]`: markdown or json  [default: markdown]
* `--help`: Show this message and exit.

### `sharur modules`

KEGG module completeness per genome, or around one protein (needs `sharur setup-kegg`).

**Usage**:

```console
$ sharur modules [OPTIONS]
```

**Options**:

* `-d, --db TEXT`: Path to DuckDB database  [default: data/sharur.duckdb]
* `-b, --bin TEXT`: Genome(s); repeatable.
* `-m, --module TEXT`: Module(s), e.g. M00175; repeatable.
* `--min-completeness FLOAT`: Report modules at or above this fraction.  [default: 0.0]
* `--around TEXT`: Protein ID: modules encoded in its neighborhood.
* `-w, --window INTEGER`: Genes on each side for --around.  [default: 10]
* `-f, --format [markdown|json]`: markdown or json  [default: markdown]
* `--help`: Show this message and exit.

### `sharur setup-kegg`

Fetch KEGG data and build the KO -> predicate map on this machine.

KEGG data is subject to KEGG's terms: KEGG REST (used by default, about 100
requests at under 3 per second) is for academic use; other users need a KEGG
license and can build from their licensed files with --inputs. The built
files stay local; Sharur ships only its own rules.

**Usage**:

```console
$ sharur setup-kegg [OPTIONS]
```

**Options**:

* `--dir PATH`: Where to build (default: $SHARUR_KEGG_DIR, else data/reference/kegg).
* `--inputs PATH`: Build from KEGG files already on disk (list/ko as ko_list.tsv, brite/<id>.json, modules.txt) instead of fetching, e.g. from a licensed KEGG copy.
* `--kofam-ko-list PATH`: KOfam ko_list (downloaded when absent).  [default: data/reference/ko_list]
* `--report PATH`: TSV of proposed pairs without evidence.
* `--help`: Show this message and exit.


## `sharur-ingest`

Run the staged ingest pipeline from the primary CLI entrypoint.

`mode=tools` is the default and mirrors the documented standard pipeline.
`mode=fast` exists only for local smoke tests and synthetic fixtures.

For the user-facing workflow, see QUICKSTART.md or src/ingest/README.md.

**Usage**:

```console
$ sharur-ingest [OPTIONS]
```

**Options**:

* `-i, --input-dir PATH`: Directory of genome FASTAs (.fna/.fa/.fasta)  [default: dummy_dataset]
* `-d, --data-dir PATH`: Dataset root where stageXX outputs will be written  [default: data]
* `-o, --output PATH`: Destination DuckDB path, usually DATASET/sharur.duckdb  [default: data/sharur.duckdb]
* `-m, --mode TEXT`: tools = run the standard staged pipeline; fast = synthetic smoke mode  [default: tools]
* `--force`: Remove existing stage outputs and rebuild
* `--skip-quast / --with-quast`: Skip optional Stage 01 by default; use --with-quast to enable it  [default: skip-quast]
* `--skip-dfast / --with-dfast`: Skip optional Stage 02 by default; use --with-dfast to enable it  [default: skip-dfast]
* `--skip-prodigal / --no-skip-prodigal`: Skip standard Stage 03 (Prodigal)  [default: no-skip-prodigal]
* `--skip-astra / --no-skip-astra`: Skip standard Stage 04 (Astra via 04_astra_scan.py)  [default: no-skip-astra]
* `--skip-gecco / --with-gecco`: Skip optional Stage 05a by default; use --with-gecco to enable it  [default: skip-gecco]
* `--skip-dbcan / --with-legacy-dbcan`: Skip deprecated Stage 05b by default; use --with-legacy-dbcan to enable it  [default: skip-dbcan]
* `--enable-cazymes`: Run the optional Stage 07 dbCAN three-tool consensus classifier
* `--pipeline-depth INTEGER RANGE`: Maximum ordered V2 transform chunks in flight during Stage 07  [default: 2; x>=1]
* `--skip-crispr / --no-skip-crispr`: Skip standard Stage 05c (CRISPR arrays via minced_crispr.py)  [default: no-skip-crispr]
* `--skip-embeddings / --no-skip-embeddings`: Skip Stage 06 post-build embeddings (06_esm2_embeddings.py)  [default: no-skip-embeddings]
* `--embedding-model TEXT`: Hugging Face protein encoder for Stage 06 (default: facebook/esm2_t6_8M_UR50D)
* `--profile TEXT`: Execution profile: auto, local (CPU), mps, or slurm  [default: auto]
* `--resume / --no-resume`: Reuse only ledger-verified stages with matching signatures and outputs  [default: resume]
* `--submit-slurm`: Submit the generated SLURM bundle; otherwise only write it
* `--run-idempotency-key TEXT`: Optional caller-stable key that deduplicates identical run creation
* `--dry-run / --no-dry-run`: Print the planned stage commands without executing them  [default: no-dry-run]
* `--help`: Show this message and exit.


## `sharur-atlas`

Plan, enqueue, and verify exhaustive genome-owned Atlas reading.

**Usage**:

```console
$ sharur-atlas [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `plan`: Build stable one-genome work units from a...
* `packet-census`: Count exact bin-scoped packets and payload...
* `verify-packet-census`: Verify the zero-model-call launch gate and...
* `enqueue`: Create an idempotent Ops campaign and one...
* `verify-coverage`: Validate every assigned genome and exact...
* `triage`: Reduce a candidate set to a short digest,...

### `sharur-atlas plan`

Build stable one-genome work units from a sealed DuckDB.

**Usage**:

```console
$ sharur-atlas plan [OPTIONS]
```

**Options**:

* `--db FILE`: [required]
* `--output-dir PATH`: [required]
* `--packet-bytes INTEGER RANGE`: Required model-facing canonical payload budget; derive it from the executor context budget.  [1024<=x<=1500000; required]
* `--seal FILE`
* `--packet-contigs INTEGER RANGE`: Diagnostic hard cap; omitted values use the byte-proportional schema ceiling.  [1<=x<=7211]
* `--packet-proteins INTEGER RANGE`: Diagnostic hard cap; omitted values use the byte-proportional schema ceiling.  [1<=x<=9375]
* `--all-annotations / --top-annotation-only`: [default: all-annotations]
* `--calibration-genomes INTEGER RANGE`: Explicit calibration size; omitted values use ceil(sqrt(live bin count)).  [x>=1]
* `--checkpoint-interval-frames INTEGER RANGE`: [default: 1; x>=1]
* `--query-result-bytes INTEGER RANGE`: [default: 2097152; x>=4097]
* `--threads INTEGER RANGE`: [default: 4; x>=1]
* `--verify-seal / --skip-seal-verification`: [default: verify-seal]
* `--help`: Show this message and exit.

### `sharur-atlas packet-census`

Count exact bin-scoped packets and payload sizes with zero model calls.

**Usage**:

```console
$ sharur-atlas packet-census [OPTIONS]
```

**Options**:

* `--plan-dir DIRECTORY`: [required]
* `--output-dir PATH`
* `--workers INTEGER RANGE`: [default: 4; x>=1]
* `--threads INTEGER RANGE`: [default: 4; x>=1]
* `--memory-limit TEXT`: [default: 16GB]
* `--temp-directory PATH`
* `--max-temp-size TEXT`: [default: 128GB]
* `--resume / --recompute`: [default: resume]
* `--verify-seal / --skip-seal-verification`: [default: verify-seal]
* `--help`: Show this message and exit.

### `sharur-atlas verify-packet-census`

Verify the zero-model-call launch gate and optional unit records.

**Usage**:

```console
$ sharur-atlas verify-packet-census [OPTIONS]
```

**Options**:

* `--plan-dir DIRECTORY`: [required]
* `--census-dir DIRECTORY`
* `--deep / --summary-only`: [default: summary-only]
* `--help`: Show this message and exit.

### `sharur-atlas enqueue`

Create an idempotent Ops campaign and one task per genome.

**Usage**:

```console
$ sharur-atlas enqueue [OPTIONS]
```

**Options**:

* `--plan-dir DIRECTORY`: [required]
* `--query-url TEXT`: [required]
* `--ops-url TEXT`: [default: http://localhost:8811]
* `--agent-id TEXT`: [default: atlas-coordinator]
* `--api-token TEXT`: [env var: SHARUR_OPS_TOKEN]
* `--priority INTEGER RANGE`: [default: 1; 0<=x<=3]
* `--max-attempts INTEGER RANGE`: [default: 5; x>=1]
* `--lease-seconds INTEGER RANGE`: [default: 900; x>=1]
* `--scan-execution-profile TEXT`: [default: atlas_scan]
* `--help`: Show this message and exit.

### `sharur-atlas verify-coverage`

Validate every assigned genome and exact contig/protein totals.

**Usage**:

```console
$ sharur-atlas verify-coverage [OPTIONS]
```

**Options**:

* `--plan-dir DIRECTORY`: [required]
* `--coverage-dir DIRECTORY`
* `--help`: Show this message and exit.

### `sharur-atlas triage`

Reduce a candidate set to a short digest, with no model calls.

Census removal, nested-cluster folding, confound flagging and interest
ranking are all deterministic. Read the digest directly, or hand it to a
model for the one step that needs judgement -- reading the evidence prose.
Feeding raw occurrences to a model instead pays per token to recompute what
a query already knows.

**Usage**:

```console
$ sharur-atlas triage [OPTIONS]
```

**Options**:

* `--ops-db FILE`: [required]
* `--top INTEGER`: Groups to display  [default: 25]
* `--evidence INTEGER`: Evidence snippets per cluster  [default: 0]
* `--min-genomes INTEGER`: Drop groups below this genome count  [default: 2]
* `--system TEXT`: Restrict to one system
* `--json-out PATH`
* `--help`: Show this message and exit.


## `sharur-review`

Operate Sharur's typed, policy-driven scientific review DAG.

**Usage**:

```console
$ sharur-review [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `policy-check`: Validate a policy and print its immutable...
* `reduce`: Reduce exact typed signatures into...
* `route`: Consume durable events and create...
* `verify`: Execute and append one bounded read-only...
* `trace`: Reconstruct a bounded...
* `status`: Report exact funnel, audit, verification,...

### `sharur-review policy-check`

Validate a policy and print its immutable execution contract.

**Usage**:

```console
$ sharur-review policy-check [OPTIONS]
```

**Options**:

* `--policy FILE`
* `--help`: Show this message and exit.

### `sharur-review reduce`

Reduce exact typed signatures into lossless versioned clusters.

**Usage**:

```console
$ sharur-review reduce [OPTIONS]
```

**Options**:

* `--campaign-id TEXT`: [required]
* `--ops-db FILE`
* `--ops-url TEXT`
* `--agent-id TEXT`: [default: review-reducer]
* `--api-token TEXT`: [env var: SHARUR_OPS_TOKEN]
* `--dataset-id TEXT`
* `--candidate-type TEXT`
* `--batch-size INTEGER RANGE`: [default: 1000; 1<=x<=10000]
* `--help`: Show this message and exit.

### `sharur-review route`

Consume durable events and create idempotent review work.

**Usage**:

```console
$ sharur-review route [OPTIONS]
```

**Options**:

* `--campaign-id TEXT`: [required]
* `--ops-db FILE`
* `--ops-url TEXT`
* `--agent-id TEXT`: [default: review-controller]
* `--api-token TEXT`: [env var: SHARUR_OPS_TOKEN]
* `--policy FILE`
* `--watch / --once`: [default: once]
* `--interval-seconds FLOAT RANGE`: [default: 2.0; 0.1<=x<=60.0]
* `--help`: Show this message and exit.

### `sharur-review verify`

Execute and append one bounded read-only DuckDB verification.

**Usage**:

```console
$ sharur-review verify [OPTIONS]
```

**Options**:

* `--review-id TEXT`: [required]
* `--claim-key TEXT`: [required]
* `--db FILE`: [required]
* `--dataset-id TEXT`: [required]
* `--specification TEXT`: Inline JSON or a path to a YAML/JSON verification spec.  [required]
* `--expected TEXT`: Inline JSON or a path to a YAML/JSON expected value.  [required]
* `--ops-db FILE`
* `--ops-url TEXT`
* `--agent-id TEXT`: [default: review-verifier]
* `--api-token TEXT`: [env var: SHARUR_OPS_TOKEN]
* `--seal FILE`
* `--verify-seal / --skip-seal-verification`: [default: verify-seal]
* `--threads INTEGER RANGE`: [default: 1; x>=1]
* `--code-commit TEXT`
* `--supersedes TEXT`
* `--help`: Show this message and exit.

### `sharur-review trace`

Reconstruct a bounded candidate-to-publication provenance graph.

**Usage**:

```console
$ sharur-review trace [OPTIONS]
```

**Options**:

* `--campaign-id TEXT`: [required]
* `--subject-kind TEXT`: candidate_cluster, finding, or unit_disposition  [required]
* `--subject-id TEXT`: [required]
* `--ops-db FILE`
* `--ops-url TEXT`
* `--agent-id TEXT`: [default: review-tracer]
* `--api-token TEXT`: [env var: SHARUR_OPS_TOKEN]
* `--help`: Show this message and exit.

### `sharur-review status`

Report exact funnel, audit, verification, and queue metrics.

**Usage**:

```console
$ sharur-review status [OPTIONS]
```

**Options**:

* `--campaign-id TEXT`: [required]
* `--ops-db FILE`
* `--ops-url TEXT`
* `--agent-id TEXT`: [default: review-observer]
* `--api-token TEXT`: [env var: SHARUR_OPS_TOKEN]
* `--help`: Show this message and exit.


## `sharur-worker`

Sharur model-worker executors

**Usage**:

```console
$ sharur-worker [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--install-completion`: Install completion for the current shell.
* `--show-completion`: Show completion for the current shell, to copy it or customize the installation.
* `--help`: Show this message and exit.

**Commands**:

* `atlas-scan`: Claim and execute `atlas_genome_read` tasks.
* `scientific-review`: Claim review-tier tasks, run their checks,...

### `sharur-worker atlas-scan`

Claim and execute `atlas_genome_read` tasks.

`--dry-run` exercises the full claim -> packet -> coverage -> disposition ->
complete path with zero model calls, which is the right smoke test before
spending any subscription budget.

**Usage**:

```console
$ sharur-worker atlas-scan [OPTIONS]
```

**Options**:

* `--ops-url TEXT`: Sharur Ops base URL  [default: http://127.0.0.1:8811]
* `--query-url TEXT`: Sharur Query base URL  [default: http://127.0.0.1:8812]
* `--agent-id TEXT`: Distinct agent identity for this worker  [required]
* `--profile TEXT`: Execution profile from the review policy  [default: atlas_scan]
* `--campaign-id TEXT`: Restrict claims to one campaign
* `--policy TEXT`: Path to a review policy YAML (default: packaged)
* `--lease-seconds INTEGER RANGE`: Lease duration per claim. A background keepalive renews it every third of this interval, so a long frame cannot outlive its lease.  [default: 2400; x>=60]
* `--stall-timeout INTEGER RANGE`: Max SILENCE tolerated from the model CLI before it is treated as a stalled connection. There is no total runtime cap: a call that keeps emitting events runs as long as it needs.  [default: 900; x>=60]
* `--max-tasks INTEGER RANGE`: Exit after N completed genomes  [x>=1]
* `--max-frames INTEGER RANGE`: Stop each genome after N frames (smoke tests)  [x>=1]
* `--idle-sleep FLOAT`: Seconds to sleep when the queue is empty  [default: 5.0]
* `--sweep-failed / --no-sweep-failed`: When the queue drains, requeue attempt-exhausted tasks whose failure looks like transport  [default: sweep-failed]
* `--max-sweeps INTEGER RANGE`: Backstop on how many sweep rounds a worker performs  [default: 20; x>=0]
* `--dry-run`: Walk packets and build coverage without calling any model (zero model calls)
* `-v, --verbose`
* `--help`: Show this message and exit.

### `sharur-worker scientific-review`

Claim review-tier tasks, run their checks, and append verified reviews.

**Usage**:

```console
$ sharur-worker scientific-review [OPTIONS]
```

**Options**:

* `-d, --db FILE`: Sealed Sharur DuckDB used for executable review checks  [required]
* `--ops-url TEXT`: Sharur Ops base URL  [default: http://127.0.0.1:8811]
* `--agent-id TEXT`: Distinct reviewer identity  [required]
* `--profile TEXT`: Scientific-review execution profile from the review policy  [required]
* `--campaign-id TEXT`: Restrict claims to one campaign
* `--policy TEXT`: Path to a review policy YAML (default: packaged)
* `--seal FILE`: Dataset seal (default: dataset.seal.json beside --db)
* `--lease-seconds INTEGER RANGE`: Lease duration; a background keepalive renews it during model and query work  [default: 1800; x>=60]
* `--stall-timeout INTEGER RANGE`: Maximum model-CLI silence before treating the connection as stalled  [default: 900; x>=60]
* `--max-members INTEGER RANGE`: Maximum candidate occurrences sampled into one review input  [default: 24; x>=1]
* `--max-input-bytes INTEGER RANGE`: Hard byte bound for the canonical model input  [default: 524288; x>=16384]
* `--verification-threads INTEGER RANGE`: DuckDB threads used by each executable check  [default: 1; x>=1]
* `--max-tasks INTEGER RANGE`: Exit after N completed reviews  [x>=1]
* `--idle-sleep FLOAT RANGE`: Seconds to wait on an empty queue  [default: 5.0; x>=0.1]
* `-v, --verbose`
* `--help`: Show this message and exit.


## `sharur-ops`

```text
usage: sharur-ops [-h] [--host HOST] [--port PORT] [--db DB]
                  [--pool-size POOL_SIZE] [--backup-dir BACKUP_DIR]
                  [--backup-interval BACKUP_INTERVAL]
                  [--allow-insecure-remote]

HTTP control plane for coordinated Sharur agents. One server process owns the
SQLite database. Remote workers communicate only through this API; this
ownership rule is especially important when the database path is on NFS or
another network filesystem.

options:
  -h, --help            show this help message and exit
  --host HOST           Bind host (default: 127.0.0.1)
  --port PORT
  --db DB
  --pool-size POOL_SIZE
  --backup-dir BACKUP_DIR
  --backup-interval BACKUP_INTERVAL
                        Seconds between online SQLite backups; 0 disables
                        scheduled backups
  --allow-insecure-remote
                        Allow remote clients without SHARUR_OPS_TOKEN (unsafe)
```


## `sharur-query`

```text
usage: sharur-query [-h] --db DB [--stage-dir STAGE_DIR] [--direct]
                    [--seal SEAL] [--reserve-gb RESERVE_GB] [--host HOST]
                    [--port PORT] [--ops-url OPS_URL] [--threads THREADS]
                    [--memory-limit MEMORY_LIMIT] [--temp-dir TEMP_DIR]
                    [--max-temp-size MAX_TEMP_SIZE]
                    [--capacity-units CAPACITY_UNITS]
                    [--heavy-weight HEAVY_WEIGHT] [--max-queue MAX_QUEUE]
                    [--light-timeout LIGHT_TIMEOUT]
                    [--heavy-timeout HEAVY_TIMEOUT]

Bounded read-only HTTP data plane for coordinated Sharur agents.

options:
  -h, --help            show this help message and exit
  --db DB
  --stage-dir STAGE_DIR
                        Campaign-local storage for the atomic immutable
                        replica
  --direct              Open the verified immutable source database in place
  --seal SEAL
  --reserve-gb RESERVE_GB
  --host HOST
  --port PORT
  --ops-url OPS_URL     Sharur Ops base URL for per-agent token introspection
  --threads THREADS
  --memory-limit MEMORY_LIMIT
  --temp-dir TEMP_DIR
  --max-temp-size MAX_TEMP_SIZE
  --capacity-units CAPACITY_UNITS
  --heavy-weight HEAVY_WEIGHT
  --max-queue MAX_QUEUE
  --light-timeout LIGHT_TIMEOUT
  --heavy-timeout HEAVY_TIMEOUT
```
