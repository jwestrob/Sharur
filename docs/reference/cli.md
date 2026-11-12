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
* `reconcile-predicates`: Bring predicates current after a...
* `predicates`: List available predicates for search.
* `preflight`: Emit one typed dataset/runtime capability...
* `compact`: Create an exact validated physical copy;...
* `seal`: Write a portable integrity seal for a...
* `migrate`: Apply pending schema/index migrations,...
* `backfill-contig-context`: Record assembly contig lengths and...
* `contig-name-topology`: Show or record which contig-name suffixes...
* `load-bgcs`: Load GECCO biosynthetic gene clusters into...
* `diff-db`: Compare every table of two DuckDB files:...
* `reverify`: Replay findings' verification triples and...
* `verify-seal`: Recompute a dataset seal and report...
* `build-vector-index`: Build mmap-ready FAISS sidecars and a...
* `doctor`: Verify external tools, reference...
* `health`: Check a dataset for problems that distort...
* `adopt-astra-hmms`: Register databases installed through Astra...
* `architecture`: Find proteins whose ordered domains match...
* `import-coverage`: Import per-sample contig coverage into the...
* `abundance`: Genome abundance per sample, or the share...
* `coverage-outliers`: Contigs whose depth departs from their...
* `import-expression`: Import per-sample gene (CDS) read counts...
* `expression`: Genome transcript share per RNA sample, a...
* `browse`: Browse a dataset in a web browser: protein...
* `describe`: What a dataset holds: annotation sources,...
* `card`: Summarize one protein: context,...
* `why`: Explain why a protein carries a predicate:...
* `modules`: KEGG module completeness per genome, or...
* `setup-kegg`: Fetch KEGG data and build the KO ->...
* `import-quality`: Fill genome completeness and contamination...
* `cas-type`: Subtype CRISPR-Cas systems with...
* `stage`: Run a writer against a staged clone of a...
* `sql`: Sharur's SQL surface: per-connection...
* `stats`: Pre-register predictions and check...
* `families`: Annotation-independent protein families...
* `remote-homology`: Remote-homology annotation: CDD per...

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

### `sharur reconcile-predicates`

Bring predicates current after a KEGG-map-only change, without a full regeneration.

Checks that only the KEGG map (or its rule files) changed since the last full generation and that
--old-kegg-map is the map that generation recorded (sha256). Regenerates exactly the proteins
carrying a KO whose predicate set changed, then records a `reconciled` generation, so
`sharur describe` / `preflight` report the predicates as current. Refuses when anything else changed.

Example:
    cp data/reference/kegg/kegg_predicates.tsv /tmp/kegg_predicates.before.tsv
    sharur setup-kegg --rejudge
    sharur reconcile-predicates --db data/DATASET/sharur.duckdb --old-kegg-map /tmp/kegg_predicates.before.tsv

**Usage**:

```console
$ sharur reconcile-predicates [OPTIONS]
```

**Options**:

* `-d, --db TEXT`: Path to DuckDB database  [default: data/sharur.duckdb]
* `--old-kegg-map PATH`: The kegg_predicates.tsv the database's predicates were generated with (copy it aside before `sharur setup-kegg --rejudge`)  [required]
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
* `--semantic-index PATH`: Optional compact semantic index directory to verify against this database.
* `-f, --format [markdown|json]`: Output format: markdown or json  [default: markdown]
* `--skip-tools`: Skip slower external binary/version probes.
* `--strict`: Exit non-zero unless every required dataset capability is available.
* `--help`: Show this message and exit.

### `sharur compact`

Create an exact validated physical copy; reseal after choosing its final location.

**Usage**:

```console
$ sharur compact [OPTIONS]
```

**Options**:

* `-d, --db PATH`: Source DuckDB, opened read-only.  [required]
* `-o, --output PATH`: Required fresh destination; source and sidecars remain in place.  [required]
* `--threads INTEGER RANGE`: CPU workers, bounded by available CPUs.  [default: 2; x>=1]
* `--memory-limit TEXT`: Bounded DuckDB buffer budget; index construction also needs memory.  [default: 8GB]
* `--storage-version TEXT`: Explicit compatibility target, e.g. v1.0.0; default uses writer format.
* `--receipt PATH`: Fresh JSON evidence file containing counts and hashes.
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

Apply pending schema/index migrations, destructive ones included, in a maintenance window.

Migrations add tables, columns and indexes, and a writable open applies those itself. v5 is
destructive: it removes defense/secretion system calls made before replicon-local validation
(raw HMM rows stay) until Stage 07 reruns the caller. Only this command applies it; a writable
open of a database holding such calls stops and names this command. Copy a pre-v5 file first
if its old calls matter.

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
Circularity comes from assembly header marks, and from contig-name suffixes
only under the conventions the database records (`sharur contig-name-topology`).
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

### `sharur contig-name-topology`

Show or record which contig-name suffixes mark circular and linear contigs in a dataset.

Suffix conventions belong to one assembler's or project's naming: srvp (`..._C`/`_RC`
circular, `..._L` linear), metamdbg (`ctg<N>_<cov>x_c|rc|l`), hifiasm (`[sp]tg<N>c|l`).
Stage 07 and `sharur backfill-contig-context` apply only the conventions a database records,
after assembly header marks. Recording leaves stored `is_circular` values as they are; the
report shows, per convention, how many contig IDs it marks and how many stored values differ.

**Usage**:

```console
$ sharur contig-name-topology [OPTIONS]
```

**Options**:

* `-d, --db PATH`: Path to a Sharur DuckDB.  [default: data/sharur.duckdb]
* `--set TEXT`: Record the conventions this dataset's contig IDs follow (comma-separated: srvp, metamdbg, hifiasm), or 'none' to clear. Without --set, only report.
* `--help`: Show this message and exit.

### `sharur load-bgcs`

Load GECCO biosynthetic gene clusters into an existing database.

Replaces every BGC locus (``loci`` rows with locus_type 'bgc', the ``bgc_loci``
view) and its member proteins. Writes the canonical database: run in a
maintenance window, then reseal.

**Usage**:

```console
$ sharur load-bgcs [OPTIONS]
```

**Options**:

* `-d, --db PATH`: Path to a writable Sharur DuckDB.  [default: data/sharur.duckdb]
* `--gecco-dir PATH`: Output root of scripts/run_gecco_bgcs.py (shard dirs with record_map.tsv and *.clusters.tsv).  [required]
* `--help`: Show this message and exit.

### `sharur diff-db`

Compare every table of two DuckDB files: schema, rows, content, keys.

**Usage**:

```console
$ sharur diff-db [OPTIONS] BEFORE AFTER
```

**Arguments**:

* `BEFORE`: Earlier database (e.g. a pre-publish backup).  [required]
* `AFTER`: Later database (e.g. a staged copy).  [required]

**Options**:

* `-t, --table TEXT`: Limit to these tables (repeatable).
* `-k, --key TEXT`: TABLE=COL[,COL] key for added/removed/changed counts (default: the table's primary or unique key).
* `--examples INTEGER`: Show up to N added/removed key values per table.  [default: 0]
* `--threads INTEGER`: [default: 4]
* `--memory-limit TEXT`: [default: 8GB]
* `-f, --format [markdown|json]`: [default: markdown]
* `--help`: Show this message and exit.

### `sharur reverify`

Replay findings' verification triples and report numbers that no longer hold.

Before publishing a refresh, `--compare -s LIVE=STAGED` shows every recorded number the
staged copy would change; after one, `--compare -b LIVE=BACKUP` does the same against a
pre-publish backup. `--reports GLOB` lists the documents that cite each affected finding.

**Usage**:

```console
$ sharur reverify [OPTIONS] FINDINGS...
```

**Arguments**:

* `FINDINGS...`: findings.jsonl files or dataset directories (data/DATASET → exploration/findings.jsonl).  [required]

**Options**:

* `-d, --db PATH`: Database for SQL triples without `attach` (default: data/DATASET/sharur.duckdb).
* `-t, --touching TEXT`: Regex: replay only findings whose queries, or the local scripts they run, match it (e.g. 'defense_systems|defensefinder').
* `--ids TEXT`: Comma-separated finding IDs to replay.
* `-s, --substitute TEXT`: OLD=NEW path substitution, repeatable; redirects queries, attach maps, duckdb.connect and ATTACH in helper scripts.
* `-b, --baseline TEXT`: OLD=NEW substitution for the compare-mode baseline run, repeatable (e.g. LIVE=BACKUP after a publish).
* `--compare`: Run each triple under the baseline (as recorded unless --baseline) and under --substitute; report every result that changes.
* `-j, --workers INTEGER`: Parallel replay processes (default: CPUs / --threads); capped so workers x (--memory-limit + 1.5 GiB) fits half of available memory.
* `--threads INTEGER`: DuckDB threads per replay connection.  [default: 2]
* `--memory-limit TEXT`: DuckDB memory limit per connection.  [default: 4GB]
* `--timeout FLOAT`: Seconds per shell triple.  [default: 1800.0]
* `--reports TEXT`: Glob (relative to the repo root) of documents to scan for citations of failing or changed findings, repeatable; e.g. 'data/srvp_joint/reports/*.md'.
* `-o, --out PATH`: Write every triple's result as TSV.
* `-f, --format [markdown|json]`: [default: markdown]
* `--strict`: Exit 1 on any FAIL/ERROR (or change, with --compare).
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

### `sharur health`

Check a dataset for problems that distort analyses: missing gene calls, unplaced genes,
annotation gaps, completeness, stale labels, seal drift.

**Usage**:

```console
$ sharur health [OPTIONS]
```

**Options**:

* `-d, --db PATH`: Dataset DuckDB (opened read-only).  [default: data/sharur.duckdb]
* `-f, --format [markdown|json]`: markdown or json  [default: markdown]
* `--strict`: Exit non-zero when any check fails.
* `--help`: Show this message and exit.

### `sharur adopt-astra-hmms`

Register databases installed through Astra in Aksha's registry.

Aksha shares Astra's database storage but keeps its own registry, so
Astra-installed databases stay on disk and unknown to `aksha search`.
This marks each one installed at its existing directory (no download);
the previous registry is kept as hmm_databases.json.bak.

**Usage**:

```console
$ sharur adopt-astra-hmms [OPTIONS]
```

**Options**:

* `--dry-run`: List what would be adopted; write nothing.
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

### `sharur import-coverage`

Import per-sample contig coverage into the abundance sidecar (the core database is unchanged).

**Usage**:

```console
$ sharur import-coverage [OPTIONS] PATHS...
```

**Arguments**:

* `PATHS...`: CoverM contig table(s) or long TSV/CSV.  [required]

**Options**:

* `-d, --db PATH`: Core DuckDB (validation and sidecar location).  [default: data/sharur.duckdb]
* `--format TEXT`: coverm (contig table) or long (sample_id, contig_id, ...).  [default: coverm]
* `--samples PATH`: TSV/CSV with sample_id and metadata columns.
* `--sidecar PATH`: Default: abundance.duckdb beside --db.
* `--allow-unknown-contigs`: Keep matching rows when some contigs are absent from the dataset.
* `--help`: Show this message and exit.

### `sharur abundance`

Genome abundance per sample, or the share of reads in genomes carrying a feature.

**Usage**:

```console
$ sharur abundance [OPTIONS]
```

**Options**:

* `-d, --db PATH`: Core DuckDB.  [default: data/sharur.duckdb]
* `--sample TEXT`: Sample(s); repeatable (default: all).
* `-b, --bin TEXT`: Genome(s); repeatable.
* `-p, --predicate TEXT`: Share of reads in genomes carrying this predicate.
* `-a, --annotation TEXT`: Share of reads in genomes carrying this accession or name.
* `-n, --top INTEGER`: Genomes shown per sample (markdown).  [default: 10]
* `--sidecar PATH`: Default: abundance.duckdb beside --db.
* `-f, --format [markdown|json]`: markdown or json  [default: markdown]
* `--help`: Show this message and exit.

### `sharur coverage-outliers`

Contigs whose depth departs from their genome's median in most samples (possible binning errors).

**Usage**:

```console
$ sharur coverage-outliers [OPTIONS] BIN_ID
```

**Arguments**:

* `BIN_ID`: Genome (bin) ID.  [required]

**Options**:

* `-d, --db PATH`: Core DuckDB.  [default: data/sharur.duckdb]
* `--min-log2 FLOAT`: Minimum |log2(contig depth / genome median)|.  [default: 1.0]
* `--sidecar PATH`: Default: abundance.duckdb beside --db.
* `-f, --format [markdown|json]`: markdown or json  [default: markdown]
* `--help`: Show this message and exit.

### `sharur import-expression`

Import per-sample gene (CDS) read counts from metatranscriptomes into the expression sidecar.

**Usage**:

```console
$ sharur import-expression [OPTIONS] PATHS...
```

**Arguments**:

* `PATHS...`: Long table (TSV/CSV/Parquet) or featureCounts table(s).  [required]

**Options**:

* `-d, --db PATH`: Core DuckDB (validation and sidecar location).  [default: data/sharur.duckdb]
* `--format TEXT`: long (sample_id, protein_id, sense|count_s1/count_s2) or featurecounts.  [default: long]
* `--orientation TEXT`: featureCounts tables: sense or antisense counts.  [default: sense]
* `--samples PATH`: TSV/CSV with sample_id and metadata columns (optional sense_orientation = s1 or s2).
* `--sidecar PATH`: Default: expression.duckdb beside --db.
* `--allow-unknown-proteins`: Keep matching rows when some proteins are absent from the dataset.
* `--gene-flags FILE`: CDS overlap table (protein_id, overlap_bp, rna_name) marking genes inside structural RNAs; flagged genes leave expression queries.
* `--min-overlap-bp INTEGER`: Overlap that flags a CDS (--gene-flags).  [default: 60]
* `--help`: Show this message and exit.

### `sharur expression`

Genome transcript share per RNA sample, a genome's top genes, or the expression of a feature.

**Usage**:

```console
$ sharur expression [OPTIONS]
```

**Options**:

* `-d, --db PATH`: Core DuckDB.  [default: data/sharur.duckdb]
* `--sample TEXT`: RNA sample(s); repeatable (default: all).
* `-b, --bin TEXT`: Most expressed genes of this genome.
* `-p, --predicate TEXT`: Expression of genes carrying a predicate.
* `-a, --annotation TEXT`: Expression of genes carrying this accession or name.
* `-n, --top INTEGER`: Rows shown per sample (markdown).  [default: 10]
* `--include-flagged`: Keep genes flagged as overlapping structural RNAs.
* `--sidecar PATH`: Default: expression.duckdb beside --db.
* `-f, --format [markdown|json]`: markdown or json  [default: markdown]
* `--help`: Show this message and exit.

### `sharur browse`

Browse a dataset in a web browser: protein cards, neighborhoods, evidence, genomes, searches.

**Usage**:

```console
$ sharur browse [OPTIONS]
```

**Options**:

* `-d, --db PATH`: Dataset DuckDB (opened read-only).  [default: data/sharur.duckdb]
* `--host TEXT`: Interface to bind; 127.0.0.1 keeps it on this machine.  [default: 127.0.0.1]
* `-p, --port INTEGER`: [default: 8800]
* `--share`: Require a random access token (printed in the link).
* `--notes PATH`: Notes and flags database (default: browser_notes.sqlite beside the dataset).
* `--assemblies PATH`: Directory of genome FASTAs named BIN_ID.fna (for CRISPR repeats and spacers); repeatable. Dataset folders are searched too.
* `--summary-cache / --no-summary-cache`: Compute startup summaries once per database version in a worker process and reuse them (keeps the server near 1 GB).  [default: summary-cache]
* `--semantic-index PATH`: Serve V2 term search and stored term rows from this compact index (its CURRENT generation), verified against the database at startup. Needs pyroaring (sharur[compact]).
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
* `--rejudge`: Re-apply the current rules to the existing build from its evidence snapshot (no download). Covers rule edits that remove pairs or touch KOs the build already maps.
* `--help`: Show this message and exit.

### `sharur import-quality`

Fill genome completeness and contamination from a CheckM, CheckM2 or GTDB table.

**Usage**:

```console
$ sharur import-quality [OPTIONS] TABLE
```

**Arguments**:

* `TABLE`: CheckM2 quality_report.tsv, CheckM --tab_table output, or GTDB ar53/bac120 metadata (.tsv or .tsv.gz).  [required]

**Options**:

* `-d, --db PATH`: Dataset DuckDB.  [default: data/sharur.duckdb]
* `--id-column TEXT`: Genome column (default: detected).
* `--completeness-column TEXT`: Completeness column, percent (default: detected; CheckM2 preferred).
* `--contamination-column TEXT`
* `--overwrite`: Replace completeness values already stored.
* `--dry-run`: Report matches and leave the dataset unchanged.
* `--help`: Show this message and exit.

### `sharur cas-type`

Subtype CRISPR-Cas systems with CRISPRCasTyper's profiles and scoring.

Searches the CCTyper Cas profiles against every genome's proteins, groups
hits into operons, scores subtypes, classifies array repeats and links
operons to arrays within 10 kb. Calls go to ``crispr_cas_systems``,
``crispr_array_types`` and ``system_proteins`` (source ``cctyper``),
replacing earlier calls from this caller.

**Usage**:

```console
$ sharur cas-type [OPTIONS]
```

**Options**:

* `-d, --db PATH`: Dataset DuckDB.  [default: data/sharur.duckdb]
* `--cctyper-db PATH`: CCTyper database directory (default: $SHARUR_CCTYPER_DB, then Aksha's CCTyper install).
* `-w, --workers INTEGER`: Parallel genome searches (default: all CPUs).  [default: 0]
* `-g, --genome TEXT`: Type only these genomes; repeatable.
* `--dry-run`: Report the calls and leave the dataset unchanged.
* `--help`: Show this message and exit.

### `sharur stage`

Run a writer against a staged clone of a dataset database, review, publish.

**Usage**:

```console
$ sharur stage [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `run`: Stage DB, run the command after `--`...
* `publish`: Replace the canonical database with a...
* `discard`: Delete a staged copy and its manifest.

#### `sharur stage run`

Stage DB, run the command after `--` against it, diff, reverify findings, optionally publish.

In the command, `{staged}` is replaced by the staged path, also exported as SHARUR_STAGED_DB:

    sharur stage run --db data/ds/sharur.duckdb --label cazy-rerun --findings data/ds -- \
        python scripts/classify_cazymes.py --db {staged} ...

**Usage**:

```console
$ sharur stage run [OPTIONS]
```

**Options**:

* `-d, --db PATH`: Canonical sharur.duckdb.  [required]
* `-l, --label TEXT`: Short name for this change (backup and staged file names).  [required]
* `--work-dir PATH`: Directory for the staged clone (default: beside the database; must share its filesystem to publish).
* `--findings PATH`: findings.jsonl files or dataset dirs to replay against the staged copy (only findings that read what changed).
* `--reports TEXT`: Globs of documents to scan for citations.
* `--report-dir PATH`: Where to write diff/reverify reports.
* `-j, --workers INTEGER`
* `--publish`: Publish when the writer succeeds.
* `--help`: Show this message and exit.

#### `sharur stage publish`

Replace the canonical database with a reviewed staged copy (backup kept, reseal).

**Usage**:

```console
$ sharur stage publish [OPTIONS] STAGED
```

**Arguments**:

* `STAGED`: The staged file `sharur stage run` printed.  [required]

**Options**:

* `--no-reseal`
* `--help`: Show this message and exit.

#### `sharur stage discard`

Delete a staged copy and its manifest.

**Usage**:

```console
$ sharur stage discard [OPTIONS] STAGED
```

**Arguments**:

* `STAGED`: [required]

**Options**:

* `--help`: Show this message and exit.

### `sharur sql`

Sharur's SQL surface: per-connection macros (neighbors, slot, feature_presence, ...).

**Usage**:

```console
$ sharur sql [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `prelude`: Print the prelude as one SQL script for a...
* `objects`: List the SQL objects installed on a...
* `query`: Run one read-only query with the prelude...

#### `sharur sql prelude`

Print the prelude as one SQL script for a generic DuckDB client.

Run it after `ATTACH 'data/DATASET/sharur.duckdb' AS ds (READ_ONLY); USE ds;`.

**Usage**:

```console
$ sharur sql prelude [OPTIONS]
```

**Options**:

* `-o, --out PATH`: Write the script here instead of stdout.
* `--help`: Show this message and exit.

#### `sharur sql objects`

List the SQL objects installed on a connection to this database, and any that failed.

**Usage**:

```console
$ sharur sql objects [OPTIONS]
```

**Options**:

* `-d, --db PATH`: Dataset DuckDB.  [required]
* `-f, --format [markdown|tsv|json]`: [default: markdown]
* `--help`: Show this message and exit.

#### `sharur sql query`

Run one read-only query with the prelude installed; results never carry sequences.

Example: sharur sql query "FROM slot_tally('pfam:PF00069', 'pfam:PF00072', split := 'K00001')" -d DB

**Usage**:

```console
$ sharur sql query [OPTIONS] SQL
```

**Arguments**:

* `SQL`: SQL to run; '-' reads it from stdin.  [required]

**Options**:

* `-d, --db PATH`: Dataset DuckDB (opened read-only).  [required]
* `-f, --format [markdown|tsv|json]`: [default: markdown]
* `-n, --limit INTEGER`: Rows to print (0 = all).  [default: 200]
* `-o, --out PATH`: Also write every row as TSV here.
* `--threads INTEGER`: DuckDB threads.  [default: 4]
* `--memory-limit TEXT`: DuckDB memory limit.  [default: 8GB]
* `--help`: Show this message and exit.

### `sharur stats`

Pre-register predictions and check statistical claims in findings.

**Usage**:

```console
$ sharur stats [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `preregister`: Freeze predictions, tests and verdict...
* `verify-prereg`: Recompute a record's hash and check its...
* `check`: Check every statistical claim as recorded:...

#### `sharur stats preregister`

Freeze predictions, tests and verdict rules before the analysis runs.

Writes <id>_<UTC time>.prereg.json (read-only, never replaced) and prints its content_sha256,
which every statistical check implementing it cites.

**Usage**:

```console
$ sharur stats preregister [OPTIONS] SPEC
```

**Arguments**:

* `SPEC`: JSON or YAML: id, title, predictions[], correction, design_inputs.  [required]

**Options**:

* `-o, --out-dir PATH`: Directory for the read-only record.  [required]
* `-d, --db PATH`: Dataset DuckDB; its dataset.seal.json fills `dataset` (read only).
* `--document PATH`: A longer pre-registration text to bind by SHA-256.
* `--help`: Show this message and exit.

#### `sharur stats verify-prereg`

Recompute a record's hash and check its document and, with --seal, its dataset.

**Usage**:

```console
$ sharur stats verify-prereg [OPTIONS] RECORD
```

**Arguments**:

* `RECORD`: A .prereg.json record.  [required]

**Options**:

* `--seal`: Also check the dataset against its seal (slow, read-only).
* `-d, --db PATH`: Current database path for --seal.
* `-f, --format TEXT`: markdown or json  [default: markdown]
* `--help`: Show this message and exit.

#### `sharur stats check`

Check every statistical claim as recorded: verdict rules, interval, controls, pre-registration,
dataset seal and input hashes. Replaying the numbers is `sharur reverify`.

**Usage**:

```console
$ sharur stats check [OPTIONS] FINDINGS...
```

**Arguments**:

* `FINDINGS...`: findings.jsonl files or dataset directories.  [required]

**Options**:

* `--root PATH`: Base for relative paths (default: repo root).
* `-f, --format TEXT`: markdown or json  [default: markdown]
* `--strict`: Exit 1 when any statistical check has issues.
* `--help`: Show this message and exit.

### `sharur families`

Annotation-independent protein families (MMseqs2): build, import, activate, assign, status.

**Usage**:

```console
$ sharur families [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `run`: Export proteins, cluster subfamilies and...
* `export`: Write proteins.faa (protein_id order,...
* `build`: Build the staging DuckDB from a finished...
* `load`: Writer: copy the clustering runs of a...
* `activate`: Stage DB, load the staged clustering runs...
* `assign`: Give families to proteins that lack one...
* `import`: Register an externally produced clustering...
* `status`: Runs, family counts, proteins without a...

#### `sharur families run`

Export proteins, cluster subfamilies and families, write a staging DuckDB (the DB is untouched).

**Usage**:

```console
$ sharur families run [OPTIONS]
```

**Options**:

* `-d, --db PATH`: Dataset sharur.duckdb. Opened read-only.  [required]
* `-o, --out-dir PATH`: Work directory (e.g. DATASET/stage09_families).  [required]
* `--subfamily-min-seq-id FLOAT`: [default: 0.5]
* `--family-min-seq-id FLOAT`: [default: 0.3]
* `--coverage FLOAT`: Coverage of both sequences (-c, cov-mode 0).  [default: 0.8]
* `-s, --sensitivity FLOAT`: [default: 6.0]
* `-t, --threads INTEGER`: MMseqs2 threads (default: $SLURM_CPUS_ON_NODE or all CPUs).
* `--duckdb-threads INTEGER`: Default: --threads.
* `--memory-limit TEXT`: DuckDB memory limit (default: 25% of the allocation).
* `--split-memory-limit TEXT`: MMseqs2 memory cap (default: 70% of the allocation).
* `--keep-fasta`
* `--mmseqs TEXT`: MMseqs2 executable.
* `--help`: Show this message and exit.

#### `sharur families export`

Write proteins.faa (protein_id order, trailing '*' removed) and export.json (count, digest).

**Usage**:

```console
$ sharur families export [OPTIONS]
```

**Options**:

* `-d, --db PATH`: Dataset sharur.duckdb. Opened read-only.  [required]
* `-o, --out-dir PATH`: [required]
* `--duckdb-threads INTEGER`: [default: 4]
* `--memory-limit TEXT`: [default: 8GB]
* `--help`: Show this message and exit.

#### `sharur families build`

Build the staging DuckDB from a finished clustering, such as one run on a cluster node.

**Usage**:

```console
$ sharur families build [OPTIONS]
```

**Options**:

* `-d, --db PATH`: Dataset sharur.duckdb. Opened read-only.  [required]
* `-c, --clustering PATH`: Directory with clustering.json, level0/1.tsv.gz and export.json.  [required]
* `-s, --staging PATH`: Staging DuckDB to write (replaced).  [required]
* `--duckdb-threads INTEGER`: [default: 4]
* `--memory-limit TEXT`: [default: 8GB]
* `--help`: Show this message and exit.

#### `sharur families load`

Writer: copy the clustering runs of a staging file into DB (same-name runs are replaced).

**Usage**:

```console
$ sharur families load [OPTIONS]
```

**Options**:

* `-d, --db PATH`: Writable database (a staged clone of a canonical one).  [required]
* `-s, --staging PATH`: [required]
* `-r, --run TEXT`: Clustering run to load (repeatable; default: all).
* `--verify-digest / --no-verify-digest`: [default: verify-digest]
* `--help`: Show this message and exit.

#### `sharur families activate`

Stage DB, load the staged clustering runs into the clone, diff, and (with --publish) swap it in.

Uses the staged-write protocol of `sharur stage` (docs/guides/staged-writes.md): the
canonical file changes only at publish, a backup is kept and the dataset is resealed.

**Usage**:

```console
$ sharur families activate [OPTIONS]
```

**Options**:

* `-d, --db PATH`: Dataset sharur.duckdb.  [required]
* `-s, --staging PATH`: Staging DuckDB from `sharur families run/build`.  [required]
* `-l, --label TEXT`: [default: protein-families]
* `--work-dir PATH`: Directory for the staged clone (default: beside DB).
* `-r, --run TEXT`: Clustering run to load (repeatable; default: all).
* `--publish`: Publish the staged copy when the load succeeds.
* `--verify-digest / --no-verify-digest`: [default: verify-digest]
* `--help`: Show this message and exit.

#### `sharur families assign`

Give families to proteins that lack one (newly added genomes) without re-clustering.

**Usage**:

```console
$ sharur families assign [OPTIONS]
```

**Options**:

* `-d, --db PATH`: Dataset sharur.duckdb.  [required]
* `--work-dir PATH`: MMseqs2 scratch (default: beside DB).
* `-t, --threads INTEGER`
* `--family-run TEXT`: Family run of the Sharur build to extend (default: the default run).
* `--in-place`: Write DB directly (a staged clone or scratch copy); default stages it.
* `-l, --label TEXT`: [default: protein-families-assign]
* `--publish`
* `--mmseqs TEXT`
* `--help`: Show this message and exit.

#### `sharur families import`

Register an externally produced clustering as a run, without re-clustering.

**Usage**:

```console
$ sharur families import [OPTIONS]
```

**Options**:

* `-d, --db PATH`: Dataset sharur.duckdb. Opened read-only.  [required]
* `--tsv PATH`: Two-column representative<TAB>member table (mmseqs createtsv).  [required]
* `-s, --staging PATH`: Staging DuckDB to add the run to (created when absent).  [required]
* `-n, --name TEXT`: Clustering run identifier, e.g. a tool and parameter tag.  [required]
* `--level TEXT`: Level label (family, subfamily or any other label).  [default: imported]
* `-p, --parameters TEXT`: The clustering command or parameters, stored verbatim.  [required]
* `--tool TEXT`: [default: mmseqs2]
* `--tool-version TEXT`
* `--default`: Make this run the default of its level (per-genome views).
* `--duckdb-threads INTEGER`: [default: 4]
* `--memory-limit TEXT`: [default: 8GB]
* `--help`: Show this message and exit.

#### `sharur families status`

Runs, family counts, proteins without a family, and whether a full re-cluster is due.

**Usage**:

```console
$ sharur families status [OPTIONS]
```

**Options**:

* `-d, --db PATH`: Dataset sharur.duckdb. Opened read-only.  [required]
* `--help`: Show this message and exit.

### `sharur remote-homology`

Remote-homology annotation: CDD per protein and HHsearch per protein family ("resembles" rows).

**Usage**:

```console
$ sharur remote-homology [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `check`: Check that the configured tools and...
* `export`: Select families (default: Pfam-dark, >= 3...
* `export-cdd`: Write every unique protein sequence (named...
* `search`: Run the searches on this host (the Slurm...
* `load-families`: Parse the family-tier results into...
* `compact-cdd`: Reduce the raw RPS-BLAST tables to...
* `load-cdd`: Load the CDD hits (per unique sequence)...
* `compact-staging`: Rewrite a staging file without free blocks...
* `activate`: Stage DB, load the remote-homology rows...
* `status`: Searches loaded, families resembling...

#### `sharur remote-homology check`

Check that the configured tools and databases are present; print tool versions.

**Usage**:

```console
$ sharur remote-homology check [OPTIONS]
```

**Options**:

* `-c, --config PATH`: Runner configuration (JSON).  [required]
* `--help`: Show this message and exit.

#### `sharur remote-homology export`

Select families (default: Pfam-dark, >= 3 genomes) and write one FASTA per family.

**Usage**:

```console
$ sharur remote-homology export [OPTIONS]
```

**Options**:

* `-d, --db PATH`: Dataset sharur.duckdb. Opened read-only.  [required]
* `-o, --out-dir PATH`: Work directory for the family tier.  [required]
* `-f, --families PATH`: Staging file with the clustering run (default: the dataset).
* `-r, --run TEXT`: Clustering run (default: the default family run).
* `--min-genomes INTEGER`: Families found in at least this many genomes.  [default: 3]
* `--include-pfam-annotated`: Also select families whose members carry Pfam hits.
* `--max-sequences INTEGER`: Unique sequences per family alignment.  [default: 250]
* `--family-id TEXT`: Export these families only (repeatable).
* `--duckdb-threads INTEGER`: [default: 4]
* `--memory-limit TEXT`: [default: 8GB]
* `--help`: Show this message and exit.

#### `sharur remote-homology export-cdd`

Write every unique protein sequence (named by MD5) in chunks for RPS-BLAST.

**Usage**:

```console
$ sharur remote-homology export-cdd [OPTIONS]
```

**Options**:

* `-d, --db PATH`: Dataset sharur.duckdb. Opened read-only.  [required]
* `-o, --out-dir PATH`: Work directory (cdd/ is created inside).  [required]
* `--chunk-size INTEGER`: Unique sequences per RPS-BLAST chunk.  [default: 2000]
* `--duckdb-threads INTEGER`: [default: 4]
* `--memory-limit TEXT`: [default: 8GB]
* `--help`: Show this message and exit.

#### `sharur remote-homology search`

Run the searches on this host (the Slurm wrapper runs the same code per shard).

**Usage**:

```console
$ sharur remote-homology search [OPTIONS] STAGE
```

**Arguments**:

* `STAGE`: families or cdd  [required]

**Options**:

* `-w, --work PATH`: [required]
* `-c, --config PATH`: [required]
* `-t, --threads INTEGER`: Default: $SLURM_CPUS_ON_NODE or all CPUs.
* `--threads-per-worker INTEGER`: [default: 4]
* `--shard INTEGER`: [default: 0]
* `--n-shards INTEGER`: [default: 1]
* `--limit INTEGER`
* `--help`: Show this message and exit.

#### `sharur remote-homology load-families`

Parse the family-tier results into family_annotations (staging file).

**Usage**:

```console
$ sharur remote-homology load-families [OPTIONS]
```

**Options**:

* `-w, --work PATH`: [required]
* `-c, --config PATH`: [required]
* `-s, --staging PATH`: Staging DuckDB (created when absent).  [required]
* `--allow-partial`: Load although some families are pending or failed.
* `--help`: Show this message and exit.

#### `sharur remote-homology compact-cdd`

Reduce the raw RPS-BLAST tables to work/cdd/compact/*.parquet (run beside them, on the cluster).

**Usage**:

```console
$ sharur remote-homology compact-cdd [OPTIONS]
```

**Options**:

* `-w, --work PATH`: [required]
* `--max-evalue FLOAT`: [default: 0.01]
* `--allow-partial`
* `-t, --threads INTEGER`: Default: $SLURM_CPUS_ON_NODE or 4.
* `--help`: Show this message and exit.

#### `sharur remote-homology load-cdd`

Load the CDD hits (per unique sequence) and models into a staging file; compacts first if needed.

**Usage**:

```console
$ sharur remote-homology load-cdd [OPTIONS]
```

**Options**:

* `-w, --work PATH`: [required]
* `-c, --config PATH`: [required]
* `-s, --staging PATH`: Staging DuckDB (created when absent).  [required]
* `--max-evalue FLOAT`: [default: 0.01]
* `--max-overlap FLOAT`: Best model per region: hits sharing more than this fraction of the shorter range belong to one region (best bit score kept).  [default: 0.5]
* `--allow-partial`
* `--help`: Show this message and exit.

#### `sharur remote-homology compact-staging`

Rewrite a staging file without free blocks or legacy indexes (checked copy, atomic replace).

**Usage**:

```console
$ sharur remote-homology compact-staging [OPTIONS]
```

**Options**:

* `-s, --staging PATH`: Remote-homology staging file to rewrite.  [required]
* `--duckdb-threads INTEGER`: [default: 4]
* `--memory-limit TEXT`: [default: 8GB]
* `--help`: Show this message and exit.

#### `sharur remote-homology activate`

Stage DB, load the remote-homology rows into the clone, diff, and (with --publish) swap it in.

Follows the staged-write protocol (docs/guides/staged-writes.md), as `sharur families
activate` does: the canonical file changes only at publish, a backup is kept and the
dataset is resealed.

**Usage**:

```console
$ sharur remote-homology activate [OPTIONS]
```

**Options**:

* `-d, --db PATH`: Dataset sharur.duckdb.  [required]
* `-s, --staging PATH`: Remote-homology staging file(s) (repeatable).  [required]
* `--families-staging PATH`: Protein-family staging file(s), loaded first in the order given (repeatable; a later file's run replaces a same-named one).
* `--family-runs TEXT`: Clustering runs taken from the families staging files: searched (those the remote-homology rows refer to) or all (every run, e.g. a dataset's own build and imported clusterings).  [default: searched]
* `-l, --label TEXT`: [default: remote-homology]
* `--work-dir PATH`: Directory for the staged clone (default: beside DB).
* `--publish`: Publish the staged copy when the load succeeds.
* `--verify-digest / --no-verify-digest`: [default: verify-digest]
* `--cdd-rows TEXT`: CDD rows written to annotations: best_per_region (the best model of each region) or all. The full hit table stays in the staging file.  [default: best_per_region]
* `--help`: Show this message and exit.

#### `sharur remote-homology status`

Searches loaded, families resembling characterized families, CDD coverage.

**Usage**:

```console
$ sharur remote-homology status [OPTIONS]
```

**Options**:

* `-d, --db PATH`: Dataset sharur.duckdb. Opened read-only.  [required]
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
* `--skip-aksha, --skip-astra / --no-skip-aksha, --no-skip-astra`: Skip standard Stage 04 (Aksha HMM search via 04_astra_scan.py)  [default: no-skip-aksha]
* `--skip-gecco / --with-gecco`: Skip optional Stage 05a by default; use --with-gecco to enable it  [default: skip-gecco]
* `--skip-dbcan / --with-legacy-dbcan`: Skip deprecated Stage 05b by default; use --with-legacy-dbcan to enable it  [default: skip-dbcan]
* `--enable-cazymes`: Run the optional Stage 07 dbCAN three-tool consensus classifier
* `--pipeline-depth INTEGER RANGE`: Maximum ordered V2 transform chunks in flight during Stage 07  [default: 2; x>=1]
* `--contig-name-topology TEXT`: Contig-name suffix conventions that mark circular/linear contigs in this dataset (comma-separated: srvp, metamdbg, hifiasm; or none); Stage 07 records them in the database. Default: what the database being replaced records, else none
* `--skip-crispr / --no-skip-crispr`: Skip standard Stage 05c (CRISPR arrays via minced_crispr.py)  [default: no-skip-crispr]
* `--skip-embeddings / --no-skip-embeddings`: Skip Stage 06 post-build embeddings (06_esm2_embeddings.py)  [default: no-skip-embeddings]
* `--reads PATH`: TSV of sample_id, read1[, read2 | interleaved]: adds optional Stage 08 per-sample coverage with CoverM (abundance.duckdb sidecar)
* `--with-families / --skip-families`: Run optional Stage 09: MMseqs2 protein families and subfamilies into stage09_families/protein_families.duckdb (activate with `sharur families activate`)  [default: skip-families]
* `--family-min-seq-id FLOAT`: Stage 09 family identity threshold  [default: 0.3]
* `--with-cdd / --skip-cdd`: Run optional Stage 10: RPS-BLAST of every unique protein against CDD into stage10_cdd/remote_homology_cdd.duckdb (needs --remote-homology-config)  [default: skip-cdd]
* `--with-family-homology / --skip-family-homology`: Run optional Stage 11: HHsearch of Pfam-dark families (>= 3 genomes) against Pfam and PDB into stage11_family_homology/ (needs --with-families and --remote-homology-config)  [default: skip-family-homology]
* `--remote-homology-config PATH`: Stages 10/11 configuration: tools and databases (JSON)
* `--family-homology-run TEXT`: Stage 11 clustering run (default: the default family run)
* `--subfamily-min-seq-id FLOAT`: Stage 09 subfamily identity threshold  [default: 0.5]
* `--embedding-model TEXT`: Hugging Face protein encoder for Stage 06 (default: facebook/esm2_t6_8M_UR50D)
* `--profile TEXT`: Execution profile: auto, local (CPU), mps, or slurm  [default: auto]
* `--slurm-cpus INTEGER RANGE`: With --profile slurm: CPUs to request for each stage job (--cpus-per-task); stages use every allocated CPU. Default: full-node jobs with no CPU request  [x>=1]
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
                    [--seal SEAL] [--semantic-index SEMANTIC_INDEX]
                    [--reserve-gb RESERVE_GB] [--host HOST] [--port PORT]
                    [--ops-url OPS_URL] [--threads THREADS]
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
  --semantic-index SEMANTIC_INDEX
                        Compact V2 term index directory for /v1/atoms/proteins
                        (needs pyroaring).
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
