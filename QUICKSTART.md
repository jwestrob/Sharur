# Sharur Quickstart: New Dataset Ingestion

**Goal:** Go from genome assembly FASTAs to an exploration-ready Sharur dataset.

## Canonical Path

The standard ingest workflow starts from nucleotide assemblies (`.fna`, `.fa`, `.fasta`) and should be run through `sharur-ingest`.

- Use `sharur-ingest` as the default interface for new dataset ingestion.
- Use the staged scripts in `src/ingest/` only when you need manual stage control, debugging, or rerunning one stage.

- Run annotation through `src/ingest/04_astra_scan.py`, which supplies Sharur's per-database Astra settings.
- Load annotation rows through `src/ingest/07_build_knowledge_base.py`, the single DuckDB writer.
- Keep `src/ingest/minced_crispr.py` in the standard pipeline; CRISPR array loci come from it.

If you only have pre-called proteins and no assemblies, see [Alternative: Protein-Only Ingest](#alternative-protein-only-ingest). That path is a special case, not the default workflow.

## Prerequisites

### Software
- Python 3.10+ with Sharur installed: `pip install -e ".[embeddings]"`
- Prodigal
- [Astra](https://github.com/Dreycey/Astra)
- MinCED
- Optional but recommended GPU access for Stage 06 embeddings

If `sharur-ingest` is not available after install, refresh the editable install:
`pip install -e ".[embeddings]"`

### Astra databases
- Standard: `PFAM`, `KOFAM`, `HydDB`, `DefenseFinder`, `dbCAN`
- Optional: `TXSScan`, `VOGdb`, `CANT-HYD`

### Reference data built or fetched once per machine
- `sharur setup-kegg` builds the KO → predicate map and KEGG module definitions locally from KEGG REST
  (academic use; `--inputs` builds from a licensed KEGG copy). KEGG-derived predicates and
  `sharur modules` use it.
- VOGdb's `vog.annotations.tsv` supplies VOG descriptions and functional categories for predicates
  and the browser. See [INSTALL.md](INSTALL.md#reference-databases).
- Optional: [CoverM](https://github.com/wwood/CoverM) for per-sample read coverage (Stage 08).

### Input data
- Genome assembly FASTAs in one directory
- Extensions supported by Stage 00: `.fna`, `.fa`, `.fasta`

## Standard Pipeline

```bash
sharur-ingest \
  --input-dir /path/to/genome_fastas \
  --data-dir data/my_dataset \
  --output data/my_dataset/sharur.duckdb \
  --pipeline-depth 2 \
  --profile auto
```

This is the primary interface for running the standard pipeline. It orchestrates the staged workflow, including the standard stages 00, 03, 04, 05c, 07, and 06 plus the internal `06i` persistent-index attempt, while exposing skip flags for optional stages when needed.

Two options shape the optional layers:

- `--embedding-model MODEL` selects the Hugging Face protein encoder for Stage 06 (ESM-2 8M by
  default). Any model loadable with `AutoModel` works; embeddings computed elsewhere load with
  `sharur build-vector-index --embeddings FILE`.
- `--reads reads.tsv` adds Stage 08, per-sample read coverage with CoverM. The table is
  tab-separated with `sample_id`, `read1`, and `read2` (paired) or `interleaved`; further columns
  become sample metadata. Coverage lands in an `abundance.duckdb` sidecar beside the dataset.

Ingest is a dependency-aware DAG, not an unconditional script list. It records runs and
stage attempts in `data/my_dataset/sharur_ops.db`. Resume is on by default: a stage is reused
only when its command, inputs, script, resource request, and dependency signatures match a
successful ledger entry and all declared outputs still match that attempt's recorded
snapshots. Use `--no-resume` or `--force` for an intentional rerun.
If any dependency executes instead of being reused, its downstream stages execute as well;
this prevents a stale derived artifact from surviving an upstream repair.

Execution profiles are explicit:

- `--profile auto`: MPS when a usable PyTorch MPS backend is detected, local CPU otherwise
- `--profile local`: bounded local CPU workers
- `--profile mps`: local CPU stages plus one exclusively locked MPS Stage 06 process
- `--profile slurm`: write a dependency-linked bundle under `data/my_dataset/slurm/`;
  add `--submit-slurm` only when ready to submit it

The default plan skips optional QUAST, DFAST, GECCO, and the deprecated legacy dbCAN
helper. Enable them deliberately with `--with-quast`, `--with-dfast`, `--with-gecco`,
or `--with-legacy-dbcan`. The distinct Stage 07 dbCAN three-tool consensus classifier is
opt-in with `--enable-cazymes`. Use `--dry-run` to inspect stage order and paths
without creating the dataset directory.

`--pipeline-depth 2` is the bounded Stage 07 default. It keeps two ordered V2
transform chunks in flight while the sole DuckDB writer commits the preceding
chunk.

## Manual Stage-by-Stage Pipeline

Use this only if you need direct control over individual stages or you are debugging a specific stage failure.

```bash
DATASET=my_dataset
INPUT=/path/to/genome/fastas

mkdir -p data/$DATASET/{source,annotations,embeddings,structures,exploration,figures,reports,survey}

# Stage 00: validate and organize assemblies
python src/ingest/00_prepare_inputs.py \
  -i $INPUT \
  -o data/$DATASET/stage00_prepared \
  --copy

# Stage 03: gene calling
python src/ingest/03_prodigal.py \
  -i data/$DATASET/stage00_prepared \
  -o data/$DATASET/stage03_prodigal \
  --max-workers 8

# Stage 04: annotation
python src/ingest/04_astra_scan.py \
  -i data/$DATASET/stage03_prodigal \
  -o data/$DATASET/stage04_astra \
  -t 12

# Stage 05c: CRISPR arrays
python src/ingest/minced_crispr.py \
  -i data/$DATASET/stage00_prepared \
  -o data/$DATASET/stage05c_crispr

# Stage 07: build DuckDB knowledge base + V2 predicates
python src/ingest/07_build_knowledge_base.py \
  -d data/$DATASET \
  -o data/$DATASET/sharur.duckdb \
  --pipeline-depth 2

# Stage 06: embeddings (run after Stage 07)
python src/ingest/06_esm2_embeddings.py \
  data/$DATASET/stage03_prodigal \
  data/$DATASET/embeddings/ \
  --device mps
```

This is the manual reference sequence behind `sharur-ingest`. `04_astra_scan.py` supplies the correct per-database flags, and `07_build_knowledge_base.py` is the loader that consolidates proteins, annotations, loci, validated systems, and V2 predicates.

Stage 07 defaults to a bounded depth-two V2 pipeline that overlaps source
reading, process-pool transformation, and ordered single-writer DuckDB commits.
Its checkpoint always describes a fully committed source prefix, so
`--resume-v2` retains exact restart semantics.

If Stage 07 reaches V2 generation and the semantic implementation changes,
reuse the completed upstream tables with:

```bash
python src/ingest/07_build_knowledge_base.py \
  -d data/$DATASET \
  -o data/$DATASET/sharur.duckdb \
  --restart-v2
```

Use `--resume-v2` for an interruption under the same semantic code/config.

## What Each Stage Does

- `00_prepare_inputs.py`: audits the complete assembly set before publishing it. It rejects
  malformed/empty records, invalid nucleotide symbols, duplicate record IDs (within or
  across files), and colliding normalized genome IDs; accepted inputs receive SHA-256
  checksums plus a processing manifest.
- `03_prodigal.py`: produces per-genome protein FASTAs and the `all_protein_symlinks/` directory used by Stage 04
- `04_astra_scan.py`: runs PFAM, KOFAM, HydDB, DefenseFinder, and dbCAN with Sharur's expected settings
- `minced_crispr.py`: finds CRISPR repeat-spacer arrays from nucleotide assemblies
- `07_build_knowledge_base.py`: loads stage outputs into `sharur.duckdb`, integrates supported validation steps, writes `semantic_atoms`/`semantic_state`, and materializes legacy-compatible predicates. Its slower dbCAN three-tool consensus path runs only with `--enable-cazymes`.
- `06_esm2_embeddings.py`: streams FASTA records through ESM2 into an atomically published
  canonical H5 without retaining the proteome in memory. Direct use builds the sidecars in a
  Torch-free child process; `sharur-ingest` records that CPU build separately as `06i`.
- `vector_index_runner.py`: produces the generation-scoped persistent FAISS sidecar, the
  disk-backed stable protein-ID map, and the atomic index manifest
- `08_coverage.py` (with `--reads`): maps each sample to the prepared assemblies with CoverM using
  every allocated CPU and imports depth, covered fraction, read counts and contig lengths into the
  abundance sidecar

## Verify the Dataset

```bash
sharur preflight --db data/my_dataset/sharur.duckdb
# Machine-readable:
sharur preflight --db data/my_dataset/sharur.duckdb --format json

# Record the completed canonical dataset state:
sharur seal --db data/my_dataset/sharur.duckdb
# Later, or after copying/archive restoration:
sharur verify-seal data/my_dataset/dataset.seal.json
```

The typed brief inspects the live dataset without mutating it. It reports
`available`, `unavailable`, `stale`, or `failed` for core tables/schema,
`annotations`-table sources, whatever structured caller resources actually exist, V2 and
compatibility coverage,
embeddings, persistent similarity index, dataset run ledger, execution profiles, and the
external toolchain. Add `--strict` for a non-zero exit when a required dataset capability is
not available; add `--skip-tools` when binary/version probes are not needed.

Inspect a structured caller result or protein without flattening raw domains
and caller-emitted names:

```bash
sharur inspect ENTITY_ID \
  --type system \
  --upstream 4 \
  --downstream 8 \
  --db data/my_dataset/sharur.duckdb
```

Run a controlled ORF-context comparison with `sharur compare-context`.
Assembly/host evidence can optionally be imported into a separate
`assembly_evidence.duckdb` sidecar. No assembly composition work runs
automatically; `sharur compute-composition-evidence` is the explicit opt-in
command. See [`docs/cases_and_evidence.md`](docs/cases_and_evidence.md).

Stage 00 is itself an integrity gate. It writes
`stage00_prepared/input_integrity.json` for both accepted and rejected input sets. A rejected
set exits non-zero and does not expose assembly links or a downstream
`processing_manifest.json`.

`sharur seal` writes `dataset.seal.json` atomically and refuses to overwrite it without
`--force`. The default structural seal fully hashes small canonical files and reads bounded,
deterministic samples from large DuckDB/H5/FAISS artifacts; it also records Stage-00 source
SHA-256 values, the live DuckDB schema and table counts, annotation sources, whatever
structured caller resources actually exist, canonical findings, and completed ingest
signatures. Use `--full` for an archival content seal that streams every discovered
canonical artifact through SHA-256. Tool/reference versions are optional provenance via
`--include-tools`; volatile operational state does not define the scientific dataset ID.
The command refuses to seal while any ingest run is active. `sharur verify-seal`
exits non-zero on canonical drift and supports `--format json`.

Stage 07 creates the final indexes and runs DuckDB `ANALYZE` before the dataset
is sealed, giving the optimizer statistics over the final table state.

## Look Around

```bash
DB=data/my_dataset/sharur.duckdb
sharur describe --db $DB                                  # sources, curated callers, predicate map state
sharur card PROTEIN_ID --db $DB                           # one protein: context, hits, evidence-backed predicates
sharur why PROTEIN_ID sam_binding --db $DB                # the evidence behind one predicate
sharur modules --db $DB --bin GENOME_ID --min-completeness 0.75
sharur architecture "TPR_* {10,}" --db $DB               # ordered-domain pattern search
sharur browse --db $DB                                    # read-only website at http://localhost:8800/
```

`sharur browse` serves the dataset as linked pages: the tree of life, functions, pathways, Pfam and
VOG families, systems, a contig viewer and protein pages with sequences. `sharur browse --share
--host 0.0.0.0` prints a link with an access token for collaborators on your network. See
[`docs/guides/browser.md`](docs/guides/browser.md).

### Coverage from existing mappings

Import CoverM `contig` tables (or a long table of `sample_id`, `contig_id`, `mean_depth`, …) and
ask about abundance:

```bash
sharur import-coverage s1.tsv s2.tsv --db $DB --samples samples.tsv
sharur abundance --db $DB --predicate nife_group1
sharur coverage-outliers GENOME_ID --db $DB
```

### Datasets built before schema 8

New ingests record each contig's assembly length and Prodigal's truncation flags. For an existing
dataset, add them once, then reseal:

```bash
sharur backfill-contig-context --db $DB --assemblies path/to/assemblies
sharur seal --db $DB --force
```

## Shared Query Service for Multi-Agent Campaigns

For a large database shared by coordinated agents, launch one bounded read-only
data plane after sealing:

```bash
pip install -e ".[ops]"

export SHARUR_OPS_URL=http://ops-host:8811

sharur-query \
  --db data/my_dataset/sharur.duckdb \
  --direct \
  --host 0.0.0.0 \
  --threads 16 \
  --memory-limit 32GB \
  --max-temp-size 256GB
```

The service verifies the dataset seal, owns one read-only DuckDB instance/cache,
and authenticates agent tokens through Sharur Ops. Direct mode avoids a
same-tier copy; `--stage-dir` remains available for a genuinely distinct
storage tier. Typed endpoints enforce queue, execution, row, request, and
result bounds. See [`docs/query_service.md`](docs/query_service.md) for
deployment, resource arithmetic, cancellation, and telemetry.

For exhaustive genome-by-genome reading, build sealed one-genome ownership
units and enqueue them through Ops:

```bash
sharur migrate --db data/my_dataset/sharur.duckdb
sharur seal --db data/my_dataset/sharur.duckdb --force

sharur-atlas plan \
  --db data/my_dataset/sharur.duckdb \
  --output-dir data/my_dataset/atlas \
  --packet-bytes 524288

sharur-atlas packet-census \
  --plan-dir data/my_dataset/atlas

sharur-atlas verify-packet-census \
  --plan-dir data/my_dataset/atlas \
  --deep

sharur-atlas enqueue \
  --plan-dir data/my_dataset/atlas \
  --ops-url http://ops-host:8811 \
  --query-url http://query-host:8812
```

Atlas packs consecutive records from exactly one bin per bounded,
sequence-free model packet. Whole contigs stay together when they fit;
oversized contigs resume by stable protein offset. The required byte target is
the model-facing budget; count ceilings are derived exactly from that byte
target and canonical serialized record floors. Real frames sampled
deterministically across genome-size quantiles provide the preliminary call
projection. The zero-model-call census must pass before enqueue, and
per-genome coverage manifests prove every frame, contig segment, and protein
total. Each scanner also emits typed candidate occurrences and one reconciled
unit disposition. Reduce and route the resulting review DAG:

```bash
sharur-review reduce --ops-url http://ops-host:8811 \
  --campaign-id CAMPAIGN_ID
sharur-review route --ops-url http://ops-host:8811 \
  --campaign-id CAMPAIGN_ID --watch
sharur-review status --ops-url http://ops-host:8811 \
  --campaign-id CAMPAIGN_ID
```

See `.claude/skills/atlas.md` and `docs/review_workflow.md`.

For a legacy H5 without sidecars:

```bash
sharur build-vector-index --db data/my_dataset/sharur.duckdb
```

Ordinary session startup discovers these artifacts but does not open H5 or FAISS. The first
similarity call opens the committed FAISS generation read-only with mmap and uses the SQLite
row-to-protein map; if sidecars are missing or stale, that call can build them.

## Start Exploring

### Claude Code

```bash
/survey
/explore --focus metabolism
```

### Python API

```python
from sharur.operators import Sharur

b = Sharur("data/my_dataset/sharur.duckdb", read_only=True)

giants = b.search_by_predicates(has=["giant", "unannotated"])
defense = b.search_by_predicates(has=["crispr_associated"])
if giants.records:
    pid = giants.records[0]["protein_id"]
    b.card(pid)                                     # bounded summary
    b.why(pid, "unannotated")                       # evidence path for one predicate
    similar = b.find_similar(pid, k=20)
b.modules(min_completeness=0.75)                    # KEGG module completeness (after setup-kegg)
b.search_architecture("Big_* {5,} . VWA")          # domain-architecture patterns
```

## Alternative: Protein-Only Ingest

If assemblies are unavailable and you only have pre-called proteins, you can bootstrap a database with:

```bash
python scripts/ingest_protein_fasta.py \
  /path/to/proteins.faa.gz \
  --output data/$DATASET/sharur.duckdb
```

Use this only when the standard assembly-based pipeline is impossible. It:
- loads bins, contigs, and proteins
- does not run Prodigal, Astra, or MinCED
- does not replace Stage 04 or Stage 07
- leaves annotation loading and downstream validation to you

## Common Mistakes

- Reconstructing the standard pipeline manually when `sharur-ingest` would do the job
- Running raw `astra search` commands instead of `src/ingest/04_astra_scan.py`
- Passing `--databases` as a space-separated list instead of repeated `-d` flags when overriding defaults
- Skipping `minced_crispr.py` and then expecting CRISPR array loci in DuckDB
- Treating legacy Stage `05b` dbCAN as the standard CAZyme path; standard ingest uses Stage 04 + Stage 07

## More Detail

- Full manual stage reference: `src/ingest/README.md`
- Tool-specific details: `docs/tools_reference.md`
- Dataset layout and archival conventions: `docs/DATA_ORGANIZATION.md`
