# Sharur Quickstart: New Dataset Ingestion

**Goal:** Go from genome assembly FASTAs to an exploration-ready Sharur dataset.

## Canonical Path

The standard ingest workflow starts from nucleotide assemblies (`.fna`, `.fa`, `.fasta`) and should be run through `sharur-ingest`.

- Use `sharur-ingest` as the default interface for new dataset ingestion.
- Use the staged scripts in `src/ingest/` only when you need manual stage control, debugging, or rerunning one stage.

- Do not call `astra search` directly. Use `src/ingest/04_astra_scan.py`.
- Do not manually insert annotation rows into DuckDB. Use `src/ingest/07_build_knowledge_base.py`.
- Do not skip `src/ingest/minced_crispr.py` in the standard pipeline.

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

## Shared Query Service for Multi-Agent Campaigns

For a large database shared by coordinated agents, launch one bounded read-only
data plane after sealing:

```bash
pip install -e ".[ops]"

export SHARUR_OPS_URL=http://ops-host:8811
export SHARUR_QUERY_STAGE_DIR=/local-nvme/sharur/my_dataset

sharur-query \
  --db data/my_dataset/sharur.duckdb \
  --host 0.0.0.0 \
  --threads 16 \
  --memory-limit 32GB \
  --max-temp-size 256GB
```

The service verifies the dataset seal, stages one atomic immutable local
replica, owns one DuckDB instance/cache, and authenticates agent tokens through
Sharur Ops. Typed endpoints enforce queue, execution, row, request, and result
bounds. See [`docs/query_service.md`](docs/query_service.md) for deployment,
resource arithmetic, cancellation, and telemetry.

For exhaustive genome-by-genome reading, build sealed one-genome ownership
units and enqueue them through Ops:

```bash
sharur migrate --db data/my_dataset/sharur.duckdb
sharur seal --db data/my_dataset/sharur.duckdb --force

sharur-atlas plan \
  --db data/my_dataset/sharur.duckdb \
  --output-dir data/my_dataset/atlas

sharur-atlas enqueue \
  --plan-dir data/my_dataset/atlas \
  --ops-url http://ops-host:8811 \
  --query-url http://query-host:8812
```

Atlas traverses every contig through bounded sequence-free packets and proves
completion with per-genome coverage manifests. See
`.claude/skills/atlas.md`.

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
    similar = b.find_similar(giants.records[0]["protein_id"], k=20)
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
- Dataset layout and archival conventions: `DATA_ORGANIZATION.md`
