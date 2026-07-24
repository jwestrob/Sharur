# Sharur Query Service

## Purpose

`sharur-query` is the shared analytical data plane for coordinated agents.
One process owns one read-only DuckDB instance, its buffer cache, worker pool,
spill directory, and memory ceiling. Agents call bounded Sharur operators over
HTTP using their existing Sharur Ops credentials.

```text
agents
  | coordination records          | typed analytical requests
  v                               v
sharur-ops                       sharur-query
  | one SQLite owner               | one DuckDB owner and cache
  v                                v
sharur_ops.db                    verified immutable database
                                (direct source or staged replica)
                                   |
                                   v
                              Sharur operators

executor / Slurm / local launcher: starts compute selected by Ops tasks
```

This keeps three concerns independent:

- Sharur Ops owns campaigns, tasks, leases, findings, and agent identity.
- Sharur Query owns shared read-only analytical access.
- Executors map logical tasks to local processes, Slurm jobs, arrays, or other
  compute runtimes.

Pipeline algorithms stay independent of scheduler packing.

## Production launch

Install the HTTP dependencies, complete all dataset writes, and seal the
canonical dataset:

```bash
pip install -e ".[ops]"
sharur migrate --db data/DATASET/sharur.duckdb
sharur seal --db data/DATASET/sharur.duckdb
sharur verify-seal data/DATASET/dataset.seal.json
```

Run migration with query services stopped. Schema migration 6 adds the
`(bin_id, contig_id)` access path required by contig pagination; migration
changes canonical database state, so existing datasets require a fresh seal
and query-service restart afterward. Staged deployments also refresh their
replica.

Launch Sharur Ops and register agent identities as described in
[`agent_ops_spec.md`](../agent_ops_spec.md). Then launch the analytical service:

```bash
export SHARUR_OPS_URL=http://ops-host:8811
sharur-query \
  --db /shared/data/DATASET/sharur.duckdb \
  --direct \
  --host 0.0.0.0 \
  --threads 16 \
  --memory-limit 32GB \
  --temp-dir /shared/data/DATASET/query_spill \
  --max-temp-size 256GB
```

Each agent uses its per-agent Ops bearer token:

```python
from sharur.query import SharurQuery

query = SharurQuery(
    "http://query-host:8812",
    api_token=worker_token,
)

overview = query.overview()
hits = query.search_predicates(
    has=["giant", "unannotated"],
    limit=100,
)
context = query.neighborhood(
    hits["raw"][0]["protein_id"],
    window=10,
    all_annotations=True,
)
contigs = query.list_contigs("genome-id", limit=25)
packet = query.genome_packet(
    "genome-id",
    max_contigs=128,
    max_proteins=500,
    max_model_payload_bytes=524288,
)
```

The same token identifies the caller in both services. Sharur Query asks
`sharur-ops /auth/whoami` to resolve agent ID and role, then caches the
token-hash-to-principal result for 30 seconds. The raw token stays out of
service state and logs. A standalone deployment can set
`SHARUR_QUERY_TOKEN`; loopback-only development can run with neither setting.
Remote binding requires Ops introspection or a query token.

## Immutable database selection

`--direct` opens the verified canonical database in place. This is the
Biotite profile because the canonical database and any replica would occupy
the same VAST tier; copying a corpus-scale database on the same filesystem
adds capacity and I/O cost while preserving the same storage performance.
One Query process remains the single read-only owner and shared cache.

Use `--stage-dir` or `SHARUR_QUERY_STAGE_DIR` when the replica target is a
genuinely distinct local or faster storage tier. Startup then performs this
sequence:

1. Recompute and verify the source `dataset.seal.json`.
2. Resolve the DuckDB canonical-artifact record from the seal.
3. Acquire a cross-process staging lock.
4. Reuse a matching dataset-ID-qualified replica when present.
5. Check free space for the database size plus a 10 GiB reserve.
6. Copy to a unique partial file while checking source size and modification
   time stability.
7. Verify the copied artifact digest, set mode `0444`, and publish it with an
   atomic rename.
8. Write a small `.stage.json` provenance sidecar.

Agent count has zero effect on replica count. The default structural seal uses
deterministic sampling for large artifacts; `sharur seal --full` supplies
full-content SHA-256 assurance when the additional streaming pass is
appropriate.

## DuckDB ownership and resource budgets

Sharur Query opens the selected immutable database with `read_only=True`, creates one root
connection, and gives each execution thread its own cursor from that root.
All cursors share the same DuckDB database instance and buffer cache. A
process-level advisory lock enforces one query-service owner per database path,
and the launcher fixes Uvicorn at one worker.

Default budgets:

| Setting | Default |
|---|---:|
| DuckDB threads | 8 |
| DuckDB memory limit | 16 GB |
| Spill limit | 128 GB |
| Admission capacity | 12 units |
| Light-query weight | 1 unit |
| Heavy-query weight | 3 units |
| Waiting queue | 128 requests |
| Queue wait limit | 30 seconds |
| Light execution limit | 30 seconds |
| Heavy execution limit | 300 seconds |
| Request body | 1 MiB |
| Serialized result | 2 MiB |
| Maximum rows per request | 500 |

The thread and memory settings belong to the shared DuckDB instance and remain
fixed as agent count grows. Weighted admission permits up to twelve
simultaneous light requests, four simultaneous heavy requests, or a FIFO mix
whose weights sum to twelve. Tune capacity, weights, memory, and spill
together for the service host.

The runtime disables external access and automatic extension installation,
then locks DuckDB configuration before serving requests. The HTTP surface
contains typed operators only; arbitrary SQL stays on the trusted local Python
surface.

## Typed HTTP surface

| Method | Path | Class |
|---|---|---|
| `GET` | `/v1/overview` | heavy |
| `POST` | `/v1/genomes/search` | light |
| `GET` | `/v1/genomes/{genome_id}` | light |
| `GET` | `/v1/genomes/{genome_id}/contigs` | light |
| `GET` | `/v1/genomes/{genome_id}/contigs/{contig_id}` | light |
| `POST` | `/v1/contigs/packet` | light |
| `POST` | `/v1/genomes/packet` | light |
| `POST` | `/v1/proteins/search` | light |
| `GET` | `/v1/proteins/{protein_id}` | light |
| `POST` | `/v1/neighborhood` | light |
| `POST` | `/v1/predicates/search` | heavy |
| `POST` | `/v1/atoms/search` | heavy |
| `POST` | `/v1/atoms/proteins` | heavy |
| `POST` | `/v1/queries/{query_id}/cancel` | control |

Requests reject unknown fields and enforce bounds on strings, lists, offsets,
windows, and row counts. Protein detail accepts verbosity 0 or 1, keeping raw
sequences outside agent-visible responses.

Contig lists and diagnostic contig packets use opaque keyset cursors. A
diagnostic packet is scoped to one contig and contains at most 250 proteins.

Atlas uses `/v1/genomes/packet`. Its cursor is scoped to the sealed dataset
ID, exact `bins.bin_id`, and exact packing-contract hash. Each packet carries
consecutive contig segments from that one bin in `contig_id`, then
`gene_index NULLS LAST, start, protein_id` order. Whole contigs stay together
when they fit a fresh packet; a contig continues by stable protein offset
when its remaining payload exceeds the contract. Defaults are 128 contigs,
500 proteins, and 512 KiB of canonical model payload. Serialized-byte,
protein, and contig limits all apply.

Each model payload nests proteins under contigs bearing the requested bin ID.
The operator also asserts the bin on every selected protein before
serialization. Both packet forms separate raw `observed_annotations` from
exact structured `named_calls` and `loci`; sequence text remains compute-only.

`sharur-atlas packet-census` traverses the real packet stream with zero model
calls, validates one-bin isolation and complete offsets, and records exact
invocation and payload-byte counts. It also checks the compact operator result
plus a service-envelope reserve against the plan's Query result limit. Atlas
enqueue requires a matching completed census.

Every request can supply `X-Sharur-Query-ID`. The Python client generates one
automatically. A caller can use that identifier with the cancellation endpoint;
workers can cancel their own requests and operators can cancel any request.
Server timeouts and client disconnects invoke DuckDB cursor interruption. The
runtime discards an interrupted cursor before its thread serves another query.

## Observability

Authenticated endpoints expose:

- `/health`: dataset identity, seal strength, resource budget, hardened DuckDB
  settings, and current admission state.
- `/metrics`: low-cardinality Prometheus counters for active/queued requests,
  admissions, rejections, queue waits, per-operator requests, errors,
  cancellations, timeouts, and duration totals.
- `/auth/whoami`: the principal resolved for the supplied credential.

Responses carry `X-Sharur-Query-ID` plus a `service` object containing the
sealed dataset ID, operator class, queue wait, and execution duration.
`SharurResult.trace.dataset_version` uses the same sealed ID, so staging paths
and file timestamps cannot change scientific trace identity.

## Access-path selection

Use `Sharur(path, read_only=True)` for one local analyst, scripts that require
trusted SQL, and modest concurrent work where each process has an explicit
thread/memory/spill budget.

Use `sharur-query` for coordinated campaigns, large shared databases, and
agent concurrency. Findings and hypotheses continue through Sharur Ops;
canonical DuckDB writes occur in ingestion or another serialized maintenance
window, followed by a fresh seal and query-service restart. Replica-backed
deployments regenerate the replica at that point.
