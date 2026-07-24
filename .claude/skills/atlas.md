# Atlas Skill

Exhaustive bottom-up reading of every genome, contig, and protein in a sealed
Sharur dataset. Atlas uses deterministic genome ownership, bounded adaptive
genome packets, retry-persistent checkpoints, and machine-verifiable coverage.

Atlas and scheduler parallelism are separate concepts:

- Atlas defines scientific work as one logical task per genome.
- Sharur Ops leases those tasks dynamically to available agents.
- Sharur Query supplies one shared read-only DuckDB owner and cache.
- An executor may use persistent workers, local processes, Slurm jobs, arrays,
  or another runtime. Scheduler packing stays outside Atlas logic.

This supports Chloroflexi-scale reading directly while retaining the same
architecture for smaller datasets and future executor strategies.

## Required references

Read these before coordinating or executing Atlas:

- `docs/biological_interpretation.md`
- `docs/subagent_guide.md`
- `docs/query_service.md`
- `agent_ops_spec.md`
- `docs/review_workflow.md`
- `.claude/skills/_validation_protocols.md`

## Invariants

1. Ownership unit: one genome.
2. Traversal order: `genome_id`, then `contig_id`, then
   `gene_index NULLS LAST, start, protein_id`.
3. Model-call unit: one bounded, sequence-free packet containing consecutive
   data from exactly one `bins.bin_id`. A packet may carry several contigs.
   Oversized contigs continue across packets by stable protein offset.
4. Checkpoint unit: a completed prefix of genome packets represented by one
   opaque packet cursor.
5. Completion proof: a per-genome coverage manifest whose frame, contig,
   segment, and protein receipts equal the sealed plan.
6. Dataset identity: every plan, task, query trace, cursor, and coverage
   manifest carries the sealed `dataset_id`.
7. Full scale remains exhaustive. Dataset size changes worker count and
   campaign duration; it never silently changes Atlas into sampling.
8. Raw nucleotide and amino-acid sequences stay outside model-visible
   prompts, reports, logs, and summaries.
9. Each unit emits typed candidates first and exactly one active unit
   disposition last. Task completion reconciles both record classes.

## Observation and naming boundary

Genome packets expose three distinct evidence classes:

- `observed_annotations`: raw per-domain observations.
- `named_calls`: exact names emitted by live structured caller resources.
- `loci`: exact normalized locus memberships emitted by live caller resources.

Only `named_calls` and `loci` support their emitted names. Domain combinations
support a hypothesis phrased as “consistent with …” and tagged `UNVERIFIED`.
Inspect whatever curated resources exist in the live dataset; treat every
static list as potentially stale and incomplete.

## Plan and enqueue

Complete all dataset writes, seal the database, and build a stable plan:

```bash
sharur migrate --db data/DATASET/sharur.duckdb
sharur seal --db data/DATASET/sharur.duckdb --force
sharur verify-seal data/DATASET/dataset.seal.json

sharur-atlas plan \
  --db data/DATASET/sharur.duckdb \
  --output-dir data/DATASET/atlas \
  --packet-contigs 128 \
  --packet-proteins 500 \
  --packet-bytes 524288 \
  --query-result-bytes 2097152 \
  --checkpoint-interval-frames 1
```

Schema migration 6 installs the `(bin_id, contig_id)` navigation index used by
large exhaustive scans. Run migration in a maintenance window with query
services stopped, then reseal and restart the verified query service. A
campaign-local replica applies when it resides on a genuinely distinct storage
tier; Biotite serves the canonical immutable database directly from VAST. The
planner fails closed when this index is absent.

The plan writes:

```text
atlas/
├── plan.json
├── units.jsonl
├── packet_census/
│   ├── census.json
│   └── units/
└── coverage/
```

`units.jsonl` contains exactly one deterministic unit per live `bins.bin_id`.
Counts come from live `contigs` and `proteins` tables. Stored `bins.n_contigs`
is retained only as a diagnostic comparison.

Before launching any model worker, enumerate the exact packet stream:

```bash
sharur-atlas packet-census \
  --plan-dir data/DATASET/atlas \
  --workers 4 \
  --threads 4

sharur-atlas verify-packet-census \
  --plan-dir data/DATASET/atlas \
  --deep
```

The census invokes the real packet builder and zero models. It records exact
model invocation count, exact canonical payload bytes, transparent 2/3/4
bytes-per-token scenarios, payload-size distributions, split-contig counts,
one-bin frame proofs, and projected compact HTTP result sizes with a reserved
service envelope. Unit records are atomic and resumable. Any target-exceeding
singleton, query-result overflow, mixed-bin frame, offset gap, duplicated
segment, or count mismatch blocks enqueue.

The census parallelizes independent genome units through thread-local cursors
over one read-only DuckDB instance, buffer cache, memory ceiling, and spill
budget. `--threads` is the shared DuckDB CPU budget and `--workers` is capped
to that value. A Slurm wrapper must pass its full allocated CPU count to both
settings; retain the remote cluster preflight rules when choosing where to run
it.

Launch `sharur-ops` and sealed `sharur-query`, inspect and approve the census,
then enqueue:

```bash
sharur-atlas enqueue \
  --plan-dir data/DATASET/atlas \
  --ops-url http://ops-host:8811 \
  --query-url http://query-host:8812 \
  --scan-execution-profile atlas_scan
```

Enqueueing is idempotent by plan ID and unit ID. Each task requires the
`atlas_reader` plus `profile:atlas_scan` capabilities and one generic CPU
slot. Add or remove workers at runtime; dynamic claiming balances genomes by
completion rate. The profile name resolves to the exact provider, model, and
effort in the review policy while Atlas retains genome-level ownership.
`enqueue` verifies the matching completed packet census and attaches its hash
to the campaign and every task.

## Worker protocol

Each worker receives only its task payload. It must read the references above,
claim one `atlas_genome_read` task, and use the same per-agent credential for
Ops and Query.

```python
from sharur.ops.client import SharurOps
from sharur.query import SharurQuery

ops = SharurOps(ops_url, agent_id=agent_id, api_token=worker_token)
query = SharurQuery(query_url, api_token=worker_token)
task = ops.claim_next_task(
    campaign_id=campaign_id,
    task_types=["atlas_genome_read"],
)
```

The claim response contains the plan unit, `query_url`, coverage path,
packet-census path/hash, checkpoint key, exact packet-packing contract, and
checkpoint interval.

### Resume

Read the latest checkpoint before traversal. A checkpoint written by an
earlier attempt remains visible to the replacement attempt.

```python
checkpoint = ops.get_task_checkpoint(
    task["id"],
    task["params"]["checkpoint_key"],
)
```

A missing checkpoint means the genome begins at its first packet. The cursor
is the next `genome_packet()` position and already carries an optional
in-contig continuation. Checkpoint payload records exact completed-frame,
contig, and protein counts plus that opaque cursor.

Expired attempts lose checkpoint and terminal-write authority. Treat a fence
error as immediate evidence that another attempt owns the task.

### Exhaustive bin-scoped packet traversal

Request adaptive packets directly from the genome endpoint:

```python
packet = query.genome_packet(
    genome_id,
    cursor=packet_cursor,
    max_contigs=task["params"]["packet_contig_limit"],
    max_proteins=task["params"]["packet_protein_limit"],
    max_model_payload_bytes=task["params"]["packet_model_payload_bytes"],
    all_annotations=task["params"]["packet_all_annotations"],
)
raw = packet["raw"]
```

For every nonempty `raw.model_payload.contigs` frame:

1. Verify `raw.dataset_id`, `raw.bin_id`, and
   `raw.packing_contract_hash` against the task.
2. Make exactly one model invocation over `raw.model_payload`.
3. Read every protein record and all three evidence classes.
4. Aggregate `raw.coverage_receipts` by contig. Protein offsets must remain
   contiguous and segment indices must increase by one.
5. Record `frame_coverage_receipt(raw)` for the final coverage manifest.
6. Copy `raw.next_cursor` exactly and continue through `raw.complete=true`.

The packet builder retains whole contigs whenever they fit a fresh frame.
Only a contig whose remaining record payload exceeds the contract is split.
Every contig and receipt carries the requested bin ID, and the operator
asserts the same bin on every selected protein before serialization.

After each configured frame interval:

```python
ops.put_task_checkpoint(
    task["id"],
    checkpoint_key,
    cursor=raw["next_cursor"],
    payload={
        "completed_frames": completed_frames,
        "completed_contigs": completed_contigs,
        "completed_proteins": completed_proteins,
        "packet_cursor": raw["next_cursor"],
    },
)
ops.heartbeat_task(task["id"])
```

The default interval is one model frame. This aligns durable progress with
paid inference and confines ordinary recovery replay to the current frame.
Tune it only after the packet census quantifies the write/replay tradeoff.

## Per-genome output

Write separate files per unit. Shared append files create contention and make
recovery ambiguous.

Recommended task-local outputs:

```text
atlas/units/{unit_id}/inventory.json
atlas/units/{unit_id}/contigs.jsonl
atlas/units/{unit_id}/findings.jsonl
atlas/coverage/{unit_id}.json
```

The inventory should include:

- genome ID, plan ID, dataset ID, and unit ID;
- exact contig/protein totals;
- observed annotation census by source;
- exact caller-emitted system/locus inventory by source table;
- unusual architectures and neighborhood hypotheses;
- annotation conflicts and unresolved candidates;
- verification queries for every specific numerical claim.

Register large outputs as content-addressed Ops artifacts. Submit discrete
candidate occurrences through Ops with `campaign_id`, `task_id`, `unit_id`,
`genome_id`, and stable idempotency keys. Each candidate keeps observed
evidence separate from caller-emitted names.

### Candidate and disposition output

Emit one record per reviewable occurrence:

```python
candidate_id = ops.create_candidate_occurrence(
    campaign_id=task["campaign_id"],
    task_id=task["id"],
    dataset_id=task["params"]["dataset_id"],
    unit_id=task["params"]["unit_id"],
    genome_id=task["params"]["genome_id"],
    candidate_type="VERSIONED_TYPED_CLASS",
    signature_schema="VERSIONED_SIGNATURE_SCHEMA",
    signature=typed_signature,
    evidence=observed_evidence,
    verification=verification_specs,
    subject_refs=stable_subject_refs,
    uncertainty=uncertainty,
    reduction_features=reduction_features,
    provenance=execution_provenance,
    idempotency_key=stable_candidate_key,
)
```

The signature is structured and type-specific. It excludes interpretive
prose. A named biological claim enters evidence only through the exact live
purpose-built caller output.

After the terminal coverage manifest exists and every candidate write has
returned, append one unit disposition:

```python
ops.record_unit_disposition(
    campaign_id=task["campaign_id"],
    task_id=task["id"],
    unit_id=task["params"]["unit_id"],
    dataset_id=task["params"]["dataset_id"],
    genome_id=task["params"]["genome_id"],
    coverage_hash=coverage["coverage_sha256"],
    candidate_count=len(candidate_ids),
    disposition="candidate" if candidate_ids else "clear",
    evidence_bundle_hash=unit_evidence_bundle_hash,
    strata=audit_strata,
    provenance=execution_provenance,
    idempotency_key=stable_disposition_key,
)
```

An audit correction appends a new disposition with
`supersedes_disposition_id`. The original version remains available for
calibration.

## Coverage completion

Build the final per-genome manifest with
`write_genome_coverage_manifest(unit, contigs, frames, output_path)` from
`sharur.atlas`. Use `frame_coverage_receipt()` on each packet. Complete the
task only when its manifest has `coverage_status="complete"`, exact expected
totals, contiguous frame and contig-segment receipts, all candidate writes,
and its active unit disposition.

After all tasks finish:

```bash
sharur-atlas verify-coverage --plan-dir data/DATASET/atlas
```

The command exits successfully only when:

- every planned unit has one manifest;
- manifest plan, dataset, unit, genome, packet, and packing identities match;
- every frame belongs to the assigned bin and every contig has a terminal
  segment;
- per-genome contig and protein totals match the assigned unit;
- campaign totals match `plan.json`;
- each manifest content hash is valid.

Coverage failure is a campaign state, not a prose caveat. Requeue the affected
units and rerun verification.

## Scientific review during reading

For each protein:

- inventory all observed domain hits;
- inspect exact predicates as retrieval aids;
- preserve coordinates and stable identifiers;
- examine local context when evidence is ambiguous or biologically notable;
- route specialized candidates to the appropriate validation skill;
- state MAG non-detection with assembly fragmentation/completeness context.

Domain presence suggests compatibility with a function. Exact functional or
system names require their authoritative caller or an explicitly labeled
hypothesis. Literature claims use the literature workflow and sourced
citations.

## Campaign synthesis

Synthesis begins after coverage verification. Aggregate unit inventories with
set-oriented queries or columnar processing. Produce:

- exact coverage accounting;
- genome-level annotation and caller-resource matrices;
- cross-genome prevalence with MAG-aware denominators;
- validated recurring patterns;
- unresolved candidate classes for specialist follow-up;
- a concise record of annotation-resource blind spots.

Atlas findings may join survey findings during synthesis. Preserve their
Atlas unit ID, genome ID, dataset ID, and verification queries so every claim
can be traced back to complete genome ownership.
