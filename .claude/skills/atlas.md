# Atlas Skill

Exhaustive bottom-up reading of every genome, contig, and protein in a sealed
Sharur dataset. Atlas uses deterministic genome ownership, bounded contig
packets, retry-persistent checkpoints, and machine-verifiable coverage.

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
3. Model-call unit: one bounded, sequence-free contig packet.
4. Checkpoint unit: a completed prefix of contigs, with optional in-contig
   packet state for unusually large contigs.
5. Completion proof: a per-genome coverage manifest whose contig and protein
   totals equal the sealed plan.
6. Dataset identity: every plan, task, query trace, cursor, and coverage
   manifest carries the sealed `dataset_id`.
7. Full scale remains exhaustive. Dataset size changes worker count and
   campaign duration; it never silently changes Atlas into sampling.
8. Raw nucleotide and amino-acid sequences stay outside model-visible
   prompts, reports, logs, and summaries.
9. Each unit emits typed candidates first and exactly one active unit
   disposition last. Task completion reconciles both record classes.

## Observation and naming boundary

Contig packets expose three distinct evidence classes:

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
  --packet-proteins 100 \
  --checkpoint-interval-contigs 25
```

Schema migration 6 installs the `(bin_id, contig_id)` navigation index used by
large exhaustive scans. Run migration in a maintenance window with query
services stopped, then reseal and restage. The planner fails closed when this
index is absent.

The plan writes:

```text
atlas/
├── plan.json
├── units.jsonl
└── coverage/
```

`units.jsonl` contains exactly one deterministic unit per live `bins.bin_id`.
Counts come from live `contigs` and `proteins` tables. Stored `bins.n_contigs`
is retained only as a diagnostic comparison.

Launch `sharur-ops` and a sealed `sharur-query` replica, then enqueue:

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
checkpoint key, packet limit, and checkpoint interval.

### Resume

Read the latest checkpoint before traversal. A checkpoint written by an
earlier attempt remains visible to the replacement attempt.

```python
checkpoint = ops.get_task_checkpoint(
    task["id"],
    task["params"]["checkpoint_key"],
)
```

A missing checkpoint means the genome begins at its first contig. The cursor
is the next `list_contigs()` page position. Checkpoint payload records exact
completed-contig and completed-protein counts plus any current large-contig
packet cursor.

Expired attempts lose checkpoint and terminal-write authority. Treat a fence
error as immediate evidence that another attempt owns the task.

### Exhaustive contig traversal

Request contig pages using `limit=checkpoint_interval_contigs`. This makes a
successful page the durable checkpoint batch and bounds recovery replay to at
most one page.

```python
page = query.list_contigs(
    genome_id,
    limit=checkpoint_interval_contigs,
    cursor=contig_cursor,
)
```

For every contig record in the page:

1. Call `query.get_contig(genome_id, contig_id)` for exact metadata.
2. Start `packet_cursor=None`.
3. Repeatedly call `query.contig_packet(...)`.
4. Read every returned protein record and every evidence class.
5. Persist a compact contig inventory in the task-local output.
6. Continue until `raw.complete` is true and `raw.next_cursor` is null.
7. Record `contig_id`, exact `protein_count`, packet count, and
   `complete=true` for the coverage manifest.

```python
packet = query.contig_packet(
    genome_id,
    contig_id,
    cursor=packet_cursor,
    limit=packet_protein_limit,
    all_annotations=True,
)
packet_cursor = packet["raw"]["next_cursor"]
```

Packet cursors are opaque and scoped to the exact dataset, genome, and contig.
Copy them exactly. For a contig spanning many packets, save current packet
state periodically under the same task checkpoint; this bounds replay while
preserving a single genome owner.

After every completed contig page:

```python
ops.put_task_checkpoint(
    task["id"],
    checkpoint_key,
    cursor=page["ref"],
    payload={
        "completed_contigs": completed_contigs,
        "completed_proteins": completed_proteins,
        "current_contig_id": None,
        "packet_cursor": None,
    },
)
ops.heartbeat_task(task["id"])
```

The default interval of 25 reduces a 10.7-million-contig campaign to roughly
428,000 central checkpoint updates while capping normal replay at 24 completed
contigs. Tune the interval from the plan when contig-size distributions justify
a different recovery/write tradeoff.

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
`write_genome_coverage_manifest()` from `sharur.atlas`. Complete the task only
when its manifest has `coverage_status="complete"`, exact expected totals, all
candidate writes, and its active unit disposition.

After all tasks finish:

```bash
sharur-atlas verify-coverage --plan-dir data/DATASET/atlas
```

The command exits successfully only when:

- every planned unit has one manifest;
- manifest plan, dataset, unit, genome, and packet identities match;
- every contig has a terminal packet;
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
