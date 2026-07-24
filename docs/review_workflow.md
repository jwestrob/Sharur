# Hierarchical Scientific Review Workflow

Sharur review converts exhaustive Atlas observations into a bounded,
auditable funnel:

```text
genome scan
  -> typed candidate occurrences + one versioned unit disposition
  -> exact typed-signature reducer
  -> deep review
  -> two blind independent reviews
  -> adjudication
  -> finding materialization
  -> canonical review + executable verification
  -> publish decision + canonical publication receipt
```

Every arrow is represented by an immutable record, a normalized relationship,
or a durable Ops event. Model contexts may disappear; the scientific lineage
remains reconstructable.

## Scope boundary

Atlas owns scientific partitioning: one logical unit per genome and exhaustive
bin-scoped packet traversal within that unit. Each model frame draws from one
exact `bins.bin_id`; its packing census and coverage receipts preserve every
contig segment and protein offset. Ops owns leases, retries, review routing,
and provenance. Executors map tasks to Codex sessions, Claude agents, local
workers, Slurm jobs, arrays, or another runtime.

This separation keeps the Chloroflexota genome campaign and future scheduler
packing as distinct concepts. Changing worker packing leaves genome ownership,
candidate identity, reduction, and review policy unchanged.

## Scanner contract

Each Atlas scanner runs under an explicit semantic execution profile such as
`atlas_scan`. The default policy resolves that profile to a provider, model,
and reasoning effort. The task records both the semantic profile and its exact
resolution.

A scanner emits:

1. zero or more `candidate_occurrences`, each carrying:
   - sealed `dataset_id`, campaign, task, unit, and genome;
   - a typed `candidate_type`;
   - a versioned `signature_schema` and structured signature;
   - observed evidence, uncertainty, subject references, and exact provenance;
   - executable `{claim, query, expected}` specifications;
   - a content hash for its evidence bundle;
2. exactly one active `unit_disposition` after all candidates:
   - `clear`, `candidate`, `incomplete`, or `failed`;
   - exact candidate count;
   - coverage and evidence-bundle hashes;
   - audit strata and execution provenance.

Atlas task completion checks this contract transactionally. A successful
completion requires an active `clear` or `candidate` disposition and exact
agreement between the task’s candidate rows and disposition count.
Task-scoped idempotency keys survive lease recovery and reassignment, so a
replacement worker resolves the original occurrence/disposition records
instead of duplicating them.

The campaign metadata pins the planned Atlas `unit_count`. Canonical
publication requires exact agreement among that plan, completed
`atlas_genome_read` tasks, active unit dispositions, and ready
`clear`/`candidate` dispositions. An indexed unfinished-task probe keeps the
per-completion readiness check bounded at 29k-genome scale.

Audit corrections append a new disposition version linked by
`supersedes_disposition_id`; the earlier record becomes historical evidence.
The controller fences any clear-unit audit targeting that historical version.

## Observations and names

Candidate evidence preserves the OBSERVED/NAMED boundary:

- raw per-domain hits are observations;
- a named system, family, pathway, or mechanism requires the exact emitted
  call from a live purpose-built resource;
- compatible domain patterns remain a hypothesis marked `UNVERIFIED`.

Review tasks require live capability discovery before a named claim. Their
input contract keeps biological sequences compute-only and model outputs
sequence-free. Finding, candidate, review, and verification-result writes
reject raw-sequence fields and sequence-like strings at the storage boundary.

## Exact reduction

Run:

```bash
sharur-review reduce \
  --ops-url http://ops-host:8811 \
  --campaign-id CAMPAIGN_ID
```

The default reducer groups only this exact key:

```text
dataset_id
+ candidate_type
+ signature_schema
+ SHA-256(canonical structured signature)
```

It performs zero prose clustering and zero model inference. Every occurrence
remains in `candidate_cluster_members`. A deterministic medoid represents the
group; scanner-authored `outlier` and `counterexample` roles remain attached.
Counts include exact occurrences, genomes, units, tasks, roles, and supplied
categorical strata.

Cluster detail responses carry exact counts and a bounded member sample.
Traverse the lossless membership table through
`GET /review/clusters/{cluster_id}/members` using
`next_after_candidate_id`; each page is capped at 1,000 rows.

New occurrences create a new cluster version under the same logical cluster
ID. The prior version becomes `superseded`, and an explicit lineage edge
connects both manifests. Type-specific adapters may validate signatures,
roles, and strata while exact typed equivalence remains the grouping boundary.
The controller fences active work targeting the superseded version, including
canonical review of findings materialized from that version. Publish decisions
and publication receipts recheck materialized source versions.

## Review tiers

The packaged policy lives at
`sharur/review/default_policy.yaml`. Validate and fingerprint it with:

```bash
sharur-review policy-check
```

The default tiers are:

| Tier | Purpose | Default quorum |
|---|---|---|
| `deepen` | Reconstruct a reduced candidate from primary evidence | one verified high-effort review |
| `independent` | Challenge the candidate from separate contexts | two verified reviews from distinct providers and reviewer identities |
| `adjudication` | Resolve agreement, disagreement, and counterevidence | one verified final review |
| `canonical` | Audit the materialized finding and every publishable claim | one verified review |

Independent tasks carry empty source-review lists and blindness flags. The
adjudication task receives the frozen source-review manifest and can compare
their discrepancies explicitly.

Every `finding_review` records:

- exact provider, model, effort, execution profile, reviewer identity, and
  producing task;
- prompt, rubric, input-bundle, policy, and sealed dataset identity;
- reconstructed observations and claim assessment;
- discrepancies, proposed evidence tasks, verdict, and confidence;
- blindness flags and optional review lineage.

The store checks task ownership and exact execution identity against the
claimed task contract. Completing a scientific-review task requires exactly
one task-owned `finding_review`; a missing or duplicate review leaves the task
active. The controller advances task-backed reviews only after the terminal
task event, regardless of whether verification records arrived before or
after completion.

## Executable verification

Each review promotion carries one or more append-only
`review_verifications`. A verification records the exact specification,
expected and actual values, sealed dataset, engine, code commit, execution
time, and status.

Verification retries create numbered revisions:

```text
pending attempt 1
  -> pass/fail attempt 2
  -> corrected attempt 3
```

Publication evaluates the latest revision for every claim. Each supporting
review needs at least one verification, and every latest result must be
`pass`.

Run a bounded read-only DuckDB verification with:

```bash
sharur-review verify \
  --ops-url http://ops-host:8811 \
  --review-id REVIEW_ID \
  --claim-key genome-count \
  --db data/DATASET/sharur.duckdb \
  --dataset-id SEALED_DATASET_ID \
  --specification verification.yaml \
  --expected '29'
```

The verifier accepts one `SELECT`/`WITH` statement, bounds rows and bytes,
disables DuckDB external access, blocks filesystem/table-function operations,
opens DuckDB read-only, verifies the dataset seal by default, and rejects
model-visible sequence fields or sequence-like results.

## Negative controls and audits

The controller deterministically samples:

- 1% of active `clear` unit dispositions, stratified by configured fields;
- 5% of rejected candidate clusters.

Sampling uses a policy seed, campaign, subject, and stratum hash. Replays
select the same subjects. Clear-unit reviews target the exact disposition
version. Rejected-case audits receive a blind task with the candidate target
and can reopen the review path when verified contrary evidence appears.

These audits estimate scanner false-negative behavior and reducer/reviewer
false-dismissal behavior from the same provenance graph as promoted cases.

## Controller and backpressure

Run one controller against the Ops HTTP owner:

```bash
sharur-review route \
  --ops-url http://ops-host:8811 \
  --campaign-id CAMPAIGN_ID \
  --watch
```

The controller consumes the durable event stream through a monotonic cursor.
Every derived task and decision has a content-derived idempotency key. A crash
between task creation and cursor advancement safely replays the event.
Cluster-supersession events cancel stale pending or leased work and invalidate
its lease token, concentrating review capacity on the current immutable
manifest.

Execution profiles define indexed pending-task ceilings. Reaching a ceiling
holds the cursor on that event and reports the constrained profile, providing
queue backpressure while preserving the work item.

The controller pins its policy hash on first use. A policy revision uses a new
`controller_id`, which creates an explicit versioned replay over historical
events and keeps decisions from different policies distinguishable.

The controller creates logical tasks only. A Codex headless worker, Claude
agent pool, or another executor registers the corresponding
`profile:<name>` capability and claims matching work.

## Finding materialization and publication

Successful adjudication creates a `materialize_finding` task. Its worker:

1. creates one immutable Ops finding;
2. links the source cluster with relation `materializes`;
3. leaves `validation_status` pending review.

Task completion verifies exactly one task-owned result finding, the typed
source-cluster link, and the expected pending validation status. The link
and terminal task events jointly create the canonical review task through an
idempotent route. A verified canonical promotion creates a `publish`
decision. The canonical JSONL writer then performs its strict locked write
and records a `canonical_publication` receipt containing the canonical URI,
record ID, and content hash. Both the decision and receipt recheck Atlas
closure and current source-cluster versions.

Ops remains the operational and review provenance graph. Canonical JSONL
remains the durable scientific archive.

## Inspection

```bash
sharur-review status \
  --ops-url http://ops-host:8811 \
  --campaign-id CAMPAIGN_ID

sharur-review trace \
  --ops-url http://ops-host:8811 \
  --campaign-id CAMPAIGN_ID \
  --subject-kind candidate_cluster \
  --subject-id CLUSTER_ID
```

`status` reports exact coverage, reduction, audit, verification, promotion,
publication, and per-profile queue counts. `trace` reconstructs a bounded
unit-, cluster-, or finding-centered lineage with reviews, verification
revisions, decisions, materialization links, and publication receipts.
