# Sharur Ops v5 — Storage, Coordination, and Review Contract

Sharur Ops is the durable control plane for multi-agent analyses. It records
campaigns, identities, logical work, scientific findings, hypotheses,
execution attempts, artifacts, and an ordered event stream in one SQLite
database.

## Ownership topology

Use one of these two modes:

1. Local, single-process work uses `OpsStore` directly.
2. Distributed work runs one `sharur-ops` HTTP server, and every agent uses
   `SharurOps`.

```python
# Local single-process mode
from sharur.ops.store import OpsStore

ops = OpsStore("data/DATASET/sharur_ops.db", agent_id="local-operator")

# Distributed mode
from sharur.ops.client import SharurOps

ops = SharurOps(
    "http://ops-host:8811",
    agent_id="worker-07",
    api_token=worker_token,
)
```

The server holds an advisory owner lock, initializes schema migrations once,
and serves requests through a bounded connection pool. This makes a
dataset-local SQLite file on NFS safe under the intended topology: every
SQLite connection originates in one server process on one host. Route all
remote clients through HTTP. Local direct clients take shared locks; the
server takes the exclusive lock, enforcing a single access topology.

SQLite remains the target backend for one service instance and tens of
low-rate agents. Select PostgreSQL when the deployment requires any of:

- multiple API replicas or automatic failover;
- database writers on multiple hosts;
- sustained lock waits or `SQLITE_BUSY` errors above 0.1% of writes;
- sustained p95 mutation latency above 100 ms after query/index tuning;
- hundreds of simultaneously active workers or a sustained write rate above
  roughly 10 mutations per second.

These are operational migration triggers, not logical-schema changes. The
HTTP API is the portability boundary.

## Coordination and execution are separate concepts

An Ops task is logical work with dependencies, requirements, ownership, and
results. Its executor can be a local process, a long-lived worker, a Slurm
job, a Slurm array element, or another scheduler adapter. The schema carries
generic capabilities and resources:

```json
{
  "required_capabilities": ["duckdb", "annotation-review"],
  "resource_request": {
    "cpu_slots": 4,
    "memory_mb": 8192,
    "accelerator_slots": 0
  }
}
```

Scheduler packing, array construction, and job submission live in executor
adapters. They do not define task semantics or pipeline stage behavior.

## Schema lifecycle

`OPS_SCHEMA_VERSION = 5`. `ensure_ops_schema()` performs an additive,
transactional migration and returns:

- `True` when it applies a migration;
- `False` when the database already has the current schema.

A current-schema open performs schema reads and connection PRAGMAs only.
Legacy JSON reference columns remain available for compatibility. New writes
also populate normalized, foreign-key-backed relationship tables.

Review-DAG output keys are task-scoped when `task_id` is present, preserving
exactly-once logical output across lease recovery and worker reassignment.
Task completion for scientific review and finding materialization validates
the typed output rows before the terminal transition.

## Primary tables

### `campaigns`

A namespace for one coordinated analysis. Fields include identity, name,
dataset path, lifecycle status, creator, timestamps, metadata, and an
idempotency key. Findings, hypotheses, tasks, runs, logs, events, and
artifacts can carry `campaign_id`.

### `agents`

Server-derived identities and dispatch capacity:

| Field | Meaning |
|---|---|
| `id` | Credential identity |
| `role` | `reader`, `worker`, `coordinator`, or `operator` |
| `token_hash`, `token_hint` | Hashed per-agent bearer credential and display hint |
| `capabilities` | JSON array of dispatch labels |
| `max_concurrent_tasks` | Active-task ceiling |
| `capacity_cpu_slots` | Generic CPU capacity |
| `capacity_memory_mb` | Generic memory capacity |
| `capacity_accelerator_slots` | Generic accelerator capacity |
| `host`, `status`, `heartbeat_ts` | Presence and drain/revoke state |
| `metadata` | Structured executor metadata |

Raw agent tokens are returned once at registration or rotation. The database
stores their SHA-256 hashes.

### `findings`

Immutable scientific records. Core fields are:

- producer, timestamp, type, domain, summary, reasoning;
- structured `evidence`;
- confidence and novelty;
- parent, campaign, and producing task;
- caller-stable idempotency key;
- record schema version and validation status;
- ID of the durable creation event.

The inline JSON ceiling is 256 KiB. Larger outputs enter `artifacts` and are
referenced by ID. An FTS5 index covers summary, reasoning, and evidence when
the SQLite build provides FTS5.

### `hypotheses`

Testable claims with status, assignment, relevant domains, source findings,
supporting findings, opposing findings, and resolution. The JSON arrays
remain the compatibility view; `hypothesis_findings` is the normalized
relationship table.

Hypothesis evidence updates acquire one writer transaction across read,
deduplication, normalized inserts, JSON-view update, and event append. This
preserves concurrent writers.

### `tasks`

The logical work queue. Important fields include:

| Field | Meaning |
|---|---|
| `status` | `pending`, `claimed`, `in_progress`, `retry_wait`, `complete`, `failed`, or `cancelled` |
| `assigned_to` | Reservation while pending; lease owner while active |
| `reserved_for` | Persistent reservation across retries |
| `priority` | 0–3 |
| `params` | Bounded task-specific JSON |
| `execution_profile` | Indexed semantic execution profile for routing/backpressure |
| `run_id`, `campaign_id` | Optional execution and analysis parents |
| `attempt_count`, `max_attempts` | Bounded retry accounting |
| `lease_seconds`, `lease_expires_ts` | Finite authority window |
| `lease_token_hash` | Hash of the attempt-specific fencing token |
| `required_capabilities` | JSON capability set |
| `resource_request` | Generic CPU/memory/accelerator request |

`task_dependencies` and `task_result_findings` are the authoritative
relationship tables. The legacy JSON arrays remain synchronized for clients
that still read them.

Preassignment creates a pending reservation. The reserved worker performs an
explicit claim and receives:

```json
{
  "lease_token": "opaque-secret",
  "lease_attempt": 2
}
```

Every heartbeat, completion, and failure presents both values. The mutation
predicate checks task ID, credential identity, attempt number, token hash,
status, and live expiration. A stale process using the same agent ID therefore
loses authority when a replacement attempt begins.

`claim_next_task()` performs recovery, retry promotion, dependency filtering,
capability matching, capacity accounting, selection, and claim inside one
`BEGIN IMMEDIATE` transaction. `available_tasks()` is a bounded read-only
queue view.

### `task_checkpoints`

One compact retry-persistent checkpoint per `(task_id, checkpoint_key)`.
Rows record the writing attempt, agent identity, opaque cursor, bounded JSON
payload, and update time. Upserts use the same active lease fence as task
completion. A replacement attempt can read the prior row and then advance it
under its new token and attempt number.

Checkpoint rows update in place and survive task retry. Workers batch updates
at scientifically safe recovery boundaries to control central write volume.

### Review-DAG tables

Ops v5 adds normalized scientific escalation state:

| Table | Authority |
|---|---|
| `unit_dispositions` | Versioned terminal result for each scanned unit |
| `candidate_occurrences` | Lossless typed observations and evidence bundles |
| `candidate_clusters` | Immutable reducer versions and exact manifests |
| `candidate_cluster_members` | Every occurrence with member/medoid/outlier/counterexample role |
| `candidate_cluster_lineage` | Supersession, split, merge, and refinement edges |
| `cluster_findings` | Materialization/support/counterexample links to findings |
| `finding_reviews` | Exact reviewer/execution identity, assessment, and verdict |
| `review_verifications` | Append-only executable claim-check attempts |
| `promotion_decisions` | Versioned policy decisions and linked reviews/tasks |
| `canonical_publications` | Canonical JSONL publication receipts |
| `review_controller_state` | Durable event cursor per policy controller/campaign |

Publication requires a matching `publish` decision, at least one verification
per supporting review, and a latest `pass` result for every recorded
claim/specification pair. Atlas campaigns also require exact closure across
the planned unit count, completed scan tasks, and active ready dispositions.
See `docs/review_workflow.md` for the end-to-end contract.

### `runs`, `run_stages`, and `run_events`

Runs describe durable executions. Stage rows are keyed by
`(run_id, stage_id, attempt)` and record signatures, commands, inputs,
outputs, resource profiles, ownership, heartbeat, and reuse origin.

The attempt number fences stale stage completion. Multiple distinct stages
can execute concurrently; run summary state retains one currently active
stage until the active set becomes empty. `run_events` provides an ordered
run-local history.

### `coordinator_log`

The coordinator's durable reasoning trail. Normalized link tables connect log
entries to findings and hypotheses. Idempotency keys make lost-response retry
safe.

### `artifacts` and `finding_artifacts`

Artifact rows register metadata; the HTTP service does not proxy bulk bytes.
Each artifact has an algorithm-qualified content hash, URI, byte size, media
type, creator, timestamps, optional campaign/task/run, and metadata. The
content hash is unique.

### `coordination_events`

The global append-only event stream. Each row has a monotonically increasing
integer ID, timestamp, actor, campaign, entity type/ID, optional task/run, and
JSON payload. Lifecycle and scientific-state mutations append their event in
the same SQLite transaction as the state change. High-rate task checkpoints
remain durable in `task_checkpoints` and use the in-process bus as a bounded
wake-up signal, limiting append-only event growth.

`GET /events?after_id=N` provides cursor replay. `GET /stream` emits SSE
records with the durable event ID and honors `Last-Event-ID`. The in-process
bus is only a wake-up signal; queue eviction cannot erase durable events.

## Authentication and authorization

`SHARUR_OPS_TOKEN` is the bootstrap operator credential. Use it to register
per-agent credentials, then give each worker its own token.

| Role | Core authority |
|---|---|
| `reader` | Read records, events, health, and metrics |
| `worker` | Write findings/hypotheses; claim, heartbeat, complete, and fail owned tasks |
| `coordinator` | Worker rights plus create/recover tasks and write coordinator logs |
| `operator` | Full administration, campaigns, runs, agents, backups, and integrity checks |

The credential establishes `agent_id` server-side. Operator credentials may
act on behalf of a named principal for bootstrap and compatibility workflows.
Other roles receive HTTP 403 for an identity mismatch.

Recovery endpoints use the server clock. Direct `OpsStore` accepts an injected
`now` only for deterministic local tests and offline administrative tooling.

## Limits, observability, and recovery

- HTTP request body default: 1 MiB.
- Inline JSON/text default: 256 KiB.
- List endpoints: bounded `limit`; time-ordered collections accept
  `before_ts`.
- Candidate-cluster detail: exact aggregate counts plus at most 1,000 sampled
  members; `/review/clusters/{cluster_id}/members` provides bounded cursor
  traversal of the lossless membership table.
- Event replay: bounded `after_id` cursor.
- SQLite request pool: four connections by default.
- Maintenance: task/run recovery every 15 seconds by default.
- Metrics: record/status counts, oldest queue age, active and expired leases,
  failed/dead-letter tasks, per-agent active resource allocation,
  database/WAL sizes, pool acquisition and SQLite writer-slot waits, rolling
  HTTP p95 latency, `SQLITE_BUSY`, SSE wake-up overflow, and backup
  age/failures.
- `GET /metrics`: Prometheus text.
- `GET /admin/integrity`: SQLite and foreign-key checks.
- `POST /admin/backup`: online SQLite backup with integrity verification.

Scheduled backups use:

```bash
export SHARUR_OPS_BACKUP_DIR=/durable/backup/path
export SHARUR_OPS_BACKUP_INTERVAL_SECONDS=3600
```

Backup filenames carry nanosecond-resolution timestamps. The online SQLite
backup API produces each copy, `PRAGMA integrity_check` verifies it, and the
service records last-success age only after verification succeeds.

## Per-agent DuckDB budgets

Analytical resource limits are explicit on each DuckDB connection:

```python
from sharur.operators import Sharur

b = Sharur(
    "data/DATASET/sharur.duckdb",
    read_only=True,
    duckdb_threads=4,
    duckdb_memory_limit="8GB",
    duckdb_temp_directory="/local-scratch/agent-07",
)
```

These settings prevent multiple analysis agents from each inheriting
host-scale DuckDB defaults. They remain independent of the scheduler used to
launch those agents.

## Operational query examples

```sql
-- Queue head with normalized dependency filtering
SELECT t.*
FROM tasks AS t
WHERE t.status = 'pending'
  AND NOT EXISTS (
      SELECT 1
      FROM task_dependencies AS d
      JOIN tasks AS parent ON parent.id = d.depends_on_task_id
      WHERE d.task_id = t.id AND parent.status <> 'complete'
  )
ORDER BY t.priority DESC, t.ts, t.id
LIMIT 50;

-- Active authority
SELECT id, assigned_to, attempt_count, lease_expires_ts
FROM tasks
WHERE status IN ('claimed', 'in_progress')
  AND lease_expires_ts > unixepoch('subsec');

-- Replay one campaign
SELECT *
FROM coordination_events
WHERE campaign_id = ? AND id > ?
ORDER BY id
LIMIT 1000;
```

Task claims and terminal transitions always use `OpsStore` or the HTTP API;
application-level token generation and transactional predicates are part of
the authority contract.
