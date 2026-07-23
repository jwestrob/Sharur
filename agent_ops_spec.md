# Sharur Ops v3 — Coordination Layer Specification

## Purpose

Sharur Ops externalizes multi-agent coordination into a durable, queryable
control plane. Agents exchange structured state through the service while
their model contexts remain independent:

```text
workers/coordinator
        |
        | HTTP + per-agent bearer identity
        v
one sharur-ops service
        |
        | bounded SQLite connection pool
        v
dataset-local sharur_ops.db
```

The service owns coordination state. `sharur-query` owns shared analytical
access to large dataset DuckDB files, and canonical research outputs remain in
their versioned finding/report artifacts.

## Scope boundary

Sharur Ops coordinates logical work. Compute launchers execute that work.

- The control plane represents campaigns, dependencies, capabilities,
  resource requests, attempts, leases, results, and events.
- An executor maps a claimed task onto a local process, persistent worker,
  Slurm job, Slurm array element, Kubernetes job, or another runtime.
- Pipeline stage algorithms and scheduler packing remain independent.

The analytical data-plane contract lives in
[`docs/query_service.md`](docs/query_service.md). Sharur Ops resolves agent
identity through `/auth/whoami`; the query service consumes that endpoint and
keeps SQLite ownership inside the Ops HTTP process.

This separation supports both large one-off analyses and general
multi-agent use while keeping scheduler-specific behavior out of the
pipeline's scientific logic.

## Deployment

Install and launch:

```bash
pip install -e ".[ops]"

export SHARUR_OPS_DB_PATH=/absolute/path/data/DATASET/sharur_ops.db
export SHARUR_OPS_TOKEN="$(python -c 'import secrets; print(secrets.token_urlsafe(32))')"
sharur-ops --host 0.0.0.0
```

Defaults:

| Setting | Default |
|---|---:|
| Bind address | `127.0.0.1:8811` |
| SQLite pool | 4 connections |
| SQLite busy timeout | 15 seconds |
| Request body limit | 1 MiB |
| Inline JSON/text limit | 256 KiB |
| Task lease | 900 seconds |
| Maintenance cadence | 15 seconds |
| Running-run recovery threshold | 300 seconds |
| Uvicorn workers | 1 |

The service acquires `.<database-name>.server.lock` for its lifetime. A second
service process targeting the same database fails startup. This single-owner
topology also governs NFS deployments: every SQLite connection belongs to one
service host, and every remote agent uses HTTP.

Direct `OpsStore` access is the local mode. Direct clients take a shared
advisory lock, while the HTTP service takes the exclusive form of that lock;
server startup fails while direct clients are active, and direct-client
startup fails while the service is active.

## Identity and roles

`SHARUR_OPS_TOKEN` authenticates a bootstrap operator. The operator registers
long-lived agent identities:

```python
from sharur.ops.client import SharurOps

operator = SharurOps(
    "http://ops-host:8811",
    agent_id="bootstrap",
    api_token=bootstrap_token,
)
registration = operator.register_agent(
    "dpann-reader-03",
    role="worker",
    capabilities=["duckdb", "synteny"],
    max_concurrent_tasks=2,
    capacity_cpu_slots=8,
    capacity_memory_mb=16384,
)
worker_token = registration["token"]  # returned once
```

The server stores a token hash and derives identity from the credential on
every request.

| Role | Authority |
|---|---|
| `reader` | Read coordination records and observability endpoints |
| `worker` | Create findings/hypotheses and mutate its fenced task attempts |
| `coordinator` | Worker authority plus task creation/recovery and reasoning logs |
| `operator` | Campaign/run/agent administration, backups, and integrity checks |

An identity mismatch returns HTTP 403. Operator credentials retain explicit
act-as authority for bootstrap and legacy client compatibility.

## Campaign model

A campaign namespaces one coordinated analysis and links:

- findings and hypotheses;
- tasks and runs;
- coordinator decisions;
- artifacts;
- durable events.

Campaign filtering keeps concurrent analyses from competing in one logical
queue while preserving one operational database.

## Task authority protocol

### Creation and reservation

Task creation accepts:

```json
{
  "task_type": "query",
  "description": "Compute a genome-level marker matrix",
  "campaign_id": "campaign-uuid",
  "depends_on": ["task-uuid"],
  "assigned_to": "worker-07",
  "required_capabilities": ["duckdb"],
  "resource_request": {
    "cpu_slots": 4,
    "memory_mb": 8192,
    "accelerator_slots": 0
  },
  "max_attempts": 3,
  "lease_seconds": 900,
  "idempotency_key": "marker-matrix:v1"
}
```

`assigned_to` creates a pending reservation. The named worker explicitly
claims it.

### Claim

`POST /tasks/claim-next` selects the highest-priority compatible task inside
one writer transaction:

1. recover expired leases;
2. promote due retries;
3. filter campaign/task type;
4. require completed normalized dependencies;
5. enforce reservations;
6. match registered capabilities;
7. account for concurrent tasks and active CPU/memory/accelerator allocation;
8. increment the attempt and mint a cryptographically random lease token.

The response carries:

```json
{
  "id": "task-uuid",
  "status": "claimed",
  "attempt_count": 2,
  "lease_attempt": 2,
  "lease_token": "opaque-attempt-secret",
  "lease_expires_ts": 1780000000.0
}
```

The database stores only the token hash.

### Heartbeat and terminal transition

Heartbeat, completion, and failure submit the token plus attempt number. The
conditional write checks:

```text
task ID
+ authenticated agent identity
+ attempt number
+ token hash
+ active status
+ server-clock lease expiry
```

This is the fencing boundary. A process from attempt 1 loses write authority
after attempt 2 begins, including when both processes use the same agent ID.

Retryable failure enters `retry_wait`. Recovery returns it to its original
reservation after the delay. Exhausted tasks enter `failed`.

HTTP recovery always uses the server clock. This removes caller-controlled
expiration.

## Resource-aware analysis

Dispatch resources are generic integer slots:

- `cpu_slots`;
- `memory_mb`;
- `accelerator_slots`.

Capabilities describe qualitative requirements. An Apple Silicon worker can
register an `mps` capability with one accelerator slot, enforcing one active
accelerated task through the same generic capacity model.

Local direct-analysis processes receive an explicit DuckDB budget:

```python
from sharur.operators import Sharur

b = Sharur(
    database_path,
    read_only=True,
    duckdb_threads=4,
    duckdb_memory_limit="8GB",
    duckdb_temp_directory="/local-scratch/worker-07",
)
```

The scheduler/executor chooses these values from the claimed task and host
allocation. Coordinated campaigns over a large shared database route typed
operators through `sharur-query`; one service-level budget then governs all
agents. Slurm array construction remains an executor concern.

## Durable writes and idempotency

Caller-stable idempotency keys cover:

- campaigns;
- findings;
- hypotheses;
- tasks;
- runs;
- coordinator-log entries.

An exact retry returns the existing ID. Reusing a key with a different
immutable payload returns a conflict. The client retries a mutation after a
transport failure only when that mutation carries an idempotency key or a
natural content identity.

Scientific references use normalized foreign-key tables:

- `task_dependencies`;
- `task_result_findings`;
- `hypothesis_findings`;
- `finding_links`;
- `coordinator_log_findings`;
- `coordinator_log_hypotheses`;
- `finding_artifacts`.

Legacy JSON arrays remain readable and synchronized.

## Events and streaming

Every coordination mutation appends a `coordination_events` row in the same
transaction as the state change. Event IDs are monotonic.

```http
GET /events?after_id=420&campaign_id=...
GET /stream?after_id=420
Last-Event-ID: 420
```

SSE records use the durable event ID:

```text
id: 421
event: finding_created
data: {"id":421,"entity_type":"finding",...}
```

The in-process event bus wakes subscribers. Clients recover any missed wake-up
from SQLite by cursor, so bounded queue eviction preserves delivery.

## Artifact boundary

Large output bytes live in file/object storage. Ops registers:

- algorithm-qualified content hash;
- URI;
- exact byte size;
- media type;
- creator and timestamp;
- campaign/task/run provenance;
- metadata.

The content hash is unique. Findings link artifacts through
`finding_artifacts`. Inline records retain a 256 KiB ceiling.

## Connection and migration lifecycle

Application startup:

1. acquires the owner lock;
2. opens the first SQLite connection;
3. applies additive schema migration v3 inside one transaction;
4. opens the remaining pool connections;
5. starts maintenance.

Each store operation checks out one pool connection. `OpsStore` and
`RunLedger` share that connection and lock, eliminating the former nested
ledger connection and duplicate schema initialization. A current v3 schema
open produces zero schema writes.

`PRAGMA foreign_keys=ON`, WAL, `synchronous=NORMAL`, a 15-second busy timeout,
and WAL autocheckpointing apply to every pool connection.

## API

### Campaigns and agents

| Method | Endpoint | Purpose |
|---|---|---|
| POST | `/campaigns` | Create an idempotent campaign |
| GET | `/campaigns` | List/filter campaigns |
| GET/PATCH | `/campaigns/{id}` | Inspect or change campaign status |
| POST | `/agents` | Register/update/rotate an agent |
| GET | `/agents` | List registered agents |
| POST | `/agents/me/heartbeat` | Agent presence/capacity heartbeat |
| PATCH | `/agents/{id}` | Drain, offline, activate, or revoke |

### Scientific state

| Method | Endpoint | Purpose |
|---|---|---|
| POST/GET | `/findings` | Create or list findings |
| GET | `/findings/search/{text}` | FTS5-backed search with fallback |
| GET | `/findings/{id}` | Fetch one finding |
| POST | `/hypotheses` | Create an idempotent hypothesis |
| GET/PATCH | `/hypotheses[/{id}]` | List or update hypotheses |
| POST/GET | `/artifacts[/{id}]` | Register or fetch artifact metadata |
| POST | `/findings/{id}/artifacts` | Attach an artifact |

### Tasks and runs

| Method | Endpoint | Purpose |
|---|---|---|
| POST/GET | `/tasks` | Create or list tasks |
| POST | `/tasks/claim-next` | Atomic compatible queue claim |
| POST | `/tasks/{id}/claim` | Atomic specific-task claim |
| POST | `/tasks/{id}/heartbeat` | Renew a fenced attempt |
| PATCH | `/tasks/{id}` | Complete, fail, or heartbeat a fenced attempt |
| POST | `/tasks/recover` | Coordinator/operator server-clock recovery |
| POST/GET | `/runs` | Create or list durable runs |
| GET | `/runs/{id}` | Fetch a run |
| POST | `/runs/{id}/start` | Start a run |
| POST | `/runs/{id}/submit` | Record scheduler handoff |
| POST | `/runs/{id}/heartbeat` | Run heartbeat |
| PATCH | `/runs/{id}` | Complete or fail a run |
| GET | `/runs/{id}/events` | Cursor-paged run events |
| GET | `/runs/{id}/stages` | Inspect attempts |
| POST | `/runs/recover` | Coordinator/operator stale-run recovery |

### Coordination and operations

| Method | Endpoint | Purpose |
|---|---|---|
| POST/GET | `/coordinator/log` | Persist/read reasoning |
| GET | `/events` | Durable global cursor replay |
| GET | `/stream` | Replayable SSE tail |
| GET | `/stats` | Structured operational metrics |
| GET | `/metrics` | Prometheus metrics |
| GET | `/health` | Schema/owner/pool health |
| POST | `/admin/backup` | Online verified backup |
| GET | `/admin/integrity` | SQLite and FK integrity |

Time-ordered list endpoints accept `before_ts`; event endpoints use
`after_id`. Every endpoint enforces a bounded `limit`.

## Recovery and backup

Maintenance runs from the service clock. It recovers expired task attempts and
stale running runs/stages.

Enable periodic online backups:

```bash
export SHARUR_OPS_BACKUP_DIR=/durable/ops-backups
export SHARUR_OPS_BACKUP_INTERVAL_SECONDS=3600
```

Each backup uses SQLite's online backup API and then runs
`PRAGMA integrity_check`. Filenames include nanosecond-resolution creation
times. Successful and failed backup attempts update service health metrics;
the scheduled path advances its last-success timestamp only after the
integrity check passes.

## Observability

`/stats` and `/metrics` expose:

- counts by record and lifecycle status;
- oldest queued-task age;
- live and expired leases;
- terminal failed-task count, attempts-exhausted count, and oldest
  dead-letter age;
- per-agent active CPU/memory/accelerator allocation;
- database, WAL, and shared-memory byte counts;
- connection-pool acquisition count/time and SQLite writer-slot wait time;
- HTTP request/error totals, `SQLITE_BUSY` count, duration histogram, and
  rolling p95 latency over the latest 4,096 requests;
- backup last-success age and failure count;
- SSE subscriber and dropped-wakeup counts. Durable cursor replay preserves
  events across wake-up eviction.

Alerts should cover:

- `expired_leases > 0` over two maintenance intervals;
- pool `wait_count` growth;
- rolling p95 HTTP latency above 100 ms;
- any `SQLITE_BUSY`;
- backup age above twice the configured interval;
- integrity failures.

## Backend scaling boundary

SQLite plus one HTTP owner is a deliberate first production tier. It provides
simple recovery, online backup, and adequate serialized-write throughput for
the intended low-rate scientific coordination workload.

Move the same API contract to PostgreSQL when deployment needs API replicas,
high availability, multi-host writers, or sustained write concurrency beyond
the thresholds in `sharur/ops/SCHEMA.md`. Bulk analytical results continue to
live outside the control-plane database in either backend.

## Validation

Focused contracts cover:

- same-ID stale-attempt fencing;
- atomic concurrent claims and claim-next selection;
- capacity/capability enforcement;
- dependency and reference integrity;
- concurrent hypothesis evidence updates;
- v0/v2 migration and current-schema zero-write opens;
- per-agent credentials, role denial, and identity spoof rejection;
- server-clock recovery;
- durable event cursor replay;
- request and inline payload limits;
- streamed/chunked request limits;
- FTS and content-addressed artifacts;
- bounded 5,000-task queue reads;
- compatible dispatch after more than 100 higher-priority incompatible tasks;
- event replay across service restart and visible SSE wake-up overflow;
- verified online backup;
- dead-letter and operational metric visibility;
- exclusive HTTP ownership;
- direct-process DuckDB budgets and shared query-service admission/resource budgets;
- concurrent run-stage summary coherence.

Run:

```bash
python -m pytest \
  tests/test_ops_store.py \
  tests/test_ops_http.py \
  tests/test_ops_control_plane.py \
  tests/test_duckdb_resource_budget.py
```
