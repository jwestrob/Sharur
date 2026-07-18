# Sharur Ops — Coordination Layer Spec

## Problem

Sharur is a multi-agent metagenomic discovery system. 10-12 AI agents run concurrently, each investigating microbial lineage databases (Omnitrophota, DPANN, Bathyarchaeia, Hinthialibacterota) for novel biology. A coordinator agent orchestrates them.

The current failure mode: all agents report back through the coordinator's context window. When many agents dump results simultaneously, the coordinator's context compacts, destroying its reasoning state. Anything not externalized is lost.

## Solution

A thin FastAPI + SQLite coordination server that externalizes all inter-agent communication into a queryable database. Agents write structured findings to the server via HTTP. The coordinator queries the server instead of holding everything in context. The coordinator also writes its own reasoning to the server, so it survives compaction.

No agent ever writes to another agent's context. All coordination flows through the ops server.

## Architecture

```
┌─────────────┐      POST /findings
│  Agent 1    │─────────────────────┐
└─────────────┘                     │
┌─────────────┐      POST /findings │     ┌──────────────────┐     ┌────────────────┐
│  Agent 2    │─────────────────────┼────▶│ sharur.ops.server│────▶│ sharur_ops.db  │
└─────────────┘                     │     │  FastAPI + SQLite │     │ (WAL mode)     │
┌─────────────┐      POST /findings │     └──────────────────┘     └────────────────┘
│  Agent N    │─────────────────────┘           ▲    │
└─────────────┘                                │    │ GET /findings
                                               │    │ GET /stream (SSE)
                                               │    ▼
                                         ┌─────────────┐
                                         │ Coordinator  │
                                         └─────────────┘
```

Agents also read from the server (checking what others have found, polling for tasks). The coordinator reads findings, assigns tasks, proposes hypotheses, and logs its own reasoning — all via HTTP.

Domain databases (DuckDB files) remain read-only reference data. Agents ATTACH them directly
for analytical queries. The ops server handles only coordination state; canonical scientific
records remain the dataset findings and hypothesis files.

## Dependencies

```
pip install -e .
```

No other dependencies. Python 3.10+.

## Running

```bash
# Start the server (creates sharur_ops.db on first run)
uvicorn sharur.ops.server:app --host 0.0.0.0 --port 8811

# Or directly
python -m sharur.ops.server
```

Works identically on laptop (localhost:8811) and cluster (node:8811). Agents just change the base URL.

The server defaults to `sharur_ops.db` in the process working directory. Pin it to the
dataset-local ledger created by `sharur-ingest` so coordination tasks and execution history
share one operational record:

```bash
export SHARUR_OPS_DB_PATH=/absolute/path/data/DATASET/sharur_ops.db
uvicorn sharur.ops.server:app --host 0.0.0.0 --port 8811
```

## Database

SQLite runs in WAL mode with a busy timeout. The HTTP server serializes writes through
an in-process lock; direct `OpsStore` instances use SQLite transactions plus a per-instance
lock. Task claim and completion paths use conditional ownership and live-lease checks.
Expired work is recovered transactionally, and run/stage heartbeats make abandoned ingest
attempts visible on the next run.

Coordination tables are `findings`, `hypotheses`, `tasks`, and `coordinator_log`. The unified
execution ledger adds `runs`, `run_events`, and `run_stages`; `ops_schema_meta` records
migration state. Existing four-table databases migrate in place when opened.

---

## Schema

### findings

The core scientific output table. Every agent writes here as it works, not just when "done."

```sql
CREATE TABLE findings (
    id TEXT PRIMARY KEY,                -- UUID, generated server-side
    agent_id TEXT NOT NULL,             -- which agent produced this
    ts REAL NOT NULL,                   -- unix timestamp, set server-side
    finding_type TEXT NOT NULL,         -- caller-defined non-empty category
    domain TEXT NOT NULL,               -- caller-defined dataset/domain slug
    summary TEXT NOT NULL,              -- one-line human-readable
    evidence TEXT NOT NULL DEFAULT '{}', -- JSON, type-specific payload
    confidence REAL NOT NULL DEFAULT 0.5, -- 0.0-1.0
    novelty INTEGER NOT NULL DEFAULT 0,   -- 0=routine, 1=interesting, 2=novel, 3=potentially_significant
    parent_finding_id TEXT,             -- nullable, links follow-ups to parent
    reasoning TEXT NOT NULL DEFAULT ''  -- agent's interpretive logic (lab notebook entry)
);

-- Indices
CREATE INDEX idx_findings_ts ON findings(ts);
CREATE INDEX idx_findings_novelty ON findings(novelty);
CREATE INDEX idx_findings_agent ON findings(agent_id);
CREATE INDEX idx_findings_type ON findings(finding_type);
CREATE INDEX idx_findings_domain ON findings(domain);
```

`finding_type` and `domain` are intentionally open strings. Use stable project slugs;
the examples below are conventions, not exhaustive enums.

**Evidence payload by finding_type:**

```json
// gene
{"gene_id": "...", "genome_id": "...", "annotation": "...", "evalue": 1e-50,
 "coordinates": {"contig": "...", "start": 100, "end": 500, "strand": "+"}}

// neighborhood / cassette
{"genes": [{"gene_id": "...", "annotation": "...", "position": 1}, ...],
 "genome_ids": ["..."], "conservation_count": 44, "phyla_observed": ["...", "..."]}

// domain_architecture
{"protein_id": "...", "length_aa": 85804,
 "domains": [{"name": "...", "start": 0, "end": 500}, ...],
 "notable_features": ["solenoid", "dockerin-cohesin"]}

// phylogenetic
{"clade": "...", "anomaly_type": "long_branch|unexpected_topology|horizontal_transfer",
 "support_value": 0.98, "taxa_involved": ["...", "..."]}

// observation
{"description": "...", "supporting_finding_ids": ["...", "..."],
 "scope": "single_genome|clade|phylum|cross_phylum"}
```

### hypotheses

Cross-agent coordination. Testable claims that need validation by other agents.

```sql
CREATE TABLE hypotheses (
    id TEXT PRIMARY KEY,
    source_agent_id TEXT NOT NULL,
    source_finding_ids TEXT NOT NULL DEFAULT '[]',  -- JSON array of finding IDs
    ts REAL NOT NULL,
    hypothesis TEXT NOT NULL,                        -- the testable claim
    status TEXT NOT NULL DEFAULT 'proposed',          -- proposed|investigating|supported|refuted|inconclusive
    assigned_agent_id TEXT,                           -- who's working on it
    domains_relevant TEXT NOT NULL DEFAULT '[]',      -- JSON array of domain names to check
    evidence_for TEXT NOT NULL DEFAULT '[]',          -- JSON array of finding IDs
    evidence_against TEXT NOT NULL DEFAULT '[]',      -- JSON array of finding IDs
    resolution_summary TEXT                           -- written when status is terminal
);

CREATE INDEX idx_hyp_status ON hypotheses(status);
```

### tasks

Work queue. Coordinator writes tasks, agents claim and complete them.

```sql
CREATE TABLE tasks (
    id TEXT PRIMARY KEY,
    created_by TEXT NOT NULL,          -- usually 'coordinator'
    assigned_to TEXT,                  -- nullable = unassigned
    ts REAL NOT NULL,
    claimed_ts REAL,
    completed_ts REAL,
    status TEXT NOT NULL DEFAULT 'pending',  -- pending|claimed|in_progress|retry_wait|complete|failed
    priority INTEGER NOT NULL DEFAULT 1,     -- 0=low, 1=normal, 2=high, 3=urgent
    task_type TEXT NOT NULL,                  -- investigate|validate_hypothesis|cross_domain_search|follow_up|survey
    description TEXT NOT NULL,
    params TEXT NOT NULL DEFAULT '{}',        -- JSON, task-specific parameters
    domain_hint TEXT,                         -- suggested domain DB
    result_finding_ids TEXT NOT NULL DEFAULT '[]', -- populated on completion
    run_id TEXT,                            -- optional parent run
    idempotency_key TEXT,                   -- unique caller-stable creation key
    dependency_ids TEXT NOT NULL DEFAULT '[]', -- JSON task IDs; all must complete
    attempt_count INTEGER NOT NULL DEFAULT 0,
    max_attempts INTEGER NOT NULL DEFAULT 3,
    lease_seconds INTEGER NOT NULL DEFAULT 900,
    lease_expires_ts REAL,
    heartbeat_ts REAL,
    retry_after_ts REAL,
    last_error TEXT
);

CREATE INDEX idx_tasks_status ON tasks(status);
CREATE INDEX idx_tasks_assigned ON tasks(assigned_to);
CREATE INDEX idx_tasks_priority ON tasks(priority);
CREATE INDEX idx_tasks_run ON tasks(run_id);
CREATE INDEX idx_tasks_lease ON tasks(status, lease_expires_ts);
CREATE UNIQUE INDEX idx_tasks_idempotency
    ON tasks(idempotency_key) WHERE idempotency_key IS NOT NULL;
```

Only dependency-satisfied tasks are claimable. Claiming atomically increments
`attempt_count` and starts a finite lease. The owning worker must heartbeat before expiry;
completion/failure is rejected after lease loss. Retryable failures enter `retry_wait`, then
return to `pending` after `retry_after_ts` until `max_attempts` is exhausted.

### coordinator_log

The coordinator's own reasoning trail. Survives context compaction — the coordinator can query this to recover its train of thought.

```sql
CREATE TABLE coordinator_log (
    id TEXT PRIMARY KEY,
    ts REAL NOT NULL,
    action_type TEXT NOT NULL,  -- synthesis|task_assignment|escalation|hypothesis_generated|checkpoint
    reasoning TEXT NOT NULL,
    referenced_finding_ids TEXT NOT NULL DEFAULT '[]',
    referenced_hypothesis_ids TEXT NOT NULL DEFAULT '[]',
    decisions_made TEXT NOT NULL DEFAULT '{}'  -- JSON, structured record of what was decided
);

CREATE INDEX idx_coordlog_ts ON coordinator_log(ts);
```

### runs, run_events, and run_stages

`runs` is one durable execution record keyed by dataset path and run type. Status is
`pending`, `submitted`, `running`, `complete`, or `failed`; `submitted` represents a scheduler
bundle that has been handed off but has no stage currently executing, including time between
dependent jobs, so queue latency is not mistaken for a stale worker. The row records
idempotency key, lifecycle timestamps, heartbeat, parent run, current stage, config, result,
and terminal error. `run_events` is the append-only event stream. `run_stages` records each
stage attempt with its signature, command, input/output snapshots, resource profile,
heartbeats, terminal error, and reuse origin.

For ingest resume, a stage is reusable only when dataset, run type, stage ID, and signature
match a prior successful/reused attempt and its required outputs match that attempt's recorded
snapshots. A reused attempt is itself recorded (`attempt=0`) with the source run and attempt.
Reuse is dependency-conservative: when any upstream stage executes, all transitive downstream
stages execute rather than reusing artifacts derived from the previous upstream attempt.

---

## API Endpoints

### Findings

| Endpoint | Action | Contract |
|--------|------|-------------|
| `POST /findings` | Create a finding | Body: `FindingIn` (see Pydantic models). Returns `{id, ts}`. |
| `GET /findings` | List findings | Query params: `since` (float), `min_novelty` (int), `finding_type`, `domain`, `agent_id`, `limit` (default 50, max 500). Returns array. |
| `GET /findings/{id}` | Get one finding | Returns full finding dict. |
| `GET /findings/search/{text}` | Full-text search | LIKE search across summary, reasoning, evidence. Query param: `limit`. |

### Hypotheses

| Endpoint | Action | Contract |
|--------|------|-------------|
| `POST /hypotheses` | Propose a hypothesis | Body: `HypothesisIn`. Returns `{id, ts}`. |
| `GET /hypotheses` | List hypotheses | Query params: `status`, `unassigned` (bool), `limit`. |
| `PATCH /hypotheses/{id}` | Update hypothesis | Body: `HypothesisUpdate`. Used to assign, add evidence, resolve. |

### Tasks

| Endpoint | Action | Contract |
|--------|------|-------------|
| `POST /tasks` | Create a task | Body: `TaskIn`. Returns `{id, ts}`. |
| `GET /tasks` | List tasks | Query params: `status`, `assigned_to`, `unassigned` (bool), `limit`. |
| `PATCH /tasks/{id}` | Update task status | Body: `TaskUpdate`. |
| `POST /tasks/{id}/claim?agent_id=X` | Atomically claim a task | Requires satisfied dependencies and remaining attempts; starts a lease. |
| `POST /tasks/{id}/heartbeat?agent_id=X` | Heartbeat owned task | Extends a still-live lease and marks it `in_progress`. |
| `POST /tasks/recover` | Recover expired tasks | Requeues attempts or marks exhausted tasks failed. |

### Runs

| Endpoint | Action | Contract |
|--------|------|-------------|
| `POST /runs` | Create/deduplicate a run | Idempotency keys reject conflicting payloads. |
| `GET /runs` | List runs | Filter by dataset path, run type, or status. |
| `GET /runs/{id}` | Get a run | Returns decoded config/result. |
| `POST /runs/{id}/start` | Start a run | Sets lifecycle timestamps and heartbeat. |
| `POST /runs/{id}/submit` | Mark scheduler handoff | Atomically moves pending → submitted. |
| `POST /runs/{id}/heartbeat` | Heartbeat a run | Updates a running record. |
| `PATCH /runs/{id}` | Complete or fail a run | Writes result or terminal error. |
| `GET /runs/{id}/events` | Read event stream | Ordered run/stage lifecycle events. |
| `GET /runs/{id}/stages` | Inspect stage attempts | Optional `stage_id` filter; decoded snapshots/resources. |
| `POST /runs/recover` | Recover stale runs | Fails expired running stage attempts and parent runs. |

### Coordinator Log

| Endpoint | Action | Contract |
|--------|------|-------------|
| `POST /coordinator/log` | Write a log entry | Body: `CoordinatorLogIn`. |
| `GET /coordinator/log` | Read recent log | Query params: `limit`, `since`. |

### Utility

| Endpoint | Action | Contract |
|--------|------|-------------|
| `GET /stream` | SSE event stream | Real-time push notifications for new findings, hypotheses, tasks. Send keepalive every 30s. |
| `GET /stats` | Dashboard stats | Counts per table, findings by novelty, tasks by status, hypotheses by status. |
| `GET /health` | Health check | Returns `{status, db, ts}`. |

---

## Pydantic Models

```python
class FindingIn(BaseModel):
    agent_id: str
    finding_type: str   # caller-defined non-empty category
    domain: str         # caller-defined dataset/domain slug
    summary: str
    evidence: dict = {}
    confidence: float = 0.5  # 0.0-1.0
    novelty: int = 0         # 0-3
    parent_finding_id: Optional[str] = None
    reasoning: str = ""

class HypothesisIn(BaseModel):
    source_agent_id: str
    source_finding_ids: list[str] = []
    hypothesis: str
    domains_relevant: list[str] = []

class HypothesisUpdate(BaseModel):
    status: str  # proposed|investigating|supported|refuted|inconclusive
    assigned_agent_id: Optional[str] = None
    evidence_for: Optional[list[str]] = None
    evidence_against: Optional[list[str]] = None
    resolution_summary: Optional[str] = None

class TaskIn(BaseModel):
    created_by: str = "coordinator"
    task_type: str  # investigate|validate_hypothesis|cross_domain_search|follow_up|survey
    description: str
    params: dict = {}
    priority: int = 1  # 0-3
    domain_hint: Optional[str] = None
    assigned_to: Optional[str] = None
    run_id: Optional[str] = None
    idempotency_key: Optional[str] = None
    depends_on: list[str] = []
    max_attempts: int = 3
    lease_seconds: int = 900

class TaskUpdate(BaseModel):
    status: str  # pending|claimed|in_progress|retry_wait|complete|failed
    agent_id: Optional[str] = None
    result_finding_ids: Optional[list[str]] = None
    error: Optional[str] = None
    retryable: bool = False
    retry_delay_seconds: int = 0

class CoordinatorLogIn(BaseModel):
    action_type: str  # synthesis|task_assignment|escalation|hypothesis_generated|checkpoint
    reasoning: str
    referenced_finding_ids: list[str] = []
    referenced_hypothesis_ids: list[str] = []
    decisions_made: dict = {}
```

---

## Client Library

`sharur/ops/client.py` provides a `SharurOps` class that wraps all endpoints. Usage:

```python
from sharur.ops.client import SharurOps

ops = SharurOps("http://localhost:8811", agent_id="dpann_surveyor")

# Post a finding
fid = ops.finding(
    finding_type="cassette",
    domain="dpann",
    summary="Novel CHAT protease defense cassette conserved across 7 DPANN phyla",
    evidence={
        "genes": [
            {"gene_id": "g001", "annotation": "GHG protease", "position": 1},
            {"gene_id": "g002", "annotation": "ELV-DEIG domain", "position": 2},
            {"gene_id": "g003", "annotation": "PTW protein", "position": 3},
        ],
        "conservation_count": 44,
        "phyla_observed": ["Nanoarchaeota", "Woesearchaeota", "Pacearchaeota",
                           "Aenigmarchaeota", "Diapherotrites", "Micrarchaeota", "Huberarchaeota"],
    },
    confidence=0.85,
    novelty=3,
    reasoning="Detected via ELSA syntenic neighborhood search at tau=0.85. "
              "GHG protease is uncharacterized but structurally similar to CRISPR-associated "
              "Lon proteases. Co-occurrence with ELV-DEIG domain protein suggests defense function. "
              "Conservation across 7 phyla indicates ancient origin pre-DPANN radiation."
)

# Propose a hypothesis for cross-domain validation
hid = ops.hypothesis(
    hypothesis="CHAT protease cassette is present in Bathyarchaeia based on shared ancestry",
    source_finding_ids=[fid],
    domains_relevant=["bathyarchaeia", "omnitrophota"],
)

# Check for assigned tasks
for task in ops.my_tasks():
    print(task["description"])

# Claim an unassigned task
available = ops.available_tasks()
if available:
    ops.claim_task(available[0]["id"])

# Check what other agents found recently
interesting = ops.recent_findings(min_novelty=2)

# Coordinator: recover reasoning after compaction
log = ops.recent_log(limit=10)

# Coordinator: log a synthesis decision
ops.log(
    action_type="synthesis",
    reasoning="Three agents independently found ELV-DEIG domains in defense contexts. "
              "Upgrading to cross-domain investigation priority.",
    referenced_finding_ids=[fid, "...", "..."],
    decisions_made={
        "action": "created cross_domain_search task",
        "rationale": "convergent findings from independent agents strongly suggest conserved system"
    },
)

# Dashboard
print(ops.stats())
```

### Client Methods Reference

| Method | Description |
|--------|-------------|
| `finding(...)` | Post a finding. Returns finding ID. |
| `recent_findings(since, min_novelty, finding_type, domain, limit)` | Query findings. |
| `get_finding(id)` | Get a single finding by ID. |
| `search_findings(text, limit)` | Full-text search across findings. |
| `hypothesis(...)` | Propose a hypothesis. Returns hypothesis ID. |
| `open_hypotheses(unassigned)` | Get hypotheses needing investigation. |
| `update_hypothesis(id, **kwargs)` | Update hypothesis status/evidence. |
| `create_task(...)` | Create a task. Returns task ID. |
| `my_tasks(status)` | Get tasks assigned to this agent. |
| `available_tasks()` | Get unassigned pending tasks. |
| `claim_task(id)` | Atomically claim a task (409 if already claimed). |
| `heartbeat_task(id, lease_seconds)` | Extend the current agent's live task lease. |
| `complete_task(id, result_finding_ids)` | Mark task complete. |
| `fail_task(id, error, retryable, retry_delay_seconds)` | Fail or schedule retry. |
| `recover_expired_tasks()` | Requeue expired attempts or fail exhausted tasks. |
| `create_run(...)`, `start_run(id)`, `submit_run(id)`, `heartbeat_run(id)` | Manage durable runs. |
| `complete_run(id, result)`, `fail_run(id, error)` | Terminate runs. |
| `list_runs(...)`, `get_run(id)`, `run_events(id)`, `run_stages(id)` | Inspect run history. |
| `recover_stale_runs(...)` | Mark abandoned run/stage attempts failed. |
| `log(...)` | Write a coordinator log entry. |
| `recent_log(limit, since)` | Read coordinator log. |
| `stats()` | Get dashboard stats. |

---

## Key Design Decisions

1. **SQLite WAL mode, not DuckDB.** DuckDB is single-writer and puts locks on the file. SQLite WAL gives concurrent reads with serialized writes. Write lock duration for one row insert is microseconds — 10 agents will almost never actually contend.

2. **HTTP, not shared files.** Works identically on laptop (localhost) and cluster (remote node). No shared filesystem, no NFS, no NDJSON append atomicity concerns.

3. **Atomic leased task claiming.** `POST /tasks/{id}/claim` uses a write transaction plus
   a conditional update. If two agents race, exactly one claims the task. Leases, heartbeats,
   bounded attempts, dependency gates, delayed retries, and idempotency keys handle worker
   death and duplicate dispatch.

4. **SSE for real-time events.** The `/stream` endpoint pushes finding/hypothesis/task creation events. The coordinator can subscribe instead of polling. Keepalive every 30s prevents timeouts.

5. **Mixed-granularity findings.** An open `finding_type` string plus flexible JSON
   `evidence` handles everything from single genes to cross-dataset observations without
   schema sprawl.

6. **Coordinator log as compaction insurance.** The coordinator writes its reasoning to `coordinator_log` before and after synthesis steps. After compaction, it queries `GET /coordinator/log?limit=10` to recover its train of thought.

7. **Parent finding chains.** `parent_finding_id` links follow-up investigations to what triggered them. Gives you a provenance graph of how discoveries were made.

8. **Flexible agent assignment.** No hardcoded domain-to-agent mapping. Any agent can work any domain. The `domain_hint` on tasks is a suggestion, not a constraint.

---

## Coordinator Loop

The intended orchestrator workflow after integration:

```
1. Query GET /findings?min_novelty=2&since={last_check_ts}
   → "What's new and interesting?"

2. Query GET /hypotheses?status=proposed&unassigned=true
   → "Any hypotheses need assignment?"

3. Synthesize across the structured results
   → This is where model capability matters — reasoning over compact structured records

4. POST /coordinator/log with synthesis reasoning
   → Externalize the train of thought BEFORE it can be compacted away

5. POST /tasks to assign follow-up investigations
   → Or POST /hypotheses to propose new cross-domain hypotheses

6. Repeat. On compaction recovery: GET /coordinator/log?limit=10
```

---

## Durability and Recovery Boundary

SQLite is the executor-of-record for Sharur task and ingest lifecycle state on a shared
filesystem. It provides transactional claims, finite leases, heartbeat recovery, bounded
retries, dependency gates, idempotent creation, and append-only run events. It does not
replace a remote scheduler: SLURM remains responsible for allocating/killing cluster jobs,
and network partitions can delay heartbeats until a lease expires. Workers must treat a lost
lease as lost authority and must not publish terminal task state afterward.

Ingest subprocesses heartbeat every 30 seconds. On the next ingest invocation, stale running
stage attempts and their parent runs are marked failed before new work starts. Resume is
signature- and output-gated; it never infers success solely from a directory's existence.
Run/stage starts are atomic as well as run creation: a duplicate launcher using the same
idempotency key cannot start the same run or stage twice. SLURM handoff uses `submitted`
status so long scheduler queue time is not crash-recovered as active worker failure.

## File Manifest

```
sharur/ops/server.py  — FastAPI server
sharur/ops/client.py  — HTTP client
sharur/ops/store.py   — direct SQLite client with matching semantics
sharur/ops/schema.py  — migration-safe shared schema
sharur/ops/ledger.py  — run, event, stage-attempt, reuse, and recovery API
agent_ops_spec.md     — schema and operating contract
```
