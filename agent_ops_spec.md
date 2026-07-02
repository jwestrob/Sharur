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
│  Agent 2    │─────────────────────┼────▶│  sharur_ops.py   │────▶│ sharur_ops.db  │
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

Domain databases (DuckDB files) remain read-only reference data. Agents ATTACH them directly for analytical queries. The ops server handles only coordination state.

## Dependencies

```
pip install fastapi uvicorn pydantic
```

No other dependencies. Python 3.10+.

## Running

```bash
# Start the server (creates sharur_ops.db on first run)
uvicorn sharur_ops:app --host 0.0.0.0 --port 8811

# Or directly
python sharur_ops.py
```

Works identically on laptop (localhost:8811) and cluster (node:8811). Agents just change the base URL.

The coordination DB defaults to `sharur_ops.db` in the process working directory. Set
`SHARUR_OPS_DB_PATH` to pin it somewhere stable — a persistent volume, or a fixed path so the
coordination state survives starting the server from a different directory:

```bash
export SHARUR_OPS_DB_PATH=/data/ops/sharur_ops.db
uvicorn sharur.ops.server:app --host 0.0.0.0 --port 8811
```

## Database

SQLite in WAL mode. Single connection, writes serialized via threading lock. This eliminates all concurrency issues — WAL gives concurrent reads, the lock serializes writes, and each write (one finding row insert) takes microseconds.

Four tables: `findings`, `hypotheses`, `tasks`, `coordinator_log`.

---

## Schema

### findings

The core scientific output table. Every agent writes here as it works, not just when "done."

```sql
CREATE TABLE findings (
    id TEXT PRIMARY KEY,                -- UUID, generated server-side
    agent_id TEXT NOT NULL,             -- which agent produced this
    ts REAL NOT NULL,                   -- unix timestamp, set server-side
    finding_type TEXT NOT NULL,         -- see enum below
    domain TEXT NOT NULL,               -- see enum below
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

**finding_type enum:** `gene`, `neighborhood`, `cassette`, `domain_architecture`, `phylogenetic`, `observation`

**domain enum:** `omnitrophota`, `dpann`, `bathyarchaeia`, `hinthialibacterota`, `cross_domain`

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
    status TEXT NOT NULL DEFAULT 'pending',  -- pending|claimed|in_progress|complete|failed
    priority INTEGER NOT NULL DEFAULT 1,     -- 0=low, 1=normal, 2=high, 3=urgent
    task_type TEXT NOT NULL,                  -- investigate|validate_hypothesis|cross_domain_search|follow_up|survey
    description TEXT NOT NULL,
    params TEXT NOT NULL DEFAULT '{}',        -- JSON, task-specific parameters
    domain_hint TEXT,                         -- suggested domain DB
    result_finding_ids TEXT NOT NULL DEFAULT '[]'  -- populated on completion
);

CREATE INDEX idx_tasks_status ON tasks(status);
CREATE INDEX idx_tasks_assigned ON tasks(assigned_to);
CREATE INDEX idx_tasks_priority ON tasks(priority);
```

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

---

## API Endpoints

### Findings

| Method | Path | Description |
|--------|------|-------------|
| `POST /findings` | Create a finding | Body: `FindingIn` (see Pydantic models). Returns `{id, ts}`. |
| `GET /findings` | List findings | Query params: `since` (float), `min_novelty` (int), `finding_type`, `domain`, `agent_id`, `limit` (default 50, max 500). Returns array. |
| `GET /findings/{id}` | Get one finding | Returns full finding dict. |
| `GET /findings/search/{text}` | Full-text search | LIKE search across summary, reasoning, evidence. Query param: `limit`. |

### Hypotheses

| Method | Path | Description |
|--------|------|-------------|
| `POST /hypotheses` | Propose a hypothesis | Body: `HypothesisIn`. Returns `{id, ts}`. |
| `GET /hypotheses` | List hypotheses | Query params: `status`, `unassigned` (bool), `limit`. |
| `PATCH /hypotheses/{id}` | Update hypothesis | Body: `HypothesisUpdate`. Used to assign, add evidence, resolve. |

### Tasks

| Method | Path | Description |
|--------|------|-------------|
| `POST /tasks` | Create a task | Body: `TaskIn`. Returns `{id, ts}`. |
| `GET /tasks` | List tasks | Query params: `status`, `assigned_to`, `unassigned` (bool), `limit`. |
| `PATCH /tasks/{id}` | Update task status | Body: `TaskUpdate`. |
| `POST /tasks/{id}/claim?agent_id=X` | Atomically claim a task | Conditional UPDATE — only succeeds if task is still `pending`. Returns 409 if already claimed. |

### Coordinator Log

| Method | Path | Description |
|--------|------|-------------|
| `POST /coordinator/log` | Write a log entry | Body: `CoordinatorLogIn`. |
| `GET /coordinator/log` | Read recent log | Query params: `limit`, `since`. |

### Utility

| Method | Path | Description |
|--------|------|-------------|
| `GET /stream` | SSE event stream | Real-time push notifications for new findings, hypotheses, tasks. Send keepalive every 30s. |
| `GET /stats` | Dashboard stats | Counts per table, findings by novelty, tasks by status, hypotheses by status. |
| `GET /health` | Health check | Returns `{status, db, ts}`. |

---

## Pydantic Models

```python
class FindingIn(BaseModel):
    agent_id: str
    finding_type: str   # gene|neighborhood|cassette|domain_architecture|phylogenetic|observation
    domain: str         # omnitrophota|dpann|bathyarchaeia|hinthialibacterota|cross_domain
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

class TaskUpdate(BaseModel):
    status: str  # pending|claimed|in_progress|complete|failed
    assigned_to: Optional[str] = None
    result_finding_ids: Optional[list[str]] = None

class CoordinatorLogIn(BaseModel):
    action_type: str  # synthesis|task_assignment|escalation|hypothesis_generated|checkpoint
    reasoning: str
    referenced_finding_ids: list[str] = []
    referenced_hypothesis_ids: list[str] = []
    decisions_made: dict = {}
```

---

## Client Library

`sharur_client.py` provides a `SharurOps` class that wraps all endpoints. Usage:

```python
from sharur_client import SharurOps

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
| `complete_task(id, result_finding_ids)` | Mark task complete. |
| `fail_task(id)` | Mark task failed. |
| `log(...)` | Write a coordinator log entry. |
| `recent_log(limit, since)` | Read coordinator log. |
| `stats()` | Get dashboard stats. |

---

## Key Design Decisions

1. **SQLite WAL mode, not DuckDB.** DuckDB is single-writer and puts locks on the file. SQLite WAL gives concurrent reads with serialized writes. Write lock duration for one row insert is microseconds — 10 agents will almost never actually contend.

2. **HTTP, not shared files.** Works identically on laptop (localhost) and cluster (remote node). No shared filesystem, no NFS, no NDJSON append atomicity concerns.

3. **Atomic task claiming.** `POST /tasks/{id}/claim` uses a conditional UPDATE (`WHERE status = 'pending'`). If two agents race, exactly one gets 200, the other gets 409. No distributed locking needed.

4. **SSE for real-time events.** The `/stream` endpoint pushes finding/hypothesis/task creation events. The coordinator can subscribe instead of polling. Keepalive every 30s prevents timeouts.

5. **Mixed-granularity findings.** The `finding_type` enum + flexible JSON `evidence` payload handles everything from single genes to cross-phylum observations without schema sprawl.

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

## File Manifest

```
sharur_ops.py      — FastAPI server (~350 lines)
sharur_client.py   — Agent client library (~180 lines)
SCHEMA.md          — Schema design document with query patterns
```

All three files are provided in the accompanying tarball.
