# Sharur Ops Layer — Schema Design

## Access Modes

Two backends expose the same API and schema:

```python
import os

# Direct SQLite (no server needed — use for local/single-machine runs)
from sharur.ops.store import OpsStore
ops = OpsStore("data/DATASET/sharur_ops.db", agent_id="my_agent")

# HTTP client (requires `pip install -e ".[ops]"` and `sharur-ops`)
from sharur.ops.client import SharurOps
ops = SharurOps(
    "http://localhost:8811",
    agent_id="my_agent",
    api_token=os.environ.get("SHARUR_OPS_TOKEN"),
)
```

Both support findings/hypotheses, leased tasks, idempotent runs, recovery, coordinator logs,
and stats. Point the HTTP server and direct clients at the same dataset-local SQLite file.
The server opens no database at import time, binds to loopback by default, requires bearer
authentication for remote use, and uses short-lived request connections. HTTP clients have
bounded timeouts and retry only idempotent reads.

## Design Principles

1. **Mixed granularity**: Findings range from single genes to multi-phylum observations.
   Common columns for everything queryable; a flexible JSON `evidence` payload for type-specific data.
2. **Flexible agents**: No hardcoded domain assignment. Agents declare capabilities;
   orchestrator routes by capability + availability.
3. **Append-heavy, update-light**: Findings are immutable once written. Hypotheses and tasks
   have status transitions but no payload edits. This keeps the write path simple.
4. **Orchestrator queries, doesn't subscribe to firehoses**: Every table has timestamp +
   novelty/priority flags so the orchestrator pulls what it needs, when it needs it.
5. **Finite ownership**: Task authority is a renewable lease. A dead worker cannot strand
   work forever, and an expired worker cannot complete a task it no longer owns.
6. **Replayable execution**: Runs have append-only events; every stage attempt records its
   signature, command, resources, inputs, outputs, heartbeat, and reuse origin.

---

## Tables

### `findings`

The core scientific output. Every agent writes here.

| Column | Type | Description |
|---|---|---|
| `id` | TEXT (uuid) | Unique finding ID, generated on insert |
| `agent_id` | TEXT | Which agent produced this |
| `ts` | REAL | Unix timestamp |
| `finding_type` | TEXT | Caller-defined non-empty category; use stable project conventions |
| `domain` | TEXT | Caller-defined dataset/domain slug |
| `summary` | TEXT | One-line human-readable description |
| `evidence` | TEXT (JSON) | Type-specific structured payload (see below) |
| `confidence` | REAL | 0.0–1.0, agent's self-assessment |
| `novelty` | INTEGER | 0=routine, 1=interesting, 2=novel, 3=potentially_significant |
| `parent_finding_id` | TEXT | Nullable. Links follow-up findings to their parent |
| `reasoning` | TEXT | Agent's interpretive logic (the lab notebook entry) |

**Evidence payload examples by finding_type:**

```json
// gene
{"gene_id": "...", "genome_id": "...", "annotation": "...", "evalue": 1e-50, 
 "coordinates": {"contig": "...", "start": 100, "end": 500, "strand": "+"}}

// neighborhood / cassette
{"genes": [{"gene_id": "...", "annotation": "...", "position": 1}, ...],
 "genome_ids": ["..."], "conservation_count": 44, "phyla_observed": ["...", "..."]}

// domain_architecture  
{"protein_id": "...", "length_aa": 85804, "domains": [{"name": "...", "start": 0, "end": 500}, ...],
 "notable_features": ["solenoid", "dockerin-cohesin"]}

// phylogenetic
{"clade": "...", "anomaly_type": "long_branch|unexpected_topology|horizontal_transfer",
 "support_value": 0.98, "taxa_involved": ["...", "..."]}

// observation
{"description": "...", "supporting_finding_ids": ["...", "..."],
 "scope": "single_genome|clade|phylum|cross_phylum"}
```

### `hypotheses`

Cross-agent coordination. Testable claims that need validation.

| Column | Type | Description |
|---|---|---|
| `id` | TEXT (uuid) | |
| `source_agent_id` | TEXT | Who proposed it |
| `source_finding_ids` | TEXT (JSON array) | Findings that motivated this hypothesis |
| `ts` | REAL | |
| `hypothesis` | TEXT | The testable claim |
| `status` | TEXT | `proposed`, `investigating`, `supported`, `refuted`, `inconclusive` |
| `assigned_agent_id` | TEXT | Nullable — who's working on it |
| `domains_relevant` | TEXT (JSON array) | Which lineage DBs should be checked |
| `evidence_for` | TEXT (JSON array) | Finding IDs supporting |
| `evidence_against` | TEXT (JSON array) | Finding IDs contradicting |
| `resolution_summary` | TEXT | Nullable — written when status is terminal |

### `tasks`

Orchestrator → agent work queue.

| Column | Type | Description |
|---|---|---|
| `id` | TEXT (uuid) | |
| `created_by` | TEXT | Usually `coordinator`, but agents can spawn sub-tasks |
| `assigned_to` | TEXT | Nullable = unassigned, agent claims by writing its ID |
| `ts` | REAL | Created timestamp |
| `claimed_ts` | REAL | Nullable — when an agent claimed it |
| `completed_ts` | REAL | Nullable |
| `status` | TEXT | `pending`, `claimed`, `in_progress`, `retry_wait`, `complete`, `failed` |
| `priority` | INTEGER | 0=low, 1=normal, 2=high, 3=urgent |
| `task_type` | TEXT | `investigate`, `validate_hypothesis`, `cross_domain_search`, `follow_up`, `survey` |
| `description` | TEXT | What to do |
| `params` | TEXT (JSON) | Task-specific parameters (query templates, gene IDs, etc.) |
| `domain_hint` | TEXT | Nullable — suggested domain DB |
| `result_finding_ids` | TEXT (JSON array) | Populated on completion |
| `run_id` | TEXT | Optional parent run |
| `idempotency_key` | TEXT | Unique caller-stable creation key |
| `dependency_ids` | TEXT (JSON array) | Tasks that must all complete before claim |
| `attempt_count`, `max_attempts` | INTEGER | Bounded attempt accounting |
| `lease_seconds` | INTEGER | Lease duration requested at claim |
| `lease_expires_ts`, `heartbeat_ts` | REAL | Live ownership boundary |
| `retry_after_ts` | REAL | Delayed retry eligibility |
| `last_error` | TEXT | Most recent worker/recovery error |

### `coordinator_log`

Orchestrator's own reasoning trail. Survives context compaction.

| Column | Type | Description |
|---|---|---|
| `id` | TEXT (uuid) | |
| `ts` | REAL | |
| `action_type` | TEXT | `synthesis`, `task_assignment`, `escalation`, `hypothesis_generated`, `checkpoint` |
| `reasoning` | TEXT | The orchestrator's thought process |
| `referenced_finding_ids` | TEXT (JSON array) | |
| `referenced_hypothesis_ids` | TEXT (JSON array) | |
| `decisions_made` | TEXT (JSON) | Structured record of what was decided and why |

### `runs`

One record per durable execution. The row includes `id`, optional unique
`idempotency_key`, `run_type`, normalized `dataset_path`, creator, status
(`pending|submitted|running|complete|failed`), lifecycle timestamps/heartbeat,
parent/current stage, canonical JSON config/result, and terminal error. `submitted` is a
scheduler handoff with no stage currently executing (before or between jobs) and is not
subject to running-heartbeat recovery.

### `run_events`

Append-only ordered lifecycle events: run ID, timestamp, event type, optional stage/attempt,
and JSON payload.

### `run_stages`

One row per `(run_id, stage_id, attempt)`. It records the deterministic stage signature,
status, timestamps/heartbeat, command, input/output snapshots, resource profile, error, and
optional source run/attempt for reuse. Reused stages use attempt zero in the new run.

### `ops_schema_meta`

Records the applied ops schema version. Opening a legacy four-table database adds task
lease/retry/dependency/idempotency columns and the run-ledger tables in place.

---

## Key Query Patterns

```sql
-- Orchestrator: "What's new and interesting?"
SELECT * FROM findings WHERE novelty >= 2 AND ts > ? ORDER BY novelty DESC, ts DESC;

-- Orchestrator: "Any hypotheses need assignment?"
SELECT * FROM hypotheses WHERE status = 'proposed' AND assigned_agent_id IS NULL;

-- Agent: "Do I have currently owned work?"
SELECT * FROM tasks
WHERE assigned_to = ?
  AND status IN ('claimed', 'in_progress')
  AND lease_expires_ts > unixepoch('subsec');

-- Orchestrator: "Cross-domain patterns?"
SELECT f1.summary, f2.summary, f1.domain, f2.domain
FROM findings f1 JOIN findings f2 
ON json_extract(f1.evidence, '$.annotation') = json_extract(f2.evidence, '$.annotation')
WHERE f1.domain != f2.domain AND f1.novelty >= 1 AND f2.novelty >= 1;

-- Agent: "Has anyone else seen this gene neighborhood?"
SELECT * FROM findings WHERE finding_type = 'neighborhood' 
AND evidence LIKE '%"GHG protease"%';

-- Orchestrator: "Recover my train of thought"
SELECT * FROM coordinator_log ORDER BY ts DESC LIMIT 10;

-- Operator: "What happened in this run?"
SELECT ts, event_type, stage_id, attempt, payload
FROM run_events WHERE run_id = ? ORDER BY id;
```

Do not claim tasks with ad-hoc SQL. Use `OpsStore.claim_task()` or
`POST /tasks/{id}/claim`; those paths transactionally recover expirations, promote eligible
retries, enforce dependencies, and increment attempts.
