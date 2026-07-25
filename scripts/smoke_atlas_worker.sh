#!/usr/bin/env bash
# End-to-end smoke test for the Atlas scan worker executor.
#
# Builds a miniature sealed dataset, plans it, runs the zero-model-call packet
# census, stands up Ops + Query, enqueues the plan, and drives one worker in
# --dry-run mode through claim -> packet walk -> coverage manifest -> candidates
# -> disposition -> complete.
#
# Makes ZERO model calls and touches no production data. Everything lives in a
# temp dir that is removed on exit.
#
#   bash scripts/smoke_atlas_worker.sh

set -euo pipefail

OPS_PORT="${OPS_PORT:-8871}"
QUERY_PORT="${QUERY_PORT:-8872}"
WORK="$(mktemp -d -t sharur-smoke-XXXXXX)"
PIDS=()

cleanup() {
    for pid in "${PIDS[@]:-}"; do kill "$pid" 2>/dev/null || true; done
    sleep 1
    rm -rf "$WORK"
}
trap cleanup EXIT

say() { printf '\n\033[1m== %s\033[0m\n' "$*"; }

say "1. build a miniature sealed dataset"
python - "$WORK" <<'PY'
import sys
from pathlib import Path
from sharur.dataset_seal import build_dataset_seal, write_dataset_seal
from sharur.storage.duckdb_store import DuckDBStore

root = Path(sys.argv[1])
db = root / "sharur.duckdb"
store = DuckDBStore(db)
store.execute(
    "INSERT INTO bins(bin_id, completeness, contamination, taxonomy, n_contigs, total_length)"
    " VALUES ('genome_a',90,2,'d__Bacteria;p__Chloroflexota',2,3000),"
    "        ('genome_b',80,3,'d__Bacteria;p__Chloroflexota',1,1500)"
)
store.execute(
    "INSERT INTO contigs(contig_id, bin_id, length) VALUES"
    " ('a_1','genome_a',1000),('a_2','genome_a',2000),('b_1','genome_b',1500)"
)
rows = []
for contig, binid, n in (("a_1", "genome_a", 8), ("a_2", "genome_a", 12), ("b_1", "genome_b", 6)):
    for i in range(1, n + 1):
        rows.append(
            f"('{contig}_p{i}','{contig}','{binid}',{i*100},{i*100+299},'+',{i},100)"
        )
store.execute(
    "INSERT INTO proteins(protein_id, contig_id, bin_id, start, end_coord, strand,"
    " gene_index, sequence_length) VALUES " + ",".join(rows)
)
store.close()
write_dataset_seal(build_dataset_seal(db, max_hash_bytes=0), root / "dataset.seal.json")
print(f"   dataset built: {db}")
PY

say "2. plan"
sharur-atlas plan --db "$WORK/sharur.duckdb" --output-dir "$WORK/atlas" \
    --packet-bytes 8192 --checkpoint-interval-frames 1 --threads 2 \
    | tail -3

say "3. packet census (zero model calls)"
sharur-atlas packet-census --plan-dir "$WORK/atlas" --workers 2 --threads 2 \
    --temp-directory "$WORK/atlas/spill" | tail -5
sharur-atlas verify-packet-census --plan-dir "$WORK/atlas" --deep | tail -3

say "4. start Ops + Query (loopback, tokenless)"
sharur-ops --db "$WORK/sharur_ops.db" --host 127.0.0.1 --port "$OPS_PORT" \
    > "$WORK/ops.log" 2>&1 &
PIDS+=($!)
sharur-query --db "$WORK/sharur.duckdb" --direct --host 127.0.0.1 --port "$QUERY_PORT" \
    --threads 2 --memory-limit 2GB > "$WORK/query.log" 2>&1 &
PIDS+=($!)

for _ in $(seq 1 30); do
    if curl -sf "http://127.0.0.1:$OPS_PORT/health" >/dev/null 2>&1 \
       && curl -sf "http://127.0.0.1:$QUERY_PORT/health" >/dev/null 2>&1; then
        break
    fi
    sleep 1
done
curl -sf "http://127.0.0.1:$OPS_PORT/health"   >/dev/null || { echo "ops failed"; tail -20 "$WORK/ops.log"; exit 1; }
curl -sf "http://127.0.0.1:$QUERY_PORT/health" >/dev/null || { echo "query failed"; tail -20 "$WORK/query.log"; exit 1; }
echo "   both services healthy"

say "5. enqueue the sealed plan"
sharur-atlas enqueue --plan-dir "$WORK/atlas" \
    --ops-url "http://127.0.0.1:$OPS_PORT" \
    --query-url "http://127.0.0.1:$QUERY_PORT" \
    --scan-execution-profile atlas_scan | tail -5

say "6. run the worker (--dry-run: zero model calls)"
sharur-worker atlas-scan \
    --ops-url "http://127.0.0.1:$OPS_PORT" \
    --query-url "http://127.0.0.1:$QUERY_PORT" \
    --agent-id smoke-worker-01 \
    --profile atlas_scan \
    --dry-run --max-tasks 2 --idle-sleep 1

say "7. verify final state"
python - "$OPS_PORT" <<'PY'
import sys, requests
base = f"http://127.0.0.1:{sys.argv[1]}"
tasks = requests.get(f"{base}/tasks", params={"limit": 50}, timeout=10).json()
rows = tasks if isinstance(tasks, list) else tasks.get("items", tasks.get("tasks", []))
by_status = {}
for t in rows:
    by_status.setdefault(t["status"], 0)
    by_status[t["status"]] += 1
print(f"   task statuses: {by_status}")
assert by_status.get("complete") == 2, f"expected 2 complete tasks, got {by_status}"
print("   PASS: both genomes claimed, covered, dispositioned, and completed")
PY

say "SMOKE TEST PASSED"
