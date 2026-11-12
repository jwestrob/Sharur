#!/usr/bin/env bash
# Submit the remote-homology searches (sharur/remote_homology/runner.py) as one Slurm array:
#   families  task i: the family-tier searches of every N-th family (largest families first)
#   cdd       task i: RPS-BLAST of every N-th chunk of unique sequences against CDD
#
# Usage:
#   scripts/submit_remote_homology_shards.sh STAGE -n N -w WORK_DIR -c CONFIG.json [-m MAX_RUNNING] \
#       [-e ENV_SCRIPT] [-s "SBATCH OPTIONS"] [-d DEPENDENCY] [-p PYTHON] [--dry-run] [-- RUNNER_ARGS...]
#
#   STAGE       families or cdd
#   WORK_DIR    the work directory from `sharur remote-homology export` / `export-cdd`
#   CONFIG      runner configuration (tools, databases with names and versions)
#   -s          site options for every task, e.g. "-A ACCOUNT -p PARTITION --qos QOS -c 24 --exclude NODES"
#   -m          array throttle: at most MAX_RUNNING tasks at once
#   -d          Slurm dependency, e.g. afterok:12345
#   -e          script sourced at the start of every task (environment setup)
#   RUNNER_ARGS extra runner.py arguments, e.g. --limit 30 for a pilot
#   Every task runs runner.py with --threads $SLURM_CPUS_ON_NODE (workers of 4 threads each).
#   Each family and chunk is skipped once finished, so resubmitting the same command continues
#   an interrupted run. The job script lands in WORK_DIR/slurm/ first; --dry-run stops there.
set -euo pipefail

usage() { sed -n '2,23p' "$0" >&2; exit 2; }
[[ $# -gt 0 ]] || usage
stage=$1; shift
[[ $stage == families || $stage == cdd ]] || usage
n="" work="" config="" max="" env_script="" opts="" dep="" dry=0 python=${PYTHON:-python3}
while [[ $# -gt 0 ]]; do
  case $1 in
    -n) n=$2; shift 2 ;;
    -w) work=$2; shift 2 ;;
    -c) config=$2; shift 2 ;;
    -m) max=$2; shift 2 ;;
    -e) env_script=$2; shift 2 ;;
    -s) opts=$2; shift 2 ;;
    -d) dep=$2; shift 2 ;;
    -p) python=$2; shift 2 ;;
    --dry-run) dry=1; shift ;;
    --) shift; break ;;
    *) usage ;;
  esac
done
[[ -n $n && -n $work && -n $config ]] || usage

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
runner=${RUNNER:-$here/../sharur/remote_homology/runner.py}
runner=$(cd "$(dirname "$runner")" && pwd)/$(basename "$runner")
[[ -f $runner ]] || { echo "runner.py not found at $runner (set RUNNER)" >&2; exit 2; }
work=$(cd "$work" && pwd)
config=$(cd "$(dirname "$config")" && pwd)/$(basename "$config")
mkdir -p "$work/slurm" "$work/logs"
args=""
if [[ $# -gt 0 ]]; then args=$(printf '%q ' "$@"); fi
array="0-$((n - 1))${max:+%$max}"
f=$work/slurm/rh_$stage.sbatch
{
  echo '#!/bin/bash'
  echo "#SBATCH -J rh_$stage"
  echo "#SBATCH --array=$array"
  echo "#SBATCH -o $work/logs/%x_%A_%a.out"
  line=""  # one #SBATCH line per option; a word without a leading dash is its value
  for o in $opts; do
    if [[ $o == -* && -n $line ]]; then echo "#SBATCH $line"; line=""; fi
    line="${line:+$line }$o"
  done
  if [[ -n $line ]]; then echo "#SBATCH $line"; fi
  echo 'set -euo pipefail'
  if [[ -n $env_script ]]; then echo "source $(printf %q "$env_script")"; fi
  echo 'echo "host $(hostname) cpus $SLURM_CPUS_ON_NODE task $SLURM_ARRAY_TASK_ID start $(date)"'
  echo "$python $(printf %q "$runner") $stage --work $(printf %q "$work") --config $(printf %q "$config") \\"
  echo "    --shard \$SLURM_ARRAY_TASK_ID --n-shards $n --threads \$SLURM_CPUS_ON_NODE $args"
  echo 'echo "done $(date)"'
} > "$f"
if (( dry )); then
  echo "wrote $f"
  exit 0
fi
job=$(sbatch --parsable ${dep:+--dependency=$dep} "$f"); job=${job%%;*}
printf '%s\t%s\t%s\n' "$stage" "$job" "$f" | tee -a "$work/jobs.tsv"
