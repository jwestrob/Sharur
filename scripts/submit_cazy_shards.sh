#!/usr/bin/env bash
# Submit a sharded dbCAN consensus run (scripts/classify_cazymes.py) as three chained Slurm jobs:
#   cazy_search  array task i: dbCAN.hmm profile part i on every protein; task 0 adds DIAMOND
#   cazy_sub     array task i: dbCAN-sub profile part i on the candidates   (afterok: cazy_search)
#   cazy_merge   consensus over every shard, written to --out               (afterok: cazy_sub)
#
# Usage:
#   scripts/submit_cazy_shards.sh -n N -w WORK_DIR [-m MAX_RUNNING] [-e ENV_SCRIPT] \
#       [-s "SBATCH OPTIONS"] [--dry-run] -- CLASSIFIER_ARGS...
#
#   CLASSIFIER_ARGS name the input, the dbCAN release and the output, for example
#     --proteins all.faa --dbcan-dir "$SHARUR_DBCAN_DIR" --results-only --out cazy.tsv
#   -s  site options for every job, e.g. "-p standard -c 24"
#   -m  array throttle: at most MAX_RUNNING tasks of each array at once
#   -e  script sourced at the start of every job (environment setup)
#   Every job gets its own TMPDIR under WORK_DIR/tmp and runs with $SLURM_CPUS_ON_NODE threads.
#   The job scripts land in WORK_DIR/slurm/ before submission; --dry-run stops there.
set -euo pipefail

usage() { sed -n '2,19p' "$0" >&2; exit 2; }
n="" work="" max="" env_script="" opts="" dry=0
while [[ $# -gt 0 ]]; do
  case $1 in
    -n) n=$2; shift 2 ;;
    -w) work=$2; shift 2 ;;
    -m) max=$2; shift 2 ;;
    -e) env_script=$2; shift 2 ;;
    -s) opts=$2; shift 2 ;;
    --dry-run) dry=1; shift ;;
    --) shift; break ;;
    *) usage ;;
  esac
done
[[ -n $n && -n $work && $# -gt 0 ]] || usage

script=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/classify_cazymes.py
mkdir -p "$work"
work=$(cd "$work" && pwd)
mkdir -p "$work/slurm" "$work/logs" "$work/tmp"
args=$(printf '%q ' "$@")
python=${PYTHON:-python}
array="0-$((n - 1))${max:+%$max}"

write_job() {  # name stage
  local name=$1 stage=$2 f=$work/slurm/$1.sbatch
  {
    echo '#!/bin/bash'
    echo "#SBATCH -J $name"
    if [[ $stage == merge ]]; then
      echo "#SBATCH -o $work/logs/%x_%j.out"
    else
      echo "#SBATCH --array=$array"
      echo "#SBATCH -o $work/logs/%x_%A_%a.out"
    fi
    local o line=""  # one #SBATCH line per option; a word without a leading dash is its value
    for o in $opts; do
      if [[ $o == -* && -n $line ]]; then echo "#SBATCH $line"; line=""; fi
      line="${line:+$line }$o"
    done
    if [[ -n $line ]]; then echo "#SBATCH $line"; fi
    echo 'set -euo pipefail'
    if [[ -n $env_script ]]; then echo "source $(printf %q "$env_script")"; fi
    echo "export TMPDIR=\$(mktemp -d $(printf %q "$work/tmp")/\${SLURM_JOB_ID}.XXXXXX)"
    echo 'trap '\''rm -rf "$TMPDIR"'\'' EXIT'
    echo 'echo "host $(hostname) cpus $SLURM_CPUS_ON_NODE task ${SLURM_ARRAY_TASK_ID:-merge} start $(date)"'
    local index=""
    if [[ $stage != merge ]]; then index=' --shard-index $SLURM_ARRAY_TASK_ID'; fi
    echo "$python $(printf %q "$script") $args--work-dir $(printf %q "$work") --n-shards $n --stage $stage$index --threads \$SLURM_CPUS_ON_NODE"
    echo 'echo "done $(date)"'
  } > "$f"
  echo "$f"
}

search=$(write_job cazy_search search)
sub=$(write_job cazy_sub sub)
merge=$(write_job cazy_merge merge)
if (( dry )); then
  printf 'wrote %s\n' "$search" "$sub" "$merge"
  exit 0
fi

j1=$(sbatch --parsable "$search"); j1=${j1%%;*}
j2=$(sbatch --parsable --dependency=afterok:"$j1" "$sub"); j2=${j2%%;*}
j3=$(sbatch --parsable --dependency=afterok:"$j2" "$merge"); j3=${j3%%;*}
printf 'search\t%s\nsub\t%s\nmerge\t%s\n' "$j1" "$j2" "$j3" | tee "$work/jobs.tsv"
