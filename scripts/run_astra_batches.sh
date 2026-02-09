#!/bin/bash
# Run Astra annotation in batches for any database
# Usage: ./run_astra_batches.sh <DB_NAME> <BATCH_SIZE> <THREADS> [EXTRA_FLAGS]
# Example: ./run_astra_batches.sh VOGdb 200 12
# Example: ./run_astra_batches.sh HydDB 0 12 --cut_ga   (0 = no batching, run all at once)
# Example: ./run_astra_batches.sh DefenseFinder 400 12

set -e

DB_NAME="${1:?Usage: $0 <DB_NAME> <BATCH_SIZE> <THREADS> [EXTRA_FLAGS]}"
BATCH_SIZE="${2:?Usage: $0 <DB_NAME> <BATCH_SIZE> <THREADS> [EXTRA_FLAGS]}"
THREADS="${3:?Usage: $0 <DB_NAME> <BATCH_SIZE> <THREADS> [EXTRA_FLAGS]}"
shift 3
EXTRA_FLAGS="$@"

SYMLINK_DIR="data/omni_production/stage03_prodigal/genomes/all_protein_symlinks"
OUTBASE="data/omni_production/stage04_astra"
DB_LOWER=$(echo "$DB_NAME" | tr '[:upper:]' '[:lower:]')
RESULT_DIR="$OUTBASE/${DB_LOWER}_results"
COMBINED_OUT="$RESULT_DIR/${DB_NAME}_hits_df.tsv"

mkdir -p "$RESULT_DIR"

# No batching mode
if [ "$BATCH_SIZE" -eq 0 ]; then
    echo "Running $DB_NAME on all genomes (no batching)..."
    astra search --installed_hmms "$DB_NAME" --threads "$THREADS" \
        --prot_in "$SYMLINK_DIR" \
        --outdir "$RESULT_DIR" \
        $EXTRA_FLAGS 2>&1 | tail -10

    if [ -f "$COMBINED_OUT" ]; then
        hits=$(($(wc -l < "$COMBINED_OUT") - 1))
        echo "$DB_NAME complete: $hits hits -> $COMBINED_OUT"
    else
        echo "ERROR: $COMBINED_OUT not created!"
        exit 1
    fi
    exit 0
fi

# Batched mode
BATCH_DIR="$OUTBASE/${DB_LOWER}_batches"

# Create batch directories if they don't exist
if [ ! -d "$BATCH_DIR" ]; then
    echo "Creating batch directories (batch size: $BATCH_SIZE)..."
    python3 -c "
import os, math
symlink_dir = '$SYMLINK_DIR'
batch_base = '$BATCH_DIR'
files = sorted(os.listdir(symlink_dir))
batch_size = $BATCH_SIZE
n_batches = math.ceil(len(files) / batch_size)
print(f'Total files: {len(files)}, Batches: {n_batches}')
for i in range(n_batches):
    batch_dir = os.path.join(batch_base, f'batch_{i:02d}')
    os.makedirs(batch_dir, exist_ok=True)
    for f in files[i*batch_size:(i+1)*batch_size]:
        src = os.path.realpath(os.path.join(symlink_dir, f))
        os.symlink(src, os.path.join(batch_dir, f))
"
else
    echo "Using existing batch directories at $BATCH_DIR"
fi

TOTAL_BATCHES=$(ls -d "$BATCH_DIR"/batch_* | wc -l | tr -d ' ')
echo "Running $DB_NAME on $TOTAL_BATCHES batches with $THREADS threads..."

for batch_dir in "$BATCH_DIR"/batch_*; do
    batch_name=$(basename "$batch_dir")
    out_dir="$OUTBASE/${DB_LOWER}_batch_results/$batch_name"

    if [ -f "$out_dir/${DB_NAME}_hits_df.tsv" ]; then
        echo "[$batch_name] Already completed, skipping"
        continue
    fi

    n_files=$(ls "$batch_dir"/*.faa 2>/dev/null | wc -l | tr -d ' ')
    echo "[$batch_name] Processing $n_files genomes..."

    astra search --installed_hmms "$DB_NAME" --threads "$THREADS" \
        --prot_in "$batch_dir" \
        --outdir "$out_dir" \
        $EXTRA_FLAGS 2>&1 | tail -3

    if [ -f "$out_dir/${DB_NAME}_hits_df.tsv" ]; then
        hits=$(wc -l < "$out_dir/${DB_NAME}_hits_df.tsv")
        echo "[$batch_name] Complete: $hits hits"
    else
        echo "[$batch_name] WARNING: No output file created!"
    fi
done

# Combine all batch results
echo ""
echo "Combining $DB_NAME batch results..."
first=true
for result_file in "$OUTBASE/${DB_LOWER}_batch_results"/batch_*/${DB_NAME}_hits_df.tsv; do
    if [ ! -f "$result_file" ]; then
        continue
    fi
    if $first; then
        cat "$result_file" > "$COMBINED_OUT"
        first=false
    else
        tail -n +2 "$result_file" >> "$COMBINED_OUT"
    fi
done

if [ -f "$COMBINED_OUT" ]; then
    total_hits=$(( $(wc -l < "$COMBINED_OUT") - 1 ))
    echo "Combined $DB_NAME results: $total_hits hits -> $COMBINED_OUT"
    # Clean up batch results
    echo "Cleaning up batch results..."
    rm -rf "$OUTBASE/${DB_LOWER}_batch_results"
    rm -rf "$BATCH_DIR"
else
    echo "ERROR: No results to combine!"
    exit 1
fi
