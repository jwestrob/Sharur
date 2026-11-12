#!/usr/bin/env bash
# Download the example genome set listed in genomes.tsv from NCBI (GenBank assemblies).
#
#   bash download.sh OUT_DIR [genomes.tsv]
#
# Writes OUT_DIR/<genome_id>.fna (uncompressed, ready for `sharur-ingest --input-dir OUT_DIR`).
# Needs curl, gzip and md5sum. Each file is checked against NCBI's md5checksums.txt.
# Files run one at a time, which keeps the load on NCBI's servers polite.
set -euo pipefail

out_dir=${1:?usage: download.sh OUT_DIR [genomes.tsv]}
table=${2:-"$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/genomes.tsv"}
base=https://ftp.ncbi.nlm.nih.gov/genomes/all
mkdir -p "$out_dir"

tail -n +2 "$table" | while IFS=$'\t' read -r genome_id assembly assembly_name _rest; do
  target="$out_dir/$genome_id.fna"
  if [[ -s "$target" ]]; then
    echo "have  $genome_id"
    continue
  fi
  digits=${assembly#GCA_}; digits=${digits%%.*}
  stem="${assembly}_${assembly_name// /_}"
  url_dir="$base/GCA/${digits:0:3}/${digits:3:3}/${digits:6:3}/$stem"
  tmp="$out_dir/.$stem.fna.gz"
  curl -fsSL --retry 5 --retry-delay 5 -o "$tmp" "$url_dir/${stem}_genomic.fna.gz"
  expected=$(curl -fsSL --retry 5 --retry-delay 5 "$url_dir/md5checksums.txt" \
    | awk -v f="./${stem}_genomic.fna.gz" '$2 == f {print $1}')
  observed=$(md5sum "$tmp" | cut -d' ' -f1)
  if [[ -z "$expected" || "$expected" != "$observed" ]]; then
    echo "checksum mismatch for $assembly (expected '${expected}', got '${observed}')" >&2
    rm -f "$tmp"
    exit 1
  fi
  gzip -dc "$tmp" > "$target.part" && mv "$target.part" "$target" && rm -f "$tmp"
  echo "got   $genome_id  ($assembly)"
done
echo "genomes in $out_dir: $(ls "$out_dir"/*.fna | wc -l)"
