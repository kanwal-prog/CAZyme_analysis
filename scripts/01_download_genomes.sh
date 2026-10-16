#!/usr/bin/env bash
# Step 1: download every genome in genomes.tsv from NCBI.
#
#   bash scripts/01_download_genomes.sh
#
# Output: <folder>/<name>.fna      e.g. bfragilis/Bfrag_NCTC9343.fna
# Genomes already downloaded are skipped.

set -euo pipefail
cd "$(dirname "$0")/.."

while IFS=$'\t' read -r name folder ncbi _ <&3; do
  fna=$folder/$name.fna
  if [[ -s $fna ]]; then
    echo "have  $fna"
    continue
  fi
  mkdir -p "$folder"
  echo "get   $fna"
  wget -q "https://ftp.ncbi.nlm.nih.gov/genomes/all/$ncbi/$(basename "$ncbi")_genomic.fna.gz" -O "$fna.gz"
  gunzip -f "$fna.gz"
done 3< <(tail -n +2 genomes.tsv)

# one line per genome: first FASTA header
for f in $(tail -n +2 genomes.tsv | awk -F'\t' '{print $2"/"$1".fna"}'); do
  echo "$f : $(head -1 "$f")"
done
