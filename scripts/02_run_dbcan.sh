#!/usr/bin/env bash
# Step 2: dbCAN on every genome in genomes.tsv (teacher's workflow).
#
#   conda activate dbcan
#   nohup nice -n 19 bash scripts/02_run_dbcan.sh > dbcan_run.log 2>&1 &
#   tail -f dbcan_run.log                          # watch; Ctrl+C stops watching only
#
# Start it ONCE (check with: pgrep -af 02_run_dbcan). Finished steps are skipped,
# so after a shutdown just start it again.
#
# Output: dbcan_teacher/<genome>/
#   1_cazyme/          CAZyme annotation               (overview.tsv)
#   2_cgc/             CAZyme gene clusters, Prodigal  (cgc_standard_out.tsv)
#   3_sub/             substrate prediction, Prodigal  (substrate_prediction.tsv)
#   4_sub_ncbi_gff/    substrate prediction with NCBI proteins + GFF  <- used for all figures
#                      (+ cgc_circos/, synteny_pdf/)

set -euo pipefail
cd "$(dirname "$0")/.."

DB=~/dbcan_db
OUT=dbcan_teacher

[[ -d $DB ]] || { echo "ERROR: database folder $DB not found (run scripts/00_install_dbcan.sh)"; exit 1; }

# a result counts as finished only if it has more than a header line
finished() { [[ -s $1 && $(wc -l < "$1") -gt 1 ]]; }

while IFS=$'\t' read -r name folder ncbi _ <&3; do
  fna=$folder/$name.fna
  dir=$OUT/$name
  mkdir -p "$dir"
  echo "=============== $name ==============="
  [[ -s $fna ]] || { echo "ERROR: $fna missing (run scripts/01_download_genomes.sh)"; exit 1; }

  # 1. CAZyme annotation
  finished "$dir/1_cazyme/overview.tsv" || \
  run_dbcan CAZyme_annotation \
    --input_raw_data "$fna" \
    --mode prok \
    --output_dir "$dir/1_cazyme" \
    --db_dir "$DB"

  # 2. CGC annotation (dbCAN predicts genes with Prodigal and makes its own GFF)
  finished "$dir/2_cgc/cgc_standard_out.tsv" || \
  run_dbcan easy_CGC \
    --input_raw_data "$fna" \
    --mode prok \
    --output_dir "$dir/2_cgc" \
    --db_dir "$DB"

  # 3. substrate prediction (Prodigal genes)
  finished "$dir/3_sub/substrate_prediction.tsv" || \
  run_dbcan easy_substrate \
    --input_raw_data "$fna" \
    --mode prok \
    --output_dir "$dir/3_sub" \
    --db_dir "$DB"

  # 4. substrate prediction with NCBI annotation.
  #    Proteins must come from NCBI too, so their IDs match the GFF's protein_id.
  gff=$dir/$name.gff
  faa=$dir/$name.faa
  url=https://ftp.ncbi.nlm.nih.gov/genomes/all/$ncbi/$(basename "$ncbi")
  if [[ ! -s $gff ]]; then
    wget -q "${url}_genomic.gff.gz" -O "$gff.gz" && gunzip -f "$gff.gz"
  fi
  if [[ ! -s $faa ]]; then
    wget -q "${url}_protein.faa.gz" -O "$faa.gz" && gunzip -f "$faa.gz"
  fi
  finished "$dir/4_sub_ncbi_gff/substrate_prediction.tsv" || \
  run_dbcan easy_substrate \
    --input_raw_data "$faa" \
    --mode protein \
    --output_dir "$dir/4_sub_ncbi_gff" \
    --db_dir "$DB" \
    --input_gff "$gff" \
    --gff_type NCBI_prok

  # 5. dbCAN's own CGC circle plot
  [[ -d $dir/4_sub_ncbi_gff/cgc_circos ]] || \
  run_dbcan cgc_circle_plot --output_dir "$dir/4_sub_ncbi_gff"
done 3< <(tail -n +2 genomes.tsv)

echo "Done. Results in $OUT/"
