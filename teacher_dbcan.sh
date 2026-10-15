#!/usr/bin/env bash
# dbCAN workflow from the teacher's example, applied to all 7 genomes.
#
#   conda activate dbcan
#   cd ~/Desktop/CAZyme_analysis
#   bash teacher_dbcan.sh
#
# Output: dbcan_teacher/<genome>/{1_cazyme, 2_cgc, 3_sub, 4_sub_ncbi_gff}

set -euo pipefail

DB=~/dbcan_db        # the database folder you downloaded (teacher's example calls it "db")
OUT=dbcan_teacher

# genome file -> NCBI assembly folder (from downloads.sh); used to get the matching GFF
declare -A ASM=(
  [ecoli/Ecoli_k12_MG1655.fna]=GCF/000/005/845/GCF_000005845.2_ASM584v2
  [ecoli/Ecoli_0157H7_Sakai.fna]=GCF/000/008/865/GCF_000008865.2_ASM886v2
  [ecoli/Ecoli_CFT073.fna]=GCF/000/007/445/GCF_000007445.1_ASM744v1
  [bfragilis/Bfrag_NCTC9343.fna]=GCF/000/025/985/GCF_000025985.1_ASM2598v1
  [bfragilis/Bfrag_YCH46.fna]=GCF/000/009/925/GCF_000009925.1_ASM992v1
  [bfragilis/Bfrag_ZY0302.fna]=GCF/045/689/415/GCF_045689415.1_ASM4568941v1
  [bfragilis/Bfrag_ZY0804.fna]=GCF/045/689/425/GCF_045689425.1_ASM4568942v1
)

[[ -d $DB ]] || { echo "ERROR: database folder $DB not found"; exit 1; }

for fna in "${!ASM[@]}"; do
  name=$(basename "$fna" .fna)
  dir=$OUT/$name
  mkdir -p "$dir"
  echo "=============== $name ==============="

  # step 1: genome .fna is already in ecoli/ or bfragilis/ (downloads.sh)

  # step 2: CAZyme annotation
  [[ -s $dir/1_cazyme/overview.tsv ]] || \
  run_dbcan CAZyme_annotation \
    --input_raw_data "$fna" \
    --mode prok \
    --output_dir "$dir/1_cazyme" \
    --db_dir "$DB"

  # step 3: CGC annotation with automatically generated GFF
  [[ -s $dir/2_cgc/cgc_standard_out.tsv ]] || \
  run_dbcan easy_CGC \
    --input_raw_data "$fna" \
    --mode prok \
    --output_dir "$dir/2_cgc" \
    --db_dir "$DB"

  # step 4: substrate prediction with automatically generated GFF
  [[ -s $dir/3_sub/substrate_prediction.tsv ]] || \
  run_dbcan easy_substrate \
    --input_raw_data "$fna" \
    --mode prok \
    --output_dir "$dir/3_sub" \
    --db_dir "$DB"

  # step 4b: download the NCBI GFF + NCBI proteins and run substrate prediction with them
  # (proteins must come from NCBI too, so their IDs match the GFF's protein_id)
  asm=${ASM[$fna]}
  gff=$dir/$name.gff
  faa=$dir/$name.faa
  if [[ ! -s $gff ]]; then
    wget -q "https://ftp.ncbi.nlm.nih.gov/genomes/all/$asm/$(basename "$asm")_genomic.gff.gz" -O "$gff.gz"
    gunzip -f "$gff.gz"
  fi
  if [[ ! -s $faa ]]; then
    wget -q "https://ftp.ncbi.nlm.nih.gov/genomes/all/$asm/$(basename "$asm")_protein.faa.gz" -O "$faa.gz"
    gunzip -f "$faa.gz"
  fi
  [[ -s $dir/4_sub_ncbi_gff/substrate_prediction.tsv ]] || \
  run_dbcan easy_substrate \
    --input_raw_data "$faa" \
    --mode protein \
    --output_dir "$dir/4_sub_ncbi_gff" \
    --db_dir "$DB" \
    --input_gff "$gff" \
    --gff_type NCBI_prok

  # step 5: CGC circle plot
  run_dbcan cgc_circle_plot --output_dir "$dir/4_sub_ncbi_gff"
done

echo "Done. Results in $OUT/"
