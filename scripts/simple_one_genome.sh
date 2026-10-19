#!/usr/bin/env bash
# dbCAN for ONE genome, written step by step like the teacher's example.
# Same commands as 02_run_dbcan.sh, without the loop over genomes.
#
#   conda activate dbcan
#   bash scripts/simple_one_genome.sh
#
# For another genome change only G, FNA and NCBI below (values are in genomes.tsv).
# Output: simple_output/<genome>_cazyme, _cgc, _sub, _sub_ncbi

set -e
cd "$(dirname "$0")/.."

G=Bfrag_NCTC9343                     # genome name
FNA=bfragilis/$G.fna                 # genome file
NCBI=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/025/985/GCF_000025985.1_ASM2598v1/GCF_000025985.1_ASM2598v1
DB=~/dbcan_db                        # dbCAN database
OUT=simple_output/$G

mkdir -p simple_output

# step 1: download the genome
[ -s $FNA ] || { wget -q ${NCBI}_genomic.fna.gz -O $FNA.gz && gunzip $FNA.gz; }

# step 2: CAZyme annotation
run_dbcan CAZyme_annotation \
  --input_raw_data $FNA \
  --mode prok \
  --output_dir ${OUT}_cazyme \
  --db_dir $DB

# step 3: CGC annotation (no --input_gff: dbCAN makes the GFF itself)
run_dbcan easy_CGC \
  --input_raw_data $FNA \
  --mode prok \
  --output_dir ${OUT}_cgc \
  --db_dir $DB

# step 4: substrate prediction
run_dbcan easy_substrate \
  --input_raw_data $FNA \
  --mode prok \
  --output_dir ${OUT}_sub \
  --db_dir $DB

# step 4b: substrate prediction with the NCBI GFF.
# Use the NCBI proteins (.faa) and --mode protein, so the protein IDs match the GFF.
wget -q ${NCBI}_genomic.gff.gz -O $OUT.gff.gz && gunzip -f $OUT.gff.gz
wget -q ${NCBI}_protein.faa.gz -O $OUT.faa.gz && gunzip -f $OUT.faa.gz
run_dbcan easy_substrate \
  --input_raw_data $OUT.faa \
  --mode protein \
  --output_dir ${OUT}_sub_ncbi \
  --db_dir $DB \
  --input_gff $OUT.gff \
  --gff_type NCBI_prok

# step 5: CGC circle plot
run_dbcan cgc_circle_plot --output_dir ${OUT}_sub_ncbi

echo "Done. Results in simple_output/"
