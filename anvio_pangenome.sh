#!/bin/bash
# Anvi'o pangenome: E. coli and B. fragilis as two SEPARATE pangenomes
# (different genera share almost no gene clusters, so never pool them).
# Run from: /home/reaz/Desktop/CAZyme_analysis
# Requires anvi'o 8:  conda activate anvio-8

set -e
THREADS=8

## ---------- 0. One-time database setup (skip if already done) ----------
# anvi-setup-scg-taxonomy
# anvi-setup-ncbi-cogs  -T $THREADS
# anvi-setup-kegg-data  -T $THREADS
# anvi-setup-cazymes                      # dbCAN HMMs, for the CAZyme focus

for SP in ecoli bfragilis; do
  mkdir -p $SP/anvio

  ## ---------- 1. Clean FASTA + build contigs databases ----------
  for f in $SP/*.fna; do
    name=$(basename "$f" .fna)            # e.g. Ecoli_CFT073, Bfrag_YCH46

    anvi-script-reformat-fasta "$f" \
      -o $SP/anvio/${name}.fa \
      --simplify-names \
      --prefix ${name} \
      --seq-type NT \
      -r $SP/anvio/${name}_rename_report.txt

    anvi-gen-contigs-database \
      -f $SP/anvio/${name}.fa \
      -o $SP/anvio/${name}.db \
      -n ${name} \
      -T $THREADS

    ## ---------- 2. Annotate ----------
    anvi-run-hmms          -c $SP/anvio/${name}.db -T $THREADS
    anvi-run-scg-taxonomy  -c $SP/anvio/${name}.db -T $THREADS
    anvi-run-ncbi-cogs     -c $SP/anvio/${name}.db -T $THREADS
    anvi-run-kegg-kofams   -c $SP/anvio/${name}.db -T $THREADS
    anvi-run-cazymes       -c $SP/anvio/${name}.db -T $THREADS
  done

  ## ---------- 3. External genomes file ----------
  anvi-script-gen-genomes-file \
    --input-dir $SP/anvio \
    -o $SP/anvio/external-genomes.txt

  # Sanity check: completeness/redundancy of each genome
  anvi-estimate-genome-completeness \
    -e $SP/anvio/external-genomes.txt \
    -o $SP/anvio/completeness.txt

  ## ---------- 4. Genomes storage ----------
  anvi-gen-genomes-storage \
    -e $SP/anvio/external-genomes.txt \
    -o $SP/anvio/${SP^^}-GENOMES.db

  ## ---------- 5. Pangenome ----------
  # mcl-inflation 10 = recommended for strains of the same species
  anvi-pan-genome \
    -g $SP/anvio/${SP^^}-GENOMES.db \
    --project-name ${SP^^}_PAN \
    --output-dir $SP/anvio/${SP^^}_PAN \
    --num-threads $THREADS \
    --minbit 0.5 \
    --mcl-inflation 10 \
    --enforce-hierarchical-clustering

  ## ---------- 6. ANI (adds ANI heatmap to the pan display) ----------
  anvi-compute-genome-similarity \
    -e $SP/anvio/external-genomes.txt \
    -o $SP/anvio/ANI \
    -p $SP/anvio/${SP^^}_PAN/${SP^^}_PAN-PAN.db \
    --program pyANI \
    --num-threads $THREADS
done

echo "Done. Visualise with the commands below."

## ---------- 7. Visualise (run one at a time, opens browser) ----------
# anvi-display-pan -g ecoli/anvio/ECOLI-GENOMES.db        -p ecoli/anvio/ECOLI_PAN/ECOLI_PAN-PAN.db
# anvi-display-pan -g bfragilis/anvio/BFRAGILIS-GENOMES.db -p bfragilis/anvio/BFRAGILIS_PAN/BFRAGILIS_PAN-PAN.db
#
# In the browser: select core / accessory / singleton bins,
# save them as a collection named "default", then:

## ---------- 8. Summarise (after saving the collection) ----------
# anvi-summarize -g ecoli/anvio/ECOLI-GENOMES.db -p ecoli/anvio/ECOLI_PAN/ECOLI_PAN-PAN.db \
#   -C default -o ecoli/anvio/ECOLI_PAN_SUMMARY
# anvi-summarize -g bfragilis/anvio/BFRAGILIS-GENOMES.db -p bfragilis/anvio/BFRAGILIS_PAN/BFRAGILIS_PAN-PAN.db \
#   -C default -o bfragilis/anvio/BFRAGILIS_PAN_SUMMARY

## ---------- 9. CAZyme presence/absence across genomes ----------
# anvi-display-functions -e bfragilis/anvio/external-genomes.txt \
#   --annotation-source CAZyme --profile-db bfragilis/anvio/BFRAGILIS_CAZYMES-PROFILE.db
# anvi-display-functions -e ecoli/anvio/external-genomes.txt \
#   --annotation-source CAZyme --profile-db ecoli/anvio/ECOLI_CAZYMES-PROFILE.db
