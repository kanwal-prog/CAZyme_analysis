#!/usr/bin/env bash
# =====================================================================
# AMR + virulence + plasmid/mobile-element layer for the pilot pangenomes
# Run SECTION BY SECTION (copy and paste), not as one script.
# Put map_hits_to_clusters.py and ani_heatmap.py in ~/Desktop/CAZyme_analysis
# =====================================================================
cd ~/Desktop/CAZyme_analysis
source "$(conda info --base)/etc/profile.d/conda.sh"


# ---------------------------------------------------------------------
# A. ANI figure (a clean PNG for a paper; the numbers come from the
#    ANI step your script already ran)
# ---------------------------------------------------------------------
conda activate anvio-9
python ani_heatmap.py ecoli
python ani_heatmap.py bfragilis
# How to see ANI INSIDE the Anvi'o pangenome display: see the notes in the chat.


# ---------------------------------------------------------------------
# B. Install the AMR / virulence tool in its OWN environment,
#    so anvio-9 is not touched
# ---------------------------------------------------------------------
conda create -n arg_tools -c conda-forge -c bioconda abricate -y
conda activate arg_tools
abricate --list          # shows which databases are installed and how many sequences
# Write the database versions (date, size) in your methods.


# ---------------------------------------------------------------------
# C. Search every genome.  Using the cleaned .fa files from anvi'o, so the
#    contig names match what anvi'o knows.
#    card / resfinder  = antibiotic resistance genes
#    vfdb              = virulence factors
#    plasmidfinder     = plasmid replicons
#    Cutoffs: >=90% identity, >=80% coverage (write them in your methods).
# ---------------------------------------------------------------------
for SP in ecoli bfragilis; do
  mkdir -p $SP/hits
  for f in $SP/anvio/*.fa; do
    g=$(basename "$f" .fa)
    for db in card resfinder vfdb plasmidfinder; do
      abricate --db $db --minid 90 --mincov 80 "$f" > $SP/hits/${g}_${db}.tsv
    done
  done
done

# Quick look: one table of the AMR genes, genomes as columns
for SP in ecoli bfragilis; do
  abricate --summary $SP/hits/*_card.tsv
  abricate --summary $SP/hits/*_vfdb.tsv | head -5
done


# ---------------------------------------------------------------------
# D. Export gene positions from anvi'o (so hits can be matched to genes)
# ---------------------------------------------------------------------
conda activate anvio-9
for SP in ecoli bfragilis; do
  mkdir -p $SP/anvio/gene_calls
  for db in $SP/anvio/*.db; do
    case "$db" in *GENOMES.db) continue;; esac
    g=$(basename "$db" .db)
    anvi-export-gene-calls -c "$db" --gene-caller pyrodigal-gv -o $SP/anvio/gene_calls/$g.txt
  done
done
# If it complains about the gene caller name, list the names:
#   anvi-export-gene-calls -c ecoli/anvio/Ecoli_k12_MG1655.db --list-gene-callers
# and use the name it shows after --gene-caller


# ---------------------------------------------------------------------
# E. Which pangenome gene cluster / bin is each hit in?
#    Writes <sp>/hits_mapped_to_clusters.tsv and <sp>/hits_presence_absence.tsv
#    and prints how many clusters in each bin look like mobile DNA.
# ---------------------------------------------------------------------
python map_hits_to_clusters.py ecoli ECOLI
python map_hits_to_clusters.py bfragilis BFRAGILIS


# ---------------------------------------------------------------------
# F. Optional, for horizontal gene transfer (heavier; add after E works)
#    - Genomic islands: upload each genome's GenBank file / accession to the
#      IslandViewer 4 web server (no install needed).
#    - Prophages and plasmids: geNomad (command-line, needs its own database
#      download of a few GB):
#         conda create -n genomad -c conda-forge -c bioconda genomad -y
#         conda activate genomad
#         genomad download-database .
#         genomad end-to-end --cleanup <genome.fa> <outdir> genomad_db
#    Then compare: do your AMR / virulence hits sit inside those regions?
#    (Compare their START/END with the island or prophage coordinates.)
# ---------------------------------------------------------------------
