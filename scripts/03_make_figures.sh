#!/usr/bin/env bash
# Step 3: all figures for the genomes in genomes.tsv.
#
#   conda activate dbcan
#   bash scripts/03_make_figures.sh
#
# Output in figures/:
#   teacher_style/fig1_cazyme_abundance   CAZymes per genome            (teacher's plot_dbcan.py)
#   teacher_style/fig2_pul_substrates     stacked PUL substrates        (teacher's plot_dbcan.py)
#   substrate_heatmap                     substrate x genome heatmap    (plot_substrates.py)
#   substrate_stacked                     stacked substrates            (plot_substrates.py)
#   circos_<genome>                       CAZyme circos plot per genome (cazyme_circos.py)
#   pul_figure                            PUL gene maps across genomes  (pul_figure.py)
# Each figure has a .tsv next to it with the numbers behind it.

set -euo pipefail
cd "$(dirname "$0")/.."

RES=dbcan_teacher
mkdir -p figures

# input folders and labels, in genomes.tsv order
dirs=() labels=()
while IFS=$'\t' read -r name _ _ label; do
  dirs+=("$RES/$name/4_sub_ncbi_gff")
  labels+=("$label")
done < <(tail -n +2 genomes.tsv)

echo "== teacher's figures (PDF + PNG)"
for fmt in pdf png; do
  python scripts/plot_dbcan.py -i "${dirs[@]}" -l "${labels[@]}" --fmt $fmt -o figures/teacher_style
done

echo "== substrate heatmap and stacked bars"
python scripts/plot_substrates.py --style heatmap -i "${dirs[@]}" -l "${labels[@]}" --fmt png -o figures/substrate_heatmap
python scripts/plot_substrates.py --style stacked -i "${dirs[@]}" -l "${labels[@]}" --fmt png -o figures/substrate_stacked

echo "== CAZyme circos plots"
python scripts/cazyme_circos.py all

echo "== PUL figure"
python scripts/pul_figure.py --fmt png
python scripts/pul_figure.py --fmt pdf

echo "Done. Figures in figures/"
