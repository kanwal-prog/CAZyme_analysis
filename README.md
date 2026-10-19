# CAZyme analysis of *B. fragilis* and *E. coli*

Carbohydrate-active enzymes (CAZymes), CAZyme gene clusters (CGCs/PULs) and predicted
substrates for 4 *Bacteroides fragilis* and 3 *Escherichia coli* genomes, using dbCAN.

## Run order

Run every command from this folder (`~/Desktop/CAZyme_analysis`).

| Step | Command | What it does | Time |
|---|---|---|---|
| 0 | `bash scripts/00_install_dbcan.sh` | Install dbCAN + database (only once) | ~1 h |
| 1 | `bash scripts/01_download_genomes.sh` | Download the genomes listed in `genomes.tsv` | minutes |
| 2 | `conda activate dbcan`<br>`nohup nice -n 19 bash scripts/02_run_dbcan.sh > dbcan_run.log 2>&1 &` | dbCAN on every genome | ~1 h per genome |
| 3 | `conda activate dbcan`<br>`bash scripts/03_make_figures.sh` | All figures | ~2 min |

- Start step 2 **only once**. Check whether it is running: `pgrep -af 02_run_dbcan`.
  Watch it: `tail -f dbcan_run.log` (Ctrl+C stops watching, not the job).
- Steps 1 and 2 skip whatever is already finished, so after a shutdown just run them again.
- **To add a genome**: add one line to `genomes.tsv`, then run steps 1, 2, 3.

## Files

```
genomes.tsv                 the genome list every script reads (name, folder, NCBI path, label)
scripts/
  00_install_dbcan.sh       step 0
  01_download_genomes.sh    step 1
  02_run_dbcan.sh           step 2  (teacher's dbCAN workflow)
  03_make_figures.sh        step 3  (runs the four Python scripts below)
  plot_dbcan.py             teacher's script: CAZyme count bars + stacked PUL substrates
  plot_substrates.py        substrate heatmap / stacked / grouped bars
  cazyme_circos.py          CAZyme circos plot per genome
  pul_figure.py             PUL gene maps across genomes
  simple_one_genome.sh      the same dbCAN steps for ONE genome, written like the teacher's example
ecoli/, bfragilis/          genomes (.fna)                      not on GitHub
dbcan_teacher/<genome>/     dbCAN results                       not on GitHub
figures/                    all figures + a .tsv of the numbers behind each
old/                        earlier versions, kept for reference - do not use
pangenome_analysis/         Anvi'o pangenome, ANI, AMR/HGT - separate repo:
                            github.com/kanwal-prog/Anvio_pangenome_analysis
```

## dbCAN results per genome (`dbcan_teacher/<genome>/`)

| Folder | Content |
|---|---|
| `1_cazyme/overview.tsv` | every CAZyme gene, its family, and how many tools found it |
| `2_cgc/cgc_standard_out.tsv` | CAZyme gene clusters (genes from Prodigal) |
| `3_sub/substrate_prediction.tsv` | predicted substrate per cluster (genes from Prodigal) |
| `4_sub_ncbi_gff/` | **same, with NCBI proteins + GFF - used for all figures**; also `cgc_circos/` and `synteny_pdf/` from dbCAN |

A CAZyme is counted when at least 2 of dbCAN's 3 tools agree (`#ofTools >= 2`).

## Figures (`figures/`)

| File | Shows |
|---|---|
| `dbcan_summary/fig1_cazyme_abundance` | number of CAZymes per genome |
| `dbcan_summary/fig2_pul_substrates` | predicted PUL substrates per genome, stacked |
| `substrate_heatmap` | substrate x genome heatmap (number of clusters) |
| `substrate_stacked` | stacked substrates (counts dbCAN-sub substrate when no PUL matched) |
| `circos_<genome>` | CAZymes by class and strand, clusters by substrate group, GC skew |
| `pul_figure` | gene maps of selected PULs (mucin, host glycan, pectin, starch) across genomes |

Single figures, e.g. other PULs or two genomes only:

```bash
python scripts/pul_figure.py --pul PUL0380 PUL0556 --fmt pdf -o figures/pul_mucin
python scripts/plot_dbcan.py -i dbcan_teacher/Bfrag_NCTC9343/4_sub_ncbi_gff dbcan_teacher/Ecoli_k12_MG1655/4_sub_ncbi_gff \
  -l "B. fragilis NCTC9343" "E. coli K12" --fmt png -o figures/two_genomes
```

## Notes

- Use `--mode protein` with the NCBI `.faa` when giving dbCAN an NCBI GFF; with `--mode prok`
  the Prodigal gene IDs do not match the GFF and no clusters are found.
- Do not pass `--input_gff gff` with `--mode prok`; dbCAN makes the GFF itself.
- Substrate predictions are homology-based (bitscore in `substrate_prediction.tsv`); weak matches
  such as the *B. fragilis* "xylan" clusters need checking before you claim them.
