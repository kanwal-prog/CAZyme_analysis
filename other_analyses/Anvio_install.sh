
anvi-pan-genome --version
which anvi-gen-contigs-database
anvi-setup-cazymes -h | head -5

conda activate anvio-8
cd /home/reaz/Desktop/CAZyme_analysis
# First time only: remove the # from the 4 anvi-setup-* lines in step 0
bash anvio_pangenome.sh


df -h ~                 # check free space first
anvi-setup-cazymes

anvi-setup-scg-taxonomy
anvi-setup-ncbi-cogs -T 4

## Install AMR tool
conda create -n arg_tools -c conda-forge -c bioconda abricate -y
conda activate arg_tools
abricate --list

##Install dbcan
mamba create -n dbcan -c conda-forge -c bioconda dbcan
mamba activate dbcan
run_dbcan database --db_dir ~/dbcan_db

#Run dbcan on each genome
cd ~/Desktop/crc_project
for f in genomes/*.fna; do
  name=$(basename $f .fna)
  run_dbcan easy_substrate --mode prok \
    --input_raw_data $f \
    --db_dir ~/dbcan_db \
    --output_dir dbcan/$name
done

#Key outputs in each folder:

#overview.tsv: every CAZyme gene. Trust rows where #ofTools ≥ 2, meaning at least 2 of the 3 methods agree.
#cgc_standard_out.tsv: CAZyme gene clusters, which are the PUL-like regions.
#substrate_prediction.tsv: what each cluster likely breaks down, such as starch or mucin.

#3. Connect CAZymes to your pangenome

dbCAN’s gene IDs won’t match Anvi’o’s, so give dbCAN Anvi’o’s own proteins:

anvi-get-sequences-for-gene-calls -c ${name}.db --get-aa-sequences -o ${name}.faa
run_dbcan CAZyme_annotation --mode protein --input_raw_data ${name}.faa \
  --db_dir ~/dbcan_db --output_dir dbcan_aa/$name
  Convert the results to Anvi’o’s format and import them. Check the column order in your overview.tsv header first.

echo -e "gene_callers_id\tsource\taccession\tfunction\te_value" > ${name}_caz.txt
awk -F'\t' 'NR>1 && $6>=2 {print $1"\tCAZyme\t"$7"\t"$7"\t0"}' \
  dbcan_aa/$name/overview.tsv >> ${name}_caz.txt
anvi-import-functions -c ${name}.db -i ${name}_caz.txt

Then rebuild the pangenome. The CAZyme families will now appear on gene clusters, and you can run anvi-compute-functional-enrichment-in-pan --annotation-source CAZyme to test pathogenic vs commensal strains.

What you’ll get
A table of CAZyme families (GH, GT, PL, CE, CBM) per strain
Which CAZymes are core and which are strain-specific
Predicted substrates of each CAZyme cluster

