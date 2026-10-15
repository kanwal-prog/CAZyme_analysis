cd ~/Desktop/CAZyme_analysis
ls ecoli/anvio ecoli/anvio/ECOLI_PAN
ls bfragilis/anvio
cat ecoli/anvio/completeness.txt
cat bfragilis/anvio/completeness.txt

## to check external genome

cd ~/Desktop/CAZyme_analysis
anvi-estimate-scg-taxonomy -e ecoli/anvio/external-genomes.txt -o ecoli_taxonomy.txt
anvi-estimate-scg-taxonomy -e bfragilis/anvio/external-genomes.txt -o bfragilis_taxonomy.txt

#ANI TABLE
cd ~/Desktop/CAZyme_analysis
ls bfragilis/anvio/ANI
cat bfragilis/anvio/ANI/ANIb_percentage_identity.txt
cat ecoli/anvio/ANI/ANIb_percentage_identity.txt

#COMPLETEMESS
cat ecoli/anvio/completeness.txt
cat bfragilis/anvio/completeness.txt

#OPen pangenome 
cd ~/Desktop/CAZyme_analysis
anvi-display-pan -g ecoli/anvio/ECOLI-GENOMES.db -p ecoli/anvio/ECOLI_PAN/ECOLI_PAN-PAN.

#Calculate cazyme i ecoli
python - <<'EOF'
import pandas as pd
df = pd.read_csv('ecoli/anvio/ECOLI_PAN_SUMMARY_final/ECOLI_PAN_gene_clusters_summary.txt.gz', sep='\t')
c = df.groupby('gene_cluster_id').agg(bin=('bin_name','first'), cazyme=('CAZyme', lambda s: s.notna().any() and (s.astype(str)!='').any()))
print(c.groupby('bin')['cazyme'].agg(['sum','count']).assign(percent=lambda d: (100*d['sum']/d['count']).round(1)))
EOF

#!/usr/bin/env bash
# =====================================================================
# Commands to run AFTER anvio_pangenome.sh has finished.
# Run them section by section (copy and paste), NOT as one script:
# the display command stays open until you press Ctrl+C.
#
# Set the species once, then run the sections in order.
# =====================================================================

# --- choose the species (E. coli shown; B. fragilis values below) ----
SP=ecoli;      N=3        # N = number of genomes in this species
# SP=bfragilis; N=4
SPU=${SP^^}               # ECOLI or BFRAGILIS
COLL=final                # name of the bin collection (use a new name if "already exists")
export SP SPU N COLL

cd ~/Desktop/CAZyme_analysis
conda activate anvio-9


# ---------------------------------------------------------------------
# 1. Genome quality: completeness and redundancy
#    Good = completeness near 100%, redundancy low (a few % at most).
# ---------------------------------------------------------------------
cat $SP/anvio/completeness.txt


# ---------------------------------------------------------------------
# 2. Species check (coarse: usually only reaches family or genus)
# ---------------------------------------------------------------------
anvi-estimate-scg-taxonomy -e $SP/anvio/external-genomes.txt -o ${SP}_taxonomy.txt
cat ${SP}_taxonomy.txt


# ---------------------------------------------------------------------
# 3. Species check by ANI (better).
#    About 95% or higher between two genomes = same species.
#    About 87% = different division/species (B. fragilis Division I vs II).
# ---------------------------------------------------------------------
ls $SP/anvio/ANI
cat $SP/anvio/ANI/ANIb_percentage_identity.txt


# ---------------------------------------------------------------------
# 4. Totals Anvi'o reported for each pangenome (keep these numbers!)
# ---------------------------------------------------------------------
grep -E "Your pangenome is ready|Num AA sequences reported" run.log


# ---------------------------------------------------------------------
# 5. Open the pangenome in the browser (Ctrl+C in the terminal to close).
#    In the browser: click Draw. Bins tab > Load bin collection > final
# ---------------------------------------------------------------------
anvi-display-pan -g $SP/anvio/${SPU}-GENOMES.db \
                 -p $SP/anvio/${SPU}_PAN/${SPU}_PAN-PAN.db


# ---------------------------------------------------------------------
# 6. Make core / accessory / singleton bins by rule (no clicking).
#    core      = cluster in all N genomes
#    accessory = in 2 to N-1 genomes
#    singleton = in 1 genome
#    Printed total must equal the pangenome total from section 4.
# ---------------------------------------------------------------------
python - <<'EOF'
import os, sqlite3, pandas as pd
sp, n = os.environ['SP'], int(os.environ['N'])
spu = os.environ['SPU']
db = f'{sp}/anvio/{spu}_PAN/{spu}_PAN-PAN.db'
t = pd.read_sql_query("select gene_cluster_id, genome_name from gene_clusters",
                      sqlite3.connect(db))
g = t.groupby('gene_cluster_id')['genome_name'].nunique()
lab = lambda k: 'core' if k == n else ('singleton' if k == 1 else 'accessory')
out = pd.DataFrame({'gene_cluster_id': g.index, 'bin_name': [lab(k) for k in g.values]})
out.to_csv(f'{sp}/anvio/{sp}_bins.tsv', sep='\t', header=False, index=False)
print(out.bin_name.value_counts())
print('total clusters:', len(out))
EOF

# Store the bins in the pangenome database, then check the count
anvi-import-collection $SP/anvio/${SP}_bins.tsv \
    -p $SP/anvio/${SPU}_PAN/${SPU}_PAN-PAN.db -C $COLL
anvi-show-collections-and-bins -p $SP/anvio/${SPU}_PAN/${SPU}_PAN-PAN.db

# Summary tables and web report for the collection
anvi-summarize -g $SP/anvio/${SPU}-GENOMES.db \
    -p $SP/anvio/${SPU}_PAN/${SPU}_PAN-PAN.db \
    -C $COLL -o $SP/anvio/${SPU}_PAN_SUMMARY_$COLL

# Open the report:
xdg-open $SP/anvio/${SPU}_PAN_SUMMARY_$COLL/index.html


# ---------------------------------------------------------------------
# 7. Share of gene clusters with a CAZyme hit, per bin
# ---------------------------------------------------------------------
python - <<'EOF'
import os, pandas as pd
sp, spu, coll = os.environ['SP'], os.environ['SPU'], os.environ['COLL']
f = f'{sp}/anvio/{spu}_PAN_SUMMARY_{coll}/{spu}_PAN_gene_clusters_summary.txt.gz'
df = pd.read_csv(f, sep='\t')
c = df.groupby('gene_cluster_id').agg(
    bin=('bin_name', 'first'),
    cazyme=('CAZyme', lambda s: s.notna().any() and (s.astype(str) != '').any()))
print(c.groupby('bin')['cazyme'].agg(['sum', 'count'])
        .assign(percent=lambda d: (100 * d['sum'] / d['count']).round(1)))
EOF


# ---------------------------------------------------------------------
# 8. Back up the results (keeps the databases and the log)
# ---------------------------------------------------------------------
zip -r ${SP}_results.zip $SP/anvio/${SPU}-GENOMES.db $SP/anvio/${SPU}_PAN \
    $SP/anvio/${SPU}_PAN_SUMMARY_$COLL $SP/anvio/completeness.txt \
    $SP/anvio/external-genomes.txt run.log