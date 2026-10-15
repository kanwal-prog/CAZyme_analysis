set -e
mkdir -p ecoli bfragilis

#___E.coli

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz \
  -O ecoli/Ecoli_k12_MG1655.fna.gz

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/865/GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2_genomic.fna.gz \
  -O ecoli/Ecoli_0157H7_Sakai.fna.gz

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/445/GCF_000007445.1_ASM744v1/GCF_000007445.1_ASM744v1_genomic.fna.gz \
  -O ecoli/Ecoli_CFT073.fna.gz

##Bacteroides fragilis (pro cancer)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/025/985/GCF_000025985.1_ASM2598v1/GCF_000025985.1_ASM2598v1_genomic.fna.gz \
  -O bfragilis/Bfrag_NCTC9343.fna.gz

#YCH46 (GCF_000009925.1)

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/925/GCF_000009925.1_ASM992v1/GCF_000009925.1_ASM992v1_genomic.fna.gz \
    -O bfragilis/Bfrag_YCH46.fna.gz


# ZY0302 (ETBF, CRC tumor)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/045/689/415/GCF_045689415.1_ASM4568941v1/GCF_045689415.1_ASM4568941v1_genomic.fna.gz \
    -O bfragilis/Bfrag_ZY0302.fna.gz

##ZY0804 (non-toxigenic, paracancerous)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/045/689/425/GCF_045689425.1_ASM4568942v1/GCF_045689425.1_ASM4568942v1_genomic.fna.gz \
    -O bfragilis/Bfrag_ZY0804.fna.gz



## ---------- Decompress ----------
gunzip -f ecoli/*.gz bfragilis/*.gz

## ---------- Check: one header per genome ----------
for f in ecoli/*.fna bfragilis/*.fna; do
  echo "$f : $(head -1 "$f")"
done

## ---------- Zip everything ----------
zip -r pilot_genomes.zip ecoli bfragilis   # sudo apt install zip, if missing

df -h ~ | tail -1