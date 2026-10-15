#!/bin/bash
# =============================================================================
# COMPLETE SECRETION SYSTEM + AMR PIPELINE
# 3 bacteria x 5 strains = 15 genomes
# Tools: TXSScan, ConjScan, TFFscan, CasFinder, RGI/CARD
# Author: your name | Supervisor: your teacher
# =============================================================================

set -e  # stop on any error
set -u  # stop if variable undefined

echo "========================================================"
echo " SECRETION SYSTEM PIPELINE — STARTING"
echo " Date: $(date)"
echo "========================================================"


# =============================================================================
# STEP 0 — SETUP FOLDER STRUCTURE
# =============================================================================
echo ""
echo "[STEP 0] Creating folder structure..."

mkdir -p raw_genomes
mkdir -p prokka_out
mkdir -p txsscan_out
mkdir -p conjscan_out
mkdir -p tffscan_out
mkdir -p casfinder_out
mkdir -p rgi_out
mkdir -p results
mkdir -p figures
mkdir -p logs

echo "Folders created."

# Define your 15 strain names — used throughout the pipeline
STRAINS=(
  "Vc_N16961"
  "Vc_O395"
  "Vc_IEM101"
  "Vc_RC385"
  "Vc_MO10"
  "Ec_Sakai"
  "Ec_E234869"
  "Ec_EC958"
  "Ec_SMS35"
  "Ec_HS"
  "Sal_LT2"
  "Sal_DT104"
  "Sal_P125109"
  "Sal_CVM19633"
  "Sal_Ty2"
)

echo "15 strains defined: ${STRAINS[@]}"


# =============================================================================
# STEP 1 — DOWNLOAD GENOMES FROM NCBI
# =============================================================================
echo ""
echo "[STEP 1] Downloading 15 genomes from NCBI RefSeq..."

# Requires: conda install -c conda-forge ncbi-datasets-cli
datasets download genome accession \
  AE003852 CP000626 CP028827 CP009967 CP001233 \
  BA000007 FM180568 HG941718 CP000970 CP000802 \
  AE006468 CP016408 AM933172 CP001127 AE014613 \
  --include genome \
  --filename raw_genomes/all15.zip \
  2>&1 | tee logs/download.log

unzip -o raw_genomes/all15.zip -d raw_genomes/
echo "Download complete."

# Rename downloaded .fna files to match our strain names
# (NCBI puts them in subfolders by accession — we copy and rename)
ACCESSIONS=(
  AE003852 CP000626 CP028827 CP009967 CP001233
  BA000007 FM180568 HG941718 CP000970 CP000802
  AE006468 CP016408 AM933172 CP001127 AE014613
)

for i in "${!ACCESSIONS[@]}"; do
  ACC="${ACCESSIONS[$i]}"
  NAME="${STRAINS[$i]}"
  FNA=$(find raw_genomes/ncbi_dataset/data/ -name "*.fna" | grep "$ACC" | head -1)
  if [ -n "$FNA" ]; then
    cp "$FNA" "raw_genomes/${NAME}.fna"
    echo "  Renamed: $ACC -> ${NAME}.fna"
  else
    echo "  WARNING: Could not find .fna for $ACC"
  fi
done

echo "All genomes renamed."


# =============================================================================
# STEP 2 — QUALITY CHECK WITH QUAST
# =============================================================================
echo ""
echo "[STEP 2] Running QUAST quality check on all genomes..."

# Requires: conda install -c bioconda quast
quast.py raw_genomes/*.fna \
  --output-dir results/quast_report \
  --threads 4 \
  2>&1 | tee logs/quast.log

echo "QUAST complete. Check results/quast_report/report.html"
echo "  Key things to verify:"
echo "  - N50 > 100000 (good assembly)"
echo "  - Number of contigs < 200 (complete genome)"
echo "  - Total length matches expected genome size"


# =============================================================================
# STEP 3 — GENOME ANNOTATION WITH PROKKA
# =============================================================================
echo ""
echo "[STEP 3] Annotating all 15 genomes with Prokka..."

# Requires: conda install -c bioconda prokka

for STRAIN in "${STRAINS[@]}"; do
  FNA="raw_genomes/${STRAIN}.fna"
  if [ ! -f "$FNA" ]; then
    echo "  SKIPPING $STRAIN — .fna not found"
    continue
  fi

  echo "  Annotating: $STRAIN"
  prokka \
    --kingdom Bacteria \
    --outdir "prokka_out/${STRAIN}" \
    --prefix "$STRAIN" \
    --force \
    --cpus 4 \
    "$FNA" \
    2>&1 >> logs/prokka.log

  echo "  Done: $STRAIN"
done

echo "Prokka annotation complete."
echo "  Output: prokka_out/<strain>/<strain>.faa (protein sequences)"
echo "  Output: prokka_out/<strain>/<strain>.gff (genome annotation)"


# =============================================================================
# STEP 4 — TXSSCAN (Type I-IX Secretion Systems + Flagellum)
# =============================================================================
echo ""
echo "[STEP 4] Running TXSScan on all 15 genomes..."

# Requires: conda install -c bioconda macsyfinder
# TXSScan models: pip install txsscan  OR  git clone from GitHub

for STRAIN in "${STRAINS[@]}"; do
  FAA="prokka_out/${STRAIN}/${STRAIN}.faa"
  if [ ! -f "$FAA" ]; then
    echo "  SKIPPING $STRAIN — .faa not found"
    continue
  fi

  echo "  Running TXSScan: $STRAIN"
  macsyfinder \
    --db-type ordered_replicon \
    --sequence-db "$FAA" \
    --models TXSScan all \
    --out-dir "txsscan_out/${STRAIN}" \
    --worker 4 \
    2>&1 >> logs/txsscan.log

  echo "  Done: $STRAIN"
done

echo "TXSScan complete."
echo "  Key output: txsscan_out/<strain>/best_solution.tsv"
echo "  Look for: T2SS, T3SS, T4SS, T5SS, T6SS, T9SS, Flagellum"


# =============================================================================
# STEP 5 — CONJSCAN (Conjugative Systems / Mobile Genetic Elements)
# =============================================================================
echo ""
echo "[STEP 5] Running ConjScan on all 15 genomes..."

# ConjScan uses same MacSyFinder engine as TXSScan
# Models: pip install conjscan  OR  from MacSyFinder model hub

for STRAIN in "${STRAINS[@]}"; do
  FAA="prokka_out/${STRAIN}/${STRAIN}.faa"
  if [ ! -f "$FAA" ]; then
    echo "  SKIPPING $STRAIN — .faa not found"
    continue
  fi

  echo "  Running ConjScan: $STRAIN"
  macsyfinder \
    --db-type ordered_replicon \
    --sequence-db "$FAA" \
    --models CONJScan all \
    --out-dir "conjscan_out/${STRAIN}" \
    --worker 4 \
    2>&1 >> logs/conjscan.log

  echo "  Done: $STRAIN"
done

echo "ConjScan complete."
echo "  Key output: conjscan_out/<strain>/best_solution.tsv"
echo "  Look for: CONJ_TypeIV, MOBF, MOBQ etc."


# =============================================================================
# STEP 6 — TFFSCAN (Type IV Filaments — Pili, TCP, MSHA)
# =============================================================================
echo ""
echo "[STEP 6] Running TFFscan on all 15 genomes..."

for STRAIN in "${STRAINS[@]}"; do
  FAA="prokka_out/${STRAIN}/${STRAIN}.faa"
  if [ ! -f "$FAA" ]; then
    echo "  SKIPPING $STRAIN — .faa not found"
    continue
  fi

  echo "  Running TFFscan: $STRAIN"
  macsyfinder \
    --db-type ordered_replicon \
    --sequence-db "$FAA" \
    --models TFFscan all \
    --out-dir "tffscan_out/${STRAIN}" \
    --worker 4 \
    2>&1 >> logs/tffscan.log

  echo "  Done: $STRAIN"
done

echo "TFFscan complete."
echo "  Key output: tffscan_out/<strain>/best_solution.tsv"
echo "  Look for: T4aP, T4bP (TCP in V. cholerae), Com pilus"


# =============================================================================
# STEP 7 — CASFINDER (CRISPR-Cas Systems)
# =============================================================================
echo ""
echo "[STEP 7] Running CasFinder on all 15 genomes..."

# Requires: conda install -c bioconda casfinder
# CasFinder takes .fna (nucleotide), not .faa

for STRAIN in "${STRAINS[@]}"; do
  FNA="raw_genomes/${STRAIN}.fna"
  if [ ! -f "$FNA" ]; then
    echo "  SKIPPING $STRAIN — .fna not found"
    continue
  fi

  echo "  Running CasFinder: $STRAIN"
  casfinder \
    -i "$FNA" \
    -o "casfinder_out/${STRAIN}" \
    -t 4 \
    2>&1 >> logs/casfinder.log

  echo "  Done: $STRAIN"
done

echo "CasFinder complete."
echo "  Key output: casfinder_out/<strain>/CRISPR_Cas.tsv"
echo "  Look for: Cas-I-E, Cas-I-F, Cas-II-A subtypes"


# =============================================================================
# STEP 8 — RGI / CARD (Antibiotic Resistance Genes)
# =============================================================================
echo ""
echo "[STEP 8] Running RGI (Resistance Gene Identifier) on all 15 genomes..."

# Requires: conda install -c bioconda rgi
# Update CARD database first: rgi load --card_json /path/to/card.json

# Download latest CARD database
echo "  Downloading CARD database..."
mkdir -p card_db
wget -q https://card.mcmaster.ca/latest/data -O card_db/card_data.tar.bz2
tar -xjf card_db/card_data.tar.bz2 -C card_db/
rgi load --card_json card_db/card.json --local 2>&1 >> logs/rgi.log

for STRAIN in "${STRAINS[@]}"; do
  FAA="prokka_out/${STRAIN}/${STRAIN}.faa"
  if [ ! -f "$FAA" ]; then
    echo "  SKIPPING $STRAIN — .faa not found"
    continue
  fi

  echo "  Running RGI: $STRAIN"
  rgi main \
    --input_sequence "$FAA" \
    --output_file "rgi_out/${STRAIN}_rgi" \
    --input_type protein \
    --alignment_tool BLAST \
    --num_threads 4 \
    --include_loose \
    --local \
    2>&1 >> logs/rgi.log

  echo "  Done: $STRAIN"
done

echo "RGI complete."
echo "  Key output: rgi_out/<strain>_rgi.txt"
echo "  Look for: Best_Hit_ARO, Drug_Class, Resistance_Mechanism"

echo ""
echo "========================================================"
echo " ALL SCANNING TOOLS COMPLETE"
echo " Now proceed to STEP 9: Post-analysis with Python"
echo "========================================================"
