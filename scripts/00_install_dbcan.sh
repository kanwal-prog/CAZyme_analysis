#!/usr/bin/env bash
# Step 0 (only once): create the "dbcan" conda environment and download the dbCAN database (~5 GB).
#
#   bash scripts/00_install_dbcan.sh

set -euo pipefail
eval "$(conda shell.bash hook)"

conda create -n dbcan -c conda-forge -c bioconda dbcan -y
conda activate dbcan
run_dbcan database --db_dir ~/dbcan_db --aws_s3
