# create conda environment and install dbcan
conda create --name dbcan python=3.10 -y
conda activate dbcan
conda install dbcan -c conda-forge -c bioconda -y

# build database from dbcan5 command line tool will take time
# 5.3 GB total in size
run_dbcan database --db_dir db --aws_s3