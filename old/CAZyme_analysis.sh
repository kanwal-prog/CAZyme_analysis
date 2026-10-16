# Genome
# Download example prokaryotic genome (E. coli K-12 MG1655)
wget -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz -O EscheriaColiK12MG1655.fna.gz
gunzip EscheriaColiK12MG1655.fna.gz

# Run CAZyme annotation
run_dbcan CAZyme_annotation \
  --input_raw_data EscheriaColiK12MG1655.fna \
  --mode prok \
  --output_dir output_EscheriaColiK12MG1655_fna \
  --db_dir db
# CGC (CAZyme Gene Cluster) Annotation)
# Run CGC annotation with automatically generated GFF
run_dbcan easy_CGC \
  --input_raw_data EscheriaColiK12MG1655.fna \
  --mode prok \
  --output_dir output_EscheriaColiK12MG1655_fna_CGC \
  --db_dir db \
  --input_gff gff \
  --gff_type prodigal

# with external gff file
# Download example GFF file for E. coli K-12 MG1655
wget -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz -O EscheriaColiK12MG1655.gff.gz
gunzip EscheriaColiK12MG1655.gff.gz

# Run CGC annotation with provided GFF
run_dbcan easy_CGC \
  --input_raw_data EscheriaColiK12MG1655.fna \
  --mode prok \
  --output_dir output_EscheriaColiK12MG1655_fna_CGC_gff \
  --db_dir db \
  --input_gff EscheriaColiK12MG1655.gff \
  --gff_type NCBI_prok

# Substrate Prediction
# Run substrate prediction
run_dbcan easy_substrate \
  --input_raw_data EscheriaColiK12MG1655.fna \
  --mode prok \
  --output_dir output_EscheriaColiK12MG1655_fna_sub \
  --db_dir db \
  --input_gff gff \
  --gff_type prodigal

# Substrate Prediction
# Run substrate prediction
run_dbcan easy_substrate \
  --input_raw_data EscheriaColiK12MG1655.fna \
  --mode prok \
  --output_dir output_EscheriaColiK12MG1655_fna_sub_2 \
  --db_dir db \
  --input_gff EscheriaColiK12MG1655.gff \
  --gff_type NCBI_prok

# plot CGC circle plot
run_dbcan cgc_circle_plot --output_dir output_EscheriaColiK12MG1655_fna_sub_2


#######################################################################################
# Proteome
# Download example prokaryotic proteome
wget -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_protein.faa.gz -O EscheriaColiK12MG1655.faa.gz
gunzip EscheriaColiK12MG1655.faa.gz


# Run CAZyme annotation (specify input format for protein sequences)
run_dbcan CAZyme_annotation \
  --input_raw_data EscheriaColiK12MG1655.faa \
  --mode protein \
  --output_dir output_EscheriaColiK12MG1655_faa \
  --db_dir db

# CGC (CAZyme Gene Cluster) Annotation)
# Run CGC annotation with automatically generated GFF
run_dbcan easy_CGC \
  --input_raw_data EscheriaColiK12MG1655.faa \
  --mode protein \
  --output_dir output_EscheriaColiK12MG1655_faa_CGC \
  --db_dir db \
  --input_gff gff \
  --gff_type prodigal

# with external gff file
# Download example GFF file for E. coli K-12 MG1655
wget -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz -O EscheriaColiK12MG1655.gff.gz
gunzip EscheriaColiK12MG1655.gff.gz

# Run CGC annotation with provided GFF
run_dbcan easy_CGC \
  --input_raw_data EscheriaColiK12MG1655.faa \
  --mode protein \
  --output_dir output_EscheriaColiK12MG1655_faa_CGC_gff \
  --db_dir db \
  --input_gff EscheriaColiK12MG1655.gff \
  --gff_type NCBI_prok

# Substrate Prediction
# Run substrate prediction
run_dbcan easy_substrate \
  --input_raw_data EscheriaColiK12MG1655.faa \
  --mode protein \
  --output_dir output_EscheriaColiK12MG1655_faa_sub \
  --db_dir db \
  --input_gff gff \
  --gff_type prodigal

# Substrate Prediction
# Run substrate prediction
run_dbcan easy_substrate \
  --input_raw_data EscheriaColiK12MG1655.faa \
  --mode protein \
  --output_dir output_EscheriaColiK12MG1655_faa_sub_2 \
  --db_dir db \
  --input_gff EscheriaColiK12MG1655.gff \
  --gff_type NCBI_prok

# plot CGC circle plot
run_dbcan cgc_circle_plot --output_dir output_EscheriaColiK12MG1655_faa_sub_1
run_dbcan cgc_circle_plot --output_dir output_EscheriaColiK12MG1655_faa_sub_2





