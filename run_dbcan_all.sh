#!/usr/bin/env bash
# Run dbCAN on every genome, import consensus CAZymes into Anvi'o, and build a summary table.
#
# Usage:
#   bash run_dbcan_all.sh            # run all steps
#   bash run_dbcan_all.sh 1          # only step 1 (or 2, 3, 4)
#
# Before the first run (once):
#   conda create -n dbcan -c conda-forge -c bioconda dbcan -y
#   conda activate dbcan
#   run_dbcan database --db_dir ~/dbcan_db --aws_s3
#
# Steps:
#   1. dbCAN easy_substrate on each genome .fna  -> CAZymes, CGCs (PULs), substrates
#   2. export Anvi'o proteins, dbCAN on them      -> CAZymes with Anvi'o gene IDs
#   3. import consensus (#ofTools >= 2) into each contigs DB as source "dbCAN_consensus"
#   4. merged tables: CAZyme family x genome counts, CGC counts per genome

set -euo pipefail

PROJ=~/Desktop/CAZyme_analysis
DB_DIR=~/dbcan_db
THREADS=8
MIN_TOOLS=2
SOURCE_NAME=dbCAN_consensus
DBCAN_ENV=dbcan
ANVIO_ENV=anvio-9
STEP=${1:-all}

source ~/miniconda3/etc/profile.d/conda.sh
cd "$PROJ"

GENOMES=(ecoli/*.fna bfragilis/*.fna)

# species folder for a genome (ecoli or bfragilis), used to find its Anvi'o contigs DB
species_of() { dirname "$1"; }

# ---------------------------------------------------------------- step 1
if [[ $STEP == all || $STEP == 1 ]]; then
  conda activate "$DBCAN_ENV"
  [[ -d $DB_DIR ]] || { echo "ERROR: $DB_DIR not found, download the dbCAN database first"; exit 1; }
  mkdir -p dbcan/genome
  for f in "${GENOMES[@]}"; do
    name=$(basename "$f" .fna)
    out=dbcan/genome/$name
    if [[ -s $out/overview.tsv ]]; then echo "[1] skip $name (done)"; continue; fi
    echo "[1] dbCAN genome mode: $name"
    run_dbcan easy_substrate --mode prok \
      --input_raw_data "$f" \
      --db_dir "$DB_DIR" \
      --input_gff gff --gff_type prodigal \
      --output_dir "$out" \
      --threads "$THREADS"
  done
  conda deactivate
fi

# ---------------------------------------------------------------- step 2
if [[ $STEP == all || $STEP == 2 ]]; then
  mkdir -p dbcan/anvio_faa dbcan/protein

  conda activate "$ANVIO_ENV"
  for f in "${GENOMES[@]}"; do
    name=$(basename "$f" .fna)
    cdb=$(species_of "$f")/anvio/$name.db
    faa=dbcan/anvio_faa/$name.faa
    [[ -s $faa ]] && continue
    echo "[2] export Anvi'o proteins: $name"
    anvi-get-sequences-for-gene-calls -c "$cdb" --get-aa-sequences -o "$faa"
  done
  conda deactivate

  conda activate "$DBCAN_ENV"
  for f in "${GENOMES[@]}"; do
    name=$(basename "$f" .fna)
    out=dbcan/protein/$name
    if [[ -s $out/overview.tsv ]]; then echo "[2] skip $name (done)"; continue; fi
    echo "[2] dbCAN protein mode: $name"
    run_dbcan CAZyme_annotation --mode protein \
      --input_raw_data dbcan/anvio_faa/$name.faa \
      --db_dir "$DB_DIR" \
      --output_dir "$out" \
      --threads "$THREADS"
  done
  conda deactivate
fi

# ---------------------------------------------------------------- step 3
if [[ $STEP == all || $STEP == 3 ]]; then
  mkdir -p dbcan/import
  conda activate "$ANVIO_ENV"
  for f in "${GENOMES[@]}"; do
    name=$(basename "$f" .fna)
    cdb=$(species_of "$f")/anvio/$name.db
    ov=dbcan/protein/$name/overview.tsv
    tbl=dbcan/import/${name}_${SOURCE_NAME}.txt
    [[ -s $ov ]] || { echo "[3] missing $ov, run step 2 first"; exit 1; }

    # find columns by header name, so a column-order change in dbCAN does not break this
    awk -F'\t' -v OFS='\t' -v min="$MIN_TOOLS" -v src="$SOURCE_NAME" '
      NR==1 {
        for (i=1; i<=NF; i++) {
          if ($i ~ /^Gene ID/)        g=i
          if ($i ~ /#ofTools/)        n=i
          if ($i ~ /^Recommend/)      r=i
        }
        if (!g || !n || !r) { print "ERROR: unexpected overview.tsv header" > "/dev/stderr"; exit 1 }
        print "gene_callers_id", "source", "accession", "function", "e_value"
        next
      }
      $n >= min && $r != "-" { print $g, src, $r, $r, 0 }
    ' "$ov" > "$tbl"

    echo "[3] import $(($(wc -l < "$tbl") - 1)) CAZymes into $cdb"
    anvi-import-functions -c "$cdb" -i "$tbl"
  done
  conda deactivate
fi

# ---------------------------------------------------------------- step 4
if [[ $STEP == all || $STEP == 4 ]]; then
  echo "[4] building summary tables in dbcan/summary/"
  mkdir -p dbcan/summary
  python3 - "$MIN_TOOLS" <<'EOF'
import csv, glob, os, re, sys
from collections import Counter

min_tools = int(sys.argv[1])
fam_counts, cgc_counts = {}, {}

for d in sorted(glob.glob("dbcan/genome/*/")):
    name = os.path.basename(d.rstrip("/"))
    ov = os.path.join(d, "overview.tsv")
    if not os.path.exists(ov):
        continue
    c = Counter()
    with open(ov) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        tools = next(k for k in rd.fieldnames if "#ofTools" in k)
        rec = next(k for k in rd.fieldnames if k.startswith("Recommend"))
        for row in rd:
            if int(row[tools]) < min_tools or row[rec] == "-":
                continue
            # "GH13_7+CBM48" -> GH13, CBM48 (subfamilies collapsed to family)
            fams = {re.sub(r"_\d+$", "", x) for x in re.split(r"[+|]", row[rec]) if x}
            c.update(fams)
    fam_counts[name] = c

    cgc = os.path.join(d, "cgc_standard_out.tsv")
    if os.path.exists(cgc):
        with open(cgc) as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            key = next((k for k in rd.fieldnames if k.lower().startswith("cgc")), None)
            cgc_counts[name] = len({row[key] for row in rd}) if key else "NA"

genomes = sorted(fam_counts)
fams = sorted({f for c in fam_counts.values() for f in c})
with open("dbcan/summary/cazyme_family_counts.tsv", "w") as out:
    out.write("family\t" + "\t".join(genomes) + "\n")
    for f in fams:
        out.write(f + "\t" + "\t".join(str(fam_counts[g][f]) for g in genomes) + "\n")

with open("dbcan/summary/per_genome_totals.tsv", "w") as out:
    out.write("genome\tCAZymes\tGH\tGT\tPL\tCE\tAA\tCBM\tCGCs\n")
    for g in genomes:
        c = fam_counts[g]
        cls = {k: sum(v for f, v in c.items() if f.startswith(k)) for k in ["GH", "GT", "PL", "CE", "AA", "CBM"]}
        out.write("\t".join([g, str(sum(c.values()))] + [str(cls[k]) for k in cls] + [str(cgc_counts.get(g, "NA"))]) + "\n")

print(f"  {len(genomes)} genomes, {len(fams)} families")
EOF
fi

cat <<'EOF'

Done. Next, so the new annotations show up in the pangenome:
  conda activate anvio-9
  cd ~/Desktop/CAZyme_analysis/bfragilis/anvio      # and the same in ecoli/anvio
  anvi-gen-genomes-storage -e external-genomes.txt -o BFRAGILIS-GENOMES.db --force-overwrite
  anvi-pan-genome -g BFRAGILIS-GENOMES.db -n BFRAGILIS_PAN ...      # same options you used before
  anvi-compute-functional-enrichment-in-pan -p BFRAGILIS_PAN/BFRAGILIS_PAN-PAN.db \
     -g BFRAGILIS-GENOMES.db --category-variable <your_category> \
     --annotation-source dbCAN_consensus -o enrichment_dbCAN.txt
EOF
