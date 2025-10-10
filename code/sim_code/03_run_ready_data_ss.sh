#!/bin/bash

# Assign input arguments to variables

# set path to code
CODE="/projects/0/vusr0599/dl-prs/code"
SYN="/projects/0/vusr0599/dl-prs/data/50k_fam_fin_files"
STATS="/projects/0/vusr0599/hapnest/synthetic_data/data/outputs/discovery_eur_100k/merged_sumstats"
OUT="/projects/0/vusr0599/dl-prs/data/sim_data_ready"
LOG_FILE="${OUT}/log_file.txt"

# Parameters for the simulation
geno="/projects/0/vusr0599/dl-prs/data/snp_data/genos.50k.snps.10k.raw"       # Genotype file

# Dynamically name the log file
log_file="${OUT}/03_run_ready_data_ss$(date +%Y%m%d_%H%M%S).log"

# Redirect all output to tee
exec > >(tee "$log_file") 2> >(tee -a "$log_file" >&2)

# Print starting message
echo "|---- Preparing data ----|"
echo ""

# Loop through files in directory and run 03_ready_data.R
for fam in $(ls $SYN/*50k.fam); do

  # get basename for all files
  tmp=$(basename "${fam}")
  base=$(basename "${tmp}" .fam)

  # make SNP file name
  tmp2=$(basename "${tmp}" .fam | sed 's/^eur_//')
  snps="${tmp2}.RData"

  # Check if output file exists
  if [ -f "${OUT}/DATA_${base}.txt" ]; then
    echo "Output file already exists for $base. Skipping..."
    echo ""
    continue
  fi

  # make GWAS summary stats file name
  base2=$(basename "${base}" _50k)
  domdev_gwas="domdev_gwas_${base2}_100k.sumstats.txt"
  linear_gwas="linear_gwas_${base2}_100k.sumstats.txt"  
  
  # print starting message
  echo "... Running 03_ready_data.R for $base ..."
  echo ""
  echo "Using files:"
  echo "  - Genotype file: $geno"
  echo "  - SNP file: $snps"
  echo "  - Dominance deviation GWAS: $domdev_gwas"
  echo "  - Linear GWAS: $linear_gwas"
  echo "  - Output directory: $OUT"
  echo ""

  # Run the R script
  Rscript "$CODE/03_ready_data_ss.R" \
    "$fam" \
    "$geno" \
    "$SYN/${snps}" \
    "${STATS}/${linear_gwas}" \
    "${STATS}/${domdev_gwas}" \
    "$OUT" \

  echo "Completed 03_ready_data.R for $base"
  echo ""

done

# Print finishing message
echo "|---- Finished preparing data for phenotypes ----|"