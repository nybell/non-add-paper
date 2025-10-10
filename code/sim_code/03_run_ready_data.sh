#!/bin/bash

# Assign input arguments to variables
pheno_h2=$1

# set path to code
CODE="/Users/nyb/phd_code/dl-prs-v2/code/nyb_sims_code_v2"
DATA="/Users/nyb/phd_code/dl-prs-v2/data/nyb_sims_data_v2"
SYN="/Users/nyb/phd_code/dl-prs-v2/data/nyb_sims_data_v2/syn_phenos"
DOMDEV="/Users/nyb/phd_code/dl-prs-v2/data/nyb_sims_data_v2/syn_phenos/domdev_sumstats_h${pheno_h2}"
LINEAR="/Users/nyb/phd_code/dl-prs-v2/data/nyb_sims_data_v2/syn_phenos/linear_sumstats_h${pheno_h2}"
OUT="/Users/nyb/phd_code/dl-prs-v2/data/nyb_sims_data_v2/ready_data"
LOG_FILE="${OUT}/log_file.txt"

# Parameters for the simulation
geno="$DATA/genos.50k.snps.10k.raw"         # Genotype file

# Dynamically name the log file
log_file="${OUT}/03_run_ready_data_$(date +%Y%m%d_%H%M%S).log"

# Redirect all output to tee
exec > >(tee "$log_file") 2> >(tee -a "$log_file" >&2)

# Print starting message
echo "|---- Preparing data for phenotypes: pheno_h2=${pheno_h2} ----|"
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
  Rscript "$CODE/03_ready_data.R" \
    "$fam" \
    "$geno" \
    "$SYN/${snps}" \
    "$LINEAR/${linear_gwas}" \
    "$DOMDEV/${domdev_gwas}" \
    "$OUT" \

  echo "Completed 03_ready_data.R for $base"
  echo ""

done

# Print finishing message
echo "|---- Finished preparing data for phenotypes: pheno_h2=${pheno_h2} ----|"