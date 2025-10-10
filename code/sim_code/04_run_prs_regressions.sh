#!/bin/bash

# set path to code
CODE="/Users/nyb/phd_code/dl-prs-v2/code/nyb_sims_code_v2"
INPUT="/Users/nyb/phd_code/dl-prs-v2/data/nyb_sims_data_v2/ready_data/sim_data_ready"
OUT="/Users/nyb/phd_code/dl-prs-v2/output/man_sims/model_out_june2025/prs_out"

# Dynamically name the log file
log_file="${OUT}/04_run_prs_regressions_$(date +%Y%m%d_%H%M%S).log"

# Redirect all output to tee
exec > >(tee "$log_file") 2> >(tee -a "$log_file" >&2)

# Print starting message
echo ""
echo "|---- Starting to run PRS regressions ----|"
echo ""

# Loop through files in directory and run 03_ready_data.R
for data_file in $(ls $INPUT/*DATA_eur_nsnps100_*.txt); do 

  # get basename for data file
  data=$(basename "${data_file}")

  # print starting message
  echo "... Running 04_prs_regression.R ..."
  echo ""
  echo "  - Using file: $data"
  echo "  - Output directory: $OUT"
  echo ""

  # Run the R script
  Rscript "$CODE/04_prs_regressions.R" \
    "$data_file" \
    "$OUT" \

  echo "Completed 04_prs_regressions.R for $data"
  echo ""

done

# Print finishing message
echo "|---- FINISHED for all files in ${INPUT} ----|"