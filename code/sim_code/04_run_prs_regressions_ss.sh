#!/bin/bash

# set path to code
CODE="/projects/0/vusr0599/dl-prs/code"
INPUT="/projects/0/vusr0599/dl-prs/data/sim_data_ready"
OUT="/projects/0/vusr0599/dl-prs/output/regression_out"

# Dynamically name the log file
log_file="${OUT}/04_run_prs_regressions_$(date +%Y%m%d_%H%M%S).log"

# Redirect all output to tee
exec > >(tee "$log_file") 2> >(tee -a "$log_file" >&2)

# Print starting message
echo ""
echo "|---- Starting to run PRS regressions ----|"
echo ""

# Loop through files in directory and run 03_ready_data.R
for data_file in $(ls $INPUT/*DATA*.txt); do 

  # get basename for data file
  data=$(basename "${data_file}")

  # Check if output file exists
  if [ -f "${OUT}/ADD_R2_OUT_${base}.txt" ]; then
    echo "Output file already exists for $base. Skipping..."
    echo ""
    continue
  fi

  # print starting message
  echo "... Running 04_prs_regression.R ..."
  echo ""
  echo "  - Using file: $data"
  echo "  - Output directory: $OUT"
  echo ""

  # Run the R script
  Rscript "$CODE/04_prs_regressions_ss.R" \
    "$data_file" \
    "$OUT" \

  echo "Completed 04_prs_regressions.R for $data"
  echo ""

done

# Print finishing message
echo "|---- FINISHED for all files in ${INPUT} ----|"