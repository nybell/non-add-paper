#!/bin/bash

# set path to code
CODE="/projects/0/vusr0599/dl-prs/code/ukb_doug_code"
INPUT="/projects/0/vusr0599/dl-prs/data/ukb_doug/ready"
OUT="/projects/0/vusr0599/dl-prs/output/ukb_doug/prs_out"

# Dynamically name the log file
log_file="${OUT}/02_run_prs_regressions_$(date +%Y%m%d_%H%M%S).log"

# Redirect all output to tee
exec > >(tee "$log_file") 2> >(tee -a "$log_file" >&2)

# Print starting message
echo ""
echo "|---- Starting to run PRS regressions ----|"
echo ""

# Loop through files in directory 
for data_file in "$INPUT"/*DATA_int*.txt; do

  # get basename for data file
  data=$(basename "${data_file}")

  # Check if output file exists
  if [[ -f "${OUT}/ADD_GENO_MODEL_OUT_${data}.RData" && -f "${OUT}/DOM_GENO_MODEL_OUT_${data}.RData" ]]; then
    echo "Output files already exists for ${data}. Skipping..."
    echo ""
    continue
  fi

  # print starting message
  echo "... Running 02_prs_regressions_UKB.R ..."
  echo ""
  echo "  - Using file: ${data}"
  echo "  - Output directory: ${OUT}"
  echo ""

  TMP_OUT="${TMPDIR}/output"
  mkdir "${TMPDIR}/output"

  # Run the R script
  Rscript "$CODE/02_prs_regressions_UKB.R" \
    "${data_file}" \
    "${TMP_OUT}" 

  echo "Completed 02_prs_regressions_UKB.R for ${data}"
  echo "Copying output to ${OUT}"
  cp "${TMP_OUT}"/* "${OUT}/"
  rm -r "${TMP_OUT}"
  echo "Output copied to ${OUT}"
  echo ""

done

# Print finishing message
echo "|---- FINISHED for all files in ${INPUT} ----|"