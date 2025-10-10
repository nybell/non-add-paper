#!/bin/bash

# set path to code
CODE="/projects/0/vusr0599/dl-prs/code/ukb_doug_code"
INPUT="/projects/0/vusr0599/dl-prs/data/ukb_doug/ready"
INPUT_DEF="/projects/0/vusr0599/dl-prs/data/ukb_doug/ready_def"

# Print starting message
echo ""
echo "|---- Re-splitting data (CLUMP 1): 70%, 15%, 15% ----|"
echo ""

# Loop through files in directory 
for data_file in "$INPUT"/*DATA_int*.txt; do

  # get basename for data file
  data=$(basename "${data_file}")

  # print starting message
  echo ""
  echo "  Re-splitting file: ${data}"
  echo "---------------------------"

  # Run the R script
  Rscript "$CODE/01c_resplit_UKB.R" \
    "${data_file}"

  echo "Completed 01c_resplit_UKB.R for ${data}"
  echo ""

done

# Print finishing message
echo "|---- FINISHED for all files in ${INPUT} ----|"

# Print starting message
echo ""
echo "|---- Re-splitting data (CLUMP 2): 70%, 15%, 15% ----|"
echo ""

# Loop through files in directory 
for data_file in "$INPUT_DEF"/*DATA_int*.txt; do

  # get basename for data file
  data=$(basename "${data_file}")

  # print starting message
  echo ""
  echo "  Re-splitting file: ${data}"
  echo "---------------------------"

  # Run the R script
  Rscript "$CODE/01c_resplit_UKB.R" \
    "${data_file}"

  echo "Completed 01c_resplit_UKB.R for ${data}"
  echo ""

done

# Print finishing message
echo "|---- FINISHED for all files in ${INPUT_DEF} ----|"