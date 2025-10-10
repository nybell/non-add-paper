#!/bin/bash

# Set base directories
CODE="/Users/nyb/phd_code/dl-prs-v2/code/nyb_sims_code_v2"
INPUT="/Users/nyb/phd_code/dl-prs-v2/data/nyb_sims_data_v2/ready_data/sim_data_ready"
OUTPUT="/Users/nyb/phd_code/dl-prs-v2/output/man_sims/model_out_june2025/dnn_out"
EPOCHS=50

# Dynamically name the log file
log_file="${OUTPUT}/04_run_dnn_models_$(date +%Y%m%d_%H%M%S).log"

# Redirect all output to tee
exec > >(tee "$log_file") 2> >(tee -a "$log_file" >&2)

# Print starting message
echo ""
echo "|---- Starting to run neural networks ----|"
echo ""

# Test file: DATA_eur_nsnps100_h0.1_a0.1_d0_50k.txt
# Loop through files in directory and run 03_ready_data.R
for input_file in $(ls $INPUT/*DATA_eur_nsnps100_*); do      # *DATA_eur_nsnps100_*
  # Extract the base name of the file (without directory and extension)
  base_name=$(basename "$input_file" .txt | sed 's/^DATA_eur_//')
  echo "base name:"
  echo $base_name
  echo ""

  # Construct the job name by adding the prefix
  job_name="DNN_${base_name}"

  # define the output file & model performance dir for test
  out_pkl="${OUTPUT}/${job_name}_epochs${EPOCHS}.pkl"
  writer_dir="${OUTPUT}/dnn_runs/${job_name}_epochs${EPOCHS}"
  # model_pkl="${OUTPUT}/dnn_models/${job_name}_epochs${EPOCHS}.pkl"


  # check if output file already exists
  exists=0

  if [[ -f "$out_pkl" ]]; then
    echo "Output file already existing: $out_pkl"
    exists=1
  fi

  # if [[ -f "$model_pkl" ]]; then
  #   echo "Models directory already existing: $models_dir"
  #   exists=1
  # fi

  if [[ -d "$writer_dir" ]]; then
    echo "Models directory already existing: $models_dir"
    exists=1
  fi

  if [[ $exists -eq 1 ]]; then
    echo "Skipping NN for $job_name"
    continue
  fi

  echo "All paths are clear. Proceeding with the script..."
  # make models output if necessary
  mkdir -p "${writer_dir}"

  # run dnn
python $CODE/run_dnn.py \
  --input_data "${input_file}" \
  --result_path "${out_pkl}" \
  --model_path "${models_dir}" \
  --writer "${writer_dir}" \
  --epochs "${EPOCHS}"

  # Log the submission for tracking
  echo "| ----- FINISHED FOR : ${job_name} ----- |"
  echo ""
done

# Print finishing message
echo ""
echo "|---- FINISHED for all files in ${INPUT} ----|"