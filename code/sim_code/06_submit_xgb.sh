#!/bin/bash

# Set base directories
BASE=/projects/0/vusr0599/dl-prs
CODE=${BASE}/code
DATA=${BASE}/data/sim_data_ready
OUTPUT=${BASE}/output
HYPER_TRIALS=100

# Iterate over all .fam files in the FAM_FILES directory
for input_file in $(ls $DATA/*DATA*.txt); do         # *DATA*.txt   
  # Extract the base name of the file (without directory and extension)
  base_name=$(basename "$input_file" .txt | sed 's/^DATA_eur_//')
  echo "base name:"
  echo $base_name
  echo ""

  # Construct the job name by adding the prefix
  job_name="XGB_${base_name}"

  # Define the output file
  out_pkl="${OUTPUT}/xgb_out/${job_name}_trials${HYPER_TRIALS}.pkl"

  # define the output file & model performance dir for test
  out_pkl="${OUTPUT}/xgb_out/${job_name}_trials${HYPER_TRIALS}.pkl"
  models_dir="${OUTPUT}/xgb_models/${job_name}_trials${HYPER_TRIALS}"

  # check if output file already exists
  exists=0

  if [[ -f "$out_pkl" ]]; then
    echo "Output file already existing: $out_pkl"
    exists=1
  fi

  # if [[ -d "$models_dir" ]]; then
  #   echo "Models directory already existing: $models_dir"
  #   exists=1
  # fi

  if [[ $exists -eq 1 ]]; then
    echo "Skipping job submission for $job_name"
    continue
  fi

  echo "All paths are clear. Proceeding with the script..."

  # Submit the job to the cluster
  sbatch -J "${job_name}" \
         -o "$BASE/logs/${job_name}_%j.o.txt" \
         -e "$BASE/logs/${job_name}_%j.e.txt" \
         $CODE/xgb.job.txt \
         $CODE \
         ${input_file} \
         ${out_pkl} \
         ${models_dir} \
         ${HYPER_TRIALS} \

  # Log the submission for tracking
  echo "Submitted job: $job_name using file: $input_file"
  echo ""
done
