#!/bin/bash

# Set base directories
BASE=/projects/0/vusr0599/dl-prs
CODE=${BASE}/code
DATA=${BASE}/data/sim_data_ready
OUTPUT=${BASE}/output
EPOCHS=100

# Iterate over all .fam files in the FAM_FILES directory
for input_file in $(ls $DATA/*DATA*.txt); do            # *DATA*.txt
  # Extract the base name of the file (without directory and extension)
  base_name=$(basename "$input_file" .txt | sed 's/^DATA_eur_//')
  echo "base name:"
  echo $base_name
  echo ""

  # Construct the job name by adding the prefix
  job_name="DNN_${base_name}"

  # # Define the output file names + paths (UNCOMMENT FOR FULL RUN)
  out_pkl="${OUTPUT}/dnn_out/${job_name}_epochs${EPOCHS}.pkl"    
  writer_dir="${OUTPUT}/model_runs/${job_name}_epochs${EPOCHS}"
  models_dir="${OUTPUT}/models/${job_name}_epochs${EPOCHS}"

  # check if output file already exists
  exists=0

  if [[ -f "$out_pkl" ]]; then
    echo "Output file already existing: $out_pkl"
    exists=1
  fi

  if [[ -d "$writer_dir" ]]; then
    echo "Writer directory already existing: $writer_dir"
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
         $CODE/dnn.job.txt \
         $CODE \
         ${input_file} \
         ${out_pkl} \
         ${models_dir} \
         ${writer_dir} \
         ${EPOCHS} \

  # Log the submission for tracking
  echo "Submitted job: $job_name using file: $input_file"
  echo ""
done
