#!/bin/bash

# Set base directories
BASE=/projects/0/vusr0599/dl-prs
CODE=${BASE}/code/ukb_doug_code
DATA=${BASE}/data/ukb_doug/ready
OUTPUT=${BASE}/output/ukb_doug
EPOCHS=100

# Iterate over all .fam files in the FAM_FILES directory
for input_file in $(ls $DATA/*DATA_int*.txt); do 
  # Extract the base name of the file (without directory and extension)
  base_name=$(basename "$input_file" .txt | sed 's/^DATA_//')
  echo "base name:"
  echo $base_name
  echo ""

  # Construct the job name by adding the prefix
  job_name="DNN_${base_name}"

  # # define the output file names + paths for test
  # geno_out_pkl="${OUTPUT}/dnn_out/GENO_${job_name}_epochs${EPOCHS}.pkl"      
  # geno_writer_dir="${OUTPUT}/dnn_model_runs/GENO_${job_name}_epochs${EPOCHS}"
  # geno_models_dir="${OUTPUT}/dnn_models/GENO_${job_name}_epochs${EPOCHS}"

  # # check if output file already exists
  # if [[ -f "$geno_out_pkl" ]]; then
  #   echo "Output or models already exist. Skipping job submission for $job_name"
  # else
  #   echo "All paths are clear. Proceeding with the script..."

  #   # Submit the job to the cluster
  #   sbatch -J "${job_name}" \
  #         -o "${BASE}/logs/ukb_doug_logs/GENO_${job_name}_%j.o.txt" \
  #         -e "${BASE}/logs/ukb_doug_logs/GENO_${job_name}_%j.e.txt" \
  #         -t 16:00:00 \
  #         "${CODE}/dnn.UKB.job.txt" \
  #         "${CODE}/03_geno_dnn_UKB.py" \
  #         "${input_file}" \
  #         "${geno_out_pkl}" \
  #         "${geno_models_dir}" \
  #         "${geno_writer_dir}" \
  #         "${EPOCHS}"

  #   # Log the submission for tracking
  #   echo "Submitted job: $job_name using file: $input_file"
  #   echo ""
  # fi

  # # define the output file names + paths for test
  # full_out_pkl="${OUTPUT}/dnn_out/FULL_${job_name}_epochs${EPOCHS}.pkl"      
  # full_writer_dir="${OUTPUT}/dnn_model_runs/FULL_${job_name}_epochs${EPOCHS}"
  # full_models_dir="${OUTPUT}/dnn_models/FULL_${job_name}_epochs${EPOCHS}"

  # # check if output file already exists
  # if [[ -f "$full_out_pkl" ]]; then
  #   echo "Output or models already exist. Skipping job submission for $job_name"
  # else
  #   echo "All paths are clear. Proceeding with the script..."

  #   # Submit the job to the cluster
  #   sbatch -J "${job_name}" \
  #         -o "${BASE}/logs/ukb_doug_logs/FULL_${job_name}_%j.o.txt" \
  #         -e "${BASE}/logs/ukb_doug_logs/FULL_${job_name}_%j.e.txt" \
  #         -t 8:00:00 \
  #         "${CODE}/dnn.UKB.job.txt" \
  #         "${CODE}/03_full_dnn_UKB.py" \
  #         "${input_file}" \
  #         "${full_out_pkl}" \
  #         "${full_models_dir}" \
  #         "${full_writer_dir}" \
  #         "${EPOCHS}"

  #   # Log the submission for tracking
  #   echo "Submitted job: $job_name using file: $input_file"
  #   echo ""
  # fi

  # define the output file names + paths for test
  covar_out_pkl="${OUTPUT}/dnn_out/COVAR_${job_name}_epochs${EPOCHS}.pkl"      
  covar_writer_dir="${OUTPUT}/dnn_model_runs/COVAR_${job_name}_epochs${EPOCHS}"
  covar_models_dir="${OUTPUT}/dnn_models/COVAR_${job_name}_epochs${EPOCHS}"

  # check if output file already exists
  if [[ -f "$covar_out_pkl" ]]; then
    echo "Output or models already exist. Skipping job submission for $job_name"
  else
    echo "All paths are clear. Proceeding with the script..."

    # Submit the job to the cluster
    sbatch -J "${job_name}" \
          -o "${BASE}/logs/ukb_doug_logs/COVAR_${job_name}_%j.o.txt" \
          -e "${BASE}/logs/ukb_doug_logs/COVAR_${job_name}_%j.e.txt" \
          -t 8:00:00 \
          "${CODE}/dnn.UKB.job.txt" \
          "${CODE}/03_covar_dnn_UKB.py" \
          "${input_file}" \
          "${covar_out_pkl}" \
          "${covar_models_dir}" \
          "${covar_writer_dir}" \
          "${EPOCHS}"

    # Log the submission for tracking
    echo "Submitted job: $job_name using file: $input_file"
    echo ""
  fi

done
