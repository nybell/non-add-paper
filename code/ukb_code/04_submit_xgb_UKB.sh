#!/bin/bash

# Set base directories
BASE=/projects/0/vusr0599/dl-prs
CODE=${BASE}/code/ukb_doug_code
DATA=${BASE}/data/ukb_doug/ready
OUTPUT=${BASE}/output/ukb_doug/xgb_out
MODEL_OUT=${BASE}/output/ukb_doug/xgb_models
HYPER_TRIALS=200

# Set model scripts 
GENO_MODEL="${CODE}/04_geno_xgb_UKB.py"
FULL_MODEL="${CODE}/04_full_xgb_UKB.py"
COVAR_MODEL="${CODE}/04_covar_xgb_UKB.py"

# Iterate over all .fam files in the FAM_FILES directory
for input_file in $(ls $DATA/*DATA_int*.txt); do         # *DATA_int*.txt   
  # Extract the base name of the file (without directory and extension)
  base_name=$(basename "$input_file" .txt | sed 's/^DATA_//')
  echo "base name:"
  echo $base_name
  echo ""

  # Construct the job name by adding the prefix
  job_name="XGB_${base_name}"

  echo ""
  echo " ---- Preparing to submit GENO-only job for $job_name ----"

  # define the output file & model performance dir for test
  geno_out_pkl="${OUTPUT}/GENO_${job_name}_trials${HYPER_TRIALS}.pkl"
  geno_models_dir="${MODEL_OUT}/GENO_${job_name}_trials${HYPER_TRIALS}.json"

  # check if output file already exists
  if [[ -f "$geno_out_pkl" || -d "$geno_models_dir" ]]; then
    echo "Output or models already exist. Skipping job submission for $job_name"
  else
    echo "All paths are clear. Proceeding with the script..."

    # Submit the job to the cluster
    sbatch -J "GENO_${job_name}" \
          -o "$BASE/logs/ukb_doug_logs/${job_name}_GENO_%j.o.txt" \
          -e "$BASE/logs/ukb_doug_logs/${job_name}_GENO_%j.e.txt" \
          -t 10:00:00 \
          "$CODE/xgb.UKB.job.txt" \
          "${GENO_MODEL}" \
          "${input_file}" \
          "${geno_out_pkl}" \
          "${geno_models_dir}" \
          "${HYPER_TRIALS}"

    # Log the submission for tracking
    echo "Submitted GENO job: $job_name using file: $input_file"
    echo "---------------------------"
    echo ""
  fi

  echo ""
  echo " ---- Preparing to submit COVAR-only job for $job_name ----"

  # define the output file & model performance dir for test
  covar_out_pkl="${OUTPUT}/COVAR_${job_name}_trials${HYPER_TRIALS}.pkl"
  covar_models_dir="${MODEL_OUT}/COVAR_${job_name}_trials${HYPER_TRIALS}.json"

  # check if output file already exists
  if [[ -f "$covar_out_pkl" || -d "$covar_models_dir" ]]; then
    echo "Output or models already exist. Skipping job submission for $job_name"
  else
    echo "All paths are clear. Proceeding with the script..."
    # Submit the job to the cluster
    sbatch -J "COVAR_${job_name}" \
          -o "$BASE/logs/ukb_doug_logs/${job_name}_COVAR_%j.o.txt" \
          -e "$BASE/logs/ukb_doug_logs/${job_name}_COVAR_%j.e.txt" \
          -t 1:00:00 \
          "${CODE}/xgb.UKB.job.txt" \
          "${COVAR_MODEL}" \
          "${input_file}" \
          "${covar_out_pkl}" \
          "${covar_models_dir}" \
          "${HYPER_TRIALS}"

    # Log the submission for tracking
    echo "Submitted COVAR job: $job_name using file: $input_file"
    echo "---------------------------"
    echo ""
  fi

  echo ""
  echo " ---- Preparing to submit FULL model job for $job_name ----"

  # define the output file & model performance dir for test
  full_out_pkl="${OUTPUT}/FULL_${job_name}_trials${HYPER_TRIALS}.pkl"
  full_models_dir="${MODEL_OUT}/FULL_${job_name}_trials${HYPER_TRIALS}.json"

  # check if output file already exists
  if [[ -f "$full_out_pkl" || -d "$full_models_dir" ]]; then
    echo "Output or models already exist. Skipping job submission for $job_name"
  else
    echo "All paths are clear. Proceeding with the script..."

    # Submit the job to the cluster
    sbatch -J "FULL_${job_name}" \
          -o "$BASE/logs/ukb_doug_logs/${job_name}_FULL_%j.o.txt" \
          -e "$BASE/logs/ukb_doug_logs/${job_name}_FULL_%j.e.txt" \
            -t 10:00:00 \
          "${CODE}/xgb.UKB.job.txt" \
          "${FULL_MODEL}" \
          "${input_file}" \
          "${full_out_pkl}" \
          "${full_models_dir}" \
          "${HYPER_TRIALS}" 

    # Log the submission for tracking
    echo "Submitted FULL model job: $job_name using file: $input_file"
    echo "---------------------------"
    echo ""
  fi

  
done
