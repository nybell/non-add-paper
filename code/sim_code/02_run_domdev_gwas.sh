#!/bin/bash

# set fam files directory based on input arg
FAM_FILES=$1

# Set base directories
BASE_DIR=/projects/0/vusr0599/hapnest/synthetic_data
DIR_50k=/projects/0/vusr0599/hapnest/synthetic_data/data/outputs/discovery_eur_50k
DIR_100k=/projects/0/vusr0599/hapnest/synthetic_data/data/outputs/discovery_eur_100k
# FAM_FILES=/projects/0/vusr0599/hapnest/synthetic_data/data/outputs/discovery_eur_100k/fam_files

# Iterate over all .fam files in the FAM_FILES directory
for fam_file in "$FAM_FILES"/*.fam; do
  # Extract the base name of the file (without directory and extension)
  base_name=$(basename "$fam_file" .fam)

  # Construct the job name by adding the prefix
  job_name="domdev_gwas_${base_name}"

  # Check if output directory exists and has all 22 .glm.linear files
  output_dir="${DIR_100k}/${job_name}"
  if [[ -d "$output_dir" ]]; then
    linear_files_count=$(ls "$output_dir"/*.glm.linear 2>/dev/null | wc -l)
    if [[ "$linear_files_count" -eq 22 ]]; then
      echo "Skipping ${job_name} – output directory exists and has all 22 .glm.linear files."
      continue
    else
      echo "Output directory exists but has only $linear_files_count .glm.linear files – resubmitting."
    fi
  fi

  # Submit the job to the cluster
  sbatch -J "$job_name" \
         -o "$BASE_DIR/logs/${job_name}.o.txt" \
         -e "$BASE_DIR/logs/${job_name}.e.txt" \
         $BASE_DIR/syn.gwas.geno.v2.job.txt \
         "$DIR_100k" \
         eur_chr \
         "$fam_file" \
         "$job_name"

  # Log the submission for tracking
  echo "Submitted job: $job_name"
  echo "using file: $fam_file"
  echo ""
done
