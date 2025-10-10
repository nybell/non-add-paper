#!/bin/bash

# set paths
LOGS="/projects/0/vusr0599/dl-prs/logs/ukb_doug_logs"
CODE="/projects/0/vusr0599/dl-prs/code/ukb_doug_code"
DATA_DIR="/projects/0/vusr0599/dl-prs/data/ukb_doug/results/clean"

# Loop over all matching files
for filepath in "$DATA_DIR"/int_*_Discovery_linear_Results_CLEAN.txt; do
    # Extract the base filename
    filename=$(basename "$filepath")
    
    # Extract the phenotype name before '_Discovery_'
    phenotype="${filename%%_Discovery_*}"

    # set job name
    job_name="CLEAN_${phenotype}"

    # Set summary statistic out paths
    add_stats="${DATA_DIR}/${phenotype}_Discovery_linear_Results_CLEAN.txt"
    dev_stats="${DATA_DIR}/${phenotype}_Discovery_genotypic_Results_CLEAN.txt"

    echo "Additive stats: $add_stats"
    echo "Dominant/Deviant stats: $dev_stats"
    echo "Phenotype: $phenotype"
    echo "Job name: $job_name"
    echo ""
    
    # Submit the job to the cluster
    sbatch -J "${job_name}" \
            -o "${LOGS}/${job_name}_%j.o.txt" \
            -e "${LOGS}/${job_name}_%j.e.txt" \
            $CODE/01_clump_clean_export.job.txt \
            ${add_stats} \
            ${dev_stats} \
            ${phenotype} 
    
    sleep 1
    
done