#!/bin/bash

# Directory where GWAS directories are located
GWAS_DIR=/projects/0/vusr0599/hapnest/synthetic_data/data/outputs/discovery_eur_100k
OUT_DIR=/projects/0/vusr0599/hapnest/synthetic_data/data/outputs/discovery_eur_100k/merged_devsim_sumstats

# Change to the gwas directory
cd $GWAS_DIR || exit

# Loop through directories containing "_gwas_eur_" in their names
for DIR in *gwas_eur_devsim*; do
    # Ensure it's a directory
    if [ -d "$DIR" ]; then
        # Extract the directory name
        DIR_NAME=$(basename "$DIR")

        # Modify the name to the desired format
        MERGED_NAME=$(echo "$DIR_NAME" | sed 's/^gwas_//' | sed 's/$/.sumstats.txt/')

        # Check if the expected output file already exists
        if [[ -f "${OUT_DIR}/${MERGED_NAME}" ]]; then                           # ${GWAS_DIR}/${DIR}/${MERGED_NAME}
            echo "Skipping ${MERGED_NAME} – output file already exists."
            continue
        fi

        # Run the merge script
        /projects/0/vusr0599/hapnest/synthetic_data/02_merge_sumstats.sh \
        "$DIR" \
        eur_chr \
        "${OUT_DIR}/${MERGED_NAME}"

    fi
done
