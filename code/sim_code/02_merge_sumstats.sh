#!/bin/bash

in_path=$1
prefix=$2
out_file=$3

# Loop through chromosomes and check for missing/empty files
for chr in {1..22}; do
  file="$in_path/${prefix}-${chr}.PHENO1.glm.linear"
  if [ -s "$file" ]; then
    ((valid_file_count++))
  else
    missing_or_empty+=("$file")
  fi
done

# Check the result
if [ "$valid_file_count" -ne 22 ]; then
  echo "Error: Expected 22 non-empty files, but found $valid_file_count."
  echo "Missing or empty files:"
  for f in "${missing_or_empty[@]}"; do
    echo " - $f"
  done
  exit 1
else
  echo "All 22 GWAS result files are present and non-empty. Proceeding with concatenation..."
fi

# First, extract the header from the first file
head -n 1 $in_path/eur_chr-1.PHENO1.glm.linear > $out_file

# Then concatenate the rest of the files, skipping the first line (header)
for i in $in_path/$prefix-{1..22}.PHENO1.glm.linear; do
  if [ -e "$i" ]; then
    # File exists, concatenate skipping the first line (header)
    tail -n +2 "$i" >> $out_file
  else
    # File does not exist, skip to the next iteration
    echo "File $i does not exist, skipping..."
    continue
  fi
done

# example usage
# ./merge_sumstats.sh $100k/poly001_no_int_gwas eur_chr $100k/poly001_no_int_gwas/100k.no.int.sumstats.txt
