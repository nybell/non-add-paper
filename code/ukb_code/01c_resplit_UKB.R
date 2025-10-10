# Make 70%, 15%, 15% split for UKB prepared data file
# Author: Nate
# Date: Jan 7, 2025

# libraries
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))

# setwd("~/Desktop/clump/causal_snps/ready_data_out/")
# data.file = "DATA_int_Alanine_aminotransferase.txt" 

args = commandArgs(trailingOnly = TRUE)

# Input arguments
data.file = args[1]                      

# # Load data
print("... Loading fam file ...")
dat = fread(data.file)

cat("Rows in file:", nrow(dat), "\n")
cat("Cols in file:", ncol(dat), "\n")

# Set seed
set.seed(2025)

# Create age bins (e.g., quintiles or deciles)
dat$age_bin = cut(dat$age, breaks = quantile(dat$age, probs = seq(0, 1, 0.2), na.rm = TRUE), include.lowest = TRUE)

# Create stratification variable
dat$strata = interaction(dat$sex, dat$age_bin, drop = TRUE)

# Initialize column
dat$split3 = NA

# Stratified assignment
for (s in unique(dat$strata)) {
  idx = which(dat$strata == s)
  n = length(idx)
  shuffled = sample(idx)
  
  n_train = floor(0.7 * n)
  n_tune  = floor(0.15 * n)
  n_test  = n - n_train - n_tune
  
  dat$split3[shuffled[1:n_train]] = "train"
  dat$split3[shuffled[(n_train + 1):(n_train + n_tune)]] = "tune"
  dat$split3[shuffled[(n_train + n_tune + 1):n]] = "test"
}

# Define the desired first column order
col_order = c(
  "FID", "IID", "array01", "sex", "age", "phenotype",
  "pop_pc1", "pop_pc2", "pop_pc3", "pop_pc4", "pop_pc5",
  "pop_pc6", "pop_pc7", "pop_pc8", "pop_pc9", "pop_pc10",
  "log_phenotype", "int_phenotype", 
  "split1", "split2", "split3",
  "additive.prs", "addditive.dom.comp", 
  "domdev.add.comp", "domdev.dom.comp", "domdev.prs",
  "age_bin", "strata"
)

# Filter to columns that actually exist in the data
col_order_existing = col_order[col_order %in% names(dat)]

# Find SNP columns: all columns not already in col_order
snp_cols = setdiff(names(dat), col_order)

# Reorder the data
dat = dat[, c(col_order_existing, snp_cols), with = FALSE]

# remove age bin and strata
dat = subset(dat, select = -c(age_bin, strata))

cat("Rows in file:", nrow(dat), "\n")
cat("Cols in file:", ncol(dat), "\n")

# save
cat("saving file to:", data.file, "\n")
fwrite(dat, file = data.file, quote = F, sep = "\t")

print("... saved ...")

