# Prepare data for (any) model input
# Author: Nate
# Date: Jan 7, 2025

# libraries
suppressMessages(library(pROC))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(caret))
suppressMessages(library(e1071))
suppressMessages(library(readxl))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(data.table))

# source functions
source("~/phd_code/dl-prs-v2/code/nyb_sims_code_v2/pheno.functions.R")

args = commandArgs(trailingOnly = TRUE)

# Input arguments
fam.file = args[1]                            # Path to fam file
geno.file = args[2]                           # Path to genotype .raw file
snps.file = args[3]                           # Path to causal snps RData file
stats.lin.file = args[4]                      # Path to linear summary statistics file
stats.dom.file = args[5]                      # Path to domdev summary statistics file
out.dir = args[6]                             # Path to out file

# # Load data
print("... Loading fam file ...")
fam = fread(fam.file)
print("... Loading geno file ...")
geno = fread(geno.file)
print("... Loading SNP file ...")
load(snps.file)
print("... Loading lin sumstats ...")
stats.linear = fread(stats.lin.file)
print("... Loading dom sumstats ...")
stats.domdev = fread(stats.dom.file)

# # Load inpputs manually (for script creation)
# fam.file = "~/phd_code/dl-prs-v2/data/nyb_sims_data_v2/syn_phenos/eur_nsnps100_h0.5_a0_d0.5_50k.fam"
# fam = fread("~/phd_code/dl-prs-v2/data/nyb_sims_data_v2/syn_phenos/eur_nsnps100_h0.5_a0_d0.5_50k.fam")
# geno = fread("~/phd_code/dl-prs-v2/data/nyb_sims_data_v2/genos.50k.snps.10k.raw")
# load("~/phd_code/dl-prs-v2/data/nyb_sims_data_v2/syn_phenos/nsnps100_h0.5_a0_d0.5_50k.RData")
# stats.domdev = fread("~/phd_code/dl-prs-v2/data/nyb_sims_data_v2/syn_phenos/domdev_sumstats_h0.5/domdev_gwas_eur_nsnps100_h0.5_a0_d0.5_100k.sumstats.txt")
# stats.linear = fread("~/phd_code/dl-prs-v2/data/nyb_sims_data_v2/syn_phenos/linear_sumstats_h0.5/linear_gwas_eur_nsnps100_h0.5_a0_d0.5_100k.sumstats.txt")
# out.dir = "/Users/nyb/Desktop"

print("... Preparing data ...")

# Adjust column names for fam file
colnames(fam) = c('FID','IID','father','mother','sex','phenotype')

# Define out name
out.name = sub(".*\\/([^\\/]+)\\.fam$", "\\1", fam.file)
out = paste(out.dir,"/DATA_",out.name,".txt", sep = "")
print(out)

# Extract SNP-h2 proportions
h2.a = as.numeric(str_extract(out.name, "(?<=_a)[0-9.]+"))
h2.d = as.numeric(str_extract(out.name, "(?<=_d)[0-9.]+"))

# make SNP id col for geno data
stats.linear = make.snp_allele.col(stats.linear, "ID", "A1", "ID2")
stats.domdev = make.snp_allele.col(stats.domdev, "ID", "A1", "ID2")

# Make splits for hybrid cross validation
fam = hybrid_cv_split(fam, test_frac = 0.1, k = 5, seed = 707, split_colname = "split1")
table(fam$split1)

# Set seed
set.seed(42)
n = nrow(fam)
k = 5
fold_ids1 = sample(rep(1:k, length.out = n))
folds1 = split(seq_len(n), fold_ids1)

# label splits
fam$split2 = NA
for (i in 1:length(folds1)) {
  fold_indices = folds1[[i]]
  fam$split2[fold_indices] = i
}

# ---- Compute PRSs ---- 
if ((h2.a != 0) & (h2.d != 0)) {
  print("... Computing both ADD & DOMDEV scores ...")
  
  # ---- Compute ADD PRS
  fam$additive.prs = compute_additive_prs(geno[, 7:ncol(geno)], stats.linear, c(causal.snp.list[[1]], causal.snp.list[[2]])) 
  # ---- Compute DOMDEV PRS
  fam$domdev.add.comp = compute_additive_prs(geno[, 7:ncol(geno)], stats.linear, causal.snp.list[[1]])
  fam$domdev.dom.comp = compute_domdev_prs(geno[, 7:ncol(geno)], stats.domdev, causal.snp.list[[2]], stats.linear)
  fam$domdev.prs = fam$domdev.add.comp + fam$domdev.dom.comp
  
} else if ((h2.a != 0) & (h2.d == 0)) {
  print("... Computing only ADD scores ...")
  fam$additive.prs = compute_additive_prs(geno[, 7:ncol(geno)], stats.linear, causal.snp.list[[1]]) 
  fam$domdev.dom.comp = 0
  fam$domdev.prs = 0
  
} else if ((h2.a == 0) & (h2.d != 0)) {
  print("... Computing ADD & DOMDEV scores for DOM only pheno ...")
  # ---- Compute ADD PRS
  fam$additive.prs = compute_additive_prs(geno[, 7:ncol(geno)], stats.linear, causal.snp.list[[1]]) 
  # ---- Compute DOMDEV PRS
  fam$domdev.add.comp = 0
  fam$domdev.dom.comp = compute_domdev_prs(geno[, 7:ncol(geno)], stats.domdev, causal.snp.list[[1]], stats.linear)
  fam$domdev.prs = fam$domdev.dom.comp + fam$domdev.add.comp
  
}

# Drop unnecessary columns from genos
if ((h2.a != 0) & (h2.d != 0)) {
  causal.snps = c(causal.snp.list[[1]], causal.snp.list[[2]])
  geno.ready = subset(geno, select = c(c("FID", "IID", causal.snps)))
} else {
  geno.ready = subset(geno, select = c(c("FID", "IID", causal.snp.list[[1]])))
}

# Merge with fam
data.ready = merge(fam, geno.ready, by = c("FID", "IID"))

# save
fwrite(data.ready, file = out, quote = F, sep = "\t")

print("... saved ...")

