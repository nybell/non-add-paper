# Prepare data for (any) model input
# Author: Nate
# Date: Jan 7, 2025

# libraries
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(caret))
suppressMessages(library(stringr))
suppressMessages(library(data.table))

# source functions
source("/projects/0/vusr0599/dl-prs/code/pheno.functions.R")

args = commandArgs(trailingOnly = TRUE)

# Input arguments
covar.file = args[1]                      
geno.file = args[2]                        
snps.file = args[3]                    
stats.lin.file = args[4]       
stats.dom.file = args[5]          
out.dir = args[6]
pheno.name = args[7]

# # Input arguments
# covar.file = "Replication_Alanine_aminotransferase_PhenoCovar.txt"                          
# geno.file = "clump5.out.raw"                         
# snps.file = "ala_phos_clump5_snps.RData"                       
# stats.lin.file = "log_Alanine_aminotransferase_Discovery_linear_Results_CLEAN.txt"            
# stats.dom.file = "log_Alanine_aminotransferase_Discovery_genotypic_Results_CLEAN.txt"             
# out.dir = "./ready_data_out/"                        

# # Load data
print("... Loading fam file ...")
covar = fread(covar.file)
print("... Loading geno file ...")
geno = fread(geno.file)
print("... Loading SNP file ...")
load(snps.file)
print("... Loading lin sumstats ...")
stats.linear = fread(stats.lin.file)
print("... Loading dom sumstats ...")
stats.domdev = fread(stats.dom.file)

print("... Preparing data ...")

# Adjust column names for covar file
colnames(covar)[c(6,17,18)] = c('phenotype', 'log_phenotype', 'int_phenotype')

# Create new filename
out.name = paste0("DATA_", pheno.name, ".txt")
out = file.path(out.dir, out.name)

# make SNP id col for geno data
stats.linear = make.snp_allele.col(stats.linear, "ID", "A1", "ID2")
stats.domdev = make.snp_allele.col(stats.domdev, "ID", "A1", "ID2")

# Set seed
set.seed(12345)

# Generate two sets of splits with 10 folds
folds1 = createFolds(covar$int_phenotype, k = 5, list = TRUE, returnTrain = FALSE)

# label splits
covar$split1 = NA
for (i in 1:length(folds1)) {
  fold_indices = folds1[[i]]
  covar$split1[fold_indices] = i
}

# Set seed
set.seed(42)

# Generate two sets of splits with 10 folds
folds2 = createFolds(covar$int_phenotype, k = 5, list = TRUE, returnTrain = FALSE)

# label splits
covar$split2 = NA
for (i in 1:length(folds2)) {
  fold_indices = folds2[[i]]
  covar$split2[fold_indices] = i
}

# ---- Compute PRSs ---- 

print("... Computing both ADD & DOMDEV scores ...")

# ---- Compute ADD PRS
covar$additive.prs = compute_additive_prs(geno[, 7:ncol(geno)], stats.linear, c(causal.snp.list[[1]], causal.snp.list[[2]])) 
covar$addditive.dom.comp = compute_additive_prs(geno[, 7:ncol(geno)], stats.linear, causal.snp.list[[2]])
# ---- Compute DOMDEV PRS
covar$domdev.add.comp = compute_additive_prs(geno[, 7:ncol(geno)], stats.linear, causal.snp.list[[1]])
covar$domdev.dom.comp = compute_domdev_prs(geno[, 7:ncol(geno)], stats.domdev, causal.snp.list[[2]], stats.linear)
covar$domdev.prs = covar$domdev.add.comp + covar$domdev.dom.comp

# ---- Drop columns & merge ----

# drop unnecessary cols
causal.snps = c(causal.snp.list[[1]], causal.snp.list[[2]])
geno.ready = subset(geno, select = c(c("FID", "IID", causal.snps)))

# merge 
data.ready = merge(covar, geno.ready, by = c("FID", "IID"))

# save
cat("saving file to:", out, "\n")
fwrite(data.ready, file = out, quote = F, sep = "\t")

print("... saved ...")

