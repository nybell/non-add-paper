# File to extract "causal" SNPs (i.e., SNPs for model input)
# Author: Nate YB
# Date: 1 May 2025

# libraries
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))

# source("~/phd_code/dl-prs-v2/code/nyb_sims_code_v2/pheno.functions.R")

args = commandArgs(trailingOnly = TRUE)

# input arguments
clump.add.file = args[1]
clump.dev.file = args[2]
stats.add.file = args[3]
stats.dev.file = args[4]
base.name = args[5]
out.dir = args[6]

# out put directory 
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

# make output file paths
snps.export.file = file.path(out.dir, paste0(base.name, ".export.txt"))
snps.causal.rdata = file.path(out.dir, paste0(base.name, "_snps.RData"))

# load files
print("loading files.")
clump.add = fread(clump.add.file)
clump.dev = fread(clump.dev.file)
stats.add = fread(stats.add.file)
stats.dev = fread(stats.dev.file)     # dev stats with only "domdev" rows in TEST col
print("finished loading files.")

# Format the column ID to make the SNP columns needed for computing PRS
make.snp_allele.col = function(data, col1, col2, new_col_name) {
  data[[new_col_name]] = paste(data[[col1]], data[[col2]], sep = "_")
  return(data)
}

# Function to extract SP2 map
extract_sp2_map = function(clump_df) {
  sp2_map = list()
  for (i in seq_len(nrow(clump_df))) {
    lead_snp = clump_df$SNP[i]
    sp2_raw = clump_df$SP2[i]
    if (is.na(sp2_raw)) {
      sp2_map[[lead_snp]] = character(0)
    } else {
      sp2_clean = gsub("\\(1\\)", "", sp2_raw)
      sp2_split = unlist(strsplit(sp2_clean, ","))
      sp2_map[[lead_snp]] = sp2_split
    }
  }
  return(sp2_map)
}

# Function to find which dominance lead SNP contains a given SNP in its SP2 list
find_lead_snp_for_sp2 = function(snp_id, sp2_map) {
  for (lead in names(sp2_map)) {
    if (snp_id %in% sp2_map[[lead]]) {
      return(lead)
    }
  }
  return(NA)
}

# Master function to compile PRS SNPs
compile_prs_snps = function(clump_lin, clump_dom) {
  linear_snps = clump_lin$SNP
  dominance_snps = clump_dom$SNP
  geno_sp2_map = extract_sp2_map(clump_dom)
  
  final_snps = character()
  designation = character()
  replaced_from_linear = character()
  
  for (snp in linear_snps) {
    if (snp %in% dominance_snps) {                        # check if linear SNP is in dominance SNPs column - if so, assign as dominance
      final_snps = c(final_snps, snp)
      designation = c(designation, "domdev")
      replaced_from_linear = c(replaced_from_linear, NA)
    } else {
      replacement = find_lead_snp_for_sp2(snp, geno_sp2_map)
      if (!is.na(replacement)) {                          # check if linear SNP is in any of the dominance SNP clumps, if so, replace and assign dominance
        final_snps = c(final_snps, replacement)
        designation = c(designation, "domdev")
        replaced_from_linear = c(replaced_from_linear, snp)
      } else {                                            # if it's not a significant dominance SNP, assign as additive
        final_snps = c(final_snps, snp)
        designation = c(designation, "additive")
        replaced_from_linear = c(replaced_from_linear, NA)
      }
    }
  }
  
  # Add dominance SNPs not yet included
  additional_dom_snps = setdiff(dominance_snps, final_snps)     # add rest of the dominance SNPs (if they are not present in linear SNP list)
  final_snps = c(final_snps, additional_dom_snps)
  designation = c(designation, rep("domdev", length(additional_dom_snps)))
  replaced_from_linear = c(replaced_from_linear, rep(NA, length(additional_dom_snps)))
  
  # Build final data.frame
  prs_snps_df = data.frame(
    SNP = final_snps,
    designation = designation,
    replaced_linear_snp = replaced_from_linear,
    stringsAsFactors = FALSE
  )
  
  prs_snps_df = prs_snps_df[!duplicated(prs_snps_df$SNP), ]
  return(prs_snps_df)
}

# Helper function to annotate a PRS data.frame with A1 column
annotate_prs_with_a1 = function(prs_df, stats_lin, stats_dom) {
  # Merge additive SNPs
  additive_snps = subset(prs_df, designation == "additive")
  domdev_snps = subset(prs_df, designation == "domdev")
  
  additive_annot = merge(additive_snps, stats_lin[, .(ID, A1)], by.x = "SNP", by.y = "ID", all.x = TRUE)
  domdev_annot = merge(domdev_snps, stats_dom[, .(ID, A1)], by.x = "SNP", by.y = "ID", all.x = TRUE)
  
  # Combine and preserve original order
  combined = rbind(additive_annot, domdev_annot)
  combined = combined[match(prs_df$SNP, combined$SNP), ]  # restore original order
  
  # Rename A1 column for clarity if desired
  # setnames(combined, "A1", "A1_from_correct_source")
  
  return(combined)
}

# Compile PRS SNPs
print("Compiling list of SNP for PRSs.")
prs.snps = compile_prs_snps(clump.add, clump.dev)
print("Compiling finished.")

# View summary
cat("PRS SNPs:", nrow(prs.snps), "\n")

# Annotate both clumping sets
print("Annotating SNPs with summary statistics.")
prs.snps.annot = annotate_prs_with_a1(prs.snps, stats.add, stats.dev)
print("Annotating finished.")

# make snp ID2 clump
print("Adding SNP ID2 column.")
prs.snps.annot = make.snp_allele.col(prs.snps.annot, "SNP", "A1", "ID2")
print("ID2 column added.")

cat("before removing SNPs with missing A1 alleles:", nrow(prs.snps.annot), "SNPs remaining.\n")
prs.snps.annot = prs.snps.annot[!is.na(prs.snps.annot$A1), ]
cat("after removing SNPs with missing A1 alleles:", nrow(prs.snps.annot), "SNPs remaining.\n")

# ---- format & save snp ---- 

# export file for plink
snps.export = subset(prs.snps.annot, select = c(SNP, A1))
fwrite(snps.export, file = snps.export.file, quote = F, sep = "\t", col.names = F)
# additive SNPs
clump.add.snps = subset(prs.snps.annot, designation == "additive")
clump.add.snps = as.vector(clump.add.snps$ID2)
# domdev SNPs
clump.dev.snps = subset(prs.snps.annot, designation == "domdev")
clump.dev.snps = as.vector(clump.dev.snps$ID2)
cat("Additive SNPs:", length(clump.add.snps), "\n")
cat("Dominance SNPs:", length(clump.dev.snps), "\n")
# causal SNP list
causal.snp.list = list(clump.add.snps, clump.dev.snps)
save(causal.snp.list, file = snps.causal.rdata)



