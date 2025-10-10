# Simulate Continuous Phenotype with Partitioned Dominance & Linear Effects
# Author: Nate
# Batch mode version

print("# ---- STARTING PHENOTYPE SIMULATION BATCH MODE ---- #")

# Load libraries
suppressMessages(library(data.table))
suppressMessages(library(jsonlite))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

# source functions
source("/Users/nyb/phd_code/dl-prs-v2/code/nyb_sims_code_v2/pheno.functions.R")

# Set file paths
DATA = "/Users/nyb/phd_code/dl-prs-v2/data/nyb_sims_data_v2"
OUT = "/Users/nyb/phd_code/dl-prs-v2/data/nyb_sims_data_v2/syn_phenos/fams_fin_snps1000"
snp.file = "/Users/nyb/phd_code/dl-prs-v2/data/nyb_sims_data_v2/syn_phenos/fams_fin_snps1000/snps.1000.RData"
config_file = "/Users/nyb/phd_code/dl-prs-v2/data/nyb_sims_data_v2/syn_phenos/fams_fin_snps1000/sim.fin.config.snps1000.csv"

# Fixed input files
geno.file1 = file.path(DATA, "genos.100k.snps.10k.raw")
geno.file2 = file.path(DATA, "genos.50k.snps.10k.raw")
fam.file1 = file.path(DATA, "eur_merged_100k.fam")
fam.file2 = file.path(DATA, "eur_merged_50k.fam")

# Load data once
print("Loading shared data ...")
load(snp.file)
genos1 = fread(geno.file1)
genos2 = fread(geno.file2)
fam1 = fread(fam.file1)
fam2 = fread(fam.file2)
set.seed(707)

# Load simulation config
grid = fread(config_file)

# Loop through simulation grid
for (i in 1:nrow(grid)) {
  total.h2 = grid$total_h2[i]
  addit.h2 = grid$additive_h2[i]
  domin.h2 = grid$dominance_h2[i]
  n.linear.snps = grid$n_causal_snps[i]
  
  h2.sum = round(addit.h2 + domin.h2, 6)
  total.round = round(total.h2, 6)
  if (h2.sum != total.round) {
    stop(paste("ERROR at row", i, ": additive + dominance != total:", addit.h2, "+", domin.h2, "!=", total.h2))
  }
  
  # Determine simulation mode
  sim.code = if (addit.h2 > 0 & domin.h2 == 0) {
    print("Simulating: ADDITIVE ONLY")
    "ADD"
  } else if (addit.h2 == 0 & domin.h2 > 0) {
    print("Simulating: NON-ADDITIVE ONLY")
    "DEV"
  } else if (addit.h2 > 0 & domin.h2 > 0) {
    print("Simulating: MIXED SNP EFFECTS")
    "MIXED"
  } else {
    stop("Invalid h2 combination.")
  }
  
  print(paste("Running sim:", total.h2, "( add:", addit.h2, ", dom:", domin.h2, ", SNPs:", n.linear.snps, ")"))
  
  # Generate output paths
  nsnps = paste0("nsnps", n.linear.snps)
  h = paste0("h", total.h2)
  a = paste0("a", addit.h2)
  d = paste0("d", domin.h2)
  s1 = paste0(nrow(genos1)/1000, "k")
  s2 = paste0(nrow(genos2)/1000, "k")
  
  # make output directory specifically for the phenotype
  h.dir = paste0(h, "_out")
  full.h.out = file.path(OUT, h.dir)
  full.sims.h.out = file.path(full.h.out, "full_sim_results")
  # check if it exists & if not make it
  if (!dir.exists(full.h.out)) {
    print("Output directory does not exist. Creating...")
    dir.create(full.h.out, showWarnings = TRUE)
  } 
  # check if it exists & if not make it
  if (!dir.exists(full.sims.h.out)) {
    print("Sims report output directory does not exist. Creating...")
    dir.create(full.sims.h.out, showWarnings = TRUE)
  } 
  
  # set up file path and names
  fam1.base = paste("eur", nsnps, h, a, d, s1, sep = "_")
  fam1.out = file.path(full.h.out, paste0(fam1.base, ".fam"))
  fam2.base = paste("eur", nsnps, h, a, d, s2, sep = "_")
  fam2.out = file.path(full.h.out, paste0(fam2.base, ".fam"))
  log.file = file.path(full.h.out, paste("log", nsnps, h, a, d, ".txt", sep = "_"))
  sim.out = file.path(full.sims.h.out, paste0("SIM_", nsnps, "_", h, "_", a, "_", d, "_", s1, ".json"))
  snps.s2.out = file.path(full.h.out, paste0(nsnps, "_", h, "_", a, "_", d, "_", s2, ".RData"))
  
  # Check if all outputs already exist
  if (file.exists(fam1.out) & file.exists(fam2.out) & file.exists(log.file) & 
      file.exists(sim.out) & file.exists(snps.s2.out)) {
    print(paste("Skipping row", i, "— all outputs already exist."))
    next
  }
  
  # Run simulation
  print("<<>> -- GENERATING PHENOYPE -- <<>>")
  pheno.result = gen.sim.pheno(snps.1000, genos1, genos2,
                               total.h2, addit.h2, domin.h2,
                               AA = 1, Aa.add = 0.5, Aa.dev = 0.25)
  print(" <<>> -- DONE -- <<>>")

  # Save models/log
  log_results(pheno.result, sim.code)
  
  # Save SNPs
  causal.snp.list = if (sim.code == "MIXED") {
    list(pheno.result$snps$add.snps, pheno.result$snps$dev.snps)
  } else if (sim.code == "ADD") {
    list(pheno.result$snps$add.snps)
  } else {
    list(pheno.result$snps$dev.snps)
  }
  save(causal.snp.list, file = snps.s2.out)

  # Save phenotypes
  fam1$V6 = pheno.result$phenotype$phenotype.s1
  fam2$V6 = pheno.result$phenotype$phenotype.s2
  fwrite(fam1, file = fam1.out, sep = "\t", col.names = FALSE)
  fwrite(fam2, file = fam2.out, sep = "\t", col.names = FALSE)

  # Save JSON summary
  pheno.result$models = lapply(pheno.result$models, function(model_list) lapply(model_list, simplify_lm))
  save_json = toJSON(pheno.result, pretty = TRUE, auto_unbox = TRUE)
  write(save_json, file = sim.out)

  print(paste("Finished sim", i, "/", nrow(grid)))
}

print("# ---- ALL SIMULATIONS COMPLETE ---- #")





