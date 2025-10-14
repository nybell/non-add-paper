# Phenotype simulation analysis functions
# Author: Nate Jeetz
# Date: Dec 14 2024

# Load libraries
library(data.table)
library(dplyr)
library(tidyr)

# ---- DEFINE FUNCTIONS ---- 

# format the column ID to be able to extract causal SNPs - separates geno SNPs into two columns
make.col.id = function(genos, snp.freqs){
  # merge .raw SNP column names with freq file for sims
  genos.snp.ids = data.frame(
    col.id = colnames(genos)[7:ncol(genos)]) %>%
    mutate(tmp.col = col.id) %>%
    separate(tmp.col, into = c("ID", "ALT"), sep = "_")
  # merge with the afreq file
  snp.freqs = merge(snp.freqs, genos.snp.ids, by = c("ID", "ALT"))
  return(snp.freqs)
}

# Format the column ID to make the SNP columns needed for computing PRS
make.snp_allele.col = function(data, col1, col2, new_col_name) {
  data[[new_col_name]] = paste(data[[col1]], data[[col2]], sep = "_")
  return(data)
}

# sample causal snps
sample.causal.snps <- function(snp.freqs, n.snps, maf.dist, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  if (!maf.dist) {
    # Uniform random sampling
    return(sample(snp.freqs$col.id, n.snps, replace = FALSE))
  }
  
  # Ensure Bin column exists
  if (!"Bin" %in% names(snp.freqs)) stop("snp.freqs must have a 'Bin' column.")
  bins <- sort(unique(snp.freqs$Bin))
  B <- length(bins)
  
  # Count available per bin
  bin_counts <- vapply(bins, function(b) sum(snp.freqs$Bin == b), integer(1))
  total_available <- sum(bin_counts)
  if (total_available < n.snps) {
    stop(sprintf("Not enough SNPs available across bins: need %d, have %d.", 
                 n.snps, total_available))
  }
  
  # Base even allocation
  base <- n.snps %/% B
  target <- pmin(base, bin_counts)  # cap by availability
  remaining <- n.snps - sum(target)
  
  # Evenly distribute the remainder, 1 at a time to bins with capacity
  # This keeps the distribution as even as possible and avoids ceiling overshoot
  while (remaining > 0) {
    capacity <- bin_counts - target
    idx <- which(capacity > 0)
    if (length(idx) == 0) stop("Allocation error: no capacity left but remaining > 0.")
    # Randomize which bins get the extra draws to avoid bias
    idx <- sample(idx)
    take <- min(remaining, length(idx))
    target[idx[seq_len(take)]] <- target[idx[seq_len(take)]] + 1
    remaining <- remaining - take
  }
  
  # Now sample per bin
  sel <- unlist(mapply(function(b, k) {
    if (k == 0) return(character(0))
    pool <- snp.freqs$col.id[snp.freqs$Bin == b]
    sample(pool, size = k, replace = FALSE)
  }, b = bins, k = target, SIMPLIFY = FALSE), use.names = FALSE)
  
  # Safety check
  if (length(sel) != n.snps) {
    stop(sprintf("Selected %d SNPs, expected %d.", length(sel), n.snps))
  }
  return(sel)
}

# model for simplfying lm results into table
simplify_lm <- function(model) {
  # Check if the model is NULL or not an lm object
  if (is.null(model) || !inherits(model, "lm")) {
    return(NULL)
  }
  
  # Simplify the lm object
  list(
    coefficients = coef(model),
    r_squared = summary(model)$r.squared,
    adj_r_squared = summary(model)$adj.r.squared
  )
}

# phenotype sim function
gen.sim.pheno = function(snps, genos1, genos2, total.h2, addit.h2, domin.h2, AA = 1, Aa.add = 0.5, Aa.dev = 0.25) {
  # get number of causal SNPs
  n.snps = length(snps)
  # Combine genotype data
  genos = rbind(genos1, genos2)
  # Calculate number of causal SNPs for each category
  n.add.snps = round((addit.h2 / total.h2) * n.snps)
  n.dev.snps = round((domin.h2 / total.h2) * n.snps)
  # get scaled additive effects if non.add.snps != 0
  if (n.add.snps != 0) {
    # Randomly sample additive causal SNPs evenly distributed across MAF bins
    causal.snps.add = sample(snps, n.add.snps)
    
    # Generate additive effects
    genos.add = genos %>% select(all_of(causal.snps.add))
    effects.mat.add = apply(genos.add, 2, function(genotype) {
      ifelse(genotype == 0, 0, ifelse(genotype == 1, Aa.add, AA))
    })
    
    # get SNP variances
    var.add = apply(effects.mat.add, 2, var)
    sum.var.add = sum(var.add)
    scaling.add = sqrt(addit.h2 / sum(var.add))
    effects.scaled.add = effects.mat.add * scaling.add
    sum.var.add.scaled = sum(apply(effects.scaled.add, 2, var))
    
  } else {
    sum.var.add.scaled = 0
  }
  
  # get scaled non-additive effects if non.add.snps != 0
  if (n.dev.snps != 0) {
    
    if (n.add.snps != 0) {
      # Get remaining SNPs that are not in the causal additive SNP list
      causal.snps.dev = setdiff(snps, causal.snps.add)
    } else if (n.add.snps == 0) {
      # sample only dev SNPs
      causal.snps.dev = sample(snps, n.dev.snps)
    }
    
    # Generate additive effects
    genos.dev = genos %>% select(all_of(causal.snps.dev))
    effects.mat.dev = apply(genos.dev, 2, function(genotype) {
      ifelse(genotype == 0, 0, ifelse(genotype == 1, Aa.dev, AA))
    })
    
    # get SNP variances
    var.dev = apply(effects.mat.dev, 2, var)
    sum.var.dev = sum(var.dev)
    scaling.dev = sqrt(domin.h2 / sum(var.dev))
    effects.scaled.dev = effects.mat.dev * scaling.dev
    sum.var.dev.scaled = sum(apply(effects.scaled.dev, 2, var))
    
  } else {
    sum.var.dev.scaled = 0
  }
  
  # Total variance
  sum.var.combined = sum.var.add.scaled + sum.var.dev.scaled
  total.pheno.var = sum.var.combined / total.h2
  
  # Noise variance
  var.noise = total.pheno.var - sum.var.combined
  sd.noise = sqrt(var.noise)
  
  # gen pheno and return
  if ((n.add.snps != 0) & (n.dev.snps != 0)) {
    # gen pheno
    phenotype = rowSums(effects.scaled.add) + rowSums(effects.scaled.dev) + rnorm(nrow(genos), mean = 0, sd = sd.noise)
    
    # Split combined data into samples
    phenotype.s1 = phenotype[1:nrow(genos1)]
    phenotype.s2 = phenotype[(nrow(genos1) + 1):nrow(genos)]
    
    effects.add.s1 = rowSums(effects.scaled.add)[1:nrow(genos1)]
    effects.add.s2 = rowSums(effects.scaled.add)[(nrow(genos1) + 1):nrow(genos)]
    
    effects.dev.s1 = rowSums(effects.scaled.dev)[1:nrow(genos1)]
    effects.dev.s2 = rowSums(effects.scaled.dev)[(nrow(genos1) + 1):nrow(genos)]
    
    # Run models for sample 1
    s1.model.mixed = lm(phenotype.s1 ~ effects.add.s1 + effects.dev.s1)
    s1.model.add = lm(phenotype.s1 ~ effects.add.s1)
    s1.model.dev = lm(phenotype.s1 ~ effects.dev.s1)
    s1.models = list(mixed = s1.model.mixed, add = s1.model.add, dev = s1.model.dev)
    
    # Run models for sample 2
    s2.model.mixed = lm(phenotype.s2 ~ effects.add.s2 + effects.dev.s2)
    s2.model.add = lm(phenotype.s2 ~ effects.add.s2)
    s2.model.dev = lm(phenotype.s2 ~ effects.dev.s2)
    s2.models = list(mixed = s2.model.mixed, add = s2.model.add, dev = s2.model.dev)
    
    # run models for combined phenotype
    combined.model.mixed = lm(phenotype ~ rowSums(effects.scaled.add) + rowSums(effects.scaled.dev))
    combined.model.add = lm(phenotype ~ rowSums(effects.scaled.add))
    combined.model.dev = lm(phenotype ~ rowSums(effects.scaled.dev))
    combined.models = list(mixed = combined.model.mixed, add = combined.model.add, dev = combined.model.dev)
    
    # make pheno list
    pheno.list = list(phenotype.s1 = phenotype.s1, phenotype.s2 = phenotype.s2)
    # effects list
    effect.list = list(effects.add = effects.scaled.add, effects.dev = effects.scaled.dev)
    # make causal SNPs list
    causal.snps.list = list(add.snps = causal.snps.add, dev.snps = causal.snps.dev)
    # models list 
    models.list = list(combined.models = combined.models, s1.models = s1.models, s2.models = s2.models)
    # make results list
    result = list(phenotype = pheno.list, snps = causal.snps.list, models = models.list) # effects = effect.list
    # return
    return(result)
    
  } else if ((n.add.snps != 0) & (n.dev.snps == 0)) {
    # gen pheno
    phenotype = rowSums(effects.scaled.add) + rnorm(nrow(genos), mean = 0, sd = sd.noise)
    
    # Split combined data into samples
    phenotype.s1 = phenotype[1:nrow(genos1)]
    phenotype.s2 = phenotype[(nrow(genos1) + 1):nrow(genos)]
    
    effects.add.s1 = rowSums(effects.scaled.add)[1:nrow(genos1)]
    effects.add.s2 = rowSums(effects.scaled.add)[(nrow(genos1) + 1):nrow(genos)]
    
    # Run models for sample 1
    s1.model.add = lm(phenotype.s1 ~ effects.add.s1)
    s1.models = list(mixed = NULL, add = s1.model.add, dev = NULL)
    
    # Run models for sample 2
    s2.model.add = lm(phenotype.s2 ~ effects.add.s2)
    s2.models = list(mixed = NULL, add = s2.model.add, dev = NULL)
    
    # run models for combined phenotype
    combined.model.add = lm(phenotype ~ rowSums(effects.scaled.add))
    combined.models = list(mixed = NULL, add = combined.model.add, dev = NULL)
    
    # make pheno list
    pheno.list = list(phenotype.s1 = phenotype.s1, phenotype.s2 = phenotype.s2)
    # effects list
    effect.list = list(effects.add = effects.scaled.add, effects.dev = NULL)
    # make causal SNPs list
    causal.snps.list = list(add.snps = causal.snps.add, dev.snps = NULL)
    # models list 
    models.list = list(combined.models = combined.models, s1.models = s1.models, s2.models = s2.models)
    # make results list
    result = list(phenotype = pheno.list, snps = causal.snps.list, models = models.list) # effects = effect.list
    # return
    return(result)
    
  } else if ((n.add.snps == 0) & (n.dev.snps != 0)) {
    # gen pheno
    phenotype = rowSums(effects.scaled.dev) + rnorm(nrow(genos), mean = 0, sd = sd.noise)
    
    # Split combined data into samples
    phenotype.s1 = phenotype[1:nrow(genos1)]
    phenotype.s2 = phenotype[(nrow(genos1) + 1):nrow(genos)]
    
    effects.dev.s1 = rowSums(effects.scaled.dev)[1:nrow(genos1)]
    effects.dev.s2 = rowSums(effects.scaled.dev)[(nrow(genos1) + 1):nrow(genos)]
    
    # Run models for sample 1
    s1.model.dev = lm(phenotype.s1 ~ effects.dev.s1)
    s1.models = list(mixed = NULL, add = NULL, dev = s1.model.dev)
    
    # Run models for sample 2
    s2.model.dev = lm(phenotype.s2 ~ effects.dev.s2)
    s2.models = list(mixed = NULL, add = NULL, dev = s2.model.dev)
    
    # run models for combined phenotype
    combined.model.dev = lm(phenotype ~ rowSums(effects.scaled.dev))
    combined.models = list(mixed = NULL, add = NULL, dev = combined.model.dev)
    
    # make pheno list
    pheno.list = list(phenotype.s1 = phenotype.s1, phenotype.s2 = phenotype.s2)
    # effects list
    effect.list = list(effects.add = NULL, effects.dev = effects.scaled.dev)
    # make causal SNPs list
    causal.snps.list = list(add.snps = NULL, dev.snps = causal.snps.dev)
    # models list 
    models.list = list(combined.models = combined.models, s1.models = s1.models, s2.models = s2.models)
    # make results list
    result = list(phenotype = pheno.list, snps = causal.snps.list, models = models.list) # effects = effect.list
    # return
    return(result)
    
  }
  
}

# phenotype sim function
# pheno = gen.rand.sim.pheno(n.snps, snp.freqs, genos1, genos2, total.h2, addit.h2, domin.h2)
gen.rand.sim.pheno = function(n.linear.snps, snp.freqs, genos1, genos2, total.h2, addit.h2, domin.h2) {
  # Combine genotype data
  genos = rbind(genos1, genos2)
  # Calculate number of causal SNPs for each category
  n.add.snps = round((addit.h2 / total.h2) * n.linear.snps)
  n.dev.snps = round((domin.h2 / total.h2) * n.linear.snps)
  print(n.add.snps)
  print(n.dev.snps)
  # get scaled additive effects if non.add.snps != 0
  if (n.add.snps != 0) {
    # Randomly sample additive causal SNPs evenly distributed across MAF bins
    causal.snps.add = sample.causal.snps(snp.freqs, n.add.snps, maf.dist)
    print("snp sampling done.")
    # Get genotypes for causal SNPs
    genos.add = genos %>% select(all_of(causal.snps.add))
    print("genotype extracting done.")
    # Generate distribution of additive effects (effect of gene substitution per Falconer, pg. 116)
    alphas.add = rnorm(n.add.snps, mean = 0, sd = 1.0)
    print("add alpha generating done.")
    # Weight genotypes by alpha values (effect size)
    # G = as.matrix(genos.add)
    # A = diag(alphas.add)
    # print("matrix converting done.")
    # effects.mat.add = G %*% A
    effects.mat.add = sweep(genos.add, 2, alphas.add, `*`)
    print("additive effects done")
    # get SNP variances
    var.add = apply(effects.mat.add, 2, var)
    sum.var.add = sum(var.add)
    scaling.add = sqrt(addit.h2 / sum(var.add))
    effects.scaled.add = effects.mat.add * scaling.add
    sum.var.add.scaled = sum(apply(effects.scaled.add, 2, var))
    
  } else {
    sum.var.add.scaled = 0
  }
  
  # get scaled non-additive effects if non.add.snps != 0
  if (n.dev.snps != 0) {
    
    if (n.add.snps != 0) {
      # subset the snp.freqs data to get those that are not additive causal SNPs
      snp.freqs.nc = snp.freqs %>%
        filter(!col.id %in% causal.snps.add)
      # Random sample non-additive causal SNPs evenly distributed across MAF bins
      causal.snps.dev = sample.causal.snps(snp.freqs.nc, n.dev.snps, maf.dist)
    } else if (n.add.snps == 0) {
      # sample only dev SNPs
      causal.snps.dev = sample.causal.snps(snp.freqs, n.dev.snps, maf.dist)
    }
    
    # Get genotypes for causal SNPs
    genos.dev = genos %>% select(all_of(causal.snps.dev))
    # Generate distribution of additive effects (effect of gene substitution per Falconer, pg. 116)
    alphas.dev = rnorm(n.dev.snps, mean = 0, sd = 1.0)
    # Generate distribution of dominance deviations: 0 < d < |alpha|
    deltas.dev = runif(
      n   = length(alphas.dev),
      min = 0,
      max = abs(alphas.dev)
    )
    
    # Initialize an empty matrix to store the effects
    effects.mat.dev = matrix(0, nrow = nrow(genos.dev), ncol = ncol(genos.dev))
    # Apply the logic column by column
    for (i in seq_len(ncol(genos.dev))) {
      geno_col = genos.dev[[i]]
      alpha = alphas.dev[i]
      delta = deltas.dev[i]
      
      effects.mat.dev[, i] = ifelse(
        geno_col == 1, alpha + delta,
        ifelse(geno_col == 2, 2 * alpha, 0)
      )
    }
    
    # get SNP variances
    var.dev = apply(effects.mat.dev, 2, var)
    sum.var.dev = sum(var.dev)
    scaling.dev = sqrt(domin.h2 / sum(var.dev))
    effects.scaled.dev = effects.mat.dev * scaling.dev
    sum.var.dev.scaled = sum(apply(effects.scaled.dev, 2, var))
    
  } else {
    sum.var.dev.scaled = 0
  }
  
  # # Genetic value per individual
  # g.add = rowSums(effects.scaled.add)
  # g.dom = rowSums(effects.scaled.dev)
  # g.tot = g.add + g.dom   # total genetic component
  # 
  # # Add environmental noise so that total h2 ≈ 0.10
  # pheno.var = var(g.tot)
  # env.var  = pheno.var * (1 - total.h2) / total.h2
  # noise = rnorm(nrow(genos), mean = 0, sd = sqrt(env.var))
  # 
  # phenotype = g.tot + noise
  # summary(lm(phenotype ~ g.tot))
  # summary(lm(phenotype ~ g.add))
  # summary(lm(phenotype ~ g.dom))
  
  # gen pheno and return
  if ((n.add.snps != 0) & (n.dev.snps != 0)) {
    # Genetic value per individual
    g.add = rowSums(effects.scaled.add)
    g.dom = rowSums(effects.scaled.dev)
    g.tot = g.add + g.dom   # total genetic component
    
    # Add environmental noise so that total h2 ≈ 0.10
    pheno.var = var(g.tot)
    env.var  = pheno.var * (1 - total.h2) / total.h2
    noise = rnorm(nrow(genos), mean = 0, sd = sqrt(env.var))
    
    phenotype = g.tot + noise
    # summary(lm(phenotype ~ g.tot))
    # summary(lm(phenotype ~ g.add))
    # summary(lm(phenotype ~ g.dom))
    
    # Split combined data into samples
    phenotype.s1 = phenotype[1:nrow(genos1)]
    phenotype.s2 = phenotype[(nrow(genos1) + 1):nrow(genos)]
    
    effects.add.s1 = rowSums(effects.scaled.add)[1:nrow(genos1)]
    effects.add.s2 = rowSums(effects.scaled.add)[(nrow(genos1) + 1):nrow(genos)]
    
    effects.dev.s1 = rowSums(effects.scaled.dev)[1:nrow(genos1)]
    effects.dev.s2 = rowSums(effects.scaled.dev)[(nrow(genos1) + 1):nrow(genos)]
    
    # Run models for sample 1
    s1.model.mixed = lm(phenotype.s1 ~ effects.add.s1 + effects.dev.s1)
    s1.model.add = lm(phenotype.s1 ~ effects.add.s1)
    s1.model.dev = lm(phenotype.s1 ~ effects.dev.s1)
    s1.models = list(mixed = s1.model.mixed, add = s1.model.add, dev = s1.model.dev)
    
    # Run models for sample 2
    s2.model.mixed = lm(phenotype.s2 ~ effects.add.s2 + effects.dev.s2)
    s2.model.add = lm(phenotype.s2 ~ effects.add.s2)
    s2.model.dev = lm(phenotype.s2 ~ effects.dev.s2)
    s2.models = list(mixed = s2.model.mixed, add = s2.model.add, dev = s2.model.dev)
    
    # run models for combined phenotype
    combined.model.mixed = lm(phenotype ~ rowSums(effects.scaled.add) + rowSums(effects.scaled.dev))
    combined.model.add = lm(phenotype ~ rowSums(effects.scaled.add))
    combined.model.dev = lm(phenotype ~ rowSums(effects.scaled.dev))
    combined.models = list(mixed = combined.model.mixed, add = combined.model.add, dev = combined.model.dev)
    
    # make pheno list
    pheno.list = list(phenotype.s1 = phenotype.s1, phenotype.s2 = phenotype.s2)
    # effects list
    effect.list = list(effects.add = effects.scaled.add, effects.dev = effects.scaled.dev)
    # make causal SNPs list
    causal.snps.list = list(add.snps = causal.snps.add, dev.snps = causal.snps.dev)
    # models list 
    models.list = list(combined.models = combined.models, s1.models = s1.models, s2.models = s2.models)
    # make results list
    result = list(phenotype = pheno.list, snps = causal.snps.list, models = models.list, effects = effect.list) # effects = effect.list
    # return
    return(result)
    
  } else if ((n.add.snps != 0) & (n.dev.snps == 0)) {
    # Genetic value per individual
    g.add = rowSums(effects.scaled.add)
    g.tot = g.add   # total genetic component
    
    # Add environmental noise so that total h2 ≈ 0.10
    pheno.var = var(g.tot)
    env.var  = pheno.var * (1 - total.h2) / total.h2
    noise = rnorm(nrow(genos), mean = 0, sd = sqrt(env.var))
    
    # Generate phenotype
    phenotype = g.tot + noise
    
    # Split combined data into samples
    phenotype.s1 = phenotype[1:nrow(genos1)]
    phenotype.s2 = phenotype[(nrow(genos1) + 1):nrow(genos)]
    
    effects.add.s1 = rowSums(effects.scaled.add)[1:nrow(genos1)]
    effects.add.s2 = rowSums(effects.scaled.add)[(nrow(genos1) + 1):nrow(genos)]
    
    # Run models for sample 1
    s1.model.add = lm(phenotype.s1 ~ effects.add.s1)
    s1.models = list(mixed = NULL, add = s1.model.add, dev = NULL)
    
    # Run models for sample 2
    s2.model.add = lm(phenotype.s2 ~ effects.add.s2)
    s2.models = list(mixed = NULL, add = s2.model.add, dev = NULL)
    
    # run models for combined phenotype
    combined.model.add = lm(phenotype ~ rowSums(effects.scaled.add))
    combined.models = list(mixed = NULL, add = combined.model.add, dev = NULL)
    
    # make pheno list
    pheno.list = list(phenotype.s1 = phenotype.s1, phenotype.s2 = phenotype.s2)
    # effects list
    effect.list = list(effects.add = effects.scaled.add, effects.dev = NULL)
    # make causal SNPs list
    causal.snps.list = list(add.snps = causal.snps.add, dev.snps = NULL)
    # models list 
    models.list = list(combined.models = combined.models, s1.models = s1.models, s2.models = s2.models)
    # make results list
    result = list(phenotype = pheno.list, snps = causal.snps.list, models = models.list, effects = effect.list) # effects = effect.list
    # return
    return(result)
    
  } else if ((n.add.snps == 0) & (n.dev.snps != 0)) {
    # Genetic value per individual
    g.dom = rowSums(effects.scaled.dev)
    g.tot = g.dom   # total genetic component
    
    # Add environmental noise so that total h2 ≈ 0.10
    pheno.var = var(g.tot)
    env.var  = pheno.var * (1 - total.h2) / total.h2
    noise = rnorm(nrow(genos), mean = 0, sd = sqrt(env.var))
    
    # Generate phenotype
    phenotype = g.tot + noise
    
    # Split combined data into samples
    phenotype.s1 = phenotype[1:nrow(genos1)]
    phenotype.s2 = phenotype[(nrow(genos1) + 1):nrow(genos)]
    
    effects.dev.s1 = rowSums(effects.scaled.dev)[1:nrow(genos1)]
    effects.dev.s2 = rowSums(effects.scaled.dev)[(nrow(genos1) + 1):nrow(genos)]
    
    # Run models for sample 1
    s1.model.dev = lm(phenotype.s1 ~ effects.dev.s1)
    s1.models = list(mixed = NULL, add = NULL, dev = s1.model.dev)
    
    # Run models for sample 2
    s2.model.dev = lm(phenotype.s2 ~ effects.dev.s2)
    s2.models = list(mixed = NULL, add = NULL, dev = s2.model.dev)
    
    # run models for combined phenotype
    combined.model.dev = lm(phenotype ~ rowSums(effects.scaled.dev))
    combined.models = list(mixed = NULL, add = NULL, dev = combined.model.dev)
    
    # make pheno list
    pheno.list = list(phenotype.s1 = phenotype.s1, phenotype.s2 = phenotype.s2)
    # effects list
    effect.list = list(effects.add = NULL, effects.dev = effects.scaled.dev)
    # make causal SNPs list
    causal.snps.list = list(add.snps = NULL, dev.snps = causal.snps.dev)
    # models list 
    models.list = list(combined.models = combined.models, s1.models = s1.models, s2.models = s2.models)
    # make results list
    result = list(phenotype = pheno.list, snps = causal.snps.list, models = models.list, effects = effect.list) # effects = effect.list
    # return
    return(result)
    
  }
  
}

# create split column for hybrid cross
hybrid_cv_split = function(df, test_frac = 0.1, k = 5, seed = 707, split_colname = "split1") {
  set.seed(seed)
  n <- nrow(df)
  test_size <- round(test_frac * n)
  
  # Randomly select test set indices
  test_indices <- sample(n, size = test_size)
  train_val_indices <- setdiff(1:n, test_indices)
  train_val_size <- length(train_val_indices)
  
  # Assign folds to the remaining data
  fold_ids <- sample(rep(1:k, length.out = train_val_size))
  folds <- split(train_val_indices, fold_ids)
  
  # Initialize split column
  df[[split_colname]] <- NA_character_
  df[[split_colname]][test_indices] <- "test"
  
  # Assign fold numbers
  for (i in 1:k) {
    df[[split_colname]][folds[[i]]] <- as.character(i)
  }
  
  return(df)
}

# # Function to compute additive PRS
compute_additive_prs = function(geno, stats, additive_snps) {
  additive_prs = numeric(nrow(geno))
  
  for (i in seq_along(additive_snps)) {
    snp = additive_snps[i]
    
    if (!snp %in% colnames(geno)) {
      print(paste("SNP not found in geno:", snp))
      next
    }
    
    snp_rows = stats[stats$ID2 == snp & stats$TEST == "ADD", ]
    if (nrow(snp_rows) == 0) {
      print(paste("Missing summary stats for SNP:", snp))
      next
    }
    
    if (nrow(snp_rows) > 1) {
      print(paste("Warning: multiple ADD rows for SNP:", snp))
    }
    
    beta = snp_rows$BETA[1]
    snp_genotypes = as.numeric(geno[[snp]])
    
    if (all(is.na(snp_genotypes))) {
      print(paste("Skipping SNP with all missing genotypes:", snp))
      next
    }
    
    if (any(is.na(snp_genotypes))) {
      print(paste("NA genotypes detected for SNP:", snp))
      print(paste("Count of NA genotypes:", sum(is.na(snp_genotypes))))
    }
    
    additive_prs = additive_prs + ifelse(is.na(snp_genotypes), 0, snp_genotypes * beta)
    
    if (i %% 100 == 0) {
      cat("Processed", i, "additive SNPs\n")
    }
  }
  
  return(additive_prs)
}

# Function to compute domdev-adjusted PRS
compute_domdev_prs = function(geno, stats, domdev_snps, stats2 = NULL) {
  domdev_prs = numeric(nrow(geno))
  
  for (i in seq_along(domdev_snps)) {
    snp = domdev_snps[i]
    
    if (!snp %in% colnames(geno)) {
      print(paste("SNP not found in geno:", snp))
      next
    }
    
    snp_rows = stats[stats$ID2 == snp, ]
    if (nrow(snp_rows) == 0) {
      print(paste("Missing summary stats for SNP:", snp))
      next
    }
    
    if (nrow(snp_rows) > 3) {
      print(paste("Warning: multiple rows for SNP:", snp))
    }
    
    beta_domdev = snp_rows[snp_rows$TEST == "DOMDEV", ]$BETA[1]
    beta_additive = snp_rows[snp_rows$TEST == "ADD", ]$BETA[1]
    
    snp_genotypes = as.numeric(geno[[snp]])
    
    if (all(is.na(snp_genotypes))) {
      print(paste("Skipping SNP with all missing genotypes:", snp))
      next
    }
    
    if (any(is.na(snp_genotypes))) {
      print(paste("NA genotypes detected for SNP:", snp))
      print(paste("Count of NA genotypes:", sum(is.na(snp_genotypes))))
    }
    
    # Handle missing beta values
    if (is.na(beta_domdev) & is.na(beta_additive)) {
      print(paste("Missing all beta values for SNP:", snp))
      if (!is.null(stats2)) {
        fallback_rows = stats2[stats2$ID2 == snp & stats2$TEST == "ADD", ]
        if (nrow(fallback_rows) > 0) {
          beta_additive = fallback_rows$BETA[1]
          score = ifelse(is.na(snp_genotypes), 0, snp_genotypes * beta_additive)
          domdev_prs = domdev_prs + score
        } else {
          print(paste("No fallback summary stats found for SNP:", snp))
        }
      }
      next
    } else if (is.na(beta_domdev)) {
      print(paste("Missing DOMDEV beta for SNP:", snp))
      next
    } else if (is.na(beta_additive)) {
      print(paste("Missing ADD beta for SNP:", snp))
      next
    }
    
    # Vectorized scoring with NA-safe logic
    score = ifelse(is.na(snp_genotypes), 0,
                   ifelse(snp_genotypes == 0, 0,
                          ifelse(snp_genotypes == 1, beta_additive + beta_domdev,
                                 ifelse(snp_genotypes == 2, 2 * beta_additive, 0))))
    
    domdev_prs = domdev_prs + score
    
    if (i %% 100 == 0) {
      cat("Processed", i, "domdev SNPs\n")
    }
  }
  
  return(domdev_prs)
}

# Function for running X-fold CV w/ linear regression
# Outputs: 1) raw object from lm(); 2) iid-level data; 3) R2 per fold
run.cv.regression = function(data, predictor.col) {
  # linear regression model results
  models = list()
  iid.probs = list()
  split.metrics = list()
  
  # loop through data and run logistic regression 10x and save results
  for (fold in 1:5) {
    
    # # 10 = test split
    # # split data (while leaving out tuning set)
    # data.train = subset(data, split1 != fold & split1 != 10)
    # # data.tune not needed
    # data.test = subset(data, split1 == 10)
    
    # split data (while leaving out tuning set)
    if (fold != 5){
      data.train = subset(data, split2 != fold & split2 != (fold+1))
      data.test = subset(data, split2 == (fold))
    } else if (fold == 5) {
      data.train = subset(data, split2 != fold & split2 != (fold-4))
      data.test = subset(data, split2 == (fold))
    }
    
    # Linear regression model
    formula = as.formula(paste("phenotype ~", predictor.col))
    lin.model = lm(formula, data = data.train)
    
    # save model
    models[[fold]] = lin.model
    
    # predict in tuning set 
    split.probs = as.data.frame(predict(lin.model, newdata=data.test, type="response")); colnames(split.probs) = 'probs'
    
    # model predictions
    split.probs$phenotype = as.numeric(as.character(data.test$phenotype))
    split.probs$ids = data.test$IID
    
    # save model probs
    iid.probs[[fold]] = split.probs
    
    # get variance explained for model predictions
    r2 = cor.test(split.probs$probs, data.test$phenotype)[['estimate']]^2
    
    # empty list for metrics
    metrics = list()
    
    # record variance explained
    metrics[["r2"]] = r2
    split.metrics[[fold]] = metrics
    
    # remove 
    rm(metrics, r2, split.probs)
    
  }
  
  # Extracting 'accuracy' values and computing mean, min, max
  r2_values = sapply(split.metrics, function(x) x$r2)
  mean_r2 = mean(r2_values)
  min_r2 = min(r2_values)
  max_r2 = max(r2_values)
  
  # Printing results
  cat("R2 - Mean:", mean_r2, "Min:", min_r2, "Max:", max_r2, "\n")
  
  # merge probs and preds from all test sets 
  merged.iid = do.call(rbind, iid.probs)
  
  # merge probs and preds from all test sets 
  models.r2 = do.call(rbind, split.metrics)
  
  # create output
  results = list(models, merged.iid, models.r2)
  names(results) = c("model.out", "iid.data", "model.2")
  
  return(results)
  
}

# function for logging results
log_results = function(sim.results, code) {
  if (code == "ADD") {
    sink(log.file)
    print("ADD PHENO ADD MODEL - Combined:")
    cat("R-squared:", summary(sim.results[["models"]][["combined.models"]][["add"]])$r.squared, '\n')
    print("ADD PHENO ADD MODEL - Sample 1:")
    cat("R-squared:", summary(sim.results[["models"]][["s1.models"]][["add"]])$r.squared, '\n')
    print("ADD PHENO ADD MODEL - Sample 2:")
    cat("R-squared:", summary(sim.results[["models"]][["s2.models"]][["add"]])$r.squared, '\n')
    sink()
  }
  
  if (code == "DEV") {
    sink(log.file)
    print("DEV PHENO DEV MODEL - Combined:")
    cat("R-squared:", summary(sim.results[["models"]][["combined.models"]][["dev"]])$r.squared, '\n')
    print("DEV PHENO DEV MODEL - Sample 1:")
    cat("R-squared:", summary(sim.results[["models"]][["s1.models"]][["dev"]])$r.squared, '\n')
    print("DEV PHENO DEV MODEL - Sample 2:")
    cat("R-squared:", summary(sim.results[["models"]][["s2.models"]][["dev"]])$r.squared, '\n')
    sink()
  }
  
  if (code == "MIXED") {
    sink(log.file)
    print("MIXED PHENO FULL MODEL - Combined:")
    cat("R-squared:", summary(sim.results[["models"]][["combined.models"]][["mixed"]])$r.squared, '\n')
    print("MIXED PHENO ADD MODEL - Combined:")
    cat("R-squared:", summary(sim.results[["models"]][["combined.models"]][["add"]])$r.squared, '\n')
    print("MIXED PHENO DEV MODEL - Combined:")
    cat("R-squared:", summary(sim.results[["models"]][["combined.models"]][["dev"]])$r.squared, '\n')
    
    print("MIXED PHENO FULL MODEL - Sample 1:")
    cat("R-squared:", summary(sim.results[["models"]][["s1.models"]][["mixed"]])$r.squared, '\n')
    print("MIXED PHENO ADD MODEL - Sample 1:")
    cat("R-squared:", summary(sim.results[["models"]][["s1.models"]][["add"]])$r.squared, '\n')
    print("MIXED PHENO DEV MODEL - Sample 1:")
    cat("R-squared:", summary(sim.results[["models"]][["s1.models"]][["dev"]])$r.squared, '\n')
    
    print("MIXED PHENO FULL MODEL - Sample 2:")
    cat("R-squared:", summary(sim.results[["models"]][["s2.models"]][["mixed"]])$r.squared, '\n')
    print("MIXED PHENO ADD MODEL - Sample 2:")
    cat("R-squared:", summary(sim.results[["models"]][["s2.models"]][["add"]])$r.squared, '\n')
    print("MIXED PHENO DEV MODEL - Sample 2:")
    cat("R-squared:", summary(sim.results[["models"]][["s2.models"]][["dev"]])$r.squared, '\n')
    sink()
  }
}

# ----
# # Function to compute additive PRS
# compute_additive_prs = function(geno, stats, additive_snps) {
#   # Initialize PRS scores
#   additive_prs = numeric(nrow(geno))
#   
#   # Compute PRS
#   for (snp in additive_snps) {
#     snp_rows = stats[stats$ID2 == snp & stats$TEST == "ADD", ]
#     if (nrow(snp_rows) > 0) {
#       print(paste("Using SNP: ", snp))
#       beta = snp_rows$BETA
#       snp_genotypes = geno[, ..snp, drop = FALSE]
#       # error checks
#       if (any(is.na(snp_genotypes))) {
#         print(paste("NA genotypes detected for SNP:", snp))
#         print(paste("Count of NA genotypes:", sum(is.na(snp_genotypes))))
#       }
#       
#       if (!snp %in% colnames(geno)) {
#         print(paste("SNP not found in geno:", snp))
#         next
#       }
#       # compute PRS
#       additive_prs = additive_prs + (snp_genotypes * beta)
#     } else {
#       print(paste("Missing SNP: ", snp))
#       print("Check that SNP ID formats match")
#     }
#   }
#   
#   return(additive_prs)
# }

# # Function to compute dominance deviation PRS
# compute_domdev_prs = function(geno, stats, domdev_snps, stats2 = NA) {
#   # Initialize PRS scores
#   domdev_prs = numeric(nrow(geno))
#   
#   # Compute PRS
#   for (snp in domdev_snps) {
#     
#     snp_rows = stats[stats$ID2 == snp, ]
#     if (nrow(snp_rows) > 0) {
#       print(paste("Using SNP: ", snp))
#       beta_domdev = snp_rows[snp_rows$TEST == "DOMDEV",]$BETA
#       beta_additive = snp_rows[snp_rows$TEST == "ADD", ]$BETA
#       snp_genotypes = as.matrix(geno[, ..snp, drop = FALSE])
#       
#       if (any(is.na(snp_genotypes))) {
#         print(paste("NA genotypes detected for SNP:", snp))
#         print(paste("Count of NA genotypes:", sum(is.na(snp_genotypes))))
#       }
#       
#       if (!snp %in% colnames(geno)) {
#         print(paste("SNP not found in geno:", snp))
#         next
#       }
#       
#       # Check for missing beta values
#       if (is.na(beta_domdev) & is.na(beta_additive)) {
#         print(paste("Missing all beta value for SNP:", snp))
#         print("Using ADD sumstats instead")
#         snp_rows = stats2[stats2$ID2 == snp, ]
#         if (nrow(snp_rows) > 0) {
#           beta_additive = snp_rows[snp_rows$TEST == "ADD", ]$BETA
#           snp_genotypes = as.matrix(geno[, ..snp, drop = FALSE])
#           additive_est = snp_genotypes * beta_additive
#           domdev_prs = domdev_prs + additive_est
#           next
#         }
#       } else if (is.na(beta_domdev)) {
#         print(paste("Missing DOMDEV beta value for SNP:", snp))
#         next
#       } else if (is.na(beta_additive)) {
#         print(paste("Missing ADD beta value for SNP:", snp))
#         next
#       }
#       
#       domdev_prs = domdev_prs + apply(snp_genotypes, 1, function(x) {
#         if (x == 0) {
#           return(0)
#         } else if (x == 1) {
#           return(beta_additive + beta_domdev)
#         } else if (x == 2) {
#           return(2 * beta_additive)
#         } else {
#           return(0) # Handle unexpected values
#         }
#       })
#     } else {
#       print(paste("SNP missing:", snp))
#     }
#   }
#   
#   return(domdev_prs)
# }
