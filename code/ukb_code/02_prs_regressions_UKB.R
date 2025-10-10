# Run regression models for PRS Scores
# Author: Nate
# Date: Jan 9, 2025

# libraries
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(caret))
suppressMessages(library(stringr))
suppressMessages(library(data.table))

# print starting message
print("", quote = F)
print("<<>> - Starting regression model testing - <<>>", quote = F)
print("... (woop woop) ...", quote = F)

args = commandArgs(trailingOnly = TRUE)

# input arguments
data.file = args[1]
out.dir = args[2]

# Function for running X-fold CV w/ linear regression
# Outputs: 1) raw object from lm(); 2) iid-level data; 3) R2 per fold
run.cv.regression.ukb = function(data, reg.formula) {
  # linear regression model results
  models = list()
  iid.probs = list()
  split.metrics = list()
  
  # loop through data and run logistic regression 10x and save results
  for (fold in 1:5) {
    
    # split data (while leaving out tuning set)
    if (fold != 5){
      data.train = subset(data, split1 != fold & split1 != (fold+1))
      data.test = subset(data, split1 == (fold+1))
    } else if (fold == 5) {
      data.train = subset(data, split1 != fold & split1 != (fold-4))
      data.test = subset(data, split1 == (fold-4))
    }
    
    # Linear regression model
    formula = as.formula(reg.formula)
    lin.model = lm(formula, data = data.train)
    summary(lin.model)
    
    # save model
    models[[fold]] = lin.model
    
    # predict in tuning set 
    split.probs = as.data.frame(predict(lin.model, newdata=data.test, type="response")); colnames(split.probs) = 'probs'
    
    # model predictions
    split.probs$int_phenotype = as.numeric(as.character(data.test$int_phenotype))
    split.probs$ids = data.test$IID
    
    # save model probs
    iid.probs[[fold]] = split.probs
    
    # get variance explained for model predictions
    r2 = cor.test(split.probs$probs, data.test$int_phenotype)[['estimate']]^2
    
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

# Function for running regression on 70%, 15%, 15% split
# Outputs: 1) raw object from lm(); 2) iid-level data; 3) R2 
run.sff.regression.ukb = function(data, reg.formula) {
  # linear regression model results
  models = list()
  iid.probs = list()
  models.r2 = list()
  
  # split data
  data.train = subset(data, split3 == "train")
  data.tune = subset(data, split3 == "tune")
  data.test = subset(data, split3 == "test")
  
  # Linear regression model
  formula = as.formula(reg.formula)
  lin.model = lm(formula, data = data.train)
  summary(lin.model)
  
  # save model
  models[["train"]] = lin.model
  
  # predict in tuneing set 
  tune.probs = as.data.frame(predict(lin.model, newdata=data.tune, type="response")); colnames(tune.probs) = 'probs'
  
  # model predictions
  tune.probs$int_phenotype = as.numeric(as.character(data.tune$int_phenotype))
  tune.probs$ids = data.tune$IID
  
  # save model probs
  iid.probs[["tune"]] = tune.probs
  
  # get variance explained for model predictions
  tune.r2 = cor.test(tune.probs$probs, data.tune$int_phenotype)[['estimate']]^2
  
  # predict in testing set 
  test.probs = as.data.frame(predict(lin.model, newdata=data.test, type="response")); colnames(test.probs) = 'probs'
  
  # model predictions
  test.probs$int_phenotype = as.numeric(as.character(data.test$int_phenotype))
  test.probs$ids = data.test$IID
  
  # save model probs
  iid.probs[["test"]] = test.probs
  
  # get variance explained for model predictions
  test.r2 = cor.test(test.probs$probs, data.test$int_phenotype)[['estimate']]^2
  
  # Printing results
  cat("Tune R2:", tune.r2, "Test R2:", test.r2, "\n")
  
  # save in metrics
  models.r2[["tune.r2"]] = tune.r2
  models.r2[["test.r2"]] = test.r2
  
  # create output
  results = list(models, iid.probs, models.r2)
  names(results) = c("model.out", "iid.data", "model.r2")
  return(results)
  
}  

# load data
print("... Loading data file ...", quote = F)
data = fread(data.file)

# get name
file.name = sub(".*DATA_(.+)\\.txt$", "\\1", data.file)

# Extracting variances for the following models
# int_phenotype ~ prs                                                 # geno
# int_phenotype ~ age + sex + array01 + pc1 ... + pc10                # covs
# int_phenotype ~ prs + age + sex + array01 + pc1 ... + pc10          # full

# -- define formulas -- 
print("... Defining model formulas ...")

# covs only 
f.covars = "int_phenotype ~ age + sex + array01 + pop_pc1 + pop_pc2 + pop_pc3 + pop_pc4 + pop_pc5 + pop_pc6 + pop_pc7 + pop_pc8 + pop_pc9 + pop_pc10"

# additive
f.geno.add = "int_phenotype ~ additive.prs"
f.full.add = "int_phenotype ~ additive.prs + age + sex + array01 + pop_pc1 + pop_pc2 + pop_pc3 + pop_pc4 + pop_pc5 + pop_pc6 + pop_pc7 + pop_pc8 + pop_pc9 + pop_pc10"

# domdev
f.geno.dom = "int_phenotype ~ domdev.prs"
f.full.dom = "int_phenotype ~ domdev.prs + age + sex + array01 + pop_pc1 + pop_pc2 + pop_pc3 + pop_pc4 + pop_pc5 + pop_pc6 + pop_pc7 + pop_pc8 + pop_pc9 + pop_pc10"

print("... Starting Part 1: Complete PRS Analyses ...")

print("... Computing COVAR only results ...", quote = F)
cat("COVAR results:", "\n")
covar.results = run.sff.regression.ukb(data, f.covars)

print("... Saving ...")
# make add output files
covar.model.out = file.path(out.dir, paste("COVAR_MODEL_OUT_", file.name, ".RData", sep=""))
covar.model.r2.out = file.path(out.dir, paste("COVAR_R2_OUT_", file.name, ".txt", sep=""))
# save geno
save(covar.results, file = covar.model.out)
fwrite(as.data.frame(covar.results[[3]]), file = covar.model.r2.out, quote = F, sep = "\t", row.names = F, col.names = T)

print("... Computing ADD results ...", quote = F)
cat("ADD PRS results:", "\n")
add.geno.results = run.sff.regression.ukb(data, f.geno.add)
cat("ADD PRS + COVARS results:", "\n")
add.full.results = run.sff.regression.ukb(data, f.full.add)

print("... Saving ...")
# make add output files
add.geno.model.out = file.path(out.dir, paste("ADD_GENO_MODEL_OUT_", file.name, ".RData", sep=""))
add.geno.model.r2.out = file.path(out.dir, paste("ADD_GENO_R2_OUT_", file.name, ".txt", sep=""))
add.full.model.out = file.path(out.dir, paste("ADD_FULL_MODEL_OUT_", file.name, ".RData", sep=""))
add.full.model.r2.out = file.path(out.dir, paste("ADD_FULL_R2_OUT_", file.name, ".txt", sep=""))
# save geno
save(add.geno.results, file = add.geno.model.out)
fwrite(as.data.frame(add.geno.results[[3]]), file = add.geno.model.r2.out, quote = F, sep = "\t", row.names = F, col.names = T)
# save full
save(add.full.results, file = add.full.model.out)
fwrite(as.data.frame(add.full.results[[3]]), file = add.full.model.r2.out, quote = F, sep = "\t", row.names = F, col.names = T)

print("... Computing DEV results ...", quote = F)
cat("DOM PRS results:", "\n")
dom.geno.results = run.sff.regression.ukb(data, f.geno.dom)
cat("DOM PRS + COVARS results:", "\n")
dom.full.results = run.sff.regression.ukb(data, f.full.dom)

print("... Saving ...")
# make dom output files
dom.geno.model.out = file.path(out.dir, paste("DOM_GENO_MODEL_OUT_", file.name, ".RData", sep=""))
dom.geno.model.r2.out = file.path(out.dir, paste("DOM_GENO_R2_OUT_", file.name, ".txt", sep=""))
dom.full.model.out = file.path(out.dir, paste("DOM_FULL_MODEL_OUT_", file.name, ".RData", sep=""))
dom.full.model.r2.out = file.path(out.dir, paste("DOM_FULL_R2_OUT_", file.name, ".txt", sep=""))
# save geno
save(dom.geno.results, file = dom.geno.model.out)
fwrite(as.data.frame(dom.geno.results[[3]]), file = dom.geno.model.r2.out, quote = F, sep = "\t", row.names = F, col.names = T)
# save full
save(dom.full.results, file = dom.full.model.out)
fwrite(as.data.frame(dom.full.results[[3]]), file = dom.full.model.r2.out, quote = F, sep = "\t", row.names = F, col.names = T)

print("... Starting Part 2: DOMDEV-only SNPs PRS Analyses ...")

# addditive                                                       # "addditive" is INTENTIONAL - typo from the 01b script.
f.geno.add.DOMCOMP = "int_phenotype ~ addditive.dom.comp"
f.full.add.DOMCOMP = "int_phenotype ~ addditive.dom.comp + age + sex + array01 + pop_pc1 + pop_pc2 + pop_pc3 + pop_pc4 + pop_pc5 + pop_pc6 + pop_pc7 + pop_pc8 + pop_pc9 + pop_pc10"

# domdev
f.geno.dom.DOMCOMP = "int_phenotype ~ domdev.dom.comp"
f.full.dom.DOMCOMP = "int_phenotype ~ domdev.dom.comp + age + sex + array01 + pop_pc1 + pop_pc2 + pop_pc3 + pop_pc4 + pop_pc5 + pop_pc6 + pop_pc7 + pop_pc8 + pop_pc9 + pop_pc10"

print("... Computing ADD results ...", quote = F)
cat("ADD PRS results:", "\n")
add.geno.results.domcomp = run.sff.regression.ukb(data, f.geno.add.DOMCOMP)
cat("ADD PRS + COVARS results:", "\n")
add.full.results.domcomp = run.sff.regression.ukb(data, f.full.add.DOMCOMP)

print("... Saving ...")
# make add output files
add.geno.model.out.domcomp = file.path(out.dir, paste("ADD_DOMCOMP_GENO_MODEL_OUT_", file.name, ".RData", sep=""))
add.geno.model.r2.out.domcomp = file.path(out.dir, paste("ADD_DOMCOMP_GENO_R2_OUT_", file.name, ".txt", sep=""))
add.full.model.out.domcomp = file.path(out.dir, paste("ADD_DOMCOMP_FULL_MODEL_OUT_", file.name, ".RData", sep=""))
add.full.model.r2.out.domcomp = file.path(out.dir, paste("ADD_DOMCOMP_FULL_R2_OUT_", file.name, ".txt", sep=""))
# save geno
save(add.geno.results.domcomp, file = add.geno.model.out.domcomp)
fwrite(as.data.frame(add.geno.results.domcomp[[3]]), file = add.geno.model.r2.out.domcomp, quote = F, sep = "\t", row.names = F, col.names = T)
# save full
save(add.full.results.domcomp, file = add.full.model.out.domcomp)
fwrite(as.data.frame(add.full.results.domcomp[[3]]), file = add.full.model.r2.out.domcomp, quote = F, sep = "\t", row.names = F, col.names = T)

print("... Computing DEV results ...", quote = F)
cat("DOM PRS results:", "\n")
dom.geno.results.domcomp = run.sff.regression.ukb(data, f.geno.dom.DOMCOMP)
cat("DOM PRS + COVARS results:", "\n")
dom.full.results.domcomp = run.sff.regression.ukb(data, f.full.dom.DOMCOMP)

print("... Saving ...")
# make dom output files
dom.geno.model.out.domcomp = file.path(out.dir, paste("DOM_DOMCOMP_GENO_MODEL_OUT_", file.name, ".RData", sep=""))
dom.geno.model.r2.out.domcomp = file.path(out.dir, paste("DOM_DOMCOMP_GENO_R2_OUT_", file.name, ".txt", sep=""))
dom.full.model.out.domcomp = file.path(out.dir, paste("DOM_DOMCOMP_FULL_MODEL_OUT_", file.name, ".RData", sep=""))
dom.full.model.r2.out.domcomp = file.path(out.dir, paste("DOM_DOMCOMP_FULL_R2_OUT_", file.name, ".txt", sep=""))
# save geno
save(dom.geno.results.domcomp, file = dom.geno.model.out.domcomp)
fwrite(as.data.frame(dom.geno.results.domcomp[[3]]), file = dom.geno.model.r2.out.domcomp, quote = F, sep = "\t", row.names = F, col.names = T)
# save full
save(dom.full.results.domcomp, file = dom.full.model.out.domcomp)
fwrite(as.data.frame(dom.full.results.domcomp[[3]]), file = dom.full.model.r2.out.domcomp, quote = F, sep = "\t", row.names = F, col.names = T)

print("<<>> - FINISHED - <<>>", quote = F)
print("", quote = F)

