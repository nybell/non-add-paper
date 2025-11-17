# Run regression models for PRS Scores
# Author: Nate
# Date: Nov 17, 2025

# CELL 1 ----

# libraries
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(caret))
suppressMessages(library(stringr))
suppressMessages(library(data.table))

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
  cat("R2 - Mean:", mean_r2, "\n")
  cat("R2 - SD:", sd(r2_values), "\n")
  
  # merge probs and preds from all test sets 
  merged.iid = do.call(rbind, iid.probs)
  
  # merge probs and preds from all test sets 
  models.r2 = do.call(rbind, split.metrics)
  
  # create output
  results = list(models, merged.iid, models.r2)
  names(results) = c("model.out", "iid.data", "model.2")
  
  return(results)
  
}

# CELL 2 ----

# input arguments
data.file = "/Users/nyb/demo_data/DATA_eur_nsnps100_h0.5_a0_d0.5_50k.txt"             # EDIT TO MATCH YOUR FILE PATH
add.pub.file = "/Users/nyb/demo_data/ADD_R2_OUT_eur_nsnps100_h0.5_a0_d0.5_50k.txt"    # EDIT TO MATCH YOUR FILE PATH
dom.pub.file = "/Users/nyb/demo_data/DOM_R2_OUT_eur_nsnps100_h0.5_a0_d0.5_50k.txt"    # EDIT TO MATCH YOUR FILE PATH

# load data
print("... Loading data file ...", quote = F)
data = fread(data.file)

# convert "test" to 10
data$split1[data$split1 == "test"] = "10"
data$split1 = as.integer(data$split1)

# get name
file.name = sub(".*DATA_(.+)\\.txt$", "\\1", data.file)

# CELL 3 ----

# Extract SNP-h2 proportions
h2.a = as.numeric(str_extract(file.name, "(?<=_a)[0-9.]+"))
h2.d = as.numeric(str_extract(file.name, "(?<=_d)[0-9.]+"))
print(paste("Running for phenotype: h2-ADD = ", h2.a, "; h2-DOM = ", h2.d), quote = F)

if (h2.d != 0) {
  print("... Running regression for ADD & DOM ...", quote = F)
  # run regression 
  print("... Computing ADD results ...", quote = F)
  add.results = run.cv.regression(data, "additive.prs")
  print("... Computing DOM results ...", quote = F)
  dom.results = run.cv.regression(data, "domdev.prs")
} else if (h2.d == 0) {
  print("... Running regression for ADD only ...", quote = F)
  # run regression 
  print("... Computing ADD results ...", quote = F)
  add.results = run.cv.regression(data, "additive.prs")
} 

# ---- CELL 4 COMPARE w/ published results ----

# additive-only PGS
add.pub = fread(add.pub.file)
cat("Additive-PGS Mean R2 (published):", (mean(add.pub$r2)), "\n")
cat("Additive-PGS SD R2 (published):", (sd(add.pub$r2)), "\n")

# additive-only PGS
dom.pub = fread(dom.pub.file)
cat("Dominance-adjusted PGS Mean R2 (published):", (mean(dom.pub$r2)), "\n")
cat("Dominance-adjusted PGS SD R2 (published):", (sd(dom.pub$r2)), "\n")









