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

# source functions
source("~/phd_code/dl-prs-v2/code/nyb_sims_code_v2/pheno.functions.R")

args = commandArgs(trailingOnly = TRUE)

# input arguments
data.file = args[1]
out.dir = args[2]

# load data
print("... Loading data file ...", quote = F)
data = fread(data.file)

# convert "test" to 10
data$split1[data$split1 == "test"] = "10"
data$split1 = as.integer(data$split1)

# get name
file.name = sub(".*DATA_(.+)\\.txt$", "\\1", data.file)

# Extract SNP-h2 proportions
h2.a = as.numeric(str_extract(file.name, "(?<=_a)[0-9.]+"))
h2.d = as.numeric(str_extract(file.name, "(?<=_d)[0-9.]+"))
print(paste("Running for phenotype: h2-ADD = ", h2.a, "; h2-DOM = ", h2.d), quote = F)

if (h2.d != 0) {
  print("... Running regression for ADD & DOM ...", quote = F)
  # make add output files
  add.iid.out = paste(out.dir, "/ADD_IID_OUT_", file.name, ".txt", sep="")
  add.model.r2.out = paste(out.dir, "/ADD_R2_OUT_", file.name, ".txt", sep="")
  # make dom output files
  dom.iid.out = paste(out.dir, "/DOM_IID_OUT_", file.name, ".txt", sep="")
  dom.model.r2.out = paste(out.dir, "/DOM_R2_OUT_", file.name, ".txt", sep="")
  # run regression 
  print("... Computing ADD results ...", quote = F)
  add.results = run.cv.regression(data, "additive.prs")
  print("... Computing DOM results ...", quote = F)
  dom.results = run.cv.regression(data, "domdev.prs")
  # save add
  fwrite(add.results[[2]], file = add.iid.out, quote = F, sep = "\t", row.names = F, col.names = T)
  fwrite(as.data.frame(add.results[[3]]), file = add.model.r2.out, quote = F, sep = "\t", row.names = F, col.names = T)
  # save dom
  fwrite(dom.results[[2]], file = dom.iid.out, quote = F, sep = "\t", row.names = F, col.names = T)
  fwrite(as.data.frame(dom.results[[3]]), file = dom.model.r2.out, quote = F, sep = "\t", row.names = F, col.names = T)
} else if (h2.d == 0) {
  print("... Running regression for ADD only ...", quote = F)
  # make add output files
  add.iid.out = paste(out.dir, "/ADD_IID_OUT_", file.name, ".txt", sep="")
  add.model.r2.out = paste(out.dir, "/ADD_R2_OUT_", file.name, ".txt", sep="")
  # run regression 
  print("... Computing ADD results ...", quote = F)
  add.results = run.cv.regression(data, "additive.prs")
  # save add
  fwrite(add.results[[2]], file = add.iid.out, quote = F, sep = "\t", row.names = F, col.names = T)
  fwrite(as.data.frame(add.results[[3]]), file = add.model.r2.out, quote = F, sep = "\t", row.names = F, col.names = T)
} 

print("<<>> - FINISHED - <<>>", quote = F)
print("", quote = F)

