# Benchmarking non-additive genetic effects on polygenic prediction and machine learning-based approaches

DOI: https://doi.org/10.1101/2025.10.10.25337750

This repository contains code, configuration files, data (simulated phenotypes), and results (from simulation analyses) from the manuscript "Benchmarking non-additive genetic effects on polygenic prediction and machine learning-based approaches". 

## Overview
The project evaluates how dominance deviations impact polygenic score (PGS) performance, and compares linear models using PGSs to XGboost and a neural network using both simulated data and real-data from the UK Biobank. All code and simulated phenotype files for the simulated analyses can be found in this repo. UK Biobank data can be accessed on their platform for approved researchers (https://www.ukbiobank.ac.uk/use-our-data/research-analysis-platform/).

Complete data for the project (raw genotypes + data prepared for ML/DL models) and UK Biobank summary statistics:     
https://zenodo.org/records/17552313     

## Software & Hardware
Local code was run on a 2021 MacBook Pro M1 Max Silicon (32GB RAM). HPC code was run on the Dutch Snellius computer, with the ML/DL models using a single NVIDIA A100 GPU.      

GWASs were run using PLINK v2.00a6LM 64-bit Intel (18 Apr 2024), which can be found here: https://www.cog-genomics.org/plink/2.0/       

## Environments Setup

### Local (MacOS M1 Pro; Sonoma 14.5)
Code to set up virtual environment with package versions used for local code (MacOS; Est. time < 5 - 10 minutes). 
```
# if pyenv is not installed
brew install pyenv pyenv-virtualenv
pyenv install 3.10.12

# if pyenv is installed, start here
pyenv virtualenv 3.10.12 myenv-3.10
pyenv activate myenv-3.10
pip install -r /non-add-paper/code/required_packages_local.txt
```

### Snellius HPC
Directions for setting up environment on HPC (Est. time = 5 - 10 minutes).      
If Miniconda is not installed (more information here: https://www.anaconda.com/docs/getting-started/miniconda/install#linux-terminal-installer)
```
cd $HOME
https://www.anaconda.com/docs/getting-started/miniconda/install#linux-terminal-installer    # download miniconda
bash miniconda.sh -b -p $HOME/miniconda3                                                    # run installer
source $HOME/miniconda3/etc/profile.d/conda.sh                                              # source & init
conda init bash
```
If needed, logout and re-login to your HPC environment and run:
```
source ~/.bashrc
```

Create virtual environment with dependencies
```
source $HOME/miniconda3/etc/profile.d/conda.sh
conda env create -f /non-add-paper/code/non-add-pgs-hpc-env.yml
conda activate ml_models
```

### Worflow Summary    
(1) Simulate phenotypes and PLINK .fam files (LOCAL)
- 01a_pheno_sim*.R
- 01b_devs_pheno_sim*.R

(2) Run additive and dominance GWAS (per chrom) and merge (HPC)
- 02_run_linear_gwas.sh
- syn.gwas.linear.v2.job.txt
- 02_run_domdev_gwas.sh
- syn.gwas.geno.v2.job.txt
- 02a_merge_all_sumstats.sh

(3) QC + compute PGSs + data formatted for model input (HPC)
- 03_ready.data.ss.job.txt

(4) Run PGS regressions
- 04_run.prs.regressions.ss.job.txt

(5 & 6) Run XGBoost and neural network models
- model_definition.py 
- dnn.job.txt
- 5_submit_dnn.sh
- xgb.job.txt
- 06_submit_xgb.sh

(7+) Visualization & figures     

Notes:     
(a) Local simulated phenotype generation will generate PLINK .fam files and a corresponding "nsnp_*.RData" file. These need to be moved to the HPC and put in the "$FAM_NSNPS" directories in order to run the GWASs.       
(b) GWASs are run with PLINK2, with one SLURM array job submitted per phenotype.      
(c) The 05 and 06 scripts submit one SLURM job per phenotype.       
(d) _The filepaths in the code have **not** been altered! If running the full pipeline, you will need to adjust to your own system. As an alternative, a test jupyter notebooks + R script can be found to run the XGBoost, and NN models locally._      

### Full pipeline
Code as run for project. 

#### Generate simulated phenotypes (local)
Executable R scripts that generate synthetic phenotypes, PLINK .fam files, and performs QC checks. Need a config file defining phenotype parameters (see '/data/fams_fin_snps100/sim.fin.config.snps100.csv' for example). Move the output files to the $FAM_NSNPS directories on the HPC to run GWASs in the next step. 
```
# phenotypes with 100 causal SNPs (Fig. 2)
Rscript $CODE/01a_pheno_sim_snps100.R

# phenotypes with 500 causal SNPs (Fig. 2)
Rscript $CODE/01a_pheno_sim_snps500.R

# phenotypes with 1000 causal SNPs (Fig. 2)
Rscript $CODE/01a_pheno_sim_snps1000.R

# phenotypes with differing k values (Fig. 3)
Rscript $CODE/01b_devs_pheno_sim_snps100.R
```

#### Run GWASs (HPC)
Scripts to submit an slurm array job to cluster per phenotype to perform linear and dominance GWASs.
```
# run linear GWASs
$HAP/02_run_linear_gwas.sh $FAM_NSNPS100 	    
$HAP/02_run_linear_gwas.sh $FAM_NSNPS500 			
$HAP/02_run_linear_gwas.sh $FAM_NSNPS1000
$HAP/02_run_linear_gwas.sh $DEVSIMS_NSNPS100 		# (k value sims)	

# run dom-dev GWASs
$HAP/02_run_domdev_gwas.sh $FAM_NSNPS100  			
$HAP/02_run_domdev_gwas.sh $FAM_NSNPS500  			
$HAP/02_run_domdev_gwas.sh $FAM_NSNPS1000
$HAP/02_run_domdev_gwas.sh $DEVSIMS_NSNPS100    # (k value sims)

# merge summary statistics per phenotype
$HAP/02a_merge_all_sumstats.sh
```

#### Prepare data for models (clean, compute PGSs, format for ML/DL)
Slurm job for cleaning and preping data for models. The 'dnn.job.txt' and 'xgb.job.txt' are run using a single NVIDIA A100 GPU. 
```
# submit ready data job (HPC)
sbatch 03_ready.data.ss.job.txt
```

#### Run models
Scripts to run models. 
```
# run PRS regressions (HPC)
sbatch 04_run.prs.regressions.ss.job.txt

# run DNN (HPC)
$DL/05_submit_dnn.sh

# run XGB (HPC)
$DL/06_submit_xgb.sh
```

Run scripts 07 and onward manually to create figures.

### Testing single phenotypes
To test models in a single phenotype, adjust and run the code below. ML/DL-ready data files are available in the Zenodo repo listed above. 

```
DATA=/path/to/data
CODE=/path/to/script
OUT=/path/to/output

# Test on single pheno - DNN
python $CODE/run_dnn.py \
--input_data $DATA/DATA_eur_nsnps100_h0.5_a0_d0.5_50k.txt \
--result_path $OUT/dnn_result.pkl \
--model_path $OUT/dnn_models \
--epochs 50 \
--writer /test

# Test on single pheno - XGB
python $CODE/run_xgb.py \
--input_data $DATA/DATA_eur_nsnps100_h0.5_a0_d0.5_50k.txt \
--result_path $OUT/xgb_result.pkl \
--model_path $OUT/xgb_models.pkl \
--max_depth 1 \
--learning_rate 0.01 0.3 

# Test polygenic scores on single pheno
# NOTE: Need to edit directory path 04_prs_regressions_ss.R (line 18)
Rscript "$CODE/04_prs_regressions_ss.R" \
    /path/to/DATA_eur_nsnps100_h0.5_a0_d0.5_50k.txt \
    /path/to/output/ \
```

