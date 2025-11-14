# Repository for "Benchmarking non-additive genetic effects on polygenic prediction and machine learning-based approaches"

DOI: https://doi.org/10.1101/2025.10.10.25337750

Repository contains code, data (simulated phenotypes), and results (from analyses using simulated phenotypes) from the pre print listed above. Code is run script by script (to confirm successful execution) and is mixed between running on a local PC (Mac Silicon M1 Pro; MacOS Sonoma 14.5) and an HPC (Snellius). UK Biobank data can be accessed on their platform for approved researchers (https://www.ukbiobank.ac.uk/use-our-data/research-analysis-platform/). The Python package versions used can be found in the '/code/required_packages_...txt' files.

Raw simulated genotype data, UK Biobank GWAS summary statistics, and data files processed for ML/DL models can be found at: 
https://zenodo.org/records/17552313

#### Setup (local)
Code to set up virtual environment with package versions used for local code (MacOS). 
```
brew install pyenv pyenv-virtualenv
pyenv install 3.10.12
pyenv virtualenv 3.10.12 myenv-3.10
pyenv activate myenv-3.10
pip install -r /path/to/required_packages_local.txt
```

#### Setup (HPC)
(1) If needed, install miniconda. More information for installing miniconda can be found here:
https://www.anaconda.com/docs/getting-started/miniconda/install#linux-terminal-installer 
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

(2) Create virtual environment with dependencies
```
source $HOME/miniconda3/etc/profile.d/conda.sh
conda env create -f /path/to/non-add-pgs-hpc-env.yml
```

#### Generate simulated phenotypes (local)
Executable R scripts that generate synthetic phenotypes, PLINK .fam files, and performs QC checks. Need a config file defining phenotype parameters (see '/data/fams_fin_snps100/sim.fin.config.snps100.csv' for example).  
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

