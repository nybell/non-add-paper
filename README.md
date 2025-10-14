# Repository for "Benchmarking non-additive genetic effects on polygenic prediction and machine learning-based approaches"

DOI: https://doi.org/10.1101/2025.10.10.25337750

Repository contains code, data (simulated phenotypes), and results (from analyses using simulated phenotypes) from the pre print listed above. Code is run script by script (to confirm successful execution) and is mixed between running on a local PC (Mac Silicon) and an HPC (Snellius). UK Biobank data can be accessed on their platform for approved researchers (https://www.ukbiobank.ac.uk/use-our-data/research-analysis-platform/). The Python package versions used can be found in the '/code/required_packages_...txt' files.

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
Slurm job for cleaning and preping data for models. The 'dnn.sh' and 'xgb.sh' are run using a single NVIDIA A100 GPU. 
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



