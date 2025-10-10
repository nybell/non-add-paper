# Script for running XGB
# Author: Nate
# Date: 12 Jan 2025

# load packages - basic
import math
import copy
import pickle
import random
import optuna
import argparse
import numpy as np
import pandas as pd
import xgboost as xgb
from tqdm import tqdm
from sklearn.metrics import r2_score
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from model_definition_UKB import xgb_objective, split_data_UKB, adjusted_r2 # type: ignore # xgb_objective

# Initialize argument parser
parser = argparse.ArgumentParser(description="Train a XGB model with configurable parameters.")
parser.add_argument("--input_data", type=str, required=True, help="Path to the input data file.")
parser.add_argument("--result_path", type=str, help="File path for saving results.")
parser.add_argument("--model_path", type=str, help="File path for saving trained models.")
parser.add_argument("--estimators", nargs="+", type=int, default=[8000,1000], help="Min and max number of estimators")
parser.add_argument("--max_depth", nargs="+", type=int, default=[1, 5], help="Min and max values for max_depth")
parser.add_argument("--learning_rate", nargs="+", type=float, default=[0.01, 0.2], help="Min and max values for learning rate")
parser.add_argument("--hyper_trials", type=int, default=200, help="Number of hyperparameter tuning trials.")
parser.add_argument("--gpu", action='store_true', help="Use GPU for training.")
args = parser.parse_args()

# Set paths
input_data = args.input_data
result_path = args.result_path
model_path = args.model_path

# Set hyperparameters
hyperparam_trials = args.hyper_trials
estimator_vals = args.estimators
max_depth_vals = args.max_depth
learning_rates = args.learning_rate
xgb_gpu = args.gpu

# Starting message
print("\n", "| ---- Starting XGB script ---- |", "\n")

print(" ... input arguments ... ")
print(" > input data:", input_data)
print(" > result_path:", result_path)
print(" > model_path:", model_path)
print(" > estimators:", args.estimators)
print(" > max_depth:", args.max_depth)
print(" > learning_rate:", args.learning_rate)

# Set your desired seed for reproducibility
seed = 42
random.seed(seed)           # Python's built-in random module
np.random.seed(seed)        # NumPy

# Load the data
print("\n", "| ... Loading data ... |")
data = pd.read_csv(input_data, sep = "\t")
# fill NAs with -1
data = data.fillna(0.0)
# drop columns that are not needed
try:
    # Try to drop all the specified columns
    data.drop(['FID', 'array01', 'sex', 'age', 'pop_pc1',
       'pop_pc2', 'pop_pc3', 'pop_pc4', 'pop_pc5', 'pop_pc6', 'pop_pc7',
       'pop_pc8', 'pop_pc9', 'pop_pc10', 'log_phenotype', 
       'phenotype', 'additive.prs', 'addditive.dom.comp',
       'domdev.add.comp', 'domdev.dom.comp', 'domdev.prs', 'split2'], axis=1, inplace=True)
except KeyError:
    # If there's an error, drop the same columns excluding 'domdev.add.comp'
    data.drop(['FID', 'array01', 'sex', 'age', 'pop_pc1',
       'pop_pc2', 'pop_pc3', 'pop_pc4', 'pop_pc5', 'pop_pc6', 'pop_pc7',
       'pop_pc8', 'pop_pc9', 'pop_pc10', 'log_phenotype', 
       'phenotype', 'additive.prs', 'domdev.add.comp',
       'domdev.dom.comp', 'domdev.prs', 'split2'], axis=1, inplace=True)

# check data
print(" ... Data info: ... ")
print(data.info())

# Step 1: Separate the first three columns (IID, split1, phenotype) and SNP data
non_snp_columns = data[['IID', 'split1', 'split3', 'int_phenotype']]  # First three columns
snp_data = data.drop(columns=['IID', 'split1', 'split3', 'int_phenotype'])  # SNP data for encoding

# Step 2: Identify binary and multi-class columns
binary_columns = snp_data.columns[snp_data.nunique() == 2]  # Columns with two unique values
multi_class_columns = snp_data.columns[snp_data.nunique() == 3]  # Columns with three unique values

# Step 3: Retain binary columns as-is
binary_data = snp_data[binary_columns]  # Keep binary columns unchanged

# Step 4: One-hot encode only the multi-class columns
# Initialize the OneHotEncoder
encoder = OneHotEncoder(sparse_output=False, categories='auto')

# Fit and transform the multi-class columns
one_hot_encoded = encoder.fit_transform(snp_data[multi_class_columns])

# Convert the one-hot encoded data to a DataFrame
one_hot_snp_data = pd.DataFrame(one_hot_encoded, index=snp_data.index)

# Combine binary columns and one-hot encoded data
encoded_snp_data = pd.concat([binary_data, one_hot_snp_data], axis=1)

# Step 5: Merge the non-SNP columns back with the SNP data
data = pd.concat([non_snp_columns, encoded_snp_data], axis=1)

# Verify final shape
print(" > Final data shape:", data.shape)

# initialize test_data
result_data = {}

## split data
X_train, y_train, X_valid, y_valid, X_test, y_test, ids_train, ids_tune, ids_test, features = split_data_UKB(data, split_col = "split3", test_split = "test", tune_split = "tune",
                                                                                            train_split = "train", drop_col = "split1", phenotype_col = "int_phenotype")

# print current split
print("| ... Running XGBoost Model ... |")
test_split = "test"

# Use Optuna to find the best hyperparameters
study = optuna.create_study(direction='maximize')
print(f" ... Hyperparameter tuning ... \n")
with tqdm(total=hyperparam_trials) as pbar:
    def wrapped_objective(trial):
        result = xgb_objective(trial, X_train, y_train, X_valid, y_valid, 
                                estimators = estimator_vals, max_depth = max_depth_vals, learning_rate = learning_rates, xgb_gpu = xgb_gpu)
        pbar.update(1)  # Update progress bar after each trial
        return result

study.optimize(wrapped_objective, n_trials=hyperparam_trials)

# Best parameters from Optuna
best_params = study.best_params
print(f"Best hyper-parameters: {best_params}")

# check for gpu 
if xgb_gpu:
    best_params.update({
        'tree_method': 'hist',   
        'device': 'cuda'})
    print(f"... Using GPU for final training ...")
else:
    best_params.update({
        'tree_method': 'hist',})  
    print(f"... Using CPU for final training ...")
# initialize
print("... Initializing model ...")

# Convert to DMatrix
dtrain = xgb.DMatrix(X_train, label=y_train)
dvalid = xgb.DMatrix(X_valid, label=y_valid)
dtest = xgb.DMatrix(X_test, label=y_test)

# Set fixed num_boost_round for final training
final_boost_rounds = max(estimator_vals)  # or any large value

# Train with early stopping
booster = xgb.train(
    params=best_params,
    dtrain=dtrain,
    num_boost_round=final_boost_rounds,
    evals=[(dvalid, 'eval')],
    early_stopping_rounds=100,
    verbose_eval=False
)

# Predict on validation and test sets
y_pred_tune = booster.predict(dvalid)
y_pred_test = booster.predict(dtest)

# Feature importance (example: gain-based)
feat_importances = booster.get_score(importance_type='gain')

# Evaluate the model on the test set
tune_r2 = r2_score(y_valid, y_pred_tune)
test_r2 = r2_score(y_test, y_pred_test)
tune_adj_r2 = adjusted_r2(y_valid, y_pred_tune, len(y_valid), X_valid.shape[1])
test_adj_r2 = adjusted_r2(y_test, y_pred_test, len(y_test), X_test.shape[1])
print(f"XGBoost Tune Set: R2 = {tune_r2}, adj-R2 = {tune_adj_r2}")
print(f"XGBoost Test Set: R2 = {test_r2}, adj-R2 = {test_adj_r2}")

# Record results
tune_dict = {"ids": ids_tune, "phenotype": y_valid, "r2": tune_r2, "adj_r2": tune_adj_r2, "feature_importances": "NA"}
test_dict = {"ids": ids_test, "phenotype": y_test, "r2": test_r2, "adj_r2": test_adj_r2, "feature_importances": feat_importances}
result_data["tune"] = tune_dict
result_data["test"] = test_dict

print("| ---- Finished XGBoost Model ---- |", "\n")

# convert results to df
print(" ... Preparing results to save ...")
result_df = pd.DataFrame.from_dict(result_data).transpose()
# add model info for plotting later
result_df['model'] = "XGB"
result_df['snps'] = "TRUE_CAUSAL"
result_df['type'] = "GENO_MODEL"
# save results
result_df.to_pickle(result_path)
booster.save_model(model_path)
# print ending messages
print(" ... Results saved to: ...")
print(result_path)
print("\n")
print(" ... GoOdByE ... \n")


