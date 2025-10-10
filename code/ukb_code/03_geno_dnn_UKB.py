# Script for running DNN
# Author: Nate
# Date: 13 May 2025

# import libraries
import os
import sys
import copy
import torch
import optuna
import pickle
import random
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
import torch.nn as nn
from pathlib import Path
import torch.optim as optim
from sklearn.metrics import r2_score
from torch.optim.lr_scheduler import StepLR, ReduceLROnPlateau
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from torch.utils.data import DataLoader, TensorDataset
from optuna.samplers import GridSampler
from torch.utils.tensorboard import SummaryWriter
from model_definition_UKB import DNN, DNNObjective, adjusted_r2, split_data_UKB, get_architecture_params

# Initialize argument parser
parser = argparse.ArgumentParser(description="Train a DNN model with configurable parameters.")
parser.add_argument("--input_data", type=str, required=True, help="Path to the input data file.")
parser.add_argument("--result_path", type=str, help="File path for saving results.")
parser.add_argument("--model_path", type=str, help="File path for saving trained models.")
parser.add_argument("--writer", type=str, default="/projects/0/vusr0599/dl-prs/output/model_runs/", help="File path for saving trained models.")
parser.add_argument("--epochs", type=int, default=100, help="Number of training epochs.")
parser.add_argument("--optuna_trials", type=int, default=50, help="Number of Optuna trials.")
parser.add_argument("--batch_size", type=int, default=32, help="Batch size for training.")
parser.add_argument("--learning_rates", nargs="+", type=float, default=[1e-5, 1e-4], help="Learning rates for the optimizer.")
# parser.add_argument("--hidden_dims", nargs="+", type=int, default=[1024, 512, 256, 128], help="List of hidden layer dimensions.")
# parser.add_argument("--use_dropout", nargs="+", type=bool, default=[False, True, False, True], help="Whether to use dropout for each hidden layer.")
parser.add_argument("--dropout_p", nargs="+", type=float, default=[0.4, 0.5], help="Dropout probabilities.")
parser.add_argument("--l2_values", nargs="+", type=float, default=[1e-6, 1e-5, 1e-4], help="L2 values.")
parser.add_argument("--use_batch_norm", type=bool, default=False, help="Whether to use batch normalization for each hidden layer.")
parser.add_argument("--search_type", type=str, default="grid", help="Random vs. grid search.")
args = parser.parse_args()

# Set paths
result_path = args.result_path
model_path = args.model_path
input_data = args.input_data
writer_path = args.writer

# Set hyperparameters
epochs = args.epochs
optuna_trials = args.optuna_trials
batch_size = args.batch_size
learning_rates = args.learning_rates
l2_values = args.l2_values
# hidden_dims = args.hidden_dims
# use_dropout = args.use_dropout
dropout_probs = args.dropout_p
use_batch_norm = args.use_batch_norm
search_type = args.search_type

# Starting message
print("\n", "| ---- Starting DNN script ---- |", "\n")

print(" ... input arguments ... ")
print(" > input data:", input_data)
print(" > result_path:", result_path)
print(" > model_path:", model_path)
print(" > epochs:", epochs)
print(" > batch_size:", batch_size)
print(" > learning_rates:", learning_rates)
# print(" > hidden_dims:", hidden_dims)
# print(" > use_dropout:", use_dropout)
print(" > dropout_p:", dropout_probs)
print(" > use_batch_norm:", use_batch_norm, "\n")
print(" > search_type:", search_type)

# create dataset class
class dataset():
    def __init__(self, x, y, ids):
        self.x = torch.tensor(x, device=device, dtype=torch.float32)        # MPSFloatType  torch.float32
        self.y = torch.tensor(y, device=device, dtype=torch.float32)        # .float32
        self.ids = ids
        self.length = self.x.shape[0]

    def __getitem__(self, idx):
        return self.x[idx], self.y[idx], self.ids[idx]

    def __len__(self):
        return self.length

# Set your desired seed for reproducibility
seed = 42
random.seed(seed)           # Python's built-in random module
np.random.seed(seed)        # NumPy
torch.manual_seed(seed)     # PyTorch

# Load the data
print("| ... Loading data ... |")
data = pd.read_csv(input_data, sep = "\t")
# fill NAs with 0.0
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
       'phenotype', 'additive.prs', 'addditive.dom.comp',
       'domdev.add.comp', 'domdev.dom.comp', 'domdev.prs', 'split2'], axis=1, inplace=True)

# check data
print(" ... Data info: ... ")
print("\n", data.info(), "\n")

print("| ... Performing one-hot encoding ... |")
# Step 1: Separate the first four columns (IID, split1, split3, int_phenotype) and SNP data
non_snp_columns = data[['IID', 'split1', 'split3', 'int_phenotype']]  # First three columns
snp_data = data.drop(columns=['IID', 'split1', 'split3', 'int_phenotype'])  # SNP data for encoding

# Step 2: One-hot encode SNPs
# Initialize the OneHotEncoder
encoder = OneHotEncoder(sparse_output=False, categories='auto')

# Fit and transform the SNP data
one_hot_encoded = encoder.fit_transform(snp_data)

# Convert the one-hot encoded data to a DataFrame (to facilitate merging)
one_hot_snp_data = pd.DataFrame(one_hot_encoded, index=data.index)

# Step 3: Merge the non-SNP columns back with the one-hot encoded SNP columns
# Concatenating the non-SNP columns (IID, split1, phenotype) with the one-hot encoded SNP data
data = pd.concat([non_snp_columns, one_hot_snp_data], axis=1)

# check data
print(" ... Checking data ... ", "\n")
print(" > Data info: ")
print(data.info(), "\n")
# Verify final shape
print(" > Final data shape:", data.shape)

# set device
print("| ... Setting device ... |")
if torch.cuda.is_available():
    device = torch.device("cuda")
    print(" > Using CUDA device:", device, " ... \n")
elif torch.has_mps == True:
    device = torch.device("mps")
    print(" > Using MPS device:", device, " ... \n")
else:
    raise RuntimeError(" !! CUDA device not available. Aborting to avoid running on CPU !!")

# Split data on split2 for hypterparameter tuning
print("\n", "| ... Preparing data for Optuna ... |")
## split data
X_train, y_train, X_valid, y_valid, X_test, y_test, ids_train, ids_tune, ids_test, features = split_data_UKB(data, split_col = "split3", test_split = "test", tune_split = "tune",
                                                                                            train_split = "train", drop_col = "split1", phenotype_col = "int_phenotype")

print(" > Train set shape:", X_train.shape, y_train.shape)
print(" > Validation set shape:", X_valid.shape, y_valid.shape)
print(" > Test set shape:", X_test.shape, y_test.shape)
# ready training set
trainset = dataset(X_train,y_train, ids_train)
trainloader = DataLoader(trainset,batch_size=batch_size,shuffle=True)
# ready tuning set
validset = dataset(X_valid, y_valid, ids_tune)
validloader = DataLoader(validset, batch_size=batch_size, shuffle=False)
# ready test set
testset = dataset(X_test, y_test, ids_test)
testloader = DataLoader(testset, batch_size=batch_size, shuffle=False)

# save test data
result_data = {}

# Define input and output dimensions
input_dim = X_train.shape[1]
output_dim = 1
print(" > Input dimension:", input_dim)

# define architecture
hidden_dims, use_dropout, use_batch_norm = get_architecture_params(input_dim, batch_norm = use_batch_norm)
print(" > Hidden layers:", hidden_dims)
print(" > Use dropout:", use_dropout)
print(" > Use batch norm:", use_batch_norm)
print("")

# Create Optuna objective instance
objective = DNNObjective(
    trainloader=trainloader,
    validloader=validloader,
    input_dim=input_dim,
    output_dim=output_dim, 
    epochs=epochs,
    hidden_dims=hidden_dims, 
    use_dropout=use_dropout, 
    learning_rate=learning_rates, 
    dropout_p=dropout_probs, 
    l2 = l2_values, 
    search_type=search_type
)

if search_type == "grid":
    print(" > Performing grid search ...")
    # define search space for learning rate and dropout probability
    search_space = {
        "learning_rate": learning_rates,
        "dropout_p": dropout_probs,
        "l2": l2_values,
    }
    sampler = GridSampler(search_space)
    study = optuna.create_study(direction="maximize", sampler=sampler)
    n_trials = len(learning_rates) * len(dropout_probs) * len(l2_values)
    study.optimize(objective, n_trials=n_trials)
else:
    print(" > Performing random search ...")
    # Optimize using Optuna
    study = optuna.create_study(direction="maximize")
    study.optimize(objective, n_trials=optuna_trials)

# Best hyperparameters
print("Best hyperparameters:", study.best_params)
print("Best validation R2:", study.best_value)

# set hyperparameters
learning_rate = study.best_params["learning_rate"]
dropout_p = study.best_params["dropout_p"]
l2 = study.best_params["l2"]

# Instantiate the model
model = DNN(
    input_dim=input_dim,
    output_dim=output_dim,
    hidden_dims=hidden_dims,
    use_dropout=use_dropout,
    dropout_p=dropout_p,
    use_batch_norm=use_batch_norm,
)

# Define loss function and optimizer
loss_fn = nn.MSELoss()
# Define optimizer
if l2 > 0:
    optimizer = torch.optim.AdamW(model.parameters(), lr=learning_rate, weight_decay=l2)
else:
    optimizer = torch.optim.AdamW(model.parameters(), lr=learning_rate)
# optimizer = torch.optim.SGD(model.parameters(), lr=1e-4, momentum=0.9)
scheduler = ReduceLROnPlateau(optimizer, mode='max', factor=0.5, patience=5, verbose=True)

# parallelize model and send to machine
model = nn.DataParallel(model)
model.to(device)

# Early stopping parameters
patience = 10  # Number of epochs with no improvement after which training will be stopped
best_loss = float('inf')
epochs_no_improve = 0

# create summary writer directory name
writer_dir = Path(writer_path)

# Initialize TensorBoard writer
print("\n", "... Initializing TensorBoard writer ...")
print(" > TensorBoard writer directory:", writer_dir, "\n")
writer = SummaryWriter(writer_dir)

# print 
print("Training DNN model ...")

# start model training
for epoch in tqdm(range(epochs)):
    
    # set model for backprop
    model.train()
    train_loss = 0.0
    
    # start batchs for training
    for i, (inputs, labels, ids) in enumerate(trainloader, 0):
        if torch.isnan(inputs).any():
            print("NaN detected in input")
            break 
        
        optimizer.zero_grad()
        outputs = model(inputs).squeeze()
        loss = loss_fn(outputs, labels.squeeze())
        loss.backward()
        optimizer.step()
        train_loss += loss.item()

    # step the scheduler
    # scheduler.step()
        
    # record training loss
    writer.add_scalar('Train Loss', train_loss / len(trainloader), epoch)
        
    # model in the validation set
    model.eval()
    val_loss = 0.0
    valid_preds = []
    valid_labels = []
    valid_ids = []
    valid_preds_prs = []
    
    with torch.no_grad():
        for inputs, labels, ids in validloader:
            outputs = model(inputs).squeeze()
            loss = loss_fn(outputs, labels.squeeze())
            val_loss += loss.item()
            valid_preds.append(outputs.cpu().numpy())
            valid_labels.append(labels.cpu().numpy())
            valid_ids.extend(ids.cpu().numpy())
            valid_preds_prs.extend(outputs.cpu().numpy().flatten())
        
    # Flatten lists of predictions and labels
    valid_preds = np.concatenate(valid_preds)
    valid_labels = np.concatenate(valid_labels)
    
    # performance metrics
    valid_n = len(valid_labels)
    valid_p = inputs.shape[1]
    valid_r2 = r2_score(valid_labels, valid_preds)
    valid_adj_r2 = adjusted_r2(valid_labels, valid_preds, valid_n, valid_p)

    # step the scheduler
    scheduler.step(valid_r2)
    
    # write to tensorboard 
    writer.add_scalar('Valid Loss', val_loss / len(validloader), epoch)
    writer.add_scalar("Valid R2", valid_r2, epoch)
    writer.add_scalar("Valid Adjusted-R2", valid_adj_r2, epoch)

    # Check for early stopping
    current_val_loss = val_loss / len(validloader)       
    if current_val_loss < best_loss:  # Compare with best_loss
        best_loss = current_val_loss  # Update best_loss
        epochs_no_improve = 0  # Reset patience counter
        best_model_state = copy.deepcopy(model.state_dict())  # Save model state
    else:
        epochs_no_improve += 1  # Increment patience counter
        if epochs_no_improve == patience:
            print("Early stopping!")
            break
            
# load best model from training
model.load_state_dict(best_model_state)

# run in test set
model.eval()
test_loss = 0.0
test_preds = []
test_labels = []
test_ids = []
test_preds_prs = []
                        
with torch.no_grad():
    for inputs, labels, ids in testloader:
        outputs = model(inputs).squeeze()
        loss = loss_fn(outputs, labels.squeeze())
        test_loss +=loss.item()
        
        # record results
        test_preds.append(outputs.cpu().numpy())
        test_labels.append(labels.cpu().numpy())
        test_ids.extend(ids.cpu().numpy())
        test_preds_prs.extend(outputs.cpu().numpy().flatten())
        
# Flatten lists of predictions and labels
test_preds = np.concatenate(test_preds)
test_labels = np.concatenate(test_labels)

# save raw pred probs
test_pred_probs = list(zip(test_ids, test_preds_prs))
        
# performance metrics
test_n = len(test_labels)
test_p = inputs.shape[1]
test_r2 = r2_score(test_labels, test_preds)
test_adj_r2 = adjusted_r2(test_labels, test_preds, test_n, test_p)

# print test performance metrics
print('> Test R2', test_r2)
print('\n') 

# set dict for recording
tune_dict = {"ids":valid_ids, "phenotype":valid_labels, "pred_probs":valid_preds_prs ,"adj_r2":valid_adj_r2, "r2":valid_r2}
test_dict = {"ids":test_ids, "phenotype":test_labels, "pred_probs":test_preds_prs ,"adj_r2":test_adj_r2, "r2":test_r2}

# save
result_data["tune"] = tune_dict
result_data["test"] = test_dict

# print
print("| ---- Finished DNN training ---- |", "\n")

# convert results to df
print(" ... Preparing results to save ...")
result_df = pd.DataFrame.from_dict(result_data).transpose()
# add model info for plotting later
result_df['model'] = "DNN"
result_df['snps'] = "TRUE_CAUSAL"
result_df['type'] = "CLUMP1"
# save results
result_df.to_pickle(result_path)
# print ending messages
print(" ... Results saved to: ...")
print(result_path)
print("\n")
print(" ... GoOdByE ... \n")