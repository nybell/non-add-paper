# Script with model definitions
# Author: Nate
# Date: 14 Jan 2025

# Load packages
import torch
import optuna
import numpy as np
import pandas as pd
import xgboost as xgb
import torch.nn as nn
from tqdm import tqdm
from sklearn.metrics import r2_score
from torch.optim.lr_scheduler import StepLR

class DNN(nn.Module):
    def __init__(self, input_dim, output_dim, hidden_dims=[128, 64], use_dropout=None, dropout_p=0.1, use_batch_norm=None):
        """
        Initialize the DNN class.

        Parameters:
        - input_dim (int): Number of input features.
        - output_dim (int): Number of output features.
        - hidden_dims (list): List specifying the number of neurons in each hidden layer.
        - use_dropout (list or None): List of booleans specifying dropout for each layer, or None to disable dropout.
        - dropout_p (float): Dropout probability (applied where use_dropout is True).
        - use_batch_norm (list or None): List of booleans specifying batch norm for each layer, or None to disable batch norm.
        """
        super(DNN, self).__init__()
        layers = []  # Initialize an empty list to hold the layers of the network
        current_dim = input_dim  # Set the current dimension to the input size

        # Validate the lengths of use_dropout and use_batch_norm
        if use_dropout is not None and len(use_dropout) != len(hidden_dims):
            raise ValueError("Length of use_dropout must match length of hidden_dims.")
        if use_batch_norm is not None and len(use_batch_norm) != len(hidden_dims):
            raise ValueError("Length of use_batch_norm must match length of hidden_dims.")

        # Create hidden layers
        for i, h_dim in enumerate(hidden_dims):
            # Add a fully connected (linear) layer
            layers.append(nn.Linear(current_dim, h_dim))

            # Optionally add batch normalization for this layer
            if use_batch_norm and use_batch_norm[i]:
                layers.append(nn.BatchNorm1d(h_dim))

            # Add activation function
            layers.append(nn.ELU())

            # Optionally add dropout for this layer
            if use_dropout and use_dropout[i]:
                layers.append(nn.Dropout(p=dropout_p))

            # Update current dimension for the next layer
            current_dim = h_dim

        # Add the output layer (final fully connected layer)
        layers.append(nn.Linear(current_dim, output_dim))

        # Combine all layers into a single sequential model
        self.network = nn.Sequential(*layers)

    def forward(self, x):
        """
        Forward pass of the model.

        Parameters:
        - x (torch.Tensor): Input tensor.

        Returns:
        - torch.Tensor: Output tensor.
        """
        return self.network(x)
    
# Define the DNN objective class    
class DNNObjective:
    def __init__(self, trainloader, validloader, input_dim, output_dim, epochs, 
                 hidden_dims=[1024, 512, 256, 128], use_dropout=[False, True, False, True], 
                 learning_rate=[1e-4, 3e-3], dropout_p=[0.1, 0.3], l2=[1e-6, 1e-5, 1e-4], search_type="random"):  
        """
        Optuna objective class for tuning a DNN model.

        Parameters:
        - trainloader: Training DataLoader object.
        - validloader: Validation DataLoader object.
        - input_dim (int): Number of input features.
        - output_dim (int): Number of output features.
        - hidden_dims (list): List specifying the number of neurons in each hidden layer.
        - use_dropout (list): Dropout settings for each layer.
        """
        self.trainloader = trainloader
        self.validloader = validloader
        self.input_dim = input_dim
        self.output_dim = output_dim
        self.epochs = epochs
        self.hidden_dims = hidden_dims
        self.use_dropout = use_dropout
        self.learning_rate = learning_rate
        self.dropout_p = dropout_p
        self.l2 = l2
        self.search_type = search_type

    def __call__(self, trial):
        # Sample hyperparameters
        if self.search_type == "grid":
            learning_rate = trial.suggest_categorical("learning_rate", self.learning_rate)
            dropout_p = trial.suggest_categorical("dropout_p", self.dropout_p)
            l2 = trial.suggest_categorical("l2", self.l2)
        else:
            # Define the search space for each hyperparameter
            learning_rate = trial.suggest_float('learning_rate', min(self.learning_rate), max(self.learning_rate), log=True)  
            dropout_p = trial.suggest_float('dropout_p', min(self.dropout_p), max(self.dropout_p))  
            l2 = trial.suggest_float('l2', min(self.l2), max(self.l2), log=True)  

        # set device
        print("\n| ... Setting device ... |")
        if torch.cuda.is_available():
            device = torch.device("cuda")
            print(" > Using CUDA device:", device, " ... ")
        elif torch.has_mps == True:
            device = torch.device("mps")
            print(" > Using MPS device:", device, " ... ")
        else:
            raise RuntimeError(" !! CUDA device not available. Aborting to avoid running on CPU !!")

        # print("... Initializing model ...")
        # Initialize the model
        model = DNN(
            input_dim=self.input_dim,
            output_dim=self.output_dim,
            hidden_dims=self.hidden_dims,
            use_dropout=self.use_dropout,
            dropout_p=dropout_p
        ).to(device)

        # Define loss function and optimizer
        loss_fn = nn.MSELoss()
        # Define optimizer
        if l2 > 0:
            optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate, weight_decay=l2)
        else:
            optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)

        # Early stopping parameters
        patience = 20
        best_loss = float('inf')
        epochs_no_improve = 0

        print("... Optimizing hyper-parameters with Optuna ...")
        epochs = self.epochs
        for epoch in tqdm(range(epochs)):
            model.train()
            train_loss = 0.0
            for inputs, labels, _ in self.trainloader:
                optimizer.zero_grad()
                outputs = model(inputs.to(device)).squeeze()
                loss = loss_fn(outputs, labels.to(device).squeeze())
                loss.backward()
                optimizer.step()
                train_loss += loss.item()

            # Evaluate the model
            model.eval()
            val_loss = 0.0
            valid_preds = []
            valid_labels = []
            with torch.no_grad():
                for inputs, labels, _ in self.validloader:
                    outputs = model(inputs.to(device)).squeeze()
                    loss = loss_fn(outputs, labels.to(device).squeeze())
                    val_loss += loss.item()
                    valid_preds.append(outputs.cpu().numpy())
                    valid_labels.append(labels.cpu().numpy())

            # Check for early stopping
            current_val_loss = val_loss / len(self.validloader)       
            if current_val_loss < best_loss:  # Compare with best_loss
                best_loss = current_val_loss  # Update best_loss
                epochs_no_improve = 0  # Reset patience counter
            else:
                epochs_no_improve += 1  # Increment patience counter
                if epochs_no_improve == patience:
                    print("Early stopping!")
                    break
        # Flatten predictions and labels
        valid_preds = np.concatenate(valid_preds)
        valid_labels = np.concatenate(valid_labels)
        # Return validation loss
        valid_r2 = r2_score(valid_labels, valid_preds)
        # return val_loss / len(self.validloader)
        return valid_r2

# function to get architecture-specific parameters
def get_architecture_params(input_dim, batch_norm=False):
    """
    Returns architecture-specific hidden_dims, use_dropout, and use_batch_norm lists
    based on the input dimensionality.
    """
    if input_dim > 10000:
        hidden_dims = [4096, 1024, 512, 64]
    elif input_dim > 5000:
        hidden_dims = [4096, 1024, 512, 64]
    elif input_dim < 100:
        hidden_dims = [128, 64, 32]
    else:
        hidden_dims = [4096, 1024, 512, 64]

    # Define dropout: apply to all but the last layer if width is large
    use_dropout = [h > 256 for h in hidden_dims]

    if batch_norm:
        # Define batch norm: enable for all but the last layer
        use_batch_norm = [True] + (len(hidden_dims)-1) * [False]
    else:
        # Disable batch norm
        use_batch_norm = [False] * len(hidden_dims)

    return hidden_dims, use_dropout, use_batch_norm


# Define function for computing adjusted R-squared
def adjusted_r2(y_true, y_pred, n, p):
    # Calculate R-squared
    r2 = r2_score(y_true, y_pred)
    # Calculate Adjusted R-squared
    adj_r2 = 1 - (1 - r2) * ((n - 1) / (n - p - 1))
    return adj_r2

def xgb_objective(trial, X_train, y_train, X_valid, y_valid,
                  estimators=[5000, 8000], max_depth=[1, 2],
                  learning_rate=[0.01, 0.2], xgb_gpu=False):

    # Suggest hyperparameters
    param = {
        'max_depth': trial.suggest_int('max_depth', min(max_depth), max(max_depth)),
        'learning_rate': trial.suggest_float('learning_rate', min(learning_rate), max(learning_rate)),
        'min_child_weight': trial.suggest_int('min_child_weight', 1, 10),
        'colsample_bytree': trial.suggest_float('colsample_bytree', 0.6, 1.0),
        'reg_lambda': trial.suggest_float('reg_lambda', 0.0, 10.0),
        'objective': 'reg:squarederror',
        'eval_metric': 'rmse',
        'random_state': 666
    }

    # GPU / CPU logic
    if xgb_gpu:
        param.update({
            'tree_method': 'hist',
            'device': 'cuda'
        })
    else:
        param.update({
            'tree_method': 'hist'
        })

    if trial.number == 0:
        print(f"XGBoost training parameters: {param}")
        print(" ... xgb_gpu = ", xgb_gpu, " ... ")
        print("... Initializing model ...")

    # Convert to DMatrix
    dtrain = xgb.DMatrix(X_train, label=y_train)
    dvalid = xgb.DMatrix(X_valid, label=y_valid)
    evals = [(dvalid, "eval")]

    # Get boosting rounds
    num_boost_round = max(estimators)        # trial.suggest_int('n_estimators', min(estimators), max(estimators), step=500)

    # Train the model with early stopping
    model = xgb.train(
        params=param,
        dtrain=dtrain,
        num_boost_round=num_boost_round,
        evals=evals,
        early_stopping_rounds=100,
        verbose_eval=False
    )

    # Predict on validation set
    y_pred = model.predict(dvalid)

    return r2_score(y_valid, y_pred)


# Define function for defining test splits on each fold
def split_data_5fold(data, split_col, test_split, drop_col, phenotype_col = "int_phenotype"):
    required_cols = {'IID', split_col, phenotype_col, drop_col}
    missing = required_cols - set(data.columns)
    if missing:
        raise ValueError(f"Missing required column(s): {missing}")
    
    # drop columns that are not needed
    data = data.drop([drop_col], axis=1)

    # test split
    test = data[data[split_col] == test_split]
    ids_test = test["IID"].to_list()
    test = test.drop(["IID", split_col], axis=1)
    
    # format for model
    X_test = test.drop(columns=[phenotype_col]).to_numpy()
    y_test = test[[phenotype_col]].to_numpy()
    
    # tuning set
    tune_split = (test_split % 5) + 1
    tune = data[data[split_col] == tune_split]
    ids_tune = tune["IID"].to_list()
    tune = tune.drop(["IID", split_col], axis=1)
    
    # format for model
    X_valid = tune.drop(columns=[phenotype_col]).to_numpy()
    y_valid = tune[[phenotype_col]].to_numpy()
    
    # training split
    train = data[(data[split_col] != test_split) & (data[split_col] != tune_split)]
    ids_train = train["IID"].to_list()
    train = train.drop(["IID", split_col], axis=1)
    
    # format for model
    X_train = train.drop(columns=[phenotype_col]).to_numpy()
    y_train = train[[phenotype_col]].to_numpy()

    return X_train, y_train, X_valid, y_valid, X_test, y_test, ids_train, ids_tune, ids_test, test_split

# Define function for defining test splits on each fold
def split_data_UKB(data, split_col, test_split = "test", tune_split = "tune", train_split = "train", drop_col = "split1", phenotype_col = "int_phenotype"):
    required_cols = {'IID', split_col, phenotype_col, drop_col}
    missing = required_cols - set(data.columns)
    if missing:
        raise ValueError(f"Missing required column(s): {missing}")
    
    # drop columns that are not needed
    data = data.drop([drop_col], axis=1)

    # test split
    test = data[data[split_col] == test_split]
    ids_test = test["IID"].to_list()
    test = test.drop(["IID", split_col], axis=1)
    
    # format for model
    X_test = test.drop(columns=[phenotype_col]).to_numpy()
    y_test = test[[phenotype_col]].to_numpy()
    
    # tuning set
    tune = data[data[split_col] == tune_split]
    ids_tune = tune["IID"].to_list()
    tune = tune.drop(["IID", split_col], axis=1)
    
    # format for model
    X_valid = tune.drop(columns=[phenotype_col]).to_numpy()
    y_valid = tune[[phenotype_col]].to_numpy()
    
    # training split
    train = data[data[split_col] == train_split]
    ids_train = train["IID"].to_list()
    train = train.drop(["IID", split_col], axis=1)
    features = train.columns.tolist()
    
    # format for model
    X_train = train.drop(columns=[phenotype_col]).to_numpy()
    y_train = train[[phenotype_col]].to_numpy()

    return X_train, y_train, X_valid, y_valid, X_test, y_test, ids_train, ids_tune, ids_test, features
