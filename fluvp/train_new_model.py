#!/usr/bin/python
# -*- coding: utf-8 -*-
# Author: lihuiru
# Created on 2023/11/18 20:29
from collections import Counter
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from joblib import dump

# Import the latest feature ranking list based on RFE
data = pd.read_csv('../data/max_AUC.csv')
mu = data.loc[:, 'features'].str.cat(sep=',')
lst = mu.split(',')
lst = [i.strip() for i in lst]

# Sort features based on their occurrence frequency in descending order
keys_sorted = sorted(Counter(lst).keys(), key=Counter(lst).get, reverse=True)

# Load your data
df = pd.read_csv('../data/mutation_matrix2.0.csv')

# Drop rows with any NaN values
df = df.dropna()

# Map string labels to integers
df.iloc[:, -1] = df.iloc[:, -1].map({'Virulent': 0, 'Avirulent': 1})

# Ensure the labels are integers
df.iloc[:, -1] = df.iloc[:, -1].astype(int)

# Separate features and labels
X = df.iloc[:, 1:-1]
y = df.iloc[:, -1]

# The number representing the top N features selected during cross-validation
num = 11

# Initialize a DataFrame to store results
results = pd.DataFrame(columns=['Seed', 'Accuracy', 'Precision', 'Recall', 'F1 Score', 'AUC'])
seed = 42  # Random seed for reproducibility

# Split data into training and test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=seed)

# Select Top-N features
X_train_selected = X_train[keys_sorted[:num]]
X_test_selected = X_test[keys_sorted[:num]]

# Initialize and train the model
model = RandomForestClassifier(random_state=seed)
model.fit(X_train_selected, y_train)

# Predict and evaluate
y_pred = model.predict(X_test_selected)
y_pred_proba = model.predict_proba(X_test_selected)[:, 1]

# Save the trained RandomForest model to disk
dump(model, '../model/random_forest_model.joblib')

# # Save the top N features used for training the model
# dump(keys_sorted[:num], 'model/top_features.joblib')
