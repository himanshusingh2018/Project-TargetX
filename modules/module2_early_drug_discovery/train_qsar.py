"""
Module 2: Early Drug Discovery
Script: train_qsar.py
Trains a Random Forest QSAR regressor using scaffold-split cross-validation.
Evaluates generalizability on out-of-scaffold chemical series.
"""

import os
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error


def run_training():
    clean_csv = "data/processed/clean_bioactivity.csv"
    X_path = "data/processed/X_features.npy"
    y_path = "data/processed/y_pIC50.npy"

    if not all(os.path.exists(p) for p in [clean_csv, X_path, y_path]):
        raise FileNotFoundError("Missing processed inputs. Make sure clean_bioactivity.csv and .npy files exist.")

    df = pd.read_csv(clean_csv)
    X = np.load(X_path)
    y = np.load(y_path)

    # Scaffold split indices from our prior step:
    # Train: MOL001, MOL002, MOL003, MOL005 -> indices [0, 1, 2, 3]
    # Test:  MOL006, MOL007                 -> indices [4, 5]
    train_idx = [0, 1, 2, 3]
    test_idx = [4, 5]

    X_train, X_test = X[train_idx], X[test_idx]
    y_train, y_test = y[train_idx], y[test_idx]

    print(f"[+] Module 2: Training QSAR Baseline")
    print(f"[+] Training Set: {len(X_train)} samples | Test Set: {len(X_test)} samples")

    # Fit Random Forest Regressor
    model = RandomForestRegressor(n_estimators=100, random_state=42, max_depth=5)
    model.fit(X_train, y_train)

    # Generate predictions
    y_train_pred = model.predict(X_train)
    y_test_pred = model.predict(X_test)

    # Calculate errors
    train_mse = mean_squared_error(y_train, y_train_pred)
    test_mse = mean_squared_error(y_test, y_test_pred)

    print("\n--- Model Performance Evaluation ---")
    print(f"Train MSE : {train_mse:.3f}")
    print(f"Test MSE  : {test_mse:.3f}")

    print("\n--- Out-of-Scaffold Test Predictions ---")
    for i, orig_idx in enumerate(test_idx):
        mol_id = df.iloc[orig_idx]["molecule_id"]
        actual = y_test[i]
        pred = y_test_pred[i]
        delta = abs(actual - pred)
        print(f"MOL: {mol_id} | Actual pIC50: {actual:.2f} | Predicted: {pred:.2f} | Error: {delta:.2f}")

    # Top Feature Importances (Bits that drive the binding prediction)
    importances = model.feature_importances_
    top_bits = np.argsort(importances)[::-1][:5]
    print(f"\n[+] Top 5 Discriminative Fingerprint Bits: {list(top_bits)}")


if __name__ == "__main__":
    run_training()