"""
Module 4: Biomarker Discovery & Patient Stratification
Script: train_classifier.py
Trains an interpretable L1-regularized classifier to predict drug response
and extract a minimal companion diagnostic (CDx) gene signature.
"""

import os
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report, roc_auc_score


def run_classifier_training():
    data_path = "data/processed/stratified_patient_cohort.csv"
    if not os.path.exists(data_path):
        raise FileNotFoundError("Missing stratified patient cohort. Run patient_clustering.py first.")

    df = pd.read_csv(data_path, index_col=0)
    gene_cols = [c for c in df.columns if c not in ["PCA_1", "PCA_2", "Subtype_Cluster"]]

    # Target: 1 = Responder (Cluster 2), 0 = Non-Responder (Clusters 0 and 1)
    df["Responder_Label"] = (df["Subtype_Cluster"] == 2).astype(int)

    X = df[gene_cols]
    y = df["Responder_Label"]

    # 70/30 Stratified Split
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.30, random_state=42, stratify=y
    )

    # Train L1-regularized Logistic Regression (sparse feature selection)
    clf = LogisticRegression(penalty="l1", solver="liblinear", C=0.5, random_state=42)
    clf.fit(X_train, y_train)

    y_pred = clf.predict(X_test)
    y_prob = clf.predict_proba(X_test)[:, 1]

    auc = roc_auc_score(y_test, y_prob)

    print("--- Companion Diagnostic Classifier Performance ---")
    print(f"Test Set ROC-AUC Score: {auc:.3f}\n")
    print("Classification Report:")
    print(classification_report(y_test, y_pred, target_names=["Non-Responder", "Responder"]))

    # Inspect the learned gene biomarker weights
    weights = pd.Series(clf.coef_[0], index=gene_cols)
    selected_biomarkers = weights[weights != 0].sort_values(ascending=False)

    print("--- Learned Sparse Biomarker Signature (Coefficients) ---")
    for gene, weight in selected_biomarkers.items():
        direction = "Positive Indicator (Sensitivity)" if weight > 0 else "Negative Indicator (Resistance)"
        print(f"Gene: {gene:<8} | Weight: {weight:>6.3f} | Role: {direction}")


if __name__ == "__main__":
    run_classifier_training()