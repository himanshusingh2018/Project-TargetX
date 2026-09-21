"""
Module 5: Clinical Trial Optimization
Script: dropout_predictor.py
Predicts high-risk patient dropout during clinical trials using baseline clinical
and operational friction features to enable proactive retention interventions.
"""

import os
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import roc_auc_score, classification_report


def generate_trial_monitoring_data(n_enrolled: int = 150, seed: int = 42) -> pd.DataFrame:
    """Generates synthetic baseline trial monitoring metrics for enrolled patients."""
    np.random.seed(seed)
    
    patient_ids = [f"ENROLL-{i:03d}" for i in range(1, n_enrolled + 1)]
    # Distance from patient residence to trial site in miles
    distance_miles = np.random.exponential(scale=25, size=n_enrolled) + 2
    # Baseline mild toxicities / adverse events (Grade 1/2 fatigue, nausea)
    baseline_ae_count = np.random.poisson(lam=1.2, size=n_enrolled)
    # Missed prior scheduled appointments in pre-screening
    missed_pre_visits = np.random.binomial(n=3, p=0.15, size=n_enrolled)
    # Patient age
    age = np.random.randint(35, 80, size=n_enrolled)
    # Caregiver support available (1 = yes, 0 = lives alone)
    has_caregiver = np.random.choice([1, 0], size=n_enrolled, p=[0.70, 0.30])

    # Dropout probability logit driven by operational friction
    logit = (
        -2.5 
        + 0.03 * distance_miles 
        + 0.45 * baseline_ae_count 
        + 0.80 * missed_pre_visits 
        - 0.60 * has_caregiver
    )
    prob_dropout = 1 / (1 + np.exp(-logit))
    dropped_out = (np.random.rand(n_enrolled) < prob_dropout).astype(int)

    df = pd.DataFrame({
        "Patient_ID": patient_ids,
        "Distance_Miles": np.round(distance_miles, 1),
        "Baseline_AE_Count": baseline_ae_count,
        "Missed_Pre_Visits": missed_pre_visits,
        "Age": age,
        "Has_Caregiver": has_caregiver,
        "Early_Dropout": dropped_out
    })
    return df


def train_dropout_model():
    df = generate_trial_monitoring_data(n_enrolled=200)
    
    feature_cols = ["Distance_Miles", "Baseline_AE_Count", "Missed_Pre_Visits", "Age", "Has_Caregiver"]
    X = df[feature_cols]
    y = df["Early_Dropout"]

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.30, random_state=42, stratify=y
    )

    model = GradientBoostingClassifier(n_estimators=60, max_depth=3, random_state=42)
    model.fit(X_train, y_train)

    y_prob = model.predict_proba(X_test)[:, 1]
    y_pred = (y_prob >= 0.50).astype(int)

    auc = roc_auc_score(y_test, y_prob)

    print("--- Clinical Trial Dropout Risk Model ---")
    print(f"Cohort Size: {len(df)} patients | Dropout Rate: {y.mean()*100:.1f}%")
    print(f"Test Set ROC-AUC Score: {auc:.3f}\n")
    print("Classification Report:")
    print(classification_report(y_test, y_pred, target_names=["Retained", "Dropped Out"]))

    # Feature Importance for Operations
    importances = pd.Series(model.feature_importances_, index=feature_cols).sort_values(ascending=False)
    print("--- Operational Drivers of Patient Dropout (Feature Importances) ---")
    for feat, imp in importances.items():
        print(f"{feat:<20}: {imp*100:.1f}%")


if __name__ == "__main__":
    os.makedirs("data/processed", exist_ok=True)
    train_dropout_model()