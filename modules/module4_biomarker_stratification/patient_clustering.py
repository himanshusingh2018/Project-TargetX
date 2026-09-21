"""
Module 4: Biomarker Discovery & Patient Stratification
Script: patient_clustering.py
Performs unsupervised transcriptomic clustering on patient cohort data
to discover molecular subtypes and predict drug response phenotypes.
"""

import os
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans


def generate_cohort_transcriptomics(n_patients: int = 120, seed: int = 42) -> pd.DataFrame:
    """Generates a biologically structured synthetic RNA-seq log2(TPM) dataset.
    Simulates three distinct cancer patient subtypes:
    - Subtype A: EGFR/ERBB2-driven (High EGFR/ERBB2, intact PTEN) -> Potential Responders
    - Subtype B: KRAS/MAPK-driven (High KRAS/BRAF/MYC)            -> Primary Non-Responders
    - Subtype C: PI3K/PTEN-loss (PTEN null, High AKT1/PIK3CA)     -> Bypass Resistant
    """
    np.random.seed(seed)
    n_per_group = n_patients // 3

    # Genes evaluated in our Module 3 Knowledge Graph
    genes = ["EGFR", "ERBB2", "GRB2", "SOS1", "KRAS", "BRAF", "MAPK1", "PIK3CA", "PTEN", "AKT1", "MTOR", "MYC"]

    # Subtype A: EGFR-driven cohort
    group_a = np.random.normal(loc=5.0, scale=0.8, size=(n_per_group, len(genes)))
    group_a[:, genes.index("EGFR")] += 3.5
    group_a[:, genes.index("ERBB2")] += 2.8
    group_a[:, genes.index("PTEN")] += 1.5

    # Subtype B: KRAS-driven cohort
    group_b = np.random.normal(loc=5.0, scale=0.8, size=(n_per_group, len(genes)))
    group_b[:, genes.index("KRAS")] += 3.8
    group_b[:, genes.index("BRAF")] += 2.5
    group_b[:, genes.index("MYC")] += 3.0

    # Subtype C: PTEN-loss / AKT survival cohort
    group_c = np.random.normal(loc=5.0, scale=0.8, size=(n_per_group, len(genes)))
    group_c[:, genes.index("PTEN")] -= 3.0
    group_c[:, genes.index("PIK3CA")] += 3.0
    group_c[:, genes.index("AKT1")] += 3.2

    X_raw = np.vstack([group_a, group_b, group_c])
    patient_ids = [f"PT-{i:03d}" for i in range(1, n_patients + 1)]

    df = pd.DataFrame(X_raw, columns=genes, index=patient_ids)
    df.index.name = "Patient_ID"
    return df


def stratify_patients(df: pd.DataFrame, n_clusters: int = 3):
    """Normalizes expression, performs PCA dimensionality reduction,
    and clusters patients into distinct molecular cohorts.
    """
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(df)

    pca = PCA(n_components=2, random_state=42)
    X_pca = pca.fit_transform(X_scaled)

    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    cluster_labels = kmeans.fit_predict(X_scaled)

    strat_df = df.copy()
    strat_df["PCA_1"] = np.round(X_pca[:, 0], 2)
    strat_df["PCA_2"] = np.round(X_pca[:, 1], 2)
    strat_df["Subtype_Cluster"] = cluster_labels

    # Determine signature biomarker profiles per cluster
    profile = strat_df.groupby("Subtype_Cluster")[df.columns].mean().round(2)

    return strat_df, profile, pca.explained_variance_ratio_


if __name__ == "__main__":
    os.makedirs("data/processed", exist_ok=True)
    cohort_df = generate_cohort_transcriptomics(n_patients=120)
    raw_cohort_path = "data/raw/patient_rnaseq_cohort.csv"
    cohort_df.to_csv(raw_cohort_path)
    print(f"[+] Cohort expression data generated -> {raw_cohort_path} ({cohort_df.shape})")

    stratified_df, cluster_profiles, var_exp = stratify_patients(cohort_df, n_clusters=3)
    strat_path = "data/processed/stratified_patient_cohort.csv"
    stratified_df.to_csv(strat_path)

    print(f"[+] PCA Explained Variance: PC1={var_exp[0]*100:.1f}%, PC2={var_exp[1]*100:.1f}%")
    print("\n--- Molecular Subtype Cluster Profiles (Mean Log2 Expression) ---")
    print(cluster_profiles[["EGFR", "ERBB2", "KRAS", "PTEN", "AKT1", "MYC"]])