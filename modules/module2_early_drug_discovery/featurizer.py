"""
Module 2: Feature Engineering & Preprocessing
Extracts Morgan Fingerprints (ECFP4 equivalent) from chemical structures
using RDKit's modern rdFingerprintGenerator API.
"""

import os
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

# Initialize modern Morgan fingerprint generator (Radius 2 = ECFP4, 1024-bit vector)
MORGAN_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)


def generate_morgan_fingerprint(smiles: str) -> np.ndarray:
    """Generates an explicit Morgan bit vector (ECFP4 equivalent) from a SMILES string.
    Returns an all-zero vector if SMILES is invalid.
    """
    if not isinstance(smiles, str) or not smiles.strip():
        return np.zeros((1024,), dtype=np.int8)

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.zeros((1024,), dtype=np.int8)

    # Directly exports to a NumPy array without deprecation warnings
    return MORGAN_GEN.GetFingerprintAsNumPy(mol)


def run_featurization():
    input_path = "data/processed/clean_bioactivity.csv"
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Missing {input_path}. Please complete Module 1 ingestion first.")

    df = pd.read_csv(input_path)
    print(f"[+] Loaded {len(df)} compounds from {input_path}")

    # Generate 1024-bit Morgan fingerprint matrix
    fingerprints = [generate_morgan_fingerprint(smiles) for smiles in df["smiles"]]

    # Feature matrix X: Shape (N_compounds, 1024)
    X = np.array(fingerprints, dtype=np.int8)

    # Target vectors
    y_reg = df["pIC50"].values.astype(np.float32)         # Continuous affinity for regression
    y_clf = df["bioactivity_class"].values.astype(np.int8) # Binary class: Active=1, Inactive=0

    os.makedirs("data/processed", exist_ok=True)
    np.save("data/processed/X_features.npy", X)
    np.save("data/processed/y_pIC50.npy", y_reg)
    np.save("data/processed/y_class.npy", y_clf)

    sparsity = (1.0 - (np.count_nonzero(X) / X.size)) * 100

    print(f"[+] Feature engineering complete.")
    print(f"    - Feature Matrix X shape : {X.shape}")
    print(f"    - Regression target y    : {y_reg.shape}")
    print(f"    - Classification target y: {y_clf.shape}")
    print(f"    - Sparsity               : {sparsity:.2f}% zeros")


if __name__ == "__main__":
    run_featurization()