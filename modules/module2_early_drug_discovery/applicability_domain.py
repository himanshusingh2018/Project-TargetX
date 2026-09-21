"""
Module 2: Early Drug Discovery
Script: applicability_domain.py
Computes chemical similarity (Tanimoto distance) to define the model's
Applicability Domain (AD) and flag out-of-distribution predictions.
"""

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

MORGAN_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)


def smiles_to_bitvect(smiles: str):
    """Generates an RDKit ExplicitBitVect for fast Tanimoto computation."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return MORGAN_GEN.GetFingerprint(mol)


def compute_applicability_domain(train_smiles: list[str], query_smiles: list[str], threshold: float = 0.40) -> pd.DataFrame:
    """Calculates max Tanimoto similarity between query molecules and training set.
    Flags compounds falling outside the Applicability Domain.
    """
    train_fps = [smiles_to_bitvect(s) for s in train_smiles if smiles_to_bitvect(s) is not None]

    results = []
    for q_smi in query_smiles:
        q_fp = smiles_to_bitvect(q_smi)
        if q_fp is None:
            results.append({"smiles": q_smi, "max_tanimoto": 0.0, "in_domain": False})
            continue

        # Compute pairwise Tanimoto similarity against every training molecule
        similarities = DataStructs.BulkTanimotoSimilarity(q_fp, train_fps)
        max_sim = max(similarities) if similarities else 0.0

        results.append({
            "smiles": q_smi,
            "max_tanimoto": round(float(max_sim), 3),
            "in_domain": bool(max_sim >= threshold)
        })

    return pd.DataFrame(results)


if __name__ == "__main__":
    clean_csv = "data/processed/clean_bioactivity.csv"
    df = pd.read_csv(clean_csv)

    # Train: MOL001, MOL002, MOL003, MOL005 (indices 0, 1, 2, 3)
    # Test:  MOL006, MOL007                 (indices 4, 5)
    train_smiles = df.iloc[[0, 1, 2, 3]]["smiles"].tolist()
    test_df = df.iloc[[4, 5]].copy()

    ad_df = compute_applicability_domain(train_smiles, test_df["smiles"].tolist(), threshold=0.40)
    test_df["max_train_tanimoto"] = ad_df["max_tanimoto"].values
    test_df["in_applicability_domain"] = ad_df["in_domain"].values

    print("--- Applicability Domain Assessment for Test Compounds ---")
    cols_to_show = ["molecule_id", "pIC50", "max_train_tanimoto", "in_applicability_domain"]
    print(test_df[cols_to_show].to_string(index=False))