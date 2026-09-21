"""
Module 2: Preprocessing & Validation Splits
Implements Bemis-Murcko Scaffold Splitting to prevent chemical data leakage.
"""

from collections import defaultdict
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold


def get_murcko_scaffold(smiles: str) -> str:
    """Generates the Bemis-Murcko scaffold SMILES for a given molecule.
    Returns empty string for acyclic molecules.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    # Strip side chains to isolate the core ring system and linkers
    scaffold = MurckoScaffold.MakeScaffoldGeneric(mol) # Generic backbone (all atoms -> C, all bonds -> single)
    scaffold_mol = MurckoScaffold.GetScaffoldForMol(mol)
    return Chem.MolToSmiles(scaffold_mol)


def scaffold_split(df: pd.DataFrame, smiles_col: str = "smiles", train_ratio: float = 0.8, seed: int = 42):
    """Splits a DataFrame such that no scaffold in the train set appears in the test set.
    """
    scaffold_dict = defaultdict(list)

    for idx, row in df.iterrows():
        scaffold = get_murcko_scaffold(row[smiles_col])
        scaffold_dict[scaffold].append(idx)

    # Sort scaffold clusters by size (largest clusters assigned first)
    scaffold_clusters = list(scaffold_dict.values())
    scaffold_clusters.sort(key=len, reverse=True)

    train_cutoff = int(len(df) * train_ratio)
    train_indices = []
    test_indices = []

    for cluster in scaffold_clusters:
        if len(train_indices) + len(cluster) <= train_cutoff:
            train_indices.extend(cluster)
        else:
            test_indices.extend(cluster)

    train_df = df.iloc[train_indices].copy().reset_index(drop=True)
    test_df = df.iloc[test_indices].copy().reset_index(drop=True)

    return train_df, test_df


if __name__ == "__main__":
    data_path = "data/processed/clean_bioactivity.csv"
    df = pd.read_csv(data_path)
    print(f"[+] Loaded {len(df)} compounds for scaffold analysis.")

    # Inspect the scaffolds present in our clean bioactivity set
    df["scaffold"] = df["smiles"].apply(get_murcko_scaffold)
    print("\n--- Identified Molecular Scaffolds ---")
    for idx, row in df.iterrows():
        scaffold_repr = row['scaffold'] if row['scaffold'] else "[Acyclic / No Ring]"
        print(f"Molecule: {row['molecule_id']:<8} | Scaffold: {scaffold_repr}")

    train_set, test_set = scaffold_split(df, smiles_col="smiles", train_ratio=0.7)
    print(f"\n[+] Scaffold Split Summary:")
    print(f"    - Training Set: {len(train_set)} compounds ({list(train_set['molecule_id'])})")
    print(f"    - Test Set    : {len(test_set)} compounds ({list(test_set['molecule_id'])})")