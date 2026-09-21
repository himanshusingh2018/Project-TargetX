"""
Module 1: Batch Ingestion & Standardization
Processes bioactivity records, handles IC50 -> pIC50 conversion,
and enforces molecular quality control filters.
"""

import os
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors


def calculate_pchembl(ic50_nm: float) -> float:
    """Converts IC50 in nanomolar (nM) to pIC50 (-log10(Molar)).
    Caps values to prevent infinite or negative values on outlier noise.
    """
    if pd.isna(ic50_nm) or ic50_nm <= 0:
        return np.nan
    # Convert nM to Molar (1 nM = 1e-9 M)
    molar = ic50_nm * 1e-9
    pIC50 = -np.log10(molar)
    return round(float(pIC50), 3)


def standardize_compound(smiles: str) -> dict:
    """Validates SMILES and computes standard physicochemical descriptors.
    Returns None if the SMILES cannot be parsed by RDKit.
    """
    if not isinstance(smiles, str) or not smiles.strip():
        return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mw = round(Descriptors.MolWt(mol), 2)
    logp = round(Descriptors.MolLogP(mol), 2)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    rotatable_bonds = Descriptors.NumRotatableBonds(mol)

    violations = sum([mw > 500, logp > 5.0, hbd > 5, hba > 10])

    return {
        "canonical_smiles": Chem.MolToSmiles(mol, canonical=True),
        "molecular_weight": mw,
        "logp": logp,
        "hbd": hbd,
        "hba": hba,
        "rotatable_bonds": rotatable_bonds,
        "ro5_violations": violations,
        "is_druglike": violations <= 1
    }


def run_pipeline():
    os.makedirs("data/raw", exist_ok=True)
    os.makedirs("data/processed", exist_ok=True)

    # Benchmark dataset: Representative EGFR small-molecule bioactivity assay
    benchmark_data = [
        {"molecule_id": "MOL001", "smiles": "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1", "ic50_nm": 0.5},   # Gefitinib
        {"molecule_id": "MOL002", "smiles": "C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1", "ic50_nm": 2.0},       # Erlotinib
        {"molecule_id": "MOL003", "smiles": "CC(=O)Nc1ccc(O)cc1", "ic50_nm": 85000.0},                             # Paracetamol (Inactive on EGFR)
        {"molecule_id": "MOL004", "smiles": "INVALID_SMILES_STRING", "ic50_nm": 12.0},                              # Corrupted Entry
        {"molecule_id": "MOL005", "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "ic50_nm": 45000.0},                 # Caffeine (Inactive on EGFR)
        {"molecule_id": "MOL006", "smiles": "CS(=O)(=O)CCNCc1ccc(o1)-c1ccc2c(c1)ncn2-c1ccc(OCc2cccc(F)c2)c(Cl)c1", "ic50_nm": 9.8}, # Lapatinib analog
        {"molecule_id": "MOL007", "smiles": "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O", "ic50_nm": 32000.0} # Atorvastatin
    ]

    raw_df = pd.DataFrame(benchmark_data)
    raw_path = "data/raw/raw_bioactivity_egfr.csv"
    raw_df.to_csv(raw_path, index=False)
    print(f"[+] Step 1: Raw data ingested with {len(raw_df)} records -> {raw_path}")

    # Process and clean
    records = []
    dropped_count = 0

    for _, row in raw_df.iterrows():
        props = standardize_compound(row["smiles"])
        if props is None:
            dropped_count += 1
            continue

        pIC50 = calculate_pchembl(row["ic50_nm"])
        # Label: active if pIC50 >= 7.0 (i.e., IC50 <= 100 nM)
        active_flag = 1 if pIC50 >= 7.0 else 0

        combined = {
            "molecule_id": row["molecule_id"],
            "smiles": props["canonical_smiles"],
            "molecular_weight": props["molecular_weight"],
            "logp": props["logp"],
            "hbd": props["hbd"],
            "hba": props["hba"],
            "rotatable_bonds": props["rotatable_bonds"],
            "ro5_violations": props["ro5_violations"],
            "is_druglike": props["is_druglike"],
            "ic50_nm": row["ic50_nm"],
            "pIC50": pIC50,
            "bioactivity_class": active_flag
        }
        records.append(combined)

    clean_df = pd.DataFrame(records)
    processed_path = "data/processed/clean_bioactivity.csv"
    clean_df.to_csv(processed_path, index=False)
    print(f"[+] Step 2: Quality control complete. Dropped {dropped_count} invalid records.")
    print(f"[+] Step 3: Saved clean bioactivity dataset with {len(clean_df)} valid compounds -> {processed_path}\n")
    print(clean_df[["molecule_id", "molecular_weight", "pIC50", "bioactivity_class", "is_druglike"]])


if __name__ == "__main__":
    run_pipeline()