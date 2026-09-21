"""
Module 1: Data Ingestion & Chemoinformatics Basics
Calculates physical descriptors and evaluates Lipinski's Rule of 5.
"""

import sys
from rdkit import Chem
from rdkit.Chem import Descriptors


def evaluate_lipinski(smiles: str) -> dict:
    """Calculates molecular properties and evaluates Lipinski's Rule of 5.
    A compound fails drug-likeness if it violates > 1 rule.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string provided: {smiles}")

    mw = round(Descriptors.MolWt(mol), 2)
    logp = round(Descriptors.MolLogP(mol), 2)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)

    # Lipinski violation checks
    violations = 0
    if mw > 500:
        violations += 1
    if logp > 5.0:
        violations += 1
    if hbd > 5:
        violations += 1
    if hba > 10:
        violations += 1

    passes_rule_of_5 = violations <= 1

    return {
        "SMILES": smiles,
        "MolecularWeight": mw,
        "LogP": logp,
        "HBD": hbd,
        "HBA": hba,
        "Violations": violations,
        "Lipinski_Pass": passes_rule_of_5
    }


if __name__ == "__main__":
    # If a SMILES is passed from the terminal, use it; otherwise default to Aspirin
    target_smiles = sys.argv[1] if len(sys.argv) > 1 else "CC(=O)Oc1ccccc1C(=O)O"
    
    results = evaluate_lipinski(target_smiles)
    print("--- Molecular Profiling & Lipinski Rule of 5 ---")
    for k, v in results.items():
        print(f"{k:>18}: {v}")

    # Analysis meterics:
    print(
        """
    LogP > 5: The drug acts like grease or oil. It precipitates out of solution in the watery bloodstream, gets trapped in body fat, or gets stuck inside cell membranes without exiting.
    LogP < 0: The drug loves water too much. It cannot cross the lipid bilayer of human cell membranes or the gut wall to enter circulation.
    1 < LogP < 4: Soluble enough to travel through blood, yet lipophilic enough to permeate cell membranes and reach intracellular targets.
    
    Molecular Weight <= 500 Da. More than 500 Da, and the drug is too large to cross cell membranes and reach intracellular targets.
    HBD <= 5 (sum of OH and HG groups): Too many hydrogen bonds hold onto water molecules, making it hard to slip through lipid membranes.
    HBA <= 10 (sum of N and O atoms): Same reason; too much hydrogen-bonding locks the molecule in aqueous solvent.
    """
    )