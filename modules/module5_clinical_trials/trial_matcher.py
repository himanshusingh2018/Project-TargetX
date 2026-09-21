"""
Module 5: Clinical Trial Optimization
Script: trial_matcher.py
Evaluates candidate patients against clinical trial Inclusion/Exclusion (I/E) criteria,
integrating clinical chemistry, ECOG performance, and companion diagnostic biomarker status.
"""

import os
import numpy as np
import pandas as pd


def generate_synthetic_ehr(n_patients: int = 50, seed: int = 42) -> pd.DataFrame:
    """Generates synthetic Electronic Health Record (EHR) screening data for trial eligibility."""
    np.random.seed(seed)
    
    patient_ids = [f"PT-{i:03d}" for i in range(1, n_patients + 1)]
    ages = np.random.randint(28, 82, size=n_patients)
    # ECOG Performance Status (0 = fully active, 1 = restricted strenuous, 2 = ambulatory/self-care only, 3+ = bedbound)
    ecog_scores = np.random.choice([0, 1, 2, 3], size=n_patients, p=[0.40, 0.35, 0.18, 0.07])
    # Serum Creatinine (normal ~0.6-1.2 mg/dL; renal impairment >1.5)
    serum_creatinine = np.round(np.random.normal(1.1, 0.35, size=n_patients), 2)
    serum_creatinine = np.clip(serum_creatinine, 0.5, 3.5)
    # AST/ALT liver enzyme ratio (multiples of Upper Limit of Normal - ULN)
    liver_alt_uln = np.round(np.random.exponential(1.2, size=n_patients), 2)
    # CDx Biomarker status (from Module 4 model)
    cdx_positive = np.random.choice([True, False], size=n_patients, p=[0.35, 0.65])
    # Prior lines of targeted therapy
    prior_lines = np.random.choice([0, 1, 2, 3], size=n_patients, p=[0.2, 0.5, 0.2, 0.1])

    df = pd.DataFrame({
        "Patient_ID": patient_ids,
        "Age": ages,
        "ECOG_Status": ecog_scores,
        "Serum_Creatinine_mg_dL": serum_creatinine,
        "ALT_xULN": liver_alt_uln,
        "Prior_Therapy_Lines": prior_lines,
        "CDx_Biomarker_Positive": cdx_positive
    })
    return df


def screen_trial_eligibility(ehr_df: pd.DataFrame) -> pd.DataFrame:
    """Applies FDA/GCP compliant Inclusion and Exclusion Protocol Filters.
    Protocol Target-X Criteria:
    - Inclusion 1: Age >= 18
    - Inclusion 2: ECOG Status <= 1 (Patient must be functional)
    - Inclusion 3: CDx Biomarker Positive == True
    - Exclusion 1: Serum Creatinine > 1.5 mg/dL (Renal failure risk)
    - Exclusion 2: Liver ALT > 3.0 x ULN (Hepatic toxicity risk)
    - Exclusion 3: Prior targeted therapy lines > 2 (Heavily pre-treated / multi-resistant)
    """
    screened_records = []

    for _, row in ehr_df.iterrows():
        disqualification_reasons = []

        if row["Age"] < 18:
            disqualification_reasons.append("Pediatric (<18)")
        if row["ECOG_Status"] > 1:
            disqualification_reasons.append(f"Poor ECOG ({row['ECOG_Status']} > 1)")
        if not row["CDx_Biomarker_Positive"]:
            disqualification_reasons.append("CDx Biomarker Negative")
        if row["Serum_Creatinine_mg_dL"] > 1.5:
            disqualification_reasons.append(f"Renal Impairment (Creatinine {row['Serum_Creatinine_mg_dL']} > 1.5)")
        if row["ALT_xULN"] > 3.0:
            disqualification_reasons.append(f"Hepatic Toxicity Risk (ALT {row['ALT_xULN']} > 3.0x ULN)")
        if row["Prior_Therapy_Lines"] > 2:
            disqualification_reasons.append("Heavily Pre-treated (>2 prior lines)")

        is_eligible = len(disqualification_reasons) == 0
        screened_records.append({
            "Patient_ID": row["Patient_ID"],
            "Eligible": is_eligible,
            "Reason_Count": len(disqualification_reasons),
            "Disqualification_Audit": "; ".join(disqualification_reasons) if not is_eligible else "Meets All Protocol Criteria"
        })

    results_df = pd.DataFrame(screened_records)
    return pd.merge(ehr_df, results_df, on="Patient_ID")


if __name__ == "__main__":
    os.makedirs("data/raw", exist_ok=True)
    os.makedirs("data/processed", exist_ok=True)

    raw_ehr = generate_synthetic_ehr(n_patients=60)
    raw_ehr.to_csv("data/raw/synthetic_clinical_ehr.csv", index=False)

    screened_cohort = screen_trial_eligibility(raw_ehr)
    out_path = "data/processed/screened_trial_cohort.csv"
    screened_cohort.to_csv(out_path, index=False)

    n_enrolled = screened_cohort["Eligible"].sum()
    total_screened = len(screened_cohort)
    enrollment_rate = (n_enrolled / total_screened) * 100

    print("--- Clinical Trial Protocol Screening Summary ---")
    print(f"Total Patients Screened : {total_screened}")
    print(f"Patients Eligible       : {n_enrolled} ({enrollment_rate:.1f}% screen success rate)")
    print(f"Patients Excluded       : {total_screened - n_enrolled}")
    
    print("\n--- Sample Screened Audit Log (First 6 Patients) ---")
    sample_cols = ["Patient_ID", "ECOG_Status", "CDx_Biomarker_Positive", "Serum_Creatinine_mg_dL", "Eligible", "Disqualification_Audit"]
    print(screened_cohort[sample_cols].head(6).to_string(index=False))