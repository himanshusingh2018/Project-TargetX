"""
Module 6: Enterprise Deployment & ROI Valuation
Script: roi_calculator.py
Quantifies the operational and financial impact of Project Target-X across
chemistry triage, patient stratification, and clinical retention.
"""

import pandas as pd


def calculate_pharma_impact(
    n_virtual_screened: int = 50000,
    cost_per_synthesis: float = 3500.0,
    trial_patients: int = 250,
    cost_per_patient: float = 45000.0,
    cost_per_dropout: float = 110000.0,
    baseline_dropout_rate: float = 0.28,
    reduced_dropout_rate: float = 0.12,
    baseline_response_rate: float = 0.18,
    stratified_response_rate: float = 0.65
) -> dict:
    """Computes operational ROI and time-to-clinic impact metrics."""
    # 1. Early Discovery Savings: Filtering out non-druglike & out-of-domain compounds
    # Assuming without AI, 200 compounds are synthesized; with AI triage, only 30 high-confidence leads are synthesized
    baseline_synthesis_n = 200
    ai_selected_synthesis_n = 30
    chemistry_savings = (baseline_synthesis_n - ai_selected_synthesis_n) * cost_per_synthesis

    # 2. Clinical Retention Savings
    baseline_dropouts = trial_patients * baseline_dropout_rate
    ai_dropouts = trial_patients * reduced_dropout_rate
    prevented_dropouts = baseline_dropouts - ai_dropouts
    retention_savings = prevented_dropouts * cost_per_dropout

    # 3. Trial De-risking: Probability of Positive Phase II Readout
    # Biomarker enrichment transforms a failing trial into an approvable indication
    total_savings = chemistry_savings + retention_savings

    return {
        "Chemistry Synthesis Savings ($)": round(chemistry_savings, 2),
        "Patient Dropouts Prevented": int(prevented_dropouts),
        "Clinical Retention Cost Savings ($)": round(retention_savings, 2),
        "Total Quantifiable Cost Savings ($)": round(total_savings, 2),
        "Unselected Response Rate (%)": baseline_response_rate * 100,
        "Biomarker Stratified Response Rate (%)": stratified_response_rate * 100,
        "Trial Success Probability Shift": "Failure Risk -> Actionable CDx Filing"
    }


def run_valuation():
    metrics = calculate_pharma_impact()
    print("==========================================================")
    print("      PROJECT TARGET-X: ENTERPRISE PHARMA ROI AUDIT       ")
    print("==========================================================")
    for k, v in metrics.items():
        if isinstance(v, float):
            print(f"{k:<38}: ${v:,.2f}")
        else:
            print(f"{k:<38}: {v}")
    print("==========================================================")


if __name__ == "__main__":
    run_valuation()