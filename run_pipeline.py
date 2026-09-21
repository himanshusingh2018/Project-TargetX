"""
Project Target-X: End-to-End Orchestrator
Executes the full pipeline across all 6 modules sequentially.
"""

import sys
import subprocess
import time
import os

# Set PYTHONPATH to the current working directory so inter-module imports work
env = os.environ.copy()
env["PYTHONPATH"] = os.getcwd()

STEPS = [
    ("Module 1: Ingestion & Standardization", "modules/module1_data_ingestion/ingest_target_data.py"),
    ("Module 2: Molecular Featurization", "modules/module2_early_drug_discovery/featurizer.py"),
    ("Module 2: QSAR Binding Prediction", "modules/module2_early_drug_discovery/train_qsar.py"),
    ("Module 2: Applicability Domain Evaluation", "modules/module2_early_drug_discovery/applicability_domain.py"),
    ("Module 3: Knowledge Graph Centrality", "modules/module3_disease_modeling/disease_knowledge_graph.py"),
    ("Module 3: Bypass Resistance Analysis", "modules/module3_disease_modeling/pathway_bypass_analyzer.py"),
    ("Module 4: Patient Cohort Stratification", "modules/module4_biomarker_stratification/patient_clustering.py"),
    ("Module 4: Companion Diagnostic Classifier", "modules/module4_biomarker_stratification/train_classifier.py"),
    ("Module 5: Clinical Protocol Matching", "modules/module5_clinical_trials/trial_matcher.py"),
    ("Module 5: Dropout Risk Prediction", "modules/module5_clinical_trials/dropout_predictor.py"),
    ("Module 6: Enterprise ROI Valuation", "modules/module6_deployment/roi_calculator.py")
]


def execute_pipeline():
    print("\n" + "="*65)
    print("       LAUNCHING FULL PIPELINE: PROJECT TARGET-X")
    print("="*65 + "\n")
    
    start_all = time.time()

    for name, script_path in STEPS:
        if not os.path.exists(script_path):
            print(f"\n[ERROR] Missing script file: {script_path}")
            sys.exit(1)

        print(f"\n[RUNNING] >>> {name} ({script_path})")
        result = subprocess.run([sys.executable, script_path], capture_output=False, env=env)
        if result.returncode != 0:
            print(f"\n[ERROR] Pipeline aborted at step: {name}")
            sys.exit(result.returncode)

    duration = time.time() - start_all
    print("\n" + "="*65)
    print(f" [+] PIPELINE COMPLETED SUCCESSFULLY IN {duration:.2f} SECONDS")
    print("="*65 + "\n")


if __name__ == "__main__":
    execute_pipeline()