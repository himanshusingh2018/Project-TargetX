# Project Target-X: An End-to-End Computational R&D Pipeline for Oncology Drug Discovery & Clinical Translation

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![Code Style: Black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

## Executive Summary
Developing a novel pharmaceutical therapeutic requires **10–15 years** and upwards of **$2.6 billion**, with over **90% of oncology clinical assets failing** in human trials. The root cause is rarely algorithmic sophistication alone; it is the structural fragmentation between early chemistry discovery, target disease biology, and clinical trial design.

**Project Target-X** is an end-to-end, modular computational pipeline demonstrating how AI and data science systematically de-risk drug development across the entire pharmaceutical value chain:
1. **Target Identification:** Biological Knowledge Graph analysis prioritizing signaling choke points over non-selective hubs to minimize off-target toxicity.
2. **Early Discovery & QSAR:** 1024-bit Morgan Fingerprint (ECFP4) featurization with Bemis-Murcko scaffold splitting (preventing chemical data leakage) and Tanimoto-based Applicability Domain gating.
3. **Precision Oncology Stratification:** Unsupervised transcriptomic subtyping (PCA/K-Means) paired with an L1-regularized Companion Diagnostic (CDx) gene signature classifier.
4. **Clinical Trial Optimization:** Automated EHR protocol Inclusion/Exclusion screening and proactive machine learning risk modeling for early patient dropout.
5. **Enterprise ROI Valuation:** Direct financial and operational quantification yielding **~$5.0M in capital efficiency** per program.

---

## Pipeline Architecture

```text
Project-TargetX/
├── data/
│   ├── raw/                             # Raw bioactivity, patient RNA-seq cohorts, EHR logs
│   └── processed/                       # Cleaned assays, graph metrics, stratified cohorts
├── modules/
│   ├── module1_data_ingestion/          # SMILES parsing, Lipinski Rule of 5, IC50 -> pIC50
│   ├── module2_early_drug_discovery/    # ECFP4 generation, Murcko scaffold split, QSAR, AD
│   ├── module3_disease_modeling/        # NetworkX Knowledge Graph, centrality, bypass routes
│   ├── module4_biomarker_stratification/# Cohort clustering, sparse CDx signature classifier
│   ├── module5_clinical_trials/         # Protocol I/E rule engine, dropout risk prediction
│   └── module6_deployment/              # Business valuation and financial ROI calculator
├── run_pipeline.py                      # Master sequential orchestrator
└── README.md


**Core Scientific Modules & Methodologies**
1. Ingestion & Quality Control (module1_data_ingestion)SMILES Syntax Validation & Parsing: 
Uses RDKit to sanitize SMILES and filter invalid structures.Physicochemical Filtering: Computes Lipinski’s Rule of 5 descriptors (Molecular Weight, LogP, HBD, HBA) to assess oral bioavailability heuristics.Affinity Standardization: Converts heterogeneous $IC_{50}$ values to continuous logarithmic binding affinity (pIC50 = -log10(IC50/10^-9)).

2. Early Discovery & QSAR Modeling (module2_early_drug_discovery)Structural Featurization: 
Encodes molecular topologies into 1024-bit circular Morgan Fingerprints (ECFP4, radius 2).Bemis-Murcko Scaffold Splitting: Groups compounds by core ring systems and side chains rather than random splits, preventing chemical data leakage and measuring true out-of-scaffold generalization.Applicability Domain Guardrail: Evaluates maximum Tanimoto similarity against the training manifold ($T_c \ge 0.40$), preventing out-of-distribution hallucinations on chemical series the model has never encountered.

3. Systems Biology & Knowledge Graphs (module3_disease_modeling)Topological Centrality: 
Constructs an oncogenic signaling network (EGFR/MAPK/PI3K/Survival axes) in NetworkX. Balances Degree Centrality against Betweenness Centrality to detect low-degree, high-bottleneck targets (e.g., AKT1, MTOR) that shut down oncogenic flow while minimizing collateral toxicity.In Silico Resistance Modeling: Simulates primary target knockout (e.g., EGFR blockade) to identify alternative bypass routes and propose rational combination therapies.

4. Patient Stratification & Companion Diagnostics (module4_biomarker_stratification)Transcriptomic Subtyping: 
Performs standard scaling, PCA, and K-Means clustering on patient RNA-seq profiles to identify receptor-dependent responders, downstream MAPK drivers, and PTEN-loss/AKT-bypass resistant populations.Sparse Feature Selection: Trains an L1-regularized Logistic Regression classifier to isolate a minimal, interpretable 8-gene Companion Diagnostic (CDx) signature suitable for affordable RT-qPCR clinical translation.

5. Clinical Trial Operations (module5_clinical_trials)Protocol Matching Engine: 
Ingests synthetic electronic health records (EHR) and deterministically screens candidates against GCP-compliant inclusion/exclusion criteria (ECOG performance status, renal/hepatic safety markers, CDx biomarker status).Dropout Prediction: Trains a Gradient Boosting classifier on operational friction metrics (distance to trial site, age, adverse events) to flag patients at risk of premature study withdrawal for decentralized clinical retention support.6. Enterprise ROI & Valuation (module6_deployment)Translates technical algorithms into biopharma business impact:Chemistry Triage: Saves ~$595,000 by eliminating non-druglike and out-of-domain compounds before bench synthesis.Clinical Retention: Prevents ~40 trial dropouts, unlocking ~$4.4M in operational retention savings.Trial De-risking: Shifts projected Phase II response rates from 18% (unselected cohort) to 65% (biomarker-selected cohort).


**Quickstart & InstallationEnvironment SetupBash# Clone repository**
git clone [https://github.com/himanshu2018/Project-TargetX.git](https://github.com/himanshu2018/Project-TargetX.git)
cd Project-TargetX

# Install dependencies using pip or uv
pip install -r pyproject.toml
# or
uv pip install rdkit pandas numpy scikit-learn networkx scipy
Running the End-to-End PipelineExecute the master orchestrator from the project root:Bashpython run_pipeline.py
Tech StackChemoinformatics: RDKitMachine Learning & Stats: scikit-learn, NumPy, Pandas, SciPyNetwork Biology: NetworkXEnvironment & Packaging: Python 3.10+, uv
---

### How to Save and Commit

Save this directly to your repository:

```bash
# Stage the updated README
git add README.md

# Commit
git commit -m "docs: add comprehensive executive README with architecture and methodology"

# Push to GitHub
git push origin main