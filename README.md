# Translational Immunology ML
Computational immunology pipeline for analyzing high-dimensional multiplex cytokine data across donors, biological matrices, and freeze-thaw conditions.

---

## Project Overview
This project analyzes multiplex Luminex cytokine assay data collected from healthy donors across EDTA plasma, heparin plasma, and serum matrices over repeated freeze-thaw cycles.

The goal was to identify:
- cytokines that remain analytically stable
- matrix-specific immune signatures
- biomarkers sensitive to sample handling conditions

The workflow combines computational biology, statistical modeling, and machine learning to study high-dimensional immune profiling data. This work contributed to a first-author publication in the Journal of Immunological Methods (Verma et al., 2026).

---

## Computational Methods

### Data Processing & QC
- Raw MFI extraction
- Limit of detection (LOD) analysis
- Detectability filtering
- Coefficient of variation (CV) analysis

### Statistical & ML Analysis
- Principal Component Analysis (PCA)
- Mixed-effects modeling (`nlme`)
- Matrix-specific slope analysis
- Random Forest classification (R and Python)
- XGBoost classification (Python)
- Feature importance analysis
- Heatmap visualization

---

## Main Analysis Scripts

| Script | Purpose |
|---|---|
| `freeze_thaw_pipeline.R` | Main statistical and computational analysis workflow |
| `panel_qc_extraction.R` | Raw MFI extraction and QC pipeline |
| `supplement_qc_table.R` | Detectability and CV summary generation |
| `matrix_prediction_model.R` | Random Forest matrix classification (R) |

---

## Python ML Notebook

A Python implementation of the matrix classification pipeline is available in `notebooks/cytokine_matrix_classifier_v2.ipynb`.

This notebook extends the original R analysis using scikit-learn and XGBoost:
- Random Forest classifier: 100% accuracy on held-out test set
- XGBoost classifier: 75% accuracy (consistent with RF's known advantage on small, high-dimensional datasets)
- Feature importance analysis identifying GCSF/CSF3 and EOTAXIN/CCL11 as top matrix-discriminating cytokines

---

## Repository Structure

```text
src/        Analysis and preprocessing scripts (R)
data/       Example processed cytokine datasets
results/    Output tables and statistical summaries
figures/    Publication-ready visualizations
notebooks/  Jupyter notebook (Python ML pipeline)
```

---

## Example Visualizations
- `PCA_MultiPanel_Analysis.png`
- `Cytokines_Heatmap_MatrixCycle.png`
- `IL6_UNIFORM.png`
- `matrix_prediction_feature_importance.png` (R)
- `feature_importance.png` (Python)

---

## Example Output Tables
- `stable_cytokines.csv`
- `decreasing_cytokines.csv`
- `matrix_effect_cytokines.csv`
- `suppl_table1.csv`

---

## Tools & Libraries

**R:** tidyverse, ggplot2, nlme, emmeans, randomForest, caret, cowplot, viridis

**Python:** scikit-learn, XGBoost, pandas, numpy, matplotlib

---

## Publication

Verma A, Sridharan K, Herschmann I, Nguyen T, Maecker HT. Stability of Luminex cytokine assays with freeze-thaw in different plasma/serum matrices. Journal of Immunological Methods. 2026. https://www.sciencedirect.com/science/article/pii/S0022175926000293
