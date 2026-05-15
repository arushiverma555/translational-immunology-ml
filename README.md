# Translational Immunology ML

Computational immunology pipeline for analyzing high-dimensional multiplex cytokine data across donors, biological matrices, and freeze-thaw conditions.

---

## Project Overview

This project analyzes multiplex Luminex cytokine assay data collected from healthy donors across EDTA plasma, heparin plasma, and serum matrices over repeated freeze-thaw cycles.

The goal was to identify:
- cytokines that remain analytically stable
- matrix-specific immune signatures
- biomarkers sensitive to sample handling conditions

The workflow combines computational biology, statistical modeling, and machine learning-oriented dimensionality reduction techniques to study high-dimensional immune profiling data.

---

## Computational Methods

### Data Processing & QC
- Raw MFI extraction
- Limit of detection (LOD) analysis
- Detectability filtering
- Coefficient of variation (CV) analysis

### Statistical & ML-Oriented Analysis
- Principal Component Analysis (PCA)
- Mixed-effects modeling (`nlme`)
- Matrix-specific slope analysis
- Cytokine classification pipelines
- Heatmap visualization
- High-dimensional immune profiling

---

## Main Analysis Scripts

| Script | Purpose |
|---|---|
| `freeze_thaw_pipeline.R` | Main statistical and computational analysis workflow |
| `panel_qc_extraction.R` | Raw MFI extraction and QC pipeline |
| `supplement_qc_table.R` | Detectability and CV summary generation |

---

## Repository Structure

```text
src/        Analysis and preprocessing scripts
data/       Example processed cytokine datasets
results/    Output tables and statistical summaries
figures/    Publication-ready visualizations
notebooks/  Exploratory analyses and workflow notes
```

---

## Example Visualizations

- `PCA_MultiPanel_Analysis.png`
- `Cytokines_Heatmap_MatrixCycle.png`
- `IL6_UNIFORM.png`

---

## Example Output Tables

- `stable_cytokines.csv`
- `decreasing_cytokines.csv`
- `matrix_effect_cytokines.csv`
- `suppl_table1.csv`

Representative outputs are available in the `figures/` and `results/` directories.

---

## Tools & Libraries

- R
- tidyverse
- ggplot2
- nlme
- emmeans
- cowplot
- viridis

---

## Future Directions

Planned extensions include:
- supervised machine learning classification of biological matrices from cytokine profiles
- feature importance analysis
- predictive modeling of cytokine stability across sample conditions
- clustering-based immune phenotype discovery

---

## Exploratory Machine Learning Extension

An exploratory random forest classifier was implemented to investigate whether cytokine slope behavior and matrix-interaction metrics could distinguish cytokine stability classes.

Feature importance analysis identified matrix interaction p-values and matrix-specific slope features as major contributors to classification performance.

Outputs include:
- random forest feature importance visualization
- cytokine stability classification framework
- matrix-specific predictive feature analysis
