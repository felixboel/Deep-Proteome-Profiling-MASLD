# Deep Proteome Profiling of Metabolic Dysfunction-Associated Steatotic Liver Disease

**Boel et al.**



## Overview

This repository contains helper functions used for data analysis in the study titled **"Deep Proteome Profiling of Metabolic Dysfunction-Associated Steatotic Liver Disease"** by Boel et al. The functions facilitate the evaluation of proteomics biomarkers through logistic regression models and cross-validation techniques, specifically utilized for generating Fig. 6 in the publication.

## Repository Contents

- **`roc_curve_utils.R`**
  - **Function:** `roc.curve`
  - **Description:** Processes ROC curves to compute mean specificities, sensitivities, AUC values, and their confidence intervals.

- **`cross_validation_utils.R`**
  - **Function:** `cross_fold_calcs_roc`
  - **Description:** Performs k-fold cross-validation with multiple repetitions, trains logistic regression models, evaluates performance metrics, and aggregates ROC curves.

- **`metric_calculations.R`**
  - **Function:** `simple_metric_calcs`
  - **Description:** Calculates basic classification metrics (Accuracy, Balanced Accuracy, Precision, Recall, F1-Score) from confusion matrices.

## Fig. 6: Performance evaluation of proteomics biomarkers for MASLD and ≥F2

Figure 6 shows the performance metrics of our proteomics-based logistic regression models after 20x five-fold cross-validation. Additionally, we compare to existing models for predicting MASLD and ≥F2.

### **Fig. 6 Details:**

- **Models:**
  - **Standard non-proteomic biomarkers:** APRI, FIB-4, NAFLD liver fat score, Hepatic steatosis index, Fatty liver index
  - **Proteomics Models:**
    - **Significant Fibrosis (≥F2):** Based on proteins from Fig. 5a (C7, COLEC11, PRAP1, CTSD, OSMR, ALDOB, PECR, A2M, ICAM1,AFM, ADIPOQ, SERPINC1, IGFBP5, APOM, C8A, IGFBP3, APOF, IGFALS, IGF2, IGF1).
    - **MASLD:** Based on proteins from Fig. 5b (ALDOB, AFM, CTSD, OAF, ICAM1, F9, SERPINF1, PRAP1, GNPTG, PRG4, ALB, TBC1D14, CETP, PON3, GSN, GPX3, ADIPOQ, SERPINC1, APOF, IGF1).
  - **Simplified Proteomics Models:**
    - **≥F2:** C7, ALDOB, ICAM1, CTSD, IGF2
    - **MASLD:** AFM, APOF, ALDOB, PRG4, IGF1

- **Metrics Covered:**
  - Receiver-Operator Characteristic (ROC) Curves
  - Balanced Accuracy
  - F1 Score

- **Cohorts:**
  - **Initial Cohort (n = 143):** Metrics a-c,e
  - **Validation Cohort (n = 41):** Metrics d,f

## Installation

Ensure you have the necessary R packages installed.

```r
install.packages(c("caret", "pROC"))
```

## Example Workflow

```r
# Load necessary libraries
library(caret)
library(pROC)

# Perform cross-validation and collect ROC curves
results <- cross_fold_calcs_roc(dataset = your_dataset, dataset2 = your_extra_dataset, k_folds = 5, reps = 20)

# Compute average ROC metrics
average_roc <- roc.curve(results$roc_curves)
average_roc_extra <- roc.curve(results$roc_curves2)

# Calculate simple metrics (APRI, NAFLD liver fat score, etc)
metrics <- simple_metric_calcs(your_confusion_matrix_dataframe)

# View Results
print(results$output_combined)
print(average_roc)
print(metrics)
```


## Citation

If you use these helper functions in your research, please cite our paper:
**Boel, et al. Deep Proteome Profiling of Metabolic Dysfunction-Associated Steatotic Liver Disease**
