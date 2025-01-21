# Reproducible Example (Fig. 6)

**Felix Boel**



## Overview

This subdirectory contains a reproducible example for Figure 6 from our study titled “Deep Proteome Profiling of Metabolic Dysfunction-Associated Steatotic Liver Disease” by Boel et al. The example demonstrates how to load clinical and proteomics data, perform classification and ROC analyses using our helper functions, and generate the corresponding plots for Figure 6.

## Repository Structure


   ```markdown
   Figure 6 example/
   ├── data/
   │   └── # (Data files will be added here upon journal acceptance)
   ├── plots/
   │   └── # Generated plots will be saved here
   ├── run_analysis.R
   └── README.md
   ```

- **`run_analysis.R`**
  - **Description:** R script that performs the entire analysis workflow to reproduce Figure 6. It includes data loading, preprocessing, model training, cross-validation, metric calculations, and plotting.

- **`data/`**
  - **Description:** Directory to store the required clinical and proteomics data files. (Note: Data files are temporarily excluded and will be added upon full acceptance by the journal.)

- **`plots/`**
  - **Description:** Directory where the generated plots for Figure 6 will be saved.

## Required R Packages

The following R packages are required:

```r
install.packages(c("caret", "pROC", "ggplot2", "dplyr", "tidyr", "forcats", "reshape2", "tibble"))
```
These packages will be installed automatically when you run the script if they are not already installed.

## The `run_analysis.R` script performs the following steps:

1. **Load Dependencies:** Installs and loads the required R packages.
2. **Source Utility Functions:** Loads helper functions from `roc_curve_utils.R`, `cross_validation_utils.R`, and `metric_calculations.R`.
3. **Load Plasma Data:** Reads and preprocesses the proteomics dataset.
4. **Load Clinical Data:** Reads and preprocesses the clinical dataset.
5. **Define Targets:** Creates binary classification targets for significant fibrosis (≥F2) and MASLD.
6. **Merge Clinical & Plasma Data:** Combines clinical data with proteomics data for both baseline and follow-up samples.
7. **Set Cross-Validation Parameters:** Configures k-fold cross-validation with specified repetitions.
8. **Perform Cross-Validated ROC Calculations:** Evaluates different models (APRI, FIB4, Proteomics) using cross-validation.
9. **Calculate Simple Metrics:** Computes basic classification metrics like Accuracy and F1-Score.
10. **Combine Performance Metrics:** Aggregates metrics from simple models and cross-validated logistic regression models.
11. **Plot Performance Metrics:** Generates performance plots saved in the `plots/` directory.
12. **Additional ROC Analysis:** Recomputes preprocesses AUC-based metrics.
13. **Plot ROC Curves:** Generates ROC curve plots saved in the `plots/` directory.

## Citation

If you use this reproducible example in your research, please cite our paper:
**Boel, et al. Deep Proteome Profiling of Metabolic Dysfunction-Associated Steatotic Liver Disease**
