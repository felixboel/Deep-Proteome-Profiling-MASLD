# Reproducible Example for Fig. 6

**Felix Boel**



## Overview

This subdirectory contains the reproducible example for Figure 6 from our study titled “Deep Proteome Profiling of Metabolic Dysfunction-Associated Steatotic Liver Disease” by Boel et al. The example demonstrates how to load clinical and proteomics data, perform classification and ROC analyses, and generate the corresponding plots for Figure 6.

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

## Example Workflow

The run_analysis.R script performs the following steps:
	1.	**Load Dependencies:** Installs and loads the required R packages.
	2.	**Source Utility Functions:** Loads helper functions from roc_curve_utils.R, cross_validation_utils.R, and metric_calculations.R.
	3.	**Load Plasma and Clinical Data:** Reads and preprocesses the proteomics and clinical datasets.
	4.	**Define Targets:** Creates binary classification targets for F2 (Fibrosis) and MASLD.
	5.	**Merge Clinical & Plasma Data:** Combines clinical data with proteomics data for both baseline and follow-up samples.
	6.	**Set Cross-Validation Parameters:** Configures k-fold cross-validation with specified repetitions.
	7.	**Perform Cross-Validated ROC Calculations:** Evaluates different models (APRI, FIB4, Proteomics) using cross-validation.
	8.	**Calculate Simple Metrics:** Computes basic classification metrics like Accuracy and F1-Score.
	9.	**Combine and Plot Metrics:** Aggregates metrics and generates performance plots saved in the plots/ directory.


## Citation

If you use this reproducible example in your research, please cite our paper:
**Boel, et al. Deep Proteome Profiling of Metabolic Dysfunction-Associated Steatotic Liver Disease**
