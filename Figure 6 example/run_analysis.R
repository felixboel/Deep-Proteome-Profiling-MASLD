###############################################################################
# File: run_analysis.R
# Project: Deep Proteome Profiling of Metabolic Dysfunction-Associated Steatotic Liver Disease
# Author:  Felix Boel
# Purpose: Demonstrates how to load clinical & proteomics data, perform
#          classification and ROC analyses, and generate plots.
###############################################################################

# Clear environment and console
rm(list = ls())
cat("\014")

#--------------------------#
#   1. Load Dependencies   #
#--------------------------#

# Define required packages
required_packages <- c("caret", "pROC", "ggplot2", "dplyr", "tidyr", 
                       "forcats", "reshape2", "tibble")

# Install and load packages
suppressMessages(suppressWarnings({
  invisible(lapply(required_packages, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }))
}))

#--------------------------#
#   2. Source Utilities    #
#--------------------------#
# These custom scripts contain helper functions for cross-validation, 
# metric calculations, and ROC curve processing.

source("../roc_curve_utils.R")       # For ROC curve averaging & confidence intervals (requirement for cross_validation_utils.R)
source("../cross_validation_utils.R") # For cross-fold calculations & performance metrics
source("../metric_calculations.R")   # For simple metric calculations (accuracy, F1, etc.)


#--------------------------#
#   3. Load Plasma Data    #
#--------------------------#

# Plasma proteomics data (Supplementary Data 3)
plasma_proteomics_data <- read.delim2(
  "data/Supplementary Data 3_plasma.txt",
  sep = "\t",
  dec = ".",
  skip = 0,
  header = TRUE,
  na.strings = c("NaN", "NA")
)

# Set protein accessions as row names and remove that column
rownames(plasma_proteomics_data) <- plasma_proteomics_data$ProteinAccessions
plasma_proteomics_data <- plasma_proteomics_data[, colnames(plasma_proteomics_data) != "ProteinAccessions"]

# Define proteins of interest for different analyses
proteins_of_interest_lists <- list(
  f2 = c(
    "P10643", "Q9BWP8", "Q96NZ9", "P07339", "Q99650", "P05062", "Q9BY49", "P01023",
    "P05362", "P43652", "Q15848", "P01008", "P24593", "O95445", "P07357",
    "P17936", "Q13790", "P35858", "P01344", "P05019"
  ),
  f2_mini = c("P10643", "P05062", "P05362", "P07339", "P01344"),
  masld = c(
    "P05062", "P43652", "P07339", "Q86UD1", "P05362", "P00740", "P36955", "Q96NZ9",
    "Q9UJJ9", "Q92954", "P02768", "Q9P2M4", "P11597", "Q15166", "P06396", "P22352",
    "Q15848", "P01008", "Q13790", "P05019"
  ),
  masld_mini = c("P43652", "Q13790", "P05062", "Q92954", "P05019")
)

# Subset plasma data by proteins of interest
plasma_f2 <- plasma_proteomics_data[rownames(plasma_proteomics_data) %in% proteins_of_interest_lists$f2, ]
plasma_mini_f2 <- plasma_proteomics_data[rownames(plasma_proteomics_data) %in% proteins_of_interest_lists$f2_mini, ]
plasma_masld <- plasma_proteomics_data[rownames(plasma_proteomics_data) %in% proteins_of_interest_lists$masld, ]
plasma_mini_masld <- plasma_proteomics_data[rownames(plasma_proteomics_data) %in% proteins_of_interest_lists$masld_mini, ]


#--------------------------#
#   4. Load Clinical Data  #
#--------------------------#

clinical_data <- read.delim2(
  "data/Supplementary Data 1_Meta.txt",
  sep = "\t",
  dec = ".",
  skip = 0,
  header = TRUE,
  na.strings = c("NaN", "NA")
)

clinical_data <- clinical_data[,c(1,3,20:24,8,13)]
clinical_data <- clinical_data[clinical_data$unique_identifier %in% colnames(plasma_proteomics_data),]
rownames(clinical_data) <- NULL
clinical_data$sample_group <- ifelse(clinical_data$sample_group == "follow_up_sample", yes = "followup", no = "initial")


#--------------------------#
#   5. Define Targets      #
#--------------------------#

# Create binary classification targets
clinical_data$F2_target <- ifelse(clinical_data$kleiner_fibrosis_grade %in% c("F2", "F3", "F4"), yes = 1, no = 0)
clinical_data$MASLD_target <- ifelse(clinical_data$saf_diagnosis %in% c("MASL", "MASH"), yes = 1, no = 0)

# Predictions based on classical scores
clinical_data$NAFLD_liver_fat_score_prediction <- ifelse(clinical_data$NAFLD_liver_fat_score > -0.640, yes = 1, no = 0)
clinical_data$Hepatic_steatosis_index_prediction <- ifelse(clinical_data$Hepatic_steatosis_index > 36, yes = 1, no = 0)
clinical_data$Fatty_liver_index_prediction <- ifelse(clinical_data$Fatty_liver_index >= 60, yes = 1, no = 0)
clinical_data$APRI_score_prediction <- ifelse(clinical_data$apri_score >= 0.7, yes = 1, no = 0)


#------------------------------------------------------#
#   6. Merge Clinical & Plasma Data (Baseline/Followup)#
#------------------------------------------------------#

# Baseline sets
f2_data_initial <- merge(
  clinical_data[clinical_data$sample_group == "initial", ], 
  t(plasma_f2), 
  by.x = "unique_identifier", 
  by.y = "row.names"
)
masld_data_initial <- merge(
  clinical_data[clinical_data$sample_group == "initial", ], 
  t(plasma_masld), 
  by.x = "unique_identifier", 
  by.y = "row.names"
)
mini_f2_data_initial <- merge(
  clinical_data[clinical_data$sample_group == "initial", ], 
  t(plasma_mini_f2), 
  by.x = "unique_identifier", 
  by.y = "row.names"
)
mini_masld_data_initial <- merge(
  clinical_data[clinical_data$sample_group == "initial", ], 
  t(plasma_mini_masld), 
  by.x = "unique_identifier", 
  by.y = "row.names"
)

# Follow-up sets
f2_data_followup <- merge(
  clinical_data[clinical_data$sample_group == "followup", ], 
  t(plasma_f2), 
  by.x = "unique_identifier", 
  by.y = "row.names"
)
masld_data_followup <- merge(
  clinical_data[clinical_data$sample_group == "followup", ], 
  t(plasma_masld), 
  by.x = "unique_identifier", 
  by.y = "row.names"
)
mini_f2_data_followup <- merge(
  clinical_data[clinical_data$sample_group == "followup", ], 
  t(plasma_mini_f2), 
  by.x = "unique_identifier", 
  by.y = "row.names"
)
mini_masld_data_followup <- merge(
  clinical_data[clinical_data$sample_group == "followup", ], 
  t(plasma_mini_masld), 
  by.x = "unique_identifier", 
  by.y = "row.names"
)


#-------------------------------#
#   7. Set Cross-Validation     #
#-------------------------------#

k_folds <- 5       # Number of folds
n_reps <- 20       # Number of repetitions


#---------------------------------------#
#   8. Cross-Validated ROC Calculations #
#---------------------------------------#
# Example metrics for F2 (Fibrosis) and MASLD tests.
# Each block calls 'cross_fold_calcs_roc()' from cross_validation_utils.R

# F2: Compare APRI, FIB4, full proteomics, and simplified proteomics models
metrics_F2_APRI <- data.frame(
  melt(as.matrix(cross_fold_calcs_roc(
    dataset  = f2_data_initial[, c(10, 4)],
    dataset2 = f2_data_followup[, c(10, 4)],
    k_folds  = k_folds,
    reps     = n_reps
  )$output_combined)),
  Test   = "F2",
  Metric = "APRI model"
)

metrics_F2_FIB4 <- data.frame(
  melt(as.matrix(cross_fold_calcs_roc(
    dataset  = f2_data_initial[, c(10, 3)],
    dataset2 = f2_data_followup[, c(10, 3)],
    k_folds  = k_folds,
    reps     = n_reps
  )$output_combined)),
  Test   = "F2",
  Metric = "FIB4 model"
)

metrics_F2_Proteomics <- data.frame(
  melt(as.matrix(cross_fold_calcs_roc(
    dataset  = f2_data_initial[, c(10, 16:ncol(f2_data_initial))],
    dataset2 = f2_data_followup[, c(10, 16:ncol(f2_data_followup))],
    k_folds  = k_folds,
    reps     = n_reps
  )$output_combined)),
  Test   = "F2",
  Metric = "Proteomics model"
)

metrics_F2_Proteomics_mini <- data.frame(
  melt(as.matrix(cross_fold_calcs_roc(
    dataset  = mini_f2_data_initial[, c(10, 16:ncol(mini_f2_data_initial))],
    dataset2 = mini_f2_data_followup[, c(10, 16:ncol(mini_f2_data_followup))],
    k_folds  = k_folds,
    reps     = n_reps
  )$output_combined)),
  Test   = "F2",
  Metric = "Proteomics simplified model"
)

# MASLD: Compare APRI, FIB4, full proteomics, and simplified proteomics models
metrics_MASLD_APRI <- data.frame(
  melt(as.matrix(cross_fold_calcs_roc(
    dataset  = masld_data_initial[, c(11, 4)],
    dataset2 = masld_data_followup[, c(11, 4)],
    k_folds  = k_folds,
    reps     = n_reps
  )$output_combined)),
  Test   = "MASLD",
  Metric = "APRI model"
)

metrics_MASLD_FIB4 <- data.frame(
  melt(as.matrix(cross_fold_calcs_roc(
    dataset  = masld_data_initial[, c(11, 3)],
    dataset2 = masld_data_followup[, c(11, 3)],
    k_folds  = k_folds,
    reps     = n_reps
  )$output_combined)),
  Test   = "MASLD",
  Metric = "FIB4 model"
)

metrics_MASLD_Proteomics <- data.frame(
  melt(as.matrix(cross_fold_calcs_roc(
    dataset  = masld_data_initial[, c(11, 16:ncol(masld_data_initial))],
    dataset2 = masld_data_followup[, c(11, 16:ncol(masld_data_followup))],
    k_folds  = k_folds,
    reps     = n_reps
  )$output_combined)),
  Test   = "MASLD",
  Metric = "Proteomics model"
)

metrics_MASLD_Proteomics_mini <- data.frame(
  melt(as.matrix(cross_fold_calcs_roc(
    dataset  = mini_masld_data_initial[, c(11, 16:ncol(mini_masld_data_initial))],
    dataset2 = mini_masld_data_followup[, c(11, 16:ncol(mini_masld_data_followup))],
    k_folds  = k_folds,
    reps     = n_reps
  )$output_combined)),
  Test   = "MASLD",
  Metric = "Proteomics simplified model"
)


#----------------------------------------------------------#
#   9. Simple Metric Calculations (APRI, HSI, FLI, etc.)   #
#----------------------------------------------------------#

# Split clinical data into baseline vs followup, removing NAs
clinical_data_initial <- clinical_data[clinical_data$sample_group == "initial", ]
clinical_data_followup <- clinical_data[clinical_data$sample_group == "followup", ]

# Example: F2 using APRI
metrics_F2_APRI_initial <- data.frame(
  melt(as.matrix(simple_metric_calcs(clinical_data_initial[, c(10, 15)]))),
  Test   = "F2",
  Metric = "APRI"
)
metrics_F2_APRI_followup <- data.frame(
  melt(as.matrix(simple_metric_calcs(clinical_data_followup[, c(10, 15)]))),
  Test   = "F2",
  Metric = "APRI"
)
metrics_F2_APRI_followup$Var1 <- paste0("EXTRA_", metrics_F2_APRI_followup$Var1)

# Example: MASLD using classical scores
metrics_MASLD_NLF_initial <- data.frame(
  melt(as.matrix(simple_metric_calcs(na.omit(clinical_data_initial[, c(11, 12)])))),
  Test   = "MASLD",
  Metric = "NLF"
)
metrics_MASLD_HSI_initial <- data.frame(
  melt(as.matrix(simple_metric_calcs(na.omit(clinical_data_initial[, c(11, 13)])))),
  Test   = "MASLD",
  Metric = "HSI"
)
metrics_MASLD_FLI_initial <- data.frame(
  melt(as.matrix(simple_metric_calcs(na.omit(clinical_data_initial[, c(11, 14)])))),
  Test   = "MASLD",
  Metric = "FLI"
)

metrics_MASLD_NLF_followup <- data.frame(
  melt(as.matrix(simple_metric_calcs(na.omit(clinical_data_followup[, c(11, 12)])))),
  Test   = "MASLD",
  Metric = "NLF"
)
metrics_MASLD_HSI_followup <- data.frame(
  melt(as.matrix(simple_metric_calcs(na.omit(clinical_data_followup[, c(11, 13)])))),
  Test   = "MASLD",
  Metric = "HSI"
)
metrics_MASLD_FLI_followup <- data.frame(
  melt(as.matrix(simple_metric_calcs(na.omit(clinical_data_followup[, c(11, 14)])))),
  Test   = "MASLD",
  Metric = "FLI"
)

# Tag follow-up rows distinctly
metrics_MASLD_NLF_followup$Var1 <- paste0("EXTRA_", metrics_MASLD_NLF_followup$Var1)
metrics_MASLD_HSI_followup$Var1 <- paste0("EXTRA_", metrics_MASLD_HSI_followup$Var1)
metrics_MASLD_FLI_followup$Var1 <- paste0("EXTRA_", metrics_MASLD_FLI_followup$Var1)


#---------------------------------------------#
#   10. Combine All Simple Metric Data        #
#---------------------------------------------#

combined_simple <- rbind(
  metrics_F2_Proteomics_mini,    metrics_F2_APRI,            metrics_F2_FIB4,           metrics_F2_Proteomics,
  metrics_F2_APRI_initial,       metrics_F2_APRI_followup,
  metrics_MASLD_Proteomics_mini, metrics_MASLD_APRI,         metrics_MASLD_FIB4,        metrics_MASLD_Proteomics,
  metrics_MASLD_NLF_initial,     metrics_MASLD_HSI_initial,  metrics_MASLD_FLI_initial,
  metrics_MASLD_NLF_followup,    metrics_MASLD_HSI_followup, metrics_MASLD_FLI_followup
)

combined_simple$Group <- ifelse(grepl("EXTRA_", combined_simple$Var1), 
                                yes = "Follow-up", 
                                no = "Test data")
combined_simple$Var1 <- gsub("EXTRA_", "", combined_simple$Var1)

combined_simple_acc <- combined_simple[combined_simple$Var1 %in% c("Balanced_accuracy", "F1_Score"), ]
combined_simple_acc$Metric <- fct_relevel(
  combined_simple_acc$Metric,
  c("APRI", "FLI", "HSI", "NLF", "APRI model", "FIB4 model", 
    "Proteomics simplified model", "Proteomics model")
)
combined_simple_acc$Group <- fct_relevel(combined_simple_acc$Group, c("Test data", "Follow-up"))
combined_simple_acc$Test <- fct_relevel(combined_simple_acc$Test, c("MASLD", "F2"))


#---------------------#
#   11. Plot Metrics  #
#---------------------#

# Example performance plots
pdf("plots/Predicting MASLD_performance metrics.pdf", width = 120/25.4, height = 40/25.4)
ggplot() +
  geom_col(
    data = combined_simple_acc[combined_simple_acc$Var2 == "Mean" & combined_simple_acc$Test == "MASLD", ],
    aes(y = Metric, x = value, fill = Metric),
    color = "black", width = 0.7
  ) +
  geom_jitter(
    data = combined_simple_acc[!(combined_simple_acc$Var2 %in% c("SD", "Mean")) & combined_simple_acc$Test == "MASLD", ],
    aes(y = Metric, x = value, fill = Metric),
    color = "black", height = 0.1, size = 0.1, alpha = 0.2, shape = 16
  ) +
  geom_errorbar(
    data = combined_simple_acc[combined_simple_acc$Var2 == "SD" & combined_simple_acc$Test == "MASLD", ],
    aes(
      xmin = -value + combined_simple_acc$value[combined_simple_acc$Var2 == "Mean" & combined_simple_acc$Test == "MASLD"],
      xmax =  value + combined_simple_acc$value[combined_simple_acc$Var2 == "Mean" & combined_simple_acc$Test == "MASLD"],
      y = Metric
    ),
    width = 0.3
  ) +
  theme_bw() + 
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    strip.background = element_rect(fill = "white"),
    axis.text = element_text(size = 5),
    strip.text = element_blank(),
    panel.spacing = unit(1, "lines"),
    axis.title = element_blank()
  ) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  facet_wrap(Group ~ Test + Var1, nrow = 2, scales = "fixed") +
  scale_fill_manual(values = c(
    "Proteomics model"            = "#860e25",
    "Proteomics simplified model" = "#B10E35",
    "FIB4 model"                  = "#324e7b",
    "APRI model"                  = "#3A86FF",
    "APRI"                        = "#FFBE0B",
    "HSI"                         = "#FFBE0B",
    "FLI"                         = "#FFBE0B",
    "NLF"                         = "#FFBE0B"
  ))
dev.off()


pdf("plots/Predicting ≥F2_performance metrics.pdf", width = 120/25.4, height = 32.5/25.4)
ggplot() +
  geom_col(
    data = combined_simple_acc[combined_simple_acc$Var2 == "Mean" & combined_simple_acc$Test == "F2", ],
    aes(y = Metric, x = value, fill = Metric),
    color = "black", width = 0.7
  ) +
  geom_jitter(
    data = combined_simple_acc[!(combined_simple_acc$Var2 %in% c("SD", "Mean")) & combined_simple_acc$Test == "F2", ],
    aes(y = Metric, x = value, fill = Metric),
    color = "black", height = 0.1, size = 0.1, alpha = 0.2, shape = 16
  ) +
  geom_errorbar(
    data = combined_simple_acc[combined_simple_acc$Var2 == "SD" & combined_simple_acc$Test == "F2", ],
    aes(
      xmin = -value + combined_simple_acc$value[combined_simple_acc$Var2 == "Mean" & combined_simple_acc$Test == "F2"],
      xmax =  value + combined_simple_acc$value[combined_simple_acc$Var2 == "Mean" & combined_simple_acc$Test == "F2"],
      y = Metric
    ),
    width = 0.3
  ) +
  theme_bw() + 
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    strip.background = element_rect(fill = "white"),
    axis.text = element_text(size = 5),
    strip.text = element_blank(),
    panel.spacing = unit(1, "lines"),
    axis.title = element_blank()
  ) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  facet_wrap(Group ~ Test + Var1, nrow = 2, scales = "fixed") +
  scale_fill_manual(values = c(
    "Proteomics model"            = "#860e25",
    "Proteomics simplified model" = "#B10E35",
    "FIB4 model"                  = "#324e7b",
    "APRI model"                  = "#3A86FF",
    "APRI"                        = "#FFBE0B",
    "HSI"                         = "#FFBE0B",
    "FLI"                         = "#FFBE0B",
    "NLF"                         = "#FFBE0B"
  ))
dev.off()


#----------------------------------------------------#
#   12. Additional ROC Analysis (Initial/Followup)   #
#----------------------------------------------------#

# Recompute AUC-based metrics for easier plotting
metrics_F2_APRI        <- data.frame(cross_fold_calcs_roc(f2_data_initial[, c(10, 4)], f2_data_followup[, c(10, 4)], k_folds, n_reps)$avg_roc,  Test = "F2",   Metric = "APRI model")
metrics_F2_FIB4        <- data.frame(cross_fold_calcs_roc(f2_data_initial[, c(10, 3)], f2_data_followup[, c(10, 3)], k_folds, n_reps)$avg_roc,  Test = "F2",   Metric = "FIB4 model")
metrics_F2_Proteomics  <- data.frame(cross_fold_calcs_roc(f2_data_initial[, c(10, 16:35)], f2_data_followup[, c(10, 16:35)], k_folds, n_reps)$avg_roc,  Test = "F2",   Metric = "Proteomics model")
metrics_F2_Prot_mini   <- data.frame(cross_fold_calcs_roc(mini_f2_data_initial[, c(10, 16:ncol(mini_f2_data_initial))], mini_f2_data_followup[, c(10, 16:ncol(mini_f2_data_followup))], k_folds, n_reps)$avg_roc, Test = "F2",   Metric = "Proteomics simplified model")

metrics_MASLD_APRI      <- data.frame(cross_fold_calcs_roc(masld_data_initial[, c(11, 4)], masld_data_followup[, c(11, 4)], k_folds, n_reps)$avg_roc,  Test = "MASLD", Metric = "APRI model")
metrics_MASLD_FIB4      <- data.frame(cross_fold_calcs_roc(masld_data_initial[, c(11, 3)], masld_data_followup[, c(11, 3)], k_folds, n_reps)$avg_roc,  Test = "MASLD", Metric = "FIB4 model")
metrics_MASLD_Proteomics <- data.frame(cross_fold_calcs_roc(masld_data_initial[, c(11, 16:35)], masld_data_followup[, c(11, 16:35)], k_folds, n_reps)$avg_roc,  Test = "MASLD", Metric = "Proteomics model")
metrics_MASLD_Prot_mini  <- data.frame(cross_fold_calcs_roc(mini_masld_data_initial[, c(11, 16:ncol(mini_masld_data_initial))], mini_masld_data_followup[, c(11, 16:ncol(mini_masld_data_followup))], k_folds, n_reps)$avg_roc, Test = "MASLD", Metric = "Proteomics simplified model")

combined_initial <- rbind(
  metrics_MASLD_Prot_mini, metrics_F2_Prot_mini,
  metrics_F2_APRI, metrics_F2_FIB4, metrics_F2_Proteomics,
  metrics_MASLD_APRI, metrics_MASLD_FIB4, metrics_MASLD_Proteomics
)
combined_initial$Cohort <- "Initial"

# Repeat for avg_roc2 (follow-up data)
metrics_F2_APRI        <- data.frame(cross_fold_calcs_roc(f2_data_initial[, c(10, 4)], f2_data_followup[, c(10, 4)], k_folds, n_reps)$avg_roc2, Test = "F2",   Metric = "APRI model")
metrics_F2_FIB4        <- data.frame(cross_fold_calcs_roc(f2_data_initial[, c(10, 3)], f2_data_followup[, c(10, 3)], k_folds, n_reps)$avg_roc2, Test = "F2",   Metric = "FIB4 model")
metrics_F2_Proteomics  <- data.frame(cross_fold_calcs_roc(f2_data_initial[, c(10, 16:35)], f2_data_followup[, c(10, 16:35)], k_folds, n_reps)$avg_roc2, Test = "F2",   Metric = "Proteomics model")
metrics_F2_Prot_mini   <- data.frame(cross_fold_calcs_roc(mini_f2_data_initial[, c(10, 16:ncol(mini_f2_data_initial))], mini_f2_data_followup[, c(10, 16:ncol(mini_f2_data_followup))], k_folds, n_reps)$avg_roc2, Test = "F2",   Metric = "Proteomics simplified model")

metrics_MASLD_APRI      <- data.frame(cross_fold_calcs_roc(masld_data_initial[, c(11, 4)], masld_data_followup[, c(11, 4)], k_folds, n_reps)$avg_roc2, Test = "MASLD", Metric = "APRI model")
metrics_MASLD_FIB4      <- data.frame(cross_fold_calcs_roc(masld_data_initial[, c(11, 3)], masld_data_followup[, c(11, 3)], k_folds, n_reps)$avg_roc2, Test = "MASLD", Metric = "FIB4 model")
metrics_MASLD_Proteomics <- data.frame(cross_fold_calcs_roc(masld_data_initial[, c(11, 16:35)], masld_data_followup[, c(11, 16:35)], k_folds, n_reps)$avg_roc2, Test = "MASLD", Metric = "Proteomics model")
metrics_MASLD_Prot_mini  <- data.frame(cross_fold_calcs_roc(mini_masld_data_initial[, c(11, 16:ncol(mini_masld_data_initial))], mini_masld_data_followup[, c(11, 16:ncol(mini_masld_data_followup))], k_folds, n_reps)$avg_roc2, Test = "MASLD", Metric = "Proteomics simplified model")

combined_followup <- rbind(
  metrics_MASLD_Prot_mini, metrics_F2_Prot_mini,
  metrics_F2_APRI, metrics_F2_FIB4, metrics_F2_Proteomics,
  metrics_MASLD_APRI, metrics_MASLD_FIB4, metrics_MASLD_Proteomics
)
combined_followup$Cohort <- "Followup"

combined_roc <- rbind(combined_initial, combined_followup)

# Reorder factor levels for plotting
combined_roc$Metric <- fct_relevel(
  combined_roc$Metric,
  c("APRI model", "FIB4 model", "Proteomics simplified model", "Proteomics model")
)
combined_roc$Cohort <- fct_relevel(combined_roc$Cohort, c("Initial", "Followup"))
combined_roc$Test   <- fct_relevel(combined_roc$Test,   c("MASLD", "F2"))


#--------------------------------------------#
#   13. Example ROC Plot for Baseline Data   #
#--------------------------------------------#

pdf("plots/Predicting MASLD (Initial cohort)_roc curves.pdf", width = 50/25.4, height = 40/25.4)
ggplot(combined_roc[combined_roc$Cohort == "Initial" & combined_roc$Test == "MASLD", ]) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray", size = 0.5) +
  geom_ribbon(
    aes(
      x = mean_specificities, ymin = ci_lower_sensitivities, ymax = ci_upper_sensitivities,
      group = Metric, fill = Metric
    ), 
    alpha = 0.25
  ) +
  geom_ribbon(
    aes(
      y = mean_sensitivities, xmin = ci_lower_specificities, xmax = ci_upper_specificities,
      group = Metric, fill = Metric
    ), 
    alpha = 0.25
  ) +
  geom_line(aes(x = mean_specificities, y = mean_sensitivities, group = Metric, color = Metric)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 5),
    axis.title = element_text(size = 5),
    panel.spacing = unit(1, "lines"),
    legend.title = element_blank(),
    strip.background = element_blank()
  ) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0), breaks = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0), breaks = c(0, 1)) +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  scale_color_manual(values = c(
    "Proteomics model"            = "#860e25",
    "Proteomics simplified model" = "#B10E35",
    "FIB4 model"                  = "#324e7b",
    "APRI model"                  = "#3A86FF"
  )) +
  scale_fill_manual(values = c(
    "Proteomics model"            = "#860e25",
    "Proteomics simplified model" = "#B10E35",
    "FIB4 model"                  = "#324e7b",
    "APRI model"                  = "#3A86FF"
  ))
dev.off()

#--------------------------------------------#
#   14. Example ROC Plot for F2 Baseline     #
#--------------------------------------------#

pdf("plots/Predicting ≥F2 (Initial cohort)_roc curves.pdf", width = 50/25.4, height = 40/25.4)
ggplot(combined_roc[combined_roc$Cohort == "Initial" & combined_roc$Test == "F2", ]) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray", size = 0.5) +
  geom_ribbon(
    aes(
      x = mean_specificities, ymin = ci_lower_sensitivities, ymax = ci_upper_sensitivities,
      group = Metric, fill = Metric
    ), 
    alpha = 0.25
  ) +
  geom_ribbon(
    aes(
      y = mean_sensitivities, xmin = ci_lower_specificities, xmax = ci_upper_specificities,
      group = Metric, fill = Metric
    ), 
    alpha = 0.25
  ) +
  geom_line(aes(x = mean_specificities, y = mean_sensitivities, group = Metric, color = Metric)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 5),
    axis.title = element_text(size = 5),
    panel.spacing = unit(1, "lines"),
    legend.title = element_blank(),
    strip.background = element_blank()
  ) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0), breaks = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0), breaks = c(0, 1)) +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  scale_color_manual(values = c(
    "Proteomics model"            = "#860e25",
    "Proteomics simplified model" = "#B10E35",
    "FIB4 model"                  = "#324e7b",
    "APRI model"                  = "#3A86FF"
  )) +
  scale_fill_manual(values = c(
    "Proteomics model"            = "#860e25",
    "Proteomics simplified model" = "#B10E35",
    "FIB4 model"                  = "#324e7b",
    "APRI model"                  = "#3A86FF"
  ))
dev.off()

print(unique(combined_roc[combined_roc$Cohort == "Initial", c(12,10,11,3)]))


###############################################################################
# End of script
###############################################################################
