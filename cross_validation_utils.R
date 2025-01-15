#' Perform Cross-Validated ROC Calculations
#'
#' This function conducts k-fold cross-validation with multiple repetitions on two datasets.
#' It trains logistic regression models, evaluates performance metrics, and aggregates ROC curves.
#'
#' @param dataset The primary dataset for cross-validation.
#' @param dataset2 An additional dataset for evaluation.
#' @param k_folds The number of folds for cross-validation.
#' @param reps The number of repetitions for cross-validation.
#'
#' @return A list containing combined performance metrics, ROC curves for both datasets,
#'         and averaged ROC metrics with confidence intervals.
#'
#' @import caret
#' @import pROC
#' @importFrom stats glm sd
#' @export
cross_fold_calcs_roc <- function(dataset, dataset2, k_folds, reps) {
  
  # Initialize datasets
  your_data <- dataset
  your_extra_data <- dataset2
  num_folds <- k_folds
  num_reps <- reps
  
  # Generate unique seed numbers for reproducibility
  seed_num <- 100:(100 + num_reps)
  
  # Initialize lists to store ROC curves
  roc_curves <- list()
  roc_curves2 <- list()
  
  # Iterate over each repetition
  for (d in 1:num_reps) {
    set.seed(seed_num[d])  # Set seed for reproducibility
    
    # Data Cleaning for primary dataset
    your_data <- your_data[!your_data[, 2] == Inf, ]
    your_data <- your_data[, colSums(is.na(your_data)) < 10]
    your_data <- na.omit(your_data)
    colnames(your_data)[1] <- "target_variable"
    if (ncol(your_data) > 2) {
      your_data[, 2:ncol(your_data)] <- lapply(your_data[, 2:ncol(your_data)], as.numeric)
    }
    
    # Data Cleaning for extra dataset
    your_extra_data <- your_extra_data[!your_extra_data[, 2] == Inf, ]
    your_extra_data <- your_extra_data[, colSums(is.na(your_extra_data)) < 10]
    your_extra_data <- na.omit(your_extra_data)
    colnames(your_extra_data)[1] <- "target_variable"
    if (ncol(your_extra_data) > 2) {
      your_extra_data[, 2:ncol(your_extra_data)] <- lapply(your_extra_data[, 2:ncol(your_extra_data)], as.numeric)
    }
    
    # Create indices for k-fold cross-validation
    folds <- createFolds(your_data$target_variable, k = num_folds, list = TRUE, returnTrain = FALSE)
    
    # Initialize vectors to store evaluation metrics for primary and extra datasets
    cv_results <- numeric(num_folds)
    roc_auc_values <- numeric(num_folds)
    precision_values <- numeric(num_folds)
    recall_values <- numeric(num_folds)
    f1_score_values <- numeric(num_folds)
    balanced_accuracy_values <- numeric(num_folds)
    
    cv_results2 <- numeric(num_folds)
    roc_auc_values2 <- numeric(num_folds)
    precision_values2 <- numeric(num_folds)
    recall_values2 <- numeric(num_folds)
    f1_score_values2 <- numeric(num_folds)
    balanced_accuracy_values2 <- numeric(num_folds)
    
    # Perform cross-validation across all folds
    for (i in 1:num_folds) {
      # Split the data into training and testing sets
      train_data <- your_data[-folds[[i]], ]
      test_data <- your_data[folds[[i]], ]
      
      # Train logistic regression model
      model <- glm(target_variable ~ ., data = train_data, family = "binomial")
      
      # Make predictions on test and extra datasets
      predictions_test <- predict(model, newdata = test_data, type = "response")
      predictions_extra <- predict(model, newdata = your_extra_data, type = "response")
      
      # Generate predicted classes based on a threshold of 0.5
      predicted_classes_test <- ifelse(predictions_test >= 0.5, 1, 0)
      predicted_classes_extra <- ifelse(predictions_extra >= 0.5, 1, 0)
      
      # Create confusion matrices for test and extra datasets
      confusion_matrix <- table(
        Actual = factor(test_data$target_variable, levels = c(0, 1)),
        Predicted = factor(predicted_classes_test, levels = c(0, 1))
      )
      
      confusion_matrix_extra <- table(
        Actual = factor(your_extra_data$target_variable, levels = c(0, 1)),
        Predicted = factor(predicted_classes_extra, levels = c(0, 1))
      )
      
      # Calculate performance metrics for test dataset
      accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
      precision <- ifelse(sum(confusion_matrix[, 2]) == 0, NA, confusion_matrix[2, 2] / sum(confusion_matrix[, 2]))
      recall <- ifelse(sum(confusion_matrix[2, ]) == 0, NA, confusion_matrix[2, 2] / sum(confusion_matrix[2, ]))
      f1_score <- ifelse(is.na(precision) | is.na(recall) | (precision + recall) == 0, NA, 2 * (precision * recall) / (precision + recall))
      balanced_accuracy <- mean(c(
        ifelse(sum(confusion_matrix[1, ]) == 0, NA, confusion_matrix[1, 1] / sum(confusion_matrix[1, ])),
        ifelse(sum(confusion_matrix[2, ]) == 0, NA, confusion_matrix[2, 2] / sum(confusion_matrix[2, ]))
      ), na.rm = TRUE)
      
      # Calculate performance metrics for extra dataset
      accuracy2 <- sum(diag(confusion_matrix_extra)) / sum(confusion_matrix_extra)
      precision2 <- ifelse(sum(confusion_matrix_extra[, 2]) == 0, NA, confusion_matrix_extra[2, 2] / sum(confusion_matrix_extra[, 2]))
      recall2 <- ifelse(sum(confusion_matrix_extra[2, ]) == 0, NA, confusion_matrix_extra[2, 2] / sum(confusion_matrix_extra[2, ]))
      f1_score2 <- ifelse(is.na(precision2) | is.na(recall2) | (precision2 + recall2) == 0, NA, 2 * (precision2 * recall2) / (precision2 + recall2))
      balanced_accuracy2 <- mean(c(
        ifelse(sum(confusion_matrix_extra[1, ]) == 0, NA, confusion_matrix_extra[1, 1] / sum(confusion_matrix_extra[1, ])),
        ifelse(sum(confusion_matrix_extra[2, ]) == 0, NA, confusion_matrix_extra[2, 2] / sum(confusion_matrix_extra[2, ]))
      ), na.rm = TRUE)
      
      # Store the calculated metrics
      cv_results[i] <- accuracy
      precision_values[i] <- precision
      recall_values[i] <- recall
      f1_score_values[i] <- f1_score
      balanced_accuracy_values[i] <- balanced_accuracy
      
      cv_results2[i] <- accuracy2
      precision_values2[i] <- precision2
      recall_values2[i] <- recall2
      f1_score_values2[i] <- f1_score2
      balanced_accuracy_values2[i] <- balanced_accuracy2
      
      # Generate ROC curves and store them
      roc_curve <- roc(test_data$target_variable, predictions_test)
      roc_curves[[length(roc_curves) + 1]] <- roc_curve
      
      roc_curve2 <- roc(your_extra_data$target_variable, predictions_extra)
      roc_curves2[[length(roc_curves2) + 1]] <- roc_curve2
    }
    
    # Compile performance metrics into dataframes
    output <- data.frame(
      Accuracy = cv_results,
      Balanced_accuracy = balanced_accuracy_values,
      Precision = precision_values,
      Recall = recall_values,
      F1_Score = f1_score_values
    )
    rownames(output) <- paste0("N", d, "_", "K", 1:nrow(output))
    output <- as.data.frame(t(output))
    
    output2 <- data.frame(
      EXTRA_Accuracy = cv_results2,
      EXTRA_Balanced_accuracy = balanced_accuracy_values2,
      EXTRA_Precision = precision_values2,
      EXTRA_Recall = recall_values2,
      EXTRA_F1_Score = f1_score_values2
    )
    rownames(output2) <- paste0("N", d, "_", "K", 1:nrow(output2))
    output2 <- as.data.frame(t(output2))
    
    # Combine primary and extra dataset metrics
    output3 <- rbind(output, output2)
    
    if (d == 1) {
      output_combined <- output3
    } else {
      output_combined <- cbind(output_combined, output3)
    }
  }
  
  # Calculate standard deviation across all repetitions and folds
  length_df <- ncol(output_combined)
  
  output_combined$Mean <- rowMeans(output_combined, na.rm = TRUE)
  output_combined$SD <- apply(output_combined[, 1:length_df], 1, sd, na.rm = TRUE)
  
  # Compute average ROC metrics with confidence intervals
  avg_roc <- roc.curve(roc_curves)
  avg_roc2 <- roc.curve(roc_curves2)
  
  # Return all aggregated results
  return(list(
    output_combined = output_combined,
    roc_curves = roc_curves,
    roc_curves2 = roc_curves2,
    avg_roc = avg_roc,
    avg_roc2 = avg_roc2
  ))
}