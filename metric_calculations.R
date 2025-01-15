#' Calculate Simple Classification Metrics
#'
#' This function computes basic classification metrics from a dataframe containing
#' actual and predicted binary values.
#'
#' @param dataframe A dataframe where the first column contains actual values
#'                  and the second column contains predicted values.
#'
#' @return A dataframe with calculated metrics: Accuracy, Balanced Accuracy,
#'         Precision, Recall, and F1-Score. The standard deviation (SD) column
#'         is set to -100 as a placeholder.
#'
#' @export
simple_metric_calcs <- function(dataframe) {
  
  # Extract actual and predicted values, ensuring they are factors with correct levels
  actual <- factor(as.numeric(unlist(dataframe[, 1])), levels = c(0, 1))
  predictor <- factor(as.numeric(unlist(dataframe[, 2])), levels = c(0, 1))
  
  # Create confusion matrix
  confusion_matrix <- table(
    Actual = actual,
    Predicted = predictor
  )
  
  # Calculate Balanced Accuracy
  balanced_accuracy <- mean(c(
    ifelse(sum(confusion_matrix[1, ]) == 0, NA, confusion_matrix[1, 1] / sum(confusion_matrix[1, ])),
    ifelse(sum(confusion_matrix[2, ]) == 0, NA, confusion_matrix[2, 2] / sum(confusion_matrix[2, ]))
  ), na.rm = TRUE)
  
  # Calculate Precision
  precision <- ifelse(sum(confusion_matrix[, 2]) == 0, NA, confusion_matrix[2, 2] / sum(confusion_matrix[, 2]))
  
  # Calculate Recall
  recall <- ifelse(sum(confusion_matrix[2, ]) == 0, NA, confusion_matrix[2, 2] / sum(confusion_matrix[2, ]))
  
  # Calculate F1-Score
  f1_score <- ifelse(is.na(precision) | is.na(recall) | (precision + recall) == 0, NA, 2 * (precision * recall) / (precision + recall))
  
  # Calculate Accuracy
  accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
  
  # Compile metrics into a dataframe
  output <- data.frame(
    Accuracy = accuracy,
    Balanced_accuracy = balanced_accuracy,
    Precision = precision,
    Recall = recall,
    F1_Score = f1_score
  )
  
  # Transpose and set column names
  output <- t(output)
  colnames(output) <- "Mean"
  
  # Convert to dataframe and add SD column as a placeholder
  output <- as.data.frame(output)
  output$SD <- -100  # Placeholder value
  
  return(output)
}