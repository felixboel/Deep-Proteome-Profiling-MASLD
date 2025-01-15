#' Compute Mean ROC Metrics with Confidence Intervals
#'
#' This function processes a list of ROC curves to calculate the mean specificities,
#' sensitivities, AUC values, and their respective standard errors and 95% confidence intervals.
#'
#' @param roc_curves A list of ROC curve objects (e.g., from pROC::roc).
#'
#' @return A list containing mean specificities, mean sensitivities, mean AUC,
#'         standard errors, and confidence intervals for specificities and sensitivities.
#'
#' @importFrom stats approx qnorm sd
#' @export
roc.curve <- function(roc_curves) {
  
  # Interpolate specificities and sensitivities to 100 points
  specs <- lapply(roc_curves, function(x) approx(rev(x$specificities), n = 100)$y)
  sens <- lapply(roc_curves, function(x) approx(rev(x$sensitivities), n = 100)$y)
  
  # Extract AUC values from each ROC curve
  aucs <- sapply(roc_curves, function(x) x$auc)
  
  # Convert lists to matrices for easier computation
  specs <- 1 - do.call(cbind, specs)      # Convert specificities to false positive rates
  sens <- do.call(cbind, sens)           # Sensitivities (true positive rates)
  
  # Calculate mean specificities, sensitivities, and AUC
  mean_specs <- apply(specs, 1, mean)
  mean_sens <- apply(sens, 1, mean)
  mean_auc <- mean(aucs)
  
  # Calculate the Standard Error of the Mean (SEM) for specificities and sensitivities
  sem_specificities <- apply(specs, 1, function(x) sd(x) / sqrt(length(x)))
  sem_sensitivities <- apply(sens, 1, function(x) sd(x) / sqrt(length(x)))
  
  # Define the confidence level (e.g., 95%)
  alpha <- 0.05
  z <- qnorm(1 - alpha / 2)  # Z-score for 95% confidence
  
  # Calculate the lower and upper bounds of the confidence intervals
  ci_lower_specificities <- mean_specs - z * sem_specificities
  ci_upper_specificities <- mean_specs + z * sem_specificities
  ci_lower_sensitivities <- mean_sens - z * sem_sensitivities
  ci_upper_sensitivities <- mean_sens + z * sem_sensitivities
  
  # Return the computed metrics as a list
  return(list(
    mean_specificities = mean_specs,
    mean_sensitivities = mean_sens,
    mean_auc = mean_auc,
    sem_specificities = sem_specificities,
    sem_sensitivities = sem_sensitivities,
    ci_lower_specificities = ci_lower_specificities,
    ci_upper_specificities = ci_upper_specificities,
    ci_lower_sensitivities = ci_lower_sensitivities,
    ci_upper_sensitivities = ci_upper_sensitivities
  ))
}