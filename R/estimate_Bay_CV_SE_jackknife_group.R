#' Jackknife estimator (groupwise deletion) of standard errors of between-group effect for 2-level latent variable model with control variables estimated with regularized Bayesian estimator
#'
#' This function performs a jackknife estimation of standard errors for the between-group effects 
#' in a two-level latent variable model using regularized Bayesian estimation. It estimates the 
#' between-group effects by deleting one group at a time and recalculating the model, providing 
#' an estimation of the standard errors through the jackknife technique.
#'
#' The model handles control variables and computes two standard errors:
#' one for the Bayesian estimator (`SE_beta_Bay_jackknife_group`) and one for the 
#' Maximum Likelihood estimator (`SE_beta_Bay_ML_jackknife_group`).
#'
#' @param data A list containing the dataset information with elements:
#'   \itemize{
#'     \item `n`: Number of individuals per group.
#'     \item `k`: Number of groups.
#'     \item `x`: Predictor variable.
#'     \item `y`: Response variable.
#'     \item `C`: Control variables matrix (if applicable).
#'     \item `kc`: Number of control variables (1 or more).
#'   }
#' @return A list containing:
#'   \itemize{
#'     \item `beta_b_Bay`: The original Bayesian estimate of the between-group effect.
#'     \item `beta_b_Bay_ML`: The original Maximum Likelihood estimate of the between-group effect.
#'     \item `SE_beta_Bay_jackknife_group`: The jackknife standard error for the Bayesian estimator.
#'     \item `SE_beta_Bay_ML_jackknife_group`: The jackknife standard error for the Maximum Likelihood estimator.
#'   }
#' @examples
#' # Assuming data is a list structured as per the requirements
#' result_Bay_CV_SE_jackknife_group <- estimate_Bay_CV_SE_jackknife_group(data)
#' print(result_Bay_CV_SE_jackknife_group)
#' @importFrom pracma sqrt
#' @export

estimate_Bay_CV_SE_jackknife_group <- function(data) {
  library(pracma)
  
  # Original estimator calculation (to get initial results)
  original_result <- estimate_Bay_CV(data)
  
  # Initialize parameters
  n <- data$n
  J <- data$k
  
  # Initialize vectors to store jackknife estimates for group deletion
  jackknife_beta_b_Bay_group <- numeric(J)
  jackknife_beta_b_Bay_ML_group <- numeric(J)
  
  # Jackknife by deleting one entire group at a time
  for (j in 1:J) {
    # Create a new dataset excluding the j-th group
    data_jackknife <- data
    data_jackknife$x <- data$x[-((n * (j - 1) + 1):(n * j))]
    data_jackknife$y <- data$y[-((n * (j - 1) + 1):(n * j))]
    if (data$kc>1){
      data_jackknife$C <- data$C[-((n * (j - 1) + 1):(n * j)), , drop = FALSE]
    } else if (data$kc==1) {
      data_jackknife$C <- data$C[-((n * (j - 1) + 1):(n * j)), drop = FALSE]
    }
    
    data_jackknife$k <- J - 1
    data_jackknife$n <- n
    data_jackknife$kn <- n*(J-1)
    
    # Recalculate estimators without the j-th group
    result_jackknife <- estimate_Bay_CV(data_jackknife)
    
    # Store jackknife estimates
    jackknife_beta_b_Bay_group[j] <- result_jackknife$beta_b_Bay
    jackknife_beta_b_Bay_ML_group[j] <- result_jackknife$beta_b_Bay_ML
  }
  
  # Calculate jackknife means for group deletion
  jackknife_mean_beta_b_Bay_group <- mean(jackknife_beta_b_Bay_group)
  jackknife_mean_beta_b_Bay_ML_group <- mean(jackknife_beta_b_Bay_ML_group)
  
  # Calculate jackknife standard errors for group deletion
  SE_beta_Bay_jackknife_group <- sqrt((J - 1) / J * sum((jackknife_beta_b_Bay_group - jackknife_mean_beta_b_Bay_group)^2))
  SE_beta_Bay_ML_jackknife_group <- sqrt((J - 1) / J * sum((jackknife_beta_b_Bay_ML_group - jackknife_mean_beta_b_Bay_ML_group)^2))
  
  # Output updated results with jackknife standard errors
  list(
    beta_b_Bay = original_result$beta_b_Bay,
    beta_b_Bay_ML = original_result$beta_b_Bay_ML,
    SE_beta_Bay_jackknife_group = SE_beta_Bay_jackknife_group,
    SE_beta_Bay_ML_jackknife_group = SE_beta_Bay_ML_jackknife_group
    #original_result = original_result
  )
}

# Example usage with your data
#result_Bay_CV_SE_jackknife_group <- estimate_Bay_CV_SE_jackknife_group(data_CV)
#print(result_Bay_CV_SE_jackknife_group)
