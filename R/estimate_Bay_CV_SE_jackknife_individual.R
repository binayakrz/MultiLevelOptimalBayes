#' Jackknife estimator (individuals deletion) of Standard Errors of  Between-Group Effect for 2-level latent variable model with control variables estimated with regularized_Bayesian_estimator
#'
#' This function performs a jackknife estimation of standard errors by individually deleting one observation
#' from each group in a two-level latent variable model with control variables. The model uses a regularized
#' Bayesian estimator to compute between-group effects. For each jackknife replication, one individual is 
#' randomly deleted from each group, and the model is re-estimated to provide an estimation of the standard errors 
#' through the jackknife technique.
#'
#' Two types of standard errors are calculated:
#'   - `SE_beta_Bay_jackknife_individual`: The jackknife standard error for the Bayesian estimator.
#'   - `SE_beta_Bay_ML_jackknife_individual`: The jackknife standard error for the Maximum Likelihood estimator.
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
#'     \item `SE_beta_Bay_jackknife_individual`: The jackknife standard error for the Bayesian estimator.
#'     \item `SE_beta_Bay_ML_jackknife_individual`: The jackknife standard error for the Maximum Likelihood estimator.
#'   }
#' @examples
#' # Assuming data is a list structured as per the requirements
#' result_Bay_CV_SE_jackknife_individual <- estimate_Bay_CV_SE_jackknife_individual(data)
#' print(result_Bay_CV_SE_jackknife_individual)
#' @importFrom pracma sqrt
#' @export

estimate_Bay_CV_SE_jackknife_individual <- function(data) {
  library(pracma)
  
  # Original estimator calculation (to get initial results)
  original_result <- estimate_Bay_CV(data)
  
  # Initialize parameters
  n <- data$n
  J <- data$k
  
  # Initialize vectors to store jackknife estimates for individual deletion
  jackknife_beta_b_Bay_individual <- numeric(n)
  jackknife_beta_b_Bay_ML_individual <- numeric(n)
  
  r = n # set up the number of replications for jackknife, r from 1 to n^J
  # Jackknife by deleting one individual from each group simultaneously
  for (i in 1:r) {
    # Generate J random indices (one for each group) from 1 to n
    random_indices <- sample(1:n, J, replace = TRUE)
    
    # Initialize a list to store indices to keep for each group
    indices_to_keep <- vector("list", J)
    
    for (j in 1:J) {
      indices_to_keep[[j]] <- (n * (j - 1) + 1):(n * j)
      indices_to_keep[[j]] <- indices_to_keep[[j]][-random_indices[j]]
    }
    
    # Create new indices after deleting one individual from each group
    indices <- unlist(indices_to_keep)
    
    # Create a new dataset excluding the selected individuals
    data_jackknife <- data
    data_jackknife$x <- data$x[indices]
    data_jackknife$y <- data$y[indices]
    if (data$kc>1){
      data_jackknife$C <- data$C[indices, , drop = FALSE]
    } else {
      data_jackknife$C <- data$C[indices, drop = FALSE]
    }
    data_jackknife$n <- n - 1
    data_jackknife$kn <- (n-1)*J
    
    # Recalculate estimators without the selected individuals
    result_jackknife <- estimate_Bay_CV(data_jackknife)
    
    # Store jackknife estimates
    jackknife_beta_b_Bay_individual[i] <- result_jackknife$beta_b_Bay
    jackknife_beta_b_Bay_ML_individual[i] <- result_jackknife$beta_b_Bay_ML
  }
  
  # Calculate jackknife means for individual deletion
  jackknife_mean_beta_b_Bay_individual <- mean(jackknife_beta_b_Bay_individual)
  jackknife_mean_beta_b_Bay_ML_individual <- mean(jackknife_beta_b_Bay_ML_individual)
  
  # Calculate jackknife standard errors for individual deletion
  SE_beta_Bay_jackknife_individual <- sqrt((n - 1) / n * sum((jackknife_beta_b_Bay_individual - jackknife_mean_beta_b_Bay_individual)^2))
  SE_beta_Bay_ML_jackknife_individual <- sqrt((n - 1) / n * sum((jackknife_beta_b_Bay_ML_individual - jackknife_mean_beta_b_Bay_ML_individual)^2))
  
  # Output updated results with jackknife standard errors
  list(
    beta_b_Bay = original_result$beta_b_Bay,
    beta_b_Bay_ML = original_result$beta_b_Bay_ML,
    SE_beta_Bay_jackknife_individual = SE_beta_Bay_jackknife_individual,
    SE_beta_Bay_ML_jackknife_individual = SE_beta_Bay_ML_jackknife_individual
    #original_result = original_result
  )
}

# Example usage with your data
#result_Bay_CV_SE_jackknife_individual <- estimate_Bay_CV_SE_jackknife_individual(data_CV)
#print(result_Bay_CV_SE_jackknife_individual)
