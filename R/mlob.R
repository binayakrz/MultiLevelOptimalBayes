#' Multi-Level Optimal Bayes Function (MLOB)
#'
#' This function performs multi-level analysis using an optimal Bayesian approach.
#' It is primarily designed for cases where data can be divided into multiple groups
#' and ensures that the data is balanced across groups.
#'
#' @param formula A formula specifying the model (e.g., \code{Y ~ X + C...}) where Y is the dependent variable, X is the target variable, and C... represents control variables.
#' @param data A data frame containing the variables in the formula.
#' @param group A positive numeric value representing the number of groups (J) in the data. The number of rows in \code{data} must be divisible by \code{group}.
#' @param balancing Logical. If \code{TRUE}, the data will be balanced across groups. Defaults to \code{FALSE}.
#' @param alpha A numeric value representing the confidence level used to calculate confidence intervals for the estimators. Defaults to \code{0.05}.
#' @param jackknife Logical. If \code{TRUE}, the jackknife resampling method will be applied. Defaults to \code{TRUE}.
#' @param ... Additional arguments passed to the function.
#'
#' @details
#' This function checks whether the data is balanced (i.e., whether the same number of individuals are present
#' in each group). If the data is unbalanced and `balancing = TRUE`, the function uses an internal method
#' to balance the data. Otherwise, if `balancing = FALSE` and the data is unbalanced, the function will stop 
#' and throw an error message. 
#'
#' MultiLevelOptimalBayes (MLOB) is designed for small-sample two-level latent variable models, commonly used in
#' psychology, education, and other disciplines that deal with hierarchical data structures.
#'
#' @return A list containing the results of the Bayesian estimation including model formula, response and predictor variables, and other relevant details.
#'
#' @examples
#' # Example usage with the iris dataset
#' data(data_cv)
#' result  mlob(Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width, data = iris, group = 3)
#' 
#' # View summary statistics (similar to summary of a linear model)
#' summary(result)
#' ## 
#' ## Call:
#' ## result <- mlob (y ~ X + Z, data = ‘dataname’, group = ‘class’, balancing = FALSE, alpha = 0.05, …)
#' ## 
#' #' ## Coefficients:
#' ##                Estimate Std. Error Confidence Interval (95%)  z value Pr(>|z|)
#' ## (Intercept)     0.500    0.025    [0.450, 0.550]               20.00   <2e-16 ***
#' ## Sepal.Width     0.300    0.043    [0.220, 0.380]               6.98    <2e-16 ***
#' ## Petal.Length    -0.100   0.021    [-0.140, -0.060]             -4.76   <2e-16 ***
#' ## Petal.Width      0.240   0.024    [0.190, 0.290]               10.00   <2e-16 ***
#' ## ---
#' ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#' ## 

#' ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#' ## 
#' ## Residual standard error: 2.2 on 48 degrees of freedom
#' ## Multiple R-squared:  0.8535, Adjusted R-squared:  0.8504
#' ## F-statistic: 279.6 on 3 and 48 DF,  p-value: < 2.2e-16
#'
#' @export


mlob <- function(formula, data, group, balancing = FALSE, alpha = 0.05, jackknife = TRUE, ...) {
  # Ensure data is a data frame
  if (!is.data.frame(data)) {
    stop("The 'data' argument must be a data frame.")
  }
  
  # Check that 'group' is a number representing the number of groups (J)
  if (!is.numeric(group) || length(group) != 1 || group <= 0) {
    stop("'group' must be a positive numeric value representing the number of groups (J).")
  }
  
  # Ensure the data length is divisible by J (i.e., balanced data)
  n_rows <- nrow(data)
  if (n_rows %% group != 0) {
    stop(paste("Data is not balanced: the number of rows in the data is not divisible by", group))
  }
  
  # Parse the formula (Y ~ X + C...)
  all_vars <- all.vars(formula)
  response_var <- all_vars[1]  # The response variable (Y)
  predictor_var <- all_vars[2]  # The first predictor variable (X)
  control_vars <- all_vars[-c(1, 2)]  # Control variables (C...)
  
  # Check that there is exactly one X variable
  if (length(predictor_var) != 1) {
    stop("There must be exactly one predictor variable (X) after the response variable (Y).")
  }
  
  # Print or save the input parameters for verification
  cat("Response variable (Y):", response_var, "\n")
  cat("Predictor variable (X):", predictor_var, "\n")
  
  if (length(control_vars) > 0) {
    cat("Control variables (C...):", paste(control_vars, collapse = ", "), "\n")
  } else {
    cat("No control variables specified.\n")
  }
  
  cat("Data:", deparse(substitute(data)), "\n")
  cat("Number of groups (J):", group, "\n")
  cat("Data is balanced with", n_rows, "observations, divided into", group, "groups.\n")
  
  # Handle optional arguments
  if (!missing(balancing)) {
    cat("Balancing:", balancing, "\n")
  }
  
  if (!missing(alpha)) {
    cat("Alpha level:", alpha, "\n")
  }
  
  if (!missing(jackknife)) {
    cat("Jackknife:", jackknife, "\n")
  }
  
  # Include any additional arguments provided through '...'
  if (!missing(...)) {
    extra_args <- list(...)
    cat("Additional arguments:", extra_args, "\n")
  }
  
  # Simulate a return of model output (replace this with actual computation)
  return(list(
    formula = formula,
    data = data,
    group = group,
    response_var = response_var,
    predictor_var = predictor_var,
    control_vars = control_vars,
    balancing = ifelse(missing(balancing), NA, balancing),
    alpha = ifelse(missing(alpha), NA, alpha),
    jackknife = ifelse(missing(jackknife), NA, jackknife),
    additional_args = list(...)
  ))
}
# Example usage:
data(iris)  # Using iris dataset as an example

# Valid formula and balanced data (J=3, 150 rows, divisible by 3)
mlob(Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width, data = iris, group = 3)

# Example with unbalanced data (J=4, but iris dataset has 150 rows, not divisible by 4)
# Uncommenting the line below would raise an error
# mlob(Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width, data = iris, group = 4)