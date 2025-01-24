% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_Bay_CV.R
\name{estimate_Bay_CV}
\alias{estimate_Bay_CV}
\title{Estimate Between-Group Effect of 2-Level Latent Variable Model with Control Variables Using Regularized Bayesian Estimator}
\usage{
estimate_Bay_CV(data)
}
\arguments{
\item{data}{A list containing the necessary data for covariance estimation. The list should include:
\describe{
  \item{k}{Number of groups or clusters.}
  \item{n}{Number of observations per group.}
  \item{kc}{Number of control variables.}
  \item{ICC_x}{Intraclass correlation coefficient for `x`.}
  \item{ICC_y}{Intraclass correlation coefficient for `y`.}
  \item{ICC_C}{Intraclass correlation coefficient for control variables.}
  \item{b0}{Parameter for intercept.}
  \item{b_w}{Parameter for within-group variance.}
  \item{b_b}{Parameter for between-group variance.}
  \item{gamma}{Parameter for the effect of control variables.}
  \item{kn}{Number of observations used for gamma estimation.}
  \item{m_x}{Parameter related to `x`.}
  \item{var_x1}{Variance component 1 for `x`.}
  \item{var_x2}{Variance component 2 for `x`.}
  \item{var_e1}{Variance component 1 for error.}
  \item{var_e2}{Variance component 2 for error.}
  \item{cov_mat}{Covariance matrix of `x`.}
  \item{cov_mat_b}{Covariance matrix of between-group variance.}
  \item{x2}{Matrix of `x` values.}
  \item{x}{Matrix of `x` values.}
  \item{e2}{Matrix of error terms.}
  \item{m_C}{Parameter related to control variables.}
  \item{var_C1}{Variance component 1 for control variables.}
  \item{var_C2}{Variance component 2 for control variables.}
  \item{C2}{Matrix of control variables (if multiple).}
  \item{C}{Matrix of control variables (if single).}
  \item{y}{Outcome variable matrix.}
  \item{x2D}{Matrix of `x2` values (if applicable).}
}}
}
\value{
A list containing the following matrices:
\describe{
  \item{Sigma_x}{Covariance matrix of `x`.}
  \item{Sigma_yx}{Covariance matrix between `x` and `y`.}
  \item{Sigma_y}{Covariance matrix of `y`.}
}
}
\description{
This function estimates covariance components and constructs covariance matrices based on input data from either R or MATLAB formats. It supports scenarios with and without control variables and computes various statistics related to the data.
}
\examples{
# Example data
data_example <- list(
  k = 5,
  n = 10,
  kc = 2,
  ICC_x = 0.5,
  ICC_y = 0.5,
  ICC_C = 0.3,
  b0 = 1,
  b_w = 0.5,
  b_b = 0.5,
  gamma = 0.2,
  kn = 50,
  m_x = 1,
  var_x1 = 1,
  var_x2 = 1,
  var_e1 = 0.5,
  var_e2 = 0.5,
  cov_mat = matrix(rnorm(25), 5, 5),
  cov_mat_b = matrix(rnorm(25), 5, 5),
  x2 = matrix(rnorm(50), 10, 5),
  x = matrix(rnorm(50), 10, 5),
  e2 = matrix(rnorm(50), 10, 5),
  m_C = 1,
  var_C1 = 0.5,
  var_C2 = 0.5,
  C2 = matrix(rnorm(50), 10, 5),
  C = matrix(rnorm(50), 10, 5),
  y = matrix(rnorm(50), 10, 5),
  x2D = matrix(rnorm(50), 10, 5)
)

result <- estimate_Bay_CV(data_example)

}
