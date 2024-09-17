#' Estimate between-group effect of 2-level latent variable model with control variables using maximum likelihood (ML) estimator
#' Maximum Likelihood Estimation of Between-Group Effect for 2-Level Latent Variable Model
#'
#' This function estimates the between-group effect in a two-level latent variable model with control variables 
#' using the Maximum Likelihood (ML) estimator. It calculates the between-group variance components for 
#' the predictor variable (`x`) and response variable (`y`), and estimates the between-group effect (`beta_b_ML`). 
#' If control variables (`C`) are provided, the function also computes the adjusted between-group effect 
#' after controlling for those variables (`beta_b_ML_tilde`).
#'
#' The function handles multiple control variables by regressing them on the predictor variable 
#' and incorporating the resulting residuals into the estimation process.
#'
#' @param data A list containing the dataset with the following elements:
#'   \itemize{
#'     \item `n`: Number of individuals per group.
#'     \item `k`: Number of groups.
#'     \item `x`: Predictor variable.
#'     \item `y`: Response variable.
#'     \item `C`: Control variables matrix (optional).
#'     \item `kc`: Number of control variables.
#'   }
#' @return A list containing:
#'   \itemize{
#'     \item `beta_b_ML_CV`: The ML estimate of the between-group effect after controlling for control variables (if applicable).
#'     \item `beta_b_ML_no_CV`: The ML estimate of the between-group effect without controlling for control variables.
#'     \item `gamma`: The estimated coefficients for the control variables.
#'     \item `tau_x2`: Between-group variance of the predictor variable.
#'     \item `sigma_x2`: Within-group variance of the predictor variable.
#'     \item `tau_yx`: Covariance between the predictor and response variables.
#'     \item `sigma_yx`: Residual covariance between the predictor and response variables.
#'   }
#' @examples
#' # Assuming data is a list structured as per the requirements
#' result <- estimate_ML_CV(data)
#' print(result)
#' @importFrom pracma sqrt
#' @export


estimate_ML_CV <- function(data) {
  
  library(pracma)
  
  ##
  
  tog <- 0 ## Matlab load toggle and change signature data - used for debug
  
  if (tog>0) {
    
    library(R.matlab)
    data <- readMat("data_CV_matlab.mat") #change file name
    #print(data)
    
    data1 <- list(
      k=(data$data[[1]])[1],
      n=(data$data[[2]])[1],
      kc=(data$data[[3]])[1],
      ICC_x=(data$data[[4]])[1],
      ICC_y=(data$data[[5]])[1],
      ICC_C =(data$data[[6]])[1],
      b0=(data$data[[7]])[1],
      b_w=(data$data[[8]])[1],
      b_b=(data$data[[9]])[1],
      gamma=(data$data[[10]])[1],
      kn=(data$data[[11]])[1],
      m_x=(data$data[[12]])[1],
      var_x1=(data$data[[13]])[1],
      var_x2=(data$data[[14]])[1],
      var_e1=(data$data[[15]])[1],
      var_e2=(data$data[[16]])[1],
      
      cov_mat =(data$data[[17]]),
      cov_mat_b =(data$data[[18]]),
      
      x2=(data$data[[19]]),
      x=(data$data[[20]]),
      e2=(data$data[[21]]),
      m_C =(data$data[[22]])[1],
      var_C1 = (data$data[[23]])[1],
      var_C2 =(data$data[[24]])[1],
      
      C2=(data$data[[25]]),
      C=(data$data[[26]]),
      y=(data$data[[27]]),
      x2D=(data$data[[28]])
    )
    
    
    data <-data1
    #print(data)
    
  }
  
  
  data$j <- data$k
  # Reshape the data.x into an n x j matrix
  x <- matrix(data$x, nrow = data$n, ncol = data$k)
  
  av_x <- mean(x)
  av_x_j <- colMeans(x)
  
  SSA <- data$n * sum(av_x_j^2) - data$n * data$k * (av_x^2)
  SSD <- sum(data$x * data$x) - data$n * sum(av_x_j * av_x_j)
  
  
  MSA <- SSA / (data$j - 1)
  MSD <- SSD / ((data$n - 1) * data$j)
  
  tau_x2 <- (MSA - MSD) / data$n
  sigma_x2 <- MSD
  
  # Reshape the data.y into an n x j matrix
  y <- matrix(data$y, nrow = data$n, ncol = data$j)
  
  av_y <- mean(y)
  av_y_j <- colMeans(y)
  
  SPA <- data$n * sum(av_y_j * av_x_j) - data$n * data$j * (av_x * av_y)
  SPD <- sum(data$x *data$y) - data$n * sum(av_y_j * av_x_j)
  
  MPA <- SPA / (data$j - 1)
  MPD <- SPD / ((data$n - 1) * data$j)
  
  tau_yx <- (MPA - MPD) / data$n
  sigma_yx <- MPD
  
  # ML estimation of beta_b
  beta_b_ML <- tau_yx / tau_x2
  
  # Multiple control variables part
  beta_b_ML_tilde <- NULL
  data$C = matrix(data$C)
  if (!is.null(data$kc) && data$kc != 0) {
    
    phi <- matrix(0, nrow = 3, ncol = data$kc)
    data$C = matrix(data$C, nrow = data$k*data$n, ncol = data$kc)
    for (i in 1:data$kc) {
      
      if (data$kc>1){
        # Reshape data$C[,i] into a matrix with n rows and J columns
        C_i_C <- matrix(data$C[, i], nrow = data$n, ncol = data$k)
      }else{
        C_i_C <- matrix(data$C, nrow = data$n, ncol = data$k)
      }
      
      av_C <- mean(C_i_C)
      av_C_j <- colMeans(C_i_C)
      
      SPA_C <- data$n * sum(av_C_j * av_x_j) - data$n * data$j* (av_x * av_C)
      SPD_C <- sum(x * data$C[, i]) - data$n * sum(av_C_j * av_x_j)
      
      MPA_C <- SPA_C / (data$j - 1)
      MPD_C <- SPD_C / ((data$n - 1) * data$j)
      
      tau_C_x <- (MPA_C - MPD_C) / data$n
      sigma_C_x <- MPD_C
      
      phi[2, i] <- tau_C_x / tau_x2
      phi[3, i] <- sigma_C_x / sigma_x2
    }
    
    C <- data$C
    x_b <- rep(av_x_j, each = data$n)
    #x_w <- as.vector(t(x)) - x_b 
    x_w <- as.vector(data$x) - x_b
    C_tilde <- C - cbind(x_b, x_w) %*% (phi[2:3, ])
    
    C_tilde_mat <- cbind(1, C_tilde)
    
    gamma_all <- solve(t(C_tilde_mat) %*% C_tilde_mat) %*% t(C_tilde_mat) %*% data$y
    gamma <- gamma_all[-1]
    
    y_tilde <- data$y - C %*% gamma
    y_tilde <- matrix(y_tilde, nrow = data$n, ncol = data$j)
    
    av_y_tilde <- mean(y_tilde)
    av_y_tilde_j <- colMeans(y_tilde)
    
    SPA_tilde <- data$n * sum(av_y_tilde_j * av_x_j) - data$n * data$j * (av_x * av_y_tilde)
    #SPD_tilde <- sum(x * as.vector(t(y_tilde))) - data$n * sum(av_y_tilde_j * av_x_j)
    SPD_tilde <- sum(as.vector(data$x) * as.vector(y_tilde)) - data$n * sum(av_y_tilde_j * av_x_j)
    
    MPA_tilde <- SPA_tilde / (data$j - 1)
    MPD_tilde <- SPD_tilde / ((data$n - 1) * data$j)
    
    tau_y_tilde_x <- (MPA_tilde - MPD_tilde) / data$n
    sigma_y_tilde_x <- MPD_tilde
    
    # ML estimation of beta_b for the model with control variables
    beta_b_ML_tilde <- tau_y_tilde_x / tau_x2
  }
  
  return(list(beta_b_ML_CV = beta_b_ML_tilde, beta_b_ML_no_CV = beta_b_ML, gamma = gamma, tau_x2 = tau_x2, sigma_x2 = sigma_x2, tau_yx = tau_yx, sigma_yx = sigma_yx))
}

# Example usage
#result <- estimate_ML_CV(data_CV)
#print(result)
