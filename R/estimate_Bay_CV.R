#' Estimate Between-Group Effect of 2-Level Latent Variable Model with Control Variables Using Regularized Bayesian Estimator
#'
#' This function estimates covariance components and constructs covariance matrices based on input data from either R or MATLAB formats. It supports scenarios with and without control variables and computes various statistics related to the data.
#'
#' @param data A list containing the necessary data for covariance estimation. The list should include:
#' \describe{
#'   \item{k}{Number of groups or clusters.}
#'   \item{n}{Number of observations per group.}
#'   \item{kc}{Number of control variables.}
#'   \item{ICC_x}{Intraclass correlation coefficient for `x`.}
#'   \item{ICC_y}{Intraclass correlation coefficient for `y`.}
#'   \item{ICC_C}{Intraclass correlation coefficient for control variables.}
#'   \item{b0}{Parameter for intercept.}
#'   \item{b_w}{Parameter for within-group variance.}
#'   \item{b_b}{Parameter for between-group variance.}
#'   \item{gamma}{Parameter for the effect of control variables.}
#'   \item{kn}{Number of observations used for gamma estimation.}
#'   \item{m_x}{Parameter related to `x`.}
#'   \item{var_x1}{Variance component 1 for `x`.}
#'   \item{var_x2}{Variance component 2 for `x`.}
#'   \item{var_e1}{Variance component 1 for error.}
#'   \item{var_e2}{Variance component 2 for error.}
#'   \item{cov_mat}{Covariance matrix of `x`.}
#'   \item{cov_mat_b}{Covariance matrix of between-group variance.}
#'   \item{x2}{Matrix of `x` values.}
#'   \item{x}{Matrix of `x` values.}
#'   \item{e2}{Matrix of error terms.}
#'   \item{m_C}{Parameter related to control variables.}
#'   \item{var_C1}{Variance component 1 for control variables.}
#'   \item{var_C2}{Variance component 2 for control variables.}
#'   \item{C2}{Matrix of control variables (if multiple).}
#'   \item{C}{Matrix of control variables (if single).}
#'   \item{y}{Outcome variable matrix.}
#'   \item{x2D}{Matrix of `x2` values (if applicable).}
#' }
#'
#' @return A list containing the following matrices:
#' \describe{
#'   \item{Sigma_x}{Covariance matrix of `x`.}
#'   \item{Sigma_yx}{Covariance matrix between `x` and `y`.}
#'   \item{Sigma_y}{Covariance matrix of `y`.}
#' }
#'
#' @examples
#' # Example data
#' data_example <- list(
#'   k = 5,
#'   n = 10,
#'   kc = 2,
#'   ICC_x = 0.5,
#'   ICC_y = 0.5,
#'   ICC_C = 0.3,
#'   b0 = 1,
#'   b_w = 0.5,
#'   b_b = 0.5,
#'   gamma = 0.2,
#'   kn = 50,
#'   m_x = 1,
#'   var_x1 = 1,
#'   var_x2 = 1,
#'   var_e1 = 0.5,
#'   var_e2 = 0.5,
#'   cov_mat = matrix(rnorm(25), 5, 5),
#'   cov_mat_b = matrix(rnorm(25), 5, 5),
#'   x2 = matrix(rnorm(50), 10, 5),
#'   x = matrix(rnorm(50), 10, 5),
#'   e2 = matrix(rnorm(50), 10, 5),
#'   m_C = 1,
#'   var_C1 = 0.5,
#'   var_C2 = 0.5,
#'   C2 = matrix(rnorm(50), 10, 5),
#'   C = matrix(rnorm(50), 10, 5),
#'   y = matrix(rnorm(50), 10, 5),
#'   x2D = matrix(rnorm(50), 10, 5)
#' )
#'
#' result <- estimate_Bay_CV(data_example)
#'
#' @export


estimate_Bay_CV <- function(data) {
  
  library(pracma)
  
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
  
  
  J <- data$k
  n <- data$n
  
  a <- rep(0, n * J + J + 1)
  
  for (i in 1:(n * J)) {
    a[i] <- -1 / (n * (n - 1) * J)
  }
  
  for (i in (n * J + 1):(n * J + J)) {
    a[i] <- (n * J - 1) / ((n - 1) * (J - 1) * J)
  }
  
  a[n * J + J + 1] <- -J / (J - 1)
  
  A <- diag(a)
  
  # B1 = [A, zeros(n*J+J+1); zeros(n*J+J+1,2*(n*J+J+1))];
  
  # B2 = [zeros(n*J+J+1), A/2; A/2, zeros(n*J+J+1)];
  
  # calculate covariance matrices:
    
  #cov(x):
    
  x <- matrix(data$x, n, J)
  
  av_x <- mean(x)
  av_x_j <- colMeans(x)
  
  SSA <- n * sum(av_x_j^2) - n * J * (av_x^2)
  SSD <- sum(data$x^2) - n * sum(av_x_j^2)
  
  MSA <- SSA / (J - 1)
  MSD <- SSD / ((n - 1) * J)
  
  tau_x2 <- (MSA - MSD) / n
  sigma_x2 <- MSD
  
  #Multiple control variables part
  # checks the existence of control variables and their number is not 0
  
  #Find 3 phi's for every of kc control variables from equation 
  # C = phi1 + phi2*X_b + phi3*X_w + eps. Our aim is phi2 and phi3  
  
  if (!is.null(data$kc) && data$kc != 0) { 
    phi <- matrix(0, 3, data$kc)
    
    C <- list()
    
    for (i in 1:data$kc) {
      C[[i]] <- list()
      
      # 2 cases for number of control variables
      if (data$kc>1){
        # Reshape data$C[,i] into a matrix with n rows and J columns
        C[[i]]$C <- matrix(data$C[, i], nrow = n, ncol = J)
      }else{
        C[[i]]$C <- matrix(data$C, nrow = n, ncol = J)
      }
      
      # Compute the average of the entire C matrix
      C[[i]]$av_C <- mean(C[[i]]$C)
      
      # Compute the average of each column in the C matrix
      C[[i]]$av_C_j <- colMeans(C[[i]]$C)
      
      # Compute SPA_C and SPD_C
      C[[i]]$SPA_C <- n * sum(C[[i]]$av_C_j * av_x_j) - n * J * (av_x * C[[i]]$av_C)
      # 2 cases for number of control variables
      if (data$kc>1){
        C[[i]]$SPD_C <- sum(data$x * data$C[, i]) - n * sum(C[[i]]$av_C_j * av_x_j)
      }else{
        C[[i]]$SPD_C <- sum(data$x * data$C) - n * sum(C[[i]]$av_C_j * av_x_j) 
      }
      # Compute MPA_C and MPD_C
      C[[i]]$MPA_C <- C[[i]]$SPA_C / (J - 1)
      C[[i]]$MPD_C <- C[[i]]$SPD_C / ((n - 1) * J)
      
      # Compute tau_C_x and sigma_C_x
      C[[i]]$tau_C_x <- (C[[i]]$MPA_C - C[[i]]$MPD_C) / n
      C[[i]]$sigma_C_x <- C[[i]]$MPD_C
      
      # Initialize phi as a vector of zeros
      C[[i]]$phi <- numeric(3)
      
      # Compute phi values
      C[[i]]$phi[2] <- C[[i]]$tau_C_x / tau_x2
      C[[i]]$phi[3] <- C[[i]]$sigma_C_x / sigma_x2
      
      # Store the phi vector in the phi matrix
      phi[, i] <- C[[i]]$phi
    }
    
    C <- data$C
    x_b <- matrix(rep(av_x_j, each = n), nrow = n * J, ncol = 1)
    #x_b <- matrix(rep(av_x_j, each = n), nrow = n, byrow = TRUE)
    x_w <- as.vector(matrix(data$x, nrow = n * J, ncol = 1)) - x_b
    C_tilde <- C - cbind(x_b, x_w) %*% phi[2:nrow(phi), ]
    C_tilde_mat <- cbind(rep(1, data$kn), C_tilde)
    
    #Estimate gamma from equation y = gamma1 + gamma*C_tilde + eps. Our
    #interest is vestor of parameters gamma without intersection gamma1
    gamma_all <- solve(t(C_tilde_mat) %*% C_tilde_mat) %*% t(C_tilde_mat) %*% data$y
    
    SE_gamma <- numeric(data$kc)
    
    for (i in 1:data$kc) {
      var_res_gamma <- sum((data$y - C_tilde_mat[, i + 1] * gamma_all[i + 1])^2) / (data$kn - 2)
      SE_gamma[i] <- sqrt(var_res_gamma / ((data$kn - 1) * var(C_tilde_mat[, i + 1])))
    }
    
    gamma <- gamma_all[-1]
    # Find y_tilde - the difference between real and estimated 
    # y (without intersection)
    if (data$kc>1){
      y_tilde <- data$y - C %*% gamma
    }else{
      y_tilde <- data$y - C * gamma
    }
    data$y <- y_tilde
  }
  #cov(x,y)
  y <- matrix(data$y, n, J)
  
  av_y <- mean(y)
  av_y_j <- colMeans(y)
  
  SPA <- n * sum(av_y_j * av_x_j) - n * J * (av_x * av_y)
  SPD <- sum(data$x * data$y) - n * sum(av_y_j * av_x_j)
  
  MPA <- SPA / (J - 1)
  MPD <- SPD / ((n - 1) * J)
  
  tau_yx <- (MPA - MPD) / n
  sigma_yx <- MPD
  
  #cov(y):
  SOA <- n * sum(av_y_j^2) - n * J * (av_y^2)
  SOD <- sum(data$y^2) - n * sum(av_y_j^2)
  
  MOA <- SOA / (J - 1)
  MOD <- SOD / ((n - 1) * J)
  
  tau_y2 <- (MOA - MOD) / n
  sigma_y2 <- MOD

  #create Sigma_x matrix of covariances of x_ij, mean(x_j), mean(x)
  
  # create Sigma blockvise with 9 blocks
  
  # block 1: cov(x_ij,x_ij): 
  a1 <- diag(sigma_y2 * rep(1, n))
  #a1 <- diag(sigma_x2, n)
  
  a2 <- matrix(tau_x2, n, n)
  
  a <- a1 + a2
  block1 <- kronecker(diag(J), a)
  # Block 2: cov(y_ij, mean(y_j))
  
  block2 <- kronecker(diag(J), (tau_x2 + sigma_x2 / n) * matrix(1, n, 1))
  #block2 <- kronecker(diag(J), (tau_y2 + sigma_y2 / n) * rep(1, n))
  
  block3 <- matrix((tau_x2 / J + sigma_x2 / (n * J)), n * J, 1)
  #block3 <- (tau_y2 / J + sigma_y2 / (n * J)) * rep(1, n * J)
  block4 <- t(block2)
  
  block5 <- (tau_x2 + sigma_x2 / n) * diag(J)
  
  block6 <- (tau_x2 / J + sigma_x2 / (n * J)) * matrix(1, J, 1)
  #block6 <- (tau_y2 / J + sigma_y2 / (n * J)) * rep(1, J
  block7 <- t(block3)
  
  block8 <- t(block6)
  
  block9 <- tau_x2 / J + sigma_x2 / (n * J)
  
  Sigma_x <- rbind(cbind(block1, block2, block3), cbind(block4, block5, block6), cbind(block7, block8, block9))
  
  a1 <- diag(sigma_yx, n)
  a2 <- matrix(tau_yx, n, n)
  
  a <- a1 + a2
  block1 <- kronecker(diag(J), a)
  
  block2 <- kronecker(diag(J), (tau_yx + sigma_yx / n) * matrix(1, n, 1))
  
  block3 <- matrix((tau_yx / J + sigma_yx / (n * J)), n * J, 1)
  
  block4 <- t(block2)
  
  block5 <- (tau_yx + sigma_yx / n) * diag(J)
  
  block6 <- (tau_yx / J + sigma_yx / (n * J)) * matrix(1, J, 1)
  
  block7 <- t(block3)
  
  block8 <- t(block6)
  
  block9 <- tau_yx / J + sigma_yx / (n * J)
  
  Sigma_yx <- rbind(cbind(block1, block2, block3), cbind(block4, block5, block6), cbind(block7, block8, block9))
  
  a1 <- diag(sigma_y2, n)
  a2 <- matrix(tau_y2, n, n)
  
  a <- a1 + a2
  block1 <- kronecker(diag(J), a)
  
  block2 <- kronecker(diag(J), (tau_y2 + sigma_y2 / n) * matrix(1, n, 1))
  
  block3 <- matrix((tau_y2 / J + sigma_y2 / (n * J)), n * J, 1)
  
  block4 <- t(block2)
  
  block5 <- (tau_y2 + sigma_y2 / n) * diag(J)
  
  block6 <- (tau_y2 / J + sigma_y2 / (n * J)) * matrix(1, J, 1)
  
  block7 <- t(block3)
  
  block8 <- t(block6)
  
  block9 <- tau_y2 / J + sigma_y2 / (n * J)
  
  Sigma_y <- rbind(cbind(block1, block2, block3), cbind(block4, block5, block6), cbind(block7, block8, block9))
  
  #spectral decomposition of Sigma_x
  
  #compute eigenvalues matrix D_x of Sigma_x
  
  # Initialize D_x
  D_x <- rep(0, n * J + J + 1)
  
  # Set values for specific indices
  D_x[(J + 2):(n * J + 1)] <- sigma_x2
  D_x[(n * J + 2):(n * J + J)] <- (n + 1) * (tau_x2 + sigma_x2 / n)
  D_x[n * J + J + 1] <- (n * J + J + 1) * (tau_x2 + sigma_x2 / n) / J
  
  # Create diagonal matrix D_x
  D_x <- diag(D_x)
  
  # Initialize V_yx matrix
  V_yx <- matrix(0, n * J + J + 1, n * J + J + 1)
  
  # lambda_i = 0 ((J+1) pieces)
  for (i in 1:J) {
    V_yx[((n * (i - 1) + 1):(n * i)), i] <- -1 / sqrt(n * (n + 1))
    V_yx[(n * J + i), i] <- sqrt(n / (n + 1))
  }
  
  V_yx[1:(n * J + J), J + 1] <- -1 / sqrt((n * J + J) * (n * J + J + 1))
  V_yx[(n * J + J + 1), J + 1] <- sqrt((n * J + J) / (n * J + J + 1))
  
  # lambda_i = sigma_yx ((n-1)J pieces)
  AU <- matrix(0, n, n - 1)
  for (i in 1:(n - 1)) {
    AU[1:i, i] <- -1 / sqrt(i * (i + 1))
    AU[(i + 1), i] <- i / sqrt(i * (i + 1))
  }
  
  for (j in 1:J) {
    V_yx[(n * (j - 1) + 1):(n * j), (J + 2 + (n - 1) * (j - 1)):(J + 1 + (n - 1) * j)] <- AU
  }
  
  # lambda_i = (n + 1)(tau_yx + sigma_yx / n) ((J - 1) pieces)
  for (i in 1:(J - 1)) {
    V_yx[1:(i * n), (n * J + 1 + i)] <- -1 / sqrt(i * (i + 1) * (n + 1))
    V_yx[(n * J + 1):(n * J + i), (n * J + 1 + i)] <- -1 / sqrt(i * (i + 1) * (n + 1))
    V_yx[(i * n + 1):((i + 1) * n), (n * J + 1 + i)] <- i / sqrt(i * (i + 1) * (n + 1))
    V_yx[(n * J + i + 1), (n * J + 1 + i)] <- i / sqrt(i * (i + 1) * (n + 1))
  }
  
  # lambda_nJ + J + 1 = (n * J + J + 1)(tau_yx + sigma_yx / n) / J (1 piece)
  V_yx[1:(n * J + J + 1), (n * J + J + 1)] <- 1 / sqrt(n * J + J + 1)
  
  # Spectral decomposition of Sigma_y
  
  # Initialize D_y
  D_y <- rep(0, n * J + J + 1)
  
  # Set values for specific indices
  D_y[(J + 2):(n * J + 1)] <- sigma_y2
  D_y[(n * J + 2):(n * J + J)] <- (n + 1) * (tau_y2 + sigma_y2 / n)
  D_y[n * J + J + 1] <- (n * J + J + 1) * (tau_y2 + sigma_y2 / n) / J
  
  # Create diagonal matrix D_y
  D_y <- diag(D_y)
  
  # Initialize V_y matrix
  V_y <- matrix(0, n * J + J + 1, n * J + J + 1)
  
  # lambda_i = 0 ((J + 1) pieces)
  for (i in 1:J) {
    V_y[((n * (i - 1) + 1):(n * i)), i] <- -1 / sqrt(n * (n + 1))
    V_y[(n * J + i), i] <- sqrt(n / (n + 1))
  }
  
  V_y[1:(n * J + J), J + 1] <- -1 / sqrt((n * J + J) * (n * J + J + 1))
  V_y[(n * J + J + 1), J + 1] <- sqrt((n * J + J) / (n * J + J + 1))
  
  # lambda_i = sigma_yx ((n - 1)J pieces)
  for (i in 1:(n - 1)) {
    AU[1:i, i] <- -1 / sqrt(i * (i + 1))
    AU[(i + 1), i] <- i / sqrt(i * (i + 1))
  }
  
  for (j in 1:J) {
    V_y[(n * (j - 1) + 1):(n * j), (J + 2 + (n - 1) * (j - 1)):(J + 1 + (n - 1) * j)] <- AU
  }
  V_x <- matrix(0, nrow = n*sum(J)+sum(J)+1, ncol = n*sum(J)+sum(J)+1)
  
  # lambda_i = (n + 1)(tau_yx + sigma_yx / n) ((J - 1) pieces)
  for (i in 1:(J - 1)) {
    V_x[1:(i * n), (n * J + 1 + i)] <- -1 / sqrt(i * (i + 1) * (n + 1))
    V_x[(n * J + 1):(n * J + i), (n * J + 1 + i)] <- -1 / sqrt(i * (i + 1) * (n + 1))
    V_x[(i * n + 1):((i + 1) * n), (n * J + 1 + i)] <- i / sqrt(i * (i + 1) * (n + 1))
    V_x[(n * J + i + 1), (n * J + 1 + i)] <- i / sqrt(i * (i + 1) * (n + 1))
  }
  
  # lambda_nJ + J + 1 = (n * J + J + 1)(tau_yx + sigma_yx / n) / J (1 piece)
  V_x[1:(n * J + J + 1), (n * J + J + 1)] <- 1 / sqrt(n * J + J + 1)
  # check the correctness of spectral decomposition:
  
  # diff_x = sum(sum((Sigma_x - V_x*D_x*V_x').^2));
  
  # disp(['difference between covariance matrix for x and its decomposition is ', num2str(diff_x)]);
  
  
  # Spectral decomposition of Sigma_yx
  
  # compute eigenvalues matrix D_yx of Sigma_yx
  D_yx <- matrix(0, n * sum(J) + sum(J) + 1, 1)
  
  for (i in ((J+2 ): ((n*J)+1))){
    D_yx[i,] <- sigma_yx
  }
  
  for (i in (((n*J)+2) : ((n*J)+J))){
    D_yx[i,] <- (n+1)*(tau_yx + sigma_yx/n)
  }
  
  D_yx[(n*J)+J+1,] <- ((n*J)+J+1)*(tau_yx + sigma_yx/n)/J
  
  D_yx <- diag(D_yx[,1])
  
  # compute eigenvectors of Sigma_yx
  
  # create matrix of eigenvectors V_yx for covariance matrix Sigma_yx
  
  V_yx <- matrix(0, nrow = n*sum(J)+sum(J)+1, ncol = n*sum(J)+sum(J)+1)
  
  
  # lambda_i = 0 ((J+1) pieces)
  
  for (i in 1:(length(seq_len(J)))){
    V_yx[((n*(i-1))+1):(n*i),i] <- -1/sqrt(n*(n+1))
    V_yx[(n*J)+i,i] <- sqrt(n/(n+1))
  }
  
  V_yx[1:((n*J)+J),J+1] <- -1/sqrt((n*J+J)*(n*J+J+1))
  V_yx[(n*J)+J+1,J+1] <- sqrt((n*J+J)/(n*J+J+1))
  
  # lambda_i = sigma_yx ((n-1)J pieces)
  # auxiliary matrix AU
  
  AU <- matrix(0, nrow = n, ncol = n - 1)
  
  for (i in 1:(n-1)){
    AU[1:i,i] <- -1/sqrt(i*(i+1))
    AU[i+1,i] <- i/sqrt(i*(i+1))
  }
  
  for (j in 1:length(seq_len(J))){
    V_y[(n*(j-1)+1):(n*j),(J+2+((n-1)*(j-1))):(J+1+((n-1)*j))] <- AU
  }
  # lambda_i = (n+1)(tau_yx + sigma_yx/n) ((J-1) pieces)
  
  for (i in 1:(J-1)){
    V_y[1:(i*n),(n*J)+1+i] <- -1/sqrt(i*(i+1)*(n+1))
    V_y[((n*J)+1):((n*J)+i),(n*J)+1+i] <- -1/sqrt(i*(i+1)*(n+1))
    V_y[(i*n+1):((i+1)*n),(n*J)+1+i] <- i/sqrt(i*(i+1)*(n+1))
    V_y[(n*J)+i+1,(n*J)+1+i] <- i/sqrt(i*(i+1)*(n+1))
  }
  
  # lambda_nJ+J+1 = (nJ+J+1)(tau_yx + sigma_yx/n)/J (1 piece)
  
  V_yx[1:((n*J)+J+1),(n*J)+J+1] <- 1/sqrt(n*J+J+1)
  
  # check the correctness of spectral decomposition:
  
  # diff_yx = sum(sum((Sigma_yx - V_yx*D_yx*V_yx').^2));
  
  # disp(['difference between covariance matrix for x and y and its decomposition is ', num2str(diff_yx)]);
  
  
  # Spectral decomposition of Sigma_y
  
  # compute eigenvalues matrix D_y of Sigma_y
  
  D_y <- matrix(0, nrow = n*sum(J)+sum(J)+1, ncol = 1)
  
  for (i in ((J+2) : ((n*J)+1))){
    D_y[i,] <- sigma_y2
  }
  
  for (i in ((n*J+2) : ((n*J)+J))){
    D_y[i,] <- (n+1)*(tau_y2 + sigma_y2/n)
  }
  
  D_y[(n*J)+J+1,] <- ((n*J)+J+1)*(tau_y2 + sigma_y2/n)/J
  
  D_y <- diag(D_y[,1])
  
  # compute eigenvectors of Sigma_y
  
  # create matrix of eigenvectors V_y for covariance matrix Sigma_y
  
  V_y <- matrix(0, nrow = n*sum(J)+sum(J)+1, ncol = n*sum(J)+sum(J)+1)
  
  # lambda_i = 0 ((J+1) pieces)
  
  for (i in 1:(length(seq_len(J)))){
    V_y[((n*(i-1))+1):(n*i),i] <- -1/sqrt(n*(n+1))
    V_y[(n*J)+i,i] <- sqrt(n/(n+1))
  }
  
  V_y[1:((n*J)+J),J+1] <- -1/sqrt((n*J+J)*(n*J+J+1))
  V_y[(n*J)+J+1,J+1] <- sqrt((n*J+J)/(n*J+J+1))
  
  # lambda_i = sigma_yx ((n-1)J pieces)
  
  # auxiliary matrix AU
  
  AU <- matrix(0, nrow = n, ncol = n - 1)
  for (i in 1:(n-1)){
    AU[1:i,i] <- -1/sqrt(i*(i+1))
    AU[i+1,i] <- i/sqrt(i*(i+1))
  }
  
  for (j in 1:J){
    V_y[(n*(j-1)+1):(n*j),(J+2+((n-1)*(j-1))):(J+1+((n-1)*j))] <- AU
  }
  
  # lambda_i = (n+1)(tau_yx + sigma_yx/n) ((J-1) pieces)
  
  for (i in 1:(J-1)){
    V_y[1:(i*n),(n*J)+1+i] <- -1/sqrt(i*(i+1)*(n+1))
    V_y[((n*J)+1):((n*J)+i),(n*J)+1+i] <- -1/sqrt(i*(i+1)*(n+1))
    V_y[(i*n+1):((i+1)*n),(n*J)+1+i] <- i/sqrt(i*(i+1)*(n+1))
    V_y[(n*J)+i+1,(n*J)+1+i] <- i/sqrt(i*(i+1)*(n+1))
  }
  
  # lambda_nJ+J+1 = (nJ+J+1)(tau_yx + sigma_yx/n)/J (1 piece)
  
  V_y[1:((n*J)+J+1),(n*J)+J+1] <- 1/sqrt(n*J+J+1)
  
  # check the correctness of spectral decomposition:
  
  # diff_y = sum(sum((Sigma_y - V_y*D_y*V_y').^2));
  
  # disp(['difference between covariance matrix for y and its decomposition is ', num2str(diff_y)]);
  
  
  # calculate the eigenvalues of transformed sum Z_x'*A*Z_x =
  # W_x'*(Sigma_x^0.5)*A*(Sigma_x^0.5)*W_x = {use spectral decomposition and
  # substitute sqrt(D) with S=sqrt(D), because D is diagonal} =
  # = W_x*V_x*S_x*V_x'*A*V_x*S_x*V_x'*W_x = W_x*C_x*W_x
  # C_x = V_x*S_x*V_x'*A*V_x*S_x*V_x' with internal part
  # L1 = S_x*V_x'*A*V_x*S_x  - diagonal -> gives eigenvalues of C_x, which we use
  # to calculate MSE:
  
  S_x <- sqrt(D_x)
  t_v_x <-t(V_x)
  L1 <- S_x %*% t_v_x %*% A %*% V_x %*% S_x 
  
  L1 <- diag(L1)
  
  L1 <- Re(L1)# could be complex insignificant tails, because of computer precision
  
  for (i in 1:(n*J+J+1)){ # delete computer precision tails
    if (abs(L1[i])<0.0000001){
      L1[i] <- 0
    }
  }
  
  # calculate the eigenvalues of transformed sum Z_yx'*A*Z_yx =
  # W_x'*(Sigma_x^0.5)*A*(Sigma_x^0.5)*W_y = {use spectral decomposition and
  # substitute sqrt(D) with S=sqrt(D), because D is diagonal} =
  # {W_x, W_y ~N(0,I)} = W_x'*V_x*S_x*V_x'*A*V_x*S_y*V_y'*W_y =
  # {H_x = V_x'*W_x ~ N(0,I), H_y = V_y'*W_y ~ N(0,I)} =
  # = H_x'*S_x*V_x'*A*V_x*S_y*H_y = H_x'*Q*H_y = {H = [H_x, H_y]} = H'*Q1*H1
  # = H' * EH^(-0.5) * EH^(0.5) * Q1 * EH^(0.5) * EH^(-0.5) *H =
  # = {EH^(-0.5) * H = H1 ~ N(0,I)} = H1' *EH^(0.5) * Q1 * EH^(0.5) * H1 =
  # {use spectral decomposition and substitute: EH = V_H*D_H*V_H' ->
  # EH^(0.5) = V_H*S_H*V_H' with with S_H=sqrt(D_H), because D_H is diagonal}
  # H1' * V_H * S_H * V_H' * Q1 * V_H * S_H * V_H' * H1 =
  # {H1_tilde = V_H' * H1 ~ N(0,I)} =
  # = H1_tilde' * S_H * V_H' * Q1 * V_H * S_H * H1_tilde, with internal part
  # L2 = S_H * V_H' * Q1 * V_H * S_H  - diagonal -> gives eigenvalues of L2,
  # which we use to calculate MSE:
  
  S_y <- sqrt(D_y)
  
  Q <- S_x %*% t_v_x %*% A %*% V_y %*% S_y ## wrong
  
  #
  upper_left <- matrix(0, nrow(Q), ncol(Q))
  lower_right <- matrix(0, nrow(Q), ncol(Q))
  
  upper_right <- Q / 2
  lower_left <- Q / 2
  
  Q1 <- rbind(cbind(upper_left, upper_right), cbind(lower_left, lower_right))
  
  #
  
  
  #
  upper_left <- diag(nrow(D_yx))
  lower_right <- diag(nrow(D_yx))
  
  pinv_S_x_S_y_D_yx <- pinv(S_x) %*% pinv(S_y) %*% D_yx
  
  
  upper_right <- pinv_S_x_S_y_D_yx
  lower_left <- pinv_S_x_S_y_D_yx
  
  EH <- rbind(cbind(upper_left, upper_right), cbind(lower_left, lower_right))
  #
  
  
  D1 <- EH[(1:(n*J+J+1)), ((n*J+J+2):(2*(n*J+J+1)))]
  D_H <- matrix(0, nrow = nrow(EH), ncol = ncol(EH))
  for (i in 1:(n*J+J+1)){
    D_H[2*i-1,2*i-1] <- 1 + D1[i,i]
    D_H[2*i,2*i] <- 1 - D1[i,i]
  }
  S_H <- ifelse(D_H>0,sqrt(D_H),0)
  rm(D1)
  
  ##-- 
  
  
  V_H <- matrix(0, nrow = nrow(D_H), ncol = ncol(D_H))
  
  for (i in 1:(n*J+J+1)){
    V_H[i,(2*i)-1] <- 1/sqrt(2)
    V_H[i, 2*i] <- 1/sqrt(2)
    V_H[(n*J)+J+1 + i, (2*i)-1] <- 1/sqrt(2)
    V_H[(n*J)+J+1 + i, 2*i] <- -1/sqrt(2)
  }
  
  L2 <- S_H %*% t(V_H) %*% Q1 %*% V_H %*% S_H
  
  L2 <- diag(L2)
  L2 <- Re(L2)# could be complex insignificant tails, because of computer precision
  
  for (i in 1:(2*(n*J+J+1))){# delete computer precision tails
    if (abs(L2[i])<0.0000001){
      L2[i] <- 0
    }
  }
  
  
  
  
  # Use L1 and L2 and ML estimator of beta_b as coefficients for computing optimal MSE
  
  beta_b_ML <- tau_yx / tau_x2
  
  K_sum1 <- sum(L1)^2/(2*(t(L1)%*%L1))# unused
  
  T_sum1 <- (t(L1)%*%L1)/sum(L1)# unused
  
  K_sum2 <- sum(L2)^2/(2*(t(L2)%*%L2))
  
  T_sum2 <- (t(L2)%*%L2)/sum(L2)
  
  beta_b <- data$b2# real value of between parameter ##check if b_b is same as b2
  
  
  
  # MSE of ML estimated parameter beta_b_ML
  
  MSE_add <- (((K_sum2 * T_sum2)/ (T_sum1*(K_sum1-1))) - beta_b)^2 +
    ((K_sum2 * T_sum2^2 * (K_sum1 + K_sum2 - 1)) / (T_sum1^2 * (K_sum1-1)^2 * (K_sum1-2)))
  
  # MSE_add = (K_sum1^2 * T_sum1/ (K_sum2 * T_sum2*(K_sum1-1)) - beta_b)^2 +...
  #    K_sum1^4 * T_sum1^2 * (K_sum1+K_sum2-1)/ (K_sum2^3 * T_sum2^2 * (K_sum1-1)^2 * (K_sum1-2)); # original
  
  
  # optimize MSE with beta_b_ML
  
  
  Tau02 <- kronecker(    matrix(1,  nrow = length(seq(0, 1, by = 0.01)), ncol = 1)   , t( matrix(seq(0.05, 10, by = 0.05), nrow = 1) )  )# restrict search interval of tau02 to 10
  
  W <-kronecker(t(matrix(seq(0, 1, by = 0.01), nrow=1)), matrix(1,nrow = length(seq(0.5, by = 0.05)), ncol=1))
  ## recheck???  
  
  MSE_ML <- matrix(0, nrow = length(W), ncol = 1)
  
  
  
  # use grid search to find w and tau02 that give smallest MSE
  for (i in 1:length(W)){
    MSE_ML[i] <- ((K_sum2 * T_sum2^2 * (K_sum2+1) * ((1-W[i])*Tau02[i] + W[i]*sum(L1)))/
                    ((((1-W[i])*Tau02[i] + W[i]*sum(L1))^2 - 2*W[i]^2*(t(L1)%*%L1)) *
                       (((1-W[i])*Tau02[i] + W[i]*sum(L1))^2 - 4*W[i]^2*(t(L1)%*%L1))) - 
                    (2*beta_b_ML *K_sum2 * T_sum2 * ((1-W[i])*Tau02[i] + W[i]*sum(L1)))/
                    ((((1-W[i])*Tau02[i] + W[i]*sum(L1))^2) - (2*W[i]^2*(t(L1)%*%L1)))+ beta_b_ML)
    
  }
  
  # optimize MSE with beta_b
  
  MSE <- matrix(0, nrow = length(W), ncol = 1)
  
  
  # use grid search to find w and tau02 that give smallest MSE
  for (i in 1:length(W)){
    MSE[i] <- ((K_sum2 * T_sum2^2 * (K_sum2+1) * ((1-W[i])*Tau02[i] + W[i]*sum(L1)))/
                    ((((1-W[i])*Tau02[i] + W[i]*sum(L1))^2 - 2*W[i]^2*(t(L1)%*%L1)) *
                       (((1-W[i])*Tau02[i] + W[i]*sum(L1))^2 - 4*W[i]^2*(t(L1)%*%L1))) - 
                    (2*data$b_b *K_sum2 * T_sum2 * ((1-W[i])*Tau02[i] + W[i]*sum(L1)))/
                    ((((1-W[i])*Tau02[i] + W[i]*sum(L1))^2) - (2*W[i]^2*(t(L1)%*%L1)))+data$b_b)
    
    
   
  }
  
  
  # use MSE to find beta_b_Bay with smallest MSE:
  
  #
  ind_MSEbetabBay <- which.min(MSE)
  MSE_beta_b_Bay <- MSE[ind_MSEbetabBay]
  #
  
  
  w_opt <- W[ind_MSEbetabBay]
  tau02_opt <- Tau02[ind_MSEbetabBay]
  
  beta_b_Bay <- tau_yx/((1-w_opt)*tau02_opt + w_opt*tau_x2)
  
  
  # use MSE_ML to find beta_b_Bay_ML with smallest MSE:
  
  #
  ind_MSEbetabBayML <- which.min(MSE_ML)
  MSE_beta_b_Bay_ML <- MSE_ML[ind_MSEbetabBayML]
  #
  
  w_opt_ML <- W[ind_MSEbetabBayML]
  tau02_opt_ML <- Tau02[ind_MSEbetabBayML]
  
  beta_b_Bay_ML <- tau_yx/((1-w_opt_ML)*tau02_opt_ML + w_opt_ML*tau_x2)
  
  
  
  # Compute standard errors of beta_b_Bay, beta_b_Bay_ML, beta_b_ML using their distributions
  
  if (w_opt != 0) {
    K_B_Bay <- (w_opt * T_sum1 * K_sum1 + (1 - w_opt) * tau02_opt)^2 / (w_opt^2 * T_sum1^2 * K_sum1)
    
    T_B_Bay <- (w_opt^2 * T_sum1^2 * K_sum1) / (w_opt * T_sum1 * K_sum1 + (1 - w_opt) * tau02_opt)
  } else {
    K_B_Bay <- 4  # to avoid negative values in SE_beta_Bay
    
    T_B_Bay <- 0.25
  }
  
  if (w_opt_ML != 0) {
    K_B_Bay_ML <- (w_opt_ML * T_sum1 * K_sum1 + (1 - w_opt_ML) * tau02_opt_ML)^2 / (w_opt_ML^2 * T_sum1^2 * K_sum1)
    
    T_B_Bay_ML <- (w_opt_ML^2 * T_sum1^2 * K_sum1) / (w_opt_ML * T_sum1 * K_sum1 + (1 - w_opt_ML) * tau02_opt_ML)
  } else {
    K_B_Bay_ML <- 4  # to avoid negative values in SE_beta_Bay_ML
    
    T_B_Bay_ML <- 0.25
  }
    
  
  
  SE_beta_ML <- abs((T_sum2 / (T_sum1 * (K_sum1 - 1))) * sqrt(abs((K_sum2 * (K_sum1 + K_sum2 - 1)) / (K_sum1 - 2))))
  
  SE_beta_Bay <- abs((T_sum2 / (T_B_Bay * (K_B_Bay - 1))) * sqrt(abs((K_sum2 * (K_B_Bay + K_sum2 - 1)) / (K_B_Bay - 2))))
  
  SE_beta_Bay_ML <- abs((T_sum2 / (T_B_Bay_ML * (K_B_Bay_ML - 1))) * sqrt(abs((K_sum2 * (K_B_Bay_ML + K_sum2 - 1)) / (K_B_Bay_ML - 2))))
  
  
  #Bay <- 0.1
  Bay <- list()
  Bay$A <- A
  #Bay$B2 = B2;
  Bay$Sigma_x <- Sigma_x
  Bay$Sigma_yx <- Sigma_yx
  Bay$Sigma_y <- Sigma_y
  Bay$tau_x2 <- tau_x2
  Bay$sigma_x2 <- sigma_x2
  Bay$tau_yx <- tau_yx
  Bay$sigma_yx <- sigma_yx
  Bay$tau_y2 <- tau_y2
  Bay$sigma_y2 <- sigma_y2
  Bay$D_x <- D_x
  Bay$V_x <- V_x
  Bay$D_yx <- D_yx
  Bay$V_yx <- V_yx
  Bay$D_y <- D_y
  Bay$V_y <- V_y
  Bay$L1 <- L1
  Bay$L2 <- L2
  Bay$MSE <- MSE
  Bay$MSE_ML <- MSE_ML
  Bay$W <- W
  Bay$Tau02 <- Tau02
  Bay$beta_b_ML <- beta_b_ML
  Bay$beta_b_Bay <- beta_b_Bay
  Bay$beta_b_Bay_ML <- beta_b_Bay_ML
  Bay$gamma <- gamma
  Bay$MSE_beta_b_Bay <- MSE_beta_b_Bay
  Bay$MSE_beta_b_Bay_ML <- MSE_beta_b_Bay_ML
  Bay$MSE_add <- MSE_add
  Bay$SE_beta_ML = SE_beta_ML
  Bay$SE_beta_Bay = SE_beta_Bay
  Bay$SE_beta_Bay_ML = SE_beta_Bay_ML
  Bay$SE_gamma = SE_gamma
  
  Bay_list <- list( 'A'=A,  'Sigma_x'=Sigma_x, 'Sigma_y'=Sigma_y, 'Sigma_yx'=Sigma_yx, 'tau_x2'=tau_x2, 
                    'sigma_x2'=sigma_x2, 'tau_y2'=tau_y2, 'sigma_y2'=sigma_y2,  'tau_yx'=tau_yx, 
                    'sigma_yx'=sigma_yx,'D_x'=D_x,  'D_y'=D_y, 'D_yx'=D_yx,'V_x'=V_x, 'V_y'=V_y, 
                    'V_yx'=V_yx,'L1'=L1, 'L2'=L2,  'MSE'=MSE,  'MSE_ML'=MSE_ML,'W'=W,'Tau02'=Tau02, 
                    'beta_b_ML'=beta_b_ML, 'beta_b_Bay'=beta_b_Bay,'beta_b_Bay_ML'=beta_b_Bay_ML, 
                    'gamma' = gamma, 'MSE_beta_b_Bay'=MSE_beta_b_Bay,'MSE_beta_b_Bay_ML'=MSE_beta_b_Bay_ML, 
                    'MSE_add'= MSE_add, 'SE_beta_ML' = SE_beta_ML, 'SE_beta_Bay' = SE_beta_Bay,
                    'SE_beta_Bay_ML' = SE_beta_Bay_ML, 'SE_gamma' = SE_gamma)
  
  #print(Bay_list)
}

#result_Bay_CV <- estimate_Bay_CV(data_CV)
#print(result_Bay_CV)
