#' this function helps calculate poisson log likelihood of our model with slack variables
#'
#' @param xi_is: the estimated slack variables from a given model
#' @param Y_it: a n by t matrix representing the observations (counts) at time t (column) and location i (row)
#' @param mu_t: a size n (location) vector representing the value of mean counts at time t
#' @return a value represent the log likelihood
#' 
#' 
#' @importFrom stats dpois
#' 
#' @keywords internal
#' 
l_likelihood = function(xi_is, Y_it, mu_t){
  
  # size of time (column)
  t = ncol(Y_it)
  
  # size of location
  n = nrow(Y_it)
  
  # vectorize Y_it and a seudo u_it
  Y_vec = matrix(c(Y_it), nrow = n*t, ncol = 1)
  
  u_vec = matrix(as.vector(mu_t), nrow = n*t, ncol = 1)
  
  # compute the likelihood
  exps_it = rep(exp(xi_is), t)
  
  u_exps = u_vec * exps_it
  
  ll = sum(dpois(Y_vec, u_exps, log = T))
  
  return(ll)
  
}


#' this function compute the EBIC measure of the Lasso penalty
#'
#' @param X: the univariate basis function matrix
#' @param BQ2: the bivariate matrix after QR decomposition, n by q matrix 
#' @param P: penalty matrix (q by q matrix)
#' @param xi_is: the estimated slack variables from a given model
#' @param Y_it: a n by t matrix representing the observations (counts) at time t (column) and location i (row) description
#' @param mu_t: a size n (location) vector representing the value of mean counts at time t description
#' @param gamma: a regulatory parameter, default as 0, note, when gamma is 0, we obtain BIC
#' @return a value represent the eBIC
#' @keywords internal
#' 
EBIC_L = function(X, BQ2, P, xi_is, Y_it, mu_t, gamma = 0){
  
  # obtain the (log) quasi_likelihood of the model
  log_likelihood = l_likelihood(xi_is, Y_it, mu_t)
  
  # number of sample
  N = nrow(Y_it) * ncol(Y_it)
  
  # number of slack
  n = nrow(xi_is)
  
  # number of total model parameters
  total_P = n
  
  # compute edf
  v = sum(xi_is != 0)
  
  # get the number of combinations (non-zero parameters) from the whole possible parameters
  log_pCv = lfactorial(total_P) - (lfactorial(v) + lfactorial(total_P-v))
  
  # compute the eBIC
  nquasi = -2 * log_likelihood
  
  v_logn = v * log(N)
  
  gamma_log_pCv = 2 * gamma * log_pCv
  
  eBIC = nquasi + v_logn + gamma_log_pCv
  
  return(eBIC)
  
}

#' this function compute the EBIC measure of the model
#'
#' @param X: the univariate basis function matrix
#' @param BQ2: the bivariate matrix after QR decomposition, n by q matrix 
#' @param P: penalty matrix (q by q matrix)
#' @param xi_is: the estimated slack variables from a given model
#' @param Y_it: a n by t matrix representing the observations (counts) at time t (column) and location i (row) description
#' @param mu_t: a size n (location) vector representing the value of mean counts at time t description
#' @param lambda_r: the penalty parameter for roughness
#' @param gamma: a regulatory parameter, default as 1, note, when gamma is 0, we obtain BIC
#' @return a value represent the eBIC
#' @keywords internal
#' 
EBIC = function(X, BQ2, P, xi_is, Y_it, mu_t, lambda_r, gamma = 1){
  
  # obtain the (log) quasi_likelihood of the model
  log_likelihood = l_likelihood(xi_is, Y_it, mu_t)
  
  # number of sample
  N = nrow(Y_it) * ncol(Y_it)
  
  # number of slack
  n = nrow(xi_is)
  
  # number of univariate parameter
  p = ncol(X)
  
  # number of bivariate paramters
  q = ncol(BQ2)
  
  # number of total model parameters
  total_P = n + p + q
  
  # compute edf
  v = compute_edf(X, BQ2, P, xi_is, lambda_r)
  
  
  # get the number of combinations (non-zero parameters) from the whole possible parameters
  log_pCv = lfactorial(total_P) - (lfactorial(v) + lfactorial(total_P-v))
  
  # compute the eBIC
  nquasi = -2 * log_likelihood
  
  v_logn = v * log(N)
  
  gamma_log_pCv = 2 * gamma * log_pCv
  
  eBIC = nquasi + v_logn + gamma_log_pCv
  
  return(eBIC)
  
}



#' this function compute the model effective degrees of freedom
#'
#' @param X: the univariate basis function matrix
#' @param BQ2: the bivariate matrix after QR decomposition, n by q matrix 
#' @param P: penalty matrix (q by q matrix)
#' @param xi_is: the estimated slack variables from a given model
#' @param lambda_r: the penalty parameter for roughness
#' @return a value represent the effective degrees of freedom of the model
#' @import Matrix
#' @keywords internal
#' 
compute_edf <- function(X, BQ2, P, xi_is, lambda_r) {
  
  N = nrow(xi_is)
  t = nrow(X) / N
  p = ncol(X)
  q = ncol(BQ2)
  
  I = Matrix(rep(t(diag(1, nrow = N)), t), ncol = N, byrow = TRUE, sparse = TRUE)
  I_sub = I[, which(xi_is > 0), drop = FALSE]
  n_sub = ncol(I_sub)
  
  H = cbind(X, BQ2, I_sub)
  
  zero1 = Matrix(0, nrow = p, ncol = p + q + n_sub, sparse = TRUE)
  zero2 = Matrix(0, nrow = q, ncol = p, sparse = TRUE)
  zero3 = Matrix(0, nrow = q, ncol = n_sub, sparse = TRUE)
  zero4 = Matrix(0, nrow = n_sub, ncol = p + q + n_sub, sparse = TRUE)
  
  R = rbind(cbind(zero2, P, zero3), zero4)
  R = rbind(zero1, R)
  
  temp1 = crossprod(H, H) + lambda_r * R
  
  # Compute the SVD of A
  svd_temp = svd(temp1)
  
  # Extract components
  U = svd_temp$u
  D = svd_temp$d  # Singular values (a vector, not a matrix)
  V = svd_temp$v
  
  # Invert the non-zero singular values
  D_pseudo = diag(c(1 / D[D > 1e-10], rep(0, sum(D <= 1e-10))))  # Create a full diagonal matrix
  
  # Construct the pseudoinverse of A
  inv_temp1 = V %*% tcrossprod(D_pseudo, U)
  
  temp2 = H %*% tcrossprod(inv_temp1, H)
  
  edf = sum(diag(temp2))
  
  return(as.numeric(edf))
}


