
#' Gradient Map Function
#'
#' Computes the gradient map for given parameters.
#'
#' @param Y A numeric matrix.
#' @param X A numeric matrix.
#' @param BQ2 A numeric matrix.
#' @param P A numeric matrix.
#' @param cur_x A numeric vector.
#' @param lambda1 A numeric value.
#' @param theta_nb A numeric value (default is 0).
#' @return A numeric vector of the gradient map.
#' 
#' @useDynLib RSTGAM, .registration = TRUE
#' @import RcppArmadillo
#' @importFrom Rcpp evalCpp
#' 
#' 
#' @keywords internal

# @importFrom Rcpp sourceCpp
gradient_map <- function(Y, X, BQ2, P, cur_x, lambda1, theta_nb = 0){
  
  family = poisson_fam()
  variance = family$variance
  linkinv = family$linkinv
  mu.eta = family$mu.eta
  
  # size of sample
  n = nrow(Y)
  t = ncol(Y)
  N = n*t
  
  # size of univariate basis coefficients parameters
  sum_pk = ncol(X)
  
  # size of bivariate basis coefficients parameters
  q = ncol(P)
  
  # total # of parameter
  m = n + sum_pk + q
  
  theta <- cur_x[1:sum_pk]
  gamma <- cur_x[(sum_pk + 1):(sum_pk + q)]
  xi <- cur_x[(sum_pk + q + 1):m]
  
  # a vector to compute the eta
  eta <- compute_eta(as.matrix(X), theta, as.matrix(BQ2), gamma, xi, t)
  
  # vectorize Y
  Y_vec <- matrix(c(Y), nrow = N, ncol = 1)
  
  # perform coordinate gradient computation
  # first compute the common term
  if(theta_nb != 0){
    temp1 = linkinv(eta)
    temp2 = -(Y_vec - temp1)*mu.eta(eta)/N
    temp3 = temp2/variance(temp1, theta_nb)
    
  }else{
    temp3 = -(Y_vec - linkinv(eta))/N
  }
  
  #for theta
  gradient_theta <- compute_gradient_theta(X, temp3)
  
  # for gamma
  gradient_gamma <- compute_gradient_gamma(as.matrix(BQ2), temp3, as.matrix(P), gamma, lambda1)
  
  # for xi
  temp3_matrix <- matrix(temp3, ncol = t)
  gradient_xi <- rowSums(temp3_matrix)
  
  # combine all three components
  gradient_total <- as.matrix(c(gradient_theta, gradient_gamma, gradient_xi))
  
  rownames(gradient_total) = NULL
  colnames(gradient_total) = NULL
  
  return(gradient_total)
}

