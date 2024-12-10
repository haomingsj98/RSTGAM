


#' this function is an intermediate step of model estimation procedure, implementing the gradient descent to get parameter estimates 
#'
#' @param Y: the response variable.
#' @param X: the covariate matrix (univariate basis).
#' @param BQ2: the bivariate matrix after QR decomposition.
#' @param P: penalty matrix.
#' @param lambda1: the tuning parameter for roughness.
#' @param lambda2: the l1 penalty parameter for slack variable.
#' @param initial_x: the initial value for the parameter, vector size of m (default as 0).
#' @param w: the weight vector for l1 penalty. 
#' @return a list include parameter estimations: 
#' \item{theta_hat}{coefficients estimated for univariate functions.}
#' \item{gamma_hat}{coefficients estimated for bivariate function.}
#' \item{xi_hat}{coefficients estimated for slack variable.}
#' 
#' @keywords internal
rgam_estimation = function(Y,X,BQ2,P,lambda1, lambda2, initial_x = NULL, w){
  
  # size of sample
  n = nrow(Y)
  t = ncol(Y)
  
  # size of univariate basis coefficients parameters
  sum_pk = ncol(X)
  
  # size of bivariate basis coefficients parameters 
  q = ncol(P)
  
  # total # of parameter
  m = n + sum_pk + q
  
  if(is.null(initial_x)){
    
    # set initial value as 0
    pre_x = matrix(0, nrow = m, ncol = 1)
  }else{
    
    pre_x = initial_x
    
  }

  # record previous y
  pre_y = pre_x
  
  # set initial step size alpha
  pre_alpha = initial_step_linesearch(pre_x, Y = Y, X = X, BQ2 = BQ2, P = P, lambda1 = lambda1, lambda2 = lambda2, w = w)
  
  # set initial hyperparatmer theta
  cur_theta = 1/3
  
  pre_gradient = gradient_map(Y, X, BQ2, P, cur_x = pre_x, lambda1 = lambda1)

  
  # compute first update of x
  cur_z = pre_x - pre_alpha * pre_gradient
  cur_x = proximal_map(cur_z, lambda = lambda2, h = pre_alpha, n = n, w = w)
  
  
  # record euclidean difference
  tol = norm(cur_x - pre_x, 'F')

  k = 1
  
  # start the algorithm
  while (tol > 10^-5 & k <= 10^5) {
    
    # the momentum computation
    cur_y = cur_x + (k-2) * (cur_x - pre_x)/(k+1) 
    
    # compute current gradient
    cur_gradient = gradient_map(Y, X, BQ2, P, cur_x = cur_y, lambda1 = lambda1)
    

    # compute the local lipschitz approximation
    diff_gradient = cur_gradient - pre_gradient
    diff_y = cur_y - pre_y
    cur_L = norm(diff_gradient, 'F')/norm(diff_y, 'F')
    
    if(is.nan(cur_L)){
      
      cur_alpha = 0
      
    }else{
      
      # compute current step size
      cur_alpha = adaptive_step_size(pre_alpha, cur_theta, cur_L)
      
    }
    
    # compute new x
    # cur_z = cur_x - cur_alpha * cur_gradient
    cur_z = cur_y - cur_alpha * cur_gradient
    next_x = proximal_map(cur_z, lambda = lambda2, h = cur_alpha, n = n, w = w)
    
    # update tolerance
    #tol = norm(next_x - cur_x, 'F')/(1 + norm(cur_x, 'F'))
    tol = sqrt(mean((next_x - cur_x)^2))

    
    # update theta
    cur_theta = cur_alpha/pre_alpha
    
    # update previous alpha
    pre_alpha = cur_alpha
    
    # update x
    pre_x = cur_x
    cur_x = next_x
    
    # update y
    pre_y = cur_y
    
    # update previous gradient
    pre_gradient = cur_gradient
    
    k = k + 1

  }

  result = list(theta_hat = cur_x[1:sum_pk, 1, drop = F], gamma_hat = cur_x[(sum_pk + 1):(sum_pk + q), 1, drop = F], xi_hat = cur_x[(sum_pk + q + 1):m, 1, drop = F])

  return(result)

  
}



