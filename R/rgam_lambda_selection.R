
#' Internal Function: Select Penalization Parameters for RST_GAM
#'
#' This internal function fits the global estimation model with selected penalization parameters.
#'
#' @param Y matrix. The response variable.
#' @param X matrix. The covariate matrix consisting of univariate basis.
#' @param BQ matrix. The matrix for bivariate basis.
#' @param P matrix. The smoothness penalty matrix.
#' @param lambda1 numeric. The roughness penalty parameter for bivariate spline estimation.
#' @param lambda2 numeric. The L1 penalty term for the slack variable.
#' @param initial_x numeric vector. The initial value for the parameter, default is 0.
#' @param w numeric vector. The weight vector for slack estimation.
#'
#' @return A list containing the following elements:
#' \item{theta_hat}{coefficients estimated for univariate functions.}
#' \item{gamma_hat}{coefficients estimated for bivariate function.}
#' \item{xi_hat}{the estimated slack parameter.}
#' \item{r_lambda}{the selected roughness penalization parameter.}
#' \item{l_lambda}{the selected L1 penalization parameter.}
#' \item{ebic}{the final ebic value}
#'
#' @keywords internal


rgam_lambda_selection = function(Y, X, BQ, P, lambda1, lambda2, initial_x, w){
  
  family = poisson_fam()
  linkinv = family$linkinv

  
  
  # sample size
  n = nrow(Y)
  t = ncol(Y)
  
  # a list to record final estimation for optimal lambda combination
  final_result = list()
  
  # length of roughness lambda
  n1 = length(lambda1)
  
  # length of l1 lambda
  n2 = length(lambda2)
  
  ebic_result = 10^8
  
  select_l = 0
  select_r = 0
  
  # select lambda smoothing parameter with 0 roughness lambda
  
  # value for plot
  ebic_l_list = c()
  l_lambda_list = c()
  
  cur_l_ebic = 10^8
  for (i in 1:n2) {
    
    # roughness lambda set 0
    r_lambda = 0
    
    # l1 lambda
    l_lambda = lambda2[n2+1-i]
    
    # run the model
    result = rgam_estimation(Y, X, BQ, P, r_lambda, l_lambda, initial_x = initial_x, w = w)
    result$mu_hat = linkinv(X %*% result$theta_hat + BQ %*% result$gamma_hat)
    
    # compute current ebic
    cur_ebic = EBIC_L(X, BQ, P, result$xi_hat, Y, result$mu_hat)

    # if ebic decreases, update 
    if(cur_ebic <= cur_l_ebic){
      
      # record the ebic and lambda value
      cur_l_ebic = cur_ebic
      select_l = l_lambda
      
      
    }
    else{
      # if ebic starts to increase, break the loop
      break
    }
    
  }
  
  # select roughness smoothing parameter with selected lambda
  
  min_ebic = 10^8
  
  for (j in 1:n1) {
    
    # selected lambda
    l_lambda = select_l
    
    # roughness lambda
    r_lambda = lambda1[j]
    
    # run the model
    result = rgam_estimation(Y, X, BQ, P, r_lambda, l_lambda, initial_x = initial_x, w = w)
    result$mu_hat = linkinv(X %*% result$theta_hat + BQ %*% result$gamma_hat)
    
    # compute current ebic
    cur_ebic = EBIC(X, BQ, P, result$xi_hat, Y, result$mu_hat, r_lambda)
    
    if(cur_ebic <= min_ebic){
      
      min_ebic = cur_ebic
      final_result = result
      select_r= r_lambda
      final_result$r_lambda = select_r
      final_result$l_lambda = select_l
      final_result$ebic = min_ebic
      
    }
    
  }
  
  return(final_result)
  
}


