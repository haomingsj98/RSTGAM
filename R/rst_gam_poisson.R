
#' Fit the Global Estimation Model (RST_GAM)
#'
#' This function fits the robust ggam estimation model to account for potentail outliers using bivariate and univariate splines with poisson family.
#'
#' @param Y matrix. The response variable.
#' @param X matrix. The covariate matrix.
#' @param BQ matrix. (Bernstein basis or similar - inferred from context).
#' @param P matrix. (Penalty matrix or similar - inferred from context).
#' @param lambda1 numeric. The roughness penalty parameter.
#' @param lambda2 numeric. The L1 penalty term.
#' @param initial_x numeric vector. Initial values.
#' @param k integer. Number of data thinning folds.
#' @param Ind vector. Indices for the subset.
#' @param theta_nb numeric. Hyperparameter for nb family (default 0).
#' @param fix_weight vector. Weights for slack estimation.
#' @param report_weight boolean. Whether to just report slack weights.
#' @param theta_nb: numeric. The hyperparameter for nb family (default is 0 - Poisson family).
#' @param fix_weight vector. The weight vector for slack estimation.
#' @param report_weight boolean. Whether to just report slack weights, used for NB model fit.
#'
#' @return A list containing the following elements:
#' \item{theta_hat}{coefficients estimated for univariate functions.}
#' \item{gamma_hat}{coefficients estimated for bivariate function.}
#' \item{xi_hat}{the estimated slack parameter.}
#' \item{mu_hat}{the estimated mean parameter.}
#' \item{df}{degree of freedom of the model.}
#' \item{r_lambda}{the selected roughness penalization parameter.}
#' \item{l_lambda}{the selected L1 penalization parameter.}
#' \item{sigma_2}{estimated pearson sigma^2.}
#'
#' @importFrom stats rmultinom 
#' 
#' 
#' @keywords internal

RST_GAM_poisson = function(Y, X, BQ, P, lambda1, lambda2, initial_x, k, Ind, theta_nb=0, fix_weight = NULL, report_weight = FALSE){
  
  if (theta_nb != 0){
    family = nb_fam()
  }
  else{
    family = poisson_fam()
  }
  variance = family$variance
  linkinv = family$linkinv
  
  # number of data sample
  t = ncol(Y)
  n = nrow(Y)
  
  ####################### stage 1 weights construction ###############################
  if(!is.null(fix_weight)){
    
    # Use pre-calculated weights (Fast Path for NB Loop)
    final_w = fix_weight
  
  }else if(k == 0){ # if no slice, just perform naive Lasso with weights 1
    
    final_w = rep(1, length(Ind))
    
  }else{
    # perform k-fold data thinning
    Y_slice = sapply(c(Y), rmultinom, n = 1, prob = rep(1/k, k))
    
    # set default weight for slack estimates
    final_w = rep(0, length(Ind))
    
    
    for (i in 1:k) {
      
      # the i-th slice
      if(k == 1){
        
        Y_i = Y
        
      }else{
        
        Y_i = matrix(Y_slice[i,], ncol = t)
        
      }
      
      # run the model
      cur_result = rgam_lambda_selection(Y_i, X, BQ, P, lambda1, lambda2, initial_x, w = rep(1, length(Ind)), theta_nb = theta_nb)
      
      # record the slack estimates
      final_w = final_w + as.vector(cur_result$xi_hat)
      
    }
    
    # obtain final weights
    final_w = 1/(final_w/k)
    
  }
  
  if (report_weight) {return(final_w)}
  
  ################################### Stage 2 Full Model Construction #########################################

  # print('Stage 2: Penalization Parameter Selection')
  
  # compute weighted estimation with Y2 to get penalization parameters
  final_result = rgam_lambda_selection(Y, X, BQ, P, lambda1, lambda2, initial_x, w = final_w, theta_nb = theta_nb)
  
  final_rlambda = final_result$r_lambda
  final_llambda = final_result$l_lambda
  
  
  # compute model effective degrees of freedom
  df = compute_edf(X, BQ, P, final_result$xi_hat, final_rlambda)
  final_result$df = df
  
  # compute sigma
  y_hat = as.vector(final_result$mu_hat * linkinv(rep(final_result$xi_hat, t)))
  
  if (theta_nb != 0){
    tt = (c(Y)-y_hat)^2/variance(y_hat, theta_nb)
  }
  else{
    tt = (c(Y)-y_hat)^2/variance(y_hat)
  }
  
  sigma_2 = 1/(n*t-final_result$df)*sum(tt)
  
  final_result$sigma_2=sigma_2

  final_result$r_lambda = final_rlambda
  final_result$l_lambda = final_llambda

  
  return(final_result)
  
}


