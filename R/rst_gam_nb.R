
#' Fit the Global Estimation Model (RST_GAM)
#'
#' This function fits the robust ggam estimation model to account for potentail outliers using bivariate and univariate splines with Negative Binomial family.
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
#' @param theta_range numeric vector. A vector of length 2 (lower, upper) specifying the search range for the NB theta parameter.
#' @param fix_weight vector. The weight vector for slack estimation.
#' 
#'
#' @return A list containing the following elements:
#' \item{theta_hat}{coefficients estimated for univariate functions.}
#' \item{gamma_hat}{coefficients estimated for bivariate function.}
#' \item{xi_hat}{the estimated slack parameter.}
#' \item{mu_hat}{the estimated mean parameter.}
#' \item{df}{degree of freedom of the model.}
#' \item{r_lambda}{the selected roughness penalization parameter.}
#' \item{l_lambda}{the selected L1 penalization parameter.}
#' \item{sigma_2}{estimated sigma^2.}
#' \item{theta_nb}{the estimated theta hyperparameter for nb family}
#'
#' 
#' 
#' @keywords internal

RST_GAM_nb = function(Y, X, BQ, P, lambda1, lambda2, initial_x, k, Ind, theta_range, fix_weight){
  
  
  # find the optimal theta value for nb family matching the pearson measure
  lower=theta_range[1]
  upper=theta_range[2]
  result1=RST_GAM_poisson(Y, X, BQ, P, lambda1, lambda2, initial_x, k, Ind, theta_nb = upper, fix_weight = fix_weight)
  if (lower == upper){
    # return result if theta is fixed
    result1$theta_nb = upper
    return(result1)
    
  }
  U_sigma=result1$sigma_2
  
  result2=RST_GAM_poisson(Y, X, BQ, P, lambda1, lambda2, initial_x, k, Ind, theta_nb = lower, fix_weight = fix_weight)
  L_sigma=result2$sigma_2
  
  
  M_sigma=2
  step=1
  
  if ( (L_sigma-1)> 0 & (U_sigma-1) > 0 ) {
    M_sigma=1
    final_result=result2
    middle=lower
    message(sprintf("NB Converged at Lower Bound: Theta=%s, Pearson=%s", round(middle,4), round(L_sigma,4)))
  }
  
  if ( (L_sigma-1)<0 & (U_sigma-1) < 0 ) {
    M_sigma=1
    final_result=result1
    middle=upper
    message(sprintf("NB Converged at Upper Bound: Theta=%s, Pearson=%s", round(middle,4), round(U_sigma,4)))
  }

  # Bisection Loop
  while(abs(M_sigma-1)>=0.01 & step <= 10){
    middle=(upper+lower)/2
    final_result=RST_GAM_poisson(Y, X, BQ, P, lambda1, lambda2, initial_x, k, Ind, theta_nb = middle, fix_weight = fix_weight)
    M_sigma=final_result$sigma_2
    
    message(sprintf("Step %d: Theta = %.4f | Pearson Measure = %.4f", step, middle, M_sigma))
    
    if((M_sigma-1)*(U_sigma-1)<0){
      lower=middle
      L_sigma=M_sigma
    } else {
      upper=middle
      U_sigma=M_sigma
    }
    step=step+1
  }
  final_result$theta_nb = middle
  return(final_result)
  
}


