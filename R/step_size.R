#' this function implementing the step size calculation from the adaptive proximal GD method by Yura Malisky 
#'
#' @param pre_step: the previous step size
#' @param pre_theta: the previous theta parameter
#' @param cur_L: the current step local Lipschitz constant estimate
#' @return a value represent the current step size
#' 
#' 
#' @keywords internal
adaptive_step_size <- function(pre_step, pre_theta, cur_L) {


  # the first bound
  bound_1 = sqrt(2/3 + pre_theta) * pre_step
  
  # the second bound
  bound_2 = pre_step/sqrt(max(c(0, 2 * pre_step^2 * cur_L^2 - 1)))
  
  # get the minimum
  cur_step = min(c(bound_1, bound_2))
  
  return(cur_step)
  
}


#' this function implementing a linesearch to find the initial step size
#'
#' @param pre_x: the previous x
#' @param initial_step: the starting step size, default as 10^-9
#' @param Y: the sample observation
#' @param X: the univariate basis function matrix
#' @param BQ2: the bivariate matrix after QR decomposition, n by q matrix 
#' @param P: penalty matrix (q by q matrix)
#' @param lambda1: the penalty parameter for roughness
#' @param lambda2: the penalty parameter for xi
#' @param w: the weight vector for l1 penalty
#' @return a value represent the initial step size
#' 
#' 
#' @keywords internal
initial_step_linesearch <- function(pre_x, initial_step = 10^-8, Y, X, BQ2, P, lambda1, lambda2, w) {
  
  cur_a = initial_step
  
  tol = 0
  
  search_seq = exp(seq(log(initial_step),log(10),length.out=200))
  
  pre_gradient = gradient_map(Y, X, BQ2, P, pre_x, lambda1)
  
  for (cur_a in search_seq) {
    
    cur_z = pre_x - cur_a * pre_gradient
    
    cur_x = proximal_map(cur_z, lambda2, cur_a, n = nrow(Y), w = w)
    
    cur_gradient = gradient_map(Y, X, BQ2, P, cur_x, lambda1)
    
    # compute L
    diff_gradient <- cur_gradient - pre_gradient
    diff_x <- cur_x - pre_x
    cur_L <- norm(diff_gradient, type = "F") / norm(diff_x, type = "F")
    
    
    if(is.nan(cur_L)){
      
      return(cur_a/100)
      break
      
    }
    
    tol = cur_L * cur_a
    
    
    if(tol >= 1/sqrt(2) & tol <= 2){
      
      return(cur_a)
      break
      
    }
    
    
  }
  
  return(cur_a/100)
  
}




