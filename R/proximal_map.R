#' this function helps compute the gradient update during each iteration with given gradient vector
#'
#' @param cur_gradient: the current step gradient
#' @param lambda: the penalty parameter for xi
#' @param h: the step size
#' @param n: size of the xi vector, same as the number of observation
#' @param w: weight vector, same size as the xi vector
#' @return a vector (same size the gradient vector) represent the proximal updates
#' @keywords internal
#' 
proximal_map = function(cur_gradient, lambda, h, n, w){

  # size of the gradient vector
  m = nrow(cur_gradient)


  # for paramters except xi, the proximal mapping is just identity map
  proximal_update = cur_gradient

  # the threshold value
  threshold = h*lambda*w

  # Compute indices for xi components
  xi_indices <- (m - n + 1):m

  # Vectorized soft thresholding for xi parameters
  cur_xi_gradient <- cur_gradient[xi_indices, 1]
  proximal_update[xi_indices, 1] <- pmax(cur_xi_gradient - threshold, 0)

  return(proximal_update)

}
