#' Summarize and Prepare RST_GAM Results for Plotting
#'
#' This function processes the output of the robust gGAM estimation model (RST_GAM) to generate summaries and prepare data for visualization.
#'
#' @param model list. The output from the RST_GAM model.
#' @param uni_sample integer. The number of evenly spaced points used to evaluate univariate components, default is 1000.
#' @param S matrix. The location matrix. Each row indicates the location of one point.
#'
#' @return A list containing the following elements:
#' \item{summary1}{A tibble summarizing estimated slack values and bivariate components.}
#' \item{X0}{A matrix of covariate values at equally spaced points, with the number of points determined by uni_sample.}
#' \item{summary2}{A matrix for each estimated univariate component, evaluated at X0.}
#'
#' 
#' @import tibble
#' 
#' @examples
#' 
#' # Refer to the detailed example in the documentation for the RST_GAM function.
#' # This function is typically used after running RST_GAM to process its output
#' # and prepare data for visualization.
#' 
#' 
#' @seealso RST_GAM
#' 
#' @export

RSTGAM_plot = function(model, uni_sample = 1000, S){
  
  # obtain the bivariate function estimation
  est_bivariate = as.vector(model$BQ %*% model$gamma_hat)
  
  # obtain the slack estimation
  est_slack = as.vector(model$xi_hat)
  
  # summarize into a dataframe
  summary1 = tibble(X = S[,1][model$Ind], Y = S[,2][model$Ind],
                    Slack = est_slack, Beta = est_bivariate)
  
  # obtain the equally spaced points for each components
  p_ind = ncol(model$X)
  X_t_exists = !is.null(model$X_t)  # Check if X_t is NULL
  p_dep = if (X_t_exists) ncol(model$X_t) else 0
  nP = p_ind + p_dep
  X0 = matrix(0, ncol = nP, nrow = uni_sample)
  
  # Time-independent components
  for (i in 1:p_ind) {
    X0[, i] = seq(min(model$X[, i]), max(model$X[, i]), length.out = uni_sample)
  }
  
  # Time-dependent components (only if X_t is not NULL)
  if (X_t_exists) {
    for (j in (p_ind + 1):nP) {
      X0[, j] = seq(min(model$X_t[, (j - p_ind)]), 
                    max(model$X_t[, (j - p_ind)]), length.out = uni_sample)
    }
  }
  
  # generate the basis matrix
  X.p=dim(X0)[2]
  XX=NULL
  p=NULL
  for(alpha in 1:X.p){
    Xalpha=Basis_generator(X0[,alpha],N=model$N,q=model$rho+1,
                           knots=model$knots[,alpha])$Bx
    XX=cbind(XX,Xalpha)
    p=c(p,dim(Xalpha)[2])
  }
  XX=as.matrix(XX)
  
  # spline estimator
  theta=model$theta_hat
  p=c(0,cumsum(p))
  mhat=NULL
  for(alpha in 1:X.p){
    m=XX[,(p[alpha]+1):(p[alpha+1])]%*%theta[(p[alpha]+1):(p[alpha+1])]
    mhat=cbind(mhat,m)
  }
  
  final_result = list(summary1 = summary1, X0 = X0, summary2 = mhat)
  
  return(final_result)
  
}