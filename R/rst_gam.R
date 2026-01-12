
#' Fit the Global Estimation Model (RST_GAM)
#'
#' This function fits the robust ggam estimation model to account for potential outliers using bivariate and univariate splines.
#'
#' @param Y matrix. The response variable.
#' @param X matrix. The covariate matrix. Each column represents one covariate.
#' @param X_t matrix. A matrix of time-changing covariates, e.g., the accumulated infected cases and human mobility, each column represents a covariate. Default is NULL.
#' @param S matrix. The location matrix. Each row indicates the location of one point.
#' @param Tr matrix. The triangulation matrix. Each row indicates one triangle.
#' @param V matrix. The vertices matrix of the triangulation.
#' @param d integer. The degree of the bivariate spline, default is 2.
#' @param r integer. The smoothness of the bivariate spline, default is 1.
#' @param rho integer. The degree of the univariate spline, default is 2.
#' @param knots matrix. The knots matrix of the univariate spline. Usually, it is a (N+2) x p matrix. The r-th column is the knots of r-th covariate. The first element and the last element are the boundary knots.
#' @param N integer. The number of interior knots, default is 4.
#' @param initial_x numeric vector. The initial value for the parameter, default is 0.
#' @param k integer. The number data thinning folds, default is 3.
#' @param lambda1 numeric. The roughness penalty parameter for bivariate spline estimation.
#' @param lambda2 numeric. The L1 penalty term for the slack variable.
#' @param negative_binomial logical. Whether to fit Negative binomial model.
#' @param theta_range vector. A size 2 vector recording the range for hyperparameter of nb family (default is NULL - Poisson Family). If negative_binomial is TRUE and theta_range is not given, default is (0.01,10).
#' @param using_fix_weight logical. Whether to use poisson weight in the NB model to allow fast computation.
#' @param slack_weight vector. A vector of weight for slack estimates, default is null and will be obtained using data-thinning.
#' @param verbose logical. Whether to report model fit message.
#'
#' @return A list containing the following elements:
#' \item{theta_hat}{coefficients estimated for univariate functions.}
#' \item{gamma_hat}{coefficients estimated for bivariate function.}
#' \item{xi_hat}{the estimated slack parameter.}
#' \item{mu_hat}{the estimated mean parameter.}
#' \item{df}{degree of freedom of the model.}
#' \item{r_lambda}{the selected roughness penalization parameter.}
#' \item{l_lambda}{the selected L1 penalization parameter.}
#' \item{weight}{the adaptive slack weight used by the model.}
#' \item{Ind}{the index of the inside observation.}
#' \item{rho}{the degree of univariate spline, same as input argument.}
#' \item{d}{the degree of bivariate spline, same as input argument.}
#' \item{r}{the smoothness of bivariate spline, same as input argument.}
#' \item{knots}{the knots matrix of univariate spline, same as input argument.}
#' \item{N}{the number of interior knots, same as input argument.}
#' \item{sigma_2}{estimated sigma^2.}
#' \item{Y}{re-ordered response variable.}
#' \item{X}{re-ordered covariate matrix.}
#' \item{X_t}{re-ordered time-changing covariate matrix.}
#' \item{U}{univariate spline basis matrix.}
#' \item{BQ}{bivariate spline basis matrix.}
#' \item{P}{energy matrix for bivariate spline.}
#' \item{theta_nb}{the theta parameter for nb family, Null if fit quasi-poisson model}
#'
#' @import BPST
#' @import mgcv
#' @import splines
#' @import bindata
#' 
#' 
#' @examples
#' \donttest{
#' # Load the example dataset
#' data(example_dataset)
#' library(BPST)
#'
#' # the lambda list for roughness
#' lambda1=0.1
#' # the lambda list for l1
#' lambda2=0.02
#' 
#' # Fit the model using RST_GAM with 0 data thinning fold
#' result <- RST_GAM(
#'   Y = example_dataset$Y,
#'   X = example_dataset$X,
#'   S = example_dataset$S,
#'   Tr = example_dataset$Tr,
#'   V = example_dataset$V,
#'   lambda1 = lambda1,
#'   lambda2 = lambda2,
#'   k = 0
#' )
#'
#'
#' library(dplyr)
#' library(ggplot2)
#' library(mgcv)
#'
#' # obtain the summarized results for visualization
#' summary_plot <- RSTGAM_plot(result, S = example_dataset$S)
#' 
#' # Plot the identified outliers with corresponding estimated slack values
#' summary_plot$summary1 %>% ggplot(aes(x=X, y=Y, color=Slack)) +
#'   geom_point(alpha = 0.6) + 
#'   ggtitle("Outlier Signal Estimation") + 
#'   scale_colour_gradient(low = "white", high = "#bd0026", name = "Outlier Signal")
#'   
#'
#' # Plot the estimated bivariate component
#' summary_plot$summary1 %>% ggplot(aes(x=X, y=Y, color=Beta)) +
#'   geom_point(alpha = 0.6) + 
#'   ggtitle("Bivariate Component Estimation") + 
#'   scale_colour_gradient(low = "#7fcdbb", high = "#bd0026", 
#'   name = "Bivariate Function Value")
#'
#' # true beta
#' m=fs.test(example_dataset$S[,1][result$Ind], 
#'           example_dataset$S[,2][result$Ind], b=1)
#'          
#' true_beta = 0.5*m
#'
#' # Compute the bivariate MISE
#' mean((as.matrix(summary_plot$summary1$Beta) - true_beta)^2)
#'
#' 
#' # Define the true coefficient functions
#'
#' beta.1 = function(x0) 1/2 * sin(2 * pi * x0) - cos(pi * x0)
#' beta.2 = function(x0) 4 * ((x0 - 0.5)^2 - 2 / 3 * (0.5)^3)
#' beta.3 = function(x0) x0 - 0.5
#'
#'
#' # the estimation from model
#' X0 = summary_plot$X0
#' mhat = summary_plot$summary2
#' 
#' 
#' # Calculate true values for given data
#' beta0.1 = beta.1(X0[, 1])
#' beta0.2 = beta.2(X0[, 2])
#' beta0.3 = beta.3(X0[, 3])
#' 
#' 
#' # Set up plotting margins and layout
#' par(mar = c(5, 9, 4, 2))
#' par(mfrow = c(1, 3))  # Arrange plots in a 2x2 grid
#'
#' # Plot for X1
#' plot(X0[, 1], mhat[, 1], xlab = 'X1', type = 'l', lwd = 3, ylab = '', 
#'      cex.axis = 1.7, cex.lab = 2, cex = 1.3, pch = 16, col = 'blue')
#' lines(X0[, 1], beta0.1, lwd = 3, col = 'red')
#' legend(0, 1, legend = c('True', 'Estimates with Slack'), 
#'        fill = c('red', 'blue'), cex = 1.6)
#'
#' # Plot for X2
#' plot(X0[, 2], mhat[, 2], xlab = 'X2', type = 'l', lwd = 3, ylab = '', 
#'      cex.axis = 1.7, cex.lab = 2, cex = 1.3, pch = 16, col = 'blue')
#' lines(X0[, 2], beta0.2, lwd = 3, col = 'red')
#'
#' # Plot for X3
#' plot(X0[, 3], mhat[, 3], xlab = 'X3', type = 'l', lwd = 3, ylab = '', 
#'      cex.axis = 1.7, cex.lab = 2, cex = 1.3, pch = 16, col = 'blue')
#' lines(X0[, 3], beta0.3, lwd = 3, col = 'red')
#'
#'
#' # the MISE for each univarite component
#' sum((beta0.1-mhat[,1])^2)/1000
#' sum((beta0.2-mhat[,2])^2)/1000
#' sum((beta0.3-mhat[,3])^2)/1000
#'}
#'

#' @export

RST_GAM = function(Y,X,X_t=NULL,S,Tr,V,d=2,r=1,rho=2,knots=NULL,N=4,initial_x = NULL, k = 3,
                       lambda1, lambda2, negative_binomial = FALSE, theta_range = NULL, using_fix_weight = TRUE, slack_weight = NULL, verbose = TRUE){
  
  if (negative_binomial & is.null(theta_range)){
    theta_range = c(0.01, 10) # default theta value for nb family
  }

  # number of data sample
  t = ncol(Y)
  
  ####################### stage 1 Spline Construction ###############################
  
  if(verbose) message('Stage 1: Construct Basis Splines')
  if (!requireNamespace("BPST", quietly = TRUE)) {
    stop("The 'BPST' package is required but not installed. Please install it using remotes::install_github('FIRST-Data-Lab/BPST').")
  }
  # generate basis function for S using bivariate spline
  ## bivariate spline basis matrix
  B0=basis(V,Tr,d,r,S)
  Bi=B0$Bi
  Q2=B0$Q2
  BQ2=as.matrix(Bi%*%Q2)
  Ind=B0$Ind
  
  stack_BQ2 = matrix(rep(t(BQ2), t), ncol = ncol(BQ2), byrow = T)
  
  ## energy matrix P
  P=B0$K
  P=t(Q2)%*%P%*%Q2
  
  # the Y matrix
  Y=Y[Ind, , drop =  F]
  
  # updated number of row
  n = nrow(Y)
  
  # generate basis function for X and X_t using univariate spline
  X=X[Ind, , drop = F]
  
  if(is.null(X_t)){
    
    X_all = matrix(rep(t(X), t), ncol = ncol(X), byrow = T)
    
  }else{
    
    X_t = X_t[((rep(0:(t-1), each = n)*nrow(S)) + rep(Ind, t)), , drop = F]
    
    X_all = cbind(matrix(rep(t(X), t), ncol = ncol(X), byrow = T), X_t)
    
  }
  
  X.p=dim(X_all)[2]
  XX=NULL
  p=NULL
  if(is.null(knots)){
    # Basis for X
    for(alpha in 1:X.p){
      UniBasis=Basis_generator(X_all[,alpha, drop = F],N=N,q=rho+1)
      Xalpha=UniBasis$Bx
      knots=cbind(knots,UniBasis$knots)
      XX=cbind(XX,Xalpha)
      p=c(p,dim(Xalpha)[2])
    }
  } else {
    # Basis for X
    for(alpha in 1:X.p){
      Xalpha=Basis_generator(X_all[,alpha, drop = F],N=N,q=rho+1,knots=knots[,alpha])$Bx
      XX=cbind(XX,Xalpha)
      p=c(p,dim(Xalpha)[2])
    }
  }
  XX_all=as.matrix(XX)
  colnames(XX_all) = NULL
  
  if(verbose) message('Stage 1 Done')
  
  ####################### stage 2 Model Fitting ###############################
  
  if(verbose) message('Stage 2: Model Fitting')
  if (negative_binomial){
    if(verbose) message(sprintf('Fit Negative Binomial Model with Theta Range: [%s, %s]', theta_range[1], theta_range[2]))
    
    # construct weights using Poisson model
    if (using_fix_weight & is.null(slack_weight)){
      detected_weights = RST_GAM_poisson(Y, XX_all, stack_BQ2, P, lambda1, lambda2, initial_x, k, Ind, fix_weight = NULL, report_weight = TRUE)
    }else if(using_fix_weight & !is.null(slack_weight)){
      detected_weights = slack_weight
    }else{
      detected_weights = NULL
    }
    # fit nb model with fixed slack weight
    final_result = RST_GAM_nb(Y, XX_all, stack_BQ2, P, lambda1, lambda2, initial_x, k, Ind, theta_range = theta_range, fix_weight = detected_weights)
  }
  else{
    if(verbose) message('Fit (Quasi)-Poisson Model')
    # fit quasi-poisson model
    final_result = RST_GAM_poisson(Y, XX_all, stack_BQ2, P, lambda1, lambda2, initial_x, k, Ind, fix_weight = slack_weight)
    final_result$theta_nb = NULL
  }
  
  if(verbose) message('Stage 2 done')
  
  
  final_result$Ind = Ind
  final_result$rho=rho
  final_result$knots=knots
  final_result$d=d
  final_result$r=r
  final_result$N=N
  final_result$Y=Y
  final_result$X=X
  final_result$X_t=X_t
  final_result$U = XX_all
  final_result$BQ = BQ2
  final_result$P = P
  
  return(final_result)
  
}


