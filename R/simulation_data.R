
#' Generate a Simulation Dataset
#'
#' This function generates a dataset for a simulation study, including the response variable, 
#' location information, covariates, and accumulated infected cases.
#'
#' @param n Integer. The number of location points.
#' @param t Integer. The number of time points.
#' @param off Integer. The number of offset observations. Default is 0.
#' @param off_value Numeric. The offset value for observations. Default is 100.
#'
#' @return A list containing:
#' \item{Y}{Matrix. The response variable.}
#' \item{S}{Matrix. The location information.}
#' \item{X}{Matrix. The covariate information.}
#' \item{I}{Vector. The accumulated infected cases.}
#' \item{Mu}{Matrix. The true mean for reponse Y}
#' \item{off}{Vector. The randomly generated outlier value vector}
#'
#' @examples
#' # Generate a simulation dataset with 100 location points, 5 time points, no offset
#' sim_data <- Simulation_Data(n = 100, t = 5)
#' 
#' @examples
#' # Generate a dataset with 100 locations, 5 time points, 
#' # 2 randomly selected outliers with value of 100
#' sim_data <- Simulation_Data(n = 100, t = 5, off = 2, off_value = 100)
#' 
#' 
#' library(dplyr)
#' library(ggplot2)
#' 
#' # plot the simulated mean parameter at each location at time 5
#' tibble(x = sim_data$S[,1], y = sim_data$S[,2], counts = sim_data$Mu[,5], 
#' offset = sim_data$off) %>%
#'   ggplot(aes(x = x, y = y, size = counts, color = offset)) +
#'   geom_point(alpha = 0.6) + 
#'   scale_radius(range = c(2,10)) +
#'   scale_colour_gradient(low = '#7fcdbb', high = '#253494')
#'
#' @import mgcv
#' @importFrom stats rpois runif
#'
#' @export

Simulation_Data=function(n, t, off = 0, off_value = 100){
  # Construct the HorseShoe boundary
  fsb=list(fs.boundary())[[1]]
  tt=1
  while(tt>0){
    # Generate gridded population
    nmax=n*3
    v=runif(nmax)*5-1
    w=runif(nmax)*2-1
    m=fs.test(v,w,b=1)
    names(fsb)=c("v","w")
    ind=inSide(fsb,x=v,y=w) # remove outsiders
    v=v[ind]
    w=w[ind]
    m=m[ind]
    index=sample((1:sum(ind)),size=n,replace=FALSE)
    x1=v[index]
    x2=w[index]
    m=m[index]


    # Simulation the univariate part
    m1=m2=m3=5
    while(abs(mean(m1)) > 0.005 | abs(mean(m2)) > 0.005 | abs(mean(m3)) > 0.005){
      z1=runif(n,0,1)
      z2=runif(n,0,1)
      z3=runif(n,0,1)
      Z=cbind(z1,z2,z3)
 
      m1=1/2*sin(2*pi*z1)-cos(pi*z1)
      m2=4*((z2-0.5)^2-2/3*(0.5)^3)
      m3=z3-0.5
    }

    # Final dataset
    Y = matrix(0, nrow = n, ncol = t)
    # assume we have 1 infected case at beginning (time 0)
    I = matrix(1, nrow = n, ncol = t)
    
    # initial true mean
    eta=m1+m2+m3+0.5*m + 0.1 * log(I[,1]) 
    mu_true=exp(eta)
    
    # randomly adding outliers
    index_off = sample(1:n, off)
    mu_off=mu_true
    mu_off[index_off] = mu_off[index_off] + off_value
    
    # simulate initial counts
    Y_true=rpois(n,mu_true)
    
    Y_off=Y_true
    Y_off[index_off]=rpois(length(index_off), mu_off[index_off])
    
    Y[,1] = Y_off
    
    Mu = matrix(0, nrow = n, ncol = t)
    Mu[,1] = mu_true
    
    # record the true slack value
    Off = log(mu_off/mu_true)
    
    # simulate for the rest time points
    if(t>1){
      for (i in 2:t) {
        
        I[,i] = rowSums(Y[, 1:(i-1), drop = F])+1
        
        eta=m1+m2+m3+0.5*m + 0.1 * log(I[,i])
        mu_true=exp(eta)
        
        Mu[,i] = mu_true
        
        mu_off=mu_true * exp(Off)

        
        Y_true=rpois(n,mu_true)
        
        Y_off=Y_true
        Y_off[index_off]=rpois(length(index_off), mu_off[index_off])
        
        Y[,i] = Y_off

      }
      
    }
    
    data=cbind(Y,x1,x2,z1,z2,z3)

    tt=sum(is.na(data))
  }

  return(list(Y=Y,S=data[,(t+1):(t+2)],X=data[,(t+3):(t+5)], Mu = Mu, I = I, off = Off))
}

