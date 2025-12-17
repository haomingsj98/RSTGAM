
#' define poisson family
#'
#' Define some function for Negative Binomial Family
#'
#' @param link: name of link function, default as log
#' @return A group of useful functions
#' @keywords internal
#' 
nb_fam=function(link='log'){
  if(link=='log') {
    linkfun=function(mu) log(mu)
    linkinv=function(eta) exp(eta)
    mu.eta=function(eta) exp(eta)
    gprime=function(mu) 1/mu
  }
  variance=function(mu, theta) mu+mu^2/theta
  initialize=function(y){
    if (any(y < 0))
      stop("negative values not allowed for the 'Negative Binomial' family")
    y+0.1
  }
  structure(list(family = "nb_bps", link=link,
                 linkfun = linkfun, linkinv = linkinv, variance = variance,
                 mu.eta = mu.eta, gprime = gprime,
                 initialize = initialize), class = "family")
}
