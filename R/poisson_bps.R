
#' define poisson family
#'
#' Define some function for Poisson Family
#'
#' @param link: name of link function, default as log
#' @return A group of useful functions
#' @keywords internal
#' 
poisson_fam=function(link='log'){
  if(link=='log') {
    linkfun=function(mu) log(mu)
    linkinv=function(eta) exp(eta)
    mu.eta=function(eta) exp(eta)
    gprime=function(mu) 1/mu
  }
  variance=function(mu) mu
  initialize=function(y){
    if (any(y < 0))
      stop("negative values not allowed for the 'quasiPoisson' family")
    y+0.1
  }
  structure(list(family = "poisson_bps", link=link,
                 linkfun = linkfun, linkinv = linkinv, variance = variance,
                 mu.eta = mu.eta, gprime = gprime,
                 initialize = initialize), class = "family")
}
