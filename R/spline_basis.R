
#' Generate constant spline bases given a vector
#'
#' This routine generates constant spline bases given a vector
#'
#' @param x Covariates.
#' @param knots An optional vector specifying the knots to be used in constructing spline bases.
#'
#' @return A list containing the following elements:
#' \item{Bx0}{Non-centered spline bases for covariates.}
#' \item{Bx}{Centered spline bases for covariates.}
#' \item{BxMean}{Means of spline bases subtracted.}
#' \item{knots}{Knots used in constructing spline bases.}
#'
#'
#' @keywords internal
#' 
cbs=function(x,knots,Boundary.knots=range(x)){
	n=length(x)
	N=length(knots)
	Ix=matrix(0L,nrow=n,ncol=N+1L)
	for (i in 1:n) {
		if (x[i]<knots[1L]) {Ix[i,1]=1L
			} else if (x[i]>=knots[N]) {Ix[i,N+1L]=1L
				} else {
					kl=max(which(x[i]>=knots))
					Ix[i,kl+1]=1L
				}
	}
	return(Ix)
}


#' Generate Basis Functions for Univariate Spline
#'
#' This function generates basis functions for univariate splines.
#'
#' @param x Covariates.
#' @param N Number of interior knots in generating spline matrix.
#' @param q Degree of polynomial spline. Default is 3.
#' @param KnotsLocation A character string naming the way for knots locations. Default is \code{"quantile"}. The only alternative is \code{"uniform"}.
#' @param knots An optional vector specifying the knots to be used in constructing spline bases.
#'
#' @return A list containing the following elements:
#' \item{Bx0}{Non-centered spline bases for covariates.}
#' \item{Bx}{Centered spline bases for covariates.}
#' \item{BxMean}{Means of spline bases subtracted.}
#' \item{knots}{Knots used in constructing spline bases.}
#'
#'
#' @importFrom splines bs
#' @importFrom stats quantile
#' 
#' 
#' 
#' @keywords internal


Basis_generator=function(x,N,q=3,KnotsLocation="quantile",knots=NULL){
	# Case I: x is a vector
	if (length(dim(x))==0) {
		# generate knots by supplied location
		if (is.null(knots)) {
			# knots generated in sample quantile
			if (KnotsLocation=="quantile")
				knots1=unique(quantile(x,probs=seq(0,1,by=1/(N+1))))
			# knots generated in uniform points
			if (KnotsLocation=="uniform")
				knots1=seq(min(x),max(x),by=(max(x)-min(x))/(N+1))
		} else knots1=knots
		N=length(knots1)-2
		# generate constant / polynomial B-spline basis
		if (q==0) {
			Bx0=cbs(x,knots1[-c(1,N+2)],Boundary.knots=knots1[c(1,N+2)])
		} else Bx0=bs(x,knots=knots1[-c(1,N+2)],degree=q,intercept=F,
					Boundary.knots=knots1[c(1,N+2)])
		Knots=knots1
	} else { # Case II: x is a matrix
		d.x=ncol(x)
		block.size=floor(sqrt(d.x))
		nblock=ceiling(d.x/block.size)
		Bx0=NULL
		Knots=NULL

		# This routine reduces the memory needed in constructing spline
		#	bases and improve the computational efficiency
		for(nj in 1:nblock) {
			Bxnj=NULL
			Knotsj=NULL
			if(nj<nblock) block.ind=(block.size*(nj-1)+1):(block.size*nj)
			if(nj==nblock) block.ind=(block.size*(nj-1)+1):d.x
			for (j in block.ind) {
				# generate knots by supplied location
				if (is.null(knots)) {
					# knots generated in sample quantile
					if (KnotsLocation=="quantile") knots1=quantile(
						x[,j],probs=seq(0,1,by=1/(N+1)))
					# knots generated in uniform points
					if (KnotsLocation=="uniform") knots1=seq(min(x[,j]),
						max(x[,j]),by=(max(x[,j])-min(x[,j]))/(N+1))
				} else knots1=knots[(N+2)*(j-1)+1:(N+2)]

				# generate constant / polynomial B-spline basis
				if (q==0) {
					bx0=cbs(x[,j],knots1[-c(1,N+2)],
						Boundary.knots=knots1[c(1,N+2)])
				} else bx0=bs(x[,j],knots=knots1[-c(1,N+2)],degree=q,
							intercept=F,Boundary.knots=knots1[c(1,N+2)])
				Bxnj=cbind(Bxnj,bx0)
				Knotsj=cbind(Knotsj,knots1)
			}
			Bx0=cbind(Bx0,Bxnj)
			Knots=cbind(Knots,Knotsj)
		}
	}
	BxMean=colMeans(Bx0)
	Bx=sweep(Bx0,2,BxMean,"-")
	list(Bx0=Bx0,Bx=Bx,BxMean=BxMean,knots=Knots)
}
