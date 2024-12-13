% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rst_gam.R
\name{RST_GAM}
\alias{RST_GAM}
\title{Fit the Global Estimation Model (RST_GAM)}
\usage{
RST_GAM(
  Y,
  X,
  X_t = NULL,
  S,
  Tr,
  V,
  d = 2,
  r = 1,
  rho = 2,
  knots = NULL,
  N = 4,
  initial_x = NULL,
  k = 3,
  lambda1,
  lambda2
)
}
\arguments{
\item{Y}{matrix. The response variable.}

\item{X}{matrix. The covariate matrix. Each column represents one covariate.}

\item{X_t}{matrix. A matrix of time-changing covariates, e.g., the accumulated infected cases and human mobility, each column represents a covariate, default as null.}

\item{S}{matrix. The location matrix. Each row indicates the location of one point.}

\item{Tr}{matrix. The triangulation matrix. Each row indicates one triangle.}

\item{V}{matrix. The vertices matrix of the triangulation.}

\item{d}{integer. The degree of the bivariate spline, default as 2.}

\item{r}{integer. The smoothness of the bivariate spline, default as 1.}

\item{rho}{integer. The degree of the univariate spline, default as 2.}

\item{knots}{matrix. The knots matrix of the univariate spline. Usually, it is a (N+2) x p matrix. The r-th column is the knots of r-th covariate. The first element and the last element are the boundary knots.}

\item{N}{integer. The number of interior knots, default as 4.}

\item{initial_x}{numeric vector. The initial value for the parameter, default is 0.}

\item{k}{The number data thinning folds, default as 3.}

\item{lambda1}{numeric. The roughness penalty parameter for bivariate spline estimation.}

\item{lambda2}{numeric. The L1 penalty term for the slack variable.}
}
\value{
A list containing the following elements:
\item{theta_hat}{coefficients estimated for univariate functions.}
\item{gamma_hat}{coefficients estimated for bivariate function.}
\item{xi_hat}{the estimated slack parameter.}
\item{mu_hat}{the estimated mean parameter.}
\item{df}{degree of freedom of the model.}
\item{r_lambda}{the selected roughness penalization parameter.}
\item{l_lambda}{the selected L1 penalization parameter.}
\item{Ind}{the index of the inside observation.}
\item{rho}{the degree of univariate spline, same as input argument.}
\item{d}{the degree of bivariate spline, same as input argument.}
\item{r}{the smoothness of bivariate spline, same as input argument.}
\item{knots}{the knots matrix of univariate spline, same as input argument.}
\item{N}{the number of interior knots, same as input argument.}
\item{sigma_2}{estimated sigma^2.}
\item{Y}{re-ordered response variable.}
\item{X}{re-ordered covariate matrix.}
\item{X_t}{re-ordered time-changing covariate matrix.}
\item{U}{univariate spline basis matrix.}
\item{BQ}{bivariate spline basis matrix.}
\item{P}{energy matrix for bivariate spline.}
}
\description{
This function fits the robust ggam estimation model to account for potentail outliers using bivariate and univariate splines.
}
\examples{
# Load the example dataset
data(example_dataset)
library(BPST)

# the lambda list for roughness
lambda1=0.1
# the lambda list for l1
lambda2=0.02

# Fit the model using RST_GAM with 0 data thinning fold
result <- RST_GAM(
  Y = example_dataset$Y,
  X = example_dataset$X,
  S = example_dataset$S,
  Tr = example_dataset$Tr,
  V = example_dataset$V,
  lambda1 = lambda1,
  lambda2 = lambda2,
  k = 0
)


library(dplyr)
library(ggplot2)
library(mgcv)

# obtain the summarized results for visualization
summary_plot <- RSTGAM_plot(result, S = example_dataset$S)

# Plot the identified outliers with corresponding estimated slack values
summary_plot$summary1 \%>\% ggplot(aes(x=X, y=Y, color=Slack)) +
  geom_point(alpha = 0.6) + 
  ggtitle("Outlier Signal Estimation") + 
  scale_colour_gradient(low = "white", high = "#bd0026", name = "Outlier Signal")
  

# Plot the estimated bivariate component
summary_plot$summary1 \%>\% ggplot(aes(x=X, y=Y, color=Beta)) +
  geom_point(alpha = 0.6) + 
  ggtitle("Bivariate Component Estimation") + 
  scale_colour_gradient(low = "#7fcdbb", high = "#bd0026", 
  name = "Bivariate Function Value")

# true beta
m=fs.test(example_dataset$S[,1][result$Ind], 
          example_dataset$S[,2][result$Ind], b=1)
         
true_beta = 0.5*m

# Compute the bivariate MISE
mean((as.matrix(summary_plot$summary1$Beta) - true_beta)^2)


# Define the true coefficient functions

beta.1 = function(x0) 1/2 * sin(2 * pi * x0) - cos(pi * x0)
beta.2 = function(x0) 4 * ((x0 - 0.5)^2 - 2 / 3 * (0.5)^3)
beta.3 = function(x0) x0 - 0.5


# the estimation from model
X0 = summary_plot$X0
mhat = summary_plot$summary2


# Calculate true values for given data
beta0.1 = beta.1(X0[, 1])
beta0.2 = beta.2(X0[, 2])
beta0.3 = beta.3(X0[, 3])


# Set up plotting margins and layout
par(mar = c(5, 9, 4, 2))
par(mfrow = c(1, 3))  # Arrange plots in a 2x2 grid

# Plot for X1
plot(X0[, 1], mhat[, 1], xlab = 'X1', type = 'l', lwd = 3, ylab = '', 
     cex.axis = 1.7, cex.lab = 2, cex = 1.3, pch = 16, col = 'blue')
lines(X0[, 1], beta0.1, lwd = 3, col = 'red')
legend(0, 1, legend = c('True', 'Estimates with Slack'), 
       fill = c('red', 'blue'), cex = 1.6)

# Plot for X2
plot(X0[, 2], mhat[, 2], xlab = 'X2', type = 'l', lwd = 3, ylab = '', 
     cex.axis = 1.7, cex.lab = 2, cex = 1.3, pch = 16, col = 'blue')
lines(X0[, 2], beta0.2, lwd = 3, col = 'red')

# Plot for X3
plot(X0[, 3], mhat[, 3], xlab = 'X3', type = 'l', lwd = 3, ylab = '', 
     cex.axis = 1.7, cex.lab = 2, cex = 1.3, pch = 16, col = 'blue')
lines(X0[, 3], beta0.3, lwd = 3, col = 'red')


# the MISE for each univarite component
sum((beta0.1-mhat[,1])^2)/1000
sum((beta0.2-mhat[,2])^2)/1000
sum((beta0.3-mhat[,3])^2)/1000


}
