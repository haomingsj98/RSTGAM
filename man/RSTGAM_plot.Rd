% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rst_gam_plot.R
\name{RSTGAM_plot}
\alias{RSTGAM_plot}
\title{Summarize and Prepare RST_GAM Results for Plotting}
\usage{
RSTGAM_plot(model, uni_sample = 1000, S)
}
\arguments{
\item{model}{list. The output from the RST_GAM model.}

\item{uni_sample}{integer. The number of evenly spaced points used to evaluate univariate components, default is 1000.}

\item{S}{matrix. The location matrix. Each row indicates the location of one point.}
}
\value{
A list containing the following elements:
\item{summary1}{A tibble summarizing estimated slack values and bivariate components.}
\item{X0}{A matrix of covariate values at equally spaced points, with the number of points determined by uni_sample.}
\item{summary2}{A matrix for each estimated univariate component, evaluated at X0.}
}
\description{
This function processes the output of the robust gGAM estimation model (RST_GAM) to generate summaries and prepare data for visualization.
}
\examples{

# Refer to the detailed example in the documentation for the RST_GAM function.
# This function is typically used after running RST_GAM to process its output
# and prepare data for visualization.


}
\seealso{
RST_GAM
}
