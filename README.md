
# Robust Spatiotemporal Epidemic Modeling with Integrated Adaptive Outlier Detection

<!-- badges: start -->
<!-- badges: end -->

In epidemic modeling, outliers can distort parameter estimation and ultimately lead to misguided public health decisions. Although there are existing robust methods that can mitigate this distortion, the ability to simultaneously detect outliers is equally vital for identifying potential disease hotspots. In this work, we introduce a robust spatiotemporal generalized additive model (RST-GAM) to address this need. We accomplish this with a mean-shift parameter to quantify and adjust for the effects of outliers and rely on adaptive Lasso regularization to model the sparsity of outlying observations. We use univariate polynomial splines and bivariate penalized splines over triangulations to estimate the functional forms and a data-thinning approach for data-adaptive weight construction. We derive a scalable proximal algorithm to estimate model parameters by minimizing a convex negative log-quasi-likelihood function. Our algorithm uses adaptive step-sizes to ensure global convergence of the resulting iterate sequence. This `RSTGAM` package provides codes to implement this robust model and some examples to illustrate its usage.

## Installation

To install the latest stable version of the package (from GitHub):

``` r
devtools::install_github("haomingsj98/RSTGAM")
```

## Example

Please check the [GitHub page](https://github.com/haomingsj98/RSTGAM/blob/main/doc/introduction.pdf) for a detailed example on how to implement the model.

