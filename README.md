
# Robust Spatiotemporal Epidemic Esimation Model with Outlier Detection

<!-- badges: start -->
<!-- badges: end -->

In epidemic modeling, the presence of outliers can distort parameter estimation and lead to misguided public health decisions. Although many existing methods achieve robustness, the ability to simultaneously detect outliers is equally vital for identifying potential disease hotspots. Here, we introduce a robust spatiotemporal generalized additive model (RST-GAM) to address this challenge by integrating a mean-shift parameter to quantify and mitigate the effects of outliers, combined with an appropriately designed adaptive Lasso regularization. Univariate polynomial splines and bivariate penalized splines over triangulations are used to estimate the functional forms, and a data-thinning approach is employed to enable adaptive weight construction. A proximal algorithm with adaptive step sizes is proposed to ensure computational efficiency for the proposed non-globally Lipschitz convex problems. This `RSTGAM` package provides codes to implement this robust model and some examples to illustrate its usage.

## Installation

To install the latest stable version of the package (from GitHub):

``` r
devtools::install_github("haomingsj98/RSTGAM")
```

## Example

Please check the [GitHub page](https://github.com/haomingsj98/RSTGAM) for a detailed example on how to implement the model.

