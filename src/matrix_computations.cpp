#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Function to compute eta
// [[Rcpp::export]]
arma::vec compute_eta(const arma::mat& X, const arma::vec& theta, const arma::mat& BQ2, const arma::vec& gamma, const arma::vec& xi, int t) {
    arma::vec eta = X * theta + BQ2 * gamma;
    eta += arma::kron(arma::ones(t, 1), xi);
    return eta;
}

// Function to compute gradient_gamma
// [[Rcpp::export]]
arma::vec compute_gradient_gamma(const arma::mat& BQ2, const arma::vec& temp3, const arma::mat& P, const arma::vec& gamma, double lambda1) {
    arma::vec gradient_gamma = BQ2.t() * temp3 + lambda1 * P * gamma;
    return gradient_gamma;
}

// Function to compute t(X) %*% temp3
// [[Rcpp::export]]
arma::vec compute_gradient_theta(const arma::mat& X, const arma::vec& temp3) {
    return X.t() * temp3;
}
