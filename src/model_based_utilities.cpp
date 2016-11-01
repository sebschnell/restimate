#include <RcppArmadillo.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

//' Model-based variance estimator for an estimated population total
//'
//' @param cov Model parameter covariance matrix
//' @param X Design matrix
//' @return The variance estimate
//'
//' @export
// [[Rcpp::export]]
double est_var_mb(const arma::mat &cov, const arma::mat &X) {
  arma::mat one_vec(X.n_rows, 1);
  one_vec.ones();
  arma::mat r = one_vec.t()*X*cov*X.t()*one_vec;
  return r(0, 0);
}

//' Return the sum of squared estimated residuals when two OLS models are
//' applied in model-based inference
//'
//' @param X Design matrix of the first model for the first-phase sample
//' @param Z Design matrix of the second model for the first-phase sample
//' @param beta Parameter estimates of the first model (column vector)
//' @param alpha Parameter estimates of the second model (column vector)
//' @param vcov_beta Variance/co-variance matrix of first-model parameters
//'   estimates
//' @return The residual sum of squares
//'
//' @export
// [[Rcpp::export]]
double res_ss_tm(const arma::mat &X,
                 const arma::mat &Z,
                 const arma::colvec beta,
                 const arma::colvec alpha,
                 const arma::mat vcov_beta) {
  arma::mat res = Z*alpha - X*beta;
  arma::mat r = res.t()*res + est_var_mb(vcov_beta, X);
  return r(0, 0);
}

//' Estimated model parameter covariance matrix for the top-level model in
//' hierarchical model-based inference setting when two models are applied
//'
//' @param Z Design matrix of the second model for the first-phase sample
//' @param res Sum of squared residuals for the top-level model
//' @return Estimated covariance matrix for model parameter estimates
//'
//' @export
// [[Rcpp::export]]
arma::mat vcov_tm(const arma::mat &Z,
               const double ss_res) {
  double mss = ss_res/(Z.n_rows - Z.n_cols - 1);
  return mss*(Z.t()*Z).i();
}

//' Return the squared estimated residuals when two OLS models are
//' applied in model-based inference and heteroscedasticity is present
//'
//' @param X Design matrix of the first model for the first-phase sample
//' @param Z Design matrix of the second model for the first-phase sample
//' @param beta Parameter estimates of the first model (column vector)
//' @param alpha Parameter estimates of the second model (column vector)
//' @param vcov_beta_hc Variance/co-variance matrix of first-model parameters
//'   estimates
//' @return The residual sum of squares
//'
//' @export
// [[Rcpp::export]]
arma::mat res_ss_hc_tm(const arma::mat &X,
                       const arma::mat &Z,
                       const arma::colvec beta,
                       const arma::colvec alpha,
                       const arma::mat vcov_beta_hc) {
  arma::colvec one_vec(X.n_rows, arma::fill::ones);
  arma::mat res = Z*alpha - X*beta;
  arma::mat r = square(res) + X*vcov_beta_hc*X.t()*one_vec;
  return r;
}

//' Estimated model parameter covariance matrix for the top-level model in
//' hierarchical model-based inference setting when two models are applied
//'
//' @param Z Design matrix of the second model for the first-phase sample
//' @param res Sum of squared residuals for the top-level model
//'
//' @details The heteroskedasticity-consistent version
//'
//' @return Estimated covariance matrix for model parameter estimates
//'
//' @export
// [[Rcpp::export]]
arma::mat vcov_tm_hc(const arma::mat &res_sq,
                     const arma::mat &X) {
  const arma::mat cross_inv = (X.t()*X).i();
  arma::mat sum(X.n_cols, X.n_cols, arma::fill::zeros);
  for (unsigned i = 0; i < X.n_rows; ++i) {
    sum += res_sq(i)*X.row(i).t()*X.row(i);
  }
  return cross_inv*sum*cross_inv;
}
