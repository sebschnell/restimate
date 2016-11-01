// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// est_var_mb
double est_var_mb(const arma::mat& cov, const arma::mat& X);
RcppExport SEXP restimate_est_var_mb(SEXP covSEXP, SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type cov(covSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(est_var_mb(cov, X));
    return rcpp_result_gen;
END_RCPP
}
// res_ss_tm
double res_ss_tm(const arma::mat& X, const arma::mat& Z, const arma::colvec beta, const arma::colvec alpha, const arma::mat vcov_beta);
RcppExport SEXP restimate_res_ss_tm(SEXP XSEXP, SEXP ZSEXP, SEXP betaSEXP, SEXP alphaSEXP, SEXP vcov_betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type vcov_beta(vcov_betaSEXP);
    rcpp_result_gen = Rcpp::wrap(res_ss_tm(X, Z, beta, alpha, vcov_beta));
    return rcpp_result_gen;
END_RCPP
}
// vcov_tm
arma::mat vcov_tm(const arma::mat& Z, const double ss_res);
RcppExport SEXP restimate_vcov_tm(SEXP ZSEXP, SEXP ss_resSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const double >::type ss_res(ss_resSEXP);
    rcpp_result_gen = Rcpp::wrap(vcov_tm(Z, ss_res));
    return rcpp_result_gen;
END_RCPP
}
// res_ss_hc_tm
arma::mat res_ss_hc_tm(const arma::mat& X, const arma::mat& Z, const arma::colvec beta, const arma::colvec alpha, const arma::mat vcov_beta_hc);
RcppExport SEXP restimate_res_ss_hc_tm(SEXP XSEXP, SEXP ZSEXP, SEXP betaSEXP, SEXP alphaSEXP, SEXP vcov_beta_hcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type vcov_beta_hc(vcov_beta_hcSEXP);
    rcpp_result_gen = Rcpp::wrap(res_ss_hc_tm(X, Z, beta, alpha, vcov_beta_hc));
    return rcpp_result_gen;
END_RCPP
}
// vcov_tm_hc
arma::mat vcov_tm_hc(const arma::mat& res_sq, const arma::mat& X);
RcppExport SEXP restimate_vcov_tm_hc(SEXP res_sqSEXP, SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type res_sq(res_sqSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(vcov_tm_hc(res_sq, X));
    return rcpp_result_gen;
END_RCPP
}