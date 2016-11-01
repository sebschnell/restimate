#' White's heteroskedasticity consistent covariance matrix estimator
#'
#' @param X Design matrix
#' @param res Model residuals
#'
#' @return Estimated parameter variance-co-variance matrix
#'
#' @references White, H. (1980). "A Heteroskedasticity-Consistent Covariance
#'   Matrix Estimator and a Direct Test for Heteroskedasticity." Econometrica
#'   48(4): 817-838.
#' @export
vcov_white <- function(X, res) {
  n = nrow(X);
  df = n - ncol(X);

  V_n = matrix(0, ncol(X), ncol(X));
  for (i in seq.int(nrow(X))) {
    V_n = V_n + n*res[i]^2/df*X[i,]%*%t(X[i,]);
  }

  return(solve(t(X)%*%X)%*%V_n%*%solve(t(X)%*%X));
}
