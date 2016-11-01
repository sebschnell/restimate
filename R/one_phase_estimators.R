#'Estimate the population total using the Horvitz-Thompson estimator
#'
#'This function estimates the population total of a study variable using the
#'Horvitz-Thompson estimator. The function also allows for domain estimation.
#'Two cases are distinguished: domian size is known and domain size is unknown.
#'
#'@param y Study variable
#'@param prob Inlcusion probabilities of the design
#'@param Nd Domain size. If provided a more precise estimator for the domain
#'  total is used.
#'@return A named vector of three elements containing the estimated total, the
#'  sample size, and the number of non-zero sample units
#'
#'@references See Särndal et al. (1992), eq. 10.3.3 and 10.3.6
#'@export
est_tot_ht <- function(y, prob, Nd = NULL) {
  n <- length(y);
  n_non_zero <- length(y[y != 0]);

  if (is.null(Nd) == TRUE) {
    t_hat <- sum(y/prob);
  } else {
    if (length(prob) == 1) {
      prob <- rep(prob, length(y));
    }
    y_hat <- sum(y/prob)/sum(1/prob);
    t_hat <- Nd*y_hat;
  }
  return(c(t_hat = t_hat, n = n, n_non_zero = n_non_zero));
}

#'Estimate the population total using the generalized regression estimator
#'(GREG)
#'
#'This function estimates the population total of a study variable using the
#'GREG estimator. Domain estimation is possible for the two cases that domain
#'size is either known or unknown. If known, the appropriate variable needs to
#'be specified.
#'
#'@param t_pred Population total of model predicitions
#'@param res Residuals between observed and predicted values in the sample
#'@param prob Inclusion probabilities of the design
#'@param Nd Domain size. If provided a more precise estimator for the domain
#'  total is used.
#'@return A named vector of three elements containing the estimated total, the
#'  sample size, and the number of non-zero sample units
#'
#'@references See Särndal et al. (1992), eq. 10.5.4 and 10.5.5
#'@export
est_tot_ma <- function(t_pred, res, prob, Nd = NULL) {
  n <- length(res);
  n_non_zero <- length(res[res != 0]);

  if (is.null(Nd) == TRUE) {
    t_hat = t_pred + sum(res/prob);
  } else {
    Nd_hat <- n/prob;
    t_hat <- t_pred + Nd/Nd_hat * sum(res/prob);
  }
  return(c(t_hat = t_hat, n = n, n_non_zero = n_non_zero));
}

#'Estimate the variance of the estimated population total under simple random
#'sampling without replacement
#'
#'This function estimates the variance of the estimated population total of a
#'study variable under simple random sampling without replacement. The function
#'also allows for domain estimation if domain size is either known or unknown.
#'
#'@param y Study variable
#'@param n Sample size. May be larger than the length of y if the aim is to
#'  estimate the variance for a certain domain.
#'@param N Size of the population
#'@param Nd Size of the domain. If provided a more precise estimator for the
#'  variance of the total is used.
#'@return A named vector of one element containing the variance estimate.
#'
#'@references Särndal et al. (1992), eq. 10.3.10 and 10.3.12
#'@export
est_var_srswor <- function(y, n, N, Nd = NULL) {
  if (is.null(Nd) == TRUE) {
    p <- length(y)/n;
    y_bar <- mean(y);
    v_hat <- N^2*(1/n - 1/N)*p*(var(y) + (1 - p)*y_bar^2);
  } else {
    nd <- length(y);
    Nd_hat <- N/n*nd;
    v_hat <- Nd^2*(1/nd - 1/Nd_hat)*var(y);
  }
  return(c(v_hat = v_hat));
}

#'Estimate the variance of the estimated population total under simple random
#'sampling without replacement followign generalised regression estimation
#'
#'This function estimates the variance of the estimated population total of a
#'study variable under simple random sampling without replacement when
#'generalised linear regression was used. The function considers domain
#'estimation.
#'
#'@param res Residuals between observed and predicted values in the sample
#'@param N Size of the population. For domain variance provide the domain size.
#'
#'@return A named vector of one element containing the variance estimate.
#'
#'@references Särndal and Hidiroglou (1989)
#'@export
est_var_srswor_ma <- function(res, N) {
  nd <- length(res);
  v_hat <- N^2*(1/nd - 1/N)*var(res);
  return(c(v_hat = v_hat));
}