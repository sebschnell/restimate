#' Internal helper function to prepare data from a two-stage sample for
#' estimation by other functions
#'
#' @param dt_samp A data.table object containing information about the PSUs and
#'   SSUs selected in the two sampling stages.
#' @param key_var Name of the variable identifying the PSUs.
#' @param y Quoted name of the variable storing the observed values of the study
#'   variable
#' @param res Quoted name of the variable under which residuals between observed
#'   and predicted values are stored. Optional in case of model-assisted
#'   estimation
#'
#' @return A data.table object containing PSU summaries needed for estimation by
#'   other functions.
#' @import data.table
get_psu_stats <- function(dt_samp, id_psu, y, res = NULL) {
  dt_psu <- dt_samp[,
                    list(sum_y = sum(eval(y), na.rm = TRUE),
                         var_y = var(eval(y), na.rm = TRUE)),
                    keyby = id_psu];
  dt_psu[is.na(var_y) == TRUE, var_y := 0];

  if (is.null(res) == FALSE) {
    dt_psu <- dt_psu[dt_samp[,
                             list(sum_res = sum(eval(res), na.rm = TRUE),
                                  var_res = var(eval(res), na.rm = TRUE)),
                             keyby = id_psu]];
    dt_psu[is.na(var_res) == TRUE, var_res := 0];
  }
  return(dt_psu);
}

#' Succesive difference estimator to estimate variance of a systematic sample
#'
#' @param y Observations
#' @param N Population size
#' @param n Sample size
#' @param w Weights in case systematic sampling of clusters of unequal size is
#'   applied. See details for how they are calculated.
#'
#' @details For the unweighted version, see the \eqn{v_4}-estimator in Wolter
#'   (2003, chapter 8.3). For the weighted version and the case of cluster of
#'   unequal size, see Nelson et al. (2008, Equation 3). The cluster weights are
#'   defined as \deqn{W_i = N_i(\sum_{i=1}^{n-1}N_i),} where \eqn{N_i} denotes
#'   the size of a cluster and \eqn{n} is the number of clusters selected into
#'   the sample.
#' @return An atomic vector with the variance estimate.
#' @references Nelson, R., Easset, T. Gobakken, G. Stahl and T. G. Gregoire
#' (2008). "Regional Forest Inventory using an Airborne Profiling LiDAR."
#' Journal of Forest Planning 13: 287-294.
#'
#' Wolter, K. (2003). Introduction to Variance Estimation, Springer.
#'
#' @export
est_var_succ_diff <- function(y, N, n, w = NULL) {
  y_b <- data.table::shift(y, n = 1, type = 'lag');
  y_a <- data.table::shift(y, n = 1, type = 'lead');
  y_succ <- (y_b - 2*y + y_a)^2;

  if (is.null(w) == TRUE) {
    v_hat <- (1 - n/N)*n/(96*(n - 2))*sum(y_succ, na.rm = TRUE);
  } else {
    w_b <- data.table::shift(w, n = 1, type = 'lag');
    w_a <- data.table::shift(w, n = 1, type = 'lead');
    w_succ <- (w_b + 2*w + w_a)^2;
    v_hat <- (1 - n/N)*n/(96*(n - 2))*sum(w_succ*y_succ, na.rm = TRUE);
  }
  return(v_hat);
}

#'Model assisted estimators for two-stage systematic sampling
#'
#'The population is separated into a set of contiguous PSUs (clusters/strips),
#'which are sampled in the first stage. SSUs are selected systematically across
#'the selected PSUs in the second stage. The estimators correspond to case C in
#'Särndal et al. (1992) adapted to systematic sampling, where auxiliaries are
#'avalaible for SSUs in the selected PSUs. Modelling is done at the element
#'level.
#'
#'@param dt_sample A data.table object containing information about the PSUs and
#'  SSUs selected in the two sampling stages.
#'@param id_psu Name of the variable identifying the PSUs. Must be present in
#'  both data.tables.
#'@param y Quoted name of the variable storing the observed values of the study
#'  variable.
#'@param res Quoted name of the variable stroing the residuals between observed
#'  and predicted values.
#'@param t_pred Quoted name of the variable storing predicted PSU totals
#'@param a_psu PSU intervall, every a-th PSU is selected into the sample.
#'@param a_ssu SSU intervall, every a-th SSU is selected into the sample.
#'@param NI Number of PSUs
#'@param nI PSU sample size
#'@param N_i Optional. PSU size in case PSUs are of different size. Needed for
#'  the successive differences variance estimator.
#'
#'@details See Sarndal et al. (1992), equations 8.9.7 for the estimation of
#'  the population total, slightly adapted to case of systematic sampling. See
#'  Ene et al. (2016, equation 13) for variance estimation using successive
#'  differences.
#'
#'@return A named list containing estimates for the total and its variance
#'  estimate
#'
#'@references
#'  Ene, L. T., E. Næsset and T. Gobakken (2016). "Simulation-based assessment
#'  of sampling strategies for large-area biomass estimation using wall-to-wall
#'  and partial coverage airborne laser scanning surveys." Remote Sensing of
#'  Environment 176: 328-340.
#'
#'  Särndal, C. E., B. Swensson and J. Wretman (1992). Model assisted survey
#'  sampling. New York, Springer-Verlag.
#'
#'@import data.table
#'@export
est_ts_ma_sys <- function(dt_samp, id_psu = 'id', y = quote(y),
                          res = quote(res), t_pred = quote(t_pred),
                          a_psu, a_ssu, NI, nI, N_i = quote(N_i)) {

  dt_psu <- get_psu_stats(dt_samp = dt_samp, y = y, id_psu = id_psu, res = res);
  dt_psu <- dt_psu[dt_samp[, list(t_pred = max(eval(t_pred)),
                                  N_i = max(eval(N_i))),
                           keyby = id_psu]];

  dt_psu[, ':='(t_res = a_ssu*sum_res,
                w_i = N_i/sum(N_i))];
  dt_psu[is.na(t_res) == TRUE, t_res := 0];
  dt_psu[, t_hat_ma := t_pred + t_res];
  dt_psu[, mean_hat_ma := t_hat_ma/N_i];

  # Equation 8.9.7 in Särndal et al. (1992), slightly adapted
  t_hat <- a_psu * dt_psu[, sum(t_hat_ma)];

  # Equation 13 in Ene et al. (2016)
  v_hat <- dt_psu[, est_var_succ_diff(y = mean_hat_ma, N = NI, n = nI, w = w_i)];

  return(list(t_hat = t_hat, v_hat = v_hat));
}

#'The population is separated into a set of contiguous PSUs (clusters/strips),
#'which are sampled in the first stage. SSUs are selected systematically across
#'the selected PSUs in the second stage. Adapted to systematic sampling.
#'
#'@param dt_sample A data.table object containing information about the PSUs and
#'  SSUs selected in the two sampling stages.
#'@param id_psu Name of the variable identifying the PSUs. Must be present in
#'  both data.tables.
#'@param y Quoted name of the variable storing the observed values of the study
#'  variable.
#'@param a_psu PSU intervall, every a-th PSU is selected into the sample.
#'@param a_ssu SSU intervall, every a-th SSU is selected into the sample.
#'@param NI Number of PSUs
#'@param nI PSU sample size
#'@param N_i Optional. PSU size in case PSUs are of different size. Needed for
#'  the successive differences variance estimator.
#'
#'@details See Särndal et al. (1992), equations 4.3.21 for the estimation of
#'  the population total, slightly adapted to case of systematic sampling. See
#'  Ene et al. (2016, equation 13) for variance estimation using successive
#'  differences.
#'
#'@return A named list containing estimates for the total and its variance
#'  estimate
#'
#'@references
#'  Ene, L. T., E. Næsset and T. Gobakken (2016). "Simulation-based assessment
#'  of sampling strategies for large-area biomass estimation using wall-to-wall
#'  and partial coverage airborne laser scanning surveys." Remote Sensing of
#'  Environment 176: 328-340.
#'
#'  Särndal, C. E., B. Swensson and J. Wretman (1992). Model assisted survey
#'  sampling. New York, Springer-Verlag.
#'
#'@import data.table
#'@export
est_ts_sys <- function(dt_samp, id_psu = 'id', y = quote(y),
                       a_psu, a_ssu, NI, nI, N_i = quote(N_i)) {

  dt_psu <- get_psu_stats(dt_samp = dt_samp, y = y, id_psu = id_psu);
  dt_psu <- dt_psu[dt_samp[, list(t_pred = max(eval(t_pred)),
                                  N_i = max(eval(N_i))),
                           keyby = id_psu]];

  dt_psu[, ':='(t_hat = a_ssu*sum_y,
                w_i = eval(N_i)/sum(eval(N_i)))];
  dt_psu[, mean_t_hat := t_hat/N_i];

  # Equation 8.9.7 in Särndal et al. (1992), slightly adapted
  t_hat <- a_psu * dt_psu[, sum(t_hat)];

  # Equation 13 in Ene et al. (2016)
  v_hat <- dt_psu[, est_var_succ_diff(y = mean_t_hat, N = NI, n = nI, w = w_i)];

  return(list(t_hat = t_hat, v_hat = v_hat));
}

#'Model-assisted estimators for two-stage simple random sampling
#'
#'The population is separated into a set of contiguous PSUs (clusters/strips),
#'which are sampled in the first stage. SSUs are selected from the selected PSUs
#'in the second stage independently. The estimators correspond to case C in
#'Särndal et al. (1992), where auxiliaries are avalaible for SSUs in the
#'selected PSUs. Modelling is done at the element level.
#'
#'@param dt_sample A data.table object containing information about the PSUs and
#'  SSUs selected in the two sampling stages.
#'@param key_var Name of the variable identifying the PSUs. Must be present in
#'  both data.tables.
#'@param y Quoted name of the variable storing the observed values of the study
#'  variable.
#'@param res Quoted name of the variable stroing the residuals between observed
#'  and predicted values.
#'@param t_pred Quoted name of the variable storing predicted PSU totals
#'@param N_i quoted name of the variable storing the size of the PSUs
#'@param n_i quoted name of the variable storing the number of selected SSUs in
#'  each PSUs
#'@param NI Number of PSUs
#'@param nI PSU sample size
#'
#'@details See Särndal et al. (1992), equation 8.9.7 for the estimation of the
#'  population total and Equation 8.9.27 for traditional variance estimation.
#'  Under some circumstances, the variance estimator may deliver negative
#'  estimates. Therefore, two alternative variance estimators are returned as
#'  well. The first is the so called prediction-based, model-assisted variance
#'  estimator and uses model-assisted estimates of PSU totals instead of
#'  Horvitz-Thompson PSU total estimates. This variance estimator should be more
#'  stable (Saarela et al. 2016). The second ignores the variation that comes
#'  from estimating PSU totals (Särndal et al., p. 154, Equation 4.6.2).
#'
#'@return A named list containing estimates for the total and its variance
#'  estimate
#'
#'@references Särndal, C. E., B. Swensson and J. Wretman (1992). Model assisted
#'  survey sampling. New York, Springer-Verlag.
#'
#'  Saarela, S.; Andersen, H.-E.; Grafstrom, A.; Schnell, S.; Gobakken, T.;
#'  Næsset, E.; Nelson, R. F.; McRoberts, R. E.; Gregoire, T. & Ståhl, G. A new
#'  prediction-based variance estimator for two-stage model-assisted surveys of
#'  forest resources. Remote Sensing of Environment, 2016
#'
#'@import data.table
#'@export
est_ts_ma_srswor <- function(dt_samp, NI, nI,
                             id_psu = 'id',
                             y = quote(y),
                             N_i = quote(N_i),
                             n_i = quote(n_i),
                             res = quote(res),
                             t_pred = quote(t_pred)) {

  dt_psu <- get_psu_stats(dt_samp = dt_samp, y = y, id_psu = id_psu, res = res);
  dt_psu <- dt_psu[dt_samp[,
                           list(t_pred = max(eval(t_pred)),
                                N_i = max(eval(N_i)),
                                n_i = max(eval(n_i))),
                           keyby = id_psu]];

  dt_psu[, t_hat := N_i/n_i*sum_y];
  dt_psu[is.nan(t_hat) == TRUE, t_hat := 0];

  dt_psu[, v_y := N_i^2*(1/n_i - 1/N_i)*var_y];
  dt_psu[is.na(v_y) == TRUE, v_y := 0];

  dt_psu[, t_res := N_i/n_i*sum_res];
  dt_psu[is.nan(t_res) == TRUE, t_res := 0];

  dt_psu[, ':='(t_hat_ma = t_pred + t_res,
                v_res = N_i^2*(1/n_i - 1/N_i)*var_res)];
  dt_psu[is.na(v_res) == TRUE, v_res := 0];

  # Equation 8.9.7  in Särndal et al. (1992)
  t_hat <- NI/nI * dt_psu[, sum(t_hat_ma)];

  # Equation 8.9.27 in Särndal et al. (1992)
  v_psu <- NI^2*(1/nI - 1/NI)*dt_psu[, var(t_hat) - 1/nI*sum(v_y)];
  v_ssu <- (NI/nI)^2*dt_psu[, sum(v_res)];
  v_ht <- v_psu + v_ssu;

  # Equation 8 in Saarela et al. (2016)
  v_pb_ma <- NI^2*(1/nI - 1/NI)*dt_psu[, var(t_hat_ma)] + NI/nI*dt_psu[, sum(v_res)];

  # Equation 4.6.2 in Särndal et al. (1992)
  v_alt <- NI^2*1/nI*dt_psu[, var(t_hat_ma)];

  return(list(t_hat = t_hat, v_ht = v_ht, v_pb_ma = v_pb_ma, v_alt = v_alt));
}

#'Horvitz-Thompson estimators for two-stage simple random sampling
#'
#'The population is separated into a set of contiguous PSUs (clusters/strips),
#'which are sampled in the first stage. SSUs are selected from the selected PSUs
#'in the second stage independently.
#'
#'@param dt_sample A data.table object containing information about the PSUs and
#'  SSUs selected in the two sampling stages.
#'@param key_var Name of the variable identifying the PSUs. Must be present in
#'  both data.tables.
#'@param y Quoted name of the variable storing the observed values of the study
#'  variable.
#'@param N_i quoted name of the variable storing the size of the PSUs
#'@param n_i quoted name of the variable storing the number of selected SSUs in
#'  each PSUs
#'@param NI Number of PSUs
#'@param nI PSU sample size
#'
#'  @details See S?rndal et al. (1992), equation 4.3.21 for the estimation of
#'  the population total and Equation 4.3.22 for variance estimation. The
#'  variance estimator may under some circumstances deliver negative estimates.
#'  Therefore, an alternative variance estimator that ignores the variation that
#'  comes from estimating PSU totals is provided as an alternative (Särndal et
#'  al., p. 154, Equation 4.6.2).
#'
#'@return A named list containing estimates for the total and its variance
#'  estimate
#'
#'@references Särndal, C. E., B. Swensson and J. Wretman (1992). Model assisted
#'  survey sampling. New York, Springer-Verlag.
#'
#'@import data.table
#'@export
est_ts_srswor <- function(dt_samp, NI, nI,
                          id_psu = 'id',
                          y = quote(y),
                          N_i = quote(N_i),
                          n_i = quote(n_i)) {

  dt_psu <- get_psu_stats(dt_samp = dt_samp, y = y, id_psu = id_psu);
  dt_psu <- dt_psu[dt_samp[,
                           list(N_i = max(N_i),
                                n_i = max(n_i)),
                           keyby = id_psu]];

  dt_psu[, t_hat := N_i/n_i * sum_y];
  dt_psu[is.nan(t_hat) == TRUE, t_hat := 0];

  dt_psu[, v_y := N_i^2*(1/n_i - 1/N_i) * var_y];
  dt_psu[is.na(v_y) == TRUE, v_y := 0];

  # Equation 4.3.21
  t_hat <- NI/nI * dt_psu[, sum(t_hat)];

  # Equation 4.3.23
  v_psu <- NI^2*(1/nI - 1/NI)*dt_psu[, var(t_hat)];
  v_ssu <- NI/nI*dt_psu[, sum(v_y)];
  v_hat <- v_psu + v_ssu;

  # Equation 4.6.2
  v_hat_alt <- NI^2*1/nI*dt_psu[, var(t_hat)];

  return(list(t_hat = t_hat, v_hat = v_hat, v_hat_alt = v_hat_alt));
}

#'Model-based inference for two-stage surveys
#'
#'The population is separated into a set of contiguous PSUs (clusters/strips),
#'which are sampled in the first stage. SSUs are selected from the selected PSUs
#'in the second stage independently.
#'
#'@param dt_sample A data.table object containing information about the PSUs and
#'  SSUs selected in the two sampling stages.
#'@param X Matrix containing the auxiliary variables (observations stored in
#'  rows, variables in columns). Can be a column vector of auxiliary averages in
#'  case linear models are used.
#'@param parm_cov Model parameter covariance matrix.
#'@param t_pred Quoted name of the variable storing predicted PSU totals
#'@param N Population size (Numbr of SSUs in the population).
#'@param NI Number of PSUs
#'@param nI PSU sample size
#'
#'@details See Ståhl et al. (2011).
#'
#'@return A named list containing estimates for the total and its variance
#'  estimate
#'
#'@references Ståhl, G., S. Holm, T. G. Gregoire, T. Gobakken, E. Næsset and R.
#'  Nelson (2011). "Model-based inference for biomass estimation in a LiDAR
#'  sample survey in Hedmark County, Norway." Canadian Journal of Forest
#'  Research 41: 96-107.
#'
#'@import data.table
#'@export
est_ts_mb <- function(dt_samp, X, parm_cov, t_pred = quote(t_pred), N, NI, nI, N_i = quote(N_i), a_psu) {
  t_hat = a_psu*dt_samp[, sum(eval(t_pred))];

  dt_samp[, ':='(w_i = eval(N_i)/sum(eval(N_i)))];
  dt_samp[, mean_pred := eval(t_pred)/eval(N_i)];

  # Equation 13 in Ene et al. (2016)
  v_pb <- dt_samp[, est_var_succ_diff(y = mean_pred, N = NI, n = nI, w = w_i)];
  v_mb <- X%*%parm_cov%*%t(X);

  v_hat <- N^2*(v_pb + v_mb);

  return(list(t_hat = t_hat, v_hat = v_hat));
}
