
#' get_power
#'
#' @param gamma, numeric vector, Multinomially distributed sample (background)
#' @param rng_gamma, numeric matrix, Multinomially distributed samples (rows)
#' @param alpha, numeric vector, Dirichlet distribution parameters (alpha)
#' @param rng_draws, integer, number of posterior draws (default 10,000)
#'
#' @return data.frame of posterior summaries
#' @export
#'
#' @examples
get_power <- function(gamma,
                      rng_gamma,
                      alpha,
                      rng_draws) {

  s <- rstan::sampling(
    object = stanmodels$dm,
    data = list(K=ncol(rng_gamma),
                B=nrow(rng_gamma),
                rng_gamma=rng_gamma,
                gamma = gamma,
                alpha = alpha),
    chains = 1,
    cores = 1,
    iter = rng_draws,
    warmup = 0,
    algorithm="Fixed_param",
    refresh = -1)

  # keep summary intervals
  probs <- c(0.005, 0.025, 0.05, 0.50, 0.95, 0.975, 0.995)

  # p_rng
  p_rng <- data.frame(rstan::summary(s, par = "p_rng", probs = probs)$summary)
  p_rng$Rhat <- NULL
  p_rng$n_eff <- NULL
  colnames(p_rng) <- paste0("p_rng_", colnames(p_rng))

  # delta_diff
  d_diff <- data.frame(rstan::summary(s, par = "d_diff", probs = probs)$summary)
  d_diff$Rhat <- NULL
  d_diff$n_eff <- NULL
  colnames(d_diff) <- paste0("diff_", colnames(d_diff))

  # delta_lor
  d_lor <- data.frame(rstan::summary(s, par = "d_lor", probs = probs)$summary)
  d_lor$Rhat <- NULL
  d_lor$n_eff <- NULL
  colnames(d_lor) <- paste0("lor_", colnames(d_lor))

  if(is.null(colnames(rng_gamma))) {
    colnames(rng_gamma) <- paste0("K", 1:ncol(rng_gamma))
  }
  d_lor$allele <- rep(x = colnames(rng_gamma), times = nrow(rng_gamma))
  d_lor$B <- rep(x = 1:nrow(rng_gamma), each = ncol(rng_gamma))
  d_lor$N <- sum(rng_gamma[1,])

  # resulting data.frame
  d <- cbind(p_rng, d_diff, d_lor)
  d$obs_theta <- rep(x = gamma/sum(gamma), times = nrow(rng_gamma))

  # return
  return(d)
}



#' get_power_run
#'
#' @param N, integer, total number of alleles in simulated samples
#' @param theta, vector of reals, allele frequencies of the simulated samples
#' @param gamma, vector of integers, allele counts of the background sample
#' @param alpha, vector of reals, Dirichlet distribution parameters (alpha)
#' @param B, integer, number of simulations to use
#' @param rng_draws, integer, number of posterior draws (default 10,000)
#' @param rng_seed integer, seed to be used by rstan for rng
#'
#' @return
#' @export
#'
#' @examples
get_power_run <- function(N,
                          theta,
                          gamma,
                          alpha,
                          B,
                          rng_draws,
                          rng_seed) {

  K <- length(theta)

  # rng
  rng_gamma <- powerHLA::get_multinomial_rng(
    theta = theta,
    B = B,
    N = N,
    rng_seed = rng_seed)
  colnames(rng_gamma) <- paste0("K", 1:ncol(rng_gamma))


  # power
  p <- powerHLA::get_power(
    gamma = gamma,
    rng_gamma = rng_gamma,
    alpha = alpha,
    rng_draws = rng_draws)

  # return
  return(p)
}



#' get_power_analysis
#'
#' @param ns, integer vector, total number of alleles in simulated samples
#' @param theta, vector of reals, allele frequencies of the simulated samples
#' @param gamma, vector of integers, allele counts of the background sample
#' @param alpha, vector of reals, Dirichlet distribution parameters (alpha)
#' @param B, integer, number of simulations to use
#' @param rng_draws, integer, number of posterior draws (default 10,000)
#' @param rng_seed integer, seed to be used by rstan for rng
#' @param cores, integer, number of cores for multicore execution
#'l
#' @return
#' @export
#'
#' @examples
get_power_analysis <- function(ns,
                               theta,
                               gamma,
                               alpha,
                               B,
                               rng_draws = 10^4,
                               rng_seed = NULL,
                               cores = 1) {

  # schedule cores
  future::plan(future::multisession, workers = cores)

  # loop over ns and execute get_power_run
  o <- future.apply::future_lapply(
    X = ns,
    FUN = get_power_run,
    theta = theta,
    gamma = gamma,
    alpha = alpha,
    B = B,
    rng_draws = rng_draws,
    rng_seed = rng_seed,
    future.seed = TRUE)

  # list names -> ns entries
  names(o) <- ns

  # return
  return(o)
}




#' get_power_summary
#'
#' @param p, object returned by function get_power_analysis
#'
#' @return data.frame
#' @export
#'
#' @examples
#' get_power_summary(p)
get_power_summary <- function(p) {

  # min(abs(x))
  get_min_abs <- function(x) {
    return(min(abs(x)))
  }

  # list to data.frame
  w <- do.call(rbind, p)

  # check if effect is detected
  w$zero_in_HDI95 <- ifelse(test = w$diff_X2.5.<=0&
                              w$diff_X97.5.>=0,
                            yes = 1, no = 0)
  w$zero_in_HDI99 <- ifelse(test = w$diff_X0.5.<=0&
                              w$diff_X99.5.>=0,
                            yes = 1, no = 0)

  # compute minimum effect
  # w$diff_min_effect_95 <- apply(
  #   X = w[,c("diff_X2.5.", "diff_X97.5.")], MARGIN = 1,
  #   FUN = get_min_abs)*w$diff_effect_HDI95

  # compute power
  s_95 <- aggregate(zero_in_HDI95~N+allele, data = w, FUN = sum)
  s_95$zero_pct_HDI95 <- s_95$zero_in_HDI95/max(w$B)*100
  # s_95$fp_HDI95 <- 100-s_95$tp_HDI95

  s_99 <- aggregate(zero_in_HDI99~N+allele, data = w, FUN = sum)
  s_99$zero_pct_HDI99 <- s_99$zero_in_HDI99/max(w$B)*100
  # s_99$fp_HDI99 <- 100-s_99$tp_HDI99


  # compute distribution of deltas
  d_95 <- aggregate(diff_mean~N+allele, data = w,
                    FUN = get_hdi, hdi_level = 0.95)
  d_95$mean_delta_L95 <- d_95$diff_mean[, 1]
  d_95$mean_delta_H95 <- d_95$diff_mean[, 2]
  d_95$diff_mean <- NULL

  s <- merge(x = s_95, y = s_99, by = c("N", "allele"))
  rm(s_95, s_99)

  s <- merge(x = s, y = d_95, by = c("N", "allele"))

  return (s)
}



#' get_hdi
#' Computes HDI for vector vec and hdi_level (e.g. 0.95). Taken (and renamed)
#' from "Doing Bayesian Analysis", section 25.2.3 R code for computing HDI of
#' a MCMC sample
#'
#' @param vec numeric vector
#' @param hdi_level number, highest density interval (HDI) (default = 0.95)
#'
#' @return vector with two elements (low, high interval of HDI)
#'
#' @examples
#'
get_hdi <- function(vec, hdi_level) {
  # Computes highest density interval from a sample of representative values,
  # estimated as shortest credible interval.
  # Arguments:
  # sampleVec
  # is a vector of representative values from a probability distribution.
  # credMass
  # is a scalar between 0 and 1, indicating the mass within the credible
  # interval that is to be estimated.
  # Value:
  # HDIlim is a vector containing the limits of the HDI
  sortedPts <- sort(vec)
  ciIdxInc <- floor(hdi_level * length(sortedPts))
  nCIs = length(sortedPts) - ciIdxInc
  ciWidth = rep(0 , nCIs)
  for (i in 1:nCIs) {
    ciWidth[i] = sortedPts[i + ciIdxInc] - sortedPts[i]
  }
  HDImin = sortedPts[which.min(ciWidth)]
  HDImax = sortedPts[which.min(ciWidth) + ciIdxInc]
  HDIlim = c(HDImin, HDImax)
  return(HDIlim)
}
