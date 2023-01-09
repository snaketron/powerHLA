
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
get_power <- function(gamma, rng_gamma, alpha, rng_draws) {
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

  # delta_or
  d_or <- data.frame(rstan::summary(s, par = "d_or", probs = probs)$summary)
  d_or$Rhat <- NULL
  d_or$n_eff <- NULL
  colnames(d_or) <- paste0("or_", colnames(d_or))

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
  d <- cbind(p_rng, d_diff, d_or, d_lor)
  d$obs_theta <- rep(x = gamma/sum(gamma), times = nrow(rng_gamma))

  # return
  return(d)
}



#' get_power_run
#'
#' @param N, integer, total number of alleles in simulated samples
#' @param theta, vector of reals, allele frequencies of the simulated samples
#' @param B, integer, number of simulations to use
#' @param y, vector of integers, allele counts of the background sample
#' @param a, vector of reals, Dirichlet distribution parameters (alpha)
#' @param rng_draws, integer, number of posterior draws (default 10,000)
#' @param cores, integer, number of cores for multicore execution
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
                          rng_draws = 10^4) {

  K <- length(theta)

  # rng
  rng_gamma <- powerHLA::get_multinomial_rng(
    theta = theta,
    B = B,
    N = N)
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
#' @param B, integer, number of simulations to use
#' @param y, vector of integers, allele counts of the background sample
#' @param a, vector of reals, Dirichlet distribution parameters (alpha)
#' @param rng_draws, integer, number of posterior draws (default 10,000)
#' @param cores, integer, number of cores for multicore execution
#'
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
                               cores = 1,
                               verbose = T) {

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
    future.seed = TRUE)

  # list names -> ns entries
  names(o) <- ns

  # return
  return(o)
}




#' get_power_summary
#'
#' @param p, object returned by function get_power_analysis
#' @param hdi_level, numeric between 0 and 1 (default=0.95), that defines the
#' highest density interval to consider when summarizing parameter posteriors
#'
#' @return data.frame
#' @export
#'
#' @examples
#' get_power_summary(p)
get_power_summary <- function(p, hdi_level = 0.95) {

  # min(abs(x))
  get_min_abs <- function(x) {
    return(min(abs(x)))
  }

  # list to data.frame
  w <- do.call(rbind, p)

  # check if effect is detected
  w$diff_effect_HDI95 <- ifelse(test = w$diff_X2.5.<=0&
                                  w$diff_X97.5.>=0,
                                yes = 0, no = 1)

  # compute minimum effect
  w$diff_min_effect_95 <- apply(
    X = w[,c("diff_X2.5.", "diff_X97.5.")], MARGIN = 1,
    FUN = get_min_abs)*w$diff_effect_HDI95

  # compute power
  s <- aggregate(diff_effect_HDI95~N+allele, data = w, FUN = sum)
  s$diff_effect_HDI95_pct <- s$diff_effect_HDI95/max(w$B)*100

  return (s)
}

