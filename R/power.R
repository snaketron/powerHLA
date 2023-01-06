
#' get_power
#'
#' @param y, numeric vector, Multinomially distributed sample (background)
#' @param y_rng, numeric matrix, Multinomially distributed samples (rows)
#' @param a, numeric vector, Dirichlet distribution parameters (alpha)
#' @param rng_draws, integer, number of posterior draws (default 10,000)
#'
#' @return data.frame of posterior summaries
#' @export
#'
#' @examples
get_power <- function(y, y_rng, a, rng_draws = 10^4) {
  s <- rstan::sampling(
    object = stanmodels$dm,
    data = list(K=ncol(y_rng),
                B=nrow(y_rng),
                y_rng=y_rng,
                y = y,
                a = a),
    chains = 1,
    cores = 1,
    iter = rng_draws,
    warmup = 0,
    algorithm="Fixed_param",
    refresh = -1)

  # keep summary intervals
  probs <- c(0.005, 0.025, 0.05, 0.25, 0.50, 0.75, 0.95, 0.975, 0.995)

  # p_rng
  p_rng <- data.frame(rstan::summary(s, par = "p_rng", probs = probs)$summary)
  colnames(p_rng) <- paste0("p_rng_", colnames(p_rng))

  # delta_diff
  d_diff <- data.frame(rstan::summary(s, par = "d_diff", probs = probs)$summary)
  colnames(d_diff) <- paste0("diff_", colnames(d_diff))

  # delta_or
  d_or <- data.frame(rstan::summary(s, par = "d_or", probs = probs)$summary)
  colnames(d_or) <- paste0("or_", colnames(d_or))

  # delta_lor
  d_lor <- data.frame(rstan::summary(s, par = "d_lor", probs = probs)$summary)
  colnames(d_lor) <- paste0("lor_", colnames(d_lor))

  if(is.null(colnames(y_rng))) {
    colnames(y_rng) <- paste0("K", 1:ncol(y_rng))
  }
  d_lor$allele <- rep(x = colnames(y_rng), times = nrow(y_rng))
  d_lor$B <- rep(x = 1:nrow(y_rng), each = ncol(y_rng))
  d_lor$N <- sum(y_rng[1,])

  # resulting data.frame
  d <- cbind(p_rng, d_diff, d_or, d_lor)

  # return
  return(d)
}





#' get_power_analysis
#'
#' @param N, integer, number of alleles in study
#' @param K, integer, number of categories
#' @param af, vector of reals, allele frequencies of the simulated samples
#' @param B, integer, number of simulations to use
#' @param y, vector of integers, allele counts of the background sample
#' @param a, vector of reals, Dirichlet distribution parameters (alpha)
#' @param rng_draws, integer, number of posterior draws (default 10,000)
#'
#' @return
#' @export
#'
#' @examples
get_power_analysis <- function(N, K, af, B, y, a,
                               rng_draws = 10^4,
                               verbose = T) {

  if(verbose) {
    base::message("1) rng")
  }

  # rng
  y_rng <- powerHLA::get_multinomial_rng(
    K = K,
    af = af,
    B = B,
    N = N)
  colnames(y_rng) <- paste0("K", 1:ncol(y_rng))

  if(verbose) {
    base::message("1) power comp.")
  }

  # power
  p <- powerHLA::get_power(
    y = y,
    y_rng = y_rng,
    a = a,
    rng_draws = rng_draws)

  # return
  return(p)
}

