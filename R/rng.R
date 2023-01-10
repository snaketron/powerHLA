
#' get_multinomial_rng
#'
#' This function simulates B samples (rows of counts of matrix) from
#' a multivariate distribution. Each count belongs to one of K categories
#' (K columns of matrix). The data is simulated based on the given
#' simplex theta and N number of total alleles.
#'
#' @param theta vector, allele frequencies, sum(theta)=1
#' @param B integer, number of simulated samples
#' @param N integer, number of alleles in simulated sample
#' @param rng_seed integer, seed to be used by rstan for rng
#'
#' @return y_rng numeric matrix, number of alleles in each simulated sample
#' @export
#'
#' @examples
#' get_multinomial_rng(theta=rep(0.1, times = 10), B=100, N=1000)
get_multinomial_rng <- function(theta, B, N, rng_seed = NULL) {

  K <- length(theta)

  if(missing(rng_seed)|is.null(rng_seed)) {
    rng_seed = sample.int(.Machine$integer.max, 1)
  }

  s <- rstan::sampling(
    object = stanmodels$rng_mult,
    data = list(K=K, N=N, theta=theta),
    seed = rng_seed,
    chains = 1,
    cores = 1,
    iter = B,
    warmup = 0,
    algorithm="Fixed_param",
    refresh = -1)

  # extract rng
  rng_gamma <- rstan::extract(
    object = s,
    permuted = FALSE,
    par = "gamma")[,1,]

  # set dimension names
  dn <- vector(mode = "list", length = 2)
  names(dn) <- c("B", "K")
  dimnames(rng_gamma) <- dn

  return(rng_gamma)
}




#' get_dirichlet_rng
#'
#' This function simulates B samples (rows of matrix) from a dirichlet
#' distribution. Each row entry (probability) belongs to one of K categories
#' (K columns of matrix). The data is simulated based on a given vector of
#' reals larger than 0.
#'
#' @param alpha vector, pseudo allele counts
#' @param B integer, number of simulated samples
#' @param rng_seed integer, seed to be used by rstan for rng
#'
#' @return y_rng numeric matrix, frequency of alleles in each simulated sample
#' @export
#'
#' @examples
#' get_dirichlet_rng(alpha=rep(0.1, times = 10), B=100)
get_dirichlet_rng <- function(alpha, B, rng_seed = NULL) {

  K <- length(alpha)

  if(missing(rng_seed)) {
    rng_seed = sample.int(.Machine$integer.max, 1)
  }

  s <- rstan::sampling(
    object = stanmodels$rng_dir,
    data = list(K=K, alpha=alpha),
    seed = rng_seed,
    chains = 1,
    cores = 1,
    iter = B,
    warmup = 0,
    algorithm="Fixed_param",
    refresh = -1)

  # extract rng
  rng_theta <- rstan::extract(
    object = s,
    permuted = FALSE,
    par = "theta")[,1,]

  # set dimension names
  dn <- vector(mode = "list", length = 2)
  names(dn) <- c("B", "K")
  dimnames(rng_theta) <- dn

  return(rng_theta)

}

