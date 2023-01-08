
#' get_multinomial_rng
#'
#' This function simulates B samples (rows of counts of matrix y) from
#' a multivariate distribution. Each count belongs to one of K categories
#' (K columns of matrix y). The data is simulated based on the given
#' simplex af and N number of total alleles.
#'
#' @param af vector, allele frequencies, sum(af)=1
#' @param B integer, number of simulated samples
#' @param N integer, number of alleles in simulated sample
#'
#' @return y_rng numeric matrix, number of alleles in each simulated sample
#' @export
#'
#' @examples
#' get_multinomial_rng(K=10, af=rep(0.1, times = 10), B=100, N=1000)
get_multinomial_rng <- function(af, B, N) {

  K <- length(af)

  s <- rstan::sampling(
    object = stanmodels$rng,
    data = list(K=K, N=N, af=af),
    chains = 1,
    cores = 1,
    iter = B,
    warmup = 0,
    algorithm="Fixed_param",
    refresh = -1)

  # extract rng
  y_rng <- rstan::extract(
    object = s,
    par = "y")$y

  # set dimension names
  dn <- vector(mode = "list", length = 2)
  names(dn) <- c("B", "K")
  dimnames(y_rng) <- dn

  return(y_rng)
}
