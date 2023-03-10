% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rng.R
\name{get_dirichlet_rng}
\alias{get_dirichlet_rng}
\title{get_dirichlet_rng}
\usage{
get_dirichlet_rng(alpha, B, N, rng_seed = NULL)
}
\arguments{
\item{alpha}{vector, pseudo allele counts}

\item{B}{integer, number of simulated samples}

\item{N}{integer, number of alleles in simulated sample}

\item{rng_seed}{integer, seed to be used by rstan for rng}
}
\value{
y_rng numeric matrix, frequency of alleles in each simulated sample
}
\description{
This function simulates B samples (rows of matrix) from a dirichlet
distribution. Each row entry (probability) belongs to one of K categories
(K columns of matrix). The data is simulated based on a given vector of
reals larger than 0.
}
\examples{
get_dirichlet_rng(alpha=rep(0.1, times = 10), B=100, N=1000)
}
