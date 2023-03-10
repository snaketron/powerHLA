% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/power.R
\name{get_power}
\alias{get_power}
\title{get_power}
\usage{
get_power(gamma, rng_gamma, alpha, rng_draws)
}
\arguments{
\item{gamma, }{numeric vector, Multinomially distributed sample (background)}

\item{rng_gamma, }{numeric matrix, Multinomially distributed samples (rows)}

\item{alpha, }{numeric vector, Dirichlet distribution parameters (alpha)}

\item{rng_draws, }{integer, number of posterior draws (default 10,000)}
}
\value{
data.frame of posterior summaries
}
\description{
get_power
}
