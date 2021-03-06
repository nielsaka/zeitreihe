% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VAR-2-estimation.R
\name{log_lik_init}
\alias{log_lik_init}
\title{Create a function to compute the log-likelihood}
\usage{
log_lik_init(Y, p)
}
\arguments{
\item{Y}{A \code{(K x N+p)} matrix carrying the data for estimation. There are
\code{N} observations for each of the \code{K} variables with \code{p} pre-sample values.}

\item{p}{An integer scalar. The lag length of the VAR(p) system.}
}
\value{
A function is returned. It takes as input a named vector \code{args}. This
vector consists of the parameters of the VAR(p) model in vectorised form.
Note that names \code{mu, a, s} are a requirement. It will return the value of
the log-likelihood at the specified parameter vector. See the example for
details.
}
\description{
Create the log-likelihood function of a VAR(p) for a particular data set.
}
\examples{
K <- 3
N <- 5E2
p <- 2

set.seed(8191)

A <- matrix(0.1, K, K * p)
Y0 <- matrix(0, K, p)
U <- matrix(rnorm(K * N), K, N)

Y <- create_varp_data(A, Y0, U)

log_lik <- log_lik_init(Y, p)

mu <- rep(0, K)
SIGMA <- diag(K)
args = c(mu = mu, a = vec(A), s = vech(SIGMA))

log_lik(args)
}
