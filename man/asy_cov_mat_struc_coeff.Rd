% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SVAR-2-estimation.R
\name{asy_cov_mat_struc_coeff}
\alias{asy_cov_mat_struc_coeff}
\title{Title}
\usage{
asy_cov_mat_struc_coeff(A, B, RA, RB)
}
\arguments{
\item{RB}{}
}
\description{
Title
}
\examples{

K <- 2

B <- matrix(c(1, 0, 0, 1), K, K)
RB <-matrix(c(1, rep(0, 6), 1), K^2, K)

asy_cov_mat_struc_coeff(B = B, RB = RB)

}
