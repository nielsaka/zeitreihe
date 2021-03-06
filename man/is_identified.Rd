% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SVAR-3-identification.R
\name{is_identified}
\alias{is_identified}
\title{Verify whether an SVAR model is identified}
\usage{
is_identified(A = NULL, B = NULL, SIGMA_U = NULL)
}
\arguments{
\item{A}{A square numeric matrix, the coefficient of contemporaneous effects
between endogenous variables. Its dimension is (K x K).}

\item{B}{A square numeric matrix, the coefficient of contemporaneous effects
of structural shocks on endogenous variables. Its dimension is (K x
K).}

\item{SIGMA_U}{A square numeric matrix, the reduced-form residual
covariances. Its dimension is (K x K). The default setting assumes unit
variance of the structural shocks.}
}
\description{
Verify whether an SVAR model is identified
}
\note{
"The default setting assumes unit variance of the structural shocks." TRUE?
}
