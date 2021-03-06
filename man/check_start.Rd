% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aux_funs.R
\name{check_start}
\alias{check_start}
\alias{check_start_all}
\title{Do all names start with "it"?}
\usage{
check_start(x, start)

check_start_all(x, start)
}
\arguments{
\item{x}{A vector. The names of the vector will be checked.}

\item{start}{A character string.}
}
\value{
\itemize{
\item check_start \cr A boolean vector of the same length as \code{x} if \code{x}
is named. Each element is either TRUE or FALSE. They are TRUE if the name
of the element starts with \code{start}. If \code{x} did not carry names, a vector of
length zero is returned.
}

\itemize{
\item \code{check_start_all} \cr A boolean scalar. TRUE if the names of \strong{all}
elements of \code{x} start with \code{start}. Otherwise FALSE.
}
}
\description{
Check whether the names of a vector all start with a certain string.
}
\examples{
x <- c(tisch = 1, teller = 2, tafel = 3)
name_start(x, "t")

x <- c(x, 4)
name_start(x, "t")

x <- 1:4
name_start(x, "t")
x <- c(tisch = 1, teller = 2, tafel = 3)
check_names(x, "t")

x <- c(x, 4)
check_names(x, "t")

x <- 1:4
check_names(x, "t")
}
