###############################################################################.
#'Create data using an AR(1)
#'
#'Given some starting values, coefficients and a sequence of error vectors,
#'\code{create_arp_data} will compute a sequence of observables through a
#'univariate autoregressive model.
#'
#'@param a    A numeric vector, specifying the lag coefficients. Its length
#'  will determine the maximum lag order `p`.
#'@param y0   A numeric vector, providing pre-sample observations. Its
#'  length must be at minimum `p` (see above).
#'@param e    A numeric vector, containing the innovations that hit the system.
#'  If unspecified, it will be drawn from a standard normal distribution.
#'@param N    An integer scalar. The sample size. Defaults to `N = 100`.
#'@param intercept A logical scalar. If the first entry of `a` is the value of
#'  the intercept, set to `TRUE`. Otherwise leave at `FALSE`.
#'@return An numeric vector of length `N`.
#'@export
create_arp_data <- function(a, y0, e = rnorm(N), N = 100, intercept = FALSE) {
  A <- t(as.matrix(a))
  Y_0 <- t(as.matrix(y0))
  EE <- t(as.matrix(e))

  as.vector(create_varp_data(A, Y_0, EE))
}
###############################################################################.
#' Create data using a reduced-form VAR
#'
#' Given some starting values, coefficients and a sequence of error vectors,
#' \code{create_varp_data} will compute a sequence of observables using a simple
#' vector autoregressive process.
#'
#' @param A   A `(K x Kp)` matrix, providing the coeffcients for lag `1` to `p`
#'   with the first row containing the coefficents of the first equation.
#'   Parameter `p` is the maximum lag length and `K` the number of variables.
#' @param Y0 A `(K x p)` matrix which will be used as starting values. The
#'   first column corresponds to the very first time period.
#' @param U  A `(K x N)` matrix, providing the sequence of error vectors.
#' @return A `(K x (N + p))` matrix holding the observations. The first `p`
#'   columns will be equal to `Y0`. Column `p + 1` will be equal to `A \\\%*\\\%
#'   Y0 + U[, 1]`, where `U` contains reduced form errors. The final observation
#'   of the `K` variables will be in column `N`.
#' @family functions for creating data
#' @section Note: For a faster implementation, see
#' [this solution](http://gallery.rcpp.org/articles/simulating-vector-autoregressive-process/)
#' by Dirk Eddelbuettel.
#' @export
create_varp_data <- function(A, Y0, U) {
  # no intercept
  K <- nrow(A)
  p <- ncol(A) / K
  Tt <- ncol(U)

  YY <- matrix(, ncol = Tt + p, nrow = K)
  rownames(YY) <- paste0("y", seq_len(K))
  Y0 <- as.matrix(Y0)
  stopifnot(dim(Y0)[1] == K && dim(Y0)[2] == p)
  YY[, 1:p] <- Y0
  for (t in 1:Tt) {
    # Alternatively, just use companion form. Indexing is even simpler.
    YY[, t + p] = A %*% as.vector(YY[, (t+p-1):t]) + U[, t]
  }
  return(YY)
}
###############################################################################.
