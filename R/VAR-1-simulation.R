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
create_arp_data <- function(a, y0, e = rnorm(N), N = 100, intercept = FALSE) {
  B <- t(as.matrix(a))
  Y_0 <- t(as.matrix(y0))
  EE <- t(as.matrix(e))

  as.vector(create_varp_data(B, Y_0, EE))
}
###############################################################################.
#' Create data using a VAR(1)
#'
#' Given some starting values, coefficients and a sequence of error vectors,
#' \code{create_var1_data} will compute a sequence of observables according to a
#' vector autoregressive model of order one.
#'
#' @param A   A (K x K) matrix, providing the coeffcients for lag one.
#' @param Y_0 A (K x 1) vector or column matrix, will be used as starting
#'   values.
#' @param EE  A (K x T) matrix, providing the sequence of error vectors.
#' @return A (K x T) matrix of observables. The first column will equal \code{A
#'   \%*\% Y_0 + EE[, 1]}. The final observation of the K variables will be in
#'   column T.
#' @section Note: For a faster implementation, see
#'   \href{http://gallery.rcpp.org/articles/simulating-vector-autoregressive-process/}{this
#'    solution} by Dirk Eddelbuettel.
#' @export
create_var1_data <- function(A, Y_0, EE) {
  YY <- matrix(, ncol = ncol(EE) + 1, nrow = nrow(EE))
  rownames(YY) <- paste0("y", seq_len(nrow(EE)))
  YY[, 1] <- Y_0
  for (t in seq_len(ncol(EE))) {
      YY[, t + 1] = A %*% YY[, t] + EE[, t]
    }
  return(YY)
}
###############################################################################.
#' Create data using a VAR(p)
#'
#' Given some starting values, coefficients and a sequence of error vectors,
#' \code{create_varp_data} will compute a sequence of observables using a simple
#' vector autoregressive process.
#'
#' @param B   A `(K x Kp)` matrix, providing the coeffcients for lag `1` to
#'   `p` with the first row containing the coefficents of the first equation.
#'   Parameter `p` is the maximum lag length and `K` the number of variables.
#' @param Z_0 A `(K x p)` matrix which will be used as starting values. The first
#'   column corresponds to the very first time period.
#' @param UU  A `(K x N)` matrix, providing the sequence of error vectors.
#' @return A `(K x (N + p))` matrix holding the observations. The first `p` columns will
#'   be equal to `Z_0`. The final observation of the `K` variables will be in
#'   column `N`.
#' @section Note: For a faster implementation, see
#' [this solution](http://gallery.rcpp.org/articles/simulating-vector-autoregressive-process/)
#' by Dirk Eddelbuettel.
#' @export
create_varp_data <- function(B, Z_0, UU) {
  # no intercept
  K <- nrow(B)
  p <- ncol(B) / K
  Tt <- ncol(UU)

  YY <- matrix(, ncol = Tt + p, nrow = K)
  rownames(YY) <- paste0("y", seq_len(K))
  Z_0 <- as.matrix(Z_0)
  stopifnot(dim(Z_0)[1] == K && dim(Z_0)[2] == p)
  YY[, 1:p] <- Z_0
  for (t in 1:Tt) {
    # Alternatively, just use companion form. Indexing is even simpler.
    YY[, t + p] = B %*% as.vector(YY[, (t+p-1):t]) + UU[, t]
  }
  return(YY)
}
###############################################################################.
