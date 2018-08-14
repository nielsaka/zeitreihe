###############################################################################.
#'Transform a matrix of time-series by adding lags
#'
#'A matrix Y is transformed into Z notation (see Lütkepohl (2005, p. 70)). This
#'facilitates estimation of VAR(p) models.
#'
#'@param Y A (K x T) matrix, containing the data of the K variables.
#'@param p An integer, the number of lags to create. Must be greater or equal 0.
#'
#'@return A ((K * (p + 1)) x (T - p)) matrix, the first K rows will contain the
#'  original data and the remaining rows will contain its lags up to order p.
#'  The columns will contain time periods, of which there are now T - p left
#'  after setting aside pre-sample values for the lags.
#'
#'@section Note: For regression purposes, the first K rows would have to be
#'  discarded as they contain Y itself; see example.
#'
#'@examples
# TODO: example not working anymore!
#' Y <- draw_data_DH08_m1(n = 5e5)
#'
#' X <- Y2Z(Y, p = 4)
#' Y <- X[seq_len(4), ] # has correct sample size
#' Z <- X[-seq_len(4), ]
#'
#' B <- Y %*% t(Z) %*% solve(Z %*% t(Z))
#'
#'@section ToDo:
#' add sanity checks... data.frame? matrix? dimension? max, min p?
#' no intercept (yet)!
#' rownames?
#'
Y2Z <- function(Y, p) {
  # TODO: add sanity checks... data.frame? matrix? dimension?
  # no intercept (yet)!
  Y <- t(as.matrix(Y))
  t(embed(Y, p + 1))
}
###############################################################################.
#' Title
#'
#' @param A A `(K x Kp)` matrix, ...
#'
#' @return
#' @export
#'
#' @examples
#'
#' A <- matrix(c(.5, .4, .1, .5, 0, .25, 0, 0), nrow = 2)
#' check_stability(A)
#'
#' A <- matrix(c(.5, .4, .1, .5, 0, .25, 1, 1), nrow = 2)
#' check_stability(A)
#'
check_stability <- function(A) {
  K <- var_length(A)
  p <- lag_length(A)

  AA <- companion_format(A)

  # some severe numeric mistakes ... floating point arithmetic?
  all(abs(eigen(AA)$values) < 1)

  # could also check solutions to characteristic polynomial
  # efficiency?
}
###############################################################################.
#' Title
#'
#' @param A
#' @param v
#'
#' @return
#' @export
#'
#' @examples
mean_process <- function(A, v) {
  K <- var_length(A)
  p <- lag_length(A)

  AA <- companion_format(A)

  solve(diag(K * p) - AA) %*% v
}
###############################################################################.
#' Title
#'
#' @param A
#' @param SIGMA
#'
#' @return
#' @export
#'
#' @examples
covariance <- function(A, SIGMA) {
# TODO
# infinite sum... convergent ... tolerance?

}
###############################################################################.
#' Convert VAR(1) companion format to VAR(p)
#'
#' @param K
#' @param p
#'
#' @return
#' @export
#'
#' @examples
#'
#' J <- selector(4, 2)
#'
selector <- function(K, p){
  out <- matrix(0, K, K * p)
  diag(out) <- 1
  out
}
###############################################################################.
#' Convert VAR(p) into VAR(1) companion format
#'
#' @param A
#'
#' @return
#' @export
#'
#' @examples
companion_format <- function(A) {
  K <- var_length(A)
  p <- lag_length(A)

  XX <- matrix(0, K * p - K, K * p)
  diag(XX) <- 1
  rbind(A, XX)
}
###############################################################################.
#' Extract number of variables
#'
#' @param A
#'
#' @return
#'
#' @examples
var_length <- function(A) {
  nrow(A)
}
###############################################################################.
#' Extract number of lags
#'
#' @param A
#'
#' @return
#'
#' @examples
#'
#' TODO Do not run:
#' A <- matrix(0, 4, 7)
#' lag_length(A)
#'
lag_length <- function(A) {
  p <- ncol(A) / nrow(A)
  stopifnot(p %% 1 == 0)
  p
}
###############################################################################.
