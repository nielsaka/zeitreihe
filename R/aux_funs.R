###############################################################################.
#'Transform a matrix of time series observations by appending lags
#'
#'A matrix `Y` is transformed into `Z` notation (see Lütkepohl (2005, p. 70)).
#'This facilitates estimation of VAR(p) models.
#'
#'@inheritParams ols_mv
#'
#'@return A `([K * (p + 1)] x [N - p])` matrix, the first K rows will contain
#'  the original data and the remaining rows will contain its lags up to order
#'  p. The columns will contain time periods, of which there are now `N - p`
#'  left after setting aside pre-sample values for the lags.
#'
#'@section Note: For regression purposes, the first K rows would have to be
#'  discarded as they contain Y itself; see example.
#'
#'@examples
#' K <- 3
#' p <- 2
#' N <- 2E4
#'
#' A <- matrix(0.1, K, K * p)
#' Y0 <- matrix(0, K, p)
#' U <- matrix(rnorm(K * N), K, N)
#' Y <- create_varp_data(A, Y0, U)
#'
#' X <- Y2Z(Y, p)
#' Y <- X[seq_len(K), ] # has correct sample size
#' Z <- X[-seq_len(K), ]
#'
#' B <- Y %*% t(Z) %*% solve(Z %*% t(Z))
#'
#'@section TODO:
#' add sanity checks... data.frame? matrix? dimension? max, min p?
#' rownames?
#'
Y2Z <- function(Y, p, const = FALSE) {
  # TODO: add sanity checks... data.frame? matrix? dimension?
  Y <- t(as.matrix(Y))
  nu <- if(const) 1 else numeric(0)
  rbind(nu, t(embed(Y, p + 1)))
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
#' A <- matrix(c(1, 0, 0, 1), nrow = 2)
#' check_stability(A)
#'
#' A <- matrix(c(1, 0, 0, 1, -0.1, 0, 0, -0.1), nrow = 2)
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
#'
#' @examples
mean_var_process <- function(A, v) {
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
#' @param h
#' @param tol
#'
#' @return
#'
#' @section TODO: TEST!!!
#'
#' @section Note: implementation is plain brute force with no regard for
#'   efficiency
#'
#' @examples
#' A <- matrix(0.0, 4, 8)
#' SIGMA <- matrix(0.5, 4, 4)
#' cov_var_process(A, SIGMA)
cov_var_process <- function(A, SIGMA, h = 0, tol = 1E-4) {
  p <- lag_length(A)

  cf_SIGMA <- big_SIGMA(SIGMA, p)
  cf_GAMMA <- cf_SIGMA

  # TODO cf_A <- companion_format(A) # TODO fix arguments
  cf_A <- big_A(A)

  cf_A_hi <- cf_A
  for (i in seq_len(h)) {
    cf_A_hi <- cf_A_hi %*% cf_A
  }

  # infinite sum... convergent ... tolerance
  # TODO better formula?

  crit <- 1
  cf_A_i <- cf_A
  repeat{
    cf_GAMMA_cand <- cf_GAMMA + cf_A_hi %*% cf_SIGMA %*% t(cf_A_i)

    crit <- sum(abs(cf_GAMMA_cand - cf_GAMMA))
    if (crit > tol) cf_GAMMA <- cf_GAMMA_cand else break

    cf_A_hi <- cf_A_hi %*% cf_A
    cf_A_i  <- cf_A_i  %*% cf_A
  }
  # TODO refactor to function?
  K <- var_length(A)
  p <- lag_length(A)
  selector(K, p)  %*% cf_GAMMA %*% t(selector(K, p))
}
###############################################################################.
#' Convert VAR(1) companion format to VAR(p)
#'
#' @param K
#' @param p
#'
#' @return
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
#'
#' @examples
companion_format <- function(Y, nu, A, U) {
  K <- var_length(A)
  p <- lag_length(A)

  list(
    YY   = big_Y(Y, p),
    AA   = big_A(A),
    nunu = big_nu(nu, p),
    UU   = big_U(U, p)
  )
}
###############################################################################.
# TODO: write code, document, test
big_A <- function(A) {
  K <- var_length(A)
  p <- lag_length(A)

  XX <- matrix(0, K * p - K, K * p)
  diag(XX) <- 1
  rbind(A, XX)
}
###############################################################################.
#' Title
#'
#' @param Y A `(K x N)` matrix, ... TODO: inherit?!
#' @param p An integer scalar, the lag length TODO: inherit?!
#'
#' @return A `(Kp x N-p+1)` matrix. The `p - 1` lags of `Y` are pasted to rows
#'  below `Y`. This leads to the loss of `p - 1` sample observations.
#'
#'@section Note: difference between big_Y and Y2Z? matrix `Z` has `K*p + 1` rows
#' whereas the output of `big_Y` has `K*p` rows. The very last lag is missing.
#'
#' @examples
#' Y <- matrix(1:15, 3, 5)
#' big_Y(Y, p = 3)
big_Y <- function(Y, p) {
  Y2Z(Y, p - 1)
}
###############################################################################.
#' Title
#'
#' @param nu
#' @param p inherit?!
#'
#' @return
#'
#' @examples
#' nu <- as.matrix(1:3)
#' big_nu(nu, p = 2)
#'
big_nu <- function(nu, p) {
  # K <- var_length(nu)
  # rbind(nu, matrix(0, K * (p - 1), 1))

  one_zeros(p) %x% nu

}
###############################################################################.
#' Convert matrix of residuals into VAR(1) companion format
#'
#' @param U A `(K x N)` matrix of residuals.
#' @param p inherit ..
#'
#' @return TODO number of columns of matrix?
#' @export
#'
#' @examples
#'
#' U <- matrix(1:21, 3, 7)
#' big_U(U, p = 3)
#'
big_U <- function(U, p) {
  # K <- var_length(U)
  # N <- obs_length(U)

  one_zeros(p) %x% U

  # more efficient (?) but less readable
  # rbind(U, matrix(0, K * (p - 1), N))[, -(seq_len(p -1))]
}
###############################################################################.
#' Convert a residual covariance matrix into VAR(1) companion form
#'
#' @param SIGMA
#' @param p
#'
#' @return
#'
#' @examples
#' SIGMA <- matrix(0.5, 4, 4)
#' big_SIGMA(SIGMA, p = 3)
#'
big_SIGMA <- function(SIGMA, p) {
  one_zeros(p, p) %x% SIGMA
}

###############################################################################.
#' Extract number of variables
#'
#' @param mat
#'
#' @return
#'
#' @examples
var_length <- function(mat) {
  # TODO: warn if not a matrix
  nrow(mat)
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
#' Title
#'
#' @param mat
#'
#' @return
#' @export
#'
#' @examples
obs_length <- function(mat) {
  ncol(mat)
}
###############################################################################.
#' Title
#'
#' @param mat
#'
#' @return
#' @export
#'
#' @examples
vech <- function(mat) {
  stopifnot(all(mat == t(mat)))
  mat[lower.tri(mat, diag = TRUE)]
}
###############################################################################.
#' Title
#'
#' @param mat
#'
#' @return
#' @export
#'
#' @examples
vec <- function(mat) {
  c(mat)
}
###############################################################################.
one_zeros <- function(r, c = 1) {
  stopifnot(r == c || c == 1)
  mat <- matrix(0, r, c)
  mat[1] <- 1
  mat
}
###############################################################################.
