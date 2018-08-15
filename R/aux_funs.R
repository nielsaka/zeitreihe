###############################################################################.
#'Create a regressor matrix for a VAR(p) model
#'
#'Matrix `Y` is transformed into `Z` notation; see Lütkepohl (2005, p. 70).
#'This transformation facilitates estimation of VAR(p) models.
#'
#'@inheritParams ols_mv
#'
#'@return A `([K * p + 1] x [N - p])` or `([K * p] x [N - p])`matrix. The
#'  dimension depends on whether a constant is included. If it is included, the
#'  first row consists of ones. The next K rows will contain the data lagged by
#'  one time period. The remaining rows will contain further lags up to order p.
#'  The columns will contain time periods of which there are now `N - p` left
#'  after setting aside pre-sample values for the lags.
#'
#'@examples
#' K <- 3
#' N <- 6
#' p <- 2
#'
#' Y <- matrix(seq_len(K*N), nrow = K)
#' Z <- Y2Z(Y, p)
Y2Z <- function(Y, p, const = TRUE) {
  K <- var_length(Y)
  nu <- if(const) 1 else numeric(0)
  Y <- t(as.matrix(Y))

  rbind(nu, t(embed(Y, p + 1))[-seq_len(K), ])
}
###############################################################################.
#' Check the stability criterion of a VAR(p) model
#'
#' This function will check whether the process implied by the coefficient
#' matrix `A` is stable. That is, whether its eigenvalues are all smaller than
#' one in absolute terms.
#'
#' @inheritParams create_varp_data
#'
#' @return Boolean scalar. `TRUE` if the criterion is satisfied, otherwise
#'   `FALSE`.
#' @export
#'
#' @examples
#' A <- matrix(c(.5, .4, .1, .5, 0, .25, 0, 0), nrow = 2)
#' check_stability(A)
#'
#' A <- matrix(c(1, 0, 0, 1), nrow = 2)
#' check_stability(A)
#'
#' A <- matrix(c(1, 0, 0, 1, -0.1, 0, 0, -0.1), nrow = 2)
#' check_stability(A)
check_stability <- function(A) {
  AA <- big_A(A)
  # some noticeable numeric mistakes ... floating point arithmetic?
  all(abs(eigen(AA)$values) < 1)
  # could also check solutions of characteristic polynomial
  # efficiency?
}
###############################################################################.
#' Compute the mean of a VAR(p) process
#'
#' @inheritParams creat_varp_data
#' @param nu A `(K x 1)` matrix or vector. The intercept parameters of the
#'   process.
#'
#' @return A `(K x 1)` matrix containing the unconditional mean of each
#'   variable.
#'
#' @examples
#'  A <- matrix(c(.5, .4, .1, .5, 0, .25, 0, 0), nrow = 2)
#'  nu <- c(2, 3)
#'  mean_var_process(A, nu)
#'
#' \dontrun{
#'
#'  A <- matrix(c(1, 0, 0, 1), nrow = 2)
#'  nu <- c(2, 3)
#'  mean_var_process(A, nu)
#'  }
mean_var_process <- function(A, nu) {
  if(check_stability(A)) {
  K <- var_length(A)
  p <- lag_length(A)

  AA <- big_A(A)
  nunu <- big_nu(nu, p)

  selector(K, p) %*% solve(diag(K * p) - AA) %*% nunu
  } else {
    stop("process not stable")
  }
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
#' Retrieve number of variables
#'
#' Given a matrix of observations, residuals or slope coefficients, retrieve the
#' number of variables `K`.
#'
#' At the moment, this function is just an aptly / intentionally named wrapper of
#' `nrow(mat)`.
#'
#' @param mat `(K x M)` matrix, with `K` being the number of variables.
#' @return An integer scalar. The number of variables `K`.
#' @examples
#' p <- 3
#' K <- 4
#' var_length(matrix(0, K, K * p)) == K
var_length <- function(mat) {
  # if mat is a vector, it will return NULL but won't complain.
  if(!is.matrix(mat)) warning("Argument 'mat' is not a matrix.")
  nrow(mat)
}
###############################################################################.
#' Retrieve number of lags
#'
#' Given a matrix of slope parameters, retrieve the number of lags.
#'
#' At the moment, this function is just an intentionally named wrapper of
#' `ncol(A) / nrow(A)`.
#'
#' @inheritParams create_varp_data
#' @return An integer scalar. The lag length `p`.
#' @examples
#' p <- 3
#' K <- 4
#' A <- matrix(0, K, K * p)
#' lag_length(A) == p
#'
#' \dontrun{
#'
#' A <- matrix(0, 4, 7)
#' lag_length(A)
#' }
lag_length <- function(A) {
  p <- ncol(A) / nrow(A)
  stopifnot(p %% 1 == 0)
  p
}
###############################################################################.
#' Retrieve number of observations
#'
#' Given a matrix of observations or residuals, retrieve the number of sample
#' observations.
#'
#' At the moment, this function is just an intentionally named wrapper of
#' `ncol(mat)`.
#'
#' @param mat A `(K x N)` matrix, with `K` being the number of variables and `N`
#'   the sample size.
#' @return An integer scalar. The sample size `N`.
#' @examples
#' K <- 4
#' N <- 1E3
#'
#' obs_length(matrix(0, K, N)) == N
obs_length <- function(mat) {
  ncol(mat)
}
###############################################################################.
#' Vectorise a symmetric matrix
#'
#' Vectorise a symmetric matrix by stacking the elements on its lower
#' triangular matrix including the diagonal.
#'
#' @param mat An `(M x M)` matrix with arbitrary dimensions `M`.
#'
#' @return A column matrix of dimension `(M^2 + M)/2 x 1`.
#'
#' @examples
#' mat <- matrix(1:8, 4, 4)
#' mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
#' vech(mat)
vech <- function(mat) {
  stopifnot(all(mat == t(mat)))
  as.matrix(mat[lower.tri(mat, diag = TRUE)])
}
###############################################################################.
#' Vectorise a matrix
#'
#' Vectorise a matrix by stacking its elements in a vector.
#'
#' @param mat An `(M x N)` matrix with arbitrary dimensions `M` and `N`.
#'
#' @return A column matrix of dimension `NM x 1`.
#'
#' @examples
#' mat <- matrix(1:15, 3, 5)
#' vec(mat)
vec <- function(mat) {
  as.matrix(c(mat))
}
###############################################################################.
one_zeros <- function(r, c = 1) {
  stopifnot(r == c || c == 1)
  mat <- matrix(0, r, c)
  mat[1] <- 1
  mat
}
###############################################################################.
