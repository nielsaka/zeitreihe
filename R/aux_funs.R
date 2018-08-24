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
#' matrix `A` is stable. That is, whether the eigenvalues of the associated
#' companion format are all smaller than one in absolute terms.
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
  stopifnot(check_stability(A))
  K <- var_length(A)
  p <- lag_length(A)

  AA <- big_A(A)
  nunu <- big_nu(nu, p)

  selector(K, p) %*% solve(diag(K * p) - AA) %*% nunu
}
###############################################################################.
#' Compute the Covariance Matrix of a VAR(p) Process
#'
#' Compute the unconditional covariance matrix of the observations
#' `y_t`.
#'
#' Computing the covariance matrix invloves an infinite sum. Computation is
#' stopped if the summed differences of the elements of two iterations of the
#' covariance matrix is less than the tolerance level `tol`.
#'
#' @inheritParams create_varp_data
#' @inheritParams big_SIGMA
#' @param h An integer scalar, the horizon at which to compute the
#'   covariances. Defaults to the contemporaneous covariances of `y_t`.
#' @param tol A numeric scalar, the tolerance level for stopping the
#'   computation. See Details.
#'
#' @return A `(K x K)` numeric matrix. It containes the covariances of `y_t` and
#'   `y_{t+h}`.
#'
#' @section Implementation: Plain brute force with no regard for
#'   efficiency.
#'
#' @examples
#' K <- 4
#' p <- 2
#'
#' A <- matrix(0.0, K, K * p)
#' SIGMA <- matrix(0.5, K, K)
#' cov_var_process(A, SIGMA)
#'
#' A <- matrix(-0.2, K, K * p); diag(A) <- 1:K / 10
#' cov_var_process(A, SIGMA)
#' cov_var_process(A, SIGMA, h = 5)
#' cov_var_process(A, SIGMA, h = 150)
# TODO better formula? frequency domain? in what ways "better"?
cov_var_process <- function(A, SIGMA, h = 0, tol = 1E-7) {
  stopifnot(check_stability(A))
  K <- var_length(A)
  p <- lag_length(A)

  cf_SIGMA <- big_SIGMA(SIGMA, p)
  cf_GAMMA <- cf_SIGMA
  cf_A <- big_A(A)

  cf_A_hi <- cf_A
  for (i in seq_len(h)) {
    cf_A_hi <- cf_A_hi %*% cf_A
  }

  crit <- 1
  cf_A_i <- cf_A
  repeat{
    cf_GAMMA_old <- cf_GAMMA
    cf_GAMMA <- cf_GAMMA + cf_A_hi %*% cf_SIGMA %*% t(cf_A_i)

    crit <- sum(abs(cf_GAMMA - cf_GAMMA_old))
    if (crit < tol) break

    cf_A_hi <- cf_A_hi %*% cf_A
    cf_A_i  <- cf_A_i  %*% cf_A
  }
  selector(K, p)  %*% cf_GAMMA %*% t(selector(K, p))
}
###############################################################################.
#' Create Selection Matrix `J`
#'
#' Create the selection matrix J for converting objects in VAR(1) companion
#' format to their original VAR(p) representation.
#'
#' @param K An integer scalar. The number of variables of the VAR(p) system.
#' @inheritParams ols_mv
#'
#' @return A `(K x Kp)` numerical matrix containing ones and zeros in selected
#'   places. The first `(K x K)` elements form an identity matrix. The rest are
#'   zeros.
#'
#' @examples
#' K <- 4
#' p <- 2
#'
#' J <- selector(K, p)
# TODO rename to `make_J` ?
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
## TODO fix arguments ?
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
#' Convert objects to VAR(1) companion format
#'
#' Convert objects such as a matrix of slope parameters `A` or intercepts `nu`,
#' a matrix of observations `Y` or residuals `U` or a covariance matrix `SIGMA`
#' such that they correspond to objects taken from a VAR(1) representation.
#'
#' @inheritParams create_varp_data
#'
#' @return * `big_A` \cr A `(Kp x Kp)` matrix with the slope parameters on the
#'   first K rows. The remaining rows carry block-diagonal unit matrices and the
#'   remainder are zeros.
#' @examples
#' K <- 4
#' N <- 7
#' p <- 2
#'
#' A <- matrix(0.1, K, K * p)
#' big_A(A)
big_A <- function(A) {
  K <- var_length(A)
  p <- lag_length(A)

  XX <- matrix(0, K * p - K, K * p)
  diag(XX) <- 1
  rbind(A, XX)
}
###############################################################################.
#' @rdname big_A
#'
#' @inheritParams ols_mv
#'
#' @return * `big_Y` \cr A `(Kp x N-p+1)` matrix. The `p-1` lags of `Y` are
#'   pasted as rows below `Y`. This leads to the loss of `p-1` sample
#'   observations.
#'
#' @section Note: difference between big_Y and Y2Z? matrix `Z` has `K*p + 1`
#'   rows whereas the output of `big_Y` has `K*p` rows. The very last lag is
#'   missing.
#'
#' @examples
#' K <- 3
#' N <- 5
#'
#' Y <- matrix(seq_len(K * N), K, N)
#' big_Y(Y, p)
big_Y <- function(Y, p) {
  Y2Z(cbind(Y, 1), p, const = FALSE)
}
###############################################################################.
#' @rdname big_A
#'
#' @inheritParams big_Y
#' @param nu A `(K x 1)` matrix of intercepts.
#'
#' @return * `big_nu` \cr A `(Kp x 1)` colum matrix. The first `K` elements
#'   correspond to `nu`. The remaining elements are zero.
#'
#' @examples
#'
#' nu <- as.matrix(seq_len(K))
#' big_nu(nu, p)
big_nu <- function(nu, p) {
  one_zeros(p) %x% nu

  # more efficient? but also more verbose
  # K <- var_length(nu)
  # rbind(nu, matrix(0, K * (p - 1), 1))
}
###############################################################################.
#' @rdname big_A
#'
#' @param U A `(K x N)` matrix of residuals.
#' @inheritParams big_Y
#'
#' @return * `big_U` \cr A `(Kp x N)` matrix. The first `K` rows contain the
#'   residuals, the last `K(p-1)` rows consist of zeros. TODO: does not really
#'   correspond to the output of `big_Y`. There, it's `N-p+1` columns; here,
#'   it's `N`.
#'
#' @examples
#'
#' U <- matrix(seq_len(K * N), K, N)
#' big_U(U, p)
big_U <- function(U, p) {
  one_zeros(p) %x% U

  # more efficient? but also more verbose
  # K <- var_length(U)
  # N <- obs_length(U)
  # rbind(U, matrix(0, K * (p - 1), N))[, -(seq_len(p -1))]
}
###############################################################################.
#' @rdname big_A
#'
#' @param SIGMA A `(K x K)` matrix of covariances. The covariance matrix of the
#'   residuals `U`.
#' @inheritParams big_Y
#'
#' @return * `big_SIGMA` \cr A `(Kp x Kp)` matrix. The upper left `(K x K)`
#'   block will contain the original covariance matrix: The remaining elements
#'   will be zero.
#'
#' @examples
#'
#' SIGMA <- matrix(0.5, K, K)
#' big_SIGMA(SIGMA, p = p)
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
#' Title
#'
#' @param K
#'
#' @return
#'
#' @examples
#' duplication_sequence(3)
duplication_sequence <- function(K) {
  indx <- numeric(K^2)
  l <- i <- 1
  for (col in seq_len(K)) {
    u <- 0
    for (row in seq_len(K)) {
      if (row >= col) {
        indx[i] <- l
        l = l + 1
      } else {
        indx[i] <- col + u
        u <- u + (K - row)
      }
      i = i + 1
    }
  }
  indx
}
###############################################################################.
#' Title
#'
#' @param K
#'
#' @return
#'
#' @examples
#'   K <- 6
#'   AA <- matrix(1:K^2, K, K)
#'   AA[upper.tri(AA)] <- t(AA)[upper.tri(AA)]
#'   AA
#'   recover_AA <- c(duplication_matrix(K) %*% vech(AA))
#'   all(c(AA) == recover_AA)
#'   all.equal(AA, matrix(recover_AA, K, K))
duplication_matrix <- function(K) {
  res <- matrix(0, K*K, K*(K+1)/2)
  res[cbind(seq_len(K^2), duplication_sequence(K))] <- 1
  res
}
###############################################################################.
# TODO compare performance and result to above
#' Title
#'
#' @param K
#'
#' @return
#'
#' @examples
# K <- 10
# all(duplication_sequence(K) == dupl_sequence_2(K))
dupl_sequence_2 <- function(K) {

  AA <- AA_sym <- matrix(seq_len(K^2), K, K)
  AA_sym[upper.tri(AA_sym)] <- t(AA)[upper.tri(AA_sym)]

  vec_AA_sym <- vec(AA_sym)
  vech_AA_sym <- vech(AA_sym)

  indx <- numeric(K^2)
  for (i in seq_len(K^2)) {
    indx[i] <- which(vec_AA_sym[i] == vech_AA_sym)
  }
  indx
}
###############################################################################.
#' Are all elements TRUE?
#'
#' Check whether all elements of a vector, matrix or array are TRUE.
#'
#' Note that numeric values are not coerced to logical. See examples.
#'
#' @param x A vector, matrix or array.
#' @return Boolean scalar. TRUE if every element of `x` is TRUE.
#' @examples
#' x <- array(TRUE, c(1, 2, 3))
#' all_true(x)
#'
#' x[1, 2, 3] <- FALSE
#' all_true(x)
#'
#' # Numeric is not coerced to logical
#' x[1, 2, 3] <- 1
#' all_true(x)
#'
#' # A list might fail
#' y <- list(a = TRUE, b = matrix(TRUE, 3, 3), c = list(TRUE, TRUE))
#' all_true(y)
all_true <- function(x) {
  res <- sapply(x, isTRUE)
  all(res) && length(res) > 0
}
###############################################################################.
#' Do all names start with "it"?
#'
#' Check whether the names of a vector all start with a certain string.
#'
#' @param x A vector. The names of the vector will be checked.
#' @param start A character string.
#'
#' @return * check_start \cr A boolean vector of the same length as `x` if `x`
#'   is named. Each element is either TRUE or FALSE. They are TRUE if the name
#'   of the element starts with `start`. If `x` did not carry names, a vector of
#'   length zero is returned.
#'
#' @examples
#' x <- c(tisch = 1, teller = 2, tafel = 3)
#' name_start(x, "t")
#'
#' x <- c(x, 4)
#' name_start(x, "t")
#'
#' x <- 1:4
#' name_start(x, "t")
check_start <- function(x, start) {
  grepl(pattern = paste0("^", start), x = names(x))
}
###############################################################################.
#' @rdname check_start
#'
#' @return * `check_start_all` \cr A boolean scalar. TRUE if the names of **all**
#' elements of `x` start with `start`. Otherwise FALSE.
#'
#' @examples
#' x <- c(tisch = 1, teller = 2, tafel = 3)
#' check_names(x, "t")
#'
#' x <- c(x, 4)
#' check_names(x, "t")
#'
#' x <- 1:4
#' check_names(x, "t")
check_start_all <- function(x, start) {
  all_true(check_start(x, start))
}
# TODO document, test(?)
# TODO memoise? another package dependency...
seq_mu <- function(K) if (K > 0) 1:K else stop("K is not positive")
seq_a  <- function(K, p) length(seq_mu(K)) + 1:(K^2 * p)
seq_s  <- function(K, p) length(c(seq_mu(K), seq_a(K, p))) + 1:(K^2 + K) / 2
###############################################################################.
#' Title
#'
#' @param Y
#' @param Nreal
#' @param Nsim
#'
#' @return
#' @export
#'
#' @examples
#' p <- 2
#' N <- 25
#'
#' Npre <- 2
#' Nest <- N
#' Noos <- 5
#'
#' Nsim <- Npre + Nest + Noos
#'
#' Y <- matrix(1:Nsim, 1, Nsim)
#' split_sample <- split_templ(Npre = Npre, Nest = Nest, Noos = Noos)
#' spl <- split_sample(Y)
#'
#' spl$est
#'
#' \dontrun{
#'split_sample_wrong <- split_templ(Npre = Npre, Nest = Nest)
#'split_sample_wrong(Y)
#' }
split_templ <- function(
  Nburn = 0, Npre = 0, Nest = 0, Ntrain = 0, Neval = 0, Noos = 0
) {
  function(Y) {
    stopifnot(Nburn + Npre + Nest + Ntrain + Neval + Noos == obs_length(Y))
    # TODO move matrix check to obs_length(Y), var_length(Y) etc?
    stopifnot(is.matrix(Y))

    K <- var_length(Y)
    slice <- function(from, to) matrix(Y[, from + seq_len(to)], nrow = K)

    list(
      burn  = slice(0, Nburn),
      pre   = slice(Nburn, Npre),
      est   = slice(Nburn + Npre, Nest),
      train = slice(Nburn + Npre + Nest, Ntrain),
      eval  = slice(Nburn + Npre + Nest + Ntrain, Neval),
      oos   = slice(Nburn + Npre + Nest + Ntrain + Neval, Noos)
    )
  }
}
