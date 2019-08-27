###############################################################################.
#'Create a regressor matrix for a VAR(p) model
#'
#'Matrix `Y` is transformed into `Z` notation; see Lütkepohl (2005, p. 70).
#'This transformation facilitates estimation of VAR(p) models.
#'
#'@inheritParams ols_mv
#'
#'@return A `([K * p + 1] x N)` or `([K * p] x N)`matrix. The
#'  dimension depends on whether a constant is included. If it is included, the
#'  first row consists of ones. The next K rows will contain the data lagged by
#'  one time period. The remaining rows will contain further lags up to order p.
#'  The columns will contain time periods of which there are now `N` left
#'  after setting aside `p` pre-sample values for as lags.
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
#' @inheritParams big_nu
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

  selector_J(K, p) %*% solve(diag(K * p) - AA) %*% nunu
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
# TODO-9 better formula? frequency domain? in what ways "better"?
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
  selector_J(K, p)  %*% cf_GAMMA %*% t(selector_J(K, p))
}
###############################################################################.
#' Create Selection Matrix `J`
#'
#' Create the selection matrix `J` for converting objects in VAR(1) companion
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
#' J <- zeitreihe:::selector_J(K, p)
selector_J <- function(K, p){
  out <- matrix(0, K, K * p)
  diag(out) <- 1
  out
}
###############################################################################.
#' Convert a VAR(p) model to VAR(1) companion format
#'
#' Convert matrices containig observations, coefficients, residuals and
#' covariances to their VAR(1) companion format.
#'
#' Convert objects such as
#' * `A` - a matrix of slope parameters
#' * `nu` -  a column vector of intercepts
#' * `Y` - a matrix of observations
#' * `U` - a matrix of  residuals
#' * `SIGMA` - a covariance matrix
#'
#'  such that they correspond to objects taken from a
#' VAR(1) representation.
#'
#' @return * `companion_format` \cr A list with five elements. Four elements
#'   `Y`, `nu`, `A` and `U` correspond to the output of `big_Y`, `big_nu`,
#'   `big_A`, and `big_U`. The fifth element is the `(Kp x N)` regressor matrix
#'   `Z` populated with lags of `Y`.
#'
#' @examples
#' set.seed(8191)
#'
#' K <- 3
#' N <- 5
#' p <- 2
#'
#' nu <- matrix(1:K, ncol = 1)
#' A  <- matrix(0.1, K, K * p); diag(A) <- 1:K / 10
#' U  <- matrix(rnorm(K * N), K, N)
#' Y0 <- matrix(0, K, p)
#' Y  <- create_varp_data(A, Y0, U)
#'
#' cf <- companion_format(Y, nu, A, U)
#'
#' cf$U
#' cf$Y - cf$A %*% cf$Z
#'
#' \dontrun{
#' # input has to be in matrix form
#'  companion_format(Y[1, ], nu[1, ], A[1, ], U[1,])
#' }
companion_format <- function(Y, nu, A, U) {
  K <- var_length(A)
  p <- lag_length(A)

  stopifnot(K == var_length(A), K == var_length(U), K == var_length(nu))
  stopifnot(obs_length(Y) == obs_length(U) + p)

  list(
    Y   = big_Y(Y, p),
    A   = big_A(A),
    nu  = big_nu(nu, p),
    U   = big_U(U, p),
    Z   = Y2Z(Y, p, const = FALSE)
  )
}
###############################################################################.
#' @rdname companion_format
#' @inheritParams create_varp_data
#'
#' @return * `big_A` \cr A `(Kp x Kp)` matrix with the slope parameters on the
#'   first K rows. The remaining rows carry block-diagonal unit matrices and the
#'   remainder are zeros.
#' @examples
#'
#' big_A(A)
big_A <- function(A) {
  K <- var_length(A)
  p <- lag_length(A)

  XX <- matrix(0, K * p - K, K * p)
  diag(XX) <- 1
  rbind(A, XX)
}
###############################################################################.
#' @rdname companion_format
#' @inheritParams ols_mv
#'
#' @return * `big_Y` \cr A `(Kp x N)` matrix. The `p-1` lags of `Y` are
#'   pasted as rows below `Y`. The `p` pre-sample observations are deleted.
#'
#' @section Note: The difference between big_Y() and [Y2Z()] is that
#'   observations in matrix `Z` are shifted back by one time period in
#'   comparison to `Y`.
#'
#' @examples
#'
#' Y <- matrix(seq_len(K * N), K, N)
#' big_Y(Y, p)
big_Y <- function(Y, p) {
  Y2Z(cbind(Y, 1), p, const = FALSE)[, -1]
}
###############################################################################.
#' @rdname companion_format
#'
#' @inheritParams big_Y
#' @param nu A `(K x 1)` matrix of intercepts.
#'
#' @return * `big_nu` \cr A `(Kp x 1)` colum matrix. The first `K` elements
#'   correspond to `nu`. The remaining elements are zero.
#'
#' @examples
#'
#' big_nu(nu, p)
big_nu <- function(nu, p) {
  stopifnot(is.matrix(nu) && ncol(nu) == 1)
  expander_e(p) %x% nu
}
###############################################################################.
#' @rdname companion_format
#'
#' @param U A `(K x N)` matrix of residuals.
#' @inheritParams big_Y
#'
#' @return * `big_U` \cr A `(Kp x N)` matrix. The first `K` rows contain the
#'   residuals, the last `K(p-1)` rows consist of zeros.
#'
#' @examples
#'
#' U <- matrix(seq_len(K * N), K, N)
#' big_U(U, p)
big_U <- function(U, p) {
  expander_e(p) %x% U

  # more efficient? but also more verbose
  # K <- var_length(U)
  # N <- obs_length(U)
  # rbind(U, matrix(0, K * (p - 1), N))[, -(seq_len(p -1))]
}
###############################################################################.
#' @rdname companion_format
#'
#' @param SIGMA A `(K x K)` matrix of covariances. The covariance matrix of the
#'   residuals `U`.
#' @inheritParams big_Y
#'
#' @return * `big_SIGMA` \cr A `(Kp x Kp)` matrix. The upper left `(K x K)`
#'   block will contain the original covariance matrix. The remaining elements
#'   will be zero.
#'
#' @examples
#'
#' SIGMA <- matrix(0.5, K, K)
#' big_SIGMA(SIGMA, p = p)
big_SIGMA <- function(SIGMA, p) {
  expander_e(p, p) %x% SIGMA
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
#'
#' \dontrun{
#'
#' var_length(1:K)
#' }
var_length <- function(mat) {
  stopifnot(is.matrix(mat))
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
#' lag_length(A[1, ])
#' }
lag_length <- function(A) {
  stopifnot(is.matrix(A))
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
#'
#' \dontrun{
#'
#' obs_length(1:K)
#' }
obs_length <- function(mat) {
  stopifnot(is.matrix(mat))
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
#' @export
#' @examples
#' mat <- matrix(1:8, 4, 4)
#' mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
#' vech(mat)
vech <- function(mat) {
  stopifnot(isSymmetric(mat, tol = 1E-10))
  as.matrix(mat[lower.tri(mat, diag = TRUE)])
}
###############################################################################.
#' Vectorise a matrix
#'
#' Vectorise a matrix by stacking its columns in a vector.
#'
#' @param mat An `(M x N)` matrix with arbitrary dimensions `M` and `N`.
#'
#' @return A column matrix of dimension `NM x 1`.
#' @export
#' @examples
#' mat <- matrix(1:15, 3, 5)
#' vec(mat)
#' all.equal(mat, matrix(vec(man), 3, 5))
vec <- function(mat) {
  as.matrix(c(mat))
}
###############################################################################.
#' Expand a Matrix
#'
#' Create a vector or matrix for expanding the dimension of another matrix by
#' adding rows and columns filled with zeros. Multiply the result of
#' `expander_e` with a matrix using the Kronecker product. See examples.
#'
#' @param r An integer scalar. It specifies the number of rows to create.
#' @param c An integer scalar. It specifies the number of columns to create.
#'
#' `r` and `c` must be equal or `c` must be a one.
#'
#' @return An `(r x c)` matrix. The very first element in the
#' top left corner is a one, the rest are zeros.
#'
#' @examples
#'
#' expander_e(3)
#' expander_e(4, 4)
#'
#' expander_e(3, 3) %x% matrix(0.1, 4, 4)
#'
#' \dontrun{
#'
#' expander_e(4, 3)
#' }
expander_e <- function(r, c = 1) {
  stopifnot(r == c || c == 1)
  mat <- matrix(0, r, c)
  mat[1] <- 1
  mat
}
###############################################################################.
#' Column indices of duplication matrix
#'
#' Create a vector containing the column indices of a duplication matrix. It's
#' assumed the output will be matched with a row index running from 1, 2, ...
#' onwards. It's helpful for creating a duplication matrix for a base matrix of
#' dimension K x K.
#'
#' @inheritParams duplication_matrix
#'
#' @return An integer vector of length K^2.
#' @details `dm_col_index` is a helper function for [`duplication_matrix`].
dm_col_index <- function(K) {
  stopifnot(K %% 1 == 0)
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
### another version, more intuitive but slower
dm_col_index_2 <- function(K) {

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
#' Indices of a duplication matrix
#'
#' Create a matrix containing row and column indices of a duplication matrix.
#' These indices specify where the duplication matrix carries ones.
#'
#' It is a helper function for [`duplication_matrix`].
#'
#' @inheritParams duplication_matrix
#'
#' @return An integer matrix of dimension 2 x K^2.  The first columns
#'   corresponds to the row index, the second column to the column index.
dm_index <- function(K) {
  cbind(seq_len(K^2), dm_col_index(K))
}
###############################################################################.
#' Duplication Matrix `D`
#'
#' Create a duplication matrix for duplicating entries of a vector. This is
#' particularly useful for converting a vector `vech(M)` to a vector `vec(M)`,
#' where `M` is a symmetric matrix. See [vec()] and [vech()] and the example
#' below.
#'
#' @inheritParams selector_J
#'
#' @return A `(K^2 x [K^2 + K]/2)` selection matrix containing ones and zeros.
#'   When multiplied with a vector, it will select entries of that vector to be
#'   duplicated.
#'
#' @examples
#' K <- 6
#'
#' AA <- matrix(1:K^2, K, K)
#' AA[upper.tri(AA)] <- t(AA)[upper.tri(AA)]
#'
#' AA
#'
#' all(vec(AA) == duplication_matrix(K) %*% vech(AA))
duplication_matrix <- function(K) {
  res <- matrix(0, K*K, K*(K+1)/2)
  res[dm_index(K)] <- 1
  res
}
###############################################################################.
#' Title
#'
#' @param r
#' @param c
#'
#' @return
#' @export
#'
#' @examples
#' r <- 3
#' c <- 2
#'
#' A <- matrix(rnorm(r * c), r, c)
#'
#' K <- commutation_matrix(r, c)
#' all(K %*% vec(A) == vec(t(A)))
commutation_matrix <- function(r, c) {
  diag(r * c)[order(rep(seq_len(r), times = c)), ]
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
###############################################################################.
#' Sample Splitting
#'
#' Create a function that will perform sample splitting. Initialise it once
#' by calling `split_templ` and then re-use where applicable.
#'
#' The function `split_templ` is a template for creating another function. This
#' other function will accept a `(K x N)` matrix holding `N` observations of `K`
#' variables and will slice the columns of that matrix according to the numbers
#' passed to `split_templ`.
#'
#' @param Nburn An integer scalar. The number of **burn-in** observations in a
#'   simulation context. These are observations that are non-existent if the
#'   data are actual observation and not synthetically created for simulation
#'   purposes.
#' @param Npre An integer scalar. The number of **pre-sample** observations for
#'   lags etc.
#' @param Nest An integer scalar. The number of observations available for
#'   **estimation**.
#' @param Ntrain An integer scalar. The number of observations for **training**,
#'   e.g. hyperparameters, forecasts.
#' @param Neval An integer scalar. The number of observations for model
#'   **evaluation**.
#' @param Noos An integer scalar. The number of **true out-of-sample**
#'   observations for simulation purposes. Can be used for evaluating the final
#'   outcome (model choice etc.).
#'
#' @return A list with six elements. Each of which is a numeric matrix
#'   corresponding to one of the six arguments passed to `split_templ`. Every
#'   matrix has `K` rows and as many columns as specified in `split_templ`. If
#'   an argument of `split_sample` is zero (the default), then the matrix will
#'   have zero columns.
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

    K <- var_length(Y)
    slice <- function(from, to) matrix(Y[, from + seq_len(to)], nrow = K)
    list(
      burn_in       = slice(0, Nburn),
      pre_sample    = slice(Nburn, Npre),
      estimation    = slice(Nburn + Npre, Nest),
      training      = slice(Nburn + Npre + Nest, Ntrain),
      evaluation    = slice(Nburn + Npre + Nest + Ntrain, Neval),
      out_of_sample = slice(Nburn + Npre + Nest + Ntrain + Neval, Noos)
    )
  }
}
