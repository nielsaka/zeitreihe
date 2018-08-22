###############################################################################.
#' Estimate Contemporaneous Structural Effects
#'
#' Estimate structural matrices `B_Y` and `B_E` by least squares or
#' maximum likelihood.
#'
#' The two matrices `B_Y` and `B_E` refer to the instantaneous effects matrices
#' in a vector autoregressive system \deqn{B_Y Y_t = A^*_1Y_{t-1} + \dots +
#' A^*_pY_{t-p} + B_E E_t.}{B_Y * Y_t = A*(L) Y_t + B_E * E_t.} Thus, `B_Y`
#' describes the contemporaneous effects between the elements of `Y`. `B_E`
#' describes the contemporaneous impact a structural error vector `E_t` has on
#' the `Y_t`.
#'
#' * `ols_cholesky` \cr Apply OLS and a Cholesky decomposition for recovering
#' `B_E` from a simple recursive system. As usual, the ordering of the variables
#' will matter for a standard Cholesky decomposition. The effect matrix `B_Y` is
#' assumed to be an identity matrix.
#'
#' @inheritParams ols_mv
#' @return * `ols_cholesky` \cr A `(K x K)` numeric lower triangular matrix.
#'
#' @examples
#' set.seed(8191)
#'
#' # number of variables, observations and lag length
#' K <- 3
#' N <- 1E6
#' p <- 2
#'
#' # prepare input
#' A <- cbind(matrix(0.1, K, K), matrix(-0.05, K, K)); diag(A) <- 0.4
#' B <- matrix(0.4, K, K); B[upper.tri(B)] <- 0
#' Y0 <-matrix(0, nrow = K, ncol = p)
#' W <- matrix(rnorm(N * K), nrow = K, ncol = N)
#'
#' # draw data and estimate B_E
#' Y <- create_svar_data(A, B, Y0, W)
#' B_E_hat <- ols_cholesky(Y, p)
#' B_E_hat
ols_cholesky <- function(Y, p) {
  ols_fit <- ols_mv(Y, p)
  chol_decomp(ols_fit$SIGMA.hat)
}
###############################################################################.
#' Initialise the Concentrated Log-Likelihood
#'
#' Create the concentrated log-likelihood function of a structural VAR(p) for a
#' particular data set. Use it for estimating the contemporaneous structural
#' parameters `BY` and `BE`.
#'
#' @inheritParams ols_mv
#' @param BY
#' @param BE
#'
#' @return A function. It takes as input a named vector `args`. This vector
#'   consists of the structural parameters of the SVAR(p) model in vectorised
#'   form.
# TODO Note that names `???` are a requirement.
#' It will return the value of the log-likelihood at the specified parameter
#' vector. See the example for details.
#'
#' @examples
# TODO factor out common code -> source modular R scripts
#' set.seed(8191)
#'
#' # number of variables, observations and lag length
#' K <- 3
#' N <- 1E6
#' p <- 2
#'
#' # prepare input
#' A <- cbind(matrix(0.1, K, K), matrix(-0.05, K, K)); diag(A) <- 0.4
#' B <- matrix(0.4, K, K); B[upper.tri(B)] <- 0
#' Y0 <-matrix(0, nrow = K, ncol = p)
#' W <- matrix(rnorm(N * K), nrow = K, ncol = N)
#'
#' # create data and matrix BE
#' Y <- create_svar_data(A, B, Y0, W)
#' init_B <- diag(K)
#' init_B_inv <- matrix(0, K, K)
#' init_B_inv[lower.tri(init_B_inv, diag = TRUE)] <- NA
#'
#' # initialise log-likelihood function and evaluate
#' log_lik <- conc_log_lik_init(Y, p, init_B, init_B_inv)
#' log_lik(rep(0.35, 6))
# TODO: use selection matrix C_B instead of specifying B and B_inv ?
# TODO check identification before creating return function
conc_log_lik_init <- function(Y, p, B, B_inv) {
  K <- var_length(Y)
  N <- obs_length(Y) - p

  ols_fit <- ols_mv(Y, p)
  # TODO take into account whether constant was computed;
  # should not be a fixed "1" below.
  # return function that accepts as params elements of A and B
  SIGMA_hat <- ols_fit$SIGMA.hat * (N - K * p - 1) / N

  constant <- - K * N * log(2 * pi)

  function(args) {
    stopifnot(!any(is.na(args)))
    # TODO how many parameters in total?
    # stopifnot(length(args) == 2 * K^2)

    # TODO only optimise over non-restricted elements!!

    # TODO move to above
    free_elems <- is.na(B)
    g1 <- sum(is.na(B))
    g2 <- sum(is.na(B_inv))

    # nothin happens if g1 is zero
    B[is.na(B)] <- args[seq_len(g1)]
    B_inv[is.na(B_inv)] <- args[g1 + seq_len(g2)]

    # TODO here notation B and B_inv is really confusing
    # maybe just use B1 and B2 from the start?
    B2 <- solve(B_inv)
    S <- t(B) %*% t(B2) %*% B2 %*% B
    # TODO log(det) using determinant()?
    .5 * (constant + N * (log(det(B)^2) - log(det(B_inv)^2) - sum(diag(S %*% SIGMA_hat))))
  }
}
###############################################################################
#' @rdname ols_cholesky
#'
#' @details * `mle_svar` \cr ...
#'
#' @inheritParams conc_log_lik_init
#'
#' @return * `mle_svar` \cr ...
#'
#' @examples
#' K <- 3
#' N <- 1E6
#' p <- 2
#'
#' A <- cbind(matrix(0.1, K, K), matrix(-0.05, K, K)); diag(A) <- 0.4
#' B <- matrix(0.4, K, K); B[upper.tri(B)] <- 0
#' Y0 <-matrix(0, nrow = K, ncol = p)
#'set.seed(8191)
#' W <- matrix(rnorm(N * K), nrow = K, ncol = N)
#'
#' Y <- create_svar_data(A, B, Y0, W)
#'
#' B <- diag(K)
#' B_inv <- matrix(0, K, K)
#' B_inv[lower.tri(B_inv, diag = TRUE)] <- NA
#'
#'# sign is not unique !
#'# WORKS!
#'mle_fit_str <- mle_svar(Y, p, B, B_inv)
#'# won't converge; determinant of SIGMA not always positive?!
#' # mle_fit_red <- mle_var(Y, p)
# TODO try constraint optimisation again (for mle_var)
# TODO solve set of equations instead!!
# TODO combine with modularised function for mle? instead of mle_svar and mle_var?
mle_svar <- function(Y, p, B, B_inv) {

  log_lik <- conc_log_lik_init(Y, p, B, B_inv)
  neg_log_lik <- function(x) -1 * log_lik(x)

  # start values
  # args <- rnorm()
  args <- rep(0.1, sum(is.na(B), is.na(B_inv)))

  mle_fit <- optim(
    args, neg_log_lik,
    method = "L-BFGS-B",
    control = list(factr = 1e-10, maxit = 50*length(args)^2)
  )
  mle_fit
}

if (FALSE) {

  set.seed(8191)
  reps <- 100
  lr_test <- numeric(reps)
  freq <- 0
  for(i in 1:reps) {
    cat(i, "\n")

    K <- 3
    N <- 1E3
    p <- 2

    A <- cbind(matrix(0.2, K, K), matrix(-0.05, K, K)); diag(A) <- 0.4
    B_true <- matrix(0.2, K, K); B_true[upper.tri(B_true)] <- 0
    diag(B_true) <- 1
    Y0 <-matrix(0, nrow = K, ncol = p)
    W <- matrix(rnorm(N * K), nrow = K, ncol = N)

    Y <- create_svar_data(A, B_true, Y0, W)

    B <- diag(K)
    B_inv <- matrix(0, K, K)
    B_inv[lower.tri(B_inv, diag = TRUE)] <- NA
    B_inv[3, 3] = 1

    ## Likelihood ratio test

    # unrestricted
    # ols_fit <- ols_mv(Y, p)
    mle_fit_red <- mle_var(Y, p)

    # restricted
    mle_fit <- mle_svar(Y, p, B = B, B_inv = B_inv)

    g1 <- sum(is.na(B))
    g2 <- sum(is.na(B_inv))
    B[is.na(B)] <- mle_fit$par[seq_len(g1)]
    B_inv[is.na(B_inv)] <- mle_fit$par[g1 + seq_len(g2)]
    B_inv
    chol_decomp(ols_fit$SIGMA.hat * (N - (K * p) - 1) / N)

    B_inv - chol_decomp(ols_fit$SIGMA.hat * (N - (K * p) - 1) / N)

    B2 <- solve(B)
    (SIGMA_r <- B2 %*% B_inv %*% t(B_inv) %*% t(B2))
    # (SIGMA_ur <- ols_fit$SIGMA.hat * (N - (K * p) - 1) / N)
    (SIGMA_ur <- mle_fit_red$SIGMA.hat)
    # mle_fit_red$SIGMA.hat - ols_fit$SIGMA.hat
    # mle_fit_red$SIGMA.hat - ols_fit$SIGMA.hat * (N - (K * p) - 1) / N
    tcrossprod(B_true)

    SIGMA_r - SIGMA_ur
    det(SIGMA_r)
    det(SIGMA_ur)
    log(det(SIGMA_r)) - log(det(SIGMA_ur))

    z_r <- determinant(SIGMA_r)
    ldSIGMA_r <- as.numeric(z_r$sign * z_r$modulus)

    z_ur <- determinant(SIGMA_ur)
    ldSIGMA_ur <- as.numeric(z_ur$sign * z_ur$modulus)

    (lr_test[i] <- N * (ldSIGMA_r - ldSIGMA_ur))
    # TODO lr test varies around 13 (without identifying restrictions). Why??
    # numerical precision; bad optimisation? better with gradient provided?
    # should really be very close to 0
    freq <- freq + (lr_test > qchisq(0.95, 1))

  }
  lr_test
}

# TODO

# check if all estimates are unbiased (BETA.hat, SIGMA.hat) from OLS and ML
# try optim with gradient specified
# try with lower convergence criterion
# try alabama::constr.Optim.nl with det(SIGMA) > 0
# compare to vars package
# work out solution to set of equations (constrained recursive system)

# Observations

# error in conc_log_lik fixed -> wrong degrees of freedom
# without overidentifying restrictions, log diff of determinants is at 5E-6. This
# will lead to lr test statistis of N * 5E-6, so at N = 1E6 it is 5 !!
# is this due to numerical inprecision??
