# first, just-identifying, recursive: ols_mv
# then decompose residual covmat using cholesky; easy.

###############################################################################.
#' Title
#'
#' @param Y
#' @param p
#'
#' @return
#'
#' @examples
#'
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
# Y <- create_svar_data(A, B, Y0, W)
# (B_hat <- ols_cholesky(Y, p))
#'
#'
ols_cholesky <- function(Y, p) {
  ols_fit <- ols_mv(Y, p)
  chol_decomp(ols_fit$SIGMA.hat)
}

# concentrated likelihood

#' @examples
#' #'
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
#' init_B <- diag(K)
#' init_B_inv <- matrix(0, K, K)
#' init_B_inv[lower.tri(init_B_inv, diag = TRUE)] <- NA
#'
#' log_lik <- conc_log_lik_init(Y, p, init_B, init_B_inv)
#'
#' log_lik(rep(0.35, 6))

# TODO: use selection matrix C_B instead of specifying B and B_inv ?
# TODO check identification before creating return function
conc_log_lik_init <- function(Y, p, B, B_inv) {
  K <- var_length(Y)
  N <- obs_length(Y) - p

  ols_fit <- ols_mv(Y, p)
  # return function that accepts as params elemnts of A and B
  SIGMA_hat <- ols_fit$SIGMA.hat * (N - K) / N

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
    .5 * (constant + N * (log(det(B)^2) - log(det(B_inv)^2) - sum(diag(S %*% SIGMA_hat))))
  }

}

# combine with modularised function for mle? instead of mle_svar and mle_var?

# second, non-recursive: maximum likelihood
# work out likelihood
# optimise over values of A, B  - concentrated likelihood
# write, document, test (example data? from book?)

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
#'# TODO try constraint optimisation again
#' # mle_fit_red <- mle_var(Y, p)
#'
#'
#'mle_fit_str$value
#'# mle_fit_red$value
#'
mle_svar <- function(Y, p, B, B_inv) {

  log_lik <- conc_log_lik_init(Y, p, B, B_inv)
  neg_log_lik <- function(x) -1 * log_lik(x)

  # start values
  # args <- rnorm()
  args <- rep(0.1, sum(is.na(B), is.na(B_inv)))

  mle_fit <- optim(
    args, neg_log_lik,
    method = "BFGS"
  )
  mle_fit
}

if (FALSE) {

  set.seed(8191)
  lr_test <- 0
  reps <- 100
  for(i in 1:reps) {

    K <- 3
    N <- 5E2
    p <- 2

    A <- cbind(matrix(0.1, K, K), matrix(-0.05, K, K)); diag(A) <- 0.4
    B_true <- matrix(200, K, K); B_true[upper.tri(B)] <- 0
    Y0 <-matrix(0, nrow = K, ncol = p)
    W <- matrix(rnorm(N * K), nrow = K, ncol = N)

    Y <- create_svar_data(A, B_true, Y0, W)

    B <- diag(K)
    B_inv <- matrix(0, K, K)
    B_inv[lower.tri(init_B_inv, diag = TRUE)] <- NA
    B_inv[3, 3] = 200


    ## Likelihood ratio test

    # unrestricted
    ols_fit <- ols_mv(Y, p)
    # mle_fit_red <- mle_var(Y, p)


    # restricted
    mle_fit <- mle_svar(Y, p, B = B, B_inv = B_inv)

    g1 <- sum(is.na(B))
    g2 <- sum(is.na(B_inv))
    B[is.na(B)] <- mle_fit$par[seq_len(g1)]
    B_inv[is.na(B_inv)] <- mle_fit$par[g1 + seq_len(g2)]
    B2 <- solve(B)
    (SIGMA_r <- B2 %*% B_inv %*% t(B_inv) %*% t(B2))
    (SIGMA_ur <- ols_fit$SIGMA.hat * (N - (K * p) - 1) / N)

    SIGMA_r - SIGMA_ur
    log(det(SIGMA_r)) - log(det(SIGMA_ur))

    (lr_test <- N * (log(det(SIGMA_r)) - log(det(SIGMA_ur))))
    # TODO lr test varies around 13 (without identifying restrictions). Why??
    # numerical precision; bad optimisation? better with gradient provided?
    # should really be very close to 0
    freq <- freq + (lr_test > qchisq(0.95, 1))

  }
  lr_test / reps
}

# TODO

# check if all estimates are unbiased (BETA.hat, SIGMA.hat) from OLS and ML
# try optim with gradient specified
# try with lower convergence criterion
# try alabama::constr.Optim.nl with det(SIGMA) > 0
