###############################################################################.
#' Estimate Contemporaneous Structural Effects
#'
#' Estimate structural matrices `By` and `Be` by least squares or
#' maximum likelihood.
#'
#' The two matrices `By` and `Be` refer to the instantaneous effects matrices
#' in a vector autoregressive system \deqn{B_y y_t = A^*_1y_{t-1} + \dots +
#' A^*_py_{t-p} + B_e e_t.}{By * y_t = A_s(L) y_t + Be * e_t.} Thus, `By`
#' describes the contemporaneous effects between the elements of `y`. `Be`
#' describes the contemporaneous impact a structural error vector `e_t` has on
#' `y_t`.
#'
#' * `ols_cholesky` \cr Apply OLS and a Cholesky decomposition for recovering
#' `Be` from a simple recursive system. As usual, the ordering of the variables
#' will matter for a standard Cholesky decomposition. The effect matrix `By` is
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
#' Be <- matrix(0.4, K, K); Be[upper.tri(Be)] <- 0
#' Y0 <-matrix(0, nrow = K, ncol = p)
#' W <- matrix(rnorm(N * K), nrow = K, ncol = N)
#'
#' # draw data and estimate Be
#' Y <- create_svar_data(A, Be, Y0, W)
#' Be_hat <- ols_cholesky(Y, p)
#' Be_hat
ols_cholesky <- function(Y, p) {
  ols_fit <- ols_mv(Y, p)
  chol_decomp(ols_fit$SIGMA.hat)
}
###############################################################################.
#' Initialise the Concentrated Log-Likelihood
#'
#' Create the concentrated log-likelihood function of a structural VAR(p) for a
#' particular data set. Use it for estimating the contemporaneous structural
#' parameters `By` and `Be`.
#'
#' @inheritParams ols_mv
#' @param By
#' @param Be
#'
#' @return A function. It takes as input a named vector `args`. This vector
#'   consists of the structural parameters of the SVAR(p) model in vectorised
#'   form.
# TODO Note that names `???` are a requirement.
#' It will return the value of the log-likelihood at the specified parameter
#' vector. See the example for details.
#'
#' @examples
#' TODO factor out common code -> source modular R scripts
#' set.seed(8191)
#'
#' # number of variables, observations and lag length
#' K <- 3
#' N <- 1E6
#' p <- 2
#'
#' # prepare input
#' A <- cbind(matrix(0.1, K, K), matrix(-0.05, K, K)); diag(A) <- 0.4
#' Be <- matrix(0.4, K, K); Be[upper.tri(Be)] <- 0
#' Y0 <-matrix(0, nrow = K, ncol = p)
#' W <- matrix(rnorm(N * K), nrow = K, ncol = N)
#'
#' # create data and matrix Be
#' Y <- create_svar_data(A, Be, Y0, W)
#' By_init <- diag(K)
#' Be_init <- matrix(0, K, K)
#' Be_init[lower.tri(Be_init, diag = TRUE)] <- NA
#'
#' # initialise log-likelihood function and evaluate
#' log_lik <- conc_log_lik_init(Y, p, By_init, Be_init)
#' log_lik(rep(0.35, 6))
# TODO: use selection matrix C_B instead of specifying By and Be ?
# TODO check identification before creating return function
conc_log_lik_init <- function(Y, p, By, Be) {
  K <- var_length(Y)
  N <- obs_length(Y) - p

  ols_fit <- ols_mv(Y, p)
  # TODO take into account whether constant was computed;
  # should not be a fixed "1" below.
  # return function that accepts as params elements of By and Be
  SIGMA_hat <- ols_fit$SIGMA.hat * (N - K * p - 1) / N

  constant <- - K * N * log(2 * pi)

  function(args) {
    stopifnot(!any(is.na(args)))
    # TODO how many parameters in total?
    # stopifnot(length(args) == 2 * K^2)

    # TODO only optimise over non-restricted elements!!

    # TODO move to above
    free_elems <- is.na(By)
    g1 <- sum(is.na(By))
    g2 <- sum(is.na(Be))

    # nothin happens if g1 is zero
    By[is.na(By)] <- args[seq_len(g1)]
    Be[is.na(Be)] <- args[g1 + seq_len(g2)]

    Beinv <- solve(Be)
    S <- t(By) %*% t(Beinv) %*% Beinv %*% By
    # TODO log(det) using determinant()?
    .5 * (constant + N * (log(det(By)^2) - log(det(Be)^2) - sum(diag(S %*% SIGMA_hat))))
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
#'
#' set.seed(8191)
#'
#' K <- 3
#' N <- 1E6
#' p <- 2
#'
#' A <- cbind(matrix(0.1, K, K), matrix(-0.05, K, K)); diag(A) <- 0.4
#' Be <- matrix(0.4, K, K); Be[upper.tri(Be)] <- 0
#' Y0 <-matrix(0, nrow = K, ncol = p)
#'set.seed(8191)
#' W <- matrix(rnorm(N * K), nrow = K, ncol = N)
#'
#' Y <- create_svar_data(A, Be, Y0, W)
#'
#' By_init <- diag(K)
#' Be_init <- matrix(0, K, K)
#' Be_init[lower.tri(Be_init, diag = TRUE)] <- NA
#'
#'# sign is not unique !
#'# WORKS!
#'mle_fit_str <- mle_svar(Y, p, By_init, Be_init)
#'# won't converge; determinant of SIGMA not always positive?!
#' # mle_fit_red <- mle_var(Y, p)
# TODO try constraint optimisation again (for mle_var)
# TODO solve set of equations instead!!
# TODO combine with modularised function for mle? instead of mle_svar and mle_var?
mle_svar <- function(Y, p, By, Be) {

  log_lik <- conc_log_lik_init(Y, p, By, Be)
  neg_log_lik <- function(x) -1 * log_lik(x)

  # start values
  # args <- rnorm()
  args <- rep(0.1, sum(is.na(By), is.na(Be)))

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

    By <- diag(K)
    Be <- matrix(0, K, K)
    Be[lower.tri(Be, diag = TRUE)] <- NA
    Be[3, 3] = 1

    ## Likelihood ratio test

    # unrestricted
    # ols_fit <- ols_mv(Y, p)
    mle_fit_red <- mle_var(Y, p)

    # restricted
    mle_fit <- mle_svar(Y, p, By = By, Be = Be)

    g1 <- sum(is.na(By))
    g2 <- sum(is.na(Be))
    By[is.na(By)] <- mle_fit$par[seq_len(g1)]
    Be[is.na(Be)] <- mle_fit$par[g1 + seq_len(g2)]
    Be
    chol_decomp(ols_fit$SIGMA.hat * (N - (K * p) - 1) / N)

    Be - chol_decomp(ols_fit$SIGMA.hat * (N - (K * p) - 1) / N)

    B2 <- solve(By)
    (SIGMA_r <- B2 %*% Be %*% t(Be) %*% t(B2))
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
