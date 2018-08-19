###############################################################################.
#' Cholesky Decomposition
#'
#' Cholesky decomposition of a positive-definite matrix using the
#' Cholesky-Banachiewicz algorithm.
#'
#' @param A A matrix, must be positive definite. If not, an error will be
#'   thrown.
#'
#' @return A matrix, the 'square root' of matrix \code{A}, as in \eqn{A =
#'   LL^{T}}{A = LL'}. IS IT L or L'?
#' @export
chol_decomp <- function(A) {
  if (!isSymmetric(A)) stop("error: matrix not symmetric")
  L <- A
  L[, ] <- 0

  n <- dim(A)[1]
  for (i in seq_len(n)) {
    for (j in seq_len(i)) {
      zz <- A[i, j]
      for (k in seq_len(j - 1)) {
        zz <- zz - L[i, k] * L[j, k]
      }
      if (i == j) {
        if (is.nan(rr <- sqrt(zz))) stop("matrix not positive definite")
        L[i, j] <- rr
      }
      else {
        L[i, j] <- zz / L[j, j]
      }
    }
  }
  L
}
###############################################################################.
# Matrix Inversion of lower triangular matrix
# using simple sequential equation solving
invLT <- function(L) {
  X <- L; X[, ] <- 0
  n <- dim(L)[1]
  I <- diag(rep(1, n))

  for (row in seq_len(n)) {
    for (col in seq_len(row)) {
      zz <- I[row, col]
      for (k in seq_len(row - 1)) {
        zz <- zz - L[row, k] * X[k, col]
      }
      X[row, col] <- zz / L[row, row]
    }
  }
  X
}
###############################################################################.
# Matrix inversion of symmetric, pos.def. matrix
invSPD <- function(A) {
  L     <- chol_decomp(A)
  L.inv <- invLT(L)
  t(L.inv) %*% L.inv
}
###############################################################################.
#' Moving Average Coefficients
#'
#' @param A A matrix.
#' @param h An integer.
#'
#' @return
#' @export
MA_coeffs <- function(A, h){
  K <- nrow(A)
  p <- ncol(A) / K
  if (p > 1 & p %% 1 == 0) {
    A.large <- rbind(A[, 1:(K * (p - 1))], diag(K * (p - 1)))
    A.large <- cbind(A.large, rbind(A[, (K * (p - 1) + 1):(K * p)],
                                    matrix(0, K * (p - 1), K)))
    J <- cbind(diag(K), matrix(0, K, K * (p - 1) ))
  } else if (p == 1) {
    A.large <- A
    J       <- diag(K)
  } else {
    stop("lag length not well defined")
  }
  PHI <- array(, dim = c(K, K, h + 1))
  A.large.iter <- diag(K*p)
  for (i in seq_len(h + 1)) {
    PHI[, , i] <- J %*% A.large.iter %*% t(J)
    A.large.iter <- A.large.iter %*% A.large
  }
  PHI
}
###############################################################################.
#' structural moving average coefficients
#'
#' @param PHI
#' @param B
#'
#' @return
#' @export
#'
#' @examples
sMA_coeffs <- function(PHI, B) {
  h <- dim(PHI)[3]
  THETA <- PHI
  for (i in seq_len(h)) {
    THETA[, , i] <- PHI[, , i] %*% B
  }
  if (!is.null(dimnames(B))) {
    dimnames(THETA) <- c(dimnames(B), list(seq_len(h) - 1))
  }
  THETA
}
###############################################################################.
#' Multivariate ordinary least squares for VARs
#'
#' Perform multivariate OLS on a set of observables. Creation of lags for the
#' VAR will be handled by the function.
#'
#' This routine applies equation-wise OLS to a VAR. It is equivalent to running
#' OLS with `p` lags for each of the variables and collecting the coefficients
#' and residuals in a matrix.
#'
#' @param Y A `(K x N)` matrix carrying the data for estimation. There are
#' `N` observations for each of the `K` variables.
#' @param p An integer scalar, the lag length used for estimation.
#' @param const A boolean scalar, indicating wether a constant should be
#'   included. Defaults to `TRUE`.
#'
#' @return A list with four elements:
#'
#'   + `BETA.hat` is a `(K x [K * p + 1])` or `(K x [K * p])` named matrix
#'   containing the coefficient estimates. Its dimension depends on wether a
#'   constant is included. Each column contains the coefficents for a particular
#'   regressor, each row corresponds to a single equation. If `Y` did not
#'   carry names, the variables will be named *y1* and counting. The names of
#'   the lagged regressors are derived by appending the variable name with a
#'   *.l1* for lag one etc.
#'   If a constant is included, its coefficients are in the first column.
#'   + `SIGMA.hat` is a `(K x K)` named matrix. It is the covariance matrix of the
#'   residuals. The columns and rows are named after the variables in `Y`.
#'   + `U.hat` is a `(K x N)` matrix of residuals. Its rows are named after the
#'   variables, too.
#'   + `std.err` is a matrix of the same dimension and naming scheme as
#'    `BETA.hat`. It has the standard errors of the coefficient estimates.
#'
#'
#' @export
ols_mv <- function(Y, p, const = TRUE) {
  K <- nrow(Y)
  Kp <- K * p
  N <- ncol(Y)
  znames <- if (const) "const" else character(0)
  ynames <- if (!is.null(rownames(Y))) rownames(Y) else paste0("y", seq_len(K))
  znames <- c(znames, paste0(ynames, rep(paste0(".l", seq_len(p)), each = K)))
  cc <- if (const) 1 else numeric(0)
  if (p > 0) {
    Z <- rbind(cc, t(embed(t(Y), p))[, 1:(N - p)])
  } else {
    Z <- rep(cc, N - p)
  }
  Y <- Y[, -seq_len(p)]
  rownames(Z) <- znames
  rownames(Y) <- ynames
  ZZ.inv <- invSPD(Z %*% t(Z))
  BETA.hat  <- Y %*% t(Z) %*% ZZ.inv
  U.hat <- Y - BETA.hat %*% Z
  SIGMA.hat <- U.hat %*% t(U.hat) / (N - p - Kp - const)
  sigma.beta.hat <- ZZ.inv %x% SIGMA.hat
  sb.hat <- matrix(, K, Kp + const)
  for (i in seq_len(K)) {
    sb.hat[i, ] <- sqrt(diag(sigma.beta.hat))[seq(i, (Kp + const) * K, K)]
  }
  dimnames(sb.hat) <- dimnames(BETA.hat)

  list(BETA.hat = BETA.hat, SIGMA.hat = SIGMA.hat, U.hat = U.hat,
       std.err = sb.hat)
}
###############################################################################.
#' Hall's percentile interval bootstrap
#'
#' @param data
#' @param p
#' @param h
#' @param alpha
#' @param reps
#'
#' @return
#' @export
#'
#' @examples
sMA_CI <- function(data, p, h, alpha, reps) {
  K     <- ncol(data)
  N.all <- nrow(data)
  N.est <- N.all - p
  Y.pre    <- t(data)[, seq_len(p)]

  VAR.hat   <- ols_mv(data, p)
  U.hat     <- VAR.hat$U.hat
  BETA.hat  <- VAR.hat$BETA.hat
  B.hat     <- chol_decomp(VAR.hat$SIGMA.hat)
  PHI.hat   <- MA_coeffs(VAR.hat$BETA.hat[, -1], h)
  THETA.hat <- sMA_coeffs(PHI.hat, B.hat)

  U.hat.demean <- (U.hat - rowMeans(U.hat))
  Y.boot <- matrix(, K, N.all)

  THETA.boot <- array(, c(K, K, h + 1, reps))

  for (r in seq_len(reps)) {
    Y.boot[, seq_len(p)] <- Y.pre
    U.boot <- U.hat.demean[, sample(1:N.est, replace = TRUE)]
    if (is.matrix(Y.pre)) {
      Z.iter <- c(1, as.vector(Y.pre[, rev(seq_len(p))]))
    } else {
      Z.iter <- c(1, Y.pre)
    }
    for (i in seq_len(N.est)) {
      Y.boot[, p + i] <- BETA.hat %*% Z.iter + U.boot[, i]
      Z.iter <- c(1, as.vector(Y.boot[, rev(i + seq_len(p))]))
    }
    VAR.boot            <- ols_mv(t(Y.boot), p)
    B.boot              <- chol_decomp(VAR.boot$SIGMA.hat)
    PHI.boot            <- MA_coeffs(VAR.boot$BETA.hat[, -1], h)
    THETA.boot[, , , r] <- sMA_coeffs(PHI.boot, B.boot) - THETA.hat
  }

  IR.CI.upper <- THETA.hat - apply(THETA.boot, c(1, 2, 3), quantile,
                                   probs = alpha/2)
  IR.CI.lower <- THETA.hat - apply(THETA.boot, c(1, 2, 3), quantile,
                                   probs = 1 - alpha/2)
  if (is.null(dimnames(IR.CI.lower))) {
    dimnames(IR.CI.upper) <- c(dimnames(B.hat), list(seq_len(h) - 1))
    dimnames(IR.CI.lower) <- c(dimnames(B.hat), list(seq_len(h) - 1))
  }
  IR.CI <- rbind(cbind(as.data.frame.table(IR.CI.upper), bound = "upper"),
                 cbind(as.data.frame.table(IR.CI.lower), bound = "lower"))
  colnames(IR.CI) <- c("Variable", "Shock", "h", "value", "bound")
  IR.CI$h <- as.numeric(levels(IR.CI$h))[IR.CI$h]
  IR.CI
}
###############################################################################.
# Forecast error variance decomposition
FEVD <- function(THETA) {
  h <- dim(THETA)[3]
  K <- dim(THETA)[1]
  THETA_outer_sum   <- matrix(0, K, K)
  THETA_elem.sq_sum <- matrix(0, K, K)
  out <- THETA
  for (j in seq_len(h)) {
    THETA_outer_sum   <- THETA_outer_sum   + THETA[, , j] %*% t(THETA[, , j])
    THETA_elem.sq_sum <- THETA_elem.sq_sum + THETA[, , j] ^ 2
    out[, , j] <- THETA_elem.sq_sum / diag(THETA_outer_sum)
  }
  out.long <- as.data.frame.table(out)
  colnames(out.long) <- c("Variable", "Shock", "h", "value")
  out.long$h <- as.numeric(levels(out.long$h))[out.long$h]
  out.long
}
###############################################################################.
#' Maximum likelihood estimation of a VAR(p)
#'
#'
#'
#'
#' @examples
#' K <- 3
#' N <- 5E2
#' p <- 2
#'
#' set.seed(8191)
#'
#' A <- matrix(0.1, K, K * p)
#' Y0 <- matrix(0, K, p)
#' U <- matrix(rnorm(K * N), K, N)
#'
#' Y <- create_varp_data(A, Y0, U)
#'
#' mu <- rep(0, K)
#' A <- matrix(0.1, K, K * p)
#' SIGMA <- matrix(0, K, K)
#' diag(SIGMA) <- 1
#'
#' args <- c(mu = mu, a = vec(A), s = vech(SIGMA))
#'
#' ols_fit <- ols_mv(Y = Y, p = p, const = TRUE)
#'
#' mle_fit <- mle_var(Y, p)
#'
#' mle_fit$BETA.hat
#' ols_fit$BETA.hat
#'
#' mle_fit$SIGMA.hat
#' ols_fit$SIGMA.hat
#'
#' mle_fit$std.err
#' ols_fit$std.err
# TODO work out gradient and feed to optim? faster? not central right now
mle_var <- function(Y, p) {

  K <- var_length(Y)
  log_lik <- log_lik_init(Y, p)

  # start values; standard? random?
  mu <- rep(0, K)
  a  <- c(rep(c(1, rep(0, K)), K)[seq_len(K^2)], rep(0, K^2 * (p - 1)))
  s  <- vech(diag(K))

  args  <- c(mu = mu, a = a, s = s)
  lower <- c(rep(-Inf, length(mu) + length(a)),
             rep(c(0, rep(-Inf, K)), K - 1), 0)

  neg_log_lik <- function(args) -1 * log_lik(args)

  mle_fit <- optim(
    args, neg_log_lik,
    method = "L-BFGS-B", lower = lower, hessian = TRUE
  )

  get_elem <- function(x, start) x[check_start(x, start)]
  get_elem_par <- function(st) get_elem(mle_fit$par, st)

  mu <- get_elem_par("mu")
  a <- get_elem_par("a")
  s <- get_elem_par("s")

  BETA.hat <- matrix(c(mu, a), K, K * p + 1)
  SIGMA.hat <- matrix(duplication_matrix(K) %*% s, K, K)
  U.hat <- Y[, -(1:p)] - mu - BETA.hat %*% Y2Z(Y, p)

  # no need to switch sign since negative was minimised above.
  vari <- diag(solve(mle_fit$hessian))
  std.err <- matrix(
    sqrt(c(get_elem(vari, "mu"), get_elem(vari, "a"))),
    nrow = K, ncol = K * p + 1
    )

  list(BETA.hat = BETA.hat,
       SIGMA.hat= SIGMA.hat,
       U.hat = U.hat,
       std.err = std.err)
  # TODO names of matrices ...
  # TODO compare to results of OLS -> test file
}
###############################################################################.
#' Create a function to compute the log-likelihood
#'
#' Create the log-likelihood function of a VAR(p) for a particular data set.
#'
#' @inheritParams ols_mv
#' @inheritParams big_Y
#'
#' @return A function is returned. It takes as input a named vector `args`. This
#'   vector consists of the parameters of the VAR(p) model in vectorised form.
#'   Note that names `mu, a, s` are a requirement. It will return the value of
#'   the log-likelihood at the specified parameter vector. See the example for
#'   details.
#'
#' @examples
#' K <- 3
#' N <- 5E2
#' p <- 2
#'
#' set.seed(8191)
#'
#' A <- matrix(0.1, K, K * p)
#' Y0 <- matrix(0, K, p)
#' U <- matrix(rnorm(K * N), K, N)
#'
#' Y <- create_varp_data(A, Y0, U)
#'
#' log_lik <- log_lik_init(Y, p)
#'
#' mu <- rep(0, K)
#' SIGMA <- diag(K)
# TODO vec and vech are not exported!
#' args = c(mu = mu, a = vec(A), s = vech(SIGMA))
#'
#' log_lik(args)
log_lik_init <- function(Y, p) {
  K <- var_length(Y)
  Z <- Y2Z(Y, p, const = FALSE)
  Y <- Y[, -seq_len(p)] # TODO return from Y2Z ??
  N <- obs_length(Z)

  constant <- - K * N * log(2 * pi) / 2

  function(args) {
    stopifnot(!any(is.na(args)))
    stopifnot(length(args) == 1.5 * K + (p + .5) * K^2)

    # unconditional means
    mu <- args[seq_len(K)]
    stopifnot(check_start_all(mu, "mu"))

    # slope parameters
    a <- args[K + seq_len(K^2 * p)]
    stopifnot(check_start_all(a, "a"))
    A  <- matrix(a, K, K * p)

    # residual covariances
    s <- args[K + K^2 * p + seq_len((K^2 + K)/2)]
    stopifnot(check_start_all(s, "s"))
    SIGMA <- matrix(duplication_matrix(K) %*% s, K, K)

    # de-meaned regressor, residuals and crossproduct
    X <- Z - rep(mu, p)
    U <- Y - mu - A %*% X
    S <- tcrossprod(U)

    # log-likelihood
    constant - .5 * N * log(det(SIGMA)) - .5 * sum(diag(solve(SIGMA) %*% S))
  }
}
