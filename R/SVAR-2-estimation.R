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
#' @return The results may not necessarily be unique. The column sign, for
#'   example, could be indeterminate.
#'* `ols_cholesky` \cr A `(K x K)` numeric lower triangular matrix.
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
#' Compute Impulse Response Functions from a Cholesky Decomposition
#'
#' @param Y
#' @param h
#' @param p
#' @param DET
#' @param CI
#' @param label_shocks
#' @param norm An integer matrix of dimension `K x K`. Specifies how to
#'   normalise the sign of the impulse response functions. A non-zero `ij`
#'   element of `norm` specifies that the sign of the `i`th response to shock
#'   `j` shall be normalised. If `norm[i, j]` is negative, the sign will be
#'   negative; otherwise positive. Leave unset if no normalisation should be imposed.
#' @param cumulate An integer vector. An index vector which specifies the
#'   response variables to be cumulated over horizon `h`. Useful for
#'   variables in (percentage) changes. Leave unset if no cumulation is required.
#'
#' @return Data frame with columns `shock`, `response`, `h`, `point` and,
#' depending on input to `CI`, further columns with quantile estimates of the
#' impulse response distribution.
#' @export
#' @section TODO: consider more general deterministic terms (breaks, trends,
#'   seasonal dummies). Pass matrix of regressors via `DET` argument. Would have
#'   to adjust in ols_mv() as well and allow for more general deterministic
#'   components.
#'
#' @examples
#'
#' data(oil, package = "zeitreihe")
#'
#' Y <- t(oil)
#' p <- 24
#' h <- 15
#' CI <- cholesky_irfs_wild_bootstrap(
#'   reps = 10,
#'   quantiles = c(0.16, 0.84),
#'   stdev = 1:2
#' )
#'
#' label_shocks <- c(
#'   "oil supply shock",
#'   "aggregate demand shock",
#'   "oil specific-demand shock"
#' )
#'
#' # normalise shocks such that response of first variable is negative on impact
#' norm <- matrix(c(0, 0, 1, 0, 0, 1, 0, 0, 1), 3, 3)
#'
#' # cumulate the response of the first variable
#' cumulate <- 1
#'
#' IRF <- cholesky_irfs(Y, p, h, CI)
#' IRF[IRF$shock == "oil supply shock", ]$value <- - IRF[IRF$shock == "oil supply shock", ]$value
#'
cholesky_irfs <- function(
  Y,
  p,
  h,
  CI,
  label_shocks = NULL,
  norm = NULL,
  cumulate = NULL,
  DET = NULL
) {

  # check inputs
  if (!is.null(DET)) warning("Deterministics beside constant not implemented.")
  if (is.null(label_shocks)) label_shocks <- rownames(Y)

  model <- ols_mv(Y, p)

  IRFs_P <- cholesky_irfs_point(model, h, DET, label_shocks, norm, cumulate)
  IRFs_CI <- CI(model, h, DET, label_shocks, norm, cumulate)

  IRF <- as.data.frame.table(IRFs_P)
  colnames(IRF) <- c("response", "shock", "h", "value")
  IRF$h <- as.numeric(levels(IRF$h))[IRF$h]
  IRF$stat = "point"

  rbind(IRF, IRFs_CI)
}

###############################################################################.
#' Title
#'
#' @param model
#' @param h
#' @param DET
#' @param label
#' @param norm
#' @param cumulate
#'
#' @return
#' @export
#'
#' @examples
cholesky_irfs_point <- function(model, h, DET, label, norm, cumulate) {
  B <- chol_decomp(model$SIGMA.hat)

  # flip sign
  if (sum(norm) != 0) {
    if (max(colSums(!norm == 0)) > 1) {
      stop("At most one non-zero element per column allowed in 'norm'.")
    }
    indx <- which(norm != 0)
    flip <- ifelse(B[indx] * norm[indx] < 0, -1, 1)
    B <- B %*% diag(flip)
  }

  PHI <- MA_coeffs(A = model$BETA.hat[, -1], h = h)
  THETA <- sMA_coeffs(PHI = PHI, B = B)

  # cumulate changes
  if (!is.null(cumulate)) {
  THETA[cumulate, , ] <-
    aperm(apply(THETA[cumulate, , , drop = FALSE], c(1, 2), cumsum), c(2, 3, 1))
  }

  colnames(THETA) <- label
  THETA
}
###############################################################################.
#' Title
#'
#' @param reps
#' @param quantiles
#' @param stdev
#' @param draw_eta
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
cholesky_irfs_wild_bootstrap <-
  function(reps, quantiles = NULL, stdev = NULL, draw_eta = rnorm, seed = 1234){

  function(model, h, DET, label, norm, cumulate){
    set.seed(seed)

    A  <- model$BETA.hat[, -1]
    nu <- model$BETA.hat[, 1]
    Y  <- model$Y
    U  <- model$U

    K <- var_length(model$U)
    N <- obs_length(model$U)
    p <- lag_length(A)

    IRFb <- array(, dim = c(K, K, h + 1, reps))
    dimnames(IRFb) <- list(rownames(Y), label, seq_len(h+1) - 1, seq_len(reps))

    for (r in seq_len(reps)) {
      if (!r %% 100) cat("bootstrap:", r, "out of ", reps, "repetitions.\n")
      # draw residuals and block of pre-sample obs
      Ur <- U * t(draw_eta(N)) %x% matrix(1, K)
      Y0 <- Y[, sample.int(N+1, 1) + 1:p - 1]

      Yr <- create_varp_data(A = A, Y0 = Y0, U = Ur, nu = nu)
      modelr <- ols_mv(Yr, p)

      IRFb[, , , r] <- cholesky_irfs_point(modelr, h, DET, label, norm, cumulate)
    }

    # summarise
    if (!is.null(quantiles)) {
      quants <- function(x) apply(IRFb, c(1, 2, 3), quantile, probs = x)
      add_stat <- function(i) IRF_qn[[i]]$stat <<- paste0("q", quantiles[i])

      IRF_qn <- lapply(lapply(quantiles, quants), as.data.frame.table)
      lapply(seq_along(quantiles), add_stat)
      IRF_qn <- do.call(rbind, IRF_qn)
    } else {
      IRF_qn <- NULL
    }
    if (!is.null(stdev)) {
      add_stdev <- function(x) IRF_P + x * IRF_SD
      sub_stdev <- function(x) IRF_P - x * IRF_SD
      add_stat <- function(i) IRF_sdp[[i]]$stat <<- paste0("sdp", stdev[i])
      sub_stat <- function(i) IRF_sdm[[i]]$stat <<- paste0("sdm", stdev[i])

      IRF_P <- cholesky_irfs_point(model, h, DET, label, norm, cumulate)
      IRF_SD <- apply(IRFb, c(1, 2, 3), sd)

      IRF_sdp <- lapply(lapply(stdev, add_stdev), as.data.frame.table)
      IRF_sdm <- lapply(lapply(stdev, sub_stdev), as.data.frame.table)

      lapply(seq_along(stdev), add_stat)
      lapply(seq_along(stdev), sub_stat)

      IRF_sd <- rbind(do.call(rbind, IRF_sdp), do.call(rbind, IRF_sdm))
    } else {
      IRF_sd <- NULL
    }
    IRF_CI <- rbind(IRF_qn, IRF_sd)

    colnames(IRF_CI) <- c("response", "shock", "h", "value", "stat")
    IRF_CI$h <- as.numeric(levels(IRF_CI$h))[IRF_CI$h]
    IRF_CI
  }
}

###############################################################################.
#' Initialise the Concentrated Log-Likelihood
#'
#' Create the concentrated log-likelihood function of a structural VAR(p) for a
#' particular data set. Maximise it for estimating the contemporaneous
#' structural parameters `By` and `Be`.
#'
#' @inheritParams ols_mv
#' @param By A `(K x K)` matrix. Determines the contemporaneous relation between
#'   the endogeneous variables `y_t`. See examples and details below.
#' @param Be A `(K x K)` matrix. Determines the contemporaneous impact of a
#'   structural shock `e_t` on the endogeneous variables `y_t`. See examples and
#'   details below.
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
#' # create data and set restrictions
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
  # Using OLS estimate here as that is unbiased
  # Alternatively, could adjust by factor (N - K * p - 1) / N
  SIGMA_hat <- ols_fit$SIGMA.hat

  constant <- - K * N * log(2 * pi)

  # TODO force arguments?!

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
#' @details * `mle_svar` \cr Estimate `By` and `Be` by maximising the
#'   concentrated log-likelihood. The reduced-form parameters are first
#'   estimated by OLS while the remaining structural parameters are estimated by
#'   MLE in a second step. The function covers recursive and non-recursive,
#'   exactly identified and over-identified systems. The necessary restrictions
#'   on matrices `By` and `Be` are imposed by assigning numeric values to the
#'   restricted elements. Unrestricted elements are left to be estimated by
#'   setting them to `NA`.
#'
#' @inheritParams conc_log_lik_init
#'
#' @return * `mle_svar` \cr A list with two elements `By` and `Be`. Both are `(K
#'   x K)` numeric matrices. The elements which previously held `NA`s were
#'   replaced with estimates.
#'
#' @examples
#' set.seed(8191)
#'
#' K <- 3
#' N <- 1E6
#' p <- 2
#'
#' A <- cbind(matrix(0.1, K, K), matrix(-0.05, K, K)); diag(A) <- 0.4
#' Be <- matrix(0.4, K, K); Be[upper.tri(Be)] <- 0
#' Y0 <-matrix(0, nrow = K, ncol = p)
#' W <- matrix(rnorm(N * K), nrow = K, ncol = N)
#'
#' Y <- create_svar_data(A, Be, Y0, W)
#'
#' By_init <- diag(K)
#' Be_init <- matrix(0, K, K)
#' Be_init[lower.tri(Be_init, diag = TRUE)] <- NA
#'
#' mle_svar(Y, p, By_init, Be_init)
# TODO won't converge; determinant of SIGMA not always positive?! --> `mle_var(Y, p)`
# TODO try constraint optimisation again (for mle_var)
# TODO solve set of equations instead!!
# TODO combine with modularised function for mle? instead of mle_svar and mle_var?
mle_svar <- function(Y, p, By, Be) {

  log_lik <- conc_log_lik_init(Y, p, By, Be)
  neg_log_lik <- function(x) -1 * log_lik(x)

  # TODO start values
  # args <- rnorm()
  args <- rep(0.1, sum(is.na(By), is.na(Be)))

  mle_fit <- optim(
    args, neg_log_lik,
    method = "L-BFGS-B",
    control = list(factr = 1e-10, maxit = 50*length(args)^2)
  )

  By[is.na(By)] <- mle_fit$par[seq_len(sum(is.na(By)))]
  Be[is.na(Be)] <- mle_fit$par[sum(is.na(By)) + seq_len(sum(is.na(Be)))]

  list(By = By, Be = Be)
}

if (FALSE) {
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
###############################################################################
#' Title
#'
#' @param A
#' @param B
#' @param RA
#' @param RB
#'
#' @return
#' @export
#'
#' @examples
#'
#' K <- 2
#'
#' B <- matrix(c(1, 0, 0, 1), K, K)
#' RB <-matrix(c(1, rep(0, 6), 1), K^2, K)
#'
#' asy_cov_mat_struc_coeff(B = B, RB = RB)
#'
asy_cov_mat_struc_coeff <- function(A, B, RA, RB) {

  # just B
  if (missing(A)) {
    K <- var_length(B)

    breadslice <- diag(K) %x% solve(B)
    Kmat <- commutation_matrix(K, K)
    I <- t(RB) %*% t(breadslice) %*% (diag(K^2) + Kmat) %*% breadslice %*% RB

    return(solve(I))
  }

  ## TODO-7 just A
  if (missing(B)) stop()

  ## TODO-7 A and B
  stop()
}
