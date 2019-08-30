context("Testing VAR reduced-form estimation")

test_that("Simple multivariate OLS succeeds", {

  # base reference

  K <- 3
  N <- 1E3
  p <- 2

  Y <- do.call("create_varp_data", prep_input_varp(K, N, p))
  expect_known_value(ols_mv(Y = Y, p), "ols_mv.rds")

  K <- 3
  N <- 1E2
  p <- 2

  Y <- do.call("create_varp_data", prep_input_varp(K, N, p))
  expect_equal(ncol(ols_mv(Y, p, const = FALSE)$BETA.hat), K * p)

  # input has no dimnames

  rownames(Y) <- NULL

  expect_identical(
    colnames(ols_mv(Y, p)$BETA.hat),
    c("const", "y1.l1", "y2.l1", "y3.l1", "y1.l2", "y2.l2", "y3.l2")
  )
  expect_identical(
    colnames(ols_mv(Y, p, const = FALSE)$BETA.hat),
    c("y1.l1", "y2.l1", "y3.l1", "y1.l2", "y2.l2", "y3.l2")
  )
  expect_identical(
    colnames(ols_mv(Y, p, const = FALSE)$SIGMA.hat),
    c("y1", "y2", "y3")
  )

  # skip_if(save_time)

  # test SIGMA.hat

  K <- 2
  N <- 1E2
  p <- 2
  reps <- 2E3

  set.seed(8191)
  input <- replicate(reps, prep_input_varp(K, N, p, seed = sample(1E5, 1)))
  Y <- lapply(seq_len(reps), function(x) do.call("create_varp_data", input[, x]))

  out <- lapply(Y, ols_mv, p = p, const = TRUE)
  (SIGMA.hat <- rowMeans(sapply(out, function(x) x$SIGMA.hat)))

  # slight bias bc of auto-regression?
  expect_equal(SIGMA.hat, c(1, 0, 0, 1), tol = 8E-2)


  # test BETA.hat

  K <- 2
  N <- 1E5
  p <- 3

  Y <- do.call("create_varp_data", prep_input_varp(K, N, p))
  out <- ols_mv(Y = Y, p)
  expect_equivalent(out$BETA.hat[, -1], make_A(K, p), tol = 1E-2)
  expect_equivalent(out$BETA.hat[,  1], make_nu(K), tol = 5E-1)
})

test_that("Simple ML estimaton of VAR succeeds", {
  K <- 3
  N <- 5E2
  p <- 2

  Y <- do.call("create_varp_data", prep_input_varp(K, N, p))

  init_mu <- rep(0, K)
  init_A  <- matrix(0.1, K, K * p)
  init_SIGMA <- diag(K)

  args <- c(mu = init_mu, a = init_A, s = vech(init_SIGMA))
  log_lik <- log_lik_init(Y, p)

  expect_equal(log_lik(args), -5934375.776099924)

  # fail without names
  args <- c(init_mu, a = init_A, s = vech(init_SIGMA))
  expect_error(log_lik(args), "check_start_all(mu, \"mu\") is not T", f = TRUE)
  args <- c(mu = init_mu, init_A, s = vech(init_SIGMA))
  expect_error(log_lik(args), "check_start_all(a, \"a\") is not T", f = TRUE)
  args <- c(mu = init_mu, a = init_A, vech(init_SIGMA))
  expect_error(log_lik(args), "check_start_all(s, \"s\") is not T", f = TRUE)
  args <- c(mu = init_mu, a = tail(c(init_A), -1), s = c(vech(init_SIGMA), 0))
  expect_error(log_lik(args), "check_start_all(a, \"a\") is not T", f = TRUE)

  # fail with NAs
  args <- c(mu = init_mu, a = init_A, s = c(tail(vech(init_SIGMA), -1), NA))
  expect_error(log_lik(args), "!any(is.na(args)) is not TRUE", f = TRUE)

  # fails with too many parameters
  args <- c(mu = init_mu, a = init_A, s = c(vech(init_SIGMA), 1))
  expect_error(log_lik(args), "length(args) == 1.5 * K + (p + 0.5) ", f = TRUE)

  #-----------
  # MLE vs OLS
  #-----------

  ols_fit <- ols_mv(Y = Y, p = p, const = TRUE)
  mle_fit <- mle_var(Y, p)

  expect_known_output(mle_fit, "mle_fit.txt", print = TRUE)

  mle_gradient_optim <- mle_var(Y, p, gradient = NULL)
  mle_gradient_analytic <- mle_var(Y, p, gradient = gradient_var_init(Y, p))

  # best solution? all correct? fastest?
  expect_equal(mle_fit, mle_gradient_optim, tolerance = 1E-6)
  expect_equal(mle_gradient_optim, mle_gradient_analytic, tolerance = 1E-4)

  expect_equivalent(mle_fit$BETA.hat, ols_fit$BETA.hat, tol = 1E-4)
  expect_equivalent(mle_fit$SIGMA.hat, ols_fit$SIGMA.hat * (N - K*p - 1) / N,
                    tol = 7E-5)
  # TODO: standard error for intercept is wrong in mle_fit
  # (currently still using standard error of mean mu instead of intercept nu)
  expect_equivalent(mle_fit$std.err, ols_fit$std.err, tol = 1E-3)
})

test_that("gradient is computed correctly", {
  set.seed(8192)

  N <- 1E4
  K <- 1
  p <- 1

  A  <- matrix(0.1, K, K * p); diag(A) <- 0.4
  Y0 <- matrix(0, K, p)
  U <- matrix(rnorm(N * K), K, N)

  Y <- create_varp_data(A = A, Y0 = Y0, U = U)

  gradient <- gradient_var_init(Y, p)

  args <- c(mu = rep(0, K), a = vec(A), s = vech(diag(K)))
  expect_known_output(gradient(args), "gradient.txt", print = TRUE)

  ols_fit <- ols_mv(Y, p, const = TRUE)
  ols_args <- c(mu = ols_fit$BETA.hat[, 1],
                a = ols_fit$BETA.hat[, -1],
                s = c(vech(ols_fit$SIGMA.hat)) * (N - 1 - K * p) / N)

  mle_fit <- mle_var(Y, p)
  mle_args <- c(mu = mle_fit$BETA.hat[, 1],
                a = mle_fit$BETA.hat[, -1],
                s = c(vech(mle_fit$SIGMA.hat)))

  expect_equal(mle_args, ols_args, tol = 1E-5)
  expect_equal(gradient(mle_args), gradient(ols_args), tol = 5E-2)

  ll <- function(...) -1 * log_lik_init(Y, p)(...)

  # TODO-2 add numDeriv to SuggestedPackages
  expect_equivalent(numDeriv::grad(ll, args), gradient(args), tol = 1E-8)
  expect_equivalent(numDeriv::grad(ll, ols_args), gradient(ols_args), tol = 1E-5)
  expect_equivalent(numDeriv::grad(ll, mle_args), gradient(mle_args), tol = 1E-5)
})

# TODO-6 comparison of estimators -> simulation with N = 1000 ..
# * ar.ols, ar.mle, ols_mv, mle_var, maxLik, vars, lm
if (FALSE) {
  # TODO-5 check
  # https://stats.stackexchange.com/questions/148722/estimation-with-mle-and-returning-the-score-gradient-qmle
  # https://www.stat.umn.edu/geyer/5931/mle/mle.pdf
  # https://www.mn.uio.no/math/tjenester/kunnskap/kompendier/num_opti_likelihoods.pdf
  ar_ols <- ar(c(Y),  aic = FALSE, order.max = p, method = "ols")
  ar_mle <- ar(c(Y),  aic = FALSE, order.max = p, method = "mle", demean = TRUE)

  ar_ols$ar
  ar_ols$x.intercept
  lm(Y[, -1] ~ Y[, -N])$coeff

  ols_fit$BETA.hat

  ar_mle$ar # closer comparison to arima0 ??
  mle_fit$BETA.hat

  ols_args
  mle_args
}
