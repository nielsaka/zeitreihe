context("Estimation routines for reduced-form VAR")

test_that("Simple multivariate OLS succeeds", {

  # base reference

  K <- 3
  N <- 1E3
  p <- 2

  Y <- do.call("create_varp_data", prep_input_varp(K, N, p))
  expect_equal_to_reference(ols_mv(Y = Y, p), "ols_mv.rds")

  # no constant

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

  skip_if(save_time)

  # test SIGMA.hat

  K <- 2
  N <- 1E1
  p <- 2
  reps <- 2E3

  set.seed(8191)
  input <- replicate(reps, prep_input_varp(K, N, p, seed = sample(1E5, 1)))
  Y <- lapply(seq_len(reps), function(x) do.call("create_varp_data", input[, x]))

  out <- lapply(Y, ols_mv, p = p, const = FALSE)
  (SIGMA.hat <- rowMeans(sapply(out, function(x) x$SIGMA.hat)))
  expect_equal(SIGMA.hat, c(1, 0, 0, 1), tol = 4E-2)
  # slight bias bc of auto-regression?

  out <- lapply(Y, ols_mv, p = p, const = TRUE)
  (SIGMA.hat <- rowMeans(sapply(out, function(x) x$SIGMA.hat)))
  expect_equal(SIGMA.hat, c(1, 0, 0, 1), tol = 2E-2)

  # test BETA.hat

  K <- 2
  N <- 2E5
  p <- 3

  Y <- do.call("create_varp_data", prep_input_varp(K, N, p))
  out <- ols_mv(Y = Y, p)
  expect_equivalent(out$BETA.hat[, -1], make_A(K, p), tol = 4E-3)
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

  expect_equal(log_lik(args), -2250.7876)

  # fail without names
  args <- c(init_mu, a = init_A, s = vech(init_SIGMA))
  expect_error(log_lik(args), "check_names(mu, \"mu\") is not TRUE", f = TRUE)
  args <- c(mu = init_mu, init_A, s = vech(init_SIGMA))
  expect_error(log_lik(args), "check_names(a, \"a\") is not TRUE", f = TRUE)
  args <- c(mu = init_mu, a = init_A, vech(init_SIGMA))
  expect_error(log_lik(args), "check_names(s, \"s\") is not TRUE", f = TRUE)
  args <- c(mu = init_mu, a = tail(c(init_A), -1), s = c(vech(init_SIGMA), 0))
  expect_error(log_lik(args), "check_names(a, \"a\") is not TRUE", f = TRUE)

  # fail with NAs
  args <- c(mu = init_mu, a = init_A, s = c(tail(vech(init_SIGMA), -1), NA))
  expect_error(log_lik(args), "!any(is.na(args)) is not TRUE", f = TRUE)

  # fails with too many parameters
  args <- c(mu = init_mu, a = init_A, s = c(vech(init_SIGMA), 1))
  expect_error(log_lik(args), "length(args) == 1.5 * K + (p + 0.5) ", f = TRUE)
})

