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

