context("Estimation routines for reduced-form VAR")

test_that("Simple multivariate OLS succeeds", {
  K <- 3
  N <- 1E3
  p <- 2

  Y <- do.call("create_varp_data", prep_input_varp(K, N, p))
  expect_equal_to_reference(ols_mv(Y = Y, p), "ols_mv.rds")

  skip_if(save_time)

  # test SIGMA.hat
  K <- 2
  N <- 1E1
  p <- 2
  reps <- 2E3

  set.seed(8191)
  input <- replicate(reps, prep_input_varp(K, N, p, seed = sample(1E5, 1)))
  Y <- lapply(seq_len(reps), function(x) do.call("create_varp_data", input[, x]))
  out <- lapply(Y, ols_mv, p = p)
  SIGMA.hat <- rowMeans(sapply(out, function(x) x$SIGMA.hat))

  expect_equal(SIGMA.hat, c(1, 0, 0, 1), tol = 2E-2)

  # test BETA.hat
  K <- 2
  N <- 2E5
  p <- 3

  Y <- do.call("create_varp_data", prep_input_varp(K, N, p))
  out <- ols_mv(Y = Y, p)
  expect_equivalent(out$BETA.hat[, -1], make_A(K, p), tol = 4E-3)
})

