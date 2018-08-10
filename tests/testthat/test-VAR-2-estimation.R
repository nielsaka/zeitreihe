context("Estimation routines for reduced-form VAR")

test_that("Simple multivariate OLS succeeds", {
  K <- 3
  N <- 1E3
  p <- 2

  Y <- do.call("create_varp_data", prep_input_varp(K, N, p))
  expect_equal_to_reference(ols_mv(Y = t(Y), p), "ols_mv.rds")

  skip_if(save_time)

  K <- 2
  N <- 1E6
  p <- 3

  Y <- do.call("create_varp_data", prep_input_varp(K, N, p))
  out <- ols_mv(Y = t(Y), p)
  expect_equivalent(out$BETA.hat[, -1], make_A(K, p), tol = 3E-3)
})

