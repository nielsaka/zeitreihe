context("Testing small auxiliary functions")

test_that("Transformation Y2Z() notation works", {

   # simple data matrix

  K <- 3
  N <- 6
  p <- 2

  Y <- matrix(seq_len(K*N), nrow = K)
  Z <- Y2Z(Y, p, const = FALSE)

  expect_equivalent(c(Y[, 2], Y[, 1]), Z[, 1])
  expect_equal(dim(Z), c(K * p, N - p))

  ###

  skip_if(save_time)

  # proper data for estimation

  K <- 2
  N <- 2E5
  p <- 3

  Y <- do.call("create_varp_data", prep_input_varp(K, N, p))
  Z <- Y2Z(Y, p, const = FALSE)
  Y <- Y[, -seq_len(p)]

  # estimate coefficients with OLS
  A <- Y %*% t(Z) %*% solve(Z %*% t(Z))

  # and compare to original values
  # deviation still quite large, but narrow down as n -> inf
  expect_equivalent(A, make_A(K, p), tolerance = 2e-3, scale = 1)
})


test_that("Stability criterion is check", {
  A <- matrix(c(.5, .4, .1, .5, 0, .25, 0, 0), nrow = 2)
  expect_true(check_stability(A))

  A <- matrix(c(1, 0, 0, 1), nrow = 2)
  expect_false(check_stability(A))

  A <- matrix(c(1, 0, 0, 1, -0.1, 0, 0, -0.1), nrow = 2)
  expect_true(check_stability(A))
})

test_that("Computing unconditional mean works", {
   A <- matrix(c(.5, .4, .1, .5, 0, .25, 0, 0), nrow = 2)
   nu <- c(2, 3)
   expect_equal(mean_var_process(A, nu), matrix(c(7.027027, 15.135135), 2, 1))

   A <- matrix(c(1, 0, 0, 1), nrow = 2)
   nu <- c(2, 3)
   expect_error(mean_var_process(A, nu), regexp = "process not stable")
})

